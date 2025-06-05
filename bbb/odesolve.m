      subroutine uedriv

*     UEDRIV is the main driver routine for the two-dimensional edge
*     plasma code.  The code solves a system of fluid equations
*     that models the edge plasma in an axisymmetric configuration.
*     The numerical procedure used is the method of lines that consist
*     of the solution of a set of coupled ODEs for the fluid variables
*     for each grid point.

      implicit none

      Use(Share)    # cutlo
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp
      Use(Math_problem_size)   # neqmx
      Use(Timing)
      Use(UEpar)    # istep,iter,isdtsfscal
      Use(Lsode)    # mmaxu,dtmax,dtinit,maxpoly,yl,yldot
      Use(Solver_work_arrays)   # liw,lrp,iwork,rwork
      Use(Jac_work_arrays)      # lwp, liwp
      Use(Timary)   # nsteps,istep_nk,nsteps_nk
      Use(Compla)   # ni,up,vy,te,ti,phi,zeff,nil,upl,tel,til,ngl,phil
      Use(Grid)     # ijac,iyld,yldmax
      Use(Ident_vars) # exmain_evals
      Use(Ynorm)    # suscal,sfscal
      Use(Oldpla)
      Use(Decomp)   # ubw,lbw
      Use(Jacaux)   # yldot1,yldot0,issfon
      Use(Err_msg_out)   # errmsgflag
      Use(Opt_input)     # inopt,iworkin,rworkin
      Use(Constraints)   # icflag,icnstr,rlx,ylprevc
      Use(Time_dep_nwt)  # ylodt,dtreal,yloext,isyloext
      Use(Indexes)       # iseqalg
      Use(UEint)         # restart
      Use(Npes_mpi)      # npes,mype,ismpion
      Use(Parallv)       # nlocal, neqg,meth,itmeth,iatol,igs,iopt,ropt,
                         # rtol_pv,atol_pv,delt_pv
      Use(Interp)        # nis,ups,tes,tis,ngs,phis,nxold,nyold
      Use(Stat)
      Use(Jacreorder)    # ireorder
      Use(Jacobian)      # nnzmx
      Use(Flags)         # iprint

c Diagnostic data
      Use(Comgeo)        # gxf,sx
      Use(RZ_grid_info)  # rm,zm,bphi,bpol
      Use(Rhsides)       # seec
      Use(Comflo)        # feex, feey
      Use(Locflux)       # conxe, floxe
      Use(Conduc)        # hcxe
      Use(Indices_domain_dcl) # ivloc2sdgl
      Use(Xpoint_indices)
      Use(Indices_domain_dcg)
      Use(Cdv)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in


      integer ifake  #forces Forthon scripts to put implicit none above here

c_mpi      integer lenrpw,lenipw,nge,ier
c_mpi      integer ii,typeneq,neqt,ionecall
c_mpi      data typeneq/51/,ionecall/0/
c_mpi      integer*4 ii4

      external ffun, jacnw, psolnw, resid, jacd2, psold, jacd1
      external rhsvd, jacvd, psolvd, rhsnk, psetnk, psolnk, jacvnk
      real vnormnk,r1mach9

c     local variables
      real tbout, dtreal_sav, initguess(neq), snesans(neq), snesminusnksol
      real fnrm, fnew, tick, tock
      external tick, tock
      integer i,ifld,lid,ilg
      #Former Aux module variables
      integer ix,iy,igsp,iv


c **- For parallel mpi case, set up preliminary mpi stuff
      if (ismpion .eq. 1) call uedriv_pll

c ... Save initial time and set accumulated times to zero.
      tstart = tick()
      ttotfe = 0.
      ttotjf = 0.
      ttimpfe = 0.
      ttimpjf = 0.
      ttmatfac = 0.
      ttmatsol = 0.
      ttjstor = 0.
      ttjrnorm = 0.
      ttjreorder = 0.

c ... Set switch to time other packages if this one is being timed.
      call sapitim (istimingon)


*  -- initialize counters --
      istep = nsteps - 1
      istep_nk = 0   # not inside if test for switching from nksol to daspk
      iter = 0

*  -- initialize the system --
      if (ismpion.eq.0) then  # Serial version
        call ueinit
      elseif (ismpion.eq.1) then  # MPI parallel version
c_mpicvode        call fpvmalloc (neqg, ts, yl, meth, itmeth, iatol,
c_mpicvode     .                  rtol_pv,atol_pv,inopt,iopt,ropt,ier)
c_mpicvode        call fcvspgmr2 (jpre, igs, maxkd, delt_pv)
      endif

*  ---------------------------------------------------------------------
*  -- continue looping until istep=nsteps, then go to resetting parameters --
   10 continue
      if (istep .ge. nsteps .or. istep_nk .ge. nsteps_nk) goto 200

*    -- set old-time values -- but only if mesh size not changing
        if(nxold == nx .and. nyold == ny) then
          do ifld = 1, nisp
            call s2copy (nx+2, ny+2, nis(0:,0:,ifld), 1, nx+2,
     .            ni0(0:,0:,ifld), 1, nx+2)
            call s2copy (nx+2, ny+2, ups(0:,0:,ifld), 1, nx+2,
     .            up0(0:,0:,ifld), 1, nx+2)
            call s2copy (nx+2, ny+2, vy(0:,0:,ifld), 1, nx+2,
     .            vy0(0:,0:,ifld), 1, nx+2)
          enddo
          do igsp = 1, ngsp
            call s2copy (nx+2, ny+2, ngs(0:,0:,igsp), 1, nx+2,
     .            ng0(0:,0:,igsp), 1, nx+2)
          enddo
          call s2copy (nx+2, ny+2, tes, 1, nx+2, te0, 1,nx+2)
          call s2copy (nx+2, ny+2, tis, 1, nx+2, ti0, 1,nx+2)
          call s2copy (nx+2, ny+2, phis,1, nx+2, phi0,1,nx+2)
        endif  #loop checking nxold, nyold with nx, ny

c...  Set present yl variables to ylodt and ylprevc
         do 13 i = 1, neq
            ylodt(i) =  (1-isyloext)*yl(i)+isyloext*yloext(i)
            ylprevc(i)= (1-isyloext)*yl(i)+isyloext*yloext(i)
 13      continue

      call idalg  # set iseqalg() to i.d. algebraic eq.

c ... Call the solver.
         if (issfon .eq. 1) then
           if (icntnunk .eq. 0 .and. isdtsfscal.eq.0) then
             if (npes <= 1) then
               call sfsetnk (neq,yl,suscal,sfscal)
             endif
           endif
         else
            call sfill (neq, 1., sfscal(1), 1)
         endif
c ... Determine the solver option to use
            call set_dt(neq, yl, yldot)  # sets dtuse for time-step models
            if (isdtsfscal.eq.1) call sfsetnk (neq, yl, suscal, sfscal)
                                         # allow dt in calc of sfscal
            iopts = 1
c ... Load optional inputs for nksol.
c     iwork(3) & iwork(4) were set in sub. allocate; allocate prob., lwp, liwp
            iwork(1) = mmaxu
            iwork(2) = 0
            iwork(3) = lwp
            iwork(4) = liwp
            iwork(5) = iprint
            iwork(6) = 0
            iwork(7) = 1-errmsgflag
            iwork(8) = itermx
            iwork(9) = incpset
            rwork(1) = stepmx
            rwork(2) = del2nksol
            rwork(3) = taunksol
*-------------------------------------------------------------------------

            call nksol(neq,yl,yldot,rhsnk,jacvnk,suscal,sfscal,ftol,
     .                 stptol,rwork,
     .                 lrw,iwork,liw,iopts,iterm,psetnk,psolnk,mfnksol,
     .                 mdif,ipflag,icflag,icnstr,rlx,epscon1,epscon2,
     .                 icntnunk,adjf1)

            if (iterm .eq. 1) exmain_evals = exmain_evals + 1
*-------------------------------------------------------------------------
              nni(1) = iwork(10)
              nli(1) = iwork(11)
              nfe(1) = iwork(12)
              nje(1) = iwork(13)
              npe(1) = iwork(14)
              nps(1) = iwork(15)
              ncfl(1) = iwork(16)

      iter = iter + 1

c...  convert solver variables back to plasma variables
      call convsr_vo1 (-1, -1, yl)  # was one call to convsr
      call convsr_vo2 (-1, -1, yl)  # was one call to convsr
      call convsr_aux (-1, -1)

c...  If nksol is in time-dependent mode, increment istep and toutlsod
      if (dtreal .lt. 1.e5) then
         istep_nk = istep_nk + 1
         istep = min(istep_nk, nsteps) - 1
         if (istep .eq. 0) then
            toutlsod(1) = dtreal
         else
            toutlsod(istep+1) = toutlsod(istep) + dtreal
         endif
         if (nsteps_nk .gt. 1) write(*,*) 'Step number = ' ,istep+1,
     .                          '   Total time = ', toutlsod(istep+1)
      endif
c...  Store the variables at the output times
      do 29 ix = 0, nx+1
         do 28 iy = 0, ny+1
            do 281 ifld=1,nisp
               nist1(istep+1,ix,iy,ifld) = ni(ix,iy,ifld)
               upst1(istep+1,ix,iy,ifld) = up(ix,iy,ifld)
 281        continue
            test1(istep+1,ix,iy) = te(ix,iy)
            tist1(istep+1,ix,iy) = ti(ix,iy)
            do 282 igsp = 1, ngsp
               ngst1(istep+1,ix,iy,igsp) = ng(ix,iy,igsp)
 282        continue
            phist1(istep+1,ix,iy) = phi(ix,iy)
 28      continue
 29   continue

c...  Diagnostic for max time-rate-of-change
      dtreal_sav = dtreal
      dtreal = 1.e20  # only temporary to get correct yldot for diag.
      call ffun(neq,tout,yl,yldot0)
      dtreal = dtreal_sav
      yldnmx(istep+1)=abs(yldot0(1)/(yl(1)+cutlo))
      iyldnmx(istep+1) = 1
      do 30 iv = 2, neq
        if (iseqalg(iv).eq.0 .and. icnstr(iv).eq.1) then # omit B.C., up, phi
           if (abs(yldot0(iv)/(yl(iv)+cutlo)) .gt. yldnmx(istep+1)) then
              yldnmx(istep+1)=abs(yldot0(iv)/(yl(iv)+cutlo))
              iyldnmx(istep+1) = iv
           endif
        endif
 30   continue

c ... Average old and new values to damp oscillations if dtdamp > 0.
      if (dtdamp > 1.e-50 .and. nil(1,1,1) > 0.) then
        fnew = 1./(1. + (dtdamp/dtreal)**itdamp)
        do ifld = 1, nisp
          do iy = 0, ny+1
            do ix = 0, nx+1
              ni(ix,iy,ifld) = fnew*ni(ix,iy,ifld) + (1.-fnew)*
     .                                            nil(ix,iy,ifld)
              up(ix,iy,ifld) = fnew*up(ix,iy,ifld) + (1.-fnew)*
     .                                            upl(ix,iy,ifld)
            enddo
          enddo
        enddo

        do iy = 0, ny+1
          do ix = 0, nx+1
             te(ix,iy) = fnew*te(ix,iy) + (1.-fnew)*tel(ix,iy)
             ti(ix,iy) = fnew*ti(ix,iy) + (1.-fnew)*til(ix,iy)
             phi(ix,iy) = fnew*phi(ix,iy) + (1.-fnew)*phil(ix,iy)
          enddo
        enddo

        do igsp = 1, ngsp
          do iy = 0, ny+1
            do ix = 0, nx+1
              ng(ix,iy,igsp) = fnew*ng(ix,iy,igsp) + (1.-fnew)*
     .                                            ngl(ix,iy,igsp)
            enddo
          enddo
        enddo
      endif

*c.... Store plasma variables as "last" or "l" quantities for possible reuse

      if (istate .ge. 0) then
         do ifld = 1, nisp
            call s2copy (nx+2, ny+2, ni(0:,0:,ifld), 1, nx+2,
     .                   nil(0:,0:,ifld), 1, nx+2)
            call s2copy (nx+2, ny+2, up(0:,0:,ifld), 1, nx+2,
     .                   upl(0:,0:,ifld), 1, nx+2)
         enddo

            call s2copy (nx+2, ny+2, te, 1, nx+2, tel, 1, nx+2)
            call s2copy (nx+2, ny+2, ti, 1, nx+2, til, 1, nx+2)
            call s2copy (nx+2, ny+2, phi, 1, nx+2, phil, 1, nx+2)

         do igsp = 1, ngsp
            call s2copy (nx+2, ny+2, ng(0:,0:,igsp), 1, nx+2,
     .                   ngl(0:,0:,igsp), 1, nx+2)
         enddo
      endif

*  ---------------------------------------------------------------------
c ... Save last dtreal timestep
      dtreal_old = dtreal

c **- Diagnostic output for parallel mpi version
*     -- increment and loop --
         istep = istep + 1
         go to 10


*  -- looping istep to nsteps is done, reset some parameters --
  200 continue

      ts = 0.
      istate = 1
      info(1) = 0
      iopts = inopt
      if (iopts .eq. 1) then
         rwork(5) = hu(istep)/tadj
         rwork(6) = 0.
         rwork(7) = 0.
         rwork(8) = 0.
         iwork(5) = 0
         iwork(6) = 0
         iwork(7) = 0
         iwork(8) = 0
         iwork(9) = 0
      endif
      tend = tock(tstart)
      if (iprinttim .eq. 1) call wtottim  # write out timing information

      return
      end
c****** end of subroutine uedriv *********************
c ------------------------------------------------------------------------
      subroutine uedriv_pll

c **- Initializes integration/solver routines for mpi parallel version

      implicit none

      Use(Dim)
      Use(Math_problem_size)
      Use(Constraints)
      Use(UEint)
      Use(UEpar)
      Use(Lsode)
      Use(Npes_mpi)
      Use(Parallv)
      Use(Flags)
C diagnostic data
      Use(Indices_domain_dcl)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in


      integer ifake  #forces Forthon scripts to put implicit none above here

CC c_mpi      include 'mpif.h'
c_mpi      integer status(MPI_STATUS_SIZE)

c     local variables
      real tbout, dtreal_sav
      integer i,ifld,lid,ier,ierr,mu,ml
      integer ii,typeneq,neqt,ionecall
      data typeneq/51/,ionecall/0/,ier/0/

c =======================================================================
c  -- initialize the system --
      restart = 1
      call ueinit

      NLOCAL = neq

      IGS = 0
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
  20    FORMAT(///' FPVINITMPI returned IER =',I5)
        STOP
      ENDIF

      MU = numvar*(nx+4)
      ML = numvar*(nx+4)

          do i=1,40
           iopt(i) = 0
           ropt(i) = 0.0
          enddo
c ==================================================================

      if (ionecall .eq. 1) then
         call pandf (-1, -1, neq, 0., yl, yldot)
 3331    continue
      return
      endif

      return
      end
c****** End of subroutine uedriv_pll ************************************
c-----------------------------------------------------------------------
      subroutine wtottim

c ... Writes out total timing data

      implicit none

      Use(Dim)        # nisp, nhsp
      Use(Timing)     # tend,tstart,ttotfe,...

         write(*,*) ' '
         write(*,900) 'Total time for last solution = ', tend
         write(*,901) 'Total full f evaluation = ', ttotfe
         write(*,902) 'Impur. part of full f evaluation = ', ttimpfe
         write(*,901) 'Total CPU Jacobian f evaluation = ', ttotjf
         write(*,901) 'Total Wall-Clock Jacobian f eval. = ', ttjstor
         write(*,902) 'Impur. part of Jacobian eval. = ', ttimpjf
         if (nisp .gt. nhsp) call wapitim
         write(*,901) 'Total Matrix factorization = ', ttmatfac
         write(*,901) 'Total Matrix backsolve = ', ttmatsol
         write(*,901) 'Total row normalization = ', ttjrnorm
         write(*,901) 'Total row and column reordering = ', ttjreorder
c         write(*,901) 'Total in other Jacobian work = ', ttjstor-ttotjf
 900     format(a36,20x,f10.4,' sec')
 901     format(a36,10x,f10.4,10x,' sec')
 902     format(a36,f10.4,20x,' sec')

         call wspltim

      return
      end

