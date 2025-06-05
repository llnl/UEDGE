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
      Use(UEpar)    # istep,iter,svrpkg,isdtsfscal
      Use(Lsode)    # mmaxu,dtmax,dtinit,maxpoly,yl,yldot
      Use(Solver_work_arrays)   # liw,lrp,iwork,rwork
      Use(Jac_work_arrays)      # lwp, liwp
      Use(Timary)   # nsteps,istep_nk,nsteps_nk
      Use(Compla)   # ni,up,vy,te,ti,phi,zeff,nil,upl,tel,til,ngl,phil
      Use(Grid)     # ig,imeth,ijac,iyld,yldmax
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


cpetsc      use PETSc_Snes_Param
cunipetsc      integer*4 uedgeComm,ii4
cunipetsc      data uedgeComm/0/

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
c...  imeth = inewton(igrid) is set in allocate & selects between
c...  time-dependent solvers(imeth=0) and Newton solvers(imeth=1)

      if (imeth .eq. 0) then
         istep = 0
         tbout = 1.
         if (nsteps .gt. 1) tbout = exp (log(trange)/(nsteps-1))
      elseif (imeth .eq. 1) then
         istep = nsteps - 1
      endif
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

      if (svrpkg.ne.'cvode') call idalg  # set iseqalg() to i.d. algebraic eq.

c ... Call one of several possible solvers.
c ... imeth=inewton(igrid) is a flag for time-dependent(=0) or Newton(=1)

      if (imeth .eq. 0) then    # Begin large imeth if test (search for imeth)

            tout = (tbout)**istep * runtim
            toutlsod(istep+1) = tout

c ... Determine the solver option to use
       if (svrpkg .eq. 'daspk') then
         ipar(1) = neq
         ipar(2) = lbw
         ipar(3) = ubw
         lid = 40
         if (info(12) .ne. 0) then  # moved from bottom of subrou. allocate
            iwork(27) = lwp
            iwork(28) = liwp
         endif
         if (info(10).eq.1 .or. info(10).eq.3) then  # use constraint checking
            lid = 40 + neq
            do i = 1, neq
               iwork(40+i) = 2*icnstr(i)
            enddo
         endif
         if (info(7).eq.1) rwork(2) = dtmax
         if (info(8).eq.1 .and. istep.eq.0) rwork(3) = dtinit
         if (info(9).eq.1) iwork(3) = maxpoly
         if (info(11).eq.1) then
            do i = 1, neq
               iwork(lid+i) = 1-2*iseqalg(i)   # iwork=-1 for alg.;=1 for diff.
            enddo
         endif
         if (info(12).eq.0) then   # direct method (banded storage only)
            iwork(1)=lbw
            iwork(2)=ubw
            if (info(13) .eq. 1) rwork(16) = epnldpk
            call sdaspk (resid,neq,ts,yl,yldot,tout,info,rtol,atol,idid,
     .                 rwork,lrw,iwork,liw,rpar,ipar,jacd1,psold,srtolpk)

         else                      # Krylov method
c        Note that iwork(27) & iwork(28) were set in subroutine allocate.
            if (info(13) .eq. 1) then
              rwork(10) = eplidpk
              rwork(16) = epnldpk
            endif
            call sdaspk (resid,neq,ts,yl,yldot,tout,info,rtol,atol,idid,
     .                 rwork,lrw,iwork,liw,rpar,ipar,jacd2,psold,srtolpk)
         endif
         if (idid.lt.0) then
            write(*,*) 'idid = ',idid
            call xerrab("")
         endif
         if (info(14).eq.1) goto 200   #init. only, jump out of time loop

       elseif (svrpkg .eq. 'vodpk') then
         ipar(1) = lbw
         ipar(2) = ubw
         iopts = inopt
         mf = 21
c        iwork(1) & iwork(2) were prev. set in sub. allocate; allocate problem
         iwork(1) = lwp
         iwork(2) = liwp
         iwork(3) = jpre
         iwork(4) = jacflg
c ... optional input for iopts=1 is stored in non-eraseable arrays
c ... iworkin(1:10) and rworkin(1:10)
         do i = 5, 10
            if (iopts .eq. 1) then
               iwork(i) = iworkin(i)
               rwork(i) = rworkin(i)
            else
               iwork(i) = 0
               rwork(i) = 0.
            endif
         enddo
         call vodpk (rhsvd, neq, yl, ts, tout, itol, rtol, atol, itask,
     .             istate, iopts, rwork, lrw, iwork, liw, jacvd,
     .             psolvd, mf, rpar, ipar, srtolpk, efacn)

       elseif (svrpkg .eq. 'cvode') then
ccc          call writeToLog('***  cvode option is not turned on ****')
ccc          call xerrab("")
c_mpicvode          CALL FCVODE (tout, ts, yl, itask, istate)

      else
         write(*,*)
         write(*,*) 'Invalid solver option, check variable svrpkg'
         write(*,*) 'Allowed values are: '
         write(*,*) '      vodpk'
         write(*,*) '      cvode'
         write(*,*) '      daspk'
         write(*,*) '      nksol'
         write(*,*) '      newton'
         write(*,*)
         call xerrab("")
       endif

       ts = tout

      elseif(imeth .eq. 1) then    # Second option for imeth if test
         if (issfon .eq. 1) then
           if (icntnunk .eq. 0 .and. isdtsfscal.eq.0) then
             if (npes <= 1) then
               call sfsetnk (neq,yl,suscal,sfscal)
             endif
           endif
         else
            call sfill (neq, 1., sfscal(1), 1)
         endif
         if((svrpkg .eq. 'nksol') .or.(svrpkg .eq. 'petsc') ) then  #not above issfon because newton poss
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

cSEK        snesdebug=1  ! TMP
*-------------------------------------------------------------------------
cpetsc      if (snesdebug. eq. 1) then
cpetsc*       Move the data to a dummy vector, yl->initguess, to run nksol on the same data
cpetsc        do i=1,neq
cpetsc          initguess(i)=yl(i)
cpetsc        end do
cpetsc        call rhsnk(neq,yl,yldot)
cpetsc        fnrm = vnormnk(neq,yldot,sfscal)
cpetsc        write(*,*)"fnrm =",fnrm
cpetsc        write(*,*)"kspatol =",(r1mach9(4)+epscon1*min(epscon2,fnrm))*fnrm
cpetsc        write(*,*)"restart =",mmaxu
cpetsc      endif

cpetsc*     Set up the stopping criteria for the nonlinear/linear solver
cpetsc*     All these values are present in the module PETSc_Snes_Param
cpetsc*      kspatol     = (r1mach9(4)+epscon1*min(epscon2,fnrm))*fnrm
cpetsc*      ksprestart = mmaxu    # from module lsode
cpetsc*      snesits    = itermx   # from module lsode

cpetsc      if ((svrpkg .eq. 'petsc') .or. (npes .gt. 1)) then
cpetsc*       Now solve the system by PETScSNES
cpetsc      flush(6)
cdb_solve      call MPI_Barrier(uedgeComm,ier)
Ccdb_solve      write(6,*)"[",mype,"] call PETSCSNES()..."
Ccdb_solve      flush(6)
cpetsc        call PETSCSNES(uedgeComm,neq,neqg,yl,yldot,iscolnorm,suscal,sfscal,nnzmx,ireorder,incpset,icflag,rlx,
cpetsc     &    int(icnstr,kind(ii4)),dtreal,iterm)
cpetsc      flush(6)

cpetsc        if (snesdebug. eq. 1) then
cpetsc          do i=1,neq
cpetsc            snesans(i)=yl(i)
cpetsc          end do
cpetsc          do i=1,neq
cpetsc            yl(i)=initguess(i)
cpetsc          end do

cpetsc*         This will let us do some nksol output I hope
cpetsc*         These variables are stored in PETSc_Snes_Param
cpetsc          nksolfnrm(1)=fnrm
cpetsc          nksoltime(1)=gettime(sec4)
cpetsc          nksollinits(1)=0
cpetsc          nksolfeval(1)=0
cpetsc        endif
cpetsc      endif

cpetsc      if ((svrpkg.eq.'nksol' .and. npes.eq.1) .or. (snesdebug. eq. 1)) then
*-------------------------------------------------------------------------

            call nksol(neq,yl,yldot,rhsnk,jacvnk,suscal,sfscal,ftol,
     .                 stptol,rwork,
     .                 lrw,iwork,liw,iopts,iterm,psetnk,psolnk,mfnksol,
     .                 mdif,ipflag,icflag,icnstr,rlx,epscon1,epscon2,
     .                 icntnunk,adjf1)

            if (iterm .eq. 1) exmain_evals = exmain_evals + 1
*-------------------------------------------------------------------------
cpetsc      endif
cpetsc      if (snesdebug. eq. 1) then
cpetsc*       Compare SNES solution and nksol solution
cpetsc        snesminusnksol=0
cpetsc        do i=1,neq
cpetsc          snesminusnksol=snesminusnksol+(yl(i)-snesans(i))*(yl(i)-snesans(i))
cpetsc        end do
cpetsc        write(*,*)"2-norm of difference SNES-NKSOL: ",snesminusnksol
cpetsc*       Hopefully this will be the output of nksol
cpetsc*       Ive hardwired this to only output the first 6 N-K iterations
cpetsc        open(UNIT=55,FILE='nkout.m',STATUS='REPLACE')
cpetsc        write(55,*)"nksol=["
cpetsc        do i=1,6
cpetsc          write(55,*)i-1,nksolfnrm(i),log10(nksolfnrm(i)),nksoltime(i)-nksoltime(1),nksollinits(i),nksolfeval(i)
cpetsc        end do
cpetsc        write(55,*)"];"
cpetsc        close(55)
cpetsc      endif

cpetsc      if (snesdebug. eq. 1) then
              nni(1) = iwork(10)
              nli(1) = iwork(11)
              nfe(1) = iwork(12)
              nje(1) = iwork(13)
              npe(1) = iwork(14)
              nps(1) = iwork(15)
              ncfl(1) = iwork(16)
cpetsc      endif
*-------------------------------------------------------------------------
         else
            call newton (ffun,neq,yl,yldot0,yldot1,npsn,
     .                   rwork,iwork,jacnw,psolnw)
            njen = ijac
         endif

      endif               # End of large imeth if test

      iter = iter + 1

c...  convert solver variables back to plasma variables
      call convsr_vo1 (-1, -1, yl)  # was one call to convsr
      call convsr_vo2 (-1, -1, yl)  # was one call to convsr
      call convsr_aux (-1, -1)

c...  If nksol is in time-dependent mode, increment istep and toutlsod
      if (svrpkg .eq. 'nksol' .and. dtreal .lt. 1.e5) then
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

c_mpi      if (mype .eq. 0) then
c_mpi         if(svrpkg .eq. 'cvode') then
c_mpi            write(*,*) 'tout, te(0,2)/1.6e-19, mype, nx =',
c_mpi     .            tout, te(0,2)/1.6e-19, mype, nx
c_mpi            write(*,*) 'dt, nst, nfe, nni, nli, mpe = ',
c_mpi     .            ropt(4),iopt(4),iopt(5),iopt(7),iopt(15),mype
c_mpi            write(*,*) 'npe, ncfn, nps, ncfl = ',
c_mpi     .            iopt(14), iopt(8), iopt(16), iopt(17)
c_mpi         elseif(svrpkg .eq. 'kinsol') then
c_mpicvode            call fkinbbdopt(lenrpw, lenipw, nge)
c_mpi            write(6,1267) iopt(4), iopt(11), iopt(5),
c_mpi     .                    iopt(12), iopt(13),  nge
c_mpi 1267 format(1x,'nni=',i4,' npe=',i4,' nfe=',i4,' nps=',i4,
c_mpi     1       ' ncfl=',i4, ' nge =',i4)
c_mpi         endif
c_mpi      endif

*     -- increment and loop --
         istep = istep + 1
         go to 10


*  -- looping istep to nsteps is done, reset some parameters --
  200 continue

      ts = 0.
      istate = 1
      info(1) = 0
      iopts = inopt
      if (svrpkg.eq."cvode" .or. svrpkg.eq."kinsol") iopts=0
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

      if (svrpkg.eq.'cvode' .and. nx.eq.6 .and. ny.eq.0) then
         write(*,*) 'IOPT(4) NST =', IOPT(4)
         write(*,*) 'IOPT(5) NFE =', IOPT(5)
         write(*,*) 'IOPT(7) NNI =', IOPT(7)
         write(*,*) ni
         write(*,*) up
         write(*,*) te
         write(*,*) ti
         write(*,*) ng
      endif
      if (svrpkg.eq.'vodpk' .and. nx.eq.6 .and. ny.eq.2) then
         write(*,*) 'iwork(11) NST =', iwork(11)
         write(*,*) 'iwork(12) NFE =', iwork(12)
         write(*,*) 'iwork(13) NPE =', iwork(13)
         write(*,*) 'iwork(20) NNI =', iwork(20)
         write(*,*) ni
         write(*,*) up
         write(*,*) te
         write(*,*) ti
         write(*,*) ng
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
cpetsc*  Added for fixing parallel implementation (Feb. 18, 2008)
cpetsc      call set_indirect_address(0)   # sets ixp1 and ixm1 on all PEs
      call ueinit

      NLOCAL = neq

c ... Check that some indices have been passed properly
ccc      write(6,*) "[",mype,"] nx+ixmnbcl+ixmxbcl, ny+iymnbcl+iymxbcl, neq:",
ccc     .              nx+ixmnbcl+ixmxbcl, ny+iymnbcl+iymxbcl, neq,
ccc     .             ";    nxg+2, nyg+2, neqg:", nxg+2, nyg+2, neqg
ccc      call flush(6)

      if (svrpkg.eq."cvode") ITASK = 0
cc   use itask = 1 for the one-step mode
cxqx      ITASK = 1
cc   use jpre=0 for no preconditioner, jpre=1 is default
cxqx      JPRE = 0
      IGS = 0

cdb_solve      IF (MYPE .EQ. 0) THEN
cdb_solve        WRITE(6,9) NLOCAL
cdb_solve 9      FORMAT('Diagonal test problem, size NLOCAL =',I7)
cdb_solve        WRITE(6,11) NEQg
cdb_solve  11    FORMAT('Diagonal test problem, size NEQg =',I7)
cdb_solve        WRITE(6,12)
cdb_solve  12    FORMAT('  yl(i = 1,...,NEQg)'/)
cdb_solve        WRITE(6,103) rtol_pv, atol_pv, delt_pv
cdb_solve  103    FORMAT('rtol_pv, atol_pv, delt_pv =',3E10.1/)
cdb_solve        WRITE(6,14)
cdb_solve  14    FORMAT('Method is BDF/NEWTON/SPGMR'/
cdb_solve     1         'Diagonal preconditioner uses approximate Jacobian'/)
cdb_solve        WRITE(6,15) NPES
cdb_solve  15    FORMAT('Number of processors =',I4)
cdb_solve      ENDIF
C
c_mpicvode      if(svrpkg .eq. 'kinsol') then
c_mpicvode        call fpvecinitmpi(nlocal, neqg, ier)
c_mpicvode      else if(svrpkg .eq. 'cvode') then
c_mpicvode        CALL FPVINITMPI (uedgeComm, NLOCAL, NEQg, IER)
c_mpicvode      endif
ctaylorkinsol for neq      call fpvecinitmpi(nlocal, neq, ier)

ctaylorkinsol for npes=size      call mpi_comm_size(uedgeComm,size,ier)
ctaylorkinsol for npes = size
ctaylorkinsol for mype     call mpi_comm_rank(uedgeComm, rank, ier)
ctaylorkinsol for mype      mype = rank
ctaylorkinsol for       baseadd = mype * nlocal

cxqxkinsol      CALL FPVINITMPI (uedgeComm, NLOCAL, NEQg, IER)

C
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

      if(svrpkg .eq. 'kinsol') then
         iopt(1) = iprint
         ropt(3) = rlx
         ropt(6) = epscon1

         do i = 1, neq
           constr(i) = float(icnstr(i))
         enddo
c_mpicvode         call fpkinmalloc(neqg, ier)
c_mpicvode         call fkinbbdinit0(maxkd, maxlrst, msbpre, mu, ml, ier)

       elseif(svrpkg .eq. 'cvode') then

cccc_mpi         CALL FPVMALLOC (neqg, ts, yl, METH, ITMETH, IATOL,
cccc_mpi     .               rtol_pv, atol_pv, INOPT, IOPT, ROPT, IER)
            IF (IER .NE. 0) THEN
              WRITE(6,300) IER
  300         FORMAT(///' FPVMALLOC returned IER =',I5)
              STOP
            ENDIF
C
c        Use either of these next two lines for kinsol with bbd
c_mpicxqx         CALL FPVBBDIN (JPRE, IGS, maxkd, 0.0D0, MU, ML, IER)
cccc_mpi         CALL FPVBBDIN (JPRE, IGS, maxkd, delt_pv, MU, ML, IER)
           IF (IER .NE. 0) THEN
             WRITE(6,35) IER
  35         FORMAT(///' FPVBBDIN returned IER =',I5)
             STOP
           ENDIF
c        FCVSPGMR2 uses our own precond.; thus comment out FPVBBDIN lines
cccc_mpi         call FCVSPGMR2 (jpre, igs, maxkd, delt_pv)
      endif

CC
c xqx start - add timer for parallel
ctdr      CALL timer_init
ctdr      CALL mpi_barrier(mpi_comm_world,ier)
ctdr      call tsecnd(wtimestep_start)
ctdr      CALL timer(timestep_start)
c xqx end
cxqx ------------------end---------------------------------------------

ctdr  mpi stuff and diagnostics
c ==================================================================

cxqx      CALL jacmap

ctdr      call MPI_BARRIER(uedgeComm, ierr)

      if (ionecall .eq. 1) then
         call pandf (-1, -1, neq, 0., yl, yldot)
ccc      do iy = 0, ny+1
ccc      write(*,*) 'uedriv vol(2,iy), iy, id =',vol(2,iy), iy, mype+1
ccc      enddo

ccc         do 3331 iy = 0, ny+1
ccc            write(*,*) 'ng(3,iy,1), iy, id', ng(3,iy,1),iy,mype+1
ccc            write(*,*) 'resng(3,iy,1)/(vol*n0g), vol(3,iy), iy, id',
ccc     .             resng(3,iy,1)/(vol(3,iy)*n0g(1)),vol(3,iy),iy,mype+1
 3331    continue
      return
      endif

      return
      end
c****** End of subroutine uedriv_pll ************************************
c-----------------------------------------------------------------------
      subroutine newton (f, neq, y, rw1, rw2, npsn, wp, iwp, jacs, jsol)

c ... Perform Newton iterations to reduce right-hand sides returned by
c     function f.

      implicit none

c ... Input arguments:
      external f, jacs, jsol
      integer neq

c ... Work-array arguments:
      real rw1(neq), rw2(neq)
      real wp(*)       # matrix elements of LU
      integer iwp(*)   # dimensions and array indices for elements of LU

c ... In-out argument:
      real y(neq)

c ... Output argument:
      integer npsn

c ... Common blocks:
      Use(Dim)       # gives dimensions of variables not used here
      Use(Math_problem_size)   # neqmx
      Use(Share)     # cutlo
      Use(Newtaux)
      Use(Indexes)   # idxu,igyl
      Use(Jacaux)    # ylold,scrit
      Use(Ynorm)     # suscal,sfscal
      Use(Constraints) # rlx,rlxv

c ... Functions:
      real res_sum_dy
      real vnormnk

c ... Local variables:
      integer i, iu, imxchng
      real tp
      integer ierr
      character*80 msgjm
      integer nrcv, ierrjm, ijmgetmr
      real sumfl2

c ... Call the RHS routine, and compute initial sumf.
      tp = 0.
      call f (neq, tp, y, rw1)
      sumf(0) = 0.
      do i = 1, neq
         sumf(0) = abs(rw1(i)*sfscal(i))+sumf(0)
      enddo
      sumf(0) = sumf(0)/float(neq)

c ... Print initial sumf.
      write(*,*) 'Entering Newton; sumf(0) = ', sumf(0)

c ... Begin iterations.
      npsn = 0
   10 npsn = npsn + 1

c ... Calculate some diagnostics.
      call nwdiagnx (neq, rw1, npsn, 1)

c ... If variables have changed significantly since the previous
c     step, compute Jacobian and its LU decomposition.  Otherwise,
c     use the LU saved in wp.
      if (npsn .eq. 1 .or. res_sum_dy (neq, y, ylold) .gt. scrit) then
         call jacs (neq, y, rw1, 1.e10, rw2, wp, iwp)
         do i = 1, neq
            ylold(i) = y(i)
         enddo
      endif

c ... Solve J*x=rw1 and return x (in rw1) as solution.
      call jsol (neq, y, rw2, wp, iwp, rw1, ierr)

c ... Set correction vector and calculate multiplicative constant
c     for updating the solution.
      saux = 1.e30
      imxchng = 1
      do i = 1, neq
         ycor(i) = -rw1(i)/suscal(i)
         if (ycor(i) .ne. 0.) then
            saux1 = abs(rlx * y(i) / ycor(i))
            iu = idxu(igyl(i,1),igyl(i,2),1)
            if (i .eq. iu) saux1 =
     .         rlxv * (abs(y(i)*suscal(i)) / abs(ycor(i)*suscal(i)) + 2.)
            if (saux1 .lt. saux .and. abs(y(i)*suscal(i)) .gt. 1.e-10) then
               saux = saux1
               imxchng = i
            endif
         endif
      enddo
      saux = min(1., saux)
      sumnew = 0.
      sumrdy = 0.
      icsum = 0
      do i = 1, neq
         sumnew = sumnew + abs(ycor(i)*suscal(i))
         if (y(i) .ne. 0.) then
            sumrdy = sumrdy + abs(ycor(i) / y(i))
         else
            icsum = icsum + 1
         endif
         y(i) = y(i) + saux * ycor(i)
      enddo

      sumnew = sumnew / float(neq)
      sumrdy = sumrdy / float(neq-icsum)

c ... Calculate and print more diagnostics.
      sumnew1(npsn) = sumnew
      sumr1dy(npsn) = sumrdy
      saux2(npsn) = saux
      write(*,*) 'sumnew1(',npsn,') = ', sumnew,
     .      ', saux2(',npsn,') = ', saux,
     .      ', sumr1dy(',npsn,') = ', sumrdy
      do i = 1, neq
         ysave(2,i) = ysave(1,i)
         ysave(1,i) = y(i)
         rw2(i) = ycor(i) / (y(i) + cutlo)
      enddo
      call nwdiagnx (neq, rw2, npsn, 2)

c ... Call the RHS routine, and compute and print sumf and maximum-
c     change diagnostics.
      call f (neq, tp, y, rw1)
      sumf(npsn) = 0.
      do i = 1, neq
         sumf(npsn) = abs(rw1(i)*sfscal(i))+sumf(npsn)
      enddo
      sumfl2 = vnormnk (neq, rw1, sfscal)
      sumf(npsn) = sumf(npsn)/float(neq)
      write(*,*) 'sumf(',npsn,') = ', sumf(npsn), ', imxchng = ',
     .     imxchng, ', ix,iy = ', igyl(imxchng,1), ',',
     .     igyl(imxchng,2), ', fnrm =', sumfl2

c ... Test if too many iterations have been tried.
      if(npsn .ge. nmaxnewt .or. npsn .gt. 101) then
         call writeToLog('***** The newton iteration has not converged')
         call xerrab("")
      endif

c ... Test if user wants to quit.
ccc      ierrjm = ijmgetmr(msgjm,80,1,nrcv)
ccc      if (ierrjm .eq. 0) then
ccc         if(msgjm(1:nrcv) .eq. 'kaboom' .or. msgjm(1:nrcv) .eq. 'k') then
ccc            call xerrab("")
ccc         endif
ccc      endif

c ... Test for convergence.
      if (sumnew .gt. rwmin) go to 10

      return
      end
c***** end of subroutine newton ***********************
c-----------------------------------------------------------------------
      subroutine nwdiagnx (neq, rw, npsn, j)

      implicit none

      integer neq
      integer npsn
      integer j
      real rw(neq)

      Use(Math_problem_size)   # neqmx (for arrays not used here)
      Use(Newtaux)   # irwd,rwdmin,rwdmax

      integer iv
      real absrw

      rwdmin(npsn,j) = 1.0001*abs(rw(1))
      rwdmax(npsn,j) = 0.

      do iv = 1, neq
         absrw = abs(rw(iv))
         if (absrw .gt. rwdmax(npsn,j)) then
            rwdmax(npsn,j) = absrw
            irwd(npsn,j) = iv
         endif
         rwdmin(npsn,j) = min(absrw, rwdmin(npsn,j))
      enddo

      return
      end
c****** end of subroutine nwdiagnx **************************
c-----------------------------------------------------------------------
      real function res_sum_dy (neq, yl, ylold)

      implicit none

c ... Return normalized sum of changes in solver variables yl.

c ... Input arguments:
      integer neq
      real yl(neq), ylold(neq)

      Use(Share)   # cutlo

c ... Local variables:
      integer iv
      real sum

      sum = 0.
      do iv = 1, neq
         sum = sum + abs((yl(iv)-ylold(iv))/(yl(iv)+cutlo))
      enddo
      res_sum_dy = sum/float(neq)

      return
      end
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
c---------------------------------------------------------------------------c

      subroutine interptrans

c ... Routine calculates flux surface averages particle and
c ... energy sources for the UEDGE interpretive mode.
c ... Note in SOL, averages are from X-point to X-point; divertor
c ... omitted. Also, in the SOL, poloidal fluxes are ignored
c ... (model inaccurate); in core, average to 0. First, read input
c ... file, then read rdfitdata_new, then execute once with
c ... issfon=0 and ftol=1e20, and finally read this file.

      implicit none
      Use(Dim)            # nx,ny
      Use(Xpoint_indices) # ixpt1,ixpt2
      Use(Aux)            # ixmp
      Use(Interprettrans) # del_wicv,del_wecv,... for interpolation
      Use(Comgeo)         # gyf
      Use(Bcond)          # curcore,pcoree,pcorei
      Use(Parallv)         # nxg,nyg
      Use(Comflo)         # fniy,fngy
      Use(Rccoef)         # recycc
      Use(Phyvar)         # qe
      Use(Compla)         # ni,te,ti
      Use(Rhsides)        # psor,seec,seic,eqp
      Use(Coefeq)         # cfloyi,cfloye
      Use(Conduc)         # eqp
      Use(Fitdata)        # taudndt
      Use(Jacaux)         # issfon
      Use(Lsode)          # ftol
      Use(UEint)          # newgeo
      Use(Math_problem_size) # neqmx
      Use(Comtra)         # difni,kye,kyi


c ... del_wicv	# dWi/dt ion radial conve heat flux [W/m**3]
c ... del_wecv	# dWe/dt elec radial conve heat flux [W/m**3]
c ... del_wicdd	# dWi/dt ion radial diffu heat flux [W/m**3]; diagnos
c ... del_wecdd	# dWe/dt elec radial diffu heat flux [W/m**3]; diagnos
c ... del_wicd	# dWi/dt ion radial diffu heat flux [W/m**3]; pow bal
c ... del_wecd	# dWe/dt elec radial diffu heat flux [W/m**3]; pow bal
c ... del_cei	# dWe/dt=-dWi/dt elec-ion coll exchange [W/m**3]
c ... del_wivdp	# dWi/dt ion p-v type pressure-work terms [W/m**3]
c ... del_wevdp	# dWi/dt elec p-v type pressure-work terms [W/m**3]

c ... Local variables
      real arat,gamp_difmin,gamei_kyimin,gamee_kyemin
      real niface,teface,tiface,rdelyf
      real fracl,fracr,gamrecycc,sycore0
      real sycore(0:ny+1),sycoregyf(0:ny+1),vtoty(0:ny+1)
      integer ncells,ifld
      #Former Aux module variables
      integer ix,iy



c ... ##########################################################
c ... # Compute geometry, array sizes, case parameters
c ... ##########################################################

      issfon = 0
      ftol = 1.e50
      call exmain

c ... Set some constants
      ncells = ixpt2(1)-ixpt1(1)
      fracl = gyf(ixmp,1)/gyf(ixmp,0)
      fracr = gyf(ixmp,ny-1)/gyf(ixmp,ny)

c ... #######################################################
c ... # Read expt radial profile coeffs, compute fits, insert
c ... # into restart/save variables nis, tes, tis
c ... #######################################################

      call fitdata2svar

c ... ##########################################################
c ... # Evaluate RHS with fit profiles before doing analysis
c ... ##########################################################

      newgeo = 0
      call exmain
      issfon = 1
      ftol = 1.e-8

c ... ##########################################################
c ... # Calculate terms in eqns; 1st particle source/radial flux
c ... ##########################################################
c ... # For testing with UEDGE solution, compute curcore(1)
      gamrecycc = 0.
      gamp(0) = 0.
      sycore(0) = 0.
      sycoregyf(0) = 0.
      do ix = ixpt1(1)+1, ixpt2(1)
        if(curcore(1) == 0.) curcore(1) = curcore(1) + fniy(ix,0,1)*qe
	gamrecycc = gamrecycc + recycc(1)*fngy(ix,0,1)
	sycore(0) = sycore(0) + sy(ix,0)
        sycoregyf(0) = sycoregyf(0) + sy(ix,0)*gyf(ix,0)
      enddo

      pfmpg(0) = gyf(ixmp,0)*sycore(0)/sycoregyf(0)
      gamp(0) = (curcore(1)/qe - gamrecycc)/sycore(0)

      do iy = 1, ny   # inclu
        vtoty(iy) = 0.
	sycore(iy) = 0.
        sycoregyf(iy) = 0.
        del_sp(iy) = 0.
	del_dndt(iy) = 0.
	del_deedt(iy) = 0.
        del_deidt(iy) = 0.
        do ix = ixpt1(1)+1, ixpt2(1)
          vtoty(iy) = vtoty(iy) + vol(ix,iy)
	  sycore(iy) = sycore(iy) + sy(ix,iy)
          sycoregyf(iy) = sycoregyf(iy) + sy(ix,iy)*gyf(ix,iy)
        enddo
        do ix = ixpt1(1)+1, ixpt2(1)
	  del_sp(iy) = del_sp(iy) + psor(ix,iy,1)/vtoty(iy)
          del_dndt(iy) = del_dndt(iy) + ni(ix,iy,1)*vol(ix,iy) /
     .                                         (taudndt(iy)*vtoty(iy))
          del_deedt(iy) = del_deedt(iy) + 1.5*te(ix,iy)*ni(ix,iy,1)*
     .                           vol(ix,iy) / (taudeedt(iy)*vtoty(iy))
          del_deidt(iy) = del_deidt(iy) + 1.5*ti(ix,iy)*ni(ix,iy,1)*
     .                           vol(ix,iy) / (taudeidt(iy)*vtoty(iy))
        enddo
        gamp(iy) = gamp(iy-1)*sycore(iy-1)/sycore(iy) +
     .               (del_sp(iy)-del_dndt(iy))*vtoty(iy)/sycore(iy)
        pfmpg(iy) = gyf(ixmp,iy)*sycore(iy)/sycoregyf(iy)
        do ix = ixpt1(1)+1, ixpt2(1)
          facgam(ix,iy) = gyf(ix,iy)*sycore(iy)/sycoregyf(iy)
          floyd(ix,iy) = cfloyi*gamp(iy)*facgam(ix,iy)*sy(ix,iy)*0.5
        enddo
      enddo

      do ix = ixpt1(1)+1, ixpt2(1)
        floyd(ix,0) = floyd(ix,1) - (floyd(ix,2)-floyd(ix,1))*fracl
      enddo

c ... #######################################################
c ... # Calculate the energy sources and fluxes
c ... #######################################################
      do iy = 1, ny
	del_wicv(iy) = 0.
	del_wecv(iy) = 0.
	del_cei(iy) = 0.
	del_wivdp(iy) = 0.
	del_wevdp(iy) = 0.
	del_wicd(iy) = 0.
	del_wecd(iy) = 0.
	del_wicdd(iy) = 0.
	del_wecdd(iy) = 0.
	do ix = ixpt1(1) + 1, ixpt2(1)
          del_wicv(iy) = del_wicv(iy) +
     .                   (floyd(ix,iy)  *(ti(ix,iy)  +ti(ix,iy+1)) -
     .                    floyd(ix,iy-1)*(ti(ix,iy-1)+ti(ix,iy)))/vtoty(iy)
          del_wecv(iy) = del_wecv(iy) +
     .                   (floyd(ix,iy)  *(te(ix,iy)  +te(ix,iy+1)) -
     .                    floyd(ix,iy-1)*(te(ix,iy-1)+te(ix,iy)))/vtoty(iy)
          del_cei(iy)  = del_cei(iy) + eqp(ix,iy)*
     .                   (te(ix,iy)-ti(ix,iy))*vol(ix,iy) / vtoty(iy)
          del_wivdp(iy) = del_wivdp(iy) + seic(ix,iy) / vtoty(iy)
          del_wevdp(iy) = del_wevdp(iy) + seec(ix,iy) / vtoty(iy)
        enddo
        del_wicd(iy)  = - del_wicv(iy) - del_cei(iy) + del_wivdp(iy)
        del_wecd(iy)  = - del_wecv(iy) + del_cei(iy) + del_wevdp(iy)
        del_wicdd(iy) = del_wicdd(iy) + del_wicd(iy)  # must sum to zero
        del_wecdd(iy) = del_wecdd(iy) + del_wecd(iy)  # must sum to zero
      enddo

c ... # Fill in guard-cell values using linear interpolation

      del_sp(0)    = del_sp(1)    - (del_sp(2)-del_sp(1))*fracl
      del_dndt(0)  = del_dndt(1)  - (del_dndt(2)-del_dndt(1))*fracl
      del_deedt(0) = del_deedt(1) - (del_deedt(2)-del_deedt(1))*fracl
      del_deidt(0) = del_deidt(1) - (del_deidt(2)-del_deidt(1))*fracl
      del_witot(0) = del_witot(1) - (del_witot(2)-del_witot(1))*fracl
      del_wetot(0) = del_wetot(1) - (del_wetot(2)-del_wetot(1))*fracl
      del_wicv(0)  = del_wicv(1)  - (del_wicv(2)-del_wicv(1))*fracl
      del_wecv(0)  = del_wecv(1)  - (del_wecv(2)-del_wecv(1))*fracl
      del_wicd(0)  = del_wicd(1)  - (del_wicd(2)-del_wicd(1))*fracl
      del_wecd(0)  = del_wecd(1)  - (del_wecd(2)-del_wecd(1))*fracl
      del_cei(0)   = del_cei(1)   - (del_cei(2)-del_cei(1))*fracl
      del_wivdp(0) = del_wivdp(1) - (del_wivdp(2)-del_wivdp(1))*fracl
      del_wevdp(0) = del_wevdp(1) - (del_wevdp(2)-del_wevdp(1))*fracl

      del_sp(ny+1)   = del_sp(ny)   + (del_sp(ny)-del_sp(ny-1))*fracr
      del_dndt(ny+1) = del_dndt(ny) + (del_dndt(ny)-del_dndt(ny-1))*fracr
      del_deedt(ny+1)= del_deedt(ny)+ (del_deedt(ny)-del_deedt(ny-1))*fracr
      del_deidt(ny+1)= del_deidt(ny)+ (del_deidt(ny)-del_deidt(ny-1))*fracr
      del_witot(ny+1)= del_witot(ny)+ (del_witot(ny)-del_witot(ny-1))*fracr
      del_wetot(ny+1)= del_wetot(ny)+ (del_wetot(ny)-del_wetot(ny-1))*fracr
      del_wicv(ny+1) = del_wicv(ny) + (del_wicv(ny)-del_wicv(ny-1))*fracr
      del_wecv(ny+1) = del_wecv(ny) + (del_wecv(ny)-del_wecv(ny-1))*fracr
      del_wicd(ny+1) = del_wicd(ny) + (del_wicd(ny)-del_wicd(ny-1))*fracr
      del_wecd(ny+1) = del_wecd(ny) + (del_wecd(ny)-del_wecd(ny-1))*fracr
      del_cei(ny+1)  = del_cei(ny)  + (del_cei(ny)-del_cei(ny-1))*fracr
      del_wivdp(ny+1)= del_wivdp(ny)+ (del_wivdp(ny)-del_wivdp(ny-1))*fracr
      del_wevdp(ny+1)= del_wevdp(ny)+ (del_wevdp(ny)-del_wevdp(ny-1))*fracr

c ... ##########################################################
c ... # Calculate energy fluxes (UEDGE must know pcoree,i)     #
c ... ##########################################################
c ... # For testing with UEDGE solution, compute pcoree,i
      if(pcorei == 0.) then
        do ix = ixpt1(1)+1, ixpt2(1)
          pcorei = pcorei + feiy(ix,0)
        enddo
      endif
      if(pcoree == 0.) then
        do ix = ixpt1(1)+1, ixpt2(1)
          pcoree = pcoree + feey(ix,0)
        enddo
      endif

      gamei(0) = pcorei/sycore(0) - cfloyi*gamp(0)*0.5*
     .                                   (ti(ixmp,0)+ti(ixmp,1))
      gamee(0) = pcoree/sycore(0) - cfloye*gamp(0)*0.5*
     .                                   (te(ixmp,0)+te(ixmp,1))

      do iy = 1, ny
        arat = sycore(iy-1)/sycore(iy)
	vyn_int(iy) = 0.                          # initialize
	vyee_int(iy) = 0.                          # initialize
	vyei_int(iy) = 0.                          # initialize
        gamei(iy)= gamei(iy-1)*arat + (del_wicd(iy)-del_deidt(iy))*
     .                                        vtoty(iy)/sycore(iy)
        gamee(iy)= gamee(iy-1)*arat + (del_wecd(iy)-del_deedt(iy))*
     .                                        vtoty(iy)/sycore(iy)
      enddo

c ... #################################################################
c ... # Calculate effective diffusion/convection coeffs from data above
c ... #################################################################

      do iy = 0, ny
        niface = 0.5*(ni(ixmp,iy+1,1) + ni(ixmp,iy,1))
        teface = 0.5*(te(ixmp,iy+1) + te(ixmp,iy))
        tiface = 0.5*(ti(ixmp,iy+1) + ti(ixmp,iy))
        dif_int(iy) = - gamp(iy)*pfmpg(iy) /
     .               ((ni(ixmp,iy+1,1) - ni(ixmp,iy,1))*gyf(ixmp,iy))
        if(dif_int(iy) < dif_use_min) then
          dif_int(iy) = dif_use_min
          gamp_difmin= -dif_int(iy)*( (ni(ixmp,iy+1,1)-ni(ixmp,iy,1))*
     .                            gyf(ixmp,iy) ) / pfmpg(iy)
          vyn_int(iy) = (gamp(iy)-gamp_difmin)/niface
        elseif(dif_int(iy) > dif_use_max) then
          dif_int(iy) = dif_use_max
          gamp_difmin= -dif_int(iy)*( (ni(ixmp,iy+1,1)-ni(ixmp,iy,1))*
     .                            gyf(ixmp,iy) ) / pfmpg(iy)
          vyn_int(iy) = (gamp(iy)-gamp_difmin)/niface
        elseif(iy > 1 .and. dif_int(iy) > maxchgdiff*dif_int(iy-1)) then
          dif_int(iy) = maxchgdiff*dif_int(iy-1)
          gamp_difmin= -dif_int(iy)*( (ni(ixmp,iy+1,1)-ni(ixmp,iy,1))*
     .                            gyf(ixmp,iy) ) / pfmpg(iy)
          vyn_int(iy) = (gamp(iy)-gamp_difmin)/niface
        elseif(iy > 1 .and. dif_int(iy) < dif_int(iy-1)/maxchgdiff) then
          dif_int(iy) = max(dif_int(iy-1)/maxchgdiff, dif_use_min)
          gamp_difmin= -dif_int(iy)*( (ni(ixmp,iy+1,1)-ni(ixmp,iy,1))*
     .                            gyf(ixmp,iy) ) / pfmpg(iy)
          vyn_int(iy) = (gamp(iy)-gamp_difmin)/niface
        endif

        kyi_int(iy) = -gamei(iy)*pfmpg(iy) /
     .          ( niface*(ti(ixmp,iy+1)-ti(ixmp,iy))*gyf(ixmp,iy) )
        kyi_int(iy) = min(kyi_int(iy), kyi_use_max)
        if (kyi_int(iy) < kyi_use_min) then
          kyi_int(iy) = kyi_use_min
          gamei_kyimin = -kyi_int(iy)*( niface*
     .                 (ti(ixmp,iy+1)-ti(ixmp,iy))*gyf(ixmp,iy) )
          if(cfloyi > 0.) vyei_int(iy) = (gamei(iy)-gamei_kyimin)/
     .                                    (cfloyi*niface*tiface)
       elseif (kyi_int(iy) > kyi_use_max) then
          kyi_int(iy) = kyi_use_max
          gamei_kyimin = -kyi_int(iy)*( niface*
     .                 (ti(ixmp,iy+1)-ti(ixmp,iy))*gyf(ixmp,iy) )
          if(cfloyi > 0.) vyei_int(iy) = (gamei(iy)-gamei_kyimin)/
     .                                    (cfloyi*niface*tiface)
        endif

        kye_int(iy) = -gamee(iy)*pfmpg(iy) /
     .          ( niface*(te(ixmp,iy+1)-te(ixmp,iy))*gyf(ixmp,iy) )
        kye_int(iy) = min(kye_int(iy), kye_use_max)
        if(kye_int(iy) < kye_use_min) then
          kye_int(iy) = kye_use_min
          gamee_kyemin = -kye_int(iy)*( niface*
     .                (te(ixmp,iy+1)-te(ixmp,iy))*gyf(ixmp,iy) )
          if(cfloye > 0.) vyee_int(iy) = (gamee(iy)-gamee_kyemin)/
     .                                  (cfloye*niface*teface)
        elseif(kye_int(iy) > kye_use_max) then
          kye_int(iy) = kye_use_max
          gamee_kyemin = -kye_int(iy)*( niface*
     .                (te(ixmp,iy+1)-te(ixmp,iy))*gyf(ixmp,iy) )
          if(cfloye > 0.) vyee_int(iy) = (gamee(iy)-gamee_kyemin)/
     .                                  (cfloye*niface*teface)
        endif
      enddo

c ... # Set guard-cell values by extrapolation
c      rdelyf = gy(ixmp,2)/gy(ixmp,1)
c      dif_int(0) = dif_int(1) - (dif_int(2)-dif_int(1))*rdelyf
c      kyi_int(0) = kyi_int(1) - (kyi_int(2)-kyi_int(1))*rdelyf
c      kye_int(0) = kye_int(1) - (kye_int(2)-kye_int(1))*rdelyf
c      vyn_int(0) = vyn_int(1) - (vyn_int(2)-vyn_int(1))*rdelyf
      dif_int(ny+1) = dif_int(ny)
      kyi_int(ny+1) = kyi_int(ny)
      kye_int(ny+1) = kye_int(ny)

      if(isdifuseinterp==1) then # use only dif_int, etc. & min/max values
        do ifld = 1, nisp
          difni(ifld) = 0.
        enddo
        kye = 0.
        kyi = 0.
c ...   Set values over SOL and PF; core reset next
      do iy = 0, ny+1
        do ix = 0, nx+1
          do ifld = 1, nisp
            if (ziin(ifld) > 0.) then
              vy_use(ix,iy,ifld) = 0.
                # now for PF and SOL regions
                if(iy<=iysptrx .and. (ix<=ixpt1(1).or.ix>=ixpt2(1)))
     .                                                       then
                  dif_use(ix,iy,ifld) = dif_int(iysptrx) + difni_pf
                elseif(iy > iysptrx) then
                  dif_use(ix,iy,ifld) = dif_int(iysptrx) + difni_sol
                endif
              endif
            enddo
              if(iy<=iysptrx .and. (ix<=ixpt1(1).or.ix>=ixpt2(1)))
     .                                                       then
                kye_use(ix,iy) = kye_int(iysptrx) + kye_pf
                kyi_use(ix,iy) = kyi_int(iysptrx) + kyi_pf
              elseif(iy > iysptrx) then
                kye_use(ix,iy) = kye_int(iysptrx) + kye_sol
                kyi_use(ix,iy) = kyi_int(iysptrx) + kyi_sol
              endif
          enddo
        enddo

c ...  Now reset values in core region from interpolation profiles

       do iy = 0, iysptrx
	 do ix = ixpt1(1)+1, ixpt2(1)
	   do ifld = 1, nisp
	     dif_use(ix,iy,ifld) = min(dif_int(iy),dif_use_max)
             dif_use(ix,iy,ifld) = max(dif_use(ix,iy,ifld),dif_use_min)
             vy_use(ix,iy,ifld) = vyn_int(iy)
           enddo
           kye_use(ix,iy) = min(kye_int(iy),kye_use_max)
           kye_use(ix,iy) = max(kye_use(ix,iy),kye_use_min)
           kyi_use(ix,iy) = min(kyi_int(iy),kyi_use_max)
           kyi_use(ix,iy) = max(kyi_use(ix,iy),kyi_use_min)
         enddo
       enddo

       endif  #test on isdifuseinterp for resetting dif_use, etc

c ... Adjust SOL profs next to core to be smooth radially
      if (isadjsolprof == 1) call adjsolprof

      return
      end
c ***  End of subroutine interptrans ---------------------------------- c

      subroutine adjsolprof

c ... Routine smoothly continues nis, tes, and tis profiles into SOL
c ... Use when tanh/spline fits for pedestal are connected to SOL soln
c ... Thus, only modifies profs for iy>iysptr and ixpt1<ix<ixpt2+1
c ... Note tis is fit only to iysptrx+1, while nis,tes are fit to iysptrx+2
c ... Thus, first extrapolate tis at iysptrx+1 (instead of pres average)
c ... Assume equal spaced mesh points here - can be done better

      implicit none
      Use(Interp)          #nis,tes,tis
      Use(Xpoint_indices)  #ixpt1,ixpt2
      Use(Dim)             #nx,iy
      Use(Interprettrans)  #denrdrop,terdrop,tirdrop

c ... Local variables
      integer ix,iy,ifld
      real niscalfac(0:nx+1),tescalfac(0:nx+1),tiscalfac(0:nx+1)
c ... Note: niscalfac only used for ixpt1(1)+1:ixpt2(1), but gnu compiler
c ...       does not accept dimensioning with ixpt1(1) and ixpt2(1)

c ... Extrapolate Ti; assume equal mesh spacing - can be more accurrate
      do ix = ixpt1(1)+1, ixpt2(1)
        tis(ix,iysptrx+2) = 2.*tis(ix,iysptrx+1) - tis(ix,iysptrx)
      enddo

      do ix = ixpt1(1)+1, ixpt2(1)
        niscalfac(ix) = denrdrop*nis(ix,iysptrx+2,1)/
     .                                          nis(ix,iysptrx+3,1)
        tescalfac(ix) = terdrop*tes(ix,iysptrx+2)/tes(ix,iysptrx+3)
        tiscalfac(ix) = tirdrop*tis(ix,iysptrx+2)/tis(ix,iysptrx+3)
      enddo

c ... Rescale SOL values in range ixpt1 < ix <= ixpt2
      do iy = iysptrx+3, ny+1
	do ix = ixpt1(1)+1, ixpt2(1)
	  nis(ix,iy,1) = niscalfac(ix)*nis(ix,iy,1)
	  tes(ix,iy) = tescalfac(ix)*tes(ix,iy)
	  tis(ix,iy) = tiscalfac(ix)*tis(ix,iy)
        enddo
      enddo

      return
      end
c ***  End of subroutine adjsolprofs ----------------------------c

c *** ----------------------------------------------------------------- c
      subroutine fitdata2svar

c ..  Evaluate profiles with Osborne fit data and store them in
c ..  nis, tes, tis variables

      implicit none
      Use(Dim)              # ny
      Use(Xpoint_indices)   # iysptrx
      Use(Interp)           # nis,tes,tis
      Use(Fitdata)          # nefit,tefit,tifit,isprof_coef
      Use(Phyvar)           # ev
      Use(Comgeo)           # gyf
      Use(Interprettrans)   # isadjsolprof

# local parameters; user may want to change iyendf
      integer ix,iy
      integer iyendf
      iyendf = iysptrx

      write(6,*) 'Now in fitdata2svar ****, isprof_coef =',isprof_coef
      call flush(6)

      if(isprof_coef == 1)  then  # compute profiles from fit coefficients

c ..    Read first set ne, Te and Ti data and evaluate via tanh or splines
        ifitset = 1
        call readne_dat('netanh1.dat')
        call readte_dat('tetanh1.dat')
        call readti_dat('tispline1.dat')

c ..    Now fit the data with tanh for ne and te & spline for ti
c ..    1-D results in tefit(0:ny+1,1), tifit(0:ny+1,1) [keV]
c ..    also, 1-D nefit(0:ny+1,1) [1/m**3]

        call fit_neteti

c ..    If a 2nd data set exists, read, and compute 2nd set of profiles
        if(isdndtfitdat == 1) then
          ifitset = 2
          call readne_dat('netanh2.dat')
          call readte_dat('tetanh2.dat')
          call readti_dat('tispline2.dat')

          call fit_neteti
        endif

      elseif(isprof_coef == 0) then  #instead interpolate from profile files

        call interp_neteti

      endif

c ..  Construct fits as specified combination of data sets 1 and 2
      do iy = 0, ny+1
         nefituse(iy) = fitfrac1*nefit(iy,1)+(1.-fitfrac1)*nefit(iy,2)
         tefituse(iy) = fitfrac1*tefit(iy,1)+(1.-fitfrac1)*tefit(iy,2)
         tifituse(iy) = fitfrac1*tifit(iy,1)+(1.-fitfrac1)*tifit(iy,2)
      enddo

c ..  Set restart variables in core region only from data set 1
      do iy = 0, iyendf
         do ix = ixpt1(1)+1, ixpt2(1)
           nis(ix,iy,1) = nefituse(iy)
           tes(ix,iy)   = 1e3*ev*tefituse(iy)
           tis(ix,iy)   = 1e3*ev*tifituse(iy)
         enddo
      enddo

c ..  Extend to cells beyond iyendf (must extrapolate ti if seperatrix)
      do iy = iyendf+1, iysptrx+2
        do ix = ixpt1(1)+1, ixpt2(1)
          nis(ix,iy,1) = nefituse(iy)
          tes(ix,iy)   =  1e3*ev*max(tefituse(iy), 1.e-3)
          if (iy <= iysptrx) then
            tis(ix,iy)   = 1e3*ev*tifituse(iy)
          elseif(iy == iysptrx+1) then
            tis(ix,iy)   = 1.e3*ev*( tifituse(iy-1) +
     .                            (tifituse(iy-1)-tifituse(iy-2))*
     .                              gyf(ix,iy-1)/gyf(ix,iy) )
          elseif (iy == iysptrx+2) then
            tis(ix,iy) = 0.5*(tis(ix,iy) + tis(ix,iy-1))
          endif
        enddo
      enddo

c ..  If second data set exists, compute taudndt, taudeedt, taudeidt
      if(isdndtfitdat == 1) then
        do iy = 0, iysptrx
          taudndt(iy) = 0.5*tim_interval_fit*
     .                             (nefit(iy,2) + nefit(iy,1))/
     .                             (nefit(iy,2) - nefit(iy,1))
          taudeedt(iy) = 0.5*tim_interval_fit*
     .       (nefit(iy,2)*tefit(iy,2) + nefit(iy,1)*tefit(iy,1))/
     .       (nefit(iy,2)*tefit(iy,2) - nefit(iy,1)*tefit(iy,1))
          taudeidt(iy) = 0.5*tim_interval_fit*
     .       (nefit(iy,2)*tifit(iy,2) + nefit(iy,1)*tifit(iy,1))/
     .       (nefit(iy,2)*tifit(iy,2) - nefit(iy,1)*tifit(iy,1))
          taudndt(iy)  = max(taudndt(iy), -tim_chng_max)
          taudndt(iy)  = min(taudndt(iy),  tim_chng_max)
          taudeedt(iy) = max(taudeedt(iy), -tim_chng_max)
          taudeedt(iy) = min(taudeedt(iy),  tim_chng_max)
          taudeidt(iy) = max(taudeidt(iy), -tim_chng_max)
          taudeidt(iy) = min(taudeidt(iy),  tim_chng_max)
        enddo
      endif

c ... Adjust SOL profs next to core to be smooth radially
      if (isadjsolprof == 1) call adjsolprof

      return
      end
c ***  End of subroutine fitdata2svar  -------------------------------c
c----------------------------------------------------------------------c

c************************************************************
c Subroutine to interpolate Osborne/Groebner DIII-D profiles
c************************************************************

      subroutine interp_neteti

      implicit none
      Use(Dim)              # ny,num_elem
      Use(Xpoint_indices)   # iysptrx
      Use(Comgeo)           # psinormc,yyc
      Use(Fitdata)          # psishift,nefit,epsi_fit,eprofile_fit,
                            # isdndtfitdat,psi_s,isprofvspsi,yyc_fit

c...  Local variables
      integer ir,iyu,idata,irs,ip
      real fracm

      do idata = 1, isdndtfitdat+1

c...    FIRST do the density profile, begin by reading data
        if (isprofvspsi == 1) then  #profiles provided vs psi
          if (idata==1) call read_exp_fit('ne_vs_psi_expt1')
          if (idata==2) call read_exp_fit('ne_vs_psi_expt2')
        else   #profiles provided vs r-r_sep
          if (idata==1) call read_exp_fit('ne_vs_r_expt1')
          if (idata==2) call read_exp_fit('ne_vs_r_expt2')
        endif

c...    Shift epsi_fit used by psishift to allow control of separatrix
        do ir = 1, num_elem
          psi_s(ir) = psishift + epsi_fit(ir)
        enddo

c...    Interpolate values from data file onto UEDGE mesh
        nefit(0:ny+1,idata) = 1.e-5  #initialize array
        irs = 1  #start location for interp intersection search
        do iyu = 0, ny+1
          do ir = irs, num_elem
            if (isprofvspsi==1 .and. psi_s(ir) >= psinormc(iyu)) then
              fracm=(psi_s(ir)-psinormc(iyu))/(psi_s(ir)-psi_s(ir-1))
              nefit(iyu,idata) = 1.e20*( fracm*eprofile_fit(ir-1) +
     &                             (1.-fracm)*eprofile_fit(ir) )
              irs = ir-1
              break
            elseif (isprofvspsi==0 .and. yyc_fit(ir) >= yyc(iyu)) then
              fracm=(yyc_fit(ir)-yyc(iyu))/(yyc_fit(ir)-yyc_fit(ir-1))
              nefit(iyu,idata) = 1.e20*( fracm*eprofile_fit(ir-1) +
     &                             (1.-fracm)*eprofile_fit(ir) )
              irs = ir-1
              break
            endif
          enddo
        enddo

c...    SECOND, do the Te profile, begin by reading data
        if (isprofvspsi == 1) then  #profiles provided vs psi
          if (idata==1) call read_exp_fit('te_vs_psi_expt1')
          if (idata==2) call read_exp_fit('te_vs_psi_expt2')
        else   #profiles provided vs r-r_sep
          if (idata==1) call read_exp_fit('te_vs_r_expt1')
          if (idata==2) call read_exp_fit('te_vs_r_expt2')
        endif

c...    Shift epsi_fit used by psishift to allow control of separatrix
        do ir = 1, num_elem
          psi_s(ir) = psishift + epsi_fit(ir)
        enddo

c...    Interpolate values from data file onto UEDGE mesh
        tefit(0:ny+1,idata) = 1.e-5   #intialize array
        irs = 1  #start location for interp intersection search
        do iyu = 0, ny+1
          do ir = irs, num_elem
            if (isprofvspsi==1 .and. psi_s(ir) >= psinormc(iyu)) then
              fracm=(psi_s(ir)-psinormc(iyu))/(psi_s(ir)-psi_s(ir-1))
              tefit(iyu,idata) = fracm*eprofile_fit(ir-1) +
     &                          (1.-fracm)*eprofile_fit(ir)
              irs = ir - 1
              break
            elseif (isprofvspsi==0 .and. yyc_fit(ir) >= yyc(iyu)) then
              fracm=(yyc_fit(ir)-yyc(iyu))/(yyc_fit(ir)-yyc_fit(ir-1))
              tefit(iyu,idata) = fracm*eprofile_fit(ir-1) +
     &                           (1.-fracm)*eprofile_fit(ir)
              irs = ir-1
              break
            endif
          enddo
        enddo

c...    THIRD, do the Ti profile, begin by reading data
        if (isprofvspsi == 1) then  #profiles provided vs psi
          if (idata==1) call read_exp_fit('ti_vs_psi_expt1')
          if (idata==2) call read_exp_fit('ti_vs_psi_expt2')
        else   #profiles provided vs r-r_sep
          if (idata==1) call read_exp_fit('ti_vs_r_expt1')
          if (idata==2) call read_exp_fit('ti_vs_r_expt2')
        endif

c...    Shift epsi_fit used by psishift to allow control of separatrix
        do ir = 1, num_elem
          psi_s(ir) = psishift + epsi_fit(ir)
        enddo

c...    Interpolate values from data file onto UEDGE mesh
        tifit(0:ny+1,idata) = 1.e-5   #intialize array
        irs = 1  #start location for interp intersection search
        do iyu = 0, ny+1
          do ir = irs, num_elem
            if (isprofvspsi==1 .and. psi_s(ir) >= psinormc(iyu)) then
              fracm=(psi_s(ir)-psinormc(iyu))/(psi_s(ir)-psi_s(ir-1))
              tifit(iyu,idata) = fracm*eprofile_fit(ir-1) +
     &                          (1.-fracm)*eprofile_fit(ir)
              irs = ir - 1
              break
            elseif (isprofvspsi==0 .and. yyc_fit(ir) >= yyc(iyu)) then
              fracm=(yyc_fit(ir)-yyc(iyu))/(yyc_fit(ir)-yyc_fit(ir-1))
              tifit(iyu,idata) = fracm*eprofile_fit(ir-1) +
     &                           (1.-fracm)*eprofile_fit(ir)
              irs = ir-1
              break
            endif
          enddo
        enddo

       enddo   # End of large loop over idata

      return
      end
c ***  End of subroutine interp_neteti --------------------------------c
c----------------------------------------------------------------------c

c *** ----------------------------------------------------------------- c
c     subroutine outputstats
c
c     This is just a dummy routine to allow nopetsc to build but still
c     call the function outputstatspetsc stored in svr/petsc-uedge.F90.
c     From within the driver script, you should call bbb.outputstats()
c     to activate this subroutine.
c
c     If you want this function to do something for the nopetsc build
c     lemme know (Mike McCourt) and I'll see what I can do. Right now it
c     only outputs results for the par and ser builds.
c *** ----------------------------------------------------------------- c
      subroutine outputstats
      implicit none

      integer :: one=1

cpetsc      if (one.eq.1) then
cpetsc        call outputstatspetsc
cpetsc      else
      write(6,*) "I am in outputstats in a nopetsc build"
      write(6,*) "There is nothing to be done"
cpetsc      end if

      return
      end
c ***  End of subroutine outputstats  -------------------------------c
