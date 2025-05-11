c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"

c-----------------------------------------------------------------------
      subroutine PandfRHS_interface(neq, time, yl, yldot)
c ... Interface for pandf rhs calculation for nksol only (added by. J.Guterl)
          implicit none
          Use(Math_problem_size) # neqmx
          Use(ParallelEval)      # ParallelPandf1
          Use(Cdv)             # pandfitertag
          integer neq
          real time, yl(neqmx),yldot(neq)

c!omp     if (ParallelPandf1.gt.0) then
c!omp         pandfitertag = "OMP iter="
c!omp         call OMPPandf1Rhs(neq, time, yl, yldot)
c!omp     else
              pandfitertag = "iter="
              call pandf(-1, -1, neq, time, yl, yldot)
c!omp     endif

      end subroutine PandfRHS_interface


      SUBROUTINE initialize_ranges(xc, yc, xlinc_loc, xrinc_loc, yinc_loc)
      IMPLICIT NONE
      integer xc, yc, ixmp2, jx, xlinc_loc, xrinc_loc, yinc_loc
      logical xccuts, xcturb
      Use(Dim)
      Use(UEpar)
      Use(Comflo)
      Use(Compla)
      Use(Selec)
      Use(Grid)
      Use(Share)
      Use(Bcond)
      Use(Indices_domain_dcl)
      Use(Xpoint_indices)
      Use(Turbulence)
      Use(Imprad)
      
c ... Set variable controlling upper limit of species loops that
c     involve ion-density sources, fluxes, and/or velocities.
      nfsp = nisp
      if (isimpon .eq. 3 .or. isimpon .eq. 4) nfsp = nhsp

      if ( (xc .lt. 0) .or. 
     .     ((0<=yc).and.(yc-yinc_loc<=0)).and.isjaccorall==1 ) then  
                                              # use full ix range near yc=0
                                              # with integrated core flux BC
         i1 = 0  # 1-ixmnbcl
         i2 = 1
         i2p = 1
         i3 = 0  # 1-ixmnbcl
         i4 = 0  # 1-ixmnbcl
         i5 = nx
         i5m = nx-1
         i6 = nx+1  # nx+ixmxbcl
         i7 = nx+1  # nx+ixmxbcl
         i8 = nx+1  # nx+ixmxbcl
      else
         i1 = max(0, xc-xlinc_loc-1)
         i2 = max(1, xc-xlinc_loc)
         i2p = max(1, xc-xrinc_loc-1)
         i3 = xc-xlinc_loc     # not used as loop indice (can be < 0)
         i4 = max(0, xc-xlinc_loc)
         i5 = min(nx, xc+xrinc_loc)
         i5m = min(nx-1, xc+xrinc_loc)
         i6 = min(nx+1, xc+xrinc_loc+1)
         i7 = xc+xrinc     # not used as loop indice (can be > nx)
         i8 = min(nx+1, xc+xrinc)
      endif
      if (yc .lt. 0) then
         j1 = 0 
         j1p = 0
         j2 = 1
         j2p = 1
         j3 = 0 
         j4 = 0 
         j5 = ny
         j5m = ny-1
         j6 = ny+1 
         j5p = ny
         j6p = ny+1
         j7 = ny+1 
         j8 = ny+1 
      else
         j1 = max(0, yc-yinc_loc-1)
         j2 = max(1, yc-yinc_loc)
         j1p = max(0, yc-yinc_loc-2)
         j2p = max(1, yc-yinc_loc-1)
         j3 = yc-yinc_loc
         j4 = max(0, yc-yinc_loc)
         j5 = min(ny, yc+yinc_loc)
         j5m = min(ny-1, yc+yinc_loc)
         j6 = min(ny+1, yc+yinc_loc)
         j5p = min(ny, yc+yinc_loc+1)
         j6p = min(ny+1, yc+yinc_loc+1)
c         j6 = min(ny+1, yc+yinc_loc+1)
         j7 = yc+yinc_loc
         j8 = min(ny+1, yc+yinc_loc)
      endif

c...  We will expand the range of possible responses when perturbing the
c...  plasma in a cell near one of the cuts.
      xccuts = .false.
      do jx = 1, nxpt
        if ( (xc-xlinc_loc<=ixpt1(jx)+1) .and. (xc+xrinc_loc+1>=ixpt1(jx)) .and.
     .       (yc-yinc_loc<=iysptrx1(jx)) .and. (iysptrx1(jx)>0) ) xccuts=.true.
        if ( (xc-xlinc_loc<=ixpt2(jx)+1) .and. (xc+xrinc_loc+1>=ixpt2(jx)) .and.
     .       (yc-yinc_loc<=iysptrx2(jx)) .and. (iysptrx2(jx)>0) ) xccuts=.true.
      enddo

c...  We must expand the range of ix in the vicinity of cells on
c...  which turbulent transport depends.
      xcturb = .false.
      do jx = 1, nxpt
         xcturb = xcturb .or. (xc.eq.ixlb(jx).and.ixmnbcl==1) .or.
     .                        (xc.eq.(ixrb(jx)+1).and.ixmxbcl==1)
      enddo
c     Set ix index for outer midplane turbulence
      if (isudsym==1) then
         ixmp2 = nxc + 2
      elseif (geometry=='dnull' .or. geometry(1:9)=="snowflake" .or.
     .        geometry=="dnXtarget" .or. geometry=="isoleg") then
         ixmp2 = ixmdp(2) + 1
      else
         ixmp2 = ixpt1(1) + nxomit + 3*(ixpt2(1)-ixpt1(1))/4
      endif
      xcturb = xcturb .or. (xc.eq.ixmp2)
      xcturb = xcturb .and. (isturbnloc.eq.1)
c...  NOTE: For a full double-null configuration, if there are 2 separatrices
c...  we use the innermost one (at iysptrx) to define the radial boundary
c...  of the turbulence.
      if (isturbcons .eq. 1) then
         xcturb = xcturb .and. yc .eq. iysptrx+1
      elseif (isturbcons .eq. 2) then
         xcturb = xcturb .and. yc .ge. iysptrx+1-diffuslimit
      else
         xcturb = xcturb .and. yc .ge. iysptrx+1
      endif

      if (xccuts .or. xcturb) then
         i1 = 0  
         i2 = 1
         i3 = 0  
         i4 = 0  
         i5 = nx
         i6 = nx+1 
         i7 = nx+1 
         i8 = nx+1 
      endif

c...  Define range for source terms to minimize calls to adpak-based routines
            ixs = i2
            ixf = i5
            iys = j2
            iyf = j5
            ixs1 = i1
            ixf6 = i6
            iys1 = j1
            iyf6 = j6
c...  Reset ioniz. and rad. indices to a point if this is a Jacobian calc.
         if (xc .ge.0 .and. yc .ge. 0) then
            ixs = xc
            ixf = xc
            iys = yc
            iyf = yc
            ixs1 = xc
            ixf6 = xc
            if (xrinc_loc .ge. 20) then
               ixs1 = 0  
               ixf6 = nx+1 
            endif
            iys1 = yc
            iyf6 = yc
            if (yinc_loc .ge. 20) then
               iys1 = 0 
               iyf6 = ny+1 
            endif
         endif

c...  Set flag that indicates wide open Jac'n "box" for subroutine bouncon.
      if (xc .lt. 0) then
         openbox = .true.
      elseif (xccuts .or. xcturb) then
         openbox = .true.
      elseif ( (0<=yc).and.(yc<=yinc_loc) ) then # for integrated core flux BC
         openbox = .true.
      else
         openbox = .false.
      endif

c...  Set flags that indicate when the Jac'n "box" overlaps left or right
c...  boundary cells of a mesh region.  Used in subroutine bouncon.
      xcnearlb = .false.
      do jx = 1, nxpt
         xcnearlb = xcnearlb .or.
     .       ( (xc-xlinc_loc.le.ixlb(jx)) .and. (xc+xrinc_loc.ge.ixlb(jx)) )
      enddo
      if (xc .lt. 0) xcnearlb = .true.
      xcnearrb = .false.
      do jx = 1, nxpt
         xcnearrb = xcnearrb .or.
     .      ( (xc-xlinc_loc.le.ixrb(jx)+1) .and. (xc+xrinc_loc.ge.ixrb(jx)) )
      enddo
      if (xc .lt. 0) xcnearrb = .true.

      END SUBROUTINE initialize_ranges


      SUBROUTINE calc_rhs(yldot)
      IMPLICIT NONE
      Use(Selec)
      Use(Dim)
      Use(UEpar)
      Use(Indexes)
      Use(Rhsides)
      Use(Comgeo)
      Use(Ynorm)
      Use(Xpoint_indices)
      Use(Indices_domain_dcl)
      real yldot(*)
      integer iy, ix, ifld, iv, jx, iv1, igsp, iv2

**********************************************************************
*  --  Equations to be solved --
**********************************************************************
      do iy = j2, j5
         do ix = i2, i5
            do ifld = 1, nisp
	       if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1-iseqalg(iv)) *
     .                        resco(ix,iy,ifld)/(vol(ix,iy)*n0(ifld))
               endif
            end do
            do ifld = 1, nusp
	       if(isuponxy(ix,iy,ifld) .eq. 1) then
                  iv = idxu(ix,iy,ifld)
                  yldot(iv) = (1-iseqalg(iv)) *
     .                        resmo(ix,iy,ifld)/(volv(ix,iy)*fnorm(ifld))
                  do jx = 1, nxpt
                     if (ix.eq.ixrb(jx) .and. ixmxbcl.eq.1) yldot(iv) =
     .                        resmo(ix,iy,ifld)/(volv(ix,iy)*fnorm(ifld))
                  enddo
               endif
            end do
            if(isteonxy(ix,iy) == 1) then
              iv =  idxte(ix,iy)
	      yldot(iv) = (1-iseqalg(iv)) *
     .                                 resee(ix,iy)/(vol(ix,iy)*ennorm)
            endif
            if(istionxy(ix,iy) == 1) then
              iv1 = idxti(ix,iy)
	      yldot(iv1) = (1-iseqalg(iv1)) *
     .                                 resei(ix,iy)/(vol(ix,iy)*ennorm)
            endif
            do igsp = 1, ngsp
	      if(isngonxy(ix,iy,igsp).eq.1) then
                iv2 = idxg(ix,iy,igsp)
                yldot(iv2) = (1-iseqalg(iv2)) *
     .                      resng(ix,iy,igsp)/(vol(ix,iy)*n0g(igsp))
              endif
	      if(istgonxy(ix,iy,igsp).eq.1) then
                iv2 = idxtg(ix,iy,igsp)
                yldot(iv2) = (1-iseqalg(iv2)) *
     .                      reseg(ix,iy,igsp)/(vol(ix,iy)*ennorm)
              endif
            end do
          end do
        end do
c ... The factor (1-iseqalg(iv)) above forces yldot=0 for algebraic
c ... equations, except up(nx,,); these yldot are subsequently set in
c ... subroutine bouncon.


c  POTEN calculates the electrostatic potential, and BOUNCON calculates the 
c  equations for the boundaries. For the vodpk solver, the B.C. are ODEs 
c  in time (rate equations).  Both bouncon and poten must be called before
c  the perturbed variables are reset below to get Jacobian correct


      END SUBROUTINE calc_rhs

      SUBROUTINE check_kaboom
      IMPLICIT NONE      
      Use(UEpar)
      character*80 msgjm
      integer nrcv, ierrjm, ijmgetmr
c
c     Check if "k" or "kaboom" has been typed to jump back to the parser
c
      if (((svrpkg.eq.'nksol') .or. (svrpkg.eq.'petsc')) .and. iskaboom.eq.1) then
                              #can only call once - preserves 's' in vodpk
        ierrjm = ijmgetmr(msgjm,80,1,nrcv)
        if (ierrjm .eq. 0) then
          if (msgjm(1:nrcv).eq.'kaboom' .or. msgjm(1:nrcv).eq.'k')then
            call xerrab("")
          endif
        endif
      endif
      END SUBROUTINE check_kaboom


      SUBROUTINE add_timestep(neq, yl, yldot)
      IMPLICIT NONE
      Use(Dim)     # nusp,nisp,ngsp
      Use(UEpar)   # svrpkg,isbcwdt,isnionxy,isuponxy,isteonxy,istionxy,
                   # isngonxy,isphionxy
      Use(Time_dep_nwt)   # nufak,dtreal,ylodt,dtuse
      Use(Indexes) # idxn,idxg,idxu,dxti,idxte,idxphi
      Use(Ynorm)   # isflxvar,isrscalf
      Use(Indices_domain_dcl) # ixmnbcl,ixmxbcl,iymnbcl,iymxbcl
      Use(Compla)  # zi
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      integer ix,iy,igsp,iv,iv1,ifld,j2l,j5l,i2l,i5l,neq
      real time, yl(neqmx),yldot(neq)

c...  Add a real timestep, dtreal, to the nksol equations 
c...  NOTE!! condition yl(neq+1).lt.0 means a call from nksol, not jac_calc

         if (isbcwdt .eq. 0) then  # omit b.c. eqns
cccMER   NOTE: what about internal guard cells (for dnbot,dnull,limiter) ???
            j2l = 1
            j5l = ny
            i2l = 1
            i5l = nx
         else                      # include b.c. eqns
            j2l = (1-iymnbcl)
            j5l = ny+1-(1-iymxbcl)
            i2l = (1-ixmnbcl)
            i5l = nx+1-(1-ixmxbcl)
         endif           
         do iy = j2l, j5l    # if j2l=j2, etc., omit the boundary equations
            do ix = i2l, i5l
              do ifld = 1, nisp
                if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = idxn(ix,iy,ifld)
                  yldot(iv) = (1.-fdtnixy(ix,iy,ifld))*yldot(iv)
                  if(zi(ifld).eq.0. .and. ineudif.eq.3) then
                    yldot(iv) = yldot(iv) - (1/n0(ifld))*
     .                          (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                  else
                    yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
                endif
              enddo
               if(ix.ne.nx+2*isbcwdt) then  
                              # nx test - for algebr. eq. unless isbcwdt=1
                  do ifld = 1, nusp
                    if(isuponxy(ix,iy,ifld).eq.1) then
                      iv = idxu(ix,iy,ifld)
                      yldot(iv) = (1.-fdtupxy(ix,iy,ifld))*yldot(iv)
                      yldot(iv) = yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                    endif
                  enddo
               endif
               if (isteonxy(ix,iy) == 1) then
                 iv =  idxte(ix,iy)
                 yldot(iv) = (1.-fdttexy(ix,iy))*yldot(iv)
                 yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv) 
               endif
               if (istionxy(ix,iy) == 1) then
                 iv1 = idxti(ix,iy)
                 yldot(iv1) = (1.-fdttixy(ix,iy))*yldot(iv1)
                 yldot(iv1)=yldot(iv1) - (yl(iv1)-ylodt(iv1))/dtuse(iv1)
               endif
               do igsp = 1, ngsp
                  if(isngonxy(ix,iy,igsp).eq.1) then
                     iv = idxg(ix,iy,igsp)
                     yldot(iv) = (1.-fdtngxy(ix,iy,igsp))*yldot(iv)
                     if(ineudif.eq.3) then
                       yldot(iv) = yldot(iv) - (1/n0g(igsp))*
     .                            (exp(yl(iv))-exp(ylodt(iv)))/dtuse(iv)
                     else
                       yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                     endif
                  endif
               enddo
               do igsp = 1, ngsp
                  if(istgonxy(ix,iy,igsp).eq.1) then
                     iv = idxtg(ix,iy,igsp)
                     yldot(iv) = (1.-fdttgxy(ix,iy,igsp))*yldot(iv)
                     yldot(iv) =yldot(iv)-(yl(iv)-ylodt(iv))/dtuse(iv)
                  endif
               enddo
               if (isphionxy(ix,iy).eq.1 .and. isbcwdt.eq.1) then
                  iv = idxphi(ix,iy)
                  yldot(iv) = (1.-fdtphixy(ix,iy))*yldot(iv)
                  yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtuse(iv)
               endif

            enddo
         enddo
      
C...  Now do an additional relaxation of the potential equations with
c...  timestep dtphi
        if (dtphi < 1e10) then
          do iy = 0, ny+1
            do ix = 0, nx+1
              if (isphionxy(ix,iy) == 1) then
                iv = idxphi(ix,iy)
                yldot(iv) = yldot(iv) - (yl(iv)-ylodt(iv))/dtphi
              endif
            enddo
          enddo
        endif

      END SUBROUTINE add_timestep

c-----------------------------------------------------------------------
      SUBROUTINE pandf (xc, yc, neq, time, yl, yldot)

c... Calculates matrix A and the right-hand side depending on the values
c... of xc, yc.
c  Definitions for argument list
c
c  Input variables:
c    xc is poloidal index of perturbed variablefor Jacobian calc, 
c       or =-1 for full RHS evaluation
c    yc is radial index for perturbed variable for Jacobian calc, 
c       or =-1 for full RHS evaluation
c    neq is the total number of variables
c    time is the present physical time; useable by VODPK but not NKSOL
c    yl is the vector of unknowns
c  Output variables:
c    yldot is the RHS of ODE solver or RHS=0 for Newton solver (NKSOL)

      IMPLICIT NONE

*  -- input arguments
      integer xc, yc, neq      
      real time, yl(*),yldot(*)
      Use(PandfTiming)
      Use(ParallelEval)
      Use(MCN_sources)
      Use(UEpar)
      Use(Ynorm)
      Use(Timing)
      Use(Selec)
      Use(Coefeq)
      Use(Time_dep_nwt)   # nufak,dtreal,ylodt,dtuse
      real tick, tock, tsfe, tsjf, tsnpg
      external tick, tock

c
c     Check if "k" or "kaboom" has been typed to jump back to the parser
      call check_kaboom

c     check if a "ctrl-c" has been type to interrupt - from basis
      call ruthere

************************************************************************
*  -- initialization --
************************************************************************
************************************************************************
*   This section is to use in the calculation of the jacobian locally.
************************************************************************

c ... Get initial value of system cpu timer.
      if(xc .lt. 0) then
         tsfe = tick()
      else
         tsjf = tick()
      endif
!     Initialize loop ranges based on xc and yc
      call initialize_ranges(xc, yc, xlinc, xrinc, yinc)

      if (xc .ge. 0 .and. yc .ge. 0) then 
          call jacobian_store_momentum(xc, yc)
          call jacobian_store_volsources(xc, yc)
      end if

      if (ismcnon .ne. 0) call initialize_ismcnon(yl(neq+1))

c... Timing of pandf components (added by J. Guterl)
        if (TimingPandf.gt.0
     .      .and. yl(neq+1) .lt. 0 .and. ParallelPandf1.eq.0
     .) then
        TimingPandfOn=1
      else
        TimingPandfOn=0
      endif
      if (TimingPandfOn.gt.0) TimePandf=tick()


************************************************************************
c... First, we convert from the 1-D vector yl to the plasma variables.
************************************************************************
      if (TimingPandfOn.gt.0) 
     .      TimeConvert0=tick()

      call convsr_vo (xc, yc, yl)  # pre 9/30/97 was one call to convsr

      if (TimingPandfOn.gt.0) 
     .      TotTimeConvert0=TotTimeConvert0+tock(TimeConvert0)
      if (TimingPandfOn.gt.0) 
     .      TimeConvert1=tick()

      call convsr_aux (xc, yc)

      if (TimingPandfOn.gt.0) 
     .      TotTimeConvert1=TotTimeConvert1+tock(TimeConvert1)
c...  TODO: gather variables calculated in convert


      call calc_plasma_diffusivities
c...  No gradients to separate out
      call initialize_driftterms


      call calc_driftterms
      if(isphion+isphiofft .eq. 1) call calc_currents
      if(isphion+isphiofft .eq. 1) call calc_fqp

c...  TODO: gather variables calculated in calc driftterms
c...        v2 needed by calc_friction
c...  TODO: Break out conditionals, move to top
      call calc_friction(xc)
************************************************************************
*     Calculate the currents fqx, fqy, fq2 and fqp, if isphion = 1
*     or if isphiofft = 1.
************************************************************************

      call calc_elec_velocities
c...  Add checks on ishosor and ispsorave: parallel only works for == 0
      call calc_volumetric_sources(xc, yc)

*****************************************************************
c In the case of neutral parallel mom, call neudif to get
c flux fngy, vy and uu, now that we have evaluated nuix etc.
*****************************************************************
      if (ineudif .eq. 1) then
         call neudif
      elseif (ineudif .eq. 2) then
c ..Timing
      if(istimingon==1) tsnpg=tick()
         call neudifpg
c ..Timing
      if(istimingon==1) ttnpg = ttnpg + tock(tsnpg)
      elseif (ineudif .eq. 3) then
         call neudifl
      else
         call neudifo
      endif

      call calc_srcmod

      call calc_plasma_viscosities

      call calc_plasma_heatconductivities

      call calc_plasma_equipartition

      call calc_gas_heatconductivities
c ... Call routine to evaluate gas energy fluxes
****************************************************************
      call engbalg

      call calc_plasma_transport
    
      call calc_fniycbo

c----------------------------------------------------------------------c
c          SCALE SOURCE TERMS FROM MONTE-CARLO-NEUTRALS MODEL
c
c     These sources are used in the residuals (resco,resmo,resee,resei)
c     so the call to scale_mcn must occur BEFORE these residuals are
c     evaluated.  Since they scale with fnix at the divertor plates,
c     the call to scale_mcn must occur AFTER fnix has been calculated.


      if (ismcnon .ne. 0) call scale_mcnsor
c----------------------------------------------------------------------c

*********************************************************************
c  Here we do the neutral gas diffusion model
c  The diffusion is flux limited using the thermal flux
**********************************************************************

      call calc_plasma_momentum(xc, yc)

c...  Compute total viscosity for nonuniform B-field; put in visvol_v,q
      if (cfvisxneov+cfvisxneoq > 0.) call upvisneo



      call calc_plasma_energy
      call calc_feeiycbo ! Nothing much to parallelize here, just do serial

      call calc_plasma_particle_residuals
      call calc_gas_continuity_residuals
      call calc_plasma_momentum_residuals()
      call calc_gas_energy_residuals
      call calc_atom_seic
c...  Requires gas energy residuals
      call calc_plasma_energy_residuals(xc, yc)

      call calc_rhs(yldot)

c  POTEN calculates the electrostatic potential, and BOUNCON calculates the 
c  equations for the boundaries. For the vodpk solver, the B.C. are ODEs 
c  in time (rate equations).  Both bouncon and poten must be called before
c  the perturbed variables are reset below to get Jacobian correct

      if (isphion.eq.1) call calc_potential_residuals (neq, yl, yldot)

      call bouncon (neq, yldot)

      if (xc .ge. 0 .and. yc .ge. 0) call jacobian_reset(xc, yc)

c ... Accumulate cpu time spent here.
      if(xc .lt. 0) then
            ttotfe = ttotfe + tock(tsfe)
      else
            ttotjf = ttotjf + tock(tsjf)
      endif
      if (TimingPandfOn.gt.0) 
     .      TotTimePandf=TotTimePandf+tock(TimePandf)


c...  ====================== BEGIN OLD PANDF1 ==========================

c...  If isflxvar=0, we use ni,v,Te,Ti,ng as variables, and the ODEs need
c...  to be modified as original equations are for d(nv)/dt, etc
c...  If isflxvar=2, variables are ni,v,nTe,nTi,ng. Boundary equations and
c...  potential equations are not reordered.

      if(isflxvar.ne.1 .and. isrscalf.eq.1) call rscalf(yl,yldot)

c
c ... Now add psuedo or real timestep for nksol method, but not both
      if (nufak.gt.1.e5 .and. dtreal.lt.1.e-5) then
         call xerrab('***Both 1/nufak and dtreal < 1.e5 - illegal***')
      endif




c...  Add a real timestep, dtreal, to the nksol equations 
c...  NOTE!! condition yl(neq+1).lt.0 means a call from nksol, not jac_calc

      if(dtreal < 1.e15) then
       if((svrpkg=='nksol' .and. yl(neq+1)<0) .or. svrpkg == 'petsc') then
        call add_timestep(neq, yl, yldot)
       endif   #if-test on svrpkg and yl(neq+1)
      endif    #if-test on dtreal



      return
      END SUBROUTINE pandf
c****** end of subroutine pandf ************

c-----------------------------------------------------------------------
      subroutine timimpfj(tick, xc)
      real tock, tick
      integer xc
      external tock
      Use(Timing)   # ttimpfe,ttimpjf

      if (xc .lt. 0) then
         ttimpfe = ttimpfe + tock(tick)
      else
         ttimpjf = ttimpjf + tock(tick)
      endif

      return
      end
c---- end of subroutine timimpfj ---------------------------------------

