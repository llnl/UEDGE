c-----------------------------------------------------------------------
c $Id: odesetup.m,v 7.26 2022/11/30 23:10:58 meyer8 Exp $
c
c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"
c-----------------------------------------------------------------------
      subroutine allocate

*     ALLOCATE is an auxiliary subroutine to allocate the variables in
*     case that the code hits an unrecoverable error that force it to
*     exit with a fatal error status.

      implicit none

      Use(Dim)      # nx,ny,nxm,nhsp,nzsp,nzspt,nusp,ngsp
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,
                               # iysptrx1,iysptrx2,iysptrx
      Use(Math_problem_size)   # neqmx,numvar,numvarbwpad
      Use(Cdv)      # ifexmain, iallcall
      Use(Share)    # nycore,nysol,nxleg,nxcore,nxomit,isgrdsym,
                    # nyomitmx,geometry
      Use(UEpar)    # isnion,isupon,isupgon,isteon,
                    # istion,isngon,isphion
      Use(UEint)    # minu,ziin,newgeo,mhdgeo,isallloc
      Use(Lsode)    # neq,jacflg,jpre,ipar,ismmaxuc,mmaxu
      Use(Solver_work_arrays)   # liw,lrw,iwork
      Use(Interp)   # nxold,nyold
      Use(Decomp)   # ubw,lbw
      Use(Jacobian)     # neqp1,nnzmx
      Use(Jacreorder)   # ireorder
      Use(Preconditioning)     # premeth,lenpfac,lenplufac,lenplumx
      Use(Nonzero_diagonals)   # ndiagmx
      Use(Model_choice)        # iondenseqn
      Use(Jacobian_part)       # nnz1mx
      Use(Jacobian_full)       # isjacmap, jacfull
      Use(Opt_input)           # inopt,iworkin(8),iworkin(24)
      Use(Constraints)
      Use(Time_dep_nwt)        # nufak
      Use(Comtra)              # kyet
      Use(Imprad)              # isimpon
      Use(Reduced_ion_interface)   # misotope,natomic,nchstate
      Use(Ynorm)               # iscolnorm
      Use(Selec)               # xlinc,xrinc,yinc
      Use(Bcond)               # nwsor,igspsori,ispfbcvsix,iswobcvsix
      Use(Jac_work_arrays)     # lwp,liwp
      Use(Coefeq)              # cngmom
cc      Use(Rccoef)
      Use(MCN_dim)             # MCN dimensions
      Use(PNC_data)			   # PNC data storage

* --  local variables
      integer lda, lenk, ngspon, nispon, nuspon, ntgspon, ifld, isor, id
      character*60 runid
      #Former Aux module variables
      integer igsp

*=======================================================================
*//computation//
      id = 1
      call gallot("Stat",0)

c...  Be sure nufak and dtreal are large if this is time-dependent run

c ... Set mesh dimensions and x-point parameters
      if (isallloc .eq. 0) then       # skip if allocating local proc
      if (newgeo .ne. 0) then         # skip if geometry is unchanged
      if (gengrid==1) then            # compute from flx and grd input
         call com_set_dims
         call gchange("Xpoint_indices",0)
c----------------------------------------------------------------------c
         if (mhdgeo == 1) then   #toroidal geo from MHD equilibrium
c----------------------------------------------------------------------c
            if (geometry=="dnull") then #snowflake external
               call set_dnull_indices
            elseif (geometry=="isoleg") then
               call set_isoleg_indices
            else
               iysptrx1(1) = nycore
               iysptrx2(1) = nycore
               iysptrx = nycore
               ixlb(1) = 0
               ixpt1(1) = nxleg(1) + nxxpt
               ixpt2(1) = nxleg(1) + nxcore(1) +
     .                                     nxcore(2) + 3*nxxpt
               ixrb(1) = ixpt2(1) + nxleg(2) + nxxpt
               if (nxomit .gt. 0) then
                  ixpt1(1) = ixpt1(1) - nxomit
                  ixpt2(1) = ixpt2(1) - nxomit
                  ixrb(1) = ixrb(1) - nxomit
               endif
            endif
c----------------------------------------------------------------------c
         elseif (mhdgeo == 2) then   #toroidal geo with circular cross-sec
c----------------------------------------------------------------------c
            iysptrx1(1) = nycore
            iysptrx2(1) = nycore
            iysptrx = nycore
            ixlb(1) = 0
            ixpt1(1) = nxleg(1)
            ixpt2(1) = nxleg(1) + nxcore(1) +
     .                                  nxcore(2)
            ixrb(1) = ixpt2(1) + nxleg(2)
c----------------------------------------------------------------------c
	 else	# cases mhdgeo=0 (cyl), mhdgeo=-1 (slab), mhdgeo=-2 (mirror)
c----------------------------------------------------------------------c
            iysptrx1(1) = nycore
            iysptrx2(1) = nycore
            iysptrx = nycore
            ixlb(1) = 0
            ixpt1(1) = -1
            ixpt2(1) = nxcore(2)
            ixrb(1) = ixpt2(1) + nxleg(2) + nxxpt
            if (isgrdsym .eq. 1) then
               ixpt1(1) = (nxm - ixpt2(1))/2
               ixpt2(1) = (nxm + ixpt2(1))/2
            endif
c----------------------------------------------------------------------c
         endif	# end if-test on mhdgeo
c----------------------------------------------------------------------c
      elseif (gengrid==0) then        # read from gridue file
         if (geometry=="dnull" .or. geometry(1:9)=="snowflake" .or.
     .       geometry=="dnXtarget" .or. geometry=="isoleg") then
            nxpt=2
         else
            nxpt=1
            call gchange("Xpoint_indices",0)  #needed to allocate ixrb
	    if (ixrb(1) == 0) then  #calc indices if not in gridue file
               iysptrx1(1) = nycore
               iysptrx2(1) = nycore
               iysptrx = nycore
               ixlb(1) = 0
               ixpt1(1) = nxleg(1) + nxxpt
               ixpt2(1) = nxleg(1) + nxcore(1) +
     .                                     nxcore(2) + 3*nxxpt
               ixrb(1) = ixpt2(1) + nxleg(2) + nxxpt
               if (nxomit .gt. 0) then
                  ixpt1(1) = ixpt1(1) - nxomit
                  ixpt2(1) = ixpt2(1) - nxomit
                  ixrb(1) = ixrb(1) - nxomit
               endif
            endif
         endif
         call gchange("Xpoint_indices",0)
         call readgridpars(trim(GridFileName),runid)  #define/redefine iysptrx1 etc
         nx = nxm - abs(nxomit)
         ny = nym - nyomitmx
      endif	# end if-test on gengrid
      endif	# end if-test on newgeo
      endif	# end if-test on isallloc


c ... Check that number neutral gas species is not larger than ngspmx
      if (ngsp .gt. ngspmx) then
         call remark('***Increase ngspmx (now <ngsp) & recompile***')
         call remark('*** And increase nuixold,psorgold in pandf ***')
         call xerrab("")
      endif

c ... Check that isupcore=2 is not being used
      do igsp = 1, ngsp
        if (isngcore(igsp) == 2) then
          write(*,*) "*** isngcore=2 option unvailable; igsp = ",igsp
          call xerrab("")
        endif
      enddo


c ... Check that a gas source and albedo is not assigned to nonexistent gas sp
      do isor = 1, nwsor
         if (igspsori(isor).gt. ngsp .or. igspsoro(isor).gt.ngsp) then
            call remark('*** igspsori,o refers to igsp > ngsp *** ')
            call xerrab("")
         endif
      enddo

c ... Check consistency of cngmom; should be zero if inertial neutrals
      if (isupgon(1).eq. 1 .and. cngmom(1).ne.0) then
         call remark('*** WARNING, likely Error: cngmom=1, isupgon=1')
ccc         call xerrab("")
      endif

c ... Check consistency of radial gradient boundary conditions
c ... Note: only ix=1 is checked for "is..." flag as they are now function of ix
      if ((isnwcono(1)==3 .or. isnwconi(1)==3) .and.
     .                                       lyni(1)+lyni(2)>1e9) then
        call remark('*** WARNING: large lyni does not give 0 flux B.C.')
      endif
      if (nx > 1000) then
        call xerrab('*** lyni allocate req nx<=1000; need recompile')
      endif
      if ((istewc==3 .or. istepfc==3) .and. lyte(1)+lyte(2) > 1e9) then
        call remark('*** WARNING: large lyte may not give 0 flux B.C.')
      endif
      if ((istiwc==3 .or. istipfc==3) .and. lyti(1)+lyti(2) > 1e9) then
        call remark('*** WARNING: large lyti may not give 0 flux B.C.')
      endif

c ... Check if old neutral energy-loss factor is used; copy and warn
      if (cgpl > 0.) then
        call remark('*** ERROR: use cgengpl, not cgpl for neut eng loss')
      endif

c ... Check consistency of impurity variables.
      nzspt = 0
      do ifld = 1, ngsp-nhgsp
         nzspt = nzspt + nzsp(ifld)
         if (nzsp(ifld) .gt. nzspmx) then
            call xerrab('*** nzspmx .gt. nzsp(ifld); reset nzspmx ***')
         endif
      enddo
      if (isimpon .ge. 4) then
         if (nzspt .lt. 1)
     .      call remark('*** Warning:  nzspt<1 for isimpon>3')
      elseif (isimpon .gt. 2) then
         if (nzspt .le. 0)
     .      call xerrab('***Error: nzspt must be >0 for isimpon=3 or 4')
         if (isnion(1) .eq. 0)
     .      call xerrab('***Error: isnion must be 1 for isimpon=3 or 4')
      else
         if (nzspt .ne. 0)
     .     call remark('***Warning: nzspt set to 0 because isimpon<3')
         nzspt = 0
      endif

c ... Check if yinc=2 for isphion=1
      if (isphion == 1 .and. yinc .ne. 2) then
	 call remark('*** Warning: yinc=2 recommended when isphion=1')
      endif

c ... Check if isnfmiy=1 when geometry is snowflake > SF15
      if (geometry=="snowflake45" .or. geometry=="snowflake75" .or.
     .    geometry=="snowflake105" .or. geometry=="snowflake135" .or.
     .    geometry=="snowflake165") then
         if (isnfmiy == 1) then
            call xerrab('*** ERROR: isnfmiy=1 not option here; set to 0')
         endif
      endif

c ... Calculate variables related to size of the problem.
      nisp = nhsp + nzspt
      if (nisp .gt. nispmx)
     .   call xerrab('**Error:  variables in Comtra limit nisp < nispmx')
      nusp = nhsp
      if (isimpon .eq. 4) nusp = nisp
      if (isimpon .eq. 5 .or. isimpon .eq. 6 .or. isimpon .eq. 7) then
         nusp = 1 + isupgon(1) + nusp_imp
         call mombal0 (nisp,nhsp,nzsp,minu,ziin,misotope,natomic,nchstate)
         call inicon
         call initmombal (misotope,natomic,nchstate)
      endif
      ngspon = 0
      ntgspon = 0
      do igsp = 1, ngsp
         ngspon = ngspon + isngon(igsp)
         ntgspon = ntgspon + istgon(igsp)
      enddo
      nispon = 0
      do ifld = 1, nisp
         nispon = nispon + isnion(ifld)
      enddo
      nuspon = 0
      do ifld = 1, nusp
         nuspon = nuspon + isupon(ifld)
      enddo
      numvar = isteon + istion + nispon + nuspon
     .                         + ngspon + ntgspon + isphion
      
c...  Allocate arrays isnionxy, isuponxy, etc in prep. for neq calc
      call gchange("UEpar",0)    # allocates isnionxy,isnioffxy,etc

cmer  The following calculations must be done in a separate subroutine
cmer  after the UEpar arrays have been allocated.
c...  Compute which equations are actually evolved from input parameters
c...  Calculate actual number of equations (some may be turned-off)
      call setonxy   #sets neq = acutal number eqns; neqmx adds 2 - usually

c...  Set PF and main-chamber boundary-condition arrays that either be uniform
c...  in the ix poloidal direction or varying by user input. The two walls are
c...  treated independtly with flags bbb.ispfbcvsix and bbb.iswobcvsix where 0
c...  yields uniform BCs versus ix and 1 utilizes array values set by user
      call setwallbcarrays
      
      neqmx = neq + 2  #1st extra flag for Jac. calc; 2nd nufak=/1psuedo dt

      neqp1 = neq + 1
      ipar(1) = neq
      ubw = (numvar+numvarbwpad)*(nx+ixpt2(nxpt)-max(0,ixpt1(1))+4)
      lbw = (numvar+numvarbwpad)*(nx+ixpt2(nxpt)-max(0,ixpt1(1))+4)
      if (isphion*isnewpot == 1) then  #phi eqn  4th order in y; larger bandw
         ubw = 4*numvar*nx
	 lbw = 4*numvar*nx
      endif
      if (xlinc.gt.20 .and. xrinc.gt.20 .and. yinc.gt.20) then
         ubw = neq       # use full matrix & search full range for Jacobian
         lbw = neq       # use full matrix & search full range for Jacobian
      endif
      lda = 2*lbw + ubw + 1
      if (jpre+jacflg+1 .eq. 0) lda = 1
c...  Increase of mmaxu with problem size is empirical
      if (ismmaxuc .eq. 1) mmaxu = neq**0.5
      if (ismmaxuc.eq.0 .and. icntnunk.eq.1) then
         call remark('WARNING: ismmaxuc=0 & icntnunk=1; Jac storage ok?')
      endif

c ... Set maximum lengths for preconditioner arrays, and allocate
c     memory for auxiliary arrays for desired storage-format option,
c     if necessary.
      if (lenpfac .lt. 9*numvar) then
        nnzmx = 9*numvar*neq
      else
        nnzmx = lenpfac*neq
      endif

c...  Init OMP/MPI/Hybrid variables after assignment of nnzmx,
c...   which is used in InitOMP/InitMPI/InitHybrid (added by  J.Guterl)
c!omp call InitParallel

      if (premeth .eq. 'ilut') then
         lwp = nnzmx + lenplufac*neq   # extra space for fill in
         lenplumx = lwp
         liwp = max(lwp+neq, nnzmx+neq+1)   # see jac_lu_decomp
c save preconditioner reordering info for icntnunk=1 case.
         if(ireorder .eq. 1) call gchange("Jacreorder",0)
      elseif (premeth .eq. 'banded') then
         lwp = lda*neq
         liwp = 3+neq
      else
         write(STDOUT,*) "*** Invalid option:  premeth=",premeth
         call xerrab("")
      endif

c ... Set lengths of work arrays for desired solver option.
cc     sizes are set for mf=1 case; opt. input iwork(1)=mmaxu needed
cc     do not reset if icntnunk=1 (saved Jac) as storage changes if mmaxu reset
         lenk = 6 + 4*neq + (neq+1)*mmaxu + 2*mmaxu*(mmaxu+1)
         lrw = 4 + 4*neq + lenk + lwp
         liw = 20 + mmaxu + liwp

c ... Set maximum lengths of arrays to hold part of Jacobian matrix,
c     and allocate space, if used.
      if (iondenseqn .eq. 'inel' .or.
     .    isimpon .eq. 3 .or. isimpon .eq. 4) then
         call xerrab("*** Ave-ion models isimpon=3,4 disabled")
ccc         nnz1mx = lenpfac*(nx+2)*(ny+2)
ccc         if (isimpon .eq. 4) nnz1mx = nzspt * nnz1mx
ccc         call gallot("Jacobian_part",0)
      endif

*  -- Allocation of the different common blocks
      if (newgeo .eq. 1) then
        if (manualgrid == 0) then
         call gallot("Comgeo",0)
         call gallot("Noggeo",0)
         call gallot("RZ_grid_info",0)
         call gallot("RZ_cell_info",0)
         call gallot("Bfield",0)
         call gallot("Linkbbb",0)
        endif
      endif
      call gchange("Bcond",0)
cc      if(iscallrccoef==1)      call gallot("Rccoef",0)
cc      call gallot("Rccoef",0)
      call gchange("Rccoef",0)
      call gchange("Outpwall",0) 
      call gchange("Timary",0)
      call gchange("Compla",0)
      call gchange("Interprettrans",0)
      call gchange("Comflo",0)
      call gchange("Cfric",0)
      call gchange("Comtra",0)
      if (kyet .gt. 1.e-20) then
         call xerrab
     .      ('kyet must = 0; turb. transport in SOL is disabled.')
         call gallot("Turbulence_diagnostics",0)
      endif
      call gchange("Poten",0)
      call gchange("Lsode",0)   # changed from gallot to switch jpre,.
      call gallot("Constraints",0)
      call gallot("Time_dep_nwt",0)
      call gchange("Ynorm",0) # preserves sfscal values for icntnunk=1
      call gchange("Selec",0) # preserves ixm1 & ixp1 values for newgeo=0
      call gallot("Indexes",0)
      call gchange("Oldpla",0)
      call gallot("Rhsides",0)
      call gchange("MCN_sources",0)
      call gallot("Jacobian_restore",0)
      call gchange("Conduc",0)   # preserves nuiz, eeli, etc for icnuiz=2
      call gallot("Locflux",0)
      call gallot("Gradients",0)
      call gallot("Condition_number",0)
      call gchange("Jacobian",0) # preserves Jacobian values for icntnunk=1
      call gchange("Jacobian_csc",0) # preserves Jacobian values for icntnunk=1
      call gchange("Jacaux",0) # preserves preconditioner data for icntnunk=1
      call gallot("Newtaux",0)
      call gallot("Wkspace",0)
      call gallot("Postproc",0)
      if (isimpon .gt. 0) call gchange("Imprad",0)
      if (isimpon .gt. 2) then
         call gallot("Impurity_source_flux",0)
         call gchange("Impurity_source",0)
         call gchange("Sources_at_walls",0)
      endif
      call gchange("Volsrc",0)
c preserves preconditioner data for icntnunk=1
      call gchange("Solver_work_arrays",0)

      call gallot("Jac_work_arrays",0)
      call gallot("Temporary_work_arrays",0)


*  -- create memory space for interpolants for the grid sequencing. --
         if (ifexmain.eq.0 .or. iallcall.eq.0) then
            nxold = nx
            nyold = ny
            call gchange("Interp",0)
         endif

c	  Allocate PNC_data group
      call gchange("PNC_data",0)

      iallcall = 1 # indicates allocate called at least once; nis allocated
      return
      end
******* end of subroutine allocate ****
***************************************
c----------------------------------------------------------------------c

      subroutine setwallbcarrays

c...  Set boundary condition flags/arrays on private-flux wall depending on
c...  input flag ispfbcvsix, where 0 gives uniform bndry conds and 1 leaves
c...  arrays as set by the user after allocation

      implicit none

      Use(Dim)      # nx,ny,nisp,ngsp
      Use(Bcond)    # wall BC arrays
#      Use(UEpar)    # isnion,isupon,isngon,isteon,istion,isphion
                    # isnionxy,isuponxy,isngonxy,isteonxy,istionxy,isphionxy
                    # isnioffxy,isupoffxy,isngoffxy,isteoffxy,istioffxy,isphioffxy
#      Use(Math_problem_size)   # neqmx
#      Use(Lsode)    # neq

      integer ix,ifld,igsp
      
c...  Need to allocate BC arrays with proper dimensions
      call gchange("Bcond",0)

      if(ispfbcvsix == 0) then  #inner private-flux wall

        do ix = 0, nx+1
          iphibcwiix(ix) = iphibcwi
          istepfcix(ix) = istepfc
          istipfcix(ix) = istipfc
          lyteix(1,ix) = lyte(1)
          lytiix(1,ix) = lyti(1)
          lyphiix(1,ix) = lyphi(1)
          do ifld = 1, nisp
            isnwconiix(ix,ifld) = isnwconi(ifld)
            isupwiix(ix,ifld) = isupwi(ifld)
            lyniix(1,ix,ifld) = lyni(1)
          enddo       
          do ifld = 1, nisp
            isupwiix(ix,ifld) = isupwi(ifld)
            lyupix(1,ix,ifld) =lyup(1)
          enddo
          do igsp = 1, ngsp
            istgpfcix(ix,igsp) = istgpfc(igsp)
          enddo  
        enddo
        
      endif   #test on ispfbcvsix for private-flux wall

c ..................
      if(iswobcvsix == 0) then  #outer main-chamber wall

        do ix = 0, nx+1
          iphibcwoix(ix) = iphibcwo
          istewcix(ix) = istewc
          istiwcix(ix) = istiwc
          lyteix(2,ix) = lyte(2)
          lytiix(2,ix) = lyti(2)
          lyphiix(2,ix) = lyphi(2)
          do ifld = 1, nisp
            lyniix(2,ix,ifld) = lyni(2)
            isnwconoix(ix,ifld) = isnwcono(ifld)
            isupwoix(ix,ifld) = isupwo(ifld)
          enddo
          do ifld = 1, nisp
            isupwoix(ix,ifld) = isupwo(ifld)
            lyupix(2,ix,ifld) =lyup(2)
          enddo
          do igsp = 1, ngsp
            istgwcix(ix,igsp) = istgwc(igsp)
          enddo
        enddo
        
      endif   #test on iswobcvsix for main-chamber wall
          
      return
      end
c----------------------------------------------------------------------c
      
      subroutine setonxy

c...  Compute which equations are evolved based on input vars isnionffxy etc.
c...  Then calculate total-number-of-eqns = neq (as some may be turned-off)

      implicit none

      Use(Dim)      # nx,ny,nisp,ngsp
      Use(UEpar)    # isnion,isupon,isngon,isteon,istion,isphion
                    # isnionxy,isuponxy,isngonxy,isteonxy,istionxy,isphionxy
                    # isnioffxy,isupoffxy,isngoffxy,isteoffxy,istioffxy,isphioffxy
      Use(Math_problem_size)   # neqmx
      Use(Lsode)    # neq

      integer ix,iy,ifld,igsp

      do ifld = 1, nisp
	do iy = 0, ny+1
	  do ix = 0, nx+1
	    isnionxy(ix,iy,ifld) = isnion(ifld)*(1-isnioffxy(ix,iy,ifld))
	    isuponxy(ix,iy,ifld) = isupon(ifld)*(1-isupoffxy(ix,iy,ifld))
	  enddo
        enddo
      enddo
      do igsp = 1, ngsp
	do iy = 0, ny+1
	  do ix = 0, nx+1
	    isngonxy(ix,iy,igsp) = isngon(igsp)*(1-isngoffxy(ix,iy,igsp))
	    istgonxy(ix,iy,igsp) = istgon(igsp)*(1-istgoffxy(ix,iy,igsp))
	  enddo
        enddo
      enddo
      do iy = 0, ny+1
        do ix = 0, nx+1
	  isteonxy(ix,iy) = isteon*(1-isteoffxy(ix,iy))
	  istionxy(ix,iy) = istion*(1-istioffxy(ix,iy))
	  isphionxy(ix,iy) = isphion*(1-isphioffxy(ix,iy))
	enddo
      enddo

      neq = 0
      do iy = 0, ny+1   #iymb,xbcl =1(0) if guard cell(no)
	do ix = 0, nx+1 #ixmb,xbcl =1(0) if guard cell(no)
	  do ifld = 1, nisp
	    neq = neq + isnionxy(ix,iy,ifld)
          enddo
	  do ifld = 1, nusp
	    neq = neq + isuponxy(ix,iy,ifld)
	  enddo
	  neq = neq+isteonxy(ix,iy)+istionxy(ix,iy)+isphionxy(ix,iy)
	  do igsp = 1, ngsp
	    neq = neq + isngonxy(ix,iy,igsp) + istgonxy(ix,iy,igsp)
	  enddo
	enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine set_dnull_indices
c     Define characteristic indices for a full double-null configuration.

      implicit none

      Use(Dim)            # nxpt
      Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb,
                          # iysptrx1,iysptrx2,iysptrx
      Use(Share)          # nycore,nxleg,nxcore,nxxpt

c     For up/down symmetric double-null:
      iysptrx1(1) = nycore
      iysptrx2(1) = nycore
      iysptrx = nycore
      iysptrx1(2) = iysptrx2(1)
      iysptrx2(2) = iysptrx1(1)
      ixlb(1)  = 0
      ixpt1(1) = nxleg(1) + nxxpt
      ixmdp(1) = ixpt1(1) + nxcore(1) - 1 + nxxpt
      ixpt2(1) = ixmdp(1) + nxcore(1) - 1 + nxxpt
      ixrb(1)  = ixpt2(1) + nxleg(1) + nxxpt
      ixlb(2)  = ixrb(1) + 2
      ixpt1(2) = ixlb(2) + nxleg(2) + nxxpt
      ixmdp(2) = ixpt1(2) + nxcore(2) - 1 + nxxpt
      ixpt2(2) = ixmdp(2) + nxcore(2) - 1 + nxxpt
      ixrb(2)  = ixpt2(2) + nxleg(2) + nxxpt

cccMER 23 Nov 1999
cccMER Although the above is not correct for general un-balanced double-
cccMER null configurations, these indices are subsequently re-set by reading
cccMER them from the 'gridue' file in subroutine nphygeo when gengrid=0.

      return
      end
c-----------------------------------------------------------------------
c----------------------------------------------------------------------c

      subroutine set_isoleg_indices
c     Define characteristic indices for a full double-null configuration.

      implicit none

      Use(Dim)            # nxpt
      Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb,
                          # iysptrx1,iysptrx2,iysptrx
      Use(Share)          # nycore,nxleg,nxcore,nxxpt

c     For mimicing an isolated X-pt:

c...  Get space for two indices
cc      nxpt = 2
cc      call gchange("Xpoint_indices",0)
cc      write(*,*) "After gchange; iysptrx1=", iysptrx1
c...  Set indices
      iysptrx1(1) = nycore
      iysptrx2(1) = nycore
      iysptrx = nycore
      iysptrx1(2) = iysptrx2(1)
      iysptrx2(2) = iysptrx1(1)
      ixlb(1)  = 0
      ixpt1(1) = nxleg(1) + nxxpt
cc      ixmdp(1) = ixpt1(1) + nxcore(1) - 1 + nxxpt
      ixpt2(1) = ixpt1(1) + nxcore(1) - 1 + nxxpt
      ixrb(1)  = ixpt2(1)
      ixlb(2)  = ixrb(1) + 2
      ixpt1(2) = ixlb(2)
cc      ixmdp(2) = ixpt1(2) + nxcore(2) - 1 + nxxpt
      ixpt2(2) = ixlb(2) + nxcore(2) - 1 + nxxpt
      ixrb(2)  = ixpt2(2) + nxleg(2) + nxxpt

      return
      end
c-----------------------------------------------------------------------
      subroutine ueinit

*     UEINIT defines the geometry, initializes the state of the
*     plasma, and defines miscellaneous other fields.

      implicit none

      Use(Dim)      # nx,ny,nhsp,nusp,nzspt,nisp,ngsp,nxpt
      Use(Share)    # nxomit,geometry,isnonog,nyomitmx
                    # nzdf,mcfilename,coronalimpfname
      Use(Multicharge)  # rtnt,rtnn,rtnsd
      Use(Comgeo)   # vol,gx,gy,dx,dy,xnrm,xvnrm,ynrm,yvnrm,sx,sy,rr,
                    # xcs,xfs,xcwi,xcwo,yyc,yyf
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
      Use(Math_problem_size)   # neqmx,numvar
      Use(UEint)    # newgeo,newaph,restart,initsol,ttbeg,
                    # tinit,tscal,ngscal,xgscal,minu,ziin,ixgb
      Use(Grid)     # ig
      Use(Compla)   # mi,zi,mg,znucl
      Use(Comflo)   # fqp,fq2,fqx,fqy
      Use(Cfric)    # frice,frici
      Use(Selec)    # ixm1,ixp1,stretcx
      Use(Phyvar)   # mp,ev
      Use(Comtra)   # kyet
      Use(Turbulence)   # isturbcons,diffusrange,diffuslimit,diffuswgts
      Use(Interp)   # isnintp,nxold,nyold,
                    # ixlbo,ixpt1o,ixpt2o,ixrbo,
                    # xnrmo,xvnrmo,xnrmox,xvnrmox,
                    # ynrmo,yvnrmo,ynrmox,yvnrmox,ixmg,iyomg,ixvmg,
                    # iyvomg,ix2g,iy2g,ixv2g,iyv2g,
                    # nis,tes,tis,phis,ups,ngs,afracs,isimesh
      Use(UEpar)    # ngbackg,,isnion,isupon,isteon,istion,isngon,isphion
                    # isnionxy,isuponxy,isteonxy,istionxy,isngonxy,isphionxy
                    # isphiofft,methg
      Use(Coefeq)   # cngtgx,cngtgy,sxgpr,xstscal
      Use(Lsode)    # neq,ires,ipar,rpar,yl,yldot,srtolpk,rtolv,icntnunk
      Use(Ynorm)    # iscolnorm,norm_cons,floor_cons,var_scale_floor,
                    # n0,temp0,nnorm,ennorm,sigbar0,vpnorm,fnorm,
                    # suscal,sfscal
      Use(Poten)    # sigma1,cfsigm
      Use(Indexes)  # isvarup
      Use(Bfield)   # b0,b0old,btot,rbfbt,rbfbt2
      Use(Bcond)    # tewalli,tiwalli,tewallo,tiwallo,tedge,isfixlb
                    # matt,matp,cion,cizb,crmb
      Use(Imprad)   # isimpon,ismctab
      Use(Impurity_source_flux)   # fnzysi,fnzyso
      Use(Variable_perturbation)  # del
      Use(Time_dep_nwt)           # nufak,inufaknk
      Use(Conduc)   # nuvl
      Use(Locflux)  # fgtdx, fgtdy
      Use(Jacobian) # isjacstnlon
      Use(Rccoef)   # feixextlb,rb;feiyexti,o
      Use(Rhsides)  # psorc, psorxr, msor, msorxr
      Use(Jacobian_restore) # psorcold, etc
      Use(Cut_indices)	# ixcut1,iycut1,ixcut2,iycut2,ixcut3,iycut3
                        # ixcut4,iycut4
      Use(Gradients) #eymask1d


*  -- external routines --
      real ssmin, s2min

*  -- local scalars --
      integer i, iu, ir, irstart, irend
      integer ifld
      integer impflag
      real crni, ctemp, cj, zn_last, proffacx, proffacy, proffacv
      integer ixmp4, jx, jy
      real diffustotal, factor
      integer ifld_fcs, ifld_lcs, igsp_lcs, jz
      #Former Aux module variables
      integer ix,iy,igsp
      real tv

*=======================================================================
*//computation//

*  ---------------------------------------------------------------------
*  preliminaries.
*  ---------------------------------------------------------------------
*  -- check nx, ny, nhsp --
      if (nx.lt.1 .or. ny.lt.1 .or. nhsp.lt.1) then
         call xerrab ('ueinit -- faulty argument nx, ny, nhsp')
      end if

*  -- Initialize some switches:
ccc      if (ngsp .eq. 1) then
ccc	 cngtgx = 0.
ccc	 cngtgy = 0.
ccc      endif

*     ------------------------------------------------------------------
*     initialize the geometry.
*     ------------------------------------------------------------------

*     The problem-dependent physics routine phygeo is called to
*     specify the geometry and the grid metric.
*     phygeo defines xcs(0:nx+1), xfs(0:nx+1), yyc(0:ny+1), yyf(0:ny+1) and
*     the (0:nx+1, 0:ny+1) subblocks of each of the arrays vol, gx, gy, sx, sy,
*     rr. xcs(ix) and yyc(iy) will contain the x- and y-coordinates of
*     the center of the (ix,iy) primary cell, for (ix,iy) in (0:nx+1,
*     0:ny+1).
*     The coordinates of the x-faces will be in xfs(0:nx+1) and those of
*     the y-faces in yyf(0:ny).
*     The volume of the (ix,iy) primary cell will be in vol(ix,iy).
*     The x-diameter and the y-diameter of that cell will be given
*     by 1/gx(ix,iy) and 1/gy(ix,iy).
*     The area of the east surface of the (ix,iy) primary cell will
*     be in sx(ix,iy), for (ix,iy) in (0:nx,0:ny+1).
*     sx(nx+1,0:ny+1) will hold 0.
*     The area of the north surface of the (ix,iy) primary cell will
*     be in sy(ix,iy), for (ix,iy) in (0:nx+1,0:ny).
*     sy(0:nx+1,ny+1) will hold 0.
*     The pitch of the field line at the center of the (ix,iy) cell
*     will be in rr(ix,iy), for (ix,iy) in (0:nx+1,0:ny+1).
*     The pitch is the ratio poloidal/parallel field.

c ... Set flag carried in yl(neq+1) to -1. meaning no pseudo time-step
c ... in the equations. Gets reset to 1.for Jac. calc. to use pseudo nufak
c ... The inverse pseudo time-step, nufak, is stored in yl(neq+2)
         yl(neq+1) = -1.
         if (inufaknk .eq. 1) then   # deter. if nufak is used in Krylov step
            yl(neq+2) = nufak
         else
            yl(neq+2) = 0.
         endif

c...  If geometry=dnbot and isudsym=0, reset isudsym = isupdown_sym flag
      if (geometry == 'dnbot' .and. isudsym == 0) then
         isudsym = 1
         call remark('*** SETTING isudsym=1 BECAUSE geometry=dndot ***')
      endif

c...  Advise to use methg=66 if isnonog=1
      if (isnonog .eq. 1 .and. methg.ne.66) then
         call remark('***********************************************')
         call remark('** CAUTION: NOT USING METHG=66 FOR ISNONOG=1 **')
         call remark('***********************************************')
      endif

c...  Check isphiofft for consistent value
      if (isphiofft.ne.0 .and. isphion.ne.0) then
         call xerrab('*** isphion.ne.0 while isphiofft.ne.0; illegal')
      endif

c...  Check for validity of the masses and atomic numbers.
c...  (mi is provided in units of mp)
      if (ssmin(nisp, mi, 1) .lt. 0) then
         call xerrab ('ueinit -- faulty input mi')
      endif
c     call sscal (nisp, mp, mi, 1)
      zn_last = 0.
      igsp = 1
      iigsp = 1
      do 1 ifld = 1, nhsp    # determine mi, zi, iigsp, and mg for hydrogen
	 mi(ifld) = minu(ifld)*mp
         zi(ifld) = ziin(ifld)
         znucl(ifld) = znuclin(ifld)
c ------------- New lines for parallel neutral momentum eq.
         if (zi(ifld).eq.0 .and. isupgon(1).eq.1) then
            iigsp = ifld
            if (mi(ifld-1) .ne. mi(ifld)) then
             call remark('**Warning: hyd ion & gas masses do not match')
            endif
         endif
         if (abs(znucl(ifld)-zn_last) .gt. 1.e-20 .and.
     .                                           zi(ifld).ne.0) then
            mg(igsp) = facmg(igsp)*mi(ifld)
            if (igsp.gt.ngsp) call xerrab ('uenit -- faulty input mg:h')
            igsp = igsp + 1
         endif
         zn_last = znucl(ifld)
         if (ifld.eq.1 .and. ishymol.eq.1) then   # for hydrogen molecules
            mg(2) = 2*mi(1)
         endif
    1 continue

      ifld_lcs = nhsp
      igsp_lcs = nhgsp
      do jz = 1, ngspmx-1    # determine mi, zi, and mg for impurities
         if (nzsp(jz)==0) break
         ifld_fcs = ifld_lcs + 1
         ifld_lcs = ifld_fcs + nzsp(jz) - 1
         do ifld = ifld_fcs, ifld_lcs
            mi(ifld) = minu(ifld)*mp
            zi(ifld) = ziin(ifld)
            znucl(ifld) = znuclin(ifld)
         enddo
         if (ngsp .gt. nhgsp) then  # impurity gas is present
            igsp_lcs = igsp_lcs + 1
            mg(igsp_lcs) = mi(ifld_lcs)
            if (igsp_lcs.gt.ngsp) call xerrab ('uenit:faulty input mg:z')
         endif
      enddo

      if (ssmin(nisp, zi, 1) .lt. 0) then
         call xerrab ('ueinit -- faulty input zi')
      endif

c...  Check for validity of nibeg and ttbeg.  Scaling for density
c...  and temperature.
c...  (ttbeg is provided in units of ev)
      if (ssmin(nisp, nibeg, 1) .lt. 0) then
         call xerrab ('ueinit -- faulty input nibeg')
      endif
      if (tinit .le. 0) then
         call xerrab ('ueinit -- faulty input ttbeg')
      endif
      ttbeg = tinit * ev

c ... Set default implicit equation-scaling array (for call sfsetnk?).
c ... This call permanently commented out.  Not needed and conflicts
c ... when performing continuation calls with NKSOL (icntnunk = 1).
c      call sfill (neq, 1., sfscal(1), 1)

c ... Set initial values of time-step arrays 
        call s2fill (nx+2, ny+2, 1.e20, dtoptx, 1, nx+2)
        call sfill (neq, 1.e20, dtoptv(:), 1)
        call sfill (neq, 1.e20, dtuse(:), 1)

c ... Set normalization constants for the yl variables seen by solvers.
c     For implicit scaling (iscolnorm .ne. 0), these settings are
c     temporary (see next 2 occurrences of iscolnorm).
      nnorm = n0(1)
      do 117 ifld = 1, nusp
         fnorm(ifld) = n0(ifld)*sqrt(mi(ifld)*temp0*ev)
  117 continue
      ennorm = 1.5*n0(1)*temp0*ev
      if (iscolnorm.eq.1 .or. iscolnorm.eq.2) then
         if (isflxvar .ne. 1) then   # ennorm,fnorm used to build floor_cons
            if (isflxvar .eq. 0) ennorm = ennorm / n0(1)
            do ifld = 1, nusp
               fnorm(ifld) = fnorm(ifld) / n0(ifld)
            enddo
         endif
      endif
      vpnorm = sqrt(temp0*ev/mi(1))
      sigbar0 = cfsigm * sigma1 * temp0**1.5

c ... Initialize array identifying up velocity variables
      do i = 1, numvar
         isvarup(i) = 0
      enddo

c ... Pack normalization constants.
      i = 0
      do ifld = 1, nisp
         if (isnion(ifld) .eq. 1) then
            i = i + 1
            norm_cons(i) = n0(ifld)
         endif
      enddo
      do ifld = 1, nusp
         if (isupon(ifld) .eq. 1) then
            i = i + 1
            norm_cons(i) = fnorm(ifld)
            isvarup(i) = 1
         endif
      enddo
      if (isteon .eq. 1) then
         i = i + 1
         norm_cons(i) = ennorm
      endif
      if (istion .eq. 1) then
         i = i + 1
         norm_cons(i) = ennorm
      endif
      do igsp = 1, ngsp
         if (isngon(igsp) .eq. 1) then
            i = i + 1
            norm_cons(i) = n0g(igsp)
         endif
	 if (istgon(igsp) .eq. 1) then
            i = i + 1
            norm_cons(i) = ennorm
         endif
      enddo
      if (isphion .eq. 1) then
         i = i + 1
         norm_cons(i) = temp0
         isvarphi(i) = 1
      endif

c ... Set floor constants.
      do i = 1, numvar
         floor_cons(i) = var_scale_floor * norm_cons(i)
         if (isvarup(i) .eq. 1) floor_cons(i) = vsf_up*norm_cons(i)
         if (isvarphi(i) .eq. 1) floor_cons(i) = vsf_phi*norm_cons(i)
                        # up velocity variables use full norm_cons
      enddo

c ... For implicit variable scaling, set some normalization constants
c     to unity.  New n0 scaling for norm_cons & in convert. may effect
c     this (3/26/96)?
      if (iscolnorm.eq.1 .or. iscolnorm.eq.2) then
         nnorm = 1.
         ennorm = 1.
         do ifld = 1, nusp
            fnorm(ifld) = 1.
         enddo
      endif

c...  Set net variable perturbation for vodpk finite-diff deriv. to del
c...  If using FORTHON (Python), use delpy to set del as del is special word
      if (delpy > 0.) del = delpy
      srtolpk = del / rtolv

c...  Set boundary conditions for Te,i on walls if arrays are zero
      if (tewalli(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tewalli(0:),1)
      if (tiwalli(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tiwalli(0:),1)
      if (tewallo(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tewallo(0:),1)
      if (tiwallo(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tiwallo(0:),1)

c...  Initialize external eng fluxes to 0 if ext_flags=0 as used anyway
      if (isextpltmod == 0) then
        do jx = 1, nxpt
          call sfill(ny+2,0.,feixextlb(0:,jx),1)
          call sfill(ny+2,0.,feixextrb(0:,jx),1)
        enddo
      endif
      if (isextwallmod == 0) then
        call sfill(nx+2,0.,feiyexti(0:),1)
        call sfill(nx+2,0.,feiyexto(0:),1)
      endif

c...  Set values of sheath potential/Te to 3.0 if values are zero
      do jx=1,nxpt
        if (kappal(ny+1,jx).lt.1.e-10) call sfill (ny+2,3.0,kappal(0:,jx),1)
        if (kappar(ny+1,jx).lt.1.e-10) call sfill (ny+2,3.0,kappar(0:,jx),1)
      enddo

c ... Initialize Multicharge rate table dimensions
      rtnt=0
      rtnn=0
      rtnsd=0
c ... Set up tables for hydrogenic atomic-physics processes.
      if (newaph == 1) call aphread
      call crumpetread
c ... Set up tables for impurity atomic-physics processes.
      if (isimpon .eq. 1) then		# obsolete option
         call xerrab ('ueinit -- option isimpon=1 is obsolete; use 2')
      elseif (isimpon .eq. 2) then	# data supplied by D. Post 1993
         call readpost(coronalimpfname)
         call splinem
         call remark('*** For isimpon=2, set afracs, not afrac ***')
      elseif ((isimpon .eq. 3) .and. (nzspt .gt. 0)) then    # avg-ion
         impflag = 1
         crni = 1.
         ctemp = 1.
         call inelinput (impflag, crni, ctemp, zi(nhsp+1), nzspt)
      elseif ((isimpon .ge. 4) .and. (isimpon .le. 6)
     .                         .and. (nzspt .gt. 0)) then    # multi-charge
         if (ismctab .eq. 1) then        # use INEL multi-charge tables
            impflag = 2
            crni = 1.
            ctemp = 1.
            call inelinput (impflag, crni, ctemp, zi(nhsp+1), nzspt)
         elseif (ismctab .eq. 2) then    # use Braams multi-charge tables
            if (newapi == 1) call readmc(nzdf,mcfilename)
         endif
      elseif (isimpon .eq. 7) then      # read both Post and Braams tables
c.....First the Post table(s):
         call readpost(coronalimpfname)
         call splinem
         call remark('*** For isimpon=7, set afracs, not afrac ***')
c.....Then the Braams table(s):
         if (ismctab .eq. 1) then        # use INEL multi-charge tables
            impflag = 2
            crni = 1.
            ctemp = 1.
            call inelinput (impflag, crni, ctemp, zi(nhsp+1), nzspt)
         elseif (ismctab .eq. 2) then    # use Braams multi-charge tables
            if (newapi == 1) call readmc(nzdf,mcfilename)
         endif
      endif

c ... If isimpon > 0 and isph_sput = 1, init. DIVIMP data for physical sputt.
      do igsp = 1, ngsp    # note only one species should have isph_sput=1
         if (isimpon.gt.0 .and. isph_sput(igsp).ge.1) then
            call syld96(matt,matp,cion,cizb,crmb)
         endif
      enddo

c ... Set up new grid geometry, if desired.
      if(newgeo .eq. 1) then
         call nphygeo
         newgeo = 0 # Disable re-reading of grid unless explicitly requested

c...  "zero out" sx at ixpt2(1) if isfixlb=2 and ixpt1(1).le.0 to prevent
c...  flux thru cut
         if (isfixlb(1).eq.2 .and. ixpt1(1).le.0 .and. ixpt2(1).ge.0) then
            do iy = 0, iysptrx1(1)
               sx(ixpt2(1),iy) = 1.e-10*sx(ixpt2(1),iy)
            enddo
         endif
         if (isfixrb(1).eq.2 .and. ixpt1(1).gt.0 .and. ixpt2(1).ge.nx) then
            do iy = 0, iysptrx1(1)
               sx(ixpt1(1),iy) = 1.e-10*sx(ixpt1(1),iy)
            enddo
         endif
c...  Calculate parallel connection length along B
         do iy = 0, ny+1  # Initialize lconi,e to large number
           do ix = 0, nx+1
             lconi(ix,iy) = 1e50
             lcone(ix,iy) = 1e50
           enddo
         enddo
         if (iysptrx >= 1 .and. nxleg(1)+nxleg(2) > 0) then 
                                          # need sep for this calc of lconi,e
           call conlen   
         endif
      endif

c...  rescale magnetic field quantities with b0
      call s2scal (nx+2, ny+2, abs(b0/b0old), btot, 1, nx+2)
      call s2scal (nx+2, ny+2, sign(1.,b0/b0old), rbfbt, 1, nx+2)
      call s2scal (nx+2, ny+2, b0old/b0, rbfbt2, 1, nx+2)
      call s2scal (nx+2, ny+2, abs(b0/b0old), rbpol, 1, nx+2)
      call s2scal (nx+2, ny+2, b0old/b0, curvrby, 1, nx+2)
      call s2scal (nx+2, ny+2, abs(b0old/b0), curvrb2, 1, nx+2)
      call s2scal (nx+2, ny+2, b0old/b0, gradby, 1, nx+2)
      call s2scal (nx+2, ny+2, abs(b0old/b0), gradb2, 1, nx+2)
         # v2cb sign change comes from rbfbt in forming uu
      b0old = b0

*  -- test vol, gx, gy, sx, sy, rr --
      if (s2min(nx+2, ny+2, vol, 1, nx+2) .le. 0) then
         call xerrab ('ueinit -- error in sign vol')
      else if (s2min(nx+2, ny+2, gx, 1, nx+2) .le. 0 .or.
     .      s2min(nx+2, ny+2, gy, 1, nx+2) .le. 0) then
         call xerrab ('ueinit -- error in sign gx, gy')
      else if (s2min(nx+1, ny+2, sx, 1, nx+2) .le. 0 .or.
     .      s2min(nx+2, ny+1, sy, 1, nx+2) .le. 0) then
         call xerrab ('ueinit -- error in sign sx, sy')
      else if (s2min(nx+2, ny+2, rr, 1, nx+2) .le. 0) then
         call xerrab ('ueinit -- error in sign rr')
      end if

c ... Compute field-line length on the SOL flux surface that is
c     half way out (in grid space).
      linelen = 0.
      iy = (iysptrx1(1) + ny) / 2
c ... MER NOTE: For a full double-null configuration, iysptrx is the last
c               closed flux surface.  For an un-balanced double-null, iy
c               could be in the region between separatrices or it could be
c               beyond the outermost separatrix.
      do ix = 1, nx
         linelen = linelen + dx(ix,iy) / rr(ix,iy)
      enddo

c ... Compute half-range and weights for digital filter of turbulent
c     diffusivity, if needed.
      if (kyet .gt. 1.e-20 .and. isturbcons .eq. 2) then
         ixmp4 = ixpt1(1) + nxomit + 3*(ixpt2(1)-ixpt1(1))/4
c        This is the approximate poloidal location of the outboard midplane
c        for a single-null magnetic configuration.
         if (nxpt==2) ixmp4 = ixmdp(2) # outboard midplane of "dnull" config
         diffuslimit = min(9,
     .                 nint(0.5 * diffusrange / dy(ixmp4,iysptrx+1)))
         do jy = diffuslimit+1, 9
           diffuswgts(jy) = 0.
           diffuswgts(-jy) = 0.
         enddo
         diffuswgts(0) = 1.
         diffustotal = diffuswgts(0)
         do jy = 1, diffuslimit
           diffuswgts(jy) = float(diffuslimit + 1 - jy) /
     .                           (diffuslimit + 1)
           diffuswgts(-jy) = diffuswgts(jy)
           diffustotal = diffustotal + 2. * diffuswgts(jy)
         enddo
         do jy = -diffuslimit, diffuslimit
           diffuswgts(jy) = diffuswgts(jy) / diffustotal
         enddo
      endif

c...  set arrays to possibly zero out the neutral diffusive velocity
c     arising from grad Ti
      call sfill(nx+2, 1., fgtdx(0:), 1)
      call sfill(ny+2, 1., fgtdy(0:), 1)
      do jx = 1, nxpt # gradT can cause BC prob.;only flux matters
        fgtdx(ixlb(jx)) = gcfacgtx
        fgtdx(ixrb(jx)) = gcfacgtx
      enddo
        fgtdy(0)  = gcfacgty
        fgtdy(ny) = gcfacgty

c...  set flux-limit arrays and account for turning-off at boundaries
      call sfill(nx+2, flalfe, flalfea(0:), 1)
      call sfill(nx+2, flalfi, flalfia(0:), 1)
      call sfill(nx+2, flalfv, flalfva(0:), 1)
      do igsp = 1, 10
        call sfill(nx+2, flalfgx(igsp), flalfgxa(0:,igsp), 1)
        call sfill(nx+2, flalfgxy(igsp), flalfgxya(0:,igsp), 1)
        call sfill(ny+2, flalfgy(igsp), flalfgya(0:,igsp), 1)
      enddo
      call sfill(nx+2, flalfvgx, flalfvgxa(0:), 1)
      call sfill(nx+2, flalfvgxy, flalfvgxya(0:), 1)
      call sfill(ny+2, flalfvgy, flalfvgya(0:), 1)
      call sfill(nx+2, flalftgx, flalftgxa(0:), 1)
ccc      call sfill(nx+2, flalftgxy, flalftgxya(0:), 1)
      call sfill(ny+2, flalftgy, flalftgya(0:), 1)

      do jx = 1, nxpt  #loop over x-points
        if (isplflxl==0) then
          flalfea(ixlb(jx)) = 1.e20
          flalfia(ixlb(jx)) = 1.e20
          flalfea(ixrb(jx)) = 1.e20
          flalfia(ixrb(jx)) = 1.e20
        endif
        if (isplflxlv==0) then   # stagger mesh ==> ixlb+1 is bndry visc
          flalfva(ixlb(jx)+1) = 1.e20
          flalfva(ixrb(jx)+1) = 1.e20
        endif
        if (isplflxlgx==0) then
          do igsp = 1, 10
            flalfgxa(ixlb(jx),igsp) = 1.e20    
            flalfgxya(ixlb(jx),igsp) = 1.e20    
            flalfgxa(ixrb(jx),igsp) = 1.e20    
            flalfgxya(ixrb(jx),igsp) = 1.e20    
          enddo
        endif
        if (isplflxlvgx==0) then
          flalfvgxa(ixlb(jx)+1) = 1.e20    
          flalfvgxya(ixlb(jx)+1) = 1.e20    
          flalfvgxa(ixrb(jx)+1) = 1.e20    
          flalfvgxya(ixrb(jx)+1) = 1.e20    
        endif
        if (isplflxltgx==0) then
          flalftgxa(ixlb(jx)) = 1.e20    
          flalftgxya(ixlb(jx)) = 1.e20    
          flalftgxa(ixrb(jx)) = 1.e20    
          flalftgxya(ixrb(jx)) = 1.e20    
        endif
      enddo  # end of loop over x-point indices (ixpt)

c...  Now set sidewall flux limit factors
      if (iswflxlgy==0) then
        do igsp = 1, 10
          flalfgya(0,igsp) = 1.e20
          flalfgya(ny,igsp) = 1.e20
        enddo
      endif
      if (iswflxlvgy==0) then
        flalfvgya(0) = 1.e20
        flalfvgya(ny) = 1.e20
      endif
      if (iswflxltgy==0) then
        flalftgya(0) = 1.e20
        flalftgya(ny) = 1.e20
      endif
c...  set wall sources
      call walsor

c...  set plate sources
      call pltsor

c...  set plate recycling coefficient profiles
      call recyprof

c...  set volume power sources if the internal Gaussian sources desired
      if (isvolsorext == 0) call volsor

c ... Set impurity sources on inner and outer walls; poss prob if nyomitmx>0
      if (isimpwallsor == 1) then  #impurity wall-flux sources
        call imp_sorc_walls (nx, nzspt, xcpf, xcwo, sy(0,0), sy(0,ny),
     .                        ixp1(0,0), ixp1(0,ny), fnzysi, fnzyso)
      endif

c ... Initialize molecular thermal equilibration array in case not computed
      do igsp = 1,ngsp
        call s2fill (nx+2, ny+2, 0.0e0, eqpg(0:,0:,igsp), 1, nx+2)
      enddo

*---  bbbbbb begin ifloop b  bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
      if (restart .eq. 0) then
*---  bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
*     ------------------------------------------------------------------
*     initialize plasma state.
*     ------------------------------------------------------------------

*  -- initialize the density and velocity
      do ifld = 1, nisp
         call s2fill (nx+2, ny+2, nibeg(ifld), ni(0:,0:,ifld), 1, nx+2)
      enddo
      do ifld = 1, nisp
         call s2fill (nx+2, ny+2, 0.0e0, uu(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0.0e0, vy(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0.0e0, up(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0.0e0, frici(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0.0e0, nuvl(0:,0:,ifld), 1, nx+2)
      enddo

*  -- initialize the gas and electron density

      do 142 iy = 0, ny+1
         do 141 ix = 0, nx+1
            ne(ix,iy) = 0.0
            nit(ix,iy) = 0.0
            nm(ix,iy,1) = 0.0
            do 140 ifld = 1, nisp  # here init. gas only for diff. model
               ng(ix,iy,1) = ngscal(igsp)*nibeg(1)*( exp(-xcs(ix)/xgscal)
     .                            + exp(-(xcs(nx+1)-xcs(ix))/xgscal) )
     .                      + ngbackg(1)
               nginit(ix,iy) = ng(ix,iy,1)
               ne(ix,iy) = ne(ix,iy) + zi(ifld)*ni(ix,iy,ifld)
               if (zi(ifld).ne.0) then
c Dont do if this is the neutral momentum equation
                  nit(ix,iy) = nit(ix,iy) + ni(ix,iy,ifld)
                  if (isimpon >= 5 .and. nusp_imp == 0)
     .                  nm(ix,iy,1)=nm(ix,iy,1)+ni(ix,iy,ifld)*mi(ifld)
               endif
               nm(ix,iy,ifld) = ni(ix,iy,ifld)*mi(ifld)
 140        continue
 141     continue
 142  continue

c...  This is redundant with above if ngsp=1, but done to set nginit
      do 145 iy = 0, ny+1
         do 144 ix = 0, nx+1
            do 143 igsp = 1, ngsp
               ng(ix,iy,igsp) = ngscal(igsp)*nibeg(1)*(
     .                                          exp(-xcs(ix)/xgscal)
     .                            + exp(-(xcs(nx+1)-xcs(ix))/xgscal) )
     .                      + ngbackg(igsp)
	       tg(ix,iy,igsp) = tscal*ttbeg
               nginit(ix,iy) = ng(ix,iy,1)
 143           continue
 144        continue
 145     continue

c...  Initialize 4th order fluxes
      do ifld = 1, nisp
        call s2fill (nx+2, ny+2, 0., fniy4ord(0:,0:,ifld), 1, nx+2)
      enddo

*  -- Initialize temperatures, potential, currents, and some nonog-fluxes.
      call s2fill (nx+2, ny+2, ttbeg, te, 1, nx+2)
      call s2fill (nx+2, ny+2, ttbeg/ev, phi, 1, nx+2)
      call s2fill (nx+2, ny+2, tscal*ttbeg, ti, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fqx, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fqy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fq2, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fqp, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., vytan, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fngxy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., feexy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., feixy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fmixy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., frice, 1, nx+2)

c ... Give some shape to these initial values
         do iy = 0, ny+1
            do ix = 0, nx+1
               if (isfixlb(1) .gt. 0) then
                 proffacx = float(nx+3-ix)/float(nx+3)
                 proffacy = float(ny+3-iy)/float(ny+3)
                 proffacv = proffacx
               else
                 proffacx = cos(float(ix-nx/2)*(pi-1.)/float(nx))
                 proffacy = float(ny+3-iy)/float(ny+3)
                 proffacv = sin(float(ix-nx/2)*(pi-1.)/float(nx))
               endif
               te(ix,iy) = ttbeg*proffacx*proffacy
               ti(ix,iy) = tscal*ttbeg*proffacx*proffacy
               do ifld = 1, nisp
                  ni(ix,iy,ifld) = nibeg(ifld)*proffacy
cc                  ni(ix,iy,ifld) = nibeg(ifld)*ttbeg/te(ix,iy)
                  up(ix,iy,ifld) = sqrt(te(0,0)/mi(1))*proffacv*
     .                             proffacy
                  up(nx+1,iy,ifld) = up(nx,iy,ifld)
               enddo
            enddo
         enddo

c ... Symmetrize the profiles in x-direction if isgrdsym=1
      if (isgrdsym.eq.1) then
         do iy = 0, ny+1
            do ix = 0, nx/2
               te(ix,iy) = te(nx+1-ix,iy)
               ti(ix,iy) = ti(nx+1-ix,iy)
               do ifld = 1, nisp
                  ni(ix,iy,ifld) = ni(nx+1-ix,iy,ifld)
                  up(ix,iy,ifld) = -up(nx-ix,iy,ifld)
                  up(nx/2,iy,ifld) = 0.
                  if (iy .le. iysptrx) then
                     up(ixpt1(1),iy,ifld) = 0.
                     up(ixpt2(1),iy,ifld) = 0.
                  endif
               enddo
            enddo
         enddo
      endif

      call convert
c...  Initializes the variables for the daspk package if this is the
c...  method chosen.
*---  bbbbbbbb ifloop b else  bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
      else
*---  bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
*---  Here the new plasma state is obtained from interpolation

*  -- test xcs, xfs, yyc and yyf --
         do 111 ix = 0, nx
            if (xcs(ix) .gt. xfs(ix) .or.
     .          xfs(ix) .gt. xcs(ix+1)) then
               call xerrab ('ueinit -- error involving xcs, xfs')
            endif
  111    continue
         do 121 iy = 0, ny
            if (yyc(iy) .gt. yyf(iy) .or.
     .          yyf(iy) .gt. yyc(iy+1)) then
               call xerrab ('ueinit -- error involving yyc, yyf')
            endif
  121    continue

*  -- Check if mesh size has changed; then must use icntnunk=0
      if (nx.ne.nxold .or. ny.ne.nyold) then
         if (icntnunk .eq. 1) then
            call xerrab ('** nx or ny changed; must set icntnunk=0 **')
         endif
      endif

*     ------------------------------------------------------------------
*     interpolate plasma state to the larger grid; either same or double
*     ------------------------------------------------------------------

      if(isnintp.eq.0 .or. isimesh.eq.1) then  # Use old interp or copy
         if (nis(nxold+1,nyold+1,1) .eq. 0.) then
            call xerrab ('variable nis=0, no saved solution for interp')
         endif
       if(nx.ne.nxold .and. ny.ne.nyold) then
         if(nx.ne.2*nxold .and. ny.ne.2*nyold) then
            call xerrab ('must double both nx and ny to interpolate')
         endif

c...  If the grid is doubled in each direction, interpolate the solution
c...  and save it for possible subsequent restarts
         call refpla

         nxold = nx
         nyold = ny
         call gchange("Interp",0)
         call gridseq

       else

c.... If the grid does not change, but restart from saved variables
          do ifld = 1, nisp
             do iy = 0, ny+1
                do ix = 0, nx+1
                   ni(ix,iy,ifld) = nis(ix,iy,ifld)
                   up(ix,iy,ifld) = ups(ix,iy,ifld)
                   if (nis(ix,iy,ifld) <= 0.) then
                      call remark('*** Error: nis <= 0 ***')
                      write (*,*) 'Error at ix=', ix,'  iy=',iy,' ifld=',ifld
                      call xerrab("")
                   endif
                enddo
             enddo
          enddo
          do iy = 0, ny+1
             do ix = 0, nx+1
                do igsp = 1, ngsp
                   ng(ix,iy,igsp) = ngs(ix,iy,igsp)
                   if (ngs(ix,iy,igsp) <= 0.) then
                      call remark('*** Error: ngs <= 0 ***')
                      write (*,*) 'Error at ix=', ix,'  iy=',iy
                      call xerrab("")
                   endif
                   tg(ix,iy,igsp) = tgs(ix,iy,igsp)
                enddo
                te(ix,iy)      = tes(ix,iy)
                ti(ix,iy)      = tis(ix,iy)
                phi(ix,iy) = phis(ix,iy)
                if (isimpon .eq. 2 .or. isimpon .eq. 7) then
                  if (afracs(1,1)+afracs(nx,ny).gt.1.e-20) then
                    afrac(ix,iy) = afracs(ix,iy)
                  else
                    afracs(ix,iy) = 1.e-20
                    afrac(ix,iy) = afracs(ix,iy)
                    call remark('***WARNING: 
     .                          setting afracs = 1.e-20; 0 is illegal')
                  endif
                endif
                if (phis(nxold-1,nyold-1) .eq. 0) phi(ix,iy) = 40.  #avoid phi=0.
             enddo
          enddo

       endif

      else          # New interpolator section for isnintp=1

c...  Calculate y-interpolated norm. poloidal grid pts, and put in xnrmox
c...  for density, etc, and in xvnrmox for poloidal velocity
c...  Likewise for y-interpolated yn values put in ynrmox and yvnrmox
c...  Also, must allocate Interp arrays with current nx,ny (diff. from nxold)

c...  Order ixstart/end of poloidal regions from new ixcut & ixlb,rb
c...  Previous values of ixsto and ixendo computed in subr gridseq
      ixst(1) = ixlb(1)
      ixend(1) = ixcut1
      if (ixlb(1) == 0 .and. ixcut1 == 0) then  # no inner leg
        ixst(2) = 0
      else
        ixst(2) = max(ixlb(1), ixcut1+1)
      endif
      ixend(2) = ixcut2
      if (nyomitmx >= nysol) then   # no inner/outer leg region
         ixst(2) = 0
         ixend(2) = nx+1
      endif
      if (nx == 1) then  #special case: 1D in radial direction
         ixend(2) = 2
      endif
c..   Now need to check if ixrb is > or < ixcut3
      ixst(3) = ixcut2+1
      if (ixcut3 > ixrb(1) .or. nxpt==1) then  #3 regions in first domain
        ixend(3) = ixrb(1)+1
      else    # 4 regions in 1st domain, end on ixcut3
        ixend(3) = ixcut3
      endif

c..   Continue ordering if double null or snowflake
      if (nxpt == 2) then
        if (ixcut3 > ixrb(1)) then  # do 3-region 2nd domain
          ixst(4) = ixlb(2)
          ixend(4) = ixcut3
          ixst(5) = ixcut3+1
        else # 4 regions in 1st domain, compl & do 2-region 2nd domain
          ixst(4) = ixcut3+1
          ixend(4) = ixrb(1)+1
          ixst(5) = ixlb(2) 
        endif  # remain indices are the same
        ixend(5) = ixcut4
        ixst(6) = ixcut4+1
        ixend(6) = ixrb(2)+1
      endif  # if-test on nxpt
              
c...  Set number of regions for interpolation; 3 for single null, 6
c...  for double null, and 1 for core-only simulations
      if (nysol <= nyomitmx) then  #core only, no divertor legs
         irstart = 2
         irend = 2
         ixst(2) = 1
         ixsto(2) = 1
         ixend(2) = ixend(2) - 1
         ixendo(2) = ixendo(2) - 1
      elseif (nxpt == 1) then  #single-null with SOL, divertor legs
         irstart = 1
         irend = 3
      else   #must be double-null with nxpt=2
          irstart = 1
          irend = 6
      endif

c...  Construct first intermediate density grid, (xnrmox,ynrmox)
      call gchange("Interp",0)
      do ir = irstart, irend
        call grdintpy(ixst(ir),ixend(ir),ixsto(ir),ixendo(ir),
     .                0,ny+1,0,nyold+1,nx,ny,nxold,nyold,
     .                xnrm,ynrm,xnrmo,ynrmo,xnrmox,ynrmox,ixmg,iyomg)
      enddo
         
c...  Construct second intermediate density grid, (xnrmnx,ynrmnx)
      do ir = irstart, irend
        call grdintpy(ixsto(ir),ixendo(ir),ixst(ir),ixend(ir),
     .                0,ny+1,0,ny+1,nxold,ny,nx,ny,
     .                xnrmox,ynrmox,xnrm,ynrm,xnrmnx,ynrmnx,ix2g,iy2g)
      enddo

c...  Construct first intermediate velocity grid (xvnrmox,yvnrmnox)
      do ir = irstart, irend
        call grdintpy(ixst(ir),ixend(ir),ixsto(ir),ixendo(ir),
     .                0,ny+1,0,nyold+1,nx,ny,nxold,nyold,
     .                xvnrm,yvnrm,xvnrmo,yvnrmo,xvnrmox,yvnrmox,
     .                ixvmg,iyvomg)
      enddo

c...  Construct second intermediate velocity grid (xvnrmnx,yvnrmnx)
      do ir = irstart, irend
        call grdintpy(ixsto(ir),ixendo(ir),ixst(ir),ixend(ir),
     .                0,ny+1,0,ny+1,nxold,ny,nx,ny,
     .                xvnrmox,yvnrmox,xvnrm,yvnrm,xvnrmnx,yvnrmnx,
     .                ixv2g,iyv2g)
      enddo

c...  Fix the special cell ixpt2(1)=ixrb(1) for geometry="isoleg"
      if (geometry == "isoleg") then
        do iy = 0, ny+1
          xnrmox(ixrbo(1)+1,iy) = 1.
          xnrmnx(ixrb(1)+1,iy) = 1.
          ynrmox(ixrbo(1)+1,iy) = ynrmox(ixrbo(1),iy)
          ynrmox(ixlb(2),iy) = ynrmox(ixlb(2)+1,iy)
          ynrmnx(ixrb(1)+1,iy) = ynrmnx(ixrb(1),iy)
          ynrmnx(ixlb(2),iy) = ynrmnx(ixlb(2)+1,iy)
          xvnrmox(ixrbo(1)+1,iy) = 1.
          xvnrmnx(ixrb(1)+1,iy) = 1.
          yvnrmox(ixrbo(1)+1,iy) = ynrmox(ixrbo(1),iy)
          yvnrmox(ixlb(2),iy) = ynrmox(ixlb(2),iy)
          yvnrmnx(ixrb(1)+1,iy) = ynrm(ixrb(1),iy)
          yvnrmnx(ixlb(2),iy) = ynrm(ixlb(2)+1,iy)
        enddo
      endif

c...  Now interpolate the plasma variables

         call intpvar (tes, te, 0, nxold, nyold)
         call intpvar (tis, ti, 0, nxold, nyold)
         call intpvar (phis, phi, 0, nxold, nyold)
c...  Interpolate the relative fraction of impurities
         if (isimpon>0 .and. afracs(1,1)+afracs(nxold,nyold)>1.e-20)
     .                        call intpvar (afracs,afrac,0,nxold,nyold)

c...  If phis(nx-1,ny-1)=0., reset to constant 40 volts
         if (phis(nxold-1,nyold-1).eq.0.)
     .                      call s2fill (nx+2, ny+2, 40., phi, 1, nx+2)

         do 610 ifld = 1, nisp
            call intpvar (nis(0:,0:,ifld), ni(0:,0:,ifld), 0, nxold, nyold)
            call intpvar (ups(0:,0:,ifld), up(0:,0:,ifld), 1, nxold, nyold)
 610     continue
         do 620 igsp = 1, ngsp
            call intpvar (ngs(0:,0:,igsp), ng(0:,0:,igsp), 0, nxold, nyold)
            call intpvar (tgs(0:,0:,igsp), tg(0:,0:,igsp), 0, nxold, nyold)
 620     continue

c...  Reset gas density to minimum if too small or negative
         do igsp = 1, ngsp
           do  iy = 0, ny+1
             do ix = 0, nx+1
               if(isngonxy(ix,iy,igsp)==1) then
                  ng(ix,iy,igsp) = max(ng(ix,iy,igsp),
     .                                 1.0e-01*ngbackg(igsp))
               endif
             enddo
           enddo
         enddo

      endif          # end of very-large 2-branch-if: (1), same mesh size
                     # or (2), index-based interp with isnintp=1 

      if (nyomitmx >= nysol+nyout) then
        call filldead_guardcells
      endif
c...  Check if any ion density is zero
      do ifld = 1, nisp
        do iy = 0, ny+1
          do ix = 0, nx+1
            if (ni(ix,iy,ifld) <= 0) then
              call remark('****** ERROR: ni <= 0 ******')
              write(*,*) 'begins at ix,iy,ifld = ',ix,iy,ifld
              call xerrab("")
            endif
          enddo
        enddo
      enddo

*  -- initialize nginit and ne to interpolated value, zero nonog-fluxes
      call s2copy (nx+2, ny+2, ng, 1, nx+2, nginit, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., ne, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., nit, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fqx, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fqy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fq2, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., fqp, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., feexy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., feixy, 1, nx+2)
      call s2fill (nx+2, ny+2, 0., frice, 1, nx+2)

      do ifld = 1, nisp
         call s2fill (nx+2, ny+2, 0., vytan(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., nm(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., psorc(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., psorxr(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., msor(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., msorxr(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., nucxi(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., nueli(0:,0:,ifld), 1, nx+2)
      enddo

      do ifld = 1, nisp   # test that s2fill does the right thing
        psorold(ifld) = 0.
        psorxrold(ifld) = 0.
        msorold(ifld) = 0.
        msorxrold(ifld) = 0.
        do iy = 0, ny+1
          do ix = 0, nx+1
            psorc(ix,iy,ifld) = 0.
            psorxr(ix,iy,ifld) = 0.
            msor(ix,iy,ifld) = 0.
            msorxr(ix,iy,ifld) = 0.
            nucxi(ix,iy,ifld) = 0.
            nueli(ix,iy,ifld) = 0.
          enddo
        enddo
      enddo

      do ifld = 1, nusp
         call s2fill (nx+2, ny+2, 0., fmixy(0:,0:,ifld), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., frici(0:,0:,ifld), 1, nx+2)
      enddo

      do igsp = 1, ngsp
         call s2fill (nx+2, ny+2, 0., fngxy(0:,0:,igsp), 1, nx+2)
         call s2fill (nx+2, ny+2, 0., fegxy(0:,0:,igsp), 1, nx+2)
      enddo

      do 713 ifld = 1, nisp
         do 712 iy = 0, ny+1
            do 711 ix = 0, nx+1
               if (zi(ifld).ne.0.) then
                  ne(ix,iy) = ne(ix,iy) + zi(ifld)*ni(ix,iy,ifld)
                  nit(ix,iy) = nit(ix,iy) + ni(ix,iy,ifld)
		  if (isimpon>=5 .and. nusp_imp==0) #note nm(ix,iy,1) initlly=0
     .                  nm(ix,iy,1)=nm(ix,iy,1)+ni(ix,iy,ifld)*mi(ifld)
               endif
               nm(ix,iy,ifld) = ni(ix,iy,ifld)*mi(ifld)
 711        continue
 712     continue
 713  continue

c...  Set boundary conditions for ni and Te,i on walls if end-element zero
      if (nwalli(nx+1).lt.1.e-10) call sfill (nx+2,nwalli(0),nwalli(0:),1)
      if (nwallo(nx+1).lt.1.e-10) call sfill (nx+2,nwallo(0),nwallo(0:),1)
      if (tewalli(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tewalli(0:),1)
      if (tiwalli(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tiwalli(0:),1)
      if (tewallo(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tewallo(0:),1)
      if (tiwallo(nx+1).lt.1.e-10) call sfill (nx+2,tedge,tiwallo(0:),1)

c...  Initialize dead pol guard cells if core-only simulation
      if (nyomitmx >= nysol+nyout) then
        call filldead_guardcells
      endif
         
      call convert
*     ------------------------------------------------------------------


*---  bbbbbb end ifloop b bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
      endif
*---  bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb

c ... Now that new indexing (e.g., ixendi, ixendo etc.) is done, set
c ... impurity wall sources that depend on these indices.
ccc     call imp_sorc_walls (nx, nzspt, xcpf, xcwo, sy(0,0), sy(0,ny),
ccc     .                        ixp1(0,0), ixp1(0,ny), fnzysi, fnzyso)

c ... Set variable-normalization array.
      call set_var_norm (iscolnorm, neq, numvar, yl, norm_cons,
     .                   floor_cons, suscal)
c...  Need sfscal initialized in jac_calc if not continuing
      if (icntnunk .eq. 0) call sfill (neq, 1., sfscal(1:), 1)

c...  set the stretching array for the poloidal coordinate used in the
c...  poloidal diffusion for the neutral gas
      do jx = 1, nxpt
      do 223 iy = 0, ny+1
         do 222 ix = ixlb(jx), ixrb(jx)+1
            factor   = 0.5*( exp(-(xcs(ix)  - xcs(ixlb(jx)))/xstscal) +
     .                       exp(-(xcs(ixrb(jx)+1)- xcs(ix))/xstscal) )
            stretcx(ix,iy) = 1 + (sxgpr**2 - 1) * factor
            if (iy .gt. max(iysptrx1(jx),iysptrx2(jx))) then
cccMER NOTE: only use sxgsol beyond outermost separatrix
               stretcx(ix,iy) = 1 + (sxgsol**2 - 1) * factor
            endif
 222     continue
 223  continue
      enddo

c ... Enable Jac stencil comp if not parallel
      if (isjacstnlon == 1) then
        call domain_dc   # comp Jacobian stencil ivl2gstnl
      endif

c...  Set eymask1d to give ey=0 in core+sep for 1d SOL pot (isphicore0=1)
      eymask1d = 1.  #2D array initialization
      if (isphicore0 == 1) then  #only solve pot eqn in SOL; phi_core const
        do jx = 1, nxpt
          do iy = 0, iysptrx
            do ix = ixpt1(jx)+1, ixpt2(jx)
              eymask1d(ix,iy) = 0.
            enddo
          enddo
        enddo
      endif

      return
      end
c***** end of subroutine ueinit ****************************************

c----------------------------------------------------------------------c

      subroutine set_indirect_address(isglobal)
c     Set indirect addressing arrays for x-direction
      implicit none

c..   Input variables

      integer isglobal  #=1, global mesh; =0, serial case or par domains

Use(Dim)                # nx,ny
Use(Share)              # geometry,nyomitmx
Use(Xpoint_indices)     # ixpt1,ixpt2,iysptrx1,iysptrx2
Use(Cut_indices)	# ixcut1,iycut1,ixcut2,iycut2,ixcut3,iycut3
                        # ixcut4,iycut4
Use(Selec)              # ixm1,ixp1
Use(Bcond)              # isfixlb,isfixrb
c     local variables --
      integer ix,iy,jx

c...  Set cut indices (duplicative for now, but used for snowflake)
      ixcut1 = ixpt1(1)
      iycut1 = iysptrx1(1)
      ixcut2 = ixpt2(1)
      iycut2 = iysptrx2(1)
      if (nxpt == 2) then
        ixcut3 = ixpt1(2)
        iycut3 = iysptrx1(2)
        ixcut4 = ixpt2(2)
        iycut4 = iysptrx2(2)
      endif

      do iy = 0, ny+1
            iym1a(0,iy) = max(0,iy-1)
            iyp1a(0,iy) = min(ny+1,iy+1)
            iym1a(nx+1,iy) = max(0,iy-1)
            iyp1a(nx+1,iy) = min(ny+1,iy+1)

c ...  First case is the default geometry=snull
            do ix = 1, nx
               iym1a(ix,iy) = max(0,iy-1)
               iyp1a(ix,iy) = min(ny+1,iy+1)
               if (iy.le.iysptrx1(1) .and. ix.eq.ixpt1(1)) then
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ixpt2(1)+1
               elseif (iy.le.iysptrx1(1) .and. ix.eq.(ixpt1(1)+1)) then
                  ixm1(ix,iy) = ixpt2(1)
                  ixp1(ix,iy) = ix+1
               elseif (iy.le.iysptrx2(1) .and. ix.eq.ixpt2(1)) then
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ixpt1(1)+1
               elseif (iy.le.iysptrx2(1) .and. ix.eq.(ixpt2(1)+1)) then
                  ixm1(ix,iy) = ixpt1(1)
                  ixp1(ix,iy) = ix+1
               else
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ix+1
               endif

               if (geometry=="dnull") then  #3 conditions for 1-cell cases
                 if (ixpt2(1)==ixpt1(1)+1.or.ixpt2(2)==ixpt1(2)+1) then
                  call xerrab("***Error: Single pol cell not supported")
                 endif
                  if (iy.le.iysptrx1(1) .and. ix.eq.ixpt1(1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixpt2(2)+1
                  elseif (iy.le.iysptrx1(1) .and. ix.eq.(ixpt1(1)+1)) then
                     ixm1(ix,iy) = ixpt2(2)
                     ixp1(ix,iy) = ix+1
                  elseif (iy.le.iysptrx2(1) .and. ix.eq.ixpt2(1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixpt1(2)+1
                  elseif (iy.le.iysptrx2(1) .and. ix.eq.(ixpt2(1)+1)) then
                     ixm1(ix,iy) = ixpt1(2)
                     ixp1(ix,iy) = ix+1
                  elseif (iy.le.iysptrx1(2) .and. ix.eq.ixpt1(2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixpt2(1)+1
                  elseif (iy.le.iysptrx1(2) .and. ix.eq.(ixpt1(2)+1)) then
                     ixm1(ix,iy) = ixpt2(1)
                     ixp1(ix,iy) = ix+1
                  elseif (iy.le.iysptrx2(2) .and. ix.eq.ixpt2(2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixpt1(1)+1
                  elseif (iy.le.iysptrx2(2) .and. ix.eq.(ixpt2(2)+1)) then
                     ixm1(ix,iy) = ixpt1(1)
                     ixp1(ix,iy) = ix+1
                  else
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif

               if (geometry=="snowflake15") then	# MER 20 Apr 2014
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut3+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ixcut3
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut2+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ix+1
                     if ((ixcut4 .eq. ixcut3+1) .and. (iy .le. iycut4)) then
                        ixp1(ix,iy) = ixcut1+1
                     endif
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     if (ixcut4 .eq. ixcut3+1) ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="snowflake15"

               if (geometry=="snowflake45" .or. geometry=="dnXtarget") then
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut2+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                     if (ixcut3 .eq. ixcut2+1) ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1
                     if ((ixcut3 .eq. ixcut2+1) .and. (iy.le.iycut2)) then
                        ixm1(ix,iy) = ixcut1
                     endif
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut3+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut3
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="snowflake45"

               if (geometry=="snowflake75") then	# MER 16 Sep 2014
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut2+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                     if ((ixcut3 .eq. ixcut2+1) .and. (iy .le. iycut3)) then
                        ixp1(ix,iy) = ixcut4+1
                     endif
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1
                     if (ixcut3 .eq. ixcut2+1) ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut3+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut3
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="snowflake75"
			   
               if (geometry=="snowflake105") then	# AK 10 DEC 2018
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                     if ((ixcut2 .eq. ixcut1+1) .and. (iy .le. iycut1)) then
                        ixp1(ix,iy) = ixcut3+1
                     endif
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     if ((ixcut2 .eq. ixcut1+1) .and. (iy .le. iycut1)) ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ixcut3+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ixcut3
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut2+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     if (ixcut4.eq.ixcut3+1) ixm1(ix,iy) = ixcut2	# never satisfied???
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="snowflake105"

               if (geometry=="snowflake135") then	# MER 24 JUL 2020
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                     if ((ixcut2 .eq. ixcut1+1) .and. (iy .le. iycut1)) then
                        ixp1(ix,iy) = ixcut3+1
                     endif
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     if ((ixcut2 .eq. ixcut1+1) .and. (iy .le. iycut1)) ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ixcut3+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ixcut3
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut2+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="snowflake135"

               if (geometry=="snowflake165") then	# MER 24 JUL 2020
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut3+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ixcut3
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut2+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ix+1
                     if ((ixcut4 .eq. ixcut3+1) .and. (iy .le. iycut4)) then
                        ixp1(ix,iy) = ixcut1+1
                     endif
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     if (ixcut4 .eq. ixcut3+1) ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="snowflake165"

               if (geometry=="isoleg") then	# TDR 03 Dec 2014
                  if ((iy .le. iycut1) .and. (ix .eq. ixcut1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ixcut4+1
                  elseif ((iy .le. iycut1) .and. (ix .eq. ixcut1+1)) then
                     ixm1(ix,iy) = ixcut4
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut2) .and. (ix .eq. ixcut2+1)) then
                     ixm1(ix,iy) = ix-1
		     ixp1(ix,iy) = ix+1  #mod
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3)) then
                     ixm1(ix,iy) = ix-1  #mod
                     ixp1(ix,iy) = ix+1
                  elseif ((iy .le. iycut3) .and. (ix .eq. ixcut3+1)) then
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                     if ((ixcut4 .eq. ixcut3+1) .and. (iy .le. iycut4)) then
                        ixp1(ix,iy) = ixcut1+1
                     endif
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4)) then
                     ixm1(ix,iy) = ix-1
                     if (ixcut4 .eq. ixcut3+1) ixm1(ix,iy) = ixcut2
                     ixp1(ix,iy) = ixcut1+1
                  elseif ((iy .le. iycut4) .and. (ix .eq. ixcut4+1)) then
                     ixm1(ix,iy) = ixcut1
                     ixp1(ix,iy) = ix+1
                  else	# when cuts do not interfere
                     ixm1(ix,iy) = ix-1
                     ixp1(ix,iy) = ix+1
                  endif
               endif	# end of geometry=="isoleg"

               if (isfixlb(1)==2.or.isfixrb(1)==2.or.iysptrx1(1)==0) then
                                         # no multiple-connections here
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ix+1
               elseif (iysptrx1(1)==0) then #phys bndry at sep
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ix+1
               endif
            enddo  # first end do-loop over ix


         do jx = 1, nxpt  #fix poloidal bdrys for multiple nulls
           ixm1(ixlb(jx),iy) = ixlb(jx)
           ixp1(ixlb(jx),iy) = ixlb(jx)+1
           ixm1(ixrb(jx)+1,iy) = ixrb(jx)
           ixp1(ixrb(jx)+1,iy) = ixrb(jx)+1
         enddo
      enddo  # end do-loop over iy

c...  Special fix for core-only cases
      if (nyomitmx >= nysol) then
         ixm1(1,ny+1) = ixpt2(1)
         ixp1(nx,ny+1) = 1
      endif

      return
      end  # end of subroutine set_indirect_address

c----------------------------------------------------------------------c
c-------------------------------------------------------------------------
      subroutine set_var_norm (job, neq, nvars, yl, norm_cons,
     .                         floor_cons, su)

c ... Set variable-normalization array su, depending on job:
c     job = 0 -- all elements of su = 1
c     job = 1 -- each element of su = inverse of the
c                normalization parameter in norm_cons
c     job = 2 -- each element of su = inverse of the
c                max of the abs of the variable value and the
c                floor parameter in floor_cons
c     job = 3 -- uses combination of global scaling with nnorm, etc
c                followed by local scaling by each yl
      implicit none

c ... Input arguments:
      integer job
      integer neq     # total number of equations (and variables)
      integer nvars   # number of variables in each grid cell
      real yl(neq)    # variables used with solver packages
      real norm_cons(nvars)    # normalization constants
      real floor_cons(nvars)   # floor constants

c ... Output argument:
      real su(neq)    # array of variable normalizations

c ... Local variables:
      integer ncells     # number of grid cells
      integer iv, i, j

      ncells = neq / nvars
      iv = 0
      if (job .eq. 0) then
         call sfill (neq, 1., su, 1)
      elseif (job .eq. 1) then
         do j = 1, ncells
            do i = 1, nvars
               iv = iv + 1
               su(iv) = 1. / norm_cons(i)
            enddo
         enddo
      elseif (job .eq. 2) then
         do j = 1, ncells
            do i = 1, nvars
               iv = iv + 1
               su(iv) = 1. / max(abs(yl(iv)), floor_cons(i))
            enddo
         enddo
      else
         do j = 1, ncells
            do i = 1, nvars
               iv = iv + 1
               su(iv) = 1. / max(abs(yl(iv)), floor_cons(i)/norm_cons(i))
            enddo
         enddo
      endif

      return
      end
c ****** End of subroutine set_var_norm ******************************
c-----------------------------------------------------------------------

      subroutine exmain_f

*
*     12/1/2019 - meyer8@llnl.gov: Changed the name of the entry point.
*     Exmain is now provided by a C source file. This allows us to 
*     trap SIGINT and provide a Basis-like debug mode in the Python
*     version of the code. There is no physics in the C source, only
*     system calls to handle the Control-C. When built with Basis the
*     C source file simply drops through to this entry point.
*
*
*
*     EXMAIN is the main subroutine for the two dimensional edge plasma
*     code UEDGE. The code solves a system of fluid equations that
*     models the edge plasma in an axisymmetric configuration.
*     The numerical procedure used is the method of lines that consist
*     of the solution of a set of coupled ODEs for the fluid variables
*     for each grid point.

******************************************************************************
*   Copyright 1994.  The Regents of the University of California.  All       *
*   rights reserved.                                                         *
*                                                                            *
*   This work was produced at the University of California, Lawrence
*   Livermore National Laboratory (UC LLNL) under contract no. W-7405-ENG-48
*   (Contract 48) between the U.S. Department of Energy (DOE) and The
*   Regents of the University of California (University) for the operation
*   of UC LLNL.  Copyright is reserved to the University for purposes of
*   controlled dissemination, commercialization through formal licensing, or
*   other disposition under terms of Contract 48; DOE policies, regulations
*   and orders; and U.S. statutes.  The rights of the Federal Government are
*   reserved under Contract 48 subject to the restrictions agreed upon by
*   the DOE and University as allowed under DOE Acquisition Letter 97-1.
*
*   			       DISCLAIMER
*
*   This software was prepared as an account of work sponsored by an agency
*   of the United States Government.  Neither the United States Government
*   nor the University of California nor any of their employees, makes any
*   warranty, express or implied, or assumes any liability or responsibility
*   for the accuracy, completeness, or usefulness of any information,
*   apparatus, product, or process disclosed, or represents that its
*   specific commercial products, process, or service by trade name,
*   trademark, manufacturer, or otherwise, does not necessarily constitute
*   or imply its endorsement, recommendation, or favoring by the United
*   States Government or the University of California.  The views and
*   opinions of authors expressed herein do not necessarily state or reflect
*   those of the United States Government or the University of California,
*   and shall not be used for advertising or product endorsement purposes.
*
******************************************************************************

      implicit none

      Use(Err_msg_out)   # errmsgflag,errunit
      Use(Cdv)      # ifexmain
      Use(Dim)      # nx,ny
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(UEpar)    # cnurn,cnuru,cnure,cnuri,cnurg,cnurp,
                    # nurlxn,nurlxu,nurlxe,nurlxi,nurlxg,nurlxp,
                    # label, istgon
      Use(Ident_vars)          # uedge_ver
      Use(Lsode)    # iterm,icntnunk
      Use(Grid)     # ngrid,nurlx,ijac,ijactot
      Use(Decomp)   # ubw,lbw
      Use(Share)    # 
      Use(Interp)   # isnintp,nxold,nyold
      Use(RZ_grid_info)  # rm,zm
      Use(UEint)               # isallloc
      Use(Rccoef)              # isoutwall
      Use(Coefeq)              # oldseec, override, cftiexclg
      Use(Flags)               # iprint
      Use(ParallelEval)
      Use(UEpar)
      Use(MCN_sources)
      Use(Imprad)

      integer ifake  #forces Forthon scripts to put implicit none above here
      integer my_pe, n_pes, ixx, icall
      data icall/0/

      call exmain_prelims

c=======================================================================
c//computation//


*     -- allocate memory for arrays --
*     -- set ifexmain=1 so that allocate knows exmain is the calling
*     -- routine.  This is necessary for grid sequencing coding.
           ifexmain = 1
           call allocate
           ifexmain = 0
	   if ((icall == 0) .and. (iprint .ne. 0)) write(*,*) 'UEDGE ',uedge_ver
           icall = 1
	   write(*,*) 'UEDGE version ',uedge_ver
           icall = 1

c TODO: Add check for inertial neutral model when Tg for atoms is on


c   Check model switches for UEDGE updates/bugs
      if (isoldalbarea .ne. 0) then
            write(*,*) "           **** WARNING ****"
            write(*,*) "Switch isoldalbarea > 0 is deprecated and should not"
            write(*,*) "be used. The option isoldalbarea > 0 will be removed" 
            write(*,*) "from future versions of UEDGE."
            write(*,*) "Please set isoldalbarea = 0 "
      endif
      if (oldseec .gt. 0) then
            write(*,*) ""
            write(*,*) ""
            write(*,*) "        **** WARNING ****"
            write(*,*) "Using old, deprecated seec model"
            write(*,*) "Set switch oldseec = 0 to use new model "
            write(*,*) "The old  model will be removed from"
            write(*,*) "future versions of UEDGE"
            write(*,*) "Please set oldseec = 0 "
            write(*,*) ""
      endif
      if (jhswitch .gt. 0) then
            write(*,*) ""
            write(*,*) ""
            write(*,*) "           **** WARNING ****"
            write(*,*) "Switch jhswitch > 0 is deprecated and should not"
            write(*,*) "be used. The option jhswitch > 0 will be removed" 
            write(*,*) "from future versions of UEDGE."
            write(*,*) "Please set jhswitch = 0 "
            write(*,*) ""
      endif
      if ((cftiexclg .gt. 1e-10) .and. (istgon(1) .eq. 1)) then
            write(*,*) ""
            write(*,*) ""
            write(*,*) "           **** WARNING ****"
            write(*,*) "The gas equation (istgon) for atoms is turned on"
            write(*,*) "while the switch cftiexclg>0, which accounts for"
            write(*,*) "the atomic energy in the ion energy equation, "
            write(*,*) "resulting in double-accounting for the atomic"
            write(*,*) "energy. "
            write(*,*) ""
            write(*,*) "Please change cftiexclg=0 when using istgon for"
            write(*,*) "atoms. (The scale factor can be changed gradually)."
            write(*,*) ""
      else if ((cftiexclg .ne. 1.0) .and. (istgon(1) .eq. 0)) then
            write(*,*) ""
            write(*,*) ""
            write(*,*) "           **** WARNING ****"
            write(*,*) "The gas equation (istgon) for atoms is turned off"
            write(*,*) "while the switch cftiexclg!=0, which accounts for"
            write(*,*) "the atomic energy in the ion energy equation, "
            write(*,*) "resulting in a discrepancy for the atomic"
            write(*,*) "energy. "
            write(*,*) ""
            write(*,*) "Please change cftiexclg=1 when not using a separate"
            write(*,*) "atom energy equation. "
            write(*,*) ""
      endif

c...    TODO: checks used to be on nigmx, a local parameter set in
c...        pandf while these checks were located there. Checks moved
c...        here to save time, nigmx obsolete. Probably need a more
c...        elegant way to control species array sizes...
      if (ngsp > 100 .or. nisp > 100) then
         call xerrab("***PANDF in oderhs.m: increase nigmx, recompile")
      endif

c... Roadblockers for  call to pandf through openmp structures (added by J.Guterl)
      if (
     .  (isimpon.gt.0 .and. 
     .      ((isimpon.ne.6) .and. (isimpon.ne.2) .and. 
     .          (isimpon.ne.7))
     .  ) .and. (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)) then
          write(*,*) "Only isimpon=0, 2, 6, or 7 is validated with"
          write(*,*) "openmp. Contact the UEDGE team to use other"
          write(*,*) "options with openmp."
        call xerrab("Invalid isimpon")
      endif

      if (
     .      (ismcnon.gt.0) .and. 
     .      (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)
     .  ) then
            call xerrab('Only ismcnon=0 is validated with openmp.
     .      Contact the UEDGE team to use other options with openmp.')
      endif

      if (
     .      (ishosor.gt.0) .and. 
     .      (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)
     .  ) then
        write(*,*) "Only ishosor=0 is validated with openmp. Contact"
        write(*,*) "the UEDGE team to use other options with openmp."
            call xerrab("Invalid ishosor")
      endif

      if (
     .      (ispsorave.gt.0) .and. 
     .      (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)
     .  ) then
        write(*,*) "Only ispsorave=0 is validated with openmp. Contact"
        write(*,*) "the UEDGE team to use other options with openmp."
            call xerrab("Invalid ispsorave")
      endif
   
      if (
     .      ((cfvisxneoq.ne.0) .or. (cfvisxneov.ne.0) ).and. 
     .      (ParallelJac.gt.0 .or. ParallelPandf1.gt.0)
     .  ) then
        write(*,*) "Only cfvisxneoq=cfvisxneov=0 is validated with "
        write(*,*) "openmp. Contact the UEDGE team to use other "
        write(*,*) "options with openmp."
            call xerrab("Invalid cfvisxneoq/cfvisxneov")
      endif






         ijac = 0
         nurlxn = cnurn*nurlx
         nurlxu = cnuru*nurlx
         nurlxe = cnure*nurlx
         nurlxi = cnuri*nurlx
         nurlxg = cnurg*nurlx
         nurlxp = cnurp*nurlx

*     -- For the continuation mode (icntnunk=1), be sure a Jacobian was
*     -- calculated on the previous step, i.e., ijac > 0
         if (icntnunk==1 .and. ijactot<=1) then
        call xerrab('**Error: need initial Jacobian-pair for icntnunk=1')
         endif

c     -- Reinitialize ijactot if icntnunk = 0; prevents ijactot=2 by 2 exmain
c     .. nksol issue
         if (icntnunk == 0) ijactot = 0

c     -- call principal driver routine --
         call uedriv

c  -- create the interpolants for the grid sequencing or restarting --
c  -- only update save variables (call gridseq) for nksol if a root
c  -- has been found, i.e., if iterm=1; other convergence failures never
c  -- reach here (kaboomed).


         #fill bndry flux arrays for export 
         if (isoutwall==1) call outwallflux  

         if (isnintp .eq. 0 ) then
            nxold = nx
            nyold = ny
            call gchange("Interp",0)
            if ( iterm .eq. 1) then
               call gridseq
               call comp_vertex_vals  # gen plasma/neut values at rm,zm(,,4)
            endif
         elseif (isnintp .eq. 1 ) then
            if ( iterm .eq. 1) then
               nxold = nx
               nyold = ny
                  call gchange("Interp",0)
                  call gridseq
                  call comp_vertex_vals  # gen plasma/neut values at rm,zm(,,4)
            endif
         endif



      return
      end
c --** End of subroutine exmain *********************************
c-----------------------------------------------------------------------
      subroutine exmain_prelims

*     EXMAIN_PRELIMS does some preliminary initialization for EXMAIN.
*     This has been separated to make it easier to duplicate EXMAIN
*     in parser code.

      implicit none

      Use(Err_msg_out)   # errmsgflag,errunit
      Use(Dim)           # nx,ny,nisp,ngsp (in UEpar)
      Use(UEpar)         # label

      character*8 cdate, ctime, rtime, rdate, rmach

c  -- set output label --

ccc  Commented out since Basis11 & Basis12 uses diff arguments
ccc      call glbheadi(cdate,ctime,rtime,rdate,rmach)
ccc      label = 'UEDGE CODE, run at '//rtime//' on '//rdate//
ccc     .   ', machine '//rmach

*  -- set flag and output unit for error messages --

      call xsetfp (errmsgflag)
      call xsetunp (errunit)

      return
      end
c------** End of subroutine exmain_prelims----------------------------
!----------------------------------------------------------c
      subroutine onedconteq
!
!  Solves 1D convection/diffusion continue eqn time-dependently
! 
!         dn/dt + div(nv -Ddn/dx) = Sp
!
! This equation is solved on the spatial domain from x=0 to x= 1 with
! boundary conditons dn/dx = 0 at x=0 and n=1 (normalized) at x=1.
! The coefficients for v and D are constructed to yield the same
! steady-state solution as for a constant D. The parameter alfz >= 1
! controls the v and D; if alfz=1, D=const, and for alf>1, there is an
! inward convection (pinch) peaking near x=1 and an enhanced diffusion
! there. Use upwinding on the convective term.

      implicit none
      Use(Convdiffeqn)

!...  Local variables
      integer ix,it
      real tim,delto,delt,delx,dtdx,vrzmax,drzmax

      call gchange("Convdiffeqn",0)

!... Compute mesh, vr, and dr
      do ix = 1, nxx
        xcz(ix) = float(ix-1)/float(nxx-1)
        xfz(ix) = xcz(ix) + 0.5/float(nxx-1)
        vrz(ix) = -vrfac*(alfz-1.)*sp*xfz(ix)**2
        drz(ix) = 1 + (alfz-1.)*xfz(ix)*(0.5*sp*(1-xfz(ix)**2)+ 1)
      enddo

!... Compute the timestep using the Courant condition
      vrzmax = 0.
      drzmax = 0.
      do ix = 1, nxx
        vrzmax = max(abs(vrz(ix)),vrzmax) + 1.e-20
        drzmax = max(abs(drz(ix)),drzmax) + 1.e-20
      enddo
      delx = 1./float(nxx-1)
      delt = courant*min(delx/vrzmax, 0.5*delx**2/drzmax)
      write(*,*) 'delt = ',delt
      dtdx = delt/delx
      delto = tendoned/float(ntim)
      ito = 1
      timo(1) = 0.

!... Set initial conditions
      do ix = 1, nxx
        dens(ix) = 1.
        nnt(ix,1) = dens(ix)
      enddo

!... Form finite volume eqn and advance in time
      do it = 1, ndtmax
        dens(1) = dens(2)  # Gives zero flux for finite volume method
cc        dens(1) = (4.*dens(2) - dens(3))/3.  # 2nd order Boundary condition
        dens(nxx) = 1.
        tim = tim + delt
        do ix = 2, nxx-1
          vrhs(ix) = vrz(ix  )*0.5*(dens(ix+1)+dens(ix  )) -
     .               vrz(ix-1)*0.5*(dens(ix  )+dens(ix-1))
          drhs(ix) = -( drz(ix  )*(dens(ix+1)-dens(ix  )) - 
     .                  drz(ix-1)*(dens(ix  )-dens(ix-1)) )/delx
          gampz(ix) = vrz(ix)*0.5*(dens(ix+1)+dens(ix)) - 
     .                drz(ix)*(dens(ix+1)-dens(ix  ))/delx
          dens(ix) = dens(ix) + (-vrhs(ix) - drhs(ix) + sp*delx)*dtdx
        enddo
        if (tim > timo(ito)+delto .and. ito<ntim) then   # store solution
          dens(1) = dens(2)  # Gives zero flux for finite volume method
cc          dens(1) = (4.*dens(2) - dens(3))/3.     # Update boundary vals
          dens(nxx) = 1.
          ito = ito + 1
          timo(ito) = tim
          do ix = 1, nxx
            nnt(ix,ito) = dens(ix)
            gampzt(ix,ito) = gampz(ix)
          enddo
        endif
        if (tim > tendoned) exit
      enddo
      return
      end

