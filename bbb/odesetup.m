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
      Use(Share)    # nycore,nysol,nxleg,nxcore,nxomit,igrid,isgrdsym,
                    # nyomitmx,geometry
      Use(UEpar)    # svrpkg,,isnion,isupon,isupgon,isteon,
                    # istion,isngon,isphion
      Use(UEint)    # minu,ziin,newgeo,mhdgeo,isallloc
      Use(Lsode)    # neq,jacflg,jpre,ipar,ismmaxuc,mmaxu
      Use(Solver_work_arrays)   # liw,lrw,iwork
      Use(Grid)     # inewton
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
      Use(Parallv)             # nxg,nyg
      Use(Jac_work_arrays)     # lwp,liwp
      Use(Indices_domain_dcg)  # isddcon
      Use(Indices_domain_dcl)  # vrsendl,visendl,ixmxbcl,ixmnbcl,iymxbcl,iymnbcl
      Use(Indices_loc_glob_map)# ivcum,ivloc2sdg,ivloc2mdg,ivl2gstnl
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
      call gallot("Grid",0)
      call gallot("Stat",0)

c...  reset inewton to unity if a Newton iteration is used
      inewton(igrid) = 0
      if((svrpkg.eq.'nksol') .or. (svrpkg.eq.'petsc') .or. 
     .                        (svrpkg.eq.'newton')) inewton(igrid) = 1
c...  Be sure nufak and dtreal are large if this is time-dependent run
      if (inewton(igrid) .eq. 0) then
         nufak = 1.e-20
         dtreal = 1.e20
      endif

c ... Set mesh dimensions and x-point parameters
      if (isallloc .eq. 0) then       # skip if allocating local proc
      nxg = nx                        # only if global allocation
      nyg = ny
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
               iysptrx1(1) = nycore(igrid)
               iysptrx2(1) = nycore(igrid)
               iysptrx = nycore(igrid)
               ixlb(1) = 0
               ixpt1(1) = nxleg(igrid,1) + nxxpt
               ixpt2(1) = nxleg(igrid,1) + nxcore(igrid,1) +
     .                                     nxcore(igrid,2) + 3*nxxpt
               ixrb(1) = ixpt2(1) + nxleg(igrid,2) + nxxpt
               if (nxomit .gt. 0) then
                  ixpt1(1) = ixpt1(1) - nxomit
                  ixpt2(1) = ixpt2(1) - nxomit
                  ixrb(1) = ixrb(1) - nxomit
               endif
               if (nyomitmx >= nysol(igrid)) then  # special case for core only
                  ixmnbcl = 0
                  ixmxbcl = 0
               endif
            endif
c----------------------------------------------------------------------c
         elseif (mhdgeo == 2) then   #toroidal geo with circular cross-sec
c----------------------------------------------------------------------c
            iysptrx1(1) = nycore(igrid)
            iysptrx2(1) = nycore(igrid)
            iysptrx = nycore(igrid)
            ixlb(1) = 0
            ixpt1(1) = nxleg(igrid,1)
            ixpt2(1) = nxleg(igrid,1) + nxcore(igrid,1) +
     .                                  nxcore(igrid,2)
            ixrb(1) = ixpt2(1) + nxleg(igrid,2)
            ixmnbcl = 0
            ixmxbcl = 0
c----------------------------------------------------------------------c
	 else	# cases mhdgeo=0 (cyl), mhdgeo=-1 (slab), mhdgeo=-2 (mirror)
c----------------------------------------------------------------------c
            iysptrx1(1) = nycore(igrid)
            iysptrx2(1) = nycore(igrid)
            iysptrx = nycore(igrid)
            ixlb(1) = 0
            ixpt1(1) = -1
            ixpt2(1) = nxcore(igrid,2)
            ixrb(1) = ixpt2(1) + nxleg(igrid,2) + nxxpt
            if (isgrdsym .eq. 1) then
               ixpt1(1) = (nxm - ixpt2(1))/2
               ixpt2(1) = (nxm + ixpt2(1))/2
            endif
            if (nyomitmx >= nysol(igrid)+nyout(1)) then
               ixmnbcl = 0
               ixmxbcl = 0
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
               iysptrx1(1) = nycore(igrid)
               iysptrx2(1) = nycore(igrid)
               iysptrx = nycore(igrid)
               ixlb(1) = 0
               ixpt1(1) = nxleg(igrid,1) + nxxpt
               ixpt2(1) = nxleg(igrid,1) + nxcore(igrid,1) +
     .                                     nxcore(igrid,2) + 3*nxxpt
               ixrb(1) = ixpt2(1) + nxleg(igrid,2) + nxxpt
               if (nxomit .gt. 0) then
                  ixpt1(1) = ixpt1(1) - nxomit
                  ixpt2(1) = ixpt2(1) - nxomit
                  ixrb(1) = ixrb(1) - nxomit
               endif
               if (nyomitmx >= nysol(igrid)) then  # special case for core only
                  ixmnbcl = 0
                  ixmxbcl = 0
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
      numvarl = numvar   # for parallel mpi version
      
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
      
      if (svrpkg .eq. "cvode" .or. svrpkg .eq. "kinsol") then
        neqmx = neq      #cvode, kinsol allow no extra storage for yl
      else
        neqmx = neq + 2  #1st extra flag for Jac. calc; 2nd nufak=/1psuedo dt
      endif

      if (svrpkg .eq. "cvode" .or. svrpkg .eq. "kinsol") neqmx = neq
                                  #cvode, kinsol allow no extra storage
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
      if (jpre+jacflg+inewton(igrid) .eq. 0) lda = 1
c...  Increase of mmaxu with problem size is empirical
      if (ismmaxuc .eq. 1) mmaxu = neq**0.5
      if (ismmaxuc.eq.0 .and. icntnunk.eq.1) then
         call remark('WARNING: ismmaxuc=0 & icntnunk=1; Jac storage ok?')
      endif
c check to make sure premeth, inewton(igrid) and svrpkg do not conflict.
c when inewton(igrid)=1 and svrpkg .ne. 'nksol', then premeth must equal
c 'banded'.

      if (inewton(igrid) .eq. 1) then
         if (svrpkg .ne. 'nksol' .and. svrpkg.ne.'petsc') then
            if ((premeth .eq. 'ilut') .or. (premeth .eq. 'inel')) then
              write(STDOUT,*) "*** Invalid option:  premeth=",premeth
              write(STDOUT,*) "When using standard Newton, premeth ",
     .                        "must equal banded."
              call xerrab("")
            endif
         endif
      endif

c ... Set maximum lengths for preconditioner arrays, and allocate
c     memory for auxiliary arrays for desired storage-format option,
c     if necessary.
      if (lenpfac .lt. 9*numvar) then
        nnzmx = 9*numvar*neq
      else
        nnzmx = lenpfac*neq
      endif
      if (premeth .eq. 'ilut') then
         lwp = nnzmx + lenplufac*neq   # extra space for fill in
         lenplumx = lwp
         liwp = max(lwp+neq, nnzmx+neq+1)   # see jac_lu_decomp
c save preconditioner reordering info for icntnunk=1 case.
         if(ireorder .eq. 1) call gchange("Jacreorder",0)
      elseif (premeth .eq. 'banded') then
         lwp = lda*neq
         liwp = 3+neq
      elseif (premeth .eq. 'inel') then
         lwp = (neq+1)*ndiagmx
         liwp = 2+ndiagmx
         call gallot("Nonzero_diagonals",0)
      else
         write(STDOUT,*) "*** Invalid option:  premeth=",premeth
         call xerrab("")
      endif

c ... Set lengths of work arrays for desired solver option.
      if (svrpkg.eq.'newton') then
         lrw = lwp
         liw = liwp
      elseif (svrpkg.eq.'vodpk') then
         lrw = 17*neq + 61 + lwp
         if (inopt.eq.1 .and. iworkin(8).gt.0) then
            lrw = lrw + (neq + iworkin(8))*(iworkin(8) + 3) -
     .                  (neq + 5)*(5 + 3)
         endif
         liw = 35 + liwp
      elseif (svrpkg.eq.'cvode') then
         # nothing to do for cvode, but this validates cvode as solver
      elseif (svrpkg.eq.'daspk') then
         if (iscolnorm .ne. 0) then
            call xerrab('**** Set iscolnorm=0 for DASPK, then restart *')
         endif
         if (info(12).eq.0) then   # direct solution method with
                                   # banded Jacobian supplied in jacd1
            lrw = 50+(2*lbw+ubw+10)*neq
            if (info(5).eq.0) lrw = lrw + 2*(neq/(lbw+ubw+1)+1)
            liw = 40+neq
         else                      # Krylov method
            lrw = 18*neq+91+lwp
            if (info(13) .eq. 1) then
               lrw = lrw + (neq + iwork(24))*(iwork(24) + 3) -
     .                     (neq + 5)*(5 + 3)
            endif
            liw = 40+liwp
         endif
         if (info(10).eq.1 .or. info(10).eq.3) liw = liw + neq
         if (info(11).eq.1 .or. info(16).eq.1) liw = liw + neq
         if (info(16).eq.1) lrw = lrw + neq
      elseif (svrpkg.eq.'nksol' .or. svrpkg.eq.'petsc') then
cc     sizes are set for mf=1 case; opt. input iwork(1)=mmaxu needed
cc     do not reset if icntnunk=1 (saved Jac) as storage changes if mmaxu reset
         lenk = 6 + 4*neq + (neq+1)*mmaxu + 2*mmaxu*(mmaxu+1)
         lrw = 4 + 4*neq + lenk + lwp
         liw = 20 + mmaxu + liwp
      else
         write(STDOUT,*) "*** Invalid option:  svrpkg=",svrpkg
         call xerrab("")
      endif

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
      if (isddcon .eq. 1) call gallot("Indices_domain_dcl",0)
      if (isallloc .eq. 0) then
         call gallot("Indices_domain_dcg",0)
         call gallot("Indices_loc_glob_map",0)
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
      call gallot("Save_terms",0)
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

c ... Load optional inputs available only in this routine.
cc      write(*,*) 'Just before nksol,iwork', '  iwork =', iwork
cc      if (svrpkg .eq. 'nksol') then
cc         iwork(3) = lwp
cc         iwork(4) = liwp
cc      elseif (svrpkg .eq. 'vodpk') then
cc         iwork(1) = lwp
cc         iwork(2) = liwp
cc      elseif (svrpkg .eq. 'daspk' .and. info(12) .ne. 0) then
cc         iwork(27) = lwp
cc         iwork(28) = liwp
cc      endif

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
      Use(Parallv)    # neq
#      Use(Indices_domain_dcl)  # ixmxbcl,ixmnbcl,iymxbcl,iymnbcl

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
      Use(Indices_domain_dcl)  # ixmxbcl,ixmnbcl,iymxbcl,iymnbcl

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
      do iy = 1-iymnbcl, ny+iymxbcl   #iymb,xbcl =1(0) if guard cell(no)
	do ix = 1-ixmnbcl, nx+ixmxbcl #ixmb,xbcl =1(0) if guard cell(no)
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
      iysptrx1(1) = nycore(igrid)
      iysptrx2(1) = nycore(igrid)
      iysptrx = nycore(igrid)
      iysptrx1(2) = iysptrx2(1)
      iysptrx2(2) = iysptrx1(1)
      ixlb(1)  = 0
      ixpt1(1) = nxleg(igrid,1) + nxxpt
      ixmdp(1) = ixpt1(1) + nxcore(igrid,1) - 1 + nxxpt
      ixpt2(1) = ixmdp(1) + nxcore(igrid,1) - 1 + nxxpt
      ixrb(1)  = ixpt2(1) + nxleg(igrid,1) + nxxpt
      ixlb(2)  = ixrb(1) + 2
      ixpt1(2) = ixlb(2) + nxleg(igrid,2) + nxxpt
      ixmdp(2) = ixpt1(2) + nxcore(igrid,2) - 1 + nxxpt
      ixpt2(2) = ixmdp(2) + nxcore(igrid,2) - 1 + nxxpt
      ixrb(2)  = ixpt2(2) + nxleg(igrid,2) + nxxpt

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
      iysptrx1(1) = nycore(igrid)
      iysptrx2(1) = nycore(igrid)
      iysptrx = nycore(igrid)
      iysptrx1(2) = iysptrx2(1)
      iysptrx2(2) = iysptrx1(1)
      ixlb(1)  = 0
      ixpt1(1) = nxleg(igrid,1) + nxxpt
cc      ixmdp(1) = ixpt1(1) + nxcore(igrid,1) - 1 + nxxpt
      ixpt2(1) = ixpt1(1) + nxcore(igrid,1) - 1 + nxxpt
      ixrb(1)  = ixpt2(1)
      ixlb(2)  = ixrb(1) + 2
      ixpt1(2) = ixlb(2)
cc      ixmdp(2) = ixpt1(2) + nxcore(igrid,2) - 1 + nxxpt
      ixpt2(2) = ixlb(2) + nxcore(igrid,2) - 1 + nxxpt
      ixrb(2)  = ixpt2(2) + nxleg(igrid,2) + nxxpt

      return
      end
c-----------------------------------------------------------------------
      subroutine ueinit

*     UEINIT defines the geometry, initializes the state of the
*     plasma, and defines miscellaneous other fields.

      implicit none

      Use(Dim)      # nx,ny,nhsp,nusp,nzspt,nisp,ngsp,nxpt
      Use(Share)    # nxomit,igrid,geometry,isnonog,nyomitmx
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
                    # svrpkg,isphiofft,methg
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
      Use(Parallv)  # nxg,nyg
      Use(Imprad)   # isimpon,ismctab
      Use(Impurity_source_flux)   # fnzysi,fnzyso
      Use(Variable_perturbation)  # del
      Use(Time_dep_nwt)           # nufak,inufaknk
      Use(Indices_domain_dcg)     # ndomain
      Use(Indices_domain_dcl)     # ndomain
      Use(Conduc)   # nuvl
      Use(Locflux)  # fgtdx, fgtdy
      Use(Jacobian) # isjacstnlon
      Use(Rccoef)   # feixextlb,rb;feiyexti,o
      Use(Rhsides)  # psorc, psorxr, msor, msorxr
      Use(Save_terms) # psorcold, etc
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
      if (svrpkg.ne.'cvode'.and.svrpkg.ne.'kinsol') then #not allowed for cvode
         yl(neq+1) = -1.
         if (inufaknk .eq. 1) then   # deter. if nufak is used in Krylov step
            yl(neq+2) = nufak
         else
            yl(neq+2) = 0.
         endif
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

c ... Set initial values of time-step arrays if svrpkg=nksol
      if ((svrpkg .eq. "nksol") .or. (svrpkg.eq."petsc")) then
        call s2fill (nx+2, ny+2, 1.e20, dtoptx, 1, nx+2)
        call sfill (neq, 1.e20, dtoptv(:), 1)
        call sfill (neq, 1.e20, dtuse(:), 1)
      endif

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
      srtolpk = del / rtolv(igrid)

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

c...  Set boundary conditions for lyni on walls; if isulynix=0, use lyni
ccTR      if (isulynix == 0) then
ccTR       do ifld = 1, nisp
ccTR          do iu = 0, nx+1
ccTR            lynix(1,iu,ifld)=lyni(1)
ccTR            lynix(2,iu,ifld)=lyni(2)
ccTR          enddo
ccTR        enddo
ccTR      endif
ccTR     if (isulytex == 0) then
ccTR       do iu = 0, nx+1
ccTR          lytex(1,iu) = lyte(1)
ccTR          lytex(2,iu) = lyte(2)
ccTR        enddo
ccTR      endif
ccTR      if (isulytix == 0) then
ccTR        do iu = 0, nx+1
ccTR          lytix(1,iu) = lyti(1)
ccTR          lytix(2,iu) = lyti(2)
ccTR        enddo
ccTR      endif
ccTR      if (isulyphix == 0) then
ccTR        do iu = 0, nx+1
ccTR          lyphix(1,iu) = lyphi(1)
ccTR          lyphix(2,iu) = lyphi(2)
ccTR        enddo
ccTR      endif

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
c ... Set up tables for impurity atomic-physics processes.
      if (isimpon .eq. 1) then		# obsolete option
         call xerrab ('ueinit -- option isimpon=1 is obsolete; use 2')
      elseif (isimpon .eq. 2) then	# data supplied by D. Post 1993
         call readpost
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
         call readpost
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
         if (iysptrx >= 1 .and. nxleg(1,1)+nxleg(1,2) > 0) then 
                                          # need sep for this calc of lconi,e
           if(ndomain <=1) call conlen    #if parall, compute only at setup
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
         if (ixmnbcl==1) fgtdx(ixlb(jx)) = gcfacgtx
	 if (ixmxbcl==1) fgtdx(ixrb(jx)) = gcfacgtx
      enddo
      if (iymnbcl==1) fgtdy(0)  = gcfacgty
      if (iymxbcl==1) fgtdy(ny) = gcfacgty

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
      if (ig .eq. 1 .and. restart .eq. 0) then
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

       if (svrpkg.eq.'daspk') then
	 call pandf1(-1,-1,0,ipar(1),tv,yl,yldot)
         do 730 ifld = 1, nhsp
 	 do 700 ix = 0, nx+1
	   if(isnionxy(ix,0,   ifld)==1) yldot(idxn(ix,0,   ifld))=0.
	   if(isnionxy(ix,ny+1,ifld)==1) yldot(idxn(ix,ny+1,ifld))=0.
	   if(isuponxy(ix,0,   ifld)==1) yldot(idxu(ix,0,   ifld))=0.
	   if(isuponxy(ix,ny+1,ifld)==1) yldot(idxu(ix,ny+1,ifld))=0.
	   if(isteonxy(ix,0   )==1) yldot(idxte(ix,0   )) = 0.
	   if(isteonxy(ix,ny+1)==1) yldot(idxte(ix,ny+1)) = 0.
	   if(istionxy(ix,0   )==1) yldot(idxti(ix,0   )) = 0.
	   if(istionxy(ix,ny+1)==1) yldot(idxti(ix,ny+1)) = 0.
           do 698 igsp = 1, ngsp
             if(isngonxy(ix,0,   igsp)==1) yldot(idxg(ix,0,   igsp))=0.
	     if(isngonxy(ix,ny+1,igsp)==1) yldot(idxg(ix,ny+1,igsp))=0.
  698      continue

  700    continue
         do 701 iy = 1, ny
	   if(isnionxy(0,iy,   ifld)==1) yldot(idxn(0,iy,   ifld)) = 0.
	   if(isnionxy(nx+1,iy,ifld)==1) yldot(idxn(nx+1,iy,ifld)) = 0.
	   if(isuponxy(0,iy,   ifld)==1) yldot(idxu(0,   iy,ifld)) = 0.
	   if(isuponxy(nx+1,iy,ifld)==1) yldot(idxu(nx+1,iy,ifld)) = 0.
	   if(isteonxy(0,   iy)==1) yldot(idxte(0,   iy)) = 0.
	   if(isteonxy(nx+1,iy)==1) yldot(idxte(nx+1,iy)) = 0.
	   if(istionxy(0,   iy)==1) yldot(idxti(0,   iy)) = 0.
	   if(istionxy(nx+1,iy)==1) yldot(idxti(nx+1,iy)) = 0.
           do 699 igsp = 1, ngsp
             if(isngonxy(0,   iy,igsp)==1) yldot(idxg(0,   iy,igsp)) = 0.
	     if(isngonxy(nx+1,iy,igsp)==1) yldot(idxg(nx+1,iy,igsp)) = 0.
  699      continue
  701    continue
  730    continue
         do ifld = 1, nzspt
           do ix = 0, nx+1
	     if (isnionxy(ix,0,nhsp+ifld)==1)
     .                             yldot(idxn(ix,0,nhsp+ifld)) = 0.
	     if (isnionxy(ix,ny+1,nhsp+ifld)==1)
     .                             yldot(idxn(ix,ny+1,nhsp+ifld)) = 0.
           enddo
           do iy = 1, ny
	     if (isnionxy(0,iy,nhsp+ifld)==1)
     .                             yldot(idxn(0,iy,nhsp+ifld)) = 0.
	     if (isnionxy(nx+1,iy,nhsp+ifld)==1)
     .                             yldot(idxn(nx+1,iy,nhsp+ifld)) = 0.
           enddo
         enddo
         do 702 iy = 0, ny+1
            do 703 ix = 0, nx+1
	       if(isphionxy(ix,iy)==1) yldot(idxphi(ix,iy)) = 0.
  703        continue
  702    continue
       endif   #ends if(svrpkg.eq.'daspk')
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
                   lng(ix,iy,igsp) = log(ng(ix,iy,igsp))
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
      if (nyomitmx >= nysol(1)) then   # no inner/outer leg region
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
      if (nysol(1) <= nyomitmx) then  #core only, no divertor legs
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
                  if (ineudif .eq. 2) then
                    lng(ix,iy,igsp) = log(ng(ix,iy,igsp))
                  endif
               endif
             enddo
           enddo
         enddo

      endif          # end of very-large 2-branch-if: (1), same mesh size
                     # or (2), index-based interp with isnintp=1 

      if (nyomitmx >= nysol(1)+nyout(1)) then
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
      if (nyomitmx >= nysol(1)+nyout(1)) then
        call filldead_guardcells
      endif
         
      call convert
      if (svrpkg.eq.'daspk') then
	 call pandf1(-1,-1,0,ipar(1),tv,yl,yldot)
 	 do 800 ix = 0, nx+1
           do ifld = 1, nisp
             if(isnionxy(ix,0,ifld)==1) yldot(idxn(ix,0,ifld)) = 0.
             if(isnionxy(ix,ny+1,ifld)==1) yldot(idxn(ix,ny+1,ifld))=0.
           enddo
           do ifld = 1, nusp
             if(isuponxy(ix,0,ifld)==1) yldot(idxu(ix,0,ifld)) = 0.
             if(isuponxy(ix,ny+1,ifld)==1) yldot(idxu(ix,ny+1,ifld))=0.
           enddo
           if(isteonxy(ix,0)==1) yldot(idxte(ix,0)) = 0.
           if(isteonxy(ix,ny+1)==1) yldot(idxte(ix,ny+1)) = 0.
           if(istionxy(ix,0)==1) yldot(idxti(ix,0)) = 0.
           if(istionxy(ix,ny+1)==1) yldot(idxti(ix,ny+1)) = 0.
           do 798 igsp = 1, ngsp
             if(isngonxy(ix,0,igsp)==1) yldot(idxg(ix,0,igsp)) = 0.
             if(isngonxy(ix,ny+1,igsp)==1) yldot(idxg(ix,ny+1,igsp))=0.
  798      continue

  800    continue
         do 801 iy = 1, ny
           do ifld = 1, nisp
             if(isnionxy(0,iy,ifld)==1) yldot(idxn(0,iy,ifld)) = 0.
             if(isnionxy(nx+1,iy,ifld)==1) yldot(idxn(nx+1,iy,ifld))=0.
           enddo
           do ifld = 1, nusp
             if(isuponxy(0,iy,ifld)==1) yldot(idxu(0,iy,ifld)) = 0.
             if(isuponxy(nx+1,iy,ifld)==1) yldot(idxu(nx+1,iy,ifld))=0.
           enddo
           if(isteonxy(0,iy)==1) yldot(idxte(0,iy)) = 0.
           if(isteonxy(nx+1,iy)==1) yldot(idxte(nx+1,iy)) = 0.
           if(istionxy(0,iy)==1) yldot(idxti(0,iy)) = 0.
           if(istionxy(nx+1,iy)==1) yldot(idxti(nx+1,iy)) = 0.
           do 799 igsp = 1, ngsp
             if(isngonxy(0,iy,igsp)==1) yldot(idxg(0,iy,igsp)) = 0.
             if(isngonxy(nx+1,iy,igsp)==1) yldot(idxg(nx+1,iy,igsp))=0.
  799      continue
  801    continue
         do 802 iy = 0, ny+1
 	    do 803 ix = 0, nx+1
 	      if(isphionxy(ix,iy)==1) yldot(idxphi(ix,iy)) = 0.
  803       continue
  802    continue
      endif

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
Use(Parallv)            # nxg,nyg
Use(Indices_domain_dcg) # ndomain,isddcon
Use(Npes_mpi)           # mype
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
         if (ndomain==1 .or. isddcon==2 .or. isglobal==1) then
                                             # no x-domain decomposition
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
               elseif (mype<=0.and.iysptrx1(1)==0) then #phys bndry at sep
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ix+1
               endif
            enddo  # first end do-loop over ix
         elseif (ndomain > 1) then    # all cells are local
           if (isddcon == 1) then			 
             do ix = 1, nx
                ixm1(ix,iy) = ix-1
                ixp1(ix,iy) = ix+1
             enddo  #second loop over ix
           elseif (isddcon == 2) then
             do ix = 1, nx
               if (iy <= iysptrxg(mype+1) .and. ix==ixpt1(1)) then
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ixpt2(1)+1
               elseif (iy <= iysptrxg(mype+1) .and. ix==ixpt1(1)+1) then
                  ixm1(ix,iy) = ixpt2(1)
                  ixp1(ix,iy) = ix+1
               elseif (iy <= iysptrxg(mype+1) .and. ix==ixpt2(1)) then
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ixpt1(1)+1
               elseif (iy <= iysptrxg(mype+1) .and. ix==ixpt2(1)+1) then
                  ixm1(ix,iy) = ixpt1(1)
                  ixp1(ix,iy) = ix+1
               else
                  ixm1(ix,iy) = ix-1
                  ixp1(ix,iy) = ix+1
               endif
             enddo  # 3rd loop over ix
           endif
         endif
         do jx = 1, nxpt  #fix poloidal bdrys for multiple nulls
           ixm1(ixlb(jx),iy) = ixlb(jx)
           ixp1(ixlb(jx),iy) = ixlb(jx)+1
           ixm1(ixrb(jx)+1,iy) = ixrb(jx)
           ixp1(ixrb(jx)+1,iy) = ixrb(jx)+1
         enddo
      enddo  # end do-loop over iy

c...  Special fix for core-only cases
      if (nyomitmx >= nysol(1)) then
         ixm1(1,ny+1) = ixpt2(1)
         ixp1(nx,ny+1) = 1
      endif

      return
      end  # end of subroutine set_indirect_address

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c

      subroutine globalvars

c ... subroutine to transfer initial plasma and gas variables to
c ... the global arrays nisg, upsg, etc. to domain decomposition

      implicit none

      Use(Dim)           # nx,ny,nhsp,nusp,nzspt,nisp,ngsp
      Use(Compla)        # ni,up,te,ti,ng
      Use(Interp)        # nxold,nyold,afracs
      Use(Global_vars)   # nisg,upsg,tesg,tisg,ngsg,afracsg
      Use(Imprad)        # isimpon

c **  local variables
      integer ifld, iy, ix, igsp

          do ifld = 1, nisp
             do iy = 0, ny+1
                do ix = 0, nx+1
                   nisg(ix,iy,ifld) = ni(ix,iy,ifld)
                   upsg(ix,iy,ifld) = up(ix,iy,ifld)
                enddo
             enddo
          enddo
          do iy = 0, ny+1
             do ix = 0, nx+1
                do igsp = 1, ngsp
                   ngsg(ix,iy,igsp) = ng(ix,iy,igsp)
                enddo
                tesg(ix,iy)      = te(ix,iy)
                tisg(ix,iy)      = ti(ix,iy)
                phisg(ix,iy) = phi(ix,iy)
                if (isimpon.gt.0.) then
                  if (afracs(1,1)+afracs(nx,ny).gt.1.e-20)
     .                                    afracsg(ix,iy) = afrac(ix,iy)
                endif
                if (phi(nx-1,ny-1) .eq. 0) phisg(ix,iy) = 40.  #avoid phi=0.
             enddo
          enddo

      return
      end
c****** end of subroutine globalvars *******************************
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
c --------------------------------------------------------------------c
      subroutine domain_dc

c ... Calculates indices and arrays needed for domain decomposition

      implicit none

      Use(Indices_domain_dcg) 	#ndleg,ndxcore,ndycore,ndysol,ndxt,ndyt,ndomain
                                #ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g,isddcon
                                #ixpt1g,ixpt2g,iysptrxg,idxpt
      Use(Parallv)              #nxg,nyg,neqg
      Use(Lsode)                #neq
      Use(Share)  		#nxleg,nxcore,nycore,nysol,isgrdsym
      Use(Dim)                  #nx,ny
      Use(Xpoint_indices)       #ixpt1,ixpt2,iysptrx
      Use(Math_problem_size)    #neqmx,numvar
      Use(Indices_loc_glob_map) #ivl2gstnl(allocated)

c ... Local variables
      integer idx,idy,idt,nysd,nysd1,nysd2,nxsd1,nxsd2,nxsd3,id,iycum,
     .        nxsd,ixcum,iv,invar,ix,iy

c ... Copy global nx and ny to nxg and nyg, so each processor can know 
      nxg = nx
      nyg = ny
      neqg = neq

c ... If isddcon=0, few variables & Jac stencil comp; no domain decomp
      if (isddcon .eq. 0) then
         ixmnbcg(1) = 1
         ixmxbcg(1) = 1
         iymnbcg(1) = 1
         iymxbcg(1) = 1
         ixmin(1) = 1
         ixmax(1) = nx
         iymin(1) = 1
         iymax(1) = ny
         ixpt1g(1) = ixpt1(1)
         ixpt2g(1) = ixpt2(1)
         iysptrxg(1) = iysptrx1(1)
         ndomain = 1
         ispwrbc(1) = 1
         neq_locgmx = neq   # to dimension ivl2gstnl
c ...    Allocate and compute Jacobian stencil
         call allocjacstnl
         return
      endif

c ... Set values
      if (isddcon == 1) then
        ndxt = ndleg(1,1)+ndleg(1,2)+ndxcore(1)   # number of domains in x
      else
        ndxt = 1                          # only 1 x-domain if issdcon=2
      endif
      ndyt = ndycore(1) + ndysol(1)       # number of domains in y
      ndomain = ndxt*ndyt                 # total number of domains
      if (ndomain > ndomainmx) then
        call xerrab("*** Exiting because ndomain > ndomainmx ***")
      endif

c ... Initialize bndry conditions to be of the interior type (domain-domain)
      do id = 1, ndomain
         ixmnbcg(id) = 0
	 ixmxbcg(id) = 0
	 iymnbcg(id) = 0
	 iymxbcg(id) = 0
         ispwrbc(id) = 0      #flag for identifying domain with core ixpt2
      enddo

c ... Number of indices (including guard cells) in each region (left leg,
c ... core, and right leg)
      if (isgrdsym .eq. 0) then
        nxsd1 =int(float(nxleg(1,1)+1)/(ndleg(1,1)+.00001) + 0.999999 )
        nxsd2 =int(float(nxcore(1,1)+nxcore(1,2))/(ndxcore(1)+.00001)+
     .             0.999999 )
        if (nxleg(1,1).eq.0 .or. nxleg(1,2).eq.0) nxsd2 = int( float(
     .       nxcore(1,1)+nxcore(1,2)+1)/(ndxcore(1)+.00001)+0.999999 )
        nxsd3 =int(float(nxleg(1,2)+1)/(ndleg(1,2)+.00001) + 0.999999 )
      elseif (isgrdsym .eq. 1) then
        nxsd1 =int(float(nxleg(1,2)/2+1)/(ndleg(1,1)+.00001) + 0.999999 )
        nxsd2 =int(float(nxcore(1,2))/(ndxcore(1)+.00001) + 0.999999 )
        nxsd3 =int(float(nxleg(1,2)/2+1)/(ndleg(1,2)+.00001) + 0.999999 )
      endif
      nysd1 = int(float(nycore(1)+1)/(ndycore(1)+.00001) + 0.99999 )
      nysd2 = int(float(nysol(1)+1)/(ndysol(1)+.00001) + 0.99999 )
      if (nycore(1) .eq. 0) then
        nysd2 = int(float(nysol(1)+2)/(ndysol(1)+.00001) + 0.99999 )
      endif

c ... Do a special case if isddcon=2: only divides into domains in y-direction

      if (isddcon .eq. 2) then
         ndomain = ndyt
         if (ndomain > ndomainmx) then
           call xerrab("*** Exiting because ndomain > ndomainmx ***")
         endif
         iycum = 0
         do idy = 1, ndycore(1)+ndysol(1)
            if (idy .le. ndycore(1)) then
              nysd = min(nysd1, iysptrx1(1)+1-iycum)
ccc              nysd = min(nysd1, iysptrx1(1)-iycum)
            else
              nysd = min(nysd2, ny-iycum)
            endif
            iycum = iycum + nysd
            ixmin(idy) = 1
            ixmax(idy) = nx
            if (idy.eq.1) then
               iymnbcg(1) = 1
               iycum = iycum - 1
               iymin(1) = 1
               iymax(1) = nysd - 1
            elseif (idy.eq.ndycore(1)+ndysol(1)) then
               iymxbcg(idy) = 1
               iymin(idy) = iycum + 1 - nysd
               iymax(idy) = ny
            else
               iymin(idy) = iycum + 1 - nysd
               iymax(idy) = iycum
            endif
            if (ixmin(idy) .eq. 1 ) ixmnbcg(idy) = 1    #true boundary
            if (ixmax(idy) .eq. nx) ixmxbcg(idy) = 1	#true boundary
            ixpt1g(idy) = ixpt1(1)
            ixpt2g(idy) = ixpt2(1)
            iysptrxg(idy) = iysptrx1(1) - (iymin(idy)-1)
            idym1g(idy) = max(0, idy-1)
            idyp1g(idy) = min(ndomain+1, idy+1)
            idxm1g(idy) = 0
            idxp1g(idy) = nx+1
         enddo
         ispwrbc(1) = 1   #for 1D radial decomp, core bdry always id=1
c ...    Allocate and compute Jacobian stencil
         call allocjacstnl
         return
      endif

c ... Setup min/max index arrays for 2D domains and boundary condition type

c ... Do the left-leg region
      iycum = 0
      do idy = 1, ndycore(1)+ndysol(1)
         if (idy .le. ndycore(1)) then
           nysd = min(nysd1, iysptrx1(1)+1-iycum)
         else
           nysd = min(nysd2, ny-iycum)
         endif
         iycum = iycum + nysd
         ixcum = 0
         do idx = 1, ndleg(1,1)
            idt = idx + ndxt*(idy-1)
            nxsd = min(nxsd1, ixpt1(1)+1-ixcum)
            ixcum = ixcum + nxsd
c ..       First do the ix indices
            if (idx.eq.1) then
               ixmnbcg(idt) = 1
               ixcum = ixcum - 1
               ixmin(idt) = 1
               ixmax(idt) = nxsd - 1
            elseif (idx.eq.ndleg(1,1)) then
               if (ndxcore(1)+ndleg(1,2).eq.0) ixmxbcg(idt) = 1
               ixmin(idt) = ixcum + 1 - nxsd
               ixmax(idt) = ixpt1(1)
            else
               ixmin(idt) = ixcum + 1 - nxsd
               ixmax(idt) = ixcum
            endif
ccc            ixmin(idt) = nxsd1*(idx-1) + 1
ccc            ixmax(idt) = min(nxsd1*idx, ixpt1(1))
c ..       Next do the iy indices
            if (idy.eq.1) then
               iymnbcg(idt) = 1
               iycum = nysd - 1
               iymin(idt) = 1
               iymax(idt) = nysd - 1
            elseif (idy.eq.ndycore(1)+ndysol(1)) then
               iymxbcg(idt) = 1
               iymin(idt) = iycum + 1 - nysd
               iymax(idt) = ny
            else
               iymin(idt) = iycum + 1 - nysd
               iymax(idt) = iycum
            endif
ccc            iymin(idt) = iycum + 1 - nysd
ccc            iymax(idt) = iycum
ccc            if (ixmin(idt) .eq. 1 ) ixmnbcg(idt) = 1 #true boundary
ccc            if (ixmax(idt) .eq. nx) ixmxbcg(idt) = 1	#true boundary
ccc            if (iymin(idt) .eq. 1 ) iymnbcg(idt) = 1	#true boundary
ccc            if (iymax(idt) .eq. ny) iymxbcg(idt) = 1	#true boundary
         enddo
      enddo

c ... Do total core region
      iycum = 0
      do idy = 1, ndycore(1)+ndysol(1)
         if (idy .le. ndycore(1)) then
           nysd = min(nysd1, iysptrx1(1)+1-iycum)
         else
           nysd = min(nysd2, ny-iycum)
         endif
         iycum = iycum + nysd
         ixcum = max(ixpt1(1), 0)
         do idx = 1, ndxcore(1)
            idt = ndleg(1,1) + idx + ndxt*(idy-1)
            nxsd = min(nxsd2, max(ixpt2(1),0)+1-ixcum)
            ixcum = ixcum + nxsd
c ..       First do the ix indices
            if (idx.eq.1 .and. ixpt1(1).le.0) then
               ixmnbcg(idt) = 1
               ixcum = ixcum - 1
               ixmin(idt) = 1
               ixmax(idt) = nxsd - 1
            elseif (idx.eq.ndxcore(1)) then
               if (ndleg(1,2).eq.0) ixmxbcg(idt) = 1
               ixmin(idt) = ixcum + 1 - nxsd
               ixmax(idt) = ixpt2(1)
               if(idy==1) ispwrbc(idt) = 1  #ident idt for core pwr BC
            else
               ixmin(idt) = ixcum + 1 - nxsd
               ixmax(idt) = ixcum
            endif
ccc            ixmin(idt) = max(ixpt1(1), 0) + nxsd2*(idx-1) + 1
ccc            ixmax(idt) = min( max(ixpt1(1),0)+nxsd2*idx, ixpt2(1) )
c ..       Next do the iy indices
            if (idy.eq.1) then
               iymnbcg(idt) = 1
               iycum = nysd - 1
               iymin(idt) = 1
               iymax(idt) = nysd - 1
            elseif (idy.eq.ndycore(1)+ndysol(1)) then
               iymxbcg(idt) = 1
               iymin(idt) = iycum + 1 - nysd
               iymax(idt) = ny
            else
               iymin(idt) = iycum + 1 - nysd
               iymax(idt) = iycum
            endif
ccc            iymin(idt) = iycum + 1 - nysd
ccc            iymax(idt) = iycum
ccc            if (ixmin(idt) .eq. 1 ) ixmnbcg(idt) = 1 #true boundary
ccc            if (ixmax(idt) .eq. nx) ixmxbcg(idt) = 1	#true boundary
ccc            if (iymin(idt) .eq. 1 ) iymnbcg(idt) = 1	#true boundary
ccc            if (iymax(idt) .eq. ny) iymxbcg(idt) = 1	#true boundary
         enddo
      enddo

c ... Do the right-leg region
      iycum = 0
      do idy = 1, ndycore(1)+ndysol(1)
         if (idy .le. ndycore(1)) then
           nysd = min(nysd1, iysptrx1(1)+1-iycum)
         else
           nysd = min(nysd2, ny-iycum)
         endif
         iycum = iycum + nysd
         if (ixpt2(1).gt.0) then
            ixcum = ixpt2(1)
         elseif (ixpt1(1).le.0 .and. ixpt2(1).le.0) then
            ixcum = 0
         endif
         do idx = 1, ndleg(1,2)
            idt = ndleg(1,1) + ndxcore(1) + idx + ndxt*(idy-1)
            nxsd = min(nxsd3, nx-ixcum)
            ixcum = ixcum + nxsd
c ..       First do the ix indices
            if (idx.eq.1 .and. ixpt2(1).le.0) then
               ixmnbcg(idt) = 1
               ixcum = ixcum - 1
               ixmin(idt) = 1
               ixmax(idt) = nxsd - 1
            elseif (idx.eq.ndleg(1,2)) then
               ixmxbcg(idt) = 1
               ixmin(idt) = ixcum + 1 - nxsd
               ixmax(idt) = nx
            else
               ixmin(idt) = ixcum + 1 - nxsd
               ixmax(idt) = ixcum
            endif
ccc            ixmin(idt) = ixpt2(1) + nxsd3*(idx-1) + 1
ccc            ixmax(idt) = min(ixpt2(1)+nxsd3*idx, nx)
            if (idy.eq.1) then
               iymnbcg(idt) = 1
               iycum = nysd - 1
               iymin(idt) = 1
               iymax(idt) = nysd - 1
            elseif (idy.eq.ndycore(1)+ndysol(1)) then
               iymxbcg(idt) = 1
               iymin(idt) = iycum + 1 - nysd
               iymax(idt) = ny
            else
               iymin(idt) = iycum + 1 - nysd
               iymax(idt) = iycum
            endif
ccc            iymin(idt) = iycum + 1 - nysd
ccc            iymax(idt) = iycum
            if (ixmin(idt) .eq. 1 ) ixmnbcg(idt) = 1 	#true boundary
            if (ixmax(idt) .eq. nx) ixmxbcg(idt) = 1	#true boundary
ccc            if (iymin(idt) .eq. 1 ) iymnbcg(idt) = 1	#true boundary
ccc            if (iymax(idt) .eq. ny) iymxbcg(idt) = 1	#true boundary
         enddo
      enddo

c ... Set up arrays for ixpt1g, etc
      do idt = 1, ndomain
         ixpt1g(idt) = ixpt1(1) - (ixmin(idt)-1)
         ixpt2g(idt) = ixpt2(1) - (ixmin(idt)-1)
         iysptrxg(idt) = iysptrx1(1) -(iymin(idt)-1)
      enddo


c ... Set up arrays to identify neighboring domains; the corner cells for
c ... idcorns follow the numbering:  3   4
c ...                                1   2
c...  Do the inner leg region (ix <= ixpt1)
      do idy = 1, ndycore(1)+ndysol(1)
         do idx = 1, ndleg(1,1)
            idt = idx + ndxt*(idy-1)
            idym1g(idt) = max(0, idt-ndxt)
            idyp1g(idt) = min(ndomain+1, idt+ndxt)
            if (idx.eq.1) then
              idxm1g(idt) = 0
            else
              idxm1g(idt) = max(0, idt-1)
            endif
            if (idy.le.ndycore(1) .and. idx.eq.ndleg(1,1)) then
              idxp1g(idt) = min(ndomain, idt+1+ndxcore(1))
            else
              idxp1g(idt) = min(ndomain+1, idt+1)
            endif
            if (idx.eq.1) then
              idcorng(idt,1) = 0
              idcorng(idt,3) = 0
            else
              idcorng(idt,1) = max(0, idxm1g(idt)-ndxt)
              idcorng(idt,3) = min(ndomain+1, idxm1g(idt)+ndxt)
            endif
            idcorng(idt,2) = max(0, idxp1g(idt)-ndxt)
            idcorng(idt,4) = min(ndomain+1, idxp1g(idt)+ndxt)
         enddo
	 if (idy.eq.ndycore(1) .and. ndleg(1,1).gt.0) then #insures 1,4 reciprocity
           idcorng(idt,4) = min(ndomain+1, idt+1+ndxt)     
         endif
      enddo

c...  Do the core region (ixpt1 < ix <= ixpt2)
      do idy = 1, ndycore(1)+ndysol(1)
         do idx = 1, ndxcore(1)
            idt = ndleg(1,1) + idx + ndxt*(idy-1)
            idym1g(idt) = max(0, idt-ndxt)
            idyp1g(idt) = min(ndomain+1, idt+ndxt)
            idxm1g(idt) = max(0, idt-1)
            if (idx.eq.1 .and. ndleg(1,1).eq.0) then
              idxm1g(idt) = 0
            elseif (idx.eq.1 .and. idy.le.ndycore(1))
              idxm1g(idt) = idt-1+ndxcore(1)
            endif
            if (idx.lt.ndxcore(1)) then
              idxp1g(idt) = min(ndomain+1, idt+1)
            elseif (idx.eq.ndxcore(1) .and. idy.le.ndycore(1) .and.
     .                                         ndleg(1,1).gt.0) then
              idxp1g(idt) = idt+1-ndxcore(1)
            elseif (idx.eq.ndxcore(1) .and. idy.gt.ndycore(1)) then
              idxp1g(idt) =  min(ndomain+1, idt+1)
            endif
            if (idx.eq.1 .and. ndleg(1,1).eq.0) then
              idcorng(idt,1) = 0
              idcorng(idt,3) = 0
            else
              idcorng(idt,1) = max(0, idxm1g(idt)-ndxt)
              idcorng(idt,3) = min(ndomain+1, idxm1g(idt)+ndxt)
            endif
            if (idx.eq.ndxcore(1) .and. ndleg(1,2).eq.0) then
              idcorng(idt,2) = ndomain + 1
              idcorng(idt,4) = ndomain + 1
            else
              idcorng(idt,2) = max(0, idxp1g(idt)-ndxt)
              idcorng(idt,4) = min(ndomain+1, idxp1g(idt)+ndxt)
            endif
            if (idx.eq. 1 .and. idy.eq.ndycore(1)) then
              idcorng(idt,3) = min(ndomain+1, idt-1+ndxt)
            endif
         enddo
	 if (idy.eq.ndycore(1)) then  # Insures reciprocity 1 <--> 4 passing
           idcorng(idt,4) = min(ndomain+1, idt+1+ndxt)
         endif
         if (ndleg(1,2).eq.0) idxp1g(idt) = ndomain + 1
      enddo

c...  Do the outer leg region
      do idy = 1, ndycore(1)+ndysol(1)
         do idx = 1, ndleg(1,2)
            idt = ndleg(1,1) + ndxcore(1) + idx + ndxt*(idy-1)
            idym1g(idt) = max(0, idt-ndxt)
            idyp1g(idt) = min(ndomain+1, idt+ndxt)
            if (idx .eq. ndleg(1,2)) then
              idxp1g(idt) = ndomain + 1
            else
              idxp1g(idt) = min(ndomain+1, idt+1)
            endif
            idxm1g(idt) = max(0, idt-1)
            if (idx.eq.1 .and. idy.le.ndycore(1))
     .                        idxm1g(idt) = idt-1-ndxcore(1)
            idcorng(idt,1) = max(0, idxm1g(idt)-ndxt)
            idcorng(idt,3) = min(ndomain+1, idxm1g(idt)+ndxt)
            if (idx .eq. ndleg(1,2)) then
              idcorng(idt,2) = ndomain + 1
              idcorng(idt,4) = ndomain + 1
            else
              idcorng(idt,2) = max(0, idxp1g(idt)-ndxt)
              idcorng(idt,4) = min(ndomain+1, idxp1g(idt)+ndxt)
            endif
            if (idx.eq.1 .and. idy.eq.ndycore(1)) then
              idcorng(idt,3) = min(ndomain+1, idt-1+ndxt)
            endif
         enddo
      enddo

c ... Compute two sets of domains/corner cells touching X-pt that need
c ... special message-passing to preserve serial diff of up staggered 
c ... poloidal mesh at X-pt for PF and core regions
      idxpt(1) = ndxt*(ndycore(1)-1) + 1
      idxpt(2) = idxpt(1) + ndxcore(1)

c ... Compute local domain neq_locg(id), find its maximum, and
c ... allocate space for ivl2gstnl Jac stencil for routine map_var_jac
c ... (must allocate outside of routine map_var_jac)

      do id = 1, ndomain
	neq_locg(id) = 0
        do iy = iymin(id)-iymnbcg(id), iymax(id)+iymxbcg(id)
          do ix = ixmin(id)-ixmnbcg(id), ixmax(id)+ixmxbcg(id)
	    do invar = 1, numvar
              iv = iv + 1
	      neq_locg(id) = neq_locg(id) + 1
	    enddo
          enddo
        enddo
      enddo

c ... Allocate and compute Jacobian stencil for isddcon=1
      call allocjacstnl

      return
      end
c ***** End of subroutine domain_dc *********************************
c ---------------------------------------------------------------------c
      subroutine allocjacstnl

c ... Allocates size for Jacobian stencil 
c ... Compute local domain neq_locg(id), find its maximum, and
c ... allocate space for ivl2gstnl Jac stencil for routine map_var_jac
c ... (must allocate outside of routine map_var_jac)

      implicit none

      Use(Indices_domain_dcg) 	#ndomain,neq_locg,iymin,iymax,iymnbcg,
                                #iymxbcg,neq_locgmx,isddcon
      Use(Math_problem_size)    #neqmx,numvar
      Use(Indices_loc_glob_map) #ivl2gstnl(allocated)

c ... Local variables
      integer iv,id,ix,iy,invar

      do id = 1, ndomain
	neq_locg(id) = 0
        do iy = iymin(id)-iymnbcg(id), iymax(id)+iymxbcg(id)
          do ix = ixmin(id)-ixmnbcg(id), ixmax(id)+ixmxbcg(id)
	    do invar = 1, numvar
              iv = iv + 1
	      neq_locg(id) = neq_locg(id) + 1
	    enddo
          enddo
        enddo
      enddo

c ... Compute max neq_locg
      neq_locgmx = neq_locg(1)
      do id = 2, ndomain
	if (neq_locg(id) > neq_locgmx) neq_locgmx = neq_locg(id)
      enddo

c ... Allocate correct space for ivl2gstnl; preserve ivloc2mdg,2sdg
      call gchange("Indices_loc_glob_map",0)

c ... Compute array maps from loc->global vars & Jac elem (ivl2gstnl)
      if (isddcon == 1) then
        call map_var_jac
      else
        call map_var_jac1d
      endif

      return
      end
c ***** End of subroutine allocjacstnl *********************************
c -------------------------------------------------------------------
      subroutine map_var_jac

c ... Calculates indices and arrays needed for domain decomposition
c ... Specifically, computes ivloc2sdg (global eqn indices for single
c ... domain ordering), ivloc2mdg (global eqn indices for looping over
c ... each of the multiple domains in order)

      implicit none

      Use(Indices_domain_dcg) 	#ndleg,ndxcore,ndycore,ndysol,ndxt,ndyt,ndomain
                                #ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g,isddcon
                                #ixpt1g,ixpt2g,iysptrxg,idxpt,neq_locg
      Use(Parallv)              #nxg,nyg,neqg
      Use(Lsode)                #neq
      Use(Share)  		#nxleg,nxcore,nycore,nysol,isgrdsym
      Use(Dim)                  #nx,ny
      Use(Xpoint_indices)       #ixpt1,ixpt2,iysptrx
      Use(Indices_loc_glob_map) #ivloc2sdg,ivloc2mdg,ivl2gstnl
      Use(Math_problem_size)    #neqmx,numvar

c ... Local variables
      integer id,iycum,nxsd,ixcum,iv,ivglobal,invar,ix,iy,ivgstart,
     .        ivgstartl,iy0,nxlu,nylu,nx2u,ii,iv1,iv2,iv3,iv4,iv5,iv6,
     .        iv7,iv8,iv9,iel,idonor,nxdu,nydu,nxl2u,nyl2u,nxd2u,nyd2u,
     .        ivsxy,ivcu,ivc

c ... Copy global nx and ny to nxg and nyg, so each processor can know 
      nxg = nx
      nyg = ny
      neqg = neq

c ... Compute indices relating local domain yl/yldot vectors, neq_locg(id) and 
c ... corresponding global indices for single domain ordering (ivloc2sdg) 
c ... and for multiple domain ordering(ivloc2mdg)

      ivgstartl = 0
      ivglobal = 0
      do id = 1, ndomain
        iv = 0
	if (ixmin(id) == 1) then
          ivgstart = ivglobal
        else
	  ivgstart = ivgstartl
        endif
        do iy = iymin(id)-iymnbcg(id), iymax(id)+iymxbcg(id)
          iy0 = iymin(id)-iymnbcg(id)
	  ivglobal = (nx+2)*numvar*(iy-iy0) + ivgstart
          do ix = ixmin(id)-ixmnbcg(id), ixmax(id)+ixmxbcg(id)
	    do invar = 1, numvar
              iv = iv + 1
              ivglobal = ivglobal + 1
              ivloc2sdg(iv,id) = ivglobal
              if (id == 1) then
	        ivloc2mdg(iv,id) = iv
	        ivcum(id) = iv
	      else
	        ivloc2mdg(iv,id) = iv + ivcum(id-1)
	        ivcum(id) = ivcum(id-1) + iv
	      endif
              iellast(iv,id) = 0   # just initialize counter array for below
	    enddo
          enddo
c ...     record last ivglobal value on bottom row of domain
          if (iy==iymin(id)-iymnbcg(id)) ivgstartl = ivglobal
        enddo
      enddo

c ... Set stencil indices (ivl2gstnl) for Jacobian
      do id = 1, ndomain  #First do case for internal donor cells: id=idonor
        idonor = id
        iv = 0
	ivsxy = 1-numvar
        nxlu = ixmax(id) - ixmin(id) + 1    # local nx
        nxl2u = nxlu + ixmnbcg(id) + ixmxbcg(id)
        nxd2u = nxl2u
        nylu = iymax(id) - iymin(id) + 1    # local ny
        do iy = 1-iymnbcg(id), nylu+iymxbcg(id)
          do ix = 1-ixmnbcg(id), nxlu+ixmxbcg(id)
            ivsxy = ivsxy + numvar #1st ivc new xy cell
 	    do invar = 1, numvar #loop over local variables
              iv = iv + 1
              iel = 0
              do ivc = 1, numvar #loop over donor variables
		ivcu = ivsxy + ivc - 1
                iv1 = ivcu - (nxd2u+1)*numvar
                if (iy > 0) then  # else iy=0 cells have no -iy stencil comp
                  if(ix>1-ixmnbcg(id) .and. iv1>0) then  # ix=0 done later 
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,id)
                  endif
                  iv2 = ivcu - nxd2u*numvar
                  if(iv2 > 0) then
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id) = ivloc2mdg(iv2,id)
                  endif
	          iv3 = ivcu - (nxd2u-1)*numvar
                  if(ix < nxlu+ixmxbcg(id) .and.iv3>0) then
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id)=ivloc2mdg(iv3,id)
                  elseif(ix < nxlu+ixmxbcg(id) .and. iv3 > 0) then
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id)=ivloc2mdg(iv3,id)
                  endif
                endif
                iv4 = ivcu - numvar
                if (ix>1-ixmnbcg(id) .and. iv4>0) then
		  iel = iel + 1
                  ivl2gstnl(iv,iel,id)=ivloc2mdg(iv4,id)
                endif
                iv5 = ivcu  #no if-test needed
                  iel = iel + 1 
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv5,id)
                iv6 = ivcu + numvar
		if (ix < nxlu+ixmxbcg(id) .and. iv6 <= neq_locg(idonor)) then 
                  iel = iel + 1
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv6,id)
                endif
                if (iy < nylu+iymxbcg(id)) then  # iy=ny+1 cells no iy+ sten compon
                  iv7 = ivcu + (nxd2u-1)*numvar
                  if(ix>1-ixmnbcg(id) .and. iv7<=neq_locg(idonor)) then
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id) = ivloc2mdg(iv7,id)
                  endif
                  iv8 = ivcu + nxd2u*numvar
                  if(iv8<=neq_locg(idonor)) then
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id) = ivloc2mdg(iv8,id)
                  endif
                  iv9 = ivcu + (nxd2u+1)*numvar
	          if(ix < nxlu+ixmxbcg(id) .and. iv9<=neq_locg(idonor)) then
                    iel = iel + 1
                    ivl2gstnl(iv,iel,id) = ivloc2mdg(iv9,id)
                  endif
                endif
              enddo
              iellast(iv,id) = iel
	    enddo
          enddo
        enddo
      enddo

c ... Now fill in boundary cell components from other domains in the 
c ... order of local eqn indices iy=0, iy=ny+1, ix=0, ix=nx+1

      do id = 1, ndomain
        nylu = iymax(id) - iymin(id) + 1                 #local ny
        nyl2u = nylu + iymnbcg(id) + iymxbcg(id)         #local eff. ny+2
        nxlu = ixmax(id) - ixmin(id) + 1                 #local nx
        nxl2u = nxlu + ixmnbcg(id) + ixmxbcg(id)         #local eff. nx+2

c ..    First do iy=0 eqn contrib to Jacobian from donor iy=ny+1 vars
        idonor = idym1g(id)
        if (idonor > 0) then 
          nydu = iymax(idonor) - iymin(idonor) + 1         #donor ny
          nyd2u = nydu + iymnbcg(idonor) + iymxbcg(idonor) #donor eff. ny+2
          nxdu = ixmax(idonor) - ixmin(idonor) + 1         #donor nx
          nxd2u = nxdu + ixmnbcg(idonor) + ixmxbcg(idonor) #donor eff. nx+2
          iv = 0                                           #local eqn number
          ivsxy = nxd2u*(nydu+iymnbcg(idonor)-1)*numvar-numvar+1 #donor left
          iv2 = iv1 + numvar                                   #donor center
          iv3 = iv2 + numvar                                   #donor right
          do ix = 1-ixmnbcg(id), nxlu+ixmxbcg(id)
            do invar = 1, numvar #loop over local variables
              iv = iv + 1
              iel = iellast(iv,id)
              do ivc = 1, numvar #loop over donor variables
                iv1 = ivsxy + ivc - 1
                iv2 = iv1 + numvar
                iv3 = iv2 + numvar
                if (ix>1-ixmnbcg(id) .and. iv1<=neq_locg(id)) then  #corn later
                  iel = iel+1
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
                endif
                if (iv2 <= neq_locg(id)) then
                  iel = iel+1
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv2,idonor)
                endif
                if (ix<nxlu+ixmxbcg(id) .and. iv3<=neq_locg(id)) then #corn later
                  iel = iel+1
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv3,idonor)
                endif
              enddo
              iellast(iv,id) = iel
            enddo
	    ivsxy = ivsxy + numvar
          enddo
        endif

c ..    Second do iy=ny+1 eqn contrib to Jacobian from donor iy=0 vars
        idonor = idyp1g(id)
        if (idonor <= ndomain) then
          ivsxy = -numvar+1 #start donor left
          iv1 = -numvar     #donor left; not used until incr by numvar below
          iv2 = 0           #donor center
          iv3 = numvar      #donor right
          iv = nxl2u*(nylu+iymnbcg(id)-1)*numvar
          do ix = 1-ixmnbcg(id), nxlu+ixmxbcg(id)
            do invar = 1, numvar #loop over local variables
              iv = iv + 1
              iel = iellast(iv,id)
              do ivc = 1, numvar #loop over donor variables
                iv1 = ivsxy + ivc - 1
                iv2 = iv1 + numvar
                iv3 = iv2 + numvar
                if (ix > 0 .and. iv1 > 0) then  #corner cell done later
                  iel = iel+1
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
                endif
                iel = iel+1
                ivl2gstnl(iv,iel,id) = ivloc2mdg(iv2,idonor)
                if (ix < nxlu+ixmxbcg(id)) then #corner cell done later
                  iel = iel+1
                  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv3,idonor)
                endif
              enddo
              iellast(iv,id) = iel
            enddo
	    ivsxy = ivsxy + numvar
          enddo
        endif

c ..    Third do ix=0 eqn contrib to Jacobian from donor ix=nx+1 vars
        idonor = idxm1g(id)
        if (idonor > 0) then
          nxdu = ixmax(idonor) - ixmin(idonor) + 1         #donor nx
          nxd2u = nxdu + ixmnbcg(idonor) + ixmxbcg(idonor) #donor eff. nx+2
          do iy = 1-iymnbcg(id), nylu+iymxbcg(id)
            ivsxy = (iy+iymnbcg(idonor)-1)*nxd2u*numvar-numvar+1 #start donor below
            iv = (iy+iymnbcg(id)-1)*nxl2u*numvar               #init iv (-1)
            do invar = 1, numvar #loop over local variables
              iv = iv + 1
              iel = iellast(iv,id)
              do ivc = 1, numvar #loop over donor variables
                iv1 = ivsxy + ivc - 1
                iv2 = iv1 + nxd2u*numvar                       #donor ctr
                iv3 = iv2 + nxd2u*numvar                       #donor above
                if(iy>1-iymnbcg(id) .and. iv1<=neq_locg(id)) then #corner later
		  iel = iel + 1
		  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
                endif
        	if (iv2 <= neq_locg(id)) then
                  iel = iel + 1
	          ivl2gstnl(iv,iel,id) = ivloc2mdg(iv2,idonor)
                endif
		if (iy<ny+iymxbcg(id) .and. iv3<=neq_locg(id)) then #corn later
                  iel = iel + 1
		  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv3,idonor)
                endif
              enddo
              iellast(iv,id) = iel
            enddo
          enddo
        endif

c ...   Fourth do ix=nx+1 eqn contrib to Jacobian from donor ix=0 vars
        idonor = idxp1g(id)
        if (idonor <= ndomain) then
          nxdu = ixmax(idonor) - ixmin(idonor) + 1         #donor nx
          nxd2u = nxdu + ixmnbcg(idonor) + ixmxbcg(idonor) #donor eff. nx+2
          do iy = 1-iymnbcg(id), nylu+iymxbcg(id)
            ivsxy = (iy+iymnbcg(idonor)-2)*nxd2u*numvar+1  #start donor below
	    iv = (iy+iymnbcg(id))*nxl2u*numvar - numvar    #init iv (-nv-1)
            do invar = 1, numvar #loop over local variables
              iv = iv + 1
              iel = iellast(iv,id)
              do ivc = 1, numvar #loop over donor variables
                iv1 = ivsxy + ivc - 1
                iv2 = iv1 + nxd2u*numvar                   #donor ctr
                iv3 = iv2 + nxd2u*numvar                   #donor above
                if(iy>1-iymnbcg(id) .and. iv1>0) then  # no corn cell iy-1,ix+1
  	          iel = iel + 1
  		  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor) #use 1000?
		endif
		if (iv2 > 0) then
	          iel = iel + 1
		  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv2,idonor) #use 1000?
                endif
		if (iy<nylu+iymxbcg(id) .and. iv3>0) then  # corn later
                  iel = iel + 1
		  ivl2gstnl(iv,iel,id) = ivloc2mdg(iv3,idonor) #use 1000?
                endif
              enddo
              iellast(iv,id) = iel
            enddo
          enddo
        elseif (idonor == ndomain+1)  then  # signifies a physical bdry;self donor
          do iy = 1-iymnbcg(id), nylu+iymxbcg(id)
	    iv = (iy+iymnbcg(id))*nxl2u*numvar - numvar    #init iv (-nv-1)
            do invar = 1, numvar  #loop over local variables
              iv = iv + 1
              do iel = 1, 9*numvar  # give bdry cell full stencil of prev cell
                ivl2gstnl(iv,iel,id) = ivl2gstnl(iv-numvar,iel,id)
              enddo
            enddo
          enddo
        endif

c ...   Lastly, do the corner cell contributions to Jacobian from donors
	  if (iymnbcg(id) == 0 .and. ixmnbcg(id) == 0) then
	     idonor = idcorng(id,1)
             nydu = iymax(idonor) - iymin(idonor) + 1         #donor ny
             nyd2u = nydu + iymnbcg(idonor) + iymxbcg(idonor) #donor eff. ny+2
             nxdu = ixmax(idonor) - ixmin(idonor) + 1         #donor nx
             nxd2u = nxdu + ixmnbcg(idonor) + ixmxbcg(idonor) #donor eff. nx+2
             iv1 = nyd2u*nxd2u*numvar - numvar
	     iv = 0
             do invar = 1, numvar  #loop over local variables
	       iv = iv + 1
               iel = iellast(iv,id)
               do ivc = 1, numvar  #loop over donor variables
                 iv1 = iv1 + 1
                 iel = iel + 1
                 ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
               enddo
	       iellast(iv,id) = iel
             enddo
          endif
          if (iymnbcg(id) == 0 .and. ixmxbcg(id) == 0) then
	     idonor = idcorng(id,2)
             nydu = iymax(idonor) - iymin(idonor) + 1         #donor ny
             nyd2u = nydu + iymnbcg(idonor) + iymxbcg(idonor) #donor eff. ny+2
             nxdu = ixmax(idonor) - ixmin(idonor) + 1         #donor nx
             nxd2u = nxdu + ixmnbcg(idonor) + ixmxbcg(idonor) #donor eff. nx+2
             iv1 = (nyd2u-1)*nxd2u*numvar
	     iv = (nxl2u-1)*numvar
             do invar = 1, numvar  #loop over local variables
	       iv = iv + 1
               iel = iellast(iv,id)
               do ivc = 1, numvar  #loop over donor variables
                 iv1 = iv1 + 1
                 iel = iel + 1
                 ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
               enddo
	       iellast(iv,id) = iel
             enddo
          endif
          if (iymxbcg(id) == 0 .and. ixmnbcg(id) == 0) then
	     idonor = idcorng(id,3)
             nydu = iymax(idonor) - iymin(idonor) + 1         #donor ny
             nyd2u = nydu + iymnbcg(idonor) + iymxbcg(idonor) #donor eff. ny+2
             iv1 = (nxd2u-1)*numvar
	     iv = (nyl2u-1)*nxl2u*numvar
             do invar = 1, numvar  #loop over local variables
	       iv = iv + 1
               iel = iellast(iv,id)
               do ivc = 1, numvar  #loop over donor variables
                 iv1 = iv1 + 1
                 iel = iel + 1
                 ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
               enddo
	       iellast(iv,id) = iel
             enddo
          endif
          if (iymxbcg(id) == 0 .and. ixmxbcg(id) == 0) then
	     idonor = idcorng(id,4)
             iv1 = 0
	     iv = nyl2u*nxl2u*numvar - numvar
             do invar = 1, numvar  #loop over local variables
	       iv = iv + 1
               iel = iellast(iv,id)
               do ivc = 1, numvar  #loop over donor variables
                 iv1 = iv1 + 1
                 iel = iel + 1
                 ivl2gstnl(iv,iel,id) = ivloc2mdg(iv1,idonor)
               enddo
	       iellast(iv,id) = iel
             enddo
           endif

      enddo   # large do-loop over id

      return
      end
c ***** End of subroutine map_var_jac *********************************
c ---------------------------------------------------------------------c
      subroutine map_var_jac1d

c ... Calculates Jacobian indices and arrays needed for 1D (radial) domain 
c ... decomposition or for no decomposition; similar to map_var_jac for 2D.
c ... Specifically, computes ivloc2sdg (global eqn indices for single
c ... domain ordering), ivloc2mdg (global eqn indices for looping over
c ... each of the multiple domains in order)

      implicit none

      Use(Indices_domain_dcg) 	#ndleg,ndxcore,ndycore,ndysol,ndxt,ndyt,ndomain
                                #ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g,isddcon
                                #ixpt1g,ixpt2g,iysptrxg,idxpt,neq_locg,
                                #neq_locgmx
      Use(Parallv)              #nxg,nyg,neqg
      Use(Lsode)                #neq
      Use(Share)  		#nxleg,nxcore,nycore,nysol,isgrdsym
      Use(Dim)                  #nx,ny
      Use(Xpoint_indices)       #ixpt1,ixpt2,iysptrx
      Use(Indices_loc_glob_map) #ivloc2sdg,ivloc2mdg,ivl2gstnl
      Use(Math_problem_size)    #neqmx,numvar
      Use(Selec)                #ixp1,ixm1
      Use(Indexes)              #ivfirst

c ... Local variables
      integer id,iycum,nxsd,ixcum,iv,ivglobal,invar,ix,iy,ivgstart,
     .        ivgstartl,iy0,nxlu,nylu,nx2u,ii,iv1,iv2,iv3,iv4,iv5,iv6,
     .        iv7,iv8,iv9,iel,idonor,nxdu,nydu,nxl2u,nyl2u,nxd2u,nyd2u,
     .        ivsxy,ivcu,ivc,iystart,iyend,iyg,iyglastid,ivdonor(9),
     .        idon,ielem

c ... Copy global nx and ny to nxg and nyg, so each processor can know 
      nxg = nx
      nyg = ny
      neqg = neq

c ... Compute indices relating local domain yl/yldot vectors, neq_locg(id) and 
c ... corresponding global indices for single domain ordering (ivloc2sdg) 
c ... and for multiple domain ordering(ivloc2mdg)

      ivgstartl = 0
      ivglobal = 0
      do id = 1, ndomain
        iv = 0 
	if (ixmin(id) == 1) then
          ivgstart = ivglobal
        else
	  ivgstart = ivgstartl
        endif
        if(isddcon == 2) then
          iystart = iymin(id)-iymnbcg(id)
          iyend = iymax(id)+iymxbcg(id)
        else
          iystart = 0
          iyend = ny+1
        endif
        do iy = iystart, iyend
          iy0 = iystart
	  ivglobal = (nx+2)*numvar*(iy-iy0) + ivgstart
          do ix = 0, nx+1
	    do invar = 1, numvar
              iv = iv + 1
              ivglobal = ivglobal + 1
              ivloc2sdg(iv,id) = ivglobal
	      ivloc2mdg(iv,id) = ivglobal
	    enddo
          enddo
c ...     record last ivglobal value on bottom row of domain
          if (iy==iystart) ivgstartl = ivglobal
        enddo
      enddo

c ... Set stencil indices (ivl2gstnl) for Jacobian
      iyglastid = 0
      do id = 1, ndomain  #First do all iy & internal ix donor cells
        idonor = id
        iv = -numvar   # not 0 as first ix=1 then starts iv=numvar+1
	ivsxy = 1-numvar
        nxlu = nx              # local nx
        nxl2u = nx+2
        nxd2u = nx+2
        if (isddcon == 2) then
          nylu = iymax(id) - iymin(id) + 1    # local ny
          iystart = 1-iymnbcg(id)
          iyend = nylu+iymxbcg(id)
        else
          nylu = ny
          iystart = 0
          iyend = ny+1
        endif
        do iy = iystart, iyend    # iy is local
          iyg = iy + iyglastid
          iv = iv + 2*numvar
	  do ix = 1, nx           # ix is local/global;ix=0,nx+1 below
            ivsxy = ivsxy + numvar          # first ivc of new xy cell
            if (iyg > 0) then
              ivdonor(1) = ivfirst(ixm1(ix,iyg-1),iyg-1)
              ivdonor(2) = ivfirst(ix,            iyg-1)
              ivdonor(3) = ivfirst(ixp1(ix,iyg-1),iyg-1)
            endif
            ivdonor(4) = ivfirst(ixm1(ix,iyg),  iyg)
            ivdonor(5) = ivfirst(ix,            iyg)
            ivdonor(6) = ivfirst(ixp1(ix,iyg),  iyg)
            if (iyg < ny+1) then
              ivdonor(7) = ivfirst(ixm1(ix,iyg+1),iyg+1)
              ivdonor(8) = ivfirst(ix,            iyg+1)
              ivdonor(9) = ivfirst(ixp1(ix,iyg+1),iyg+1)
            endif
	    do invar = 1, numvar #loop over local variables
              iv = iv + 1
              ielem = 0
              do idon = 1, 9  #loop over donor cells (not variables)
                do ivc = 1, numvar #loop over donor variables
 		  if (iyg > 0 .and. idon <= 3) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  elseif (idon > 3 .and. idon <= 6) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  elseif (iyg < ny+1 .and. idon > 6) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  endif
                enddo   # ivc loop
              enddo     # idon loop
	    enddo       # invar loop
          enddo         # ix loop
        enddo           # iy loop
      iyglastid = iyg
      enddo             # id loop

c ... Now do ix=0 boundary cells
      iyglastid = 0
      do id = 1, ndomain  #First do case for internal donor cells: id=idonor
        idonor = id
        iv = -(nx+1)*numvar
	ivsxy = 1-numvar
        nxlu = nx              # local nx
        nxl2u = nx+2
        nxd2u = nx+2
        if (isddcon == 2) then
          nylu = iymax(id) - iymin(id) + 1    # local ny
          iystart = 1-iymnbcg(id)
          iyend = nylu+iymxbcg(id)
        else
          nylu = ny
          iystart = 0
          iyend = ny+1
        endif
        do iy = iystart, iyend    # iy is local
          iv = iv + (nx+1)*numvar
          iyg = iy + iyglastid
	  ix = 0
            ivsxy = ivsxy + numvar          # first ivc of new xy cell
            if (iyg > 0) then
              ivdonor(2) = ivfirst(ix,            iyg-1)
              ivdonor(3) = ivfirst(ixp1(ix,iyg-1),iyg-1)
            endif
            ivdonor(5) = ivfirst(ix,            iyg)
            ivdonor(6) = ivfirst(ixp1(ix,iyg),  iyg)
            if (iyg < ny+1) then
              ivdonor(8) = ivfirst(ix,            iyg+1)
              ivdonor(9) = ivfirst(ixp1(ix,iyg+1),iyg+1)
            endif
	    do invar = 1, numvar #loop over local ix=0 variables
              iv = iv + 1 
              ielem = 0
              do idon = 1, 9  #loop over donor cells (not variables)
                do ivc = 1, numvar #loop over donor variables
 		  if (iyg > 0 .and. idon > 1 .and. idon <= 3) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  elseif (idon > 4 .and. idon <= 6) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  elseif (iyg < ny+1 .and. idon > 7) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  endif
                enddo
              enddo
	    enddo
        enddo
        iyglastid = iyg
      enddo

c ... Now do ix=nx+1 boundary cells
      iyglastid = 0
      do id = 1, ndomain  
        idonor = id
        iv = 0
	ivsxy = 1-numvar
        nxlu = nx              # local nx
        nxl2u = nx+2
        nxd2u = nx+2
        if (isddcon == 2) then
          nylu = iymax(id) - iymin(id) + 1    # local ny
          iystart = 1-iymnbcg(id)
          iyend = nylu+iymxbcg(id)
        else
          nylu = ny
          iystart = 0
          iyend = ny+1
        endif
        do iy = iystart, iyend    # iy is local
          iv = iv + (nx+1)*numvar
          iyg = iy + iyglastid
	  ix = nx+1
            ivsxy = ivsxy + numvar          # first ivc of new xy cell
            if (iyg > 0) then
              ivdonor(1) = ivfirst(ixm1(ix-1,iyg-1),iyg-1)
              ivdonor(2) = ivfirst(ixm1(ix,  iyg-1),iyg-1)
              ivdonor(3) = ivfirst(ix,              iyg-1)
            endif
            ivdonor(4) = ivfirst(ixm1(ix-1,iyg),  iyg)
            ivdonor(5) = ivfirst(ixm1(ix,  iyg),  iyg)
            ivdonor(6) = ivfirst(ix,              iyg)
            if (iy < nylu+1) then
              ivdonor(7) = ivfirst(ixm1(ix-1,iyg+1),iyg+1)
              ivdonor(8) = ivfirst(ixm1(ix,  iyg+1),iyg+1)
              ivdonor(9) = ivfirst(ix,              iyg+1)
            endif
	    do invar = 1, numvar #loop over local variables
              iv = iv + 1
              ielem = 0
              do idon = 1, 9  #loop over donor cells (not variables)
                do ivc = 1, numvar #loop over donor variables
 		  if (iyg > 0 .and. idon <= 3) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  elseif (idon > 3 .and. idon <= 6) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  elseif (iyg < ny+1 .and. idon > 6) then
                    ielem = ielem + 1
                    ivl2gstnl(iv,ielem,id) = ivdonor(idon) + ivc-1
                  endif
                enddo
              enddo
	    enddo
        enddo
        iyglastid = iyg
      enddo

      return
      end
c ***** End of subroutine map_var_jac1d *********************************
c-----------------------------------------------------------------------

      subroutine jacstnlout

      implicit none
c_mpi      include 'mpif.h'

c ... Output a map of the Jacobian stensil fro  map_var_jac1d

c ... All information is passed through common for convenience in
c     calling this routine from the UEDGE> prompt.

c ... Common blocks:
      Use(Math_problem_size)   # neqmx,numvar
      Use(Lsode)           # neq
      Use(Jacobian_full)   # jacfull
      Use(Indices_domain_dcg)  #neq_locgmx
      Use(Indices_loc_glob_map) #ivloc2sdg,ivloc2mdg,ivl2gstnl

c ... Local variables:
      integer ierr
      integer us,iv,ielem
c_mpi      integer my_pe
      character*24 filename

c ... Allocate full Jacobian for jacmap; warning of size
      call remark("*** CAUTION: allocating large jacfull(neq,neq)***")
      call gallot("Jacobian_full",0)
      write (STDOUT,*) '*** Full Jacobian size is neq**2 = ', neq*neq

c ... Fill the full matrix with 1 wherever an elements might exist
      do iv = 1, neq
        do ielem = 1, 9*numvar
          if (ivl2gstnl(iv,ielem,1) > 1.e-50) then
            jacfull(iv,ivl2gstnl(iv,ielem,1)) = 1.
          endif
        enddo
      enddo

c ... Open a file, and output the map.
      call freeus (us)
      filename = 'Jac_stencil.dat'
c_mpi      if(MY_PE().eq.0) then
c_mpi        us = 59
c_mpi        filename = 'Jacobian_map.dat0'
c_mpi      else
c_mpi        us = 69
c_mpi        filename = 'Jacobian_map.dat1'
c_mpi      endif
      open(unit=us, file=filename, status='unknown')
      call jmap (neq, jacfull, us)

c ... Close file, and report file name.
      close(us)
      write (STDOUT,*) ' Jacobian map in data file:  ', filename

      return
      end
c-----------------------------------------------------------------------

      subroutine PackSend_dc_ind(iv_toti)

*     Computes iv_toti(), num. of elem. in visend to be used by send_dc_ind()

      implicit none

      Use(Indices_domain_dcg) 	# ndomain,ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g
                                #vrsend,visend

c ... output variables
      integer iv_toti(128)  #number of elem. in visend sent to proc(id-1)

c ... local variables
      integer id

c ... Pack indices for domain decomposition to each processor
      if (15 > nvisend) then
	  call xerrab('**ERROR send_dc_ind: iv_totc>nvisend; reset nvisend')
        endif
      do id = 1, 128
	if (id <= ndomain) then
          iv_toti(id) = 15
	else
	  iv_toti(id) = 0
        endif
      enddo

      return
      end
c **** End of subroutine PackSend_dc_ind **************
c ---------------------------------------------------------------------c
      subroutine send_dc_ind(iv_toti,visendl_mype)

*     This subroutine packs and sends the indice data to different processesor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcg) 	# ndomain,ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g
                                #vrsend,visend,neq_locg
      Use(Npes_mpi)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typeind, ierr
c_mpi      data typeind/1/

c ... output variables
      integer iv_toti(128)          #number of elem. in visend sent to proc(id-1)
      integer visendl_mype(nvisend)

c ... local variables
      integer id, iv, ic

c ... send indices for domain decomposition to each processor
      do id = 1, ndomain
        iv = 1
        visend(iv) = ixmax(id) - ixmin(id) + 1   # local nx
        iv = iv + 1
        visend(iv) = iymax(id) - iymin(id) + 1   # local ny

        iv = iv + 1
        visend(iv) = ixmnbcg(id)
        iv = iv + 1
        visend(iv) = ixmxbcg(id)
        iv = iv + 1
        visend(iv) = iymnbcg(id)
        iv = iv + 1
        visend(iv) = iymxbcg(id)
        iv = iv + 1
        visend(iv) = idxp1g(id)
        iv = iv + 1
        visend(iv) = idxm1g(id)
        iv = iv + 1
        visend(iv) = idyp1g(id)
        iv = iv + 1
        visend(iv) = idym1g(id)
        iv = iv + 1
	visend(iv) = neq_locg(id)
        do ic = 1,4
          iv = iv + 1
          visend(iv) = idcorng(id,ic)
        enddo
    
cpetsc*        MPI_SEND() to mype hangs. Copy visend to the temp array visendl_mype.
cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv_toti(id)
cpetsc            visendl_mype(ic) = visend(ic)
cpetsc          enddo
cdb_senddc          write(6,*) "[",mype,"]  send_dc_ind() copy to visendl, iv=",iv,"visendl_mype((2)=",visendl_mype(2)
cdb_senddc          call flush(6)
cpetsc        else
cdb_senddc          write(6,*) " [",mype,"] MPI_SEND() iv=",iv," to [",id-1,"] typeind ",typeind
cdb_senddc          call flush(6)
c_mpi         call MPI_SEND(visend, iv_toti(id), MPI_UE_INT, id-1, typeind,
c_mpi     .                   uedgeComm, ierr)
cpetsc        endif
      enddo

      return
      end
c **** End of subroutine send_dc_ind **************
c ---------------------------------------------------------------------c
      subroutine recv_dc_ind(iv_toti, id,visendl_mype)

*     This subroutine receives the index data on different processor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,visendl
      Use(Npes_mpi)             #npes,mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer ierr, typeind, status(MPI_STATUS_SIZE)
c_mpi      data typeind/1/

c ... input variables
      integer id, iv_toti(128)  #domain number (my_pe+1), number elems. 
                                #in visend sent to this process
c      integer visendl_mype(nvisendl)
C     Above works because we compile with no promotion flags
c     To be certain, do
      integer*4 visendl_mype(nvisendl)

c ... local variables
      integer iv, ic

c ... Receive and unpack indices for domain decomposition to each processor
cpetsc      if (mype .ne. 0) then
c_mpi         call MPI_RECV(visendl, iv_toti(id), MPI_UE_INT, 0, typeind,
c_mpi     .                 uedgeComm, status, ierr)
cpetsc      else
cpetsc        do ic = 1,iv_toti(id)
cpetsc          visendl(ic) = visendl_mype(ic)
cpetsc        end do
cpetsc      end if

      iv = 1
      nx_loc = visendl(iv)
      iv = iv + 1
      ny_loc = visendl(iv)
      iv = iv + 1
      ixmnbcl = visendl(iv)
      iv = iv + 1
      ixmxbcl = visendl(iv)
      iv = iv + 1
      iymnbcl = visendl(iv)
      iv = iv + 1
      iymxbcl = visendl(iv)
      iv = iv + 1
      idxp1 = visendl(iv)
      iv = iv + 1
      idxm1 = visendl(iv)
      iv = iv + 1
      idyp1 = visendl(iv)
      iv = iv + 1
      idym1 = visendl(iv) 
      iv = iv + 1
      neq_locl = visendl(iv)
      do ic = 1,4
        iv = iv + 1
        idcorn(ic) = visendl(iv)
      enddo

      return
      end
c **** End of subroutine recv_dc_ind **************
c ---------------------------------------------------------------------c
      subroutine send_yl_map

*     This subroutine packs and sends the indice data to different processesor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcg) 	# ndomain,ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g
                                #vrsend,visend,ivlocsdg,ivlocmdg,neq_locg
      Use(Indices_loc_glob_map) #ivloc2sdg,ivloc2mdg,ivl2gstnl
      Use(Math_problem_size)    #neqmx,numvar
      Use(Indices_domain_dcl) 	#ivlocsdgl,ivlocmdgl,neq_locl	
      Use(Npes_mpi)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typeylind, ierr
c_mpi      data typeylind/4/

c ... output variables
      integer visendl_mype(nvisend)

c ... local variables
	integer id,iv

c ... send indices for domain decomposition to each processor
    
      do id = 1, ndomain
	 do iv = 1, neq_locg(id)
	    visend(iv) = ivloc2sdg(iv,id)
	    visend(iv+neq_locg(id)) = ivloc2mdg(iv,id)
         enddo

c_mpi        if (id-1 > mype) then
cdb_ylmap         write(6,*) "[",mype,"] MPI_SEND yl_map to ",id-1
cdb_ylmap         call flush(6)
c_mpi          call MPI_SEND(visend, 2*neq_locg(id), MPI_UE_INT, id-1, typeylind,
c_mpi     .                   uedgeComm, ierr)
c_mpi        else   #Set PE0 arrays not easily passed via mpi_send
c_mpi          neq_locl = neq_locg(1)
c_mpi          do iv = 1, neq_locl
c_mpi	         ivloc2sdgl(iv) = ivloc2sdg(iv,1)
c_mpi	         ivloc2mdgl(iv) = ivloc2mdg(iv,1)
c_mpi          enddo
c_mpi        endif

      enddo

      return
      end
c **** End of subroutine send_yl_map **************
c ---------------------------------------------------------------------c
      subroutine recv_yl_map(id)

*     This subroutine receives the index data on different processor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,visendl,
                                #ivlocsdgl,ivlocmdgl,neq_locl	
      Use(Npes_mpi)             #npes,mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer ierr, typeylind, status(MPI_STATUS_SIZE)
c_mpi      data typeylind/4/

c ... input variables
      integer id        

c ... local variables
      integer iv

c ... Receive and unpack indices for domain decomposition to each processor
cpetsc      if (mype .ne. 0) then
c_mpi          call MPI_RECV(visendl,2*neq_locl, MPI_UE_INT, 0, typeylind,
c_mpi     .                 uedgeComm, status, ierr)
cdb_ylmap          write(6,*) "[",mype,"] MPI_RECV yl_map from 0"
cdb_ylmap          call flush(6)
c_mpi         do iv = 1, neq_locl
c_mpi	        ivloc2sdgl(iv) = visendl(iv)
c_mpi	        ivloc2mdgl(iv) = visendl(iv+neq_locl)
c_mpi         enddo
cpetsc      endif

      return
      end
c **** End of subroutine recv_yl_map **************
c ---------------------------------------------------------------------c
      subroutine SendRecv_dc_ind(iv_toti)

*     Processor[0] sends the indice data to all other processesors for
*     domain decomposition
*     Prosessor[rank>0] receives the indice data

      implicit none

      Use(Indices_domain_dcg) 	# ndomain,ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g
                                #vrsend,visend
      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,visendl
      Use(Npes_mpi)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typeind, ierr, status(MPI_STATUS_SIZE)
c_mpi      data typeind/1/

c ... output variables
      integer iv_toti(128)  #number of elem. in visendl sent to proc(id-1)

c ... local variables
      integer id, iv, ic

c ... send indices for domain decomposition to each processor
      if (mype .le. 0) then
        do id = ndomain,1,-1
        iv = 1
        visendl(iv) = ixmax(id) - ixmin(id) + 1   # local nx
        iv = iv + 1
        visendl(iv) = iymax(id) - iymin(id) + 1   # local ny

        iv = iv + 1
        visendl(iv) = ixmnbcg(id)
        iv = iv + 1
        visendl(iv) = ixmxbcg(id)
        iv = iv + 1
        visendl(iv) = iymnbcg(id)
        iv = iv + 1
        visendl(iv) = iymxbcg(id)
        iv = iv + 1
        visendl(iv) = idxp1g(id)
        iv = iv + 1
        visendl(iv) = idxm1g(id)
        iv = iv + 1
        visendl(iv) = idyp1g(id)
        iv = iv + 1
        visendl(iv) = idym1g(id)
        do ic = 1,4
          iv = iv + 1
          visendl(iv) = idcorng(id,ic)
        enddo
    
c_mpi      if (id .ne. 1) then
cdb_senddc          write(6,*) " [",mype,"] MPI_SEND() iv=",iv," to [",id-1,"] typeind ",typeind
cdb_senddc          call flush(6)
c_mpi         call MPI_SEND(visendl, iv_toti(id), MPI_UE_INT, id-1, typeind,
c_mpi     .                   uedgeComm, ierr)
c_mpi      endif
      enddo
      endif

c ... Receive and unpack indices for domain decomposition to each processor
c_mpi      if (mype .ne. 0) then
c_mpi         call MPI_RECV(visendl, iv_toti(mype+1), MPI_UE_INT, 0, typeind,
c_mpi     .                 uedgeComm, status, ierr)
c_mpi      end if

      iv = 1
      nx_loc = visendl(iv)
      iv = iv + 1
      ny_loc = visendl(iv)
      iv = iv + 1
      ixmnbcl = visendl(iv)
      iv = iv + 1
      ixmxbcl = visendl(iv)
      iv = iv + 1
      iymnbcl = visendl(iv)
      iv = iv + 1
      iymxbcl = visendl(iv)
      iv = iv + 1
      idxp1 = visendl(iv)
      iv = iv + 1
      idxm1 = visendl(iv)
      iv = iv + 1
      idyp1 = visendl(iv)
      iv = iv + 1
      idym1 = visendl(iv)
      do ic = 1,4
        iv = iv + 1
        idcorn(ic) = visendl(iv)
      enddo

      return
      end
c **** End of subroutine SendRecv_dc_ind **************

c ---------------------------------------------------------------------c
c ---------------------------------------------------------------------c
      subroutine iSendRecv_dc_ind(iv_toti)

*     Processor[0] sends the indice data to all other processesors for
*     domain decomposition Prosessor[rank>0] receives the indice data
*     This is done with a nonblocking send and receive.

      implicit none

      Use(Indices_domain_dcg) 	# ndomain,ixmin,ixmax,iymin,iymax
                                #ixmnbcg,ixmxbcg,iymnbcg,iymxbcg
                                #idxp1g,idxm1g,idyp1g,idym1g
                                #vrsend,visend
      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,visendl
      Use(Npes_mpi)
      Use(Dim)                  #nx
      Use(Bcond)                #matwallo,matwalli
      Use(Parallv)              #nxg,nyg
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
    
      integer ifake  #forces Forthon scripts to put implicit none above here
CC c_mpi      integer typeind,ierr,status(MPI_STATUS_SIZE),req
c_mpi      integer typeind
c_mpi      integer(int4) :: ierr, status(MPI_STATUS_SIZE), req
c_mpi      data typeind/1/
           
c ... output variables
      integer iv_toti(128)  #number of elem. in visendl sent to proc(id-1)

c ... local variables
      integer id, iv, ic, ixu

c ... Post receive indices for this processor domain decomposition
c_mpi      if (mype .ne. 0) then
c_mpi         call MPI_IRECV(visendl,iv_toti(mype+1),MPI_UE_INT,0,typeind,
c_mpi     .                 uedgeComm,req,ierr)
c_mpi      end if

c ... Proc[0] sends indices for domain decomposition to each processor
      if (mype .eq. 0) then
        do id = ndomain,1,-1
            iv = 1
            visendl(iv) = ixmax(id) - ixmin(id) + 1   # local nx
            iv = iv + 1
            visendl(iv) = iymax(id) - iymin(id) + 1   # local ny

            iv = iv + 1
            visendl(iv) = ixmnbcg(id)
            iv = iv + 1
            visendl(iv) = ixmxbcg(id)
            iv = iv + 1
            visendl(iv) = iymnbcg(id)
            iv = iv + 1
            visendl(iv) = iymxbcg(id)
            iv = iv + 1
            visendl(iv) = idxp1g(id)
            iv = iv + 1
            visendl(iv) = idxm1g(id)
            iv = iv + 1
            visendl(iv) = idyp1g(id)
            iv = iv + 1
            visendl(iv) = idym1g(id)
            iv = iv + 1
            visendl(iv) = neq_locg(id)
            do ic = 1,4
              iv = iv + 1
              visendl(iv) = idcorng(id,ic)
            enddo
            do ixu = 0, nx+1
              iv = iv+1
              vrsendl(iv) = matwallog(ixu)
              iv = iv+1
              vrsendl(iv) = matwallig(ixu)
            enddo

c_mpi      if (id .ne. 1) then
cdb_senddc          write(6,*) " [",mype,"] MPI_ISEND() iv=",iv," to [",id-1,"] typeind ",typeind
cdb_senddc          call flush(6)
c_mpi         call MPI_ISEND(visendl, iv_toti(id), MPI_UE_INT, id-1, typeind,
c_mpi     .                  uedgeComm, req, ierr)
c ...         wait until this msg is sent, otherwise visendl cannot be reused!
c_mpi         call MPI_WAIT(req,status,ierr)
cdb_senddc          write(6,*) " [",mype,"] MPI_WAIT() done for id=",id
cdb_senddc          call flush(6)
c_mpi      endif
        enddo
      endif

c_mpi      if (mype .ne. 0) then
c_mpi        call MPI_WAIT(req, status, ierr)
cdb_senddc          write(6,*) " [",mype,"] MPI_WAITALL() done"
cdb_senddc          call flush(6)
c_mpi      endif
    
c ... Unpack indices for this processor
      iv = 1
      nx_loc = visendl(iv)
      iv = iv + 1
      ny_loc = visendl(iv)
      iv = iv + 1
      ixmnbcl = visendl(iv)
      iv = iv + 1
      ixmxbcl = visendl(iv)
      iv = iv + 1
      iymnbcl = visendl(iv)
      iv = iv + 1
      iymxbcl = visendl(iv)
      iv = iv + 1
      idxp1 = visendl(iv)
      iv = iv + 1
      idxm1 = visendl(iv)
      iv = iv + 1
      idyp1 = visendl(iv)
      iv = iv + 1
      idym1 = visendl(iv)
      iv = iv + 1
      neq_locl = visendl(iv)
      do ic = 1,4
        iv = iv + 1
        idcorn(ic) = visendl(iv)
      enddo
      do ixu = 0, nx+1
        iv = iv+1
        matwallo(ixu) = vrsendl(iv)
        iv = iv+1
        matwalli(ixu) = vrsendl(iv)
      enddo

      return
      end
c **** End of subroutine iSendRecv_dc_ind *************

c ---------------------------------------------------------------------c
c ---------------------------------------------------------------------c
      subroutine PackSendglobal(iv_totv, iv_totrz, vrsendlv_mype, vrsendlz_mype)

*     Pack iv_totv, iv_totrz, number of the global data to different processesor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcg)   #ndleg,ndxcore,ndycore,ndysol,ndxt,ndyt,ndomain
                                #ixmin,ixmax,iymin,iymax
                                #vrsend,visend
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Interp)               #nxold,nyold
      Use(RZ_grid_global)       #rmg,zmg,psig,brg,bzg,bpolg,bphig,bg
      Use(Global_vars)          #nisg,upsg,tesg,tisg,ngsg,phisg,afracsg
      Use(Selec)                #ixm1,ixp1
      Use(Xpoint_indices)       #iysptrx
      Use(Npes_mpi)
      Use(Comgeo_g)             #lcong,lconeg,lconig
      Use(Bcond)                #matwallog,matwallig
      Use(Parallv)              #nxg,nyg
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr
c_mpi      data typev/2/, typerz/3/

c ... output variables    #number of elems. in vrsend for vars & rz mesh
      integer iv_totv(128), iv_totrz(128)
      real    vrsendlv_mype(nvrsend),vrsendlz_mype(nvrsend)

c ... local variables
      integer id, ix, iy, iv, ifld, igsp, ic, ixu

c ... Initialize arrays
      do iv = 1, 128
	iv_totv(iv) = 0
	iv_totrz(iv) = 0
      enddo
      do iv = 1, nvrsend
	vrsendlv_mype(iv) = 0.0e0
	vrsendlz_mype(iv) = 0.0e0
      enddo

c ... Pack and send global plasma variables to different processor
      do id = 1, ndomain
        iv = 0
        do iy = iymin(id)-1, iymax(id)+1
          do ix = ixmin(id)-1, ixmax(id)+1
            if (ix .eq. ixmin(id)-1)  then
               ixu = ixm1(ixmin(id),iy)
            elseif(ix .eq. ixmax(id)+1) then
               ixu = ixp1(ixmax(id),iy)
            else
               ixu = ix
            endif
            do ifld = 1, nisp
              iv = iv+1
              vrsend(iv) = nisg(ixu,iy,ifld)
            enddo

            do ifld = 1, nusp
              iv = iv+1
              vrsend(iv) = upsg(ixu,iy,ifld)
            enddo

            iv = iv+1
            vrsend(iv) = tesg(ixu,iy)
            iv = iv+1
            vrsend(iv) = tisg(ixu,iy)

            do igsp = 1, ngsp
              iv = iv+1
              vrsend(iv) = ngsg(ixu,iy,igsp)
            enddo

            iv = iv+1
            vrsend(iv) = phisg(ixu,iy)
            iv = iv+1
            vrsend(iv) = afracsg(ixu,iy)
          enddo
        enddo
        iv_totv(id) = iv
	if (iv > nvrsend) then
	  call xerrab('**ERROR sendglobal: iv_totc>nvrsend; reset nvrsend')
        endif

cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv
cpetsc            vrsendlv_mype(ic) = vrsend(ic)
cpetsc          enddo
cpetsc        endif
      enddo

c ... Pack and send geometry and magnetic values
      do id = 1, ndomain
        iv = 0
        do iy = iymin(id)-1, iymax(id)+1
          do ix = ixmin(id)-1, ixmax(id)+1
            if (ix .eq. ixmin(id)-1)  then
               ixu = ixm1(ixmin(id),iy)
            elseif(ix .eq. ixmax(id)+1) then
               ixu = ixp1(ixmax(id),iy)
            else
               ixu = ix
            endif
            do ic = 0, 4
              iv = iv+1
              vrsend(iv) = rmg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = zmg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = psig(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = brg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bzg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bpolg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bphig(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bg(ixu,iy,ic)
            enddo
            iv = iv+1
            vrsend(iv) = lcong(ixu,iy)
            iv = iv+1
            vrsend(iv) = lconig(ixu,iy)
            iv = iv+1
            vrsend(iv) = lconig(ixu,iy)
          enddo
        enddo
        iv_totrz(id) = iv
	if (iv > nvrsend) then
	  call xerrab('**ERROR sendglobal: iv_totcz>nvrsend; reset nvrsend')
        endif

cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv
cpetsc            vrsendlz_mype(ic) = vrsend(ic)
cpetsc          enddo
cpetsc        endif
      enddo

      return
      end
c **** End of subroutine PackSendglobal **************
c ---------------------------------------------------------------------c
      subroutine sendglobal(iv_totv, iv_totrz, vrsendlv_mype, vrsendlz_mype)

*     This subroutine sends the global data to different processesor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcg) 	#ndleg,ndxcore,ndycore,ndysol,ndxt,ndyt,ndomain
                                #ixmin,ixmax,iymin,iymax
                                #vrsend,visend
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Interp)               #nxold,nyold
      Use(RZ_grid_global)       #rmg,zmg,psig,brg,bzg,bpolg,bphig,bg
      Use(Global_vars)          #nisg,upsg,tesg,tisg,ngsg,phisg,afracsg
      Use(Selec)                #ixm1,ixp1
      Use(Xpoint_indices)       #iysptrx
      Use(Npes_mpi)
      Use(Comgeo_g)             #lcong,lconeg,lconig
      Use(Bcond)                #matwallog,matwallig
      Use(Parallv)              #nxg,nyg
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr
c_mpi      data typev/2/, typerz/3/

c ... output variables    #number of elems. in vrsend for vars & rz mesh
      integer iv_totv(128), iv_totrz(128)
      real    vrsendlv_mype(nvrsend),vrsendlz_mype(nvrsend)

c ... local variables
      integer id, ix, iy, iv, ifld, igsp, ic, ixu

c ... Pack and send global plasma variables to different processor
      do id = 1, ndomain
        iv = 0
        do iy = iymin(id)-1, iymax(id)+1
          do ix = ixmin(id)-1, ixmax(id)+1
            if (ix .eq. ixmin(id)-1)  then
               ixu = ixm1(ixmin(id),iy)
            elseif(ix .eq. ixmax(id)+1) then
               ixu = ixp1(ixmax(id),iy)
            else
               ixu = ix
            endif

            do ifld = 1, nisp
              iv = iv+1
              vrsend(iv) = nisg(ixu,iy,ifld)
            enddo

            do ifld = 1, nusp
              iv = iv+1
              vrsend(iv) = upsg(ixu,iy,ifld)
            enddo

            iv = iv+1
            vrsend(iv) = tesg(ixu,iy)
            iv = iv+1
            vrsend(iv) = tisg(ixu,iy)

            do igsp = 1, ngsp
              iv = iv+1
              vrsend(iv) = ngsg(ixu,iy,igsp)
            enddo

            iv = iv+1
            vrsend(iv) = phisg(ixu,iy)
            iv = iv+1
            vrsend(iv) = afracsg(ixu,iy)
          enddo
        enddo
*        iv_totv(id) = iv
	if (iv > nvrsend) then
	  call xerrab('**ERROR sendglobal: iv_totc>nvrsend; reset nvrsend')
        endif

cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv
cpetsc            vrsendlv_mype(ic) = vrsend(ic)
cpetsc          enddo
cpetsc        else
cdb_sendglobal          write(6,*) " [",mype,"] sendglobal MPI_SEND() to [",id-1,"] iv ",iv," typev ",typev
cdb_sendglobal          call flush(6)
c_mpi           call MPI_SEND(vrsend, iv, MPI_DOUBLE_PRECISION, id-1, typev,
c_mpi     .                   uedgeComm, ierr)
cpetsc        endif
      enddo

c ... Pack and send geometry and magnetic values
      do id = 1, ndomain
        iv = 0
        do iy = iymin(id)-1, iymax(id)+1
          do ix = ixmin(id)-1, ixmax(id)+1
            if (ix .eq. ixmin(id)-1)  then
               ixu = ixm1(ixmin(id),iy)
            elseif(ix .eq. ixmax(id)+1) then
               ixu = ixp1(ixmax(id),iy)
            else
               ixu = ix
            endif
            do ic = 0, 4
              iv = iv+1
              vrsend(iv) = rmg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = zmg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = psig(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = brg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bzg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bpolg(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bphig(ixu,iy,ic)
              iv = iv+1
              vrsend(iv) = bg(ixu,iy,ic)
            enddo
            iv = iv+1
            vrsend(iv) = lcong(ixu,iy)
            iv = iv+1
            vrsend(iv) = lconig(ixu,iy)
            iv = iv+1
            vrsend(iv) = lconig(ixu,iy)
          enddo
        enddo
*        iv_totrz(id) = iv
	if (iv > nvrsend) then
	  call xerrab('**ERROR sendglobal: iv_totcz>nvrsend; reset nvrsend')
        endif

cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv
cpetsc            vrsendlz_mype(ic) = vrsend(ic)
cpetsc          enddo
cpetsc        else
cdb_sendglobal          write(6,*) " [",mype,"] sendglobal MPI_SEND() to [",id-1,"] iv ",iv," typerz ",typerz
cdb_sendglobal          call flush(6)
c_mpi           call MPI_SEND(vrsend, iv, MPI_DOUBLE_PRECISION, id-1, typerz,
c_mpi     .                   uedgeComm, ierr)
cpetsc        endif
      enddo

      return
      end
c **** End of subroutine sendglobal **************
c ---------------------------------------------------------------------c
      subroutine recvglobal(iv_totv, iv_totrz, id, vrsendlv_mype, vrsendlz_mype)

*     This subroutine receives the global data on different processor for
*     domaini decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,vrsendl
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Interp)               #nis,ups,tes,tis,ngs,phis,afracs,nxold,nyold
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Npes_mpi)             #npes,mype
      Use(Comgeo)               #lcon,lcone,lconi
      Use(Bcond)                #matwallo,matwalli
      Use(Parallv)              #nxg,nyg
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr, status(MPI_STATUS_SIZE)
c_mpi      data typev/2/, typerz/3/

c ... input variables
      integer id                          #domain number (my_pe+1)
      integer iv_totv(128), iv_totrz(128) # elems. in vrsend (for mpi)
      real    vrsendlv_mype(nvrsendl),vrsendlz_mype(nvrsendl)

c ... local variables
      integer ix, iy, iv, ifld, igsp, ic

c ... Receive and unpack plasma variables
cpetsc      if (mype .ne. 0) then
c_mpi         call MPI_RECV(vrsendl, iv_totv(id), MPI_DOUBLE_PRECISION, 0,
c_mpi     .                typev, uedgeComm, status, ierr)
cpetsc      else
cpetsc        do ic = 1,iv_totv(id)
cpetsc          vrsendl(ic) = vrsendlv_mype(ic)
cpetsc        end do
cpetsc      end if

      iv = 0
      do iy = 0, ny_loc+1
        do ix = 0, nx_loc+1

           do ifld = 1, nisp
            iv = iv+1
            nis(ix,iy,ifld) = vrsendl(iv)
          enddo

          do ifld = 1, nusp
            iv = iv+1
            ups(ix,iy,ifld) = vrsendl(iv)
          enddo

          iv = iv+1
          tes(ix,iy) = vrsendl(iv)
          iv = iv+1
          tis(ix,iy) = vrsendl(iv)

          do igsp = 1, ngsp
            iv = iv+1
            ngs(ix,iy,igsp) = vrsendl(iv)
           enddo

          iv = iv+1
          phis(ix,iy) = vrsendl(iv)
          iv = iv + 1
          afracs(ix,iy) = vrsendl(iv)
        enddo
      enddo

c ... Receive and unpack geometry and magnetic values
cpetsc      if (mype .ne. 0) then
c_mpi         call MPI_RECV(vrsendl, iv_totrz(id), MPI_DOUBLE_PRECISION, 0,
c_mpi     .                typerz, uedgeComm, status, ierr)
cpetsc      else
cpetsc        do ic = 1,iv_totrz(id)
cpetsc          vrsendl(ic) = vrsendlz_mype(ic)
cpetsc        end do
cpetsc      end if
      iv = 0
      do iy = 0, ny_loc+1
        do ix = 0, nx_loc+1
          do ic = 0, 4
            iv = iv+1
            rm(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            zm(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            psi(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            br(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            bz(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            bpol(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            bphi(ix,iy,ic) = vrsendl(iv)
            iv = iv+1
            b(ix,iy,ic) = vrsendl(iv)
          enddo
          iv = iv+1
          lcon(ix,iy) = vrsendl(iv)
          iv = iv+1
          lcone(ix,iy) = vrsendl(iv)
          iv = iv+1
          lconi(ix,iy) = vrsendl(iv)
        enddo
      enddo

      return
      end
c **** End of subroutine recvglobal **************
c ---------------------------------------------------------------------c
      subroutine sendglobal_xpt(iv_totv, iv_totrz, 
     .                                    vrsendlv_mype, vrsendlz_mype) 

*     This subroutine sends the global data to different processesor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcg) 	#ndleg,ndxcore,ndycore,ndysol,ndxt,ndyt,ndomain
                                #ixmin,ixmax,iymin,iymax
                                #vrsend,visend
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Interp)               #nxold,nyold
      Use(RZ_grid_global)       #rmg,zmg,psig,brg,bzg,bpolg,bphig,bg
      Use(Global_vars)          #nisg,upsg,tesg,tisg,ngsg,phisg,afracsg
      Use(Selec)                #ixm1,ixp1
      Use(Xpoint_indices)       #iysptrx
      Use(Npes_mpi)
      Use(Comgeo_g)             #lcong,lconeg,lconig
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr
c_mpi      data typev/2/, typerz/3/

c ... output variables    #number of elems. in vrsend for vars & rz mesh
      integer iv_totv(128), iv_totrz(128)
      real    vrsendlv_mype(nvrsend),vrsendlz_mype(nvrsend)

c ... local variables
      integer id, idu, ix, iy, iv, ifld, igsp, ic, ixu

c ... Pack and send global plasma variables to different processor
      do idu = 1, 2
        id = idxpt(idu)     #gives X-pt domain being sent to
        iv = 0
        iy = iysptrx + 1   #global ix,iy sent to id up-right corner
        if (idu == 1) then
          ix = ixpt2(1) + 1
        else
          ix = ixpt1(1) + 1
        endif

        do ifld = 1, nisp
          iv = iv+1
          vrsend(iv) = nisg(ix,iy,ifld)
        enddo
        do ifld = 1, nusp
          iv = iv+1
          vrsend(iv) = upsg(ix,iy,ifld)
        enddo
        iv = iv+1
        vrsend(iv) = tesg(ix,iy)
        iv = iv+1
        vrsend(iv) = tisg(ix,iy)
        do igsp = 1, ngsp
          iv = iv+1
          vrsend(iv) = ngsg(ix,iy,igsp)
        enddo

        iv = iv+1
        vrsend(iv) = phisg(ix,iy)
        iv = iv+1
        vrsend(iv) = afracsg(ix,iy)
*        iv_totv(id) = iv

	if (iv > nvrsend) then
	  call xerrab('**ERROR sendglobal: iv_totc>nvrsend; reset nvrsend')
        endif

cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv
cpetsc            vrsendlv_mype(ic) = vrsend(ic)
cpetsc          enddo
cpetsc        else
c_mpi           call MPI_SEND(vrsend, iv, MPI_DOUBLE_PRECISION, id-1, typev,
c_mpi     .                   uedgeComm, ierr)
cdb_sendglobal_xpt           write(6,*) " [",mype,"] sendglobal_xpt MPI_SEND() to [",id-1,"] typev ",typev
cdb_sendglobal_xpt           call flush(6)
cpetsc      endif
      enddo

c ... Pack and send geometry and magnetic values
      do idu = 1, 2
        id = idxpt(idu)     #gives X-pt domain being sent to
        iv = 0
        iy = iysptrx + 1   #global ix,iy sent to id up-right corner
        if (idu == 1) then
          ix = ixpt2(1) + 1
        else
          ix = ixpt1(1) + 1
        endif
        do ic = 0, 4
          iv = iv+1
          vrsend(iv) = rmg(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = zmg(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = psig(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = brg(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = bzg(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = bpolg(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = bphig(ix,iy,ic)
          iv = iv+1
          vrsend(iv) = bg(ix,iy,ic)
        enddo
        iv = iv+1
        vrsend(iv) = lcong(ix,iy)
        iv = iv+1
        vrsend(iv) = lconig(ix,iy)
        iv = iv+1
        vrsend(iv) = lconig(ix,iy)
c        iv_totrz(id) = iv

        if (iv > nvrsend) then
          call xerrab('**ERROR sendglobal:iv_totcz>nvrsend;reset nvrsend')
        endif

cpetsc        if (id-1 .eq. mype) then
cpetsc          do ic = 1,iv
cpetsc            vrsendlz_mype(ic) = vrsend(ic)
cpetsc          enddo
cpetsc        else
c_mpi           call MPI_SEND(vrsend, iv, MPI_DOUBLE_PRECISION, id-1, typerz,
c_mpi     .                   uedgeComm, ierr)
cdb_sendglobal_xpt           write(6,*) " [",mype,"] sendglobal_xpt MPI_SEND() to [",id-1,"] typerz ",typerz
cdb_sendglobal_xpt           call flush(6)
cpetsc       endif

      enddo

      return
      end
c **** End of subroutine sendglobal_xpt **************
c ---------------------------------------------------------------------c
      subroutine recvglobal_xpt(iv_totv, iv_totrz, id,
     .                           vrsendlv_mype, vrsendlz_mype)

*     This subroutine receives the global data on different processor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,vrsendl
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Interp)               #nis,ups,tes,tis,ngs,phis,afracs,nxold,nyold
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Npes_mpi)             #npes,mype
      Use(Comgeo)               #lcon,lcone,lconi
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr, status(MPI_STATUS_SIZE)
c_mpi      data typev/2/, typerz/3/

c ... input variables
      integer id          #receiving domain number
      integer iv_totv(128), iv_totrz(128) # elems. in vrsend (for mpi)
      real    vrsendlv_mype(nvrsendl),vrsendlz_mype(nvrsendl)

c ... local variables
      integer ix, iy, iv, ifld, igsp, ic

c ... Receive and unpack plasma variables
c_mpi      if (mype .ne. 0) then
c_mpi         call MPI_RECV(vrsendl, iv_totv(id), MPI_DOUBLE_PRECISION, 0,
c_mpi     .                typev, uedgeComm, status, ierr)
cdb_sendglobal_xpt         write(6,*) " [",mype,"] recvglobal_xpt MPI_RECV() from [",0,"] typev ",typev
cdb_sendglobal_xpt         call flush(6)
c_mpi      else    # PE0 cannot send to itself; vrsendlv_mpye holds data
c_mpi        do ic = 1,iv_totv(id)
c_mpi          vrsendl(ic) = vrsendlv_mype(ic)
c_mpi        enddo
c_mpi      endif

      iv = 0
c ... Fill upper-right corner ghost cell only
      iy = ny_loc+1
      ix = nx_loc+1

        do ifld = 1, nisp
          iv = iv+1
          nis(ix,iy,ifld) = vrsendl(iv)
        enddo

        do ifld = 1, nusp
          iv = iv+1
          ups(ix,iy,ifld) = vrsendl(iv)
        enddo

        iv = iv+1
        tes(ix,iy) = vrsendl(iv)
        iv = iv+1
        tis(ix,iy) = vrsendl(iv)

        do igsp = 1, ngsp
          iv = iv+1
          ngs(ix,iy,igsp) = vrsendl(iv)
        enddo

        iv = iv+1
        phis(ix,iy) = vrsendl(iv)
        iv = iv + 1
        afracs(ix,iy) = vrsendl(iv)

c_mpi      if (mype .ne. 0) then
c_mpi         call MPI_RECV(vrsendl, iv_totrz(id), MPI_DOUBLE_PRECISION, 0,
c_mpi     .                typerz, uedgeComm, status, ierr)
cdb_sendglobal_xpt         write(6,*) " [",mype,"] recvglobal_xpt MPI_RECV() from [",0,"] typerz ",typerz
cdb_sendglobal_xpt         call flush(6)
c_mpi      else    # PE0 cannot send to itself; vrsendlz_mpye holds data
c_mpi        do ic = 1,iv_totrz(id)
c_mpi          vrsendl(ic) = vrsendlz_mype(ic)
c_mpi        enddo
c_mpi      endif

      iv = 0
c ... Fill upper-right corner ghost cell only
      iy = ny_loc+1
      ix = nx_loc+1
      do ic = 0, 4
        iv = iv+1
        rm(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        zm(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        psi(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        br(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        bz(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        bpol(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        bphi(ix,iy,ic) = vrsendl(iv)
        iv = iv+1
        b(ix,iy,ic) = vrsendl(iv)
      enddo
      iv = iv+1
      lcon(ix,iy) = vrsendl(iv)
      iv = iv+1
      lcone(ix,iy) = vrsendl(iv)
      iv = iv+1
      lconi(ix,iy) = vrsendl(iv)

c ... Received and unpack geometry and magnetic values
c_mpi*      write(6,*) " [",mype,"] recvglobal_xpt MPI_RECV() from [",0,"] iv ",iv," typerz ",typerz

      return
      end
c **** End of subroutine recvglobal_xpt **************
c ---------------------------------------------------------------------c
      subroutine sendbdry_xpt(iv_totvxpt,pedestin)

*     This subroutine sends local corner data for two X-pt processesors
*     domain decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,vrsendl,
                                #typebdy,typecn,iv_totbdy
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Imprad)               #afrac
      Use(Npes_mpi)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr
c_mpi      data typev/2/, typerz/3/

c ... input variables
      integer pedestin      #destination processor (domain-1)

c ... output variables    #number of elems. in vrsend for vars & rz mesh
      integer iv_totv(128)
      integer iv_totvxpt  #number plasma variables sent in mpi message

c ... local variables
      integer ix, iy, iv, ifld, igsp

c ... Pack and send global plasma variables to processor ids
        iv = 0
        iy = 1             #local ix,iy sent to id up-right corner
        ix = 1

        do ifld = 1, nisp
          iv = iv+1
          vrsendl(iv) = ni(ix,iy,ifld)
        enddo
        do ifld = 1, nusp
          iv = iv+1
          vrsendl(iv) = up(ix,iy,ifld)
        enddo
        iv = iv+1
        vrsendl(iv) = te(ix,iy)
        iv = iv+1
        vrsendl(iv) = ti(ix,iy)
        do igsp = 1, ngsp
          iv = iv+1
          vrsendl(iv) = ng(ix,iy,igsp)
        enddo

        iv = iv+1
        vrsendl(iv) = phi(ix,iy)
        iv = iv+1
        vrsendl(iv) = afrac(ix,iy)
*        iv_totv(ids) = iv
        iv_totvxpt = iv

	if (iv > nvrsendl) then
	  call xerrab('**ERROR sendglobal_xpt: iv_totc>nvrsendl; reset nvrsendl')
        endif

ccc        write(6,*) " [",mype,"] sendbdry_xpt MPI_SEND() to [",pedestin,"] iv ",iv," typev ",typev
ccc        call flush(6)
c_mpi           call MPI_SEND(vrsendl, iv, MPI_DOUBLE_PRECISION, pedestin,
c_mpi     .                    typev, uedgeComm, ierr)

      return
      end
c **** End of subroutine sendbdry_xpt **************
c ---------------------------------------------------------------------c
      subroutine recvbdry_xpt(iv_totvxpt, pedonor)

*     This subroutine receives the global data on different processor for
*     domaini decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,vrsendl
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Imprad)               #afrac
      Use(Npes_mpi)             #npes,mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev, typerz, ierr, status(MPI_STATUS_SIZE)
c_mpi      data typev/2/, typerz/3/

c ... input variables
      integer pedonor      #donor processor (domain-1)
      integer iv_totv(128), iv_totrz(128) # elems. in vrsend (for mpi)
      integer iv_totvxpt  #number plasma variables sent in mpi message

c ... local variables
      integer ix, iy, iv, ifld, igsp

c ... Receive and unpack plasma variables
c_mpi         call MPI_RECV(vrsendl, iv_totvxpt, MPI_DOUBLE_PRECISION,
c_mpi     .               pedonor, typev, uedgeComm, status, ierr)

      iv = 0
c ... Fill upper-right corner ghost cell only
      iy = ny_loc+1
      ix = nx_loc+1

        do ifld = 1, nisp
          iv = iv+1
          ni(ix,iy,ifld) = vrsendl(iv)
        enddo

        do ifld = 1, nusp
          iv = iv+1
          up(ix,iy,ifld) = vrsendl(iv)
        enddo

        iv = iv+1
        te(ix,iy) = vrsendl(iv)
        iv = iv+1
        ti(ix,iy) = vrsendl(iv)

        do igsp = 1, ngsp
          iv = iv+1
          ng(ix,iy,igsp) = vrsendl(iv)
        enddo

        iv = iv+1
        phi(ix,iy) = vrsendl(iv)
        iv = iv + 1
        afrac(ix,iy) = vrsendl(iv)

      return
      end
c **** End of subroutine recvbdry_xpt **************
c ---------------------------------------------------------------------c
      subroutine sendbdry(id,vrsendl_mype)

*     This subroutine sends boundary data from one processor to another

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,vrsendl,
                                #typebdy,typecn,iv_totbdy
      Use(Indices_domain_dcg) 	#idym1g,idyp1g,idxm1g,idxp1g,idcorng
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Npes_mpi)             #mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer ierr

c ... input variables
      integer id                #domain number (my_pe+1)
C$$$      real    vrsendl_mype(nvrsendl) #rm later!

      REAL, ALLOCATABLE :: vrsendl_mype(:)

c ... local variables
      integer ix, iy, iv, ifld, igsp, ii

c ... Pack and send plasma boundary variables along x-direction
      do iy = 1, ny, ny-1
        if (iy.eq.1  .and. iymnbcl.eq.1) goto 10   #bndry is exter bndry, skip
        if (iy.eq.ny .and. iymxbcl.eq.1) goto 10   #bndry is exter bndry, skip
        iv = 0
        do ix = 1-ixmnbcl, nx+ixmxbcl
          do ifld = 1, nisp
          iv = iv+1
          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
          vrsendl(iv) = te(ix,iy)
          iv = iv+1
          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
          vrsendl(iv) = phi(ix,iy)
        enddo

c_mpi        if (iy.eq.1) then
cdb_sendbdry           write(6,*) " [",mype,"] send to [",idym1-1,"] in sendbdry() "
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idym1-1) then
cdb_sendbdry             write(6,*) " [",mype,"] **** ERROR: "
cdb_sendbdry             call flush(6)
c_mpi             do ii = 1,iv
c_mpi               vrsendl_mype = vrsendl(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_SEND(vrsendl, iv, MPI_DOUBLE_PRECISION,
c_mpi     .              idym1-1, typebdy(1), uedgeComm, ierr)
c_mpi           endif
c_mpi        else
cdb_sendbdry           write(6,*) " [",mype,"] send to [",idyp1-1,"] in sendbdry() "
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idyp1-1) then
cdb_sendbdry             write(6,*) " [",mype,"] **** ERROR: "
cdb_sendbdry             call flush(6)
c_mpi             do ii = 1,iv
c_mpi               vrsendl_mype = vrsendl(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_SEND(vrsendl, iv, MPI_DOUBLE_PRECISION,
c_mpi     .              idyp1-1, typebdy(2), uedgeComm, ierr)
c_mpi           endif
c_mpi        endif

  10  continue
      enddo


c ... Pack and send plasma boundary variables along y-direction
      do ix = 1, nx, nx-1
        if (ix.eq.1  .and. ixmnbcl.eq.1) goto 20   #bndry is exter bndry, skip
        if (ix.eq.nx .and. ixmxbcl.eq.1) goto 20   #bndry is exter bndry, skip
        iv = 0
        do iy = 1-iymnbcl, ny+iymxbcl
          do ifld = 1, nisp
          iv = iv+1
          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
          vrsendl(iv) = te(ix,iy)
          iv = iv+1
          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
          vrsendl(iv) = phi(ix,iy)
        enddo

c_mpi        if (ix.eq.1) then
cdb_sendbdry           write(6,*) " [",mype,"] send to idxm1-1 [",idxm1-1,"] in sendbdry() "
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idxm1-1) then
cdb_sendbdry             write(6,*) " [",mype,"] **** ERROR: copy to visendl_mype "
cdb_sendbdry             call flush(6)
c_mpi             do ii = 1,iv
c_mpi               vrsendl_mype = vrsendl(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_SEND(vrsendl, iv, MPI_DOUBLE_PRECISION,
c_mpi     .             idxm1-1, typebdy(3), uedgeComm, ierr)
c_mpi           endif
c_mpi        else
cdb_sendbdry           write(6,*) " [",mype,"] send to idxp1 [",idxp1-1,"] in sendbdry() "
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idxp1-1) then
cdb_sendbdry             write(6,*) " [",mype,"] **** ERROR: copy to visendl_mype"
cdb_sendbdry             call flush(6)
c_mpi             do ii = 1,iv
c_mpi               vrsendl_mype = vrsendl(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_SEND(vrsendl, iv, MPI_DOUBLE_PRECISION,
c_mpi     .             idxp1-1, typebdy(4), uedgeComm, ierr)
c_mpi           endif
c_mpi        endif

  20  continue
      enddo

c ... Pack and send corner data (ix=nx, iy=ny)
      do ix = 1, nx, nx-1
        if (ix.eq.1  .and. ixmnbcl.eq.1) goto 40   #bndry is exter bndry, skip
        if (ix.eq.nx .and. ixmxbcl.eq.1) goto 40   #bndry is exter bndry, skip
        do iy = 1, ny, ny-1
          iv = 0
          if (iy.eq.1  .and. iymnbcl.eq.1) goto 30 #bndry is exter bndry, skip
          if (iy.eq.ny .and. iymxbcl.eq.1) goto 30 #bndry is exter bndry, skip
          do ifld = 1, nisp
          iv = iv+1
          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
          vrsendl(iv) = te(ix,iy)
          iv = iv+1
          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
          vrsendl(iv) = phi(ix,iy)

          if (ix.eq.1) then
            if(iy.eq.1) then
              ii = 1
            else
              ii = 3
            endif
          else
            if (iy.eq.1) then
              ii = 2
            else
              ii = 4
            endif
          endif
cdb_sendbdry           write(6,*) " [",mype,"] send to  idcorn(ii)-1 [",idcorn(ii)-1,"] in sendbdry() "
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idcorn(ii)-1) then
             write(6,*) " [",mype,"] **** ERROR: "
             call flush(6)
c_mpi             do ii = 1,iv
c_mpi               vrsendl_mype = vrsendl(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_SEND(vrsendl, iv, MPI_DOUBLE_PRECISION,
c_mpi     .         idcorn(ii)-1, typecn(ii), uedgeComm, ierr)
c_mpi           endif
  30    continue
        enddo
  40  continue
      enddo

cdb_sendbdry      write(6,*)"[",mype,"] exit sendbdry"
cdb_sendbdry      call flush(6)

      return
      end
c **** End of subroutine sendbdry **************
c-----------------------------------------------------------------------
      subroutine recvbdry(id,vrsendl_mype)

*     This subroutine receives the global data on different processor for
*     domain decomposition

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,vrsendl
                                #typebdy,typecn,iv_totbdy
      Use(Indices_domain_dcg) 	#idym1g,idyp1g,idxm1g,idxp1g,idcorng
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Npes_mpi)             #mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer status(MPI_STATUS_SIZE)
c_mpi      integer(int4) ierr

c ... input variables
      integer id
C$$$      real    vrsendl_mype(nvrsendl) #rm later!

      REAL, ALLOCATABLE :: vrsendl_mype(:)

c ... local variables
      integer ix, iy, iv, ifld, igsp, iis, iir,ii

c ... Receive and unpack plasma boundary variables along x-direction
      do iy = 0, ny+1, ny+1
        if (iy.eq.0  .and. iymnbcl.eq.1) goto 10   #bndry is exter bndry, skip
        if (iy.eq.ny+1 .and. iymxbcl.eq.1) goto 10 #bndry is exter bndry, skip

c_mpi        if (iy.eq.0) then
cdb_sendbdry           write(6,*) " [",mype,"] recv from idym1-1 [",idym1-1,"] in recvbdry()"
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idym1-1) then
c_mpi             do ii = 1,iv_totbdy(1)
c_mpi               vrsendl(ii) = vrsendl_mype(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_RECV(vrsendl,iv_totbdy(1),MPI_DOUBLE_PRECISION,
c_mpi     .           idym1-1, typebdy(2), uedgeComm, status, ierr)
c_mpi           endif
c_mpi        else
cdb_sendbdry           write(6,*) " [",mype,"] recv from idyp1-1 [",idyp1-1,"] in recvbdry()"
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idyp1-1) then
c_mpi             do ii = 1,iv_totbdy(2)
c_mpi               vrsendl(ii) = vrsendl_mype(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_RECV(vrsendl,iv_totbdy(2),MPI_DOUBLE_PRECISION,
c_mpi     .          idyp1-1, typebdy(1), uedgeComm, status, ierr)
c_mpi           endif
c_mpi        endif

        iv = 0
        do ix = 1-ixmnbcl, nx+ixmxbcl

          do ifld = 1, nisp
          iv = iv+1
          ni(ix,iy,ifld) = vrsendl(iv)
          enddo

          do ifld = 1, nusp
            iv = iv+1
            up(ix,iy,ifld) = vrsendl(iv)
          enddo

          iv = iv+1
          te(ix,iy) = vrsendl(iv)
          iv = iv+1
          ti(ix,iy) = vrsendl(iv)

          do igsp = 1, ngsp
            iv = iv+1
            ng(ix,iy,igsp) = vrsendl(iv)
          enddo

          iv = iv+1
          phi(ix,iy) = vrsendl(iv)
        enddo
  10  continue
      enddo


c ... Receive and unpack plasma boundary variables along y-direction
      do ix = 0, nx+1, nx+1
        if (ix.eq.0  .and. ixmnbcl.eq.1) goto 20   #bndry is exter bndry, skip
        if (ix.eq.nx+1 .and. ixmxbcl.eq.1) goto 20 #bndry is exter bndry, skip

c_mpi        if (ix.eq.0) then
cdb_sendbdry           write(6,*) " [",mype,"] recv from idxm1-1 [",idxm1-1,"] in recvbdry()"
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idxm1-1) then
c_mpi             do ii = 1,iv_totbdy(3)
c_mpi               vrsendl(ii) = vrsendl_mype(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_RECV(vrsendl,iv_totbdy(3),MPI_DOUBLE_PRECISION,
c_mpi     .           idxm1-1, typebdy(4), uedgeComm, status, ierr)
c_mpi           endif
c_mpi        else
cdb_sendbdry           write(6,*) " [",mype,"] recv from idxp1-1 [",idxp1-1,"] in recvbdry() of size [", iv_totbdy(4), "]"
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idxp1-1) then
c_mpi             do ii = 1,iv_totbdy(4)
c_mpi               vrsendl(ii) = vrsendl_mype(ii)
c_mpi             enddo
c_mpi           else
c_mpi             call MPI_RECV(vrsendl,iv_totbdy(4),MPI_DOUBLE_PRECISION,
c_mpi     .          idxp1-1, typebdy(3), uedgeComm, status, ierr)
c_mpi           endif
c_mpi        endif

        iv = 0
        do iy = 1-iymnbcl, ny+iymxbcl
          do ifld = 1, nisp
          iv = iv+1
          ni(ix,iy,ifld) = vrsendl(iv)
          enddo

          do ifld = 1, nusp
            iv = iv+1
            up(ix,iy,ifld) = vrsendl(iv)
          enddo

          iv = iv+1
          te(ix,iy) = vrsendl(iv)
          iv = iv+1
          ti(ix,iy) = vrsendl(iv)

          do igsp = 1, ngsp
            iv = iv+1
            ng(ix,iy,igsp) = vrsendl(iv)
          enddo

          iv = iv+1
          phi(ix,iy) = vrsendl(iv)
        enddo
  20  continue
      enddo

c ... Receive and unpack corner data (ix=nx, iy=ny)
      do ix = 0, nx+1, nx+1
        if (ix.eq.0  .and. ixmnbcl.eq.1) goto 40   #bndry is exter bndry, skip
        if (ix.eq.nx+1 .and. ixmxbcl.eq.1) goto 40 #bndry is exter bndry, skip
        do iy = 0, ny+1, ny+1
          iv = 0
          if (iy.eq.0  .and. iymnbcl.eq.1) goto 30  #bndry is exter bndry, skip
          if (iy.eq.ny+1 .and. iymxbcl.eq.1) goto 30#bndry is exter bndry, skip

          if (ix.eq.0) then
            if(iy.eq.0) then
              iir = 1          # note for passing  4 <--> 1
              iis = 4
ccc              iir = 4          # note for passing  4 <--> 1
ccc              iis = 1
            else
              iir = 3          # note for passing  2 <--> 3
              iis = 2
ccc              iir = 2          # note for passing  2 <--> 3
ccc              iis = 3
            endif
          else
            if (iy.eq.0) then
              iir = 2          # note for passing  3 <--> 2
              iis = 3
ccc              iir = 3          # note for passing  3 <--> 2
ccc              iis = 2
            else
              iir = 4          # note for passing  1 <--> 4
              iis = 1
ccc              iir = 1          # note for passing  1 <--> 4
ccc              iis = 4
            endif
          endif

cdb_sendbdry          write(6,*) " [",mype,"] recv from idcorn(iir)-1 [",idcorn(iir)-1,"] in recvbdry() of size [",iv_totbdy(4+iir),"]"
cdb_sendbdry          call flush(6)
c_mpi          if (mype .eq. idcorn(iir)-1) then
c_mpi            do ii = 1,iv_totbdy(4+iir)
c_mpi              vrsendl(ii) = vrsendl_mype(ii)
c_mpi            enddo
c_mpi          else
c_mpi            call MPI_RECV(vrsendl,iv_totbdy(4+iir),
c_mpi     .                  MPI_DOUBLE_PRECISION, idcorn(iir)-1,
c_mpi     .                  typecn(iis), uedgeComm, status, ierr)
c_mpi          endif
          do ifld = 1, nisp
          iv = iv+1
          ni(ix,iy,ifld) = vrsendl(iv)
          enddo

          do ifld = 1, nusp
            iv = iv+1
            up(ix,iy,ifld) = vrsendl(iv)
          enddo

          iv = iv+1
          te(ix,iy) = vrsendl(iv)
          iv = iv+1
          ti(ix,iy) = vrsendl(iv)

          do igsp = 1, ngsp
            iv = iv+1
            ng(ix,iy,igsp) = vrsendl(iv)
          enddo

          iv = iv+1
          phi(ix,iy) = vrsendl(iv)

  30    continue
        enddo
  40  continue
      enddo
c_mpi      call MPI_BARRIER(uedgeComm, ierr)
cdb_sendbdry      write(6,*)"[",mype,"] exit recvbdry"
cdb_sendbdry      call flush(6)

      return
      end
c **** End of subroutine recvbdry **************
c ---------------------------------------------------------------------c
      subroutine sendbdry_asz(id,myiv)

*     This subroutine sends the array size, iv_totbdy(pe,8) for the vector
*     vrsendl used by sendbdry; a nonzero value indicates a message will
*     be sent

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,
                                #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl,vrsendl
                                #typebdyi,typecni,iv_totbdy
      Use(Indices_domain_dcg) 	#idym1g,idyp1g,idxm1g,idxp1g,idcorng
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Npes_mpi)             #mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer ierr,req

c ... input variables
      integer id                #domain number (my_pe+1)
      integer myiv              #iv previously sent to mype

c ... local variables
      integer ix, iy, iv, ifld, igsp, ii

c ... Pack and send plasma boundary variables along x-direction
c_mpi      if (ny .eq. 1) then
c_mpi        write(6,*) '**** ERROR: only one iy-cell for PE =', mype
c_mpi        call flush(6)
c_mpi      endif
c_mpi      myiv = -1

      do iy = 1, ny, ny-1
        if (iy.eq.1  .and. iymnbcl.eq.1) goto 10   #bndry is exter bndry, skip
        if (iy.eq.ny .and. iymxbcl.eq.1) goto 10   #bndry is exter bndry, skip
        iv = 0
        do ix = 1-ixmnbcl, nx+ixmxbcl
          do ifld = 1, nisp
          iv = iv+1
ccc          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
ccc           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = te(ix,iy)
          iv = iv+1
ccc          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
ccc            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = phi(ix,iy)
        enddo

        write(6,*) " [",mype,"] idym1 idyp1:",idym1,idyp1," idxm1, idxp1:",idxm1,idxp1
        call flush(6)

c_mpi        if (iy.eq.1) then
           write(6,*) " [",mype,"] send to [",idym1-1,"] typebdyi(1) ",typebdyi(1)
           call flush(6)
c_mpi           if (mype .eq. idym1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself myiv",iv," - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             call MPI_SEND(iv, 1, MPI_UE_INT,
c_mpi     .              idym1-1, typebdyi(1), uedgeComm, ierr)
c_mpi           endif
c_mpi        else
cdb_sendbdry           write(6,*) " [",mype,"] send to [",idyp1-1,"] typebdyi(2) ",typebdyi(2)
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idyp1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             call MPI_SEND(iv, 1, MPI_UE_INT,
c_mpi     .              idyp1-1, typebdyi(2), uedgeComm, ierr)
c_mpi           endif
c_mpi        endif

  10  continue
      enddo


c ... Pack and send plasma boundary variables along y-direction
c_mpi      if (nx .eq. 1) then
c_mpi        write(6,*) '**** ERROR: only one ix-cell for PE =', mype
c_mpi        call flush(6)
c_mpi      endif
      do ix = 1, nx, nx-1
        if (ix.eq.1  .and. ixmnbcl.eq.1) goto 20   #bndry is exter bndry, skip
        if (ix.eq.nx .and. ixmxbcl.eq.1) goto 20   #bndry is exter bndry, skip
        iv = 0
        do iy = 1-iymnbcl, ny+iymxbcl
          do ifld = 1, nisp
          iv = iv+1
ccc          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
ccc           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = te(ix,iy)
          iv = iv+1
ccc          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
ccc            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = phi(ix,iy)
        enddo

c_mpi        if (ix.eq.1) then
cdb_sendbdry           write(6,*) " [",mype,"] send to [",idxm1-1,"] typebdyi(3) ",typebdyi(3)
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idxm1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself myiv",iv,"- ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             call MPI_SEND(iv, 1, MPI_UE_INT,
c_mpi     .             idxm1-1, typebdyi(3), uedgeComm, ierr)
c_mpi           endif
c_mpi        else
cdb_sendbdry           write(6,*) " [",mype,"] send to [",idxp1-1,"] typebdyi(4) ",typebdyi(4)
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idxp1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             call MPI_SEND(iv, 1, MPI_UE_INT,
c_mpi     .             idxp1-1, typebdyi(4), uedgeComm, ierr)
c_mpi           endif
c_mpi        endif

  20  continue
      enddo

c ... Pack and send corner data (ix=nx, iy=ny)
      do ix = 1, nx, nx-1
        if (ix.eq.1  .and. ixmnbcl.eq.1) goto 40   #bndry is exter bndry, skip
        if (ix.eq.nx .and. ixmxbcl.eq.1) goto 40   #bndry is exter bndry, skip
        do iy = 1, ny, ny-1
          iv = 0
          if (iy.eq.1  .and. iymnbcl.eq.1) goto 30 #bndry is exter bndry, skip
          if (iy.eq.ny .and. iymxbcl.eq.1) goto 30 #bndry is exter bndry, skip
          do ifld = 1, nisp
          iv = iv+1
ccc          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
ccc           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = te(ix,iy)
          iv = iv+1
ccc          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
ccc            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = phi(ix,iy)

          if (ix.eq.1) then
            if(iy.eq.1) then
              ii = 1
            else
              ii = 3
            endif
          else
            if (iy.eq.1) then
              ii = 2
            else
              ii = 4
            endif
          endif
cdb_sendbdry           write(6,*) " [",mype,"] send to [",idcorn(ii)-1,"] corner typecni(ii)",typecni(ii)
cdb_sendbdry           call flush(6)
c_mpi           if (mype .eq. idcorn(ii)-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi            call MPI_SEND(iv, 1, MPI_UE_INT, idcorn(ii)-1,
c_mpi     .                  typecni(ii), uedgeComm, ierr)
c_mpi           endif

  30    continue
        enddo
  40  continue
      enddo

      return
      end
c **** End of subroutine sendbdry_asz **************
c-----------------------------------------------------------------------
      subroutine recvbdry_asz(id,myiv)

*     This subroutine receives the array iv_totbdy on the proper processor
*     so that recvbdry knows how many data values will be coming on each

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,vrsendl
                                #typebdyi,typecni,iv_totbdy
      Use(Indices_domain_dcg)       #idym1g,idyp1g,idxm1g,idxp1g,idcorng
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Npes_mpi)             #mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer ierr, status(MPI_STATUS_SIZE),req

c ... input varibles
      integer id
      integer myiv

c ... local variables
      integer ix, iy, iv, ifld, igsp, iis, iir

c ... Receive and unpack array sizes along x-direction (fixed iy = 0 or ny+1)
      do iy = 0, ny+1, ny+1
        if (iy.eq.0  .and. iymnbcl.eq.1) goto 10   #bndry is exter bndry, skip
        if (iy.eq.ny+1 .and. iymxbcl.eq.1) goto 10 #bndry is exter bndry, skip

c_mpi        if (iy.eq.0) then
c_mpi           if (mype .eq. idym1-1) then
c_mpi             iv_totbdy(1) = myiv
        WRITE (*,*) "iv_totbdy(1) set to [", iv_totbdy(1), "]"
c_mpi           else
c_mpi             call MPI_RECV(iv_totbdy(1), 1, MPI_UE_INT, idym1-1,
c_mpi     .               typebdyi(2), uedgeComm, status, ierr)
c_mpi           endif
           write(6,*) " [",mype,"] recv from [",idym1-1,"] of value [", iv_totbdy(1), "] into iv_totbdy(1)"
           call flush(6)
c_mpi        else
c_mpi           if (mype .eq. idyp1-1) then
c_mpi             iv_totbdy(2) = myiv
           WRITE (*,*) "iv_totbdy(2) set to [", iv_totbdy(2), "]"
c_mpi           else
c_mpi             call MPI_RECV(iv_totbdy(2), 1, MPI_UE_INT, idyp1-1,
c_mpi     .               typebdyi(1), uedgeComm, status, ierr)
c_mpi           endif
           write(6,*) " [",mype,"] recv from [",idyp1-1,"] of value [", iv_totbdy(2), "] into iv_totbdy(2)"
           call flush(6)
c_mpi        endif

  10     continue
      enddo


c ... Receive and unpack array sizes along y-direction (fixed ix = 0 or nx+1)
      do ix = 0, nx+1, nx+1
        if (ix.eq.0  .and. ixmnbcl.eq.1) goto 20   #bndry is exter bndry, skip
        if (ix.eq.nx+1 .and. ixmxbcl.eq.1) goto 20 #bndry is exter bndry, skip

c_mpi        if (ix.eq.0) then
c_mpi           if (mype .eq. idxm1-1) then
c_mpi             iv_totbdy(3) = myiv
        WRITE (*,*) "iv_totbdy(3) set to [", iv_totbdy(3), "]"
c_mpi           else
c_mpi             call MPI_RECV(iv_totbdy(3), 1, MPI_UE_INT, idxm1-1,
c_mpi     .                typebdyi(4), uedgeComm, status, ierr)
c_mpi           endif
           write(6,*) " [",mype,"] recv from [",idxm1-1,"] of value [", iv_totbdy(3), "] into iv_totbdy(3)"
           call flush(6)
c_mpi        else
c_mpi           if (mype .eq. idxp1-1) then
c_mpi             iv_totbdy(4) = myiv
           WRITE (*,*) "iv_totbdy(4) set to [", iv_totbdy(4), "]"
c_mpi           else
c_mpi             call MPI_RECV(iv_totbdy(4), 1, MPI_UE_INT, idxp1-1,
c_mpi     .                typebdyi(3), uedgeComm, status, ierr)
c_mpi           endif
           write(6,*) " [",mype,"] recv from [",idxp1-1,"] of value [", iv_totbdy(4), "] into iv_totbdy(4)"
           call flush(6)
c_mpi        endif

  20     continue
      enddo

c ... Receive and unpack corner data (ix=nx, iy=ny)
      do ix = 0, nx+1, nx+1
        if (ix.eq.0  .and. ixmnbcl.eq.1) goto 40   #bndry is exter bndry, skip
        if (ix.eq.nx+1 .and. ixmxbcl.eq.1) goto 40 #bndry is exter bndry, skip
        do iy = 0, ny+1, ny+1
          iv = 0
          if (iy.eq.0  .and. iymnbcl.eq.1) goto 30  #bndry is exter bndry, skip
          if (iy.eq.ny+1 .and. iymxbcl.eq.1) goto 30#bndry is exter bndry, skip

          if (ix.eq.0) then
            if(iy.eq.0) then
              iir = 1          # note for passing  4 <--> 1
              iis = 4
ccc              iir = 4          # note for passing  4 <--> 1
ccc              iis = 1
            else
              iir = 3          # note for passing  2 <--> 3
              iis = 2
ccc              iir = 2          # note for passing  2 <--> 3
ccc              iis = 3
            endif
          else
            if (iy.eq.0) then
              iir = 2          # note for passing  3 <--> 2
              iis = 3
ccc              iir = 3          # note for passing  3 <--> 2
ccc              iis = 2
            else
              iir = 4          # note for passing  1 <--> 4
              iis = 1
ccc              iir = 1          # note for passing  1 <--> 4
ccc              iis = 4
            endif
          endif

c_mpi          call MPI_RECV(iv_totbdy(4+iir), 1, MPI_UE_INT, idcorn(iir)-1,
c_mpi     .                  typecni(iis), uedgeComm, status, ierr)
           write(6,*) " [",mype,"] recv from [",idcorn(iir)-1,"] of value [", iv_totbdy(4+iir), "] into iv_totbdy(",4+iir,")"
           call flush(6)

  30    continue
        enddo
  40    continue
      enddo

      return
      end
c **** End of subroutine recvbdry_asz **************
c-----------------------------------------------------------------------
      subroutine isendrecv_bdry_asz(id,myiv,a_sz)

*     This subroutine sends the array size, iv_totbdy(pe,8) for the vector
*     vrsendl used by sendbdry; a nonzero value indicates a message will
*     be sent.  Then it receives the array iv_totbdy on the proper processor
*     so that recvbdry knows how many data values will be coming on each

      implicit none

      Use(Indices_domain_dcl) 	#nx_loc,ny_loc,idxp1,idxm1,idyp1,idcorn,vrsendl
                                #typebdyi,typecni,iv_totbdy
      Use(Indices_domain_dcg) 	#idym1g,idyp1g,idxm1g,idxp1g,idcorng
      Use(Dim)                  #nx,ny,nxm,nisp,nusp,ngsp
      Use(RZ_grid_info)         #rm,zm,psi,br,bz,bpol,bphi,b
      Use(Compla)               #ni,up,te,ti,ng,phi
      Use(Npes_mpi)             #mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
CC c_mpi      integer ierr,status(MPI_STATUS_SIZE),req(8),reqSEND,to_recv(8)
c_mpi      integer(int4) :: ierr, status(MPI_STATUS_SIZE),req(8),
c_mpi     .  reqSEND, to_recv(8), foo

c ... input varibles
      integer id                #domain number (my_pe+1)
      integer myiv              #iv previously sent to mype
      integer :: a_sz(8), b_sz(8)           #space for the nonblocking sends and recvs
ccc      integer(mpiue) :: a_sz(8), b_sz(8)           #space for the nonblocking sends and recvs

c ... local variables
      integer ix, iy, iv, ifld, igsp, ii, iis, iir

c_mpi      do ii=1,8
c_mpi        to_recv(ii)=0
c_mpi      end do

c ... Pack and send plasma boundary variables along x-direction
c_mpi      if (ny .eq. 1) then
c_mpi        write(6,*) '**** ERROR: only one iy-cell for PE =', mype
c_mpi        call flush(6)
c_mpi      endif
c_mpi      myiv = -1

      do iy = 1, ny, ny-1
        if (iy.eq.1  .and. iymnbcl.eq.1) goto 11   #bndry is exter bndry, skip
        if (iy.eq.ny .and. iymxbcl.eq.1) goto 11   #bndry is exter bndry, skip
        iv = 0
        do ix = 1-ixmnbcl, nx+ixmxbcl
          do ifld = 1, nisp
          iv = iv+1
ccc          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
ccc           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = te(ix,iy)
          iv = iv+1
ccc          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
ccc            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = phi(ix,iy)
        enddo

c_mpi        if (iy.eq.1) then
c_mpi           if (mype .eq. idym1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself myiv",iv," - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             a_sz(1)=iv
cdb_sendbdry             write(6,*) "[",mype,"] -> [",idym1-1,"] | [",typebdyi(1),"] : [",a_sz(1),"] @ [ 1 ]"
cdb_sendbdry             call flush(6)
c_mpi             call MPI_ISEND(a_sz(1),1,MPI_UE_INT,idym1-1,
c_mpi     .               typebdyi(1),uedgeComm,reqSEND,ierr)
c_mpi             iv_totbdy(1)=iv
c_mpi           endif
c_mpi        else
c_mpi           if (mype .eq. idyp1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             a_sz(2)=iv
cdb_sendbdry             write(6,*) "[",mype,"] -> [",idyp1-1,"] | [",typebdyi(2),"] : [",a_sz(2),"] @ [ 1 ]"
cdb_sendbdry             call flush(6)
c_mpi             call MPI_ISEND(a_sz(2),1,MPI_UE_INT,idyp1-1,
c_mpi     .               typebdyi(2),uedgeComm,reqSEND,ierr)
c_mpi             iv_totbdy(2)=iv
c_mpi           endif
c_mpi        endif

  11  continue
      enddo
      
c ... Pack and send plasma boundary variables along y-direction
c_mpi      if (nx .eq. 1) then
c_mpi        write(6,*) '**** ERROR: only one ix-cell for PE =', mype
c_mpi        call flush(6)
c_mpi      endif
      do ix = 1, nx, nx-1
        if (ix.eq.1  .and. ixmnbcl.eq.1) goto 21   #bndry is exter bndry, skip
        if (ix.eq.nx .and. ixmxbcl.eq.1) goto 21   #bndry is exter bndry, skip
        iv = 0
        do iy = 1-iymnbcl, ny+iymxbcl
          do ifld = 1, nisp
          iv = iv+1
ccc          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
ccc           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = te(ix,iy)
          iv = iv+1
ccc          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
ccc            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = phi(ix,iy)
        enddo

c_mpi        if (ix.eq.1) then
c_mpi           if (mype .eq. idxm1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself myiv",iv,"- ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             a_sz(3)=iv
cdb_sendbdry             write(6,*) "[",mype,"] -> [",idxm1-1,"] | [",typebdyi(3),"] : [",a_sz(3),"] @ [ 1 ]"
cdb_sendbdry             call flush(6)
c_mpi             call MPI_ISEND(a_sz(3),1,MPI_UE_INT,idxm1-1,
c_mpi     .               typebdyi(3),uedgeComm,reqSEND,ierr)
c_mpi             iv_totbdy(3)=iv
c_mpi           endif
c_mpi        else
c_mpi           if (mype .eq. idxp1-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi             a_sz(4)=iv
cdb_sendbdry             write(6,*) "[",mype,"] -> [",idxp1-1,"] | [",typebdyi(4),"] : [",a_sz(4),"] @ [ 1 ]"
cdb_sendbdry             call flush(6)
c_mpi             call MPI_ISEND(a_sz(4),1,MPI_UE_INT,idxp1-1,
c_mpi     .               typebdyi(4),uedgeComm,reqSEND,ierr)
c_mpi             iv_totbdy(4)=iv
c_mpi           endif
c_mpi        endif

  21  continue
      enddo

c ... Pack and send corner data (ix=nx, iy=ny)
      do ix = 1, nx, nx-1
        if (ix.eq.1  .and. ixmnbcl.eq.1) goto 41   #bndry is exter bndry, skip
        if (ix.eq.nx .and. ixmxbcl.eq.1) goto 41   #bndry is exter bndry, skip
        do iy = 1, ny, ny-1
          iv = 0
          if (iy.eq.1  .and. iymnbcl.eq.1) goto 31 #bndry is exter bndry, skip
          if (iy.eq.ny .and. iymxbcl.eq.1) goto 31 #bndry is exter bndry, skip
          do ifld = 1, nisp
          iv = iv+1
ccc          vrsendl(iv) = ni(ix,iy,ifld)
          enddo

          do ifld = 1, nusp
            iv = iv+1
ccc           vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = te(ix,iy)
          iv = iv+1
ccc          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
ccc            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
ccc          vrsendl(iv) = phi(ix,iy)

          if (ix.eq.1) then
            if(iy.eq.1) then
              ii = 1
            else
              ii = 3
            endif
          else
            if (iy.eq.1) then
              ii = 2
            else
              ii = 4
            endif
          endif
c_mpi           if (mype .eq. idcorn(ii)-1) then
c_mpi             myiv = iv
cdb_sendbdry             write(6,*) " [",mype,"] send to itself - ERROR in sendbdry_asz()"
cdb_sendbdry             call flush(6)
c_mpi           else
c_mpi            a_sz(4+ii)=iv
cdb_sendbdry             write(6,*) "[",mype,"] -> [",idcorn(ii)-1,"] | [",typecni(ii),"] : [",a_sz(4+ii),"] @ [ 1 ]"
cdb_sendbdry             call flush(6)
c_mpi            call MPI_ISEND(a_sz(4+ii),1,MPI_UE_INT,idcorn(ii)-1,
c_mpi     .              typecni(ii),uedgeComm,reqSEND,ierr)
c_mpi             iv_totbdy(4+ii)=iv
c_mpi           endif

  31    continue
        enddo
  41  continue
      enddo

cc_mpi      call MPI_BARRIER(uedgeComm,ierr)

c ... Receive and unpack array sizes along x-direction (fixed iy = 0 or ny+1)
      do iy = 0, ny+1, ny+1
        if (iy.eq.0  .and. iymnbcl.eq.1) goto 10   #bndry is exter bndry, skip
        if (iy.eq.ny+1 .and. iymxbcl.eq.1) goto 10 #bndry is exter bndry, skip

c_mpi        if (iy.eq.0) then
c_mpi           if (mype .eq. idym1-1) then
c_mpi             iv_totbdy(1) = myiv
        WRITE (*,*) "iv_totbdy(1) is set to [", iv_totbdy(1), "]"
c_mpi           else
c_mpi             call MPI_IRECV(b_sz(1),1,MPI_UE_INT,idym1-1,
c_mpi     .               typebdyi(2),uedgeComm,req(1),ierr)
cdb_mpass        write(6,*) "[",mype,"] <- [",idym1-1,"] | [",typebdyi(2),"] : [*",req(1),"*] @ [ 1 ]"
cdb_mpass        call flush(6)
c_mpi             to_recv(1)=1
c_mpi           endif
c_mpi        else
c_mpi           if (mype .eq. idyp1-1) then
c_mpi             iv_totbdy(2) = myiv
c_mpi           else
c_mpi             call MPI_IRECV(b_sz(2),1,MPI_UE_INT,idyp1-1,
c_mpi     .               typebdyi(1),uedgeComm,req(2),ierr)
cdb_mpass        write(6,*) "[",mype,"] <- [",idyp1-1,"] | [",typebdyi(1),"] : [*",req(2),"*] @ [ 1 ]"
cdb_mpass        call flush(6)
c_mpi             to_recv(2)=1
c_mpi           endif
c_mpi        endif

  10     continue
      enddo
      
c ... Receive and unpack array sizes along y-direction (fixed ix = 0 or nx+1)
      do ix = 0, nx+1, nx+1
        if (ix.eq.0  .and. ixmnbcl.eq.1) goto 20   #bndry is exter bndry, skip
        if (ix.eq.nx+1 .and. ixmxbcl.eq.1) goto 20 #bndry is exter bndry, skip

c_mpi        if (ix.eq.0) then
c_mpi           if (mype .eq. idxm1-1) then
c_mpi             iv_totbdy(3) = myiv
c_mpi        WRITE (*,*) "iv_totbdy(3) is set to [", iv_totbdy(3), "]"
c_mpi           else
c_mpi             call MPI_IRECV(b_sz(3),1,MPI_UE_INT,idxm1-1,
c_mpi     .               typebdyi(4),uedgeComm,req(3),ierr)
cdb_mpass        write(6,*) "[",mype,"] <- [",idxm1-1,"] | [",typebdyi(4),"] : [*",req(3),"*] @ [ 1 ]"
cdb_mpass        call flush(6)
c_mpi             to_recv(3)=1
c_mpi           endif
c_mpi        else
c_mpi           if (mype .eq. idxp1-1) then
c_mpi             iv_totbdy(4) = myiv
c_mpi           else
c_mpi             call MPI_IRECV(b_sz(4),1,MPI_UE_INT,idxp1-1,
c_mpi     .               typebdyi(3),uedgeComm,req(4),ierr)
cdb_mpass        write(6,*) "[",mype,"] <- [",idxp1-1,"] | [",typebdyi(3),"] : [*",req(4),"*] @ [ 1 ]"
cdb_mpass        call flush(6)
c_mpi             to_recv(4)=1
c_mpi           endif
c_mpi        endif

  20     continue
      enddo

c ... Receive and unpack corner data (ix=nx, iy=ny)
      do ix = 0, nx+1, nx+1
        if (ix.eq.0  .and. ixmnbcl.eq.1) goto 40   #bndry is exter bndry, skip
        if (ix.eq.nx+1 .and. ixmxbcl.eq.1) goto 40 #bndry is exter bndry, skip
        do iy = 0, ny+1, ny+1
          iv = 0
          if (iy.eq.0  .and. iymnbcl.eq.1) goto 30  #bndry is exter bndry, skip
          if (iy.eq.ny+1 .and. iymxbcl.eq.1) goto 30#bndry is exter bndry, skip

          if (ix.eq.0) then
            if(iy.eq.0) then
              iir = 1          # note for passing  4 <--> 1
              iis = 4
ccc              iir = 4          # note for passing  4 <--> 1
ccc              iis = 1
            else
              iir = 3          # note for passing  2 <--> 3
              iis = 2
ccc              iir = 2          # note for passing  2 <--> 3
ccc              iis = 3
            endif
          else
            if (iy.eq.0) then
              iir = 2          # note for passing  3 <--> 2
              iis = 3
ccc              iir = 3          # note for passing  3 <--> 2
ccc              iis = 2
            else
              iir = 4          # note for passing  1 <--> 4
              iis = 1
ccc              iir = 1          # note for passing  1 <--> 4
ccc              iis = 4
            endif
          endif

c_mpi         call MPI_IRECV(b_sz(4+iir),1,MPI_UE_INT,idcorn(iir)-1,
c_mpi     .           typecni(iis),uedgeComm,req(4+iir),ierr)
cdb_mpass        write(6,*) "[",mype,"] <- [",idcorn(iir)-1,"] | [",typecni(iis),"] : [*",req(4+iir),"*] @ [ 1 ]"
cdb_mpass        call flush(6)
c_mpi         to_recv(4+iir)=1

  30    continue
        enddo
  40    continue
      enddo

c_mpi       do ii=1,8
c_mpi         if (to_recv(ii).eq.1) then
c_mpi           foo = req(ii)
c_mpi           call MPI_WAIT(req(ii),status,ierr)
c_mpi           iv_totbdy(ii) = b_sz(ii)
cdb_mpass      write(6,*) "[",mype,"] <| [*",foo,"*] : [",iv_totbdy(ii),"]"
c_mpi         end if
c_mpi       end do

      return
      end
c **** End of subroutine isendrecv_bdry_asz **************
c-----------------------------------------------------------------------
      subroutine sendloc_glob(id)

*     This subroutine sends the local data to the master processor (PE0)
*     in order to construct the global solution

      implicit none

      Use(Indices_domain_dcl) 	#ixmnbcl,ixmxbcl,iymnbcl,iymxbcl
      Use(Dim)                  #nisp,nusp,ngsp
      Use(Interp)               #nxoldg,nyoldg,afracs
      Use(Compla)               #ni,up,te,ti,ng,phi,afracs
      Use(Indices_domain_dcg)       #vrsend
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev2, typeiv2, typev2_start, typeiv2_start, ierr
c_mpi      data typev2_start/10/, typeiv2_start/138/  #here 138=128+10; 128 max proc

c ... input/output variables    # domain number id (= mype+1)
      integer id

c ... local variables
      integer ix, iy, iv, ivmx, ifld, igsp, iv_totv(128)

c ... Pack and send global plasma variables to different processor
      iv = 0
      do iy = 1-iymnbcl, ny+iymxbcl
        do ix = 1-ixmnbcl, nx+ixmxbcl

          do ifld = 1, nisp
            iv = iv+1
            vrsendl(iv) = ni(ix,iy,ifld)
          enddo
           do ifld = 1, nusp
            iv = iv+1
            vrsendl(iv) = up(ix,iy,ifld)
          enddo

          iv = iv+1
          vrsendl(iv) = te(ix,iy)
          iv = iv+1
          vrsendl(iv) = ti(ix,iy)

          do igsp = 1, ngsp
            iv = iv+1
            vrsendl(iv) = ng(ix,iy,igsp)
          enddo

          iv = iv+1
          vrsendl(iv) = phi(ix,iy)
          iv = iv+1
          vrsendl(iv) = afracs(ix,iy)
        enddo
      enddo
      ivmx = iv
      if (ivmx > nvrsendl) then
      call xerrab('**ERROR sendloc_glob: iv_totc>nvrsendl; reset nvrsendl')
      endif
c_mpi      if (id > 1) then #if id=1, no send/recv to self; put in vrsend
c_mpi         call MPI_SEND(ivmx, 1, MPI_UE_INT,
c_mpi     .             0, typeiv2_start+id, uedgeComm, ierr)

c_mpi         call MPI_SEND(vrsendl, ivmx, MPI_DOUBLE_PRECISION,
c_mpi     .             0, typev2_start+id, uedgeComm, ierr)
c_mpi      else
c_mpi        do iv = 1, ivmx
c_mpi          vrsend(iv) = vrsendl(iv)
c_mpi        enddo
c_mpi      endif

      return
      end
c **** End of subroutine sendloc_glob **************
c ---------------------------------------------------------------------c
      subroutine recvloc_glob(id)

*     This subroutine receives the local data from different processor
*     to construct the global solution

      implicit none

      Use(Indices_domain_dcg) 	#ixmnbcg,ixmxbcg,iymnbcg,iymxbcg,
                                #ixmin,ixmax,iymin,iymax
                                #vrsend,visend
      Use(Dim)                  #nisp,nusp,ngsp
      Use(Interp)               #nxoldg,nyoldg
      Use(Global_vars)          #nisg,upsg,tesg,tisg,ngsg,phisg,afracsg
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
c_mpi      integer typev2, typeiv2, typev2_start, typeiv2_start, ierr,
c_mpi     . status(MPI_STATUS_SIZE)
c_mpi      data typev2_start/10/, typeiv2_start/138/  #here 138=128+10; 128 max proc

c ... input variables
      integer id                          #domain number (my_pe+1)

c ... local variables
      integer ix, iy, iv, ifld, igsp, iv_totv(128)

c_mpi      if (id > 1) then #if id = 1, no mpi send/recv to self; in vrsend
c ... Receive iv_totv, the length of vrsend
c_mpi         call MPI_RECV(iv_totv(id), 1, MPI_UE_INT, id-1,
c_mpi     .               typeiv2_start+id, uedgeComm, status, ierr)

c ... Receive and unpack plasma variables
c_mpi         call MPI_RECV(vrsend, iv_totv(id), MPI_DOUBLE_PRECISION, id-1,
c_mpi     .                typev2_start+id, uedgeComm, status, ierr)
c_mpi      endif

      iv = 0
      do iy = iymin(id)-iymnbcg(id), iymax(id)+iymxbcg(id)
        do ix = ixmin(id)-ixmnbcg(id), ixmax(id)+ixmxbcg(id)
           do ifld = 1, nisp
            iv = iv+1
            nisg(ix,iy,ifld) = vrsend(iv)
          enddo

          do ifld = 1, nusp
            iv = iv+1
            upsg(ix,iy,ifld) = vrsend(iv)
          enddo

          iv = iv+1
          tesg(ix,iy) = vrsend(iv)
          iv = iv+1
          tisg(ix,iy) = vrsend(iv)

          do igsp = 1, ngsp
            iv = iv+1
            ngsg(ix,iy,igsp) = vrsend(iv)
          enddo

          iv = iv+1
          phisg(ix,iy) = vrsend(iv)
          iv = iv + 1
          afracsg(ix,iy) = vrsend(iv)
        enddo
      enddo

      return
      end
c **** End of subroutine recvloc_glob **************
c-----------------------------------------------------------------------
      subroutine init_pll

c **- This rountine sets up the mesh and passes necessary information
c     to the various processors when UEDGE is run on a parallel machine
c     with MPI

      implicit none

      Use(Dim)
      Use(Xpoint_indices)
      Use(Indices_domain_dcg)
      Use(Share)
      Use(Indices_domain_dcl)
      Use(Grid)
      Use(Interp)
      Use(Cdv)
      Use(Npes_mpi)
      Use(UEint)
ccc      Use(PETSc_lib)
      Use(Comgeo)            #area_core
      Use(Bcond)             #fngysog,fngysig,albedoog,mawallog,...
      Use(Parallv)           #nxg,nyg
      Use(Comgeo_g)          #lcong,lconeg,lconig
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here

      integer a_sz(8),igl,ixu,iyu
c_mpi      integer id,iv
c_mpi      integer ierr, status(MPI_STATUS_SIZE),
c_mpi     .        iv_toti(128), iv_totv(128), iv_totrz(128), nd_tot

c_mpi*  tmp arrays need to be optimized later!
C$$$      integer   visendl_mype(nvisendl),myiv
C$$$      real      vrsendlv_mype(nvrsendl),vrsendlz_mype(nvrsendl)

      INTEGER, ALLOCATABLE :: visendl_mype(:)
      REAL, ALLOCATABLE :: vrsendlv_mype(:)
      REAL, ALLOCATABLE :: vrsendlz_mype(:)

      INTEGER :: myiv

      ALLOCATE(visendl_mype(nvisendl))
      ALLOCATE(vrsendlv_mype(nvrsendl))
      ALLOCATE(vrsendlz_mype(nvrsendl))

      if (mype .le. 0) then
         ifexmain = 1
            call allocate
            call gallot("Comgeo_g",0)
         ifexmain = 0
         ig = 1
         call ueinit
         call s2copy (nx+2,ny+2,lcon,1,nx+2,lcong,1,nx+2)
         call s2copy (nx+2,ny+2,lconi,1,nx+2,lconig,1,nx+2)
         call s2copy (nx+2,ny+2,lcone,1,nx+2,lconeg,1,nx+2)
         nxoldg = nx
         nyoldg = ny

c ...    Section for interpolating to new mesh size
         if (isumesh2.eq.1) then
            nxold = nx
            nyold = ny
            call gchange("Interp",0)
            call gridseq
            nxleg(1,1) = nxleg(2,1)
            nxleg(1,2) = nxleg(2,2)
            nxcore(1,1) = nxcore(2,1)
            nxcore(1,2) = nxcore(2,2)
            nysol(1) = nysol(2)
            nycore(1) = nycore(2)
ccc            nxomit = ??
            restart = 1
            ifexmain = 1
               call allocate
               call gallot("Comgeo_g",0)
            ifexmain = 0
            ig = 1
            call ueinit
            isimesh = 1
            nxoldg = nx
            nyoldg = ny
         endif
c ...    End section for interpolating to new mesh

         call gchange("Global_vars",0)
         call gallot("Indices_domain_dcg",0)
         call globalvars
         call domain_dc
         if (ndomain .ne. npes) then
            write(6,*) '****[',mype,'] ERROR, ndomain',ndomain,'.ne.npes',npes
ctdr            call flush(6)
            stop
         endif

c_mpi         call PackSend_dc_ind(iv_toti)
c_mpi         call PackSendglobal(iv_totv, iv_totrz,vrsendlv_mype,vrsendlz_mype)
      endif  # end of ifloop for mype .le. 0

       call gchange("Indices_domain_dcl",0)

c_mpi       call MPI_BCAST(ndomain, 1, MPI_UE_INT, 0, uedgeComm, ierr)
c_mpi       call MPI_BCAST(nxg, 1, MPI_UE_INT, 0, uedgeComm, ierr)
c_mpi       call MPI_BCAST(nyg, 1, MPI_UE_INT, 0, uedgeComm, ierr)
c_mpi       call MPI_BCAST(neqg, 1, MPI_UE_INT, 0, uedgeComm, ierr)
c_mpi       call MPI_BCAST(idxpt, 2, MPI_UE_INT, 0, uedgeComm, ierr)

c_mpi       call MPI_BCAST(iysptrxg, 128, MPI_UE_INT, 0, uedgeComm, ierr)
c_mpi       call MPI_BCAST(ixpt1g, 128, MPI_UE_INT, 0, uedgeComm, ierr)
c_mpi       call MPI_BCAST(ixpt2g, 128, MPI_UE_INT, 0,uedgeComm, ierr)
c_mpi       call MPI_BCAST(ispwrbc, ndomainmx, MPI_UE_INT, 0,uedgeComm, ierr)

c_mpi       call MPI_BCAST(ixmin, ndomainmx, MPI_UE_INT, 0,uedgeComm, ierr)
c_mpi       call MPI_BCAST(iymin, ndomainmx, MPI_UE_INT, 0,uedgeComm, ierr)
c_mpi       call MPI_BCAST(ixmax, ndomainmx, MPI_UE_INT, 0,uedgeComm, ierr)
c_mpi       call MPI_BCAST(iymax, ndomainmx, MPI_UE_INT, 0,uedgeComm, ierr)

c ... Send global wall source information
c_mpi      call gchange("Bcond",0)
c_mpi      call MPI_BCAST(fngysog, (nxg+2)*ngsp, MPI_DOUBLE_PRECISION, 0, 
c_mpi     .               uedgeComm, ierr)
c_mpi      call MPI_BCAST(fngysig, (nxg+2)*ngsp, MPI_DOUBLE_PRECISION, 0, 
c_mpi     .               uedgeComm, ierr)
c_mpi      call MPI_BCAST(albedoog, (nxg+2)*ngsp, MPI_DOUBLE_PRECISION, 0, 
c_mpi     .               uedgeComm, ierr)
c_mpi      call MPI_BCAST(albedoig, (nxg+2)*ngsp, MPI_DOUBLE_PRECISION, 0, 
c_mpi     .               uedgeComm, ierr)
ccc Problem using integer BCAST for static build only - reverted direct send
ccc_mpi      call MPI_BCAST(matwallog, nxg, MPI_UE_INT, 0, 
ccc_mpi     .               uedgeComm, ierr)
ccc_mpi      call MPI_BCAST(matwallig, nxg, MPI_UE_INT, 0,
ccc_mpi     .               uedgeComm, ierr)

ccc       if (mype .eq. 1) then
ccc          ixpt1(1) = ixpt1(1) - 4
ccc          ixpt2(1) = ixpt2(1) - 4
ccc       endif

ccc       call MPI_BCAST(iv_toti, npes, MPI_UE_INT, 0,
c_mpi       call MPI_BCAST(iv_toti, 128, MPI_UE_INT, 0,uedgeComm, ierr)
cc ... Send global variables and mesh
ccc       call MPI_BCAST(iv_totv, npes, MPI_UE_INT, 0,
c_mpi       call MPI_BCAST(iv_totv,  128, MPI_UE_INT, 0,uedgeComm, ierr)
ccc       call MPI_BCAST(iv_totrz, npes, MPI_UE_INT, 0,
c_mpi       call MPI_BCAST(iv_totrz, 128, MPI_UE_INT, 0,uedgeComm, ierr)

c ...  proc[0] sends indice information;  proc[>0] receive. 
c_mpi       call iSendRecv_dc_ind(iv_toti)

c_mpi      if (mype .le. 0) then
c_mpi         call send_yl_map
c_mpi      endif
c_mpi      call recv_yl_map(mype+1)

c ...  With indice info, can now do allocation for each processor
c         if (mype .eq. 0) then
cc_mpi            nxg = nx
cc_mpi            nyg = ny
c         endif

         nx = nx_loc
         nxm = nx
         ny = ny_loc
	 ixlb(1) = 0    #limited to single null geometry
	 ixrb(1) = nx
         isallloc = 1

c_mpi      if (mype .eq. 0) then
c_mpi         call sendglobal(iv_totv, iv_totrz,vrsendlv_mype,vrsendlz_mype)
c_mpi      endif
         call allocate

c ... Copy portion of global wall sources to local arrays; all guard cell
c ... values are left empty
      if (isddcon > 0) then  # only needed if isddcon=1, but do anyway
cc        do ixu = 1, nx  #scopy doesnt work for static build - integer issue?
cc          matwallo(ixu) = matwallog(ixmin(mype+1)-1+ixu)
cc          matwalli(ixu) = matwallig(ixmin(mype+1)-1+ixu)
cc        enddo
        do igl = 1, ngsp
          call scopy (nx,fngysog(ixmin(mype+1),igl),1,fngyso(1,igl),1)
          call scopy (nx,fngysig(ixmin(mype+1),igl),1,fngysi(1,igl),1)
          call scopy (nx,albedoog(ixmin(mype+1),igl),1,albedoo(1,igl),1)
          call scopy (nx,albedoig(ixmin(mype+1),igl),1,albedoi(1,igl),1)
        enddo
      endif

cdb_initpll         write(6,*) 'nx_loc, ny_loc, mype =', nx,ny,mype
cdb_initpll         write(6,*) 'ixmnbcl, ixmxbcl, mype =', ixmnbcl,ixmxbcl,mype
cdb_initpll         write(6,*) 'iymnbcl, iymxbcl, mype =', iymnbcl,iymxbcl,mype
cdb_initpll         call flush(6)

c ...  Receive variables and grid info on each processor
c_mpi       call recvglobal(iv_totv, iv_totrz, mype+1, vrsendlv_mype,vrsendlz_mype)

c ... Fix global data at two X-point cells to get correct up eqn; also
c ... fill local pe0 arrays with yl map
c_mpi       if (mype == 0) then
c_mpi         if(isddcon == 1) call sendglobal_xpt(iv_totv,iv_totrz,
c_mpi     .                                vrsendlv_mype,vrsendlz_mype)
ccc         do iv = 1, neq_locg(1)
ccc            ivloc2sdgl(iv) = ivloc2sdg(iv,1)
ccc            ivloc2mdgl(iv) = ivloc2mdg(iv,1)
ccc         enddo
c_mpi       endif
c_mpi       if (mype+1 == idxpt(1) .and. isddcon == 1) then
c_mpi        call recvglobal_xpt(iv_totv,iv_totrz,idxpt(1),
c_mpi     .                                    vrsendlv_mype,vrsendlz_mype)
c_mpi       elseif (mype+1 == idxpt(2) .and. isddcon == 1) then
c_mpi          call recvglobal_xpt(iv_totv,iv_totrz,idxpt(2),
c_mpi     .                                    vrsendlv_mype,vrsendlz_mype)
c_mpi       endif

       iysptrx1(1) = max(iysptrxg(mype+1),0)
       iysptrx1(1) = min(iysptrx1(1), ny+1)
       iysptrx2(1) = iysptrx1(1)   #LIKELY ONLY GOOD FOR SINGLE NULL
       iysptrx = iysptrx1(1)
       ixpt1(1) = max(0, ixpt1g(mype+1))
       ixpt1(1) = min(ixpt1(1), nx)
       ixpt2(1) = min(nx+1, ixpt2g(mype+1))
       ixpt2(1) = max(0, ixpt2(1))

ccc         do iy = 0, ny+1
ccc         write(6,*) 'aft rm(3,iy,0), iy, id =',rm(3,iy,0), iy, mype+1
ccc         enddo

c ...  Do initial call to sendbdry_asz to get each processor contribution
c ...  to the array sizes iv_totbdy; then recvbdry_asz receives the
c ...  appropriate data on each processor
c ...  The communications are nonblocking in the isendrecv subroutine
c       call sendbdry_asz(mype+1,myiv)
cdb_initpll       write(6,*) ' Just send sendbdry_asz from mype =', mype,'myiv ',myiv
cdb_initpll       call flush(6)
c       call recvbdry_asz(mype+1,myiv)
       call isendrecv_bdry_asz(mype+1,myiv,a_sz)

cdb_initpll      call MPI_BARRIER(uedgeComm, ierr)
cdb_initpll      write(6,*) "[",mype,"] exit init_pll ***** "
cdb_initpll      call flush(6)

       DEALLOCATE(visendl_mype)
       DEALLOCATE(vrsendlv_mype)
       DEALLOCATE(vrsendlz_mype)

      return
      end
c ****** End of subroutine init_pll **********************************
c-----------------------------------------------------------------------
      subroutine init_par_meshg

c **- This rountine for parallel cases builds the R,Z global mesh on
c     PE0 and then sends the resulting core area (com.area_core) to the
c     other processors

      implicit none

      Use(Dim)            # nx,ny
      Use(Comgeo)         # area_com
      Use(Npes_mpi)       # mype
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
       
c_mpi      integer(int4) :: ierr

      if (mype == 0) then
        call globalmesh
      endif

c ... Now broadcast results of com.area_core to all processors
c_mpi       call MPI_BCAST(area_core,1,MPI_DOUBLE_PRECISION,0,uedgeComm,ierr)

      return
      end
c ***** End of subroutine init_par_meshg *****************************      
c ---------------------------------------------------------------------c

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
                    # label,svrpkg, istgon
      Use(Ident_vars)          # uedge_ver
      Use(Lsode)    # iterm,icntnunk
      Use(Grid)     # ngrid,inewton,imeth,nurlx,ijac,ijactot
      Use(Decomp)   # ubw,lbw
      Use(Share)    # igrid
      Use(Interp)   # isnintp,nxold,nyold
      Use(RZ_grid_info)  # rm,zm
      Use(Indices_domain_dcl)  # nx_loc,ny_loc
      Use(Indices_domain_dcg)  # ndomain
      Use(Npes_mpi)            # npes,mype,ispmion
      Use(UEint)               # isallloc
      Use(Rccoef)              # isoutwall
      Use(Coefeq)              # oldseec, override, cftiexclg
      Use(Flags)               # iprint
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer ifake  #forces Forthon scripts to put implicit none above here
      integer my_pe, n_pes, ixx, icall
CC c_mpi      integer ierr, ireturn
c_mpi      integer ireturn
c_mpi      integer(int4) :: ierr, myfoo
      data icall/0/

      call exmain_prelims

c=======================================================================
c//computation//
c_mpi         call MPI_BARRIER(uedgeComm, ierr)

      do 100 igrid = 1, ngrid

*     -- allocate memory for arrays --
*     -- set ifexmain=1 so that allocate knows exmain is the calling
*     -- routine.  This is necessary for grid sequencing coding.
         if (ismpion.eq.0) then   # MPI parallel version from serial version
           ifexmain = 1
           call allocate
           ifexmain = 0
	   if ((icall == 0) .and. (iprint .ne. 0)) write(*,*) 'UEDGE ',uedge_ver
           icall = 1
         elseif (ismpion.eq.1 .and. icall==0) then
           call init_pll
	   write(*,*) 'UEDGE version ',uedge_ver
           icall = 1
         endif

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





c_mpi         call MPI_BARRIER(uedgeComm, myfoo)

         imeth = inewton(igrid)
         ijac(igrid) = 0
         nurlxn = cnurn*nurlx
         nurlxu = cnuru*nurlx
         nurlxe = cnure*nurlx
         nurlxi = cnuri*nurlx
         nurlxg = cnurg*nurlx
         nurlxp = cnurp*nurlx
c_mpi         call MPI_BARRIER(uedgeComm, myfoo)

*     -- For the continuation mode (icntnunk=1), be sure a Jacobian was
*     -- calculated on the previous step, i.e., ijac > 0
         if (icntnunk==1 .and. ijactot<=1 .and. svrpkg=='nksol') then
            call xerrab('**Error: need initial Jacobian-pair for icntnunk=1')
         endif

c     -- Reinitialize ijactot if icntnunk = 0; prevents ijactot=2 by 2 exmain
c     .. nksol issue
         if (icntnunk == 0) ijactot = 0

c     -- call principal driver routine --
c_mpi         call MPI_BARRIER(uedgeComm, ierr)
         call uedriv

c  -- create the interpolants for the grid sequencing or restarting --
c  -- only update save variables (call gridseq) for nksol if a root
c  -- has been found, i.e., if iterm=1; other convergence failures never
c  -- reach here (kaboomed).

c ...    If a parallel run, send and gather data to PE0 first
         if (ismpion == 1) then
           if (isparmultdt == 0) call build_global_soln
           call gridseq
         endif

         #fill bndry flux arrays for export 
         if (ismpion==0 .and. isoutwall==1) call outwallflux  

         if (isnintp .eq. 0 .and. mype <= 0) then
            nxold = nx
            nyold = ny
            call gchange("Interp",0)
            if (((svrpkg .ne. 'nksol') .and. (svrpkg.ne."petsc")).or. iterm .eq. 1) then
               call gridseq
               call comp_vertex_vals  # gen plasma/neut values at rm,zm(,,4)
            endif
         elseif (isnintp .eq. 1 .and. mype <= 0) then
            if (((svrpkg .ne. 'nksol') .and. (svrpkg.ne."petsc")) .or. iterm .eq. 1) then
               nxold = nx
               nyold = ny
               if(ismpion == 0) then
                  call gchange("Interp",0)
                  call gridseq
                  call comp_vertex_vals  # gen plasma/neut values at rm,zm(,,4)
               endif
            endif
        if (iprint .ne. 0) write(6,*) "Interpolants created; mype =", mype
         endif

  100 continue

      igrid = min(igrid, ngrid)  # prevents igrid>1 problem for ngrid=1

      return
      end
c --** End of subroutine exmain *********************************

c-----------------------------------------------------------------------
      subroutine uedge_petscInsertOpts
c     subroutine simply calls PetscInsertOptions in uedge/svr/petsc-uedge.F90
c     This is just a device to allow calling from the Python or CXX wrapper
cpetsc    
cpetsc      call PetscInsertOpts()
      return
      end
c --** End of subroutine uedge_petscFinal
c
c-----------------------------------------------------------------------
      subroutine uedge_petscFinal

c     Subroutine simply calls PetscFinalizeWrap in uedge/svr/petsc-uedge.F90
c     This is just a device to allow calling from the Python or CXX wrapper

cpetsc      call PetscFinalizeWrap()
      return
      end
c --** End of subroutine uedge_petscFinal

c-----------------------------------------------------------------------
      subroutine uedge_mpiFinal

c     Subroutine simply calls MPI_Finalize()
c     This is just a device to allow calling from the Python or CXX wrapper
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

c_mpi      integer(int4) ierr
c_mpi      call MPI_Barrier(uedgeComm, ierr)
c_mpi      call MPI_Finalize(ierr)
      return
      end
c --** End of subroutine uedge_mpiFinal

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
c-----------------------------------------------------------------------
      subroutine gather_pll_soln

c     GATHER_SOLN combines primary variables from multiple processors
c     into nisg etc. and then copies to resized ni etc. on PE0

      implicit none

      Use(Err_msg_out)         # errmsgflag,errunit
      Use(Dim)                 # nx,ny,nxm,nym
      Use(Npes_mpi)            # mype
      Use(Indices_domain_dcg)  # ndomain,ixpt1g,ixpt2g,iysptrxg
      Use(Indices_domain_dcl)  # ixpt1l,ixpt2l,iysptrxl,ixlbl,ixrbl
      Use(Global_vars)         # nisg,upsg,...
      Use(Parallv)             # nxg,nyg
      Use(Compla)              # ni,up,...
      Use(Imprad)              # isimpon
      Use(Interp)              # nis,ups,...,nxold,nyold
      Use(RZ_grid_info)        # rm,zm
      Use(RZ_grid_global)      # rmg,zmg
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx1,iysptrx2,iysptrx
      Use(Comgeo)              # lcon,lcone,lconi
      Use(Comgeo_g)            # lcong,lconeg,lconig

c ... Local variables
      integer id,ifld,igsp,igmh

c ... Received and manipulate data on PE0
      do id = 1, ndomain
         call recvloc_glob(id)
      enddo

c ... Resize arrays and copy data to appropriate arrays
      nx = nxg
      ny = nyg
      nxm = nx
      nym = ny
      call gallot("RZ_grid_info",0)
      call gchange("Compla",0)
      call gchange("Imprad",0)
      call gchange("Comgeo",0)
      nxold = nxg
      nyold = nyg
      ndomain_orig = ndomain
      ndomain = 1
      call gchange("Interp",0)
c ... Copy nisg -> ni, upsg -> up, etc
      do  ifld = 1, nisp
          call s2copy (nx+2, ny+2, nisg(0:,0:,ifld), 1, nx+2,
     .            ni(0:,0:,ifld), 1, nx+2)
          call s2copy (nx+2, ny+2, nisg(0:,0:,ifld), 1, nx+2,
     .            nis(0:,0:,ifld), 1, nx+2)
      enddo

      do ifld = 1, nusp
         call s2copy (nx+2, ny+2, upsg(0:,0:,ifld), 1, nx+2,
     .            up(0:,0:,ifld), 1, nx+2)
         call s2copy (nx+2, ny+2, upsg(0:,0:,ifld), 1, nx+2,
     .            ups(0:,0:,ifld), 1, nx+2)
      enddo

      do igsp = 1, ngsp
         call s2copy (nx+2, ny+2, ngsg(0:,0:,igsp), 1, nx+2,
     .            ng(0:,0:,igsp), 1, nx+2)
         call s2copy (nx+2, ny+2, ngsg(0:,0:,igsp), 1, nx+2,
     .            ngs(0:,0:,igsp), 1, nx+2)
      enddo

      call s2copy (nx+2, ny+2, tesg, 1, nx+2, te, 1, nx+2)
      call s2copy (nx+2, ny+2, tisg, 1, nx+2, ti, 1, nx+2)
      call s2copy (nx+2, ny+2, phisg, 1, nx+2, phi, 1, nx+2)
      call s2copy (nx+2, ny+2, tesg, 1, nx+2, tes, 1, nx+2)
      call s2copy (nx+2, ny+2, tisg, 1, nx+2, tis, 1, nx+2)
      call s2copy (nx+2, ny+2, phisg, 1, nx+2, phis, 1, nx+2)
      if (isimpon.gt.0) then
         call s2copy (nx+2, ny+2, afracsg, 1, nx+2, afrac, 1, nx+2)
         call s2copy (nx+2, ny+2, afracsg, 1, nx+2, afracs, 1, nx+2)
      endif

      do igmh = 0, 4  # mesh (R,Z) coordinates set to global
         call s2copy (nx+2, ny+2, rmg(0:,0:,igmh), 1, nx+2,
     .            rm(0:,0:,igmh), 1, nx+2)
         call s2copy (nx+2, ny+2, zmg(0:,0:,igmh), 1, nx+2,
     .            zm(0:,0:,igmh), 1, nx+2)
         call s2copy (nx+2, ny+2, psig(0:,0:,igmh), 1, nx+2,
     .            psi(0:,0:,igmh), 1, nx+2)
      enddo
      call s2copy (nx+2,ny+2,lcong,1,nx+2,lcon,1,nx+2)
      call s2copy (nx+2,ny+2,lconeg,1,nx+2,lcone,1,nx+2)
      call s2copy (nx+2,ny+2,lconig,1,nx+2,lconi,1,nx+2)

c ... Store some local values of indices for sub reset_pe0_vars
      ixpt1l = ixpt1(1)
      ixpt2l = ixpt2(1)
      iysptrx1l = iysptrx1(1)
      ixlbl = ixlb(1)
      ixrbl = ixrb(1)

c ... Reset some indices to global values for PE0 (mype=0)
      ixpt1(1) = ixpt1g(1)
      ixpt2(1) = ixpt2g(1)
      iysptrx1(1) = iysptrxg(1)  # limited to single-null
      iysptrx2(1) = iysptrxg(1)
      iysptrx = iysptrx1(1)
      ixlb(1) = 0      # next two lines limited to single-null geo
      ixrb(1) = nx

      return
      end
c------------** End of subroutine gather_pll_soln ----------------------
c-----------------------------------------------------------------------
      subroutine reset_pe0_vars

c     Reset variables on PE0 back to the original domain-decomposed values
c     that were overwritten by gather_pll_soln in order to do a dump of the
c     global data to a hdf5 file; thus, the reverse of sub gather_pll_soln

      implicit none

      Use(Err_msg_out)         # errmsgflag,errunit
      Use(Dim)                 # nx,ny,nxm,nym
      Use(Npes_mpi)            # mype
      Use(Indices_domain_dcg)  # ndomain
      Use(Indices_domain_dcl)  # nx_loc,ny_loc,ixpt1l,ixpt2l,iysptrxl,
                               # ixlbl,ixrbl
      Use(Global_vars)         # nisg,upsg,...
      Use(Parallv)             # nxg,nyg
      Use(Compla)              # ni,up,...
      Use(Imprad)              # isimpon
      Use(Interp)              # nis,ups,...,nxold,nyold
      Use(RZ_grid_info)        # rm,zm
      Use(RZ_grid_global)      # rmg,zmg
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx1,iysptrx2,iysptrx
      Use(Comgeo)              # lcon,lcone,lconi
      Use(Comgeo_g)            # lcong,lconeg,lconig

c ... Local variables
      integer id,ifld,igsp,igmh

c ... Resize arrays and copy data to appropriate arrays
      nx = nx_loc
      ny = ny_loc
      nxm = nx
      nym = ny
      call gallot("RZ_grid_info",0)
      call gchange("Compla",0)
      call gchange("Imprad",0)
      call gchange("Comgeo",0)
      nxold = nx
      nyold = ny
      ndomain = ndomain_orig
      call gchange("Interp",0)
c ... Copy nisg -> ni, upsg -> up, etc
      do  ifld = 1, nisp
          call s2copy (nx+2, ny+2, nisg(0:,0:,ifld), 1, nx+2,
     .            ni(0:,0:,ifld), 1, nx+2)
          call s2copy (nx+2, ny+2, nisg(0:,0:,ifld), 1, nx+2,
     .            nis(0:,0:,ifld), 1, nx+2)
      enddo

      do ifld = 1, nusp
         call s2copy (nx+2, ny+2, upsg(0:,0:,ifld), 1, nx+2,
     .            up(0:,0:,ifld), 1, nx+2)
         call s2copy (nx+2, ny+2, upsg(0:,0:,ifld), 1, nx+2,
     .            ups(0:,0:,ifld), 1, nx+2)
      enddo

      do igsp = 1, ngsp
         call s2copy (nx+2, ny+2, ngsg(0:,0:,igsp), 1, nx+2,
     .            ng(0:,0:,igsp), 1, nx+2)
         call s2copy (nx+2, ny+2, ngsg(0:,0:,igsp), 1, nx+2,
     .            ngs(0:,0:,igsp), 1, nx+2)
      enddo

      call s2copy (nx+2, ny+2, tesg, 1, nx+2, te, 1, nx+2)
      call s2copy (nx+2, ny+2, tisg, 1, nx+2, ti, 1, nx+2)
      call s2copy (nx+2, ny+2, phisg, 1, nx+2, phi, 1, nx+2)
      call s2copy (nx+2, ny+2, tesg, 1, nx+2, tes, 1, nx+2)
      call s2copy (nx+2, ny+2, tisg, 1, nx+2, tis, 1, nx+2)
      call s2copy (nx+2, ny+2, phisg, 1, nx+2, phis, 1, nx+2)
      if (isimpon.gt.0) then
         call s2copy (nx+2, ny+2, afracsg, 1, nx+2, afrac, 1, nx+2)
         call s2copy (nx+2, ny+2, afracsg, 1, nx+2, afracs, 1, nx+2)
      endif

      do igmh = 0, 4  # mesh (R,Z) coordinates set to global
         call s2copy (nx+2, ny+2, rmg(0:,0:,igmh), 1, nx+2,
     .            rm(0:,0:,igmh), 1, nx+2)
         call s2copy (nx+2, ny+2, zmg(0:,0:,igmh), 1, nx+2,
     .            zm(0:,0:,igmh), 1, nx+2)
         call s2copy (nx+2, ny+2, psig(0:,0:,igmh), 1, nx+2,
     .            psi(0:,0:,igmh), 1, nx+2)
      enddo
      call s2copy (nx+2,ny+2,lcong,1,nx+2,lcon,1,nx+2)
      call s2copy (nx+2,ny+2,lconeg,1,nx+2,lcone,1,nx+2)
      call s2copy (nx+2,ny+2,lconig,1,nx+2,lconi,1,nx+2)

c ... Reset some indices to local values for PE0 (mype=0)
      ixpt1(1) = ixpt1l
      ixpt2(1) = ixpt2l
      iysptrx1(1) = iysptrx1l  # limited to single-null
      iysptrx2(1) = iysptrx1l
      iysptrx = iysptrx1l
      ixlb(1) = ixlbl      # next two lines limited to single-null geo
      ixrb(1) = ixrbl

      return
      end
c------------** End of subroutine reset_pe0_vars ----------------------
c-----------------------------------------------------------------------
      subroutine build_global_soln

c     Subroutine sends all PE plasma data (ni, up, ...) to PE0 and 
c     then these are copied into save-variables nis, ups, etc.

      implicit none

      Use(Npes_mpi)            # mype

      call sendloc_glob(mype+1)   # send data to PE0
      if (mype == 0) then           
        call gather_pll_soln      # recieve data put in ni, up, ...
        call gridseq              # copy datat to nis, ups, ...
      endif

      return
      end
c--------------** End of subroutine build_global_soln ----------------

       SUBROUTINE read_zag
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     READS INPUT DATA for Zagorski 2-D edge code
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      Use (Dim)
      Use (Zag_output)
      REAL XP,XS,XK,YP,YK,PI,EC,XEM,XMI,ZM,AMU,AMZ,SING,RT,RP,RS,BT,
     &   Q,DIFA,VINW,GAMAD,GAMAW,DIFH,VINH,ETAY0,CONDE,CONDI,DIFE,
     &   AREA,TME
      INTEGER ILZER,LNST1,I,J,ILN,L

      call gallot("Zag_output",0)


C
      OPEN(30,FILE= 'zag_file',STATUS='OLD')
C
ccc       CALL PLOTDATA

      READ(30,*) XP,XS,XK,YP,YK,IMX,IMY,PI,EC,XEM,XMI,ZM,AMU,AMZ,
     &           SING,RT,RP,RS,BT,Q,DIFA,VINW,GAMAD,GAMAW,DIFH,VINH,
     &           ETAY0,CONDE,CONDI,DIFE,AREA,ILZER,LNST,WX,WY,RMACH,
     &           RECYC,ZRECYC,SNZ,SN,SQE,SQI,FE0,FI0,FLIME,FLIMI,
     &           INH,INZ,IHF,LST
      LNST1=LNST+1
      READ(30,100) (Xzag(I),I=1,IMX+2)
      READ(30,100) (Yzag(J),J=1,IMY+2)
      READ(30,100) ((NEzag(I,J),I=1,IMX+1),J=1,IMY+1)
      READ(30,100) ((NIzag(I,J),I=1,IMX+1),J=1,IMY+1)
      READ(30,100) ((TIzag(I,J),I=1,IMX+1),J=1,IMY+1)
      READ(30,100) ((TEzag(I,J),I=1,IMX+1),J=1,IMY+1)
      READ(30,100) ((VIPARzag(I,J),I=1,IMX),J=1,IMY)
      READ(30,100) ((UIzag(I,J),I=1,IMX),J=1,IMY)
      READ(30,100) ((VEPARzag(I,J),I=1,IMX),J=1,IMY)
      READ(30,100) ((UEzag(I,J),I=1,IMX),J=1,IMY)
      READ(30,100) ((CURXzag(I,J),I=1,IMX),J=1,IMY)
      READ(30,100) ((CURYzag(I,J),I=1,IMX),J=1,IMY)
      READ(30,100) (((N0zag(I,J,ILN),I=1,IMX+1),J=1,IMY+1),ILN=1,3),
     &             (((V0zag(I,J,ILN),I=1,IMX),J=1,IMY),ILN=1,3),
     &             (((U0zag(I,J,ILN),I=1,IMX),J=1,IMY),ILN=1,3),
     &             (((T0zag(I,J,ILN),I=1,IMX+1),J=1,IMY+1),ILN=1,3)
      IF(INZ.EQ.0) THEN
      READ(30,100) (((NZzag(I,J,L),L=1,LNST1),J=1,IMY+1),I=1,IMX+1)
      READ(30,100) (((VZPARzag(I,J,L),L=1,LNST1),I=1,IMX),J=1,IMY)
      READ(30,100) (((UZzag(I,J,L),L=1,LNST1),I=1,IMX),J=1,IMY)
      READ(30,100) ((ZEFFzag(I,J),I=1,IMX+1),J=1,IMY+1)
      READ(30,100) ((TZ0zag(I,J),I=1,IMX+1),J=1,IMY+1)
      ENDIF
      READ(30,*) TME
C
      IF(INZ.EQ.0) THEN
      READ(30,100) (YIELH(J),(YIELZ(J,L),L=1,LST+1),J=1,IMY+1)
      ENDIF
  100 FORMAT(1P,10E18.8)
  110 FORMAT(E18.8)
      CLOSE(30)
C
C     SLAB GEOMETRY
C     Bp/Btot=0.1
C
      DO 5 I=1,IMX+2
      DO 6 J=1,IMY+2
      RXzag(I,J)=1.
      RYzag(I,J)=1.
      BRATIOzag(I,J)=0.1
      GGzag(I,J)=1.
   6  CONTINUE
   5  CONTINUE
      DO 138 J=1,IMY
      DO 38 I=1,IMX
      VIzag(I,J)=VIPARzag(I,J)*BRATIOzag(I,J)
      VEzag(I,J)=VEPARzag(I,J)*BRATIOzag(I,J)
      DO 48 L=1,LNST+1
      VZzag(I,J,L)=VZPARzag(I,J,L)*BRATIOzag(I,J)
   48 CONTINUE
   38 CONTINUE
  138 CONTINUE

ccc      CALL CONST
ccc      CALL CONDCOEF
ccc      CALL EFIELD

      RETURN
      END
c ** End of read_zag subroutine **** ------------------------------------
c -----------------------------------------------------------------------
      subroutine write_profs

*     Writes a file out_ue_s of plasma profiles to be read by the parallel
*     UEDGE on the T3E

      implicit none

      Use(Dim)           # nx,ny,nisp,ngsp
      Use(Compla)        # ni,up,te,ti,ng,phi

c     Local variables
      integer ix,iy,ifld,igsp,ip,basopen

      external basopen,basclose

*  -- set output label --

      ip = basopen("out_ue_s","w")

      do iy = 0, ny+1
         do ix = 0, nx+1
            do ifld = 1, nisp
               if (abs(up(ix,iy,ifld)).lt.1e-99) up(ix,iy,ifld)=1e-99
            enddo
         enddo
      enddo

      write(ip,131) (((ni(ix,iy,ifld),ix=0,nx+1),iy=0,ny+1),
     .                                            ifld=1,nisp)
      write(ip,131) (((up(ix,iy,ifld),ix=0,nx+1),iy=0,ny+1),
     .                                            ifld=1,nisp)
      write(ip,131) ((te(ix,iy),ix=0,nx+1),iy=0,ny+1)
      write(ip,131) ((ti(ix,iy),ix=0,nx+1),iy=0,ny+1)
      write(ip,131) (((ng(ix,iy,igsp),ix=0,nx+1),iy=0,ny+1),
     .                                            igsp=1,ngsp)
      write(ip,131) ((phi(ix,iy),ix=0,nx+1),iy=0,ny+1)

 131  format(1P,8e14.5)

      call basclose(ip)

      return
      end
c ** End of subroutine write_profs **** ----------------------------------
c ------------------------------------------------------------------------
      subroutine read_profs

*     Reads output file out_ue from the T3E for plasma profiles

      implicit none

      Use(Dim)           # nx,ny,nisp,ngsp
      Use(Interp)        # nis,ups,tes,tis,ngs,phis

c     Local variables
      integer ix,iy,ifld,igsp,ip,basopen

      external basopen,basclose

*  -- set output label --

      ip = basopen("out_ue","r")

      read(ip,131) (((nis(ix,iy,ifld),ix=0,nx+1),iy=0,ny+1),
     .                                            ifld=1,nisp)
      read(ip,131) (((ups(ix,iy,ifld),ix=0,nx+1),iy=0,ny+1),
     .                                            ifld=1,nisp)
      read(ip,131) ((tes(ix,iy),ix=0,nx+1),iy=0,ny+1)
      read(ip,131) ((tis(ix,iy),ix=0,nx+1),iy=0,ny+1)
      read(ip,131) (((ngs(ix,iy,igsp),ix=0,nx+1),iy=0,ny+1),
     .                                            igsp=1,ngsp)
      read(ip,131) ((phis(ix,iy),ix=0,nx+1),iy=0,ny+1)

 131  format(1P,8e14.5)

      call basclose(ip)

      return
      end
c ** End of subroutine read_profs **** ----------------------------------
c *** -------------------------------------------------------------------
      subroutine write_profs_boris(filename)

*     Writes a file out_ue_s of plasma profiles to be read by the parallel
*     UEDGE on the T3E

      implicit none

      Use(Dim)           # nx,ny,nisp,ngsp
      Use(Compla)        # ni,up,te,ti,ng,phi
      Use(Phyvar)        # ev
      Use(Share)         # nxomit,igrid,geometry,isnonog,nyomitmx
      Use(RZ_grid_info)  # rm,zm

c     Local variables
      integer ix,iy,ifld,igsp,ip,basopen
      character filename*32

      external basopen,basclose

*  -- set output label --

      ip = basopen(filename,"w")

      do iy = 0, ny+1
         do ix = 0, nx+1
            do ifld = 1, nisp
               if (abs(up(ix,iy,ifld)).lt.1e-99) up(ix,iy,ifld)=1e-99
            enddo
         enddo
      enddo

      write(ip,132)
cc      write(ip, FMT = '(6F11.6, 7ES13.5)') (( 0., 0., 0.,
cc     .                 rm(nxomit+ix,iy,0), zm(nxomit+ix,iy,2),
cc     .                 zm(nxomit+ix,iy,0),
cc     .                 te(ix,iy)/ev, ti(ix,iy)/ev,
cc     .                 up(ix,iy,1), 0., 0.,
cc     .                 ni(ix,iy,1), ng(ix,iy,1),
cc     .                                      ix=0,nx+1),iy=0,ny+1)
      write(ip, FMT = '(6F11.6, 10ES13.5)') (( 0., 0., 0.,
     .                 rm(nxomit+ix,iy,0), zm(nxomit+ix,iy,2),
     .                 zm(nxomit+ix,iy,0),
     .                 te(ix,iy)/ev, ti(ix,iy)/ev,
     .                 up(ix,iy,1),  ni(ix,iy,1), 0., 0.,
     .                 ng(ix,iy,1), 0., 0., 0.,
     .                                      ix=0,nx+1),iy=0,ny+1)

 132  format(6x,"s   ",6x,"theta ",6x,"phi",8x,"r",9x,"z_up",7x,"z_nT",8x,
     .       "Te",11x,"Ti",11x,"up",11x,"ni",11x,"v1",11x,"v2",11x,"n0",
     .        11x,"upg",10x,"vg1",10x,"vg2")

      call basclose(ip)

      return
      end
c ** End of subroutine write_profs_boris **** ----------------------------
c ------------------------------------------------------------------------
      subroutine read_profs_boris(filename,ivers)

*     Reads output file out_ue from the T3E for plasma profiles

      implicit none

      Use(Dim)           # nx,ny,nisp,ngsp
      Use(Interp)        # nis,ups,tes,tis,ngs,phis
      Use(Phyvar)        # ev
      Use(Share)         # nxomit,igrid,geometry,isnonog,nyomitmx
      Use(RZ_grid_info)  # rm,zm

c     Input variables
      character filename*32
      integer ivers          # =0 for old case, =1 for new ng case

c     Local variables
      integer ix,iy,ifld,igsp,ip,basopen
      real dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9
      character cdum*16

      external basopen,basclose

*  -- set output label --

      ip = basopen(filename,"r")

cc  First read the header line, but ignore this data
      read(ip, FMT = '(A)') cdum
      if (ivers == 0) then
         read(ip, FMT = '(6F11.6, 7ES13.5)') (( dum1, dum2, dum3,
     .                 rm(nxomit+ix,iy,0),
     .                 zm(nxomit+ix,iy,2), zm(nxomit+ix,iy,0),
     .                 tes(ix,iy), tis(ix,iy),
     .                 ups(ix,iy,1), dum5, dum6,
     .                 nis(ix,iy,1), ngs(ix,iy,1),
     .                                      ix=0,nx+1),iy=0,ny+1)

      else
         read(ip, FMT = '(6F11.6, 10ES13.5)') (( dum1, dum2, dum3,
     .                 rm(nxomit+ix,iy,0),
     .                 zm(nxomit+ix,iy,2), zm(nxomit+ix,iy,0),
     .                 tes(ix,iy), tis(ix,iy),
     .                 ups(ix,iy,1), nis(ix,iy,1), dum5, dum6,
     .                 ngs(ix,iy,1), dum7, dum8, dum9,
     .                                      ix=0,nx+1),iy=0,ny+1)
      endif
      call basclose(ip)

      return
      end
c ** End of subroutine read_profs_boris **** ----------------------------------
c
c ** Next are some subroutines needed for the parallel version only **
c ...........................................................
      SUBROUTINE GLOCFN (NLOC, T, YLOC, GLOC)
C Routine to define local approximate function g, here the same as f.
      REAL NLOC,T,YLOC,GLOC
      DIMENSION YLOC(*), GLOC(*)
C
      Use(Indices_domain_dcg)

      isddcon = 0
      CALL FFUN(NLOC, T, YLOC, GLOC)
      isddcon = 1
C
      RETURN
      END

c ...........................................................
      SUBROUTINE COMMFN (NLOC, T, YLOC)
C Routine to perform communication required for evaluation of g.
      REAL NLOC,T,YLOC
      RETURN
      END

c ...........................................................
      SUBROUTINE GCOMMFN (NLOC, T, YLOC)
C Routine to perform communication required for evaluation of g.
      REAL NLOC,T,YLOC
      RETURN
      END

c ...........................................................
      subroutine set_uedgeComm(COMM)
c...  Sets the communicator
      Use(Npes_mpi)            # hascomm
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      integer, intent(inout) :: COMM
c_mpi      integer(int4) :: ireturn
      integer my_pe, n_pes
c_mpi      uedgeComm = COMM
      hascomm=1
      npes = N_PES()
      mype = MY_PE()
      if (mype .eq. 0) then
!         write (*, *) "set_uedgeComm: UEDGE has communicator, ", uedgeComm,
!     . ", and ", npes, " processors on that."
      endif
cpetsc      call PetscInitWrap()
      return
      end

c ...........................................................
      subroutine uedge_mpiInit
c...  Sets the communicator
      Use(Npes_mpi)            # hascomm
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

c****
c_mpi      integer(int4) :: ireturn
      integer :: lcomm
!      Initialize MPI, Get NPES and MYPE
       if (hascomm==0) then
! c_mpi         write(*,*) "Initializing MPI on rank ", mype, "."
c_mpi         call MPI_Init(ireturn)
! c_mpi         if(ireturn/=0) write(*,*) "ireturn=", ireturn
         lcomm = 0
c_mpi         lcomm = mpi_comm_world
         call set_uedgeComm(lcomm)
! c_mpi         write(*,*) "MPI initialized on rank ", mype, "."
       endif
      return
      end

c ...........................................................
      subroutine uedge_petscInit
c...  Sets the communicator
c****
cpetsc      call PetscInitWrap()
      return
      end

c ...........................................................
      INTEGER FUNCTION MY_PE()
c...  Returns processor number for calling PE.  0 if serial.
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
      Use(Npes_mpi)
      ! integer me, ireturn
c_mpi      integer(int4) ireturn
      integer me
      me = 0
c_mpi      call MPI_Comm_rank(uedgeComm, me, ireturn)
      MY_PE = me
      return
      end

c ...........................................................
      INTEGER FUNCTION N_PES()
c...  Returns total number of processors.  1 if not parallel.
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
      Use(Npes_mpi)
c_mpi      integer(int4) :: ireturn
      integer :: mynpes
      mynpes = 1
c_mpi      call MPI_Comm_size(uedgeComm,mynpes,ireturn)
      N_PES=mynpes
      return
      end

c ------------------------------------------------------------ c
      subroutine setLogFile(filename)
c...  sets the name of the logfile to which to write, also sets up the logger proper
      implicit none

      Use(Logging)
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      character*(*), intent(in) :: filename
      character(7) :: rankstr
      integer :: my_pe, rank

      rank = MY_PE()
      IF (rank>9999) THEN
        WRITE(rankstr,fmt='(i5.5)') rank
      ELSE IF (rank>999) THEN
        WRITE(rankstr,fmt='(i4.4)') rank
      ELSE IF (rank>99) THEN
        WRITE(rankstr,fmt='(i3.3)') rank
      ELSE IF (rank>9) THEN
        WRITE(rankstr,fmt='(i2.2)') rank
      ELSE
        WRITE(rankstr,fmt='(i1.1)') rank
      ENDIF

      logfname = TRIM(filename) // "_" // TRIM(rankstr) // ".log"
      OPEN(UNIT=6, FILE=logfname, STATUS="UNKNOWN")
        return
      end
c... ** End of subroutine setLogFile
c -------------------------------------------------------- c
      subroutine writeToLog(message)
c...  writes the message to the logfile
      implicit none
      character*(*), intent(in) :: message
        WRITE(6, fmt='(a)') trim(message)
c        WRITE(6, fmt='(...I3...O4)') "Last message had length = ",
c     .    len(message), ", with last character = ",
c     .    message(len(message):len(message))
c        WRITE(6, fmt='(O4)') message(len(message):len(message))
	return
      end
c... ** End of subroutine writeToLog ---------------------- c

        subroutine uedge_reset
      implicit none
      Use(Dim)      # nx,ny
      Use(Interp)   # nis,ngs,...
      Use(Oldpla)   # ni0,ng0, ...
c     resets Uedge variables to their state prior to last call to ueAdvance.
c     returns time at beginning of last call to ueAdvance
      nis=ni0
      ups=up0
      ngs=ng0
      tes=te0
      tis=ti0
      return
      end
c... ** End of subroutine uedge_reset -------------------- c
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

      subroutine getixiyloc(ixg,iyg,ixl,iyl,iownit)
c     Determines if the local processor owns the global index (ixg,iyg); if yes, sets iownit to 1
c     and determines the local index pair (ixl,iyl).  
c     iownit is set to 1 only if the data point is in the physical domain, i.e. in an interior cell or 
c     in an active (physical boundary) ghost cell.  Hence only 1 pe will own any point.
      implicit none
      Use(Indices_domain_dcg)       #ixmin,ixmax,iymin,iymax
      Use(Indices_domain_dcl)       #ixmnbcl,ixmxbcl,iymnbcl,iymxbcl
      Use(Npes_mpi)  # mype, ismpion
      integer ixg,iyg,ixl,iyl,iownit
      integer ixltest,ixutest,iyltest,iyutest

      iownit = 0
      if (ismpion .eq. 0) then
         # If not mpi, then yes the one processor owns it, and local index = global index
         iownit = 1
         ixl = ixg
         iyl = iyg
         return
      endif
c     If still here, we aare running mpi.
c     construct global bounds of active cells in this domain,  recalling that
c     the range of cells on this domain is from ixmin to ixmax and iymin to iymax,
c     but the first and/or last cells in this range are inactive for interior cells
c     and active if the cell touches a physical boundary (i.e. for those ghost cells the
c     data is copied from a neighboring cell).   The test for active or not is in the 
c     variables ixmnbcl, ixmxbcl, iymnbcl, iymxbcl
      ixltest = ixmin(mype)+1-ixmnbcl
      iyltest = iymin(mype)+1-iymnbcl
      ixutest = ixmax(mype)-1+ixmxbcl
      iyutest = iymax(mype)-1+iymxbcl
      if ((ixg .ge. ixltest) .and. (ixg .le. ixutest) .and. (iyg .ge. iyltest) 
     1 .and. (iyg .le. iyutest)) then
         iownit = 1
         ixl = ixg-ixmin(mype)
         iyl = iyg-iymin(mype)
      endif
      return

      end

      subroutine getbdyindexlims
c     calculates values of running index corresponding to start and end of various
c     portions of the edge boundary
      implicit none
      Use(Bdy_indexlims) # ib_idiv, ie_idiv, ib_comwall, ie_comwall, ib_odiv, ie_odiv
                                          # ib_opfwall,ie_opfwall, ib_ipfewall, ie_ipfwall
      Use(Parallv)              #nxg,nyg
      Use(Npes_mpi)  # mype, ismpion
      Use(Share)   # geometry
      Use(Dim)           # nx,ny,nisp,ngsp
      Use(Xpoint_indices)  # ixpt1,ixpt2

      integer nx1,ny1,nyglob,nlpfi1,nlpfo1
c     sets limits of running index on various segments of boundary
      if (ismpion .eq. 1) then
         nx1=nxg-1
         ny1=nyg-1
         nyglob=nyg
      else
         nx1=nx-1
         ny1=ny-1
         nyglob=ny
      endif

c     lower private flux is divided into two parts in uedge grid
      nlpfi1=ixpt1(1)-1    # number phys. cells in inner portion of lpf wall-1
      nlpfo1=nx1-ixpt2(1)    # number phys. cells in outer part of lpf wall-1
      if ((geometry == 'snull') .or. (geometry == 'uppersn')) then
c        carry indices for both beginning and end of surfaces for clarity
         ib_idiv = 1   # begin inner divertor
         ie_idiv = nyglob    # end inner divertor
         ib_comwall = ie_idiv+1  # begin common flux wall
         ie_comwall = ib_comwall+nx1 # end common flux wall
         ib_odiv = ie_comwall + 1  # begin outer divertor
         ie_odiv = ib_odiv + ny1  # end outer divertor
         ib_opfwall = ie_odiv+1 # begin outer portion of pf wall
         ie_opfwall = ib_opfwall+nlpfi1 # end outer part pf wall
         ib_ipfwall = ie_opfwall+1  # begin inner part of pf wall
         ie_ipfwall = ib_ipfwall+nlpfo1 # end inner part of pf wall
      endif
      return
      end

      subroutine getixiybdy(lindex,ix,iy,surfacename)
c     returns poloidal and radial indices corresponding to single running
c     index that moves around the UEDGE boundary.  Also returns the name
c     of the bounding surface ("innerdiv", "outerdiv","outerwall","privwall")
c     Single null, index starts at junction of inner divertor and private flux wall,
c      runs across inner divertor, around outer wall, across outer divertor, and around
c      private flux wall to starting point.
c     Upper single null: uedge handles this as mirror image, same index
c      progression as lower single null
c     Double null, index starts in lower inner private flux, progresses
c      across lower inner divertor, up inner wall, across upper inner
c      divertor, across upper private flux, across upper outer divertor,
c      down outer wall, across lower outer divertor, across lower private
c      flux; NOT YET IMPLEMENTED
c     Half double null: not implemented
c
c    Must call getbdyindexlims prior to calling getixiybdy

      implicit none
      Use(Dim)      # nxpt
      Use(Share)    # geometry
      Use(Xpoint_indices) # ixpt1
      Use(Parallv)             # nxg,nyg
      Use(Bdy_indexlims) # ib_idiv, ie_idiv, ib_comwall, ie_comwall, ib_odiv, ie_odiv
                                          # ib_opfwall,ie_opfwall, ib_ipfewall, ie_ipfwall
      integer lindex,ix,iy
      character*(*), intent(out):: surfacename

      if (geometry == 'snull' .or. geometry == 'uppersn') then
c     (Lower) single null, index starts in inner private flux, progressses
c     across inner divertor, up inner wall and around to outer wall,
c     across outer divertor, then around private flux boundary.

         if ((lindex .lt. ib_idiv) .or. (lindex .gt. ie_ipfwall)) then
            call xerrab('*** Index out of bounds ***')
            endif
         if (lindex .le. ie_idiv) then   # inboard divertor
            ix = 0
            iy = lindex  # Recall iy=0 is a thin ghost cell
            surfacename = "innerdiv"
            return
         endif
         if (lindex .le. ie_comwall) then # common flux outer wall
            ix = lindex-ie_idiv
            iy = nyg   # outer radial ghost cell value
            surfacename = "outerwall"
            return
         endif
         if (lindex .le. ie_odiv) then   # outer divertor 
            ix = nxg+1
            iy = ie_odiv-lindex+1
            surfacename = "outerdiv" 
            return
         endif
         if (lindex .le. ie_opfwall) then  # Outer portion private flux wall
            ix = nxg+ib_opfwall-lindex
            iy = 0
            surfacename = "privwall" 
            return
         else   # Inner portion private flux wall
            ix = ixpt1(1)+ib_ipfwall-lindex
            iy = 0
            surfacename = "privwall"
            return
         endif
      endif
      if (geometry .eq. 'dnull') then
         call xerrab ("*** getixiybdy double null not implemented ***")
      endif
      return
      end

      subroutine set2dat2dpoint(darray,ix,iy,val)
c     Sets value of 2D array "darray" at global index point ix,iy to value "val".
c     Assumes that darray is dimensioned (0:nx+1,0:ny+1).
      
      implicit none

      Use(Dim) # nx, ny
      Use(Npes_mpi)  # mype, ismpion

      real darray(0:nx+1,0:ny+1)
      real val
      integer ix,iy,ixl,iyl,iownit
c     If parallel, determine if the pcurernt processor owns ix,iy and if so get 
c      the local indices, and set
      if (ismpion .eq. 1) then
         call getixiyloc(ix,iy,ixl,iyl,iownit)
         if (iownit == 1) then
            darray(ixl,iyl) = val
         endif
      else     #serial
         darray(ix,iy) = val
      endif
      return
      end

      subroutine set1dat1dpoint(darray,lindex,val)
c     Sets value of 1D array "darray" which runs over a boundary side, to value "val"
c     The darray is assumed to run over x values or y values based on the string returned
c     by surfacename

      implicit none

      real darray(1),val
      character*10 surfacename
      integer lindex,ix,iy,ixl,iyl,iownit

      call getixiybdy(lindex,ix,iy,surfacename)
      call getixiyloc(ix,iy,ixl,iyl,iownit)
      if (iownit == 1) then
         if ((surfacename == "innerdiv") .or. (surfacename == "outerdiv")) then
c           running index is y
            darray(iyl+1) = val   # this should work in serial or parallel
         endif
         if ((surfacename == "privwall") .or. (surfacename == "outerwall")) then
c           running index is x
            darray(ixl+1) = val   # this should work in serial or parallel
         endif
      endif
      return
      end

c ------------------------------------------------------------------ c
      real function getat2dpoint(darray,ix,iy)
c     Returns value of 2D array "darray" at global index point ix,iy.
c     Assumes that darray is dimensioned (0:nx+1,0:ny+1).
      implicit none

      Use(Dim) # nx, ny
      Use(Npes_mpi)  # mype, ismpion
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      real darray(0:nx+1,0:ny+1)
      real val
      integer ix,iy,ixl,iyl,iownit
c_mpi      integer typeval, ierr, status(MPI_STATUS_SIZE),pedonor
c_mpi      data typeval/1/

c     If parallel, determine if the current processor owns ix,iy and if so get 
c      the local indices, and set
      if (ismpion .eq. 1) then
         call getixiyloc(ix,iy,ixl,iyl,iownit)
         if (iownit == 1) then
            val = darray(ixl,iyl)
c .....     Now broadcast results of val to all processors
c_mpi       call MPI_BCAST(val,1,MPI_DOUBLE_PRECISION,0,uedgeComm,ierr)
	 else
c_mpi         call MPI_RECV(val, 1, MPI_DOUBLE_PRECISION, pedonor,
c_mpi     .               typeval, uedgeComm, status, ierr)

         endif
      else     #serial
         val = darray(ix,iy)
      endif
      getat2dpoint = val
      return
      end

      real function getat1dpoint(darray,lindex)
c     Gets value of 2D array "darray" assumed dimensioned 0:nx+1,0:ny+1, at value
c     of index "lindex" that runs around periphery
      implicit none

      Use(Dim) # nx, ny

      real darray(0:nx+1,0:ny+1)
      real getat2dpoint
      character*10 surfacename
      integer lindex,ix,iy
      external getat2dpoint

      call getixiybdy(lindex,ix,iy,surfacename)
      getat1dpoint =  getat2dpoint(darray,ix,iy)
      return
      end
!-------------------------------------------------------- c
      integer function ru_active(amumass,znucleus,charge)

!*** Checks if an ion/neutral with a given mass in AMU,  nucleus charge 
!*** in units of the fundamental (-electron) charge, and atomic charge.
!*** Returns 1 if yes and 0 if no.

      implicit none

      Use(Dim)    #nfsp,ngsp
      Use(UEint)  #minu, ziin, znuclin
      Use(Compla) #mg
      Use(Phyvar) #mp
 
      integer amumass,charge,znucleus,ij

      do ij = 1, nfsp
        if ( amumass == minu(ij) .and. znucleus == znuclin(ij) .and.
     .       charge == ziin(ij) ) then
          ru_active = 1   # species is present
          return
        endif
      enddo
      
!*** Check neutral gas (diffusive-only model)
      if (charge == 0) then
        do ij = 1, ngsp
          if (abs(mg(ij)/mp - amumass) < 5.e-2*amumass) then
            ru_active = 1
            return
          endif
        enddo
      endif
 
      ru_active = 0   # no match in do loop, so species not present

      return
      end
!... ** End of function ru_active  -------------------- c
