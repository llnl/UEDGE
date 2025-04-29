c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"


c-----------------------------------------------------------------------
c***************************************************
c Subroutine to initialize an external neutrals model
c***************************************************
      subroutine init_neutrals()

      implicit none

      Use(Ext_neutrals) #extneutopt

      if (extneutopt .eq. 1) then      	#EIRENE
        call init_eirene()
      else if (extneutopt .eq. 2) then 	#DEGAS2
        call init_degas2()
      end if
      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to initialize DEGAS2 neutrals model
c***************************************************
      subroutine init_degas2()

      implicit none

c ... Common blocks:
      Use(Ext_neutrals) #gecmd, geufile

c ... Local variables:
      character*256 :: cmd
      logical :: isechocmdon =.false.  #uses "echo cmd" for debugging purposes

c ... Set up atomic data files
      cmd="./datasetup"
      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
      call system(trim(cmd))

c ... Set problem specific files
      cmd="./problemsetup"
      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
      call system(trim(cmd))

c ... Process geometry file
      cmd=trim(gecmd)//" "//trim(geufile)
      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
      call system(trim(cmd))

c ... Set options for MC run
      cmd="./tallysetup"
      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
      call system(trim(cmd))

c ... If desired, process background file before flighttest


      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to initialize EIRENE neutrals model
c***************************************************
      subroutine init_eirene()

      character*80 :: cmd

      cmd="echo "//"Hello eirene!"
      call system(trim(cmd))

      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to call an external neutrals model
c***************************************************
      subroutine run_neutrals()

      implicit none

c ... Common blocks:
      Use(Ext_neutrals) #extneutopt

      if (extneutopt .eq. -1) then      #spawn external UEDGE
        call run_uedge
      else if (extneutopt .eq. 1) then  #spawn external EIRENE
        call run_eirene
      else if (extneutopt .eq. 2) then 	#spawn external DEGAS2
        call run_degas2
      else 							 	#internal UEDGE
        call uedge_neutrals
      end if

      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to call UEDGE from shell
c***************************************************
      subroutine run_uedge()

      implicit none

c ... Common blocks:
      Use(Ext_neutrals) #uedgesave, uedgefile, uedgescript, uedgecmd

c ... Local variables:
      integer :: s
      character*256 :: cmd
      character*1 :: backslash="\\"  !pgf90, ifort, gfortran
      logical :: isechocmdon =.false.  #uses "echo cmd" for debugging purposes

c ... External functions:
      integer uedge_save_pdb
c      integer uedge_read_pdb

c ... Save current state
      s=uedge_save_pdb(trim(uedgesave))

c ... Execute uedgecmd, read uedgescript, and exit
c      write(*,*) '*************************************'
      write(*,*) "------------------------------------------------"
      write(*,*) "Running a new UEDGE session"
      
c ... Run a new UEDGE session
 	  cmd = trim(uedgecmd//" ""read "//uedgescript//";quit""")
      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
	  call system(trim(cmd))

c ... Read output data from uedgefile
c      s=uedge_read_pdb(trim(uedgefile))

      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to call DEGAS2 neutrals model
c***************************************************
      subroutine run_degas2()

      implicit none

c ... Common blocks:
      Use(Dim)      	#nisp, ngsp, ...
      Use(MCN_dim)      #nstra
      Use(Ext_neutrals)

c ... Local variables:
      integer :: istra
      character*256 :: cmd
	  character*16 :: npextstr, nflightstr(nstra)
      character*1 :: backslash="\\"  !pgf90, ifort, gfortran
      logical :: isechocmdon =.false.  #uses "echo cmd" for debugging purposes

      if(ext_verbose) then
c       write(*,*) '*************************************'
        write(*,*) "------------------------------------------------"
        write(*,*) "Running DEGAS2 Monte Carlo neutrals"
        write(*,*) "... until steady state"	
      endif

c ... Allocate MNC data storage ???
c     call gchange("MCN_sources",0) #wsor, sni, ... see

c ... Write uedge plasma data
      if(isechocmdon) then
        cmd="call writemcnfile("""//trim(bkufile)//""","""//trim(runid)//""")"
        print *,trim(cmd)
      else
	    call writemcnfile(bkufile,runid)	  	  	 
      end if

c ... Process plasma data into degas2 background file format 
	  cmd = trim(bkcmd)//" "//trim(bkufile)
c      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
	  call system(trim(cmd))

c ... Set parameters in degas2 background file
      do istra=1,nstra
	    write(nflightstr(istra),'(i0)') mcnflights(istra)
      enddo
	  cmd = trim(ncsetcmd)//" "//trim(bkdfile)//" "//trim(ncsetvar)
	  if (isechocmdon) then
       cmd=trim(cmd)//" "//backslash//"'"
      else
       cmd=trim(cmd)//" '"
      endif
      cmd = trim(cmd)//trim(nflightstr(1))
      do istra=2,nstra
        cmd = trim(cmd)//" , "//trim(nflightstr(istra))
      enddo
	  if (isechocmdon) cmd=trim(cmd)//backslash
      cmd = trim(cmd)//"'"
c      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *, trim(cmd)
	  call system(trim(cmd))

c ... Run degas2 MC neutrals code
	  if(ismpicmdon) then
	    write(npextstr,'(i0)') npext 
 	    cmd = trim(degas2mpi)
		cmd=trim(mpicmd)//" "//trim(npopt)//" "//trim(npextstr)//" "//trim(cmd)
      else
 	    cmd = trim(degas2cmd)
	  end if
      if (istimecmdon) cmd="time "//trim(cmd)
      if (isechocmdon) cmd="echo "//trim(cmd)
      if (ext_verbose) print *,trim(cmd)
	  call system(trim(cmd))

c ... Read data and convert
      if(get_neutral_sources) then 
c ... Read source data
        if(isechocmdon) then
          cmd="call readmcnsor("""//trim(degas2file)//""")"
          print *,trim(cmd)
        else
          call readmcnsor(degas2file)
        end if
c ...  Convert source data
        if(isechocmdon) then
          cmd="call convertmcnsor"
          print *,trim(cmd)
        else
          call convertmcnsor
        end if
      end if 

      if(get_neutral_moments) then
c ... Generate moment data using outputbrowser
        cmd = trim(degas2outcmd) //" "//trim(degas2outscript)        
c        if (istimecmdon) cmd="time "//trim(cmd)
        if (isechocmdon) cmd="echo "//trim(cmd)
        if (ext_verbose) print *,trim(cmd)
        call system(trim(cmd))
c ... Filter moment data
        cmd="cd "//trim(neut_output_dir)
        cmd=trim(cmd)//"; ../"//trim(degas2outsh)
        cmd=trim(cmd)//"; cd .."
        if (ext_verbose) print *,trim(cmd)
        call system(trim(cmd))
c ... Read moment data
        if(isechocmdon) then
          cmd="call readmcnmoments("""//trim(neut_output_dir)//""")"
          print *,trim(cmd)
        else
          call readmcnmoments(trim(neut_output_dir))
        end if
c ... Convert moment data
        if(isechocmdon) then
          cmd="call convertmcnmoments"
          print *,trim(cmd)
        else
          call convertmcnmoments
        end if
      end if 

      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to call EIRENE neutrals model
c***************************************************
      subroutine run_eirene()

      implicit none

      character*80 :: cmd

      cmd="echo "//"Run eirene!"
      call system(trim(cmd))

      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to solve UEDGE fluid plasma model
c***************************************************

      subroutine uedge_plasma()

      implicit none

c ... Common blocks:
      Use(Dim)      	#nisp, ngsp, nhsp ...
      Use(UEpar)    	#isnion, isupon, ... 
	  Use(Math_problem_size) #neqmx
      Use(Time_dep_nwt) #dtreal
	  Use(PNC_params) 	#dtplasma

c ... Local variables:
      integer isnionold(nispmx),isuponold(nispmx),istionold,isteonold,isphionold
      integer isngonold(ngspmx), isupgonold(ngspmx), nhspold, ngspold
      real dtrealold 

      if(pnc_verbose) then
c        write(*,*) "************************************************"
        write(*,*) "------------------------------------------------"
        write(*,*) "Solving UEDGE plasma model without neutral gas"
        write(*,*) "dtplasma=",dtplasma
      endif

c ... Save original settings
      dtrealold=dtreal
      isnionold=isnion
      isuponold=isupon
      istionold=istion
      isteonold=isteon
      isphionold=isphion
      isngonold=isngon
      isupgonold=isupgon
      nhspold=nhsp
      ngspold=ngsp

c ... Turn on fluid plasma model & turn off neutrals model
c ... The settings below correspond to the original model
c ... Now, assume the desired settings are specified by the user
c     isnion=1			 
c      isupon=1			 
c      istion=1		 
c      isteon=1			 
c      isphion=1	 
c      if (nhsp .gt. 1) then 
c        isnion(iigsp)=0
c        isupon(iigsp)=0
c      endif

c ... Solve fluid plasma model 
      dtreal=dtplasma
      call exmain

c ... Restore original settings  
      dtreal=dtrealold    
      isnion=isnionold
      isupon=isuponold
      istion=istionold
      isteon=isteonold
      isphion=isphionold
      isngon=isngonold
      isupgon=isupgonold
      nhsp=nhspold
      ngsp=ngspold

      end
c-----------------------------------------------------------------------
c***************************************************
c Subroutine to solve UEDGE fluid neutrals model
c***************************************************

      subroutine uedge_neutrals()
c ... This subroutine can be used with explicit PNC, but not with implicit JFNK

      implicit none

c ... Common blocks:
      Use(Dim)      	#nisp, ngsp, ...
      Use(UEint)        #ziin
      Use(UEpar)    	#isnion, isupon, iigsp, ... 
	  Use(Math_problem_size) #neqmx
      Use(Time_dep_nwt)  #dtreal
      Use(MCN_dim)
      Use(MCN_sources)  #ismcnon
      Use(PNC_params) 	#dtneut


c ... Local variables:
      integer isnionold(nispmx),isuponold(nispmx),istionold,isteonold,isphionold
      integer ismcnonold, nhspold, ngspold, ziinold(nispmx)
      integer isngonold(ngspmx), isupgonold(ngspmx)
      real dtrealold 

      if(pnc_verbose) then
c        write(*,*) "************************************************"
        write(*,*) "------------------------------------------------"
        write(*,*) "Solving UEDGE neutral gas model for ng, upg"
        write(*,*) "dtneut=", dtneut
      endif

c ... Save original settings
      dtrealold=dtreal
      isnionold=isnion
      isuponold=isupon
      istionold=istion
      isteonold=isteon
      isphionold=isphion
      ismcnonold=ismcnon

      isupgonold=isupgon
      isngonold=isngon
      ngspold=ngsp
      nhspold=nhsp
      ziinold=ziin

c ... Turn on fluid neutrals model & turn off fluid plasma model
      dtreal=dtneut
      isnion=0
      isupon=0
      istion=0
      isteon=0
      isphion=0
      isnion(iigsp)=1
      isupon(iigsp)=1
      ismcnon=0

      ngsp=1
      nhsp=2
      ziin(nhsp)=0
      isupgon(1)=1	#neutrals model with parallel flow
      isngon=0

c ... Solve fluid neutrals model 
      call exmain

c ... Restore original settings  
      dtreal=dtrealold    
      isnion=isnionold
      isupon=isuponold
      istion=istionold
      isteon=isteonold
      isphion=isphionold
      ismcnon=ismcnonold

      isupgon=isupgonold
      isngon=isngonold
      ngsp=ngspold
      nhsp=nhspold
      ziin=ziinold

      end
c-----------------------------------------------------------------------
c******************************************************
c Subroutine to convert DEGAS2 sources to UEDGE sources
c******************************************************
      subroutine convertmcnsor()
# This function converts the plasma source data in file sources.out
# from DEGAS2 to the proper form for the UEDGE plasma model.

      implicit none

c ... Common blocks:
      Use(Dim)
      Use(MCN_dim)
      Use(MCN_sources) #uesor_ni, ..., mcnsor_ni ...
      Use(MCN_bkgd) 	#v2c,vyc,upc
      Use(Compla) 		#upi,mi,v2,vy
      Use(Bfield) 		#rbfbt
      Use(Comgeo) 		#rr
      Use(RZ_grid_info) #br,bz,bpol
      Use(Selec)  		#ixm1,ixp1


c ... Local Variables:
      integer ix,iy,ifld,istra,ix2
      real smox(0:nx+1,0:ny+1,1:nfl,1:nstra)
      real smop(0:nx+1,0:ny+1,1:nfl,1:nstra)
      real smov(0:nx+1,0:ny+1,1:nfl,1:nstra)
      real seit(0:nx+1,0:ny+1,1:nstra)
      real sinb1, cosb1
c     real sinbp(0:nx+1,0:ny+1),cosbp(0:nx+1,0:ny+1) 

c ... Allocate MCN_bkgd group
      call gchange("MCN_bkgd",0) 

c      cosb1=br(:,:,0)/bpol(:,:,0)
c      sinb1=bz(:,:,0)/bpol(:,:,0)


# compute UEDGE flow velocity at cell centers:
      do ifld=1,nisp
        do ix=1,nx
          do iy=1,ny
            ### 3 orthogonal components of the velocity are:
            v2c(ix,iy,ifld)=(v2(ixm1(ix,iy),iy,ifld)+v2(ix,iy,ifld))/2.
	        vyc(ix,iy,ifld)=(vy(ix,iy-1,ifld)+vy(ix,iy,ifld))/2.
            upc(ix,iy,ifld)=(upi(ixm1(ix,iy),iy,ifld)+upi(ix,iy,ifld))/2.
          enddo
       enddo
      enddo
c      v2c=(v2)/2.
c      vyc(1:nx,1:ny,1:nisp)=(vy(1:nx,0:ny-1,1:nisp)+vy(1:nx,1:ny,1:nisp))/2.
c      upc=(upi)/2.

# Calculate parallel momentum source from components and interpolate
# onto UEDGE velocity mesh; convert from total ion energy source to
# ion thermal energy source:

      smox=0.;
      smop=0.;
      smov=0.;
      seit=sei;
      do istra=1,nstra
        do ifld=1,nfl
          do iy=1,ny
            do ix=1,nx
              cosb1=br(ix,iy,0)/bpol(ix,iy,0)
              sinb1=bz(ix,iy,0)/bpol(ix,iy,0)
              ix2=ixp1(ix,iy)
      	      # poloidal component (on UEDGE density mesh) is:
      	      smox(ix,iy,ifld,istra) = smor(ix,iy,ifld,istra) * cosb1
     .                                  + smoz(ix,iy,ifld,istra) * sinb1
              # parallel component (on UEDGE density mesh) is:
              smop(ix,iy,ifld,istra) = smox(ix,iy,ifld,istra) * rr(ix,iy)
     .                                + smophi(ix,iy,ifld,istra) * rbfbt(ix,iy)
              # interpolate parallel component (on UEDGE velocity mesh) is:
              smov(ix,iy,ifld,istra) = 
     .                      0.5*(smop(ix,iy,ifld,istra)+smop(ix2,iy,ifld,istra))
              # contribution to ion thermal energy
              seit(ix,iy,istra) = seit(ix,iy,istra)
     .                           -upc(ix,iy,ifld)*smop(ix,iy,ifld,istra)
     .                  +(0.5*mi(ifld)*upc(ix,iy,ifld)**2)*sni(ix,iy,ifld,istra)
            enddo
          enddo
        enddo
      enddo

#Fill in Border Data
      smov(0,:,:,:)=0.5*smop(1,:,:,:)
      smov(nx+1,:,:,:)=0.5*smop(nx,:,:,:)


# Put data in generic source arrays
      mcnsor_ni=sni
      mcnsor_up=smov 
      mcnsor_te=see
      mcnsor_ti=seit

# Set reference source strengths for each stratum:
### NOTE: This needs to be generalized for multi-species (nfl>1) ###
      mcncurr = sum(sum(mcnsor_ni(:,:,1,:),dim=2),dim=1)

# Initialize plasma source terms; sum over strata:
      if (nisp==2*nfl) then
        uesor_ni(:,:,1:nfl)=sum(mcnsor_ni,dim=4)
        uesor_up(:,:,1:nfl)=sum(mcnsor_up,dim=4)
        uesor_ni(:,:,(nfl+1):2*nfl)=-uesor_ni(:,:,1:nfl)
        uesor_up(:,:,(nfl+1):2*nfl)=-uesor_up(:,:,1:nfl)
      else
        uesor_ni=sum(mcnsor_ni,dim=4)
        uesor_up=sum(mcnsor_up,dim=4)
      endif
      uesor_ti=sum(mcnsor_ti,dim=3)
      uesor_te=sum(mcnsor_te,dim=3)

      end

c----------------------------------------------------------------------c

c******************************************************
c Subroutine to convert components of DEGAS2 vector to UEDGE vector
c******************************************************
c     [Vx]            [   c  0   s   ] [Vr] 
c     [Vy] = R.V    = [  -s  0   s   ] [Vt] 
c     [Vp]            [bp*c  bt  s*bp] [Vz]
c
c     [Vr]            [  c    -s    0   ] [Vr] 
c     [Vt] = R^-1.V = [-bp/bt  0   1/bt ] [Vt] 
c     [Vz]            [  s     c    0   ] [Vz]

      subroutine convertmcnvec(mcvar,uevar,mcvar_rsd,uevar_rsd,sgn)

      implicit none

c ... Common blocks:
      Use(Dim)          #nx,ny
      Use(MCN_dim)      #nfl
      Use(Bfield) 		#rbfbt
      Use(Comgeo) 		#rr
      Use(RZ_grid_info) #br,bz,bpol

c ... Input Variables:
      integer, intent(in) :: sgn     # sign of transformation: +1 for vector, -1 for transpose(vector)
      real, dimension(0:nx+1,0:ny+1,nfl,3), intent(in) :: mcvar,mcvar_rsd

c ... Output Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3), intent(out) :: uevar,uevar_rsd

c ... Local Variables:
      integer :: ifld
      real, dimension(0:nx+1,0:ny+1) :: sinb,cosb
      real, dimension(0:nx+1,0:ny+1,nfl,3) :: mcvar_sd

# Initialize
      mcvar_sd=mcvar_rsd*mcvar
      cosb=br(:,:,0)/bpol(:,:,0)
      sinb=bz(:,:,0)/bpol(:,:,0)

# Calculate field-aligned UEDGE components:
      if (sgn .ge. 0) then
        do ifld=1,nfl
          uevar(:,:,ifld,1) =   cosb*mcvar(:,:,ifld,1) +  sinb*mcvar(:,:,ifld,3)
          uevar(:,:,ifld,2) = - sinb*mcvar(:,:,ifld,1) +  cosb*mcvar(:,:,ifld,3)
          uevar(:,:,ifld,3) =     rr*uevar(:,:,ifld,1) + rbfbt*mcvar(:,:,ifld,2)

          uevar_rsd(:,:,ifld,1) = sqrt((cosb*mcvar_sd(:,:,ifld,1))**2  + (sinb*mcvar_sd(:,:,ifld,3))**2)  
          uevar_rsd(:,:,ifld,2) = sqrt((sinb*mcvar_sd(:,:,ifld,1))**2  + (cosb*mcvar_sd(:,:,ifld,3))**2)
          uevar_rsd(:,:,ifld,3) = sqrt( (rr*uevar_rsd(:,:,ifld,1))**2 + (rbfbt*mcvar_sd(:,:,ifld,2))**2 ) 
        enddo
      else
        do ifld=1,nfl
          uevar(:,:,ifld,1) =   cosb*mcvar(:,:,ifld,1) -  sinb*mcvar(:,:,ifld,2)
          uevar(:,:,ifld,3) =   sinb*mcvar(:,:,ifld,1) +  cosb*mcvar(:,:,ifld,2)
          uevar(:,:,ifld,2) =  (- rr*uevar(:,:,ifld,1) +       mcvar(:,:,ifld,3))/rbfbt

          uevar_rsd(:,:,ifld,1) = sqrt((cosb*mcvar_sd(:,:,ifld,1))**2  + (sinb*mcvar_sd(:,:,ifld,2))**2)  
          uevar_rsd(:,:,ifld,3) = sqrt((sinb*mcvar_sd(:,:,ifld,1))**2  + (cosb*mcvar_sd(:,:,ifld,2))**2)
          uevar_rsd(:,:,ifld,2) = sqrt( (rr*uevar_rsd(:,:,ifld,1))**2 + mcvar_sd(:,:,ifld,3)**2)/abs(rbfbt) 
        enddo
      endif

      where (uevar .ne. 0) uevar_rsd=uevar_rsd/abs(uevar)
      where (uevar .eq. 0) uevar_rsd=1

      end

c----------------------------------------------------------------------c

c******************************************************
c Subroutine to interpolate DEGAS2 vector to UEDGE vector
c******************************************************

      subroutine interpmcnvec(mcvar,uevar,mcvar_rsd,uevar_rsd)
# This function interpolates a vector from the centered grid to the velocity grid

      implicit none

c ... Common blocks:
      Use(Dim)          #nx,ny
      Use(MCN_dim)      #nfl
      Use(Comgeo) 		#dx,dy
      Use(Selec)  		#ixm1,ixp1

c ... Input Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3), intent(in) :: mcvar,mcvar_rsd

c ... Output Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3), intent(out) :: uevar,uevar_rsd

c ... Local Variables:
      integer ix,iy,ifld,ix2,iy2
      real dx2,dy2

# Interpolate to UEDGE velocity mesh:
      do ifld=1,nfl
        do iy=1,ny
          do ix=1,nx
            ix2=ixp1(ix,iy) 
            dx2=dx(ix,iy)+dx(ix2,iy)
            uevar(ix,iy,ifld,1) =     
     .        (mcvar(ix,iy,ifld,1)*dx(ix2,iy)+mcvar(ix2,iy,ifld,1)*dx(ix,iy))/dx2
            uevar_rsd(ix,iy,ifld,1) = 
     .        (mcvar_rsd(ix,iy,ifld,1)*dx(ix2,iy)+mcvar_rsd(ix2,iy,ifld,1)*dx(ix,iy))/dx2

            iy2 = iy+1
            dy2 = dy(ix,iy)+dy(ix,iy2)
            uevar(ix,iy,ifld,2) =
     .        (mcvar(ix,iy,ifld,2)*dy(ix,iy2) + mcvar(ix,iy2,ifld,2)*dy(ix,iy))/dy2
            uevar_rsd(ix,iy,ifld,2) =  
     .        (mcvar_rsd(ix,iy,ifld,2)*dy(ix,iy2) + mcvar_rsd(ix,iy2,ifld,2)*dy(ix,iy))/dy2
          enddo
        enddo
      enddo

      end
c----------------------------------------------------------------------c

c******************************************************************
c Subroutine to convert DEGAS2 vector to UEDGE vector on UEDGE mesh
c******************************************************************

      subroutine convertmcnvector(mcvar,uevar,mcvar_rsd,uevar_rsd)

      implicit none

c ... Common blocks:
      Use(Dim)          #nx,ny
      Use(MCN_dim)      #nfl

c ... Input Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3) :: mcvar,mcvar_rsd

c ... Output Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3) :: uevar,uevar_rsd

c ... Local Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3) :: uevec,uevec_rsd

c      write(*,*) "convertmcnvector: vec"
      call convertmcnvec(mcvar,uevec,mcvar_rsd,uevec_rsd,1)
c      write(*,*) "convertmcnvector: interp" 
      call  interpmcnvec(uevec,uevar,uevec_rsd,uevar_rsd)
c     write(*,*) "convertmcnvector: done"
      end

c----------------------------------------------------------------------c

c******************************************************
c Subroutine to convert DEGAS2 tensor to UEDGE tensor
c******************************************************

      subroutine convertmcntensor(mcvar,uevar,mcvar_rsd,uevar_rsd)
# This function converts the neutral data generated by outputbrowser
# from DEGAS2 to the proper form for the UEDGE plasma model.

      implicit none

c ... Common blocks:
      Use(Dim)
      Use(MCN_dim)

c ... Input Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3,3), intent(in) :: mcvar,mcvar_rsd

c ... Output Variables:
      real, dimension(0:nx+1,0:ny+1,nfl,3,3), intent(out) :: uevar,uevar_rsd

c ... Local Variables:
      integer :: idim 
      real, dimension(0:nx+1,0:ny+1,nfl,3,3) :: uevec,uevec_rsd
      real, dimension(0:nx+1,0:ny+1,nfl,3) :: vec,vec_rsd,var,var_rsd

c ... Calculate field-aligned UEDGE components:
c ...     [Vx] = R.v = [ cos, sin] [Vr]
c ...     [Vy]         [-sin, cos] [Vz]
c ...     ...   etc. 
c ...     [P]_xy = R . P_[rz] . R^T

c ... First multiply by R matrix
      do idim=1,3
        var=mcvar(:,:,:,idim,:)
        var_rsd=mcvar_rsd(:,:,:,idim,:)
        call convertmcnvec(var,vec,var_rsd,vec_rsd,1)
        uevec(:,:,:,idim,:)=vec
        uevec_rsd(:,:,:,idim,:)=vec_rsd
      end do

c ... Then multiply by Transpose(R) matrix
      do idim=1,3
        vec=uevec(:,:,:,idim,:)
        vec_rsd=uevec_rsd(:,:,:,idim,:)
        call convertmcnvec(vec,var,vec_rsd,var_rsd,-1)
        uevar(:,:,:,idim,:)=var
        uevar_rsd(:,:,:,idim,:)=var_rsd
      end do
c      write(*,*) "convertmcntensor: done"

      end

c----------------------------------------------------------------------c
      subroutine scale_mcnsor
      implicit none

      Use(Dim)	# nx,ny,nisp,nusp
      Use(UEpar)      # isupgon,iigsp
      Use(Comflo)	# fnix
      Use(MCN_dim)
      Use(MCN_sources)
      Use(Interp) # nis
      Use(Compla) # ni

      integer ix,iy,ifld,istra

      external remark
c 	  write(*,*) 'TEST SCALE_MCN START: ismcnvar=',ismcnvar

c     This subroutine scales plasma source terms obtained from the 
c     Monte-Carlo-Neutrals code.  These sources are assumed to scale
c     with the neutral source currents at the divertor plates.
c     The flag that controls the scaling is ismcnvar :
c         ismcnvar=0  -->  sources are fixed
c         ismcnvar=1  -->  sources scale with current
c
c     In the terminology of the EIRENE code, a stratum is a surface
c     or volume element where neutral particles originate.  Here
c     we consider only two possible strata: the inboard and outboard
c     divertor plates.  We neglect strata associated with gas puffing
c     or recombination.
c     
c     Plasma source terms from EIRENE are :
c         mcnsor_ni
c         mcnsor_up
c         mcnsor_te
c         mcnsor_ti
c     where mcnsor_ni(ix,iy,ifld,istra) is the particle source
c     for ion fluid 'ifld' at cell (ix,iy) due to neutral source
c     stratum 'istra'.
c
c     mcncurr(istra) is the incident ion current that
c     normalizes the plasma sources due to neutrals
c     that originate at stratum 'istra'.
c
c     The scaled plasma source terms for UEDGE are :
c         uesor_ni
c         uesor_up
c         uesor_te
c         uesor_ti
c     where uesor_ni is obtained by summing the scaled mcnsor_ni
c     over all strata, and similarly for up,te,ti.

ccc=================================================================
ccc     NOTE: For the cmod-box problem there is only ONE strata !!!
ccc     (in general, we need to devise a test to identify strata)
ccc=================================================================

      if (ismcnvar==0) return # Already set by convertmcnsor
c=====else

c     Compute scale factors for each strata :
      strascal=1.	# default is no scaling
      uesor_ni=0.
      uesor_up=0.
      uesor_ti=0.
      uesor_te=0.

      if (ismcnvar==1) then	# sources scale with current
c		 call remark("ismcnvar=1")
c         write(*,*) 'ismcnvar=1' 
         do istra=1,nstra
            uecurr(istra)=0.
            if (istra==1) then	# for east target plate only
               do iy=1,ny
                  uecurr(istra)=uecurr(istra)+fnix(nx,iy,1)
               enddo
c            else
c               call remark("***")
c               call remark("***    subroutine scale_mcn    ***")
c               call remark("***  not defined for nstra > 1  ***")
c               call remark("***")
            endif
            if(mcncurr(istra) > 0) then
               strascal(istra)=uecurr(istra)/mcncurr(istra)
            endif
         enddo
c       Scale source terms and sum over strata :
        do istra=1,nstra
           do iy=0,ny+1
              do ix=0,nx+1
                 do ifld=1,nisp
                    uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*strascal(istra)
                 enddo
                 do ifld=1,nusp
                    uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*strascal(istra)
                 enddo
                 uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*strascal(istra)
                 uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*strascal(istra)
              enddo
           enddo
        enddo
      else if (ismcnvar==2) then	# sources scale with electron density
c	     call remark("ismcnvar=2")
c         write(*,*) 'ismcnvar=2' 
         # Scale all sources to electron (ion) density 
         do istra=1,nstra
            do iy=0,ny+1
               do ix=0,nx+1
                  do ifld=1,nisp
                     uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
                  enddo
                  do ifld=1,nusp
                     uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
                  enddo
                  uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
                  uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
               enddo
            enddo
         enddo
      else if (ismcnvar==3) then	# sources scale with old UEDGE current
c        call remark("ismcnvar=3")
c        write(*,*) 'ismcnvar=3' 
         do istra=1,nstra
            uecurr(istra)=0.
            if (istra==1) then	# for east target plate only
               do iy=1,ny
                  uecurr(istra)=uecurr(istra)+fnix(nx,iy,1)
               enddo
c            else
c               call remark("***")
c               call remark("***    subroutine scale_mcn    ***")
c               call remark("***  not defined for nstra > 1  ***")
c               call remark("***")
            endif
            if(mcncurr(istra) > 0) then
               strascal(istra)=uecurr(istra)/olduecurr(istra)
            endif
         enddo
c       Scale source terms and sum over strata :
        do istra=1,nstra
           do iy=0,ny+1
              do ix=0,nx+1
                 do ifld=1,nisp
                    uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*strascal(istra)
                 enddo
                 do ifld=1,nusp
                    uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*strascal(istra)
                 enddo
                 uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*strascal(istra)
                 uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*strascal(istra)
              enddo
           enddo
        enddo
      else if (ismcnvar==4) then	
         # Sources scale with electron density and target plate current (neutral density)
         # This is a combination of ismcnvar = 2 and 3
c	     call remark("ismcnvar=4")
c         write(*,*) 'ismcnvar=4' 
         # Scale all sources electron (ion) density 
         do iy=0,ny+1
            do ix=0,nx+1
               do ifld=1,nisp
                  uesor_ni(ix,iy,ifld)=0.
               enddo
               do ifld=1,nusp
                  uesor_up(ix,iy,ifld)=0.
               enddo
               uesor_te(ix,iy)=0.
               uesor_ti(ix,iy)=0.
            enddo
         enddo
         do istra=1,nstra
            do iy=0,ny+1
               do ix=0,nx+1
                  do ifld=1,nisp
                     uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
                  enddo
                  do ifld=1,nusp
                     uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
                  enddo
                  uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
                  uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*ni(ix,iy,ifld)/nis(ix,iy,ifld)
               enddo
            enddo
         enddo
      # Scale sources linearly with old UEDGE current to target plate 
      # This option corresponds to ismcnvar=3
         do istra=1,nstra
            uecurr(istra)=0.
            if (istra==1) then	# for east target plate only
               do iy=1,ny
                  uecurr(istra)=uecurr(istra)+fnix(nx,iy,1)
               enddo
c            else
c               call remark("***")
c               call remark("***    subroutine scale_mcn    ***")
c               call remark("***  not defined for nstra > 1  ***")
c               call remark("***")
            endif
            if(mcncurr(istra) > 0) then
               strascal(istra)=uecurr(istra)/olduecurr(istra)
            endif
         enddo
c       Scale source terms and sum over strata :
        do istra=1,nstra
           do iy=0,ny+1
              do ix=0,nx+1
                 do ifld=1,nisp
                    uesor_ni(ix,iy,ifld)=uesor_ni(ix,iy,ifld)
     &                    +mcnsor_ni(ix,iy,ifld,istra)*strascal(istra)
                 enddo
                 do ifld=1,nusp
                    uesor_up(ix,iy,ifld)=uesor_up(ix,iy,ifld)
     &                    +mcnsor_up(ix,iy,ifld,istra)*strascal(istra)
                 enddo
                 uesor_te(ix,iy)=uesor_te(ix,iy)
     &                    +mcnsor_te(ix,iy,istra)*strascal(istra)
                 uesor_ti(ix,iy)=uesor_ti(ix,iy)
     &                    +mcnsor_ti(ix,iy,istra)*strascal(istra)
              enddo
           enddo
        enddo
      endif


      return
      end

c-----------------------------------------------------------------------

c***************************************************
c Subroutine to run coupled plasma + neutrals models
c***************************************************
      subroutine run_pnc()

      implicit none

c ... Common blocks:
      Use(Dim)  		# needed for PNC
      Use(Ext_neutrals) # extneutopt
      Use(MCN_dim)  	# needed for PNC
      Use(PNC_params)	# pnc_opt

      if (extneutopt .eq. 1) then      		#EIRENE
c        call uedge_eirene
      else if (extneutopt .eq. 2) then 		#DEGAS2
        call uedge_degas2
      else 							 		#UEDGE neutral fluid
        call uedge_uedge
      endif

      end

c-----------------------------------------------------------------------

c********************************************************
c Subroutine to run coupled uedge plasma + uedge neutrals
c********************************************************
      subroutine uedge_uedge()

      implicit none

c ... Common blocks:
      Use(Dim) 			# nx, ny, nisp
      Use(Compla) 		# ni, up, ...
      Use(Ext_neutrals)	# extneutopt
      Use(PNC_params)	# pnc_step
      Use(Math_problem_size)
      Use(MCN_dim)		# needed for PNC_data
      Use(PNC_data)		# ni_pnc, ... sni_pnc, ...
      Use(Time_dep_nwt) # dtreal

c ... External functions:
      integer pnc_save_pdb

c ... Local variables:
      character*8 stepstr
      character*256 pnc_save_name

      write(*,*) '*************************************'
      write(*,*) 'UEDGE plasma + UEDGE neutral coupling'

c ... Allocate PNC_data group
      call gchange("PNC_data",0) 

      do pnc_step=pnc_step+1,pnc_maxstep
        pnc_time = pnc_time + dtplasma
        write(*,*) '*************************************'
        write(*,*) 'Step: ',pnc_step
        write(*,*) 'Time: ',pnc_time

        #Neutral Step      Use(PNC_params)	#dtneut
		call store_neutrals
        if (extneutopt .ge. 0) then
          dtreal=dtneut
          call uedge_neutrals
        else
          call run_uedge
        endif
		call update_neutrals

        #Plasma Step
        dtreal=dtplasma
        call store_plasma	
        call uedge_plasma	
c        call exmain  			#should use exmain to be consistent with other pnc routines
        call update_plasma

        call pnc_diagnostics    

        #Save PNC data
        if(mod(pnc_step,pnc_nsave)==0) then
          write(stepstr,'(i0)') pnc_step
          pnc_save_name=trim(pnc_savefile)//trim(stepstr)//'.pdb'
          if(.not.(pnc_save_pdb(pnc_save_name)==0)) 
     .      write(*,*) 'uedge_uedge: error writing '//pnc_save_name
        endif
      enddo
      pnc_step=pnc_step-1  #return value to pnc_maxstep

      end

c-----------------------------------------------------------------------

c*********************************************************************
c Subroutine to run coupled uedge plasma + degas2 monte carlo neutrals
c*********************************************************************
      subroutine uedge_degas2()

      implicit none

c ... Common blocks:
      Use(Dim) 			# nx, ny, nisp
      Use(Compla) 		# ni, up, ...
      Use(Math_problem_size)
      Use(MCN_dim)  	# needed for PNC
      Use(PNC_params)	# pnc_step, dtneut, dtplasma
      Use(Time_dep_nwt) # dtreal

c ... External functions:
      integer pnc_save_pdb

c ... Local variables:
      character*256 pnc_save_name
      character*8 stepstr
    

c      write(*,*) '**************************************'
c        write(*,*) "************************************************"
      
      write(*,*) "------------------------------------------------"
      write(*,*) 'UEDGE plasma + DEGAS2 neutral coupling'

c ... Allocate PNC_data group
      call gchange("PNC_data",0) 

      do pnc_step=pnc_step+1,pnc_maxstep
        pnc_time = pnc_time + dtplasma
        write(*,*) '**************************************'
        write(*,*) 'Step: ',pnc_step
        write(*,*) 'Time: ',pnc_time

        #Neutral step
        dtreal=dtneut
		call store_neutrals
        call run_neutrals
		call update_neutrals

        #Plasma step
		dtreal=dtplasma
        call store_plasma
c        write(*,*) "Solving UEDGE plasma model without neutral gas"
c        write(*,*) "dtreal=",dtplasma	
c        call exmain  			#not uedge_plasma to be consistent with other pnc routines
        call uedge_plasma
        call update_plasma
        call pnc_diagnostics    

        #Save PNC data
        if(mod(pnc_step,pnc_nsave)==0) then
          write(stepstr,'(i0)') pnc_step
          pnc_save_name=trim(pnc_savefile)//trim(stepstr)//'.pdb'
          if(.not.(pnc_save_pdb(pnc_save_name)==0)) 
     .      write(*,*) 'uedge_degas2: error writing '//pnc_save_name
        endif

      enddo
      pnc_step=pnc_step-1 #return value to pnc_maxstep

      end

c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to get fnrm
c*******************************************************************
      real function get_fnrm(dtreal_try)

      implicit none

c ... Common blocks:
      Use(Dim)			  # neq 
      Use(Math_problem_size) # neqmx
      Use(Lsode)    	  # yl, yldot
      Use(Time_dep_nwt)   # dtreal
      Use(Ynorm)		  # sfscal

c ... Local variables:
      real dtreal_try, dtreal_save

      dtreal_save=dtreal
      dtreal=dtreal_try
      dtuse=dtreal_try
c      call convert() #convert ni,up, ti... to yl
      call pandf(-1, -1, neq, 1., yl, yldot)
      get_fnrm=sqrt(sum((yldot(1:neq)*sfscal(1:neq))**2))
      dtreal=dtreal_save
      
      end function get_fnrm

c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to get fnrm from pandf (infinite time step)
c   WARNING! there is often a rescaling b/w pandf and pandf1 results
c*******************************************************************
      real function get_fnrm_pandf()

      implicit none

c ... Common blocks:
      Use(Dim)			  # neq 
      Use(Math_problem_size) # neqmx
      Use(Lsode)    	  # yl, yldot
      Use(Time_dep_nwt)   # dtreal
      Use(Ynorm)		  # sfscal

c ... Local variables:
      real dtreal_save

      dtreal_save=dtreal
      dtreal=1.0d20
c      call convert() #convert ni,up, ti... to yl
      call pandf(-1, -1, neq, 1., yl, yldot)
      get_fnrm_pandf=sqrt(sum((yldot(1:neq)*sfscal(1:neq))**2))
      dtreal=dtreal_save
      
      end function get_fnrm_pandf


c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to save neutrals during coupled plasma + neutrals run
c*******************************************************************
      subroutine pnc_diagnostics()

      implicit none

c ... Common blocks:
      Use(Dim) 		# nx, ny, nisp
      Use(Ext_neutrals) #extneutopt
      Use(MCN_dim)  # needed for PNC
      Use(PNC_params) # pnc_step, ...
      Use(PNC_data) # ni_pnc, ..., sni_pnc
      Use(Phyvar)   # ev


c ... Local variables:
      integer ioerr
      real fnrm_new
      character*128 cmd
      logical :: isechocmdon =.false.  #uses "echo cmd" for debugging purposes

c ... External functions:
      real get_fnrm,get_fnrm_pandf 
c     basis subroutine dobalance must be predefined within basis parser

      fnrm_new=get_fnrm(1.0d20)

      if (pnc_step == 1) then	# Write header information into file
        open(unit=pnc_fp,file=trim(pnc_histfile),action="write", 
     .        status="replace",iostat=ioerr)
        if(ioerr .ne. 0) then 
          write(*,*) "Trouble opening ",trim(pnc_histfile)," iostat =",ioerr
        endif
        write(pnc_fp,*) "***************************************"
        write(pnc_fp,*) "TIME HISTORY OF RESIDUALS AND SUMMARIES"
        write(pnc_fp,*) "***************************************"
        write(pnc_fp,*) "relax_p  = ", relax_p
        write(pnc_fp,*) "relax_g  = ", relax_g
      endif

      open(unit=pnc_fp,file=trim(pnc_histfile),action="write",status="old",  
     .      position="append",iostat=ioerr)
      if(ioerr .ne. 0) then 
        write(*,*) "Trouble writing ",trim(pnc_histfile)," iostat =",ioerr
      endif

      write(pnc_fp,*) "***************************************"
      write(pnc_fp,*) 'step = ', pnc_step
      write(pnc_fp,*) 'time = ', pnc_time
      write(pnc_fp,*) "***************************************"
      write(pnc_fp,*) 'fnrm      = ', fnrm_new
      write(pnc_fp,*) 'res_ni    = ', res_ni
      write(pnc_fp,*) 'res_up    = ', res_up
      write(pnc_fp,*) 'res_ti    = ', res_ti
      write(pnc_fp,*) 'res_te    = ', res_te
      write(pnc_fp,*) 'res_phi   = ', res_phi

      if (extneutopt==0) then 
        write(pnc_fp,*) 'res_ng    = ', res_ng
        write(pnc_fp,*) 'res_upg   = ', res_upg
      else
        if (get_neutral_sources) then
          write(pnc_fp,*) 'res_sni   = ', res_sni
          write(pnc_fp,*) 'res_smor  = ', res_smor
          write(pnc_fp,*) 'res_smophi= ', res_smophi
          write(pnc_fp,*) 'res_smoz  = ', res_smoz
          write(pnc_fp,*) 'res_sei   = ', res_sei
          write(pnc_fp,*) 'res_see   = ', res_see
        endif
        if (get_neutral_moments) then
          write(pnc_fp,*) 'res_ng   = ', res_ng
          write(pnc_fp,*) 'res_tg   = ', res_tg
          write(pnc_fp,*) 'res_sng   = ', res_sng
          write(pnc_fp,*) 'res_seg   = ', res_seg
        endif
      endif

      close(unit=pnc_fp,iostat=ioerr)
      if(ioerr .ne. 0) then 
        write(*,*) "Trouble closing ",trim(pnc_histfile)," iostat =",ioerr
      endif

      if(pnc_dobalance) then
c       call parsestr('rundiagnostics()')
        call parsestr('dobalance(pnc_balancefile)')
c        call sleep(0.005)  #Need time lag before touching files
c       cmd="touch "//trim(pnc_histfile)
c       call system(trim(cmd))        
        cmd='cat '//trim(pnc_balancefile)//' >> '//trim(pnc_histfile)        
        if (isechocmdon) cmd="echo "//trim(cmd)
        if (ext_verbose) print *,trim(cmd)
        call system(trim(cmd)) 
      endif

      end

c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to save neutrals during coupled plasma + neutrals run
c*******************************************************************
      subroutine store_neutrals()

      implicit none

c ... Common blocks:
      Use(Dim) 			# nx, ny, nisp
      Use(Ext_neutrals) # extneutopt, get_neutral_moments
      Use(Interp)		# nis, ups
      Use(MCN_dim)
      Use(MCN_sources) 	# sni, ..., see, ng_mc,... tg_mc
      Use(PNC_params)	# pnc_step
      Use(PNC_data)		# sni_pnc, ...

      if (extneutopt==0) then 
        ni_pnc=nis
        up_pnc=ups
      else 
        if (get_neutral_sources) then
          sni_pnc   =sni
          smor_pnc  =smor
          smophi_pnc=smophi
          smoz_pnc  =smoz
          sei_pnc   =sei
          see_pnc   =see
        endif
        if (get_neutral_moments) then 
          ng_pnc  = ng_ue
          upg_pnc = upg_ue
          tg_pnc  = tg_ue
        endif
      endif

      end

c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to update neutrals during coupled plasma + neutrals run
c*******************************************************************
      subroutine update_neutrals()

      implicit none

c ... Common blocks:
      Use(Dim) 			# nx, ny, nisp
      Use(Compla)		# ni, up
      Use(Interp)		# nis, ups
      Use(Ext_neutrals) # extneutopt
      Use(MCN_dim)
      Use(MCN_sources) 	# sni, ..., see
      Use(PNC_params)	# pnc_step, ...
      Use(PNC_data)     # sni_pnc, ...
      Use(Phyvar)   	# ev
      Use(UEpar)		# iigsp

c ... Local Variables
      integer ifld

      if(extneutopt==0) then 
        ni(:,:,iigsp)=relax_p*ni(:,:,iigsp)+(1-relax_p)*ni_pnc(:,:,iigsp)
        up(:,:,iigsp)=relax_p*up(:,:,iigsp)+(1-relax_p)*up_pnc(:,:,iigsp)

        nis(:,:,iigsp)=ni(:,:,iigsp)
        ups(:,:,iigsp)=up(:,:,iigsp)

        res_ng    = sqrt(sum((ni(:,:,iigsp)-ni_pnc(:,:,iigsp))**2/real(nx*ny)))
        res_upg   = sqrt(sum((up(:,:,iigsp)-up_pnc(:,:,iigsp))**2/real(nx*ny)))
        
        del_ng    = maxval(abs(ni(:,:,iigsp)-ni_pnc(:,:,iigsp)))
        del_upg   = maxval(abs(up(:,:,iigsp)-up_pnc(:,:,iigsp)))

        if (pnc_print_norm==1) then #normalize to maximum value
          res_ng     = res_ng/maxval(abs(ni(:,:,iigsp)))
          res_upg     = res_upg/maxval(abs(up(:,:,iigsp)))

          del_ng     = del_ng/maxval(abs(ni(:,:,iigsp)))
          del_upg     = del_upg/maxval(abs(up(:,:,iigsp)))
        else if (pnc_print_norm==2) then #normalize to previous maximum value
          res_ng     = res_ng/maxval(abs(ni_pnc))
          res_upg     = res_upg/maxval(abs(up_pnc))

          del_ng     = del_ng/maxval(abs(ni_pnc))
          del_upg     = del_upg/maxval(abs(up_pnc))
        end if
        if(pnc_verbose) then
          write(*,*) 'Variable           Max Res                        Std Dev Res'
          write(*,*) 'ng        ', del_ng,  res_ng
          write(*,*) 'upg       ', del_upg, res_upg
        endif
      endif
      if (get_neutral_sources) then
        sni   =relax_g*sni   +(1-relax_g)*sni_pnc
        smor  =relax_g*smor  +(1-relax_g)*smor_pnc
        smophi=relax_g*smophi+(1-relax_g)*smophi_pnc
        smoz  =relax_g*smoz  +(1-relax_g)*smoz_pnc
        sei   =relax_g*sei   +(1-relax_g)*sei_pnc
        see   =relax_g*see   +(1-relax_g)*see_pnc

        res_sni    = sqrt(sum((sni-sni_pnc)**2)/real(nx*ny*nfl*nstra))
        res_smor   = sqrt(sum((smor-smor_pnc)**2)/real(nx*ny*nfl*nstra))
        res_smophi = sqrt(sum((smophi-smophi_pnc)**2)/real(nx*ny*nfl*nstra))
        res_smoz   = sqrt(sum((smoz-smoz_pnc)**2)/real(nx*ny*nfl*nstra))
        res_sei    = sqrt(sum((sei-sei_pnc)**2)/real(nx*ny*nstra))
        res_see    = sqrt(sum((see-see_pnc)**2)/real(nx*ny*nstra))

        del_sni    = maxval(abs(sni-sni_pnc))
        del_smor   = maxval(abs(smor-smor_pnc))
        del_smophi = maxval(abs(smophi-smophi_pnc))
        del_smoz   = maxval(abs(smoz-smoz_pnc))
        del_sei    = maxval(abs(sei-sei_pnc))
        del_see    = maxval(abs(see-see_pnc))

        if (pnc_print_norm==1 .or. (pnc_step<=1 .and. pnc_print_norm==2) ) then
          res_sni    = res_sni/maxval(abs(sni))
          res_smor   = res_smor/maxval(abs(smor))
          res_smophi = res_smophi/maxval(abs(smophi))
          res_smoz   = res_smoz/maxval(abs(smoz))
          res_sei    = res_sei/maxval(abs(sei))
          res_see    = res_see/maxval(abs(see))

          del_sni    = del_sni/maxval(abs(sni))
          del_smor   = del_smor/maxval(abs(smor))
          del_smophi = del_smophi/maxval(abs(smophi))
          del_smoz   = del_smoz/maxval(abs(smoz))
          del_sei    = del_sei/maxval(abs(sei))
          del_see    = del_see/maxval(abs(see))
        else if (pnc_print_norm==2) then #normalize to previous maximum value
          res_sni    = res_sni/maxval(abs(sni_pnc))
          res_smor   = res_smor/maxval(abs(smor_pnc))
          res_smophi = res_smophi/maxval(abs(smophi_pnc))
          res_smoz   = res_smoz/maxval(abs(smoz_pnc))
          res_sei    = res_sei/maxval(abs(sei_pnc))
          res_see    = res_see/maxval(abs(see_pnc))

          del_sni    = del_sni/maxval(abs(sni_pnc))
          del_smor   = del_smor/maxval(abs(smor_pnc))
          del_smophi = del_smophi/maxval(abs(smophi_pnc))
          del_smoz   = del_smoz/maxval(abs(smoz_pnc))
          del_sei    = del_sei/maxval(abs(sei_pnc))
          del_see    = del_see/maxval(abs(see_pnc))
        end if

         if (pnc_verbose .and. pnc_step .gt. 1) then
          write(*,*) 'Variable           Max Res                         Std Dev Res'
          write(*,*) 'sni     ', del_sni,    res_sni
          write(*,*) 'smor    ', del_smor,   res_smor
          write(*,*) 'smophi  ', del_smophi, res_smophi
          write(*,*) 'smoz    ', del_smoz,   res_smoz
          write(*,*) 'sei     ', del_sei,    res_sei
          write(*,*) 'see     ', del_see,    res_see
        endif
      endif
      if (get_neutral_moments) then
        if(pnc_ngs_mc) then
          ngs = relax_g*ng_ue + (1-relax_g)*ng_pnc
        endif
        if(pnc_upgs_mc) then
            ups(:,:,iigsp) = relax_g*upg_ue(:,:,1) + (1-relax_g)*upg_pnc(:,:,1)
        endif
        if(pnc_tgs_mc) then
c         tgs = relax_g*tg_ue(:,:,1) + (1-relax_g)*tg_pnc(:,:,1)
        endif

        res_ng = sqrt(sum((ng_ue-ng_pnc)**2)/real(nx*ny*nfl*nstra))
        res_tg = sqrt(sum((tg_ue-tg_pnc)**2)/real(nx*ny*nfl*nstra))

        del_ng = maxval(abs(ng_ue-ng_pnc))
        del_tg = maxval(abs(tg_ue-tg_pnc))

        if (pnc_print_norm==1 .or. (pnc_step<=1 .and. pnc_print_norm==2) ) then       
          res_ng = res_ng/maxval(abs(ng_ue))
          res_tg = res_tg/maxval(abs(tg_ue))
          del_ng = del_ng/maxval(abs(ng_ue))
          del_tg = del_tg/maxval(abs(tg_ue))

        else if (pnc_print_norm==2) then #normalize to previous maximum value
          res_ng = res_ng/maxval(abs(ng_pnc))
          res_tg = res_tg/maxval(abs(tg_pnc))
          del_ng = del_ng/maxval(abs(ng_pnc))
          del_tg = del_tg/maxval(abs(tg_pnc))

        endif

         if (pnc_verbose .and. pnc_step .gt. 1) then
          if (.not. get_neutral_sources) then 
            write(*,*) 'Variable           Max Res                         Std Dev Res'
          endif
          write(*,*) 'ng     ', del_ng, res_ng
          write(*,*) 'tg     ', del_tg, res_tg
        endif
      endif
      if(pnc_verbose) then 
c        write(*,*) "************************************************"
        write(*,*) "------------------------------------------------"
      endif
      end

c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to store plasma data during coupled plasma + neutrals run
c*******************************************************************
      subroutine store_plasma()

      implicit none

c ... Common blocks:
      Use(Dim) 		# nx, ny, nisp
      Use(Compla) 	# ni, up, ...
      Use(Ext_neutrals)		# get_neutral_moments
      Use(MCN_dim)  # needed for PNC_data
      Use(MCN_sources)  # sng_ue, ...
      Use(PNC_params) # pnc_step, ...
      Use(PNC_data) # ni_pnc, ...

      ni_pnc=ni
      up_pnc=up
      ti_pnc=ti
      te_pnc=te
      phi_pnc=phi

c ... Store neutral sources calculated by pandf
      if (get_neutral_moments) then 
 	      sng_pnc = sng_ue
		  seg_pnc = seg_ue
      endif

      end

c-----------------------------------------------------------------------

c*******************************************************************
c Subroutine to update plasma during coupled plasma + neutrals run
c*******************************************************************
      subroutine update_plasma()

      implicit none

c ... Common blocks:
      Use(Dim) 		# nx, ny, nisp
      Use(Compla) 	# ni, up, ...
      Use(Ext_neutrals)		# get_neutral_moments
      Use(Interp) 	# nis, ups, ...
      Use(MCN_dim)  # needed for PNC
      Use(MCN_sources)  # needed for PNC
      Use(PNC_params) # pnc_step, ...
      Use(PNC_data) # ni_pnc, ..., sni_pnc
      Use(Phyvar)   # ev

      #Perform update
      ni =relax_p*ni +(1-relax_p)*ni_pnc
      up =relax_p*up +(1-relax_p)*up_pnc
      ti =relax_p*ti +(1-relax_p)*ti_pnc
      te =relax_p*te +(1-relax_p)*te_pnc
      phi=relax_p*phi+(1-relax_p)*phi_pnc

      nis  = ni
      ups  = up 
      tis  = ti
      tes  = te
      phis = phi

      #Compute residuals
      res_ni   = sqrt(sum((ni-ni_pnc)**2)/real(nx*ny*nisp))
      res_up   = sqrt(sum((up-up_pnc)**2)/real(nx*ny*nisp))
      res_ti   = sqrt(sum((ti-ti_pnc)**2)/real(nx*ny))
      res_te   = sqrt(sum((te-te_pnc)**2)/real(nx*ny))
      res_phi  = sqrt(sum((phi-phi_pnc)**2)/real(nx*ny))
       
      del_ni   = maxval(abs(ni-ni_pnc))
      del_up   = maxval(abs(up-up_pnc))
      del_ti   = maxval(abs(ti-ti_pnc))
      del_te   = maxval(abs(te-te_pnc))
      del_phi  = maxval(abs(phi-phi_pnc))

c ... Calculate residuals for neutral sources calculated by pandf
      if (get_neutral_moments) then 
        res_sng = sqrt(sum((sng_ue-sng_pnc)**2)/real(nx*ny*nfl*nstra))
        res_seg = sqrt(sum((seg_ue-seg_pnc)**2)/real(nx*ny*nfl*nstra))
        del_sng = maxval(abs(sng_ue-sng_pnc))
        del_seg = maxval(abs(seg_ue-seg_pnc))

c       For debugging calculation
c		res_sng = maxval(abs(sng_ue))
c		res_seg = maxval(abs(seg_ue))
c		del_sng = maxval(abs(fngx_ue))
c	    del_seg = maxval(abs(fegx_ue))
      endif

      if (pnc_print_norm==1 .or. (pnc_step<=1 .and. pnc_print_norm==2)) then #normalize to maximum value
        res_ni     = res_ni/maxval(abs(ni))
        res_up     = res_up/maxval(abs(up))
        res_ti     = res_ti/maxval(abs(ti))
        res_te     = res_te/maxval(abs(te))
        res_phi    = res_phi/maxval(abs(phi))

        del_ni     = del_ni/maxval(abs(ni))
        del_up     = del_up/maxval(abs(up))
        del_ti     = del_ti/maxval(abs(ti))
        del_te     = del_te/maxval(abs(te))
        del_phi    = del_phi/maxval(abs(phi))

        if (get_neutral_moments) then
         res_sng = res_sng/maxval(abs(sng_ue))
         res_seg = res_seg/maxval(abs(seg_ue))
         del_sng = del_sng/maxval(abs(sng_ue))
         del_seg = del_seg/maxval(abs(seg_ue))
        endif  
      else if (pnc_print_norm==2) then #normalize to previous maximum value
        res_ni     = res_ni/maxval(abs(ni_pnc))
        res_up     = res_up/maxval(abs(up_pnc))
        res_ti     = res_ti/maxval(abs(ti_pnc))
        res_te     = res_te/maxval(abs(te_pnc))
        res_phi    = res_phi/maxval(abs(phi_pnc))

        del_ni     = del_ni/maxval(abs(ni_pnc))
        del_up     = del_up/maxval(abs(up_pnc))
        del_ti     = del_ti/maxval(abs(ti_pnc))
        del_te     = del_te/maxval(abs(te_pnc))
        del_phi    = del_phi/maxval(abs(phi_pnc))

        if (get_neutral_moments) then
          res_sng = res_sng/maxval(abs(sng_pnc))
          res_seg = res_seg/maxval(abs(seg_pnc))
          del_sng = del_sng/maxval(abs(sng_pnc))
          del_seg = del_seg/maxval(abs(seg_pnc))
        endif  
      end if
      if(pnc_verbose) then
        write(*,*) 'Variable           Max Res                         Std Dev Res'
          if (get_neutral_moments) then
            write(*,*) 'sng     ', del_sng, res_sng
            write(*,*) 'seg     ', del_seg, res_seg
          endif  

        write(*,*) 'ni      ', del_ni,  res_ni
        write(*,*) 'up      ', del_up,  res_up
        write(*,*) 'ti      ', del_ti,  res_ti
        write(*,*) 'te      ', del_te,  res_te
        write(*,*) 'phi     ', del_phi, res_phi
      endif
c      write(*,*) '************************************************'
       write(*,*) '------------------------------------------------'
      end

c-----------------------------------------------------------------------

c******************************************************
c Function for saving UEDGE data to PDB file
c******************************************************
      integer function uedge_save_pdb(filename)

      implicit none

c ... Common blocks:
      Use(Dim)
      Use(Interp) 		#nis
      Use(PNC_params) 	#pnc_savefile

c ... Local variables:
      character*256 filename,varstr
      character*8 nxstr,nystr,nispstr,ngspstr
      integer fileid, flen

c ... External functions: 
      integer pfopen, pfwrta, pfclos, pfgerr # PDB API

      flen=len(trim(filename))
      fileid = pfopen(flen, trim(filename), 'w')
      if (fileid .eq. 0)
     .  write(*,*) pfgerr
cc     .   write(*,*) 'PDB error opening ', filename #call errproc
c     write(*,*) 'PDB open ',trim(filename)

      write(nxstr,'(i0)') nx+1
      write(nystr,'(i0)') ny+1
      write(nispstr,'(i0)') nisp
      write(ngspstr,'(i0)') ngsp

      varstr='ngs(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(ngspstr)//')'
c      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', ngs) == 0)
     .  write(*,*) 'PDB error writing ngs', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='nis(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nispstr)//')'
c      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', nis) == 0)
     .  write(*,*) 'PDB error writing nis', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='ups(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nispstr)//')'
c      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', ups) == 0)
     .  write(*,*) 'PDB error writing ups', pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='tis(0:'//trim(nxstr)//',0:'//trim(nystr)//')'
c      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', tis) == 0)
     .  write(*,*) 'PDB error writing tis', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='tes(0:'//trim(nxstr)//',0:'//trim(nystr)//')'
c      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', tes) == 0)
     .  write(*,*) 'PDB error writing tes', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='phis(0:'//trim(nxstr)//',0:'//trim(nystr)//')'
c      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', phis) == 0)
     .  write(*,*) 'PDB error writing phis', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      if (pfclos(fileid)==0) write(*,*) 'Error closing file', pfgerr
c     .   write(*,*) 'PDB error closing ', filename #call errproc
cc      write(*,*) 'PDB close ',trim(filename)

      uedge_save_pdb=fileid
      end function uedge_save_pdb

c-----------------------------------------------------------------------

c******************************************************
c Function for saving UEDGE variables to default PDB filename
c******************************************************
      integer function uedge_save()

      implicit none

c ... Common blocks:
      Use(Dim) 		# for Interp
      Use(Interp) 	# uedge_savefile

c ... External functions:
      integer uedge_save_pdb
     
      uedge_save=uedge_save_pdb(trim(uedge_savefile))

      end function uedge_save

c-----------------------------------------------------------------------

c******************************************************
c Function for reading UEDGE data from PDB file
c WARNING: Not working yet!!!
c******************************************************
      integer function uedge_read_pdb(filename)

      implicit none

c ... Common blocks:
      Use(Dim)
      Use(Interp) 		#nis
      Use(PNC_params) 	#pnc_savefile

c ... Local variables:
      character*256 filename,varstr
      character*8 nxstr,nystr,nispstr,ngspstr
      integer fileid, flen

c ... External functions: 
      integer pfopen, pfread, pfclos, pfgerr # PDB API

      flen=len(trim(filename))
      fileid = pfopen(flen, trim(filename), 'r')
      if (fileid .eq. 0)
     .  write(*,*) 'PDB error opening file', pfgerr
cc     .   write(*,*) 'PDB error opening ', filename #call errproc
cc      write(*,*) 'PDB open ',trim(filename)

      varstr='ngs(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(ngspstr)//')'
c     write(*,*) trim(varstr), len(trim(varstr))
      if (pfread(fileid, len(trim(varstr)), trim(varstr), ngs)==0)
     .  write(*,*) 'PDB error reading ngs', pfgerr
cc     .   write(*,*) 'PDB error reading ', filename #call errproc

      write(nxstr,'(i0)') nx+1
      write(nystr,'(i0)') ny+1
      write(nispstr,'(i0)') nisp
      write(ngspstr,'(i0)') ngsp

      varstr='nis(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nispstr)//')'
c     write(*,*) trim(varstr), len(trim(varstr))
      if (pfread(fileid, len(trim(varstr)), trim(varstr), nis)==0)
     .  write(*,*) 'PDB error reading nis', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='ups(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nispstr)//')'
c     write(*,*) trim(varstr), len(trim(varstr))
      if (pfread(fileid, len(trim(varstr)), trim(varstr), ups)==0)
     .  write(*,*) 'PDB error reading ups', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='tis(0:'//trim(nxstr)//',0:'//trim(nystr)//')'
c     write(*,*) trim(varstr), len(trim(varstr))
      if (pfread(fileid, len(trim(varstr)), trim(varstr), tis)==0)
     .  write(*,*) 'PDB error reading tis', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='tes(0:'//trim(nxstr)//',0:'//trim(nystr)//')'
c     write(*,*) trim(varstr), len(trim(varstr))
      if (pfread(fileid, len(trim(varstr)), trim(varstr),tes)==0)
     .  write(*,*) 'PDB error reading tes', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='phis(0:'//trim(nxstr)//',0:'//trim(nystr)//')'
c     write(*,*) trim(varstr), len(trim(varstr))
      if (pfread(fileid, len(trim(varstr)), trim(varstr), phis)==0)
     .  write(*,*) 'PDB error reading phis', pfgerr
cc     .   write(*,*) 'PDB error writing ', filename #call errproc

      if (pfclos(fileid)==0) write(*,*) 'Error closing file', pfgerr
c     .   write(*,*) 'PDB error closing ', filename #call errproc
c     write(*,*) 'PDB close ',trim(filename)

      uedge_read_pdb=fileid
      end function uedge_read_pdb


c-----------------------------------------------------------------------

c******************************************************
c Function for reading UEDGE variables from default PDB filename
c******************************************************
      integer function uedge_read()

      implicit none

c ... Common blocks:
      Use(Dim) 		# for Interp
      Use(Interp) 	# uedge_savefile

c ... External functions:
      integer uedge_read_pdb
     
      uedge_read=uedge_read_pdb(trim(uedge_savefile))

      end function uedge_read

c-----------------------------------------------------------------------

c******************************************************
c Function for saving neutral source data to PDB file
c******************************************************
      integer function mcnsor_save_pdb(filename)

      implicit none

c ... Common blocks: 
      Use(Dim)
      Use(MCN_dim)
      Use(MCN_sources)  #sni

c ... Local variables:
      character*256 filename,varstr
      character*8 nxstr,nystr,nflstr,nstrastr
      integer fileid, flen

c ... External functions: 
      integer pfopen, pfwrta, pfclos, pfgerr # PDB API

      flen=len(trim(filename))

      fileid = pfopen(flen, trim(filename), 'w')
      if (fileid .eq. 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error opening ', filename #call errproc
      write(*,*) 'PDB open ',trim(filename)

      write(nxstr,'(i0)') nx+1
      write(nystr,'(i0)') ny+1
      write(nflstr,'(i0)') nfl
      write(nstrastr,'(i0)') nstra

      varstr='sni(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', sni) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='smor(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', smor) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='smophi(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', smophi) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='smoz(0:'//trim(nxstr)//',0:'//trim(nystr)//
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', smoz) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='sei(0:'//trim(nxstr)//',0:'//trim(nystr)//
     .  ',1:'//trim(nstrastr)//')pnc_save_pdb'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', sei) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='see(0:'//trim(nxstr)//',0:'//trim(nystr)//
     .  ',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', see) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      if (pfclos(fileid)==0) write(*,*) pfgerr
c     .   write(*,*) 'PDB error closing ', filename #call errproc
      write(*,*) 'PDB close ',trim(filename)

      mcnsor_save_pdb=fileid
      end function mcnsor_save_pdb

c-----------------------------------------------------------------------

c******************************************************
c Function for appending neutral source data to PDB file
c******************************************************
      integer function mcnsor_append_pdb(filename)

      implicit none

c ... Common blocks: 
      Use(Dim)
      Use(MCN_dim)
      Use(MCN_sources)  #sni
      Use(PNC_params) 	#pnc_savefile

c ... Local variables:
      character*256 filename,varstr
      character*8 nxstr,nystr,nflstr,nstrastr
      integer fileid, flen

c ... External functions: 
      integer pfopen, pfwrta, pfclos, pfgerr # PDB API

      flen=len(trim(filename))
      fileid = pfopen(flen, trim(filename), 'a')
      if (fileid .eq. 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error opening ', filename #call errproc
      write(*,*) 'PDB open ',trim(filename)

      write(nxstr,'(i0)') nx+1
      write(nystr,'(i0)') ny+1
      write(nflstr,'(i0)') nfl
      write(nstrastr,'(i0)') nstra

      varstr='sni(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', sni) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='smor(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', smor) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='smophi(0:'//trim(nxstr)//',0:'//trim(nystr)// 
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', smophi) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='smoz(0:'//trim(nxstr)//',0:'//trim(nystr)//
     .  ',1:'//trim(nflstr)//',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', smoz) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='sei(0:'//trim(nxstr)//',0:'//trim(nystr)//
     .  ',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', sei) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      varstr='see(0:'//trim(nxstr)//',0:'//trim(nystr)//
     .  ',1:'//trim(nstrastr)//')'
      write(*,*) trim(varstr), len(trim(varstr))
      if (pfwrta(fileid, len(trim(varstr)), trim(varstr), 
     .   len('double'), 'double', see) == 0)
     .  write(*,*) pfgerr
c     .   write(*,*) 'PDB error writing ', filename #call errproc

      if (pfclos(fileid)==0) write(*,*) pfgerr
c     .   write(*,*) 'PDB error closing ', filename #call errproc
      write(*,*) 'PDB close ',trim(filename)

      mcnsor_append_pdb=fileid
      end function mcnsor_append_pdb

c-----------------------------------------------------------------------

c******************************************************
c Function for saving UEDGE variables and neutral sources to default PDB filename
c******************************************************
      integer function pnc_save_pdb(filename)

      implicit none

c ... Common blocks:
      Use(Ext_neutrals) # extneutopt

c ... External functions:
      integer uedge_save_pdb, mcnsor_append_pdb

c ... Local variables:
      character*256 filename
     
      pnc_save_pdb=uedge_save_pdb(filename)

      if (extneutopt .gt. 0) then
        pnc_save_pdb=mcnsor_append_pdb(filename)
      endif

      end function pnc_save_pdb

c----------------------------------------------------------------------c

      subroutine readmcndens(fname)
      implicit none
      character*(*) fname

      Use(Dim)		# nx,ny
      Use(MCN_dim)
      Use(MCN_sources)
      Use(Ext_neutrals) 		 # ext_verbose
      integer nunit,ix,iy,ifl,idum
      real rdum

c     Read data from DEGAS2 code output data file 'density.out':

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

c     Plasma density and pressure:
         do ifl=1,1 #ngsp
           do ix=1,nx
             do iy=1,ny
               read (nunit,*) idum, idum,
     .              ng_mc(ix,iy,ifl), ng_mc_rsd(ix,iy,ifl),
     .              pg_mc(ix,iy,ifl), pg_mc_rsd(ix,iy,ifl),
     .              pxz_mc(ix,iy,ifl), pxz_mc_rsd(ix,iy,ifl)
             enddo
           enddo
         enddo

      close (nunit)
      if (ext_verbose)
     . call remark(" *** neutral density read from DEGAS2 file "
     .                                                  //fname//" ***")

      end

c----------------------------------------------------------------------c

      subroutine readmcnoutput(fname,hskip,var,rsd)
      implicit none

      Use(Dim)		# nx,ny
      Use(MCN_dim)
      Use(MCN_sources)
      Use(Ext_neutrals) 		 # ext_verbose

c ... Input Variables
      character*(*) fname
      integer hskip

c ... Output Variables
      real var(0:nx+1,0:ny+1,1:ngsp), rsd(0:nx+1,0:ny+1,1:ngsp)

c ... Local Variables
      integer nunit,ix,iy,ifl,idum,tskip,il
      real rdum 

c     Read data from DEGAS2 code outputbrowser data files such as 'neutral_density.dat':
      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      var=0
      rsd=0

      tskip=5
c     Skip header
      do il=1,hskip
        read(nunit,*)
      enddo
c     Read variable and rel. std. dev.
      do ifl=1,ngsp
c       write(*,*) "ifl =",ifl
        do ix=1,nx
          do iy=1,ny
            read (nunit,*) idum, idum,
     .        var(ix,iy,ifl), rsd(ix,iy,ifl)
c           if (abs(var(ix,iy,ifl)) .gt. 0)
c    .        write(*,*) ix, iy, ifl, var(ix,iy,ifl), rsd(ix,iy,ifl)
           enddo
        enddo
c       Skip trailing data
        do il=1,tskip
          read(nunit,*)
        enddo
      enddo


      close (nunit)
      if (ext_verbose)
     . call remark(" *** output read from DEGAS2 file "
     .                                                  //fname//" ***")

      end

c----------------------------------------------------------------------c

      subroutine readmcnmoments(dname)
      implicit none

      Use(Dim)				# nx,ny
      Use(MCN_dim)
      Use(MCN_sources)
      Use(Ext_neutrals) 	# ext_verbose, moment filenames

c ... Input Variables
      character*(*) dname

c ... Local Variables      
      logical verbose
      character*512 fname
      integer sskip, vskip, idim
      real var1(0:nx+1,0:ny+1,nfl), rsd1(0:nx+1,0:ny+1,nfl)

c ... External Functions
      external readmcnoutput, mcndivide

c ... Read Data
      sskip=4
      vskip=5

      verbose=ext_verbose
      ext_verbose=.false.

c     Read data from DEGAS2 code outputbrowser data files such as 'neutral_density.dat':
      fname=dname//"/"//neut_ng_file
      call readmcnoutput(trim(fname),sskip,var1,rsd1)
      ng_mc=var1
      ng_mc_rsd=rsd1
      fname=dname//"/"//neut_pg_file
      call readmcnoutput(trim(fname),sskip,var1,rsd1)
      pg_mc=var1
      pg_mc_rsd=rsd1

      fname=dname//"/"//neut_jng1_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      jng_mc(:,:,:,1)=var1
      jng_mc_rsd(:,:,:,1)=rsd1
      fname=dname//"/"//neut_jng2_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      jng_mc(:,:,:,2)=var1
      jng_mc_rsd(:,:,:,2)=rsd1
      fname=dname//"/"//neut_jng3_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      jng_mc(:,:,:,3)=var1
      jng_mc_rsd(:,:,:,3)=rsd1

      fname=dname//"/"//neut_pg11_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      stressg_mc(:,:,:,1,1)=var1
      stressg_mc_rsd(:,:,:,1,1)=rsd1
      fname=dname//"/"//neut_pg22_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      stressg_mc(:,:,:,2,2)=var1
      stressg_mc_rsd(:,:,:,2,2)=rsd1
      fname=dname//"/"//neut_pg33_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      stressg_mc(:,:,:,3,3)=var1
      stressg_mc_rsd(:,:,:,3,3)=rsd1

      fname=dname//"/"//neut_pg23_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      stressg_mc(:,:,:,2,3)=var1
      stressg_mc_rsd(:,:,:,2,3)=rsd1
      stressg_mc(:,:,:,3,2)=var1
      stressg_mc_rsd(:,:,:,3,2)=rsd1
      fname=dname//"/"//neut_pg31_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      stressg_mc(:,:,:,3,1)=var1
      stressg_mc_rsd(:,:,:,3,1)=rsd1
      stressg_mc(:,:,:,1,3)=var1
      stressg_mc_rsd(:,:,:,1,3)=rsd1
      fname=dname//"/"//neut_pg12_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      stressg_mc(:,:,:,1,2)=var1
      stressg_mc_rsd(:,:,:,1,2)=rsd1
      stressg_mc(:,:,:,2,1)=var1
      stressg_mc_rsd(:,:,:,1,2)=rsd1

      fname=dname//"/"//neut_jeg1_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      jeg_mc(:,:,:,1)=var1
      jeg_mc_rsd(:,:,:,1)=rsd1
      fname=dname//"/"//neut_jeg2_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      jeg_mc(:,:,:,2)=var1
      jeg_mc_rsd(:,:,:,2)=rsd1
      fname=dname//"/"//neut_jeg3_file
      call readmcnoutput(trim(fname),vskip,var1,rsd1)
      jeg_mc(:,:,:,3)=var1
      jeg_mc_rsd(:,:,:,3)=rsd1

      ext_verbose=verbose
       if (ext_verbose)
     . call remark(" *** neutral moments read from DEGAS2 directory "
     .                                                  //dname//" ***")     

      end

c----------------------------------------------------------------------c

      subroutine convertmcnmoments()
      implicit none

      Use(Dim)				# nx,ny
      Use(Comflo) 			# fngx,fngy,...
      Use(Comgeo)			# sx,sy
      Use(Compla)			# ng,tg
      Use(Interp)			# ngs,tis
      Use(MCN_dim)
      Use(MCN_sources)      # ng_mc,ng_mc_rsd,... 
      Use(Ext_neutrals) 	# ext_verbose
      Use(PNC_params)		# pnc moment filenames

c ... Local Variables
      integer idim,ifld

c ... Scalars
c     write(*,*) "convertmcnmoments: scalar"
      call mcnblend(ng_ue,ngs,ng_mc,ng_ue_rsd,ng_mc_rsd,mcalpha_ng)
      call mcnblend(pg_ue,ngs(:,:,1)*tis,pg_mc,pg_ue_rsd,pg_mc_rsd,mcalpha_pg)
      call mcuedivide(tg_ue,pg_ue,ng_ue,tg_ue_rsd,pg_ue_rsd,ng_ue_rsd)

c ... Vectors
c ... 1. Convert vectors
c ... 2. Multiply current density by area in order to get current
c ... 3. Blend kinetic fluxes with fluid fluxes
c ... 4. Calculate velocity

c     write(*,*) "convertmcnmoments: vector"
      call convertmcnvector(jng_mc,jng_ue,jng_mc_rsd,jng_ue_rsd)
      do ifld=1,nfl
        fngx_mc(:,:,ifld)=sx*jng_ue(:,:,ifld,1) 
        fngy_mc(:,:,ifld)=sy*jng_ue(:,:,ifld,2)
      enddo
      fngx_mc_rsd=jng_ue_rsd(:,:,:,1)
      fngy_mc_rsd=jng_ue_rsd(:,:,:,2)
      call mcnblend(fngx_ue,fngx,fngx_mc,fngx_ue_rsd,fngx_mc_rsd,mcalpha_fng)
      call mcnblend(fngy_ue,fngy,fngy_mc,fngy_ue_rsd,fngy_mc_rsd,mcalpha_fng)
      fngx_ue=fngx
      fngy_ue=fngy

      do idim=1,3
      do ifld=1,nfl
        call mcuedivide(vg_mc(:,:,ifld,idim),jng_mc(:,:,ifld,idim),ng_ue,
     .                  vg_mc_rsd(:,:,ifld,idim),jng_mc_rsd(:,:,ifld,idim),ng_ue_rsd)
        call mcuedivide(vg_ue(:,:,ifld,idim),jng_ue(:,:,ifld,idim),ng_ue,
     .                  vg_ue_rsd(:,:,ifld,idim),jng_ue_rsd(:,:,ifld,idim),ng_ue_rsd)
      enddo
      enddo
      upg_ue = vg_ue(:,:,:,3)
      upg_ue_rsd = vg_ue_rsd(:,:,:,3)
     
      call convertmcnvector(jeg_mc,jeg_ue,jeg_mc_rsd,jeg_ue_rsd)
      do ifld=1,nfl
        fegx(:,:,ifld)=0.0
        fegy(:,:,ifld)=0.0
        fegx_mc(:,:,ifld)=sx*jeg_ue(:,:,ifld,1) 
        fegy_mc(:,:,ifld)=sy*jeg_ue(:,:,ifld,2)
      enddo
      fegx_mc_rsd=jeg_ue_rsd(:,:,:,1)
      fegy_mc_rsd=jeg_ue_rsd(:,:,:,2)
      call mcnblend(fegx_ue,fegx,fegx_mc,fegx_ue_rsd,fegx_mc_rsd,mcalpha_feg)
      call mcnblend(fegy_ue,fegy,fegy_mc,fegy_ue_rsd,fegy_mc_rsd,mcalpha_feg)
 
c ... Tensors
c ... 1. convertmcntensor
c ... 2. find parallel component, perpendicular curl, etc.
c ... 3. Blend kinetic and fluid results

c     write(*,*) "convertmcnmoments: tensor"
      call convertmcntensor(stressg_mc,stressg_ue,stressg_mc_rsd,stressg_ue_rsd)
      do ifld=1,nfl
        fmgx(:,:,ifld)=fmix(:,:,2*ifld)
        fmgy(:,:,ifld)=fmiy(:,:,2*ifld)
        fmgx_mc(:,:,ifld)=sx*stressg_ue(:,:,ifld,1,3) 
        fmgy_mc(:,:,ifld)=sy*stressg_ue(:,:,ifld,2,3)
      end do
      fmgx_mc_rsd=stressg_ue_rsd(:,:,:,1,3)
      fmgy_mc_rsd=stressg_ue_rsd(:,:,:,2,3)
      call mcnblend(fmgx_ue,fmgx,fmgx_mc,fmgx_ue_rsd,fmgx_mc_rsd,mcalpha_fmg)
      call mcnblend(fmgy_ue,fmgy,fmgy_mc,fmgy_ue_rsd,fmgy_mc_rsd,mcalpha_fmg)

c ... Finally, find contribution to vorticity ... is this already handled via friction forces?

      end

c----------------------------------------------------------------------

      subroutine mcndivide(out,var,dens,out_rsd,var_rsd,dens_rsd)
      implicit none

      Use(Dim)				# nx,ny
      Use(MCN_dim)          # nfl

c ... Input Variables
      real var(0:nx+1,0:ny+1,1:ngsp),dens(0:nx+1,0:ny+1,1:nfl)
      real var_rsd(0:nx+1,0:ny+1,1:ngsp),dens_rsd(0:nx+1,0:ny+1,1:nfl)

c ... Output Variables
      real out(0:nx+1,0:ny+1,1:nfl),out_rsd(0:nx+1,0:ny+1,1:nfl)     

      out=0.;out_rsd=1.;
      where(dens .gt. 0.0) out=var/dens
      where(dens .gt. 0.0) out_rsd=sqrt(var_rsd**2+dens_rsd**2) 
     
      end

c----------------------------------------------------------------------

      subroutine mcuedivide(out,var,dens,out_rsd,var_rsd,dens_rsd)
      implicit none

      Use(Dim)				# nx,ny
      Use(MCN_dim)          # nfl

c ... Input Variables
      real var(0:nx+1,0:ny+1,1:nfl),dens(0:nx+1,0:ny+1,1:nfl)
      real var_rsd(0:nx+1,0:ny+1,1:nfl),dens_rsd(0:nx+1,0:ny+1,1:nfl)

c ... Output Variables
      real out(0:nx+1,0:ny+1,1:nfl),out_rsd(0:nx+1,0:ny+1,1:nfl)     

      out=0.;out_rsd=1.;
      where(dens .ne. 0.0) out=var/dens
      where(dens .ne. 0.0) out_rsd=sqrt(var_rsd**2+dens_rsd**2)
     
      end

c----------------------------------------------------------------------

      subroutine mcnrsdfix(mcrsd)
      implicit none

      Use(Dim)				# nx,ny,ngsp
      Use(MCN_dim)          # nfl

c ... Input/Output Variables
      real mcrsd(0:nx+1,0:ny+1,1:nfl)

      where(mcrsd .le. 0 .or. mcrsd .gt. 1) mcrsd=1

      end

c----------------------------------------------------------------------

      subroutine mcnblend(out,uevar,mcvar,outrsd,mcrsd,alpha)
      implicit none

      Use(Dim)				# nx,ny,ngsp
      Use(MCN_dim)          # nfl

c ... Input Variables
      real, intent(in), dimension(0:nx+1,0:ny+1,1:nfl) :: uevar,mcvar,mcrsd 
      real, intent(in) :: alpha

c ... Output Variables
      real, intent(inout), dimension(0:nx+1,0:ny+1,1:nfl) :: outrsd  
      real, intent(out), dimension(0:nx+1,0:ny+1,1:nfl) :: out   

c ... Local Variables
      real, dimension(0:nx+1,0:ny+1,1:nfl) :: mcfrac 

      call mcnrsdfix(mcrsd)
      mcfrac=(1-mcrsd**2)**alpha
c      mcfrac=(1-mcrsd**alpha)
      out = (1-mcfrac)*uevar + mcfrac*mcvar 
      where(out .ne. 0) outrsd = mcrsd*mcfrac*mcvar/out
      where(out .eq. 1) outrsd = 1

c      write(*,*) "mcnblend ",out
      end

c-----------------------------------------------------------------------

      function mult23(var2,var3,n3)
      implicit none

      Use(Dim)				# nx,ny
      Use(MCN_dim)

c ... Input Variables
      integer n3
      real var2(0:nx+1,0:ny+1)
      real var3(0:nx+1,0:ny+1,n3)

c ... Output Variables
      real mult23(0:nx+1,0:ny+1,n3)    

c ... Local Variables
      integer i

      do i=1,n3
        mult23(:,:,i)=var3(:,:,i)*var2
      end do
      
      end function mult23

c-----------------------------------------------------------------------

      function mult24(var2,var4,n3,n4)
      implicit none

      Use(Dim)				# nx,ny

c ... Input Variables
      integer n3,n4
      real var2(0:nx+1,0:ny+1)
      real var4(0:nx+1,0:ny+1,n3,n4)

c ... Output Variables
      real mult24(0:nx+1,0:ny+1,n3,n4)   

c ... Local Variables
      integer i,j

      do i=1,n3
      do j=1,n4
        mult24(:,:,i,j)=var4(:,:,i,j)*var2
      end do
      end do

      end function mult24

c-----------------------------------------------------------------------

      function mult34(var3,var4,n3,n4)
      implicit none

      Use(Dim)				# nx,ny

c ... Input Variables
      integer n3, n4
      real var3(0:nx+1,0:ny+1,n3)
      real var4(0:nx+1,0:ny+1,n3,n4)

c ... Output Variables
      real mult34(0:nx+1,0:ny+1,n3,n4)    

c ... Local Variables
      integer i

      do i=1,n4
        mult34(:,:,:,i)=var4(:,:,:,i)*var3
      end do

      end function mult34

c-----------------------------------------------------------------------

c******************************************************
c Subroutines for testing external types in VDF file
c******************************************************

      subroutine testfunfun(f,x,y)
      implicit none
	  external f
      real x,y,f
      y=f(x)
      print *,'y=f(x)'
      end

      real function testfun(x)
      implicit none
      real x
      testfun=2*x
      end function testfun

c      real function norm(v)
c      implicit none
c      real :: v 
c      norm=sum(v**2)
c      norm=sqrt(norm)
c      end function norm
    
      subroutine test_opt(optarg)
c ... This fails when argument not present: parxfcn: incorrect number of argument 
      implicit none

      logical test

c      real opt, optarg
c      opt=3.14
c      if(optarg .ne. real(DEFAULT)) opt=optarg

      character*128 opt, optarg
      opt='test_string'
      test=(optarg==' y')
      write(*,*) "test =", test,", optarg=",optarg,", opt=",opt
      if(.not.(optarg==' y')) opt=trim(optarg)

      call remark(trim(opt))
      end

      subroutine test_parser(optarg)
c ... This fails when argument not present: parxfcn: incorrect number of argument 
      implicit none

c     Use(Dim)
      Use(PNC_params)

      character*128 opt, optarg

      opt=trim(adjustl(optarg))//char(0)
c     opt='remark 1'
c     opt='dobalance(pnc_balancefile)'

      write(*,*) "opt=",opt
      call parsestr(opt)

      end
