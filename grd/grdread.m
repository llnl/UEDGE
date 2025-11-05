c!include "../mppl.h"
c     ------------------------------------------------------------------

      subroutine copyflx
      implicit none
Use(Dimflxgrd)
Use(Comflxgrd)
Use(Curves)
      integer j,n

      do j=1,jdim
        npointg(j)=npoint(j)
        do n=1,npoint(j)
          xcurveg(n,j)=xcurve(n,j)
          ycurveg(n,j)=ycurve(n,j)
        enddo
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine readflx
      implicit none
Use(Share)	#nycore,nysol,igrid
Use(Dimflxgrd)	#jdim,noregs,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)
Use(Aeqflxgrd)
Use(Dimensions)
Use(Curves)
Use(Linkco)
Use(Transfm)
Use(Spline)
Use(Transit)
Use(System)
Use(Mmod)
Use(Argfc)
      integer iunit, ios
      external remark, xerrab, gallot, gchange
      external rdflx1, rdflx2, rdflx3

c **************** read the output of the flx package ******************

      data iunit /55/
      open (iunit, file='flx-grd', form='unformatted', iostat=ios,
     &      status='old')
      if (ios .ne. 0) then
         call xerrab("**** flx-grd file not found")
      endif

c     read dimensioning parameters and allot storage space for arrays --

      read(iunit) jdim,noregs,npts,nycore(igrid),nysol(igrid)
      read(iunit) jmin,jmax,jsptrx,jaxis
      call gallot("Dimensions",0)
      call gallot("Curves",0)
      call rdflx1(iunit)

      read(iunit) nxefit,nyefit
c     set length for 2-d spline workspace --
      nwork = nxefit*nyefit +
     .         2*max(kxord*(nxefit+1),kyord*(nyefit+1))
      call gallot("Comflxgrd",0)
      call rdflx2(iunit)

      read(iunit) nlim
      call gchange("Comflxgrd",0)
      call rdflx3(iunit)

      read(iunit) eshot,etime,rseps,zseps,
     &             rvsin,zvsin,rvsout,zvsout

      read(iunit) xlbnd,xubnd,ylbnd,yubnd
      read(iunit) runid
      read(iunit) geometry

      close (iunit)

      call setidim    # sets angle-like indices and dimensioning parameters
                      # for Inmesh, Linkco and Mmod groups
      call gchange("Inmesh",0)
      call gallot("Linkco",0)
      call gallot("Transfm",0)
      call gallot("Spline",0)
      call gallot("Transit",0)
      call gallot("System",0)
      call gchange("Mmod",0)
c     Arrays for SLATEC spline routine FC --
      ndata=npts
      nbkpt=npts
      nwdim=8*npts+11	# >=5*(nbkpt-1)+2*max(ndata,nbkpt)+nbkpt+16
      niwdim=2*npts-6	# >=2*(nbkpt-3)
      call gchange("Argfc",0)

      return
      end

c     ------------------------------------------------------------------

      subroutine rdflx1(iunit)
      implicit none
Use(Dimflxgrd)	#jdim,npts
Use(Dimensions)
Use(Curves)
      integer iunit

      read(iunit) npointg,xcurveg,ycurveg

      return
      end

c     ------------------------------------------------------------------

      subroutine rdflx2(iunit)
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      integer iunit

      read(iunit) xold,yold,fold
      read(iunit) fpol
      read(iunit) bcentr,rcentr,rmagx,zmagx,simagx,sibdry,
     &             rgrid1,xdim,zdim,zmid

      return
      end

c     ------------------------------------------------------------------

      subroutine rdflx3(iunit)
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      integer iunit

      read(iunit) xlim,ylim

      return
      end

c     ------------------------------------------------------------------

      subroutine readgridpars(fname, runid)
      implicit none
Use(Share)            # geometry
Use(Dim)              # nxm,nym
Use(Xpoint_indices)   # ixlb,ixpt1,ixmdp,ixpt2,ixrb,iysptrx1,iysptrx2
Use(Flags)            # iprint
      character*(*) fname, runid
      integer nuno,ios
      real simagxs_tmp, sibdrys_tmp
      external freeus,remark,xerrab,gallot,rdgrid

      if (isgriduehdf5 .eq. 1) then
        call parsestr('import uedge.gridue as gue;gue.read_gridpars()')
      else
c     Read mesh parameters from a UEDGE code grid data file
      simagxs_tmp=0
      sibdrys_tmp=0
      call freeus (nuno)
      if (iprint .ne. 0) write(*,*) 'Reading grid from file:',trim(fname)
      open (nuno, file=trim(fname), form='formatted', iostat=ios,
     &      status='old')
      if (ios .ne. 0) then
         call xerrab("**** requested grid data file not found")
      endif

      if (geometry=="dnull" .or. geometry=="snowflake15" .or.
     .    geometry=="snowflake45" .or. geometry=="snowflake75" .or.
     .    geometry=="snowflake105" .or. geometry=="snowflake135" .or.
     .    geometry=="snowflake165" .or.
     .    geometry=="dnXtarget" .or. geometry=="isoleg") then
         read(nuno,1999) nxm,nym
         read(nuno,1999) iysptrx1(1),iysptrx2(1)
         read(nuno,1999) ixlb(1),ixpt1(1),ixmdp(1),ixpt2(1),ixrb(1)
         read(nuno,1999) iysptrx1(2),iysptrx2(2)
         read(nuno,1999) ixlb(2),ixpt1(2),ixmdp(2),ixpt2(2),ixrb(2)
	 if (geometry=="dnXtarget") nxc = ixmdp(1)
      else
         read(nuno,'( 5i4, :, f16.10, f16.10 )') nxm,nym,ixpt1(1),ixpt2(1),iysptrx1(1), simagxs_tmp, sibdrys_tmp
         simagxs = simagxs_tmp
         sibdrys = sibdrys_tmp
         ixlb(1)=0
         ixrb(1)=nxm
         iysptrx2(1)=iysptrx1(1)
      endif

      close (nuno)

 1999 format(5i4)
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine readgrid(fname, runid)
      implicit none
Use(Share)            # geometry, isgriduehdf5
Use(Dim)              # nxm,nym
Use(Xpoint_indices)   # ixlb,ixpt1,ixmdp,ixpt2,ixrb,iysptrx1,iysptrx2
Use(Flags)            # iprint
      character*(*) fname, runid
      integer nuno,ios
      real simagxs_tmp, sibdrys_tmp
      external freeus,remark,xerrab,gallot,rdgrid

      if (isgriduehdf5 .eq. 1) then
        call parsestr('import uedge.gridue as gue;gue.read_gridue()')
      else

c     Read a UEDGE code grid data file
      simagxs_tmp=0
      sibdrys_tmp=0
      call freeus (nuno)
      if (iprint .ne. 0) write(*,*) 'Reading grid from file:',trim(fname)
      open (nuno, file=trim(fname), form='formatted', iostat=ios,
     &      status='old')
      if (ios .ne. 0) then
         call xerrab("**** requested grid data file not found")
      endif

      if (geometry=="dnull" .or. geometry=="snowflake15" .or.
     .    geometry=="snowflake45" .or. geometry=="snowflake75" .or.
     .    geometry=="snowflake105" .or. geometry=="snowflake135" .or.
     .    geometry=="snowflake165" .or.	  
     .    geometry=="dnXtarget" .or. geometry=="isoleg") then
         read(nuno,1999) nxm,nym
         read(nuno,1999) iysptrx1(1),iysptrx2(1)
         read(nuno,1999) ixlb(1),ixpt1(1),ixmdp(1),ixpt2(1),ixrb(1)
         read(nuno,1999) iysptrx1(2),iysptrx2(2)
         read(nuno,1999) ixlb(2),ixpt1(2),ixmdp(2),ixpt2(2),ixrb(2)
	 if (geometry=="dnXtarget") nxc = ixmdp(1)
      else
         read(nuno,'( 5i4, :, f16.10, f16.10 )' ) nxm,nym,ixpt1(1),ixpt2(1),iysptrx1(1), simagxs_tmp, sibdrys_tmp
         simagxs = simagxs_tmp
         sibdrys = sibdrys_tmp
         ixlb(1)=0
         ixrb(1)=nxm
         iysptrx2(1)=iysptrx1(1)
      endif
 1999 format(5i4)
      call gallot("RZ_grid_info",0)
      call rdgrid(nuno, runid)

      close (nuno)
      endif # end isgriduehdf5 check

      return
      end

c     ------------------------------------------------------------------

      subroutine rdgrid(nuno, runid)
      implicit none
      integer nuno
      character*(*) runid

Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b

      integer ix,iy,n

      read(nuno,2000)
      read(nuno,2001) (((rm(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((zm(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((psi(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((br(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((bz(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((bpol(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((bphi(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2000)
      read(nuno,2001) (((b(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      read(nuno,2002) runid
 2000 format()
ifelse([WORDSIZE],64,\
 2001 format(1p3e23.15)
,\
 2001 format(1p3d23.15)
)\
 2002 format(a60)

      return
      end

#----------------------------------------------------------------------#

      subroutine setidim
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Dimensions)
Use(Share)         # islimon
Use(Linkco)
Use(Inmesh)
Use(Mmod)
Use(Xmesh)
      integer region
      external gchange

c     Set angle-like parameters and allocate space for arrays --

      if ( (geometry .eq. "dnbot") .or. (geometry .eq. "dnull") .or.
     .     (geometry == "isoleg") .or. (islimon .ne. 0) ) then
         nxuse(1) = max(0,nxcore(igrid,1)-1)  # nxcore includes guard cells
         nxuse(2) = nxcore(igrid,2) - 1       # at the internal mesh boundaries
      else
         nxuse(1) = nxcore(igrid,1)           # nxcore specifies the number
         nxuse(2) = nxcore(igrid,2)           # of finite-size cells
      endif

c     Set some angle-surface parameters --
      idim=0
      do region=1,noregs
         ixpoint(1,region) = nxuse(region) + 1
         ixpoint(2,region) = ixpoint(1,region) + 1
         ixpoint(3,region) = ixpoint(2,region) + 1
         ilmax(region) = ixpoint(3,region) + nxleg(igrid,region)
         idim = max( idim, ilmax(region) )
      enddo

c     Allocate space for angle-dependent arrays --
         call gchange("Linkco",0)
         call gchange("Inmesh",0)
         call gchange("Mmod",0)

c     and for poloidal mesh distribution data --
         call gchange("Xmesh",0)

      return
      end

