c----------------------------------------------------------------------c

      subroutine readpost
      implicit none
      character*(500) fname
      logical fileExists
Use(P93dat)     # atn,atw,nt,nr,nn
Use(Impdata)    # coronalimppath, coronalimpfilename

c     local variables --
      integer ios, nget

c     procedures --
      external xerrab, remark, gallot, readpost1

c----------------------------------------------------------------------c
c     Set up path to mist.dat file
c----------------------------------------------------------------------c


       fname=TRIM(coronalimppath) // '/'//TRIM(coronalimpfname)
       INQUIRE(FILE=TRIM(fname),EXIST=fileExists)
       if (fileExists .neqv. .TRUE.) then
         fname=TRIM(coronalimpfname)
         INQUIRE(FILE=TRIM(fname),EXIST=fileExists)
         if (fileExists .neqv. .TRUE.) then
            write(*,*) "**** Coronal Equilibrium data file not found! ****"
            write(*,*) ""
            write(*,*) "Cannot find "//TRIM(coronalimpfname)//" in:"
            write(*,*) TRIM(coronalimppath)
            write(*,*) " or current directory."
            write(*,*) ""
            write(*,*) "Specify the file name using coronalimpfname"
            write(*,*) "and its path using coronalimppath"
            call xerrab("")
         endif
       endif


c----------------------------------------------------------------------c
c     Read impurity emissivity and charge state from POST93 data files
c----------------------------------------------------------------------c

      data nget /55/
      open (nget, file=fname, form='formatted', iostat=ios,
     .     status='old')
      if (ios .ne. 0) then
         call remark("**** data file mist.dat not found --")
         call xerrab("")
      endif

c     read array dimensions
      read(nget,9000)
      read(nget,9000)
      read(nget,9000)
      read(nget,9000)
      read(nget,9001) atn, atw, nt, nr, nn

c     allocate storage for arrays
      call gallot("P93dat",0)

c     read data arrays
      call readpost1 (nget)

      close (nget)

 9000 format()
 9001 format(5(1x,i2/))

      return
      end

c----------------------------------------------------------------------c

      subroutine readpost1 (nget)
      implicit none
      integer nget
Use(P93dat)   # nt,nr,nn,tdatm,rdatm,ndatm,emdatm,z1datm,z2datm
Use(Physical_constants2)  # ev2

c     local variables --
      integer jn, jr, jt

      do jn=1,nn
         do jr=1,nr
            do jt=1,nt
               read(nget,9010) tdatm(jt,jr,jn), rdatm(jt,jr,jn), 
     &                         ndatm(jt,jr,jn), emdatm(jt,jr,jn), 
     &                         z1datm(jt,jr,jn), z2datm(jt,jr,jn)
            enddo
         enddo
      enddo

c     convert to SI units:
      do jt=1,nt
         do jr=1,nr
            do jn=1,nn
               tdatm(jt,jr,jn)=tdatm(jt,jr,jn)*ev2
               ndatm(jt,jr,jn)=ndatm(jt,jr,jn)*1.0e+06
               emdatm(jt,jr,jn)=emdatm(jt,jr,jn)*1.0e-06
            enddo
         enddo
      enddo

 9010 format(6(1x,e12.5))

      return
      end

c----------------------------------------------------------------------c

      subroutine splinem
      implicit none
Use(P93dat)    # nt,nr,nn
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata
      external gallot, splinem1

c     Construct 3-dimensional B-spline representation for impurity
c     emissivity and charge state versus e-temperature, ng/ne and
c     n*tau (data from POST '93 tables)

c     Allocate arrays for spline fitting --
      nxdata_api=nt		# temperature
      nydata_api=nr		# density ratio
      nzdata=nn		# n*tau
      nwork2 = kyords_api*kzords + 3*max(kxords_api,kyords_api,kzords) + kzords + 2
      nwork3 = nxdata_api*nydata_api*nzdata + 2*max( kxords_api*(nxdata_api+1),
     &           kyords_api*(nydata_api+1), kzords*(nzdata+1) )
      call gallot("Imslwrk",0)

      call splinem1

      return
      end

c----------------------------------------------------------------------c

      subroutine splinem1
      implicit none
Use(P93dat)    # tdatm,rdatm,ndatm,emdatm,z1datm,z2datm
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,fdata_api,ldf_api,mdf,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,emcoef,
               # z1coef,z2coef
      integer i,j,k
      external B3INT

c     Define data arrays --
      do i=1,nxdata_api
        xdata_api(i)=log10( tdatm(i,1,1) )
      enddo
      do j=1,nydata_api
        ydata_api(j)=log10( rdatm(1,j,1) )
      enddo
      do k=1,nzdata
        zdata(k)=log10( ndatm(1,1,k) )
      enddo

c     Define the order of the spline fit
c      kxords_api=4		# cubic in x=log10(temperature)
c      kyords_api=4		# cubic in y=log10(density ratio)
c      kzords=4		# cubic in z=log10(n*tau)

      ldf_api=nxdata_api
      mdf=nydata_api
      iflagi = 1	# let B3INT choose knots

c     Compute the coefficients --
c     first, for emissivity:
      do i=1,nxdata_api
        do j=1,nydata_api
          do k=1,nzdata
            fdata_api(i,j,k)=log10( emdatm(i,j,k) )
            emcoef(i,j,k)=fdata_api(i,j,k)
          enddo
        enddo
      enddo
      call B3INT (xdata_api, nxdata_api, ydata_api, nydata_api, zdata, nzdata,
     &            kxords_api, kyords_api, kzords, xknots_api, yknots_api, zknots,
     &            emcoef, ldf_api, mdf, work3, iflagi)
c     next, for average Z:
      do i=1,nxdata_api
        do j=1,nydata_api
          do k=1,nzdata
            fdata_api(i,j,k)=z1datm(i,j,k)
            z1coef(i,j,k)=fdata_api(i,j,k)
          enddo
        enddo
      enddo
      call B3INT (xdata_api, nxdata_api, ydata_api, nydata_api, zdata, nzdata,
     &            kxords_api, kyords_api, kzords, xknots_api, yknots_api, zknots,
     &            z1coef, ldf_api, mdf, work3, iflagi)
c     next, for average Z**2:
      do i=1,nxdata_api
        do j=1,nydata_api
          do k=1,nzdata
            fdata_api(i,j,k)=z2datm(i,j,k)
            z2coef(i,j,k)=fdata_api(i,j,k)
          enddo
        enddo
      enddo
      call B3INT (xdata_api, nxdata_api, ydata_api, nydata_api, zdata, nzdata,
     &            kxords_api, kyords_api, kzords, xknots_api, yknots_api, zknots,
     &            z2coef, ldf_api, mdf, work3, iflagi)

      return
      end

c-----------------------------------------------------------------------

      real function emissbs (vte,vnr,vnt)
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,emcoef

      integer nxcoef,nycoef,nzcoef
      real vlogw,w,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     emissivity [Watts-m**3] versus :
c     		vte = e-temperature [J]
c     		vnr = ng/ne density ratio
c     		vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      icont = 0
      vlogw=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords, 
     &             emcoef, ldf_api, mdf, icont, iworki, work2, iflagi)
      w=10**vlogw
      emissbs=w

      return
      end

c----------------------------------------------------------------------c

      real function z1avgbs (vte,vnr,vnt)
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,z1coef

      integer nxcoef,nycoef,nzcoef
      real vz1,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     charge state Z versus :
c     		vte = e-temperature [J]
c     		vnr = ng/ne density ratio
c     		vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      icont = 0
      vz1=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords, 
     &             z1coef, ldf_api, mdf, icont, iworki, work2, iflagi)
      z1avgbs=vz1

      return
      end

c-----------------------------------------------------------------------

      real function z2avgbs (vte,vnr,vnt)
      implicit none
      real vnt,vnr,vte
Use(Imslwrk)   # nxdata_api,nydata_api,nzdata,xdata_api,ydata_api,zdata,
               # kxords_api,kyords_api,kzords,xknots_api,yknots_api,zknots,z2coef

      integer nxcoef,nycoef,nzcoef
      real vz2,xuse,yuse,zuse

      real B3VAL
      external B3VAL

c     Evaluate 3-dimensional B-spline representation for impurity
c     charge state Z**2 versus :
c     		vte = e-temperature [J]
c     		vnr = ng/ne density ratio
c     		vnt = ne*tau [sec/m**3].

      xuse=min(max(xdata_api(1),log10(vte)),xdata_api(nxdata_api))
      yuse=min(max(ydata_api(1),log10(vnr)),ydata_api(nydata_api))
      zuse=min(max(zdata(1),log10(vnt)),zdata(nzdata))

      nxcoef=nxdata_api
      nycoef=nydata_api
      nzcoef=nzdata

      icont = 0
      vz2=B3VAL (xuse, yuse, zuse, 0, 0, 0, xknots_api, yknots_api, zknots,
     &             nxcoef, nycoef, nzcoef, kxords_api, kyords_api, kzords, 
     &             z2coef, ldf_api, mdf, icont, iworki, work2, iflagi)
      z2avgbs=vz2

      return
      end

