c!include "../mppl.h"
c----------------------------------------------------------------------c

      subroutine ingrd
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)
Use(Dimensions)
Use(Curves)
Use(Inmesh)	#isspnew,rstrike,zstrike
		#istpnew,rtpnew,ztpnew
		#x0g,y0g,xlast,ylast
Use(Linkco)
Use(Mmod)	# delmax
Use(Xmesh)
Use(Share)

c     local variables --
      integer n
      real dxefit,dyefit,rmin,rmax

c     Set input parameters for the grd package --

      dxefit = xdim/(nxefit-1)
      dyefit = zdim/(nyefit-1)
      if (dsmin .eq. 0) then  # use this value if user does not supply one
         dsmin = 0.25*dxefit
      endif
      dsminx = 2.0*dsmin
      dyjump = 1.0*dyefit
      delmax = 0.01*dxefit

c     Set data for first and last seedpoints on the separatrix
c     for inboard and outboard regions of the mesh --

c     First seedpoints --
c     The uppermost point for each region
      if (istpnew .eq. 0) then		# use default specification
         if ((geometry .eq. 'snull').or.(geometry .eq. "uppersn")) then
c           choose radial position of magnetic axis
            x0g(1) = rmagx
            x0g(2) = rmagx
         else	# geometry="dnbot"
c           assume up/down symmetric double-null configuration
c           choose innermost and outermost points on the separatrix
            rmin=rgrid1+xdim
            do n=1,npointg(jsptrx(1))
               if (xcurveg(n,jsptrx(1)) .lt. rmin) then
                  rmin=xcurveg(n,jsptrx(1))
               else
                  break
               endif
            enddo
            rmax=rgrid1
            do n=1,npointg(jsptrx(2))
               if (xcurveg(n,jsptrx(2)) .gt. rmax) then
                  rmax=xcurveg(n,jsptrx(2))
               else
                  break
               endif
            enddo
            x0g(1)=rmin
            x0g(2)=rmax
         endif	# end if-test on "geometry"
      elseif (istpnew .eq. 1) then	# user-specified values
         x0g(1) = rtpnew(1)
         y0g(1) = ztpnew(1)
         x0g(2) = rtpnew(2)
         y0g(2) = ztpnew(2)
      endif

c     Last seedpoints --
c     at the strike point for each region
      if (isspnew .eq. 0) then		# use eqdsk values
ccc Temporary code for FRC mesh generation --
         if (nxleg(igrid,1)==0) then
            xlast(1) = rseps
            ylast(1) = zseps
         else
            xlast(1) = rvsin
            ylast(1) = zvsin
         endif
         xlast(2) = rvsout
         ylast(2) = zvsout
      elseif (isspnew .eq. 1) then	# user-specified values
         xlast(1) = rstrike(1)
         ylast(1) = zstrike(1)
         xlast(2) = rstrike(2)
         ylast(2) = zstrike(2)
      endif

c     Definition of "conjugate magnetic axis" (in private region) --
      rtanpl=(rvsin+rvsout)/2
      ztanpl=(zvsin+zvsout)/2

      return
      end

c     ------------------------------------------------------------------

      subroutine idealgrd

*  -- Subroutine for initializing the arrays rm, zm, etc. in group
*  -- RZ_grid_info for cylindrical IDEAL. Re-allocation does not take
*  -- effect in the same routine where nxm and nym are changed; thus, the
*  -- computations are done in a following routine, called idlcomp.

      implicit none

      Use(Share)      # igrid,nxomit

      external com_set_dims, gallot

c ... Compute dimension variables nx, nym, and nxm.
      igrid = 1
      nxomit = 0
      call com_set_dims

c ... Allocate space, then compute.
      call gallot("RZ_grid_info",0)
      call idlcomp

      return
      end

c     --------------------------------------------------------------------

      subroutine idlcomp

*  -- Subroutine for computing the vertex locations of the cylindrical
*  -- or rectangular grid (formerly, the IDEAL grid)

      implicit none

      Use(Dim)              # nxm,nym
      Use(RZ_grid_info)     # rm,zm,psi,br,bz,bpol,bphi,b
      Use(Xpoint_indices)   # ixpt2,iysptrx1
      Use(Analgrd)          # alfyt,alfxt,tctr,za0,zaxpt,zax,radm,radx,rad0,
                            # btfix,bpolfix,isgdistort,ixdstar,iynod,rnod
                            # agsidx,agsrsp,sratiopf,rscalcore,tnoty,isadjalfxt
      Use(Share)            # nxcore,nxleg,isgrdsym

*  -- local variables
      real dznl, dznr, dznu, t1x, t2x, t1y, t2y, rixm, adrt, 
     .     xshift, radmu, ta_r,drn1,drn2,drnu1,drnu2
      real btorfix
      character*60 runidl
      integer ix, iy, n, ixm, ixpt2_init
      real lpi
      data lpi/3.141592654/

      runidl = 'ideal geometry'

      ixpt2_init = nxcore(1,2)  # define to allow isgrdsym=1 case to work
      if (iysptrx1(1) .eq. 0) rad0 = radm  #allows old cases without rad0 to run
      if (ixpt2_init .gt. 0) dznu = (zaxpt - za0) / float(ixpt2_init)
      if (tctr .gt. 0.0000001) dznl = tctr*(zax - zaxpt - za0) / 
     .                               (1 - exp(+alfxt*tctr))
      if (tctr .lt. 0.9999999) dznr = (1 - tctr)*(zax - zaxpt - za0) / 
     .                               (1 - exp(-alfxt*(1-tctr)))
cc      write(*,*) 'dznu, dznr =',dznu,dznr
      if (isadjalfxt > 0) then  #change alfxt slightly so dznr=dznu
        alfxt = -log(1. - (1.-tctr)*(zax - zaxpt - za0)*
     .                      (1.-exp(-alfxt/float(nxm-ixpt2_init)))/dznu)
     ,                                                      /(1.-tctr)
        write(*,*) '*****'
        write(*,*) 'alfxt adjusted, smooth dx at ixpt2: alfxt=',alfxt
        write(*,*) '*****'
      endif
      
      rixm = tctr * float(nxm-ixpt2_init) + .5
      ixm = int(rixm)
      
      do 20 iy = nym, 1, -1
         if (iy .gt. iysptrx1(1)) then
            t1y = float(iy-1-iysptrx1(1)) / float(nym-iysptrx1(1))
            t2y = float(iy-iysptrx1(1)) / float(nym-iysptrx1(1))
            radmu = rad0
            if (tnoty .eq. 0) then
              drn1=(radx - rad0)*(1-exp(-alfyt*t1y))/(1 - exp(-alfyt))
              drn2=(radx - rad0)*(1-exp(-alfyt*t2y))/(1 - exp(-alfyt))
            else  # new (6/99) generalized form; note, make alfyt=-alfyt
              drn1=(radx-rad0)*
     .                 (tanh(alfyt*tnoty)+tanh(alfyt*(t1y-tnoty)))/
     .                   (tanh(alfyt*tnoty)+tanh(alfyt*(1-tnoty)))
              drn2=(radx-rad0)*
     .                 (tanh(alfyt*tnoty)+tanh(alfyt*(t2y-tnoty)))/
     .                   (tanh(alfyt*tnoty)+tanh(alfyt*(1-tnoty)))
            endif
         else
            if (sratiopf .eq. 0.) sratiopf = (rad0-radm)/(radx-rad0)   
                                      #gives same expan. rate around sep.
            radmu = rad0
            if (tnoty.eq.0) then
              t1y = -sratiopf*float(iy-1-iysptrx1(1)) / float(iysptrx1(1))
              t2y = -sratiopf*float(iy-iysptrx1(1)) / float(iysptrx1(1))
              drn1 = (radm - rad0)*(1-exp(-alfyt*t1y)) / 
     .                                       (1 - exp(-sratiopf*alfyt))
              drn2 = (radm - rad0)*(1-exp(-alfyt*t2y)) / 
     .                                       (1 - exp(-sratiopf*alfyt))
            else  # new (6/99) generalized form; note, not sratiopf
              drn1 = (rm(1,iysptrx1(1)+1,2)-rm(1,2*iysptrx1(1)+1-iy,3))
              drn1 = -drn1*(rad0-radm)/
     .                    (rm(1,iysptrx1(1)+1,2)-rm(1,2*iysptrx1(1),3))
              drn2 = (rm(1,iysptrx1(1)+1,2)-rm(1,2*iysptrx1(1)+1-iy,2))
              drn2 = -drn2*(rad0-radm)/
     .                    (rm(1,iysptrx1(1)+1,2)-rm(1,2*iysptrx1(1),3))
             endif
         endif          

         do 10 ix = 1, nxm
            if (ix.le.ixpt2_init .and. iy.le.iysptrx1(1)) then
               drnu1 = rscalcore*drn1  # give diff radial range to core region
               drnu2 = rscalcore*drn2  
            else
               drnu1 = drn1
               drnu2 = drn2
            endif
            rm(ix,iy,1) = radmu + drnu1
            rm(ix,iy,2) = rm(ix,iy,1)
            rm(ix,iy,3) = radmu + drnu2
            rm(ix,iy,4) = rm(ix,iy,3)
            rm(ix,iy,0) = 0.25 * ( rm(ix,iy,1) + rm(ix,iy,2) +
     .                             rm(ix,iy,3) + rm(ix,iy,4) )

            if (ix .le. ixpt2_init) then      # uniform grid region
               zm(ix,iy,1) = za0 + dznu * (ix-1)
               zm(ix,iy,2) = za0 + dznu * ix
            elseif (ix .gt. ixm) then    # nonuniform grid with decreasing dx
               t1x = float(ix-ixpt2_init-1) / float(nxm-ixpt2_init)
               t2x = float(ix-ixpt2_init) / float(nxm-ixpt2_init)
               zm(ix,iy,1) = zaxpt + dznr * (1-exp(-alfxt*(t1x-tctr))) 

               zm(ix,iy,2) = zaxpt + dznr * (1-exp(-alfxt*(t2x-tctr))) 
            else                         # nonuniform grid with increasing dx
               t1x = float(ix-ixpt2_init-1) / float(nxm-ixpt2_init)
               t2x = float(ix-ixpt2_init) / float(nxm-ixpt2_init)
               zm(ix,iy,1) = zaxpt + dznl * (1-exp(+alfxt*t1x)) 
               zm(ix,iy,2) = zaxpt + dznl * (1-exp(+alfxt*t2x)) 
            endif
            zm(ix,iy,3) = zm(ix,iy,1)
            zm(ix,iy,4) = zm(ix,iy,2)

            zm(ix,iy,0) = 0.25 * ( zm(ix,iy,1) + zm(ix,iy,2) +
     .                             zm(ix,iy,3) + zm(ix,iy,4) )
  10     continue
  20  continue

ccc      write(*,*) 'zm(,1,1)]',zm(1:nxm,1,1)

c...  Tilt the mesh in the poloidal direction from ixsnog if tiltang.ne.0
      if (tiltang.ne.0.) then
         do iy = 1, nym
            do ix = ixsnog, nx
               ta_r = (zm(ix,1,2)-zm(ixsnog,1,1))*tiltang*(lpi/180.) /
     .                (zm(nx,1,2)-zm(ixsnog,1,1))
               zm(ix,iy,4) = zm(ix,1,2) + (rm(ix,iy,3)-rm(ix,1,1)) *
     .                                                       tan(ta_r)
               if (ix.gt.1) zm(ix,iy,3) = zm(ix-1,iy,4)
               if (iy.gt.1) then
                  zm(ix,iy,1) = zm(ix,iy-1,3)
                  zm(ix,iy,2) = zm(ix,iy-1,4)
               endif
               zm(ix,iy,0) = 0.25 * ( zm(ix,iy,1) + zm(ix,iy,2) +
     .                                zm(ix,iy,3) + zm(ix,iy,4) )
            enddo
         enddo
      endif

c...  Distort mesh poloidally to test nonorthogonal coding/performance
      if (isgdistort .eq. 1) then   # only zm changes, rm the same
         do iy = 1, nym
            adrt = agsindx*(iy-iynod)/nym + 
     .                            agsrsp*(rm(1,iy,3)-rnod)/rm(1,nym,3)
            do ix = ixdstar, nx
               xshift = adrt*sin((ix-ixdstar+1)*lpi/(nx-ixdstar+1))
               zm(ix,iy,4) = zm(ix,iy,4) + (zm(ix,iy,4)-zm(ix,iy,3))*
     .                                     0.1*float(nx-ixdstar)*xshift
            # 0.1*(nx-ixdstar) factor keeps distortion ~same as nx increases
               if (ix.gt.1) zm(ix,iy,3) = zm(ix-1,iy,4)
               if (iy.gt.1) then 
                  zm(ix,iy,1) = zm(ix,iy-1,3)
                  zm(ix,iy,2) = zm(ix,iy-1,4)
               elseif (iy.eq.1) then   # need to do bottom face as iy=0
                  xshift = ( agsindx*(-iynod)/nym + 
     .                           agsrsp*(rm(1,1,1)-rnod)/rm(1,nym,3) ) *
     .                        sin((ix-ixdstar+1)*lpi/(nx-ixdstar+1)) 
                  zm(ix,1,2) = zm(ix,1,2) + (zm(ix,1,2)-zm(ix,1,1))*
     .                                      0.1*float(nx-ixdstar)*xshift
            # 0.1*(nx-ixdstar) factor keeps distortion ~same as nx increases
                  if (ix.gt.1) zm(ix,1,1) = zm(ix-1,1,2)
               endif
               zm(ix,iy,0) = 0.25 * ( zm(ix,iy,1) + zm(ix,iy,2) +
     .                                zm(ix,iy,3) + zm(ix,iy,4) )
            enddo
         enddo
      endif

c...  Make a symmetric grid in zm if isgrdsym=1; redistribute zm's and
      if (isgrdsym .eq. 1) then
         if ((nxm-1)/2.eq.nxm/2 .or. 
     .                          (ixpt2_init-1)/2.eq.ixpt2_init/2) then
            call xerrab('*** nxm or ixpt2_init odd, cannot use isgrdsym ***')
         endif
         do iy = 1, nym  #don't change points ,2 and ,4 until second ix loop
            do ix = 1, nx/2
               zm(ix,iy,1) = -zm((nx-2*(ix-1)),iy,2)
               zm(ix,iy,3) = -zm((nx-2*(ix-1)),iy,4)
               if (ix .eq. 1) then
                  zm(ix+nx/2,iy,1) = 0.
                  zm(ix+nx/2,iy,3) = 0.
               else
                  zm(ix+nx/2,iy,1) = zm(2*(ix-1),iy,2)
                  zm(ix+nx/2,iy,3) = zm(2*(ix-1),iy,4)
               endif
            enddo
            do ix = 1, nx
               if (ix .ne. nx) then   # ix=nx case satisfied from orgin. grid
                  zm(ix,iy,2) = zm(ix+1,iy,1)
                  zm(ix,iy,4) = zm(ix+1,iy,3)
               endif
               zm(ix,iy,0) = 0.25 * ( zm(ix,iy,1) + zm(ix,iy,2) +
     .                                zm(ix,iy,3) + zm(ix,iy,4) )
            enddo
         enddo
      endif          

c...  Reset x-point index ixpt2 if the zxpt_reset is positive; used for
c...  one continuous mesh (zaxpt=0, ncore=0) with no uniform dz region
      if (zxpt_reset .ge. 1.e-10) then
         do ix = 1, nx
            if (zm(ix,1,1).lt.zxpt_reset .and. 
     .                                zm(ix,1,2).ge.zxpt_reset) then
               ixpt2(1) = ix
            endif
         enddo
      endif      

c...  Set constant values for B-fields

      do 40 iy = 1, nym
         do 30 ix = 1, nxm
            do 25 n = 0, 4
# MVU 6-jan-2022
               btorfix=sqrt(btfix**2-bpolfix**2)
               bphi(ix,iy,n) = btorfix*(rm(ix,iy,n)/rmajfix)**sigma_btor
               bpol(ix,iy,n) = bpolfix*(rm(ix,iy,n)/rmajfix)**sigma_bpol
               b(ix,iy,n)    = sqrt(bphi(ix,iy,n)**2 + bpol(ix,iy,n)**2)

               ##-Note: this is poloidal flux per radian, correct only for the cylindrical case
               psi(ix,iy,n)  = ((bpolfix*rmajfix**2)/(sigma_bpol+2)) *
     &               (rm(ix,iy,n)/rmajfix)**(sigma_bpol+2)
		 
# IJ 2016/10/04: radial and vertical components must be set for MCN conversions
               br  (ix,iy,n) = 0.    
               bz  (ix,iy,n) = - bpol (ix,iy,n) 
 25         continue
 30      continue
 40   continue
 
c ... Uncomment the following line to write out the data.
ccc      call writedata ('grididl', runidl)

      return
      end

c     ------------------------------------------------------------------

      subroutine grdgen

      implicit none

Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)	#ylbnd
Use(Dimensions)
Use(Share)       #ishalfm
Use(Inmesh)
Use(Linkco)

c     local variables --
      integer region,rbegin

c     procedures --
      external prune,extend,splfit,sow,meshgen,meshfin

      call prune

      if (yextend .lt. ylbnd) call extend

      if (dxleft .gt. 0) call exleft		# extrapolate to smaller R

      call splfit

      call sow

      rbegin = 1
      if (ishalfm==1) rbegin = 2
      do 90 region=rbegin,noregs
         call meshgen(region)
   90 continue


      call meshfin
ccc   call remark("***** subroutine grdgen completed")

      return
      end

c     ------------------------------------------------------------------
c     ------------------------------------------------------------------

      subroutine torangrd

*  -- Subroutine for initializing the arrays rm, zm, etc. in group
*  -- RZ_grid_info for a toroidal annulus
*  -- effect in the same routine where nxm and nym are changed; thus, the
*  -- computations are done in a following routine, called torancomp.

      implicit none

      Use(Share)      # igrid,nxomit

      external com_set_dims, gallot

c ... Compute dimension variables nx, nym, and nxm.
      igrid = 1
      nxomit = 0
      call com_set_dims

c ... Allocate space, then compute.
      call gallot("RZ_grid_info",0)
      call gallot("Torannulus",0)
      call torancomp

      return
      end

c     --------------------------------------------------------------------
c     ------------------------------------------------------------------
      subroutine torancomp

*  -- Subroutine for computing the vertex locations of a simple
*  -- toroidal annulus

      implicit none

      Use(Dim)              # nxm,nym
      Use(RZ_grid_info)     # rm,zm,psi,br,bz,bpol,bphi,b
      Use(Xpoint_indices)   # ixpt2,iysptrx1
      Use(Torannulus)       # acore,edgewid,dthlim,rm0,bpol0,btor0,
                            # radf,thpf
      Use(Share)            # nxcore,nysol,islimon

*  -- local variables
      real dely, delthp, wshf
      character*60 runtoran
      integer ix, iy, ij
      real pi  # local pi to avoid using full Phyvar group
      data pi/3.14159265358979323/
      runtoran = 'Circular toroidal annulus'

c ... Compute the poloidal angle and radial position of cell faces
      dely = edgewid/(float(nycore(1)))
      if (islimon==0) then   # bottom of annulus is theta=pi/2; clockwise
        delthp = 2*pi/(float(nxcore(1,1)+nxcore(1,2)))
	thpf(1,1) = 0.5*pi  
	thpf(1,2) = 0.5*pi + delthp
	thpf(1,3) = 0.5*pi  
	thpf(1,4) = 0.5*pi + delthp
	thpf(1,0) = 0.5*pi + 0.5*delthp
      else  
        delthp = 2*(pi-dthlim)/(float(nxcore(1,1)+nxcore(1,2)-2)) 
        wshf = 2*float(nxcore(1,2)-1)/float(nxcore(1,1)+nxcore(1,2)-2)           
	thpf(1,1) = 0.5*pi + wshf*dthlim 
	thpf(1,2) = 0.5*pi + wshf*dthlim + delthp
	thpf(1,3) = 0.5*pi + wshf*dthlim 
	thpf(1,4) = 0.5*pi + wshf*dthlim + delthp
        thpf(1,0) = 0.5*(thpf(1,1)+thpf(1,2))
      endif
      do iy = 0, ny
        radf(iy,1) = acore + dely*float(max(iy-1,0))
	radf(iy,2) = acore + dely*float(max(iy-1,0))
        radf(iy,3) = acore + dely*float(iy  )
        radf(iy,4) = acore + dely*float(iy  )
        radf(iy,0) = 0.5*(radf(iy,1)+radf(iy,3))
      enddo

c ... Compute poloidal angles
      do ix = 2, nxcore(1,1)-islimon   # ix=1 set above
        thpf(ix,1) = thpf(ix-1,2)
        thpf(ix,2) = thpf(ix-1,2) + delthp
        thpf(ix,3) = thpf(ix-1,4)
        thpf(ix,4) = thpf(ix-1,4) + delthp
        thpf(ix,0) = 0.5*(thpf(ix,1)+thpf(ix,2))
      enddo
c ... Insert limiter angles if switched on
      if (islimon==1) then
        do ix = nxcore(1,1), nxcore(1,1)+1
          thpf(ix,1) = thpf(ix-1,2) 
          thpf(ix,2) = thpf(ix-1,2) + dthlim
          thpf(ix,3) = thpf(ix-1,4)
          thpf(ix,4) = thpf(ix-1,4) + dthlim
          thpf(ix,0) = 0.5*(thpf(ix,1)+thpf(ix,2))
        enddo
      endif
c ... Finish angles
      do ix = nxcore(1,1)+1+islimon, nxm   # ix=nx set above
        thpf(ix,1) = thpf(ix-1,2)
        thpf(ix,2) = thpf(ix-1,2) + delthp
        thpf(ix,3) = thpf(ix-1,4)
        thpf(ix,4) = thpf(ix-1,4) + delthp
        thpf(ix,0) = 0.5*(thpf(ix,1)+thpf(ix,2))
      enddo

c ... Compute coordinates (rm,zm) for mesh vertices and centers
      do iy = 1, nym
        do ix = 1, nxm
	  do ij = 0, 4
	    rm(ix,iy,ij) = rm0 + radf(iy,ij)*cos(thpf(ix,ij))
	    zm(ix,iy,ij) = -radf(iy,ij)*sin(thpf(ix,ij))
          enddo
        enddo
      enddo

c...  Set B-fields
      do iy = 1, nym
        do ix = 1, nxm
          do ij= 0, 4
            bpol(ix,iy,ij) = bpol0*(rm0/rm(ix,iy,ij))**ibpmodel
            psi (ix,iy,ij) = bpol0*rm0*radf(iy,ij)*
     .                  (1-(radf(iy,ij)*cos(thpf(ix,ij))/rm0)**ibpmodel)
            br  (ix,iy,ij) = -bpol(ix,iy,ij)*sin(thpf(ix,ij))
            bz  (ix,iy,ij) = -bpol(ix,iy,ij)*cos(thpf(ix,ij))
            bphi(ix,iy,ij) = btor0*rm0/rm(ix,iy,ij)
            b   (ix,iy,ij) = sqrt(bphi(ix,iy,ij)**2 + bpol(ix,iy,ij)**2)
          enddo
        enddo
      enddo
 
c ... Uncomment the following line to write out the data.
      call writedata ('gridue', runtoran)

      return
      end
c     ------------------------------------------------------------------

      subroutine mirrorgrd

*  -- Subroutine for initializing the arrays rm, zm, etc. in group
*  -- RZ_grid_info for the annulus of a mirror B-field initially for an FRC
*  -- effect in the same routine where nxm and nym are changed; thus, the
*  -- computations are done in a following routine, called mirrorcomp.

      implicit none

      Use(Share)      # igrid,nxomit
      Use(Dim)        # nx,ny,nxm,nym
      Use(Magmirror)  # zu,ru,bru,bzu,bmag

      external com_set_dims, gallot

c ... Compute dimension variables nx, nym, and nxm.
      igrid = 1
      nxomit = 0
      nxm = nzc
      nym = nrc
      nx = nxm
      ny = nym
##      call com_set_dims

c ... Allocate space, then compute.
      call gallot("RZ_grid_info",0)
      call mirrorcomp

      return
      end
c     --------------------------------------------------------------------

      subroutine mirrorcomp

*  -- Subroutine for computing the vertex locations of a mirror field

      implicit none

      Use(Dim)              # nxm,nym
      Use(RZ_grid_info)     # rm,zm,psi,br,bz,bpol,bphi,b
      Use(Magmirror)        # zu,ru,bru,bzu,bmag

*  -- local variables
      character*60 runmirror
      integer ix, iy, ij
      runmirror = 'Magnetic mirror (FRC-annulus)'

c ... Compute coordinates (rm,zm) for mesh vertices and centers
      do iy = 1, nym
        do ix = 1, nxm
	  do ij = 0, 4
	    rm(ix,iy,ij) = ru(ix,iy,ij)
	    zm(ix,iy,ij) = zu(ix,iy,ij)
          enddo
        enddo
      enddo

c...  Set B-fields
      do iy = 1, nym
        do ix = 1, nxm
          do ij= 0, 4
            br  (ix,iy,ij) = bru(ix,iy,ij)
            bz  (ix,iy,ij) = bzu(ix,iy,ij)
            bpol(ix,iy,ij) = bmag(ix,iy,ij)
            b   (ix,iy,ij) = bmag(ix,iy,ij)
            bphi(ix,iy,ij) = 1.e-20*b(ix,iy,ij)
          enddo
        enddo
      enddo
 
c ... Uncomment the following line to write out the data.
      call writedata ('gridue', runmirror)

      return
      end

c     --------------------------------------------------------------------
      subroutine clear
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Spline)
Use(System)
Use(Transfm)
Use(Transit)

c     local variables --
      integer region,j,k,l

c     This subroutine clears the output and working arrays used in          
c     construction of the mesh

c     loop over all regions --
      do 800 region=1,noregs

      do 2208 j=jmin(region),jmax(region)
         do 2208 l=1,mseg
            do 2208 k=1,npts
               splcoef(k,l,j)=0.0
               xknts(k,l,j)=0.0
 2208 continue
c
      do 2212 j=jmin(region),jmax(region)
         do 2212 l=1,npts
            isegment(l,j)=0
 2212 continue
c
      do 2209 j=jmin(region),jmax(region)
         nseg(j)=0
         do 2209 l=1,mseg
            isys(l,j)=0
            istartg(l,j)=0
            iendg(l,j)=0
            ncap7(l,j)=0
            m(l,j)=0
 2209 continue
c
      do 2210 k=1,npts
         wg(k)=1.0
         xtrans(k)=0.0
         ytrans(k)=0.0
 2210 continue
  800 continue          # end of loop over regions

      return
      end

c     ------------------------------------------------------------------

      subroutine codsys(j,icood,iseg,is,dy,region,alpha1)
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)
Use(Transfm)
Use(Spline)
Use(Inmesh)   # dalpha
Use(Transit)
Use(System)
Use(Argfc)
Use(Share)    # cutlo
      integer j		# flux surface index for this contour
      integer icood	# gives the initial "direction" of the contour
			# trajectory at input data point is
      integer iseg	# on input, this is the index of the previous segment
			# of this contour
			# on output, this is the index of the current segment
			# of this contour
      integer is	# on input, this is the index of the first data point
			# on the current segment of this contour
			# on output, this is the index of the first data point
			# on the next segment of this contour
      real dy		# the "jump" criterion
      integer region	# the inboard/outboard designator
      real alpha1	# the "direction" criterion
c
      logical jump
      integer ibeg, irun, istart1, iend1, irange, i, k, nkint
      real rlimit, z1, z2, tau, s, t, dy1, beta
      external xerrab

c     icood=1, flux curve goes right, transform to system beta=0.
c     icood=2, flux curve goes up, transform to system beta=Pi/2.
c     icood=3, flux curve goes left, transform to system beta=Pi.
c     icood=4, flux curve goes down, transform to system, beta=3*Pi/2.

c     a small real number whose reciprocal does not overflow --
      rlimit = cutlo

c     find the index of the last data point for this spline segment --

         jump=.false.
         ibeg=is
         irun=is
 22      continue
         irun=irun+1
         z1=ycurveg(irun,j)-ycurveg(irun-1,j)
         z2=xcurveg(irun,j)-xcurveg(irun-1,j)

         if(abs(z1) .gt. dy) jump=.true.

         if ((j .eq. jsptrx(region)) .and.
     &             (irun .eq. ixpointc(2,region))) then
            jump=.true.
            z2=1.
         endif

         if (abs(z2) .le. abs(z1)*rlimit) then
            tau=90.
         else
            s=z1/z2
            t=atan(s)
            tau=(180./Pi)*t
         endif
c
         if (jump .or. (irun .gt. npointg(j)))  then
            if (icood.eq.1) tau=1000.
            if (icood.eq.2) tau=0.
            if (icood.eq.3) tau=1000.
            if (icood.eq.4) tau=0.
         endif
      if ((icood .eq. 1) .or. (icood .eq. 3)) then
         if (abs(tau) .le. alpha1+dalpha) goto 22
      endif
      if ((icood .eq. 2) .or. (icood .eq. 4)) then
         if(abs(tau) .gt. alpha1-dalpha) goto 22
      endif
      if ((j .eq. jsptrx(region)) .and.
     &    (irun .eq. ixpointc(1,region)))  goto 22

c     When we leave loop 22, irun-1 is the index of the last point
c     for the spline segment being constructed --

c     Set the starting index for the next spline segment --
         is=irun
         if (jump .and. (irun .le. npointg(j))) then
            ijump(j)=irun-1
            if (j .eq. jsptrx(region)) then
               is=irun+3        # to allow inclusion of data at is-1
                                # and is-2 for next spline segment
            else
               is=irun+2        # to allow inclusion of data at is-1
                                # and is-2 for next spline segment
            endif
         endif

c     Temporary definition of end points for spline data
      istart1 = max(1,ibeg-2)   # extra data points at ibeg-1 and ibeg-2 are
                                # included (for overlap with previous segment)
				# unless this is the first segment of the contour 
      if (jump .or. (irun .gt. npointg(j))) then
         iend1 = irun - 1       # jump occurs at irun-1
                                # no additional data points are included
      else
         dy1=ycurveg(irun+1,j)-ycurveg(irun,j)
         if (abs(dy1) .gt. dy) then
            iend1 = irun        # jump occurs at irun
                                # data point at irun is included
         else
            iend1 = irun + 1    # data is continuous
                                # data points at irun and irun+1 are included
         endif
      endif

c     Transform into system icood (xcurveg,ycurveg) ---> (xtrans,ytrans)
      beta = (icood-1)*90.
      call transys(beta,istart1,iend1,j)
      irange = iend1 - istart1 -1
      if (irange .lt. 5) goto 1000  # insufficient range for spline fit

c     Test for monotonicity of end-point data
      if ( (xtrans(1) .gt. xtrans(2)) .or.
     &                    (xtrans(2) .gt. xtrans(3)) ) then
         istart1=istart1+1
         call transys(beta,istart1,iend1,j)
      endif
      if (xtrans(1) .gt. xtrans(2)) then
         istart1=istart1+1
         call transys(beta,istart1,iend1,j)
      endif
      irange = iend1 - istart1 + 1
      if ( (xtrans(irange) .lt. xtrans(irange-1)) .or.
     &           (xtrans(irange-1) .lt. xtrans(irange-2)) ) then
         iend1 = iend1-1
         irange = irange-1
         call transys(beta,istart1,iend1,j)
      endif
      if (xtrans(irange) .lt. xtrans(irange-1)) then
         iend1 = iend1-1
         irange = irange-1
         call transys(beta,istart1,iend1,j)
      endif
      if (irange .lt. 5) goto 1000  # insufficient range for spline fit

c     Definitions for this spline segment --
       iseg=iseg+1
       if (iseg .gt. mseg) then
          write (STDOUT,901) mseg,j
  901     format ("more than mseg = ",i2," segments on",
     &            " flux contour j = ",i2)
          write (STDOUT,*) "Try changing the value of alpha1 or"
          write (STDOUT,*) "increase psi0sep or increase psi0min1."
          call xerrab("")
       endif
       isys(iseg,j)=icood
       istartg(iseg,j)=istart1
       iendg(iseg,j)=iend1
       do 10 i=ibeg,irun-1
          isegment(i,j)=iseg
  10   continue
       m(iseg,j)=iendg(iseg,j)-istartg(iseg,j)+1

c     Call SLATEC subroutine FC to do least squares fit --
      ndata=m(iseg,j)
      do i=1,ndata
         xdatag(i)=xtrans(i)
         ydatag(i)=ytrans(i)
      enddo
      do i=1,nord			# end knots
         bkpt(i)=xtrans(1)
      enddo
c     Define interior knots --
      nkint=0
      do i=3,m(iseg,j)-2,2		# Note stride=2 for definition of knots
					# Also note that 2 data points on each
					# end are excluded as interior knots
         nkint=nkint+1
         bkpt(nkint+nord)=xtrans(i)
      enddo
      do i=nkint+nord+1,nkint+2*nord	# end knots
         bkpt(i)=xtrans(m(iseg,j))
      enddo
      nbkpt=nkint+2*nord		# total number of knots
      mode=1
      iwsla(1)=nwdim
      iwsla(2)=niwdim
      call fch(ndata,xdatag,ydatag,sddata,nord,nbkpt,bkpt,
     &        nconst,xconst,yconst,nderiv,mode,coeff,wsla,iwsla)
      if (mode .ne. 0) then
         call remark("  *** subroutine codsys ***")
         call remark("error from SLATEC routine FC")
         call xerrab("")
      endif
c     Store the spline data --
      ncap7(iseg,j)=nbkpt
      do k=1,nbkpt
         xknts(k,iseg,j)=bkpt(k)
         splcoef(k,iseg,j)=coeff(k)
      enddo

 1000 continue

      return
      end

c     ------------------------------------------------------------------

      subroutine evalspln(iseg,j,x,y)
      implicit none
Use(Dimflxgrd)	#,jdim,npts
Use(Dimensions)
Use(Spline)
      integer iseg, j
      real x, y(4)
      integer k, ideriv, inbv, iflagv, ncoef
      real workv(12)
      real B1VAhL
      external B1VAhL

      external xerrab

c     evaluates the (iseg,j) spline at x :
c     returns array y(k) = k-1 th derivative for k=1,4

      ncoef = ncap7(iseg,j)-4
      inbv = 1
      iflagv = 0
      do k=1,4
         ideriv = k-1
         y(k) = B1VAhL (x, ideriv, xknts(1:npts,iseg,j), ncoef, 4,
     &                  splcoef(1:npts,iseg,j), inbv, workv, iflagv)
      enddo
      if (iflagv .ne. 0) then
         write (STDOUT,901) iseg,j
  901    format ("spline evaluation error on segment ",i2," of flux",
     &           " contour ",i2)
         write (STDOUT,902) x
  902    format ("   x = ",1p1e14.6," in rotated coordinate system")
         call xerrab("")
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine extend
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)
Use(Inmesh)
Use(Linkco)

c     local variables --
      integer region,j,nex

c     procedures --
      external remark,xerrab

c     This subroutine ensures that each contour extends to (or below)
c     a lower boundary at y = yextend by linearly extrapolating if
c     necessary.

c     loop over all regions --
      do 800 region=1,noregs

      do 10 j=jmin(region),jmax(region)
      nex=npointg(j)-1
  11  continue
      nex=nex+1
      if (nex .ge. npts) then
         call remark("***** error in subroutine extend")
         call remark("***** number of data points exceeds npts")
         call xerrab("")
      endif
      npointg(j)=nex
      ycurveg(nex+1,j)=2.*ycurveg(nex,j)-ycurveg(nex-1,j)
      xcurveg(nex+1,j)=2.*xcurveg(nex,j)-xcurveg(nex-1,j)
      if ((ycurveg(nex+1,j) .gt. yextend) .and.
     &    (xcurveg(nex+1,j) .gt. xlbnd) .and.
     &    (xcurveg(nex+1,j) .lt. xubnd)) goto 11
  10  continue

  800 continue          # end of loop over regions

      return
      end

c     ------------------------------------------------------------------

      subroutine findalph(nsys,iseg,j,xob,yob,alphab)
      implicit none
Use(Dimflxgrd)	#jdim,npts
Use(Dimensions)
Use(Transfm)
Use(Spline)
      integer nsys, iseg, j
      real xob, yob, alphab
      real alphaf,sinalpf,cosalpf,tanalp,alpha
      real xo, yo(4)

c     this subroutine finds the rotation angle of the moving system
c     relative to the basic system, given the coordinates of the
c     point (on the spline) in the basic system.

c     find the rotation angle of the fixed system in which the
c     spline was fitted.
      alphaf=alphasys(nsys)

c     transform the x of the point into the fixed system.
      sinalpf=sin(alphaf)
      cosalpf=cos(alphaf)
      xo=xob*cosalpf+yob*sinalpf

c     find the slope of the spline at xo.
c     check for possible out-of-range error condition:
         if (xo .lt. xknts(1,iseg,j)) then
            call remark("*** error from s.r. findalph")
            write (STDOUT,901) iseg,j,xob,yob
            call xerrab("")
         elseif (xo .gt. xknts(ncap7(iseg,j),iseg,j)) then
            call remark("*** error from s.r. findalph")
            write (STDOUT,901) iseg,j,xob,yob
            call xerrab("")
         endif
 901     format("iseg,j,xob,yob = ",i3,i3,f11.4,f11.4)

      call evalspln(iseg,j,xo,yo)
      tanalp=yo(2)

c     transform back to the moving system.
      alpha=atan(tanalp)
      alphab=alpha+alphaf

      return
      end

c     ------------------------------------------------------------     c

      subroutine intersect (x,y,i1,i2,xfp,yfp,thetafp,xc,yc,ic)
      implicit none
      integer i1,i2,ic
      real x(*),y(*),xfp,yfp,thetafp,xc,yc

c     Find the intersection of the curve (x(i),y(i)) i=i1,i2
c     with the straight line y = yfp + (x-xfp)*tan(thetafp);
c     Return the intersection point (xc,yc) and the node index ic
c     such that the intersection point lies between nodes ic and ic+1.

      integer i
      real bignum,tanfp,tanic
      data bignum /1.0e+6/
      tanfp=min(tan(thetafp),+bignum)
      tanfp=max(tanfp,-bignum)

      ic=0		# return this 'error' flag value if no intersection

      do i=i1,i2-1	# search forward from x-point; should also consider
			# backward search from end of original mesh ?
         if ( (y(i)-yfp-(x(i)-xfp)*tanfp)*\
                  (y(i+1)-yfp-(x(i+1)-xfp)*tanfp) .le. 0) then
c     line falls between point i and point i+1
            ic=i
            tanic=(y(i+1)-y(i))/(x(i+1)-x(i))
            xc=((y(i)-x(i)*tanic)-(yfp-xfp*tanfp))/(tanfp-tanic)
            yc=yfp+(xc-xfp)*tanfp
            break
         endif
      enddo

      return
      end

c     ------------------------------------------------------------     c

      subroutine intrcept(iseg,j,xo,yo,alpha,xm,sguess,s,isegalt)
      implicit none
      integer iseg, j, noiter, iter, isegalt
      real xo, yo, alpha, xm, sguess, s, sinalp, cosalp, xf, yf
      real x1, y1, slope, dx, xlimit, xfold, sfold, s1, xfnew
      real ss(4),sfnew(4)
Use(Dimflxgrd)	#jdim,npts
Use(Dimensions)
Use(Spline)

c     *** this routine finds the intercept of the spline and the straight
c     *** line through the point (xm,0) parallel to the y axis of the
c     *** moving system.

ccc   If a potential out-of-range spline evaluation error is encountered,
ccc   then isegalt suggests an alternative value for iseg to avoid the error.

c
c     *** transform xm and sguess into the fixed system (of the spline).
         sinalp=sin(alpha)
         cosalp=cos(alpha)
c
         xf=xo+xm*cosalp-sguess*sinalp
         yf=yo+xm*sinalp+sguess*cosalp
c     *** transform (xm,0) into the fixed system.
         x1=xo+xm*cosalp
         y1=yo+xm*sinalp
c
c     ***
      if (abs(alpha).gt.0.001) then
         slope=-1.0/tan(alpha)
         dx=0.01*abs(sguess)
         xlimit=0.1*dx
c
c     *** obtain the solution of the equation s-y1-slope*(x-x1)=0.
         xfold=xf
         sfold=yf
         noiter=10
c
      isegalt=iseg
      do 10 iter=1,noiter

c     find the derivative of s at xfold
ccc   begin error check
         if (xfold .lt. xknts(1,iseg,j)) then
            isegalt=iseg-1	# switch to previous spline segment
            goto 12
         elseif (xfold .gt. xknts(ncap7(iseg,j),iseg,j)) then
            isegalt=iseg+1	# switch to next spline segment
            goto 12
         endif
ccc   end error check
         call evalspln(iseg,j,xfold,ss)
         s1=ss(2)

         xfnew=xfold-(sfold-y1-slope*(xfold-x1))/(s1-slope)
ccc   begin error check
         if (xfnew .lt. xknts(1,iseg,j)) then
            isegalt=iseg-1	# switch to previous spline segment
            goto 12
         elseif (xfnew .gt. xknts(ncap7(iseg,j),iseg,j)) then
            isegalt=iseg+1	# switch to next spline segment
            goto 12
         endif
ccc   end error check
         call evalspln(iseg,j,xfnew,sfnew)

         if(abs(xfnew-xfold).lt.xlimit) goto 11
            xfold=xfnew
            sfold=sfnew(1)
c
   10 continue
   11 continue
c
         else
            xfnew=x1
ccc   begin error check
            if (xfnew .lt. xknts(1,iseg,j)) then
               isegalt=iseg-1	# switch to previous spline segment
               goto 12
            elseif (xfnew .gt. xknts(ncap7(iseg,j),iseg,j)) then
               isegalt=iseg+1	# switch to next spline segment
               goto 12
            endif
ccc   end error check
            call evalspln(iseg,j,xfnew,sfnew)
         endif
c
c     *** transform the answer into the moving system.
         s=-(xfnew-xo)*sinalp+(sfnew(1)-yo)*cosalp
c
   12 continue
      return
      end

c     ------------------------------------------------------------------

      subroutine meshfin
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)
Use(Linkco)
Use(Inmesh)
Use(Mmod)
Use(Share)   # ismmon,ishalfm
      integer ix1,ix2,i,irange,i1,i2,region
      integer j,jrange,j1,j2,jj, rbegin
      real xavg,yavg,cmx,cmy

      external meshmod1,meshmod2,meshmod3

      if ( (geometry .eq. "dnbot") .or. (geometry .eq. "dnull") .or.
     .    (geometry == "isoleg") ) then
         if (istpnew==0) then	# 1st seed point is NOT user-specified
c     Move the first angle-like surface to the exact magnetic midplane
            do j=jmin(1),jmax(2)
               cmeshy(1,j) = zmagx
            enddo
         endif
      endif

c     Corrections at cut interfaces where two halves of the mesh
c     are joined together:
      if (ishalfm==0) then

c     Correct the curve going from the x-point to the center.
      ix1=ixpoint(1,1)-1
      ix2=ixpoint(1,2)-1
      irange=ixpoint(2,1)-ixpoint(1,1)+1
      jrange=jmax(1)-jsptrx(1)
      do i=1,irange
         i1=ix1+i
         i2=ix2+i
         do j=0,jrange
            j1=jsptrx(1)+j
            j2=jsptrx(2)-j
            xavg=0.5*(cmeshx(i1,j1)+cmeshx(i2,j2))
            yavg=0.5*(cmeshy(i1,j1)+cmeshy(i2,j2))
            cmeshx(i1,j1)=xavg
            cmeshy(i1,j1)=yavg
            cmeshx(i2,j2)=xavg
            cmeshy(i2,j2)=yavg
         enddo
      enddo

c     Correct the curve going from the x-point to the private region.
      jrange=jmax(1)-jsptrx(1)
      i1=ixpoint(3,1)
      i2=ixpoint(3,2)
      do j=0,jrange
         j1=jsptrx(1)+j
         j2=jsptrx(2)-j
         xavg=0.5*(cmeshx(i1,j1)+cmeshx(i2,j2))
         yavg=0.5*(cmeshy(i1,j1)+cmeshy(i2,j2))
         cmeshx(i1,j1)=xavg
         cmeshy(i1,j1)=yavg
         cmeshx(i2,j2)=xavg
         cmeshy(i2,j2)=yavg
      enddo

c     Correct the curve through the top of the mesh.
      if ( ((geometry .eq. "snull").or.(geometry .eq. "uppersn"))
     .                      .and. (islimon==0) ) then
         do jj = jmin(1),jmax(1)
            cmx = 0.5 * (cmeshx(1,jj)+cmeshx(1,jmax(2)-jj+jmin(1)))
            cmy = 0.5 * (cmeshy(1,jj)+cmeshy(1,jmax(2)-jj+jmin(1)))
            cmeshx(1,jj) = cmx
            cmeshx(1,jmax(2)-jj+jmin(1)) = cmx
            cmeshy(1,jj) = cmy
            cmeshy(1,jmax(2)-jj+jmin(1)) = cmy
         enddo
      endif

      endif   # end corrections on cut interfaces

# Options for creating non-orthogonal mesh:
      if (ismmon .ne. 0) then		# save copy of orthogonal mesh
         call s2copy(idim,jdim,cmeshx,1,idim,cmeshx0,1,idim)
         call s2copy(idim,jdim,cmeshy,1,idim,cmeshy0,1,idim)
      endif
      if (ismmon .eq. 1) then
c..      set flag to pass to nphygeo in package bbb via the Share group
         isnonog = 1
c        define upstream and plate arrays, if necessary
         if (istream .eq. 0) call getu	# default upstream definition
         if (iplate  .eq. 0) call getp	# default orthogonal definition
         do region=1,noregs
            call meshmod1(region)
         enddo
c        smooth angle-like surfaces, if desired
         call smoother
      elseif (ismmon .eq. 2) then
c..      set flag to pass to nphygeo in package bbb via the Share group
         isnonog = 1
c        define upstream and plate arrays, if necessary
         if (istream .eq. 0) call getu	# default upstream definition
         if (iplate  .eq. 0) call getp	# default orthogonal definition
         do region=1,noregs
            call meshmod2(region)
         enddo
c        smooth angle-like surfaces, if desired
         call smoother
      elseif (ismmon .eq. 3) then
c..      set flag to pass to nphygeo in package bbb via the Share group
         isnonog = 1
         call gett			# top-of-mesh definition
         call getd			# downstream definition
c        define upstream and plate arrays, if necessary
         if (istream .eq. 0) call getu	# default upstream definition
         if (iplate  .eq. 0) call getp	# default orthogonal definition
         do region=1,noregs
            call meshmod3(region)
         enddo
c        smooth angle-like surfaces, if desired
         call smoother
      endif

# Option for including limiter surface:
      if ((islimon .ne. 0) .and. ((geometry .eq. 'snull')
     .        .or. (geometry .eq. 'uppersn'))) then
         call s2copy(idim,jdim,cmeshx,1,idim,cmeshx0,1,idim)
         call s2copy(idim,jdim,cmeshy,1,idim,cmeshy0,1,idim)
         isnonog = 1   # non-orthog mesh: pass this flag to bbb.nphygeo
         call getlim   # define limiter interface with in/outboard regions
         do region=1,noregs
            call meshlim(region)
         enddo
c        do smoothing on angle-like surfaces ????
         call setlimindex      # before or after smoothing ???
      endif

# Assign values to cmesh for the "axis" flux surface :
#
      j = jaxis
#
# above the x-point (magnetic axis) :
      i1=1
      i2=ixpoint(2,1)
      do i=i1,i2
            cmeshx(i,j) = rmagx
            cmeshy(i,j) = zmagx
      enddo
#
# below the x-point (conjugate axis) :
      i1=ixpoint(3,1)
      i2=max(ilmax(1),ilmax(2))
      do i=i1,i2
         cmeshx(i,j) = 0.5 * ( cmeshx(ilmax(1),jmax(1))
     &                      +  cmeshx(ilmax(2),jmin(2)) )
         cmeshy(i,j) = 0.5 * ( cmeshy(ilmax(1),jmax(1))
     &                      +  cmeshy(ilmax(2),jmin(2)) )
      enddo

# assign values to cmesh when ilmax(1) and ilmax(2) are different by
# creating zero-volume cells at the plates :
	if (ilmax(1) .lt. ilmax(2)) then
	   do i=ilmax(1)+1,ilmax(2)
	      do j=jmin(1),jmax(1)
	         cmeshx(i,j) = cmeshx(ilmax(1),j)
	         cmeshy(i,j) = cmeshy(ilmax(1),j)
	      enddo
	   enddo
	elseif (ilmax(2) .lt. ilmax(1)) then
	   do i=ilmax(2)+1,ilmax(1)
	      do j=jmin(2),jmax(2)
	         cmeshx(i,j) = cmeshx(ilmax(2),j)
	         cmeshy(i,j) = cmeshy(ilmax(2),j)
	      enddo
	   enddo
	endif

      if (mincorecellwidth .eq. 1) then
      rbegin = 1
      if (ishalfm==1) rbegin = 2
      do region = rbegin, noregs
        call corexpointfix(region)
      end do
      endif

      return
      end

c----------------------------------------------------------------------c

      subroutine meshgen (region)
      implicit none
      integer region
Use(Share)        #isfrc
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)
Use(Dimensions)
Use(Curves)
Use(Linkco)
Use(Transfm)
Use(Inmesh)

c     local variables --
      integer i,ii,iseg,isptrx,j,jj,j1,j2,nsys
      real alphab,alphab1,alphab2,alphabsx,distance,distmin,
     &     xob,xosptrx,yob,yosptrx

c     procedures ---
      external findalph, xerrab, orthogrd, orthogx, remark

c     initialization --
      data alphab /0./
      data alphab1 /0./
      data alphab2 /0./
      data alphabsx /0./
      data distance /0./
      data distmin /0./
      data xob /0./
      data yob /0./
      data xosptrx /0./
      data yosptrx /0./

c     This subroutine generates a mesh by finding an orthogonal set
c     of curves to (xcurveg,ycurveg)

         jj=jsptrx(region)

c     loop over the seed points --
      do 20 i=1,ilmax(region)
         xosptrx=cmeshx(i,jj)
         yosptrx=cmeshy(i,jj)


      if((i.lt.ixpoint(1,region)).or.(i.gt.ixpoint(3,region))) then
c     *** find which segment this point is on
         distmin=(xubnd-xlbnd)+(yubnd-ylbnd)
         do 30 ii=1,npointg(jj)
            distance = sqrt((xcurveg(ii,jj)-xosptrx)**2
     .                     +(ycurveg(ii,jj)-yosptrx)**2)
            if(distance.gt.distmin) goto 30
            distmin=distance
            isptrx=ii
   30    continue
c     *** we now know the nearest point on the curve.
c     *** find the segment and the corresponding system.
         iseg=isegment(isptrx,jj)
         if (iseg .eq. 0) then
            call remark ("*** error in subroutine meshgen")
            write (STDOUT,899) isptrx,jj
  899       format ("data point i=",i3," on contour j=",i3," is not",
     &              " assigned to a spline segment")
            call xerrab("")
         endif
         nsys=isys(iseg,jj)
c     *** find the rotation angle of the moving system.
         call findalph(nsys,iseg,jj,xosptrx,yosptrx,alphabsx)
         alphab=alphabsx
      endif


c     *** next go in the direction of smaller j and find the nearest
c     *** point on the next curve satisfying the orthogonality
c     *** condition
c     Finding the first point off the separatrix requires special 
c     treatment when starting from the x-point because the slope there 
c     is not uniquely defined => use s.r. orthogx vs orthogrd
      if (jj .gt. jmin(region)) then
        xob=xosptrx
        yob=yosptrx
        j1=jj
        j2=jj-1
        j=j2
        if (i.lt.ixpoint(1,region)) then
           call orthogrd(ityp(1,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
        elseif (i.eq.ixpoint(1,region)) then
           call orthogx(ityp(1,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
           alphab1=alphab
        elseif (i.eq.ixpoint(2,region)) then
           cmeshx(i,j)=cmeshx(i-1,j)
           cmeshy(i,j)=cmeshy(i-1,j)
           xob=cmeshx(i,j)
           yob=cmeshy(i,j)
           alphab=alphab1
ccc Modification for field-reversed configuration
ccc        elseif (i.eq.ixpoint(3,region)) then
        elseif (i.eq.ixpoint(3,region) .and. isfrc==0) then
           call orthogx(ityp(3,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
ccc Modification for field-reversed configuration
ccc	elseif (i.gt.ixpoint(3,region)) then
	elseif (i.gt.ixpoint(3,region) .and. isfrc==0) then
           call orthogrd(ityp(3,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
        endif

        do 50 j=jj-2,jmin(region),-1
           j1=j+1
           j2=j
	   if (i.le.ixpoint(2,region)) then
              call orthogrd(ityp(1,region),i,j1,j2,xob,yob,alphab)
              cmeshx(i,j)=xob
              cmeshy(i,j)=yob
ccc Modification for field-reversed configuration
ccc	   elseif (i.ge.ixpoint(3,region)) then
	   elseif (i.ge.ixpoint(3,region) .and. isfrc==0) then
              call orthogrd(ityp(3,region),i,j1,j2,xob,yob,alphab)
              cmeshx(i,j)=xob
              cmeshy(i,j)=yob
	   endif
   50   continue
      endif


c     *** repeat the process for j greater than jsptrx.
      if (jj .lt. jmax(region)) then
        xob=xosptrx
        yob=yosptrx
        alphab=alphabsx
        j1=jj
        j2=jj+1
        j=j2
        if (i.lt.ixpoint(1,region)) then
           call orthogrd(ityp(4,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
        elseif (i.eq.ixpoint(1,region)) then
           call orthogx(ityp(4,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
           alphab2=alphab
        elseif (i.eq.ixpoint(2,region)) then
           cmeshx(i,j)=cmeshx(i-1,j)
           cmeshy(i,j)=cmeshy(i-1,j)
           xob=cmeshx(i,j)
           yob=cmeshy(i,j)
           alphab=alphab2
        elseif (i.eq.ixpoint(3,region)) then
           call orthogx(ityp(6,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
	elseif (i.gt.ixpoint(3,region)) then
           call orthogrd(ityp(6,region),i,j1,j2,xob,yob,alphab)
           cmeshx(i,j)=xob
           cmeshy(i,j)=yob
        endif

        do 60 j=jj+2,jmax(region),+1
           j1=j-1
           j2=j
	   if (i.le.ixpoint(2,region)) then
              call orthogrd(ityp(4,region),i,j1,j2,xob,yob,alphab)
              cmeshx(i,j)=xob
              cmeshy(i,j)=yob
	   elseif (i.ge.ixpoint(3,region)) then
              call orthogrd(ityp(6,region),i,j1,j2,xob,yob,alphab)
              cmeshx(i,j)=xob
              cmeshy(i,j)=yob
	   endif
   60   continue
      endif


c     *** close the loop along the seed points.
   20 continue

      return
      end

c----------------------------------------------------------------------c

      subroutine meshmod1 (region)
      implicit none
      integer region
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)
Use(Dimensions)
Use(Curves)
Use(Inmesh)
Use(Linkco)
Use(Mmod)
Use(Transfm)
      external intersect2, remark, xerrab

c     Local variables --
      integer i,istart,iend,ii,ii1,ii2,j,j1,j2,k,ierr
      integer i1crv,i2crv,i1msh,i2msh,imu,kmu,icu,kcu,icp,kcp
      real xmu,ymu,xcu,ycu,xcp,ycp
      real dsmp0,dscp,frac,eps

c     This subroutine modifies the original orthogonal mesh by
c     expanding or compressing the grid along each flux surface between
c     the divertor plate and some upstream reference surface.

      data eps/1.0e-06/		# extrapolation parameter

      j1=jmin(region)
      j2=jmax(region)
      iend=ilmax(region)

c     Define the upstream reference surface that stays fixed:
      if (region .eq. 1) then
         nupstream=nupstream1
         if (nupstream .gt. nptmp) then
            call xerrab("dimension of r,zupstream1 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream1(k)
            zupstream(k)=zupstream1(k)
         enddo
      elseif (region .eq. 2) then
         nupstream=nupstream2
         if (nupstream .gt. nptmp) then
            call xerrab("dimension of r,zupstream2 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream2(k)
            zupstream(k)=zupstream2(k)
         enddo
      endif

c     Define the surface that coincides with the new divertor plate:
      if (region .eq. 1) then
         nplate=nplate1
         if (nplate .gt. nptmp) then
            call xerrab("dimension of r,zplate1 exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate1(k)
            zplate(k)=zplate1(k)
         enddo
      elseif (region .eq. 2) then
         nplate=nplate2
         if (nplate .gt. nptmp) then
            call xerrab("dimension of r,zplate2 exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate2(k)
            zplate(k)=zplate2(k)
         enddo
      endif

c     Define the surface that coincides with the original divertor plate:
      nplate0=j2-j1+1
         if (nplate0 .gt. nptmp) then
            call xerrab("dimension of r,zplate0 exceeds nptmp")
         endif
      do k=j1,j2
         rplate0(k-j1+1)=cmeshx(iend,k)
         zplate0(k-j1+1)=cmeshy(iend,k)
      enddo
c     extrapolate end points slightly beyond the the original plate
      rplate0(1)=2*rplate0(1)-rplate0(2)
      zplate0(1)=2*zplate0(1)-zplate0(2)
      rplate0(nplate0)=2*rplate0(nplate0)-rplate0(nplate0-1)
      zplate0(nplate0)=2*zplate0(nplate0)-zplate0(nplate0-1)

c----------------------------------------------------------------------c
c     Loop over all flux surfaces:
      do j=j1,j2

c     Define some temporary variables for (xcurveg,ycurveg) surface data
         i1crv=ijump(j)+1	# start index
         i2crv=npointg(j)	# end index
         if (j .eq. jsptrx(region)) i1crv=1
         do i=i1crv,i2crv
            xcrv(i)=xcurveg(i,j)
            ycrv(i)=ycurveg(i,j)
         enddo

c     Define some temporary variables for (cmeshx,y) surface data
         i1msh=1			# start at top of SOL
         i2msh=ilmax(region)		# end index
         if ( (j .gt. jsptrx(1)) .and. (j .lt. jsptrx(2)) ) then
            i1msh=ixpoint(3,region)	# start at p.f. cut
         endif
         do i=i1msh,i2msh
            xmsh(i)=cmeshx(i,j)
            ymsh(i)=cmeshy(i,j)
         enddo
c     extrapolate upstream slightly from the initial point
         xmsh(i1msh)=xmsh(i1msh)+eps*(xmsh(i1msh)-xmsh(i1msh+1))
         ymsh(i1msh)=ymsh(i1msh)+eps*(ymsh(i1msh)-ymsh(i1msh+1))

c     Intersection of (cmeshx,y) with upstream reference surface
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xmsh(i1msh:i2msh),ymsh(i1msh:i2msh),i1msh,i2msh,
     &                   xmu,ymu,kmu,imu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,885) j
 885        format ("no upstream intersection for cmeshx,y on j=",i3)
            call xerrab("")
         endif

c        ensure that separatrix points are not modified above x-point:
         if ( (j .eq. jsptrx(region)) .and. 
     &        (imu .lt. ixpoint(3,region)) ) then
            imu=ixpoint(3,region)
            xmu=rseps
            ymu=zseps
         endif

         istart=imu+1	# index of first point below reference surface
         iend=i2msh	# index of last point below reference surface

c     Distance from upstream reference to each (cmeshx,y) point
         dsm(istart-1)=0.
         dsm(istart)=sqrt( (xmsh(istart)-xmu)**2 +
     &                        (ymsh(istart)-ymu)**2 )
         do i=istart+1,iend
            dsm(i)=dsm(i-1) +
     &                   sqrt( (xmsh(i)-xmsh(i-1))**2 +
     &                         (ymsh(i)-ymsh(i-1))**2 )
         enddo
c     Total distance from upstream reference to old divertor plate
c     along the cmeshx,y representation of the flux surface
         dsmp0=dsm(iend)
c     Normalize
         do i=istart,iend
            dsm(i)=dsm(i)/dsmp0
         enddo

c     Intersection of (x,ycurveg) with upstream reference surface
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcu,ycu,kcu,icu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j
 886        format ("no upstream intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif

c        ensure that separatrix points are not modified above x-point:
         if ( (j .eq. jsptrx(region)) .and. 
     &        (icu .lt. ijump(j)) ) then
            icu=ijump(j)
            xcu=rseps
            ycu=zseps
         endif

         ii1=icu+1	# index of first point below reference surface
c     Re-define the data point at index ii1-1 to lie exactly on the
c     reference surface:
         xcrv(ii1-1)=xcu
         ycrv(ii1-1)=ycu

c     Intersection of (x,ycurveg) with new divertor plate
         call intersect2(rplate,zplate,1,nplate,
     &                    xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                    xcp,ycp,kcp,icp,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,887) j
 887        format ("no new plate intersection for x,ycurveg on j=",i3)
            write (STDOUT,889)
 889        format ("*** You may have to extend the plate end points")
            call xerrab("")
         endif
         ii2=icp+1	# index of first point below new divertor plate

c     Distance from upstream reference to each (x,ycurveg) point
         dsc(ii1-1)=0.
         do ii=ii1,i2crv
            dsc(ii)=dsc(ii-1) +
     &                   sqrt( (xcrv(ii)-xcrv(ii-1))**2 +
     &                         (ycrv(ii)-ycrv(ii-1))**2 )
         enddo

c     Total distance from upstream reference to new divertor plate
         dscp=dsc(ii2)-sqrt( (xcrv(ii2)-xcp)**2 +
     &                               (ycrv(ii2)-ycp)**2 )


c     Compute new arc lengths --
         do i=istart,iend
            dsmesh(i)=dsm(i)*dscp
         enddo

c     Modified cmeshx,cmeshy values --
            do i=istart,iend
               do ii=ii1,ii2
                  if (dsmesh(i) .le. dsc(ii)) then
c     modified point i lies between (x,ycurveg) points ii-1 and ii
                     frac=(dsmesh(i)-dsc(ii-1))/
     &                    (dsc(ii)-dsc(ii-1))
                     cmeshx(i,j)=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                     cmeshy(i,j)=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                     break
                  endif
               enddo
            enddo

      enddo	# end of loop over all j flux surfaces
c----------------------------------------------------------------------c

      return
      end

c----------------------------------------------------------------------c

      subroutine meshmod2 (region)
      implicit none
      integer region
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)	# nix,idim
Use(Curves)	# xcurveg,ycurveg
Use(Inmesh)	# ilmax
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# 
Use(Transfm)	# ijump
      external intersect2, remark, xerrab

c     local variables --
      integer i,ii,istart,iend,j,jj,k
      integer i1crv,i2crv,kcu,icu,ierr,ii1,kcp,icp,ii2
      real xcu,ycu,xcp,ycp
      real dstotal,dscp,frac

c     This subroutine places mesh points on flux surfaces
c     using the same normalized distribution as the seed points
c     on the separatrix.

c----------------------------------------------------------------------c
c     Define the upstream reference surface:
      if (region .eq. 1) then
         nupstream=nupstream1
         if (nupstream .gt. nptmp) then
            call xerrab("meshmod2: dim'n of r,zupstream1 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream1(k)
            zupstream(k)=zupstream1(k)
         enddo
      elseif (region .eq. 2) then
         nupstream=nupstream2
         if (nupstream .gt. nptmp) then
            call xerrab("meshmod2: dim'n of r,zupstream2 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream2(k)
            zupstream(k)=zupstream2(k)
         enddo
      endif

c----------------------------------------------------------------------c
c     Define the divertor plate surface:
      if (region .eq. 1) then
         nplate=nplate1
         if (nplate .gt. nptmp) then
            call xerrab("meshmod2: dimension of r,zplate exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate1(k)
            zplate(k)=zplate1(k)
         enddo
      elseif (region .eq. 2) then
         nplate=nplate2
         if (nplate .gt. nptmp) then
            call xerrab("meshmod2: dimension of r,zplate exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate2(k)
            zplate(k)=zplate2(k)
         enddo
      endif

c----------------------------------------------------------------------c
c     Loop over all flux surfaces:
      do j=jmin(region),jmax(region)

c     Define some temporary variables for (xcurveg,ycurveg) surface data
         i1crv=ijump(j)+1	# start index
         i2crv=npointg(j)	# end index
         if (j .eq. jsptrx(region)) i1crv=1
         do i=i1crv,i2crv
            xcrv(i)=xcurveg(i,j)
            ycrv(i)=ycurveg(i,j)
         enddo

c     Intersection of reference surface with x,ycurveg
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcu,ycu,kcu,icu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j
 886        format 
     & ("meshmod2: no reference intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
         ii1=icu+1	# index of first point below reference surface
c     Re-define the data point at index ii1-1 to lie exactly on the
c     reference surface:
         xcrv(ii1-1)=xcu
         ycrv(ii1-1)=ycu

c     Intersection of divertor plate surface with x,ycurveg)
         call intersect2(rplate,zplate,1,nplate,
     &                    xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                    xcp,ycp,kcp,icp,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,887) j
 887        format 
     & ("meshmod2: no plate intersection for x,ycurveg on j=",i3)
            write (STDOUT,889)
 889        format ("*** You may have to extend the plate end points")
            call xerrab("")
         endif
         ii2=icp+1	# index of first point below divertor plate

c     Distance from reference surface to each (x,ycurveg) point
         dsc(ii1-1)=0.
         do ii=ii1,i2crv
            dsc(ii)=dsc(ii-1) +
     &                   sqrt( (xcrv(ii)-xcrv(ii-1))**2 +
     &                         (ycrv(ii)-ycrv(ii-1))**2 )
         enddo

c     Total distance from reference surface to divertor plate
         dscp=dsc(ii2)-sqrt( (xcrv(ii2)-xcp)**2 +
     &                               (ycrv(ii2)-ycp)**2 )

c     Normalized distance along the separatrix:

ccc   NOTE: The following assumes that the upstream reference surface 
ccc         coincides with the i=1 surface in the SOL and the 
ccc         i=ixpoint(3,region) in the private flux region.  This 
ccc         may NOT be consistent with the upstream data that the 
ccc         user supplies !!!

         if ( (j .ge. jsptrx(1)) .and. (j .le. jsptrx(2)) ) then
            istart=ixpoint(3,region)	# for psi<=1 flux surfaces
         else
            istart=1			# for psi>1 flux surfaces
         endif
         iend=ilmax(region)
         jj=jsptrx(region)
         dsnorm(istart)=0.
         do i=istart+1,iend
            dsnorm(i)=dsnorm(i-1) +
     &                sqrt( (cmeshx(i,jj)-cmeshx(i-1,jj))**2 +
     &                      (cmeshy(i,jj)-cmeshy(i-1,jj))**2 )
         enddo
         dstotal=dsnorm(iend)
         do i=istart,iend
            dsnorm(i)=dsnorm(i)/dstotal
         enddo

c     Find cmeshx,y by interpolation along x,ycurveg
         do i=istart+1,iend
            dsmesh(i)=dsnorm(i)*dscp
            do ii=ii1,ii2
               if (dsmesh(i) .le. dsc(ii)) then
c           modified point i lies between (x,ycurveg) points ii-1 and ii
                  frac=(dsmesh(i)-dsc(ii-1))/
     &                 (dsc(ii)-dsc(ii-1))
                  cmeshx(i,j)=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                  cmeshy(i,j)=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                  break
               endif
            enddo
         enddo

      enddo		# end of loop over flux surfaces

      return
      end

c----------------------------------------------------------------------c

      subroutine meshmod3 (region)
      implicit none
      integer region
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)	# rseps,zseps
Use(Dimensions) # nptmp
Use(Curves)     # npointg,xcurveg,ycurveg
Use(Inmesh)     # ilmax
Use(Linkco)     # ixpoint,cmeshx,cmeshy
Use(Mmod)       # ntop,ntop1,ntop2,rtop,rtop1,rtop2,ztop,ztop1,ztop2,
                # nupstream,nupstream1,nupstream2,rupstream,rupstream1,
                # rupstream2,zupstream,zupstream1,zupstream2,ndnstream,
                # ndnstream1,ndnstream2,rdnstream,rdnstream1,rdnstream2,
                # zdnstream,zdnstream1,zdnstream2,nplate,nplate1,
                # nplate2,rplate,rplate1,rplate2,zplate,zplate1,zplate2,
                # nplate0,rplate0,zplate0,cmeshx0,cmeshy0,fuzzm,xmsh,
                # ymsh,dss,dssleg,xcrv,ycrv,dsm,dsc,dmix0,delmax,
                # dsmesh0,dsmesh1,dsmesh2,dsmesh,wtmesh1
Use(Transfm)    # ijump
      external intersect2, remark, xerrab, twixt
      logical twixt

c     Local variables --
      integer i,iend,ii,j,j1,j2,k,ierr,i0,i0leg
      integer i1crv,i2crv,i1msh,i2msh,imu,kmu,icu,kcu,icp,kcp
      integer ict,kct,imt,icp0,kcp0,icd,kcd,imd,im12d
      integer na,ica,kca
      real xmu,ymu,xcu,ycu,xcp,ycp,xct,yct,xmt,ymt,xcp0,ycp0,dscp0
      real dsmu,dsmp0,dscu,dscp,frac,eps,dssp0,dsslegp0,dssu,fracu
      real xcd,ycd,xmd,ymd,dscd,len,dscmi
      real ra(3),za(3),xca,yca

c     This subroutine modifies the original orthogonal mesh by
c     re-distributing the mesh points along each flux surface.
c     It combines features of the ismmon=0,1,2 options via
c     weight factor wtmesh1 and normalized mixing length dmix0.

      data eps/1.0e-06/		# extrapolation parameter

      j1=jmin(region)
      j2=jmax(region)
      iend=ilmax(region)

c     Define the top-of-mesh surface :
      if (region .eq. 1) then
         ntop=ntop1
         if (ntop .gt. nptmp) then
            call xerrab("dimension of r,ztop1 exceeds nptmp")
         endif
         do k=1,ntop
            rtop(k)=rtop1(k)
            ztop(k)=ztop1(k)
         enddo
      elseif (region .eq. 2) then
         ntop=ntop2
         if (ntop .gt. nptmp) then
            call xerrab("dimension of r,ztop2 exceeds nptmp")
         endif
         do k=1,ntop
            rtop(k)=rtop2(k)
            ztop(k)=ztop2(k)
         enddo
      endif

c     Define the upstream reference surface :
      if (region .eq. 1) then
         nupstream=nupstream1
         if (nupstream .gt. nptmp) then
            call xerrab("dimension of r,zupstream1 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream1(k)
            zupstream(k)=zupstream1(k)
         enddo
      elseif (region .eq. 2) then
         nupstream=nupstream2
         if (nupstream .gt. nptmp) then
            call xerrab("dimension of r,zupstream2 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream2(k)
            zupstream(k)=zupstream2(k)
         enddo
      endif

c     Define the downstream reference surface :
      if (region .eq. 1) then
         ndnstream=ndnstream1
         if (ndnstream .gt. nptmp) then
            call xerrab("dimension of r,zdnstream1 exceeds nptmp")
         endif
         do k=1,ndnstream
            rdnstream(k)=rdnstream1(k)
            zdnstream(k)=zdnstream1(k)
         enddo
      elseif (region .eq. 2) then
         ndnstream=ndnstream2
         if (ndnstream .gt. nptmp) then
            call xerrab("dimension of r,zdnstream2 exceeds nptmp")
         endif
         do k=1,ndnstream
            rdnstream(k)=rdnstream2(k)
            zdnstream(k)=zdnstream2(k)
         enddo
      endif

c     Define the new divertor plate surface :
      if (region .eq. 1) then
         nplate=nplate1
         if (nplate .gt. nptmp) then
            call xerrab("dimension of r,zplate1 exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate1(k)
            zplate(k)=zplate1(k)
         enddo
      elseif (region .eq. 2) then
         nplate=nplate2
         if (nplate .gt. nptmp) then
            call xerrab("dimension of r,zplate2 exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate2(k)
            zplate(k)=zplate2(k)
         enddo
      endif

c     Define the old (orthogonal) divertor plate surface :
      nplate0=j2-j1+1
         if (nplate0 .gt. nptmp) then
            call xerrab("dimension of r,zplate0 exceeds nptmp")
         endif
      do k=j1,j2
         rplate0(k-j1+1)=cmeshx0(iend,k)
         zplate0(k-j1+1)=cmeshy0(iend,k)
      enddo
c     extrapolate end points slightly beyond the the original plate
      rplate0(1)=2*rplate0(1)-rplate0(2)
      zplate0(1)=2*zplate0(1)-zplate0(2)
      rplate0(nplate0)=2*rplate0(nplate0)-rplate0(nplate0-1)
      zplate0(nplate0)=2*zplate0(nplate0)-zplate0(nplate0-1)

c----------------------------------------------------------------------c
c     Compute original meshpoint distribution along separatrix
c     for use with ismmon=2 option.

      j=jsptrx(region)

c     Define some temporary variables for (cmeshx,y) surface data
         i1msh=1			# start index
         i2msh=ilmax(region)		# end index
         do i=i1msh,i2msh
            xmsh(i)=cmeshx0(i,j)
            ymsh(i)=cmeshy0(i,j)
         enddo

c     Distance from top-of-mesh to each downstream meshpoint along separatrix
         i0=1				# index of top-of-mesh surface
					# for psi > 1
         dss(i0)=0.
         do i=i0+1,i2msh
            dss(i)=dss(i-1) +
     &                   sqrt( (xmsh(i)-xmsh(i-1))**2 +
     &                         (ymsh(i)-ymsh(i-1))**2 )
         enddo
         i0leg=ixpoint(3,region)	# index of top-of-mesh surface
					# for psi < 1
         dssleg(i0leg)=0.
         do i=i0leg+1,i2msh
            dssleg(i)=dssleg(i-1) +
     &                   sqrt( (xmsh(i)-xmsh(i-1))**2 +
     &                         (ymsh(i)-ymsh(i-1))**2 )
         enddo

c     Distance from top-of-mesh to orthogonal divertor plate
         dssp0=dss(i2msh)
         dsslegp0=dssleg(i2msh)

c----------------------------------------------------------------------c
c     Loop over all flux surfaces:
      do j=j1,j2

c     Define some temporary variables for (xcurveg,ycurveg) surface data
         i1crv=ijump(j)+1	# start index
         i2crv=npointg(j)	# end index
         if (j .eq. jsptrx(region)) i1crv=1
         do i=i1crv,i2crv
            xcrv(i)=xcurveg(i,j)
            ycrv(i)=ycurveg(i,j)
         enddo

c     Define some temporary variables for (cmeshx,y) surface data
         i1msh=1			# start at top of SOL
         i2msh=ilmax(region)		# end index
         if ( (j .gt. jsptrx(1)) .and. (j .lt. jsptrx(2)) ) then
            i1msh=ixpoint(3,region)	# start at p.f. cut
         endif
         do i=i1msh,i2msh
            xmsh(i)=cmeshx0(i,j)
            ymsh(i)=cmeshy0(i,j)
         enddo
c     extrapolate upstream slightly from the initial point
         xmsh(i1msh)=xmsh(i1msh)+eps*(xmsh(i1msh)-xmsh(i1msh+1))
         ymsh(i1msh)=ymsh(i1msh)+eps*(ymsh(i1msh)-ymsh(i1msh+1))

c     Top-of-mesh for cmeshx,y:
         if ( (jsptrx(1) .le. j) .and. (j .le. jsptrx(2)) ) then
            imt=ixpoint(3,region)
         else
            imt=1
         endif
         xmt=xmsh(imt)
         ymt=ymsh(imt)
c        index of first point below top-of-mesh is imt+1

c     Intersection of (cmeshx,y) with upstream reference surface
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xmsh(i1msh:i2msh),ymsh(i1msh:i2msh),i1msh,i2msh,
     &                   xmu,ymu,kmu,imu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,880) j
 880        format ("no upstream intersection for cmeshx,y on j=",i3)
            call xerrab("")
         endif
         if (imu .lt. imt) then
            imu=imt
            xmu=xmt
            ymu=ymt
         endif

c     Downstream reference surface for cmeshx,y:
      imd=ixpoint(3,region)
      xmd=xmsh(imd)
      ymd=ymsh(imd)

c     Distance from top-of-mesh to each meshpoint along cmeshx,y
         dsm(imt)=0.
         do i=imt+1,i2msh
            dsm(i)=dsm(i-1) +
     &                   sqrt( (xmsh(i)-xmsh(i-1))**2 +
     &                         (ymsh(i)-ymsh(i-1))**2 )
         enddo

c     Distance from top-of-mesh to upstream reference surface
         dsmu=dsm(imu)+sqrt( (xmu-xmsh(imu))**2 +
     &                       (ymu-ymsh(imu))**2 )

c     Distance from top-of-mesh to orthogonal divertor plate
         dsmp0=dsm(i2msh)

c----------------------------------------------------------------------c

c     Intersection of (x,ycurveg) with top-of-mesh surface
         call intersect2(rtop,ztop,1,ntop,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xct,yct,kct,ict,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,884) j
 884        format ("no top-of-mesh intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
         if (j .eq. jsptrx(region)) then	# top is at x-point
            ict=ijump(j)+2	# there are 3 identical data points here
            xct=rseps
            yct=zseps
         endif
c        index of first point below top-of-mesh is ict+1

c     Intersection of (x,ycurveg) with upstream reference surface
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcu,ycu,kcu,icu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,885) j
 885        format ("no upstream intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below reference surface is icu+1
c        Ensure that upstream surface lies "below" top-of-mesh surface
         if (icu .lt. ict) then
            icu=ict
            xcu=xct
            ycu=yct
         endif

c     Intersection of (x,ycurveg) with downstream reference surface
         call intersect2(rdnstream,zdnstream,1,ndnstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcd,ycd,kcd,icd,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j
 886        format ("no dnstream intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below reference surface is icd+1
c        Ensure that downstream surface lies "below" upstream surface
         if (icd .lt. icu) then
            icd=icu
            xcd=xcu
            ycd=ycu
         endif

c     Intersection of (x,ycurveg) with new divertor plate
c     NOTE: start the search at the top-of-mesh index, ict, to avoid
c           a possible problem with multiple intersections (MER 95/09/05)
         call intersect2(rplate,zplate,1,nplate,
     &                    xcrv(ict:i2crv),ycrv(ict:i2crv),ict,i2crv,
     &                    xcp,ycp,kcp,icp,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,888) j
 888        format ("no new plate intersection for x,ycurveg on j=",i3)
            write (STDOUT,889)
 889        format ("*** You may have to extend the plate end points")
            call xerrab("")
         endif
c        index of first point below new divertor plate is icp+1

c     Intersection of (x,ycurveg) with old divertor plate
         call intersect2(rplate0,zplate0,1,nplate0,
     &                    xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                    xcp0,ycp0,kcp0,icp0,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,890) j
 890        format ("no old plate intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below old divertor plate is icp0+1

c     Distance from top-of-mesh to each data point along x,ycurveg
         xcrv(ict)=xct
         ycrv(ict)=yct
         dsc(ict)=0.
         do ii=ict+1,i2crv
            dsc(ii)=dsc(ii-1) +
     &                   sqrt( (xcrv(ii)-xcrv(ii-1))**2 +
     &                         (ycrv(ii)-ycrv(ii-1))**2 )
         enddo

c     Distance from top-of-mesh to upstream reference
         dscu=dsc(icu)+sqrt( (xcu-xcrv(icu))**2 +
     &                       (ycu-ycrv(icu))**2 )

c     Distance from top-of-mesh to downstream reference
         dscd=dsc(icd)+sqrt( (xcd-xcrv(icd))**2 +
     &                       (ycd-ycrv(icd))**2 )

c     Distance from top-of-mesh to new divertor plate
         dscp=dsc(icp)+sqrt( (xcp-xcrv(icp))**2 +
     &                       (ycp-ycrv(icp))**2 )

c     Distance from top-of-mesh to old divertor plate
         dscp0=dsc(icp0)+sqrt( (xcp0-xcrv(icp0))**2 +
     &                         (ycp0-ycrv(icp0))**2 )

c     Ensure that downstream reference is above new divertor plate
         if (dscp .lt. dscd) then
            icd=icp
            xcd=xcp
            ycd=ycp
            dscd=dscp
         endif

c     Use dmix0 to define length of mesh transition region
         len = dmix0*(dscd-dscu)
         dscd = dscu + len
         # find new icd,xcd,ycd
         do i=icu,icp
            if ( (dsc(i) .le. dscd) .and. (dscd .lt. dsc(i+1)) ) then
               icd=i
               frac=(dscd-dsc(i))/(dsc(i+1)-dsc(i))
               xcd=xcrv(i)+frac*(xcrv(i+1)-xcrv(i))
               ycd=ycrv(i)+frac*(ycrv(i+1)-ycrv(i))
               break
            endif
         enddo
         # find new imd (index of first cmeshx,y point upstream of x,ycd)
         # Strategy: compute the distance along x,ycurveg to ith meshpoint
         # and when this distance exceeds dscd we have identified imd+1
         do i=imu+1,i2msh
            dscmi=0.
            do ii=ict,i2crv
               if ( twixt(xcrv(ii),xmsh(i),xcrv(ii+1),delmax) .and.
     &              twixt(ycrv(ii),ymsh(i),ycrv(ii+1),delmax) ) then
                  dscmi=dsc(ii)+sqrt( (xmsh(i)-xcrv(ii))**2
     &                               +(ymsh(i)-ycrv(ii))**2 )
                  break
               endif
            enddo
            if (dscmi .gt. dscd) then
               imd=i-1
               break
            endif
         enddo

c----------------------------------------------------------------------c
c     Compute new meshpoints (along x,ycurveg) for each grid option
c----------------------------------------------------------------------c
c     First, the ismmon=0 option:
         dsmesh0(imt)=0.
         do i=imt+1,i2msh
c     Intersection of x,ycurveg with ith angle (orthogonal) surface
            na=3
            if (j .eq. jmin(region)) then
               ra(1)=2*cmeshx0(i,j)-cmeshx0(i,j+1)
               za(1)=2*cmeshy0(i,j)-cmeshy0(i,j+1)
            else
               ra(1)=cmeshx0(i,j-1)
               za(1)=cmeshy0(i,j-1)
            endif
            ra(2)=cmeshx0(i,j)
            za(2)=cmeshy0(i,j)
            if (j .eq. jmax(region)) then
               ra(3)=2*cmeshx0(i,j)-cmeshx0(i,j-1)
               za(3)=2*cmeshy0(i,j)-cmeshy0(i,j-1)
            else
               ra(3)=cmeshx0(i,j+1)
               za(3)=cmeshy0(i,j+1)
            endif
            call intersect2(ra,za,1,na,
     &                    xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                    xca,yca,kca,ica,fuzzm,ierr)
            if (ierr .ne. 0) then
               write (STDOUT,891) j,i
 891           format ("no angle intersection for x,ycurveg on j=",i3,
     &                 " at i=",i3)
               call xerrab("")
            endif
c     Index of first (x,ycrv) point below angle surface is ica+1
            dsmesh0(i)=dsc(ica)+sqrt( (xca-xcrv(ica))**2
     &                               +(yca-ycrv(ica))**2 )
         enddo
c----------------------------------------------------------------------c
c     Then, the ismmon=1 option:
         do i=imt,i2msh
            if (i .le. imu) then
               dsmesh1(i)=dsmesh0(i)
            else
               dsmesh1(i)=dscu+(dscp-dscu)*
     &                     (dsmesh0(i)-dscu)/(dscp0-dscu)
            endif
         enddo
c----------------------------------------------------------------------c
c     Then, the ismmon=2 option:
         fracu=(dsmu-dsm(imu))/(dsm(imu+1)-dsm(imu))
         dssu=dss(imu)+fracu*(dss(imu+1)-dss(imu))
         do i=imt,i2msh
            if (i .le. imu) then
               dsmesh2(i)=dsmesh0(i)
            else
               if ( (jsptrx(1) .le. j) .and. (j .le. jsptrx(2)) ) then
                  dsmesh2(i)=dscu+(dscp-dscu)*
     &              (dssleg(i)-dssleg(imu))/(dssleg(i2msh)-dssleg(imu))
               else
                  dsmesh2(i)=dscu+(dscp-dscu)*
     &              (dss(i)-dssu)/(dss(i2msh)-dssu)
               endif
            endif
         enddo
c----------------------------------------------------------------------c

c     Combine the 1 & 2 meshes according to weight factor wtmesh1:
	    do i=imt,i2msh
	       dsmesh(i)=wtmesh1*dsmesh1(i)+(1-wtmesh1)*dsmesh2(i)
	    enddo
c     Find the index of the mesh point just above the downstream surface
            do i=imt+1,i2msh
               im12d=i-1
               if ( (dsmesh(i)   .ge. dscd) .and. 
     &              (dsmesh(i-1) .le. dscd) ) break
            enddo

c     Mix in the 0 mesh, starting with weight factor 0 at the downstream
c     position and increasing linearly to weight factor 1 at the upstream
c     position:

	 do i=imt,i2msh
            if (i .le. imu) then
               dsmesh(i)=dsmesh0(i)
            elseif ( (imu .lt. i) .and. (i .le. im12d) ) then
               dsmesh(i)=dscu+(dscd-dscu)*
     &            (dsmesh0(i)-dscu)/(dsmesh0(i)-dscu+dscd-dsmesh(i))
            else
		# use dsmesh as is
            endif
	 enddo


c     Modified cmeshx,cmeshy values --
            do i=imt+1,i2msh
               do ii=ict+1,i2crv
                  if (dsmesh(i) .le. dsc(ii)) then
c     modified point i lies between (x,ycurveg) points ii-1 and ii
                     frac=(dsmesh(i)-dsc(ii-1))/
     &                    (dsc(ii)  -dsc(ii-1))
                     cmeshx(i,j)=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                     cmeshy(i,j)=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                     break
                  endif
               enddo
            enddo

      enddo	# end of loop over all j flux surfaces


      return
      end

c----------------------------------------------------------------------c
      function pointdist(x1, y1, x2, y2) result(dist)
        IMPLICIT NONE
        REAL, INTENT(IN):: x1, x2, y1, y2
        REAL dist
        dist = ((x1-x2)**2 + (y1-y2)**2)**0.5
      end function pointdist

      subroutine corexpointfix(region)
        IMPLICIT NONE
        INTEGER, INTENT(IN):: region
        Use(Dimflxgrd)	#jdim,npts,noregs
        Use(Comflxgrd)
        Use(Dimensions)
        Use(Curves)	# npointg,xcurveg,ycurveg
        Use(Linkco)	# cmeshx,cmeshy
        Use(Mmod)	# nupstream,rupstream,zupstream as temporary arrays
                        # ndnstream,rdnstream,zdnstream as temporary arrays
                # xcrv,ycrv,dsc,fuzzm,wtold
        Use(Transfm)	# ijump
        Use(Inmesh)     # ilmax
        INTEGER:: j, i, ii, jmx, jmn, imin, kcu, icu, ierr, dir
        REAL:: xxpt, yxpt, pointdist, xline(2), yline(2)
        REAL:: xintersect, yintersect, xline2(2), yline2(2), dist
        REAL:: deltadist, totdist

        ! Set up indices and bounds depending on the region being worked on    
        if (region .eq. 1) then
            jmn = jsptrx(1) + 1
            jmx = jmax(1)
            dir = 1
        else
            jmn = jsptrx(2) - 1
            jmx = jmax(1) + 2
            dir = -1
        endif
        ! Set up totdist to account for expansion towards sep
        deltadist = (MIN(maxdist(region), 0.4) - MIN(mindist(region),0.2)) / (ABS(jmn-jmx))
        totdist = MIN(maxdist(region),0.4) + deltadist
        ! Loop through all flux surfaces from sep to core
        do j = jmn, jmx, dir
            ! Work arrays
            xcrv = xcurveg(:,j)
            ycrv = ycurveg(:,j)
            ! Add curvature to core cells
            totdist = totdist - deltadist
            ! Get midplane index
            if (region .eq. 1) then
                imin = MINLOC(cmeshx(:ixpoint(1, region),j),1)
            else
                imin = MAXLOC(cmeshx(:ixpoint(1,region),j),1)
            endif
            ! Loop through poloidal cells between midplane and X-point 
            do i = ixpoint(1, region) - 1, imin, -1
                ! Get distance to the neighboring grid cell
                dist = pointdist(
     .                      cmeshx(i,j), cmeshy(i,j), 
     .                      cmeshx(i+1,j), cmeshy(i+1,j)
     .          )
                ! If the seed point is within the minimum distance
                ! or on the _wrong_ side of the neighboring seed point, 
                ! shift it to the minimum distance
                if ( 
     .              (dist < totdist) .or. 
     .              (cmeshx(i,j)*dir > dir*cmeshx(i+1,j) ) 
     .          ) then
                    ! Define vertical line that is mindist from the X-point
                    if (region .eq. 1) then
                        xline = cmeshx(i+1,j) - totdist
                    else
                        xline = cmeshx(i+1,j) + totdist
                    endif
                    yline(1) = cmeshy(i+1,j)-0.1
                    yline(2) = cmeshy(i+1,j)+0.1
                    ! Find intersect between vertical line and flux surface
                    call intersect2(    xline, yline, 1, 2,
     .                          xcrv, ycrv, 1 , ijump(j),
     .                          xintersect, yintersect, kcu, icu,
     .                          fuzzm,ierr)
                    ! Set new seed point to intersect with mindist line
                    cmeshx(i,j) = xintersect
                    cmeshy(i,j) = yintersect
            endif
        end do
        end do


      end subroutine corexpointfix


      subroutine smooth(i,j1,j2)
      implicit none
      integer i,j1,j2
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)	# npointg,xcurveg,ycurveg
Use(Linkco)	# cmeshx,cmeshy
Use(Mmod)	# nupstream,rupstream,zupstream as temporary arrays
                # ndnstream,rdnstream,zdnstream as temporary arrays
		# xcrv,ycrv,dsc,fuzzm,wtold
Use(Transfm)	# ijump
      external intersect2, remark, xerrab

c     Local variables --
      integer ii,j,i1crv,i2crv,icu,kcu,icp,kcp,ierr,region,istart,iend
      real xcu,ycu,xcp,ycp,xjm1,yjm1,distnew,frac

c     Smooth spatial irregularities in angle-like surfaces of the mesh

      xjm1=cmeshx(i,j1)
      yjm1=cmeshy(i,j1)

      do j=j1+1,j2-1

	 if (j .lt. jaxis) then
	    region=1
	 elseif (j .gt. jaxis) then
	    region=2
	 else
	    call xerrab("smooth: region not defined")
	 endif

c     Temporary variables for (xcurveg,ycurveg) surface data
         i1crv=1		# start index
         i2crv=npointg(j)	# end index
         do ii=i1crv,i2crv
            xcrv(ii)=xcurveg(ii,j)
            ycrv(ii)=ycurveg(ii,j)
         enddo
   
c     Straight line between j-1 and j+1 points on angle-like surface
	 nupstream=2
	 rupstream(1)=cmeshx(i,j-1)
	 zupstream(1)=cmeshy(i,j-1)
	 rupstream(2)=cmeshx(i,j+1)
	 zupstream(2)=cmeshy(i,j+1)

c     Intersection of straight line with x,ycurveg is (xcu,ycu):
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcu,ycu,kcu,icu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j,i
 886        format 
     & ("smooth: no intersection for x,ycurveg on j=",i3," with i=",i3)
            call xerrab("")
         endif
c     The point (x,ycu) lies between icu and icu+1 on x,ycrv.

c     Segmented line through j-1, j and j+1 points on angle-like surface
         ndnstream=3
         rdnstream(1)=cmeshx(i,j-1)
         zdnstream(1)=cmeshy(i,j-1)
         rdnstream(2)=cmeshx(i,j)
         zdnstream(2)=cmeshy(i,j)
         rdnstream(3)=cmeshx(i,j+1)
         zdnstream(3)=cmeshy(i,j+1)

c     Intersection of segmented line with x,ycurveg is (xcp,ycp):
         call intersect2(rdnstream,zdnstream,1,ndnstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcp,ycp,kcp,icp,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,887) j,i
 887        format 
     & ("smooth: no intersection for x,ycurveg on j=",i3," with i=",i3)
            call xerrab("")
         endif
c     The point (x,ycp) lies between icp and icp+1 on x,ycrv.

c     Update meshpoint at j-1, now that we are done using it:
         cmeshx(i,j-1)=xjm1
         cmeshy(i,j-1)=yjm1

         if (icp .lt. icu) then		# x,ycp is  upstream  from x,ycu
            istart=icp
            xcrv(istart)=xcp
            ycrv(istart)=ycp
            iend=icu+1
            xcrv(iend)=xcu
            ycrv(iend)=ycu
            dsc(istart)=0.
            do ii=istart+1,iend
               dsc(ii)=dsc(ii-1)
     &                     + sqrt( (xcrv(ii)-xcrv(ii-1))**2
     &                            +(ycrv(ii)-ycrv(ii-1))**2 )
            enddo
            distnew=(1-wtold)*dsc(iend)
            do ii=istart+1,iend
               if ( distnew .le. dsc(ii) ) then
c        new meshpoint lies between ii-1 and ii
                  frac = (distnew-dsc(ii-1))
     &                         /(dsc(ii)-dsc(ii-1))
                  xjm1=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                  yjm1=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                  break
               endif
            enddo
         elseif (icp .gt. icu) then	# x,ycp is downstream from x,ycu
            istart=icu
            xcrv(istart)=xcu
            ycrv(istart)=ycu
            iend=icp+1
            xcrv(iend)=xcp
            ycrv(iend)=ycp
            dsc(istart)=0.
            do ii=istart+1,iend
               dsc(ii)=dsc(ii-1)
     &                     + sqrt( (xcrv(ii)-xcrv(ii-1))**2
     &                            +(ycrv(ii)-ycrv(ii-1))**2 )
            enddo
            distnew=wtold*dsc(iend)
            do ii=istart+1,iend
               if (distnew .le. dsc(ii)) then
c        new meshpoint lies between ii-1 and ii
                  frac = (distnew-dsc(ii-1))
     &                         /(dsc(ii)-dsc(ii-1))
                  xjm1=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                  yjm1=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                  break
               endif
            enddo
         else	# both x,ycu and x,ycp lie on same segment of x,ycrv
            xjm1=(1-wtold)*xcu+wtold*xcp
            yjm1=(1-wtold)*ycu+wtold*ycp
         endif

      enddo

      cmeshx(i,j2-1)=xjm1
      cmeshy(i,j2-1)=yjm1

      return
      end

c----------------------------------------------------------------------c

      subroutine orthogrd(ixtyp,i,j0,j,xob,yob,alphab)
      implicit none
      integer ixtyp, i, j0, j, ii, ipoint, iseg, nsys, noiter, iter
      integer isegalt,numerr,numerrmx
      real xob, yob, alphab, distmin, distance, alphaf, alpha
      real cosalpf, sinalpf, xo, yo, sinalpb, cosalpb, xb, yb, xm
      real ym, xoldv, sold, xlimit, dx, dxsq, xpdx, xmdx, s, sp
      real sm, s1, s2, xnew, tanbet, beta
Use(Dimflxgrd)	#jdim,npts
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)
Use(Transfm)
Use(Linkco)
      data numerrmx /2/
c     NOTE:
c     On entry, (xob,yob,alphab) give the (R,Z) coordinates and local
c     rotation angle for a mesh point on the j0th flux surface.
c
c     On exit, (xob,yob,alphab) give the (R,Z) coordinates and local
c     rotation angle for a mesh point on the jth flux surface.
c
c     These points define a segment of the ith orthogonal mesh surface
c     connecting the j0 and j flux surfaces.
c
c     Procedure:
c
c     *** find the nearest point to (xob,yob), on the xycurveg( ,j) :
c     *** ixtyp is a flag that tells us whether to search the upper,
c     *** lower, or entire range of flux surface data.
c
      distmin=(xubnd-xlbnd)+(yubnd-ylbnd)
c
      if (ixtyp.eq.0) then
      do 10 ii=1,npointg(j)
         distance=sqrt((xcurveg(ii,j)-xob)**2+(ycurveg(ii,j)-yob)**2)
         if(distance.gt.distmin) goto 10
         distmin=distance
         ipoint=ii
   10 continue
c
      elseif (ixtyp.eq.1) then
      do 30 ii=1,ijump(j)
         distance=sqrt((xcurveg(ii,j)-xob)**2+(ycurveg(ii,j)-yob)**2)
         if(distance.gt.distmin) goto 30
         distmin=distance
         ipoint=ii
   30 continue
c
      elseif (ixtyp.eq.2) then
      do 40 ii=ijump(j)+1,npointg(j)
         distance=sqrt((xcurveg(ii,j)-xob)**2+(ycurveg(ii,j)-yob)**2)
         if(distance.gt.distmin) goto 40
         distmin=distance
         ipoint=ii
   40 continue
c
      endif
c
c     *** determine the segment on which this point lies.
         iseg=isegment(ipoint,j)
         isegalt=iseg
         numerr=0
   44 continue	# possible re-start from here if out-of-range spline
		# evaluation error occurs in loop 20 from s.r. intrcept
         if (isegalt .ne. iseg) then
            numerr=numerr+1
            if (numerr .gt. numerrmx) then
               call remark ("*** error in subroutine orthogrd")
               call remark ("Too many out-of-range spline errors for")
               write (STDOUT,901) i,j
  901          format ("orthogonal surface i = ",i3,
     &                 " at flux surface j = ",i3)
               call xerrab("")
            endif
            iseg=isegalt
         endif
         if (iseg .eq. 0) then
            call remark ("*** error in subroutine orthogrd")
            write (STDOUT,899) ipoint,j
  899       format ("data point i=",i3," on contour j=",i3," is not",
     &              " assigned to a spline segment")
            call xerrab("")
         endif
         nsys=isys(iseg,j)
         alphaf=alphasys(nsys)
c     *** transform the moving system origin into the fixed system.
         alpha=alphab-alphaf
         cosalpf=cos(alphaf)
         sinalpf=sin(alphaf)
         xo= xob*cosalpf+yob*sinalpf
         yo=-xob*sinalpf+yob*cosalpf
c     *** transform the coordinates of the nearest point into the
c     *** moving system (from the basic system).
         sinalpb=sin(alphab)
         cosalpb=cos(alphab)
c
         xb=xcurveg(ipoint,j)
         yb=ycurveg(ipoint,j)
         xm=+(xb-xob)*cosalpb+(yb-yob)*sinalpb
         ym=-(xb-xob)*sinalpb+(yb-yob)*cosalpb
c
c     *** obtain the solution of the orthogonality equation using xm and ym
c     *** as the starting point.
         xoldv=xm
         sold=ym
         xlimit=0.001*distmin
         dx=0.01*distmin
         dxsq=dx*dx
c
c     *** start the iteration process.
         noiter=20
c
      do 20 iter=1,noiter
c     *** we need the first and the second derivatives of the spline.
c
         xpdx=xoldv+dx
         xmdx=xoldv-dx
      call intrcept(iseg,j,xo,yo,alpha,xoldv,sold,s,isegalt)
      if (isegalt .ne. iseg) then
         goto 44	# re-start with alternative spline representation
      endif
      call intrcept(iseg,j,xo,yo,alpha,xpdx,sold,sp,isegalt)
      if (isegalt .ne. iseg) then
         goto 44	# re-start with alternative spline representation
      endif
      call intrcept(iseg,j,xo,yo,alpha,xmdx,sold,sm,isegalt)
      if (isegalt .ne. iseg) then
         goto 44	# re-start with alternative spline representation
      endif
c
c     *** form the derivatives of s.
         s1=0.5*(sp-sm)/dx
         s2=(sp-2.0*s+sm)/dxsq
c
         xnew=xoldv-(2.*s*xoldv+(s*s-xoldv*xoldv)*s1)/
     &        (2.*s*(1.+s1*s1)+(s*s-xoldv*xoldv)*s2)
c
      if(abs(xnew-xoldv).lt.xlimit) goto 21
         xoldv=xnew
         sold=s
c
   20 continue
c
      goto 22
   21 continue
c     *** xnew and s are the desired solutions in the moving system.
c     *** transform into the basic system, and make this point the new
c     *** origin of the moving system.
         xob=xob+xnew*cosalpb-s*sinalpb
         yob=yob+xnew*sinalpb+s*cosalpb
c     *** compute the angle of the new moving coordinate system, as viewed
c     *** from the old moving coordinate system.
         tanbet=s1
         beta=atan(tanbet)
         alpha=alpha+beta
c     *** transform into the basic system.
         alphab=alpha+alphaf
c
      return
   22 continue
c     *** iterations have failed to converge
         write(STDOUT,900) i,j0,j
  900    format("orthogrd not converged for i=",i2," surface",
     &          " between j0=",i2," and j=",i2," flux contours")
c
      return
      end

c     ------------------------------------------------------------------

      subroutine orthogx(ixtyp,i,j0,j,xob,yob,alphab)
      implicit none
      integer ixtyp, i, j0, j, ii, ipoint, iseg, nsys, noiter, iter
      integer isegalt, numerr, numerrmx
      real xob, yob, alphab, distmin, distance, alphaf, alpha
      real cosalpf, sinalpf, xo, yo, sinalpb, cosalpb, xb, yb, xm
      real ym, xoldv, sold, xlimit, dx, dxsq, xpdx, xmdx, s, sp
      real sm, s1, s2, xnew, tanbet, beta
Use(Dimflxgrd)	#jdim,npts
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)
Use(Transfm)
Use(Linkco)
      data numerrmx /2/
c     NOTE:
c     On entry, (xob,yob) give the (R,Z) coordinates for a point
c     (the x-point) on the j0th flux surface (the separatrix).
c
c     On exit,  (xob,yob,alphab) give the (R,Z) coordinates and local
c     rotation angle for a mesh point on the jth flux surface.
c     
c     *** finds the orthogonal set of curves.
c
c     *** find the nearest point to (xob,yob), on the xycurveg( ,j).
c
         distmin=(xubnd-xlbnd)+(yubnd-ylbnd)
c
c
      if(ixtyp.eq.0) then
      do 10 ii=1,npointg(j)
         distance=sqrt((xcurveg(ii,j)-xob)**2+(ycurveg(ii,j)-yob)**2)
      if(distance.gt.distmin) goto 10
         distmin=distance
         ipoint=ii
   10 continue
c
      else if(ixtyp.eq.1) then
      do 30 ii=1,ijump(j)
         distance=sqrt((xcurveg(ii,j)-xob)**2+(ycurveg(ii,j)-yob)**2)
      if(distance.gt.distmin) goto 30
         distmin=distance
         ipoint=ii
   30 continue
c
      else if(ixtyp.eq.2) then
      do 40 ii=ijump(j)+1,npointg(j)
         distance=sqrt((xcurveg(ii,j)-xob)**2+(ycurveg(ii,j)-yob)**2)
      if(distance.gt.distmin) goto 40
         distmin=distance
         ipoint=ii
   40 continue
c
      endif
c
         distmin=sqrt(distmin)
c
c     *** determine the segment on which this point lies.
         iseg=isegment(ipoint,j)
         isegalt=iseg
         numerr=0
   44 continue	# possible re-start from here if out-of-range spline
		# evaluation error occurs in loop 20 from s.r. intrcept
         if (isegalt .ne. iseg) then
            numerr=numerr+1
            if (numerr .gt. numerrmx) then
               call remark ("*** error in subroutine orthogx")
               call remark ("Too many out-of-range spline errors for")
               write (STDOUT,901) i,j
  901          format ("orthogonal surface i = ",i3,
     &                 " at flux surface j = ",i3)
               call xerrab("")
            endif
            iseg=isegalt
         endif
         if (iseg .eq. 0) then
            call remark ("*** error in subroutine orthogx")
            write (STDOUT,899) ipoint,j
  899       format ("data point i=",i3," on contour j=",i3," is not",
     &              " assigned to a spline segment")
            call xerrab("")
         endif
         nsys=isys(iseg,j)
         alphaf=alphasys(nsys)
c
c     *** orient the moving system the same way as the fixed system.
         alphab=alphaf
c
c     *** transform the moving system origin into the fixed system.
         alpha=alphab-alphaf
         cosalpf=cos(alphaf)
         sinalpf=sin(alphaf)
         xo= xob*cosalpf+yob*sinalpf
         yo=-xob*sinalpf+yob*cosalpf
c     *** transform the coordinates of the nearest point into the
c     *** moving system (from the basic system).
         sinalpb=sin(alphab)
         cosalpb=cos(alphab)
c
         xb=xcurveg(ipoint,j)
         yb=ycurveg(ipoint,j)
         xm=+(xb-xob)*cosalpb+(yb-yob)*sinalpb
         ym=-(xb-xob)*sinalpb+(yb-yob)*cosalpb
c
c     *** obtain the solution of the orthogonality equation using xm and ym
c     *** as the starting point.
         xoldv=xm
         sold=ym
         xlimit=0.001*distmin
         dx=0.01*distmin
         dxsq=dx*dx
c
c     *** start the iteration process.
         noiter=20
c
      do 20 iter=1,noiter
c     *** we need the first and the second derivatives of the spline.
c
         xpdx=xoldv+dx
         xmdx=xoldv-dx
      call intrcept(iseg,j,xo,yo,alpha,xoldv,sold,s,isegalt)
      if (isegalt .ne. iseg) then
         goto 44	# re-start with alternative spline representation
      endif
      call intrcept(iseg,j,xo,yo,alpha,xpdx,sold,sp,isegalt)
      if (isegalt .ne. iseg) then
         goto 44	# re-start with alternative spline representation
      endif
      call intrcept(iseg,j,xo,yo,alpha,xmdx,sold,sm,isegalt)
      if (isegalt .ne. iseg) then
         goto 44	# re-start with alternative spline representation
      endif
c
c     *** form the derivatives of s.
         s1=0.5*(sp-sm)/dx
         s2=(sp-2.0*s+sm)/dxsq
c
         xnew=xoldv-(xoldv+s*s1)/(1.+s1*s1+s*s2)
c
      if(abs(xnew-xoldv).lt.xlimit) goto 21
         xoldv=xnew
         sold=s
c
   20 continue
c
      goto 22
   21 continue
c     *** xnew and s are the desired solutions in the moving system.
c     *** transform into the basic system, and make this point the new
c     *** origin of the moving system.
         xob=xob+xnew*cosalpb-s*sinalpb
         yob=yob+xnew*sinalpb+s*cosalpb
c     *** compute the angle of the new moving coordinate system, as viewed
c     *** from the old moving coordinate system.
         tanbet=s1
         beta=atan(tanbet)
         alpha=alpha+beta
c     *** transform into the basic system.
         alphab=alpha+alphaf
c
      return
   22 continue
c     *** iterations have failed to converge
         write(STDOUT,900) i,j0,j
  900    format("orthogx not converged for i=",i2," surface",
     &          " between j0=",i2," and j=",i2," flux contours")

      return
      end

c     ------------------------------------------------------------------

      subroutine prune
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)
Use(Dimensions)
Use(Curves)
Use(Inmesh)
Use(Linkco)
Use(System)

c     local variables --
      real ds,dstest
      integer i,j,k,ii,k1,ixpointb,numpoint,region

c     loop over all regions --
      do 800 region=1,noregs

c     Find the minimum distance between x-point and separatrix contour
      j=jsptrx(region)
      dstest=(xubnd-xlbnd)+(yubnd-ylbnd)
      do 10 k=1,npointg(j)
         ds=sqrt((xcurveg(k,j)-rseps)**2+(ycurveg(k,j)-zseps)**2)
         dstest=min(dstest,ds)
   10 continue

      dsminx = 2*max(dstest,dsmin)

c     First, prune the points which are within dsminx of the x-point
c     on the separatrix --

      j=jsptrx(region)
         i=0
         k=0
         ixpointb=0
   11 continue
      if ((k .eq. i+1) .and. (ixpointb .eq. 0)) then
c     *** insert the first x-point into the sequence.
         i=i+1
         ixpointb=i
         xwork(i)=rseps
         ywork(i)=zseps
      endif
c
         k=k+1
      if (k .gt. npointg(j)) goto 12
         ds=sqrt((xcurveg(k,j)-rseps)**2+(ycurveg(k,j)-zseps)**2)
      if (ds .lt. dsminx) goto 11
         i=i+1
         xwork(i)=xcurveg(k,j)
         ywork(i)=ycurveg(k,j)
      goto 11
   12 continue
         numpoint=i

c     Re-define the data for the separatrix contour --
      do 20 ii=1,npointg(j)
         xcurveg(ii,j)=0.0
         ycurveg(ii,j)=0.0
   20 continue
      do 30 ii=1,numpoint
         xcurveg(ii,j)=xwork(ii)
         ycurveg(ii,j)=ywork(ii)
   30 continue
      npointg(j)=numpoint

c     Now, for each contour, prune all points which are separated from
c     their neighbors by less than dsmin --

      do 100 j=jmin(region),jmax(region)
         i=1
         k=1
         xwork(i)=xcurveg(k,j)
         ywork(i)=ycurveg(k,j)
   40 continue
         k1=k
   50 continue
         k=k+1
      if (k .gt. npointg(j)) goto 60
         ds=sqrt((xcurveg(k,j)-xcurveg(k1,j))**2 +
     &           (ycurveg(k,j)-ycurveg(k1,j))**2)
      if (ds .lt. dsmin) goto 50
      if ((j .eq. jsptrx(region)) .and. (k1 .eq. ixpointb)) then
         ixpointc(1,region)=i
         i=i+1
         ixpointc(2,region)=i
         xwork(i)=xcurveg(k1,j)
         ywork(i)=ycurveg(k1,j)
         i=i+1
         ixpointc(3,region)=i
         xwork(i)=xcurveg(k1,j)
         ywork(i)=ycurveg(k1,j)
      endif
         i=i+1
         xwork(i)=xcurveg(k,j)
         ywork(i)=ycurveg(k,j)
      goto 40
   60 continue
         numpoint=i

c     Re-define the data for the jth flux contour --
         do 70 ii=1,npointg(j)
            xcurveg(ii,j)=0.0
            ycurveg(ii,j)=0.0
   70    continue
         do 80 ii=1,numpoint
            xcurveg(ii,j)=xwork(ii)
            ycurveg(ii,j)=ywork(ii)
   80    continue
         npointg(j)=numpoint

  100 continue

  800 continue          # end of loop over regions

      return
      end

c     ------------------------------------------------------------------

      subroutine seed0(i,j,x0,y0,xcurveg,ycurveg,cmeshx,cmeshy,
     &                       nptsg,idimg,jdimg,lim1,lim2,km,kp,testflag)
      implicit none
      integer i, j, nptsg, idimg, jdimg, lim1, lim2, k, km, kp, testflag
      real x0, y0, dx1, dx2, dy1, dy2, delta
      real xcurveg(nptsg,jdimg),ycurveg(nptsg,jdimg)
      real cmeshx(idimg,jdimg),cmeshy(idimg,jdimg)

      external twixt,remark,xerrab
      logical twixt,xtest,ytest

c     This subroutine finds the indicies (km,kp) of the (x,ycurveg)
c     data points that straddle the top-of-mesh position (x0,y0)
c     on the jth flux contour.  For testflag=0, only the horizontal
c     coordinate is tested; for testflag=2, only the vertical
c     coordinate is tested; for testflag=1, a more stringent test
c     on both horizontal and vertical coordinates is used to
c     locate the end-of-mesh.  This is useful when there are multiple
c     solutions for the test on horizontal position only.  Finally, 
c     the coordinates of the top-of-mesh data point (x0,y0)
c     are put into the mesh arrays (cmeshx,y) for angle index i and
c     flux index j.

      data delta/0/	# treat data as exact in twixt tests

c     In the following,
c        xtest is TRUE if the k and k+1 data points straddle x0
c        ytest is TRUE if the k and k+1 data points straddle y0
c     Otherwise, xtest and ytest are FALSE.

      km=0
      kp=0
      if (testflag==0) then		# test on x only
         do k=lim1,lim2-1
            xtest=twixt(xcurveg(k,j), x0, xcurveg(k+1,j), delta)
            if (xtest) then
               km=k
               kp=k+1
               break
            endif
         enddo
         if (km==0) then			# test not satisfied
            call remark("***     failure in subroutine seed0      ***")
            call remark("*** first point not found on separatrix ***")
            call xerrab("")
         else				# interpolate to find y
            dx1=abs(xcurveg(km,j)-x0)
            dx2=abs(xcurveg(kp,j)-x0)
            cmeshx(i,j)=x0
            cmeshy(i,j)=(dx1*ycurveg(kp,j)+dx2*ycurveg(km,j))/(dx1+dx2)
         endif

      elseif (testflag==1) then		# test on x and y
         do k=lim1,lim2-1
            xtest=twixt(xcurveg(k,j), x0, xcurveg(k+1,j), delta)
            ytest=twixt(ycurveg(k,j), y0, ycurveg(k+1,j), delta)
            if (xtest .and. ytest) then
               km=k
               kp=k+1
               break
            endif
         enddo
         if (km==0) then			# test not satisfied
            call remark("***     failure in subroutine seed0      ***")
            call remark("*** first point not found on separatrix ***")
            call xerrab("")
         else				# put data in cmeshx,y
            cmeshx(i,j)=x0
            cmeshy(i,j)=y0
         endif

      elseif (testflag==2) then		# test on y only
         do k=lim1,lim2-1
            ytest=twixt(ycurveg(k,j), y0, ycurveg(k+1,j), delta)
            if (ytest) then
               km=k
               kp=k+1
               break
            endif
         enddo
         if (km==0) then			# test not satisfied
            call remark("***     failure in subroutine seed0      ***")
            call remark("*** first point not found on separatrix ***")
            call xerrab("")
         else				# interpolate to find x
            dy1=abs(ycurveg(km,j)-y0)
            dy2=abs(ycurveg(kp,j)-y0)
            cmeshx(i,j)=(dy1*xcurveg(kp,j)+dy2*xcurveg(km,j))/(dy1+dy2)
            cmeshy(i,j)=y0
         endif

      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine seedl(i,j,xl,yl,xcurveg,ycurveg,cmeshx,cmeshy,
     &                       nptsg,idimg,jdimg,lim1,lim2,km,kp,testflag)
      implicit none
      integer i, j, nptsg, idimg, jdimg, lim1, lim2, k, km, kp, testflag
      real xl, yl, dx1, dx2, dy1, dy2, delta
      real xcurveg(nptsg,jdimg),ycurveg(nptsg,jdimg)
      real cmeshx(idimg,jdimg),cmeshy(idimg,jdimg)

      external twixt,remark,xerrab
      logical twixt,xtest,ytest

c     This subroutine finds the indicies (km,kp) of the (x,ycurveg)
c     data points that straddle the end-of-mesh position (xl,yl)
c     on the jth flux contour.  For testflag=0, only the horizontal
c     coordinate is tested; for testflag=2, only the vertical
c     coordinate is tested; for testflag=1, a more stringent test
c     on both horizontal and vertical coordinates is used to
c     locate the end-of-mesh.  This is useful when there are multiple
c     solutions for the test on horizontal position only.  Finally, 
c     the coordinates of the end-of-mesh data point (xl,yl)
c     are put into the mesh arrays (cmeshx,y) for angle index i and
c     flux index j.

      data delta/0/	# treat data as exact in twixt tests

c     In the following,
c        xtest is TRUE if the k and k+1 data points straddle xl
c        ytest is TRUE if the k and k+1 data points straddle yl
c     Otherwise, xtest and ytest are FALSE.

      km=0
      kp=0
      if (testflag==0) then		# test on x only
         do k=lim1,lim2-1
            xtest=twixt(xcurveg(k,j), xl, xcurveg(k+1,j), delta)
            if (xtest) then
               km=k
               kp=k+1
               break
            endif
         enddo
         if (km==0) then			# test not satisfied
            call remark("***     failure in subroutine seedl      ***")
            call remark("*** strike point not found on separatrix ***")
            call xerrab("")
         else				# interpolate to find y
            dx1=abs(xcurveg(km,j)-xl)
            dx2=abs(xcurveg(kp,j)-xl)
            cmeshx(i,j)=xl
            cmeshy(i,j)=(dx1*ycurveg(kp,j)+dx2*ycurveg(km,j))/(dx1+dx2)
         endif

      elseif (testflag==1) then		# test on x and y
         do k=lim1,lim2-1
            xtest=twixt(xcurveg(k,j), xl, xcurveg(k+1,j), delta)
            ytest=twixt(ycurveg(k,j), yl, ycurveg(k+1,j), delta)
            if (xtest .and. ytest) then
               km=k
               kp=k+1
               break
            endif
         enddo
         if (km==0) then			# test not satisfied
            call remark("***     failure in subroutine seedl      ***")
            call remark("*** strike point not found on separatrix ***")
            call xerrab("")
         else				# put data in cmeshx,y
            cmeshx(i,j)=xl
            cmeshy(i,j)=yl
         endif

      elseif (testflag==2) then		# test on y only
         do k=lim1,lim2-1
            ytest=twixt(ycurveg(k,j), yl, ycurveg(k+1,j), delta)
            if (ytest) then
               km=k
               kp=k+1
               break
            endif
         enddo
         if (km==0) then			# test not satisfied
            call remark("***     failure in subroutine seedl      ***")
            call remark("*** strike point not found on separatrix ***")
            call xerrab("")
         else				# interpolate to find x
            dy1=abs(ycurveg(km,j)-yl)
            dy2=abs(ycurveg(kp,j)-yl)
            cmeshx(i,j)=(dy1*xcurveg(kp,j)+dy2*xcurveg(km,j))/(dy1+dy2)
            cmeshy(i,j)=yl
         endif

      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine sow
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)	# rseps,zseps
Use(Dimensions)
Use(Curves)
Use(Inmesh)	# x0g,xlast,ylast,isztest
Use(Linkco)
Use(System)
Use(Xmesh)
Use(Share)      # islimon,theta_split,ishalfm
Use(Mmod)

c     local variables --
      integer region,rbegin,ict,kct,ierr
      integer i,i0,i1,ic,ii,im,ip,ix,ixtop,j,lim1,lim2,nxtotal
      real x1,x2,xr0,xrl,xtop,y1,y2,yr0,yrl,ds,ds1,ds2,t,xct,yct

c     procedures --
      real xfcn,xfcn2,xfcn3,xfcn4
      external seedl,xfcn,xfcn2,xfcn3,xfcn4,intersect2

c     This subroutine generates the seed points along the separatrix.

c     Loop over regions --
      rbegin = 1
      if (ishalfm==1) rbegin = 2
      do 10 region=rbegin,noregs
         j=jsptrx(region)

c----------------------------------------------------------------------c
c     Insert the first seed point --

         i0=1
         xr0=x0g(region)
         yr0=y0g(region)
         lim1=1
         lim2=ixpointc(1,region)

         if (istpnew==1) then	# user-defined first-seed-point
            call seed0(i0,j,xr0,yr0,xcurveg,ycurveg,cmeshx,cmeshy,
     &              npts,idim,jdim,lim1,lim2,im,ip,istptest(region))
         else
c           Define top-of-mesh surface
	    ntop=2
	    if ((geometry .eq. 'snull').or.(geometry .eq. "uppersn")) then
                                                # straight line outward
                                                # from magnetic axis
                                                # along theta_split
               rtop(1)=rmagx
               ztop(1)=zmagx
               rtop(2)=rtop(1)+max(xdim,zdim)*cos(theta_split)
               ztop(2)=ztop(1)+max(xdim,zdim)*sin(theta_split)
            else
						# horizontal line
						# through magnetic axis
               rtop(1)=0*rmagx
               ztop(1)=zmagx
               rtop(2)=2*rmagx
               ztop(2)=zmagx
	    endif
c           Intersection of (x,ycurveg) with top-of-mesh surface
            call intersect2(rtop,ztop,1,ntop,
     &                   xcurveg(lim1:lim2,j),ycurveg(lim1:lim2,j),lim1,lim2,
     &                   xct,yct,kct,ict,fuzzm,ierr)
            if (ierr .ne. 0) then
               write (STDOUT,884) j
 884           format ("sow: no top-of-mesh for x,ycurveg on j=",i3)
               call xerrab("")
            endif
c           index of first point below top-of-mesh is ict+1
            cmeshx(i0,j)=xct
            cmeshy(i0,j)=yct
            im=ict
            ip=ict+1
         endif

c        Compute the distance along the separatrix as measured from
c        the first seed point
         dissep(1,region)=0.
         x1=cmeshx(i0,j)
         y1=cmeshy(i0,j)
         ii=1
         istartc(region)=ip	# from subroutine seed0 above
      do 30 i=istartc(region),npointg(j)
         ii=ii+1
         x2=xcurveg(i,j)
         y2=ycurveg(i,j)
         dissep(ii,region)=dissep(ii-1,region)+
     &                     sqrt((x2-x1)**2+(y2-y1)**2)
         x1=x2
         y1=y2
         if(i.eq.ixpointc(1,region)) ix=ii
   30 continue

c     Distance from first seed point to x-point is
         distxp(region)=dissep(ix,region)

c----------------------------------------------------------------------c
c     Insert the last seed point --

      if (nxleg(igrid,region)==0) then
c        Modification for FRC mesh generation
         cmeshx(ilmax(region),j) = rseps
         cmeshy(ilmax(region),j) = zseps
         distxpxl(region)=0.
      else
         i1=ilmax(region)
         xrl=xlast(region)
         yrl=ylast(region)
         lim1=ixpointc(3,region)
         lim2=npointg(j)
         call seedl(i1,j,xrl,yrl,xcurveg,ycurveg,cmeshx,cmeshy,
     &              npts,idim,jdim,lim1,lim2,im,ip,isztest(region))

c     Compute distance from x-point to last seed point --
         ds=sqrt((cmeshx(ilmax(region),j)-xcurveg(im,j))**2
     &          +(cmeshy(ilmax(region),j)-ycurveg(im,j))**2)
         iendc(region)=im+1
         dissep(iendc(region)-istartc(region)+2,region)=
     &               dissep(im-istartc(region)+2,region)+ds
         distxpxl(region)=
     &               dissep(iendc(region)-istartc(region)+2,region)
     &                 - distxp(region)
      endif

c     End of loop over regions
   10 continue



c     Compute data for xfcn(t) --
      ndat=5
      xdat(1)=0.
      xdat(2)=xdat(1)+distxpxl(1)
      xdat(3)=xdat(2)+distxp(1)
      xdat(4)=xdat(3)+distxp(2)
      xdat(5)=xdat(4)+distxpxl(2)
      nxtotal=nxleg(igrid,1)+nxuse(1)
     &                      +nxuse(2)+nxleg(igrid,2)
      tdat(1)=0.
      tdat(2)=tdat(1)+float(nxleg(igrid,1))/float(nxtotal)
      tdat(3)=tdat(2)+float(nxuse(1))/float(nxtotal)
      tdat(4)=tdat(3)+float(nxuse(2))/float(nxtotal)
      tdat(5)=tdat(4)+float(nxleg(igrid,2))/float(nxtotal)

      if (kxmesh .eq. 0) then           # use manual x-mesh as for PLANET code

      rbegin = 1
      if (ishalfm==1) rbegin = 2
      do 53 region=rbegin,noregs
      do 51 i=1,ixpoint(1,region)
         seed(i,region)=0.01*seedxp(i,region)*distxp(region)
   51 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 52 i=ixpoint(3,region),ilmax(region)
         ii=i-ixpoint(2,region)
         seed(i,region)=seed(ixpoint(2,region),region)+
     &                  0.01*seedxpxl(ii,region)*distxpxl(region)
   52 continue
   53 continue

      elseif (kxmesh .eq. 1) then       # use analytic x-mesh from xfcn(t)

      ixtop = nxuse(1) + nxleg(igrid,1)      # B2 code poloidal index
                                        # ix at PLANET code index i=1
      xtop = distxp(1) + distxpxl(1)    # B2 code poloidal distance x at
                                        # PLANET code index i=1
      if (ishalfm==0) then
         region=1
         do 61 i=1,ixpoint(1,region)
            ix = ixtop - (i-1)
            t = float(ix)/float(nxtotal)
            seed(i,region) = xtop - xfcn(t)
   61    continue
         seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
         do 62 i=ixpoint(3,region),ilmax(region)
            ix = ixtop - (i-3)
            t = float(ix)/float(nxtotal)
            seed(i,region) = xtop - xfcn(t)
   62    continue
      endif      # end if-test on ishalfm

      region=2
      do 63 i=1,ixpoint(1,region)
         ix = ixtop + (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn(t) - xtop
   63 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 64 i=ixpoint(3,region),ilmax(region)
         ix = ixtop + (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn(t) - xtop
   64 continue

      elseif (kxmesh .eq. 2) then       # use analytic x-mesh from xfcn2(t)

      ixtop = nxuse(1) + nxleg(igrid,1)      # B2 code poloidal index
                                        # ix at PLANET code index i=1
      xtop = distxp(1) + distxpxl(1)    # B2 code poloidal distance x at
                                        # PLANET code index i=1
      region=1
      do 71 i=1,ixpoint(1,region)
         ix = ixtop - (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xtop - xfcn2(t)
   71 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 72 i=ixpoint(3,region),ilmax(region)
         ix = ixtop - (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xtop - xfcn2(t)
   72 continue

      region=2
      do 73 i=1,ixpoint(1,region)
         ix = ixtop + (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn2(t) - xtop
   73 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 74 i=ixpoint(3,region),ilmax(region)
         ix = ixtop + (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn2(t) - xtop
   74 continue

      elseif (kxmesh .eq. 3) then       # use analytic x-mesh from xfcn3(t)
c     insert extra data points --
         ndat=7
         ndatp2=ndat+2
         tdat(7)=tdat(5)
         xdat(7)=xdat(5)
         tdat(6)=tdat(7)-dt2
         xdat(6)=xdat(7)-dx2
         do i=5,3,-1
           tdat(i)=tdat(i-1)
           xdat(i)=xdat(i-1)
         enddo
         tdat(2)=tdat(1)+dt1
         xdat(2)=xdat(1)+dx1
c     set derivative values at endpoints --
         dleft  = dxgas(1) * nxtotal
         dright = dxgas(2) * nxtotal
c     compute spline coefficients for xfcn3(t) --
         call xcscoef

      ixtop = nxuse(1) + nxleg(igrid,1)       # B2 code poloidal index ix at
                                        # PLANET code index i=1
      xtop = distxp(1) + distxpxl(1)    # B2 code poloidal distance x at
                                        # PLANET code index i=1
      region=1
      do 81 i=1,ixpoint(1,region)
         ix = ixtop - (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xtop - xfcn3(t)
   81 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 82 i=ixpoint(3,region),ilmax(region)
         ix = ixtop - (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xtop - xfcn3(t)
   82 continue

      region=2
      do 83 i=1,ixpoint(1,region)
         ix = ixtop + (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn3(t) - xtop
   83 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 84 i=ixpoint(3,region),ilmax(region)
         ix = ixtop + (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn3(t) - xtop
   84 continue

      elseif (kxmesh .eq. 4) then       # use analytic x-mesh from xfcn4(t)
         ndat=5
         ndatp2=ndat+2
c     re-define endpoints of data for spline
         alfx(1)=max(.0001,alfx(1))
         alfx(2)=max(.0001,alfx(2))
         tdat(1)=float(nxgas(1))/float(nxtotal)
         xdat(1)=dxgas(1)*(exp(alfx(1)*nxgas(1))-1.)
     &                   /(exp(alfx(1))-1.)
         tdat(5)=tdat(5)-float(nxgas(2))/float(nxtotal)
         xdat(5)=xdat(5)-dxgas(2)*(exp(alfx(2)*nxgas(2))-1.)
     &                           /(exp(alfx(2))-1.)
c     set derivative values at endpoints --
         ileft = 1
         dleft = dxgas(1)*alfx(1)*nxtotal*exp(alfx(1)*nxgas(1))
     &                                  /(exp(alfx(1))-1.)
         iright = 1
         dright = dxgas(2)*alfx(2)*nxtotal*exp(alfx(2)*nxgas(2))
     &                                   /(exp(alfx(2))-1.)
c     compute spline coefficients for xfcn4(t) --
         call xcscoef

      ixtop = nxuse(1) + nxleg(igrid,1)       # B2 code poloidal index ix at
                                        # PLANET code index i=1
      xtop = distxp(1) + distxpxl(1)    # B2 code poloidal distance x at
                                        # PLANET code index i=1
      region=1
      do 85 i=1,ixpoint(1,region)
         ix = ixtop - (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xtop - xfcn4(t,nxtotal)
   85 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 86 i=ixpoint(3,region),ilmax(region)
         ix = ixtop - (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xtop - xfcn4(t,nxtotal)
   86 continue

      region=2
      do 87 i=1,ixpoint(1,region)
         ix = ixtop + (i-1)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn4(t,nxtotal) - xtop
   87 continue
      seed(ixpoint(2,region),region) = seed(ixpoint(1,region),region)
      do 88 i=ixpoint(3,region),ilmax(region)
         ix = ixtop + (i-3)
         t = float(ix)/float(nxtotal)
         seed(i,region) = xfcn4(t,nxtotal) - xtop
   88 continue

      endif
c
c     *** finally obtain the coordinates of the seed points by interpolation.
c
c     Loop over regions --
      rbegin = 1
      if (ishalfm==1) rbegin = 2
      do 102 region=rbegin,noregs
      j=jsptrx(region)

      ii=1
      do 92 i=2,ilmax(region)-1
      ii=ii-1
   91 continue
      ii=ii+1
      if(dissep(ii,region).gt.seed(i,region)) then
         ds1=seed(i,region)-dissep(ii-1,region)
         ds2=dissep(ii,region)-seed(i,region)
         ds=ds1+ds2
         ic=ii+istartc(region)-2
         if(ic.eq.istartc(region)) then
            cmeshx(i,j)=(ds1*xcurveg(ic,j)+ds2*cmeshx(1,j))/ds
            cmeshy(i,j)=(ds1*ycurveg(ic,j)+ds2*cmeshy(1,j))/ds
         else if(ic.eq.iendc(region)) then
            cmeshx(i,j)=(ds1*cmeshx(ilmax(region),j)
     .                    +ds2*xcurveg(ic-1,j))/ds
            cmeshy(i,j)=(ds1*cmeshy(ilmax(region),j)
     .                    +ds2*ycurveg(ic-1,j))/ds
         else
            cmeshx(i,j)=(ds1*xcurveg(ic,j)+ds2*xcurveg(ic-1,j))/ds
            cmeshy(i,j)=(ds1*ycurveg(ic,j)+ds2*ycurveg(ic-1,j))/ds
         endif
      else
         goto 91
      endif
   92 continue

c     End loop over regions --
  102 continue

      return
      end

c     ------------------------------------------------------------------

      subroutine exponseed

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculate seed points using exponentials from plates and from X-pts  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      Use(Inmesh)          #seedxp,seedxpxl
      Use(Dimflxgrd)       #noregs,npts for use in Inmesh arrays
      Use(Dimensions)      #idim for use in Inmesh arrays
      Use(Share)           #geometry,nxcore
      Use(Comgeo)          #xfs
      Use(Xpoint_indices)  #ixpt1,ixpt2
      Use(UEgrid)          #ixtop
      Use(Dim)             #nxm
      Use(Expseed)         #fraclplt,alfxdiv,alfxcore,rleg2core_len,
                           #shift_seed_leg,shift_seed_core,nxlplt,nxlixpt,
                           #fcorenunif

      integer im,ir,immax,isunifm,nxcore_use(2)
      real seedunifm,rleg2core_len(2)

      if (geometry=="snull") then
        nxcore_use(1) = nxcore(1,1)
        nxcore_use(2) = nxcore(1,2)
      elseif (geometry == "dnbot") then
        nxcore_use(1) = nxcore(1,1) - 1
        nxcore_use(2) = nxcore(1,2) - 1
      endif

c.. Compute seedxpxl in divertor legs
      rleg2core_len(1) = xfs(ixpt1(1))/(xfs(ixtop) - xfs(ixpt1(1)))
      rleg2core_len(2) = (xfs(nxm)-xfs(ixpt2(1)))/(xfs(ixpt2(1)) - 
     .                                                   xfs(ixtop))

c..  First do mesh near the plate (seeds large, start at 100)
c  ########################################################### 
      immax = nxlplt(1) + nxlxpt(1)
      seedxpxl(immax+1,1) = 100.
      do im = immax, nxlxpt(1), -1
        ir = immax - im
        seedxpxl(im,1) = 100. - 100.*fraclplt(1)*
     .                             (exp(alfxdiv(1)*(ir+1)) - 1) / 
     .                             (exp(alfxdiv(1)*nxlplt(1)) - 1)
      enddo
      immax = nxlplt(2) + nxlxpt(2)
      seedxpxl(immax+1,2) = 100.
      do im = immax, nxlxpt(2), -1
        ir = immax - im
        seedxpxl(im,2) = 100. - 100.*fraclplt(2)*
     .                             (exp(alfxdiv(2)*(ir+1)) - 1) / 
     .                             (exp(alfxdiv(2)*nxlplt(2)) - 1)
      enddo

c..  Do divertor legs near X-point (seeds small, start at 0)
c  #############################################################
      seedxpxl(1,1) = 0.
      do im = 2, nxlxpt(1)
        seedxpxl(im,1) = shift_seed_leg(1) + \
                     (100.*(1.-fraclplt(1)) - shift_seed_leg(1))* 
     .                              (exp(alfxdiv(1)*(im-1)) - 1) / 
     .                              (exp(alfxdiv(1)*nxlxpt(1)) - 1)
      enddo

      seedxpxl(1,2) =0.
      do im = 2, nxlxpt(2)
        seedxpxl(im,2) = shift_seed_leg(2) + \
                     (100.*(1.-fraclplt(2)) - shift_seed_leg(2))* 
     .                              (exp(alfxdiv(2)*(im-1)) - 1) / 
     .                              (exp(alfxdiv(2)*nxlxpt(2)) - 1)
      enddo

c  Last, do in/out core regions; large seeds near x-point
c  ###########################################################

      seedxp(nxcore_use(1)+1,1) = 100.
      seedxp(nxcore_use(1),1) = 100.-rleg2core_len(1)*seedxpxl(2,1)-
     .                                      shift_seed_core(1)
      isunifm = 0
      do im = nxcore_use(1)-1, 2, -1
        if (isunifm == 0) then
          seedxp(im,1) = 100. - (100. - seedxp(im+1,1)) - 
     .                          rleg2core_len(1)*seedxpxl(2,1)* 
     .                         exp(alfxcore(1)*float(nxcore_use(1)-im))
          seedunifm = seedxp(im,1)/float(im-1)
          if (seedxp(im+1,1)-seedxp(im,1)>fcorenunif*seedunifm) isunifm=1
        elseif (isunifm == 1) then
          seedxp(im,1) = seedxp(im+1,1) - seedunifm
        endif
        seedxp(1,1) = 0.
      enddo

      seedxp(nxcore_use(2)+1,2) = 100.
      seedxp(nxcore_use(2),2) = 100. - rleg2core_len(2)*seedxpxl(2,2) - 
     .                                      shift_seed_core(2)
      isunifm = 0
      do im = nxcore_use(2)-1, 2, -1
        if (isunifm == 0) then
          seedxp(im,2) = 100. - (100. - seedxp(im+1,2)) - 
     .                          rleg2core_len(2)*seedxpxl(2,2)* 
     .                         exp(alfxcore(2)*float(nxcore_use(2)-im))
          seedunifm = seedxp(im,2)/float(im-1)
          if (seedxp(im+1,1)-seedxp(im,1)>fcorenunif*seedunifm) isunifm=1
        elseif (isunifm == 1) then
          seedxp(im,2) = seedxp(im+1,2) - seedunifm
        endif
        seedxp(1,2) = 0.
      enddo

      return
      end
c  end of subroutine exponseed
c     ------------------------------------------------------------------

      subroutine splfit
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Curves)
Use(Inmesh)
Use(Linkco)
Use(System)
Use(Transfm)
Use(Share)    # cutlo,ishalfm

c     local variables --
      integer region,rbegin,j,iseg,is,icood,ii
      real z1,z2,s,t,alpha,rlimit

c     procedures --
      external clear,codsys

c     a small real number whose reciprocal does not overflow --
      rlimit = cutlo

c     initialization -- define rotation angles of fixed coordinate systems
      alphasys(1) = 0
      alphasys(2) = Pi/2
      alphasys(3) = Pi
      alphasys(4) = 3*Pi/2

c     clear the workspace and output arrays --
      call clear

c     loop over all regions --
      rbegin = 1
      if (ishalfm==1) rbegin = 2
      do 9696 region=rbegin,noregs

c     loop over all flux contours --
      do 1 j=jmin(region),jmax(region)

      iseg=0
      is=1

    2 continue          # top of loop over spline-fit segments

      if(is.gt.npointg(j)) goto 1000    # jump out of loop 2

c     icood=1, flux curve goes right, transform to system beta=0.
c     icood=2, flux curve goes up, transform to system beta=Pi/2.
c     icood=3, flux curve goes left, transform to system beta=Pi.
c     icood=4, flux curve goes down, transform to system, beta=3*Pi/2.

c     check the initial slope --
         z1=ycurveg(is,j)-ycurveg(is+1,j)
         z2=xcurveg(is,j)-xcurveg(is+1,j)
         if (abs(z2) .le. abs(z1)*rlimit) then
            alpha=90.
         else
            s=z1/z2
            t=atan(s)
            alpha=(180./Pi)*t
         endif
         if(abs(alpha) .le. alpha1) then
            if(xcurveg(is+1,j) .lt. xcurveg(is,j)) then
               icood=3
            else
               icood=1
            endif
         else
            if (ycurveg(is+1,j) .lt. ycurveg(is,j)) then
               icood=4
            else
               icood=2
            endif
         endif

      call codsys(j,icood,iseg,is,dyjump,region,alpha1)

      goto 2            # end of loop over spline-fit segments


 1000 continue
      nseg(j)=iseg
      if(j.eq.jsptrx(region)) then
         ii=ixpointc(2,region)
         isegment(ii-2,j)=isegment(ii-3,j)
         isegment(ii-1,j)=isegment(ii-3,j)
         isegment(ii+1,j)=isegment(ii+3,j)
         isegment(ii+2,j)=isegment(ii+3,j)
      endif
      ii=npointg(j)
      isegment(ii,j)=isegment(ii-2,j)
      isegment(ii-1,j)=isegment(ii-2,j)

    1 continue          # end of loop over flux contours

 9696 continue          # end of loop over regions

      return
      end

c     ------------------------------------------------------------------

      subroutine transys(gamma,k1,k2,j)
      implicit none
      integer k1, k2, j, k, l
      real gamma, rad, a, b
Use(Dimflxgrd)	#jdim,npts
Use(Dimensions)
Use(Curves)
Use(Transit)
c
      do 233 k=1,npts
         xtrans(k)=0.0
         ytrans(k)=0.0
 233  continue
c
      rad=(Pi/180.)*gamma
c
      a=cos(rad)
      b=sin(rad)
c
      do 9299 l=k1,k2
         xtrans(l-k1+1)=xcurveg(l,j)*a+ycurveg(l,j)*b
         ytrans(l-k1+1)=-xcurveg(l,j)*b+ycurveg(l,j)*a
 9299 continue
c
      return
      end

c     ------------------------------------------------------------------

      real function xfcn(t)
      implicit none
      real t
Use(Xmesh)
      real t1,t2,t3,t4,t5,x1,x2,x3,x4,x5,
     &     aa,beta2,beta4,x2p,x3p,x4p

      t1=tdat(1)
      t2=tdat(2)
      t3=tdat(3)
      t4=tdat(4)
      t5=tdat(5)
      x1=xdat(1)
      x2=xdat(2)
      x3=xdat(3)
      x4=xdat(4)
      x5=xdat(5)

c     Compute parameters for analytic representation of the x-mesh;
c     e.g., on t1 < t < t2 use the rational function
c
c           x(t) = x1 + (x2-x1)*(t-t1)*(beta2-t2)/((t2-t1)*(beta2-t))
c
c     and on t2 < t < t4  use a cubic that passes through points
c     (t2,x2), (t3,x3), and (t4,x4) with adjustable slope at (t3,x3)
c
c           x(t) =   x2*(t-t3)*(t-t4)/((t2-t3)*(t2-t4))
c                  + x3*(t-t2)*(t-t4)/((t3-t2)*(t3-t4))
c                  + x4*(t-t2)*(t-t3)/((t4-t2)*(t4-t3))
c                  + aa*(t-t2)*(t-t3)*(t-t4)
c
c     Choose beta2 so that x'(t) is continuous at t2.
c     And similarly for the outboard leg of the divertor,
c     on the interval t4 < t < t5 use the representation
c
c           x(t) = x5 + (x4-x5)*(t-t5)*(beta4-t4)/((t4-t5)*(beta4-t))
c
c     Choose beta4 so that x'(t) is continuous at t4.


c     The slope at t3 with only the quadratic terms included is

      x3p =      x2*(t3-t4)/((t2-t3)*(t2-t4))
     &    + x3*(2*t3-t2-t4)/((t3-t2)*(t3-t4))
     &         + x4*(t3-t2)/((t4-t2)*(t4-t3))

c     The cubic term enhances this slope by a factor slpxt if we choose

      aa = (slpxt-1)*x3p/((t3-t2)*(t3-t4))

c     The slopes at t2 and t4 are then given by

      x2p = x2*(2*t2-t3-t4)/((t2-t3)*(t2-t4))
     &         + x3*(t2-t4)/((t3-t2)*(t3-t4))
     &         + x4*(t2-t3)/((t4-t2)*(t4-t3)) + aa*(t2-t3)*(t2-t4)

      x4p =      x2*(t4-t3)/((t2-t4)*(t2-t3))
     &         + x3*(t4-t2)/((t3-t4)*(t3-t2))
     &    + x4*(2*t4-t3-t2)/((t4-t3)*(t4-t2)) + aa*(t4-t2)*(t4-t3)

c     For continuity at t2 and t4 we choose

      if (t2==t1) then
         beta2 = 999999999.   # dummy value -- beta2 will not be used
      else
         beta2 = (x2p*t2*(t2-t1)-(x2-x1)*t1)/(x2p*(t2-t1)-(x2-x1))
      endif
      if (t5==t4) then
         beta4 = 999999999.   # dummy value -- beta4 will not be used
      else
         beta4 = (x4p*t4*(t4-t5)-(x4-x5)*t5)/(x4p*(t4-t5)-(x4-x5))
      endif

c     Finally, x(t) is given by

      if (t .lt. t1) then
         xfcn=x1
      elseif ((t1 .le. t) .and. (t .lt. t2)) then
         xfcn=x1+(x2-x1)*(t-t1)*(beta2-t2)/((t2-t1)*(beta2-t))
      elseif ((t2 .le. t) .and. (t .le. t4)) then
         xfcn=x2*(t-t3)*(t-t4)/((t2-t3)*(t2-t4))
     &          + x3*(t-t2)*(t-t4)/((t3-t2)*(t3-t4))
     &             + x4*(t-t2)*(t-t3)/((t4-t2)*(t4-t3))
     &                + aa*(t-t2)*(t-t3)*(t-t4)
      elseif ((t4 .lt. t) .and. (t .le. t5)) then
         xfcn=x5+(x4-x5)*(t-t5)*(beta4-t4)/((t4-t5)*(beta4-t))
      elseif (t5 .lt. t) then
         xfcn=x5
      endif

      return
      end

c     ------------------------------------------------------------------

      real function xfcn2(t)
      implicit none
      real t
Use(Xmesh)
      real t1,t2,t3,t4,t5,x1,x2,x3,x4,x5,
     &     aa,beta2,beta4,x2p,x3p,x4p

      t1=tdat(1)
      t2=tdat(2)
      t3=tdat(3)
      t4=tdat(4)
      t5=tdat(5)
      x1=xdat(1)
      x2=xdat(2)
      x3=xdat(3)
      x4=xdat(4)
      x5=xdat(5)

c     Compute parameters for analytic representation of the x-mesh;
c     e.g., on t1 < t < t2 use the exponential function
c
c           x(t) = x1 + (x2-x1)*(t-t1)*exp(beta2*(t-t2))/(t2-t1)
c
c     and on t2 < t < t4  use a cubic that passes through points
c     (t2,x2), (t3,x3), and (t4,x4) with adjustable slope at (t3,x3)
c
c           x(t) =   x2*(t-t3)*(t-t4)/((t2-t3)*(t2-t4))
c                  + x3*(t-t2)*(t-t4)/((t3-t2)*(t3-t4))
c                  + x4*(t-t2)*(t-t3)/((t4-t2)*(t4-t3))
c                  + aa*(t-t2)*(t-t3)*(t-t4)
c
c     Choose beta2 so that x'(t) is continuous at t2.
c     And similarly for the outboard leg of the divertor,
c     on the interval t4 < t < t5 use the representation
c
c           x(t) = x5 + (x4-x5)*(t-t5)*exp(beta4*(t-t4))/(t4-t5)
c
c     Choose beta4 so that x'(t) is continuous at t4.

c     The slope at t3 with only the quadratic terms included is

      x3p =      x2*(t3-t4)/((t2-t3)*(t2-t4))
     &    + x3*(2*t3-t2-t4)/((t3-t2)*(t3-t4))
     &         + x4*(t3-t2)/((t4-t2)*(t4-t3))

c     The cubic term enhances this slope by a factor slpxt if we choose

      aa = (slpxt-1)*x3p/((t3-t2)*(t3-t4))

c     The slopes at t2 and t4 are then given by

      x2p = x2*(2*t2-t3-t4)/((t2-t3)*(t2-t4))
     &         + x3*(t2-t4)/((t3-t2)*(t3-t4))
     &         + x4*(t2-t3)/((t4-t2)*(t4-t3)) + aa*(t2-t3)*(t2-t4)

      x4p =      x2*(t4-t3)/((t2-t4)*(t2-t3))
     &         + x3*(t4-t2)/((t3-t4)*(t3-t2))
     &    + x4*(2*t4-t3-t2)/((t4-t3)*(t4-t2)) + aa*(t4-t2)*(t4-t3)

c     For continuity at t2 and t4 we choose

      beta2 = (x2p*(t2-t1)-(x2-x1))/((t2-t1)*(x2-x1))
      beta4 = (x4p*(t4-t5)-(x4-x5))/((t4-t5)*(x4-x5))

c     Finally, x(t) is given by

      if (t .lt. t1) then
         xfcn2=x1
      elseif ((t1 .le. t) .and. (t .lt. t2)) then
         xfcn2=x1+(x2-x1)*(t-t1)*exp(beta2*(t-t2))/(t2-t1)
      elseif ((t2 .le. t) .and. (t .le. t4)) then
         xfcn2=x2*(t-t3)*(t-t4)/((t2-t3)*(t2-t4))
     &          + x3*(t-t2)*(t-t4)/((t3-t2)*(t3-t4))
     &             + x4*(t-t2)*(t-t3)/((t4-t2)*(t4-t3))
     &                + aa*(t-t2)*(t-t3)*(t-t4)
      elseif ((t4 .lt. t) .and. (t .le. t5)) then
         xfcn2=x5+(x4-x5)*(t-t5)*exp(beta4*(t-t4))/(t4-t5)
      elseif (t5 .lt. t) then
         xfcn2=x5
      endif

      return
      end

#----------------------------------------------------------------------#
      subroutine xcscoef
      implicit none
Use(Xmesh)
c     compute spline interpolant for x(t)

      call BINTh4 (tdat, xdat, ndat, ileft, iright, dleft, dright,
     &            kntopt, tknt, z1cscoef, ndatp2, kord, z1work)

      return
      end

#----------------------------------------------------------------------#
      real function xfcn3(t)
      implicit none
      real t
Use(Xmesh)
      real B1VAhL
      external B1VAhL
      integer inbv

c     evaluate the spline for x(t)

*
      inbv = 1		# initialization parameter for B1VAhL
         xfcn3 = B1VAhL (t, 0, tknt, ndatp2, kord, z1cscoef, inbv,
     &                  wrk1, iflag1)

      return
      end


#----------------------------------------------------------------------#
      real function xfcn4(t,nxtotal)
      implicit none
      integer nxtotal
      real t
Use(Xmesh)
      real B1VAhL
      external B1VAhL
      integer inbv

      if ((0. .le. t) .and. (t .lt. tdat(1))) then
c     use exponential form --
         xfcn4 = dxgas(1)*(exp(alfx(1)*t*nxtotal)-1.)
     &                       /(exp(alfx(1))-1.)

      elseif ((tdat(1) .le. t) .and. (t .lt. tdat(5))) then
c     evaluate the spline for x(t) --
         inbv = 1
         xfcn4 = B1VAhL (t, 0, tknt, ndatp2, kord, z1cscoef, inbv,
     &                  wrk1, iflag1)

      elseif ((tdat(5) .le. t) .and. (t .le. 1.)) then
c     use exponential form --
         xfcn4 = xdat(5) + dxgas(2)*(exp(alfx(2)*nxgas(2))-
     &            exp(alfx(2)*(1.-t)*nxtotal))   /(exp(alfx(2))-1.)

      endif

      return
      end

#----------------------------------------------------------------------#
      subroutine intersect2 (x1,y1,i1min,i1max,x2,y2,i2min,i2max,
     .                     xc,yc,i1c,i2c,fuzz,ierr)
      implicit none
      integer i1min,i1max,i2min,i2max,i1c,i2c,ierr
      real x1(i1min:i1max),y1(i1min:i1max)
      real x2(i2min:i2max),y2(i2min:i2max),xc,yc,fuzz

c     Find the intersection of the two segmented curves :
c        (x1(i),y1(i)) i=i1min,i1max and (x2(i),y2(i)) i=i2min,i2max
c     Return the intersection point (xc,yc) and the node indices i1c
c     and i2c such that the intersection point lies between nodes 
c     i1c and i1c+1 of curve 1 and nodes ic2 and ic2+1 of curve 2.
c     Uncertainty in data points is specified by input parameter fuzz.
c     Return error flag ierr=1 if no intersection is found

      external remark,twixt
      logical twixt

c     local variables --
      integer i1,i2
      real tan1,tan2

c     initialize error flag --
      ierr=1

      do i1=i1min,i1max-1
         i1c=i1
         if (x1(i1+1) .ne. x1(i1)) then
c           # non-vertical segment of curve 1
	    tan1=(y1(i1+1)-y1(i1))/(x1(i1+1)-x1(i1))
	    do i2=i2min,i2max-1
	       i2c=i2
               if (x2(i2+1) .ne. x2(i2)) then
c                 # non-vertical segment of curve 2
		  tan2=(y2(i2+1)-y2(i2))/(x2(i2+1)-x2(i2))
                  if (tan2 .ne. tan1) then
c                    # compute intersection of extended line segments
		     xc=((y1(i1)-x1(i1)*tan1)-(y2(i2)-x2(i2)*tan2))
     .                  /(tan2-tan1)
		     yc=y2(i2)+(xc-x2(i2))*tan2
c                    # test for actual intersection
                     if ( twixt(x2(i2),xc,x2(i2+1),fuzz) .and.
     .                    twixt(x1(i1),xc,x1(i1+1),fuzz) .and.
     .                    twixt(y2(i2),yc,y2(i2+1),fuzz) .and.
     .                    twixt(y1(i1),yc,y1(i1+1),fuzz) ) then
                        ierr=0
                        break (2) # exit do_i1
                     endif
                  endif
               else
c                 # vertical segment of curve 2
c                 # compute intersection of extended line segments
                  xc=x2(i2)
                  yc=y1(i1)+(xc-x1(i1))*tan1
c                 # test for actual intersection
                  if ( twixt(x2(i2),xc,x2(i2+1),fuzz) .and.
     .                 twixt(x1(i1),xc,x1(i1+1),fuzz) .and.
     .                 twixt(y2(i2),yc,y2(i2+1),fuzz) .and.
     .                 twixt(y1(i1),yc,y1(i1+1),fuzz) ) then
                     ierr=0
                     break (2) # exit do_i1
                  endif
               endif
	    enddo
         else
c           # vertical segment of curve 1
	    do i2=i2min,i2max-1
	       i2c=i2
               if (x2(i2+1) .ne. x2(i2)) then
c                 # non-vertical segment of curve 2
		  tan2=(y2(i2+1)-y2(i2))/(x2(i2+1)-x2(i2))
c                 # compute intersection of extended line segments
		  xc=x1(i1)
		  yc=y2(i2)+(xc-x2(i2))*tan2
c                 # test for actual intersection
                  if ( twixt(x2(i2),xc,x2(i2+1),fuzz) .and.
     .                 twixt(x1(i1),xc,x1(i1+1),fuzz) .and.
     .                 twixt(y2(i2),yc,y2(i2+1),fuzz) .and.
     .                 twixt(y1(i1),yc,y1(i1+1),fuzz) ) then
                     ierr=0
                     break (2) # exit do_i1
                  endif
               endif
	    enddo
         endif
      enddo # do_i1

      return
      end
 
#----------------------------------------------------------------------#

      logical function twixt(x1,xc,x2,dx)
      implicit none
      real x1,xc,x2,dx

c     test for whether xc lies between x1 and x2

c     dx = magnitude of possible uncertainty in values of x1,xc,x2

      twixt=.false.

      if (x1 .le. x2) then
         if ( (x1-dx .le. xc) .and. (xc .le. x2+dx) ) twixt=.true.
      elseif (x1 .gt. x2) then
         if ( (x2-dx .le. xc) .and. (xc .le. x1+dx) ) twixt=.true.
      endif

      return
      end

c----------------------------------------------------------------------c

      subroutine getu
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Mmod)
Use(Share)      # ismmon

      external gchange

# Define the default upstream reference surface arrays

      if (ismmon==2) then			# top+cut
c for the inboard half of the mesh:
	 nupstream1=jmax(1)-jmin(1)+3
	 call gchange("Mmod",0)
         call getu21
c for the outboard half of the mesh:
	 nupstream2=jmax(2)-jmin(2)+3
	 call gchange("Mmod",0)
         call getu22
      else
         if (isupstreamx==0) then               # midplane+cut
c for the inboard half of the mesh:
            nupstream1=jmax(1)-jsptrx(1)+3
            call gchange("Mmod",0)
            call getu11
c for the outboard half of the mesh:
	    nupstream2=jsptrx(2)-jmin(2)+3
	    call gchange("Mmod",0)
            call getu12
         elseif (isupstreamx==1) then		# surface thru x-point
c for the inboard half of the mesh:
	    nupstream1=jmax(1)-jmin(1)+1
	    call gchange("Mmod",0)
            call getu41
c for the outboard half of the mesh:
	    nupstream2=jmax(2)-jmin(2)+1
	    call gchange("Mmod",0)
            call getu42
	 endif
      endif
      return
      end

c----------------------------------------------------------------------c

      subroutine getu11
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)upstream1

      integer iref1,j,k
      real eps
      data eps /1./

      iref1=ixpoint(3,1)			# the cut in the p.f. region
      k=0
      do j=jmax(1),jsptrx(1),-1
         k=k+1
         rupstream1(k)=cmeshx(iref1,j)
         zupstream1(k)=cmeshy(iref1,j)
      enddo
      k=k+1
      zupstream1(k)=zmagx			# the magnetic axis
      rupstream1(k)=rmagx
      k=k+1
      zupstream1(k)=zmagx			# the inboard midplane
      rupstream1(k)=0.
c extrapolate end point slightly:
      rupstream1(1)=rupstream1(1)+eps*(rupstream1(1)-rupstream1(2))
      zupstream1(1)=zupstream1(1)+eps*(zupstream1(1)-zupstream1(2))

      return
      end

c----------------------------------------------------------------------c

      subroutine getu12
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)upstream2

      integer iref2,j,k
      real eps
      data eps /1./

      iref2=ixpoint(3,2)		# the cut in the p.f. region
      k=0
      do j=jmin(2),jsptrx(2)
         k=k+1
         rupstream2(k)=cmeshx(iref2,j)
         zupstream2(k)=cmeshy(iref2,j)
      enddo
      k=k+1
      zupstream2(k)=zmagx		# the magnetic axis
      rupstream2(k)=rmagx
      k=k+1
      zupstream2(k)=zmagx		# the outboard midplane
      rupstream2(k)=2*rmagx
c extrapolate end point slightly:
      rupstream2(1)=rupstream2(1)+eps*(rupstream2(1)-rupstream2(2))
      zupstream2(1)=zupstream2(1)+eps*(zupstream2(1)-zupstream2(2))

      return
      end

c----------------------------------------------------------------------c

      subroutine getu41
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# nupstream1,(r,z)upstream1

      integer iref1,j,k,m
      real eps
      data eps /1./

      iref1=ixpoint(3,1)		# the inboard PF cut surface
      k=0
      do j=jmax(1),jmin(1),-1
         k=k+1
         rupstream1(k)=cmeshx(iref1,j)
         zupstream1(k)=cmeshy(iref1,j)
      enddo
c extrapolate end points slightly:
      rupstream1(1)=rupstream1(1)+eps*(rupstream1(1)-rupstream1(2))
      zupstream1(1)=zupstream1(1)+eps*(zupstream1(1)-zupstream1(2))
      m=nupstream1
      rupstream1(m)=rupstream1(m)+eps*(rupstream1(m)-rupstream1(m-1))
      zupstream1(m)=zupstream1(m)+eps*(zupstream1(m)-zupstream1(m-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getu42
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# nupstream2,(r,z)upstream2

      integer iref2,j,k,m
      real eps
      data eps /1./

      iref2=ixpoint(3,2)		# the outboard PF cut surface
      k=0
      do j=jmin(2),jmax(2)
         k=k+1
         rupstream2(k)=cmeshx(iref2,j)
         zupstream2(k)=cmeshy(iref2,j)
      enddo
c extrapolate end points slightly:
      rupstream2(1)=rupstream2(1)+eps*(rupstream2(1)-rupstream2(2))
      zupstream2(1)=zupstream2(1)+eps*(zupstream2(1)-zupstream2(2))
      m=nupstream2
      rupstream2(m)=rupstream2(m)+eps*(rupstream2(m)-rupstream2(m-1))
      zupstream2(m)=zupstream2(m)+eps*(zupstream2(m)-zupstream2(m-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getu21
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)upstream1

      integer iref1,j,k
      real eps
      data eps /1./

      iref1=ixpoint(3,1)		# the cut in the p.f. region
      k=0
      do j=jmax(1),jsptrx(1),-1
         k=k+1
         rupstream1(k)=cmeshx(iref1,j)
         zupstream1(k)=cmeshy(iref1,j)
      enddo
      k=k+1
      zupstream1(k)=zmagx		# the magnetic axis
      rupstream1(k)=rmagx
      iref1=1				# the top of the mesh
      do j=jsptrx(1),jmin(1),-1
         k=k+1
         rupstream1(k)=cmeshx(iref1,j)
         zupstream1(k)=cmeshy(iref1,j)
      enddo
c extrapolate end points slightly:
      rupstream1(1)=rupstream1(1)+eps*(rupstream1(1)-rupstream1(2))
      zupstream1(1)=zupstream1(1)+eps*(zupstream1(1)-zupstream1(2))
      rupstream1(nupstream1)=rupstream1(nupstream1)
     &          +eps*(rupstream1(nupstream1)-rupstream1(nupstream1-1))
      zupstream1(nupstream1)=zupstream1(nupstream1)
     &          +eps*(zupstream1(nupstream1)-zupstream1(nupstream1-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getu22
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)upstream2

      integer iref2,j,k
      real eps
      data eps /1./

      iref2=ixpoint(3,2)		# the cut in the p.f. region
      k=0
      do j=jmin(2),jsptrx(2)
         k=k+1
         rupstream2(k)=cmeshx(iref2,j)
         zupstream2(k)=cmeshy(iref2,j)
      enddo
      k=k+1
      zupstream2(k)=zmagx		# the magnetic axis
      rupstream2(k)=rmagx
      iref2=1				# the top of the mesh
      do j=jsptrx(2),jmax(2)
         k=k+1
         rupstream2(k)=cmeshx(iref2,j)
         zupstream2(k)=cmeshy(iref2,j)
      enddo
c extrapolate end points slightly:
      rupstream2(1)=rupstream2(1)+eps*(rupstream2(1)-rupstream2(2))
      zupstream2(1)=zupstream2(1)+eps*(zupstream2(1)-zupstream2(2))
      rupstream2(nupstream2)=rupstream2(nupstream2)\
                +eps*(rupstream2(nupstream2)-rupstream2(nupstream2-1))
      zupstream2(nupstream2)=zupstream2(nupstream2)\
                +eps*(zupstream2(nupstream2)-zupstream2(nupstream2-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getp
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Mmod)

      external gchange

# Define the default (orthogonal) divertor plate arrays


c for the inboard half of the mesh:
      nplate1=jmax(1)-jmin(1)+1
      call gchange("Mmod",0)
      call getp1

c for the outboard half of the mesh:
      nplate2=jmax(2)-jmin(2)+1
      call gchange("Mmod",0)
      call getp2

      return
      end

c----------------------------------------------------------------------c

      subroutine getp1
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)	# ilmax
Use(Linkco)	# cmeshx,y
Use(Mmod)	# (r,z)plate1

      integer j,k,iref1
      real eps
      data eps /1./

      iref1=ilmax(1)		# the inboard-orthogonal-mesh plate

      k=0
      do j=jmax(1),jmin(1),-1
	 k=k+1
	 rplate1(k)=cmeshx(iref1,j)
	 zplate1(k)=cmeshy(iref1,j)
      enddo
c extrapolate end points slightly:
      rplate1(1)=rplate1(1)+eps*(rplate1(1)-rplate1(2))
      zplate1(1)=zplate1(1)+eps*(zplate1(1)-zplate1(2))
      rplate1(nplate1)=rplate1(nplate1)+eps*
     &                      (rplate1(nplate1)-rplate1(nplate1-1))
      zplate1(nplate1)=zplate1(nplate1)+eps*
     &                      (zplate1(nplate1)-zplate1(nplate1-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getp2
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)	# ilmax
Use(Linkco)	# cmeshx,y
Use(Mmod)	# (r,z)plate2

      integer j,k,iref2
      real eps
      data eps /1./

      iref2=ilmax(2)		# the outboard-orthogonal-mesh plate

      k=0
      do j=jmin(2),jmax(2)
	 k=k+1
	 rplate2(k)=cmeshx(iref2,j)
	 zplate2(k)=cmeshy(iref2,j)
      enddo
c extrapolate end points slightly:
      rplate2(1)=rplate2(1)+eps*(rplate2(1)-rplate2(2))
      zplate2(1)=zplate2(1)+eps*(zplate2(1)-zplate2(2))
      rplate2(nplate2)=rplate2(nplate2)+eps*
     &                      (rplate2(nplate2)-rplate2(nplate2-1))
      zplate2(nplate2)=zplate2(nplate2)+eps*
     &                      (zplate2(nplate2)-zplate2(nplate2-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine smoother
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)	# ilmax
Use(Linkco)	# ixpoint
Use(Mmod)
      external smooth

      integer n,i

# Smooth all angle-like surfaces in the mesh
#	wtold   = weight parameter in smoothing algorithm
#	nsmooth = number of passes

      do n=1,nsmooth

c Above the x-point (SOL and core, inboard and outboard):
	 do i=2,ixpoint(1,1)-1
	    call smooth(i,jmin(1),jmax(1))
	 enddo
	 do i=2,ixpoint(1,2)-1
	    call smooth(i,jmin(2),jmax(2))
	 enddo

c Through the x-point (SOL only, inboard and outboard):
         call smooth(ixpoint(1,1),jmin(1),jsptrx(1))
         call smooth(ixpoint(2,1),jmin(1),jsptrx(1))
         call smooth(ixpoint(3,1),jmin(1),jsptrx(1))
         call smooth(ixpoint(1,2),jsptrx(2),jmax(2))
         call smooth(ixpoint(2,2),jsptrx(2),jmax(2))
         call smooth(ixpoint(3,2),jsptrx(2),jmax(2))

c Below the x-point (SOL and p.f., inboard and outboard):
	 do i=ixpoint(3,1)+1,ilmax(1)
	    call smooth(i,jmin(1),jmax(1))
	 enddo
	 do i=ixpoint(3,2)+1,ilmax(2)
	    call smooth(i,jmin(2),jmax(2))
	 enddo

      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine smoother2
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)	# ilmax
Use(Linkco)	# ixpoint
Use(Mmod)
      external smooth

      integer n,i

# Smooth all angle-like surfaces in the outer half of the mesh only
#	wtold   = weight parameter in smoothing algorithm
#	nsmooth = number of passes

      do n=1,nsmooth

c Above the x-point (SOL and core, inboard and outboard):
	 do i=2,ixpoint(1,2)-1
	    call smooth(i,jmin(2),jmax(2))
	 enddo

c Through the x-point (SOL only, inboard and outboard):
         call smooth(ixpoint(1,2),jsptrx(2),jmax(2))
         call smooth(ixpoint(2,2),jsptrx(2),jmax(2))
         call smooth(ixpoint(3,2),jsptrx(2),jmax(2))

c Below the x-point (SOL and p.f., inboard and outboard):
	 do i=ixpoint(3,2)+1,ilmax(2)
	    call smooth(i,jmin(2),jmax(2))
	 enddo

      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine gett
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Mmod)

      external gchange

# Define the top-of-mesh reference surface arrays

c for the inboard half of the mesh:
	 ntop1=jmax(1)-jmin(1)+3
	 call gchange("Mmod",0)
         call gett1

c for the outboard half of the mesh:
	 ntop2=jmax(2)-jmin(2)+3
	 call gchange("Mmod",0)
         call gett2

      return
      end

c----------------------------------------------------------------------c

      subroutine gett1
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)top1

      integer iref1,j,k
      real eps
      data eps /1./

      iref1=ixpoint(3,1)		# the cut in the p.f. region
					# for psi < 1
      k=0
      do j=jmax(1),jsptrx(1),-1
         k=k+1
         rtop1(k)=cmeshx(iref1,j)
         ztop1(k)=cmeshy(iref1,j)
      enddo
      k=k+1
      ztop1(k)=zmagx			# the magnetic axis
      rtop1(k)=rmagx
      iref1=1				# the top of the mesh
					# for psi > 1
      do j=jsptrx(1),jmin(1),-1
         k=k+1
         rtop1(k)=cmeshx(iref1,j)
         ztop1(k)=cmeshy(iref1,j)
      enddo
c extrapolate end points slightly:
      rtop1(1)=rtop1(1)+eps*(rtop1(1)-rtop1(2))
      ztop1(1)=ztop1(1)+eps*(ztop1(1)-ztop1(2))
      rtop1(ntop1)=rtop1(ntop1)
     &          +eps*(rtop1(ntop1)-rtop1(ntop1-1))
      ztop1(ntop1)=ztop1(ntop1)
     &          +eps*(ztop1(ntop1)-ztop1(ntop1-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine gett2
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)top2

      integer iref2,j,k
      real eps
      data eps /1./

      iref2=ixpoint(3,2)		# the cut in the p.f. region
      k=0
      do j=jmin(2),jsptrx(2)
         k=k+1
         rtop2(k)=cmeshx(iref2,j)
         ztop2(k)=cmeshy(iref2,j)
      enddo
      k=k+1
      ztop2(k)=zmagx			# the magnetic axis
      rtop2(k)=rmagx
      iref2=1				# the top of the mesh
      do j=jsptrx(2),jmax(2)
         k=k+1
         rtop2(k)=cmeshx(iref2,j)
         ztop2(k)=cmeshy(iref2,j)
      enddo
c extrapolate end points slightly:
      rtop2(1)=rtop2(1)+eps*(rtop2(1)-rtop2(2))
      ztop2(1)=ztop2(1)+eps*(ztop2(1)-ztop2(2))
      rtop2(ntop2)=rtop2(ntop2)\
                +eps*(rtop2(ntop2)-rtop2(ntop2-1))
      ztop2(ntop2)=ztop2(ntop2)\
                +eps*(ztop2(ntop2)-ztop2(ntop2-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getd
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Mmod)

      external gchange

# Define the downstream reference surface arrays

c for the inboard half of the mesh:
	 ndnstream1=jmax(1)-jmin(1)+1
	 call gchange("Mmod",0)
         call getd1

c for the outboard half of the mesh:
	 ndnstream2=jmax(2)-jmin(2)+1
	 call gchange("Mmod",0)
         call getd2

      return
      end

c----------------------------------------------------------------------c

      subroutine getd1
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)dnstream1

      integer iref1,j,k
      real eps
      data eps /1./

      iref1=ixpoint(3,1)
      k=0
      do j=jmax(1),jmin(1),-1
         k=k+1
         rdnstream1(k)=cmeshx(iref1,j)
         zdnstream1(k)=cmeshy(iref1,j)
      enddo
c extrapolate end points slightly:
      rdnstream1(1)=rdnstream1(1)+eps*(rdnstream1(1)-rdnstream1(2))
      zdnstream1(1)=zdnstream1(1)+eps*(zdnstream1(1)-zdnstream1(2))
      rdnstream1(ndnstream1)=rdnstream1(ndnstream1)
     &          +eps*(rdnstream1(ndnstream1)-rdnstream1(ndnstream1-1))
      zdnstream1(ndnstream1)=zdnstream1(ndnstream1)
     &          +eps*(zdnstream1(ndnstream1)-zdnstream1(ndnstream1-1))

      return
      end

c----------------------------------------------------------------------c

      subroutine getd2
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)	# ixpoint,cmeshx,cmeshy
Use(Mmod)	# (r,z)dnstream2

      integer iref2,j,k
      real eps
      data eps /1./

      iref2=ixpoint(3,2)
      k=0
      do j=jmin(2),jmax(2)
         k=k+1
         rdnstream2(k)=cmeshx(iref2,j)
         zdnstream2(k)=cmeshy(iref2,j)
      enddo
c extrapolate end points slightly:
      rdnstream2(1)=rdnstream2(1)+eps*(rdnstream2(1)-rdnstream2(2))
      zdnstream2(1)=zdnstream2(1)+eps*(zdnstream2(1)-zdnstream2(2))
      rdnstream2(ndnstream2)=rdnstream2(ndnstream2)\
                +eps*(rdnstream2(ndnstream2)-rdnstream2(ndnstream2-1))
      zdnstream2(ndnstream2)=zdnstream2(ndnstream2)\
                +eps*(zdnstream2(ndnstream2)-zdnstream2(ndnstream2-1))

      return
      end

c     ------------------------------------------------------------------

      subroutine exleft
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)	#
Use(Curves)	# xcurveg,ycurveg,npointg
Use(Inmesh)
Use(Linkco)	# dxleft,ndxleft

c     local variables --
      integer n,j,npt,j1,npt1
      real dx,dy,rleft,slope,yc1

      external remark,xerrab

c     This subroutine extrapolates flux contours beyond the left
c     boundary of the EFIT computational domain, to permit construction
c     of orthogonal meshes for high-triangularity TCV plasmas.
c
c     NDXLEFT data points are added to each flux surface in radial
c     steps of size DXLEFT.
c
c     Only the inboard flux surfaces are extrapolated at present.
c     In the future this may be modified to include outboard surfaces.

      rleft=rgrid1		# initialize left radial limit

      do n=1,ndxleft
         rleft=rleft-dxleft	# increment left radial limit

c     Extrapolate the separatrix: 
         j=jsptrx(1)
         npt=npointg(j)
         dx=rleft-xcurveg(npt,j)
         slope=(ycurveg(npt,j)-ycurveg(npt-1,j))/
     .            (xcurveg(npt,j)-xcurveg(npt-1,j))
         dy=slope*dx
         xcurveg(npt+1,j)=xcurveg(npt,j)+dx
         ycurveg(npt+1,j)=ycurveg(npt,j)+dy
         if (npointg(j) .lt. npts) then
            npointg(j)=npointg(j)+1
         else
            call remark("***  error in subroutine exleft  ***")
            call remark("*** npointg(j) exceeds npts limit ***")
            call xerrab("")
         endif

c     Extrapolate SOL flux surfaces:
         do j=jsptrx(1)-1,jmin(1),-1
            npt=npointg(j)
            dx=rleft-xcurveg(npt,j)
            slope=(ycurveg(npt,j)-ycurveg(npt-1,j))/
     .            (xcurveg(npt,j)-xcurveg(npt-1,j))
            dy=slope*dx
            xcurveg(npt+1,j)=xcurveg(npt,j)+dx
            ycurveg(npt+1,j)=ycurveg(npt,j)+dy
c     Prevent flux surface intersection
            j1=j+1
            npt1=npointg(j1)
            yc1=ycurveg(npt1,j1)
            ycurveg(npt+1,j)=max(ycurveg(npt+1,j),yc1)
	    if (npointg(j) .lt. npts) then
	       npointg(j)=npointg(j)+1
	    else
	       call remark("***  error in subroutine exleft  ***")
	       call remark("*** npointg(j) exceeds npts limit ***")
	       call xerrab("")
	    endif
         enddo

c     Extrapolate private flux surfaces:
         do j=jsptrx(1)+1,jmax(1)
            npt=npointg(j)
            dx=rleft-xcurveg(npt,j)
            slope=(ycurveg(npt,j)-ycurveg(npt-1,j))/
     .            (xcurveg(npt,j)-xcurveg(npt-1,j))
            dy=slope*dx
            xcurveg(npt+1,j)=xcurveg(npt,j)+dx
            ycurveg(npt+1,j)=ycurveg(npt,j)+dy
c     Prevent flux surface intersection
            j1=j-1
            npt1=npointg(j1)
            yc1=ycurveg(npt1,j1)
            ycurveg(npt+1,j)=min(ycurveg(npt+1,j),yc1)
	    if (npointg(j) .lt. npts) then
	       npointg(j)=npointg(j)+1
	    else
	       call remark("***  error in subroutine exleft  ***")
	       call remark("*** npointg(j) exceeds npts limit ***")
	       call xerrab("")
	    endif
         enddo

      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine meshff (region)
      implicit none
      integer region
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)	# rseps,zseps
Use(Dimensions) # nptmp
Use(Curves)     # npointg,xcurveg,ycurveg
Use(Inmesh)     # ilmax
Use(Linkco)     # ixpoint,cmeshx,cmeshy
Use(Mmod)       # ntop,ntop1,ntop2,rtop,rtop1,rtop2,ztop,ztop1,ztop2,
                # nupstream,nupstream1,nupstream2,rupstream,rupstream1,
                # rupstream2,zupstream,zupstream1,zupstream2,
                # nplate,nplate1,nplate2,
                # rplate,rplate1,rplate2,zplate,zplate1,zplate2,
                # fuzzm,xmsh,ymsh,xcrv,ycrv,dsm,dsc,dsmesh,
                # nff,nff1,nff2,rff,rff1,rff2,zff,zff1,zff2,
                # cmeshx3,cmeshy3,dsmesh3,dsmeshff,
                # isxtform,iswtform,cwtffu,cwtffd,
                # wtff1,slpxff1,slpxffu1,slpxffd1,nxdff1,
                # wtff2,slpxff2,slpxffu2,slpxffd2,nxdff2
Use(Transfm)    # ijump

      real xtform, wtform
      external intersect2, xtform, wtform
      external remark, xerrab

c     Local variables --
      integer i,iend,ii,j,j1,j2,k,ierr
      integer i1crv,i2crv,i1msh,i2msh,icu,kcu,icp,kcp
      integer ict,kct,imt,icff,kcff,imff,kmff,nxdff
      integer na,ica,kca
      real xcu,ycu,xcp,ycp,xct,yct,xmt,ymt
      real dscp,frac,eps
      real xcff,ycff,xmff,ymff,dscff
      real wtff,slpxff,slpxffu,slpxffd
      real ra(3),za(3),xca,yca
      real t,t1,t2,t3,x1,x2,x3,wt

c     This subroutine modifies the mesh stored in (cmeshx3,cmeshy3)
c     by re-distributing mesh points along each flux surface so
c     that the "flame front" region has better spatial resolution.

      data eps/1.0e-06/		# extrapolation parameter

      j1=jmin(region)
      j2=jmax(region)
      iend=ilmax(region)

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
ccc
ccc   preparation procedures for subroutine meshff
ccc
c----------------------------------------------------------------------c

c     Define the top-of-mesh surface :
      if (region .eq. 1) then
         ntop=ntop1
         if (ntop .gt. nptmp) then
            call xerrab("dimension of r,ztop1 exceeds nptmp")
         endif
         do k=1,ntop
            rtop(k)=rtop1(k)
            ztop(k)=ztop1(k)
         enddo
      elseif (region .eq. 2) then
         ntop=ntop2
         if (ntop .gt. nptmp) then
            call xerrab("dimension of r,ztop2 exceeds nptmp")
         endif
         do k=1,ntop
            rtop(k)=rtop2(k)
            ztop(k)=ztop2(k)
         enddo
      endif

c     Define the upstream reference surface :
      if (region .eq. 1) then
         nupstream=nupstream1
         if (nupstream .gt. nptmp) then
            call xerrab("dimension of r,zupstream1 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream1(k)
            zupstream(k)=zupstream1(k)
         enddo
      elseif (region .eq. 2) then
         nupstream=nupstream2
         if (nupstream .gt. nptmp) then
            call xerrab("dimension of r,zupstream2 exceeds nptmp")
         endif
         do k=1,nupstream
            rupstream(k)=rupstream2(k)
            zupstream(k)=zupstream2(k)
         enddo
      endif

c     Define the flamefront surface :
      if (region .eq. 1) then
         nff=nff1
         if (nff .gt. nptmp) then
            call xerrab("dimension of r,zff1 exceeds nptmp")
         endif
         do k=1,nff
            rff(k)=rff1(k)
            zff(k)=zff1(k)
         enddo
      elseif (region .eq. 2) then
         nff=nff2
         if (nff .gt. nptmp) then
            call xerrab("dimension of r,zff2 exceeds nptmp")
         endif
         do k=1,nff
            rff(k)=rff2(k)
            zff(k)=zff2(k)
         enddo
      endif

c     Define some flamefront mesh control parameters:
      if (region .eq. 1) then
         nxdff=nxdff1
         wtff=wtff1
         slpxff=slpxff1
         slpxffu=slpxffu1
         slpxffd=slpxffd1
      elseif (region .eq. 2) then
         wtff=wtff2
         slpxff=slpxff2
         slpxffu=slpxffu2
         slpxffd=slpxffd2
         nxdff=nxdff2
      endif

c     Define the divertor plate surface :
      if (region .eq. 1) then
         nplate=nplate1
         if (nplate .gt. nptmp) then
            call xerrab("dimension of r,zplate1 exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate1(k)
            zplate(k)=zplate1(k)
         enddo
      elseif (region .eq. 2) then
         nplate=nplate2
         if (nplate .gt. nptmp) then
            call xerrab("dimension of r,zplate2 exceeds nptmp")
         endif
         do k=1,nplate
            rplate(k)=rplate2(k)
            zplate(k)=zplate2(k)
         enddo
      endif

c----------------------------------------------------------------------c

c     Determine flamefront position along the separatrix :
      j=jsptrx(region)
c     Define some temporary variables for (cmeshx,y) surface data
      i1msh=ixpoint(3,region)+1	# start index (just below x-point)
      i2msh=ilmax(region)	# end index (at plate)
      do i=i1msh,i2msh
         xmsh(i)=cmeshx3(i,j)
         ymsh(i)=cmeshy3(i,j)
      enddo
c     Intersection of (cmeshx,y) with flamefront surface
      call intersect2(rff,zff,1,nff,
     &                xmsh(i1msh:i2msh),ymsh(i1msh:i2msh),i1msh,i2msh,
     &                xmff,ymff,kmff,imff,fuzzm,ierr)
      if (ierr .ne. 0) then
         write (STDOUT,801) j
 801     format ("no flamefront intersection for cmeshx,y on j=",i3)
         call xerrab("")
      endif
c     Allow the user to specify the number of cells downstream
c     between the flamefront and divertor plate
      if (nxdff .gt. 0) then
         imff=max(ixpoint(3,region)+1,ilmax(region)-nxdff)
      endif

c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c     Loop over all flux surfaces:
      do j=j1,j2

c     Define some temporary variables for (xcurveg,ycurveg) surface data
         i1crv=ijump(j)+1	# start index
         i2crv=npointg(j)	# end index
         if (j .eq. jsptrx(region)) i1crv=1
         do i=i1crv,i2crv
            xcrv(i)=xcurveg(i,j)
            ycrv(i)=ycurveg(i,j)
         enddo

c     Define some temporary variables for (cmeshx,y) surface data
         i1msh=1			# start at top of SOL
         i2msh=ilmax(region)		# end index
         if ( ((region .eq. 1) .and. (j .ge. jsptrx(1)))   .or.
     &        ((region .eq. 2) .and. (j .le. jsptrx(2))) ) then
            i1msh=ixpoint(3,region)	# start at p.f. cut
         endif
         do i=i1msh,i2msh
            xmsh(i)=cmeshx3(i,j)
            ymsh(i)=cmeshy3(i,j)
         enddo
c     extrapolate upstream slightly from the initial point
         xmsh(i1msh)=xmsh(i1msh)+eps*(xmsh(i1msh)-xmsh(i1msh+1))
         ymsh(i1msh)=ymsh(i1msh)+eps*(ymsh(i1msh)-ymsh(i1msh+1))

c     Top-of-mesh for cmeshx,y:
         if ( ((region .eq. 1) .and. (j .ge. jsptrx(1)))   .or.
     &        ((region .eq. 2) .and. (j .le. jsptrx(2))) ) then
            imt=ixpoint(3,region)
         else
            imt=1
         endif
         xmt=xmsh(imt)
         ymt=ymsh(imt)
c        index of first point below top-of-mesh is imt+1

c     Intersection of (x,ycurveg) with top-of-mesh surface
         call intersect2(rtop,ztop,1,ntop,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xct,yct,kct,ict,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,884) j
 884        format ("no top-of-mesh intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
         if (j .eq. jsptrx(region)) then	# top is at x-point
            ict=ijump(j)+2	# there are 3 identical data points here
            xct=rseps
            yct=zseps
         endif
c        index of first point below top-of-mesh is ict+1

c     Intersection of (x,ycurveg) with upstream reference surface
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcu,ycu,kcu,icu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,885) j
 885        format ("no upstream intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below reference surface is icu+1
c        Ensure that upstream surface lies "below" top-of-mesh surface
         if (icu .lt. ict) then
            icu=ict
            xcu=xct
            ycu=yct
         endif

c     Intersection of (x,ycurveg) with flamefront surface
         call intersect2(rff,zff,1,nff,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcff,ycff,kcff,icff,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j
 886        format ("no flamefront intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below flamefront surface is icff+1
c        Ensure that flamefront surface lies "below" upstream surface
         if (icff .lt. icu) then
            icff=icu
            xcff=xcu
            ycff=ycu
         endif

c     Intersection of (x,ycurveg) with divertor plate
c     NOTE: start the search at the top-of-mesh index, ict, to avoid
c     a possible problem with multiple intersections (MER 95/09/05)
         call intersect2(rplate,zplate,1,nplate,
     &                    xcrv(ict:i2crv),ycrv(ict:i2crv),ict,i2crv,
     &                    xcp,ycp,kcp,icp,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,888) j
 888        format ("no plate intersection for x,ycurveg on j=",i3)
            write (STDOUT,889)
 889        format ("*** You may have to extend the plate end points")
            call xerrab("")
         endif
c        index of first point below divertor plate is icp+1

c     Distance from top-of-mesh to each data point along x,ycurveg
         xcrv(ict)=xct
         ycrv(ict)=yct
         dsc(ict)=0.
         do ii=ict+1,i2crv
            dsc(ii)=dsc(ii-1) +
     &                   sqrt( (xcrv(ii)-xcrv(ii-1))**2 +
     &                         (ycrv(ii)-ycrv(ii-1))**2 )
         enddo

c     Distance from top-of-mesh to flamefront
         dscff=dsc(icff)+sqrt( (xcff-xcrv(icff))**2 +
     &                       (ycff-ycrv(icff))**2 )

c     Distance from top-of-mesh to divertor plate
         dscp=dsc(icp)+sqrt( (xcp-xcrv(icp))**2 +
     &                       (ycp-ycrv(icp))**2 )

c     Ensure that flamefront is above divertor plate
         if (dscp .lt. dscff) then
            icff=icp
            xcff=xcp
            ycff=ycp
            dscff=dscp
         endif

c----------------------------------------------------------------------c
c     Compute new meshpoints (along x,ycurveg)
c----------------------------------------------------------------------c
c     First, the original distribution:
         dsmesh3(imt)=0.
         do i=imt+1,i2msh
c     Intersection of x,ycurveg with ith angle (orthogonal) surface
            na=3
            if (j .eq. jmin(region)) then
               ra(1)=2*cmeshx3(i,j)-cmeshx3(i,j+1)
               za(1)=2*cmeshy3(i,j)-cmeshy3(i,j+1)
            else
               ra(1)=cmeshx3(i,j-1)
               za(1)=cmeshy3(i,j-1)
            endif
            ra(2)=cmeshx3(i,j)
            za(2)=cmeshy3(i,j)
            if (j .eq. jmax(region)) then
               ra(3)=2*cmeshx3(i,j)-cmeshx3(i,j-1)
               za(3)=2*cmeshy3(i,j)-cmeshy3(i,j-1)
            else
               ra(3)=cmeshx3(i,j+1)
               za(3)=cmeshy3(i,j+1)
            endif
            call intersect2(ra,za,1,na,
     &                    xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                    xca,yca,kca,ica,fuzzm,ierr)
            if (ierr .ne. 0) then
               write (STDOUT,891) j,i
 891           format ("no angle intersection for x,ycurveg on j=",i3,
     &                 " at i=",i3)
               call xerrab("")
            endif
c     Index of first (x,ycrv) point below angle surface is ica+1
            dsmesh3(i)=dsc(ica)+sqrt( (xca-xcrv(ica))**2
     &                               +(yca-ycrv(ica))**2 )
         enddo
c----------------------------------------------------------------------c
c     Then, the flamefront distribution:
         t1=ixpoint(3,region)
         x1=dsmesh3(ixpoint(3,region))
         t2=imff
         x2=dscff
         t3=i2msh
         x3=dscp
         do i=imt,i2msh
            if (i .le. ixpoint(3,region)) then
               dsmeshff(i)=dsmesh3(i)
            else
               t=i
               dsmeshff(i)=xtform(t,t1,t2,t3,x1,x2,x3,
     &                              slpxffu,slpxff,slpxffd,isxtform)
            endif
         enddo
c----------------------------------------------------------------------c

c     Combine the 3 & ff meshes according to weight factor wtff:
	 do i=imt,i2msh
            t=i
            wt=wtform(t,t1,t2,t3,wtff,cwtffu,cwtffd,iswtform)
	    dsmesh(i)=(1-wt)*dsmesh3(i)+wt*dsmeshff(i)
	 enddo

c     Modified cmeshx,cmeshy values --
         do i=imt+1,i2msh
            do ii=ict+1,i2crv
               if (dsmesh(i) .le. dsc(ii)) then
c     modified point i lies between (x,ycurveg) points ii-1 and ii
                  frac=(dsmesh(i)-dsc(ii-1))/
     &                 (dsc(ii)  -dsc(ii-1))
                  cmeshx(i,j)=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                  cmeshy(i,j)=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                  break
               endif
            enddo
         enddo

      enddo	# end of loop over all j flux surfaces


      return
      end

c----------------------------------------------------------------------c

      real function wtform(t,t1,t2,t3,wtmax,gamma12,gamma23,iflag)
      implicit none
      integer iflag
      real t,t1,t2,t3,wtmax,gamma12,gamma23

      external remark,xerrab

      if (iflag .eq. 0) then
         wtform = wtmax
      elseif (iflag .eq. 1) then
         if (t .lt. t1) then
            wtform = 0.
         elseif ((t .ge. t1) .and. (t .le. t2)) then
            wtform = wtmax * ((t-t1)/(t2-t1))**gamma12
         elseif (t .gt. t2) then
            wtform = wtmax * ((t3-t)/(t3-t2))**gamma23
         endif
      else
         call xerrab("*** wtform: unknown iflag option ***")
      endif

      return
      end

c----------------------------------------------------------------------c

      real function xtform(t,t1,t2,t3,x1,x2,x3,
     &                       slopefac1,slopefac2,slopefac3,iflag)
      implicit none
      integer iflag
      real t,t1,t2,t3,x1,x2,x3,slopefac1,slopefac2,slopefac3

      real xtform1,xtform2,xtform3
      external xtform1,xtform2,xtform3,remark,xerrab

      if (iflag .eq. 1) then
         xtform = xtform1(t,t1,t2,t3,x1,x2,x3,
     &                      slopefac2)
      elseif (iflag .eq. 2) then
         xtform = xtform2(t,t1,t2,t3,x1,x2,x3,
     &                      slopefac1,slopefac2)
      elseif (iflag .eq. 3) then
         xtform = xtform3(t,t1,t2,t3,x1,x2,x3,
     &                      slopefac1,slopefac2,slopefac3)
      else
         call xerrab("*** xtform: unknown iflag option ***")
      endif

      return
      end

c----------------------------------------------------------------------c

      real function xtform1(t,t1,t2,t3,x1,x2,x3,slopefac2)
      implicit none
      real t,t1,t2,t3,x1,x2,x3,slopefac2

c     local variables --
      real alpha,beta,x2p

c     Compute parameters for analytic representation of the x-mesh;
c     On t1 < t < t2 use the rational function
c
c           x(t) = x1 + (x2-x1)*(t-t1)/((t2-t1)+alpha*(t2-t)
c
c     and on t2 < t < t3 use the rational function
c
c           x(t) = x2 + (x3-x2)*(t-t2)/((t3-t2)+beta*(t3-t))
c
c     where alpha and beta are chosen to give a specified slope
c     at t2.

c     Let the slope at t2 be some multiple of the average slope

      x2p = slopefac2*(x3-x1)/(t3-t1)

c     Then alpha and beta are given by

      alpha = (x2p*(t2-t1))/(x2-x1) - 1

      beta  = (x3-x2)/((t3-t2)*x2p) - 1

c     and x(t) is given by

      if (t .lt. t1) then
         xtform1=x1
      elseif ((t1 .le. t) .and. (t .lt. t2)) then
         xtform1=x1+(x2-x1)*(t-t1)/((t2-t1)+alpha*(t2-t))
      elseif ((t2 .le. t) .and. (t .le. t3)) then
         xtform1=x2+(x3-x2)*(t-t2)/((t3-t2)+beta*(t3-t))
      elseif (t3 .lt. t) then
         xtform1=x3
      endif

      return
      end

c----------------------------------------------------------------------c

      real function xtform2(t,t1,t2,t3,x1,x2,x3,slopefac1,slopefac2)
      implicit none
      real t,t1,t2,t3,x1,x2,x3,slopefac1,slopefac2

c     local variables --
      real alpha,beta,gamma,x1p,x2p

c     Compute parameters for analytic representation of the x-mesh;
c     On t1 < t < t2 use the rational function
c
c           x(t) = x1 + (x2-x1)*(t-t1)/((t2-t1)+alpha*(t2-t)+gamma*(t2-t)**2)
c
c     and on t2 < t < t3 use the rational function
c
c           x(t) = x2 + (x3-x2)*(t-t2)/((t3-t2)+beta*(t3-t))
c
c     where alpha, beta and gamma are chosen to give specified slopes
c     at t1 and t2.

c     Let the slope at t2 be some multiple of the average slope

      x2p = slopefac2*(x3-x1)/(t3-t1)

c     Then alpha and beta are given by

      alpha = (x2p*(t2-t1))/(x2-x1) - 1

      beta  = (x3-x2)/((t3-t2)*x2p) - 1

c     Let the slope at t1 be some multiple of the average slope

      x1p = slopefac1*(x3-x1)/(t3-t1)

c     Then gamma is given by

      gamma = (((x2-x1)/((t2-t1)*x1p))-1-alpha)/(t2-t1)

c     Finally, x(t) is given by

      if (t .lt. t1) then
         xtform2=x1
      elseif ((t1 .le. t) .and. (t .lt. t2)) then
         xtform2=x1+(x2-x1)*(t-t1)/((t2-t1)+alpha*(t2-t)+gamma*(t2-t)**2)
      elseif ((t2 .le. t) .and. (t .le. t3)) then
         xtform2=x2+(x3-x2)*(t-t2)/((t3-t2)+beta*(t3-t))
      elseif (t3 .lt. t) then
         xtform2=x3
      endif

      return
      end

c----------------------------------------------------------------------c

      real function xtform3(t,t1,t2,t3,x1,x2,x3,
     &                     slopefac1,slopefac2,slopefac3)
      implicit none
      real t,t1,t2,t3,x1,x2,x3,slopefac1,slopefac2,slopefac3

c     local variables --
      real alpha,beta,gamma,delta,x1p,x2p,x3p

c     Compute parameters for analytic representation of the x-mesh;
c     On t1 < t < t2 use the rational function
c
c           x(t) = x1 + (x2-x1)*(t-t1)/((t2-t1)+alpha*(t2-t)+gamma*(t2-t)**2)
c
c     and on t2 < t < t3 use the rational function
c
c           x(t) = x3 + (x2-x3)*(t-t3)/((t2-t3)+beta*(t2-t)+delta*(t2-t)**2)
c
c     where alpha, beta, gamma and delta are chosen to give specified
c     slopes at t1, t2 and t3.

c     Let the slope at t2 be some multiple of the average slope

      x2p = slopefac2*(x3-x1)/(t3-t1)

c     Then alpha and beta are given by

      alpha = (x2p*(t2-t1))/(x2-x1) - 1

      beta  = (x2p*(t2-t3))/(x2-x3) - 1

c     Let the slope at t1 be some multiple of the average slope

      x1p = slopefac1*(x3-x1)/(t3-t1)

c     Then gamma is given by

      gamma = (((x2-x1)/((t2-t1)*x1p))-1-alpha)/(t2-t1)

c     Let the slope at t3 be some multiple of the average slope

      x3p = slopefac3*(x3-x1)/(t3-t1)

c     Then delta is given by

      delta = (((x2-x3)/((t2-t3)*x3p))-1-beta)/(t2-t3)

c     Finally, x(t) is given by

      if (t .lt. t1) then
         xtform3=x1
      elseif ((t1 .le. t) .and. (t .lt. t2)) then
         xtform3=x1+(x2-x1)*(t-t1)/((t2-t1)+alpha*(t2-t)+gamma*(t2-t)**2)
      elseif ((t2 .le. t) .and. (t .le. t3)) then
         xtform3=x3+(x2-x3)*(t-t3)/((t2-t3)+beta*(t2-t)+delta*(t2-t)**2)
      elseif (t3 .lt. t) then
         xtform3=x3
      endif

      return
      end

c----------------------------------------------------------------------c

      subroutine getlim
      implicit none
Use(Limiter)    # nlimu,nptnma,nsplit1,nsplit2

c     Define the limiter surfaces that interface with the inboard and
c     outboard regions of the mesh

      if (nptnma .eq. 0) then
         call remark("***")
         call remark("getlim: limiter point nma not defined")
         call remark("***")
         call xerrab("")
      endif

c     For inboard region of the mesh:
      nsplit1=nptnma+1
      call gchange("Limiter",0)
      call getlim1

c     For outboard region of the mesh:
      nsplit2=nlimu-nptnma+2
      call gchange("Limiter",0)
      call getlim2

      return
      end

c----------------------------------------------------------------------c

      subroutine getlim1
      implicit none
Use(Limiter)    # nlimu,rlimu,zlimu,nptnma,rptnma,zptnma
                # rsplit1,zsplit1
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)

c     Define the limiter surface that interfaces with the inboard
c     region of the mesh

      integer n
      do n=1,nptnma
         rsplit1(n)=rlimu(n)
         zsplit1(n)=zlimu(n)
      enddo
      rsplit1(nptnma+1)=rmagx
      zsplit1(nptnma+1)=zmagx

      return
      end

c----------------------------------------------------------------------c

      subroutine getlim2
      implicit none
Use(Limiter)    # nlimu,rlimu,zlimu,nptnma,rptnma,zptnma
                # rsplit2,zsplit2
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)

c     Define the limiter surface that interfaces with the outboard
c     region of the mesh

      integer n
      do n=nlimu,nptnma,-1
         rsplit2(nlimu-n+1)=rlimu(n)
         zsplit2(nlimu-n+1)=zlimu(n)
      enddo
      rsplit2(nlimu-nptnma+2)=rmagx
      zsplit2(nlimu-nptnma+2)=zmagx

      return
      end

c----------------------------------------------------------------------c

      subroutine meshlim (region)
      implicit none
      integer region
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions) # nptmp
Use(Curves)     # npointg,xcurveg,ycurveg
Use(Inmesh)
Use(Linkco)     # ixpoint,cmeshx,cmeshy
Use(Mmod)       # ntop,rtop,ztop,ntop0,rtop0,ztop0,nbot,rbot,zbot,
                # cmeshx0,cmeshy0,fuzzm,
                # xmsh,ymsh,xcrv,ycrv,dsm,dsc,
                # dsmesh0,dsmesh1,dsmesh,wtm1
Use(Transfm)    # ijump
Use(Limiter)	# nsplit1,rsplit1,zsplit1,nsplit2,rsplit2,zsplit2
      external intersect2, remark, xerrab

c     Local variables --
      integer i,ii,j,j1,j2,k,ierr
      integer i1crv,i2crv,i1msh,i2msh,itop,ibot
      integer ict,kct,ict0,kct0,icb,kcb,imb
      integer na,ica,kca
      real xct,yct,xct0,yct0,xcb,ycb,xmb,ymb
      real frac,dsct,dsct0,dsmt0
      real ra(3),za(3),xca,yca

c     This subroutine modifies the original mesh by
c     re-distributing the mesh points along each flux surface.

      j1=jmin(region)
      j2=jmax(region)

c     Define the new (limiter interface) top-of-mesh surface :
      if (region .eq. 1) then
         ntop=nsplit1
         if (ntop .gt. nptmp) then
            call xerrab("meshlim:dimension of r,zsplit1 exceeds nptmp")
         endif
         do k=1,ntop
            rtop(k)=rsplit1(k)
            ztop(k)=zsplit1(k)
         enddo
      elseif (region .eq. 2) then
         ntop=nsplit2
         if (ntop .gt. nptmp) then
            call xerrab("meshlim:dimension of r,zsplit2 exceeds nptmp")
         endif
         do k=1,ntop
            rtop(k)=rsplit2(k)
            ztop(k)=zsplit2(k)
         enddo
      endif

c     Define the original top-of-mesh surface :
      itop=1
      ntop0=j2-j1+1
         if (ntop0 .gt. nptmp) then
            call xerrab("meshlim:dimension of r,ztop0 exceeds nptmp")
         endif
      do k=j1,j2
         rtop0(k-j1+1)=cmeshx0(itop,k)
         ztop0(k-j1+1)=cmeshy0(itop,k)
      enddo
c     extrapolate end points slightly beyond the the original mesh
      rtop0(1)=2*rtop0(1)-rtop0(2)
      ztop0(1)=2*ztop0(1)-ztop0(2)
      rtop0(ntop0)=2*rtop0(ntop0)-rtop0(ntop0-1)
      ztop0(ntop0)=2*ztop0(ntop0)-ztop0(ntop0-1)

c----------------------------------------------------------------------c
c     Define the lower boundary separating the modified and unmodified
c     regions of the mesh; usually, the surface at i=ixpoint(1,region)
      ibot=ixpoint(1,region)
      nbot=j2-j1+1
         if (nbot .gt. nptmp) then
            call xerrab("meshlim:dimension of r,zbot exceeds nptmp")
         endif
      do k=j1,j2
         rbot(k-j1+1)=cmeshx0(ibot,k)
         zbot(k-j1+1)=cmeshy0(ibot,k)
      enddo
c     extrapolate end points slightly beyond the the original mesh
      rbot(1)=2*rbot(1)-rbot(2)
      zbot(1)=2*zbot(1)-zbot(2)
      rbot(nbot)=2*rbot(nbot)-rbot(nbot-1)
      zbot(nbot)=2*zbot(nbot)-zbot(nbot-1)

c----------------------------------------------------------------------c
c     Loop over all flux surfaces:
      do j=j1,j2

c     Define some temporary variables for (xcurveg,ycurveg) surface data
         i1crv=1			# start index
         if (ijump(j) .ne. 0) then
            i2crv=ijump(j)		# end at last core data point
         else
            i2crv=npointg(j)		# end at last SOL data point
         endif
         if (j .eq. jsptrx(region)) i2crv=npointg(j)
         do i=i1crv,i2crv
            xcrv(i)=xcurveg(i,j)
            ycrv(i)=ycurveg(i,j)
         enddo

c     Define some temporary variables for (cmeshx,y) surface data
         i1msh=1			# start at top of SOL
         i2msh=ixpoint(1,region)	# end index (at x-point)
         do i=i1msh,i2msh
            xmsh(i)=cmeshx0(i,j)
            ymsh(i)=cmeshy0(i,j)
         enddo

c     bottom-of-mesh for cmeshx,y:
         imb=ibot
         xmb=xmsh(imb)
         ymb=ymsh(imb)

c     Distance from bottom to each meshpoint along cmeshx,y
         dsm(imb)=0.
         do i=imb-1,i1msh,-1
            dsm(i)=dsm(i+1) +
     &                   sqrt( (xmsh(i)-xmsh(i+1))**2 +
     &                         (ymsh(i)-ymsh(i+1))**2 )
         enddo

c     Distance from bottom to old top-of-mesh
         dsmt0=dsm(i1msh)

c----------------------------------------------------------------------c

c     Intersection of (x,ycurveg) with new top-of-mesh surface
         call intersect2(rtop,ztop,1,ntop,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xct,yct,kct,ict,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,884) j
 884        format ("no top-of-mesh intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below top-of-mesh is ict+1

c     Intersection of (x,ycurveg) with old top-of-mesh surface
         call intersect2(rtop0,ztop0,1,ntop0,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xct0,yct0,kct0,ict0,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,885) j
 885        format ("no top-of-mesh intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below top-of-mesh is ict0+1

c     Intersection of (x,ycurveg) with bottom surface
         call intersect2(rbot,zbot,1,nbot,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcb,ycb,kcb,icb,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j
 886        format ("no bot intersection for x,ycurveg on j=",i3)
            call xerrab("")
         endif
c        index of first point below bottom surface is icb+1


c     Distance from bottom-of-mesh to each data point along x,ycurveg
         xcrv(icb+1)=xcb
         ycrv(icb+1)=ycb
         dsc(icb+1)=0.
         do ii=icb,i1crv,-1
            dsc(ii)=dsc(ii+1) +
     &                   sqrt( (xcrv(ii)-xcrv(ii+1))**2 +
     &                         (ycrv(ii)-ycrv(ii+1))**2 )
         enddo

c     Distance from bottom to new top-of-mesh
         dsct=dsc(ict)-sqrt( (xct-xcrv(ict))**2 +
     &                       (yct-ycrv(ict))**2 )

c     Distance from bottom to old top-of-mesh
         dsct0=dsc(ict0)-sqrt( (xct0-xcrv(ict0))**2 +
     &                         (yct0-ycrv(ict0))**2 )

c----------------------------------------------------------------------c
c     Compute new meshpoints (along x,ycurveg) for each grid option
c----------------------------------------------------------------------c
c     First, the islimon=0 option:
         dsmesh0(imb)=0.
         do i=imb-1,i1msh,-1
c     Intersection of x,ycurveg with ith angle (orthogonal) surface
            na=3
            if (j .eq. jmin(region)) then
               ra(1)=2*cmeshx0(i,j)-cmeshx0(i,j+1)
               za(1)=2*cmeshy0(i,j)-cmeshy0(i,j+1)
            else
               ra(1)=cmeshx0(i,j-1)
               za(1)=cmeshy0(i,j-1)
            endif
            ra(2)=cmeshx0(i,j)
            za(2)=cmeshy0(i,j)
            if (j .eq. jmax(region)) then
               ra(3)=2*cmeshx0(i,j)-cmeshx0(i,j-1)
               za(3)=2*cmeshy0(i,j)-cmeshy0(i,j-1)
            else
               ra(3)=cmeshx0(i,j+1)
               za(3)=cmeshy0(i,j+1)
            endif
            call intersect2(ra,za,1,na,
     &                    xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                    xca,yca,kca,ica,fuzzm,ierr)
            if (ierr .ne. 0) then
               write (STDOUT,891) j,i
 891           format ("no angle intersection for x,ycurveg on j=",i3,
     &                 " at i=",i3)
               call xerrab("")
            endif
c     Index of first (x,ycrv) point below angle surface is ica+1
            dsmesh0(i)=dsc(ica)-sqrt( (xca-xcrv(ica))**2
     &                               +(yca-ycrv(ica))**2 )
         enddo
c----------------------------------------------------------------------c
c     Then, the non-zero islimon option:
         do i=imb,i1msh,-1
            dsmesh1(i)=dsmesh0(i)*dsct/dsct0
         enddo
c----------------------------------------------------------------------c
c     Combine options:  weighting will be implemented later
         do i=imb,i1msh,-1
            wtm1(i)=1.
            dsmesh(i)=(1-wtm1(i))*dsmesh0(i) + wtm1(i)*dsmesh1(i)
         enddo
c----------------------------------------------------------------------c

c     Modified cmeshx,cmeshy values --
            do i=imb-1,i1msh,-1
               do ii=icb,i1crv,-1
                  if (dsmesh(i) .le. dsc(ii)) then
c     modified point i lies between (x,ycurveg) points ii and ii+1
                     frac=(dsmesh(i)-dsc(ii+1))/
     &                    (dsc(ii)  -dsc(ii+1))
                     cmeshx(i,j)=xcrv(ii+1) + frac*(xcrv(ii)-xcrv(ii+1))
                     cmeshy(i,j)=ycrv(ii+1) + frac*(ycrv(ii)-ycrv(ii+1))
                     break
                  endif
               enddo
            enddo

      enddo	# end of loop over all j flux surfaces

      return
      end

c----------------------------------------------------------------------c

      subroutine setlimindex
      implicit none
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Linkco)  # cmeshx,y
Use(Inmesh)
Use(Share)   # nxleg,nxcore,nxomit,nxxpt,ix_lim,iy_lims
Use(Limiter) # dslims
      integer j,jj,j1
      real ds

# Set the poloidal index, ix_lim, of the limiter interface:
      ix_lim = nxleg(1,1) + nxcore(1,1) + 2*nxxpt - max(0, nxomit)

# Starting at the innermost core flux surface, search radially
# outward for the index of the first flux surface that has a poloidal
# gap at the inboard/outboard mesh interface.  This defines the 
# starting index iy=iy_lims for limiter recycling boundary conditions.
      do j=jmin(2),jmax(2)
         jj = j - jmin(2)
         j1 = jmax(1) - jj
         ds = sqrt( (cmeshx(1,j)-cmeshx(1,j1))**2
     .             +(cmeshy(1,j)-cmeshy(1,j1))**2 )
         if (ds .ge. dslims) then
            iy_lims = jj
            break
         endif
      enddo

      return
      end

c----------------------------------------------------------------------c

