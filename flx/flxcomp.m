c!include "../mppl.h"
c     ------------------------------------------------------------------
      subroutine flxgen
      implicit none
Use(Share)    #isfrc,ishalfm
Use(Dim)
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim
Use(Comflxgrd)	#
Use(Dimflx)	#nsearch
Use(Flxin)
Use(Inpf0)
Use(Inpf)
Use(Polflx)

      integer ns,kk,ll,j,n
      real dxf,dyf
      external refine, contours, remark

c     insert a sub-grid into the EFIT mesh --
      call refine

c     *** clear the working arrays
      do 20 kk=1,jdim
         ilast(kk)=0
         npoint(kk)=0
      do 20 ll=1,npts
         xcurve(ll,kk)=0.0
         ycurve(ll,kk)=0.0
   20 continue

c     initialize the contour indexing variables
      ncmin=0
      ncmax=0

c     For each search region, call the contouring routines
      do 500 ns=1,nsearch
         if (nsearch==4) then  # nycore is non-zero
            if (isfrc==1 .and. ns==2) then
               continue  # skip inboard half of p.f.
            elseif (isfrc==1 .and. ns==4) then
               continue  # skip outboard half of p.f.
            elseif (ishalfm==1 .and. (ns==1 .or. ns==2)) then
               continue  # skip inboard half
            else
               call contours(ns)
            endif
         elseif (nsearch==2) then  # nycore=0
            if (ishalfm==1 .and. ns==1) then
               continue  # skip inboard half
            else
               call contours(ns)
            endif
         endif
  500 continue

c     For core/p.f. contours, find the "jump" index:
      if (nycore .gt. 0) then
         do j=jsptrx(1)+1,jsptrx(2)-1
            ijumpf(j)=0
            do n=1,npoint(j)-1
               dxf=xcurve(n+1,j)-xcurve(n,j)
               dyf=ycurve(n+1,j)-ycurve(n,j)
               if (sqrt(dxf**2+dyf**2) .gt. dsjumpf) then
                  ijumpf(j)=n
                  break
               endif
            enddo
         enddo
      endif

ccc   call remark("***** subroutine flxgen completed")

      return
      end

c     ------------------------------------------------------------------

      subroutine refine
      implicit none
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Polflx)
      integer i,ix,iy,n
      real dx,dy,frac
      real B2VAhL
      external B2VAhL

#     This subroutine inserts a subgrid into the EFIT mesh


      frac=1/float(mrfac)
c******* first refine the mesh in x-direction *******

      ix=0
      do i=1,nxefit-1
         dx=xold(i+1)-xold(i)
         ix=ix+1
         x(ix)=xold(i)
         do n=1,mrfac-1
            ix=ix+1
            x(ix)=xold(i)+n*frac*dx
         enddo
      enddo
      x(nx4)=xold(nxefit)

c******* now refine the mesh in y-direction *******

      iy=0
      do i=1,nyefit-1
         dy=yold(i+1)-yold(i)
         iy=iy+1
         y(iy)=yold(i)
         do n=1,mrfac-1
            iy=iy+1
            y(iy)=yold(i)+n*frac*dy
         enddo
      enddo
      y(ny4)=yold(nyefit)

#     Evaluate the spline on the refined grid --
      do ix=1,nx4
         do iy=1,ny4
            f(ix,iy) = B2VAhL (x(ix), y(iy), 0, 0, xknot, yknot, nxefit,
     .                nyefit, kxord, kyord, bscoef, ldf, work, iflag)
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine contours(ns)
      implicit none
      integer k, nc, i, j, ii, n, ndown, iprev, nup, nleft, nright
      real tstval,tstvalh,tstvalv
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Dimflx)	#nsearch
Use(Inpf0)
Use(Inpf)
Use(Polflx)
      integer ns
      integer right,down,left,up
      parameter(right=1,down=2,left=3,up=4)

c.....Get parameters for contour search path ns
      ncmin=ncmin0(ns)
      ncmax=ncmax0(ns)
      ncontr=ncmax-ncmin+1
      do k=1,ncontr
         plflux(k)=plflux0(k,ns)
      enddo
      xminf=xminf0(ns)
      xmaxf=xmaxf0(ns)
      yminf=yminf0(ns)
      ymaxf=ymaxf0(ns)
c.....Start contour search at same point for all flux values
      do j=1,jdim
         iserch(j)=iserch0(ns)
         jserch(j)=jserch0(ns)
      enddo
      istepf=istepf0(ns)
      jstepf=jstepf0(ns)


      do 17 nc=ncmin,ncmax
         pfr=plflux(nc-ncmin+1)

c........Find starting point (iserch,jserch,leadir) on each flux contour

         if (istepf==0) then                 # vertical search (in p.f.)
            i=iserch(nc)
            j=jserch(nc)-jstepf
 91         continue
               j=j+jstepf
               if ( (j+jstepf .lt. jmins) .or. (jmaxs .lt. j+jstepf) ) then
                  write(*,*) "Error in loop 91 of subroutine contours:"
                  write(*,*) "  contour psi0 = ",(pfr-simagx)/(sibdry-simagx),
     .                 " not found in region ",ns
                  call xerrab("")
               endif
               tstvalv=(f(i,j+jstepf)-pfr)*(f(i,j)-pfr)
               if (tstvalv > 0) go to 91     # continue vertical stepping
# contour intersects vertical face between vertices (i,j) and (i,j+jstepf)
            iserch(nc)=i
            jserch(nc)=min(j,j+jstepf)  # min => go up, so leadir= 1 or 3
            if ( (ns==2) .and. nsearch==4 ) then
               # inboard half of p.f. region
               leadir(nc)=1      # search up; contour right, then left
            elseif (ns==4) then  # outboard half of p.f. region
               leadir(nc)=3      # search up; contour left, then right
            endif


         elseif (jstepf==0) then             # horizontal search (in SOL)
            i=iserch(nc)-istepf
            j=jserch(nc)
 92         continue
               i=i+istepf
               if ( (i+istepf .lt. imin ) .or. (imax  .lt. i+istepf) ) then
                  write(*,*) "Error in loop 92 of subroutine contours:"
                  write(*,*) "  contour psi0 = ",(pfr-simagx)/(sibdry-simagx),
     .                 " not found in region ",ns
                  call xerrab("")
               endif
               tstvalh=(f(i+istepf,j)-pfr)*(f(i,j)-pfr)
               if (tstvalh > 0) go to 92     # continue horizontal stepping
# contour intersects horizontal face between vertices (i+istepf,j) and (i,j)
            iserch(nc)=min(i,i+istepf)  # min => go right, so leadir= 2 or 4
            jserch(nc)=j
            if (ns==1) then      # inboard SOL
               leadir(nc)=4      # search right; contour up, then down
            elseif ( ns==3 .or. ((ns==2).and.(nsearch==2)) ) then
               # outboard SOL
               leadir(nc)=4      # search right; contour up, then down
            endif


         elseif (istepf*jstepf .ne. 0) then  # diagonal search (in p.f.)
            i=iserch(nc)-istepf
            j=jserch(nc)-jstepf
 93         continue
               i=i+istepf
               j=j+jstepf
               if ( (i+istepf .lt. imin ) .or. (imax  .lt. i+istepf) .or.
     .              (j+jstepf .lt. jmins) .or. (jmaxs .lt. j+jstepf) ) then
                  write(*,*) "Error in loop 93 of subroutine contours:"
                  write(*,*) "  contour psi0 = ",(pfr-simagx)/(sibdry-simagx),
     .                 " not found in region ",ns
                  call xerrab("")
               endif
               tstval=(f(i+istepf,j+jstepf)-pfr)*(f(i,j)-pfr)
               if (tstval .gt. 0.) go to 93  # continue diagonal stepping
# contour intersects diagonal between vertices (i,j) and (i+istepf,j+jstepf)
            tstvalh=(f(i+istepf,j)-pfr)*(f(i,j)-pfr)
            tstvalv=(f(i,j+jstepf)-pfr)*(f(i,j)-pfr)
            if (tstvalh .lt. 0.) then
# contour intersects horizontal face between vertices (i,j) and (i+istepf,j)
               iserch(nc)=min(i,i+istepf) # => go right, so leadir= 2 or 4
               jserch(nc)=j
               if ( (ns==2) .and. nsearch==4 ) then
                  # inboard half of p.f. region
                  leadir(nc)=2      # search right; contour down, then up
               elseif (ns==4) then  # outboard half of p.f. region
                  leadir(nc)=4      # search right; contour up, then down
               endif
            elseif (tstvalv .lt. 0) then
# contour intersects vertical face between vertices (i,j) and (i,j+jstepf)
               iserch(nc)=i
               jserch(nc)=min(j,j+jstepf) # => go up, so leadir= 1 or 3
               if ( (ns==2) .and. nsearch==4 ) then
                  # inboard half of p.f. region
                  leadir(nc)=1      # search up; contour right, then left
               elseif (ns==4) then  # outboard half of p.f. region
                  leadir(nc)=3      # search up; contour left, then right
               endif
            else
# contour intersects vertical face at i+istepf between j and j+jstepf
# and horizontal face at j+jstepf between i and i+istepf
               iserch(nc)=i+istepf
               jserch(nc)=min(j,j+jstepf) # => go up, so leadir= 1 or 3
               if ( (ns==2) .and. nsearch==4 ) then
                  # inboard half of p.f. region
                  leadir(nc)=1      # search up; contour right, then left
               elseif (ns==4) then  # outboard half of p.f. region
                  leadir(nc)=3      # search up; contour left, then right
               endif
            endif
         endif  # end if-test on istepf and jstepf


c........Go to appropriate section for search and subsequent contouring:
         go to (140,110,130,120,150,160) leadir(nc)

 110     continue
c     ** Next is horizontal search from left to right
c     ** and contouring, down then up.
c     ** Find the first point on the "pfr" contour
         i=iserch(nc)-1
         j=jserch(nc)
 100     continue  # top of loop that steps horizontally to the right
            i=i+1
            if (i+1 > imax) then
               write(*,*) "Error in section 110 of subroutine contours:"
               write(*,*) "flux contour not found on jserch ",
     .                       "between iserch and imax"
               write(*,*) "psi0 =   ",(pfr-simagx)/(sibdry-simagx)
               write(*,*) "jserch = ",jserch(nc)
               write(*,*) "iserch = ",iserch(nc)
               write(*,*) "imax =   ",imax
               call xerrab("")
            endif
            tstval=(f(i+1,j)-pfr)*(f(i,j)-pfr)
            if (tstval > 0.0) go to 100
c     ** We know that pfr is on horiz. face between vertices (i,j) and (i+1,j)
c     ** going down
         jstart=j
         istart=i
         call go(down,n,nc,ns)
         ndown=n
         npoint(nc)=npoint(nc)+ndown
         if (n > 0) then
c     ** copy xcn(1) to xcn(n) into xcurve in reverse order
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(n-ii+1)
               ycurve(iprev+ii,nc)=ycn(n-ii+1)
               xcn(n-ii+1)=0.
               ycn(n-ii+1)=0.
            enddo
            ilast(nc)=ilast(nc)+n
         endif
c     ** going up
         jstart=j-1
         istart=i
         call go(up,n,nc,ns)
         npoint(nc)=npoint(nc)+n
         if (n > 0) then
c     ** copy xcn(1) to xcn(n) into xcurve
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(ii)
               ycurve(iprev+ii,nc)=ycn(ii)
               xcn(ii)=0.
               ycn(ii)=0.
            enddo
            ilast(nc)=ilast(nc)+n
         endif
         go to 17

 120     continue
c     ** Next is horizontal search from left to right
c     ** and contouring, up then down
c     ** Find the first point on the "pfr" contour
         i=iserch(nc)-1
         j=jserch(nc)
 200     continue  # top of loop that steps horizontally to the right
            i=i+1
            if (i+1 > imax) then
               write(*,*) "Error in section 120 of subroutine contours:"
               write(*,*) "flux contour not found on jserch",
     .                       "between iserch and imax"
               write(*,*) "psi0 =   ",(pfr-simagx)/(sibdry-simagx)
               write(*,*) "jserch = ",jserch(nc)
               write(*,*) "iserch = ",iserch(nc)
               write(*,*) "imax =   ",imax
               call xerrab("")
            endif
            tstval=(f(i+1,j)-pfr)*(f(i,j)-pfr)
            if (tstval > 0.0) go to 200
c     ** We know that pfr is on horiz. face between vertices (i,j) and (i+1,j)
c     ** going up
         jstart=j-1
         istart=i
         call go(up,n,nc,ns)
         nup=n
         npoint(nc)=npoint(nc)+nup
         if (n > 0) then
c     ** copy xcn(1) to xcn(n) into xcurve in reverse order
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(n-ii+1)
               ycurve(iprev+ii,nc)=ycn(n-ii+1)
               xcn(n-ii+1)=0.
               ycn(n-ii+1)=0.
            enddo
            ilast(nc)=ilast(nc)+n
         endif
c     ** going down
         jstart=j
         istart=i
         call go(down,n,nc,ns)
         npoint(nc)=npoint(nc)+n
         if (n > 0) then
c     ** copy xcn(1) to xcn(n) into xcurve
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(ii)
               ycurve(iprev+ii,nc)=ycn(ii)
               xcn(ii)=0.
               ycn(ii)=0.
            enddo
            ilast(nc)=ilast(nc)+n
         endif
         go to 17

 130     continue
c     ** Next is vertical search upward
c     ** and contouring, left then right
c     ** Find the first point on the "pfr" contour
         i=iserch(nc)
         j=jserch(nc)-1
 300     continue  # top of loop that steps vertically upward
            j=j+1
            if (j+1 > jmaxs) then
               write(*,*) "Error in section 130 of subroutine contours:"
               write(*,*) "flux contour not found on iserch ",
     .                       "between jserch and jmaxs"
               write(*,*) "psi0 =   ",(pfr-simagx)/(sibdry-simagx)
               write(*,*) "iserch = ",iserch(nc)
               write(*,*) "jserch = ",jserch(nc)
               write(*,*) "jmaxs =  ",jmaxs
               call xerrab("")
            endif
            tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
            if (tstval > 0.0) go to 300
c     ** We know that pfr is on vert. face between vertices (i,j) and (i,j+1)
c     ** going left
         istart=i
         jstart=j
         call go(left,n,nc,ns)
         nleft=n
         npoint(nc)=npoint(nc)+nleft
         if (n > 0) then
c     ** copy xcn(1) to xcn(n) into xcurve in reverse order
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(n-ii+1)
               ycurve(iprev+ii,nc)=ycn(n-ii+1)
               xcn(n-ii+1)=0.0
               ycn(n-ii+1)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
c     ** going right
         istart=i-1
         jstart=j
         call go(right,n,nc,ns)
         npoint(nc)=npoint(nc)+n
         if (n > 0) then
c     ** copy xcn and ycn into xcurve and ycurve
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(ii)
               ycurve(iprev+ii,nc)=ycn(ii)
               xcn(ii)=0.0
               ycn(ii)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
         go to 17

 140     continue
c     ** Next is vertical search upward
c     ** and contouring, right then left
c     ** Find the first point on the "pfr" contour
         i=iserch(nc)
         j=jserch(nc)-1
 400     continue  # top of loop that steps vertically upward
            j=j+1
            if (j+1 > jmaxs) then
               write(*,*) "Error in section 140 of subroutine contours:"
               write(*,*) "flux contour not found on iserch ",
     .                       "between jserch and jmaxs"
               write(*,*) "psi0 =   ",(pfr-simagx)/(sibdry-simagx)
               write(*,*) "iserch = ",iserch(nc)
               write(*,*) "jserch = ",jserch(nc)
               write(*,*) "jmaxs =  ",jmaxs
               call xerrab("")
            endif
            tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
            if (tstval > 0.0) go to 400
c     ** We know that pfr is on vert. face between vertices (i,j) and (i,j+1)
c     ** going right
         istart=i-1
         jstart=j
         call go(right,n,nc,ns)
         nright=n
         npoint(nc)=npoint(nc)+nright
         if (n > 0) then
c     ** copy xcn and ycn into xcurve and ycurve in reverse order
            iprev=ilast(nc)
            do ii=1,n
               xcurve( iprev+ii,nc)=xcn(n-ii+1)
               ycurve( iprev+ii,nc)=ycn(n-ii+1)
               xcn(n-ii+1)=0.0
               ycn(n-ii+1)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
c     ** going left
         istart=i
c        jstart=j
         call go(left,n,nc,ns)
         npoint(nc)=npoint(nc)+n
         if (n > 0) then
c     ** copy xcn(1) to xcn(n) into xcurve
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(ii)
               ycurve(iprev+ii,nc)=ycn(ii)
               xcn(ii)=0.0
               ycn(ii)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
         go to 17

 150     continue
c     ** Next is vertical search downward
c     ** and contouring, left then right
c     ** Find the first point on the "pfr" contour
         i=iserch(nc)
         j=jserch(nc)
 500     continue  # top of loop that steps vertically downward
            j=j-1
            if (j < jmins) then
               write(*,*) "Error in section 150 of subroutine contours:"
               write(*,*) "flux contour not found on iserch ",
     .                       "between jmins and jserch"
               write(*,*) "psi0 =   ",(pfr-simagx)/(sibdry-simagx)
               write(*,*) "iserch = ",iserch(nc)
               write(*,*) "jmins =  ",jmins
               write(*,*) "jserch = ",jserch(nc)
               call xerrab("")
            endif
            tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
            if (tstval > 0.0) go to 500
c     ** We know that pfr is on vert. face between vertices (i,j) and (i,j+1)
c     ** going left
         istart=i
         jstart=j
         call go(left,n,nc,ns)
         nleft=n
         npoint(nc)=npoint(nc)+nleft
         if (n > 0) then
c     ** copy xcn and ycn into xcurve and ycurve in reverse order
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(n-ii+1)
               ycurve(iprev+ii,nc)=ycn(n-ii+1)
               xcn(n-ii+1)=0.0
               ycn(n-ii+1)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
c     ** going right
         istart=i-1
         jstart=j
         call go(right,n,nc,ns)
         npoint(nc)=npoint(nc)+n
         if (n > 0) then
c     ** copy xcn and ycn into xcurve and ycurve
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(ii)
               ycurve(iprev+ii,nc)=ycn(ii)
               xcn(ii)=0.0
               ycn(ii)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
         go to 17

 160     continue
c     ** Next is vertical search downward
c     ** and contouring, right then left
c     ** Find the first point on the "pfr" contour
         i=iserch(nc)
         j=jserch(nc)
 600     continue  # top of loop that steps vertically downward
         j=j-1
            if (j < jmins) then
               write(*,*) "Error in section 160 of subroutine contours:"
               write(*,*) "flux contour not found on iserch ",
     .                       "between jmins and jserch"
               write(*,*) "psi0 =   ",(pfr-simagx)/(sibdry-simagx)
               write(*,*) "iserch = ",iserch(nc)
               write(*,*) "jmins =  ",jmins
               write(*,*) "jserch = ",jserch(nc)
               call xerrab("")
            endif
            tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
            if (tstval > 0.0) go to 600
c     ** We know that pfr is on vert. face between vertices (i,j) and (i,j+1)
c     ** going right
         istart=i-1
         jstart=j
         call go(right,n,nc,ns)
         nright=n
         npoint(nc)=npoint(nc)+nright
         if (n > 0) then
c    *** copy xcn and ycn into xcurve and ycurve in reverse order
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(n-ii+1)
               ycurve(iprev+ii,nc)=ycn(n-ii+1)
               xcn(n-ii+1)=0.0
               ycn(n-ii+1)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
c     ** going left
         istart=i
         jstart=j
         call go(left,n,nc,ns)
         npoint(nc)=npoint(nc)+n
         if (n > 0) then
c     ** copy xcn and ycn into xcurve and ycurve
            iprev=ilast(nc)
            do ii=1,n
               xcurve(iprev+ii,nc)=xcn(ii)
               ycurve(iprev+ii,nc)=ycn(ii)
               xcn(ii)=0.0
               ycn(ii)=0.0
            enddo
            ilast(nc)=ilast(nc)+n
         endif
         go to 17

 17   continue  # end of do-loop over flux contours in region ns

      return
      end

c     ------------------------------------------------------------------

      subroutine go(idir,n,jj,ns)
      implicit none
      integer idir, n, jj, ns, i, j, nmx, ncod
      real xnew, dx, dy, tstval, dy1, rat, dx1
Use(Dim)
Use(Dimflxgrd)	#jdim,npts
Use(Comflxgrd)	#
Use(Flxin)	#istchkon
Use(Polflx)
      external theta_ok
      logical theta_ok

c     This subroutine computes data points on the jjth flux contour
c     in subregion ns, and stores them in the (xcn,ycn) arrays.
c
c     On entry, common block variables (istart,jstart) are the indices of
c     a vertex adjacent to the jjth flux contour, pfr, in subregion ns;
c     (istart,jstart) refer to (i,j) indices on the refined-EFIT rectangular
c     mesh, so x(i) and y(j) are the (R,Z) coord's of this vertex.

c     On entry, idir is the flag that indicates which face of the (i,j)
c     "cell" intersects the pfr contour:
c     idir=1    # go right; pfr is on vert face between (i+1,j) and (i+1,j+1)
c     idir=2    # go down; pfr is on horiz face between (i,j) and (i+1,j) 
c     idir=3    # go left; pfr is on vert face between (i,j) and (i,j+1)
c     idir=4    # go up; pfr is on horiz face between (i,j+1) and (i+1,j+1)

c     On exit, the n data points are stored in the first n entries of the
c     common block arrays (xcn,ycn).

      i=istart
      j=jstart
      n=0
      nmx=npts-npoint(jj)
      xnew=x(i)

      go to(1,2,3,4)idir
c
    1 continue
c     *** going right
      if (n.ge.nmx) then
         write (STDOUT,901) npts, jj
         return
      endif
         ncod=1
         i=i+1
         if (i .ge. imax) return        # MER 12/13/90
         dx=x(i+1)-x(i)
         dy=y(j+1)-y(j)
c     *** check on upright
         tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
         if(tstval.gt.0.0) go to 11
c     *** find coordinates
c         xold=xnew
         xnew=x(i)
c     if(xnew.gt.xold) return
         n=n+1
         dy1=dy*(pfr-f(i,j))/(f(i,j+1)-f(i,j))
         xcn(n)=xnew
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if(i.eq.imax) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 12
c
   11 continue
         write(STDOUT,900) ncod
      return
c
c     ***
   12 continue
c     *** check diagonal
         tstval=(f(i+1,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 4
c     *** find coordinates
         rat=(pfr-f(i,j))/(f(i+1,j+1)-f(i,j))
         dx1=dx*rat
         dy1=dy*rat
         n=n+1
         xcn(n)=x(i)+dx1
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 14
c
c     ***
   14 continue
c     *** check horizontal
         tstval=(f(i+1,j)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 1
c     *** find coordinates
         dx1=dx*(pfr-f(i,j))/(f(i+1,j)-f(i,j))
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.lt.xold) return
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if(j.eq.jmins) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 2
c
c     ***************************
    2 continue
c     *** going down
      if (n.ge.nmx) then
         write (STDOUT,901) npts, jj
         return
      endif
         ncod=2
         j=j-1
         if (j .le. jmins) return        # MER 12/13/90
         dx=x(i+1)-x(i)
         dy=y(j+1)-y(j)
c     *** check diagonal
         tstval=(f(i+1,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 22
c     *** find coordinates
         rat=(pfr-f(i,j))/(f(i+1,j+1)-f(i,j))
         dx1=dx*rat
         dy1=dy*rat
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.gt.xold) return
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 21
c
c     ***
   21 continue
c     *** check horizontal
         tstval=(f(i+1,j)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 1
c     *** find coordinates
         dx1=dx*(pfr-f(i,j))/(f(i+1,j)-f(i,j))
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.gt.xold) return
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if(j.eq.jmins) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 2
c
c     ***
   22 continue
c     *** check vertical
         tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) return
c     *** find coordinates
c         xold=xnew
         xnew=x(i)
c     if(xnew.gt.xold) return
         n=n+1
         dy1=dy*(pfr-f(i,j))/(f(i,j+1)-f(i,j))
         xcn(n)=xnew
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if(i.eq.imin) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 3
c
c     ******************************
    3 continue
c     *** going left
      if (n.ge.nmx) then
         write (STDOUT,901) npts, jj
         return
      endif
         ncod=3
         i=i-1
         if (i .le. imin) return        # MER 12/13/90
         dx=x(i+1)-x(i)
         dy=y(j+1)-y(j)
c     *** check diagonal
         tstval=(f(i+1,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 32
c     *** find coordinates
         rat=(pfr-f(i,j))/(f(i+1,j+1)-f(i,j))
         dx1=dx*rat
         dy1=dy*rat
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.gt.xold) return
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 31
c
c     ***
   31 continue
c     *** check vertical
         tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 4
c     *** find coordinates
c         xold=xnew
         xnew=x(i)
         dy1=dy*(pfr-f(i,j))/(f(i,j+1)-f(i,j))
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if(i.eq.imin) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 3
c
c     ***
   32 continue
c     *** check horizontal
         tstval=(f(i+1,j)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) return
c     *** ind coordinates
         dx1=dx*(pfr-f(i,j))/(f(i+1,j)-f(i,j))
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.gt.xold) return
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if(j.eq.jmins) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 2
c
c     ***************************
    4 continue
c     *** going up
      if (n.ge.nmx) then
         write (STDOUT,901) npts, jj
         return
      endif
         ncod=4
         j=j+1
         if (j .ge. jmaxs) return        # MER 12/13/90
         dx=x(i+1)-x(i)
         dy=y(j+1)-y(j)
c     *** check horizontal
         tstval=(f(i+1,j)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) return
c     *** find coordinates
         dx1=dx*(pfr-f(i,j))/(f(i+1,j)-f(i,j))
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.gt.xold) return
         n=n+1
         xcn(n)=xnew
         ycn(n)=y(j)
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if(j.eq.jmaxs) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 41
c
c     ***
   41 continue
c     *** check diagonal
         tstval=(f(i+1,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 1
c     *** find coordinates
         rat=(pfr-f(i,j))/(f(i+1,j+1)-f(i,j))
         dx1=dx*rat
         dy1=dy*rat
c         xold=xnew
         xnew=x(i)+dx1
c     if(xnew.gt.xold) return
         n=n+1
         xcn(n)=x(i)+dx1
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 42
c
c     ***
   42 continue
c     *** check vertical
         tstval=(f(i,j+1)-pfr)*(f(i,j)-pfr)
      if(tstval.gt.0.0) go to 4
c         xold=xnew
         xnew=x(i)
c     if(xnew.gt.xold) return
         dy1=dy*(pfr-f(i,j))/(f(i,j+1)-f(i,j))
         n=n+1
         xcn(n)=x(i)
         ycn(n)=y(j)+dy1
      if((xcn(n).gt.xmaxf).or.(xcn(n).lt.xminf)) return
      if((ycn(n).gt.ymaxf).or.(ycn(n).lt.yminf)) return
      if(i.eq.imin) return
      if (istchkon .eq. 1) then                     # MER 97/03/24
         if ( .not. theta_ok(xcn(n),ycn(n),ns) ) return
      endif
      go to 3
c
  900 format(' go - section = ',i2)
  901 format(' more than npts = ',i4,
     &       ' data points on contour segment ',i3)
c     ***
  999 continue
      return
      end


c----------------------------------------------------------------------c

      logical function theta_ok (r, z, n)
      implicit none
      integer n 
      real r, z

c     This function checks to see that the point (r,z) lies within
c     an acceptable polar angle for various flux-contouring regions.
c     The origin of the coordinate sytem for calculating the polar
c     angle is the magnetic axis.  By setting limits on the polar
c     angle we avoid tracing endless loops around closed core plasma
c     flux contours.

Use(Dim)
Use(Dimflxgrd)	#jdim,nxefit,nyefit,nlim
Use(Comflxgrd)	#
Use(Aeqflxgrd)	#
Use(Dimflx)	#nsearch
Use(Flxin)      #theta1fac,theta2fac

c     local variables
      real theta, thetap, thetamaxp, theta1, theta2, pi, twopi

c     initialize
      pi = 4.*atan2(1.,1.)
      twopi = 2*pi
      theta_ok = .false.
      theta = atan2 (z-zmagx, r-rmagx)    # -pi < theta < +pi

c----------------------------------------------------------------------c
      if (nsearch .eq. 2) then     # SOL only; no core region
c----------------------------------------------------------------------c
         if (n .eq. 1) then        # region 1: inboard SOL
            # in rotated frame where thetamin(1) = 0 (see subr inflx1)
            thetap = theta - thetamin(1)
            if (thetap .lt. 0) thetap = thetap + twopi
            if (thetap .gt. twopi) thetap = thetap - twopi
            thetamaxp = thetamax(1) - thetamin(1)
            if (thetamaxp .lt. 0) thetamaxp = thetamaxp + twopi
            if (thetamaxp .gt. twopi) thetamaxp = thetamaxp - twopi
            if (thetap .lt. thetamaxp) theta_ok = .true.
         elseif (n .eq. 2) then    # region 2: outboard SOL
            # in rotated frame where thetamin(2) = 0 (see subr inflx1)
            thetap = theta - thetamin(2)
            if (thetap .lt. 0) thetap = thetap + twopi
            if (thetap .gt. twopi) thetap = thetap - twopi
            thetamaxp = thetamax(2) - thetamin(2)
            if (thetamaxp .lt. 0) thetamaxp = thetamaxp + twopi
            if (thetamaxp .gt. twopi) thetamaxp = thetamaxp - twopi
            if (thetap .lt. thetamaxp) theta_ok = .true.
         else
            call remark("*** ")
            call remark("*** function theta_ok: illegal argument n")
            call remark("*** ")
            call xerrab("")
         endif
c----------------------------------------------------------------------c
      elseif (nsearch .eq. 4) then # SOL + finite core region
c----------------------------------------------------------------------c
         if (n .eq. 1) then        # region 1: inboard SOL & core
            # in rotated frame where thetamin(1) = 0 (see subr inflx1)
            thetap = theta - thetamin(1)
            if (thetap .lt. 0) thetap = thetap + twopi
            if (thetap .gt. twopi) thetap = thetap - twopi
            thetamaxp = thetamax(1) - thetamin(1)
            if (thetamaxp .lt. 0) thetamaxp = thetamaxp + twopi
            if (thetamaxp .gt. twopi) thetamaxp = thetamaxp - twopi
            if (thetap .lt. thetamaxp) theta_ok = .true.
         elseif (n .eq. 2) then    # region 2: inboard private flux
# MER (09 Feb 2015) Avoid discontinuity in theta=atan2(y,x) near +/-pi
# by transforming all theta's such that 0 <= theta < +twopi
            theta1 = theta1fac*pi
            theta2 =  twopi + thetax + dtheta_overlap_pf(1)
	    if (theta .lt. 0) theta = theta + twopi
            if ((theta1 .lt. theta) .and. (theta .lt. theta2)) 
     .         theta_ok = .true.
         elseif (n .eq. 3) then    # region 3: outboard SOL & core
            # in rotated frame where thetamin(2) = 0 (see subr inflx1)
            thetap = theta - thetamin(2)
            if (thetap .lt. 0) thetap = thetap + twopi
            if (thetap .gt. twopi) thetap = thetap - twopi
            thetamaxp = thetamax(2) - thetamin(2)
            if (thetamaxp .lt. 0) thetamaxp = thetamaxp + twopi
            if (thetamaxp .gt. twopi) thetamaxp = thetamaxp - twopi
            if (thetap .lt. thetamaxp) theta_ok = .true.
         elseif (n .eq. 4) then    # region 4: outboard private flux
            theta1 =  thetax - dtheta_overlap_pf(2)
# MER (09 Feb 2015) Insert user-defined multiplicative factor in theta2:
            theta2 = theta2fac*pi
            if ((theta1 .lt. theta) .and. (theta .lt. theta2)) 
     .         theta_ok = .true.
         else
            call remark("*** ")
            call remark("*** function theta_ok: illegal argument n")
            call remark("*** ")
            call xerrab("")
         endif
c----------------------------------------------------------------------c
      else
c----------------------------------------------------------------------c
         call remark("*** ")
         call remark("*** function theta_ok: nsearch must be 2 or 4")
         call remark("*** ")
         call xerrab("")
c----------------------------------------------------------------------c
      endif
c----------------------------------------------------------------------c

      return
      end




