c!include "../mppl.h"
c     ------------------------------------------------------------------

      subroutine inflx
      implicit none
Use(Share)            # geometry,nycore,nysol,nyout,isfrc
Use(Dim)              # nxpt
Use(Xpoint_indices)
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Dimflx)	#nsearch
Use(Polflx)           # mrfac
Use(Flxin)            # psi0sep1,psi0sep2,iseqdskr

c     This subroutine allocates arrays, sets some indices and reads the
c     eqdsk files.

c     Set the number of x-points in the mesh:
      if (geometry=="dnull") then
         nxpt = 2
      else
         nxpt = 1
      endif
      call gchange("Xpoint_indices",0)
 
# Read EFIT data from aeqdsk and neqdsk files
      if (mdsefit .eq. 0) then 
         call readefit
      endif
# Set contour data array dimensions
      npts = 4*mrfac*(nxefit+nyefit)
      if (nycore(igrid) .eq. 0) then
         nsearch = 2
      else
         nsearch = 4
      endif
      if (kymesh==0 .or. geometry=="dnull") then
         nym = nycore(igrid) + nysol(igrid) + nyout(igrid)
      else
         nym = nycore(igrid) + nysol(igrid)
      endif
      jdim = 2*(nym+1) + 1       # includes magnetic axis point
# Allot storage space for arrays
      call gchange("Comflxgrd",0)
      call gallot("Inpf0",0)
      call gallot("Inpf",0)
      call gallot("Polflx",0)
      call gchange("Flxin",0)    # save existing psitop,psibot when kymesh=0

# Set characteristic indices in radial direction
      jmin(1) = 1
      if (kymesh==0 .or. geometry=="dnull") then
c
c     psi0sep1 refers to the lower x-point
c     psi0sep2 refers to the upper x-point
cccMER NOTE: psi0sep1 and psi0sep2 must be supplied by the user; these
cccMER       both default to zero, in which case we get the usual
cccMER       results based on nycore and nysol (nyout=0 is default).
c
         if (iseqdskr==0) then  # we are doing the lower half mesh
            if (psi0sep1 .gt. psi0sep2) then
               jsptrx(1) = jmin(1) + nyout(igrid)
               jmax(1) = jsptrx(1) + nysol(igrid) + nycore(igrid)
               jaxis = jmax(1) + 1
               jmin(2) = jaxis + 1
               jsptrx(2) = jmin(2) + nycore(igrid) + nysol(igrid)
               jmax(2) = jsptrx(2) + nyout(igrid)
            else
               jsptrx(1) = jmin(1) + nyout(igrid) + nysol(igrid)
               jmax(1) = jsptrx(1) + nycore(igrid)
               jaxis = jmax(1) + 1
               jmin(2) = jaxis + 1
               jsptrx(2) = jmin(2) + nycore(igrid)
               jmax(2) = jsptrx(2) + nysol(igrid) + nyout(igrid)
            endif
         else                   # we are doing the upper half mesh
            if (psi0sep2 .gt. psi0sep1) then
               jsptrx(1) = jmin(1) + nyout(igrid)
               jmax(1) = jsptrx(1) + nysol(igrid) + nycore(igrid)
               jaxis = jmax(1) + 1
               jmin(2) = jaxis + 1
               jsptrx(2) = jmin(2) + nycore(igrid) + nysol(igrid)
               jmax(2) = jsptrx(2) + nyout(igrid)
            else
               jsptrx(1) = jmin(1) + nyout(igrid) + nysol(igrid)
               jmax(1) = jsptrx(1) + nycore(igrid)
               jaxis = jmax(1) + 1
               jmin(2) = jaxis + 1
               jsptrx(2) = jmin(2) + nycore(igrid)
               jmax(2) = jsptrx(2) + nysol(igrid) + nyout(igrid)
            endif
         endif  # end if-test on iseqdskr
      else  # use only nycore and nysol for kymesh > 0
         jsptrx(1) = jmin(1) + nysol(igrid)
         jmax(1) = jsptrx(1) + nycore(igrid)
         jaxis = jmax(1) + 1	# There is a dummy surface
				# between jmax(1) and jmin(2)
				# that represents the magnetic axis.
         jmin(2) = jaxis + 1
         jsptrx(2) = jmin(2) + nycore(igrid)
         jmax(2) = jsptrx(2) + nysol(igrid)
      endif  # end if-test on kymesh and geometry

# Special coding for field-reversed configuration
      if (isfrc==1) jmin(1)=jsptrx(1)

# Define flux contour values
      call inflx1

      return
      end

c     ------------------------------------------------------------------

      subroutine readefit
      implicit none

c     read in aeqdsk file from EFIT code --
      call aeqdsk

c     read in neqdsk file from EFIT code --
      call neqdsk

      call procefit

      return
      end

c     ------------------------------------------------------------------

      subroutine procefit
      implicit none
      integer i, j
      real dxdim, dydim
Use(Dim)
Use(Dimflxgrd)	#jdim,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Share)
Use(Aeqflxgrd)
Use(Flxin)    # iseqdskr
      real B2VAhL
      external B2INhT, B2VAhL

      if (geometry=="uppersn" .or. iseqdskr==1) call convert_eqdsk

c...  Now can pass simagxs and sibdrys to Share group for use in other pkgs
      simagxs = simagx
      sibdrys = sibdry

c   rcentr = major radius of the point where btoroidal = bcentr
c   rgrid1 = distance between major axis and the left side of
c   EFIT computational grid

c     Reconstruct the EFIT mesh -- with shift of vertical coordinate
      dxdim=xdim/(nxefit-1)
      dydim=zdim/(nyefit-1)
      do 10 i=1,nxefit
         xold(i)=(i-1)*dxdim+rgrid1
  10  continue
      do 20 j=1,nyefit
         yold(j)=(j-1)*dydim
 20   continue

c     Construct the 2-d spline interpolant for the flux
      call s2copy (nxefit,nyefit,fold,1,nxefit,bscoef,1,nxefit) 
      ldf = nxefit
      iflag = 1	# let B2INhT define the knots
      call B2INhT (xold, nxefit, yold, nyefit, kxord, kyord,
     .             xknot, yknot, bscoef, ldf, work, iflag)

c     shift vertical coordinate data from EFIT so as to make z=0 at the
c     bottom of the mesh (as in definition of yold above) --
      zshift = zdim/2. - zmid
      do 8 i=1,nbdry
         zbdry(i) = zbdry(i) + zshift
    8 continue
      do 9 i=1,nlim
         ylim(i)  = ylim(i)  + zshift
    9 continue
      zmid    = zmid   + zshift
      zmagx   = zmagx  + zshift
      zseps   = zseps  + zshift
      zseps2  = zseps2 + zshift
      zvsin   = zvsin  + zshift
      zvsout  = zvsout + zshift

      rseps1 = rseps
      zseps1 = zseps
c     compute the flux value, sibdry1, at the lower x-point --
      if (rseps1 .gt. 0) then
         sibdry1 = B2VAhL(rseps1, zseps1, 0, 0, xknot, yknot, nxefit,
     .              nyefit, kxord, kyord, bscoef, ldf, work, iflag)
      endif

c     compute the flux value, sibdry2, at the upper x-point --
      if (rseps2 .gt. 0) then
         sibdry2 = B2VAhL(rseps2, zseps2, 0, 0, xknot, yknot, nxefit,
     .              nyefit, kxord, kyord, bscoef, ldf, work, iflag)
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine aeqdsk
      implicit none
      integer ktime1, jflag, lflag
Use(Aeqflxgrd)   # vmonth,vday,vyear
              # eshot,etime
              # mco2v,mco2r
              # rseps,zseps,rseps2,zseps2
              # rvsin,zvsin,rvsout,zvsout
              # nsilop,magpri,nfcoil,nesum

# The array dimensioning parameters nsilop, magpri, nfcoil and nesum
# are defaulted to DIII-D values as of 10/15/92.
# Change to magpri=29 for old aeqdsk files (e.g., APS '88 work)
# Change to nsilop=26, magpri=26, nfcoil=13 and nesum=3 for C-MOD
#
# The format for the aeqdsk and neqdsk files is given in the file
# efit:weqdsk.for on the GA VAX's.  One also needs the parameter
# and dimensioning declarations contained in files efit:parmd.for
# and efit:ecomd1.for and efit:ecomd2.for.
#
      character*10 uday
      character*3 limloc
      character*3 qmflag
      integer iunit, ios
      integer efitvers
      external efitvers

c ***************** read the output of the EFIT code *****************
      data iunit /55/
      open (iunit, file=aeqdskfname, form='formatted', iostat=ios,
     .      status='old')
      if (ios .ne. 0) then
         call xerrab("**** aeqdsk file not found")
      endif

        read (iunit,1056) uday,vmonth,vday,vyear
c       Convert from 2-digit to 4-digit year for old EFIT --
        if ( (75 .le. vyear) .and. (vyear .le. 99) ) then
           vyear = 1900 + vyear
        endif
c       Special case is the widely-used ITER aeqdsk file, aeqibm1,
c       which contains a faulty day and date format --
        if (uday=="06.1./94  ") then
           vday = 12
           vmonth = 6
           vyear = 1994
        endif
        read (iunit,*) eshot,ktime1
        read (iunit,1040) etime
c     The following allows for the possibility of reading 2 additional
c     flags from newer versions of EFIT, but the read statement has
c     not yet been changed here (MER 97/06/03).
      if (efitvers(vmonth,vday,vyear) .ge. 1) then
        read (iunit,1060) etime,jflag,lflag,limloc,mco2v,mco2r,qmflag
      else
        read (iunit,1060) etime,jflag,lflag,limloc,mco2v,mco2r,qmflag
      endif

        call gchange("Aeqflxgrd",0)

        call aeqdsk1(iunit)

      close (iunit)

c     convert coordinate units from [cm] to [m] --
      rseps  = .01 * rseps
      zseps  = .01 * zseps
      rseps2 = .01 * rseps2
      zseps2 = .01 * zseps2
      rvsin  = .01 * rvsin
      zvsin  = .01 * zvsin
      rvsout = .01 * rvsout
      zvsout = .01 * zvsout

c     NOTE : simagx and sibdry flux values from aeqdsk and neqdsk
c            have opposite signs; the neqdsk values are consistent
c            with the EFIT flux array, fold(xold,yold).

 1040 format (1x,4e16.9)
 1050 format (1x,i5,11x,i5)
 1056 format (1x,a10,i2,1x,i2,1x,i4)
 1060 format (1x,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3)

      return
      end

c     ------------------------------------------------------------------

      subroutine aeqdsk1 (iunit)
      implicit none
      integer iunit, k
      real tsaisq, rcencm, bcentr, pasmat, cpasma, rout, zout
      real aout, eout, doutu, doutl, vout, rcurrt, zcurrt, qsta
      real betat, betap, ali, oleft, oright, otop, obott, qpsib
      real vertn, shearb, bpolav, s1, s2, s3, qout, olefs, orighs
      real otops, sibdry, areao, wplasm, terror, elongm, qqmagx
      real cdflux, alpha, rttt, psiref, xndnt, sepexp, obots
      real btaxp, btaxv, aaq1, aaq2, aaq3, seplim, rmagx, zmagx
      real simagx, taumhd, betapd, betatd, wplasmd, fluxx, vloopt
      real taudia, qmerci, tavem, pbinj, vsurfa, wpdot, wbdot
      real slantu, slantl, zuperts, chipre, cjor95, pp95, xdum
      real xxxxxx, xxx
Use(Aeqflxgrd)   # rseps,zseps,rseps2,zseps2
              # rvsin,zvsin,rvsout,zvsout
              # nsilop,magpri,nfcoil,nesum
      character*42 header
      integer efitvers
      external efitvers

        read (iunit,1040) tsaisq,rcencm,bcentr,pasmat
        read (iunit,1040) cpasma,rout,zout,aout
        read (iunit,1040) eout,doutu,doutl,vout
        read (iunit,1040) rcurrt,zcurrt,qsta,betat
        read (iunit,1040) betap,ali,oleft,oright
        read (iunit,1040) otop,obott,qpsib,vertn
        read (iunit,1040) (rco2v(k),k=1,mco2v)
        read (iunit,1040) (dco2v(k),k=1,mco2v)
        read (iunit,1040) (rco2r(k),k=1,mco2r)
        read (iunit,1040) (dco2r(k),k=1,mco2r)
        read (iunit,1040) shearb,bpolav,s1,s2
        read (iunit,1040) s3,qout,olefs,orighs
        read (iunit,1040) otops,sibdry,areao,wplasm
        read (iunit,1040) terror,elongm,qqmagx,cdflux
        read (iunit,1040) alpha,rttt,psiref,xndnt
        read (iunit,1040) rseps,zseps,rseps2,zseps2
        read (iunit,1040) sepexp,obots,btaxp,btaxv
        read (iunit,1040) aaq1,aaq2,aaq3,seplim
        read (iunit,1040) rmagx,zmagx,simagx,taumhd
        read (iunit,1040) betapd,betatd,wplasmd,fluxx
        read (iunit,1040) vloopt,taudia,qmerci,tavem
        if (efitvers(vmonth,vday,vyear) .ge. 1) 
     .     read (iunit,1044) nsilop,magpri,nfcoil,nesum
        call gchange("Aeqflxgrd",0)
        call aeqdsk2 (iunit)
        read (iunit,1040) pbinj,rvsin,zvsin,rvsout
        read (iunit,1040) zvsout,vsurfa,wpdot,wbdot
        read (iunit,1040) slantu,slantl,zuperts,chipre
        read (iunit,1040) cjor95,pp95,xdum,xdum
        read (iunit,1040) xxxxxx,xxx,xxx,xxx
        read (iunit,1042) header

 1040 format (1x,4e16.9)
 1042 format (1x,a42)
 1044 format (1x,4i5)

      return
      end

c     ------------------------------------------------------------------

      subroutine aeqdsk2 (iunit)
      implicit none
      integer iunit,k
Use(Aeqflxgrd)

        read (iunit,1040) (csilop(k),k=1,nsilop),(cmpr2(k),k=1,magpri)
        read (iunit,1040) (ccbrsp(k),k=1,nfcoil)
        read (iunit,1040) (eccurt(k),k=1,nesum)

 1040 format (1x,4e16.9)

      return
      end

c     ------------------------------------------------------------------

      integer function efitvers (month, day, year)
      implicit none
      integer month, day, year

      efitvers = 0

c     New version of EFIT on 05/24/97 writes aeqdsk that includes
c     data values for parameters nsilop,magpri,nfcoil and nesum.
      if ( ((year .eq. 1997).and.(month .eq. 05).and.(day .ge. 24))
     .      .or. ((year .eq. 1997).and.(month .gt. 05))
     .      .or. (year .gt. 1997) ) efitvers = 1

      return
      end

c     ------------------------------------------------------------------

      subroutine neqdsk
      implicit none
      integer i, idum
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#geqdskfname
Use(Polflx)
      integer iunit, ios
      external gallot, xerrab, remark
      character*8 label(6)

c *************** read the output of the EFIT code ******************
      data iunit /55/
      open (iunit, file=geqdskfname, form='formatted', iostat=ios,
     .      status='old')
      if (ios .ne. 0) then
         call xerrab("**** geqdsk or neqdsk file not found")
      endif

      read(iunit,2000) (label(i),i=1,6),idum,nxefit,nyefit
      runid=label(1)//label(2)//label(3)//label(4)//label(5)//label(6)

c     set length for 2-d spline workspace --
      nwork = nxefit*nyefit + 2*max(kxord*(nxefit+1),kyord*(nyefit+1))
      call gallot("Comflxgrd",0)
      nx4 = mrfac*(nxefit-1) + 1
      ny4 = mrfac*(nyefit-1) + 1
      call gallot("Polflx",0)
      call rdgrp1(iunit)

      read(iunit,2022) nbdry,nlim
      call gchange("Comflxgrd",0)
      call rdgrp2(iunit)

      close (iunit)

 2000 format(6a8,3i4)
 2022 format(2i5)

      return
      end

c     ------------------------------------------------------------------

      subroutine rdgrp1(iunit)
      implicit none
      integer iunit
      real xdum
Use(Dimflxgrd)	#nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
      integer i, j
      read(iunit,2020) xdim,zdim,rcentr,rgrid1,zmid
      read(iunit,2020) rmagx,zmagx,simagx,sibdry,bcentr
      read(iunit,2020) cpasma,simagx,xdum,rmagx,xdum
      read(iunit,2020) zmagx,xdum,sibdry,xdum,xdum
      read(iunit,2020) (fpol(i),i=1,nxefit)
      read(iunit,2020) (pres(i),i=1,nxefit)
      read(iunit,2020) (workk(i),i=1,nxefit)
      read(iunit,2020) (workk(i),i=1,nxefit)
      read(iunit,2020) ((fold(i,j),i=1,nxefit),j=1,nyefit)
      read(iunit,2020) (qpsi(i),i=1,nxefit)
 2020 format(5e16.9)

      return
      end

c     ------------------------------------------------------------------

      subroutine rdgrp2(iunit)
      integer iunit
Use(Dimflxgrd)	#nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
      integer i

      read(iunit,2020) (rbdry(i),zbdry(i),i=1,nbdry)
      read(iunit,2020) (xlim(i),ylim(i),i=1,nlim)
 2020 format(5e16.9)

      return
      end

c     ------------------------------------------------------------------

      subroutine convert_eqdsk
      implicit none
Use(Share)    # geometry
Use(Dimflxgrd)	#nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Aeqflxgrd)
      integer i,j
      real temp,tempr,tempz

c     This subroutine converts eqdsk data for an upper-single-null
c     configuration to make it look like a lower-single-null.
c     Effectively, the observer is rotated by 180 degrees about
c     the midplane of the EFIT domain, so (R,Z)-->(R,-Z) and the
c     signs of the toroidal current and magnetic field are reversed.

      zseps  = 2*zmid - zseps
      zseps2 = 2*zmid - zseps2
      zmagx  = 2*zmid - zmagx
      zvsin  = 2*zmid - zvsin
      zvsout = 2*zmid - zvsout
      do i=1,nlim
         ylim(i) = 2*zmid - ylim(i)
      enddo
      do i=1,nbdry
         zbdry(i) = 2*zmid - zbdry(i)
      enddo

C     Change signs of toroidal current and magnetic field:
      cpasma = - cpasma
      bcentr = - bcentr

c     Change poloidal magnetic flux data:
c        (1) interchange top and bottom
c        (2) change sign of psi
      do i=1,nxefit
         do j=1,nyefit/2
            temp = fold(i,j)
            fold(i,j) = - fold(i,nyefit-j+1)
            fold(i,nyefit-j+1) = - temp
         enddo
c     if nyefit is odd, middle row must be treated separately:
         if (mod(nyefit,2) .eq. 1) then
            fold(i,nyefit/2+1) = - fold(i,nyefit/2+1)
         endif
      enddo
      simagx = - simagx

c     Interchange upper and lower x-point data:
      if ( (geometry=="dnbot") .or. (geometry=="dnull") .or.
     .     (geometry=="isoleg") ) then
c     We're doing the upper half of a double-null configuration
         sibdry = - sibdry2
      else
c     We're doing an upper single-null configuration
         sibdry = - sibdry
      endif

      tempr=rseps
      tempz=zseps
      rseps=rseps2
      zseps=zseps2
      rseps2=tempr
      zseps2=tempz

      return
      end

c     ------------------------------------------------------------------

      subroutine inflx1
      implicit none
Use(Share)      # geometry,islimon,theta_split
Use(Dim)
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Dimflx)	#nsearch
Use(Inpf0)	# istcvon,altsearch,isetpath
Use(Inpf)
Use(Aeqflxgrd)
Use(Polflx)
Use(Flxin)      # psi0lim,sfaclim,kymesh,xcutoff1,icutoff1,
                # ycutoff1,jcutoff1,ymax1fac,ymax2fac
Use(Limiter)    # nlimu,rlimu,zlimu,nptnma,rptnma,zptnma

c     local variables --
      integer i, j, n, ns, nmax, nc2, nc3
      real dxefit, dyefit, alfcy_outer, alfold, alfnew, b
      real eps, dpsidt_inner, dpsidt_outer
      real t1, t2, t3, r1, r2, r3, r2p, dx, dy, pi, twopi

      external remark, xerrab, rho1

      data nmax /30/	# maximum number of newton iterations for alfcy
      data eps /1.0e-10/	# convergence criterion for alfcy

      psi0max_outer = psi0max
      alfcy_outer   = alfcy
# For single-null, ensure that SOL flux contours have same
# values on inboard and outboard halves of mesh:
      if ((geometry .eq. 'snull').or.(geometry .eq. "uppersn")) then
         psi0max_inner = psi0max_outer
         alfcy_inner = alfcy_outer
      endif

# Check psi0max to avoid possible upper x-point for single nulls
      if (((geometry .eq. 'snull').or.(geometry .eq. "uppersn"))
     &          .and. (rseps2 .gt. 0)) then
         psi0sep2 = (sibdry2-simagx)/(sibdry-simagx)
         if (psi0max .gt. psi0sep2) then
	    psi0max = 0.999 * psi0sep2
	    call remark(" ")
	    call remark("***    a second x-point exists    ***")
	    call remark("***  psi0max has been re-defined  ***")
	    call remark(" ")
	    write (STDOUT,801) psi0max
 801        format (5x,"psi0max = ",1p1e14.6)
	    call remark(" ")
            psi0max_inner = psi0max
            psi0max_outer = psi0max
	 endif
      endif

c     Flux contour values for radial mesh:
c----------------------------------------------------------------------c
      if (kymesh==0) then   # use existing psitop and psibot arrays
c----------------------------------------------------------------------c
         continue
c----------------------------------------------------------------------c
      else                  # compute new psitop and psibot arrays
c----------------------------------------------------------------------c
# Define distribution of radial mesh points via analytic function
# rho(t) where t is mesh index and rho is normalized flux.
      do i=0,nym
         tflx(i) = i
      enddo
      t1 = tflx(0)
      t2 = tflx(nycore(igrid))
      t3 = tflx(nycore(igrid)+nysol(igrid))

c     Set some parameters for various mesh forms:
      if (kymesh==1) then # Set SOL exponential factors

# First, to avoid overflow from 1/sinh(alfcy):
        alfcy_outer = max (alfcy_outer, 1.0e-07)
        alfcy_inner = max (alfcy_inner, 1.0e-07)
# When psi0max_inner differs from psi0max_outer we ensure that the core 
# (and  private-flux) distributions on the two halves of the mesh are
# the same by requiring that d(psi)/d(t) at the separatrix has the same
# value on the inboard and outboard halves of the mesh.
        dpsidt_outer = ((psi0max_outer - psi0sep)/(t3-t2))
     &               * alfcy_outer / sinh (alfcy_outer)
        dpsidt_inner = ((psi0max_inner - psi0sep)/(t3-t2))
     &               * alfcy_inner / sinh (alfcy_inner)
        if (psi0max_outer .eq. psi0max_inner) then	# use alfcy_outer
         alfcy_inner = alfcy_outer			# for both halves
         if ( (geometry .eq. 'dnbot').or.(geometry .eq. 'dnull') .or.
     .        (geometry == 'isoleg') ) then
	    call remark(" ")
            call remark("***  alfcy_inner is set to alfcy_outer  ***")
	    call remark(" ")
            write (STDOUT,804) alfcy_inner
 804        format (5x,"alfcy_inner = ",1p1e14.6)
	    call remark(" ")
         endif
        elseif (psi0max_outer .lt. psi0max_inner) then	# re-compute alfcy_inner
	 b = ((psi0max_outer - psi0sep)/(psi0max_inner - psi0sep))
     &       * alfcy_outer / sinh (alfcy_outer)
	 alfold = 1.1 * log ((1+sqrt(1-b**2))/b)	# this initial guess
							# ensures the non-zero
							# root is obtained
	 do n=1,nmax
	    alfnew = alfold - (alfold-b*sinh(alfold))/(1.-b*cosh(alfold))
	    if (abs((alfnew-alfold)/alfcy_outer) .le. eps) goto 555
	    alfold = alfnew
	 enddo
  555    alfcy_inner = alfnew     
	 call remark(" ")
         call remark("***  alfcy_inner has been re-defined  ***")
	 call remark(" ")
         write (STDOUT,802) alfcy_inner
 802     format (5x,"alfcy_inner = ",1p1e14.6)
	 call remark(" ")
        elseif (psi0max_inner .lt. psi0max_outer) then	# re-compute alfcy_outer
	 b = ((psi0max_inner - psi0sep)/(psi0max_outer - psi0sep))
     &       * alfcy_inner / sinh (alfcy_inner)
	 alfold = 1.1 * log ((1+sqrt(1-b**2))/b)	# this initial guess
							# ensures the non-zero
							# root is obtained
	 do n=1,nmax
	    alfnew = alfold - (alfold-b*sinh(alfold))/(1.-b*cosh(alfold))
	    if (abs((alfnew-alfold)/alfcy_inner) .le. eps) goto 666
	    alfold = alfnew
	 enddo
  666    alfcy_outer = alfnew
         alfcy = alfcy_outer     
	 call remark(" ")
         call remark("***  alfcy has been re-defined  ***")
	 call remark(" ")
         write (STDOUT,803) alfcy
 803     format (5x,"alfcy = ",1p1e14.6)
	 call remark(" ")
        endif

      elseif (kymesh==2) then   # Set radial mesh slope factor at separatrix
        r2p = (psi0sep-psi0min1)/float(nycore(igrid))   # rho'(t2)
      endif     # end if-test on kymesh for setting mesh parameters

# For inboard half of mesh:
# In core region + SOL
      t1 = tflx(0)
      t2 = tflx(nycore(igrid))
      t3 = tflx(nycore(igrid)+nysol(igrid))
      r1 = psi0min1
      r2 = psi0sep
      r3 = psi0max_inner
      if (kymesh==1) then    # rational + exponential function rho1(t)
         call rho1(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,alfcy_inner)
      elseif (kymesh==2) then # uniform radial mesh in core 
         call rho5(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,r2p)
      endif
# Core region may be modified for limiter configuration --
      if ((islimon .ne. 0) .and. (sfaclim .ne. 1.) .and.
     .    (psi0lim .le. psi0sep) .and. (psi0lim .ge. psi0min1)) then
         nc2 = int(((psi0lim-psi0min1)/(psi0sep-psi0min1))*nycore(igrid))
         nc3 = nycore(igrid)
         t1 = tflx(0)
         t2 = tflx(nc2)
         t3 = tflx(nc3)
         r1 = psi0min1
         r2 = psi0lim
         r3 = psi0sep
         call rho4(tflx,rho,nc3+1,t1,t2,t3,r1,r2,r3,sfaclim)
      endif
      do i=0,nym
	 psitop(i+1) = simagx+(sibdry-simagx)*rho(nym-i)
      enddo
      psitop(nym+2) = simagx

# In private region + SOL
      t1 = tflx(0)
      t2 = tflx(nycore(igrid))
      t3 = tflx(nycore(igrid)+nysol(igrid))
      r1 = psi0min2
      r2 = psi0sep
      r3 = psi0max_inner
      if (kymesh==1) then    # rational + exponential function rho1(t)
         call rho1(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,alfcy_inner)
      elseif (kymesh==2) then # uniform radial mesh in core 
         call rho5(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,r2p)
      endif
      do i=0,nym
	 psibot(i+1) = simagx+(sibdry-simagx)*rho(nym-i)
      enddo
      psibot(nym+2) = simagx

# For outboard half of mesh:
# In core region + SOL
      t1 = tflx(0)
      t2 = tflx(nycore(igrid))
      t3 = tflx(nycore(igrid)+nysol(igrid))
      r1 = psi0min1
      r2 = psi0sep
      r3 = psi0max_outer
      if (kymesh==1) then    # rational + exponential function rho1(t)
         call rho1(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,alfcy_outer)
      elseif (kymesh==2) then # uniform radial mesh in core 
         call rho5(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,r2p)
      endif
# Core region may be modified for limiter configuration --
      if ((islimon .ne. 0) .and. (sfaclim .ne. 1.) .and.
     .    (psi0lim .le. psi0sep) .and. (psi0lim .ge. psi0min1)) then
         nc2 = int(((psi0lim-psi0min1)/(psi0sep-psi0min1))*nycore(igrid))
         nc3 = nycore(igrid)
         t1 = tflx(0)
         t2 = tflx(nc2)
         t3 = tflx(nc3)
         r1 = psi0min1
         r2 = psi0lim
         r3 = psi0sep
         call rho4(tflx,rho,nc3+1,t1,t2,t3,r1,r2,r3,sfaclim)
      endif
      do i=0,nym
	 psitop(jdim-i) = simagx+(sibdry-simagx)*rho(nym-i)
      enddo

# In private region + SOL
      t1 = tflx(0)
      t2 = tflx(nycore(igrid))
      t3 = tflx(nycore(igrid)+nysol(igrid))
      r1 = psi0min2
      r2 = psi0sep
      r3 = psi0max_outer
      if (kymesh==1) then    # rational + exponential function rho1(t)
         call rho1(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,alfcy_outer)
      elseif (kymesh==2) then # uniform radial mesh in core 
         call rho5(tflx,rho,nym+1,t1,t2,t3,r1,r2,r3,r2p)
      endif
      do i=0,nym
	 psibot(jdim-i) = simagx+(sibdry-simagx)*rho(nym-i)
      enddo
c----------------------------------------------------------------------c
      endif # end if-test on kymesh
c----------------------------------------------------------------------c

# define some search bounds for contouring -
      xlbnd = max(rgrid1,xcutoff1)
      xubnd = rgrid1+xdim
      ylbnd = max(0.,ycutoff1)
      yubnd = zdim
      pi = 4*atan2(1.,1.)
      twopi = 2*pi
      if (islimon .ne. 0) then  # compute theta_split
         istchkon=1
         if (nlimu .lt. 2) then
            call remark("***")
            call remark("*** inflx1: limiter surface r,zlimu not defined")
            call remark("***")
            call xerrab("")
         endif
         # find point on limiter that is nearest the magnetic axis:
         call findptnma(nlimu,rlimu,zlimu,rmagx,zmagx,
     .                                    nptnma,rptnma,zptnma)
         # define poloidal angle to limiter contact point:
         theta_split=atan2 (zptnma-zmagx, rptnma-rmagx)
      endif
      if (istchkon .gt. 0) then
         # do not limit radial search in regions 1 and 2:
         xoverlap(1)=1000
         xoverlap(2)=1000
         # compute angular bounds:
         thetax = atan2 (zseps-zmagx, rseps-rmagx)
         if ( (geometry=="dnbot" .or. geometry=="dnull" .or.
     .         geometry=="isoleg") .and. isthmmxn==1) then
            thetamin(1)=pi-dtheta_overlap_sol(1)
            thetamax(2)=twopi+dtheta_overlap_sol(2)
         else
            thetamin(1)=theta_split-dtheta_overlap_sol(1)
            thetamax(2)=theta_split+dtheta_overlap_sol(2)
         endif
         thetamax(1)=thetamin(1)+twopi-dtheta_exclude(1)
         thetamin(2)=thetamax(2)-twopi+dtheta_exclude(2)
      endif

# minimum and maximum cell indices on the refined EFIT grid --
      imin = 1
      imax = mrfac*(nxefit-1)+1
      jmins = 1
      jmaxs = mrfac*(nyefit-1)+1

# define indices of cells containing the magnetic axis and the x-point :
      dxefit = xdim/(nxefit-1)
      dyefit = zdim/(nyefit-1)
      dx = dxefit/float(mrfac)
      dy = dyefit/float(mrfac)
      imagx = (rmagx-rgrid1)/dx+1
      jmagx = (zmagx)/dy+1
      iseps = (rseps-rgrid1)/dx+1
      jseps = (zseps)/dy+1
      icutoff1 = max(0.,xcutoff1-rgrid1)/dx+1
      jcutoff1 = max(0.,ycutoff1)/dy+1


      ns=0    # temporary counter for number of search regions
      ns=ns+1
# Search region 1 (inboard half of scrape-off layer and core plasma) :
      ncmin0(ns) = jmin(1)
      ncmax0(ns) = jmax(1)
      do n=1,ncmax0(ns)-ncmin0(ns)+1
         plflux0(n,ns) = psitop(ncmin0(ns)+n-1)
      enddo
      xminf0(ns) = xlbnd
      xmaxf0(ns) = min( xubnd, max(rmagx,rseps)+xoverlap(1)*dxefit )
      yminf0(ns) = ylbnd
      if ((geometry .eq. 'snull').or.(geometry .eq. "uppersn")) then
         ymaxf0(ns) = yubnd
      else
         ymaxf0(ns) = ymax1fac*zmagx + 5*dyefit
      endif
# Set initial search position and direction:
      if (isetpath==0) then   # search radially outward from left boundary
      iserch0(ns) = icutoff1
      jserch0(ns) = jmagx
      istepf0(ns) = 1
      jstepf0(ns) = 0
      else                    # search radially inward from magnetic axis
      iserch0(ns) = imagx
      jserch0(ns) = jmagx
      istepf0(ns) = -1
      jstepf0(ns) = 0
      endif

      if (nycore(igrid) .gt. 0) then
         ns=ns+1
# Search region 2 (inboard half of divertor triangle) :
	 ncmin0(ns) = jsptrx(1) + 1
	 ncmax0(ns) = jmax(1)
         do n=1,ncmax0(ns)-ncmin0(ns)+1
            plflux0(n,ns) = psibot(ncmin0(ns)+n-1)
         enddo
         xminf0(ns) = xlbnd
         xmaxf0(ns) = min( xubnd, rseps+xoverlap(1)*dxefit )
	 yminf0(ns)=ylbnd
	 ymaxf0(ns) = ymax2fac*zseps
cccMER (28 Mar 2002) fix for Ignitor config'n with inner strike pt above xpt
ccc	 ymaxf0(ns)=zseps
cccTDR (24 Feb 2008) fix for snowflake
ccc	 ymaxf0(ns)=zseps*1.1
# Set initial search position and direction:
         if (istcvon==1) then
            write(*,*) "**WARNING: istcvon=1 OBSOLETE; setting altsearch=2"
            altsearch = 2
         endif
         if (altsearch==0) then       # search upward toward x-point
            iserch0(ns) = iseps
            jserch0(ns) = jcutoff1
            istepf0(ns) = 0
            jstepf0(ns) = 1
         elseif (altsearch==1) then   # search downward from x-point
            iserch0(ns) = iseps
            jserch0(ns) = jseps
            istepf0(ns) = 0
            jstepf0(ns) = -1
         elseif (altsearch==2) then   # search diagonally from x-point
            iserch0(ns) = iseps
            jserch0(ns) = jseps
	    istepf0(ns) = -2          # MER (28 Aug 2012) fix for NSTX
            jstepf0(ns) = -1
         elseif (altsearch==3) then   #	MVU (7-Dec-2023) fix for NT configuration
            iserch0(ns) = iseps
            jserch0(ns) = jseps
	    istepf0(ns) =  2
            jstepf0(ns) = -1
         else
	    print *, "Invalid altsearch value ", altsearch
         endif
      endif  # end if-test on nycore for region 2

      ns=ns+1
# Search region 3 (outboard half of core plasma and scrape-off layer) :
      ncmin0(ns) = jmin(2)
      ncmax0(ns) = jmax(2)
      do n=1,ncmax0(ns)-ncmin0(ns)+1
         plflux0(n,ns) = psitop(ncmin0(ns)+n-1)
      enddo
      xmaxf0(ns) = xubnd
      xminf0(ns) = max( xlbnd, min(rmagx,rseps)-xoverlap(2)*dxefit )
      yminf0(ns) = ylbnd
      if ((geometry .eq. 'snull').or.(geometry .eq. "uppersn")) then
         ymaxf0(ns) = yubnd
      else
         ymaxf0(ns) = zmagx + 5*dyefit
      endif
# Set initial search position and direction:
      iserch0(ns) = imagx+1
      jserch0(ns) = jmagx
      istepf0(ns) = 1
      jstepf0(ns) = 0

      if (nycore(igrid) .gt. 0) then
         ns=ns+1
# Search region 4 (outboard half of divertor triangle) 
	 ncmin0(ns) = jmin(2)
	 ncmax0(ns) = jsptrx(2) - 1
         do n=1,ncmax0(ns)-ncmin0(ns)+1
            plflux0(n,ns) = psibot(ncmin0(ns)+n-1)
         enddo
         xminf0(ns) = max( xlbnd, rseps-xoverlap(2)*dxefit )
         xmaxf0(ns) = xubnd
	 yminf0(ns) = ylbnd
	 ymaxf0(ns) = ymax4fac*zseps #-MVU 19-dec-2023, adding ymax4fac	/1.0/
# Set initial search position and direction:
         if (istcvon==1) then
            write(*,*) "**WARNING: istcvon=1 OBSOLETE; setting altsearch=2"
            altsearch = 2
         endif
         if (altsearch==0) then       # search upward toward x-point
            iserch0(ns) = iseps
            jserch0(ns) = jcutoff1
            istepf0(ns) = 0
            jstepf0(ns) = 1
         elseif (altsearch==1) then   # search downward from x-point
            iserch0(ns) = iseps
            jserch0(ns) = jseps
            istepf0(ns) = 0
            jstepf0(ns) = -1
         elseif (altsearch==2) then   # search diagonally from x-point
            iserch0(ns) = iseps
            jserch0(ns) = jseps
            istepf0(ns) = -2          # MER (28 Aug 2012) fix for NSTX
            jstepf0(ns) = -1
	 elseif (altsearch==3) then   # MVU (7-Dec-2023) fix for NT configuration
            iserch0(ns) = iseps
            jserch0(ns) = jseps
            istepf0(ns) =  2
            jstepf0(ns) = -1
         else
            print *, "Invalid altsearch value ", altsearch
         endif
      endif  # end if-test on nycore for region 4

      if (ns .ne. nsearch) then
	 call remark(" ")
	 call remark("*** inflx1 : search region indexing error")
	 call remark(" ")
	 call xerrab("")
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine rho1(t,rho,nt,t1,t2,t3,r1,r2,r3,alf)
      implicit none

c     input/output arguments --
      integer nt
      real rho(nt),r1,r2,r3,alf
      real t(nt),t1,t2,t3

c     local variables --
      integer n
      real epslon
      real a,b
      real r2p

c     procedures --
      external xerrab

c     initialization --
      data epslon /1.0e-7/

c     rho(t) is a rational function of the form (f+gt)/(p+qt) on
c     the interval t1 < t < t2 and has an exponential form on the 
c     interval t2 < t < t3 with continuous first derivative at t2.

      if (t3 .le. t2) then
         write (STDOUT,900)
         call xerrab("")
      endif
  900 format("error in subroutine rho1; bad input parameters: t3 .le. t2")

      alf = max (epslon,alf)

      r2p = ((r3-r2)/(t3-t2)) * (2*alf)/(exp(alf)-exp(-alf))
      a = r2p * (t2-t1)**2 / (r2-r1)**2
      b = ((r2-r1) - r2p*(t2-t1)) / (r2-r1)**2

      do n = 1, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r1 + (t(n)-t1) / (a + b*(t(n)-t1))
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + r2p * ((t3-t2)/(2*alf)) *
     &           ( exp(+alf*(t(n)-t2)/(t3-t2)) -
     &             exp(-alf*(t(n)-t2)/(t3-t2)) )
         elseif (t3 .le. t(n)) then
            rho(n) = r3
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho1dn(t,rho,nt,t1,t2,t3,t4,r1,r2,r3,r4,alf)
      implicit none

c     input/output arguments --
      integer nt
      real rho(0:nt),r1,r2,r3,r4,alf
      real t(0:nt),t1,t2,t3,t4

c     local variables --
      integer n
      real epslon
      real a,b,c,d
      real r2p,r3p

c     initialization --
      data epslon /1.0e-10/

c     rho(t) is a rational function of the form (f+gt)/(p+qt) on
c     the interval t1 < t < t2 and also on the interval t2 < t < t3,
c     and has an exponential form on the interval t3 < t < t4
c     with rho(t) and rho'(t) continuous at t2 and t3.

      if (t4 .le. t3) then
         call remark("error in subroutine rho1dn; bad input parameters:")
         write (STDOUT,901) t4, t3
         call xerrab("")
      endif
  901 format("t4=",f4.1," .le. t3=",f4.1)
      if (r2 .le. r1) then
         call remark("error in subroutine rho1dn; bad input parameters:")
         write (STDOUT,902) r2, r1
         call xerrab("")
      endif
  902 format("r2=",f9.6," .le. r1=",f9.6)
      if (r3 .le. r2) then
         call remark("error in subroutine rho1dn; bad input parameters:")
         write (STDOUT,903) r3, r2
         call xerrab("")
      endif
  903 format("r3=",f9.6," .le. r2=",f9.6)

      alf = max (epslon,alf)

c     Assume rho(t) has the form:
c		rho(t)=rho(t1)+(t-t1)/(a+b*(t-t1))
c     on the interval t1 < t < t2 (in the core), and
c		rho(t)=rho(t2)+(t-t2)/(c+d*(t-t2))
c     on the interval t2 < t < t3 (between separatrices), and
c		rho(t)=rho(t3)+k*sinh(alf*(t-t3)/(t4-t3))
c     on the interval t3 < t < t4 (in outer SOL)  where
c		k=(rho(t4)-rho(t3))/sinh(alf)
c
c     Choose coefficients (a,b) and (c,d) so that rho(t) and rho'(t)
c     are continuous at t2 and t3 :
      r3p = ((r4-r3)/(t4-t3)) * alf/sinh(alf)	# rho'(t3)
      d = ((r3-r2) - r3p*(t3-t2)) / (r3-r2)**2
      c = (1-d*(r3-r2))*(t3-t2)/(r3-r2)
      r2p=1/c					# rho'(t2)
      b = ((r2-r1) - r2p*(t2-t1)) / (r2-r1)**2
      a = (1-b*(r2-r1))*(t2-t1)/(r2-r1)

      do n = 0, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r1 + (t(n)-t1) / (a + b*(t(n)-t1))
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + (t(n)-t2) / (c + d*(t(n)-t2))
         elseif (t3 .le. t(n) .and. t(n) .lt. t4) then
            rho(n) = r3 + r3p*(t4-t3)*sinh(alf*(t(n)-t3)/(t4-t3))/alf
         elseif (t4 .le. t(n)) then
            rho(n) = r4
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho2dn(t,rho,nt,t1,t2,t3,t4,r1,r2,r3,r4,fac)
      implicit none

c     input/output arguments --
      integer nt
      real rho(0:nt),r1,r2,r3,r4,fac
      real t(0:nt),t1,t2,t3,t4

c     local variables --
      integer n
      real a,b,c,d,e,f
      real r2p,r3p

c     rho(t) is a piece-wise rational function of the form (f+gt)/(p+qt)
c     on the intervals t1 < t < t2, t2 < t < t3 and t3 < t < t4 with
c     rho(t) and rho'(t) continuous at t2 and t3.
c     The slopes rho'(t2) and rho'(t3) scale with input argument fac
c     according to the relations:
c                     rho'(t2) = ((r3-r2)/(t3-t2))*fac
c     and
c                     rho'(t3) = ((r3-r2)/(t3-t2))/fac
c     
c     NOTE: For fac=1, rho(t) is linear on the interval t2 < t < t3.

      if (t2 .le. t1) then
         call remark("error in subroutine rho2dn; bad input parameters:")
         write (STDOUT,899) t2, t1
         call xerrab("")
      endif
  899 format("t2=",f4.1," .le. t1=",f4.1)
      if (t3 .le. t2) then
         call remark("error in subroutine rho2dn; bad input parameters:")
         write (STDOUT,900) t3, t2
         call xerrab("")
      endif
  900 format("t3=",f4.1," .le. t2=",f4.1)
      if (t4 .le. t3) then
         call remark("error in subroutine rho2dn; bad input parameters:")
         write (STDOUT,901) t4, t3
         call xerrab("")
      endif
  901 format("t4=",f4.1," .le. t3=",f4.1)
      if (r2 .le. r1) then
         call remark("error in subroutine rho2dn; bad input parameters:")
         write (STDOUT,902) r2, r1
         call xerrab("")
      endif
  902 format("r2=",f9.6," .le. r1=",f9.6)
      if (r3 .le. r2) then
         call remark("error in subroutine rho2dn; bad input parameters:")
         write (STDOUT,903) r3, r2
         call xerrab("")
      endif
  903 format("r3=",f9.6," .le. r2=",f9.6)
      if (r4 .le. r3) then
         call remark("error in subroutine rho2dn; bad input parameters:")
         write (STDOUT,904) r4, r3
         call xerrab("")
      endif
  904 format("r4=",f9.6," .le. r3=",f9.6)

c     Assume rho(t) has the form:
c		rho(t)=rho(t1)+(t-t1)/(a+b*(t-t1))
c     on the interval t1 < t < t2 (in the core), and
c		rho(t)=rho(t2)+(t-t2)/(c+d*(t-t2))
c     on the interval t2 < t < t3 (between separatrices), and
c		rho(t)=rho(t3)+(t-t3)/(e+f*(t-t3))
c     on the interval t3 < t < t4 (beyond outer separatrix).
c
c     Choose coefficients (a,b),(c,d) and (e,f) so that rho(t) and rho'(t)
c     are continuous at t2 and t3 :
      r2p = ((r3-r2)/(t3-t2))*fac               # rho'(t2)
      r3p = ((r3-r2)/(t3-t2))/fac               # rho'(t3)

      b = ((r2-r1) - r2p*(t2-t1)) / (r2-r1)**2
      a = (1-b*(r2-r1))*(t2-t1)/(r2-r1)

      c = 1/r2p
      d = (1-c*(r3-r2)/(t3-t2))/(r3-r2)

      e = 1/r3p
      f = (1-e*(r4-r3)/(t4-t3))/(r4-r3)

      do n = 0, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r1 + (t(n)-t1) / (a + b*(t(n)-t1))
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + (t(n)-t2) / (c + d*(t(n)-t2))
         elseif (t3 .le. t(n) .and. t(n) .lt. t4) then
            rho(n) = r3 + (t(n)-t3) / (e + f*(t(n)-t3))
         elseif (t4 .le. t(n)) then
            rho(n) = r4
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho3dn(t,rho,nt,t1,t2,t3,t4,r1,r2,r3,r4,slp2fac,slp3fac,
     .                         r2p,r3p)
      implicit none

c     input/output arguments --
      integer nt
      real rho(0:nt),r1,r2,r3,r4,slp2fac,slp3fac,r2p,r3p
      real t(0:nt),t1,t2,t3,t4

c     local variables --
      integer n
      real a,b,c,d,e,f,g

c     rho(t) is a cubic function on the interval t2 < t < t3 and a
c     piece-wise rational function of the form (f+gt)/(p+qt)
c     on the intervals t1 < t < t2 and t3 < t < t4 with
c     rho(t) and rho'(t) continuous at t2 and t3.
c     The slopes rho'(t2) and rho'(t3) scale with input arguments
c     slp2fac and slp3fac according to the relations:
c                     rho'(t2) = ((r3-r2)/(t3-t2))*slp2fac
c     and
c                     rho'(t3) = ((r3-r2)/(t3-t2))*slp3fac
c     
c     so for slp2fac=slp3fac=1, rho(t) is strictly linear on t2 < t < t3.

      if (t2 .le. t1) then
         call remark("error in subroutine rho3dn; bad input parameters:")
         write (STDOUT,899) t2, t1
         call xerrab("")
      endif
  899 format("t2=",f4.1," .le. t1=",f4.1)
      if (t3 .le. t2) then
         call remark("error in subroutine rho3dn; bad input parameters:")
         write (STDOUT,900) t3, t2
         call xerrab("")
      endif
  900 format("t3=",f4.1," .le. t2=",f4.1)
      if (t4 .le. t3) then
         call remark("error in subroutine rho3dn; bad input parameters:")
         write (STDOUT,901) t4, t3
         call xerrab("")
      endif
  901 format("t4=",f4.1," .le. t3=",f4.1)
      if (r2 .le. r1) then
         call remark("error in subroutine rho3dn; bad input parameters:")
         write (STDOUT,902) r2, r1
         call xerrab("")
      endif
  902 format("r2=",f9.6," .le. r1=",f9.6)
      if (r3 .le. r2) then
         call remark("error in subroutine rho3dn; bad input parameters:")
         write (STDOUT,903) r3, r2
         call xerrab("")
      endif
  903 format("r3=",f9.6," .le. r2=",f9.6)
      if (r4 .le. r3) then
         call remark("error in subroutine rho3dn; bad input parameters:")
         write (STDOUT,904) r4, r3
         call xerrab("")
      endif
  904 format("r4=",f9.6," .le. r3=",f9.6)

c     Assume rho(t) has the form:
c		rho(t)=rho(t1)+(t-t1)/(a+b*(t-t1))
c     on the interval t1 < t < t2 (in the core), and
c		rho(t)=rho(t2)+c*(t-t2)+d*(t-t2)**2+g*(t-t2)**3
c     on the interval t2 < t < t3 (between separatrices), and
c		rho(t)=rho(t3)+(t-t3)/(e+f*(t-t3))
c     on the interval t3 < t < t4 (beyond outer separatrix).
c
c     Choose coefficients (a,b),(c,d,g) and (e,f) so that rho(t) and rho'(t)
c     are continuous at t2 and t3 with slopes given by:
      r2p = ((r3-r2)/(t3-t2))*slp2fac               # rho'(t2)
      r3p = ((r3-r2)/(t3-t2))*slp3fac               # rho'(t3)

      b = ((r2-r1) - r2p*(t2-t1)) / (r2-r1)**2
      a = (1-b*(r2-r1))*(t2-t1)/(r2-r1)

      c = r2p
      d = (-2*c - r3p + 3*(r3-r2)/(t3-t2))/(t3-t2)
      g = (c + r3p -2*(r3-r2)/(t3-t2))/(t3-t2)**2

      e = 1/r3p
      f = (1-e*(r4-r3)/(t4-t3))/(r4-r3)

      do n = 0, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r1 + (t(n)-t1) / (a + b*(t(n)-t1))
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + c*(t(n)-t2) + d*(t(n)-t2)**2 + g*(t(n)-t2)**3
         elseif (t3 .le. t(n) .and. t(n) .lt. t4) then
            rho(n) = r3 + (t(n)-t3) / (e + f*(t(n)-t3))
         elseif (t4 .le. t(n)) then
            rho(n) = r4
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho1l(t,rho,nt,t1,t2,r1,r2,r1p)
      implicit none

c     input/output arguments --
      integer nt
      real rho(0:nt),r1,r2,r1p
      real t(0:nt),t1,t2

c     local variables --
      integer n
      real a,b

c     procedures --
      external xerrab

c     rho1l(t) is a rational function of the form (f+gt)/(p+qt) on
c     the interval t1 < t < t2 with fixed slope at t1

      if (t2 .le. t1) then
         write (STDOUT,900)
         call xerrab("")
      endif
  900 format("error in subroutine rho1l; bad input parameters: t2 .le. t1")

      a = 1/r1p
      b = -a/(t2-t1) + 1./(r2-r1)

      do n = 0, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r1 + (t(n)-t1) / (a + b*(t(n)-t1))
         elseif (t2 .le. t(n)) then
            rho(n) = r2
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho1r(t,rho,nt,t1,t2,r1,r2,r2p)
      implicit none

c     input/output arguments --
      integer nt
      real rho(0:nt),r1,r2,r2p
      real t(0:nt),t1,t2

c     local variables --
      integer n
      real a,b

c     procedures --
      external xerrab

c     rho1r(t) is a rational function of the form (f+gt)/(p+qt) on
c     the interval t1 < t < t2 with fixed slope at t2

      if (t2 .le. t1) then
         write (STDOUT,900)
         call xerrab("")
      endif
  900 format("error in subroutine rho1r; bad input parameters: t2 .le. t1")

      a = 1/r2p
      b = -a/(t1-t2) + 1./(r1-r2)

      do n = 0, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r2 + (t(n)-t2) / (a + b*(t(n)-t2))
         elseif (t2 .le. t(n)) then
            rho(n) = r2
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho2(t,rho,nt,t1,t2,t3,r1,r2,r3)
      implicit none

c     input/output arguments --
      integer nt
      real rho(nt),r1,r2,r3
      real t(nt),t1,t2,t3

c     local variables --
      integer n
      real a,b,c

c     procedures --
      external xerrab

c     rho(t) is quadratic on the interval t1 < t < t2
c     and linear on the interval t2 < t < t3 with continuous
c     first derivative at t2.

      if (t3 .le. t2) then
         write (STDOUT,900)
         call xerrab("")
      endif
  900 format("error in subroutine rho2; bad input parameters: t3 .le. t2")

      do n = 1, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            a = r1
            b = 2*(r2-r1)/(t2-t1) - (r3-r2)/(t3-t2)
            c = ((r3-r2)/(t3-t2) - (r2-r1)/(t2-t1))/(t2-t1)
            rho(n) = a + b*(t(n)-t1) + c*(t(n)-t1)**2
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + (r3-r2)*(t(n)-t2)/(t3-t2)
         elseif (t3 .le. t(n)) then
            rho(n) = r3
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho3(t,rho,nt,t1,t2,t3,r1,r2,r3)
      implicit none

c     input/output arguments --
      integer nt
      real rho(nt),r1,r2,r3
      real t(nt),t1,t2,t3

c     local variables --
      integer n
      real a,b,c,d

c     procedures --
      external xerrab

c     rho(t) is cubic on the interval t1 < t < t2 and linear on
c     the interval t2 < t < t3 with continuous first and second
c     derivatives at t2.

      if (t3 .le. t2) then
         write (STDOUT,900)
         call xerrab("")
      endif
  900 format("error in subroutine rho3; bad input parameters: t3 .le. t2")

      do n = 1, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            a = r1
            b = 3*(r2-r1)/(t2-t1) - 2*(r3-r2)/(t3-t2)
            c = -3*((r2-r1)/(t2-t1) - (r3-r2)/(t3-t2))/(t2-t1)
            d = ((r2-r1)/(t2-t1) - (r3-r2)/(t3-t2))/(t2-t1)**2
            rho(n) = a + b*(t(n)-t1) + c*(t(n)-t1)**2 + d*(t(n)-t1)**3
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + (r3-r2)*(t(n)-t2)/(t3-t2)
         elseif (t3 .le. t(n)) then
            rho(n) = r3
         endif
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine rho4(t,rho,nt,t1,t2,t3,r1,r2,r3,s2)
      implicit none

c     input/output arguments --
      integer nt
      real rho(nt),r1,r2,r3,s2
      real t(nt),t1,t2,t3

c     local variables --
      integer n
      real r2p_quadratic,d

c     procedures --
      external xerrab

c     rho(t) is cubic on the interval t1 < t < t3 with the slope
c     at t2 enhanced by a factor s2 over that obtained for a quadratic
c     form.

      if (t3 .le. t2) then
         write (STDOUT,901)
  901 format("error in subroutine rho4; bad input parameters: t3 .le. t2")
         call xerrab("")
      elseif (t2 .le. t1) then
         write (STDOUT,902)
  902 format("error in subroutine rho4; bad input parameters: t2 .le. t1")
         call xerrab("")
      endif

c     For a quadratic form, the slope at t2 would be
      r2p_quadratic = r1*(t2-t3)/((t1-t3)*(t1-t2))
     .              + r2*(2*t2-t1-t3)/((t2-t1)*(t2-t3))
     .              + r3*(t2-t1)/((t3-t1)*(t3-t2))

c     A cubic term, + d*(t-t1)*(t-t2)*(t-t3), adds to the slope at t=t2:
c     r2p = r2p_quadratic + d*(t2-t1)*(t2-t3), so in order to get
c     r2p = s2*r2p_quadratic we should choose

      d = (s2-1.)*r2p_quadratic/((t2-t1)*(t2-t3))

      do n = 1, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r1 *((t(n)-t2)*(t(n)-t3))/((t1-t2)*(t1-t3))
     .             + r2 *((t(n)-t1)*(t(n)-t3))/((t2-t1)*(t2-t3))
     .             + r3 *((t(n)-t1)*(t(n)-t2))/((t3-t1)*(t3-t2))
     .             + d*(t(n)-t1)*(t(n)-t2)*(t(n)-t3)
         elseif (t3 .le. t(n)) then
            rho(n) = r3
         endif
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine rho5(t,rho,nt,t1,t2,t3,r1,r2,r3,r2p)
      implicit none

c     input/output arguments --
      integer nt
      real rho(0:nt),r1,r2,r3,r2p
      real t(0:nt),t1,t2,t3

c     local variables --
      integer n
      real a,b,c,d

c     procedures --
      external xerrab

c     This defines rho(t) as a piece-wise rational function of the form
c          rho = r1 + (t-t1)/(a + b*(t-t1))
c     on the interval t1 < t < t2 and
c          rho = r2 + (t-t2)/(c + d*(t-t2))
c     on the interval t2 < t < t3 
c     with a specified derivative rho'(t2) = r2p

      if (t2 .le. t1) then
         write (STDOUT,900)
         call xerrab("")
      endif
  900 format("error in subroutine rho5; bad input parameters: t2 .le. t1")
      if (t3 .le. t2) then
         write (STDOUT,901)
         call xerrab("")
      endif
  901 format("error in subroutine rho5; bad input parameters: t3 .le. t2")

      a = r2p*((t2-t1)/(r2-r1))**2
      b = (-a+(t2-t1)/(r2-r1))/(t2-t1)
      c = 1/r2p
      d = (-c+(t3-t2)/(r3-r2))/(t3-t2)

      do n = 0, nt
         if (t(n) .lt. t1) then
            rho(n) = r1
         elseif (t1 .le. t(n) .and. t(n) .lt. t2) then
            rho(n) = r1 + (t(n)-t1) / (a + b*(t(n)-t1))
         elseif (t2 .le. t(n) .and. t(n) .lt. t3) then
            rho(n) = r2 + (t(n)-t2) / (c + d*(t(n)-t2))
         elseif (t3 .le. t(n)) then
            rho(n) = r3
         endif
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine findptnma(nlim,rlim,zlim,rmagx,zmagx,
     .                                    nptnma,rptnma,zptnma)
      implicit none
      integer nlim,nptnma
      real rlim(1:nlim),zlim(1:nlim),rmagx,zmagx,rptnma,zptnma

# find the coordinates of the point on the limiter surface that is
# nearest the magnetic axis
      integer n
      real ds, dsmin

      dsmin=sqrt((rlim(1)-rmagx)**2 +(zlim(1)-zmagx)**2)
      nptnma=1
      do n=2,nlim
         ds=sqrt((rlim(n)-rmagx)**2 +(zlim(n)-zmagx)**2)
         if (ds .lt. dsmin) then
            dsmin=ds
            nptnma=n
         endif
      enddo

      rptnma=rlim(nptnma)
      zptnma=zlim(nptnma)

      return
      end

c----------------------------------------------------------------------c

      subroutine findstrike(js,rs,zs)
      implicit none

c     input/output arguments --
      integer js	# flux surface index (input)
      real rs,zs	# strike point coordinates (output)

c     common blocks --
Use(Dimflxgrd)	#jdim,npts,nxefit,nyefit,nlim,nwork
Use(Comflxgrd)	#
Use(Polflx)

c     procedures --
      external intersect2

c     local variables --
      integer icurv,ilim,ierr
      real fuzzf

      data fuzzf /1.0e-06/
 
c     This procedure finds the first intersection of a poloidal
c     flux surface with the limiter surface from EFIT.  The flux 
c     surface index label is js; the flux surface coordinates are
c     xcurve(,js) and ycurve(,js); the limiter surface coordinates
c     are xlim and ylim.

      ierr=0
      call intersect2(xcurve(1,js),ycurve(1,js),1,npoint(js),
     .                xlim,ylim,1,nlim,rs,zs,icurv,ilim,fuzzf,ierr)

c     Also put output into rs_com and zs_com for python version
      rs_com = rs
      zs_com = zs

      if (ierr .ne. 0) then
         call remark("*************************")
         call remark("no intersection was found")
         call remark("*************************")
      endif

      return
      end

