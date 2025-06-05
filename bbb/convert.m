c-----------------------------------------------------------------------
      subroutine convert

*     CONVERT changes from fluid variables inteo variables by solvers.
*     It also calculates the tolerances and the indexes for the cuts.

      implicit none

      Use(Dim)         # nx,ny,nisp,nusp,ngsp
      Use(Math_problem_size)   # neqmx
      Use(UEpar)       # cniatol,cngatol,cupatol,cteatol,ctiatol,cphiatol,
                       # tolbf,isnionxy,isuponxy,isteonxy,istionxy,isngonxy,
                       # isphionxy
      Use(Lsode)       # rtolv,rtol,atol,yl
      Use(Compla)
      Use(Indexes)     # igyl
      Use(Selec)       # ixm1,ixp1
      Use(Ynorm)       # isflxvar,nnorm,ennorm,fnorm,temp0,n0,n0g
      Use(Share)       # 
      Use(Coefeq)      # cngtgx,cngtgy
      Use(Constraints) # icnstr
      Use(Indices_domain_dcl)   #ixmxbcl,ixmnbcl,iymxbcl,iymnbcl
      Use(Phyvar)      # ev

c...  local variables
      real bfac, ntemp
      integer ifld,isfirstvar
      #Former Aux module variables
      integer ix,iy,igsp,iv,ix1,ix2

      isfirstvar = 0   #flag signals first var at ix,iy found
      iv = 0
      do 8 iy = 1-iymnbcl, ny+iymxbcl
         do 6 ix = 1-ixmnbcl, nx+ixmxbcl
            bfac = 1.
            if(ix.eq.0 .or. iy.eq.0 .or. iy.eq.ny+1) bfac = tolbf
            ix1 = ixp1(ix,iy)
            isfirstvar = 0   #flag signals first var at ix,iy found
            do ifld = 1, nisp
	       if(isnionxy(ix,iy,ifld) .eq. 1) then
                  iv = iv + 1
                  yl(iv) = ni(ix,iy,ifld)/n0(ifld)
                  rtol(iv) = rtolv*bfac
                  atol(iv) = cniatol*rtol(iv)*bfac*abs(yl(iv))
                  idxn(ix,iy,ifld) = iv
                  if (isfirstvar==0) ivfirst(ix,iy) = iv
                  isfirstvar = 1
                  igyl(iv,1) = ix
                  igyl(iv,2) = iy
                  icnstr(iv) = 1
               endif
            enddo
               do ifld = 1, nusp
		  if(isuponxy(ix,iy,ifld) .eq. 1) then
                  iv = iv + 1
                  ntemp =0.5*( nm(ix,iy,ifld) + nm(ix1,iy,ifld) )
                  if(isflxvar.eq.0 .or. isflxvar.eq.2) ntemp = mi(ifld)*
     .                                                           n0(ifld)
                  yl(iv) = up(ix,iy,ifld)*ntemp/fnorm(ifld)
                  rtol(iv) = rtolv*bfac
                  atol(iv) = cupatol*rtol(iv)*bfac*sqrt(te(ix,iy)/
     .                                     mi(ifld))*ntemp/fnorm(ifld)
                  idxu(ix,iy,ifld) = iv 
                  if (isfirstvar==0) ivfirst(ix,iy) = iv
                  isfirstvar = 1
                  igyl(iv,1) = ix
                  igyl(iv,2) = iy
                  icnstr(iv) = 0
                  endif
               enddo
	 if(isteonxy(ix,iy) .eq. 1) then
            iv = iv + 1
            ntemp = ne(ix,iy)
            if(isflxvar .eq. 0) ntemp = nnorm
            yl(iv) = 1.5*ntemp*te(ix,iy)/ennorm
            rtol(iv) = rtolv*bfac
            atol(iv) = cteatol*rtol(iv)*bfac*abs(yl(iv))
            idxte(ix,iy) = iv
            if (isfirstvar==0) ivfirst(ix,iy) = iv
            isfirstvar = 1
            igyl(iv,1) = ix
            igyl(iv,2) = iy
            icnstr(iv) = 1
         endif
	 if(istionxy(ix,iy) .eq. 1) then
            iv = iv + 1
            ntemp = nit(ix,iy) + cngtgx(1)*ng(ix,iy,1)
            if(isflxvar .eq. 0) ntemp = nnorm
            yl(iv) = 1.5*ntemp*ti(ix,iy)/ennorm
            rtol(iv) = rtolv*bfac
            atol(iv) = ctiatol*rtol(iv)*bfac*abs(yl(iv))
            idxti(ix,iy) = iv
            if (isfirstvar==0) ivfirst(ix,iy) = iv
            isfirstvar = 1
            igyl(iv,1) = ix
            igyl(iv,2) = iy
            icnstr(iv) = 1
c...  Omit constraint check on x-boundaries for Ti - ckinfl problem
            if (ix.eq.0 .or. ix.eq.nx+1) icnstr(iv) = 0
         endif
         do igsp = 1, ngsp
	  if(isngonxy(ix,iy,igsp) .eq. 1) then
            iv = iv + 1
            yl(iv) = ng(ix,iy,igsp)/n0g(igsp)
            rtol(iv) = rtolv*bfac
            atol(iv) = cngatol*rtol(iv)*bfac*abs(yl(iv))
            idxg(ix,iy,igsp) = iv
            if (isfirstvar==0) ivfirst(ix,iy) = iv
            isfirstvar = 1
            igyl(iv,1) = ix
            igyl(iv,2) = iy
            icnstr(iv) = 1
          endif
         enddo
         do igsp = 1, ngsp
	  if(istgonxy(ix,iy,igsp) .eq. 1) then
            istgcon(igsp) = -1. 
            iv = iv + 1
            ntemp = ng(ix,iy,igsp)
            if(isflxvar == 0) ntemp=n0g(igsp)
	    yl(iv) = 1.5*ntemp*tg(ix,iy,igsp)/ennorm
            rtol(iv) = rtolv*bfac
            atol(iv) = cngatol*rtol(iv)*bfac*abs(yl(iv))
            idxtg(ix,iy,igsp) = iv
            if (isfirstvar==0) ivfirst(ix,iy) = iv
            isfirstvar = 1
            igyl(iv,1) = ix
            igyl(iv,2) = iy
            icnstr(iv) = 1
          endif
         enddo
	 if(isphionxy(ix,iy) .eq. 1) then
            iv = iv + 1
            yl(iv) = phi(ix,iy)/temp0
            rtol(iv) = rtolv*bfac
            atol(iv) = cphiatol*rtol(iv)*bfac*abs(yl(iv))
            idxphi(ix,iy) = iv
            if (isfirstvar==0) ivfirst(ix,iy) = iv
            isfirstvar = 1
            igyl(iv,1) = ix
            igyl(iv,2) = iy
            icnstr(iv) = 0
         endif
  6      continue
 8    continue

      return
      end
c ***** end of subroutine convert ***********
c-----------------------------------------------------------------------
      subroutine convsr_vo (ixl, iyl, yl)
      implicit none
      integer ixl, iyl
      real yl(*)
    
      call convsr_vo1 (ixl, iyl, yl)
      call convsr_vo2 (ixl, iyl, yl)
      end subroutine convsr_vo



      subroutine convsr_vo1 (ixl, iyl, yl)

c ... Converts the yl variables into plasma variables over a restricted
c ... range of indices
c ... This routine only unloads the yl variables into ni, up, etc.,
c ... The other variables are added in the routine convr_auxo


      implicit none

*  -- input arguments
      integer ixl, iyl, inc
      real yl(*)

*  -- local variables
      real ntemp
      integer is, ie, js, je
      integer ifld, id
      integer inegt, inegng, inegni, ixneg, iyneg, ifldneg, igspneg
      #Former Aux module variables
      integer ix,iy,igsp,ix1,ix2
      real t1,t2

c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
      Use(Dim)                 # nx,ny,nhsp,nzsp,nisp,nusp,ngsp
      Use(Xpoint_indices)      # ixpt1,ixpt2
      Use(Math_problem_size)   # neqmx
      Use(Compla)      # ,zi,zeff,zimpc
      Use(Ynorm)       # isflxvar,temp0,nnorm,ennorm,fnorm,n0,n0g
      Use(UEpar)       # itrap_negt,itrap_negng,
                       # isnionxy,isuponxy,isteonxy,istionxy,
                       # isngonxy,isphionxy
      Use(Selec)       # yinc,ixm1,ixp1
      Use(Indexes)
      Use(Comgeo)
      Use(Noggeo)      # fxm,fx0,fxp,fxmy,fxpy
      Use(Gradients)
      Use(Phyvar)      # pi,ev
      Use(Coefeq)      # cngtgx
      Use(Imprad)      # isimpon
      Use(Indices_domain_dcl)   #ixmxbcl,ixmnbcl,iymxbcl,iymnbcl
                                #typebdy,typecn,iv_totbdy
      Use(Indices_domain_dcg)   #isddcon
      Use(Npes_mpi)             #mype
      Use(ParallelEval)
      Use(Share)                #isnonog
 
      integer ifake  #forces Forthon scripts to put implicit none above here

c ... Set mpi indices, etc
CC c_mpi      include 'mpif.h'
c_mpi      integer status(MPI_STATUS_SIZE)
c_mpi      integer ierr

      id = 1
      if(ixl .lt. 0 .or. yinc .ge. 6) then
         is = 1-ixmnbcl
         ie = nx+ixmxbcl
      else
         is = ixl
         ie = ixl
      endif

      if(iyl .lt. 0 .or. yinc .ge. 6) then
         js = 1-iymnbcl
         je = ny+iymxbcl
      else
         js = iyl
         je = iyl
      endif

c... Added the following for OMPPandf1rhs call (added by .J.Guterl)
      if(ixl .lt. 0 .and. iyl.ge.0) then
         js = max(1-iymnbcl,iyl-yinc)
         je = min(ny+iymxbcl,iyl+yinc)
      endif

      if (ParallelPandfCall.gt.0) then
        is = i4
        ie = i8
        js = j4
        je = j8
      end if

        do 20 iy = js, je
          do 19 ix = is, ie  # was nx+1
            ne(ix,iy) = 0.0e0
            nit(ix,iy) = 0.0e0
            nm(ix,iy,1) = 0.0e0
            nz2(ix,iy) = 0.0e0
   19     continue
   20   continue

      inegt = 0    # not used since te & ti have temin eV floor
      inegng = 0
      inegni = 0

        do 50 ifld = 1, nisp  # use full loop even if some isnionxyy=0;
         do 40 iy = js, je
            do 30 ix = is, ie  # was nx+1
           if(isnionxy(ix,iy,ifld).eq.1) then
                 ni(ix,iy,ifld) =  yl(idxn(ix,iy,ifld))*n0(ifld)
                 if (ni(ix,iy,ifld) < 0) then
                   inegni = 1
                   ixneg = ix
                   iyneg = iy
                   ifldneg = ifld
                 endif
               endif
               ne(ix,iy) = ne(ix,iy) + zi(ifld)*ni(ix,iy,ifld)
               if (isupgon(1).eq.1 .and. zi(ifld).eq.0) then
                  ng(ix,iy,1) = ni(ix,iy,ifld)
               else
                  nit(ix,iy) = nit(ix,iy) + ni(ix,iy,ifld)
                  if (isimpon.ge.5 .and. nusp_imp.eq.0) 
     .                 nm(ix,iy,1)=nm(ix,iy,1)+ni(ix,iy,ifld)*mi(ifld) 
                  nz2(ix,iy) = nz2(ix,iy) + ni(ix,iy,ifld)*zi(ifld)**2
               endif
               nm(ix,iy,ifld) = ni(ix,iy,ifld)*mi(ifld)
!               if (isnonog.eq.1) 
!     .                  logni(ix,iy,ifld) = log(ni(ix,iy,ifld))
   30       continue
   40    continue
   50   continue

      do 7 iy = js, je
         do 6 ix = is, ie
            ntemp = ne(ix,iy)
            if(isflxvar .eq. 0) ntemp = nnorm
        if(isteonxy(ix,iy) .eq. 1) then
               te(ix,iy)=yl(idxte(ix,iy))*ennorm/(1.5*ntemp)
               te(ix,iy) = max(te(ix,iy), temin*ev)  #NEW Feb4,2018
               if (isnonog.eq.1) logte(ix,iy) = log(te(ix,iy))
            endif
            do 65 igsp =1, ngsp
             if(isngonxy(ix,iy,igsp) .eq. 1) then
               ng(ix,iy,igsp) = yl(idxg(ix,iy,igsp))*n0g(igsp)
               if (ng(ix,iy,igsp) < 0) then
                     inegng = 1
                     ixneg = ix
                     iyneg = iy
                     igspneg = igsp
               endif
              endif
               if (isnonog.eq.1) logng(ix,iy,igsp) = log(ng(ix,iy,igsp))
               if(istgonxy(ix,iy,igsp) .eq. 1) then
                 ntemp = ng(ix,iy,igsp)
                 if(isflxvar == 0) ntemp = n0g(igsp)
                 tg(ix,iy,igsp) = yl(idxtg(ix,iy,igsp))*ennorm/
     .                                                  (1.5*ntemp)
                 tg(ix,iy,igsp) = max(tg(ix,iy,igsp), tgmin*ev) 
                 if (isnonog.eq.1) logtg(ix,iy,igsp) = log(tg(ix,iy,igsp))
               endif
 65         continue
            ntemp = nit(ix,iy) + cngtgx(1)*ng(ix,iy,1)
            if(isflxvar .eq. 0) ntemp = nnorm
            if(istionxy(ix,iy) .eq. 1) then
               ti(ix,iy)=yl(idxti(ix,iy))*ennorm/(1.5*ntemp)
               ti(ix,iy) = max(ti(ix,iy), temin*ev)
               if (isnonog.eq.1) logti(ix,iy) = log(ti(ix,iy))
            endif
            if(isphionxy(ix,iy) .eq. 1)
     .                           phi(ix,iy) = yl(idxphi(ix,iy))*temp0
    6    continue
    7 continue

      if (inegni .gt. 0 .and. itrap_negni.eq.1) then 
         call remark("***  ni is negative - calculation stopped")
           write(*,*) 'At  ix =', ixneg, ' iy =', iyneg, ' ifld =', ifldneg
         call xerrab("")
      endif
      if (inegng .gt. 0 .and. itrap_negng.eq.1) then 
         call remark("***  ng is negative - calculation stopped")
          write(*,*) 'At  ix =', ixneg, ' iy =', iyneg, ' igsp =', igspneg
         call xerrab("")
      endif
cc Since Te and Ti have temin eV floors, this not used
cc      if (inegt .gt. 0 .and. itrap_negt.eq.1) then 
cc         call xerrab("***  Te or Ti is negative - calculation stopped")
cc      endif

C the message passing is done twice here to get nm for up - very inefficient
c *** Mpi message passing if this is a parallel calculation - only need for 
c *** isflxvar.ne.0
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (isddcon .ge. 1 .and. ixl .lt. 0) then
        call sendbdry(mype+1)
        call recvbdry(mype+1)
      endif
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end


      subroutine convsr_vo2 (ixl, iyl, yl)

c ... Converts the yl variables into plasma variables over a restricted
c ... range of indices
c ... This routine only unloads the yl variables into ni, up, etc.,
c ... The other variables are added in the routine convr_auxo


      implicit none

*  -- input arguments
      integer ixl, iyl, inc
      real yl(*)

*  -- local variables
      real ntemp
      integer is, ie, js, je
      integer ifld, id
      integer inegt, inegng, inegni, ixneg, iyneg, ifldneg, igspneg
      #Former Aux module variables
      integer ix,iy,igsp,ix1,ix2
      real t1,t2

c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
      Use(Dim)                 # nx,ny,nhsp,nzsp,nisp,nusp,ngsp
      Use(Xpoint_indices)      # ixpt1,ixpt2
      Use(Math_problem_size)   # neqmx
      Use(Compla)      # ,zi,zeff,zimpc
      Use(Ynorm)       # isflxvar,temp0,nnorm,ennorm,fnorm,n0,n0g
      Use(UEpar)       # itrap_negt,itrap_negng,
                       # isnionxy,isuponxy,isteonxy,istionxy,
                       # isngonxy,isphionxy
      Use(Selec)       # yinc,ixm1,ixp1
      Use(Indexes)
      Use(Comgeo)
      Use(Noggeo)      # fxm,fx0,fxp,fxmy,fxpy
      Use(Gradients)
      Use(Phyvar)      # pi,ev
      Use(Coefeq)      # cngtgx
      Use(Imprad)      # isimpon
      Use(Indices_domain_dcl)   #ixmxbcl,ixmnbcl,iymxbcl,iymnbcl
                                #typebdy,typecn,iv_totbdy
      Use(Indices_domain_dcg)   #isddcon
      Use(Npes_mpi)             #mype
      Use(ParallelEval)
 
      integer ifake  #forces Forthon scripts to put implicit none above here

c ... Set mpi indices, etc
CC c_mpi      include 'mpif.h'
c_mpi      integer status(MPI_STATUS_SIZE)
c_mpi      integer ierr

      id = 1
      if(ixl .lt. 0 .or. yinc .ge. 6) then
         is = 1-ixmnbcl
         ie = nx+ixmxbcl
      else
         is = ixl
         ie = ixl
      endif

      if(iyl .lt. 0 .or. yinc .ge. 6) then
         js = 1-iymnbcl
         je = ny+iymxbcl
      else
         js = iyl
         je = iyl
      endif

c... Added the following for OMPPandf1rhs call (added by .J.Guterl)
      if(ixl .lt. 0 .and. iyl.ge.0) then
         js = max(1-iymnbcl,iyl-yinc)
         je = min(ny+iymxbcl,iyl+yinc)
      endif

      if (ParallelPandfCall.gt.0) then
        is = i4
        ie = i8
        js = j4
        je = j8
      end if

      do 10 ifld = 1, nusp
         do 9 iy = js, je
            do 8 ix = is, ie
              if(isuponxy(ix,iy,ifld) .eq. 1) then
                ix1 = ixp1(ix,iy)
                ix2 = max(1-ixmnbcl, ixm1(ix,iy))
                t1 = 0.5*( nm(ix2,iy,ifld)+nm(ix, iy,ifld) )
                t2 = 0.5*( nm(ix, iy,ifld)+nm(ix1,iy,ifld) )
                if(isflxvar .eq. 0 .or. isflxvar .eq. 2) then
                  t1 = mi(ifld)*n0(ifld)
                  t2 = mi(ifld)*n0(ifld)
                endif
                up(ix2,iy,ifld) = yl(idxu(ix2,iy,ifld))*fnorm(ifld)/t1
                up(ix,iy,ifld) = yl(idxu(ix,iy,ifld))*fnorm(ifld)/t2
                if(isup1up2==1) then #temp model for up2 if isupon(2)=0
                  up(ix2,iy,2) = rup21*up(ix2,iy,1)
                  up(ix,iy,2) = rup21*up(ix,iy,1)
                endif
              endif  
 8          continue
 9       continue
 10   continue

c *** Mpi message passing if this is a parallel calculation; should only do up here
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (isddcon .ge. 1 .and. ixl .lt. 0) then
        call sendbdry(mype+1)
        call recvbdry(mype+1)
      endif
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      return
      end

c ***** end of subroutines convsr_vo ********



      subroutine convsr_aux (ixl, iyl)
      implicit none
      integer ixl, iyl

      call convsr_aux1(ixl, iyl)
      call convsr_aux2(ixl, iyl)
    
      end subroutine convsr_aux

      subroutine convsr_aux1 (ixl, iyl)

c...  Calculates various plasmas quantities used repeatedly in pandf

      implicit none

*  -- input arguments
      integer ixl, iyl, inc
*  -- local variables
      integer is, ie, js, je
      integer ifld
      #Former Aux module variables
      integer ix,iy,igsp,ix1,ix2
*  -- external functions

      Use(Dim)                 # nx,ny,nhsp,nzsp,nisp,nusp,ngsp
      Use(UEpar)
      Use(Selec)
      Use(Coefeq)
      Use(Compla)      # ,zi,zeff,zimpc
      Use(Phyvar)      # pi,ev
      Use(Imprad)      # isimpon
      Use(ParallelEval)

      if(ixl .lt. 0 .or. yinc .ge. 6) then
         is = 0
         ie = nx+1
      else
         is = ixl
         ie = ixl
      endif

      if(iyl .lt. 0 .or. yinc .ge. 6) then
         js = 0
         je = ny+1
      else
         js = iyl
         je = iyl
      endif

c... Added the following for OMPPandf1rhs call (added by .J.Guterl)
      if(ixl .lt. 0 .and. iyl.ge.0) then
         js = max(0,iyl-yinc)
         je = min(ny+1,iyl+yinc)
      endif

      if (ParallelPandfCall.gt.0) then
        is = i4
        ie = i8
        js = j4
        je = j8
      end if

      do iy = js, je
        do ix = is, ie
          pr(ix,iy) = 0.0e0
        enddo
      enddo

      do ifld = 1, nisp
         do iy = js, je
            do ix = is, ie
               pri(ix,iy,ifld) = ni(ix,iy,ifld) * ti(ix,iy)
               if (ifld.eq.iigsp .and. istgon(1).eq.1) 
     .                  pri(ix,iy,ifld) = 
     .                                     ni(ix,iy,ifld) * tg(ix,iy,1)
               if (zi(ifld).ne.0.) then
                  pr(ix,iy) = pr(ix,iy) + pri(ix,iy,ifld)
               endif
            enddo
          enddo
        enddo

c ... Replace values of ne calculated above by values consistent
c     with those of the INEL average-ion impurity model, if it is
c     being used.  Note that zi(2) passed to rnec will not be used.
      if (isimpon .eq. 3) then
         call xerrab("**** Ave-ion model with isimpon=3 disabled")
ccc         impflag = 1
ccc         do iy = js, je
ccc            do ix = is, ie
ccc               zimpc(ix,iy) = zimp(te(ix,iy))
ccc               ne(ix,iy) = rnec (ni(ix,iy,1), nzsp, ni(ix,iy,nhsp+1),
ccc     .                           te(ix,iy), zi(2), impflag,
ccc     .                           zimpc(ix,iy))
ccc            enddo
ccc         enddo
      endif

      do iy = js, je
        do ix = is, ie
	        pre(ix,iy) = ne(ix,iy) * te(ix,iy)
            pr(ix,iy) = pr(ix,iy) + pre(ix,iy)
            do igsp = 1, ngsp
               if(istgcon(igsp) > -1.e-20) 
     .              tg(ix,iy,igsp) = 
     .                       (1-istgcon(igsp))*rtg2ti(igsp)*ti(ix,iy) + 
     .                          istgcon(igsp)*tgas(igsp)*ev
	       pg(ix,iy,igsp) = ng(ix,iy,igsp)*tg(ix,iy,igsp)
           if(mod(methg,10).eq.6) logpg(ix,iy,igsp) = log(pg(ix,iy,igsp))
           enddo
        enddo
      enddo

      return
      end subroutine convsr_aux1
c-----------------------------------------------------------------------

      subroutine convsr_aux2 (ixl, iyl)

c...  Calculates various plasmas quantities used repeatedly in pandf

      implicit none

*  -- input arguments
      integer ixl, iyl, inc
*  -- local variables
      real ntemp
      integer is, ie, js, je, k, l
      integer ifld, id
      integer impflag, inegt, inegng
      integer jx,ixlp1,ixlp2,ixrm1
      #Former Aux module variables
      integer ix,iy,igsp,ix1,ix2
      real t1,t2
*  -- external functions
      real zimp, rnec, zeffc, intpnog

      Use(Dim)                 # nx,ny,nhsp,nzsp,nisp,nusp,ngsp
      Use(Xpoint_indices)      # ixpt1,ixpt2,iysptrx1,iysptrx2,ixlb,ixrb
      Use(Math_problem_size)   # neqmx
      Use(Compla)      # ,zi,zeff,zimpc
      Use(Ynorm)       # isflxvar,temp0,nnorm,ennorm,fnorm,n0,n0g
      Use(UEpar)       # itrap_negt,itrap_negng
      Use(Selec)       # yinc,ixm1,ixp1
      Use(Indexes)
      Use(Comgeo)      # xcs, yyc
      Use(Noggeo)      # fxm,fx0,fxp,fxmy,fxpy
      Use(Gradients)
      Use(Phyvar)      # pi,ev
      Use(Coefeq)      # cngtgx
      Use(Imprad)      # isimpon
      Use(Interp)      # nis,tis,phis,nxold,nyold
      Use(Share)       # nysol,nyomitmx
      Use(RZ_grid_info)   # rm,zm 
      Use(Volsrc)      # pondpot     
      Use(ParallelEval)

*  -- procedures --
      real interpte,interpti,interpphi,interpni,interppri,interpng,
     .     interpnis,interptis,interpphis,interppg,interptg
c                               # interpolate 2-D array with a 5-point stencil
      interpte(ix,iy,k) = fxm (ix,iy,k)*te(ixm1(ix,iy+k)  ,iy+k  ) + 
     .                    fx0 (ix,iy,k)*te(ix             ,iy+k  ) +
     .                    fxp (ix,iy,k)*te(ixp1(ix,iy+k)  ,iy+k  ) +
     .                    fxmy(ix,iy,k)*te(ixm1(ix,iy+1-k),iy+1-k) +
     .                    fxpy(ix,iy,k)*te(ixp1(ix,iy+1-k),iy+1-k) 
      interpti(ix,iy,k) = fxm (ix,iy,k)*ti(ixm1(ix,iy+k)  ,iy+k  ) + 
     .                    fx0 (ix,iy,k)*ti(ix             ,iy+k  ) +
     .                    fxp (ix,iy,k)*ti(ixp1(ix,iy+k)  ,iy+k  ) +
     .                    fxmy(ix,iy,k)*ti(ixm1(ix,iy+1-k),iy+1-k) +
     .                    fxpy(ix,iy,k)*ti(ixp1(ix,iy+1-k),iy+1-k) 
      interptis(ix,iy,k)=fxm (ix,iy,k)*tis(ixm1(ix,iy+k)  ,iy+k  ) + 
     .                   fx0 (ix,iy,k)*tis(ix             ,iy+k  ) +
     .                   fxp (ix,iy,k)*tis(ixp1(ix,iy+k)  ,iy+k  ) +
     .                   fxmy(ix,iy,k)*tis(ixm1(ix,iy+1-k),iy+1-k) +
     .                   fxpy(ix,iy,k)*tis(ixp1(ix,iy+1-k),iy+1-k) 
      interptg(ix,iy,k,l) = (
     .                   fxm (ix,iy,k)*tg(ixm1(ix,iy+k)  ,iy+k  ,l) + 
     .                   fx0 (ix,iy,k)*tg(ix             ,iy+k  ,l) +
     .                   fxp (ix,iy,k)*tg(ixp1(ix,iy+k)  ,iy+k  ,l) +
     .                   fxmy(ix,iy,k)*tg(ixm1(ix,iy+1-k),iy+1-k,l) +
     .                   fxpy(ix,iy,k)*tg(ixp1(ix,iy+1-k),iy+1-k,l) )
      interpphi(ix,iy,k)= fxm (ix,iy,k)*phi(ixm1(ix,iy+k)  ,iy+k  ) + 
     .                    fx0 (ix,iy,k)*phi(ix             ,iy+k  ) +
     .                    fxp (ix,iy,k)*phi(ixp1(ix,iy+k)  ,iy+k  ) +
     .                    fxmy(ix,iy,k)*phi(ixm1(ix,iy+1-k),iy+1-k) +
     .                    fxpy(ix,iy,k)*phi(ixp1(ix,iy+1-k),iy+1-k) 
      interpphis(ix,iy,k)=fxm (ix,iy,k)*phis(ixm1(ix,iy+k)  ,iy+k  ) + 
     .                    fx0 (ix,iy,k)*phis(ix             ,iy+k  ) +
     .                    fxp (ix,iy,k)*phis(ixp1(ix,iy+k)  ,iy+k  ) +
     .                    fxmy(ix,iy,k)*phis(ixm1(ix,iy+1-k),iy+1-k) +
     .                    fxpy(ix,iy,k)*phis(ixp1(ix,iy+1-k),iy+1-k) 
      interpni(ix,iy,k,l) = exp (
     .                  fxm (ix,iy,k)*log(ni(ixm1(ix,iy+k)  ,iy+k  ,l)) + 
     .                  fx0 (ix,iy,k)*log(ni(ix             ,iy+k  ,l)) +
     .                  fxp (ix,iy,k)*log(ni(ixp1(ix,iy+k)  ,iy+k  ,l)) +
     .                  fxmy(ix,iy,k)*log(ni(ixm1(ix,iy+1-k),iy+1-k,l)) +
     .                  fxpy(ix,iy,k)*log(ni(ixp1(ix,iy+1-k),iy+1-k,l)) )
      interpnis(ix,iy,k,l) = exp (
     .                 fxm (ix,iy,k)*log(nis(ixm1(ix,iy+k)  ,iy+k  ,l)) + 
     .                 fx0 (ix,iy,k)*log(nis(ix             ,iy+k  ,l)) +
     .                 fxp (ix,iy,k)*log(nis(ixp1(ix,iy+k)  ,iy+k  ,l)) +
     .                 fxmy(ix,iy,k)*log(nis(ixm1(ix,iy+1-k),iy+1-k,l)) +
     .                 fxpy(ix,iy,k)*log(nis(ixp1(ix,iy+1-k),iy+1-k,l)) )
      interppri(ix,iy,k,l) = exp (
     .                 fxm (ix,iy,k)*log(pri(ixm1(ix,iy+k)  ,iy+k  ,l)) + 
     .                 fx0 (ix,iy,k)*log(pri(ix             ,iy+k  ,l)) +
     .                 fxp (ix,iy,k)*log(pri(ixp1(ix,iy+k)  ,iy+k  ,l)) +
     .                 fxmy(ix,iy,k)*log(pri(ixm1(ix,iy+1-k),iy+1-k,l)) +
     .                 fxpy(ix,iy,k)*log(pri(ixp1(ix,iy+1-k),iy+1-k,l)) )
      interpng(ix,iy,k,l) = exp (
     .                  fxm (ix,iy,k)*log(ng(ixm1(ix,iy+k)  ,iy+k  ,l)) + 
     .                  fx0 (ix,iy,k)*log(ng(ix             ,iy+k  ,l)) +
     .                  fxp (ix,iy,k)*log(ng(ixp1(ix,iy+k)  ,iy+k  ,l)) +
     .                  fxmy(ix,iy,k)*log(ng(ixm1(ix,iy+1-k),iy+1-k,l)) +
     .                  fxpy(ix,iy,k)*log(ng(ixp1(ix,iy+1-k),iy+1-k,l)) )
      interppg(ix,iy,k,l) = exp (
     .                  fxm (ix,iy,k)*log(pg(ixm1(ix,iy+k)  ,iy+k  ,l)) + 
     .                  fx0 (ix,iy,k)*log(pg(ix             ,iy+k  ,l)) +
     .                  fxp (ix,iy,k)*log(pg(ixp1(ix,iy+k)  ,iy+k  ,l)) +
     .                  fxmy(ix,iy,k)*log(pg(ixm1(ix,iy+1-k),iy+1-k,l)) +
     .                  fxpy(ix,iy,k)*log(pg(ixp1(ix,iy+1-k),iy+1-k,l)) )
 
      id = 1
      if(ixl .lt. 0 .or. yinc .ge. 6) then
         is = 0
         ie = nx+1
      else
         is = ixl
         ie = ixl
      endif

      if(iyl .lt. 0 .or. yinc .ge. 6) then
         js = 0
         je = ny+1
      else
         js = iyl
         je = iyl
      endif

c... Added the following for OMPPandf1rhs call (added by .J.Guterl)
      if(ixl .lt. 0 .and. iyl.ge.0) then
         js = max(0,iyl-yinc)
         je = min(ny+1,iyl+yinc)
      endif

      if (ParallelPandfCall.gt.0) then
        is = i4
        ie = i8
        js = j4
        je = j8
      end if

      do iy = js, je
        do ix = is, ie
          zeff(ix,iy) = 0.0e0
        enddo
      enddo



      do ifld = 1, nisp
         do iy = js, je
            do ix = is, ie
               if (zi(ifld).ne.0.) then
                  zeff(ix,iy)=zeff(ix,iy)+zi(ifld)**2*ni(ix,iy,ifld)
               endif
            enddo
          enddo
        enddo

c ... Replace values of ne calculated above by values consistent
c     with those of the INEL average-ion impurity model, if it is
c     being used.  Note that zi(2) passed to rnec will not be used.
      if (isimpon .eq. 3) then
         call xerrab("**** Ave-ion model with isimpon=3 disabled")
ccc         impflag = 1
ccc         do iy = js, je
ccc            do ix = is, ie
ccc               zimpc(ix,iy) = zimp(te(ix,iy))
ccc               ne(ix,iy) = rnec (ni(ix,iy,1), nzsp, ni(ix,iy,nhsp+1),
ccc     .                           te(ix,iy), zi(2), impflag,
ccc     .                           zimpc(ix,iy))
ccc            enddo
ccc         enddo
      endif

      do iy = js, je
        do ix = is, ie
            zeff(ix,iy) = zeff(ix,iy) / ne(ix,iy)
            znot(ix,iy) = ne(ix,iy)*zeff(ix,iy)/ni(ix,iy,1) - 1
        enddo
      enddo



c ... Set zeff=constant if iszeffcon=1
      if (iszeffcon == 1) then
        do iy = js, je
	  do ix = is, ie
	    zeff(ix,iy) = zeffcon
	  enddo
        enddo
      endif

cccc ... Replace values of zeff calculated above by values consistent
cccc     with those of the INEL average-ion impurity model, if it is
cccc     being used.  Note that zi(2) passed to zeffc will not be used.
ccc      if (isimpon .eq. 3) then
ccc         impflag = 1
ccc         do iy = js, je
ccc            do ix = is, ie
ccc               zeff(ix,iy) = zeffc (nzsp, ni(ix,iy,1), ni(ix,iy,nhsp+1),
ccc     .                              te(ix,iy), zi(2), impflag,
ccc     .                              zimpc(ix,iy))
ccc            enddo
ccc         enddo
ccc      endif

      do 18 iy = js, je
	 inc = isign(max(1,iabs(ie-ixm1(ie,iy))),ie-ixm1(ie,iy))
	 do 17 ix = ixm1(is,iy), min(nx,ie), inc
	    gprx(ix,iy) = 0.0
   17    continue
   18 continue
     
c Tom:  add comments here to explain the indices used on do 20 and 19
      do 20 iy = max(js-1,0), min(ny,je)
	 inc = isign(max(1,iabs(ie-ixm1(ie,js))),ie-ixm1(ie,js))
         do 19 ix = ixm1(is,js), min(nx,ie), inc
            ney0(ix,iy) = 0.0
            ney1(ix,iy) = 0.0
            nity0(ix,iy) = 0.0
            nity1(ix,iy) = 0.0
	    gpry(ix,iy) = 0.0
   19    continue
         ix = ixp1(ie,iy)
         ney0(ix,iy) = 0.0
         ney1(ix,iy) = 0.0
         nity0(ix,iy) = 0.0
         nity1(ix,iy) = 0.0
         gpry(ix,iy) = 0.0
   20 continue

c Tom:  add comments here to explain the indices used on do 21
      do 23 ifld = 1, nisp
         do 22 iy = js, je
	    inc = isign(max(1,iabs(ie-ixm1(ie,iy))),ie-ixm1(ie,iy))
            do 21 ix = ixm1(is,iy), min(nx,ie), inc
	       ix1 = ixp1(ix,iy)
	       gpix(ix,iy,ifld) =
     .               (pri(ix1,iy,ifld)-pri(ix,iy,ifld))*gxf(ix,iy)
              if (zi(ifld).ne.0.) 
     .          gprx(ix,iy) = gprx(ix,iy) +
     .             gpix(ix,iy,ifld)
   21       continue
   22    continue
c...  fix the corners for smooth vycd diamagnetic drifts
cc      gpix(0,0,ifld) = gpix(1,0,ifld)
cc      gpix(0,ny+1,ifld) = gpix(1,ny+1,ifld)
cc      gpix(nx+1,0,ifld) = gpix(nx,0,ifld)
cc      gpix(nx+1,ny+1,ifld) = gpix(nx,ny+1,ifld)
   23 continue

c Tom:  add comments here to explain the indices used on do 25 and 24
      do 26 ifld = 1, nisp
         do 25 iy = max(js-1,0), min(je,ny)
	    inc = isign(max(1,iabs(ie-ixm1(ie,js))),ie-ixm1(ie,js))
            do 24 ix = ixm1(is,js), min(nx,ie), inc
               niy0(ix,iy,ifld) = interpni(ix,iy,0,ifld)
               niy1(ix,iy,ifld) = interpni(ix,iy,1,ifld)
               if (nx==nxold .and. ny==nyold .and. 
     .                                       nis(1,1,1).ne.0.) then
                  niy0s(ix,iy,ifld) = interpnis(ix,iy,0,ifld)
                  niy1s(ix,iy,ifld) = interpnis(ix,iy,1,ifld)
               endif
               nity0(ix,iy) = nity0(ix,iy) + niy0(ix,iy,ifld)
               nity1(ix,iy) = nity1(ix,iy) + niy1(ix,iy,ifld)
               ney0(ix,iy) = ney0(ix,iy) + zi(ifld)*niy0(ix,iy,ifld)
               ney1(ix,iy) = ney1(ix,iy) + zi(ifld)*niy1(ix,iy,ifld)
               priy0(ix,iy,ifld) = interppri(ix,iy,0,ifld)
               priy1(ix,iy,ifld) = interppri(ix,iy,1,ifld) 
               gpiy(ix,iy,ifld) = (priy1(ix,iy,ifld) -
     .                             priy0(ix,iy,ifld))/dynog(ix,iy)
               if (zi(ifld).ne.0.) gpry(ix,iy) = gpry(ix,iy) +
     .            gpiy(ix,iy,ifld)
   24       continue
            ix = ixp1(ie,iy)
            niy0(ix,iy,ifld) = interpni(ix,iy,0,ifld)
            niy1(ix,iy,ifld) = interpni(ix,iy,1,ifld) 
            if (nx==nxold .and. ny==nyold .and. 
     .                                     nis(1,1,1).ne.0.) then
               niy0s(ix,iy,ifld) = interpnis(ix,iy,0,ifld)
               niy1s(ix,iy,ifld) = interpnis(ix,iy,1,ifld) 
            endif
            nity0(ix,iy) = nity0(ix,iy) + niy0(ix,iy,ifld)
            nity1(ix,iy) = nity1(ix,iy) + niy1(ix,iy,ifld)
            ney0(ix,iy) = ney0(ix,iy) + zi(ifld)*niy0(ix,iy,ifld)
            ney1(ix,iy) = ney1(ix,iy) + zi(ifld)*niy1(ix,iy,ifld)
            priy0(ix,iy,ifld) = interppri(ix,iy,0,ifld)
            priy1(ix,iy,ifld) = interppri(ix,iy,1,ifld) 
            gpiy(ix,iy,ifld) = (priy1(ix,iy,ifld) -
     .                          priy0(ix,iy,ifld))/dynog(ix,iy)
            if (zi(ifld).ne.0.) 
     .              gpry(ix,iy) = gpry(ix,iy) +
     .                  gpiy(ix,iy,ifld)
   25    continue
   26 continue

c Tom:  add comments here to explain the indices used on do 264 and 263
      do 264 iy = max(0,js-1), min(je,ny)
         inc = isign(max(1,iabs(ie-ixm1(ie,js))),ie-ixm1(ie,js))
         do 263 ix = ixm1(is,js), min(nx,ie), inc
            tey0(ix,iy) = interpte(ix,iy,0)
            tey1(ix,iy) = interpte(ix,iy,1)
            tiy0(ix,iy) = interpti(ix,iy,0)
            tiy1(ix,iy) = interpti(ix,iy,1)
            phiy0(ix,iy) =interpphi(ix,iy,0) 
            phiy1(ix,iy) =interpphi(ix,iy,1) 
            if (nx==nxold .and. ny==nyold) then
              tiy0s(ix,iy) =interptis(ix,iy,0) 
              tiy1s(ix,iy) =interptis(ix,iy,1) 
              phiy0s(ix,iy) =interpphis(ix,iy,0) 
              phiy1s(ix,iy) =interpphis(ix,iy,1) 
            endif
  263   continue
         ix = ixp1(ie,iy)
         tey0(ix,iy) = interpte(ix,iy,0)
         tey1(ix,iy) = interpte(ix,iy,1)
         tiy0(ix,iy) = interpti(ix,iy,0)
         tiy1(ix,iy) = interpti(ix,iy,1)
         phiy0(ix,iy) =interpphi(ix,iy,0)
         phiy1(ix,iy) =interpphi(ix,iy,1)
         if (nx==nxold .and. ny==nyold) then
           tiy0s(ix,iy) = interptis(ix,iy,0)
           tiy1s(ix,iy) = interptis(ix,iy,1)
           phiy0s(ix,iy) =interpphis(ix,iy,0)
           phiy1s(ix,iy) =interpphis(ix,iy,1)
         endif
  264 continue

c Tom:  add comments here to explain the indices used on do 266 and 265
      do 267 igsp = 1, ngsp
         do 266 iy = max(js-1,0), min(je,ny)
            inc = isign(max(1,iabs(ie-ixm1(ie,js))),ie-ixm1(ie,js))
	    do 265 ix = ixm1(is,js), min(nx,ie), inc
               ngy0(ix,iy,igsp) = interpng(ix,iy,0,igsp)
               ngy1(ix,iy,igsp) = interpng(ix,iy,1,igsp)
               tgy0(ix,iy,igsp) = interptg(ix,iy,0,igsp)
               tgy1(ix,iy,igsp) = interptg(ix,iy,1,igsp)
  265       continue
            ix = ixp1(ie,iy)
            ngy0(ix,iy,igsp) = interpng(ix,iy,0,igsp)
            ngy1(ix,iy,igsp) = interpng(ix,iy,1,igsp)
            tgy0(ix,iy,igsp) = interptg(ix,iy,0,igsp)
            tgy1(ix,iy,igsp) = interptg(ix,iy,1,igsp)
  266    continue
  267 continue
      
C ... Calculate pgy0,1 
        do igsp = 1, ngsp
          do iy = max(js-1,0), min(je,ny)
            inc = isign(max(1,iabs(ie-ixm1(ie,js))),ie-ixm1(ie,js))
            do ix = ixm1(is,js), min(nx,ie), inc
              pgy0(ix,iy,igsp) = interppg(ix,iy,0,igsp)
              pgy1(ix,iy,igsp) = interppg(ix,iy,1,igsp)
            enddo
            ix = ixp1(ie,iy)
            pgy0(ix,iy,igsp) = interppg(ix,iy,0,igsp)
            pgy1(ix,iy,igsp) = interppg(ix,iy,1,igsp)
          enddo
        enddo

      do iy = js, je    # inc index gives stride for jac perturb at cuts
         inc = isign(max(1,iabs(ie-ixm1(ie,iy))),ie-ixm1(ie,iy))
	 do ix = ixm1(is,iy), min(nx,ie), inc
	    ix1 = ixp1(ix,iy)
	    gpex(ix,iy) = (pre(ix1,iy)-pre(ix,iy))*gxf(ix,iy)
	    gtex(ix,iy) = (te(ix1,iy)-te(ix,iy))*gxf(ix,iy)
            gtix(ix,iy) = (ti(ix1,iy)-ti(ix,iy))*gxf(ix,iy)
            gprx(ix,iy) = gprx(ix,iy) + gpex(ix,iy)
            gpondpotx(ix,iy) = (pondpot(ix1,iy)-pondpot(ix,iy))*
     .                                                     gxf(ix,iy)
            if (isphion+isphiofft .eq. 1) then  #ponderom eff in phi if on
               ex(ix,iy) = (phi(ix,iy)-phi(ix1,iy))*gxf(ix,iy)
            endif
         enddo
      ##   reset ex(ixlb,), ex(ixrb,) to avoid using sheath phi(0,) & phi(nx+1,)
	 if (iysptrx < ny) then  # Otherwise this is core-only case w/o plates
	   do jx = 1, nxpt 
             ex(ixlb(jx),iy) = ex(ixlb(jx)+1,iy)
             ex(ixrb(jx),iy) = ex(ixrb(jx)-1,iy)
           enddo
         endif
	 if (islimon.ne.0 .and. iy>=iy_lims) then  # limiter like plate
           ex(ix_lim,iy) = 0.
           ex(ix_lim-1,iy) = ex(ix_lim-2,iy)
           ex(ix_lim+1,iy) = ex(ix_lim+2,iy)
         endif
      enddo

c Tom:  add comments here to explain the indices used on do 30 and 29
      do 30 iy = max(js-1,0), min(ny,je)
         inc = isign(max(1,iabs(ie-ixm1(ie,js))),ie-ixm1(ie,js))
	 do 29 ix = ixm1(is,js), min(nx,ie), inc
            gpey(ix,iy) = (ney1(ix,iy)*tey1(ix,iy) - 
     .                     ney0(ix,iy)*tey0(ix,iy)) / dynog(ix,iy)
            gtey(ix,iy) = (tey1(ix,iy) - tey0(ix,iy)) / dynog(ix,iy)
            gtiy(ix,iy) = (tiy1(ix,iy) - tiy0(ix,iy)) / dynog(ix,iy)
            ey(ix,iy) = - eymask1d(ix,iy)*
     .                      (phiy1(ix,iy) - phiy0(ix,iy)) / dynog(ix,iy)
            gpry(ix,iy) = gpry(ix,iy) + gpey(ix,iy)
   29    continue
         ix = ixp1(ie,iy)
         gpey(ix,iy) = (ney1(ix,iy)*tey1(ix,iy) - 
     .                  ney0(ix,iy)*tey0(ix,iy)) / dynog(ix,iy)
         gtey(ix,iy) = (tey1(ix,iy) - tey0(ix,iy)) / dynog(ix,iy)
         gtiy(ix,iy) = (tiy1(ix,iy) - tiy0(ix,iy)) / dynog(ix,iy)
         ey(ix,iy) = -eymask1d(ix,iy)*
     .                  (phiy1(ix,iy) - phiy0(ix,iy)) / dynog(ix,iy)
         gpry(ix,iy) = gpry(ix,iy) + gpey(ix,iy)
 30   continue

c.... Define vertex values using linear interpolation

c,,,  Note that here we used ixm1(ix,iy) and not ixm1(ix,js) as above
c...  when the iy-loop starts at js-1; seems to work, but should check

         do 32 iy = max(0,js-1), min(ny,je)
	    inc = isign(max(1,iabs(ie-ixm1(ie,iy))),ie-ixm1(ie,iy))
            do 31 ix = ixm1(is,iy), min(nx,ie), inc
               ix1 = ixp1(ix,iy)
               ix2 = ixp1(ix,iy+1)
               phiv(ix,iy) = 0.25*( phi(ix,iy) + phi(ix1,iy) +
     .                              phi(ix,iy+1) + phi(ix2,iy+1) )
               tiv(ix,iy) = 0.25*( ti(ix,iy) + ti(ix1,iy) +
     .                             ti(ix,iy+1) + ti(ix2,iy+1) )
               tev(ix,iy) = 0.25*( te(ix,iy) + te(ix1,iy) +
     .                             te(ix,iy+1) + te(ix2,iy+1) )
c...  add electron contribution to prtv; ion contribution added below
               prev(ix,iy) = 0.25*( pre(ix,iy) + pre(ix1,iy) +
     .                              pre(ix,iy+1) + pre(ix2,iy+1) )

               prtv(ix,iy) = prev(ix,iy)
  31       continue
 32      continue

      do 40 ifld = 1, nisp
         do 39 iy = max(0,js-1), min(ny,je)
	    inc = isign(max(1,iabs(ie-ixm1(ie,iy))),ie-ixm1(ie,iy))
            do 38 ix = ixm1(is,iy), min(nx,ie), inc
               ix1 = ixp1(ix,iy)
               ix2 = ixp1(ix,iy+1)
               priv(ix,iy,ifld) = 0.25*( pri(ix,iy,ifld) + 
     .                           pri(ix1,iy,ifld) + pri(ix,iy+1,ifld) + 
     .                           pri(ix2,iy+1,ifld) )
               if (zi(ifld).ne.0.) 
     .                  prtv(ix,iy) = prtv(ix,iy) + 
     .                                      priv(ix,iy,ifld)

 38         continue
 39      continue
 40   continue

c.... reset the x-point value(s) all the time as it is easier and perhaps
c.... cheaper than checking

      if (nyomitmx < nysol+nyout) then  # otherwise core only - skip

      do jx = 1, nxpt  # loop over mesh regions
         is = ixpt1(jx)
         js = iysptrx1(jx)
         if (jx==1) then
            ie = ixpt2(nxpt) # adjacent cells are in mesh region jx=nxpt
         else
            ie = ixpt2(jx-1) # adjacent cells are in the previous mesh region
         endif
      if (is.lt.0 .or. ie.lt.0 .or. ie.gt.nx) goto 45
c ... Last test (ie.gt.nx) to fix parallel version with mpi - check
      phiv(is,js) = 0.125*( 
     .              phi(is,js)+phi(is+1,js)+phi(is,js+1)+phi(is+1,js+1)+
     .              phi(ie,js)+phi(ie+1,js)+phi(ie,js+1)+phi(ie+1,js+1) )
      phiv(ie,js) = phiv(is,js)
      tiv(is,js) = 0.125*(
     .              ti(is,js)+ti(is+1,js)+ti(is,js+1)+ti(is+1,js+1)+ 
     .              ti(ie,js)+ti(ie+1,js)+ti(ie,js+1)+ti(ie+1,js+1) )
      tiv(ie,js) = tiv(is,js)
      tev(is,js) = 0.125*(
     .              te(is,js)+te(is+1,js)+te(is,js+1)+te(is+1,js+1)+ 
     .              te(ie,js)+te(ie+1,js)+te(ie,js+1)+te(ie+1,js+1) )
      tev(ie,js) = tev(is,js)
      prev(is,js) = 0.125*(
     .              pre(is,js)+pre(is+1,js)+pre(is,js+1)+pre(is+1,js+1)+ 
     .              pre(ie,js)+pre(ie+1,js)+pre(ie,js+1)+pre(ie+1,js+1) )
      prev(ie,js) = prev(is,js)
      prtv(is,js) = prev(is,js)

      do 43 ifld = 1, nisp
      priv(is,js,ifld) = 0.125*(pri(is,js  ,ifld) + pri(is+1,js  ,ifld) 
     .                        + pri(is,js+1,ifld) + pri(is+1,js+1,ifld)
     .                        + pri(ie,js  ,ifld) + pri(ie+1,js  ,ifld)
     .                        + pri(ie,js+1,ifld) + pri(ie+1,js+1,ifld) )
      priv(ie,js,ifld) = priv(is,js,ifld)
      if (zi(ifld).ne.0.) 
     .  prtv(is,js) = prtv(is,js) + priv(is,js,ifld)
 43   continue
      prtv(ie,js) = prtv(is,js)
      enddo  # end do-loop over nxpt mesh regions

      endif  # test on nyomit at top of do loop just above

 45   continue

      return
      end
c ***** end of subroutine convsr_aux ********
c ----------------------------------------------------------------------

      subroutine comp_vertex_vals

c...  Calculates plasmas quantities at cell vertices as diagnostic
c...  Simple averages are used

      implicit none

*  -- local variables
      integer is,ie,js,je,jx,ifld
      #Former Aux module variables
      integer ix,iy,igsp,ix1,ix2

      Use(Dim)            # nx,ny,nhsp,nzsp,nisp,nusp,ngsp
      Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1
      Use(Compla)         # ni,up,..,niv,upv,
      Use(Selec)          # ixp1
      Use(Share)          # nysol,nyomitmx
      
c.. Do all interior cells as 4-pt ave to upper vertex; reset X-point below
      do ix = 1, nx   
        do iy = 1, ny
          ix1 = ixp1(ix,iy)
          ix2 = ixp1(ix,iy+1)
          tev(ix,iy) = 0.25*( te(ix,iy  ) + te(ix1,iy  ) +
     .                        te(ix,iy+1) + te(ix2,iy+1) )
          tiv(ix,iy) = 0.25*( ti(ix,iy  ) + ti(ix1,iy  ) +
     .                        ti(ix,iy+1) + ti(ix2,iy+1) )
          phiv(ix,iy) = 0.25*( phi(ix,iy  ) + phi(ix1,iy  ) +
     .                         phi(ix,iy+1) + phi(ix2,iy+1) )

          do ifld = 1, nisp
            niv(ix,iy,ifld) =0.25*(ni(ix,iy,  ifld) + ni(ix1,iy,  ifld) +
     .                             ni(ix,iy+1,ifld) + ni(ix2,iy+1,ifld) )
            upv(ix,iy,ifld) =0.5*( up(ix,iy,ifld) + up(ix,iy+1,ifld) )
          enddo

          do igsp = 1, ngsp
            ngv(ix,iy,igsp) =0.25*(ng(ix,iy,  igsp) + ng(ix1,iy,  igsp) +
     .                             ng(ix,iy+1,igsp) + ng(ix2,iy+1,igsp) )
          enddo

        enddo
      enddo

c.. Do all x-bdry cells as 2-pt y-ave to upper vertex
      do ix = 0, nx+1, nx+1   
        do iy = 1, ny    # note: corner cells relegated to y-bdry here
          tev(ix,iy) = 0.5*( te(ix,iy) + te(ix,iy+1) )
          tiv(ix,iy) = 0.5*( ti(ix,iy) + ti(ix,iy+1) )
          phiv(ix,iy) = 0.5*( phi(ix,iy) + phi(ix,iy+1) )
          do ifld = 1, nisp
	    niv(ix,iy,ifld) = 0.5*( ni(ix,iy,ifld) + ni(ix,iy+1,ifld) )
	    upv(ix,iy,ifld) = 0.5*( up(ix,iy,ifld) + up(ix,iy+1,ifld) )
          enddo
          do igsp = 1, ngsp
	    ngv(ix,iy,igsp) = 0.5*( ng(ix,iy,igsp) + ng(ix,iy+1,igsp) )
          enddo
        enddo
      enddo

c.. Do all y-bdry cells as 2-pt x-ave to upper vertex
      do ix = 0, nx+1  
        do iy = 0, ny+1, ny+1
          ix1 = ixp1(ix,iy)
          tev(ix,iy) = 0.5*( te(ix,iy) + te(ix1,iy) )
          tiv(ix,iy) = 0.5*( ti(ix,iy) + ti(ix1,iy) )
	  phiv(ix,iy) = 0.5*( phi(ix,iy) + phi(ix1,iy) )
          do ifld = 1, nisp
	    niv(ix,iy,ifld) = 0.5*( ni(ix,iy,ifld) + ni(ix1,iy,ifld) )
            upv(ix,iy,ifld) = up(ix,iy,ifld)
          enddo
          do igsp = 1, ngsp
	    ngv(ix,iy,igsp) = 0.5*( ng(ix,iy,igsp) + ng(ix1,iy,igsp) )
          enddo
        enddo
      enddo

c.. Now reset x-point values; mostly 8-pt ave
      if (nyomitmx < nysol) then  # otherwise core only - skip

        do jx = 1, nxpt  # loop over mesh regions
          is = ixpt1(jx)
          js = iysptrx1(jx)
          if (jx==1) then
            ie = ixpt2(nxpt) # adjacent cells are in mesh region jx=nxpt
          else
            ie = ixpt2(jx-1) # adjacent cells are in the previous mesh region
          endif
          if (is.lt.0 .or. ie.lt.0 .or. ie.gt.nx) return
c ... Last test (ie.gt.nx) to fix parallel version with mpi - check
          tev(is,js) = 0.125*(
     .              te(is,js)+te(is+1,js)+te(is,js+1)+te(is+1,js+1)+ 
     .              te(ie,js)+te(ie+1,js)+te(ie,js+1)+te(ie+1,js+1) )
          tev(ie,js) = tev(is,js)
          tiv(is,js) = 0.125*(
     .              ti(is,js)+ti(is+1,js)+ti(is,js+1)+ti(is+1,js+1)+ 
     .              ti(ie,js)+ti(ie+1,js)+ti(ie,js+1)+ti(ie+1,js+1) )
          tiv(ie,js) = tiv(is,js)
          phiv(is,js) = 0.125*( 
     .              phi(is,js)+phi(is+1,js)+phi(is,js+1)+phi(is+1,js+1)+
     .              phi(ie,js)+phi(ie+1,js)+phi(ie,js+1)+phi(ie+1,js+1) )
          phiv(ie,js) = phiv(is,js)
         
          do ifld = 1, nisp
            niv(is,js,ifld) = 0.125*( ni(is,js,ifld)+ni(is+1,js,ifld)+
     .                                ni(is,js+1,ifld)+ni(is+1,js+1,ifld)+ 
     .                                ni(ie,js,ifld)+ni(ie+1,js,ifld)+
     .                                ni(ie,js+1,ifld)+ni(ie+1,js+1,ifld) )
            niv(ie,js,ifld) = niv(is,js,ifld)
            upv(is,js,ifld) = 0.25*( up(is,js,ifld)+up(is,js+1,ifld)+
     .                               up(ie,js,ifld)+up(ie+1,js+1,ifld) )
          enddo

          do igsp = 1, ngsp
            ngv(is,js,igsp) = 0.125*( ng(is,js,igsp)+ng(is+1,js,igsp)+
     .                                ng(is,js+1,igsp)+ng(is+1,js+1,igsp)+ 
     .                                ng(ie,js,igsp)+ng(ie+1,js,igsp)+
     .                                ng(ie,js+1,igsp)+ng(ie+1,js+1,igsp) )
            ngv(ie,js,igsp) = ngv(is,js,igsp)
          enddo
            
        enddo
      endif

      return
      end
c ***** end of subroutine comp_vertex_vals **************
