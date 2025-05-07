c-----------------------------------------------------------------------
      subroutine curr (iy)

      implicit none

*  -- input arguments
      integer iy

      Use(Dim)      # nx 
      Use(Phyvar)   # qe
      Use(Compla)   # zi
      Use(Comflo)   # fnix,fqx
      Use(Poten)    # bctype,capx
      Use(Coefeq)   # fac2sp

*  -- local variables
      real fqx1, fqx2, t0

*-----------------------------------------------------------------------

*  -- Set the current flow fqx in a flux tube such that
*
*                       fqx1 < fqx < fqx2
*
*     where fqx1 and fqx2 are the ion saturation currents at the
*     left and right boundaries of the flux tube.  We further
*     assume that no current flows across flux surfaces so fqx(ix,iy)
*     is independent of ix on open flux tubes.

      fqx1 = qe * zi(1) * fac2sp*fnix(0,iy,1)
      fqx2 = qe * zi(1) * fac2sp*fnix(nx,iy,1)
      t0 = bctype(iy) * 0.5 * ((fqx1+fqx2) + (fqx1-fqx2)*tanh(capx(iy)))

      call sfill(nx+2, t0, fqx(0,iy), 1)

      return
      end
c-----------------------------------------------------------------------
      subroutine calc_currents
      implicit none

*  -- local variables
      real  nbarx, nbary, sigbarx, sigbary, zfac, temp1, utm, ut0, utp,
     .      ut0n, fqp_old, omgci, zfac0, zfac1, tiy1d, tiy0d,
     .      niy1d, niy0d, phiy1d, phiy0d, fp1, fp2, nzvibtot
      integer iy1, ifld
      #Former Aux module variables
      integer ix,iy,ix1,ix2,ix3,ix4
      real t0,t1
      integer jx,ixl,ixlp1,ixlp2,ixr,ixrm1,ixrm2,iyp2

      Use(Dim)               # nx,ny,nxpt
      Use(Share)             # geometry,nxc,islimon,isudsym
      Use(Xpoint_indices)    # ixlb,ixpt1,ixpt2,ixrb,iysptrx
      Use(Math_problem_size) # neqmx(for arrays not used here)
      Use(Phyvar)            # ev,qe
      Use(Coefeq)            # cfjp2,cfjpy,cftnm,cfvycf,cfqybf,cfq2bf,cfqydt
      Use(Selec)             # i1,i5,i6,j1,j5,j6,ixm1,ixp1,j1p,j2p,j5p,j6p,i3
      Use(Comgeo)            # gx,gy,gxf,gyf,gxc,gyc,sx,sy,rr,isxptx,isxpty,
                             # isixcore
      Use(Compla)            # te,prtv,phi,ne,zeff,up,vycf,netap,niy1,tiy1,
                             # phiy1
      Use(Comtra)            # difutm
      Use(Comflo)            # fqp,fq2,fqx,fqy,fqya,fmity,fnxg,fngxy,
                             # fqyd,fq2d,fqydt
      Use(Ynorm)             # sigbar0
      Use(Poten)             # cthe,sigma1,rsigpl,rsigplcore,cfsigm
      Use(Bfield)            # btot,rbfbt,rbfbt2,bfacyrozh,bfacxrozh
      Use(Gradients)         # gprx,gpry,ey
      Use(RZ_grid_info)      # rm
      Use(UEpar)             # isnewpot,r0slab,isfqpave,rrmin,frfqpn
      Use(Imprad)            # isimpon
      Use(Conduc)            # nucx
      Use(Bcond)             # isexunif,fqpsatlb,fqpsatrb,isextrnp
      Use(Parallv)           # nxg,nyg
      Use(Time_dep_nwt)      # dtreal
      Use(Interp)            # nxold,nyold
      Use(Indices_domain_dcl)# ixmnbcl,ixmxbcl
      Use(Volsrc)            # pondpot

      ifld = 1

************************************************************************
*     Calculate fqp, the parallel current if isimpon.ne.5
************************************************************************
      if (isimpon.eq.5) goto 255

      do 25 iy = j1p, j6p
         do 24 ix = i1, i5
	    iy1 = max(0,iy-1)
            ix1 = ixp1(ix,iy)
            ix3 = ixp1(ix,iy1)
	    t0 = max(te(ix1,iy),temin*ev)
	    t1 = max(te(ix,iy),temin*ev)
            zfac0 = 1. / ( zeff(ix1,iy) * (1.193 - 0.2205*zeff(ix1,iy)
     .                                + 0.0275*zeff(ix1,iy)**2) )
            zfac1 = 1. / ( zeff(ix,iy) * (1.193 - 0.2205*zeff(ix,iy)
     .                                + 0.0275*zeff(ix,iy)**2) )
            if (isfqpave .eq. 0) then   # use linear interpolation for fqp
               zfac = (zfac0*gx(ix1,iy) + zfac1*gx(ix,iy))/
     .                                      (gx(ix1,iy) + gx(ix,iy))
               nbarx = (ne(ix1,iy)*gx(ix1,iy) + ne(ix,iy)*gx(ix,iy))/
     .                                      (gx(ix1,iy) + gx(ix,iy))
               sigbarx = zfac * cfsigm * sigma1 *
     .                   (rr(ix1,iy)*t0**1.5*gx(ix1,iy)
     .                   +rr(ix,iy)*t1**1.5*gx(ix,iy))
     .                   / ((gx(ix1,iy)+gx(ix,iy))*ev**1.5)
            else   # use simple average for fqp components
               zfac = 0.5*(zfac0 + zfac1)
               nbarx = 0.5*(ne(ix1,iy) + ne(ix,iy))
               sigbarx = zfac * cfsigm * sigma1 * rrv(ix,iy) *
     .                                 ( 0.5*(t0+t1)/ev )**1.5
            endif
            netap(ix,iy) = nbarx/sigbarx   # used for frice in pandf
c           temp1 = (gpry(ix,iy) + gpry(ix,iy1) +
c    .               gpry(ix1,iy) + gpry(ix3,iy1))
c           if(isxptx(ix,iy).eq.0) temp1 =
c    .              4.0*(prtv(ix,iy) - prtv(ix,iy1)) * gyc(ix,iy)
            fqp(ix,iy) = (rrv(ix,iy)*sx(ix,iy)*sigbarx*gxf(ix,iy)/qe)*
     .                       ( (pre(ix1,iy) - pre(ix,iy))/nbarx
     .                  - qe * (phi(ix1,iy) - phi(ix,iy))
     .                  + qe * (pondpot(ix1,iy) - pondpot(ix,iy))
     .                  + pondomfpare_use(ix,iy)/(rrv(ix,iy)*nbarx)
     .               + cthe * ( te(ix1,iy) -  te(ix,iy)) )

c           Special calculation for ix-boundary cells; able to retrieve
c           old case with frfqpn=0 ...
            do jx = 1, nxpt  # loop over all mesh regions
               ixl   = ixlb(jx)     # analog of ix=0
               ixlp1 = ixlb(jx)+1   # analog of ix=1
               ixlp2 = ixlb(jx)+2   # analog of ix=2
               ixr   = ixrb(jx)+1   # analog of ix=nx+1
               ixrm1 = ixrb(jx)     # analog of ix=nx
               ixrm2 = ixrb(jx)-1   # analog of ix=nx-1
               if (ix==ixl .and. ixmnbcl==1) then  # at left boundary
                  fqp_old = fqp(ix,iy)
                  nbarx = ne(ixlp1,iy)
                  sigbarx = zfac*cfsigm*sigma1*rrv(ixlp1,iy)*
     .                                         (te(ixlp1,iy)/ev)**1.5
                  fqp(ix,iy) =
     .                  (rrv(ixlp1,iy)*sx(ixlp1,iy)*sigbarx*gxf(ixlp1,iy)/qe)*
     .                       (  (pre(ixlp2,iy)-pre(ixlp1,iy))/nbarx -
     .               qe*(phi(ixlp1,iy)-phi(ixl,iy))*gxf(ixl,iy)/gxf(ixlp1,iy) +
     .                  cthe*(te(ixlp2,iy)- te(ixlp1,iy))  )
                  fqp(ix,iy) = (1.-frfqpn)*fqp_old + frfqpn*fqp(ix,iy)
                  fqpsatlb(iy,jx) = -qe*isfdiax*( ne(ixl,iy)*v2ce(ixl,iy,1)*
     .                         rbfbt(ixl,iy)*sx(ixl,iy) + fdiaxlb(iy,jx) )
                  do ifld = 1, nusp   # note fqp,fqpsat are poloidal proj. of || curr
                     fqpsatlb(iy,jx)=fqpsatlb(iy,jx)-qe*zi(ifld)*ni(ixl,iy,ifld)*
     .                         up(ixl,iy,ifld)*sx(ixl,iy)*rrv(ixl,iy)
                  enddo
                  if (fqp(ixl,iy) < 0.) then #limit to saturation current
		     fp1 = fqp(ixl,iy)
                     fp2 = cffqpsat*fqpsatlb(iy,jx)
                     fqp(ixl,iy)=-( abs(fp1*fp2)**exjbdry/
     .                             (abs(fp1)**exjbdry +
     .                              abs(fp2)**exjbdry) )**(1/exjbdry)
                  endif
               elseif (ix==ixrm1 .and. ixmxbcl==1) then  # at right boundary
                  fqp_old = fqp(ix,iy)
                  nbarx = ne(ixrm1,iy)
                  sigbarx = zfac*cfsigm*sigma1*rrv(ixrm2,iy)*
     .                                         (te(ixrm1,iy)/ev)**1.5
                  fqp(ix,iy) =
     .                  (rrv(ixrm2,iy)*sx(ixrm2,iy)*sigbarx*gxf(ixrm2,iy)/qe)*
     .                       (  (pre(ixrm1,iy)-pre(ixrm2,iy))/nbarx -
     .               qe*(phi(ixr,iy)-phi(ixrm1,iy))*gxf(ixrm1,iy)/gxf(ixrm2,iy) +
     .                  cthe*(te(ixrm1,iy)-te(ixrm2,iy))  )
                  fqp(ix,iy) = (1.-frfqpn)*fqp_old +  frfqpn*fqp(ix,iy)
                  fqpsatrb(iy,jx) = qe*isfdiax*( ne(ixr,iy)*v2ce(ixrm1,iy,1)*
     .                      rbfbt(ixr,iy)*sx(ixrm1,iy) + fdiaxrb(iy,jx) )
                  do ifld = 1, nusp   # note fqp,fqpsat are poloidal proj. of || curr
                     fqpsatrb(iy,jx)=fqpsatrb(iy,jx)+qe*zi(ifld)*ni(ixr,iy,ifld)*
     .                        up(ixrm1,iy,ifld)*sx(ixrm1,iy)*rrv(ixrm1,iy)
                  enddo
                  if (fqp(ixrm1,iy) > 0.) then #limit to saturation current
		     fp1 = fqp(ixrm1,iy)
                     fp2 = cffqpsat*fqpsatrb(iy,jx)
                     fqp(ixrm1,iy)= ( abs(fp1*fp2)**exjbdry/
     .                               (abs(fp1)**exjbdry +
     .                                abs(fp2)**exjbdry) )**(1/exjbdry)
                  endif
               endif  # end if-test for ix boundary cells
            enddo  # end do-loop on number of mesh regions

            if(isudsym==1 .and. ix==nxc)  fqp(ix,iy) = 0.

            if (iy .eq. 1) then
               dphi_iy1(ix) = -fqp(ix,iy)/(rrv(ix,iy)*sx(ix,iy)*
     .                         sigbarx*gxf(ix,iy)) +
     .                         (pre(ix1,iy)-pre(ix,iy))/(qe*nbarx)
            endif
c    .        + cfjpx * sx(ix,iy) * 0.125 * temp1
c    .               * (rbfbt2(ix,iy) + rbfbt2(ix1,iy))
ccc         if (iy.eq.1 .and. isnewpot.eq.1) fqp(ix,1) = fqp(ix,2)
   24    continue
   25 continue

 255  continue   # jump here if isimpon=5
************************************************************************
*     Calculate fq2, the 2-current
************************************************************************
      do 37 iy = j1p, j6p
         do 36 ix = i1, i5
	    iy1 = max(0,iy-1)
            ix1 = ixp1(ix,iy)
            ix2 = ixp1(ix,iy1)
c            temp1 =   (gpry(ix,iy) + gpry(ix,iy1) +
c     .                 gpry(ix1,iy) + gpry(ix2,iy1))
c... sknam: grad P from priv
            temp1 = 4.0*(prtv(ix,iy) - prtv(ix,iy1)) * gyc(ix,iy)
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
            if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .          temp1 = 4.0*(prtv(ix,iy) - prtv(ix,iy-1)) * gyc(ix,iy)
            fq2d(ix,iy) = sx(ix,iy) * 0.25 * temp1 *
     .                         (rbfbt(ix1,iy) + rbfbt(ix,iy)) /
     .                           (btot(ix,iy) + btot(ix1,iy))
            fq2(ix,iy) = cfjp2 * bfacxrozh(ix,iy) * fq2d(ix,iy)
   36    continue
   37 continue

************************************************************************
*     Calculate radial current fqy from toroidal momentum balance eqn
************************************************************************

      do 27 iy = j1p, j5p
         do 26 ix = i1, i6
            ix3 = ixm1(ix,iy)
            ix4 = ixm1(ix,iy+1)
  	    t0 = max(te(ix,iy+1),temin*ev)
            t1 = max(te(ix,iy),temin*ev)
            nbary = (ne(ix,iy+1)*gy(ix,iy+1) + ne(ix,iy)*gy(ix,iy)) /
     .                       (gy(ix,iy+1) + gy(ix,iy))
            zfac = 1. / ( zeff(ix,iy) * (1.193 - 0.2205*zeff(ix,iy)
     .                                + 0.0275*zeff(ix,iy)**2) )
            sigbary = zfac * rsigpl * sigbar0
            if (iy < iysptrx .and. ix > ixpt1(1) .and. ix < ixpt2(1)+1)
     .                        sigbary = sigbary+zfac*rsigplcore*sigbar0
            if (iy.eq.0) sigbary=0.
cc       sigbary = zfac * rsigpl * cfsigm * sigma1 *
cc   .              (rr(ix,iy+1)**2*t0**1.5*gy(ix,iy+1)
c...  the gx in the next line is a bug, but the code likes it - it makes
c...  sigbary very small for iy=0
cc   .              +rr(ix,iy)**2*t1**1.5*gx(ix,iy))
cc   .              / ((gy(ix,iy+1)+gy(ix,iy))*ev**1.5)
c            temp1 = (gprx(ix,iy) + gprx(ix3,iy) +
c     .               gprx(ix,iy+1) + gprx(ix4,iy+1))
c... sknam: grad P from priv
            temp1 = 4.0*(prtv(ix,iy) - prtv(ix3,iy)) * gxc(ix,iy)
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
            if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .          temp1 = 4.0*(prtv(ix,iy) - prtv(ix3,iy)) * gxc(ix,iy)
            fqyae(ix,iy) = ( sy(ix,iy)*sigbary/(dynog(ix,iy)*qe) ) * (
     .        + (ney1(ix,iy)*tey1(ix,iy) - ney0(ix,iy)*tey0(ix,iy))/nbary
     .        - qe * (phiy1(ix,iy) - phiy0(ix,iy)) )
            fqyai(ix,iy) = -(sy(ix,iy)*sigbary/(dynog(ix,iy)*qe*zi(1)) ) * (
     .        + (niy1(ix,iy,1)*tiy1(ix,iy) - niy0(ix,iy,1)*tiy0(ix,iy))/nbary
     .        + qe*zi(1) * (phiy1(ix,iy) - phiy0(ix,iy)) )  # off by default
            fqyao(ix,iy) = cfqyao*(cfqyae*fqyae(ix,iy) + cfqyai*fqyai(ix,iy))
            fqyd(ix,iy) = - sy(ix,iy) * 0.125 * temp1
     .                         * (rbfbt2(ix,iy) + rbfbt2(ix,iy+1))
            nzvibtot = 0.
            do ifld = 1, nisp  #note:vycb=0 if zi(ifld)=0, so ok to include
              nzvibtot = nzvibtot + 0.5*zi(ifld)*(
     .            niy0(ix,iy,ifld)+niy1(ix,iy,ifld) )*vycb(ix,iy,ifld)
            enddo
            fqyb(ix,iy) = qe*sy(ix,iy)*( nzvibtot -
     .                     0.5*(ney0(ix,iy)+ney1(ix,iy))*veycb(ix,iy) )
ccc            fqy(ix,iy) = fqya(ix,iy) + cfydd*fqyd(ix,iy)
   26    continue
   27 continue

      do ifld = 1, nisp
       if (zi(ifld) > 1.e-10) then
         do iy = max(j1p, 1), min(j5p, ny)  #min(j5p, ny-1)
           iyp2 = min(iy+2,ny+1)
           do ix = i1, i6
             ix3 = ixm1(ix,iy)
             ix4 = ixm1(ix,iy-1)
             utm = (4/(btot(ix,iy-1)+btot(ix,iy))**2)*
     .             ( ey(ix,iy-1) - 2*cfgpijr*gpiy(ix,iy-1,ifld)/
     .           (qe*zi(ifld)*(niy1(ix,iy-1,ifld)+niy0(ix,iy-1,ifld))) )

             ix3 = ixm1(ix,iy+1)
             ix4 = ixm1(ix,iy)
             ut0 = (4/(btot(ix,iy)+btot(ix,iy+1))**2)*
     .            ( ey(ix,iy) - 2*cfgpijr*gpiy(ix,iy,ifld)/
     .            (qe*zi(ifld)*(niy1(ix,iy,ifld)+niy0(ix,iy,ifld))) )

             ix3 = ixm1(ix,iyp2)
             ix4 = ixm1(ix,iy+1)
             if(iy < ny) then
               utp = (4/(btot(ix,iy+1)+btot(ix,iyp2))**2)*
     .            ( ey(ix,iy+1) - 2*cfgpijr*gpiy(ix,iy+1,ifld)/
     .           (qe*zi(ifld)*(niy1(ix,iy+1,ifld)+niy0(ix,iy+1,ifld))) )
             else
               utp = 0
             endif
             fmity(ix,iy,ifld) = -0.25*mi(ifld)*
     .                          (difutm(ifld)+dutm_use(ix,iy,ifld))*
     .                         ( (niy1(ix,iy,ifld)+niy0(ix,iy,ifld))*
     .                       (2*r0slab+rm(ix,iy,0)+rm(ix,iy+1,0))*ut0 -
     .                         (niy1(ix,iy-1,ifld)+niy0(ix,iy-1,ifld))*
     .                      (2*r0slab+rm(ix,iy-1,0)+rm(ix,iy,0))*utm )*
     .                                                         gy(ix,iy)
c... Use gy not 1/dynog as diff is btwn 2 cells, not interp pts

             fmity(ix,iy+1,ifld) = -0.25*mi(ifld)*
     .                            (difutm(ifld)+dutm_use(ix,iy+1,ifld))*
     .                        ( (niy1(ix,iy+1,ifld)+niy0(ix,iy+1,ifld))*
     .                      (2*r0slab+rm(ix,iy+1,0)+rm(ix,iyp2,0))*utp -
     .                           (niy1(ix,iy,ifld)+niy0(ix,iy,ifld))*
     .                       (2*r0slab+rm(ix,iy,0)+rm(ix,iy+1,0))*ut0 )*
     .                                                      gy(ix,iy+1)

             omgci = qe*zi(ifld)*b(ix,iy,0)/mi(ifld)
             ix3 = ixm1(ix,iy+1)
             ix4 = ixm1(ix,iy)

          # here fqym is from inertia; fqyn is calc in calc_curr_cx, must
          # follow update of nucx.
             fqymi(ix,iy,ifld)= qe*0.5*
     ,                            (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))*
     .                            (vyce(ix,iy,ifld)+vycp(ix,iy,ifld))*(
     .                    -0.5*( (btot(ix,iy+1)+btot(ix,iyp2))*utp -
     .                           (btot(ix,iy-1)+btot(ix,iy))*utm ) ) *
     .                                  0.5*gyf(ix,iy)*sy(ix,iy)/omgci

          # add fqydt from time-dep inertia term
            tiy1d = (tiy1(ix,iy) - tiy1s(ix,iy))/dtreal
            tiy0d = (tiy0(ix,iy) - tiy0s(ix,iy))/dtreal
            niy1d = (niy1(ix,iy,ifld) - niy1s(ix,iy,ifld))/dtreal
            niy0d = (niy0(ix,iy,ifld) - niy0s(ix,iy,ifld))/dtreal
            phiy1d = (phiy1(ix,iy) - phiy1s(ix,iy))/dtreal
            phiy0d = (phiy0(ix,iy) - phiy0s(ix,iy))/dtreal

c... Next diffs btwn interp pts (niy1&niy0), thus use 1/dynog, not gyf
            fqydti(ix,iy,ifld) = ( -qe*0.5*
     .                             (niy1(ix,iy,ifld)+niy0(ix,iy,ifld))/
     .                   (mi(ifld)*omgci**2*dynog(ix,iy)) )*sy(ix,iy)* (
     .              (tiy1d-tiy0d) +
     .              (tiy1d+tiy0d)*(niy1(ix,iy,ifld) - niy0(ix,iy,ifld))/
     .                 (niy1(ix,iy,ifld)+niy0(ix,iy,ifld)) +
     .                 0.5*(tiy1(ix,iy)+tiy0(ix,iy))*
     .              (niy1d/niy1(ix,iy,ifld) - niy0d/niy0(ix,iy,ifld)) +
     .                 qe*(phiy1d - phiy0d) )

           enddo  #loop over ix
         enddo    #loop over iy
       endif      #test zi>1.e-10
      enddo       #loop over ifld


c ... Sum ion species contributions to fqya, fqym, fqydti
      do iy = max(j1p, 1), min(j5p, ny) ## min(j5p, ny-1)
        do ix = i1, i6
          fqya(ix,iy) = 0.
          fqym(ix,iy) = 0.
          fqydt(ix,iy) = 0.
          do ifld = 1, nisp
            if (zi(ifld) > 1.e-10) then
              fqya(ix,iy)=fqya(ix,iy) + (2/(rm(ix,iy,0)+rm(ix,iy+1,0)))*
     .                        ( (fmity(ix,iy+1,ifld)-fmity(ix,iy,ifld))*
     .                                        gyf(ix,iy)*sy(ix,iy) )
              fqym(ix,iy) = fqym(ix,iy) + fqymi(ix,iy,ifld)
              fqydt(ix,iy) = fqydt(ix,iy) + fqydti(ix,iy,ifld)
            endif
          enddo
        enddo
      enddo

c ... Zero fqya on iy=0,1 and iy=ny,1 as unimportant (use isnewpot=0 bc)
      do ix = i1, i6
         if (isixcore(ix)==1) then
            do iy = 0, nfqya0core
              fqya(ix,iy) = 0.
            enddo
         else
            do iy = 0, nfqya0pf
              fqya(ix,iy) = 0.
            enddo
         endif
         do iy = ny, ny+1-nfqya0ow, -1
           fqya(ix,iy) = 0.
         enddo
      enddo

c ... Sum contributions for fqy; ave old fqyao & fqya with rnewpot
      do iy = j1p, j5p
         do ix = i1, i6
            fqy(ix,iy) = (1.-rnewpot)*fqyao(ix,iy) +
     .                     rnewpot*fqya(ix,iy) + cfqybf*fqyb(ix,iy) +
     .          cfqym*fqym(ix,iy)+cfjpy*bfacyrozh(ix,iy)*fqyd(ix,iy)
            fqygp(ix,iy) = (1.-rnewpot)*fqyao(ix,iy) +
     .                     rnewpot*fqya(ix,iy) +
     .                     cfqym*fqym(ix,iy)+fqyd(ix,iy)

            if (nx==nxold .and. ny==nyold) then
               fqy(ix,iy) = fqy(ix,iy) + cfqydt*fqydt(ix,iy)
            endif
            if (cfvycf .ne. 0.) then  # use for classical Braginskii model
              fqy(ix,iy) = qe*0.5*(niy1(ix,iy,1)+niy0(ix,iy,1))*
     .                      vycf(ix,iy)
            endif
         enddo
      enddo

************************************************************************
*     Calculate fqx, the poloidal current
************************************************************************
      do 47 iy = j1p, j6p
         do 46 ix = i1, i5
            ix1 = ixp1(ix,iy)
            nzvibtot = 0.
            do ifld = 1, nisp  #note:v2cb=0 if zi(ifld)=0, so ok to include
              nzvibtot = nzvibtot + 0.5*zi(ifld)*(
     .                ni(ix,iy,ifld)+ni(ix1,iy,ifld) )*v2cb(ix,iy,ifld)
            enddo
            fqxb(ix,iy) = qe*sx(ix,iy)*( nzvibtot -
     .                  0.5*(ne(ix,iy)+ne(ix1,iy))*ve2cb(ix,iy) )*
     .                         0.5*(rbfbt(ix1,iy)+rbfbt(ix,iy))
	    fqx(ix,iy) =  fqp(ix,iy) + fq2(ix,iy) + cfq2bf*fqxb(ix,iy)
c ...     Force boundary fqx to be uniform; these fqx only for phi B.C.
            do jx = 1, nxpt
               if (ix==ixlb(jx)+1 .and. isexunif*ixmnbcl==1)
     .                     fqx(ixlb(jx),iy) = fqx(ixlb(jx)+1,iy)
               if (ix==ixrb(jx) .and. isexunif*ixmxbcl==1)
     .                     fqx(ixrb(jx),iy) = fqx(ixrb(jx)-1,iy)
            enddo
   46    continue
         if ((isudsym==1.or.(geometry.eq.'dnXtarget'))
     &                           .and. nxc.gt.0) fqx(nxc,iy) = 0.
         if (islimon.ne.0 .and. iy.ge.iy_lims) fqx(ix_lim,iy)=0.
   47 continue

c ... Add anomalous perp vis vy using calc_currents result - awkward,change
      if (cfvyavis > 0.) then
        do ifld = 1, 1  # nfsp  # only good for ifld=1
          do iy = max(j1,2), min(j5,ny-1)
            do ix = max(i1,2), min(i6,nx-1)
              vyavis(ix,iy,ifld) = fqya(ix,iy)*2/(
     .                  qe*(niy1(ix,iy,1)+niy0(ix,iy,1))*sy(ix,iy) )
              vy(ix,iy,ifld) = vy(ix,iy,ifld) + cfvyavis*vyavis(ix,iy,ifld)
            enddo
          enddo
        enddo
      endif          


      return
      end
c ***  End of subroutine calc_currents  **********


c --------------------------------------------------------------------c
      subroutine calc_curr_cx
      implicit none

*  -- local variables
      real  nbarx, nbary, sigbarx, sigbary, zfac, temp1, utm, ut0, utp,
     .      ut0n, fqp_old, omgci
      integer iy1
      #Former Aux module variables
      integer ix,iy,ix1,ix2,ix3,ix4

      Use(Dim)               # ny
      Use(Phyvar)            # ev,qe
      Use(Coefeq)            # cfqyn
      Use(Selec)             # i1,i6,j1p,j5p,j6p,ixm1
      Use(Comgeo)            # gx,gy,gxf,gyf,gxc,gyc,sx,sy
      Use(Compla)            # ngy0,niy0,v2ce,b2cd,mi
      Use(Comflo)            # fqy,fqyn
      Use(RZ_grid_info)      # b
      Use(Conduc)            # nucx,nuneo
      Use(UEpar)             # isnewpot

      do iy = max(j1p, 1), min(j5p, ny-1)
         do ix = i1, i6
            omgci = qe*b(ix,iy,0)/mi(1)
            ix3 = ixm1(ix,iy+1)
            ix4 = ixm1(ix,iy)
            fqyn(ix,iy) = qe*0.125*(
     .              (ngy0(ix,iy,1)+ngy1(ix,iy,1))*nucx(ix,iy,1) +
     .              (niy0(ix,iy,1)+niy1(ix,iy,1))*nuneo )*sy(ix,iy)* (
     .                             v2ce(ix ,iy  ,1)+v2cd(ix ,iy  ,1) +
     .                             v2ce(ix ,iy+1,1)+v2cd(ix ,iy+1,1) +
     .                             v2ce(ix4,iy  ,1)+v2cd(ix4,iy  ,1) +
     .                             v2ce(ix3,iy+1,1)+v2cd(ix3,iy+1,1) )/
     .                                                           omgci
         enddo
      enddo

cc      endif   # for tests on isnewpot

c ... Add contributions to fqy
      do iy = j1p, j5p
         do ix = i1, i6
            fqy(ix,iy) = fqy(ix,iy) + cfqyn*fqyn(ix,iy)
         enddo
      enddo

      return
      end
c ***  End of subroutine calc_cur_cx  ************
c-----------------------------------------------------------------------

      subroutine calc_potential_residuals (neq, yl, yldot)

c...  Calculates the potential equation.  It can not be called alone,
c...  It requires that pandf be called first to evaluate the range
c...  of the calculation in the grid.

      implicit none

*  -- input arguments
      integer neq
      real yl(neq),yldot(neq)

*  -- local variables
      integer jx
      #Former Aux module variables
      integer ix,iy,iv3,ix1,ix2
      real dtuse
      logical isgc,isgc1

      Use(Dim)      # nx,ny,[nhsp,nisp,ngsp(for arrays in Rhsides not used)]
                    # nxpt
      Use(Xpoint_indices)      # ixlb,ixrb,ixpt1,ixpt2,iysptrx
      Use(Math_problem_size)   # neqmx
      Use(UEpar)    # nurlxp
      Use(Selec)    # i2,i5,j2,j5,ixm1,ixp1
      Use(Comgeo)   # vol
      Use(Comflo)   # fqx,fqy
      Use(Rhsides)  # resphi
      Use(Indexes)  # idxphi
      Use(Ynorm)    # temp0,dx0,sigbar0
      Use(Bcond)    # isexunif
      Use(Parallv)  # nxg,nyg
      Use(Compla)   # phi
      Use(Volsrc)   # voljcsor
      Use(Aux)      # ixmp

***********************************************************************
*     Set up the equation for the electrostatic potential phi
***********************************************************************

      do 162 iy = j2p, j5p
         do 161 ix = i2, i5
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            isgc = .false.
            do jx = 1, nxpt
               if ((ix==ixlb(jx)) .or. (ix==(ixrb(jx)+1))) isgc=.true.
            enddo
            if ( isgc ) then # set resphi=0 for guard cells
               resphi(ix,iy) = 0.
            else
               resphi(ix,iy) =  ( nurlxp*dx0**2/sigbar0 ) *
     .                            ( fqx(ix1,iy  ) - fqx(ix,iy)
     .                            + fqy(ix ,iy-1) - fqy(ix,iy)
     .                            + voljcsor(ix,iy) )
            endif
  161    continue
  162 continue

*****************************************************************
*   Potential equation to be solved
*****************************************************************
      do 270 iy = j2p, j5p
         do 260 ix = i2, i5
            iv3 = idxphi(ix,iy)
            isgc = .false.
            isgc1 = .false.
            do jx = 1, nxpt
               if (ix==ixlb(jx) .or. ix==(ixrb(jx)+1)) isgc=.true.
               if (ix==(ixlb(jx)+1) .or. ix==ixrb(jx)) isgc1=.true.
            enddo
            if (isexunif==0) then
               if (.not. isgc) then
                  yldot(iv3) = resphi(ix,iy)/(vol(ix,iy)*temp0)
               endif
            else
               if ((.not. isgc) .and. (.not. isgc1)) then
                  yldot(iv3) = resphi(ix,iy)/(vol(ix,iy)*temp0)
               endif
            endif
  260    continue
ccc         yldot(idxphi(1,iy)) = -nurlxp*(phi(1,iy) - phi(0,iy))/temp0
ccc         yldot(idxphi(nx,iy)) = -nurlxp*(phi(nx,iy) -
ccc     .                                           phi(nx+1,iy))/temp0
  270 continue
cc    If isphicore0=1, eset core potential everywhere to midplane pot
cc    just outside separatrix (phi(ixmp,iysptrx+1) 
      if (isphicore0 == 1) then
        do jx = 1, nxpt
          do iy = 0, iysptrx  #iy=0 & 1 set by BCs
            do ix = ixpt1(1)+1, ixpt2(1)
              iv3 = idxphi(ix,iy)
              yldot(iv3) = -nurlxp*(phi(ix,iy)-phi(ixmp,iysptrx+1))/
     .                                                          temp0
            enddo
          enddo
        enddo
      endif

      return
      end subroutine calc_potential_residuals
c-----------------------------------------------------------------------
      subroutine potent_1dsol

      implicit none

*  -- input arguments
      integer iy

      Use(Dim)         # nx
      Use(Phyvar)      # qe
      Use(Comgeo)      # gx
      Use(Compla)      # te,phi
      Use(Poten)       # kappal,kappar
      Use(Gradients)   # ex
      Use(Selec)       # ixp1
      Use(Share)       # nxc
      Use(Xpoint_indices) # iysptrx,ixpt1,ixpt2

*  -- local variables
      integer ix, ix1, ix2
      real dxf

*-----------------------------------------------------------------------

*  -- This subroutine computes the potential along a flux tube by
*     integrating the poloidal electric field, starting from the
*     left boundary where the potential is defined to be zero for snull.
*     For up/down symm cases (isudsym=1), do inner half from inner plate
*     and then outer half from outer plate.
*     For the core region we have previously set the
*     electric field to zero, so we can integrate over all ix, even
*     for flux tubes in the private flux region.

      if (isudsym == 0) then  # single-null
        do iy = j2, j5
          phi(0,iy) = kappal(iy,1) * te(0,iy) / qe
          do ix = i2, i5
            ix1 = ixp1(ix,iy)
            dxf= 0.5 * (gx(ix,iy)+gx(ix1,iy))/(gx(ix,iy)*gx(ix1,iy))
            phi(ix1,iy) = phi(ix,iy) - ex(ix,iy) * dxf
          enddo
        enddo
      else  # for lower half of up/down sym, sep phi from inner/outer plates
        do iy = j2, j5
          phi(0,iy) = kappal(iy,1) * te(0,iy) / qe
          do ix = 0, nxc-1
            ix1 = ixp1(ix,iy)
            dxf= 0.5 * (gx(ix,iy)+gx(ix1,iy))/(gx(ix,iy)*gx(ix1,iy))
            phi(ix1,iy) = phi(ix,iy) - ex(ix,iy) * dxf
          enddo
        enddo
        do iy = j2, j5
          phi(nx+1,iy) = kappar(iy,1) * te(nx+1,iy) / qe
          do ix = nx, nxc+1, -1
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            dxf= 0.5 * (gx(ix,iy)+gx(ix1,iy))/(gx(ix,iy)*gx(ix1,iy))
            phi(ix,iy) = phi(ix2,iy) + ex(ix,iy) * dxf
          enddo
        enddo
      endif

c ... Impose zero y-gradient for iy < iysptrx in core region and 
c ... for all ny+1; could generalize to double null
      do ix = i1, i6
        if (ix > ixpt1(1) .and. ix < ixpt2(1)+1) then
          do iy = j1, min(iysptrx,j5)
            phi(ix,iy) = phi(ix,iysptrx+1)
          enddo
        endif
        phi(ix,ny+1) = phi(ix,ny)
      enddo

      return
      end
c-----------------------------------------------------------------------

      SUBROUTINE initialize_driftterms
      IMPLICIT NONE
      Use(Selec)
      Use(Compla)
      Use(Phyvar)
      Use(UEpar)
      Use(Conduc)
      Use(Bfield)
      integer iy, ix, ix1
      real teev, nexface

c ... Compute log_lambda
      do iy = j1, j6
        do ix = i1, i6
          ix1 = ixp1(ix,iy)
          teev = 0.5*(te(ix,iy)+te(ix1,iy))/ev
          nexface = 0.5*(ne(ix,iy)+ne(ix1,iy))
          if (islnlamcon == 1) then
            loglambda(ix,iy) = lnlam  # set to constant
          elseif (teev < 50.) then    # Braginskii formula, teev<50
            loglambda(ix,iy) = 23.4-1.15*log10(1.e-6*nexface)+
     .                              3.45*log10(teev)
          else                        #teev > 50
            loglambda(ix,iy) = 25.3-1.15*log10(1.e-6*nexface)+
     .              2.33167537087122D+00*log10(teev)
          endif
          ctaui(ix,iy,1) = 2.1e13*sqrt(mi(1)/mp)/ loglambda(ix,iy)
          ctaue(ix,iy,1) = 3.5e11/loglambda(ix,iy) #both for zi=1
        enddo
      enddo

c ... Calculate collis. factors eta1 and rtaue for the simple Braginski model
      do iy = j1, j6
        do ix = i1, i6
           eta1(ix,iy) = cfeta1*0.3*nm(ix,iy,1)*ti(ix,iy)*
     .                   (1/(qe*btot(ix,iy))) / omgci_taui
           rtaue(ix,iy) = cfrtaue*(1/(qe*btot(ix,iy))) / omgce_taue
           dclass_i(ix,iy) = cfcl_i*eta1(ix,iy)/(0.3*nm(ix,iy,1))
           dclass_e(ix,iy) = cfcl_e*te(ix,iy)*rtaue(ix,iy)
        enddo
      enddo


      END SUBROUTINE initialize_driftterms


      SUBROUTINE calc_driftterms
      IMPLICIT NONE
      Use(Selec)
      Use(Compla)
      Use(Phyvar)
      Use(UEpar)
      Use(Dim)
      Use(Coefeq)
      Use(Bcond)
      Use(Conduc)
      Use(Share)
      Use(Comtra)
      Use(Comgeo)
      Use(Bfield)
      Use(Comflo)
      Use(Noggeo)
      Use(Gradients)
      Use(Indices_domain_dcl)
      Use(Xpoint_indices)
      integer iy, ix, ix1, ifld, iy1, ix6, ix4, ix3, ix2, jx, iyp1, iym1,
     .iy2
      real teev, nexface, t0, t1, dgeyy1, difnimix, dgeyy0, geyy0, geyym,
     .geyyp, grdnv, qion, lambd_ci, lambd_ce, temp1, temp2, temp3, temp4,
     .v2dia
      real ave, etaper
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)
      etaper(ix,iy) = 3.234e-9*loglambda(ix,iy)/(max(te(ix,iy),temin*ev)
     .                                          /(1000.*ev))**(1.5)


************************************************************************
*  Transverse Drifts in y-direction and in 2-direction 
*  (normal to B and y)
************************************************************************
*  ---------------------------------------------------------------------
*  compute drifts
*  ---------------------------------------------------------------------
*  -- loop over species number --

      do ifld = 1, nfsp
c --- If this is the neutral species (zi(ifld).eq.0)) we dont want velocities
        if(zi(ifld) > 1.e-10) then  # if not, skip to end of 100 loop
         qion = zi(ifld)*qe
         do iy = j1, j5
            iyp1 = min(iy+1,ny+1)
            iym1 = max(iy-1,0)
            do ix = i1, i6
              ix3 = ixm1(ix,iy)
              ix4 = ixm1(ix,iy+1)
              temp1 = 
     .                ( ex(ix ,iy) + ex(ix ,iy+1) +
     .                  ex(ix3,iy) + ex(ix4,iy+1) )
c... sknam: grad P from priv, prev
              temp1 = (-4.0)*(phiv(ix,iy) - phiv(ix3,iy))*gxc(ix,iy)
              temp2 = 4.0*(priv(ix,iy,ifld) - priv(ix3,iy,ifld))*gxc(ix,iy) 
              temp3 = 4.0*(prev(ix,iy) - prev(ix3,iy))*gxc(ix,iy) 
 
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
              if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 ) then
                temp1 = (-4.0)*(phiv(ix,iy) - phiv(ix3,iy))*gxc(ix,iy)
                temp2 = 4.0*(priv(ix,iy,ifld) - priv(ix3,iy,ifld))*gxc(ix,iy) 
                temp3 = 4.0*(prev(ix,iy) - prev(ix3,iy))*gxc(ix,iy) 
              endif
c...    Calc collisionality factors nu_s/(1 + nu_s) = 1/(1 + lambda_s)
              lambd_ci = 1e16*(ti(ix,iy)/ev)**2/nit(ix,iy)  # approx
              lambd_ce = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)   # approx
              coll_fi(ix,iy) = cfnus_i/(cfnus_i + (lambd_ci/(lconi(ix,iy))))
              coll_fe(ix,iy) = cfnus_e/(cfnus_e + (lambd_ce/(lcone(ix,iy))))
              vyce(ix,iy,ifld) = 0.125 * temp1
     .                          * ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) )
              vycb(ix,iy,ifld) = ( cfcurv*( 0.5*(ti(ix,iy)+ti(ix,iyp1)) +
     .                mi(ifld)*(0.25*(up(ix,iy,ifld)+up(ix,iyp1,ifld)+
     .                         up(ix3,iy,ifld)+up(ix4,iyp1,ifld)))**2 )*
     .                                       curvrby(ix,iy)/qion +
     .                   cfgradb*0.5*( ti(ix,iy)+ti(ix,iyp1) )*
     .                               gradby(ix,iy)/qion )*coll_fi(ix,iy)
              veycb(ix,iy) = ( -cfcurv*0.5*(te(ix,iy)+te(ix,iyp1))*
     .                                               curvrby(ix,iy)/qe -
     .                   cfgradb*0.5*( te(ix,iy)+te(ix,iyp1) )*
     .                                  gradby(ix,iy)/qe )*coll_fe(ix,iy)

              vycp(ix,iy,ifld) = -0.25 * temp2
     .               * (rbfbt2(ix,iy)+rbfbt2(ix,iy+1)) /
     .            (qion*(niy0(ix,iy,ifld)+niy1(ix,iy,ifld)))
              veycp(ix,iy) =  0.25 * temp3
     .               * (rbfbt2(ix,iy)+rbfbt2(ix,iy+1)) /
     .                    (qe*(ney0(ix,iy)+ney1(ix,iy)))
c...   zero the vy-diamagnetic velocity on the y guard-cell faces
              vycp(ix,0,ifld) = 0.
              vycp(ix,ny,ifld) = 0.
              veycp(ix,0) = 0.
              veycp(ix,ny) = 0.

c...  Precompute radial velocities from fixed BOUT turbulence fluxes
              vy_cft(ix,iy,ifld) = 2*fniyos_use(ix,iy,ifld)/
     .                            (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))
              vyte_cft(ix,iy) = 2*feeyosn_use(ix,iy)/
     .                            (tey0(ix,iy)+tey1(ix,iy))
              vyti_cft(ix,iy) = 2*feiyosn_use(ix,iy)/
     .                            (tiy0(ix,iy)+tiy1(ix,iy))
              vyrd(ix,iy,ifld) = - 2. * gpry(ix,iy) /
     .                       (btot(ix,iy)**2/etaper(ix,iy) +
     .                        btot(ix,iy+1)**2/etaper(ix,iy+1) )
              vydd(ix,iy,ifld) = vcony(ifld) + vy_use(ix,iy,ifld) +
     .                                         vy_cft(ix,iy,ifld) -
     .                         (difpr(ifld) + difp_use(ix,iy,ifld)) * 
     .                    ( 2*gpry(ix,iy)/(pr(ix,iy+1) + pr(ix,iy)) -
     .                      3.0*gtey(ix,iy)/(tey1(ix,iy)+tey0(ix,iy)) )
c ...   Note that the density grad. term for vydd added below
           if (cfrtaue.ne.0.) then  #special classical mom. transfer term
              vycr(ix,iy) = -0.5*(rtaue(ix,iy)+rtaue(ix,iyp1)) * (
     .                        (gpiy(ix,iy,1) + gpey(ix,iy))/
     .                        (0.5*(niy1(ix,iy,1)+niy0(ix,iy,1))) -
     .                          1.5*gtey(ix,iy) )
           endif
           if (cfeta1.ne.0. .and. iy.le.ny-1 .and. iy.gt.0) then  
                                             #special classical vis. term
              geyym = 2*gpiy(ix,iym1,1)/(ney1(ix,iym1)+ney0(ix,iym1))-
     .                qe*ey(ix,iym1)
              geyy0 = 2*gpiy(ix,iy,1)/(ney1(ix,iy)+ney0(ix,iy)) -
     .                qe*ey(ix,iy)
              geyyp = 2*gpiy(ix,iyp1,1)/(ney1(ix,iyp1)+ney0(ix,iyp1))-
     .                qe*ey(ix,iyp1)
              dgeyy0 = (geyy0-geyym)*eta1(ix,iy)*gy(ix,iy)
              dgeyy1 = (geyyp-geyy0)*eta1(ix,iyp1)*gy(ix,iyp1)
              vycf(ix,iy) = 2*(dgeyy1-dgeyy0)*gy(ix,iy) /
     .                                 ( (ney1(ix,iy)+ney0(ix,iy))*
     .                     (qe*0.5*(btot(ix,iy)+btot(ix,iym1)))**2 )
           endif

              diffusivwrk(ix,iy)=fcdif*difni(ifld)+dif_use(ix,iy,ifld)
          end do
        end do

c
c ... Compute diffusive part of radial velocity.
c .. Needs further cleaning; no turbulence model used now TDR 9/1/15
         do iy = j1, j5
            do ix = i1, i6
              difnimix = diffusivwrk(ix,iy)

c ... Alter diffusivity in the SOL by mixing fixed diffusivity
c     with anomalous diffusivity computed in subroutine turb_diffus but
c     reduced by the factor difnit(ifld).  The mixing ratio is given by
c     cdifnit.  Diffusivity is unaltered if difnit(ifld) = 0.
c...MER NOTE: For a full double-null configuration, the SOL is defined to
c...  be the region outside the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
cc              if (difnit(ifld) .gt. 1.e-20 .and. zi(ifld) .eq. 1.
cc     .                                 .and. iy .gt. iysptrx) then
cc                 difnimix = (1. - cdifnit) * 
cc     .                      (fcdif*difni(ifld) + dif_use(ix,iy,ifld)) +
cc     .                               cdifnit * difnit(ifld) * difnimix
cc              endif

              vydd(ix,iy,ifld) = vydd(ix,iy,ifld) 
     .           -1. * difnimix * (
     .            2*(1-isvylog)*( (niy1(ix,iy,ifld) - niy0(ix,iy,ifld)) /
     .              dynog(ix,iy) ) / (niy1(ix,iy,ifld)+niy0(ix,iy,ifld))+
     .              isvylog*(log(niy1(ix,iy,ifld)) - 
     .                            log(niy0(ix,iy,ifld))) /dynog(ix,iy) )

c ... Compute total radial velocity.
              vy(ix,iy,ifld) = cfydd *bfacyrozh(ix,iy) *
     .                                 vycp(ix,iy,ifld) + 
     .                         cfrd  * vyrd(ix,iy,ifld) +
     .                                 vydd(ix,iy,ifld) + 
     .                         cfyef * vyce(ix,iy,ifld) +
     .                         cfybf * vycb(ix,iy,ifld) +
     .                        cfvycf * vycf(ix,iy) + 
     .                        cfvycr * vycr(ix,iy)
c ... Compute radial vel v_grad_P eng eqn terms;cfydd+cfybf=1 or 0
              vygp(ix,iy,ifld) = (cfydd+cfybf)*bfacyrozh(ix,iy) *
     .                                         vycp(ix,iy,ifld) + 
     .                                 cfrd  * vyrd(ix,iy,ifld) +
     .                                         vydd(ix,iy,ifld) + 
     .                                 cfyef * vyce(ix,iy,ifld) +
     .                                cfvycf * vycf(ix,iy) + 
     .                                cfvycr * vycr(ix,iy)
              if (isybdrywd == 1) then  #make vy diffusive in wall cells
                 if (iy==0 .and. matwalli(ix) > 0) then
                    vy(ix,iy,ifld) = vydd(ix,iy,ifld)
                 elseif (iy==ny .and. matwallo(ix) > 0) then
                    vy(ix,iy,ifld) = vydd(ix,iy,ifld)
                 endif
              endif
            enddo  #loop over iy
         enddo     #loop over ix

	    do iy = j1, j6
	      do ix = i1, i6
	      iy1 = max(0,iy-1)            # does iy=0 properly
              iy2 = min(ny+1,iy+1) # use ex*fqx since phi(0,) may be large 
	      ix2 = ixp1(ix,iy)
	      ix4 = ixp1(ix,iy1)
              ix6 = ixp1(ix,iy2)
              do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))*gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                 elseif (ix==ixrb(jx) .and. ixmxbcl==1) then
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))*gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                 else  # not a boundary
                    temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                    temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*
     .                                                       gyc(ix,iy)
                    temp3 = 4.*(prev(ix,iy) - prev(ix,iy1))* gyc(ix,iy)
                    temp4 = 2.5*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
                endif
              enddo  #vis end do-loop over nxpt mesh regions
c...MER NOTE: For a full double-null configuration, the following test will
c...  use the radial index of the innermost separatrix (see iysptrx definition
c...  in subroutine nphygeo)
              if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 ) then
                 temp1 =  (-4.)* (phiv(ix,iy) - phiv(ix,iy1))*gyc(ix,iy)
                 temp2 = 4.*(priv(ix,iy,ifld) - priv(ix,iy1,ifld))*
     .                                                       gyc(ix,iy)
              endif

              v2ce(ix,iy,ifld) = - 0.5 * temp1
     .             / ( btot(ix,iy) + btot(ix2,iy) )
              v2cb(ix,iy,ifld) =(cfcurv*( 0.5*(tiv(ix,iy)+tiv(ix,iy1)) +
     .                 mi(ifld)*up(ix,iy,ifld)**2 )*curvrb2(ix,iy) +
     .                     cfgradb*0.5*( tiv(ix,iy)+tiv(ix,iy1) )*
     .                            gradb2(ix,iy))/qion
              ve2cb(ix,iy) = -(cfcurv*0.5*(tev(ix,iy)+tev(ix,iy1))*
     .                                                 curvrb2(ix,iy) +
     .                          cfgradb*0.5*(tev(ix,iy)+tev(ix,iy1))*
     .                            gradb2(ix,iy))/qe
              v2cd(ix,iy,ifld) = temp2
     .                    / ((btot(ix,iy)+btot(ix2,iy))*qion*
     .                           (ni(ix,iy,ifld)+ni(ix2,iy,ifld)))
              ve2cd(ix,iy,1) = -temp3
     .                    / ((btot(ix,iy)+btot(ix2,iy))*qe*
     .                           (ni(ix,iy,ifld)+ni(ix2,iy,ifld)))
              q2cd(ix,iy,ifld) = (priv(ix,iy,ifld)+priv(ix,iy1,ifld))*temp4
     .                    / ( (btot(ix,iy)+btot(ix2,iy))*qion )

c...  Calculate plate electr diamag flux used to find sheath potential
              do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                 # use ix=ixlb+1 values to avoid BC variations
                    v2dia = -0.5*( gpey(ixlb(jx)+1,iy)+gpey(ixlb(jx)+1,iy1) ) /
     .                          ( btot(ixlb(jx)+1,iy)*qe*ne(ixlb(jx)+1,iy) )
                    fdiaxlb(iy,jx) = ne(ixlb(jx)+1,iy) * sx(ixlb(jx),iy) *
     .                              v2dia * rbfbt(ixlb(jx)+1,iy)
                 endif
                 if (ix==ixrb(jx) .and. ixmxbcl==1) then
                    v2dia = -0.5*( gpey(ixrb(jx),iy)+gpey(ixrb(jx),iy1) ) /
     .                          ( btot(ixrb(jx),iy)*qe*ne(ixrb(jx),iy) )
                    fdiaxrb(iy,jx) = ne(ixrb(jx),iy) * sx(ixrb(jx),iy) *
     .                              v2dia * rbfbt(ixrb(jx),iy)
                 endif
              enddo  # end do-loop over nxpt mesh regions

              v2rd(ix,iy,ifld) = - 2. * gprx(ix,iy) /
     .           ( btot(ix,iy)/(etaper(ix,iy)*rbfbt2(ix,iy)) +
     .             btot(ix2,iy)/(etaper(ix2,iy)*rbfbt2(ix2,iy)) )
              v2dd(ix,iy,ifld) = - 2. * difpr2(ifld) * gprx(ix,iy) /
     .                                ( pr(ix2,iy)/rbfbt(ix2,iy) +
     .                                  pr(ix,iy)/rbfbt(ix,iy) ) -
     .            2. * (fcdif*difni2(ifld) + dif2_use(ix,iy,ifld)) * 
     .                               (ni(ix2,iy,ifld)-ni(ix,iy,ifld)) /
     .                      (ni(ix2,iy,ifld)/(rbfbt(ix2,iy)*gx(ix2,iy))+
     .                          ni(ix,iy,ifld)/(rbfbt(ix,iy)*gx(ix,iy)))
              v2(ix,iy,ifld) = cf2dd * bfacxrozh(ix,iy) *
     .                                 v2cd(ix,iy,ifld) + 
     .                         cfrd  * v2rd(ix,iy,ifld) +
     .                                 v2dd(ix,iy,ifld) + 
     .                         cf2ef * v2ce(ix,iy,ifld) +
     .                         cf2bf * v2cb(ix,iy,ifld)
c ...         Compute v2 for v2x_gradx_P eng terms; cf2dd+cf2bf=1 or 0
              v2xgp(ix,iy,ifld) =  0.5*(rbfbt(ix,iy)+rbfbt(ix2,iy)) * (
     .                 (cf2dd+cf2bf) * bfacxrozh(ix,iy) *
     .                                 v2cd(ix,iy,ifld) + 
     .                         cfrd  * v2rd(ix,iy,ifld) +
     .                                 v2dd(ix,iy,ifld) + 
     .                         cf2ef * v2ce(ix,iy,ifld) )
         if (isnonog.eq.1 .and. iy.le.ny) then
c            grdnv = ( 1/( fym (ix,iy,1)/ni(ix2,iy1,ifld) +
c     .                    fy0 (ix,iy,1)/ni(ix2,iy ,ifld) +
c     .                    fyp (ix,iy,1)/ni(ix2,iy2,ifld) +
c     .                    fymx(ix,iy,1)/ni(ix ,iy1,ifld) +
c     .                    fypx(ix,iy,1)/ni(ix, iy2,ifld) ) -
c     .                1/( fym (ix,iy,0)/ni(ix ,iy1,ifld) +
c     .                    fy0 (ix,iy,0)/ni(ix ,iy ,ifld) +
c     .                    fyp (ix,iy,0)/ni(ix ,iy2,ifld) + 
c     .                    fymx(ix,iy,0)/ni(ix4,iy1,ifld) +
c     .                    fypx(ix,iy,0)/ni(ix6,iy2,ifld) ) )
c     .                                                 / dxnog(ix,iy)
cc            grdnv = ( exp( fym (ix,iy,1)*log(ni(ix2,iy1,ifld)) + 
cc     .                     fy0 (ix,iy,1)*log(ni(ix2,iy ,ifld)) +
cc     .                     fyp (ix,iy,1)*log(ni(ix2,iy2,ifld)) +
cc     .                     fymx(ix,iy,1)*log(ni(ix ,iy1,ifld)) +
cc     .                     fypx(ix,iy,1)*log(ni(ix, iy2,ifld)) ) 
cc     .               -exp( fym (ix,iy,0)*log(ni(ix ,iy1,ifld)) + 
cc     .                     fy0 (ix,iy,0)*log(ni(ix ,iy ,ifld)) +
cc     .                     fyp (ix,iy,0)*log(ni(ix ,iy2,ifld)) +
cc     .                     fymx(ix,iy,0)*log(ni(ix4,iy1,ifld)) +
cc     .                     fypx(ix,iy,0)*log(ni(ix6,iy2,ifld)) ) ) /
cc     .                                                    dxnog(ix,iy)
            grdnv = (    ( fym (ix,iy,1)*log(ni(ix2,iy1,ifld)) + 
     .                     fy0 (ix,iy,1)*log(ni(ix2,iy ,ifld)) +
     .                     fyp (ix,iy,1)*log(ni(ix2,iy2,ifld)) +
     .                     fymx(ix,iy,1)*log(ni(ix ,iy1,ifld)) +
     .                     fypx(ix,iy,1)*log(ni(ix, iy2,ifld)) ) 
     .                  -( fym (ix,iy,0)*log(ni(ix ,iy1,ifld)) + 
     .                     fy0 (ix,iy,0)*log(ni(ix ,iy ,ifld)) +
     .                     fyp (ix,iy,0)*log(ni(ix ,iy2,ifld)) +
     .                     fymx(ix,iy,0)*log(ni(ix4,iy1,ifld)) +
     .                     fypx(ix,iy,0)*log(ni(ix6,iy2,ifld)) ) ) /
     .                                                      dxnog(ix,iy)
            vytan(ix,iy,ifld)=(fcdif*difni(ifld) + dif_use(ix,iy,ifld)) *
     .                                      (grdnv/cos(angfx(ix,iy)) - 
     .                       (log(ni(ix2,iy,ifld)) - log(ni(ix,iy,ifld)))
     .                                                 * gxf(ix,iy) )
            if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
              vytan(ix,iy,ifld) = 0.
            endif
            if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
c             non-physical interface between upper target plates for dnull
              vytan(ix,iy,ifld) = 0.
            endif
         endif

         end do
        end do
         do ix = i1, i6
            vy(ix,ny+1,ifld) = 0.0   
         end do
        else    # test on zi > 1.e-10 to skip whole loop
        endif
        end do  # Giant loop over ifld (species)


c ... Need to calculate new currents (fqp) after saving old & before frice,i
      if(isphion+isphiofft .ne. 1)  then
c ... Add anomalous perp vis vy using calc_currents result - awkward,change
          if (cfvyavis > 0.) then
            do ifld = 1, 1  # nfsp  # only good for ifld=1
              do iy = max(j1,2), min(j5,ny-1)
                do ix = max(i1,2), min(i6,nx-1)
                  vyavis(ix,iy,ifld) = fqya(ix,iy)*2/(
     .                  qe*(niy1(ix,iy,1)+niy0(ix,iy,1))*sy(ix,iy) )
                  vy(ix,iy,ifld) = vy(ix,iy,ifld) + cfvyavis*vyavis(ix,iy,ifld)
                enddo
              enddo
            enddo
          endif          
      endif


      END SUBROUTINE calc_driftterms


