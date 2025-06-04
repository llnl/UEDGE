c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"


      SUBROUTINE calc_plasma_momentum_coeffs
      IMPLICIT NONE
      Use(UEpar)
      Use(Compla)
      Use(Dim)
      Use(Selec)
      Use(Locflux)
      Use(Imprad)
      Use(Comgeo)
      Use(Bfield)
      Use(Coefeq)
      Use(Conduc)
      Use(Xpoint_indices)
      Use(Indices_domain_dcg)
      Use(Npes_mpi)
      Use(Comflo)
      Use(Share)
      Use(Noggeo)
      Use(Phyvar)
      Use(Comtra)
      Use(Rhsides)
      integer ifld, iy, ix, ix1, ix2, ix4, iy1, ix3, ix5, iysepu, k, k1, 
     .k2, ixpt1u, ixpt2u
      real uuv, grdnv, vtn, b_ctr, dbds_m, dbds_p, eta_h0, eta_hm,
     .  eta_hp, drag_1, drag_2, nu_ii, drag_3, mf_path, frac_col, 
     .  t0, t1
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

*****************************************************************
*  Here starts the old MOMBAL_B2
*****************************************************************


*  ---------------------------------------------------------------------
*  loop over all species.
*  ---------------------------------------------------------------------

      do ifld = 1, nusp
      if(isupon(ifld) .ne. 0) then
*     ------------------------------------------------------------------
*     compute the residual.
*     ------------------------------------------------------------------

*  -- evaluate flox and conx --

         do iy = j4omp, j8omp
            flox(0,iy,ifld) = 0.0e0
            conx(0,iy,ifld) = 0.0e0
            do ix = i2omp, i6pomp
               ix1 = ixm1(ix,iy)
               if (isimpon.ge.5 .and. ifld.eq.1) then
                   #up(,,1) is total mass vel, whereas uu(,,i) for each ion
                  uuv =0.5*( (up(ix1,iy,ifld)*rrv(ix1,iy)+
     .                        up(ix,iy,ifld)*rrv(ix,iy)) +
     .                     (v2(ix1,iy,ifld)+v2(ix,iy,ifld))*rbfbt(ix,iy)-
     .                     (vytan(ix1,iy,ifld)+vytan(ix,iy,ifld)) )
               else
                  uuv = 0.5 * (uu(ix1,iy,ifld)+uu(ix,iy,ifld))
               endif
               flox(ix,iy,ifld) = cmfx * nm(ix,iy,ifld) * uuv *
     .                          vol(ix,iy) * gx(ix,iy)
ccc Distance between veloc. cell centers:
               if (isgxvon .eq. 0) then     # dx(ix)=1/gx(ix)
                 conx(ix,iy,ifld) = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .                            * gx(ix,iy)
               elseif (isgxvon .eq. 1) then # dx(ix)=.5/gxf(ix-1) + .5/gxf(ix)
                 conx(ix,iy,ifld) = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .               * 2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
               endif
            end do
         end do

*  -- evaluate floy and cony without averaging over two ix cells --

         do  iy = j1omp1, j5omp
            if (nxpt == 1 .or. iy <= iysptrx1(1)) then
              iysepu = iysptrx1(1)
              if (ndomain > 1) iysepu = iysptrxg(mype+1)  # and ixpt1,2u??
              ixpt1u = ixpt1(1)
              ixpt2u = ixpt2(1)
            else  # nxpt=2 and iy > iysptrx1(1), use second separatrix
              iysepu = iysptrx1(2)
              ixpt1u = ixpt1(2)
              ixpt2u = ixpt2(2)
            endif
            do ix = i4omp, i8omp
               ix2 = ixp1(ix,iy)
               ix4 = ixp1(ix,iy+1)
               cony(ix,iy,ifld) = .5 * syv(ix,iy) *
     .                       (ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                           visy(ix,iy+1,ifld)*gy(ix,iy+1)) +
     .                        ave(visy(ix2,iy,ifld)*gy(ix2,iy),
     .                            visy(ix4,iy+1,ifld) * gy(ix4,iy+1)))
               if (iy==iysepu .and. (ix==ixpt1u .or. ix==ixpt2u)) then
                 cony(ix,iy,ifld) = syv(ix,iy) *
     .                       ( ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                             visy(ix,iy+1,ifld)*gy(ix,iy+1)) )
                 floy(ix,iy,ifld) = (cmfy/2) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) )*
     .                                         vy(ix,iy,ifld)
                 if (ifld==1) then  # add user-specified convection
                   floy(ix,iy,ifld) = floy(ix,iy,ifld)+(cmfy/2)*syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) )*
     .                                                 vyup_use(ix,iy)
                 endif
               elseif(isugfm1side == 1 .and. zi(ifld) == 0.) then
                 floy(ix,iy,ifld) = (cmfy/4) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vy(ix,iy,ifld) + vy(ix,iy,ifld) )
               else
                 floy(ix,iy,ifld) = (cmfy/4) * syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vy(ix,iy,ifld) + vy(ix2,iy,ifld) )
                 if (ifld==1) then  # add user-specified convection
                    floy(ix,iy,ifld) = floy(ix,iy,ifld)+(cmfy/4)*syv(ix,iy) *(
     .                    ave(nm(ix,iy,ifld),nm(ix,iy+1,ifld)) +
     .                    ave(nm(ix2,iy,ifld),nm(ix4,iy+1,ifld))) *
     .                       ( vyup_use(ix,iy) + vyup_use(ix2,iy) )
                 endif
               endif
	     if(ishavisy == 1) then
               cony(ix,iy,ifld) = .5 * syv(ix,iy) *
     .                       (ave(visy(ix,iy,ifld)*gy(ix,iy),
     .                           visy(ix,iy+1,ifld)*gy(ix,iy+1)) +
     .                        ave(visy(ix2,iy,ifld)*gy(ix2,iy),
     .                            visy(ix4,iy+1,ifld) * gy(ix4,iy+1)))
             else
               cony(ix,iy,ifld) = .25 * cfaccony*syv(ix,iy) *
     .                       ( visy(ix,iy,ifld)*gy(ix,iy) +
     .                         visy(ix,iy+1,ifld)*gy(ix,iy+1) +
     .                         visy(ix2,iy,ifld)*gy(ix2,iy) +
     .                         visy(ix4,iy+1,ifld)*gy(ix4,iy+1) )
             endif

            end do
        end do
        end if
        end do


      END SUBROUTINE calc_plasma_momentum_coeffs


      SUBROUTINE calc_plasma_momentum(xc, yc)
      IMPLICIT NONE
      Use(UEpar)
      Use(Compla)
      Use(Dim)
      Use(Selec)
      Use(Locflux)
      Use(Imprad)
      Use(Comgeo)
      Use(Bfield)
      Use(Coefeq)
      Use(Conduc)
      Use(Xpoint_indices)
      Use(Indices_domain_dcg)
      Use(Npes_mpi)
      Use(Comflo)
      Use(Share)
      Use(Noggeo)
      Use(Phyvar)
      Use(Comtra)
      Use(Rhsides)
      integer xc, yc
      integer ifld, iy, ix, ix1, ix2, ix4, iy1, ix3, ix5, iysepu, k, k1, 
     .k2, ixpt1u, ixpt2u
      real uuv, grdnv, vtn, b_ctr, dbds_m, dbds_p, eta_h0, eta_hm,
     .  eta_hp, drag_1, drag_2, nu_ii, drag_3, mf_path, frac_col, 
     .  t0, t1
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

*****************************************************************
*  Here starts the old MOMBAL_B2
*****************************************************************


*  ---------------------------------------------------------------------
*  loop over all species.
*  ---------------------------------------------------------------------

      do ifld = 1, nusp
      if(isupon(ifld) .ne. 0) then
*     ------------------------------------------------------------------
*     compute the residual.
*     ------------------------------------------------------------------
*  -- compute the momentum transport --

         call fd2tra (nx,ny,flox(:,:,ifld),floy(:,:,ifld),conx(:,:,ifld),cony(:,:,ifld),
     .                up(0:nx+1,0:ny+1,ifld),fmix(0:nx+1,0:ny+1,ifld),
     .                fmiy(0:nx+1,0:ny+1,ifld),1, methu)

      if (isnonog .eq. 1) then

c     Compute y-component fmixy of nonorthogonal diffusive momentum flux.
c     The convective component is already already added through uu(ix,iy).
c     Average fym, etc in ix to get staggered velocity-grid values fymv, etc.
c     The density-stencil dxnog has to be averaged as well.
         do iy = j2omp, j5omp
            iy1 = max(iy-1,0)
            do ix = i2omp, i5omp+1    # ixp1(i5,iy)
               ix1 = ixm1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               grdnv = ( 
     .                  fymv (ix,iy,1)*up(ix ,iy1 ,ifld)+
     .                  fy0v (ix,iy,1)*up(ix ,iy  ,ifld)+
     .                  fypv (ix,iy,1)*up(ix ,iy+1,ifld)+
     .                  fymxv(ix,iy,1)*up(ix3,iy1 ,ifld)+
     .                  fypxv(ix,iy,1)*up(ix5,iy+1,ifld)-
     .                  fymv (ix,iy,0)*up(ix3,iy1 ,ifld)-
     .                  fy0v (ix,iy,0)*up(ix1,iy  ,ifld)-
     .                  fypv (ix,iy,0)*up(ix5,iy+1,ifld)-
     .                  fymxv(ix,iy,0)*up(ix ,iy1 ,ifld)-
     .                  fypxv(ix,iy,0)*up(ix ,iy+1,ifld) )*2/
     .                              (dxnog(ix,iy)+dxnog(ix1,iy))
               if (isgxvon .eq. 0) then
                  fmixy(ix,iy,ifld) = cfvisxy(ifld)*visy(ix,iy,ifld) *
     .              ( grdnv/cos(0.5*(angfx(ix1,iy)+angfx(ix,iy))) - 
     .               (up(ix,iy,ifld) - up(ix1,iy,ifld))*gx(ix,iy) ) *
     .              0.5*(sx(ix1,iy)+sx(ix,iy))
               elseif (isgxvon .eq. 1) then
                  fmixy(ix,iy,ifld) = cfvisxy(ifld)*visy(ix,iy,ifld) *
     .              ( grdnv/cos(0.5*(angfx(ix1,iy)+angfx(ix,iy))) - 
     .               (up(ix,iy,ifld) - up(ix1,iy,ifld))*
     .                     ( 2*gxf(ix,iy)*gxf(ix1,iy) /
     .                        (gxf(ix,iy)+gxf(ix1,iy)) ) ) *
     .                     0.5*(sx(ix1,iy)+sx(ix,iy))
               endif
c...  Now flux limit with flalfvgxy if ifld=2
               if (ifld==2) then
                 t0 = max(tg(ix,iy,1),tgmin*ev)
                 vtn = sqrt(t0/mg(1))
                 qfl = flalfvgxya(ix)*0.5*(sx(ix,iy)+sx(ix1,iy))*vtn**2*
     .                                        nm(ix,iy,ifld) + cutlo
                 fmixy(ix,iy,ifld) = fmixy(ix,iy,ifld) /
     .                             sqrt(1+(fmixy(ix,iy,ifld)/qfl)**2)
               endif

            end do
        end do
      endif

c...  Compute viscous drag from nonuniform B-field, then add to smoc
      if (isupdrag .eq. 1 .and. ifld .eq. 1) then
        do iy = j2omp, j5omp
          do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
c ...   First, the short mfp drag 
            b_ctr = 0.5*(btot(ix,iy)+btot(ix2,iy)) 
                    # derviatives dbds_m and dbds_p are one-sided
            dbds_m = (btot(ix,iy) - btot(ix1,iy))*
     .                                      gxf(ix1,iy)*rrv(ix1,iy)
            dbds_p = (btot(ix2,iy) - btot(ix,iy))*
     .                                        gxf(ix,iy)*rrv(ix,iy)
            eta_h0 = visx(ix,iy,1)/b_ctr**2.5
            eta_hm = 0.5*(visx(ix,iy,1)+visx(ix1,iy,1))/
     .                                              btot(ix,iy)**2.5
            eta_hp = 0.5*(visx(ix2,iy,1)+visx(ix,iy,1))/
     .                                              btot(ix2,iy)**2.5
            drag_1 = -2*eta_h0*0.5*(up(ix2,iy,1)-up(ix1,iy,1))*
     .                            (btot(ix2,iy)-btot(ix,iy))*
     .                                   (gxf(ix,iy)*rrv(ix,iy))**2
            drag_2 = up(ix,iy,1)*(eta_hp*dbds_p - eta_hm*dbds_m)*
     .                                        gxf(ix,iy)*rrv(ix,iy)
c ...   now for the trapped particle drag (sloppy)
            nu_ii = ni(ix,iy,1)*(2*mp/mi(1))**0.5/
     .                                   (3e12*(ti(ix,iy)*ev)**1.5)
            drag_3 = -mi(1)*ni(ix,iy,1)*up(ix,iy,1)*nu_ii*frac_pt
            mf_path = (2*ti(ix,iy)/mi(1))**0.5 / nu_ii
            frac_col = 1 / (1 + (mf_path/con_leng)**2)
            smoc(ix,iy,1) = smoc(ix,iy,1) + (0.6666667*frac_col*
     .                      b_ctr**2.5*(drag_1 + drag_2) +
     .                      (1 - frac_col)*drag_3) * volv(ix,iy)
          enddo
        enddo
      endif
      
c...  Now fix the fluxes touching the x-point(s):


         do k = 1, nxpt   # loop over all x-points
           k1 = k      # region argument of ixpt1 that touches this x-point
           k2 = k-1    # region argument of ixpt2 that touches this x-point
           if (k==1) k2 = nxpt

           if (nxpt==2) then      # set ghxpt,gvxpt,sxyxpt for full double null
              if (k==1) then      # this is the lower x-point
                 ghxpt = ghxpt_lower
                 gvxpt = gvxpt_lower
                 sxyxpt = sxyxpt_lower
              elseif (k==2) then  # this is the upper x-point
                 ghxpt = ghxpt_upper
                 gvxpt = gvxpt_upper
                 sxyxpt = sxyxpt_upper
              endif
           endif

	   if( ((2*(yc-iysptrx1(k1))-1)/4 .le. 1) .or. j1 == 0 ) then
           if( ((2*(xc-ixpt1(k1))-1)/4)*((2*(xc-ixpt2(k2))-1)/4).eq.0 .or. 
     .                                                        i1.eq.0 ) then
           if(isnfmiy .eq. 1) then

           fmiy(ixpt1(k1),iysptrx1(k1),ifld) = 0.
           fmiy(ixpt2(k2),iysptrx2(k2),ifld) = 0.
           nixpt(ifld,k1) = 0.125 * ( 
     .           ni(ixpt1(k1),iysptrx1(k1)  ,ifld) + ni(ixpt1(k1)+1,iysptrx1(k1)  ,ifld)
     .         + ni(ixpt1(k1),iysptrx1(k1)+1,ifld) + ni(ixpt1(k1)+1,iysptrx1(k1)+1,ifld)
     .         + ni(ixpt2(k2),iysptrx2(k2)  ,ifld) + ni(ixpt2(k2)+1,iysptrx2(k2)  ,ifld)
     .         + ni(ixpt2(k2),iysptrx2(k2)+1,ifld) + ni(ixpt2(k2)+1,iysptrx2(k2)+1,ifld) )
           visyxpt(ifld,k1) = 0.125 * ( 
     .        visy(ixpt1(k1),iysptrx1(k1)  ,ifld) + visy(ixpt1(k1)+1,iysptrx1(k1)  ,ifld)
     .      + visy(ixpt1(k1),iysptrx1(k1)+1,ifld) + visy(ixpt1(k1)+1,iysptrx1(k1)+1,ifld)
     .      + visy(ixpt2(k2),iysptrx2(k2)  ,ifld) + visy(ixpt2(k2)+1,iysptrx2(k2)  ,ifld)
     .      + visy(ixpt2(k2),iysptrx2(k2)+1,ifld) + visy(ixpt2(k2)+1,iysptrx2(k2)+1,ifld) )
           upxpt(ifld,k1) = 0.25 * (
     .           up(ixpt1(k1),iysptrx1(k1)  ,ifld) + up(ixpt2(k2),iysptrx2(k2)  ,ifld)
     .         + up(ixpt1(k1),iysptrx1(k1)+1,ifld) + up(ixpt2(k2),iysptrx2(k2)+1,ifld) )
           vyvxpt(ifld,k1) = (0.707*0.25) * (
     .           vy(ixpt1(k1)  ,iysptrx1(k1),ifld) - vy(ixpt2(k2)  ,iysptrx2(k2),ifld)
     .         - vy(ixpt1(k1)+1,iysptrx1(k1),ifld) + vy(ixpt2(k2)+1,iysptrx2(k2),ifld) )
           vyhxpt(ifld,k1) = (0.707*0.25) * (
     .         - vy(ixpt1(k1)  ,iysptrx1(k1),ifld) + vy(ixpt2(k2)  ,iysptrx2(k2),ifld)
     .         - vy(ixpt1(k1)+1,iysptrx1(k1),ifld) + vy(ixpt2(k2)+1,iysptrx2(k2),ifld) )
cccMER The convective contributions to fmihxpt and fmivxpt seem to have an
cccMER erroneous multiplicative factor -1/2 (from original code) ???
           fmihxpt(ifld,k1) = cfnfmiy*( - cmfy*(mi(ifld)/2)*sxyxpt*nixpt(ifld,k1)*
     .                     vyhxpt(ifld,k1)*upxpt(ifld,k1)
     .                   - sxyxpt*visyxpt(ifld,k1)*(up(ixpt2(k2),iysptrx2(k2)+1,ifld)
     .                     - up(ixpt1(k1),iysptrx1(k1)+1,ifld))*ghxpt )
           fmivxpt(ifld,k1) = cfnfmiy*(- cmfy*(mi(ifld)/2)*sxyxpt*nixpt(ifld,k1)*
     .                     vyvxpt(ifld,k1)*upxpt(ifld,k1)
     .                   - sxyxpt*visyxpt(ifld,k1)*(up(ixpt2(k2),iysptrx2(k2),ifld)
     .                     - up(ixpt1(k1),iysptrx1(k1),ifld))*gvxpt )
           smoc(ixpt1(k1),iysptrx1(k1)+1,ifld) = smoc(ixpt1(k1),iysptrx1(k1)+1,ifld) 
     .                                - fmihxpt(ifld,k1)
           smoc(ixpt2(k2),iysptrx2(k2)+1,ifld) = smoc(ixpt2(k2),iysptrx2(k2)+1,ifld) 
     .                                + fmihxpt(ifld,k1)
           smoc(ixpt1(k1),iysptrx1(k1)  ,ifld) = smoc(ixpt1(k1),iysptrx1(k1)  ,ifld) 
     .                                - fmivxpt(ifld,k1)
           smoc(ixpt2(k2),iysptrx2(k2)  ,ifld) = smoc(ixpt2(k2),iysptrx2(k2)  ,ifld) 
     .                                + fmivxpt(ifld,k1)
         
           endif # end if-test on isnfmiy
           endif # end if-test on xc
           endif # end if-test on yc

         enddo # end do-loop over nxpt x-points
        end if
        end do


      END SUBROUTINE calc_plasma_momentum

      SUBROUTINE calc_plasma_momentum_residuals
      IMPLICIT NONE
      Use(Selec)
      Use(Compla)
      Use(UEpar)
      Use(Coefeq)
      Use(Rhsides)
      Use(MCN_sources)
      Use(Comgeo)
      Use(Conduc)
      Use(Volsrc)
      Use(Dim)
      Use(Share)
      Use(Comflo)
      Use(Comtra)
      Use(Cfric)
      Use(Wkspace)
      Use(Phyvar)
      integer iy, ix, ix2, ifld, jfld 
      real awoll, t0, t1, awll, tv


      do ifld = 1, nusp
      if(isupon(ifld) .ne. 0) then
*  -- source term and pressure gradient --

         do iy = j2omp, j5omp
            do ix = i2omp, i5omp
               ix2 = ixp1(ix,iy)
               if (zi(ifld) .ne. 0) then  # additions only for charged ions
                  dp1 =  cngmom(ifld)*(1/fac2sp)*
     .                      ( ng(ix2,iy,1)*tg(ix2,iy,1)-
     .                        ng(ix ,iy,1)*tg(ix ,iy,1) ) 
                  resmo(ix,iy,ifld) = 0.
                  resmo(ix,iy,ifld) = 
     .                  smoc(ix,iy,ifld)
     .                + smov(ix,iy,ifld) * up(ix,iy,ifld)
     .                - cfneut * cfneutsor_mi * sx(ix,iy) * rrv(ix,iy) * dp1
     .                - cfneut * cfneutsor_mi * cmwall(ifld)*0.5*(ng(ix,iy,1)+ng(ix2,iy,1))
     .                      * mi(ifld)*up(ix,iy,ifld)*0.5
     .                      *(nucx(ix,iy,1)+nucx(ix2,iy,1))*volv(ix,iy)
     .                + cmneut * cmneutsor_mi * uesor_up(ix,iy,ifld)
     .                + cfmsor*(msor(ix,iy,ifld) + msorxr(ix,iy,ifld)) #### IJ 2017: needs *cfneut for multi-charge state ions & MC neutrals?
     .                + volmsor(ix,iy,ifld)
     .                + cfvisxneov*visvol_v(ix,iy,ifld)
     .                + cfvisxneoq*visvol_q(ix,iy,ifld)
c  Add drag with cold, stationary impurity neutrals
                  if (ifld > nhsp) then
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld) - cfupimpg*
     .                    0.25*mi(ifld)*(ni(ix,iy,ifld)+ni(ix2,iy,ifld))*
     .                     ( nucxi(ix,iy,ifld)+nucxi(ix2,iy,ifld)+
     .                       nueli(ix,iy,ifld)+nueli(ix2,iy,ifld) )*
     .                       up(ix,iy,ifld)*volv(ix,iy)
                  endif
               endif

               if (isupgon(1) .eq. 1) then

c     If we are solving the parallel neutral mom eq. we need different/addtnl
c     source terms. Beware that cngmom and cmwall should be zero so that the
c     main ions do not get coupled twice to the neutrals!
c     The CX, ionization, recomb. friction for the parallel momentum eqs
c     for the neutrals and main ions are included here.
c     Assumes the neutrals have index iigsp and corresponding ions index 1

                  if (ifld .eq. 1) then
c     The main ions, momentum coupling:
                     resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                   + cfneut * cfneutsor_mi * cfupcx*0.25*volv(ix,iy)*
     .                       (nucx(ix,iy,1)+nucx(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                       (up(ix,iy,iigsp)-up(ix,iy,1))
     .                   + cfneut * cfneutsor_mi * 0.25*volv(ix,iy)*
     .                       (  (nuiz(ix,iy,1)+nuiz(ix2,iy,1))*
     .                          (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                          up(ix,iy,iigsp)
     .                        - (nurc(ix,iy,1)+nurc(ix2,iy,1))*
     .                          (nm(ix,iy,1)+nm(ix2,iy,1))*up(ix,iy,1)
     .                       )
                  elseif ((isupgon(1) .eq. 1) .and. ifld .eq. iigsp) then
c     The neutral species, momentum coupling AND other source terms:
                      resmo(ix,iy,iigsp) =   # TR resmo(ix,iy,ifld) #IJ 2016
     .                    - cmneut * cmneutsor_mi * uesor_up(ix,iy,1) 
     .                    -sx(ix,iy) * rrv(ix,iy) * 
     .                       cpgx*( cftiexclg*(ni(ix2,iy,iigsp)*ti(ix2,iy)-
     .                                ni(ix,iy,iigsp)*ti(ix,iy))+
     .                              (1.0-cftiexclg)*
     .                               (ni(ix2,iy,iigsp)*tg(ix2,iy,1)-
     .                                ni(ix,iy,iigsp)*tg(ix,iy,1)) ) 
     .                    -cfupcx*0.25*volv(ix,iy)*
     .                       (nucx(ix,iy,1)+nucx(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                       (up(ix,iy,iigsp)-up(ix,iy,1))
     .                    -0.25*volv(ix,iy)*(
     .                       (nuiz(ix,iy,1)+nuiz(ix2,iy,1))*
     .                       (nm(ix,iy,iigsp)+nm(ix2,iy,iigsp))*
     .                                            up(ix,iy,iigsp)
     .                       -(nurc(ix,iy,1)+nurc(ix2,iy,1))*
     .                       (nm(ix,iy,1)+nm(ix2,iy,1))*up(ix,iy,1) )
                  endif
               endif
            end do
        end do

*  -- divergence of momentum flow --

         if (isnonog.eq.1) then
            do iy = j2omp, j5omp
               do ix = i2omp, i5omp
                  ix2 = ixp1(ix,iy)
c ... IJ 2016/10/10 use cfneutdiv_fmg multiplier for neutrals 
c                 if (ifld .ne. iigsp) then
	          if(zi(ifld) > 1.e-20) then # IJ 2016; depends if ion or neut
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                 + (fmixy(ix2,iy,ifld) - fmixy(ix,iy,ifld))
                  else
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                 + cfneutdiv*cfneutdiv_fmg*(fmixy(ix2,iy,ifld) - fmixy(ix,iy,ifld))
c***	IJ 2017/09/21: Need to add similar fmgxy calculation for MC neutrals on nonorthogonal mesh *** 
                  endif
                end do
            end do
         endif

         do iy = j2omp, j5omp
            do ix = i2omp, i5omp
               ix2 = ixp1(ix,iy)
c IJ 2016/10/10 add cfneutdiv_fmg multiplier for neutrals to control fraction of momentum to add 
c               if (ifld .ne. iigsp) then
c ... IJ 2016 resmo contrib changes if ion or neut
	       if(zi(ifld) > 1.e-20) then
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                                 - (fmix(ix2,iy,ifld) - fmix(ix,iy  ,ifld)
     .                        + fluxfacy*(fmiy(ix ,iy,ifld) - fmiy(ix,iy-1,ifld)) )
               else
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                      - cfneutdiv*cfneutdiv_fmg*(fmix(ix2,iy,ifld) - fmix(ix,iy  ,ifld)
     .                        + fluxfacy*(fmiy(ix ,iy,ifld) - fmiy(ix,iy-1,ifld)) )
                 if(cmneutdiv_fmg .ne. 0.0) then 
                    jfld=1
                    resmo(ix,iy,ifld) = resmo(ix,iy,ifld)
     .                    - cmneutdiv*cmneutdiv_fmg*( (fmgx_ue(ix2,iy,jfld) - fmgx_ue(ix,iy  ,jfld))
     .                        + fluxfacy*(fmgy_ue(ix ,iy,jfld) - fmgy_ue(ix,iy-1,jfld)) )
     .                        * (ni(ix,iy,ifld)*ti(ix,iy))/(ni(ix,iy,ifld)*ti(ix,iy))
                 endif
              endif
            end do
        end do

c  -- Include frictional drag in parallel flow here if isofric=1; otherwise
c  -- it is included in frici from mombal or mombalni

        if (isofric.eq.1 .and. nusp .gt.1) then
*  -- w0 now accumulates friction coefficient --
*     -- set w2 = vol*ti**(-1.5) --
         do iy = j1, j6
           do ix = i1omp, i6omp
             fricnrl(ix,iy,ifld) = 0.  #diagnostic ~ ni*mi*nu*(up1-up2)
             w0(ix,iy) = 0.0e0
             w2(ix,iy) = vol(ix,iy) / (ti(ix,iy)*sqrt(ti(ix,iy)))
           enddo
         enddo

*  -- consider all other species --

         do jfld = 1, nusp
           if (jfld .ne. ifld) then
*     -- common factor in collision frequency --
             awoll = zi(ifld)**2 * zi(jfld)**2 *
     .            (qe**4/(12*pi**2*eps0**2)) *
     .            sqrt (2*pi*mi(ifld)*mi(jfld)/(mi(ifld)+mi(jfld)))

*     -- frictional coupling --
             do iy = j1, j6
               do ix = i1omp, i5omp
                 ix2 = ixp1(ix,iy)
                 t0 = ni(ix,iy,ifld) * ni(ix,iy,jfld) * w2(ix,iy)
                 t1 = ni(ix2,iy,ifld)*ni(ix2,iy,jfld)*w2(ix2,iy)
                 awll = awoll*loglambda(ix,iy)
                 tv  = awll*(t0+t1)/2
                 resmo(ix,iy,ifld) = resmo(ix,iy,ifld) +
     .                        tv * (up(ix,iy,jfld)-up(ix,iy,ifld))
                 fricnrl(ix,iy,ifld) = fricnrl(ix,iy,ifld) +
     .                tv*(up(ix,iy,jfld)-up(ix,iy,ifld))/vol(ix,iy)
                 w0(ix,iy) = w0(ix,iy) - tv
               enddo
             enddo
           endif
         enddo
        endif   # if test on isofric.eq.1
        end if
        end do



      END SUBROUTINE calc_plasma_momentum_residuals


      SUBROUTINE calc_friction(xc)
      IMPLICIT NONE
      Use(Imprad)
      Use(Selec)
      Use(Compla)
      Use(Comgeo)
      Use(Gradients)
      Use(Share)
      Use(Phyvar)
      Use(Comtra)
      Use(Cfric)
      Use(Poten)
      Use(Coefeq)
      Use(Comflo)
      Use(Dim)
      Use(UEpar)
      Use(Xpoint_indices)
      Use(Indices_domain_dcl)
      Use(Timing)
      Use(Bfield)
      Use(Bcond)
      Use(Jacobian_restore)
      Use(Conduc)
      integer xc, yc
      integer iy, ix, ix2, ix1, jx, ifld, ixt0, ixt, ixt1
      real nbarx, ltmax, lmfpe, flxlimf, nexface, 
     .cutlo3, tsimp, argx, ueb
      real tick
      cutlo3 = cutlo**0.3
    
 
c ... Save values returned by Hirshman mombal for Jacobian calc. to
c ... minimize calls - restore the "m" or ix-1 values at the end of pandf
c ... The Jacobian ix loop can then be reduced to only include ix-1 and ix
c ... Suffix "o" refers to "old" value at ix, and suffix "om" means "old" 
c ... value at ix-1.

c ... Calc friction forces from Braginskii; no individ chg-states;isimpon < 5.

      if (isimpon < 5) then
         do iy = j1omp1, j6omp    #iys1, iyf6
            do ix = i1omp, i6omp
               ix2 = ixp1(ix,iy)
               nbarx = 0.5*(ne(ix,iy)+ne(ix2,iy))
               ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
               lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
               flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
               frice(ix,iy) = -cthe*flxlimf*nbarx*rrv(ix,iy)*gtex(ix,iy) +
     .                  cfnetap*qe*netap(ix,iy)*fqp(ix,iy)/sx(ix,iy)
               frici(ix,iy,1) = - frice(ix,iy)
               if (fac2sp .gt. 1.1 .and. nusp .eq. 2) then
                  frici(ix,iy,1) = - frice(ix,iy)/fac2sp
                  frici(ix,iy,2) = - frice(ix,iy)/fac2sp
               endif
            enddo
         enddo
      endif

c ... For use within subroutine mombal, the poloidal electric field is
c     calculated from || Ohms law if isphion = 0.  This field is not
c     intended for use in computing cross-field drifts, so a test of
c     cfyef is included. Both isphiofft=0 or 1 cases included in one loop

      if (isphion .eq. 0) then   # ex calc here assumes no parallel current
         do iy = iys1, iyf6
            do ix = i1omp, i6omp
               ix1 = ix
               do jx = 1, nxpt
                 if (ix==ixlb(jx) .and. ixmnbcl==1) then
                   ix1 = ixlb(jx) + 1
                 elseif (ix==ixrb(jx) .and. ixmxbcl==1) then
                   ix1 = ixrb(jx) - 1 
                 endif
               enddo
               ix2 = ixp1(ix1,iy)
               ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
               lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
               flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
               nexface = 0.5*(ne(ix2,iy)+ne(ix1,iy))
               ex(ix,iy) = (1 - isphiofft) * (
     .                                     -( gpex(ix1,iy)/nexface +
     .                                cthe*flxlimf*gtex(ix1,iy) )/qe -
     .                                              gpondpotx(ix,iy) +
     .                 pondomfpare_use(ix,iy)/(qe*rrv(ix,iy)*nexface) ) + 
     .                      isphiofft * (
     .                            (phi(ix1,iy)-phi(ix2,iy))*gxf(ix1,iy) )
            enddo
         enddo
      endif

c ... Loop over cells (restricted to poloidal slice of box if doing
c     Jacobian), calling mombal if it is providing parallel flow, and
c     taking poloidal projection of parallel flow to get poloidal flow.
c     Unperturbed values of the parallel-flow contribution to uu are
c     saved here so they can be restored below.



c...  TODO: remove xc dependency from here
      do iy = iys1, iyf6
         if (xc .gt. 0) then
            ix = xc
            ix1 = ixm1(ix,iy)
            if (isimpon .eq. 5) then   # Hirshmans reduced-ion approx.
               if (istimingon .eq. 1) tsimp = tick()
               call mombal (ix1,ix,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            elseif(isimpon .eq. 6 .or. isimpon .eq. 7) then # Force balance without inertia
               if (istimingon .eq. 1) tsimp = tick()
               call mombalni (ix1,ix,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            endif
            do ifld = 1, nfsp
               if (ifld .le. nusp) then
                 upi(ix1,iy,ifld) = up(ix1,iy,ifld)
               else
                 do jx = 1, nxpt
                    if ( (ix1==ixlb(jx).and.ixmnbcl==1) .or.
     .                   (ix1==ixrb(jx).and.ixmxbcl==1) ) then
                       # constrain boundary velocity
                       if (zi(ifld) .gt. 1.e-10) then
                          argx = abs((2-2*upi(ix1,iy,ifld)/
     .                                     (upi(ix1,iy,1)+cutlo3))**3)
                          argx = min(20., argx)
                          upi(ix1,iy,ifld) = upi(ix1,iy,1) + 
     .                       (upi(ix1,iy,ifld)-upi(ix1,iy,1))*exp(-argx)
                       endif  # end if-test on zi
                    endif  # end if-test on ix
                 enddo  # end do-loop over nxpt mesh regions
               endif
               uup(ix1,iy,ifld) = rrv(ix1,iy)*upi(ix1,iy,ifld)
            enddo
         endif
         do ix = ixs1, min(ixf6, nx+1-ixmxbcl)
            ix2 = ixp1(ix,iy)
            if (isimpon .eq. 5) then
               if (istimingon .eq. 1) tsimp = tick()
               call mombal (ix,ix2,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            elseif(isimpon .eq. 6 .or. isimpon .eq. 7) then # Force balance without inertia
               if (istimingon .eq. 1) tsimp = tick()
               call mombalni (ix,ix2,iy)
               if (istimingon .eq. 1) call timimpfj (tsimp, xc)
            endif
            do ifld = 1, nfsp
               if (ifld .le. nusp) then
                 upi(ix,iy,ifld) = up(ix,iy,ifld)
               else
                 do jx = 1, nxpt
                    if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                   (ix==ixrb(jx).and.ixmxbcl==1) ) then
                       # constrain boundary velocity
                       if (zi(ifld) .gt. 1.e-10) then
                          argx = abs((2-2*upi(ix,iy,ifld)/
     .                                     (upi(ix,iy,1)+cutlo3))**3)
                          argx = min(20., argx)
                          upi(ix,iy,ifld) = upi(ix,iy,1) + 
     .                         (upi(ix,iy,ifld)-upi(ix,iy,1))*exp(-argx)
                       endif  # end if-test on zi
                    endif  # end if-test on ix
                 enddo  # end do-loop over nxpt mesh regions
               endif
               uup(ix,iy,ifld) = rrv(ix,iy)*upi(ix,iy,ifld)
            enddo
         enddo
      enddo

c ... Add contributions to poloidal velocity from cross-field drifts
c     to those from parallel flow.

      do ifld = 1, nfsp
         do iy = j1omp1, j6omp
            if (i1 .gt. 0) then  # il is initial ix; here for uu(ixm1(i1-1,,)
               ix = i1
               ix1 = ixm1(ix,iy)
               uu(ix1,iy,ifld) = uup(ix1,iy,ifld) +
     .                           0.5 * (rbfbt(ix,iy) + rbfbt(ix1,iy)) *
     .                           v2(ix1,iy,ifld) - vytan(ix1,iy,ifld) -
     .                         difax(ifld) * 0.5 * ( ( 0.5*(
     .                         ni(ix1,iy,ifld)/ni(ix,iy,ifld) +
     .                         ni(ix,iy,ifld)/ni(ix1,iy,ifld)) -1)**2 ) *
     .                        (ni(ix,iy,ifld)-ni(ix1,iy,ifld))*gxf(ix1,iy)
     .                       /(ni(ix,iy,ifld)+ni(ix1,iy,ifld))
               uz(ix1,iy,ifld) = -uup(ix1,iy,ifld)/rrv(ix1,iy)*
     .          0.5*(rbfbt(ix,iy) + rbfbt(ix1,iy)) + sign(1.,b02d(ix,iy))*
     .               (cftef*v2ce(ix1,iy,ifld)+cftdd*v2cd(ix1,iy,ifld))*
     .                                                       rrv(ix1,iy)
            endif
            do ix = i1momp, i6pomp    #now the remainder of the uu(ix,,)
               ix2 = ixp1(ix,iy)
               uu(ix,iy,ifld) = uup(ix,iy,ifld) +
     .                          0.5 * (rbfbt(ix,iy) + rbfbt(ix2,iy)) *
     .                          v2(ix,iy,ifld) - vytan(ix,iy,ifld) -
     .                         difax(ifld) * 0.5 * ( ( 0.5*(
     .                         ni(ix,iy,ifld)/ni(ix2,iy,ifld) +
     .                         ni(ix2,iy,ifld)/ni(ix,iy,ifld)) -1)**2 ) *
     .                        (ni(ix2,iy,ifld)-ni(ix,iy,ifld))*gxf(ix,iy)
     .                       /(ni(ix2,iy,ifld)+ni(ix,iy,ifld))
               uz(ix,iy,ifld) = -uup(ix,iy,ifld)/rrv(ix,iy)*
     .          0.5*(rbfbt(ix,iy) + rbfbt(ix2,iy)) + sign(1.,b02d(ix,iy))*
     .               (cftef*v2ce(ix,iy,ifld)+cftdd*v2cd(ix,iy,ifld))*
     .                                                       rrv(ix,iy)
            enddo
         enddo
      enddo

c...  If upi not from full ||mom eq (e.g.,isimpon=6), set impurity
c...  uu(ixrb,,) & upi(ixrb,,) via generalized Bohm cond.
      if(isimpon > 0) then
        do jx = 1, nxpt
	  ixt0 = ixlb(jx)
          ixt = ixrb(jx)+1 
          ixt1 = ixrb(jx)
	  do ifld = nhsp+1, nfsp
            if(ifld > nusp) then  #species without full mom eqn
	      do iy = j1omp1, j6omp
c ..          first left plate(s)
                if(isfixlb(jx) == 0) then # set upi for left plate
                  cs = csfacrb(ifld,jx)*sqrt( (te(ixt0,iy) +
     .                            csfacti*ti(ixt0,iy))/mi(ifld) )
                  ueb = cfueb*( cf2ef*v2ce(ixt0,iy,ifld)*rbfbt(ixt0,iy) -
     .                            vytan(ixt0,iy,ifld) )/rrv(ixt0,iy)
	          uu(ixt0,iy,ifld) = -rrv(ixt0,iy)*cs
                  upi(ixt0,iy,ifld) = -(cs - ueb)
                endif
c ..          switch to right plate(s)
                if(isfixrb(jx) == 0) then
                  cs = csfacrb(ifld,jx)*sqrt( (te(ixt1,iy) +
     .                            csfacti*ti(ixt1,iy))/mi(ifld) )
                  ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .                            vytan(ixt1,iy,ifld) )/rrv(ixt1,iy)
	          uu(ixt1,iy,ifld) = rrv(ixt1,iy)*cs
	          uu(ixt,iy,ifld) = uu(ixt1,iy,ifld)
                  upi(ixt1,iy,ifld) = cs - ueb
                  upi(ixt,iy,ifld) = upi(ixt1,iy,ifld)
                endif
              enddo
            endif   #checks if ifld > nusp
          enddo
        enddo
      endif         # checks if isimpon > 0


      END SUBROUTINE calc_friction


c-----------------------------------------------------------------------
      subroutine mombal0 (nisp, nhsp, nzsp, minu, ziin,
     .                                      misotope, natomic, nchstate)
c ... Compute 'misotope', 'nchstate', and 'natomic', and allocate memory
c     for arrays used in subroutine mombal.
 
      implicit none

c ... Input arguments:
      integer nisp      # total number of ion species
      integer nhsp      # total number of hydrogenic ion species
      integer nzsp(ngspmx-1)   # number of charge states for each imp isotope
      real minu(nisp)   # mass (in amu) of ion species
      real ziin(nisp)   # charge (in units of e) of ion species

c ... Output arguments:
      integer misotope  # total number of isotopes (including electrons)
      integer natomic(*)   # maximum charge state of each isotope
      integer nchstate     # maximum charge state among all isotopes

c ... Local variables:
      integer misa, ifld, jz
	  
c ... Loop over ion species, looking for change to a new isotope, and
c     finding maximum charge state.
      natomic(1) = 1   # electrons are "isotope 1"
      nchstate = 0
      misa = 2
      do ifld = 1, nhsp
         natomic(misa) = max(nint(ziin(ifld)), 1)   # must be .ge. 1
         nchstate = max(nchstate, natomic(misa)) 
         if (ifld .eq. nhsp) go to 50
         if (minu(ifld+1) .ne. minu(ifld)) misa = misa + 1
      enddo
 50   misotope = misa
      do jz = 1, ngspmx-1
         if (nzsp(jz)==0) break
         misotope = misotope + 1
	 if (misotope .gt. MXMISO) then
           call remark("subroutine mombal0 error: ")
           call remark("To avoid write out-of-bounds for array natomic")
           call remark("increase the value of MXMISO and recompile.")
	       call xerrab("")
         endif
         natomic(misotope) = nzsp(jz)
         nchstate = max(nchstate, natomic(misotope)) 
      enddo

c ... Allocate memory for arrays used in subroutine mombal.
      call gallot("Reduced_ion_interface", 0)

      return
      end
c****** end of subroutine mombal0 ************
c-----------------------------------------------------------------------
      subroutine mombal (ix,ix1,iy)
c ... Prepare information needed to call Steve Hirshman reduced-ion
c     momentum-balance routine, and distribute results from it.
c     Note that results are computed at a poloidal density-cell face.
c     The parallel flow velocities are returned in arrays up and upe.

      implicit none

c ... Input arguments:
      integer ix    # poloidal index of density cell to left of face
      integer ix1   # poloidal index of density cell to right of face
      integer iy    # radial index of density cell

c ... Common blocks:
      Use(Dim)          # nx,ny,nisp
      Use(Selec)        # isupgoon
      Use(Comgeo)       # sx,rrv
      Use(UEpar)        # lnlam,isupgon
      Use(Cfric)        # frice,frici
      Use(Compla)       # ni,up,upe,upi,te,ti,ng,ne
      Use(Gradients)    # gpix,gtix,ex,gtex,gpex
      Use(Comflo)       # fqp
      Use(UEint)        # minu,ziin
      Use(Reduced_ion_interface) # misotope,nchstate,natomic,
                                 # amu,tempa,qneut,uneut,den,gradp,
                                 # gradt,friction,nuion,nurec,qcond,ucond
                                 # friccomp
      Use(Imprad)       # ismctab
      Use(Phyvar)       # ev
      Use(Comtra)       # fricflf,sigvi_floor
      Use(Share)        # cutlo

c ... External functions:
      real rra, rsa

c ... Local variables:
      integer ifld, misa, nz
      integer ldir
      real rdum, dloglam, epar, parcurrent, umass, lmfpe, lmfpi,
     .     ltmax, tif, flxlimf, umassni, massni

c ... Determine a flux-limit factor for all input terms by finding the min
c ... scale length. First consider the electrons
      den(1,1) = 0.5 * (ne(ix,iy) + ne(ix1,iy))
      tempa(1) = 0.5 * (te(ix,iy) + te(ix1,iy))
      tif = 0.5 * (ti(ix,iy) + ti(ix1,iy))
      ltmax = min( abs(tempa(1)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .             abs(tif/(rrv(ix,iy)*gtix(ix,iy) + cutlo)),
     .             abs(den(1,1)*tempa(1)/
     .                          (rrv(ix,iy)*gpex(ix,iy) + cutlo)) )
      lmfpe = 1.e16*((tempa(1)/ev)**2/den(1,1)) #Approx Coulomb e-mean-free-path
      lmfpi = 1.e16*(tif/ev)**2/den(1,1)  # Approx. Coulomb i-mean-free-path
c ... Now check ion pressure gradient scale-lengths
      do nz = 1, 1   # previously nisp - only do majority species now
         if (zi(nz) .gt. 1.e-10) then  # omits neutrals species
            ltmax = min( ltmax, abs( 0.5*(pri(ix,iy,nz)+pri(ix1,iy,nz))/
     .                           (rrv(ix,iy)*gpix(ix,iy,nz) + cutlo) ) )
         endif
      enddo
      flxlimf = 1 / (1 + fricflf*((lmfpe+lmfpi)/ltmax)**2)

c ... Store electron density and temperature, gradients of electron
c     pressure and temperature, (hydrogenic) gas density, ionization rate
c     of gas, and a few other (physically meaningless) quantities for which
c     there are array locations.
      amu(1) = 5.45e-4
      gradp(1,1) = flxlimf*rrv(ix,iy) * gpex(ix,iy)
      gradt(1,1) = flxlimf*rrv(ix,iy) * gtex(ix,iy)
      den(2,0) = 0.5 * (ng(ix,iy,1) + ng(ix1,iy,1))
      nuion(2,0) = den(1,1) * (rsa(tempa(1), den(1,1), 0., 0)
     .                         + sigvi_floor)
      den(1,0) = 0.
      qneut(1) = 0.
      uneut(1) = up(ix,iy,1)  # netural velocity. Cant be zero?; was umass
      nuion(1,0) = 0.
      nuion(1,1) = 0.
      nurec(1,1) = 0.

c ... Loop over isotopes.
      ifld = 1
      do misa = 2, misotope
         amu(misa) = minu(ifld)   # Store mass of this isotope

c ... Store ionization rate of neutral if this isotope is an impurity.
         if (natomic(misa) .gt. 1) then
            den(misa,0) = 0.              # impurity neutral density
            if (ismctab .eq. 1) then
               call imprates(tempa(1), 0, natomic(misa),
     .                       nuion(misa,0), rdum, rdum) 
               nuion(misa,0) = den(1,1)*nuion(misa,0) #sigv-->ne*sigv
            elseif (ismctab .eq. 2) then
               call mcrates(den(1,1), tempa(1), 0., 0, natomic(misa),
     .                      znucl(ifld), nuion(misa,0), rdum, rdum)
               nuion(misa,0) = den(1,1)*nuion(misa,0) #sigv-->ne*sigv
            endif
         endif

c ...    Loop over charged states, storing ni and parallel gradients
c        of pressure and Ti at poloidal density-cell faces.
         do nz = 1, natomic(misa)
            den(misa,nz) = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            gradp(misa,nz) = flxlimf*rrv(ix,iy) * gpix(ix,iy,ifld)
            tempa(misa) = 0.5 * (ti(ix,iy) + ti(ix1,iy))
            gradt(misa,nz) = flxlimf*rrv(ix,iy) * gtix(ix,iy) 
c .......   Get ionization and recombination rates.
c           Note that nuion has no meaning for the fully-stripped state,
c           but space is available to store zero returned by imprates.
            if (isupimpap .eq. 0) then  # omit atomic physics coupling in up
               nuion(misa,nz) = 0.
               nurec(misa,nz) = 0.
            elseif (natomic(misa) .eq. 1) then   # hydrogenic isotope
               nuion(misa,1) = 0.
               nurec(misa,1) = den(1,1) *
     .                         rra(tempa(1), den(1,1), 0., 1)
            else                             # impurity isotope
               if (ismctab .eq. 1) then
                  call imprates(tempa(1), nz, natomic(misa),
     .                          nuion(misa,nz), nurec(misa,nz), rdum)
                  nuion(misa,nz) = den(1,1)*nuion(misa,nz) #sigv-->ne*sigv
                  nurec(misa,nz) = den(1,1)*nurec(misa,nz)
               elseif (ismctab .eq. 2) then
                  call mcrates(den(1,1), tempa(1), 0., nz, natomic(misa),
     .                znucl(ifld), nuion(misa,nz), nurec(misa,nz), rdum)
                  nuion(misa,nz) = den(1,1)*nuion(misa,nz) #sigv-->ne*sigv
                  nurec(misa,nz) = den(1,1)*nurec(misa,nz)
               endif
            endif
            ifld = ifld + 1
            if (natomic(misa).eq.1.and.isupgon(ifld-1).eq.1) ifld=ifld+1
         enddo    # end of loop over charge states with index nz
      enddo    # end of loop over isotopes with index misa

c ... Set up other inputs for fmombal, including flow velocities.
 50   epar = flxlimf*rrv(ix,iy) * ex(ix,iy)
      umassni = 0.
      massni = 0.
      do ifld = 1, nusp
        if(zi(ifld) .gt. 1.e-20) then
          umassni = umassni + mi(ifld)*ni(ix,iy,ifld)*up(ix,iy,ifld)
          massni = massni + mi(ifld)*ni(ix,iy,ifld)
        endif
      enddo
      umass = umassni/massni
      do misa = 2, misotope
         qneut(misa) = 0.
         uneut(misa) = up(ix,iy,1)  # use hydr ion as default (orig)
      enddo
      ldir = 2
      dloglam = loglambda(ix,iy)

c ... Call Steve Hirshman reduced-ion momentum-balance routine.
      call fmombal(amu,den,dloglam,epar,friction,gradp,gradt,
     >	 nuion,nurec,qcond,qneut,ucond,uneut,umass,
     >   parcurrent,tempa,natomic,misotope,nchstate,ldir,friccomp)

c ... Distribute results into arrays used in pandf.  Note that we do
c     nothing with qcond.
      fqp(ix,iy) = cfparcur*parcurrent * rrv(ix,iy)*sx(ix,iy)
      frice(ix,iy) = friction(1,1)
      upe(ix,iy) = ucond(1,1)
      ifld = 0
      do misa = 2, misotope
         do nz = 1, natomic(misa)
 60         ifld = ifld + 1
            if (ziin(ifld) .lt. 1.e-10) goto 60    #omits neutral species
            frici(ix,iy,ifld) = friction(misa,nz)
            upi(ix,iy,ifld) = ucond(misa,nz)
            upifmb(ix,iy,ifld) = ucond(misa,nz)    #diagnostic only
         enddo
      enddo

      return
      end
c****** end of subroutine mombal ************
c-----------------------------------------------------------------------
      subroutine mombalni (ix,ix1,iy)
c ... Use force balance for the impurity momentum equation, neglecting
c     inertia, viscosity, and atomic-physics coupling ala Knoll, Campbell.
c     We also use Keilhacker, et al., Nucl. Fusion., Vol. 31, 537 (1991)
c     which differs somewhat from Campbell: alfi gets divided by zeffv
c     but we retain the extra term 0.6*... in betai from Campbell.
c     Note that results are computed at a poloidal density-cell face.
c     The parallel flow velocities are returned in arrays up and upe.

      implicit none

c ... Input arguments:
      integer ix    # poloidal index of density cell to left of face
      integer ix1   # poloidal index of density cell to right of face
      integer iy    # radial index of density cell

c ... Common blocks:
      Use(Dim)          # nx,ny,nisp,nhsp,nusp
      Use(Selec)        # isupgon
      Use(Comgeo)       # sx,rrv,vol
      Use(UEpar)        # lnlam,isupgon,isofric,is_z0_imp_const,z0_imp_const
      Use(Cfric)        # frice,frici,cfgte,cfgti,cftaud,alfe,betai
      Use(Compla)       # ni,up,upe,upi,te,ti,ng,ne,zeff,netap
      Use(Gradients)    # gpix,gtix,ex,gtex,gpex
      Use(Comflo)       # fqp
      Use(UEint)        # minu,ziin
      Use(Reduced_ion_interface) # misotope,nchstate,natomic,dztot,
                                 # amu,tempa,qneut,uneut,den,gradp,
                                 # gradt,friction,nuion,nurec,qcond,ucond
      Use(Imprad)       # ismctab
      Use(Phyvar)       # ev,qe
      Use(Comtra)       # fricflf,fupe_cur
      Use(Share)        # cutlo
      Use(Coefeq)       # cfnetap
      Use(Volsrc)       # volmsor
      Use(Npes_mpi)
      Use(Conduc)       # pondomfpari_use


c ... Local variables:
      integer ifld, misa, nz
      integer ldir
      real rdum, tdum, dloglam, epar, parcurrent, umass, lmfpe, lmfpi,
     .     zeffv, z0, taud, taudeff, ltmax, tif, flxlimf,
     .     dzj, dzz2tot


c ... Determine a flux-limit factor for all input terms by finding the min
c ... scale length. First consider the electrons
      den(1,1) = 0.5 * (ne(ix,iy) + ne(ix1,iy))
      tempa(1) = 0.5 * (te(ix,iy) + te(ix1,iy))
      tif = 0.5 * (ti(ix,iy) + ti(ix1,iy))
      ltmax = min( abs(tempa(1)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .             abs(tif/(rrv(ix,iy)*gtix(ix,iy) + cutlo)),
     .             abs(den(1,1)*tempa(1)/
     .                          (rrv(ix,iy)*gpex(ix,iy) + cutlo)) )
      lmfpe = 1.e16*((tempa(1)/ev)**2/den(1,1)) #Approx Coulomb e-mean-free-path
      lmfpi = 1.e16*((tif/ev)**2/den(1,1))  # Approx. Coulomb i-mean-free-path
c ... Now check ion pressure gradient scale-lengths
      do nz = 1, 1      # previously nisp - only do majority species now
         if (zi(nz) .gt. 1.e-10) then  # omits neutrals species
            ltmax = min( ltmax, abs( 0.5*(pri(ix,iy,nz)+pri(ix1,iy,nz))/
     .                           (rrv(ix,iy)*gpix(ix,iy,nz) + cutlo) ) )
         endif
      enddo
      flxlimf = 1 / (1 + fricflf*((lmfpe+lmfpi)/ltmax)**2)

c ... Store electron density and temperature, gradients of electron
c     pressure and temperature, (hydrogenic) gas density, ionization rate
c     of gas, and a few other (physically meaningless) quantities for which
c     there are array locations.
      gradp(1,1) = rrv(ix,iy) * gpex(ix,iy)
      gradt(1,1) = rrv(ix,iy) * gtex(ix,iy)
      frice(ix,iy) = -0.71*flxlimf*den(1,1)*gradt(1,1) +
     .                cfnetap*qe*netap(ix,iy)*fqp(ix,iy)/sx(ix,iy)

c ... Loop over charge states to get total impurity density
      ifld = nhsp
      dztot = 0.
      dzz2tot = 0.
      do misa = 3, misotope
         do nz = 1, natomic(misa)
            ifld = ifld + 1
            dzj = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            dztot = dztot + dzj
            dzz2tot = dzz2tot + dzj*zi(ifld)**2
         enddo
      enddo

c ... Set the hydrogen values based on electrons (trace limit)
ccc      frici(ix,iy,1) = - frice(ix,iy)   # needed for hydrogen
ccc   For arbitrary impurity concentration use the following frici:
      frici(ix,iy,1) = - frice(ix,iy) *
     .                   ni(ix,iy,1)*zi(1)**2/(ne(ix,iy)*zeff(ix,iy))
      upe(ix,iy) = up(ix,iy,1) - fupe_cur*fqp(ix,iy)/( sx(ix,iy)*qe*
     .                         rrv(ix,iy)*0.5*(ne(ix,iy)+ne(ix1,iy)) )
      upi(ix,iy,1) = up(ix,iy,1)
      den(2,1) = 0.5 * (ni(ix,iy,1) + ni(ix1,iy,1))
           
c ... Loop over isotopes for friction coefficients
      ifld = 1
      do misa = 3, misotope      # only executed if impurities are present

c ... Loop over charged states, storing ni and parallel gradients
c ... of pressure and Ti at poloidal density-cell faces.
         do nz = 1, natomic(misa)  # note: hydrogen friction set above
            ifld = ifld + 1
            if (ziin(ifld) .lt. 1.e-10) ifld = ifld+1 #skip gas index
            den(misa,nz) = 0.5 * (ni(ix,iy,ifld) + ni(ix1,iy,ifld))
            gradp(misa,nz) = rrv(ix,iy) * gpix(ix,iy,ifld) -
     .                                      pondomfpari_use(ix,iy,ifld)
            tempa(misa) = 0.5 * (ti(ix,iy) + ti(ix1,iy))
            gradt(misa,nz) = rrv(ix,iy) * gtix(ix,iy) 
            zeffv = 0.5*(zeff(ix,iy)+zeff(ix1,iy))
            if (is_z0_imp_const == 0) then
              z0 = den(1,1)*zeffv/den(2,1) - 1.
            else # fixed user input
              z0 = z0_imp_const
            endif
            if (isbetaicalc(ifld) == 1) then
              betai(ifld)=cfgti*1.56*zi(ifld)**2*(1+1.414*z0)*(1+.52*z0)/
     .                       ( (1+2.65*z0)*(1+.285*z0)*( z0 + 
     .                         sqrt( 0.5*(mi(1)+mi(ifld))/mi(ifld)) ) )
     .                       + 0.6*(zi(ifld)**2*den(misa,nz)/dzz2tot - 1.)
            endif
            if (isalfecalc(ifld) == 1) then
              alfe(ifld) = cfgte*2.2*zi(ifld)**2*(1+.52*zeffv) /
     .                      ( (1+2.65*zeffv)*(1+.285*zeffv)*zeffv )
            endif
c... NOTE:next coefficient 12*pi*sqrt(pi/2)*epsilon**2/e**4 = 5.624e54 in mks
            taud =cftaud*5.624e54*sqrt(mi(1))*mi(ifld)*tempa(misa)**1.5 /
     .             ( loglambda(ix,iy)*den(misa,nz)*zi(ifld)**2*
     .                                                 (mi(1)+mi(ifld)) )
            taudeff = flxlimf*taud*den(misa,nz)*(1+2.65*z0)*(1+.285*z0) /
     .                         ( den(1,1)*(1+.24*z0)*(1+.93*z0) )
            upi(ix,iy,ifld) = up(ix,iy,1) + (taudeff/mi(1)) * (
     .                         - gradp(misa,nz)/den(misa,nz)
     .                         + alfe(ifld)*gradt(1,1) 
     .                         + betai(ifld)*gradt(misa,nz) 
     .                         + qe*zi(ifld)*rrv(ix,iy)*ex(ix,iy) 
     .                         + volmsor(ix,iy,ifld)/
     .                                       (den(misa,nz)*vol(ix,iy)) )
c ...       For force balance, frici just balances E-field and pressure
c ...       No flxlimf for 1st option; it only enhances (1/taudeff)
            if (nusp-isupgon(1) .eq. 1) then  #only frici(,,1) used here
              frici(ix,iy,ifld) =-qe*zi(ifld)*den(misa,nz)*
     .                           rrv(ix,iy)*ex(ix,iy) + gradp(misa,nz)  
            else # multi ion mom eqns; drag calc elsewhere (w0) if isofric=1
              frici(ix,iy,ifld) =flxlimf*den(misa,nz)*( 
     .                                          alfe(ifld)*gradt(1,1) + 
     .                                     betai(ifld)*gradt(misa,nz) + 
     .                                               (1-isofric)*mi(1)*
     .                           (up(ix,iy,1)-up(ix,iy,ifld))/ taudeff )
            endif
         enddo
      enddo

      return
      end
c****** end of subroutine mombalni ************

      SUBROUTINE calc_plasma_diffusivities
      IMPLICIT NONE
      Use(Dim)
      Use(Comtra)
      Use(Selec)
      Use(Conduc)
      Use(Compla)
      Use(Phyvar)
      Use(Bfield)
      Use(Coefeq)
      Use(Aux)
      Use(Xpoint_indices)
      Use(RZ_grid_info)
      Use(Share)
      integer ifld, iy, iyp1, ix, ix1
      real bpfac, bpolmin, bscalf, t0, t1
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)
c ... Calculate the Bohm diffusion rates (units are m**2/s)
      do ifld = 1, nisp
       if (facbni+facbup+facbee+facbei>0 .and. isbohmcalc>0) then
         do iy = j1omp1, j6omp
            iyp1 = min(ny+1, iy+1)
            do ix = i1omp, i6omp
               ix1 = ixp1(ix,iy)
               kybohm(ix,iy) = (te(ix,iy)+te(ix,iyp1)) /
     .                        (16*ev*(btot(ix,iy)+btot(ix,iyp1)))
               kxbohm(ix,iy) = (te(ix,iy)+te(ix1,iy)) /
     .                        (16*ev*(btot(ix,iy)+btot(ix1,iy)))
               dif_use(ix,iy,ifld)  = facbni*kybohm(ix,iy)
               dif2_use(ix,iy,ifld) = facbni2*kxbohm(ix,iy)
               tray_use(ix,iy,ifld)  = facbup*kybohm(ix,iy)
               kye_use(ix,iy)  = facbee*kybohm(ix,iy)
               kyi_use(ix,iy)  = facbei*kybohm(ix,iy)
	       dutm_use(ix,iy,ifld) = kybohm(ix,iy)
            enddo
         enddo
         if (isbohmcalc.eq.2) then  # calc. recip. average with const D
           fcdif = 0.   # used to zero constant diff. if recip. ave used
           do iy = j1omp1, j6omp
             do ix = i1omp, i6omp
               dif_use(ix,iy,ifld)  = 0.5*ave(difni(ifld),  dif_use(ix,iy,ifld))
               dif2_use(ix,iy,ifld) = 0.5*ave(difni2(ifld), dif2_use(ix,iy,ifld))
               tray_use(ix,iy,ifld)  = 0.5*ave(travis(ifld), tray_use(ix,iy,ifld))
               kye_use(ix,iy)  = 0.5*ave(kye, kye_use(ix,iy))
               kyi_use(ix,iy)  = 0.5*ave(kyi, kyi_use(ix,iy))
               dutm_use(ix,iy,ifld)  = 0.5*ave(difutm(1), dutm_use(ix,iy,ifld))
             enddo
           enddo
         endif
       endif

c ... If isbohmcalc=3, then give (B0/B)**inbdif profile to diff
       if (isbohmcalc==3) then  # use inbtdif, inbpdif for btot, bpol scaling
         bpolmin = bpol(ixpt2(1)-ixbpmin,iysptrx,0)
         do iy = j1omp1, j6omp
            do ix = i1omp, i6omp
              ix1 = ixp1(ix,iy)
	      bscalf=((.5*(btot(ixmp,iysptrx)/btot(ix,iy)) +
     .               .5*(btot(ixmp,iysptrx)/btot(ix1,iy)))**inbtdif)
     .       * ((bpol(ixmp,iysptrx,3)+bpol(ixmp,iysptrx,4))/
     .          (bpol(ix,iy,3)+bpol(ix,iy,4)+bpolmin))**inbpdif
	      dif_use(ix,iy,ifld)  = difniv(iy,ifld)*bscalf
	      difp_use(ix,iy,ifld) = difprv(iy,ifld)*bscalf          
              dif2_use(ix,iy,ifld) = difniv2(iy,ifld)*bscalf
              tray_use(ix,iy,ifld)  = travisv(iy,ifld)*bscalf
              kye_use(ix,iy)  = kyev(iy)*bscalf
              kyi_use(ix,iy)  = kyiv(iy)*bscalf
              dutm_use(ix,iy,ifld) = difutmv(iy,ifld)*bscalf
	      vy_use(ix,iy,ifld) = vconyv(iy,ifld)*bscalf
            enddo
         enddo
       endif
      enddo  # loop over species lfld

c ,,, Add diffusion propto betap**iexpbp and (B0/B)**inbdif (as for isbohmcalc=3)
      if (isdifbetap == 1) then
       do ifld = 1, nisp
         if(zi(ifld) > 0.) then
           bpolmin = bpol(ixpt2(1)-ixbpmin,iysptrx,0)
           do iy = j1omp1, j6omp
             do ix = i1omp, i6omp
               ix1 = ixp1(ix,iy)
               betap(ix,iy) = 8.*pi*1.e-7*pr(ix,iy)/bpol(ix,iy,0)**2
               bpfac = betap(ix,iy)**iexpbp
 	       bscalf = ((.5*(btot(ixmp,iysptrx)/btot(ix,iy)) +
     .                   .5*(btot(ixmp,iysptrx)/btot(ix1,iy)))**inbtdif)
     .                  * ((bpol(ixmp,iysptrx,3)+bpol(ixmp,iysptrx,4))/
     .                   (bpol(ix,iy,3)+bpol(ix,iy,4)+bpolmin))**inbpdif
	       dif_use(ix,iy,ifld)  = difniv(iy,ifld)*bscalf +
     .                                                    dfacbp*bpfac
	       difp_use(ix,iy,ifld) = difprv(iy,ifld)*bscalf          
               dif2_use(ix,iy,ifld) = difniv2(iy,ifld)*bscalf + 
     .                                                    dfacbp*bpfac
               tray_use(ix,iy,ifld)  = travisv(iy,ifld)*bscalf + 
     .                                                   trfacbp*bpfac
               trax_use(ix,iy,ifld) = trfacbp*bpfac
               kye_use(ix,iy)  = kyev(iy)*bscalf +  kefacbp*bpfac
               kyi_use(ix,iy)  = kyiv(iy)*bscalf + kifacbp*bpfac
               kxe_use(ix,iy) = kefacbp*bpfac
               kxi_use(ix,iy) = kifacbp*bpfac
               dutm_use(ix,iy,ifld) = difutmv(iy,ifld)*bscalf
	       vy_use(ix,iy,ifld) = vconyv(iy,ifld)*bscalf
             enddo
           enddo
         endif
       enddo
      endif   # test on isdifbetap


      END SUBROUTINE calc_plasma_diffusivities


      SUBROUTINE calc_elec_velocities
      IMPLICIT NONE
      Use(Selec)
      Use(Compla)
      Use(Imprad)
      Use(Dim)
      Use(Comtra)
      Use(Comflo)
      Use(Comgeo)
      Use(Phyvar)
      Use(Coefeq)
      Use(Bfield)
      Use(UEpar)
      Use(Bcond)
      integer iy, ix, ifld, ix1
      real afqp
************************************************************************
*     Calculate the electron velocities, vex, upe, ve2, vey
************************************************************************
       do iy = j1omp1, j6omp
         do ix = i1momp, i6omp
            vex(ix,iy) = 0.
            vey(ix,iy) = 0.
          end do
        end do

      if (isimpon.ne.5) then    # have upe from mombal

      do iy = j1omp1, j6omp    #iys1, iyf6
         do ix = i1momp, i6omp
            upe(ix,iy) = 0.
         enddo
      enddo

      do ifld = 1, nfsp
         do iy = j1omp1, j6omp    #iys1, iyf6
	    do ix = i1momp, i6omp
               ix1 = ixp1(ix,iy)
	       upe(ix,iy) = upe(ix,iy) + upi(ix,iy,ifld)*zi(ifld)*0.5*
     .                      ( ni(ix,iy,ifld)+ni(ix1,iy,ifld) )
            enddo
         enddo
        end do
      afqp = 1.
      if (isimpon.eq.6 .or. isimpon.eq.7) afqp = fupe_cur  #allows gradual fix for old cases
      do iy = j1omp1, j6omp    #iys1, iyf6
         do ix = i1momp, i6omp
            ix1 = ixp1(ix,iy)
	    upe(ix,iy) = (upe(ix,iy) -afqp*fqp(ix,iy)/
     .                               (rrv(ix,iy)*sx(ix,iy)*qe))/
     .                             (0.5*( ne(ix,iy)+ne(ix1,iy) ))
         enddo
      enddo

      end if

      do iy = j1omp1, j6omp   # ExB same all species;if cf2dd=1, no imp yet
	    do ix = i1momp, i6omp
            ix1 = ixp1(ix,iy)
            vex(ix,iy) = upe(ix,iy)*rrv(ix,iy) + 
     .                   (cf2ef*v2ce(ix,iy,1) + cf2bf*ve2cb(ix,iy) + 
     .                         cf2dd*bfacxrozh(ix,iy)*ve2cd(ix,iy,1) ) *
     .                           0.5*(rbfbt(ix,iy) + rbfbt(ix1,iy)) -
     .                                               vytan(ix,iy,1) 
     
        end do
       end do

      do ifld = 1, nfsp
	 do iy = j1omp1, j5omp
	    do  ix = i1momp, i6omp   # grad_B will be ok as next fqy is subtr.
	       vey(ix,iy) = vey(ix,iy) + vy(ix,iy,ifld)*zi(ifld)*0.5*
     .                      ( niy0(ix,iy,ifld)+niy1(ix,iy,ifld) )
          end do
         end do
        end do

        do iy = j1omp1, j5omp
	     do ix = i1momp, i6omp
	      vey(ix,iy) = (vey(ix,iy)-cfjve*fqy(ix,iy)/(sy(ix,iy)*qe))/
     .                    (0.5*( ney0(ix,iy)+ney1(ix,iy) ))
         end do
        end do
	 
c ... if isnewpot=0, vey(,0) needs to be redone since fqy(,0)=0
      if (isnewpot==1) then
        do ix = i1momp, i6omp  # ExB vyce same all species
          vey(ix,0) = cfybf*veycb(ix,0) + vydd(ix,0,1) +
     .                cfyef*vyce(ix,0,1)
        enddo
      endif

c ... If isybdrywd = 1, make vey diffusive, just like vy
      if (isybdrywd == 1) then  #make vy diffusive in wall cells
        do ix = i1momp, i6omp
          if (matwalli(ix) > 0) vey(ix,0)  = vydd(ix,0,1)
          if (matwallo(ix) > 0) vey(ix,ny) = vydd(ix,ny,1)
        enddo
      endif
      


      END SUBROUTINE calc_elec_velocities



      SUBROUTINE calc_plasma_viscosities
      IMPLICIT NONE
      Use(Dim)
      Use(UEpar)
      Use(Compla)
      Use(Selec)
      Use(Phyvar)
      Use(Comtra)
      Use(Coefeq)
      Use(Comgeo)
      Use(Conduc)
      Use(Share)
      Use(Wkspace)
      integer ifld, iy, iyp1, ix, ix1, ix2, ix3, jfld
      real vtn, lmfpn, lmfppar, lmfpperp, rrfac, tgupyface, nmxface, 
     .  ngupyface, n1upyface, ng2upyface, tv, lambda, tv2, 
     .  epstmp, visxtmp, rt2nus, t0, a
*****************************************************************
*  Other physics coefficients. (old PHYVIS)
*****************************************************************

* -- loop over species number --

      do ifld = 1, nfsp

c
c     neutral viscosity for isupgon=1
c
c  Logic here is specialized; only ifld=2 (igsp=1) has viscosity calc.
c  If more neutral species have full parallel mom eqn, need to redo loops
         if(isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0)then
            do iy = j1,j6
               iyp1 = min(iy+1,ny+1)
               do ix = i1,i6
c
                  ix1 = ixm1(ix,iy)
                  vtn = sqrt(max(tg(ix,iy,1),tgmin*ev)/mi(ifld))
 		  qfl = flalfvgxa(ix)*nm(ix,iy,ifld)*vtn**2
                  if(isvisxn_old == 1) then
                    lmfpn = 1./(sigcx * 
     .                          (ni(ix,iy,1) + rnn2cx*ni(ix,iy,ifld)))
                  elseif(isvisxn_old==0 .and. ishymol==0) then
                    lmfppar = vtn/(kelhihg*ni(ix,iy,1) +
     .                                         kelhghg*ni(ix,iy,ifld))
                    lmfpperp = vtn/( vtn*sigcx*ni(ix,iy,1) + 
     .                      kelhihg*ni(ix,iy,1)+kelhghg*ni(ix,iy,ifld) )
                    rrfac = rr(ix,iy)*rr(ix,iy)
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  else   # (isvisxn_old=0 .and. ishymol=1) then #with mols
                    lmfppar = vtn/(kelhihg*ni(ix,iy,1) +
     .                     kelhghg*ni(ix,iy,ifld) + kelhmhg*ng(ix,iy,2))
                    lmfpperp = vtn/( vtn*sigcx*ni(ix,iy,1) + 
     .                     kelhihg*ni(ix,iy,1) +kelhghg*ni(ix,iy,ifld) +
     .                     kelhmhg*ng(ix,iy,2) )
                    rrfac = rr(ix,iy)*rr(ix,iy)
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  endif
                  csh = lmfpn*nm(ix,iy,ifld)*vtn*
     .                                      lgvmax/(lgvmax + lmfpn)  
                  if (isgxvon .eq. 0) then 
                    qsh = csh * (up(ix1,iy,ifld)-up(ix,iy,ifld))
     .                                         *gx(ix,iy)
                  elseif (isgxvon .eq. 1) then
                    qsh = csh * (up(ix1,iy,ifld)-up(ix,iy,ifld))
     .                *2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
                  endif
                  visx(ix,iy,ifld)= cfvisxn*csh/ 
     .               (1 + (abs(qsh/(qfl+cutlo))**flgamvg))**(1./flgamvg)
     .               + cfanomvisxg*travis(ifld)*nm(ix,iy,ifld)

c    Now do y-direction; use ni on up y-face
                  ix2 = ixp1(ix,iy)
                  ix3 = ixp1(ix,iyp1)
                  tgupyface = 0.25*( tg(ix,iy,1)+
     .                         tg(ix,iyp1,1)+tg(ix2,iy,1)+
     .                         tg(ix3,iyp1,1) )
                  vtn = sqrt(max(tgupyface,tgmin*ev)/mi(ifld))
                  nmxface = 0.5*(nm(ix,iy,ifld)+nm(ix2,iy,ifld))
                  ngupyface = 0.25*( ni(ix,iy,ifld)+
     .                         ni(ix,iyp1,ifld)+ni(ix2,iy,ifld)+
     .                         ni(ix3,iyp1,ifld) )
                  n1upyface = 0.25*( ni(ix,iy,1)+
     .                         ni(ix,iyp1,1)+ni(ix2,iy,1)+
     .                         ni(ix3,iyp1,1) )
                  if(ishymol == 0) then
                    lmfppar = vtn/(kelhihg*n1upyface +
     .                                         kelhghg*ngupyface)
                    lmfpperp = vtn/( vtn*sigcx*n1upyface + 
     .                      kelhihg*n1upyface+kelhghg*ngupyface )
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  else  # ishymol=1
                    ng2upyface = 0.25*( ng(ix,iy,2)+
     .                         ng(ix,iyp1,2)+ng(ix2,iy,2)+
     .                         ng(ix3,iyp1,2) )
                    lmfppar = vtn/(kelhihg*n1upyface +
     .                     kelhghg*n1upyface + kelhmhg*ng2upyface)
                    lmfpperp = vtn/( vtn*sigcx*n1upyface + 
     .                     kelhihg*n1upyface +kelhghg*ngupyface +
     .                     kelhmhg*ng2upyface )
                    lmfpn = lmfppar*rrfac + lmfpperp*(1-rrfac)
                  endif
        
                  csh = lmfpn*ngupyface*mi(ifld)*vtn*
     .                                      lgvmax/(lgvmax + lmfpn)  

 		  qfl = flalfvgya(iy)*ngupyface*mi(ifld)*vtn**2
                  qsh = csh * (up(ix,iy,ifld)-up(ix,iyp1,ifld)) * 
     .                                                        gyf(ix,iy)
                  visy(ix,iy,ifld)= cfvisyn*csh / 
     .               (1 + (abs(qsh/(qfl+cutlo))**flgamvg))**(1./flgamvg)
     .               + cfanomvisyg*travis(ifld)*nmxface
c
            end do
        end do
         endif
c
c
       if(zi(ifld) > 1.e-20) then
         do iy = j1, j6
            do ix = i1, i6
               w(ix,iy) = 0.0e0
            end do
        end do

         do jfld = 1, nisp
            tv = zi(jfld)**2 / sqrt((mi(ifld)+mi(jfld))/(2*mp))
            do iy = j1, j6
               do ix = i1, i6
                  w(ix,iy) = w(ix,iy) + tv*ni(ix,iy,jfld)
            end do
          end do
        end do

         do iy = j1, j6
            do ix = i1, i6
	       ctaui(ix,iy,ifld) = 2.1e13/(loglambda(ix,iy)*zi(ifld)**2) # mass fac?
               tv2 = ctaui(ix,iy,ifld)/(ev*sqrt(ev))
               if (convis .eq. 0) then
                  a = max (ti(ix,iy), temin*ev)
               else
                  a = afix*ev
               endif
	       epstmp = max(epsneo(ix,iy), 1.e-50)
               visxtmp = tv2 * coef * rr(ix,iy) * rr(ix,iy) *
     .                            a*a*sqrt(a) * ni(ix,iy,ifld)/w(ix,iy)
               visx(ix,iy,ifld) = parvis(ifld)*visxtmp +
     .                            trax_use(ix,iy,ifld)*nm(ix,iy,ifld)
               nuii(ix,iy,ifld) = w(ix,iy)/(tv2*a*sqrt(a))
               nuiistar(ix,iy,ifld) = ( lconneo(ix,iy)*nuii(ix,iy,ifld)/
     .                                 (epstmp**1.5*(2*a/mi(ifld))**0.5)
     .                                           + 1.e-50 )
               visxneo(ix,iy,ifld) = visxtmp*
     .                 (1./(1.+epstmp**(-1.5)/nuiistar(ix,iy,ifld)))*
     .                               (1./(1.+1./nuiistar(ix,iy,ifld)))
               rt2nus = 1.414*nuiistar(ix,iy,ifld)
               ktneo(ix,iy,ifld) = (-0.17 + 1.05*rt2nus**.5 + 
     .                     2.7*rt2nus**2*epstmp**3) / ( 1.+
     .                0.7*rt2nus**.5 + rt2nus**2*epstmp**3 )
               alfneo(ix,iy,ifld) = (8./15.)*(ktneo(ix,iy,ifld) - 1.)*
     .                 (1./(1.+epstmp**(-1.5)/nuiistar(ix,iy,ifld)))*
     .                            ( 1./(1.+1./nuiistar(ix,iy,ifld)) )
               k2neo(ix,iy,ifld) =(.66 + 1.88*epstmp**.5 - 1.54*epstmp)/ 
     .                            (1. + 1.03*rt2nus**.5 + 0.31*rt2nus) +
     .             1.17*epstmp**3*rt2nus/(1. + 0.74*epstmp**1.5*rt2nus)
c...  flux limit the viscosity; beware of using visx(0,iy) and 
c...  visx(nx+1,iy) as they are meaningless when flux limited
               ix1 = ixm1(ix,iy)
               t0 = max (ti(ix,iy), temin*ev)
               vtn = sqrt(t0/mi(ifld))
               mfl = flalfv * nm(ix,iy,ifld) * rr(ix,iy) *
     .               vol(ix,iy) * gx(ix,iy) * (t0/mi(ifld)) 
ccc  Distance between veloc. cell centers:
               if (isgxvon .eq. 0) then     # dx(ix)=1/gx(ix)
                 csh = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .                            * gx(ix,iy)
               elseif (isgxvon .eq. 1) then # dx(ix)=.5/gxf(ix-1)+.5/gxf(ix)
                 csh = visx(ix,iy,ifld) * vol(ix,iy) * gx(ix,iy)
     .               * 2*gxf(ix,iy)*gxf(ix1,iy)/(gxf(ix,iy)+gxf(ix1,iy))
               endif
ccc 
               msh = abs( csh*(upi(ix1,iy,ifld) - upi(ix,iy,ifld)) )
               visx(ix,iy,ifld) = visx(ix,iy,ifld)
     .               / (1 + (msh/(mfl+1.e-20*msh))**flgamv )**(1/flgamv)
               visy(ix,iy,ifld)=(fcdif*travis(ifld)+ tray_use(ix,iy,ifld))*
     .                              nm(ix,iy,ifld) +  4*eta1(ix,iy)
            end do
        end do
       endif      # test if zi(ifld) > 1.e-20
        end do    # large loop for ifld = 1, nfsp


      END SUBROUTINE calc_plasma_viscosities


