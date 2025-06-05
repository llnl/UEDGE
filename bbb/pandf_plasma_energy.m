c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"


      SUBROUTINE calc_plasma_energy(xc, yc)
      IMPLICIT NONE
      Use(Selec)
      Use(Locflux)
      Use(Comflo)
      Use(Wkspace)
      Use(Rhsides)
      Use(Compla)
      Use(UEpar)
      Use(Phyvar)
      Use(Dim)
      Use(Xpoint_indices)
      Use(Indices_domain_dcl)
      Use(Comtra)
      Use(Comgeo)
      Use(Conduc)
      Use(Share)
      Use(Gradients)
      Use(Coefeq)
      Use(Poten)
      Use(MCN_sources)
      Use(Bcond)
      Use(Rccoef)
      Use(Bfield)
      Use(Volsrc)
      Use(Noggeo)
      Use(Jacobian_restore)
      Use(Ext_neutrals)
      Use(Imprad)
      Use(Timing)
      Use(ParallelEval)
      integer xc, yc
      integer iy, ix, iy1, ix1, ix2, ix3, ix4, ix5, ix6, jx, jfld, jz, 
     .  igsp, iy_min, iy_max, j2pwr, j5pwr, i2pwr, i5pwr, ifld, nsm1, zn,
     .  znuc, zmax, iixt, iym1, iyp1, iyp2
      real grdnv, fcd, t0, t1, vttn, vttp, isfe, l0, feexflr, feixflr, 
     .  ne_sgvi, dene, rdum, radmc, radz(0:1), erl1, erl2, 
     .  up1cc, upgcc, vycc, v2cc, tsimp, tick, emissbs, radneq,  
     .  argth, fac_rad, radimpmc, wj, cfwj, telim, thetacc, 
     .  dupdx, dupdy, upxavep1, upxave0, upxavem1, upf0, upfm1,
     .  denz(0:1), sv_crumpet, vt1, wallfac, lxtec, qr, lxtic,
     .  ltmax, lmfpe, flxlimf, fniy_recy, temp1, dtdym1, dtdy0,
     .  d2tdy20, d2tdy2p1, d3tdy3, dtdyp1, vt0
      external sv_crumpet, radmc

*****************************************************************
*****************************************************************
*  Here starts the old ENEBAL
*****************************************************************
*  ---------------------------------------------------------------------
*  compute temperature conductances.
*  ---------------------------------------------------------------------
*  -- initialize to 0 --

        floxe = 0.0e0
        floxi = 0.0e0
        floye = 0.0e0
        floyi = 0.0e0
        feiycbo = 0.0e0
        feeycbo = 0.0e0
        w0 = 0.0e0
        w1 = 0.0e0
        wvh = 0.0e0    #wvh is 3D vector (ix,iy,nusp)

*  -- compute conxe and conxi --

*     (The computation of conxe involves a flux limit)

      do iy = j4omp, j8omp
         do ix = i1momp, i5omp
            ix2 = ixp1(ix,iy)
            t0 = max (te(ix,iy), temin*ev)
            t1 = max (te(ix2,iy), temin*ev)
            vt0 = sqrt(t0/me)
            vt1 = sqrt(t1/me)
            wallfac = 1.
            do jx = 1, nxpt
               if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                (ix==ixrb(jx).and.ixmxbcl==1) )
     .              .and. (isplflxl==0) ) wallfac = flalfepl/flalfe
            enddo
            qfl = wallfac*flalfe * sx(ix,iy) * rrv(ix,iy) *
     .          (ne(ix,iy)*vt0*t0 + ne(ix2,iy)*vt1*t1) / 2
            csh = sx(ix,iy) * hcxe(ix,iy) * gxf(ix,iy)
	    lxtec = 0.5*(te(ix,iy)+te(ix2,iy)) /
     .             ( abs(te(ix,iy)-te(ix2,iy))*gxf(ix,iy) + 100.*cutlo )
	    qsh = csh * (te(ix,iy)-te(ix2,iy)) * (1. + lxtec/lxtemax)
            qr = (1-isflxlde)*abs(qsh/qfl)
            conxe(ix,iy) = (1-isflxlde)*csh / (1 + qr)**2 +
     .                isflxlde*csh /(1 + abs(qsh/qfl)**flgam)**(1/flgam)
            floxe(ix,iy) = floxe(ix,iy) + (sign(qr*qr,qsh)/(1 + qr)**2)
     .                  *flalfea(ix) * sx(ix,iy) * ( ne(ix,iy)*
     .                   rr(ix,iy)*vt0 + ne(ix2,iy)*rr(ix2,iy)*vt1 ) / 2
c.... Now do the ions (hcxi is flux-limited previously when it is built)
          if (isflxldi .ne. 2) then    # Else flux limit done on hcxij
            t0 = max (ti(ix,iy), temin*ev)
            t1 = max (ti(ix2,iy), temin*ev)
            vt0 = sqrt(t0/mi(1))
            vt1 = sqrt(t1/mi(1))
            wallfac = 1.
            do jx = 1, nxpt
               if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                (ix==ixrb(jx).and.ixmxbcl==1) )
     .              .and. (isplflxl==0) ) wallfac = flalfipl/flalfi
            enddo
            qfl = wallfac*flalfia(ix) * sx(ix,iy) * rrv(ix,iy) *
     .                    (ne(ix,iy)*vt0*t0 + ne(ix2,iy)*vt1*t1) / 2
            csh = sx(ix,iy) * hcxi(ix,iy) * gxf(ix,iy)
            lxtic = 0.5*(ti(ix,iy)+ti(ix2,iy)) /
     .             ( abs(ti(ix,iy)-ti(ix2,iy))*gxf(ix,iy) + 100.*cutlo)
	    qsh = csh * (ti(ix,iy)-ti(ix2,iy)) * (1. + lxtic/lxtimax)
            qipar(ix,iy,1) = qsh/(rrv(ix,iy)*sx(ix,iy))
            qr = (1-isflxldi)*abs(qsh/qfl)
            conxi(ix,iy) = (1-isflxldi)*csh / (1 + qr)**2 +
     .                isflxldi*csh / (1 + abs(qsh/qfl)**flgam)**(1/flgam)
            floxi(ix,iy) = floxi(ix,iy) + (sign(qr*qr,qsh)/(1 + qr)**2)
     .                  *flalfia(ix) * sx(ix,iy) *( ne(ix,iy)*rr(ix,iy)*
     .                   vt0 + ne(ix2,iy)*rr(ix2,iy)*vt1 ) / 2 
          else
            conxi(ix,iy) = sx(ix,iy) * hcxi(ix,iy) * gxf(ix,iy)
          endif
        end do
         conxe(nx+1,iy) = 0
         conxi(nx+1,iy) = 0
        end do

*  -- compute conye and conyi --

      do iy = j1omp1, j5omp
         do ix = i4omp, i8omp
            conye(ix,iy) = sy(ix,iy) * hcye(ix,iy) / dynog(ix,iy)
            conyi(ix,iy) = sy(ix,iy) * hcyi(ix,iy) / dynog(ix,iy)
          end do
      end do

      do ix = i1omp, i6omp
         conye(ix,ny+1) = 0.0e0
         conyi(ix,ny+1) = 0.0e0
      end do

*  ---------------------------------------------------------------------
*  compute the strength of convection for te and ti.
*  ---------------------------------------------------------------------
*  includes a correction to account for e-velocity .ne. i-velocity
*  (ccn term), and also a conduction contribution to electron heat flux
*  due to friction (cthe term).  (M.E. Rensink, 07/20/90)
*  ---------------------------------------------------------------------
*     floxe, floxi  contain the cross-derivative terms now
*                        JLM      5/3/90
*  ---------------------------------------------------------------------

      do iy = j4omp, j8omp
         do ix = i1momp, i5omp  
            ix1 = ixp1(ix,iy)
            ltmax = min( abs(te(ix,iy)/(rrv(ix,iy)*gtex(ix,iy) + cutlo)),
     .                   lcone(ix,iy) )
            lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)
            flxlimf = flalftf*ltmax/(flalftf*ltmax + lmfpe)
            floxe(ix,iy) = floxe(ix,iy) + cfcvte*1.25*
     .                  (ne(ix,iy)+ne(ix1,iy))*vex(ix,iy)*sx(ix,iy)
     .                   - cthe*flxlimf*cfjhf*fqp(ix,iy)/ev
         end do
         floxe(nx+1,iy) = 0.0e0
        end do

c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add 
      do ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then  #neutrals
            do iy = j4omp, j8omp
               do ix = i1omp, i5omp
                  floxi(ix,iy) = floxi(ix,iy) +
     .                 cftiexclg*cfcvti*2.5*cfneut*cfneutsor_ei*fnix(ix,iy,ifld) 
               end do # next correct for incoming neut pwr = 0
               do jx = 1, nxpt  #if at plate, sub (1-cfloxiplt)*neut-contrib
                 if(ixmnbcl==1) then  #real plate-need for parallel UEDGE
                   iixt = ixlb(jx) #left plate
                   if(fnix(iixt,iy,ifld) > 0.) then
                     floxi(iixt,iy) = floxi(iixt,iy) - (1.-cfloxiplt)*
     .                 cftiexclg*cfcvti*2.5*cfneut*cfneutsor_ei*fnix(iixt,iy,ifld)
                   endif
                 endif
                 if(ixmxbcl==1) then #real plate-need for parallel UEDGE
                   iixt = ixrb(jx) # right plate
                   if(fnix(iixt,iy,ifld) < 0.) then
                     floxi(iixt,iy) = floxi(iixt,iy) - (1.-cfloxiplt)*
     .                 cftiexclg*cfcvti*2.5*cfneut*cfneutsor_ei*fnix(iixt,iy,ifld)
                   endif
                   floxi(ixrb(jx)+1,iy) = 0.0e0  #cosmetic
                 endif
               enddo
            end do
         else  #ions
            do iy = j4omp, j8omp
               do ix = i1momp, i5omp
                  floxi(ix,iy) = floxi(ix,iy) +
     .                           cfcvti*2.5*fnix(ix,iy,ifld)
                end do 
               floxi(nx+1,iy) = 0.0e0
            end do
         endif
        end do

*  -- compute floye and floyi --

      do iy = j1omp1, j5omp    # note: cfloye usually = 2.5 or 1.5 (ExB turb)
         do ix = i4omp, i8omp
            floye(ix,iy) = floye(ix,iy) + (cfloye/2.)*
     .                    (ney0(ix,iy)+ney1(ix,iy))*vey(ix,iy)*sy(ix,iy)
     .                + (vyte_use(ix,iy)+vyte_cft(ix,iy))*0.5*sy(ix,iy)*
     .                     (ney0(ix,iy)+ney1(ix,iy))
        end do
      end do

      do ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then
            do iy = j1omp1, j5omp
               do ix = i4omp, i8omp
                  floyi(ix,iy) = floyi(ix,iy)
     .                 + cftiexclg*cfneut * cfneutsor_ei * 2.5 * fniy(ix,iy,ifld)
               enddo
            enddo
c ...       Make correction at walls to prevent recyc neutrals injecting pwr
            do ix = i4omp, i8omp
              if (matwallo(ix) > 0 .and. recycwot(ix,1)>0.) then
                fniy_recy = max(recycwot(ix,1)*fac2sp*fniy(ix,ny,1), 0.)
                floyi(ix,ny) = floyi(ix,ny) + 
     .                 cftiexclg*cfneut*cfneutsor_ei*2.5*(1.-cfloygwall)*fniy_recy
              endif
              if (matwalli(ix) > 0 .and. recycwit(ix,1,1)>0.) then
                fniy_recy = min(recycwit(ix,1,1)*fac2sp*fniy(ix,0,1), 0.)
                floyi(ix,0) = floyi(ix,0) +
     .                 cftiexclg*cfneut*cfneutsor_ei*2.5*(1.-cfloygwall)*fniy_recy
              endif
            enddo 

         else
            do iy = j1omp1, j5omp # note: cfloyi usually = 2.5 or 1.5 (ExB turb)
               do ix = i4omp, i8omp
                  floyi(ix,iy) = floyi(ix,iy)
     .                            + cfloyi * fniy(ix,iy,ifld)
     .                            + (vyti_use(ix,iy)+vyti_cft(ix,iy))*
     .                                                  0.5*sy(ix,iy)*
     .                              (niy0(ix,iy,ifld)+niy1(ix,iy,ifld))
               end do
             end do
         endif
        end do

c...  Next B x grad(T), first for the ions
      if(abs(cfbgt) .gt. 0 .or. cfeexdbo+cfeixdbo > 0.) then

      do ifld = 1, nfsp
        do iy = j4omp, j8omp
           do ix = i1omp, i5omp
	     iy1 = max(0,iy-1)
             ix1 = ixp1(ix,iy)
             ix5 = ixp1(ix,iy1)
	     if (iy.ne.0 .or. iy.ne.ny+1) then
c... sknam: grad T from tiv
             temp1 = 4.0*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .                 temp1 = 4.0*(tiv(ix,iy) - tiv(ix,iy1))*gyc(ix,iy)

             if (zi(ifld) > 1.e-10) then
                floxibgt(ix,iy,ifld)=(5*sx(ix,iy)/(32*qe*zi(ifld) )) *
     .                           ( ni(ix,iy,ifld) + ni(ix1,iy,ifld) ) *
     .                           ( rbfbt2(ix,iy) + rbfbt2(ix1,iy) ) *
     .                            temp1
             endif
             floxi(ix,iy) = floxi(ix,iy) + cfbgt*floxibgt(ix,iy,ifld)
            endif
            end do
          end do
        end do
 
      do ifld = 1, nfsp
        do iy = j1omp1, j5omp
           do ix = i4omp, i8omp
             ix3 = ixm1(ix,iy)
             ix4 = ixm1(ix,iy+1)
             do jx = 1, nxpt
                if ( (ix.ne.ixlb(jx).and.ixmnbcl.ne.1) .or.
     .               (ix.ne.ixrb(jx)+1.and.ixmxbcl.ne.1) ) then
c... sknam: grad T from tiv
             temp1 = 4.0*(tiv(ix,iy) - tiv(ix3,iy))*gxc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .              temp1 = 4.0*(tiv(ix,iy) - tiv(ix3,iy))*gxc(ix,iy)
             if (zi(ifld) > 1.e-10) then
                floyi(ix,iy) = floyi(ix,iy) - 
     .                      cfbgt*( 5*sy(ix,iy) / (32*qe*zi(ifld) )) * 
     .                          ( ni(ix,iy,ifld) + ni(ix,iy+1,ifld) ) *
     .                          ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) ) *
     .                           temp1
             endif
             endif
             enddo
            end do
          end do
        end do

c...  Now B x grad(T) for the electrons

      do iy = j4omp, j8omp
         do ix = i1omp, i5omp
	     iy1 = max(0,iy-1)
             ix1 = ixp1(ix,iy)
             ix5 = ixp1(ix,iy1)
	     if (iy.ne.0 .or. iy.ne.ny+1) then
c             temp1 = (gtey(ix,iy ) + gtey(ix1,iy ) +
c     .                gtey(ix,iy1) + gtey(ix5,iy1))
c... sknam: grad T from tev
             temp1 = 4.0*(tev(ix,iy) - tev(ix,iy1))*gyc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
             if ( isxptx(ix,iy)==0 .and. iysptrx.gt.0 )
     .                    temp1 = 4.0*(tev(ix,iy) - tev(ix,iy1))*gyc(ix,iy)
             floxebgt(ix,iy) = ( 5*sx(ix,iy) / (32*qe) ) *
     .                           ( ne(ix,iy) + ne(ix1,iy) ) *
     .                          ( rbfbt2(ix,iy) + rbfbt2(ix1,iy) ) *
     .                           temp1
             floxe(ix,iy) = floxe(ix,iy) - cfbgt*floxebgt(ix,iy)
            endif
          end do
        end do
 
      do iy = j1omp1, j5omp
	   do ix = i4omp, i8omp
	    ix3 = ixm1(ix,iy)
	    ix4 = ixm1(ix,iy+1)
            do jx = 1, nxpt
               if ( (ix.ne.ixlb(jx).and.ixmnbcl.ne.1) .or.
     .              (ix.ne.ixrb(jx)+1.and.ixmxbcl.ne.1) ) then
c... sknam: grad T from tev
            temp1 = 4.0*(tev(ix,iy) - tev(ix3,iy))*gxc(ix,iy)
cccMER For full double-null configuration, iysptrx is last closed flux surface.
            if ( isxpty(ix,iy)==0 .and. iysptrx.gt.0 )
     .               temp1 = 4.0*(tev(ix,iy) - tev(ix3,iy))*gxc(ix,iy)
            floye(ix,iy) = floye(ix,iy) + 
     .                       cfbgt*( 5*sy(ix,iy) / (32*qe) ) * 
     .                         ( ne(ix,iy) + ne(ix,iy+1) ) *
     .                         ( rbfbt2(ix,iy) + rbfbt2(ix,iy+1) ) *
     .                          temp1
            endif
            enddo
          end do
        end do
      endif

c...Add the charge-exhange neutral contributions to ion+neutral temp eq.


         do iy = j4omp, j8omp
            do ix = i1omp, i5omp
               floxi(ix,iy) = floxi(ix,iy) +
     .          cftiexclg*cfneut*cfneutsor_ei*cngtgx(1)*cfcvti*2.5*fngx(ix,iy,1)
            end do
            floxi(nx+1,iy) = 0.0e0
        end do
*  --Adds to floyi --

         do iy = j1omp1, j5omp
            do ix = i4omp, i8omp
               floyi(ix,iy) = floyi(ix,iy)
     .             + cftiexclg*cfneut*cfneutsor_ei*cngtgy(1)*2.5*fngy(ix,iy,1)
            end do
        end do

*  ---------------------------------------------------------------------
*  compute the electron and the ion energy flow.
*  ---------------------------------------------------------------------

*  -- compute the electron energy flow --
      if(isteon .eq. 1) call fd2tra (nx,ny,floxe,floye,conxe,conye,
     .                               te,feex,feey,0,methe)

*  -- compute the ion thermal energy flow --
      if(istion .eq. 1) call fd2tra (nx,ny,floxi,floyi,conxi,conyi,
     .                               ti,feix,feiy,0,methi)

c  -- Add rad flux of 4th order diff operator; damp grid-scale oscillations
      if (abs(kye4order)>1.e-50 .or. abs(kyi4order)>1.e-50) then
        do iy = j2p, j5m   #range to iy=1:ny-1 for feey4ord,feiy4ord
          iym1 = max(iy-1,0)
          iyp1 = min(iy+1,ny+1)
          iyp2 = min(iy+2,ny+1)
          do ix = i4omp, i8omp
            dtdym1 = (te(ix,iy)-te(ix,iym1))*gyf(ix,iym1)
            dtdy0 = (te(ix,iyp1)-te(ix,iy))*gyf(ix,iy)
            dtdyp1 = (te(ix,iyp2)-te(ix,iyp1))*gyf(ix,iyp1)
            d2tdy20 = (dtdy0 - dtdym1)*gy(ix,iy)
            d2tdy2p1 = (dtdyp1 - dtdy0)*gy(ix,iyp1)
            d3tdy3 = (d2tdy2p1 - d2tdy20)*gyf(ix,iy)
            feey4ord(ix,iy) = kye4order*d3tdy3*ney1(ix,iy)*
     .                                     sy(ix,iy)/gyf(ix,iy)**2
            feey(ix,iy) = feey(ix,iy) + feey4ord(ix,iy)
          enddo
          do ix = i4omp, i8omp
            dtdym1 = (ti(ix,iy)-ti(ix,iym1))*gyf(ix,iym1)
            dtdy0 = (ti(ix,iyp1)-ti(ix,iy))*gyf(ix,iy)
            dtdyp1 = (ti(ix,iyp2)-ti(ix,iyp1))*gyf(ix,iyp1)
            d2tdy20 = (dtdy0 - dtdym1)*gy(ix,iy)
            d2tdy2p1 = (dtdyp1 - dtdy0)*gy(ix,iyp1)
            d3tdy3 = (d2tdy2p1 - d2tdy20)*gyf(ix,iy)
            feiy4ord(ix,iy) = kyi4order*d3tdy3*niy1(ix,iy,1)*
     .                                     sy(ix,iy)/gyf(ix,iy)**2
            feiy(ix,iy) = feiy(ix,iy) + feiy4ord(ix,iy)
          enddo
        enddo
      endif

c...  Add y-component of nonorthogonal diffusive flux; convective component 
c...  already added to uu(ix,iy)
      if (isnonog .eq. 1) then

         do iy = j1omp1, j6omp
            if (iy .le. ny) then 
            iy1 = max(iy-1,0)
            do ix = i1momp, i6omp 
c...  First do the Te equation
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1)
                 grdnv=(    ( fym (ix,iy,1)*logte(ix2,iy1 ) +  
     .                        fy0 (ix,iy,1)*logte(ix2,iy  ) +
     .                        fyp (ix,iy,1)*logte(ix2,iy+1) +  
     .                        fymx(ix,iy,1)*logte(ix ,iy1 ) +
     .                        fypx(ix,iy,1)*logte(ix ,iy+1) ) 
     .                     -( fym (ix,iy,0)*logte(ix ,iy1 ) +
     .                        fy0 (ix,iy,0)*logte(ix ,iy  ) +
     .                        fyp (ix,iy,0)*logte(ix ,iy+1) +
     .                        fymx(ix,iy,0)*logte(ix4,iy1 ) +  
     .                        fypx(ix,iy,0)*logte(ix6,iy+1) ) ) / 
     .                                                   dxnog(ix,iy)  
               feexy(ix,iy) = exp( 0.5*
     .                         (logte(ix2,iy) + logte(ix,iy)) )* 
     .                               (fcdif*kye+kye_use(ix,iy))*0.5*
     .                                       (ne(ix2,iy)+ne(ix,iy))*
     .                                     (grdnv/cosangfx(ix,iy) - 
     .                         (logte(ix2,iy) - logte(ix,iy))* 
     .                                         gxf(ix,iy))*sx(ix,iy)

c...  Now do the Ti equation.
c --- If we are using the parallel neutral momentum equation, we automatically
c --- change to a combined neutral+ion energy equation. We thus need to
c --- include the neutral heat conductivity. Since it is is isotropic
c --- we could use hcxn though we take the radial derivative; but this is
c --- only true if we dont flux limit.  Thus, we use 4-pt average of hcyn.
c --- Note: this four-point average results in not getting the full Jac. for
c --- a nonorthogonal mesh because of niy1,0 - see def. of hcyn

                 grdnv =(    ( fym (ix,iy,1)*logti(ix2,iy1 ) +  
     .                         fy0 (ix,iy,1)*logti(ix2,iy  ) +
     .                         fyp (ix,iy,1)*logti(ix2,iy+1) +  
     .                         fymx(ix,iy,1)*logti(ix ,iy1 ) +
     .                         fypx(ix,iy,1)*logti(ix ,iy+1) ) 
     .                      -( fym (ix,iy,0)*logti(ix ,iy1 ) +
     .                         fy0 (ix,iy,0)*logti(ix ,iy  ) +
     .                         fyp (ix,iy,0)*logti(ix ,iy+1) +
     .                         fymx(ix,iy,0)*logti(ix4,iy1 ) +  
     .                         fypx(ix,iy,0)*logti(ix6,iy+1) ) ) / 
     .                                                   dxnog(ix,iy)  
               feixy(ix,iy) = exp( 0.5*
     .                       (logti(ix2,iy) + logti(ix,iy)) )*
     .                           ( (fcdif*kyi+kyi_use(ix,iy))*0.5*
     .                                     (nit(ix2,iy)+nit(ix,iy))
     .          + cftiexclg*cfneut*cfneutsor_ei*0.25*(hcyn(ix ,iy)+hcyn(ix ,iy1)
     .                              +hcyn(ix2,iy)+hcyn(ix4,iy1)) ) *
     .                                 (  grdnv/cosangfx(ix,iy)
     .                         - (logti(ix2,iy) - logti(ix,iy))*
     .                                        gxf(ix,iy) )*sx(ix,iy)
c...  Flux limit with flalftxt even though hcys have parallel FL built in
               t0 = max(ti(ix,iy),temin*ev)
               t1 = max(ti(ix2,iy),temin*ev)
               vttn = t0*sqrt( t0/mi(1) )
               vttp = t1*sqrt( t1/mi(1) )
               qfl = flalftxy * (cftiexclg*0.125+(1.-cftiexclg)*0.25) * sx(ix,iy) * (vttn + vttp) * 
     .               (ni(ix,iy,1)+cftiexclg*ng(ix,iy,1)+ni(ix2,iy,1)+cftiexclg*ng(ix2,iy,1))
               feixy(ix,iy) = feixy(ix,iy) /
     .                              sqrt(1. + (feixy(ix,iy)/qfl)**2)

            end do
          endif
        end do

c...  Fix the fluxes with the same indice range as in fd2tra
         do iy = j4omp, j8omp
            do ix = i1momp, i5omp
               feex(ix,iy) = feex(ix,iy) - feexy(ix,iy)
               feix(ix,iy) = feix(ix,iy) - feixy(ix,iy)
            enddo
         enddo

      endif
c...  Finished with nonorthogonal mesh part

c ... Demand that net feex cannot be out of the plates
      if (isfeexpl0 == 1) then
        do iy = j4omp, j8omp
          do jx = 1, nxpt
            if(feex(ixlb(jx),iy) > 0. .and. ixmnbcl==1) then
              feexflr = ni(ixlb(jx),iy,1)*1.e4*ev*sx(ixlb(jx),iy)
              feex(ixlb(jx),iy) = feex(ixlb(jx),iy)/
     .                (1.+ (feex(ixlb(jx),iy)/feexflr)**4)
            endif
            if(feex(ixrb(jx),iy) < 0. .and. ixmxbcl==1) then
              feexflr = ni(ixrb(jx),iy,1)*1.e4*ev*sx(ixrb(jx),iy)
              feex(ixrb(jx),iy) = feex(ixrb(jx),iy)/
     .                (1.+ (feex(ixrb(jx),iy)/feexflr)**4)
            endif
          enddo
        enddo
      endif

      if (isfeixpl0 == 1) then
        do iy = j4omp, j8omp
          do jx = 1, nxpt
            if(feix(ixlb(jx),iy) > 0.) then
              feixflr = ni(ixlb(jx),iy,1)*1.e3*ev*sx(ixlb(jx),iy)
              feix(ixlb(jx),iy) = feix(ixlb(jx),iy)/
     .                (1.+ (feix(ixlb(jx),iy)/feixflr)**4)
            endif
            if(feix(ixrb(jx),iy) < 0.) then
              feixflr = ni(ixrb(jx),iy,1)*1.e3*ev*sx(ixrb(jx),iy)
              feix(ixrb(jx),iy) = feix(ixrb(jx),iy)/
     .                (1.+ (feix(ixrb(jx),iy)/feixflr)**4)
            endif
          enddo
        enddo
      endif

      do iy = j2omp, j5omp
	    if((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
            feex(nxc-1,iy) = 0.
            feix(nxc-1,iy) = 0.
            feex(nxc,iy) = 0.
            feix(nxc,iy) = 0.
            feex(nxc+1,iy) = 0.
            feix(nxc+1,iy) = 0.
         endif
         if (islimon .ne. 0 .and. iy .ge. iy_lims) then
            feex(ix_lim,iy) = 0.
            feix(ix_lim,iy) = 0.
         endif
         if (nxpt==2 .and. ixmxbcl==1) then
            feex(ixrb(1)+1,iy) = 0.
            feix(ixrb(1)+1,iy) = 0.
         endif
         do ix = i2omp, i5omp
c ... ## IJ 2016/10/19 add MC neutral flux
           if(get_neutral_moments .and. cmneutdiv_feg .ne. 0.0) then   
              jfld=1
              seg_ue(ix,iy,jfld)=-( (fegx_ue(ix,iy,jfld)-fegx_ue(ix1,iy,  jfld))
     .                   + fluxfacy*(fegy_ue(ix,iy,jfld)-fegy_ue(ix, iy-1,jfld)) )
     .                  *( (ni(ix,iy,jfld)*ti(ix,iy))/(ni(ix,iy,jfld)*ti(ix,iy)) )
            endif
        end do
      end do

*  -- total energy residual and equipartition --

c...  Electron radiation loss -- ionization and recombination
            do iy = iys1, iyf6  #iys,iyf
               do ix = ixs1, ixf6  #iys, iyf
                  erlizold = erliz(ix,iy)
                  erlrcold = erlrc(ix,iy)
                  ne_sgvi = ne(ix,iy)
                  if (ifxnsgi.eq.1) ne_sgvi = cne_sgvi  # fix density dependence
                  if (istabon==16) then      # compute from b2frates data file
                     zmax=1
                     znuc=1
                     denz(0)=ng(ix,iy,1)     # use ngbackg as below ?
                     denz(1)=ni(ix,iy,1)     # use fac2sp  as below ?
                     dene=ne_sgvi
                     rdum=radmc(zmax,znuc,te(ix,iy),dene,denz,radz)
                     erliz(ix,iy)=chradi*radz(0)*vol(ix,iy)
                     if (isrecmon .ne. 0) erlrc(ix,iy)=chradr*radz(1)*vol(ix,iy)
                  else                       # compute from other data files
                     erliz(ix,iy) = chradi *
     .                           erl1(te(ix,iy),ne_sgvi,rtau(ix,iy))
     .                                  * (ng(ix,iy,1)-ngbackg(1)* 
     .                    (0.9+0.1*(ngbackg(1)/ng(ix,iy,1))**ingb) ) * 
     .                                                       vol(ix,iy)
                     if (isrecmon .ne. 0) erlrc(ix,iy) = chradr *
     .                               erl2(te(ix,iy),ne_sgvi,rtau(ix,iy))
     .                             * fac2sp*ni(ix,iy,1) * vol(ix,iy)
                  endif
                  eeliold = eeli(ix,iy)
                  if (icnuiz.le.1 .and. psor(ix,iy,1).ne.0.) 
     .                                           eeli(ix,iy) = 13.6*ev + 
     .                               erliz(ix,iy)/(fac2sp*psor(ix,iy,1))

                   edisse(ix,iy)=-(1-ishymol*ismolcrm)*ediss*ev*(-0.5*psordis(ix,iy,2)) +
     .                               ishymol*ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                               sv_crumpet(te(ix,iy), ne(ix,iy), 20)

                  pradhyd(ix,iy)= ( (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1)+
     .                                         erlrc(ix,iy) )/vol(ix,iy)
          end do    
        end do

      do iy = iys1, iyf6  #j2, j5
        do ix = ixs1, ixf6  #i2, i5
          vsoreec(ix,iy) =
     .          - cfneut*cfneutsor_ee*cnsor*13.6*ev*fac2sp*psorc(ix,iy,1)
     .          + cfneut*cfneutsor_ee*cnsor*13.6*ev*fac2sp*psorrgc(ix,iy,1)
     .          - cfneut*cfneutsor_ee*cnsor*erliz(ix,iy)
     .          - cfneut*cfneutsor_ee*cnsor*erlrc(ix,iy)
     .          - 2*(0.5-ishymol) * cfneut*cfneutsor_ee*cnsor*cmesore*edisse(ix,iy)
        enddo
      enddo

ccc         if (ishosor.eq.1) then  #full RHS eval
ccc
ccc           if (svrpkg.eq."cvode") then    # cannot access yl(neq+1)
ccc            call xerrab('*** svrpkg=cvode not allowed for ishosor=1 **')
ccc           endif 
ccc
ccc           if (yl(neq+1).lt.0) then  #full RHS eval
ccc
cccc ...    integ. source over cells (but not for Jac) for higher-order accuracy
ccc
ccc             call volave(nx, ny, j2, j5, i2, i5, ixp1(0,0), ixm1(0,0),
ccc     .                         fsprd, vol(0,0), psor_tmpov(0,0), vsoree)
ccc       
ccc           endif   # end of if (yl(neq+1).lt.0) test
ccc         endif    # end of integrating over sources and ishosor test

c*************************************************************
c   Perform 5pt average of source terms as volume integral
c*************************************************************
         if (iseesorave == 0.) then  #use only single-cell value
           do iy = iys1, iyf6
             do ix = ixs1, ixf6
               vsoree(ix,iy) = vsoreec(ix,iy)
             enddo
           enddo

         elseif (iseesorave > 0.) 

            if (xc < 0) then  #full RHS eval
              j2pwr = j2
              j5pwr = j5
            else  # Jacobian
              j2pwr = max(1, yc-1)
              j5pwr = min(ny, yc+1)
            endif 
              if (ParallelPandfCall.gt.0) then
                j2pwr = j2omp
                j5pwr = j5omp
              end if
            do iy = j2pwr, j5pwr
              if (xc < 0) then #full RHS eval
                i2pwr = i2
                i5pwr = i5
              else  #Jacobian eval
                i2pwr = max(1,ixm1(xc,yc))
                i5pwr = min(nx, ixp1(xc,yc))
              endif
              if (ParallelPandfCall.gt.0) then
                i2pwr = i2omp
                i5pwr = i5omp
              end if
              do ix = i2pwr, i5pwr
                ix1 = ixm1(ix,iy)
                ix2 = ixp1(ix,iy)
                vsoree(ix,iy) = (1.-iseesorave*0.5)*
     .                                  vsoreec(ix,iy)+ 
     .                               0.125*iseesorave*vol(ix,iy)*
     .                           ( vsoreec(ix,iy-1)/vol(ix,iy-1) + 
     .                             vsoreec(ix,iy+1)/vol(ix,iy+1) +
     .                             vsoreec(ix1,iy)/vol(ix1,iy)   + 
     .                             vsoreec(ix2,iy)/vol(ix2,iy) )
              enddo
            enddo

         endif    # end of integrating over sources and iseesorave test


      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)
c ... Energy density change due to molecular dissociation ("Franck-Condon")
            emolia(ix,iy,1) = 0
            emolia(ix,iy,2) = 0
            if (ishymol .eq. 1) then
                emolia(ix, iy, 1) = emolia(ix,iy,1)
     .              + ishymol*ismolcrm*ng(ix,iy,2)*vol(ix,iy)* 
     .              sv_crumpet(te(ix,iy), ne(ix,iy), 21) 
     .              * ( psordis(ix,iy,1)/(-2*psordisg(ix,iy,2)))
                emolia(ix, iy, 2) = emolia(ix,iy,2)
     .              + ishymol*ismolcrm*ng(ix,iy,2)*vol(ix,iy)* 
     .              sv_crumpet(te(ix,iy), ne(ix,iy), 21) 
     .              * ( psordis(ix,iy,2)/(-2*psordisg(ix,iy,2)))
            endif

            if (isupgon(1).eq.1) then
c             Set up helper arrays for velocities
              up1cc = 0.5*(up(ix,iy,1)+up(ix1,iy,1))
              upgcc = 0.5*(up(ix,iy,iigsp)+up(ix1,iy,iigsp))
              vycc = (cfnidhgy**0.5)*0.5*(vy(ix,iy,iigsp)+vy(ix1,iy,iigsp))
              v2cc = (cfnidhg2**0.5)*0.5*(v2(ix,iy,iigsp)+v2(ix1,iy,iigsp))

c             IONS
c             -------------------------------------------------------------
c             Ion rate from CX
              psicx(ix,iy) = cfticx*nucx(ix,iy,1)*ng(ix,iy,1)*vol(ix,iy)

c             Ion internal energy sink/source from ioniz & recom
              seit(ix,iy) = 1.5*( 
     .                tg(ix,iy,1) * (psor(ix,iy,1) + psicx(ix,iy))
     .              - ti(ix,iy) * (psorrg(ix,iy,1) +psicx(ix,iy))
     .        )

c             Ion energy source/sink from ioniz & recom
              seik(ix,iy) = cfneut * cfneutsor_ei * cfnidh * 
     .              0.5*mi(1) * ( (up1cc-upgcc)**2 + vycc**2 + v2cc**2 ) * 
     .              (  psor(ix,iy,1) + cftiexclg*psorrg(ix,iy,1)
     .              + (1 + cftiexclg) * psicx(ix,iy) )

c             Ion energy source from mol. diss ("Franck Condon")
              seid(ix,iy) = cftiexclg * cfneut * cfneutsor_ei 
     .              * 2*(0.5 - ishymol) 
     .              * cnsor*eion*(1-ishymol*ismolcrm)*ev*psordis(ix,iy,2) 
     .              + emolia(ix,iy,1) + cftiexclg*emolia(ix,iy,2) # CRM FC

c             Ion energy source from drift heating 
              seidh(ix,iy) = cfnidh2* 
     .              ( -mi(1)*up1cc*upgcc*(psor(ix,iy,1)+psicx(ix,iy))
     .              + 0.5*mi(1)*up1cc**2
     .              * (psor(ix,iy,1)+psorrg(ix,iy,1)+2*psicx(ix,iy)) )

c             ATOMS
c             -------------------------------------------------------------
c             Atom kinetic energy source from recom & CX
              seak(ix,iy) = 0.5*mg(1) * ( (up1cc-upgcc)**2 + vycc**2 + v2cc**2 )
     .              * (psorrg(ix,iy,1)+psicx(ix,iy))

c              Atom kinetic energy source from mol. diss
               sead(ix,iy) = (1-ishymol*ismolcrm)*2*(0.5-ishymol)*(
     .              eion*ev*psordis(ix,iy,2) 
     .          ) + emolia(ix,iy,2) 

c              Atom energy source from drift heating 
               seadh(ix,iy) = cfnidh2* ( -mg(1) *up1cc*upgcc 
     .              * (psorrg(ix,iy,1)+psicx(ix,iy)) 
     .              + 0.5*mg(1) * (upgcc**2 + vycc**2 + v2cc**2)
     .              * (psor(ix,iy,1)+psorrg(ix,iy,1)+2*psicx(ix,iy)) )
            end if
          end do
        end do

*  -- impurity radiation --

      if (isimpon .ge. 2) then
         if (istimingon .eq. 1) tsimp = tick()

         do iy = iys1, iyf6  #if Jacobian, only 1 cell done - local sor
            do ix = ixs1, ixf6
               ntau(ix,iy) = atau(ix,iy) * ne(ix,iy)
               nratio(ix,iy) = ng(ix,iy,1)/ne(ix,iy)
               pradold = pwrzec(ix,iy)
               if (isimpon .eq. 2) then   # fixed-fraction model
                  na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
                  pradcff(ix,iy) = na(ix,iy) * ne(ix,iy) *
     .               emissbs (te(ix,iy), nratio(ix,iy), ntau(ix,iy))
                  pradc(ix,iy) = pradcff(ix,iy)
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (isimpon .eq. 3) then   # average-ion model
                  na(ix,iy) = ni(ix,iy,nhsp+1)
                  pradc(ix,iy) = na(ix,iy) * ne(ix,iy) *
     .               radneq (te(ix,iy), nratio(ix,iy))
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (nzspt .eq. 0) then
                  # Multicharge group not allocated, so avoid radimpmc
                  pradc(ix,iy) = 0.
                  pwrzec(ix,iy) = pradc(ix,iy)
               elseif (isimpon .ge. 4) then  # multi-charge model
                  pradc(ix,iy) = 0.
                  pwrzec(ix,iy) = 0.
                  nsm1 = nhsp 
                  do igsp = nhgsp+1, ngsp  # loop over diff imp species
                     jz = igsp - nhgsp
                     zn = znucl(nsm1+nzsp(jz))
                     if (ngsp .ge. nhgsp) nzloc(0) = ng(ix,iy,igsp)
                     do ifld = 1, nzsp(jz)
                        nzloc(ifld) = ni(ix,iy,nsm1+ifld)
                     enddo
                     nsm1 = nsm1 + nzsp(jz)   # setup for next igsp
                     argth = (te(ix,iy)-1.*ev)/(del_te_ro*ev)
                     fac_rad = 1.
                     if(del_te_ro.lt. 100.) fac_rad=0.5*(1+tanh(argth))
                     if (ismctab .eq. 1) then
                        pwrzec(ix,iy)= pwrzec(ix,iy) + fac_rad*
     .                                radimpmc (nzsp(jz), te(ix,iy), 
     .                                    ne(ix,iy), nzloc, impradloc)
                     elseif (ismctab .eq. 2) then
                        pwrzec(ix,iy)= pwrzec(ix,iy) + fac_rad*
     .                                radmc(nzsp(jz), zn, te(ix,iy), 
     .                                   ne(ix,iy), nzloc, impradloc)
                     endif
                     
                     do ifld = 0, nzsp(jz)
                        pradzc(ix,iy,ifld,jz) = impradloc(ifld)
                        pradc(ix,iy) = pradc(ix,iy)+impradloc(ifld)
                     enddo

                  enddo

		  if (isimpon .eq. 7) then  # add fixed-fraction contrib
                     na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
                     pradcff(ix,iy) = na(ix,iy)* ne(ix,iy)*
     .                     emissbs(te(ix,iy), nratio(ix,iy), ntau(ix,iy))
                     pradc(ix,iy) = pradc(ix,iy) + pradcff(ix,iy)
                     pwrzec(ix,iy) = pwrzec(ix,iy) + pradcff(ix,iy)
                  endif

               endif
            end do
        end do

c*************************************************************
c   Perform 5pt average of source terms as volume integral
c*************************************************************
cc         if (ishosor.eq.0) then  #use only single-cell value
         if (iseesorave.eq.0.) then  #use only single-cell value
           do iy = iys1, iyf6
             do ix = ixs1, ixf6
               pwrze(ix,iy) = pwrzec(ix,iy)
               prad(ix,iy) = pradc(ix,iy)
               do igsp = nhgsp+1, ngsp
                 jz = igsp - nhgsp
                 do ifld = 0, nzsp(jz)
                   pradz(ix,iy,ifld,jz) = pradzc(ix,iy,ifld,jz)
                 enddo
               enddo
             enddo
           enddo

cc         elseif (ishosor .ne. 0) 
         elseif (iseesorave > 0.) 

            if (xc < 0) then  #full RHS eval
              j2pwr = j2
              j5pwr = j5
            else  # Jacobian
              j2pwr = max(1, yc-1)
              j5pwr = min(ny, yc+1)
            endif 
              if (ParallelPandfCall.gt.0) then
                j2pwr = j2omp
                j5pwr = j5omp
              end if
            do iy = j2pwr, j5pwr
              if (xc < 0) then #full RHS eval
                i2pwr = i2
                i5pwr = i5
              else  #Jacobian eval
                i2pwr = max(1,ixm1(xc,yc))
                i5pwr = min(nx, ixp1(xc,yc))
              endif
              if (ParallelPandfCall.gt.0) then
                i2pwr = i2omp
                i5pwr = i5omp
              end if
              do ix = i2pwr, i5pwr
                ix1 = ixm1(ix,iy)
                ix2 = ixp1(ix,iy)
                pwrze(ix,iy) = (1.-iseesorave*0.5)*pwrzec(ix,iy) + 
     .                                            0.125*iseesorave*
     .                         ( pwrzec(ix,iy-1)+ pwrzec(ix,iy+1)+
     .                           pwrzec(ix1,iy) + pwrzec(ix2,iy) )
                if (isimpon < 4) prad(ix,iy) = pwrze(ix,iy)
                if (isimpon >= 4) then  #prad, pradz only diagnostic
                  prad(ix,iy) = (1.-iseesorave*0.5)*pradc(ix,iy) + 
     .                                         0.125*iseesorave*
     .                          ( pradc(ix,iy-1)+ pradc(ix,iy+1)+
     .                            pradc(ix1,iy) + pradc(ix2,iy) )
                  do igsp = nhgsp+1, ngsp
                    jz = igsp - nhgsp
                    do ifld = 0, nzsp(jz)
                      pradz(ix,iy,ifld,jz) = (1.-iseesorave*0.5)*
     .                                        pradzc(ix,iy,ifld,jz) + 
     .                                            0.125*iseesorave*
     .              ( pradzc(ix,iy-1,ifld,jz)+ pradzc(ix,iy+1,ifld,jz)+
     .                pradzc(ix1,iy,ifld,jz) + pradzc(ix2,iy,ifld,jz) )
                    enddo
                  enddo
                endif
              enddo
            enddo

         endif    # end of integrating over sources and iseesorave test

c*******************************************************************
c ... Define a background elec energy source to prevent very low Te
c******************************************************************
      do iy = iys, iyf  #j2, j5
        do ix = ixs, ixf  #i2, i5
          pwrebkgold = pwrebkg(ix,iy)
          if (isimpon == 0) then
            pwrebkg(ix,iy) = (tebg*ev/te(ix,iy))**iteb*pwrbkg_c
          else  #add impurity rad loss
            pwrebkg(ix,iy) = (tebg*ev/te(ix,iy))**iteb*pwrbkg_c
          endif
        enddo
      enddo

         if (istimingon .eq. 1) call timimpfj (tsimp, xc)
      endif  #loop for isimpon==2
 

c*******************************************************************
c ... Define a background ion energy source to prevent very low Ti
c******************************************************************
      do iy = iys, iyf  #j2, j5
        do ix = ixs, ixf  #i2, i5
          pwribkgold = pwribkg(ix,iy)
          pwribkg(ix,iy) = (tibg*ev/ti(ix,iy))**iteb*pwribkg_c
        enddo
      enddo

      END SUBROUTINE calc_plasma_energy



      SUBROUTINE calc_feeiycbo
      IMPLICIT NONE
      Use(Comflo)
      Use(Compla)
      Use(Coefeq)
      Use(Comgeo)
      Use(UEpar)
      Use(Dim)
      Use(Bcond)
      integer  ifld

*  ---------------------------------------------------------------------
*  compute the electron and the ion energy flow.
*  --------------------------------------------------------
      feeycbo =  cfloye*( ne(:,0)*te(:,0)*sy(:,0) ) *
     .      ( (1-cfeeybbo)*cfybf*veycb(:,0) - cfeeydbo*(1-cfydd)*veycp(:,0) )


      do ifld = 1, nfsp
         if ((isupgon(1) .eq. 1) .and. (ifld .eq. iigsp)) then
         else
               feiycbo = feiycbo(:) + cfloyi*fniycbo(:,ifld)*ti(:,0)
         end if
      end do


      END SUBROUTINE calc_feeiycbo

      SUBROUTINE calc_plasma_energy_residuals(xc, yc)
      IMPLICIT NONE
      Use(Selec)
      Use(Rhsides)
      Use(Compla)
      Use(Volsrc)
      Use(Coefeq)
      Use(MCN_sources)
      Use(Conduc)
      Use(Comgeo)
      Use(Poten)
      Use(Share)
      Use(Dim)
      Use(Noggeo)
      Use(Comflo)
      Use(Comtra)
      Use(UEpar)
      Use(Phyvar)
      Use(Xpoint_indices)
      Use(Jacobian_restore)
      Use(Indices_domain_dcl)
      Use(Ext_neutrals)
      Use(Wkspace)
      Use(Imprad)
      Use(Timing)
      Use(Gradients)
      integer xc, yc
      integer iy, ix, iy1, ix1, ix2, ix3, ix4, ix5, ix6, jx, jfld, jz, 
     .  igsp, iy_min, iy_max, j2pwr, j5pwr, i2pwr, i5pwr, ifld, nsm1, zn,
     .  znuc, zmax
      real grdnv, fcd, t0, t1, vttn, vttp, isfe, l0, feexflr, feixflr, 
     .  ne_sgvi, dene, rdum, radmc, radz(0:1), erl1, erl2, 
     .  up1cc, upgcc, vycc, v2cc, tsimp, tick, emissbs, radneq,  
     .  argth, fac_rad, radimpmc, wj, cfwj, telim, thetacc, 
     .  dupdx, dupdy, upxavep1, upxave0, upxavem1, upf0, upfm1,
     .  denz(0:1), sv_crumpet
      external sv_crumpet, radmc
*  ---------------------------------------------------------------------
*  compute the energy residuals.
*  ---------------------------------------------------------------------

*  -- source terms --

      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            resee(ix,iy) = 
     .             seec(ix,iy) + seev(ix,iy) * te(ix,iy)
     .           + pwrsore(ix,iy)
     .           + cmneut * cmneutsor_ee * uesor_te(ix,iy)
     .           - nuvl(ix,iy,1)*vol(ix,iy)*bcee*ne(ix,iy)*te(ix,iy) 
            resei(ix,iy) = 
     .             seic(ix,iy) + seiv(ix,iy) * ti(ix,iy)
     .           + pwrsori(ix,iy)
     .           + cmneut * cmneutsor_ei * uesor_ti(ix,iy)
     .           - nuvl(ix,iy,1)*vol(ix,iy)*bcei*ne(ix,iy)*ti(ix,iy) 
        end do 
      end do  


      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)
            resee(ix,iy) = resee(ix,iy)
     .                  - ( feex(ix,iy) - feex(ix1,iy)
     .          + fluxfacy*(feey(ix,iy) - feey(ix,iy-1)) )
c ... ## IJ 2017 cfneutsor_ei flags above control neutral contrib.
            resei(ix,iy) = resei(ix,iy)					
     .                  - ( feix(ix,iy) - feix(ix1,iy)
     .          + fluxfacy*(feiy(ix,iy) - feiy(ix,iy-1)) )

c ... ## IJ 2016/10/19 add MC neutral flux
           if(get_neutral_moments .and. cmneutdiv_feg .ne. 0.0) then   
              jfld=1
              resei(ix,iy) = resei(ix,iy) +
     .                    cftiexclg*cmneutdiv*cmneutdiv_feg*seg_ue(ix,iy,jfld)
c     .                             cmneutdiv*cmneutdiv_feg*seg_ue(ix,iy,jfld)
            endif
        end do
      end do

      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)
            w0(ix,iy) = vol(ix,iy) * eqp(ix,iy) * (te(ix,iy)-ti(ix,iy))
            resee(ix,iy) = resee(ix,iy) - w0(ix,iy) + vsoree(ix,iy)
            if (isupgon(1).eq.1) then
c             Set up helper arrays for velocities
              up1cc = 0.5*(up(ix,iy,1)+up(ix1,iy,1))
              upgcc = 0.5*(up(ix,iy,iigsp)+up(ix1,iy,iigsp))
              vycc = (cfnidhgy**0.5)*0.5*(vy(ix,iy,iigsp)+vy(ix1,iy,iigsp))
              v2cc = (cfnidhg2**0.5)*0.5*(v2(ix,iy,iigsp)+v2(ix1,iy,iigsp))

c             IONS
c             -------------------------------------------------------------
c             Ion rate from CX
              resei(ix,iy) = resei(ix,iy) + w0(ix,iy)
     .              + seik(ix,iy) 
     .              + (1.0-cftiexclg) * seit(ix,iy)
     .              + seid(ix,iy)
     .              + seidh(ix,iy)

                if (ishymol .eq. 0) then

c                   Ion energy source from mol. drift heating
                    resei(ix,iy) = resei(ix,iy)
     .                  - cftiexclg * cfneut * cfneutsor_ei * cnsor 
     .                  * cfnidhdis * 0.5*mg(1)
     .                  * (upgcc**2 + vycc**2 + v2cc**2) 
     .                  * ( 
     .                      (1-ishymol*ismolcrm)*psordis(ix,iy,2) 
     .                  +   ishymol*ismolcrm * psordis(ix,iy,1)
     .                  )
                endif
            else
               resei(ix,iy) = resei(ix,iy) + w0(ix,iy)
     .             + cfneut * cfneutsor_ei * ctsor*1.25e-1*mi(1)*
     .                    (upi(ix,iy,1)+upi(ix1,iy,1))**2*
     .                    fac2sp*psor(ix,iy,1)
     .             + cfneut * cfneutsor_ei * ceisor*(
     .                  + cmesori*(emolia(ix,iy,1)+emolia(ix,iy,2)) 
     .              )
     .             - cfneut * cfneutsor_ei * ccoldsor*ng(ix,iy,1)*nucx(ix,iy,1)*
     .                    (  1.5*ti(ix,iy)
     .                     - 0.125*mi(1)*(upi(ix,iy,1)+upi(ix1,iy,1))**2
     .                     - eion*ev  ) * vol(ix,iy)
            endif
          end do
        end do


*  -- Energy transfer to impurity neutrals at tg(,,igsp)
      if (ngsp >= 2) then   # for now, specialized to igsp=2 only
        do ifld = nhsp+1, nisp
          do iy = j2omp, j5omp    # iys,iyf limits dont seem to work(?)
            do ix = i2omp, i5omp
              resei(ix,iy) =resei(ix,iy) -cftiimpg*1.5*ni(ix,iy,ifld)*
     .                      (nucxi(ix,iy,ifld)+nueli(ix,iy,ifld))*
     .                      (ti(ix,iy) - tg(ix,iy,2))*vol(ix,iy)
            enddo
          enddo
        enddo
      endif

*  -- impurity radiation --

      if (isimpon .ge. 2) then
         if (istimingon .eq. 1) tsimp = tick()
c******************************************************************
c...  Update resee over whole "box" because initially set to zero 
c******************************************************************
         do iy = j2omp, j5omp
            do ix = i2omp, i5omp           
               resee(ix,iy) = resee(ix,iy) -
     .                            cnimp*pwrze(ix,iy)*vol(ix,iy) +
     .                                pwrebkg(ix,iy)*vol(ix,iy)
            end do
        end do

         if (istimingon .eq. 1) call timimpfj (tsimp, xc)
      endif  #loop for isimpon==2
  
*  -- joule heating --

      if (jhswitch > 0) then  # relies on div(J)=0, so omit iy=1 & ny
         if (isnewpot .eq. 1) then
            iy_min = 2
            iy_max = ny-1
         else
            iy_min = 1
            iy_max = ny
         endif
         if (jhswitch == 1) then   # div(J)=0 gives -grad(phi).J=-div(phi.J)
           do iy = max(iy_min, j2omp), min(iy_max, j5omp)
             do ix = i2, i5   
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               wjdote(ix,iy) = 
     .                          - 0.5*(fqp(ix,iy)+fq2(ix,iy))*
     .                                (phi(ix2,iy)+phi(ix,iy))
     .                          + 0.5*(fqp(ix1,iy)+fq2(ix1,iy))*
     .                                (phi(ix,iy)+phi(ix1,iy))
     .                          - 0.5*fqygp(ix,iy)*
     .                                (phi(ix,iy+1)+phi(ix,iy))
     .                          + 0.5*fqygp(ix,iy-1)*
     .                                (phi(ix,iy)+phi(ix,iy-1))
               resee(ix,iy) = resee(ix,iy) + wjdote(ix,iy) / ( 1. +
     .                             cfwjdotelim*(te(ix,iy)/tebg)**iteb )
             enddo
           enddo
         else  # for jhswitch > 1
           do iy = max(iy_min, j2omp), min(iy_max, j5omp)
             do ix = i2, i5    # use ex*fqx since phi(0,) may be large 
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               wjdote(ix,iy) = 
     .                       0.5*( ex(ix1,iy) *fqx(ix1,iy) +
     .                             ex(ix, iy) *fqx(ix, iy) )/gx(ix,iy)
     .                     + 0.5*( ey(ix, iy) *fqy(ix, iy) +
     .                             ey(ix,iy-1)*fqy(ix,iy-1) )/gy(ix,iy) 
               resee(ix,iy) = resee(ix,iy) + wjdote(ix,iy)
             enddo
           enddo
         endif
      endif

*  -- Now we introduce the viscous heating; one-side derviatives are used
*  -- on either side of the x-point where isxpty = 0

      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            do ifld = 1, nusp  # if nusp --> nfsp, problems from y-term
               ix1 = ixm1(ix,iy)
               ix2 = ixm1(ix,iy+1)
               ix3 = ixm1(ix,iy-1)
	           thetacc = 0.5*(angfx(ix1,iy) + angfx(ix,iy))
	           dupdx = gx(ix,iy)*(upi(ix,iy,ifld)-upi(ix1,iy,ifld))
               wvh(ix,iy,ifld) = cfvcsx(ifld)*cfvisx*cos(thetacc)*
     .                                    visx(ix,iy,ifld)*dupdx**2
               if ( isxpty(ix,iy)==0 ) then  #1-sided deriv down in y
                 dupdy = 0.5*( upi(ix,iy,  ifld)+upi(ix1,iy  ,ifld) -
     .                         upi(ix,iy-1,ifld)-upi(ix3,iy-1,ifld) )*
     .                                                    gyf(ix,iy-1)
               elseif (isxpty(ix,iy)== -1) then #1-sided up in y
                 dupdy = 0.5*( upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld) -
     .                         upi(ix,iy  ,ifld)-upi(ix1,iy  ,ifld) )*
     .                                                    gyf(ix,iy)
               elseif (isxpty(ix,iy)==1.and.isvhyha==1) then
                                 #use harm y-ave for up face-values
                                 #take abs() to avoid near-zero denomin;
                                 #small err in wvh because up then small
                 upxavep1 = 0.5*(upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld))
                 upxave0 =  0.5*(upi(ix,iy  ,ifld)+upi(ix1,iy  ,ifld))
                 upxavem1 = 0.5*(upi(ix,iy-1,ifld)+upi(ix3,iy-1,ifld))
                 upf0  = 2.*upxavep1*upxave0*(upxavep1+upxave0) /
     .                           ( (upxavep1+upxave0)**2 + upvhflr**2 )
                 upfm1 = 2.*upxave0*upxavem1*(upxave0+upxavem1) /
     .                           ( (upxave0+upxavem1)**2 + upvhflr**2 )
                 dupdy = (upf0 - upfm1)*gy(ix,iy)
               else	#V7.08.04 option - linear ave in y-direction
		 dupdy = 0.25*( (upi(ix,iy+1,ifld)+upi(ix2,iy+1,ifld) - 
     .                           upi(ix,iy  ,ifld)-upi(ix1,iy  ,ifld))*
     .                                                     gyf(ix,iy) +
     .                          (upi(ix,iy  ,ifld)+upi(ix1,iy  ,ifld) -
     .                           upi(ix,iy-1,ifld)-upi(ix3,iy-1,ifld))*
     .                                                     gyf(ix,iy-1) )
               endif
               wvh(ix,iy,ifld) = wvh(ix,iy,ifld) + cfvcsy(ifld)*cfvisy*
     .                                   visy(ix,iy,ifld)*dupdy**2
	           wvh(ix,iy,ifld) = wvh(ix,iy,ifld) -
     .                             sin(thetacc)*cfvcsy(ifld)*cfvisy*
     .                                   visy(ix,iy,ifld)*dupdx*dupdy
            if (zi(ifld)==0.0 .and. ifld.eq.iigsp) then 
              resei(ix,iy) = resei(ix,iy) + cftiexclg*wvh(ix,iy,ifld)*vol(ix,iy)
            else
              resei(ix,iy) = resei(ix,iy) + wvh(ix,iy,ifld)*vol(ix,iy)
            endif
            end do   # loop over up species ifld
          end do
        end do
        


c*******************************************************************
c ... Define a background ion energy source to prevent very low Ti
c******************************************************************
      do iy = j2omp, j5omp
        do ix = i2omp, i5omp
          resei(ix,iy) = resei(ix,iy) + pwribkg(ix,iy)*vol(ix,iy)
        enddo
      enddo


      END SUBROUTINE calc_plasma_energy_residuals


      SUBROUTINE calc_plasma_heatconductivities
      IMPLICIT NONE
      Use(Selec)
      Use(Conduc)
      Use(Dim)
      Use(Compla)
      Use(Wkspace)
      Use(Comgeo)
      Use(UEpar)
      Use(Comtra)
      Use(Phyvar)
      Use(Xpoint_indices)
      Use(Indices_domain_dcl)
      Use(Share)
      Use(Comflo)
      Use(Coefeq)
      Use(MCN_sources)
      integer iy, ix, ifld, jfld, ix1, iyp1, jx, iy1
      real tv, a, fxet, fxit, niavex, niavey, kyemix, fcd, 
     .  kyimix, tiave, lmfpi, wallfac, qflx, cshx, lxtic, qshx, 
     .  teave, zeffave, lmfpe, neavex, tgavex, tgavey, noavex, 
     .  noavey, lmfpn, qfly, cshy, qshy

*****************************************************************
*****************************************************************
*  Heat Conduction. (old PHYTHC)
*****************************************************************
*  ---------------------------------------------------------------------
*  compute conductivities on cell faces
*  ---------------------------------------------------------------------

*  -- initialize to 0 --

      do iy = j1omp1, j6omp
         do ix = i1momp, i6omp
            hcxe(ix,iy) = 0.0e0
            hcxi(ix,iy) = 0.0e0
            hcxineo(ix,iy) = 0.0e0
            hcye(ix,iy) = 0.0e0
            hcyi(ix,iy) = 0.0e0
            do ifld = 1, nisp
               hcxij(ix,iy,ifld) = 0.0e0
               hcyij(ix,iy,ifld) = 0.0e0
            enddo
        end do
      end do

*  -- loop over species number --

      do ifld = 1, nisp
c -- Skip this if these are the neutrals (zi(ifld).eq.0)
         if (zi(ifld) .ne. 0.0e0) then

c...  Initialize w1 and w2 for each species
       w1 = 0.0e0
       w2 = 0.0e0

*     -- conductivities --
*        The poloidal conductivities  are initially computed without
*        the factor rr**2 * tv**2.5)

         do jfld = 1, nisp
            tv = zi(jfld)**2
            a = zi(jfld)**2 *
     .            sqrt(2*mi(ifld)*mi(jfld)/(mi(ifld)+mi(jfld)))
            do  iy = j1omp1, j6omp
               do ix = i1, i6
                  ix1 = ixp1(ix,iy)
                  w1(ix,iy) = w1(ix,iy) + tv*(ni(ix,iy,jfld)*gx(ix,iy) +
     .                                     ni(ix1,iy,jfld)*gx(ix1,iy)) / 
     .                                        (gx(ix,iy) + gx(ix1,iy))
                  w2(ix,iy) = w2(ix,iy) + a*(ni(ix,iy,jfld)*gx(ix,iy) +
     .                                     ni(ix1,iy,jfld)*gx(ix1,iy)) / 
     .                                        (gx(ix,iy) + gx(ix1,iy))
                end do
            end do
        end do

         do iy = j1omp1, j6omp
            do ix = i1momp, i6omp
               ix1 = ixp1(ix,iy)
               iyp1 = min(ny+1, iy+1)
               ctaue(ix,iy,ifld) = 3.5e11*zi(ifld)/loglambda(ix,iy)
               ctaui(ix,iy,ifld) =2.1e13/(loglambda(ix,iy)*zi(ifld)**2)
               fxe = kxe * ce * ctaue(ix,iy,ifld) / (me*ev*sqrt(ev))
               fxi = kxi * ci * ctaui(ix,iy,ifld) / (ev*sqrt(ev*mp))
               fxet = fxe
               fxit = fxi
               do jx = 1, nxpt  #reduce kxe inside sep by rkxecore fac
                  if ( (iy.le.iysptrx) .and. 
     .                    ix.gt.ixpt1(jx) .and. ix.le.ixpt2(jx) ) then
                     fxet = fxe/( 1. + (rkxecore-1.)*
     .                              (yyf(iy)/(yyf(0)+4.e-50))**inkxc )
                     fxit = kxicore * fxi
                  endif
               enddo
               niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                    ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                     (gx(ix,iy) + gx(ix1,iy))
               niavey = ( niy0(ix,iy,ifld)*gy(ix,iy) +
     .                    niy1(ix,iy,ifld)*gy(ix,iyp1)) /
     .                     (gy(ix,iy) + gy(ix,iyp1))
               hcxe(ix,iy) = hcxe(ix,iy)+fxet*niavex/w1(ix,iy)

c ... Use fixed diffusivity inside the separatrix, anomalous outside,
c     if anomalous-diffusivity multiplier is nonzero.
               kyemix = fcdif*kye + kye_use(ix,iy)
cccMER NOTE: when there are multiple x-points, as in 'dnull' configuration,
ccc          iysptrx is the last closed flux surface (see S.R. nphygeo)
               if(kyet .gt. 1.e-20 .and. iy .gt. iysptrx) then
                  kyemix = (1. - ckyet) * kyemix +
     .               ckyet * kyet * diffusivwrk(ix,iy)
               endif
               hcye(ix,iy) = hcye(ix,iy) + ( kyemix +
     .                       2.33*(dclass_e(ix,iy)+dclass_e(ix,iyp1)) )*
     .                                               zi(ifld) * niavey

               hcxij(ix,iy,ifld) = fxit*niavex/w2(ix,iy)
               kyimix = fcdif*kyi + kyi_use(ix,iy)
cccMER NOTE: when there are multiple x-points, as in 'dnull' configuration,
ccc          iysptrx is the last closed flux surface (see S.R. nphygeo)
               if(kyit .gt. 1.e-20 .and. iy .gt. iysptrx) then
                  kyimix = (1. - ckyit) * kyimix +
     .               ckyit * kyit * diffusivwrk(ix,iy)
               endif
               hcyij(ix,iy,ifld) = hcyij(ix,iy,ifld) + ( kyimix +
     .                           (dclass_i(ix,iy)+dclass_i(ix,iyp1)) )*
     .                                                         niavey
              end do  
            end do
          end if 
        end do

  
c ... Add ion temp. dep. for pol. terms, flux limit, & build total ion hcx,yi
        do ifld = 1, nisp
         if (zi(ifld) .ne. 0.e0) then
         do iy = j1omp1, j6omp
            do ix = i1momp, i6omp
               ix1 = ixp1(ix,iy)
               if (concap .eq. 0) then
                  tiave = (ti(ix,iy)*gx(ix,iy) + ti(ix1,iy)*gx(ix1,iy)) /
     .                                         (gx(ix,iy) + gx(ix1,iy))
                  do jx = 1, nxpt
                    if (ix==ixlb(jx).and.ixmnbcl==1) tiave=ti(ixlb(jx)+1,iy)
                    if (ix==ixrb(jx).and.ixmxbcl==1) tiave=ti(ixrb(jx),iy)
                  enddo
                  a = max (tiave, temin*ev)
               else
                  a = afix*ev
               endif
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld)*rrv(ix,iy)*
     .                                            rrv(ix,iy)*a*a*sqrt(a)
c...  reduce hcxij if ti very flat; prevents large conduction for high ti
c...  or if lmfpi exceeds a mean-free path limit, lmfplim
               lmfpi = 1.e16*(tiave/ev)**2/ni(ix,iy,1) # i-mean-free-path
               niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                    ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                                     (gx(ix,iy) + gx(ix1,iy))
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld)/(1.+lmfpi/lmfplim)
               hcxij(ix,iy,ifld) = hcxij(ix,iy,ifld) *
     .                         ( cutlo + (ti(ix,iy)-ti(ix1,iy))**2 ) /
     .                         ( cutlo + (ti(ix,iy)-ti(ix1,iy))**2 +
     .                      (0.5*alfkxi*(ti(ix,iy)+ti(ix1,iy)))**2 ) +
     .                         kxi_use(ix,iy)*niavex

c ... Flux limit individ. hcxij in poloidal direction if isflxldi=2
               if (isflxldi .eq. 2) then
                  niavex = ( ni(ix ,iy,ifld)*gx(ix ,iy) +
     .                       ni(ix1,iy,ifld)*gx(ix1,iy)) /
     .                                        (gx(ix,iy) + gx(ix1,iy))
                  wallfac = 1.
                  do jx = 1, nxpt
                     if ( ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                      (ix==ixrb(jx).and.ixmxbcl==1) )
     .                    .and. (isplflxl==0) ) wallfac = flalfipl/flalfi
                  enddo
                  qflx = wallfac*flalfi * rrv(ix,iy) *
     .                                    sqrt(a/mi(ifld)) * niavex * a
                  cshx = hcxij(ix,iy,ifld)
	          lxtic = 0.5*(ti(ix,iy)+ti(ix1,iy)) /
     .               (abs(ti(ix,iy)-ti(ix1,iy))*gxf(ix,iy) + 100.*cutlo)
	          qshx = cshx * (ti(ix,iy)-ti(ix1,iy)) * gxf(ix,iy) *
     .                                              (1. + lxtic/lxtimax)
                  hcxij(ix,iy,ifld) = cshx  / (1 + abs(qshx/qflx))
               endif
               hcxi(ix,iy) = hcxi(ix,iy) + hcxij(ix,iy,ifld)
               hcxineo(ix,iy) = hcxineo(ix,iy) + hcxij(ix,iy,ifld)*
     .                      1.5676*epsneo(ix,iy)**1.5/k2neo(ix,iy,ifld)
               hcyi(ix,iy) = hcyi(ix,iy) + hcyij(ix,iy,ifld)
               qipar(ix,iy,ifld) = hcxij(ix,iy,ifld)*gxf(ix,iy)*
     .                                (ti(ix,iy)-ti(ix1,iy))/rrv(ix,iy)
            enddo
         enddo
        end if
        end do

c...  Now include elec. temp and other dep. in poloidal terms + diff. neut.
      do iy = j1omp1, j6omp
         do ix = i1momp, i6omp
            ix1 = ixp1(ix,iy)
            iyp1 = min(ny+1, iy+1)
            if (concap .eq. 0) then
               teave = (te(ix,iy)*gx(ix,iy) + te(ix1,iy)*gx(ix1,iy)) /
     .                                       (gx(ix,iy) + gx(ix1,iy))
               do jx = 1, nxpt
                 if(ix==ixlb(jx).and.ixmnbcl==1) teave=te(ixlb(jx)+1,iy)
                 if(ix==ixrb(jx).and.ixmxbcl==1) teave=te(ixrb(jx),iy)
               enddo
               a = max (teave, temin*ev)
            else
               a = afix*ev
            endif
            zeffave = (zeff(ix,iy)*gx(ix,iy) + zeff(ix1,iy)*gx(ix1,iy)) /
     .                                         (gx(ix,iy) + gx(ix1,iy))
            zcoef = 0.308 + 0.767*zeffave - 0.075*zeffave**2
            hcxe(ix,iy) = hcxe(ix,iy)*rrv(ix,iy)*rrv(ix,iy)*a*a*sqrt(a)
     .                    *zcoef

c...  reduce hcxe if te very flat; prevents very large conduction for high te
c...  or if the mean-free path exceeds lmfplim
            lmfpe = 2e16*(te(ix,iy)/ev)**2/ne(ix,iy)  #mfp for elec [m]
            neavex = ( ne(ix ,iy)*gx(ix ,iy) +
     .                 ne(ix1,iy)*gx(ix1,iy))/(gx(ix,iy) + gx(ix1,iy))
            hcxe(ix,iy) = hcxe(ix,iy) *
     .                         ( cutlo + (te(ix,iy)-te(ix1,iy))**2 ) /
     .                         ( cutlo + (te(ix,iy)-te(ix1,iy))**2 +
     .                     (0.5*alfkxe*(te(ix,iy)+te(ix1,iy)))**2 ) +
     .                        kxe_use(ix,iy)*neavex
            hcxe(ix,iy) = hcxe(ix,iy) /( (1. + lmfpe/lmfplim) *
     .                  (1+hcxe(ix,iy)*gx(ix,iy)**2*tdiflim/ne(ix,iy)) )
            if (isupgon(1).eq.0) then   # add diff. gas cx contrib. to hcxi
               hcxn(ix,iy) = 0.
               hcyn(ix,iy) = 0.
c..1dn0802
c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
               hcxi(ix,iy) = hcxi(ix,iy)
     .           + cftiexclg*cfneut*cfneutsor_ei*kxn*( ng(ix ,iy,1)*ti(ix ,iy)
     .                              +ng(ix1,iy,1)*ti(ix1,iy) ) /
     .                  (mi(1)*(nucx(ix,iy,1) + nucx(ix1,iy,1)))
               hcyi(ix,iy) = hcyi(ix,iy)
     .           + cftiexclg*cfneut*cfneutsor_ei*kyn*( ngy0(ix,iy,1)*tiy0(ix,iy)
     .                              +ngy1(ix,iy,1)*tiy1(ix,iy) ) /
     .                  (mi(1)*(nucx(ix,iy,1) + nucx(ix,iyp1,1)))
            endif
          end do
        end do

c
c
      if (isupgon(1).eq.1) then
c
c ----- Section for the inertial neutral fluid; we need to do different
c ----- things than for the ions. Note third index=iigsp is neutral species
c ----- The inertial neutrals coeff. are flux-limited and add to total here
         do iy = j1omp1, j6omp
            iy1 = min(iy,ny)   #dont use j5 because hcx also in loop (not imp.)
            do ix = i1momp, i6omp
               ix1 = ixp1(ix,iy)
               tgavex = max(0.5*(tg(ix,iy,1) + tg(ix1,iy,1)), temin*ev)
               tgavey= max(0.5*(tgy0(ix,iy,1)+tgy1(ix,iy,1)), temin*ev)
               niavex = 0.5*(ni(ix,iy,1) + ni(ix1,iy,1)) #only for coll. term
               niavey = 0.5*(niy0(ix,iy1,1) + niy1(ix,iy1,1)) #only coll. term
               noavex = ( ni(ix ,iy,iigsp)*gx(ix ,iy) +
     .                    ni(ix1,iy,iigsp)*gx(ix1,iy)) /
     .                     (gx(ix,iy) + gx(ix1,iy))
               noavey = 0.5*(niy0(ix,iy1,iigsp) + niy1(ix,iy1,iigsp))

c          Set up flux-limit variables (no rrv here) 
c          First limit the poloidal coeff, then radial
c IJ 2016/10/10	add cfneutsor_ei multiplier to control fraction of neutral energy to add
               qflx = flalftgxa(ix) * sqrt(tgavex/mi(iigsp)) * noavex *
     .                                                     tgavex
               lmfpn = 1./(sigcx * (niavex + rnn2cx*noavex))
               cshx = lmfpn*sqrt(tgavex/mi(iigsp))*noavex * 
     .                         lgtmax(iigsp)/(lmfpn + lgtmax(iigsp))
               qshx = cshx * (tg(ix,iy,1)-tg(ix1,iy,1)) * gxf(ix,iy)
	       hcxn(ix,iy) = cshx  / 
     .                      (1 + (abs(qshx/qflx))**flgamtg)**(1./flgamtg)
               hcxi(ix,iy) = hcxi(ix,iy) + 
     .                          cftiexclg*cfneut*cfneutsor_ei*hcxn(ix,iy)
c          Now for the radial flux limit - good for nonorthog grid too
               qfly = flalftgya(iy) * sqrt(tgavey/mi(iigsp)) * noavey *
     .                                                     tgavey
               lmfpn = 1./(sigcx * (niavey + rnn2cx*noavey))
               cshy = lmfpn*sqrt(tgavey/mi(iigsp))*noavey * 
     .                         lgtmax(iigsp)/(lmfpn + lgtmax(iigsp))
               qshy = cshy * (tgy0(ix,iy1,1)-tgy1(ix,iy1,1))/dynog(ix,iy)
               hcyn(ix,iy) = cshy / 
     .                      (1 + (abs(qshy/qfly))**flgamtg)**(1./flgamtg)
               hcyi(ix,iy) = hcyi(ix,iy) + 
     .                          cftiexclg*cfneut*cfneutsor_ei*hcyn(ix,iy)
c     
          end do
        end do
      endif
      END SUBROUTINE calc_plasma_heatconductivities

      SUBROUTINE calc_plasma_equipartition
      IMPLICIT NONE
      Use(Selec)
      Use(Wkspace)
      Use(Dim)
      Use(Compla)
      Use(UEpar)
      Use(Phyvar)
      Use(Comtra)
      Use(Conduc)
      Use(Share)
      integer iy, ix, ifld, ix2
      real tv, a, loglmcc
*  Equipartition (old PHYEQP)
*****************************************************************
*  ---------------------------------------------------------------------
*  compute equipartition.
*  ---------------------------------------------------------------------

*     -- initialize w3 --
      w3=0


*  -- compute equipartition --
ccc In detail, coef1 = qe**4*sqrt(me)*lnlam / ((2*pi)**1.5*eps0**2)
      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            do ifld = 1, nisp
              tv = zi(ifld)**2/mi(ifld)
              w3(ix,iy) = w3(ix,iy) + tv*ni(ix,iy,ifld)
            end do
            ix2 = ixm1(ix,iy)
            a = max (te(ix,iy), temin*ev)
            loglmcc = 0.5*(loglambda(ix,iy)+loglambda(ix2,iy))
            coef1 = feqp*4.8e-15*loglmcc*sqrt(ev)*ev*mp
            eqp(ix,iy) = coef1 * w3(ix,iy) * ne(ix,iy) / (a*sqrt(a))
c...       reduce eqp when (te-ti)/(te+ti) << 1
            eqp(ix,iy) = eqp(ix,iy) * (a-ti(ix,iy))**2 / ( cutlo +
     .                (a-ti(ix,iy))**2 + (alfeqp*(a+ti(ix,iy)))**2 )
        end do
      end do
      END SUBROUTINE calc_plasma_equipartition


