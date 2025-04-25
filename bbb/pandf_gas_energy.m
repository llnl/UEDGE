c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"


c --------------------------------------------------------------------------
c SUBROUTINE TO SET UP ENERGY EQUATION FOR NEUTRAL GAS
c --------------------------------------------------------------------------

      subroutine engbalg

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,ffyi,ffyo
      real vt0,vt1,wallfac,lxtgc,dupdx,dupdy,fniy_recy,thetacc
      real vttn,vttp
      integer ifld,iixt,iy1, methgx, methgy, iy2, jx
      #Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t0,t1,t2,tv,a
      real uuxgcc, vygcc, v2gcc, upgcc, vycc, v2cc
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx1
      Use(Share)    # geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy,nlimiy,rld2dx,rld2dy

      Use(Coefeq)   # cngfx,cngfy,cngmom,cmwall,cdifg,rld2dxg,rld2dyg
      Use(Bcond)    # albedoo,albedoi
      Use(Parallv)  # nxg,nyg
      Use(Rccoef)   # recylb,recyrb,recycw,recycz,sputtr
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,i8,j1,j2,j3,j4,j5,j6,j7,j8
                    # xlinc,xrinc,yinc,ixm1,ixp1,stretcx
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv,eiamoldiss
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
      Use(Wkspace)  # w0,w1,etc
      Use(MCN_dim)  #
      Use(MCN_sources)   # cfneutsor_ei

*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

*  -- Compute sources terms v_grad_Pg; first initialize array over range
c  -- This v_grad_Pg term first added by MZhao

      do igsp = 1, ngsp
        do iy = j2, j5
          do ix = i2, i5
            segc(ix,iy,igsp) = 0.0
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
        if(istgon(igsp) == 1) then 
          do iy = j2, j5
            do ix = i2, i5
              ix1 = ixm1(ix,iy)
              ix2 = ixp1(ix,iy)
              iy1 = max(0,iy-1)
              tv = (pg(ix2,iy,igsp) - pg(ix ,iy,igsp))
              t1 = (pg(ix ,iy,igsp) - pg(ix1,iy,igsp))
              segc(ix,iy,igsp) = 0.5*cvgpg*( 
     ,                 uuxg(ix, iy,igsp)*ave(gx(ix2,iy),gx(ix, iy))*tv +
     .                 uuxg(ix1,iy,igsp)*ave(gx(ix ,iy),gx(ix1,iy))*t1 )*
     .                                                      vol(ix,iy)
              t2 = cvgpg*0.5*( vyg(ix, iy,igsp)*dynog(ix, iy)*
     .                       (pgy1(ix, iy,igsp)-pgy0(ix, iy,igsp)) +
     .                   vyg(ix,iy1,igsp)*dynog(ix,iy1)*
     .                       (pgy1(ix,iy1,igsp)-pgy0(ix,iy1,igsp)) )
              segc(ix,iy,igsp)=segc(ix,iy,igsp) + cvgpg*t2*vol(ix,iy)
            enddo
          enddo
        endif
      enddo



*  -- Compute flux terms; initialize some arrays to 0 --

      do igsp = 1, ngsp
        do iy = j1, j6
          do ix = i1, i6
            floxge(ix,iy,igsp) = 0.0e0
            floyge(ix,iy,igsp) = 0.0e0
            conxge(ix,iy,igsp) = 0.0e0
            conyge(ix,iy,igsp) = 0.0e0
          enddo
        enddo
      enddo

*  ---------------------------------------------------------------------
*  Compute thermal conductances
*  ---------------------------------------------------------------------
c ... Compute poloidal conduction
      do igsp = 1,ngsp
        do iy = j4, j8
          do ix = i1, i5
            ix2 = ixp1(ix,iy)

            t0 = max (tg(ix,iy,igsp), temin*ev)
            t1 = max (tg(ix2,iy,igsp), temin*ev)
            vt0 = sqrt(t0/mg(igsp))
            vt1 = sqrt(t1/mg(igsp))
c... flux-limit occurs in building hcxg - do not flux-limit 2nd time
            conxge(ix,iy,igsp) = sx(ix,iy) * hcxg(ix,iy,igsp) * gxf(ix,iy) 
          enddo
          conxge(nx+1,iy,igsp) = 0
        enddo
      enddo

*  -- compute radial conduction conyge
      do igsp = 1, ngsp
        do iy = j1, j5
          do ix = i4, i8
            conyge(ix,iy,igsp) = sy(ix,iy)*hcyg(ix,iy,igsp)/dynog(ix,iy)
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
        do ix = i1, i6
          conyge(ix,ny+1,igsp) = 0.0e0
        enddo
      enddo
*  ---------------------------------------------------------------------
*  compute convective flow of tg
*  ---------------------------------------------------------------------

*  -- compute floxge --

      do igsp = 1, ngsp
        do iy = j4, j8
          do ix = i1, i5
            floxge(ix,iy,igsp) = cfcvtg*2.5*fngx(ix,iy,igsp)
          enddo
          floxge(nx+1,iy,igsp) = 0.
        enddo
      enddo

*  -- Correct bdry:remove any inward power from plates; ok in parallel
      do igsp = 1, ngsp
        do iy = j4, j8
          do jx = 1, nxpt  #if at plate, sub (1-cfloxiplt)*neut-contrib
            if(ixmnbcl==1) then  #real div plt -need for parallel UEDGE
              iixt = ixlb(jx) #left plate
              if(fngx(iixt,iy,igsp) > 0.) then
                floxge(iixt,iy,igsp) = floxge(iixt,iy,igsp) -
     .                    (1.-cfloxiplt)*cfcvti*2.5*fngx(iixt,iy,igsp)
              endif
            endif
            if(ixmxbcl==1) then #real div plt -need for parallel UEDGE
              iixt = ixrb(jx) # right plate
              if(fngx(iixt,iy,igsp) < 0.) then
                floxge(iixt,iy,igsp) = floxge(iixt,iy,igsp) -
     .                   (1.-cfloxiplt)*cfcvti*2.5*fngx(iixt,iy,igsp)
              endif
              floxge(ixrb(jx)+1,iy,igsp) = 0.0e0 #cosmetic
            endif
          enddo
        enddo
      enddo

*  -- compute floyge --

      do igsp = 1, ngsp
        do iy = j1, j5
          do ix = i4, i8
            floyge(ix,iy,igsp) = cfcvtg*2.5*fngy(ix,iy,igsp)
          enddo
        enddo
      enddo
 
*  -- Correct bdry:remove any inward power from walls; ok in parallel
      do igsp = 1, ngsp
        do ix = i4, i8
          do jx = 1, nxpt  #if on PF wall, sub (1-cfloygw)*neut-contrib
            if(iymnbcl==1) then  #real PFw-need for parallel UEDGE?
              if(ix <= ixpt1(jx) .or. ix > ixpt2(jx)) then
                if(fngy(ix,0,igsp) > 0.) then
                  floyge(ix,0,igsp) = floyge(ix,0,igsp) -
     .                       (1.-cfloygwall)*cfcvtg*2.5*fngy(ix,0,igsp)
                endif
              endif
            endif
          enddo
          if(iymxbcl==1) then #real outer-w-need for parallel UEDGE?
            if(fngy(ix,ny,igsp) < 0.) then
              floyge(ix,ny,igsp) = floyge(ix,ny,igsp) -
     .                    (1.-cfloygwall)*cfcvtg*2.5*fngy(ix,ny,igsp)
            endif
          endif
          floyge(ix,ny+1,igsp) = 0.0e0 #cosmetic
        enddo   #ix loop
      enddo

        fegx = 0
        fegxy = 0
*  -- Combine conduction/convection to compute thermal energy flow --
      do igsp = 1,ngsp
        if(istgon(igsp) == 1) then
          call fd2tra (nx,ny,floxge(0:nx+1,0:ny+1,igsp), 
     .          floyge(0:nx+1,0:ny+1,igsp), conxge(0:nx+1,0:ny+1,igsp),
     .          conyge(0:nx+1,0:ny+1,igsp),tg(0:nx+1,0:ny+1,igsp),
     .          fegx(0:nx+1,0:ny+1,igsp),fegy(0:nx+1,0:ny+1,igsp),
     .          0,methi)
        endif
      enddo

c...  Add y-component of nonorthogonal diffusive flux; convective component 
c...  already added to uug(ix,iy,igsp)
      if (isnonog == 1) then
        do igsp = 1, ngsp
          if(istgon(igsp) == 0) cycle
          do iy = j1, j6
            if (iy .gt. ny) cycle
            iy1 = max(iy-1,0)
            do ix = i1, i6
              ix1 = ixm1(ix,iy)
              ix2 = ixp1(ix,iy)
              ix3 = ixm1(ix,iy1)
              ix4 = ixp1(ix,iy1)
              ix5 = ixm1(ix,iy+1)
              ix6 = ixp1(ix,iy+1)
              t0 = max(tg(ix,iy,igsp),tgmin*ev) 
              t1 = max(tg(ix2,iy,igsp),tgmin*ev)
              vtn = sqrt( t0/mg(igsp) )
              vtnp = sqrt( t1/mg(igsp) )
              nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
              nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)

c --- Note: this four-point average results in not getting the full Jac. for
c --- a nonorthogonal mesh because of ngy1,0 - see def. of hcyn

                grdnv =( ( fym (ix,iy,1)*log(tg(ix2,iy1 ,igsp)) +  
     .                     fy0 (ix,iy,1)*log(tg(ix2,iy  ,igsp)) +
     .                     fyp (ix,iy,1)*log(tg(ix2,iy+1,igsp)) +  
     .                     fymx(ix,iy,1)*log(tg(ix ,iy1 ,igsp)) +
     .                     fypx(ix,iy,1)*log(tg(ix ,iy+1,igsp)) ) 
     .                  -( fym (ix,iy,0)*log(tg(ix ,iy1 ,igsp)) +
     .                     fy0 (ix,iy,0)*log(tg(ix ,iy  ,igsp)) +
     .                     fyp (ix,iy,0)*log(tg(ix ,iy+1,igsp)) +
     .                     fymx(ix,iy,0)*log(tg(ix4,iy1 ,igsp)) +  
     .                     fypx(ix,iy,0)*log(tg(ix6,iy+1,igsp)) ) ) / 
     .                                                  dxnog(ix,iy)  
               difgx2 = ave( tg(ix ,iy,igsp)/nu1,
     .                       tg(ix2,iy,igsp)/nu2 )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))

               fegxy(ix,iy,igsp) = cfegxy*exp( 0.5*
     .                 (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                      difgx2*ave(ng(ix2,iy,igsp),ng(ix,iy,igsp))*
     .                                 ( grdnv/cos(angfx(ix,iy))
     .                  - (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))*
     .                                        gxf(ix,iy) )*sx(ix,iy)
c...  Flux limit with flalftxt even though hcys have parallel FL built in
               t0 = max(tg(ix,iy,igsp),tgmin*ev)
               t1 = max(tg(ix2,iy,igsp),tgmin*ev)
               vttn = t0*sqrt( t0/mg(igsp) )
               vttp = t1*sqrt( t1/mg(igsp) )
       if(isfegxyqflave == 0) then
               qfl = flalftgxy(igsp)*0.25*sx(ix,iy) * (vttn+vttp) * 
     .                              (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
       else  #use harmonic average of T*vt and ng to face
               qfl = flalftgxy(igsp) * sx(ix,iy) * ave(vttn,vttp) * 
     .                            ave(ng(ix,iy,igsp),ng(ix2,iy,igsp))
       endif
               fegxy(ix,iy,igsp) = fegxy(ix,iy,igsp) /
     .                          sqrt(1. + (fegxy(ix,iy,igsp)/qfl)**2)
               fegx(ix,iy,igsp) = fegx(ix,iy,igsp) - fegxy(ix,iy,igsp)
             enddo  #loop over ix
           enddo    #loop over iy
        enddo       #loop over igsp
      endif         #if-test on isnonog

        fegx(nx+1,:,:) = 0
        fegxy(nx+1,:,:) = 0
*  ---------------------------------------------------------------------
*  compute the energy residuals.
*  ---------------------------------------------------------------------

*  -- total energy residual and equipartition --

      do igsp = 1, ngsp
        do iy = j2, j5
          iy1 = max(0,iy-1)
          do ix = i2, i5
            ix1 = ixm1(ix,iy)

*           Compute thermal equipartition rate with ion-atom fluid
*           ------------------------------------------------------------
            if (nisp >= 2) then   # uses ni(,,2), so must have atoms
                eqpg(ix,iy,igsp) = cftgeqp*(
     .              ng(ix,iy,igsp)*ni(ix,iy,1)*keligig(igsp)
     .              + cftiexclg*ng(ix,iy,igsp)*ni(ix,iy,2)*keligig(igsp))
            endif

*           ------------------------------------------------------------
*                                GAS-ION/ATOM TERMS
*           ------------------------------------------------------------

*           Divergence of gaseous flows & v-grad-P heating
*           ------------------------------------------------------------
            reseg(ix,iy,igsp)= -( fegx(ix,iy,igsp)-fegx(ix1,iy,  igsp)+
     .                            fegy(ix,iy,igsp)-fegy(ix, iy1,igsp) )
     .                                                + segc(ix,iy,igsp)

*           ------------------------------------------------------------
*                                GAS-GAS TERMS
*           ------------------------------------------------------------
            if (igsp.eq.1) then  #..for D0, we should include D+ and D0 in Ti
*               Thermal equipartition coupling of atoms and ions
*               --------------------------------------------------------
c               Should scale with cftiexclg to conserve energy when 
c               transitioning with 0<cftiexclg<1
                reseg(ix,iy,1)= reseg(ix,iy,1) 
     .              + vol(ix,iy)*eqpg(ix,iy,1)*(ti(ix,iy)-tg(ix,iy,1))
                seic(ix,iy) = seic(ix,iy)
     .              - (1.0-cftiexclg)*vol(ix,iy)*eqpg(ix,iy,1)
     .              * (ti(ix,iy)-tg(ix,iy,1))
            else
*               Thermal equipartition coupling of ions and gas
*               --------------------------------------------------------
                reseg(ix,iy,igsp)= reseg(ix,iy,igsp) 
     .              + vol(ix,iy)*eqpg(ix,iy,igsp)*(ti(ix,iy)-tg(ix,iy,igsp))
                seic(ix,iy) = seic(ix,iy)
     .              - vol(ix,iy)*eqpg(ix,iy,igsp)*(ti(ix,iy)-tg(ix,iy,igsp))

*               Thermal equipartition coupling of atoms and gas
*               --------------------------------------------------------
                reseg(ix,iy,igsp) = reseg(ix,iy,igsp)
     .              + cftgeqp*(1.0-cftiexclg)*vol(ix,iy)*kelighg(igsp)
     .              * ng(ix,iy,igsp)*ng(ix,iy,1)*(tg(ix,iy,1)-tg(ix,iy,igsp))
                reseg(ix,iy,1) = reseg(ix,iy,1) 
     .              - cftgeqp*vol(ix,iy)*kelighg(igsp)
     .              * ng(ix,iy,igsp)*ng(ix,iy,1)*(tg(ix,iy,1)-tg(ix,iy,igsp))

              if (ishymol.eq.1 .and. igsp.eq.2) then  #..D2 dissociation

*               Internal ion/atom energy source due to dissociation
*               ----------------------------------------------------
*               IONS
                eiamoldiss(ix,iy,1)=ishymol*ismolcrm*(3/4)*tg(ix,iy,2)*(-psordis(ix,iy,1)/(-2*psordisg(ix,iy,2)))

*               ATOMS
                eiamoldiss(ix,iy,2) = ishymol*ismolcrm * (
     .              (3/4)*tg(ix,iy,2)*(-psordis(ix,iy,2)/(-2*psordisg(ix,iy,2)))
     .          )

*               INCLUDE IN RESIDUALS
                reseg(ix,iy,1) = reseg(ix,iy,1)
     .              + (1-cftiexclg)*eiamoldiss(ix,iy,2)

                seic(ix,iy) = seic(ix,iy)
     .              + eiamoldiss(ix,iy,1) + cftiexclg*eiamoldiss(ix,iy,2)

*               Internal molecular energy sink due to dissociation
*               ----------------------------------------------------
                reseg(ix,iy,2) = reseg(ix,iy,2) + psorg(ix,iy,2)*1.5*tg(ix,iy,2)


*               Drift heating energy source for molecules
*               ----------------------------------------------------
c               The below is adapted from the original implementation,
c               where the first term assumes v_m = 0 when molecular 
c               dissociation is implicitly assumed. The remaining terms
c               are corrections for (v_m - v_a)**2. However, the original
c               implementation seems to mix the poloidal and parallel 
c               velocities arbitrarily, not actually completing the 
c               square.

c               The switches are mixed: cfnidhdis for v_a but cfnidhmol
c               for the molecular terms: use cfnidhmol for all?

c               Only apply drift heating for inertial atoms?
*               ----------------------------------------------------
                if (1.eq.1) then

                    upgcc = 0.5*(up(ix,iy,iigsp)+up(ix1,iy,iigsp))
                    vycc = (cfnidhgy**0.5)*0.5*(vy(ix,iy,iigsp)+vy(ix1,iy,iigsp))
                    v2cc = (cfnidhg2**0.5)*0.5*(v2(ix,iy,iigsp)+v2(ix1,iy,iigsp))

                    reseg(ix,iy,1) = reseg(ix,iy,1) 
     .                  - cfnidh*cfnidhdis*0.5*mg(1)* (upgcc**2 + vycc**2 + v2cc**2)
     .                  * psordis(ix,iy,2)
                    seic(ix,iy) = seic(ix,iy)
     .                  - cftiexclg * cfneut * cfneutsor_ei * cnsor * cfnidhdis
     .                  * 0.5*mg(1)*(upgcc**2 + vycc**2 + v2cc**2) 
     .                  * ( (1-ishymol*ismolcrm)*psordis(ix,iy,2) + ishymol*ismolcrm*psordis(ix,iy,1) )

                    uuxgcc = (cfnidhmol**0.5)*0.5*(uuxg(ix,iy,2)+uuxg(ix1,iy,2))
                    vygcc = (cfnidhmol**0.5)*0.5*(vyg(ix,iy,2)+vyg(ix1,iy,2))
                    v2gcc = 0. #.. molecule v in the tol direction, it seems to be assumed as 0 in neudifpg?

                    reseg(ix,iy,1) = reseg(ix,iy,1) 
     .                  - cfnidhdis*0.5*mg(1)*(uuxgcc**2 + vygcc**2 + v2gcc**2 )*psordis(ix,iy,2)
                    seic(ix,iy) = seic(ix,iy) 
     .                  - cftiexclg*cfnidhdis*0.5*mg(1)*(uuxgcc**2 + vygcc**2 + v2gcc**2 )
     .                  * ( (1-ishymol*ismolcrm)*psordis(ix,iy,2) + ishymol*ismolcrm*psordis(ix,iy,1) )

                    uuxgcc = cfnidhmol*0.25*(uuxg(ix,iy,2)+uuxg(ix1,iy,2))
     .                      *(uuxg(ix,iy,1)+uuxg(ix1,iy,1))
                    vygcc = cfnidhmol*0.25*(vyg(ix,iy,2)+vyg(ix1,iy,2))
     .                      *(vyg(ix,iy,1)+vyg(ix1,iy,1))
                    v2gcc = 0.

                    reseg(ix,iy,1) = reseg(ix,iy,1) 
     .                  + cfnidhdis*mg(1)*(uuxgcc + vygcc + v2gcc)*psordis(ix,iy,2)
                    seic(ix,iy) = seic(ix,iy) 
     .                  + cftiexclg*cfnidhdis*mg(1)*(uuxgcc + vygcc + v2gcc)
     .                  * ( (1-ishymol*ismolcrm)*psordis(ix,iy,2) + ishymol*ismolcrm*psordis(ix,iy,1) )

                else
                    upgcc = 0.5*(up(ix,iy,iigsp)+up(ix1,iy,iigsp))
                    vycc = (cfnidhgy**0.5)*0.5*(vy(ix,iy,iigsp)+vy(ix1,iy,iigsp))
                    v2cc = (cfnidhg2**0.5)*0.5*(v2(ix,iy,iigsp)+v2(ix1,iy,iigsp))

                    uuxgcc = (cfnidhmol**0.5)*0.5*(uuxg(ix,iy,2)+uuxg(ix1,iy,2))
                    vygcc = (cfnidhmol**0.5)*0.5*(vyg(ix,iy,2)+vyg(ix1,iy,2))
                    v2gcc = 0. #.. molecule v in the tol direction, it seems to be assumed as 0 in neudifpg?

                    reseg(ix,iy,1) = reseg(ix,iy,1) 
     .                  - cfnidhdis*0.5*mg(1)
     .                  * ((uuxgcc-upgcc)**2 + (vygcc-vycc)**2 + (v2gcc-v2cc)**2 )
     .                  * psordis(ix,iy,2)

                    seic(ix,iy) = seic(ix,iy) 
     .                  - cftiexclg*cfnidhdis*0.5*mg(1)
     .                  * ((uuxgcc-upgcc)**2 + (vygcc-vycc)**2 + (v2gcc-v2cc)**2 )
     .                  * ( (1-ishymol*ismolcrm)*psordis(ix,iy,2) + ishymol*ismolcrm*psordis(ix,iy,1) )
                endif



              endif
            endif
	    #..zml place holder for neutral-neutral collision,
	    #..    not included above?
          enddo
        enddo
      enddo

      if((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
        do igsp = 1, ngsp
          do iy = j2, j5
            fegx(nxc-1,iy,igsp) = 0. 
            fegx(nxc  ,iy,igsp) = 0. 
            fegx(nxc+1,iy,igsp) = 0. 
          enddo
        enddo
      endif

*  -- Energy transfer to impurity neutrals at tg(,,igsp)
      if (ngsp >= 2) then   # for now, specialized to igsp=2 only
        do ifld = nhsp+1, nisp
          do iy = j2, j5    # iys,iyf limits dont seem to work(?)
            do ix = i2, i5
              #      possible bugs here? resei is replaced with seic
	      #      since resei is not defined before this subroutine called
	      #      more needs to be done here.. e.g. for segc
              seic(ix,iy) =seic(ix,iy) -cftiimpg*1.5*ni(ix,iy,ifld)*
     .                      (nucxi(ix,iy,ifld)+nueli(ix,iy,ifld))*
     .                      (ti(ix,iy) - tg(ix,iy,2))*vol(ix,iy)
              #..zml place holder for things remain to be done for segc
	      #..    to include neutral-impurity ion collisions
            enddo
          enddo
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
c  END subroutine engbalg - THE NEUTRAL GAS ENERGY EQUATION
c-----------------------------------------------------------------------


      SUBROUTINE calc_gas_heatconductivities
      IMPLICIT NONE
      Use(Dim)
      Use(Selec)
      Use(Compla)
      Use(Comgeo)
      Use(UEpar)
      Use(Phyvar)
      Use(Comtra)
      Use(Coefeq)
      Use(Conduc)
      integer igsp, iy, iy1, ix, ix1
      real tgavex, tgavey, noavex, niavex, naavex, noavey, niavey, 
     .  naavey, nuelmolx, qflx, cshx, qshx, nuelmoly, qfly, cshy, qshy

*********************************************************
c ... Gas thermal coefficients, initially for molecules *
*********************************************************
*
c ... Gas thermal conductivity coeffs - from self-collisions
*****************************************************************

      if (nisp >= 2) then  # uses ni(,,2), so must have atoms
       do igsp = 1, ngsp
        do iy = j1, j6
        iy1 = min(iy,ny)
          do ix = i1, i6
            ix1 = ixp1(ix,iy)
            tgavex = max( (tg(ix,iy,igsp)*gx(ix,iy) +
     .                              tg(ix1,iy,igsp)*gx(ix1,iy)) /
     .                             (gx(ix,iy) + gx(ix1,iy)), temin*ev )
            tgavey=max(0.5*(tgy0(ix,iy,igsp)+tgy1(ix,iy,igsp)),temin*ev)
            noavex = ( ng(ix,iy,igsp)*gx(ix,iy) +
     .                                   ng(ix1,iy,igsp)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            niavex = ( ni(ix,iy,1)*gx(ix,iy) +
     .                                   ni(ix1,iy,1)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            naavex = ( ni(ix,iy,2)*gx(ix,iy) +
     .                                   ni(ix1,iy,2)*gx(ix1,iy)) /
     .                                   (gx(ix,iy) + gx(ix1,iy))
            noavey = 0.5*(ngy0(ix,iy1,igsp) + ngy1(ix,iy1,igsp))
            niavey = 0.5*(niy0(ix,iy1,1) + niy1(ix,iy1,1))
            naavey = 0.5*(niy0(ix,iy1,2) + niy1(ix,iy1,2))
            nuelmolx = noavex*kelhmhm + niavex*kelhmhg + 
     .                 naavex*kelhmhg
            qflx = flalftmx*sqrt(tgavex/mg(igsp))*noavex*tgavex
            cshx = cftgcond*noavex*tgavex/(mg(igsp)*nuelmolx)  #assume K not fcn Tg
            qshx = cshx * (tg(ix,iy,igsp)-tg(ix1,iy,igsp)) * gxf(ix,iy)
            hcxg(ix,iy,igsp) = cshx / 
     .                     (1.+ (abs(qshx/qflx))**flgamtg)**(1./flgamtg)
            hcxg(ix,iy,igsp)=(1.-cfhcxgc(igsp))*hcxg(ix,iy,igsp)+
     .                     cfhcxgc(igsp)*noavex*kxg_use(ix,iy,igsp)
c..   Now radial direction
            nuelmoly = noavey*kelhmhm + niavey*kelhmhg + 
     .                 naavey*kelhmhg
            qfly = flalftmy*sqrt(tgavey/mg(igsp))*noavey*tgavey
            cshy = cftgcond*noavey*tgavey/(mg(igsp)*nuelmoly)  #assume Kel_s not fcn Tg
            qshy = cshy*(tgy0(ix,iy1,igsp)-tgy1(ix,iy1,igsp))/
     .                                                  dynog(ix,iy)
            hcyg(ix,iy,igsp) = cshy / 
     .                     (1 + (abs(qshy/qfly))**flgamtg)**(1./flgamtg)
            hcyg(ix,iy,igsp)=(1-cfhcygc(igsp))*hcyg(ix,iy,igsp)+
     .                     cfhcygc(igsp)*noavey*kyg_use(ix,iy,igsp)
          enddo
        enddo
        if (igsp.eq.1 .and. isupgon(igsp).eq.1) then 
          hcxg(:,:,igsp) = hcxn(:,:)
          hcyg(:,:,igsp) = hcyn(:,:)
        endif
       enddo
      endif

      END SUBROUTINE calc_gas_heatconductivities


