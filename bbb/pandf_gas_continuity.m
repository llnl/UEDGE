c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"



c ======================================================================
c  Here we do the neutral gas diffusion model, if isngon=1.
c  The diffusion is flux limited using the thermal flux.
c
c  If we are solving for the neutral parallel momentum (isupgon=1)
c  we only want the perpendicular flux of neutrals
c  and the corresponding perp velocity.
c  Franck-Condon neutrals complicate the picture with parallel momentum.
c  Not yet resolved in the coding.
c ======================================================================
c
      subroutine neudif

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,ffyi,ffyo
      integer iy1, methgx, methgy, iy2, jx
      logical isxyfl
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
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
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
	  
*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10
	  
      do igsp = 1, ngsp

c.... First the flux in the x-direction

      do iy = j4, j8
         do ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
            qfl = flalfgxa(ix,igsp) * sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
            csh = (1-isgasdc) * cdifg(igsp) * sx(ix,iy) * gxf(ix,iy) *
     .                ave( stretcx(ix,iy)*vtn**2/nu1,
     .                     stretcx(ix2,iy)*vtnp**2/nu2 ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng +
     .                        rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                          0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
            qtgf = cngfx(igsp) * fgtdx(ix) * sx(ix,iy) *
     .            ave( stretcx(ix,iy)*gx(ix,iy)/nu1 ,
     .                 stretcx(ix2,iy)*gx(ix2,iy)/nu2 )
     .                     * (vtn**2 - vtnp**2)
            if (isupgon(igsp).eq.1) then # only 2->pol project.;up has rest
               csh = csh*(1 - rrv(ix,iy)**2)
               qtgf = qtgf*(1 - rrv(ix,iy)**2)
            endif
            vygtan(ix,iy,igsp) = 0.  # vygtan is grad(T) rad vel on x-face
            if (isnonog .eq. 1 .and. iy .le. ny) then
                  grdnv =(    ( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                       -( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               vygtan(ix,iy,igsp) = exp( 0.5*
     .                     (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                      ( cngfx(igsp) / (mg(igsp)*0.5*(nu1+nu2)) ) *
     .                                     ( grdnv/cos(angfx(ix,iy)) - 
     .                       (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))
     .                                                 * gxf(ix,iy) ) 
             if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
               vygtan(ix,iy,igsp) = 0.
             endif
             if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
               vygtan(ix,iy,igsp) = 0.
             endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - vygtan(ix,iy,igsp)*sx(ix,iy)
            if (isupgon(igsp) .eq. 1) then
              qtgf = qtgf + rrv(ix,iy)*up(ix,iy,iigsp)*sx(ix,iy)
            endif
            nconv = 2.0*(ng(ix,iy,igsp)*ng(ix2,iy,igsp)) /
     .                  (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
c...   Use upwind for "convective" grad T term if methgx .ne. 2
            if(methgx.ne.2) nconv =
     .                         ng(ix ,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ng(ix2,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ng(ix,iy,igsp)-ng(ix2,iy,igsp)) + qtgf * nconv
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  bndry face - should not matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floxg
            floxg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = 
     .              floxg(ix,iy) + cngflox(1)*sx(ix,iy)*uu(ix,iy,1)

        end do
         conxg(nx+1,iy) = 0
      end do

c.... Now the flux in the y-direction

      do iy = j1, j5
         do ix = i4, i8
	    t0 = max(tg(ix,iy,igsp),temin*ev)
	    t1 = max(tg(ix,iy+1,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix,iy+1,igsp) + vtnp/lgmax(igsp)
            qfl = flalfgya(iy,igsp) * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .                              ( ngy0(ix,iy,igsp)*gy(ix,iy) + 
     .                                ngy1(ix,iy,igsp)*gy(ix,iy+1) ) / 
     .                                     (8*(gy(ix,iy)+gy(ix,iy+1)))
             csh = (1-isgasdc) * cdifg(igsp) *sy(ix,iy)/(dynog(ix,iy)) *
     .                                  ave(vtn**2/nu1, vtnp**2/nu2) +
     .            isgasdc * sy(ix,iy) * gyf(ix,iy) * difcng +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*(1/gyf(ix,iy))*
     .                       0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))
            qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                            ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2)
     .                                      * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0 ) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) )
              endif
              qtgf = cngfy(igsp) * fgtdy(iy)* sy(ix,iy) * 
     .                           ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ngy0(ix,iy,igsp)-ngy1(ix,iy,igsp)) + qtgf*nconv
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floyg
	    floyg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) =  
     .             floyg(ix,iy)+cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)

        end do
      end do

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             ng(0:nx+1,0:ny+1,igsp),fngx(0:nx+1,0:ny+1,igsp),
     .             fngy(0:nx+1,0:ny+1,igsp),0,methg)

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then

         do iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
	       t0 = max(tg(ix ,iy,igsp),temin*ev)
	       t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn =  sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               nu1 = nuix(ix ,iy,igsp) + vtn/lgmax(igsp)
               nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
ccc            MER: Set flag to apply xy flux limit except at target plates
               isxyfl = .true.
               do jx = 1, nxpt
                  if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                 (ix==ixrb(jx).and.ixmxbcl==1) ) isxyfl = .false.
               enddo
               if (methgx .eq. 6) then  # log interpolation
               grdnv =(   ( fym (ix,iy,1)*log(ng(ix2,iy1 ,igsp)) +
     .                      fy0 (ix,iy,1)*log(ng(ix2,iy  ,igsp)) +
     .                      fyp (ix,iy,1)*log(ng(ix2,iy+1,igsp)) + 
     .                      fymx(ix,iy,1)*log(ng(ix ,iy1 ,igsp)) +
     .                      fypx(ix,iy,1)*log(ng(ix, iy+1,igsp)) )
     .                   -( fym (ix,iy,0)*log(ng(ix ,iy1 ,igsp)) +
     .                      fy0 (ix,iy,0)*log(ng(ix ,iy  ,igsp)) +
     .                      fyp (ix,iy,0)*log(ng(ix ,iy+1,igsp)) +
     .                      fymx(ix,iy,0)*log(ng(ix4,iy1 ,igsp)) + 
     .                      fypx(ix,iy,0)*log(ng(ix6,iy+1,igsp)) ) )/ 
     .                                                  dxnog(ix,iy)
               elseif (methgx .eq. 7) then  # inverse interpolation
               grdnv =( 1/(fym (ix,iy,1)/ng(ix2,iy1 ,igsp) + 
     .                     fy0 (ix,iy,1)/ng(ix2,iy  ,igsp) +
     .                     fyp (ix,iy,1)/ng(ix2,iy+1,igsp) + 
     .                     fymx(ix,iy,1)/ng(ix ,iy1 ,igsp) +
     .                     fypx(ix,iy,1)/ng(ix, iy+1,igsp))
     .                - 1/(fym (ix,iy,0)/ng(ix ,iy1 ,igsp) +
     .                     fy0 (ix,iy,0)/ng(ix ,iy  ,igsp) +
     .                     fyp (ix,iy,0)/ng(ix ,iy+1,igsp) +
     .                     fymx(ix,iy,0)/ng(ix4,iy1 ,igsp) + 
     .                     fypx(ix,iy,0)/ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               else                   # linear interpolation
               grdnv =( (fym (ix,iy,1)*ng(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*ng(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*ng(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*ng(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*ng(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*ng(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*ng(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*ng(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*ng(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               endif
               difgx2 = ave( tg(ix ,iy,igsp)/nu1,
     .                       tg(ix2,iy,igsp)/nu2 )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
              if (methgx .eq. 6) then
               fngxy(ix,iy,igsp) =  exp( 0.5*
     .                     (log(ng(ix2,iy,igsp))+log(ng(ix,iy,igsp))) )*
     .                               difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                     (log(ng(ix2,iy,igsp)) - log(ng(ix,iy,igsp)))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              else
               fngxy(ix,iy,igsp) = difgx2*( grdnv/cos(angfx(ix,iy)) -
     .                             (ng(ix2,iy,igsp) - ng(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              endif
c...  Now flux limit with flalfgxy
               t0 = max(tg(ix,iy,igsp),temin*ev)
               t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn = sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               qfl = flalfgxya(ix,igsp)*sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
ccc   MER NOTE:  no xy flux limit for ix=0 or ix=nx in original code
               if (isxyfl) then
                  fngxy(ix,iy,igsp) = fngxy(ix,iy,igsp) /
     .                           sqrt( 1 + (fngxy(ix,iy,igsp)/qfl)**2 )
               endif
            end do
        end do

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Flux-limit the total poloidal flux here
            do iy = j4, j8
               do ix = i1, i5
                  ix2 = ixp1(ix,iy)
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)-fngxy(ix,iy,igsp)
                  t0 = max(tg(ix,iy,igsp),temin*ev)
                  t1 = max(tg(ix2,iy,igsp),temin*ev)
                  vtn = sqrt( t0/mg(igsp) )
                  vtnp = sqrt( t1/mg(igsp) )
                  qfl = flalfgnx * sx(ix,iy)*(vtn + vtnp)*rt8opi*
     .                             (ng(ix,iy,igsp)+ng(ix2,iy,igsp))/16
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/
     .                              sqrt(1 + (fngx(ix,iy,igsp)/qfl)**2)
c ...          adjust fluxes to prevent pumpout
                    fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/( 1 - 2*nlimgx +
     .                          nlimgx*(ng(ix2,iy,igsp)/ng(ix,iy,igsp) +
     .                                  ng(ix,iy,igsp)/ng(ix2,iy,igsp)) )
               enddo
            enddo
c ...   adjust y-fluxes to prevent pumpout
            do iy = j1, j5    # same loop ranges as for fngy in fd2tra
               do ix = i4, i8
                     fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/( 1 - 2*nlimgy +
     .                         nlimgy*(ng(ix,iy+1,igsp)/ng(ix,iy,igsp)+
     .                                 ng(ix,iy,igsp)/ng(ix,iy+1,igsp)) )
	            t0 = max(tg(ix,iy,igsp),temin*ev)
	            t1 = max(tg(ix,iy+1,igsp),temin*ev)
                    vtn = sqrt( t0/mg(igsp) )
                    vtnp = sqrt( t1/mg(igsp) )
                    qfl = flalfgny*sy(ix,iy)*(vtn+vtnp)*rt8opi*
     .                        (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp))/16
                    fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/
     .                            sqrt(1 + (fngy(ix,iy,igsp)/qfl)**2)
               enddo
            enddo     

      endif
c...  Finished with nonorthogonal mesh part

c ... Calculate the neutral flow velocity from v = flux/ng; these are
c ... diagnostic if isupgon=0, but used for uu and vy of the inertial
c ... gas if isupgon=1.  However, even when isupgon=1, the particle
c ... fluxes are fngx -> fnix and fngy -> fniy in pandf, i.e., methg
c ... determines the differencing for the inertial particle fluxes, not 
c ... methn
      do iy = j1, j5
         do ix = i1,i5
            ix1 = ixp1(ix,iy)
            if (1.-rrv(ix,iy) > 1.e-4) then #combine binormal & par components
              uug(ix,iy,igsp) = fngx(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                      *sx(ix,iy) )
            else   # binormal component negligable small
               uug(ix,iy,igsp) = up(ix,iy,iigsp)
            endif
            vyg(ix,iy,igsp) = fngy(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
     .                                                      *sy(ix,iy) )
c --------------- transfer inertial gas velocities to neutral ion species
            if (isupgon(igsp).eq.1) then
               vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
cc               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            end if
         enddo
cfw ----   If doing only the outer half we want this boundary condition:
         if (iy.le.iysptrx2(1) .and. isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp) = 0
      enddo

c **- loop for uu just as in the previous version - needed for correct Jac?
      if (isupgon(igsp) .eq. 1) then
         do iy = j4, j6
            do ix = i1, i6
               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            enddo
         enddo
      endif

c.... Calculate the residual for the gas equation for diffusive neutral case

      if (isupgon(igsp).eq.0) then
         do iy = j2, j5
	      if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
	      fngx(nxc-1,iy,igsp) = 0.
	      fngx(nxc,  iy,igsp) = 0.
	      fngx(nxc+1,iy,igsp) = 0.
	    endif
            if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
            if (nxpt==2) fngx(ixrb(1)+1,iy,igsp)=0.
       	    do ix = i2, i5
               ix1 = ixm1(ix,iy)
               resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .                         psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .                       + volpsorg(ix,iy,igsp)                    
     .                       - fngx(ix,iy,igsp) + fngx(ix1,iy  ,igsp)
     .             - fluxfacy*(fngy(ix,iy,igsp) - fngy(ix ,iy-1,igsp))
     .                       + psgov_use(ix,iy,igsp)*vol(ix,iy)
               if (igsp.eq.1 .and. ishymol.eq.1) resng(ix,iy,igsp) =
     .                  resng(ix,iy,igsp)+psordis(ix,iy,2)
            end do
        end do
      endif

       end do   # end of igsp loop from the beginning of subroutine

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                         nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudif
c --------------------------------------------------------------------------

c ======================================================================
c  Below neudifpg similar to neudif, except that the neutral pressure, ng*tg, 
c  is the dependent variable differenced, rather than ng and tg separately
c  Here we do the neutral gas diffusion model, if isngon=1.
c  The diffusion is flux limited using the thermal flux.
c
c  If we are solving for the neutral parallel momentum (isupgon=1)
c  we only want the perpendicular flux of neutrals
c  and the corresponding perp velocity.
c  Franck-Condon neutrals complicate the picture with parallel momentum.
c  Not yet resolved in the coding.
c ======================================================================
c
      subroutine neudifpg

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,nu1,nu2,tgf,ffyi,ffyo
      real tsngxlog,tsngylog,tsngfd2,tsngfxy
      real dndym1,dndy0,dndyp1,d2ndy20,d2ndy2p1,d3ndy3
      real dndxm1,dndx0,dndxp1,d2ndx20,d2ndx2p1,d3ndx3
      real flalfgx_adj, flalfgy_adj, flalfgxy_adj, ngxface, ngyface
      integer iy1, methgx, methgy, iy2, jx, jfld, ifld
      integer iym1,iyp1,iyp2,ixm1b,ixp1b,ixp2b
      logical isxyfl
      # Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a,tick,tock
      external tick, tock

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
                    # xlinc,xrinc,yinc,ixm1,ixp1
      Use(Comgeo)   # vol, gx, gy, sx ,xy
      Use(Noggeo)   # fym,fy0,fyp,fymx,fypx,angfx
      Use(Compla)   # mi, zi, ni, uu, up, v2, v2ce, vygtan, mg
      Use(Comflo)   # fngx,fngy,fngxy,fnix,fniy
      Use(Conduc)   # nuiz, nucx, nuix
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
      Use(Timing)   # ttngxlog,ttngylog,ttngfd2,ttngfxy

      Use(Ext_neutrals) # get_neutral_moments, ...
      Use(MCN_dim)      # ngsp, ...
      Use(MCN_sources)  # cfneut_sng, cfneutdiv_fng, ... mcfngx, mcfngy, ...
      Use(Interp)		# ngs, tgs 
      Use(Bfield)   # rbfbt 

*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10
	  
      do igsp = 1, ngsp

c *********************************************
c.... First the flux in the x-direction
c *********************************************

c ..Timing;initialize 
      if(istimingon==1) tsngxlog = tick()

      do 888 iy = j4omp, j8omp
         do 887 ix = i1momp, i5pomp
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            ngxface = 0.5*(ng(ix,iy,igsp)+ng(ix2,iy,igsp))
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
	    tgf = 0.5*(tg(ix,iy,igsp)+tg(ix2,iy,igsp))
            flalfgx_adj = flalfgxa(ix,igsp)*( 1. +
     .                    (cflbg*ngbackg(igsp)/ngxface)**inflbg )
            qfl = flalfgx_adj * sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
            csh = (1-isgasdc) * cdifg(igsp) * sx(ix,iy) * gxf(ix,iy) *
     .                (1/mg(igsp))* ave( 1./nu1,1./nu2 ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng / tgf +
     .                        rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                    0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))/tgf
            qtgf = alftng * fgtdx(ix) * sx(ix,iy) *
     .            ave( gx(ix,iy)/nu1 ,gx(ix2,iy)/nu2 )
     .                     * (vtn**2 - vtnp**2)
            if (isupgon(igsp).eq.1) then # only 2->pol project.;up has rest
               csh = csh*(1 - rrv(ix,iy)**2)
               qtgf = qtgf*(1 - rrv(ix,iy)**2)
            endif
            vygtan(ix,iy,igsp) = 0.  # vygtan is grad(T) rad vel on x-face
                                     # vygtan) only from thermal force
            if (isnonog .eq. 1 .and. iy .le. ny) then
                  grdnv =(    ( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                       -( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               vygtan(ix,iy,igsp) = exp( 0.5*
     .                     (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                      ( alftng / (mg(igsp)*0.5*(nu1+nu2)) ) *
     .                                     ( grdnv/cos(angfx(ix,iy)) - 
     .                       (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))
     .                                                 * gxf(ix,iy) ) 
             if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
               vygtan(ix,iy,igsp) = 0.
             endif
             if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
               vygtan(ix,iy,igsp) = 0.
             endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - vygtan(ix,iy,igsp)*sx(ix,iy)
            if (isupgon(igsp) .eq. 1) then
              qtgf = qtgf + rrv(ix,iy)*up(ix,iy,iigsp)*sx(ix,iy)
            endif
            nconv = 2.0*(ng(ix,iy,igsp)*ng(ix2,iy,igsp)) /
     .                  (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
c...   Use upwind for "convective" grad T term if methgx .ne. 2
            if(methgx.ne.2) nconv =
     .                         ng(ix ,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ng(ix2,iy,igsp)*0.5*(1-sign(1.,qtgf))

            qsh = csh * (pg(ix,iy,igsp)-pg(ix2,iy,igsp)) + qtgf * nconv
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  boundary face - shouldnot matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  themal force temperature gradient term is included in floxg
	    floxg(ix,iy) = (qtgf/tgf) / (1 + qr**flgamg)**(1/flgamg)
c...  now add the ion convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = floxg(ix,iy) +
     .                            cngflox(1)*sx(ix,iy)*uu(ix,iy,1)/tgf
c...  For one impurity, add convect vel for elastic scattering with ions
         if(igsp == 2 .and. ishymol == 0) then #Caution; need to weight uu
           do ifld = nhsp+1, nisp
             floxg(ix,iy) = floxg(ix,iy) +
     .                cngniflox(ifld,igsp)*sx(ix,iy)*uu(ix,iy,ifld)/tgf
           enddo
         endif  

  887    continue
         conxg(nx+1,iy) = 0
  888 continue

c ..Timing; add info if timing is on
      if(istimingon==1) ttngxlog = ttngxlog + tock(tsngxlog)

c *******************************************************
c.... Now the flux in the y-direction
c *******************************************************

c ..Timing; initiate time for y-direction calc
      if(istimingon==1) tsngylog = tick()

      do iy = j1omp1, j5omp
         do ix = i4momp, i8pomp
            ngyface = 0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
	    t0 = max(tg(ix,iy,igsp),tgmin*ev)
	    t1 = max(tg(ix,iy+1,igsp),tgmin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            nu1 = nuix(ix,iy,igsp) + vtn/lgmax(igsp)
            nu2 = nuix(ix,iy+1,igsp) + vtnp/lgmax(igsp)
	    tgf = 0.5*(tg(ix,iy,igsp)+tg(ix,iy+1,igsp))
            flalfgy_adj = flalfgya(iy,igsp)*( 1. +
     .                   (cflbg*ngbackg(igsp)/ngyface)**inflbg )
            qfl = flalfgy_adj * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .                              ( ngy0(ix,iy,igsp)*gy(ix,iy) + 
     .                                ngy1(ix,iy,igsp)*gy(ix,iy+1) ) / 
     .                                     (8*(gy(ix,iy)+gy(ix,iy+1)))
            if (iy==0) then  #at bdry, ng ave to avoid ng->0 prob
              qfl = flalfgy_adj * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .              (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) / 8.
            endif
            csh = (1-isgasdc) * (cdifg(igsp) *sy(ix,iy)/dynog(ix,iy)) *
     .                          (1/mg(igsp))* ave(1./nu1, 1./nu2) +
     .            isgasdc * sy(ix,iy) * difcng /(dynog(ix,iy)*tgf) +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*dynog(ix,iy)*
     .                     0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))/tgf

            qtgf = alftng * fgtdy(iy) * sy(ix,iy) * 
     .                            ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2)
     .                                      * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0 ) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) )
              endif
              qtgf = alftng * fgtdy(iy)* sy(ix,iy) * 
     .                           ave(gy(ix,iy)/nu1, gy(ix,iy+1)/nu2) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (pgy0(ix,iy,igsp)-pgy1(ix,iy,igsp)) + qtgf*nconv
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  thermal force temperature gradient term is included in floyg
	    floyg(ix,iy) = (qtgf/tgf) / (1 + qr**flgamg)**(1/flgamg)
c...  now add ion convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) = floyg(ix,iy) +
     .                            cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)/tgf
c...  For impurities, add convect vel for elastic scattering with ions
         if(igsp == 2 .and. ishymol == 0) then #Caution; need to weight uu 
           do ifld = nhsp+1, nisp
             floyg(ix,iy) = floyg(ix,iy) +
     .                cngnifloy(ifld,igsp)*sy(ix,iy)*vy(ix,iy,ifld)/tgf
           enddo
         endif
        end do
      end do

c ..Timing; add increment if timing is on
      if(istimingon==1) ttngylog = ttngylog + tock(tsngylog)

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------
c ..Timing
      if(istimingon==1) tsngfd2 = tick()

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             pg(0:nx+1,0:ny+1,igsp),fngx(0:nx+1,0:ny+1,igsp),
     .             fngy(0:nx+1,0:ny+1,igsp),0,methg)
c ..Timing
      if(istimingon==1) ttngfd2 = ttngfd2 + tock(tsngfd2)

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then
c ..Timing
      if(istimingon==1) tsngfxy = tick()

         do iy = j1omp1, min(j6omp, ny)
            iy1 = max(iy-1,0)
            do ix = i1momp, min(i6pomp, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
	       t0 = max(tg(ix ,iy,igsp),temin*ev)
	       t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn =  sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               nu1 = nuix(ix ,iy,igsp) + vtn/lgmax(igsp)
               nu2 = nuix(ix2,iy,igsp) + vtnp/lgmax(igsp)
ccc            MER: Set flag to apply xy flux limit except at target plates
               isxyfl = .true.
               do jx = 1, nxpt
                 if( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .               (ix==ixrb(jx).and.ixmxbcl==1) ) isxyfl = .false.
               enddo
               if (methgx .eq. 6) then  # log interpolation
               grdnv =(   ( fym (ix,iy,1)*log(pg(ix2,iy1 ,igsp)) + 
     .                      fy0 (ix,iy,1)*log(pg(ix2,iy  ,igsp)) +
     .                      fyp (ix,iy,1)*log(pg(ix2,iy+1,igsp)) + 
     .                      fymx(ix,iy,1)*log(pg(ix ,iy1 ,igsp)) +
     .                      fypx(ix,iy,1)*log(pg(ix, iy+1,igsp)) )
     .                   -( fym (ix,iy,0)*log(pg(ix ,iy1 ,igsp)) +
     .                      fy0 (ix,iy,0)*log(pg(ix ,iy  ,igsp)) +
     .                      fyp (ix,iy,0)*log(pg(ix ,iy+1,igsp)) +
     .                      fymx(ix,iy,0)*log(pg(ix4,iy1 ,igsp)) + 
     .                      fypx(ix,iy,0)*log(pg(ix6,iy+1,igsp)) ) )/ 
     .                                                  dxnog(ix,iy)
               elseif (methgx .eq. 7) then # inverse interpolation
               grdnv =( 1/(fym (ix,iy,1)/pg(ix2,iy1 ,igsp) + 
     .                     fy0 (ix,iy,1)/pg(ix2,iy  ,igsp) +
     .                     fyp (ix,iy,1)/pg(ix2,iy+1,igsp) + 
     .                     fymx(ix,iy,1)/pg(ix ,iy1 ,igsp) +
     .                     fypx(ix,iy,1)/pg(ix, iy+1,igsp))
     .                - 1/(fym (ix,iy,0)/pg(ix ,iy1 ,igsp) +
     .                     fy0 (ix,iy,0)/pg(ix ,iy  ,igsp) +
     .                     fyp (ix,iy,0)/pg(ix ,iy+1,igsp) +
     .                     fymx(ix,iy,0)/pg(ix4,iy1 ,igsp) + 
     .                     fypx(ix,iy,0)/pg(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               else                   # linear interpolation
               grdnv =( (fym (ix,iy,1)*pg(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*pg(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*pg(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*pg(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*pg(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*pg(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*pg(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*pg(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*pg(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*pg(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               endif
               difgx2 = ave( 1./nu1,
     .                       1./nu2 )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
              if (methgx .eq. 6) then
               fngxy(ix,iy,igsp) =  exp( 0.5*
     .                     (log(pg(ix2,iy,igsp))+log(pg(ix,iy,igsp))) )*
     .                               difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                     (log(pg(ix2,iy,igsp)) - log(pg(ix,iy,igsp)))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              else
               fngxy(ix,iy,igsp) = difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                             (pg(ix2,iy,igsp) - pg(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
              endif
c...  Now flux limit with flalfgxy
               ngxface = 0.5*(ng(ix,iy,igsp)+ng(ix2,iy,igsp))
               t0 = max(tg(ix,iy,igsp),temin*ev)
               t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn = sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               flalfgxy_adj = flalfgxya(ix,igsp)*( 1. +
     .                     (cflbg*ngbackg(igsp)/ngxface)**inflbg )
               qfl = flalfgxy_adj*sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
ccc   MER NOTE:  no xy flux limit for ix=0 or ix=nx in original code
               if (isxyfl) then
                  fngxy(ix,iy,igsp) = fngxy(ix,iy,igsp) /
     .                           sqrt( 1 + (fngxy(ix,iy,igsp)/qfl)**2 )
               endif
            end do
        end do
c ..Timing
      if(istimingon==1) ttngfxy = ttngfxy + tock(tsngfxy)

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Flux-limit the total poloidal flux here
            do iy = j4omp, j8omp
               do ix = i1momp, i5pomp
                  ix2 = ixp1(ix,iy)
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)-fngxy(ix,iy,igsp)
                  t0 = max(tg(ix,iy,igsp),temin*ev)
                  t1 = max(tg(ix2,iy,igsp),temin*ev)
                  vtn = sqrt( t0/mg(igsp) )
                  vtnp = sqrt( t1/mg(igsp) )
                  qfl = flalfgnx * sx(ix,iy)*(vtn + vtnp)*rt8opi*
     .                             (ng(ix,iy,igsp)+ng(ix2,iy,igsp))/16
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/
     .                              sqrt(1 + (fngx(ix,iy,igsp)/qfl)**2)
c ...          adjust fluxes to prevent pumpout
                    fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/( 1 - 2*nlimgx +
     .                          nlimgx*(ng(ix2,iy,igsp)/ng(ix,iy,igsp) +
     .                                  ng(ix,iy,igsp)/ng(ix2,iy,igsp)) )
               enddo
            enddo
c ...   adjust y-fluxes to prevent pumpout
            do iy = j1, j5    # same loop ranges as for fngy in fd2tra
               do ix = i4omp, i8omp
                     fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/( 1-2*nlimgy +
     .                         nlimgy*(ng(ix,iy+1,igsp)/ng(ix,iy,igsp)+
     .                                ng(ix,iy,igsp)/ng(ix,iy+1,igsp)) )
	            t0 = max(tg(ix,iy,igsp),temin*ev)
	            t1 = max(tg(ix,iy+1,igsp),temin*ev)
                    vtn = sqrt( t0/mg(igsp) )
                    vtnp = sqrt( t1/mg(igsp) )
                    qfl = flalfgny*sy(ix,iy)*(vtn+vtnp)*rt8opi*
     .                        (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp))/16
                    fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/
     .                            sqrt(1 + (fngy(ix,iy,igsp)/qfl)**2)
               enddo
            enddo     

      endif
c...  Finished with nonorthogonal mesh part

c...  Add 4th order radial diffusion op; damp grid-scale oscill
        if (abs(difgy4order(igsp)) > 1.e-50 .and. isngon(igsp)==1) then
          do iy = j2p, j5m   #limits to range iy=1:ny-1 for fngy4ord
            iym1 = max(iy-1,0)
            iyp1 = min(iy+1,ny+1)
            iyp2 = min(iy+2,ny+1)
            do ix = i4omp, i8omp
              dndym1 = (ng(ix,iy,igsp)-ng(ix,iym1,igsp))*gyf(ix,iym1)
              dndy0 = (ng(ix,iyp1,igsp)-ng(ix,iy,igsp))*gyf(ix,iy)
              dndyp1 = (ng(ix,iyp2,igsp)-ng(ix,iyp1,igsp))*gyf(ix,iyp1)
              d2ndy20 = (dndy0 - dndym1)*gy(ix,iy)
              d2ndy2p1 = (dndyp1 - dndy0)*gy(ix,iyp1)
              d3ndy3 = (d2ndy2p1 - d2ndy20)*gyf(ix,iy)
              fngy4ord(ix,iy,igsp) = difgy4order(igsp)*d3ndy3*sy(ix,iy)/
     .                                                  gyf(ix,iy)**2
              fngy(ix,iy,igsp) = fngy(ix,iy,igsp) + fngy4ord(ix,iy,igsp)
            enddo
          enddo
        endif

c...  Add 4th order poloidal diffusion op; damp grid-scale oscill
C...  NOTE: PRESENTLY ONLY CODED FOR SIMPLY-CONNECTED DOMAIN
        if (abs(difgx4order(igsp)) > 1.e-50 .and. isngon(igsp)==1) then
          do ix = i2p, i5m   #limits to range ix=1:nx-1 for fngx4ord
            ixm1b = max(ix-1,0)
            ixp1b = min(ix+1,nx+1)
            ixp2b = min(ix+2,nx+1)
            do iy = j4omp, j8omp
              dndxm1 = (ng(ix,iy,igsp)-ng(ixm1b,iy,igsp))*gxf(ixm1b,iy)
              dndx0 = (ng(ixp1b,iy,igsp)-ng(ix,iy,igsp))*gxf(ix,iy)
              dndxp1 = (ng(ixp2b,iy,igsp)-ng(ixp1b,iy,igsp))*gxf(ixp1b,iy)
              d2ndx20 = (dndx0 - dndxm1)*gx(ix,iy)
              d2ndx2p1 = (dndxp1 - dndx0)*gx(ixp1b,iy)
              d3ndx3 = (d2ndx2p1 - d2ndx20)*gxf(ix,iy)
              fngx4ord(ix,iy,igsp) = difgx4order(igsp)*d3ndx3*sx(ix,iy)/
     .                                                  gxf(ix,iy)**2
              fngx(ix,iy,igsp) = fngx(ix,iy,igsp) + fngx4ord(ix,iy,igsp)
            enddo
          enddo
        endif

c ... Calculate the neutral flow velocity from v = flux/ng; these are
c ... diagnostic if isupgon=0, but used for uu and vy of the inertial
c ... gas if isupgon=1.  However, even when isupgon=1, the particle
c ... fluxes are fngx -> fnix and fngy -> fngy in pandf, i.e., methg
c ... determines the differencing for the inertial particle fluxes, not 
c ... methn
      do iy = j1omp1, j5omp
        do ix = i1momp,i5pomp
          ix1 = ixp1(ix,iy)
          if (1.-rrv(ix,iy) > 1.e-4 .or. isupgon(igsp)==0) then 
                           #combine binormal/par comps or x only diffusive
             uug(ix,iy,igsp) = fngx(ix,iy,igsp) / (
     .                      0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                    *sx(ix,iy) )
c ...    remove nonorthogonal vy component from uug for seic contrib
             uuxg(ix,iy,igsp) = (fngx(ix,iy,igsp)+fngxy(ix,iy,igsp))/ (
     .                      0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                    *sx(ix,iy) )
          else   # binormal component negligable small
             uug(ix,iy,igsp) = up(ix,iy,iigsp)
          endif
          vyg(ix,iy,igsp) = fngy(ix,iy,igsp) / (
     .                      0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
     .                                                      *sy(ix,iy) )
c --------------- transfer inertial gas velocities to neutral ion species
          if (isupgon(igsp).eq.1) then
             vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
          endif
        enddo
cfw ----   If doing only the outer half we want this boundary condition:
        if(iy.le.iysptrx2(1).and.isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp)=0
      enddo

c **- loop for uu just as in the previous version - needed for correct Jac?
      if (isupgon(igsp) .eq. 1) then
         do iy = j4omp, j6omp
            do ix = i1momp, i6pomp
               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
               v2(ix,iy,iigsp) = ( uuxg(ix,iy,igsp) 
     .                            - up(ix,iy,iigsp)*rrv(ix,iy) )
     .                       /(rbfbt(ix,iy) + rbfbt(ixp1(ix,iy),iy))*2.
            enddo
         enddo
      endif

      end do   # end of igsp loop from the beginning of subroutine
      return
      end
c
c --------------------------------------------------------------------------
c END subroutine neudifpg
c --------------------------------------------------------------------------

      SUBROUTINE calc_gas_continuity_residuals
      IMPLICIT NONE
      Use(UEpar)
      Use(Selec)
      Use(Share)
      Use(Comflo)
      Use(Compla)
      Use(Dim)
      Use(Indices_domain_dcl)
      Use(Xpoint_indices)
      Use(Rhsides)
      Use(Coefeq)
      Use(Volsrc)
      Use(Comgeo)
      Use(MCN_sources)
      Use(Comtra)
      Use(Ext_neutrals)
      Use(Conduc)
      Use(Rccoef)
      integer igsp, iy, ix, ix1, jfld
      real tnuiz, lmfp, ngnot

c.... Calculate the residual for the gas equation for diffusive neutral case

      do igsp = 1, ngsp
      if (isupgon(igsp).eq.0) then
        do iy = j2omp, j5omp
          if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then 
	     fngx(nxc-1,iy,igsp) = 0.
	     fngx(nxc,  iy,igsp) = 0.
	     fngx(nxc+1,iy,igsp) = 0.
	  endif
          if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
          if (nxpt==2.and.ixmxbcl==1) fngx(ixrb(1)+1,iy,igsp)=0.
          do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)

c ... 2016/09/16 IJ: coding to blend MC neutral flux !!! here ***
c ... is it correct to use ng instead of ni??? i.e. will ng enter jacobian?
            resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .               psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .             + volpsorg(ix,iy,igsp)
     .             + psgov_use(ix,iy,igsp)*vol(ix,iy)
            if (igsp.eq.1 .and. ishymol.eq.1)
     .          resng(ix,iy,igsp) = resng(ix,iy,igsp)+psordis(ix,iy,2)
            resng(ix,iy,igsp) = resng(ix,iy,igsp) - cfneutdiv*
     .          cfneutdiv_fng*((fngx(ix,iy,igsp) - fngx(ix1,iy, igsp)) +
     .          fluxfacy*(fngy(ix,iy,igsp) - fngy(ix,iy-1,igsp)) )

c ... IJ 2016/10/19 add MC neut flux if flags set
             if (get_neutral_moments .and. cmneutdiv_fng .ne. 0.0) then  
               jfld=1  
               sng_ue(ix,iy,jfld) = - ( 
     .                     (fngx_ue(ix,iy,jfld)-fngx_ue(ix1,iy, jfld))
     .           +fluxfacy*(fngy_ue(ix,iy,jfld)-fngy_ue(ix,iy-1,jfld)) )
     .           *( (ng(ix,iy,jfld)*ti(ix,iy))/
     .                                  (ng(ix,iy,jfld)*ti(ix,iy)) )
               resng(ix,iy,igsp) = resng(ix,iy,igsp) + 
     .                     cmneutdiv*cmneutdiv_fng*sng_ue(ix,iy,igsp)
             endif

            end do 
        end do
      endif

      end do   # end of igsp loop from the beginning of subroutine

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
            lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                         nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
            lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
               lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
               lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif


      END SUBROUTINE calc_gas_continuity_residuals


c --------------------------------------------------------------------------
c   Below subroutine neudifl is just like subroutine neudif, except that the
c   log of the gas density is used, and then converted back to give the 
c   physically meaningful gas variables (flux, velocity, etc)
c --------------------------------------------------------------------------

      subroutine neudifl

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1,nlmt,ffyo,ffyi
      integer iy1, methgx, methgy, iy2, jx
      logical isxyfl
      # Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
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
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
	  
*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10
	  
      do igsp = 1, ngsp

c.... First the flux in the x-direction

      do iy = j4, j8
         do ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgxa(ix,igsp) * sx(ix,iy) * (vtn + vtnp)*rt8opi/8
            csh = (1-isgasdc) * cdifg(igsp) * sx(ix,iy) * gxf(ix,iy) *
     .                ave( stretcx(ix,iy)*vtn**2/nuix(ix,iy,igsp),
     .                     stretcx(ix2,iy)*vtnp**2/nuix(ix2,iy,igsp) ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng +
     .                        rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                          0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
            qtgf = cngfx(igsp) * fgtdx(ix) * sx(ix,iy) *
     .            ave( stretcx(ix,iy)*gx(ix,iy)/nuix(ix,iy,igsp) ,
     .                 stretcx(ix2,iy)*gx(ix2,iy)/nuix(ix2,iy,igsp) )
     .                     * (vtn**2 - vtnp**2)
            if (isupgon(igsp).eq.1) then # only 2->pol project.;up has rest
               csh = csh*(1 - rrv(ix,iy)**2)
               qtgf = qtgf*(1 - rrv(ix,iy)**2)
            endif
            vygtan(ix,iy,igsp) = 0.  # vygtan is grad(T) rad vel on x-face
            if (isnonog .eq. 1 .and. iy .le. ny) then
                  grdnv =(    ( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                       -( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               vygtan(ix,iy,igsp) = exp( 0.5*
     .                     (log(tg(ix2,iy,igsp))+log(tg(ix,iy,igsp))) )*
     .                                  ( cngfx(igsp) / (mg(igsp)*0.5*
     .                         (nuix(ix,iy,igsp)+nuix(ix2,iy,igsp))) ) *
     .                             ( grdnv/cos(angfx(ix,iy)) - 
     .                       (log(tg(ix2,iy,igsp)) - log(tg(ix,iy,igsp)))
     .                                                 * gxf(ix,iy) ) 
               if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
                  vygtan(ix,iy,igsp) = 0.
               endif
               if (nxpt==2 .and. ix==ixrb(1)+1) then
                  vygtan(ix,iy,igsp) = 0.
               endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - vygtan(ix,iy,igsp)*sx(ix,iy)
            if (isupgon(igsp) .eq. 1) then
              qtgf = qtgf + rrv(ix,iy)*up(ix,iy,iigsp)*sx(ix,iy)
            endif
            qsh = csh * (lng(ix,iy,igsp)-lng(ix2,iy,igsp)) + qtgf
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  boundary face - shouldnt matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .              (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floxg
            floxg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = 
     .              floxg(ix,iy) + cngflox(1)*sx(ix,iy)*uu(ix,iy,1)
         floxg(ix,iy) = floxg(ix,iy)*2/(lng(ix,iy,igsp)+lng(ix2,iy,igsp))

        end do
         conxg(nx+1,iy) = 0
      end do

c.... Now the flux in the y-direction

      do iy = j1, j5
         do ix = i4, i8
	    t0 = max(tg(ix,iy,igsp),temin*ev)
	    t1 = max(tg(ix,iy+1,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgya(iy,igsp) * sy(ix,iy) * (vtn + vtnp)*rt8opi/8
            csh = (1-isgasdc) * (cdifg(igsp) *sy(ix,iy)/dynog(ix,iy)) *
     .                            ave( vtn**2/nuix(ix,iy,igsp) ,
     .                                 vtnp**2/nuix(ix,iy+1,igsp) ) +
     .            isgasdc * sy(ix,iy) * difcng / dynog(ix,iy) +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*dynog(ix,iy)*
     .                       0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))
c               csh = sy(ix,iy) * gyf(ix,iy) * ( (vtn**2+vtnp**2)/
c     .                 (nuix(ix,iy,igsp)+nuix(ix,iy+1,igsp)) )
            qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                     ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                          gy(ix,iy+1)/nuix(ix,iy+1,igsp) )
     .                    * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0 ) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp )) )
              endif
              qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                      ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                           gy(ix,iy+1)/nuix(ix,iy+1,igsp) ) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (lng(ix,iy,igsp)-lng(ix,iy+1,igsp)) + qtgf
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floyg
	    floyg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) =  
     .             floyg(ix,iy)+cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)

         floyg(ix,iy)=floyg(ix,iy)*2/(lng(ix,iy,igsp)+lng(ix,iy+1,igsp))
        end do
      end do

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             lng(0:nx+1,0:ny+1,igsp),flngx(0:nx+1,0:ny+1,igsp),
     .             flngy(0:nx+1,0:ny+1,igsp),0,methg)

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then

         do iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
ccc            MER: Set flag to apply xy flux limit except at target plates
               isxyfl = .true.
               do jx = 1, nxpt
                  if ( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .                 (ix==ixrb(jx).and.ixmxbcl==1) )isxyfl = .false.
               enddo
               grdnv =( (fym (ix,iy,1)*lng(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*lng(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*lng(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*lng(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*lng(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*lng(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*lng(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*lng(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*lng(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*lng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)

               difgx2 = ave( tg(ix ,iy,igsp)/nuix(ix ,iy,igsp),
     .                       tg(ix2,iy,igsp)/nuix(ix2,iy,igsp) )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))

               flngxy(ix,iy,igsp) = difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                             (lng(ix2,iy,igsp) - lng(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)

c...  Now flux limit with flalfgxy
               t0 = max(tg(ix,iy,igsp),temin*ev)
               t1 = max(tg(ix2,iy,igsp),temin*ev)
               vtn = sqrt( t0/mg(igsp) )
               vtnp = sqrt( t1/mg(igsp) )
               qfl = flalfgxya(ix,igsp)*sx(ix,iy)* (vtn + vtnp)*rt8opi/8
ccc   MER NOTE:  no xy flux limit for ix=0 or ix=nx in original code
               if (isxyfl) then
                  fngxy(ix,iy,igsp) = fngxy(ix,iy,igsp) /
     .                           sqrt( 1 + (fngxy(ix,iy,igsp)/qfl)**2 )
               endif
            end do
        end do

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Flux-limit the total poloidal flux here
            do iy = j4, j8
              do ix = i1, i5
                flngx(ix,iy,igsp)= flngx(ix,iy,igsp)-flngxy(ix,iy,igsp)
              enddo
            enddo

      endif
c...  Finished with nonorthogonal mesh part

c ... Calculate the neutral flow velocity from v = flux/ng; these are
c ... diagnostic if isupgon=0, but used for uu and vy of the inertial
c ... gas if isupgon=1.  However, even when isupgon=1, the particle
c ... fluxes are fngx -> fnix and fngy -> fngy in pandf, i.e., methg
c ... determines the differencing for the inertial particle fluxes, not 
c ... methn
      do iy = j1, j5
         do ix = i1,i5
            ix1 = ixp1(ix,iy)
            uug(ix,iy,igsp) = flngx(ix,iy,igsp) / sx(ix,iy)
            vyg(ix,iy,igsp) = flngy(ix,iy,igsp) / sy(ix,iy)
c --------------- transfer inertial gas velocities to neutral ion species
            if (isupgon(igsp).eq.1) then
               vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
cc               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            end if
         enddo
cfw ----   If doing only the outer half we want this boundary condition:
         if (iy.le.iysptrx1(1).and.isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp) = 0
      enddo

c **- loop for uu just as in the previous version - needed for correct Jac?
      if (isupgon(igsp) .eq. 1) then
         do iy = j4, j6
            do ix = i1, i6
               uu(ix,iy,iigsp) = uug(ix,iy,igsp)
            enddo
         enddo
      endif

c.... Calculate the particle flux, fnix,y, from flux of lng, i.e., flngx,y
      do iy = j4, j8
         do ix = i1, i5
            ix2 = ixp1(ix,iy)
            fngx(ix,iy,igsp) = flngx(ix,iy,igsp) *
     .                     exp(0.5*(lng(ix,iy,igsp)+lng(ix2,iy,igsp)))
         enddo
      enddo
c ...   now do fniy
      do iy = j1, j5    # same loop ranges as for fngy in fd2tra
         do ix = i4, i8
            fngy(ix,iy,igsp) = flngy(ix,iy,igsp) *
     .                     exp(0.5*(lng(ix,iy,igsp)+lng(ix,iy+1,igsp)))
         enddo
      enddo

c.... Calculate the residual for the gas equation for diffusive neutral case

      if (isupgon(igsp).eq.0) then
         do iy = j2, j5
            if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
	      fngx(nxc-1,iy,igsp) = 0.
	      fngx(nxc,  iy,igsp) = 0.
	      fngx(nxc+1,iy,igsp) = 0.
	    endif
            if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
            if (nxpt==2.and.ixmxbcl==1) fngx(ixrb(1)+1,iy,igsp)=0.
                   do ix = i2, i5
               ix1 = ixm1(ix,iy)
               resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .                         psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .                       + volpsorg(ix,iy,igsp)                     
     .                       - fngx(ix,iy,igsp) + fngx(ix1,iy  ,igsp)
     .             - fluxfacy*(fngy(ix,iy,igsp) - fngy(ix ,iy-1,igsp))
     .                       + psgov_use(ix,iy,igsp)*vol(ix,iy)

               if (igsp.eq.1 .and. ishymol.eq.1)  
     .              resng(ix,iy,igsp) = resng(ix,iy,igsp)+psordis(ix,iy,2)
            end do
        end do
      endif

      end do # end of igsp loop from the beginning of subroutine

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
	   lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                         nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
	   lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
	      lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
	      lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudifl
c --------------------------------------------------------------------------
c --------------------------------------------------------------------------
      subroutine neudifo

c ..  Older version of neudif where the gas velocities are deduced from
c ..  the gas fluxes and then used to form fnix if isupgon(igsp)=1

      implicit none

*  -- local variables
      real vtn, vtnp, qr, qtgf, nconv, grdnv, difgx2
      real tnuiz,ngnot,lmfp,ty0,ty1
      integer iy1, methgx, methgy, iy2, jx

      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nxpt
      Use(Xpoint_indices)      # ixlb,ixpt1,ixpt2,ixrb,iysptrx
      Use(Share)    # nxpt,geometry,nxc,isnonog,cutlo,islimon,ix_lim,iy_lims
      Use(Phyvar)   # pi,me,mp,ev,qe,rt8opi
      Use(UEpar)    # methg,
                    # qfl,csh,qsh,cs,
                    # isupgon,iigsp,nlimgx,nlimgy
       # Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,ix1,ix2,ix3,ix4,ix5,ix6
      real t,t0,t1,t2,a
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
      Use(Rhsides)  # resng,psor,psorg,psorrg,sniv
      Use(Comtra)   # flalfgx,flalfgy
      Use(Locflux)  # floxg,floyg,conxg,conyg
      Use(Indices_domain_dcl)    # iymnbcl,iymxbcl
      Use(Volsrc)   # volpsorg
	  
*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

c ------------------
      methgx = mod(methg, 10)
      methgy = methg/10

      do igsp = 1, ngsp

c.... First the flux in the x-direction

      do iy = j4, j8
         do ix = i1, i5
            iy1 = max(0,iy-1)
            iy2 = min(ny+1,iy+1)
            ix2 = ixp1(ix,iy)
            ix4 = ixp1(ix,iy1)
            ix6 = ixp1(ix,iy2)
            t0 = max(tg(ix,iy,igsp),temin*ev)
            t1 = max(tg(ix2,iy,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgxa(ix,igsp) * sx(ix,iy) * (vtn + vtnp)*rt8opi*
     .           (ng(ix,iy,igsp)*gx(ix,iy) + ng(ix2,iy,igsp)*gx(ix2,iy))
     .                                      / (8*(gx(ix,iy)+gx(ix2,iy)))
c            csh = sx(ix,iy) * gxf(ix,iy) *
c     .                 (stretcx(ix,iy)*vtn**2+stretcx(ix2,iy)*vtnp**2)/
c     .                 (nuix(ix,iy,igsp)+nuix(ix2,iy,igsp))
            csh = (1-isgasdc) * cdifg(igsp)* sx(ix,iy) * gxf(ix,iy) *
     .                ave( stretcx(ix,iy)*vtn**2/nuix(ix,iy,igsp),
     .                     stretcx(ix2,iy)*vtnp**2/nuix(ix2,iy,igsp) ) +
     .              isgasdc * sx(ix,iy) * gxf(ix,iy) * difcng +
     .                       rld2dxg(igsp)**2*sx(ix,iy)*(1/gxf(ix,iy))*
     .                          0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
            qtgf = cngfx(igsp) * fgtdx(ix) * sx(ix,iy) *
     .            ave( stretcx(ix,iy)*gx(ix,iy)/nuix(ix,iy,igsp) ,
     .                 stretcx(ix2,iy)*gx(ix2,iy)/nuix(ix2,iy,igsp) )
     .                     * (vtn**2 - vtnp**2)
            vygtan(ix,iy,igsp) = 0.
            if (isnonog .eq. 1 .and. iy .le. ny) then
               if (isintlog .eq. 0) then
                  grdnv =( fym (ix,iy,1)*tg(ix2,iy1,igsp) +
     .                     fy0 (ix,iy,1)*tg(ix2,iy ,igsp) +
     .                     fyp (ix,iy,1)*tg(ix2,iy2,igsp) +
     .                     fymx(ix,iy,1)*tg(ix ,iy1,igsp) +
     .                     fypx(ix,iy,1)*tg(ix, iy2,igsp) -
     .                     fym (ix,iy,0)*tg(ix ,iy1,igsp) -
     .                     fy0 (ix,iy,0)*tg(ix ,iy ,igsp) -
     .                     fyp (ix,iy,0)*tg(ix ,iy2,igsp) - 
     .                     fymx(ix,iy,0)*tg(ix4,iy1,igsp) -
     .                     fypx(ix,iy,0)*tg(ix6,iy2,igsp) )/dxnog(ix,iy)
               elseif (isintlog .eq. 1) then
                  grdnv =( exp( fym (ix,iy,1)*log(tg(ix2,iy1,igsp)) +
     .                          fy0 (ix,iy,1)*log(tg(ix2,iy ,igsp)) +
     .                          fyp (ix,iy,1)*log(tg(ix2,iy2,igsp)) +
     .                          fymx(ix,iy,1)*log(tg(ix ,iy1,igsp)) +
     .                          fypx(ix,iy,1)*log(tg(ix, iy2,igsp)) )
     .                    -exp( fym (ix,iy,0)*log(tg(ix ,iy1,igsp)) +
     .                          fy0 (ix,iy,0)*log(tg(ix ,iy ,igsp)) +
     .                          fyp (ix,iy,0)*log(tg(ix ,iy2,igsp)) + 
     .                          fymx(ix,iy,0)*log(tg(ix4,iy1,igsp)) +
     .                          fypx(ix,iy,0)*log(tg(ix6,iy2,igsp)) ) )/ 
     .                                                      dxnog(ix,iy)
               endif
               vygtan(ix,iy,igsp) = ( cngfx(igsp) / (mg(igsp)*0.5*
     .                         (nuix(ix,iy,igsp)+nuix(ix2,iy,igsp))) ) *
     .                             ( grdnv/cos(angfx(ix,iy)) - 
     .                             (tg(ix2,iy,igsp) - tg(ix,iy,igsp))
     .                                                 * gxf(ix,iy) ) 
               if (islimon.eq.1.and. ix.eq.ix_lim.and. iy.ge.iy_lims) then
                  vygtan(ix,iy,igsp) = 0.
               endif
               if (nxpt==2 .and. ix==ixrb(1)+1 .and. ixmxbcl==1) then
                  vygtan(ix,iy,igsp) = 0.
               endif
            endif
c --- In the neutral momentum case we are after the 2-direction flux,
c --- which has no non-orthogonal part.
            qtgf = qtgf - (1-isupgon(igsp))*vygtan(ix,iy,igsp)*sx(ix,iy)
            nconv = 2.0*(ng(ix,iy,igsp)*ng(ix2,iy,igsp)) /
     .                  (ng(ix,iy,igsp)+ng(ix2,iy,igsp))
c...   Use upwind for "convective" grad T term if methgx .ne. 2
            if(methgx.ne.2) nconv =
     .                         ng(ix ,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ng(ix2,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ng(ix,iy,igsp)-ng(ix2,iy,igsp)) + qtgf * nconv
            qr = abs(qsh/qfl)
c...  Because guard-cell values may be distorted from B.C., possibly omit terms on
c...  boundary face - shouldnt matter(just set BC) except for guard-cell values
            do jx = 1, nxpt
               if( (ix==ixlb(jx).and.ixmnbcl==1) .or.
     .             (ix==ixrb(jx).and.ixmxbcl==1) ) then
                  qr = gcfacgx*qr
                  qtgf = gcfacgx*qtgf
               endif
            enddo
            conxg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifxg_aug .eq. 1) conxg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floxg
            floxg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)
c...  now add the convective velocity for charge-exchange neutrals
         if(igsp .eq. 1) floxg(ix,iy) = 
     .             floxg(ix,iy) + cngflox(1)*sx(ix,iy)*uu(ix,iy,1)

         end do
         conxg(nx+1,iy) = 0
       end do

c.... Now the flux in the y-direction

      do iy = j1, j5
         do ix = i4, i8
	    t0 = max(tg(ix,iy,igsp),temin*ev)
	    t1 = max(tg(ix,iy+1,igsp),temin*ev)
            vtn = sqrt( t0/mg(igsp) )
            vtnp = sqrt( t1/mg(igsp) )
            qfl = flalfgya(iy,igsp) * sy(ix,iy) * (vtn + vtnp)*rt8opi*
     .                              ( ngy0(ix,iy,igsp)*gy(ix,iy) + 
     .                                ngy1(ix,iy,igsp)*gy(ix,iy+1) ) / 
     .                                     (8*(gy(ix,iy)+gy(ix,iy+1)))
            csh = (1-isgasdc) * (cdifg(igsp) *sy(ix,iy) /dynog(ix,iy)) *
     .                            ave( vtn**2/nuix(ix,iy,igsp) ,
     .                                 vtnp**2/nuix(ix,iy+1,igsp) ) +
     .            isgasdc * sy(ix,iy) * difcng / dynog(ix,iy) +
     .                      rld2dyg(igsp)**2*sy(ix,iy)*dynog(ix,iy)*
     .                       0.5*(nuiz(ix,iy,igsp)+nuiz(ix,iy+1,igsp))
c               csh = sy(ix,iy) * ( ((vtn**2+vtnp**2)/ dynog(ix,iy)) /
c     .                 (nuix(ix,iy,igsp)+nuix(ix,iy+1,igsp)) )
            qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                     ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                          gy(ix,iy+1)/nuix(ix,iy+1,igsp) )
     .                    * (vtn**2 - vtnp**2)
            if (isnonog.eq.1 .and. iy.le.ny) then
              if (isintlog .eq. 0) then
                ty0 = fxm (ix,iy,0)*tg(ixm1(ix,iy)  ,iy  ,igsp) + 
     .                fx0 (ix,iy,0)*tg(ix           ,iy  ,igsp) +
     .                fxp (ix,iy,0)*tg(ixp1(ix,iy)  ,iy  ,igsp) +
     .                fxmy(ix,iy,0)*tg(ixm1(ix,iy+1),iy+1,igsp) +
     .                fxpy(ix,iy,0)*tg(ixp1(ix,iy+1),iy+1,igsp)
                ty1 = fxm (ix,iy,1)*tg(ixm1(ix,iy+1),iy+1,igsp) + 
     .                fx0 (ix,iy,1)*tg(ix           ,iy+1,igsp) +
     .                fxp (ix,iy,1)*tg(ixp1(ix,iy+1),iy+1,igsp) +
     .                fxmy(ix,iy,1)*tg(ixm1(ix,iy)  ,iy  ,igsp) +
     .                fxpy(ix,iy,1)*tg(ixp1(ix,iy)  ,iy  ,igsp)
              elseif (isintlog .eq. 1) then
                ty0=exp(fxm (ix,iy,0)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) + 
     .                  fx0 (ix,iy,0)*log(tg(ix           ,iy  ,igsp)) +
     .                  fxp (ix,iy,0)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) +
     .                  fxmy(ix,iy,0)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) +
     .                  fxpy(ix,iy,0)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) )
                ty1=exp(fxm (ix,iy,1)*log(tg(ixm1(ix,iy+1),iy+1,igsp)) + 
     .                  fx0 (ix,iy,1)*log(tg(ix           ,iy+1,igsp)) +
     .                  fxp (ix,iy,1)*log(tg(ixp1(ix,iy+1),iy+1,igsp)) +
     .                  fxmy(ix,iy,1)*log(tg(ixm1(ix,iy)  ,iy  ,igsp)) +
     .                  fxpy(ix,iy,1)*log(tg(ixp1(ix,iy)  ,iy  ,igsp)) )
              endif
              qtgf = cngfy(igsp) * fgtdy(iy) * sy(ix,iy) * 
     .                      ave( gy(ix,iy)/nuix(ix,iy,igsp) ,
     .                           gy(ix,iy+1)/nuix(ix,iy+1,igsp) ) *
     .                                            (ty0 - ty1)/mg(igsp)
            endif       # Better interpolation of nuix could be done here
            nconv = 2.0*(ngy0(ix,iy,igsp)*ngy1(ix,iy,igsp)) /
     .                  (ngy0(ix,iy,igsp)+ngy1(ix,iy,igsp)) 
c...   Use upwind for "convective" grad T term if methgy .ne. 2
            if(methgy.ne.2) nconv =
     .                         ngy0(ix,iy,igsp)*0.5*(1+sign(1.,qtgf)) +
     .                         ngy1(ix,iy,igsp)*0.5*(1-sign(1.,qtgf))
            qsh = csh * (ngy0(ix,iy,igsp)-ngy1(ix,iy,igsp)) + qtgf*nconv
            qr = abs(qsh/qfl)
            if(iy.eq.0 .and. iymnbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            if(iy.eq.ny .and. iymxbcl.eq.1) then
               qr = gcfacgy*qr
               qtgf = gcfacgy*qtgf
            endif
            conyg(ix,iy) = csh / (1 + qr**flgamg)**(1/flgamg)
            if (isdifyg_aug .eq. 1) conyg(ix,iy) = csh*(1+qr) #augment diffusion

c...  the temperature gradient term is included in floyg
	    floyg(ix,iy) = qtgf / (1 + qr**flgamg)**(1/flgamg)  
c...  now add the convective velocity for the charge-exchange species
         if(igsp .eq. 1) floyg(ix,iy) = 
     .               floyg(ix,iy)+cngfloy(1)*sy(ix,iy)*vy(ix,iy,1)
         end do
        end do

*  --------------------------------------------------------------------
*  compute the neutral particle flow
*  --------------------------------------------------------------------

      call fd2tra (nx,ny,floxg,floyg,conxg,conyg,
     .             ng(0:nx+1,0:ny+1,igsp),fngx(0:nx+1,0:ny+1,igsp),
     .             fngy(0:nx+1,0:ny+1,igsp),0,methg)

c ... Calculate the neutral flow velocity from v = flux/ng
      do iy = j1, j5
         do ix = i1,i5
            ix1 = ixp1(ix,iy)
            uug(ix,iy,igsp) = fngx(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix1,iy,igsp))
     .                                                      *sx(ix,iy) )
            vyg(ix,iy,igsp) = fngy(ix,iy,igsp) / (
     .                        0.5*(ng(ix,iy,igsp)+ng(ix,iy+1,igsp))
     .                                                      *sy(ix,iy) )
            if (isupgon(igsp).eq.1) then
c --------------- We need to transfer the diffusive radial neutral
c --------------- velocity to the "ion" species containing the neutrals
               vy(ix,iy,iigsp) = vyg(ix,iy,igsp)
            end if
          end do
cfw ----   If doing only the outer half we want this boundary condition:
         if (iy.le.iysptrx1(1).and.isfixlb(1).eq.2) uug(ixpt2(1),iy,igsp) = 0
        end do

c ... For nonorthogonal mesh and diffusive neutrals, limit fngy for pump out
      if (isnonog.eq.1 .and. isupgon(igsp).eq.0) then
         do iy = j1, j5
            do ix = i4, i8
               if (nlimgy*fngy(ix,iy,igsp)*
     .                 (ng(ix,iy,igsp)-ng(ix,iy+1,igsp)) .lt. 0) then
                  fngy(ix,iy,igsp) = fngy(ix,iy,igsp)/( 1 - 2*nlimgy +
     .                      nlimgy*(ng(ix,iy+1,igsp)/ng(ix,iy,igsp)+
     .                              ng(ix,iy,igsp)/ng(ix,iy+1,igsp)) )
               endif
            enddo
         enddo     
      endif            

c...  Addition for nonorthogonal mesh
      if (isnonog .eq. 1) then

         do iy = j1, min(j6, ny)
            iy1 = max(iy-1,0)
            do ix = i1, min(i6, nx)
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               ix3 = ixm1(ix,iy1)
               ix4 = ixp1(ix,iy1)
               ix5 = ixm1(ix,iy+1)
               ix6 = ixp1(ix,iy+1) 
               if (methgx .eq. 6) then  # log interpolation
               grdnv =( exp(fym (ix,iy,1)*log(ng(ix2,iy1 ,igsp)) + 
     .                      fy0 (ix,iy,1)*log(ng(ix2,iy  ,igsp)) +
     .                      fyp (ix,iy,1)*log(ng(ix2,iy+1,igsp)) + 
     .                      fymx(ix,iy,1)*log(ng(ix ,iy1 ,igsp)) +
     .                      fypx(ix,iy,1)*log(ng(ix, iy+1,igsp)))
     .                - exp(fym (ix,iy,0)*log(ng(ix ,iy1 ,igsp)) +
     .                      fy0 (ix,iy,0)*log(ng(ix ,iy  ,igsp)) +
     .                      fyp (ix,iy,0)*log(ng(ix ,iy+1,igsp)) +
     .                      fymx(ix,iy,0)*log(ng(ix4,iy1 ,igsp)) + 
     .                      fypx(ix,iy,0)*log(ng(ix6,iy+1,igsp))) ) / 
     .                                                  dxnog(ix,iy)
               elseif (methgx .eq. 7) then  # inverse interpolation
               grdnv =( 1/(fym (ix,iy,1)/ng(ix2,iy1 ,igsp) + 
     .                     fy0 (ix,iy,1)/ng(ix2,iy  ,igsp) +
     .                     fyp (ix,iy,1)/ng(ix2,iy+1,igsp) + 
     .                     fymx(ix,iy,1)/ng(ix ,iy1 ,igsp) +
     .                     fypx(ix,iy,1)/ng(ix, iy+1,igsp))
     .                - 1/(fym (ix,iy,0)/ng(ix ,iy1 ,igsp) +
     .                     fy0 (ix,iy,0)/ng(ix ,iy  ,igsp) +
     .                     fyp (ix,iy,0)/ng(ix ,iy+1,igsp) +
     .                     fymx(ix,iy,0)/ng(ix4,iy1 ,igsp) + 
     .                     fypx(ix,iy,0)/ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               else                   # linear interpolation
               grdnv =( (fym (ix,iy,1)*ng(ix2,iy1 ,igsp) + 
     .                   fy0 (ix,iy,1)*ng(ix2,iy  ,igsp) +
     .                   fyp (ix,iy,1)*ng(ix2,iy+1,igsp) + 
     .                   fymx(ix,iy,1)*ng(ix ,iy1 ,igsp) +
     .                   fypx(ix,iy,1)*ng(ix, iy+1,igsp))
     .                - (fym (ix,iy,0)*ng(ix ,iy1 ,igsp) +
     .                   fy0 (ix,iy,0)*ng(ix ,iy  ,igsp) +
     .                   fyp (ix,iy,0)*ng(ix ,iy+1,igsp) +
     .                   fymx(ix,iy,0)*ng(ix4,iy1 ,igsp) + 
     .                   fypx(ix,iy,0)*ng(ix6,iy+1,igsp)) ) / 
     .                                                  dxnog(ix,iy)
               endif
               difgx2 = ave( tg(ix ,iy,igsp)/nuix(ix ,iy,igsp),
     .                       tg(ix2,iy,igsp)/nuix(ix2,iy,igsp) )/mg(igsp)
     .                         + rld2dxg(igsp)**2*(1/gxf(ix,iy)**2)*
     .                           0.5*(nuiz(ix,iy,igsp)+nuiz(ix2,iy,igsp))
               fngxy(ix,iy,igsp) = difgx2*(grdnv/cos(angfx(ix,iy)) -
     .                             (ng(ix2,iy,igsp) - ng(ix,iy,igsp))*
     .                                 gxf(ix,iy) ) * sx(ix,iy)
            end do
        end do

c...  Fix the total fluxes; note the loop indices same as fd2tra
c...  Dont bother if  solving the neutral momentum equation
         if (isupgon(igsp).eq.0) then
            do iy = j4, j8
               do ix = i1, i5
                  ix2 = ixp1(ix,iy)
                  fngx(ix,iy,igsp) = fngx(ix,iy,igsp) - fngxy(ix,iy,igsp)
c ...          adjust fluxes to prevent pump out
                  if (nlimgx*fngx(ix,iy,igsp)
     .                   *(ng(ix,iy,igsp)-ng(ix2,iy,igsp)) .lt. 0.) then
                    fngx(ix,iy,igsp) = fngx(ix,iy,igsp)/( 1 - 2*nlimgx +
     .                          nlimgx*(ng(ix2,iy,igsp)/ng(ix,iy,igsp) +
     .                                  ng(ix,iy,igsp)/ng(ix2,iy,igsp)) )
                  endif
               enddo
            enddo
         endif

      endif
c...  Finished with nonorthogonal mesh part

c.... Calculate the residual or right-hand-side for the gas equation

      if (isupgon(igsp).eq.0) then
         do iy = j2, j5
            if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
	      fngx(nxc-1,iy,igsp) = 0.
	      fngx(nxc,  iy,igsp) = 0.
	      fngx(nxc+1,iy,igsp) = 0.
	    endif
            if (islimon.ne.0.and.iy.ge.iy_lims) fngx(ix_lim,iy,igsp)=0.
            if (nxpt==2.and.ixmxbcl==1) fngx(ixrb(1)+1,iy,igsp)=0.
                   do ix = i2, i5
               ix1 = ixm1(ix,iy)
               resng(ix,iy,igsp) = cngsor * (psorg(ix,iy,igsp) +
     .                         psorcxg(ix,iy,igsp) + psorrg(ix,iy,igsp))
     .                       + volpsorg(ix,iy,igsp)                     
     .                       - fngx(ix,iy,igsp) + fngx(ix1,iy  ,igsp)
     .                       - fngy(ix,iy,igsp) + fngy(ix ,iy-1,igsp)
     .                       + psgov_use(ix,iy,igsp)*vol(ix,iy)
               if (igsp.eq.1 .and. ishymol.eq.1) 
     .              resng(ix,iy,igsp) = resng(ix,iy,igsp)-psordis(ix,iy,2)
            end do
        end do
      else if (isupgon(igsp).eq.1) then

c --- Form the poloidal velocity uu(,,2) from
c --- a) the projection of the neutral parallel velocity up(,,2),
c --- b) the projection of the 2-direction gas velocity and
c --- c) vygtan, the grad(T) part of the non-orth diffusive radial velocity.
c --- d) fngxy, the grad(n) part of the non-orth diffusive radial velocity.
c --- The fngxy contribution could have been kept separate and added to
c --- fnix in PARBAL, but we include it here so that it automatically gets
c --- taken into account in PARBAL and MOMBAL_B2.
c --- By multiplying uu(,,2) with sx*ng in PARBAL we get the 
c --- TOTAL neutral particle flux out of the poloidal face.
c --- By multiplying uu(,,2) with sx*ng*mi(1)*up(,,iigsp) in MOMBAL_B2
c --- we get the TOTAL parallel momentum flux out of the poloidal face.
c --- Note that rrv=Bpol/B is defined at a velocity point
c --- In order to get Bt/B at a VELOCITY point we cannot use rbfbt,
c --- use (Bt/B)**2=1-rrv**2 instead.
         do iy = j4, j6
            do ix = i1, i6
               ix2 = ixp1(ix,iy)
               uu(ix,iy,iigsp) = rrv(ix,iy)*up(ix,iy,iigsp) +
     .              (1.-rrv(ix,iy)*rrv(ix,iy))*uug(ix,iy,igsp)
               if (isnonog .eq. 1) uu(ix,iy,iigsp) =
     .              uu(ix,iy,iigsp) - vygtan(ix,iy,igsp) -
     .              fngxy(ix,iy,igsp) / (sx(ix,iy)*
     .                    0.5*(ng(ix,iy,igsp)+ng(ix2,iy,igsp)))
            end do
        end do
      end if

      end do

c...  Special coding for the 1-D gas-box model
      if (is1D_gbx.eq.1) then
         tnuiz = 0.
         do ix = 0, ixgb-1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp + agdc*(ix-ixgb))*
     .                                          nuiz(ix,1,1)/gx(ix,1)
         enddo        
         do ix = ixgb, nx+1
	    lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
            tnuiz = tnuiz + exp(-pcolwid/lmfp)*nuiz(ix,1,1)/gx(ix,1)
         enddo
         ngnot = max( recycp(1)*fnix(nx,1,1)/(sx(nx,1)*tnuiz),
     .                                                    ngbackg(1) )
         if (i5+1 .gt. ixgb) then  # In gas-box region; i5+1 to be safe
            do ix = ixgb, nx+1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1) = -nurlxg*vol(ix,1)*( ng(ix,1,1) -
     .                                     ngnot*exp(-pcolwid/lmfp) )
            enddo
         endif
         if (i2 .lt. ixgb) then  # In gas-free region, make ng=10*ngbackg
                                 # -ng(ixgb)*exp(agdc*(ix-ixgb))
            do ix = i2, ixgb-1
	       lmfp = sqrt(2*tg(ix,1,1)/(mi(1)*nuiz(ix,1,1)*nucx(ix,1,1)))
               resng(ix,1,1)= -nurlxg*vol(ix,1)*( ng(ix,1,1)-ngnot*
     .                            exp(-pcolwid/lmfp + agdc*(ix-ixgb)) )
            enddo
         endif
      endif

      return

      end
c
c --------------------------------------------------------------------------
c END subroutine neudifo - old version of neudif
c --------------------------------------------------------------------------


