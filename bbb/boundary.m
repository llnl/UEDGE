c!include "../mppl.h"
c!include "../sptodp.h"
c-----------------------------------------------------------------------
      subroutine bouncon(neq, yl, yldot)

*   Bouncon provides the evaluation of the equations for the boundaries.

      implicit none

      integer neq
      real yl(neq), yldot(neq)

      Use(Dim)      # nx,ny,nhsp,nzspt,nzsp,nisp,ngsp,nusp,nxpt
      Use(Share)    # nxpt,nxc,geometry,cutlo,islimon,ix_lim,iy_lims
		    # isudsym
      Use(Xpoint_indices) # ixlb,ixpt1,ixpt2,ixrb,iysptrx1,iysptrx2
      Use(Math_problem_size)   # neqmx 
      Use(Phyvar)
      Use(UEpar)    # isnewpot,r0slab,cslim,dcslim,csfaclb,csfacrb,csfacti,
                    # isnion,isupon,isteon,istion,isngon,isnionxy,isuponxy,
                    # isteonxy,istionxy,isngonxy,isphionxy
      Use(Aux)      # ixmp
      Use(Coefeq)   # fac2sp,cf2ef,exjbdry
      Use(Bcond)    # iflux,ncore,tcoree,tcorei,tbmin,nbmin,ngbmin,
                    # tepltl,tipltl,tepltr,tipltr,
                    # istewc,istiwc,istepfc,istipfc,
                    # tewalli,tiwalli,tewallo,tiwallo,isextrnp,
                    # isextrnpf,isextrtpf,isextrngc,isextrnw,isextrtw,
                    # iflcore,pcoree,pcorei,ifluxni,ckinfl,isupss,
                    # isnwconiix,isnwconoix,nwalli,nwallo,iscpli,iscplo,
                    # fngysi,fngyso,albedoo,albedoi,matwallo,matwalli,
                    # sinphi,isfixlb,nib,teb,tib,nibprof,tebprof,tibprof,
                    # engbsr,epsbs,rlimiter,ngcore,isngcore,isutcore,
                    # ixfixnc,incixc,isupcore,isfixrb,chemsputi,chemsputo
                    # islbcn,islbcu,islbce,islbci,islbcg,isexunif
                    # fchemygwi,fchemylb,fphysylb,fchemygwo,fchemyrb,fphysyrb
                    # xcnearlb,xcnearrb,openbox,fqpsatlb,fqpsatrb
                    # cfueb,ikapmod,cfvytanbc
      Use(Parallv)  # nxg,nyg
      Use(Selec)    # i1,i2,i3,i4,i5,i6,i7,j1,j2,j3,j4,j5,j6,j7,xlinc
      Use(Comgeo)   # gx,gy,gyf,sx,sy,xcwi,xcwo,yylb,rrv,sygytotc,isixcore
      Use(Compla)   # mi, mg
      Use(Comflo)   # fqx,fqy,fnix,fniy,feex,feey,feix,feiy,fngx,fngy
                    # fdiaxlb, fdiaxrb
      Use(Conduc)   # visx
      Use(Indexes)
      Use(Ynorm)    # temp0,nnorm,ennorm
      Use(Poten)    # newbcl,newbcr,bcee,bcei,rsigpl,bcel,bcer,bcil,bcir,
                    # kappal,kappar,bctype,phi0l,phi0r,isfdiax
      Use(Rccoef)   # recylb,recyrb,alblb,albrb,recycw,sputtr,
                    # recycm,recyce,recycmlb,recycmrb,recyllim,recyrlim
      Use(Bfield)   # rbfbt,btot
      Use(Imprad)   # isimpon
      Use(Impurity_source_flux)   # fnzysi,fnzyso
      Use(Gradients)  # ey
      Use(RZ_grid_info)      # rm
      Use(Indices_domain_dcl)   #ixmxbcl,ixmnbcl,iymxbcl,iymnbcl,ispwrbcl
      Use(Interp)
      Use(Jacaux)   # yldot_diag
      Use(Npes_mpi) # npes
      Use(Indices_domain_dcg) # ndomain,ispwrbc
c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in

      Use(MCN_dim)
      Use(MCN_sources) # edisspl, edisspr, cmntgpl, cmntgpl

c...  local scalars
      real totfeix, totfeex, kfeix, vyn, cosphi,
     .     ueb, nbound, tbound, vxn, ut0, sumb, feeytotc, feiytotc,
     .     r_major, fniytotc, fng_chem, vbound, eng_sput, flx_incid,
     .     yld_chm, t0p, zflux_chm, fqytotc, flux_inc,
     .     totfnex, totfnix, fqpsate, qpfac, aq, expkmx, arglgphi, faceel,
     .     faceel2, csfac, lambdae, uztotc, uztotc1, uztotc2,
     .     fngytotc, fmiytotc, sytotc, f_cgpld, vparn, sfeeytotc, sfeiytotc
      integer ii,isphion2, nzsp_rt, jz
      real hflux, zflux
      integer ifld, ihyd, iimp, ix_fl_bc, ixc1, igsp2
      real dif_imp_flux, fng_alb, fngyw, nharmave
      real upbdry, upbdry1, upbdry2, uugoal, fniy_recy, lengg, xtotc
      integer ixt, ixt1, ixt2, ixt3, jx, ixc, ierr
      integer ixtl, ixtl1, ixtr,ixtr1
      #Former Aux module variables
      integer ix,iy,igsp,iv,iv1,iv2,iv3,iv4,ix1,ix2,ix3,ix4
      real t0,t1

*  -- external procedures --
      real sdot, yld96, kappa
      external sdot

*  -- procedures --
      real ave
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)

*****************************************************************
*  --  Here we write the equations for the boundaries.
*****************************************************************
c...  now we reset the boundary conditions around the edge

c...  Initialization for constant
      expkmx = exp(-kappamx)

c ====================================================================
c ======================== The iy=0 boundary =========================
c ====================================================================

c...  do the iy = 0 boundary
c...  if extrapolation b.c. on p.f. region, isextrpf=1, otherwise isextrpf=0
      if (iymnbcl .eq. 0) goto 1100   # interior domain boundary; no bdry eqn
      ixc1 = max(0, ixpt1(1)+1)       # 1st  core cell;used for core flux BC
      if ((geometry(1:9)=="snowflake" .and. geometry(10:11).ne."15")
     &                     .or. geometry=="dnXtarget") then
        ix_fl_bc = min(ixpt2(1), nx)
      else
        ix_fl_bc = min(ixpt2(nxpt), nx) # last core cell;use for flux BC
      endif 

c  Note: j3 is local range index for iy passed from pandf in oderhs.m
      if (j3 .le. isextrnpf .or. j3 .le. isextrtpf .or. 
     .                                         j3 .le. isextrngc) then
              # this is a very large if-loop; check for 'j3 index' to find end

      do 274 ifld = 1, nisp
         do 272 ix = i4+1-ixmnbcl, i8-1+ixmxbcl # ix-loop over ion dens
           if (isnionxy(ix,0,ifld)==1) then
             iv1 = idxn(ix,0,ifld)

c ...  First check if this is a neutral species (zi=0) & set BC for it
c ---  Different boundary conditions if neutral momentum equation use;
c ---  typically hydrogen only, so DIVIMP chem sputt not used here.
c ---  Below, isupon(1) & ngbackg(1) used, so implies hydrogen
             if (isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0.0) then
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                 if (isngcore(1) .eq. 0) then
                   t0 = max(tg(ix,0,1),tgmin*ev)
                   vyn = sqrt( 0.5*t0/(pi*mi(ifld)) )
                   nharmave = 2.*(ni(ix,0,ifld)*ni(ix,1,ifld)) /
     .                           (ni(ix,0,ifld)+ni(ix,1,ifld))
                   fng_alb = (1-albedoc(1))*
     .                                  nharmave*vyn*sy(ix,0)
                   yldot(iv1) = -nurlxg*(fniy(ix,0,ifld) + fng_alb)/
     .                                     (vpnorm*sy(ix,0)*n0(ifld))
                 elseif (isngcore(1) .eq. 1) then
                   yldot(iv1) = nurlxn*(ngcore(1)-ni(ix,0,ifld))/
     .                                                     n0(ifld)
                 elseif (isngcore(1) .eq. 2) then
                   call xerrab("*** isngcore=2 option unvailable  ***")
                   lengg = sqrt(ti(ix,0)/
     .                              (mi(1)*nuix(ix,0,1)*nuiz(ix,0,1)))
                   yldot(iv1) = -nurlxn*( ni(ix,0,ifld) -
     .                           max( ni(ix,1,ifld) -
     .                        0.5*(ni(ix,1,ifld) + ni(ix,0,ifld))/
     .                        (gyf(ix,0)*lengg), 0.1*ngbackg(1) ) )
     .                                                       /n0(ifld)
                 elseif (isngcore(1) .eq. 3) then
                   nbound =  ng(ix,1,1) - gyf(ix,1)*
     .                          (ng(ix,2,1)-ng(ix,1,1))/gyf(ix,0)
                   nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                           ng(ix,1,1)-1) ) ) + 0.2*ng(ix,1,1)
                   yldot(iv1) = nurlxn *(nbound - ng(ix,0,1))/n0(ifld)
                 else  # old condition before 2/21/99
                   yldot(iv1) = nurlxn*(ni(ix,1,ifld)-ni(ix,0,ifld))/
     .                                                    n0(ifld)
                 endif
               else   # ix is not part of the core boundary
                  if (iscpli(ix) .eq. 1) call wsmodi(1)
		  t0 = max(tg(ix,0,1),tgmin*ev)
                  vyn = sqrt( 0.5*t0/(pi*mi(1)) )
                  fng_chem = 0.
		  do ii = 1, ngsp  #chem sputt of hydrogen - strange = 0
                     flx_incid = ng(ix,1,ii)*0.25*sqrt(8*t0/(pi*mg(ii)))
                     fng_chem= fng_chem + chemsputi(1,ii)*flx_incid*
     .                                                        sy(ix,0)
                  enddo
                  nharmave = 2.*(ni(ix,0,ifld)*ni(ix,1,ifld)) /
     .                          (ni(ix,0,ifld)+ni(ix,1,ifld))
                  fng_alb = (1-albedoi(ix,1))*nharmave*vyn*sy(ix,0)
                  yldot(iv1) = -nurlxg*( fniy(ix,0,ifld) + fng_alb -
     .                                   fng_chem ) / 
     .                                        (vyn*sy(ix,0)* n0(ifld))
cc                  if (fng_chem .ne. 0.) yldot(iv1) = -nurlxg*(
cc     .               fniy(ix,0,ifld) - fng_chem )/(vyn*sy(ix,0)*n0(ifld))
c...   Caution: the wall source models assume gas species 1 only is inertial
                  if(matwalli(ix) .gt. 0) then
                    if (recycwit(ix,1,1) .gt. 0.) then  
                      fniy_recy = recycwit(ix,1,1)*fac2sp*fniy(ix,0,1)
                      if (isrefluxclip==1) fniy_recy=min(fniy_recy,0.)
                      yldot(iv1)=-nurlxg*( fniy(ix,0,ifld) + fniy_recy -
     .                            fngyi_use(ix,1) - fngysi(ix,1) + fng_alb - 
     .                            fng_chem ) / (vyn*n0(ifld)*sy(ix,0))
                    elseif (recycwit(ix,1,1) < -1) then
                      yldot(iv1)=nurlxg*(ngbackg(1)-ni(ix,0,ifld))/n0(ifld)
                    elseif (recycwit(ix,1,1) .le. 0.) then  # treat recycwit as albedo
                      nharmave = 2.*(ni(ix,0,ifld)*ni(ix,1,ifld)) /
     .                              (ni(ix,0,ifld)+ni(ix,1,ifld))
                      yldot(iv1) = -nurlxg*( fniy(ix,0,ifld) +
     .                   (1+recycwit(ix,1,1))*nharmave*vyn*sy(ix,0) )/
     .                                             (vyn*n0(ifld)*sy(ix,0))
                    endif                  
                  endif
                  if(fngysi(ix,1)+fngyi_use(ix,1) .ne. 0. 
     .                                           .and. matwalli(ix)==0.)
     .               yldot(iv1)= -nurlxg*(fniy(ix,0,ifld) - fngysi(ix,1) -
     .                            fngyi_use(ix,1) ) / (vyn*sy(ix,0)* n0(ifld))
               endif   # end if-test for core and p.f. boundaries
             else  # ifld is NOT inertial neutrals but still hydrogen ion
c
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  if (isnicore(ifld) .eq. 1) then # set uniform, fixed density
                     yldot(iv1) = nurlxn*(ncore(ifld)-ni(ix,0,ifld))/
     .                                                          n0(ifld)
                  elseif (isnicore(ifld) .eq. 0) then # set flux to curcore
                     yldot(iv1) = - nurlxn* ( qe*
     .                       (fniy(ix,0,ifld)-fniycbo(ix,ifld))/sy(ix,0)
     .                            - curcore(ifld)*gyf(ix,0)/sygytotc ) /
     .                                              (qe*vpnorm*n0(ifld))
                  elseif (isnicore(ifld) .eq. 2) then # set flux & ni over range
                     yldot(iv1) = - nurlxn* ( qe*
     .                       (fniy(ix,0,ifld)-fniycbo(ix,ifld))/sy(ix,0)
     .                            - curcore(ifld)*gyf(ix,0)/sygytotc ) /
     .                                              (qe*vpnorm*n0(ifld))
                     if (ix.ge.ixfixnc .and. ix.le.ixfixnc+incixc) then
                        yldot(iv1) = nurlxn*(ncore(ifld)-ni(ix,0,ifld))/
     .                                                           n0(ifld)
                     endif                       
		  elseif (isnicore(ifld) .eq. 3) then # const ni; flux set to
                                                      # curcore-recycc*fngy
                     yldot(iv1) = -nurlxn*( ni(ix,0,ifld) -
     .                                ni(ixp1(ix,0),0,ifld) ) / n0(ifld)
                     if (ix.eq.ix_fl_bc) then
                        ii=max(0,ixpt1(1)+1)
                        fniytotc=fniy(ii,0,ifld)-fniycbo(ii,ifld)
                        fngytotc = fngy(ii,0,1)  # needs generalization
                        do
                          ii=ixp1(ii,0)
                          fniytotc = fniytotc + fniy(ii,0,ifld) -
     .                                          fniycbo(ii,ifld)
                          fngytotc = fngytotc + fngy(ii,0,1) #generalize
                          if (ii==ix_fl_bc) break
                        enddo
			if (ifld > 1) fngytotc = 0.  # crude fix for imp.
                        yldot(iv1)= -nurlxn*(qe*fniytotc +
     .                                       qe*recycc(1)*fngytotc -
     .                               curcore(ifld))/(qe*vpnorm*n0(ifld))
                     endif
                  elseif (isnicore(ifld)==5) then # set dni/dy=-ni/lynicore
                                                  # at midpl; the pol const ni
                     if (ix .ne. ixmp) then  #  flux func for ni
                       yldot(iv1) = -nurlxn*( ni(ix,0,ifld) -
     .                                 ni(ixp1(ix,0),0,ifld) ) / n0(ifld)
                     else  # close by setting midp gradient 1/lynicore
                       yldot(iv1)=-nurlxn*( ni(ix,0,ifld) - ni(ix,1,ifld)*
     .                                 (2*gyf(ixmp,0)*lynicore(ifld)+1)/
     .                                 (2*gyf(ixmp,0)*lynicore(ifld)-1)- 
     .                                     ncoremin(ifld) ) / n0(ifld)
                     endif
                  else
		            call xerrab ('** isnicore value not valid option **')
                  endif
               elseif (isnwconiix(ix,ifld) .eq. 0) then
               # ix not on core boundary; set zero gradient (or flux)
                  yldot(iv1) = nurlxn * ( (1-ifluxni)*
     .                      (niy1(ix,0,ifld) - niy0(ix,0,ifld))
     .                 - ifluxni*( fniy(ix,0,ifld)/(sy(ix,0)*vpnorm)    
     .                 - 0.001*ni(ix,1,ifld)*vy(ix,0,ifld)/vpnorm 
     .                 ))/n0(ifld)
c...  the last term going as 0.001 is to prevent very small densities
               elseif (isnwconiix(ix,ifld) .eq. 1) then
               # ix is not part of the core boundary; set fixed density
                  yldot(iv1) = nurlxn*(nwalli(ix) - ni(ix,0,ifld))/
     .                                                        n0(ifld)
               elseif (isnwconiix(ix,ifld) .eq. 2) then
               # ix is not part of the core boundary; set extrapolation bc
                  nbound =  ni(ix,1,ifld) - gyf(ix,1)*
     .                         (ni(ix,2,ifld)-ni(ix,1,ifld))/gyf(ix,0)
                  nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                          ni(ix,1,ifld)-1) ) ) + 0.2*ni(ix,1,ifld)
                  yldot(iv1) = nurlxn *(nbound - ni(ix,0,ifld))/n0(ifld)
               elseif (isnwconiix(ix,ifld) .eq. 3) then  #spec. gradient
                    yldot(iv1) =-nurlxn*( niy0(ix,0,ifld) -
     .                niy1(ix,0,ifld)*(2*gyf(ix,0)*lyniix(1,ix,ifld)-1)/
     .                                (2*gyf(ix,0)*lyniix(1,ix,ifld)+1) -
     .                                  nwimin(ifld) ) / n0(ifld)
               endif   # end if-test for core and p.f. boundaries

             endif           # endif for if (isupgon.eq.1 ...
           endif             # end for if-test on isnionxy
  272    continue  # ix-loop for ion dens

c...  Reset density at corners
         do jx = 1, nxpt
           if(isfixlb(jx).ne.2 .and. ixmnbcl*iymnbcl.eq.1 .and.
     .                       isnionxy(ixlb(jx),0,ifld)==1) then
	     yldot(idxn(ixlb(jx),0,ifld)) = nurlxn * 
     .            ( ave(ni(ixlb(jx),1,ifld),ni(ixlb(jx)+1,0,ifld))
     .                        - ni(ixlb(jx),0,ifld) ) / n0(ifld)
           endif
           if (isfixrb(jx).ne.2 .and. ixmxbcl*iymnbcl.eq.1 .and.
     .                      isnionxy(ixrb(jx)+1,0,ifld)==1) then
	      yldot(idxn(ixrb(jx)+1,0,ifld)) = nurlxn * 
     .           ( ave(ni(ixrb(jx)+1,1,ifld),ni(ixrb(jx),0,ifld))
     .                       - ni(ixrb(jx)+1,0,ifld) ) / n0(ifld)
           endif
         enddo
      
  274 continue  # ifld loop over ion species

c...  Now do the parallel velocity BC at iy = 0 for full mom eqn species
      do ifld = 1, nusp
         do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  #large ix-loop for up
            if (isuponxy(ix,0,ifld) .eq. 1) then
               ix2 = ixp1(ix,0)
               iv2 = idxu(ix,0,ifld)
               if (isixcore(ix)==1) then   # ix is part of the core boundary
                  upbdry = up(ix,0,ifld)
                  upbdry1 = up(ix,1,ifld)
                  upbdry2 = up(ix,2,ifld)
	          if (isupcore(ifld) == 0) then
		    yldot(iv2) = nurlxu * (upcore(ifld) - upbdry)/vpnorm
		  elseif (isupcore(ifld) ==1) then # d(up)/dy=0
                    yldot(iv2) = nurlxu * (upbdry1 - upbdry) / vpnorm
		  elseif(isupcore(ifld) ==2) then  # d^2(up)/dy^2 = 0
		    yldot(iv2) = nurlxu * ( (upbdry1-upbdry )*gy(ix,1) -
     .                                    (upbdry2-upbdry1)*gy(ix,2) )/
     .                                                 (gy(ix,1)*vpnorm)
		  elseif(isupcore(ifld)==3) then # set fmiy=0 on core bdry
                     yldot(iv2) = -nurlxu * fmiy(ix,0,ifld) / 
     .                                    (vpnorm*sy(ix,0)*fnorm(ifld))
		  elseif (isupcore(ifld)==4 .or. isupcore(ifld)==5) then 
                      # n*uz/rm=const & for =4- radial flux=lzflux; 
                      # or, for =5- ave uz=utorave 
                    yldot(iv2)=-nurlxe*( uz(ix,0,ifld)-uz(ixp1(ix,0),0,ifld)*
     .                           rm(ix,0,2)*ni(ixp1(ix,0),0,ifld)/
     .                          (rm(ixp1(ix,0),0,2)*ni(ix,0,ifld)) )/vpnorm
                    if (ix.eq.ix_fl_bc) then # not all Jac elems included
                      # sum core parallel momentum flux:
                      ii = max(0, ixpt1(1)+1)
                      fmiytotc = rm(ii,0,0)*fmiy(ii,0,ifld)
                      sytotc = sy(ii,0)
                      uztotc = uz(ii,0,ifld)/gxf(ii,0)
                      xtotc = dx(ii,0)
                      do    # loop over ii as changed by ii=ixp1 statement
                        ii = ixp1(ii,0)
                        if (isudsym==1 .and. (ii==nxc .or. ii==nxc+1)) then
                          continue  # skip m.p. guard cells
                        else
                          fmiytotc=fmiytotc+ rm(ii,0,0)*fmiy(ii,0,ifld)
                          sytotc = sytotc + sy(ii,0)
                          uztotc = uztotc + uz(ii,0,ifld)/gxf(ii,0)
                          xtotc = xtotc + dx(ii,0)
                        endif
                        if (ii==ix_fl_bc) break
                      enddo
                      uztotc = uztotc/xtotc   # divide by total core distance
                      if (isupcore(ifld)==4) then
                         yldot(iv2)= -nurlxe*( fmiytotc/sytotc -
     .                                                 lzflux(ifld) )/
     .                                (vpnorm*fnorm(ifld)*rm(ixmp,0,0))
                      elseif (isupcore(ifld)==5) then
                         yldot(iv2)= -nurlxe*(uztotc - utorave(ifld))/vpnorm
                      endif
                    endif
                  else   # isupcore does not correspond to cases above
                    call xerrab ("*** Illegal setting of isupcore ***")
                  endif

               elseif (isupwiix(ix,ifld)==1) then   # isixcore if-test; PF bdry
                  yldot(iv2) = -nurlxu * fmiy(ix,0,ifld) / 
     .                                    (vpnorm*sy(ix,0)*fnorm(ifld))
               elseif (isupwiix(ix,ifld)==2) then   # PF bdry
                  yldot(iv2) = nurlxu * nm(ix,0,ifld) / fnorm(ifld) *
     .                               (up(ix,1,ifld) - up(ix,0,ifld))
               elseif (isupwiix(ix,ifld)==3) then   # PF bdry
                  yldot(iv2) = -nurlxu * nm(ix,0,ifld) / fnorm(ifld) *
     .                         ( up(ix,0,ifld) - up(ix,1,ifld)*
     .                            (2*gyf(ix,0)*lyup(1)-1)/
     .                                     (2*gyf(ix,0)*lyup(1)+1) )
               else                          # PF bdry
                  yldot(iv2) = nurlxu * nm(ix,0,ifld) / fnorm(ifld) *
     .                                      (0. - up(ix,0,ifld))
               endif   # end if-test for core and p.f. boundaries
            endif      # endif for if (isuponxy .eq. 1)
	 enddo         # ix-loop for up eqns
      enddo

c ... Do boundary conditions for impurities along iy = 0.
c     Force corner cells (0,0) and (nx+1,0) to relax to an average of
c     the adjacent cells.
      do ifld = 1, nzspt
         iimp = nhsp + ifld
         if (isimpon .ge. 3 .and. isimpon .le. 7) then 
           do ix = i2, i5   # ix-loop for ion dens
             if (isnionxy(ix,0,iimp)==1) then
               iv = idxn(ix,0,iimp)
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  if (isnicore(iimp) .eq. 1) then # set uniform, fixed density
                     yldot(iv) = nurlxn*(ncore(iimp)-ni(ix,0,iimp))/
     .                                                         n0(iimp)
                  elseif (isnicore(iimp) .eq. 0) then # set flux to curcore
                     yldot(iv) = -nurlxn * ( qe*
     .                           (fniy(ix,0,iimp)-fniycbo(ix,iimp))/sy(ix,0)
     .                             - curcore(iimp)*gyf(ix,0)/sygytotc ) /
     .                                               (qe*vpnorm*n0(iimp))
                  elseif (isnicore(iimp) .eq. 2) then # set flux & ni over range
                     yldot(iv) = -nurlxn * ( qe*
     .                           (fniy(ix,0,iimp)-fniycbo(ix,iimp))/sy(ix,0)
     .                             - curcore(iimp)*gyf(ix,0)/sygytotc ) /
     .                                               (qe*vpnorm*n0(iimp))
                     if (ix.ge.ixfixnc .and. ix.le.ixfixnc+incixc) then
                        yldot(iv) = nurlxn*(ncore(iimp)-ni(ix,0,iimp))/
     .                                                         n0(iimp)
                     endif
                  elseif (isnicore(iimp) .eq. 3) then # set integr. flux, const ni 
                     yldot(iv1) = nurlxn*( ni(ix,0,iimp) -
     .                                 ni(ixp1(ix,0),0,iimp) ) / n0(iimp)
                     if (ix.eq.ix_fl_bc) then
                        ii=max(0,ixpt1(1)+1)
                        fniytotc=fniy(ii,0,iimp)-fniycbo(ii,iimp)
                        do
                          ii=ixp1(ii,0)
                          fniytotc=fniytotc+fniy(ii,0,iimp)-fniycbo(ii,iimp)
                          if (ii==ix_fl_bc) break
                        enddo
                        yldot(iv1)= -nurlxn*(qe*fniytotc-
     .                               curcore(iimp))/(qe*vpnorm*n0(iimp))
                     endif
                  elseif (isnicore(iimp) .eq. 4) then  # use bndry impur source
                     yldot(iv) = -nurlxn *
     .                           (fniy(ix,0,iimp) - fnzysi(ix,ifld)) /
     .                           (sy(ix,0) * n0(iimp) * vpnorm)
                  elseif (isnicore(iimp) == 5) then # set dni/dy=0
                     yldot(iv1) = nurlxn*( ni(ix,1,iimp) -
     .                                       ni(ix,0,iimp) ) / n0(ifld)
                  endif
               else   # ix is not part of the core boundary:
                  if (isnwconiix(ix,iimp) .eq. 0) then
                     yldot(iv) = -nurlxn *
     .                           (fniy(ix,0,iimp) - fnzysi(ix,ifld)) /
     .                           (sy(ix,0) * n0(iimp) * vpnorm)
                  elseif (isnwconiix(ix,iimp) .eq. 1) then
                     yldot(iv) = nurlxn*(nwalli(ix) - ni(ix,0,iimp))/
     .                                                        n0(iimp)
                  elseif (isnwconiix(ix,iimp) .eq. 2) then  #extrapolation bc
                     nbound =  ni(ix,1,iimp) - gyf(ix,1)*
     .                         (ni(ix,2,iimp)-ni(ix,1,iimp))/gyf(ix,0)
                     nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                          ni(ix,1,iimp)-1) ) ) + 0.2*ni(ix,1,iimp)
                     yldot(iv) = nurlxn *(nbound - ni(ix,0,iimp))/n0(iimp)
                  elseif (isnwconiix(ix,iimp) .eq. 3) then  #spec. gradient
                     yldot(iv) = -nurlxn*( niy0(ix,0,iimp) -
     .                 niy1(ix,0,iimp)*(2*gyf(ix,0)*lyniix(1,ix,iimp)-1)/
     .                                 (2*gyf(ix,0)*lyniix(1,ix,iimp)+1) -
     .                                     nwimin(iimp) ) / n0(iimp)
                  endif
               endif   # end if-test for core and p.f. boundaries
             endif     # end if-test for isnionxy
	   enddo       # end for do ix = i2, i5 loop for ion dens
c...   set the corner values here
            if (ixmnbcl*iymnbcl.eq.1) then
               do jx = 1, nxpt
		  if (isnionxy(ixlb(jx),0,iimp)==1) then
                    iv = idxn(ixlb(jx),0,iimp)
                    yldot(iv) = nurlxn *
     .                 ( ave(ni(ixlb(jx),1,iimp),ni(ixlb(jx)+1,0,iimp))
     .                           - ni(ixlb(jx),0,iimp) ) / n0(iimp)
                  endif
               enddo
            endif
            if (ixmxbcl*iymnbcl.eq.1) then
               do jx = 1, nxpt
	         if (isnionxy(ixrb(jx)+1,0,iimp)==1) then
                   iv = idxn(ixrb(jx)+1,0,iimp)
                   yldot(iv) = nurlxn *
     .                 ( ave(ni(ixrb(jx)+1,1,iimp),ni(ixrb(jx),0,iimp))
     .                           - ni(ixrb(jx)+1,0,iimp) ) / n0(iimp)
                 endif
               enddo
            endif

         endif
      enddo          # end of impurities along iy = 0

CCC  Do Te, Ti & later diffusive neutrals, phi along iy=0 boundary
ccc  -------------------------------------------------------------------
c
c    if flux energy condition, precompute feeytotc & feiytotc for use
      if (iflcore==1) then  #sum e-i core pwr fee,iytotc; unnecess if te,i frozen
        feeytotc = 0.
        feiytotc = 0.
        do ix = i4+1-ixmnbcl, i8-1+ixmxbcl
          if (iymnbcl==1 .and. isixcore(ix)==1) then  # domain part of core bdry
            ii = max(0, ixpt1(1)+1)
            feeytotc = feey(ii,0)-feeycbo(ii)
            feiytotc = feiy(ii,0)-feiycbo(ii)
            do    # loop over ii as changed by ii=ixp1 statement
              ii = ixp1(ii,0)
              if (isudsym==1 .and.(ii==nxc .or. ii==nxc+iexclnxc1)) then
                 continue  # skip m.p. guard cells
              else
                 feeytotc = feeytotc+feey(ii,0)-feeycbo(ii)
                 feiytotc = feiytotc+feiy(ii,0)-feiycbo(ii)
              endif
              if (ii==ix_fl_bc) break
            enddo
          endif
        enddo
      endif
c ... Sum for serial case
      ispwrbcl = 0
      if (ndomain == 1) then
        sfeeytotc = feeytotc
        sfeiytotc = feiytotc
        ispwrbcl = 1
      else
        if (ispwrbc(mype+1) == 1) ispwrbcl = 1
      endif

c     Sum for parallel case, sum/distribute to all processors
c_mpi      call mpi_allreduce(feeytotc, sfeeytotc, 1,
c_mpi     .      MPI_DOUBLE_PRECISION, MPI_SUM, uedgeComm, ierr)
c_mpi      call mpi_allreduce(feiytotc, sfeiytotc, 1,
c_mpi     .      MPI_DOUBLE_PRECISION, MPI_SUM, uedgeComm, ierr)

c ... Set Te and Ti BCs
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl #long ix-loop for Te, Ti
         if (isteonxy(ix,0) .eq. 1) then
           iv1 = idxte(ix,0)
            if (isixcore(ix)==1) then   # ix is part of the core boundary:
              yldot(iv1)=nurlxe*(tcoree*ev-te(ix,0))*1.5*ne(ix,0)/ennorm
              if (iflcore .eq. 1) then
                 yldot(iv1) = -nurlxe*(te(ix,0)-te(ixp1(ix,0),0))*n0(1)/ennorm
                 if (ix.eq.ix_fl_bc) then # not all Jac elems included
                    # sum core power:
                    ii = max(0, ixpt1(1)+1)
                    feeytotc = feey(ii,0)-feeycbo(ii)
                    do    # loop over ii as changed by ii=ixp1 statement
                       ii = ixp1(ii,0)
                       if (isudsym==1 .and. (ii==nxc .or. ii==nxc+iexclnxc1)) then
                          continue  # skip m.p. guard cells
                       else
                          feeytotc = feeytotc+feey(ii,0)-feeycbo(ii)
                       endif
                       if (ii==ix_fl_bc) break
                    enddo
                    yldot(iv1)= -nurlxe*(feeytotc-pcoree)/(vpnorm*ennorm)
                 endif
              elseif (iflcore .eq. -1) then
                 yldot(iv1) = -nurlxe*(te(ix,0)-te(ix,1))*n0(1)/ennorm
              endif
            else   # ix is not part of the core boundary; various PF cases:
               if (istepfcix(ix) .eq. 0) then          # zero electron energy flux
                 yldot(iv1) = -nurlxe*(feey(ix,0)/(n0(1)*vpnorm*sy(ix,0)))
     .                                                       / (temp0*ev) 
               elseif (istepfcix(ix) .eq. 1) then      # fixed Te
                 yldot(iv1) =nurlxe*(tewalli(ix)*ev-te(ix,0))/(temp0*ev)
               elseif (istepfcix(ix) .eq. 2) then       # extrapolation
                 tbound = te(ix,1) - gyf(ix,1)*
     .                                     (te(ix,2)-te(ix,1))/gyf(ix,0)
                 tbound = max(tbound, tbmin*ev)
                 yldot(iv1) = nurlxe *(tbound - te(ix,0))/(temp0*ev)
               elseif (istepfcix(ix) .eq. 3) then         # specified gradient
                 yldot(iv1) = nurlxe*( (te(ix,1) - te(ix,0)) -
     .                             0.5*(te(ix,1) + te(ix,0))/
     .                              (gyf(ix,0)*lyteix(1,ix)) )/(temp0*ev)
               elseif (istepfcix(ix) .eq. 4) then      
                 yldot(iv1) = -nurlxe*(feey(ix,0) - bceew*ne(ix,0)*
     .                               vey(ix,0)*sy(ix,0)*te(ix,0)) /
     .                              (vpnorm*ennorm*sy(ix,0))
               endif
            endif   # end if-test for core and p.f. boundaries
         endif   # end if-test on isteonxy

         if (istionxy(ix,0) .eq. 1) then
            iv2 = idxti(ix,0)
            if (isixcore(ix)==1) then   # ix is part of the core boundary:
              yldot(iv2)=nurlxi*(tcorei*ev-ti(ix,0))*1.5*ne(ix,0)/ennorm
              if (iflcore .eq. 1) then
                 yldot(iv2) = -nurlxi*(ti(ix,0)-ti(ixp1(ix,0),0))*
     .                                                     n0(1) /ennorm
                 if (ix.eq.ix_fl_bc) then # not all Jac elems included
                    # sum core power:
                    ii = max(0, ixpt1(1)+1)
                    feiytotc = feiy(ii,0)-feiycbo(ii)
                    do
                       ii = ixp1(ii,0)
                       if (isudsym==1 .and. (ii==nxc .or. ii==nxc+iexclnxc1)) then
                          continue  # skip m.p. guard cells
                       else
                          feiytotc = feiytotc+feiy(ii,0)-feiycbo(ii)
                       endif
                       if (ii==ix_fl_bc) break
                    enddo
                    yldot(iv2)= -nurlxi*(feiytotc-pcorei)/(vpnorm*ennorm)
                 endif
              elseif (iflcore .eq. -1) then
                 yldot(iv2) = -nurlxe*(ti(ix,0)-ti(ix,1))*n0(1)/ennorm
              endif
ccc            elseif (matwalli(ix) .gt. 0) then
ccc               bcen = matwalli(ix)
ccc               yldot(iv2) = -nurlxi*(feiy(ix,0) -
ccc     .             (bcei*fniy(ix,0,1) + bcen*fniy(ix,0,iigsp))*ti(ix,0)) /
ccc     .             (vpnorm*ennorm*sy(ix,0))
            else                                 # various PF cases
               if (istipfcix(ix) .eq. 0) then          # zero ion energy flux
                 yldot(iv2) = -nurlxi*(feiy(ix,0)/(n0(1)*vpnorm*sy(ix,0)))
     .                                                       / (temp0*ev) 
               elseif (istipfcix(ix) .eq. 1) then      # fixed Ti
                 yldot(iv2) =nurlxi*(tiwalli(ix)*ev-ti(ix,0))/(temp0*ev)
               elseif (istipfcix(ix) .eq. 2) then      # extrapolation
                 tbound = ti(ix,1) - gyf(ix,1)*
     .                                     (ti(ix,2)-ti(ix,1))/gyf(ix,0)
                 tbound = max(tbound, tbmin*ev)
                 yldot(iv2) = nurlxi *(tbound - ti(ix,0))/(temp0*ev)
               elseif (istipfcix(ix) .eq. 3) then      # specified gradient
                 yldot(iv2) = nurlxi*( (ti(ix,1) - ti(ix,0)) -
     .                           0.5*(ti(ix,1) + ti(ix,0))/
     .                             (gyf(ix,0)*lytiix(1,ix)) )/(temp0*ev)
               elseif (istipfcix(ix) .eq. 4) then      # sheath-like condition
                 t0 = max(tg(ix,0,1), temin*ev)
                 fngyw=0.25*sqrt(8*t0/(pi*mg(1)))*ng(ix,1,1)*sy(ix,0)
                 yldot(iv2) = -nurlxi*( feiy(ix,0) - bceiw*fniy(ix,0,1)*
     .                                                        ti(ix,0) -
     .                     cftiexclg*( bcenw*fniy(ix,0,iigsp)*tg(ix,0,1) -
     .                                cgengw*2.*tg(ix,0,iigsp)*fngyw ) ) /
     .                                          (vpnorm*ennorm*sy(ix,0))
               endif
            endif   # end if-test for core and p.f. boundaries
         endif   # end if-test on istionxy
      enddo  # long ix-loop for Te and Ti BCs

ccc  Diffusive neutrals on iy=0 boundary
ccc  - - - - - - - - - - - - -
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl #long ix-loop for diff neut dens
         nzsp_rt = nhsp
         do 275 igsp = 1, ngsp
	   jz = max(igsp - nhgsp, 1)  # identify impurity index
	   if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fniy
           if (isngonxy(ix,0,igsp) .eq. 1) then  # ends just before 275 continue statem.
            if (iscpli(ix) .eq. 1) call wsmodi(igsp)
            iv = idxg(ix,0,igsp)   
	    t0 = max(cdifg(igsp)*tg(ix,0,igsp),tgmin*ev)
            vyn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
            t1 = engbsr * max(tg(ix,0,1),tgmin*ev)
c ... prepare impurity ion flux for possible recycling
            zflux = 0.
            if (igsp .gt. nhgsp) then
              do iimp = 1, nzsp(jz)
		zflux = zflux + min(fniy(ix,0,nzsp_rt+iimp), 0.)
              enddo
	    endif
c ... Set core BC
            if (isixcore(ix)==1) then   # ix is part of the core boundary:
               if (isngcore(igsp) .eq. 0) then
                 nharmave = 2.*(ng(ix,0,igsp)*ng(ix,1,igsp)) /
     .                         (ng(ix,0,igsp)+ng(ix,1,igsp))
                 fng_alb = (1-albedoc(igsp))*
     .                                 nharmave*vyn*sy(ix,0)
  
                 yldot(iv) = -nurlxg*(fngy(ix,0,igsp) + fng_alb)/
     .                                    (vyn*sy(ix,0)*n0g(igsp))
               elseif (isngcore(igsp).eq.1) then
                 yldot(iv) = nurlxg*(ngcore(igsp)-ng(ix,0,igsp)) /
     .                                                       n0g(igsp)
               elseif (isngcore(igsp) .eq. 2) then
                 lengg = sqrt(tg(ix,0,igsp)/
     .                  (mg(igsp)*(nuix(ix,0,igsp)*nuiz(ix,0,igsp))) )
                 yldot(iv) = nurlxn*(
     .                         (ng(ix,1,igsp) - ng(ix,0,igsp)) -
     .                        0.5*(ng(ix,1,igsp) + ng(ix,0,igsp))/
     .                                 (gyf(ix,0)*lengg) )/n0g(igsp)
               elseif (isngcore(igsp) .eq. 3) then
                 nbound =  ng(ix,1,igsp) - gyf(ix,1)*
     .                        (ng(ix,2,igsp)-ng(ix,1,igsp))/gyf(ix,0)
                 nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                         ng(ix,1,igsp)-1) ) ) + 0.2*ng(ix,1,igsp)
                 yldot(iv) = nurlxn *(nbound - ng(ix,0,igsp))/n0g(igsp)
               elseif (isngcore(igsp) .eq. 4) then
                 yldot(iv) = nurlxn *(ng(ix,1,igsp) - ng(ix,0,igsp))/
     .                                                       n0g(igsp)
               endif

###  End of core bounary conditions 

            else  # cell is not part of core boundary
### wall BCs. calc diff neut sputt fluxes
               fng_chem = 0.
               do igsp2 = 1, 1+ishymol  #include hydrogen neut sputt
                 t0 = max(tg(ix,0,igsp2),tgmin*ev)
                 flx_incid=ng(ix,1,igsp2)*0.25*sqrt(8*t0/(pi*mg(igsp2)))
		 if (isch_sput(igsp).ne.0) then
                   call sputchem (isch_sput(igsp), t0/ev, tvwalli(ix),
     .                                       flx_incid, yld_carbi(ix))
                   fng_chem = fng_chem + fchemygwi(igsp2)*yld_carbi(ix)*
     .                                                flx_incid*sy(ix,0)
                 else
                   fng_chem = fng_chem + chemsputi(igsp,igsp2)* 
     .                                                flx_incid*sy(ix,0)
                 endif
               enddo

c ... Include gas BC from sputtering by ions
               sputflxpf(ix,igsp) = 0.
               if (isi_sputpf(igsp).ge.1 .and. igsp>1) then # incl phys sput from ions
                 do ifld = ipsputt_s, ipsputt_e 
                   eng_sput = ( 2*ti(ix,0) + zi(ifld)*
     .                                           3*te(ix,0) )/ev
                   if (zi(ifld) > 0.) sputflxpf(ix,igsp) = sputflxpf(ix,igsp) +
     .                                min(fniy(ix,0,ifld),0.)*fphysyiwi(igsp)*
     .                                      yld96(matp,matt,eng_sput)
                 enddo
                  if (isi_sputpf(igsp) .ge. 2) then # add chem sput from ions
                    do ifld = 1,1  # (ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                      eng_sput = ( 2*ti(ix,0) + zi(ifld)*
     .                                             3*te(ix,0) )/ev
		      flx_incid = -min(fniy(ix,0,ifld)/sy(ix,0),0.)
                      call sputchem (isch_sput(igsp), eng_sput, tvwalli(ix),
     .                                            flx_incid, yld_chm)
                      sputflxpf(ix,igsp) = sputflxpf(ix,igsp)+fchemyiwi(igsp)*
     .                                        min(fniy(ix,0,ifld),0.)*yld_chm
                    enddo
                  endif
               endif

 
               nharmave = 2.*(ng(ix,0,igsp)*ng(ix,1,igsp)) /
     .                       (ng(ix,0,igsp)+ng(ix,1,igsp))
               fng_alb = (1-albedoi(ix,igsp))*nharmave*vyn*sy(ix,0) 
               yldot(iv) = -nurlxg*( fngy(ix,0,igsp) + fng_alb -
     .                                   fng_chem + sputflxpf(ix,igsp) ) / 
     .                                        (vyn*sy(ix,0)*n0g(igsp))
               if(matwalli(ix) .gt. 0) then
                 if (recycwit(ix,igsp,1) .gt. 0.) then
                   fniy_recy = fac2sp*fniy(ix,0,1)
                   if(isrefluxclip==1) fniy_recy=min(fniy_recy,0.)
                   if (igsp .gt. nhgsp) fniy_recy = zflux
                   if (ishymol.eq.1 .and. igsp.eq.2) then # 2 atoms per molecule
                     if (isupgon(1) .eq. 1) then
		       fniy_recy = 0.5*( fniy(ix,0,1) + fniy(ix,0,2) )
                     else
                       fniy_recy = 0.5*( fniy(ix,0,1) + fngy(ix,0,1) )
                     endif
                     if(isrefluxclip==1) fniy_recy=min(fniy_recy,0.)
                   endif
                   yldot(iv) = -nurlxg*( fngy(ix,0,igsp) + fniy_recy*
     .                          recycwit(ix,igsp,1) - fngyi_use(ix,igsp) -
     .                                        fngysi(ix,igsp) + fng_alb -
     .                             fng_chem + sputflxpf(ix,igsp) ) / 
     .                                        (vyn*n0g(igsp)*sy(ix,0))
                 elseif (recycwit(ix,igsp,1) < -1) then
                   yldot(iv)=nurlxg*(ngbackg(igsp)-ng(ix,0,igsp))/
     .                                                          n0g(igsp)
                elseif (recycwit(ix,igsp,1) .le. 0.) then # treat recycwit as albedo
                   nharmave = 2.*(ng(ix,0,igsp)*ng(ix,1,igsp)) /
     .                           (ng(ix,0,igsp)+ng(ix,1,igsp))
                   yldot(iv) = -nurlxg*( fngy(ix,0,igsp) +
     .                     (1+recycwit(ix,igsp,1))*nharmave*vyn*
     .                            sy(ix,0) )/(vyn*n0g(igsp)*sy(ix,0))
                 endif 
               endif
               if(fngysi(ix,igsp)+fngyi_use(ix,igsp) .ne.0. 
     .                                            .and. matwalli(ix).eq.0.) 
     .                                            yldot(iv) = -nurlxg*
     .                             ( fngy(ix,0,igsp) - fngysi(ix,igsp) )
     .                                         / (vyn*sy(ix,0)*n0g(igsp))
            endif   # end if-test for core and p.f. boundaries
           endif    # End for if(isngon(igsp).eq.1) after do 275
 275     continue   # igsp loop over species
      enddo  # end of long ix-loop for diff neuts

c... BC for neutral gas temperature/energy at iy=0
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl # ix-loop for diff neut temp 
        do igsp = 1, ngsp
          if (istgonxy(ix,0,igsp) == 1) then
            iv = idxtg(ix,0,igsp)
            if (isixcore(ix)==1) then
              if (istgcore(igsp) == 0) then  #set to ti
                yldot(iv)=nurlxg*(ti(ix,0)*cftgticore(igsp) 
     .                            -tg(ix,0,igsp))/(temp0*ev)
              elseif (istgcore(igsp) == 1) then   # set to tgcore
                yldot(iv)=nurlxg*(tgcore(igsp)*ev-tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              elseif (istgcore(igsp) == 2) then 
                #if (isupgon(igsp)==1) then
                  t0 = max(tg(ix,0,igsp),tgmin*ev)
                  vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
                  nharmave = 2.*(ng(ix,0,igsp)*ng(ix,1,igsp)) /
     .                          (ng(ix,0,igsp)+ng(ix,1,igsp))
                  fng_alb = (1-albedoc(igsp))*
     .                                 nharmave*vyn*sy(ix,0)
                  yldot(iv)=-nurlxg*(fegy(ix,0,igsp) + cfalbedo*fng_alb*t0)
     .                                      /(vpnorm*ennorm*sy(ix,0))
                #else
                #  write(*,*) 'Istgcore',igsp,'==1 is not available'
                #endif
              else # all others, set to zero y-gradient
                yldot(iv)=nurlxg*(tg(ix,1,igsp)-tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              endif
            else   # PF wall
              if (istgpfcix(ix,igsp) == 0) then   # fixed Tg
                yldot(iv) = nurlxg*(tgwall(igsp)*ev-tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              elseif (istgpfcix(ix,igsp) == 1)    # extrapolation
                tbound = tg(ix,1,igsp) - gyf(ix,1)*
     .                      (tg(ix,2,igsp)-tg(ix,1,igsp))/gyf(ix,0)
                tbound = max(tbound, 0.25*tbmin*ev)  #tbmin=.1 eV
                yldot(iv) = nurlxi *(tbound - tg(ix,0,igsp))/(temp0*ev)
              elseif (istgpfcix(ix,igsp) == 2)    # specified gradient
                yldot(iv) = nurlxi*( (tg(ix,1,igsp) - tg(ix,0,igsp)) -
     .                           0.5*(tg(ix,1,igsp) + tg(ix,0,igsp))/
     .                           (gyf(ix,0)*lytg(1,igsp)) )/(temp0*ev)
              elseif (istgpfcix(ix,igsp) == 3)  #Maxwell thermal flux to wall
                t0 = max(cdifg(igsp)*tg(ix,1,igsp), temin*ev)
                vyn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                yldot(iv) =  -nurlxg*( fegy(ix,0,igsp) + 2*cgengmw*
     .                               ng(ix,1,igsp)*vyn*t0*sy(ix,0) )/
     .                                     (sy(ix,0)*vpnorm*ennorm)
              elseif (istgpfcix(ix,igsp) == 4) 
        	t0 = max(tg(ix,0,igsp),tgmin*ev)
                vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
                nharmave = 2.*(ng(ix,0,igsp)*ng(ix,1,igsp)) /
     .                        (ng(ix,0,igsp)+ng(ix,1,igsp))
                fng_alb = (1-albedoi(ix,1))*
     .                               nharmave*vyn*sy(ix,0)
                fng_chem = 0.  #..double check
                yldot(iv)=-nurlxg*(fegy(ix,0,igsp) + cfalbedo*fng_alb*t0
     .                                             - 2.*fng_chem*t0)
     .                                     /(vpnorm*ennorm*sy(ix,0))
                fniy_recy = 0.
                if (matwalli(ix) .gt. 0) then
                  if (recycwit(ix,igsp,1) .gt. 0) then
                    fniy_recy = recycwit(ix,igsp,1)*fac2sp*fniy(ix,0,1)
                    if (isrefluxclip==1) fniy_recy=min(fniy_recy,0.)
                    yldot(iv)=-nurlxg*(fegy(ix,0,igsp) + cfalbedo*fng_alb*t0
     .                                             - 2.*fng_chem*t0
     .                                             + fniy_recy*(1.-cfdiss)
     .                                              *cfalbedo*recycwe
     .                                              *ti(ix,0) )
     .                                     /(vpnorm*ennorm*sy(ix,0))
                  endif
                endif
              elseif (istgpfc(igsp) == 5) then   # tg=ti*cftgtipfc
                yldot(iv) = nurlxg*(ti(ix,0)*cftgtipfc(igsp)
     .                                        -tg(ix,0,igsp))/
     .                                                      (temp0*ev)
              elseif (istgpfc(igsp) > 5)
                 call xerrab("***Input error: invalid istgpfc ***")
              endif

            endif           
          endif
        enddo  
      enddo  # end ix-loop for neutral temps

ccc  Now do potential along iy=0 boundary if isnewpot=-1 or -2; others next page
ccc  - - - - - - - - - - - - - - 
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl #long ix-loop for potential
         if (isphionxy(ix,0) .eq. 1) then
            iv3 = idxphi(ix,0)
            yldot(iv3) = nurlxp*( (phi(ix,1) - phi(ix,0)) -
     .                             0.5*(phi(ix,1) + phi(ix,0))/
     .                              (gyf(ix,0)*lyphiix(1,ix)) )/temp0
            if (isnewpot==-1) then
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  ix4 = ixm1(ix,0)
                  ut0 = (rr(ix,0)+rr(ix,1))/(btot(ix,0)+btot(ix,1))*
     .                  ( ey(ix,0) - gpiy(ix,0,1)/
     .                         (qe*zi(1)*niy0(ix,0,1)) ) +
     .                        0.25*(rbfbt(ix,0)+rbfbt(ix,1))*
     .                              (up(ix,0,1)+up(ix4,0,1))
                  yldot(iv3) = -nurlxp*ut0/vpnorm
               endif   # end if-test for core and p.f. boundaries
            endif
            if (isnewpot==-2) then  # phi is constant on core boundary
               if (isixcore(ix)==1) then   # ix is part of the core boundary:
                  if (ix .ne. ix_fl_bc) then
                     yldot(iv3) = -nurlxp*(phi(ix,0)-phi(ixp1(ix,0),0))/temp0
                  elseif (ix==ix_fl_bc) then # total core current = icoreelec
                     fqytotc = 0.
                     ii=ixc1
                     fqytotc=fqyao(ii,0)
                     do
                        ii=ixp1(ii,0)
                        fqytotc = fqytotc + fqyao(ii,0)
                        if (ii==ix_fl_bc) break
                     enddo
                     yldot(iv3) = -nurlxp*(fqytotc - icoreelec)/
     .                                     (qe*n0(1)*vpnorm*sy(ixc1,0))
                  endif   # end if-test for integral BC on core
               endif   # end if-test for core and p.f. boundaries
            endif   # end if-test for isnewpot==-2
         endif

      enddo # long ix loop for potential
ccc  -------------------------------------------------------------------

c...  Now do the corners: the ion density corners are done elsewhere
       if ((xcnearlb.or.openbox) .and. ixmnbcl*iymnbcl.eq.1) then  
         do ifld = 1, nusp
            do jx = 1, nxpt
	     if(isuponxy(ixlb(jx),0,ifld)==1) then
               yldot(idxu(ixlb(jx),0,ifld)) = - nurlxu*
     .                      ( up(ixlb(jx),0,ifld) -
     .             0.5*(up(ixlb(jx),1,ifld)+up(ixlb(jx)+1,0,ifld)) )
     .                                         /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
	     if (isteonxy(ixlb(jx),0)==1) then
               yldot(idxte(ixlb(jx),0)) = nurlxe*
     .             ( 0.5*(te(ixlb(jx)+1,0)+te(ixlb(jx),1))
     .                               - te(ixlb(jx),0) ) / (temp0*ev)
             endif
         enddo
         do jx = 1, nxpt
	    if (istionxy(ixlb(jx),0)==1) then
               yldot(idxti(ixlb(jx),0)) = nurlxi*
     .             ( 0.5*(ti(ixlb(jx)+1,0)+ti(ixlb(jx),1))
     .                               - ti(ixlb(jx),0) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if (isngonxy(ixlb(jx),0,igsp)==1) then
               yldot(idxg(ixlb(jx),0,igsp)) = nurlxg*
     .           (ng(ixlb(jx)+1,0,igsp) - ng(ixlb(jx),0,igsp))/n0g(igsp)
             endif
           enddo
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if(istgonxy(ixlb(jx),0,igsp)==1) then
               yldot(idxtg(ixlb(jx),0,igsp)) = nurlxg*
     .          (tg(ixlb(jx)+1,0,igsp) - tg(ixlb(jx),0,igsp))/(temp0*ev)
             endif
           enddo
         enddo
       endif
       if ((xcnearrb.or.openbox) .and. ixmxbcl*iymnbcl.eq.1) then
         do ifld = 1, nusp
           do jx = 1, nxpt
	     if(isuponxy(ixrb(jx),0,ifld)==1) then
               yldot(idxu(ixrb(jx),0,ifld)) = - nurlxu*
     .                      ( up(ixrb(jx),0,ifld) -
     .             0.5*(up(ixrb(jx)-1,0,ifld)+up(ixrb(jx),1,ifld)) )
     .                                         /vpnorm
               yldot(idxu(ixrb(jx)+1,0,ifld)) = - nurlxu*
     .                   ( up(ixrb(jx)+1,0,ifld)-up(ixrb(jx),0,ifld) )
     .                                            /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
	    if(isteonxy(ixrb(jx)+1,0)==1) then
               yldot(idxte(ixrb(jx)+1,0)) = nurlxe*
     .             ( 0.5*(te(ixrb(jx)+1,1)+te(ixrb(jx),0))
     .                               - te(ixrb(jx)+1,0) ) / (temp0*ev)
            endif
         enddo
         do jx = 1, nxpt
	    if(istionxy(ixrb(jx)+1,0)==1) then
               yldot(idxti(ixrb(jx)+1,0)) = nurlxi*
     .             ( 0.5*(ti(ixrb(jx)+1,1)+ti(ixrb(jx),0))
     .                               - ti(ixrb(jx)+1,0) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if(isngonxy(ixrb(jx),0,igsp)==1) then
               yldot(idxg(ixrb(jx)+1,0,igsp)) = nurlxg*
     .           (ng(ixrb(jx),0,igsp) - ng(ixrb(jx)+1,0,igsp))/n0g(igsp)
             endif
           enddo
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
	     if(istgonxy(ixrb(jx),0,igsp)==1) then
               yldot(idxtg(ixrb(jx)+1,0,igsp)) = nurlxg*
     .          (tg(ixrb(jx),0,igsp) - tg(ixrb(jx)+1,0,igsp))/(temp0*ev)
             endif
           enddo
         enddo
       endif

      endif           # end of iy = 0 boundary conditions except for isnewpot=1

c  ####################################################################
c  ###  Now do the potential for isnewpot=1 case along iy=0 boundary
c  ####################################################################
c...  If isnewpot=1, phi boundary involves two eqns at iy=0 and 1
c...  Note: j3 is local range index for iy passed from pandf in oderhs.m
      if (isnewpot*isphion.eq.1 .and. j3.le.3) then
        do ix = min(i4+1-ixmnbcl,ixpt1(1)+1), max(i8-1+ixmxbcl,ixpt2(nxpt))
          if(isphionxy(ix,0)*isphionxy(ix,1)==1) then
            iv  = idxphi(ix,0)
            iv1 = idxphi(ix,1)

c  ################ do BC on core bdry, then on PF bdry #################
            if (isixcore(ix)==1 .and. isphicore0==0) then # ix is core bdry
               ix3 = ixm1(ix,1)
               ix4 = ixm1(ix,0)
               r_major = 0.5*(rm(ix,0,0) + rm(ix,1,0)) + r0slab
c         Set two B.C. (int. cur. and int. ut0 + pol. uniform phi(,0), te(,1)
c         Set phi poloidally constant at iy=0 on core bdry
                 yldot(iv) = -nurlxp*(phi(ix,0)-phi(ixp1(ix,0),0))/temp0

c         Next set phi at iy=1 depending on iphibcc flag

c  ######################################################################
c  ### First core phi BC at iy=0 & 1 over full ixcore range; then ix=ixmp
c  ######################################################################
cc     First check that iphibcc=1,2, or 3; othewise abort with message
               if (iphibcc < 1 .or. iphibcc > 3) then
                 call xerrab("**INPUT ERROR: only iphibcc=1,2,3 available")
               endif

                 if (iphibcc==1) then
                   yldot(iv1) = -nurlxp*( (ey(ix,1)-ey(ix,0))*gy(ix,1) -
     .                                  (ey(ix,2)-ey(ix,1))*gy(ix,2) )/
     .                                                  (gy(ix,1)*temp0)
                 elseif (iphibcc==2) then
			yldot(iv1)=-nurlxp*(te(ix,1)-te(ixp1(ix,1),1))/
     .                                                          (ev*temp0)
                 elseif (iphibcc==3) then
                   yldot(iv1)= -nurlxp*(phi(ix,1)-phi(ixp1(ix,1),1))/temp0
                 else  # no longer allowed; need to reconsider
                   yldot(iv1)= -nurlxp*( phi(ix,1)-phi(ixp1(ix,1),1) +
     .                                      dphi_iy1(ix) ) / temp0
                   if (isteon.eq.0) yldot(iv1)=nurlxp*
     .                                     (tes(ix,1)-te(ix,1))/temp0
                 endif

c  ####################################################################
c  ### Reset phi BC at iy=0,1 at ix=ixmp only using iphibcc & isutcore
c  ####################################################################
               if (ix == ixmp) then  #redefine eqn for midplane ix
                 fqytotc = 0.
                 uztotc = 0.
	         uztotc1 = 0.
                 uztotc2 = 0.
c NOTE: total poloidal distance around the core:
                 ii=ixc1
                 xtotc = dx(ii,0)
                 do
                    ii=ixp1(ii,1)
                    xtotc = xtotc + dx(ii,0)
                    if (ii==ix_fl_bc) break
                 enddo
                 ii=ixc1
c NOTE: bndry surface flux should exclude fqyb and include fqyd
                 fqytotc = fqya(ii,1)+ cfqyn*fqyn(ii,1) + cfqym*fqym(ii,1) 
     .                     + cfqybbo*fqyb(ii,1) + cfqydbo*fqyd(ii,1)
                 uztotc  = uz(ii,0,1) / ( gxf(ii,0)*xtotc )
                 uztotc1 = uz(ii,1,1) / ( gxf(ii,1)*xtotc )
                 uztotc2 = uz(ii,2,1) / ( gxf(ii,2)*xtotc )
                 do
                    ii=ixp1(ii,1)
                    fqytotc = fqytotc +
     .                        fqya(ii,1) + cfqyn*fqyn(ii,1) + cfqym*fqym(ii,1) 
     .                           + cfqybbo*fqyb(ii,1) + cfqydbo*fqyd(ii,1)
                    uztotc  = uztotc  + uz(ii,0,1) / ( gxf(ii,0)*xtotc )
                    uztotc1 = uztotc1 + uz(ii,1,1) / ( gxf(ii,1)*xtotc )
                    uztotc2 = uztotc2 + uz(ii,2,1) / ( gxf(ii,2)*xtotc )
                    if (ii==ix_fl_bc) break
                 enddo
c ... Midplane phi(ixmp,0) BC is always to set rad current to icoreelec
                 yldot(iv)=-nurlxp*(fqytotc - icoreelec)/
     .                                     (qe*n0(1)*vpnorm*sy(ixc1,0))

c ... Midplane phi(ixmp,1) BC set by iphibcc if =1,2,3; if iphibbc otherwise,
c ... isutcore sets phi(ixmp,1) BC
		 if (iphibcc.eq.1) then  # set d^2(ey)/dy^2=0
                   yldot(iv1) = -nurlxp*( (ey(ix,1)-ey(ix,0))*gy(ix,1) -
     .                                   (ey(ix,2)-ey(ix,1))*gy(ix,2) )/
     .                                                  (gy(ix,1)*temp0)
                 elseif (iphibcc==2 .or. iphibcc==3) then
                   yldot(iv1)=-nurlxp*(ey(ix,0)-eycore)/(gyf(ix,0)*temp0)
                 else  # any iphibcc > 3 uses isutcore BC at ix=ixmp
                   if (isutcore .eq. 0) then # flux-ave R*n*m*uz = lzcore
                     yldot(iv1)=-nurlxp*(uztotc -
     .                        lzcore(1)/(mi(1)*n0(1)*r_major)) / vpnorm
                   elseif (isutcore == 1) then  # zero radial deriv. of ave uz
                     yldot(iv1) = nurlxp*(uztotc-uztotc1)/vpnorm
		   else     # second deriv. of ey at the midplane
                     yldot(iv1) = -nurlxp*( (ey(ix,1)-ey(ix,0))*gy(ix,1) -
     .                                   (ey(ix,2)-ey(ix,1))*gy(ix,2) )/
     .                                                  (gy(ix,1)*temp0)
                   endif
                 endif
               endif  #end outer midplane cells ix=ixmp potential BCs
c  ####################################################################
c  ... Completed phi BC on core bdry
c  ####################################################################
c  ... Special case rarely used
               if (cfvycf .gt. 1e-20) then  # preserves PET97 cases
                 yldot(iv) = -nurlxp*( ey(ix,0) - gpiy(ix,0,1)/
     .                     (qe*zi(1)*niy0(ix,0,1)) )/(btot(ix,0)*vpnorm)
                 yldot(iv1) =nurlxp*(fqy(ix,1)-(fqx(ix,1)-fqx(ix3,1))) /
     .                              (rrv(ix,0)*sy(ix,0)*vpnorm*ev*n0(1))
               endif

c  ###################################################################
c  ### Finally set private-flux (PF) phi BC at iy=0 & 1
c  ###################################################################
            else   # ix is not part of the core boundary:
               if(iphibcwiix(ix) == 0) then
                 yldot(iv) = nurlxp*( phi(ix,1) - phi(ix,0) )/temp0
               elseif(iphibcwiix(ix) == 1) then
                 yldot(iv) = nurlxp*
     .                     (phintewi*te(ix,0)/ev - phi(ix,0))/temp0
	       elseif(iphibcwiix(ix) == 3) then
                 yldot(iv) = nurlxp*( (phi(ix,1) - phi(ix,0)) -
     .                             0.5*(phi(ix,1) + phi(ix,0))/
     .                               (gyf(ix,0)*lyphiix(1,ix)) )/temp0
               elseif(iphibcwiix(ix) == 4) then
                 yldot(iv) = nurlxp*
     .                     (phiwi(ix) - phi(ix,0))/temp0
               endif
            endif   # end if-test for core and p.f. boundaries
          endif     # end if-test on isphionxy for ix=0,1
        enddo       # end loop over ix
      endif         # large if testing isnewpot*isphion=1, j3
 1100 continue   # jump to here if iymnbcl = 0; interior bdry

c ====================================================================
c ======================== The iy=ny+1 boundary ======================
c ====================================================================

c...  now do the iy = ny+1 boundary
c...  if extrapolation b.c.on outer wall, isextrw=1, otherwise isextrw=0
      if (iymxbcl .eq. 0) goto 1200  #skip setting eqn because interior bdry
      if (j7 .ge. (ny+1)-isextrnw .or. j7 .ge. (ny+1)-isextrtw) then
      do 280 ifld = 1 , nisp
        do 278 ix = i4+1-ixmnbcl, i8-1+ixmxbcl
          if (isnionxy(ix,ny+1,ifld)==1) then
            iv1 = idxn(ix,ny+1,ifld)

c ...  First check if this is neutral species (zi=0) & do BC
c ---  Different boundary conditions if neutral momentum equation use;
c ---  typically hydrogen only, so DIVIMP chem sputt not used here
            if (isupgon(1) .eq. 1 .and. zi(ifld) .eq. 0.0) then
                if (iscplo(ix) .eq. 1) call wsmodo(1)
                t0 = max(tg(ix,ny+1,1),tgmin*ev)
                vyn = sqrt( 0.5*t0/(pi*mi(ifld)) )
                fng_chem = 0.
                do ii = 1, ngsp
                   flx_incid =ng(ix,ny+1,ii)*0.25*sqrt(8*t0/(pi*mg(ii)))
                   fng_chem = fng_chem + chemsputo(1,ii)*flx_incid*
     .                                                        sy(ix,ny)
                enddo
                nharmave = 2.*(ni(ix,ny,ifld)*ni(ix,ny+1,ifld)) /
     ,                        (ni(ix,ny,ifld)+ni(ix,ny+1,ifld))
                fng_alb = (1-albedoo(ix,1))*nharmave*vyn*sy(ix,ny)
                yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) - fng_alb + 
     .                             fng_chem ) / (vyn*sy(ix,ny)* n0(ifld))

c...   Caution: the wall source models assume gas species 1 only is inertial
                if(matwallo(ix) .gt. 0) then
                  if (recycwot(ix,1) .gt. 0.) then
                    fniy_recy = recycwot(ix,1)*fac2sp*fniy(ix,ny,1)
                    if(isrefluxclip==1) fniy_recy=max(fniy_recy,0.)
                    yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) + fniy_recy +
     .                 fngyo_use(ix,1)+fngyso(ix,1)-fng_alb + fng_chem ) / 
     .                                          (vyn*n0(ifld)*sy(ix,ny))
                  elseif (recycwot(ix,1) < -1) then
                    yldot(iv1)=nurlxg*(ngbackg(1)-ni(ix,ny+1,ifld))/
     .                                                           n0(ifld)
                  elseif (recycwot(ix,1) .le. 0.) then  # treat recycwot as albedo
                    yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) -
     .                     (1+recycwot(ix,1))*ni(ix,ny+1,ifld)*vyn*sy(ix,ny) )/
     .                                            (vyn*n0(ifld)*sy(ix,ny))
                  endif 
                endif
                if(fngyso(ix,1)+fngyo_use(ix,1).ne.0. .and. matwallo(ix)==0.)
     .                        yldot(iv1) = nurlxg*( fniy(ix,ny,ifld) + 
     .                                    fngyo_use(ix,1) + fngyso(ix,1) ) 
     .                                        / (vyn*sy(ix,ny)* n0(ifld))
            else

              if(isnwconoix(ix,ifld) .eq. 0) then
                 yldot(iv1) = nurlxn * ( (1-ifluxni)*
     .                         (niy0(ix,ny,ifld) - niy1(ix,ny,ifld))
     .            + ifluxni*( fniy(ix,ny,ifld)/(sy(ix,ny)*vpnorm) 
     .           - 0.001*ni(ix,ny,ifld)*vy(ix,ny,ifld)/vpnorm ) )
     .                                                      /n0(ifld)
c...  the last term going as 0.001 is to prevent very small densities
c...  This next if is a recalc. if constant ni wanted - could be done better
              elseif (isnwconoix(ix,ifld) .eq. 1) then
                 yldot(iv1) =nurlxn*(nwallo(ix) - ni(ix,ny+1,ifld))/
     .                                                      n0(ifld)
              elseif (isnwconoix(ix,ifld) .eq. 2) then  #extrapolation
                 nbound = ni(ix,ny,ifld) + gyf(ix,ny-1)*
     .                      (ni(ix,ny,ifld)-ni(ix,ny-1,ifld))/gyf(ix,ny)
                 nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                         ni(ix,ny,ifld)-1) ) ) + 0.2*ni(ix,ny,ifld)
ccc                 nbound = max(nbound, 0.3*ni(ix,ny,ifld))
                 yldot(iv1) = nurlxn*(nbound - ni(ix,ny+1,ifld))/n0(ifld)
              elseif (isnwconoix(ix,ifld) .eq. 3) then   #spec. gradient 
                  yldot(iv1) = -nurlxn*( niy1(ix,ny,ifld) -
     .              niy0(ix,ny,ifld)*(2*gyf(ix,ny)*lyniix(2,ix,ifld)-1)/
     .                               (2*gyf(ix,ny)*lyniix(2,ix,ifld)+1) -
     .                         nwomin(ifld) ) / n0(ifld)

              endif            # endif for if (isnwconoix.eq.1 .and. ..
            endif              # endif for neut. ni; i.e., if (isupgon .eq. 1
          endif                # endif for isnionxy
  278   continue               # end loop over ix for ni BC
c...  Reset density at corners
          if (ixmnbcl*iymxbcl.eq.1) then
            do jx = 1, nxpt
              if(isnionxy(ixlb(jx),ny+1,ifld)==1) then
	         yldot(idxn(ixlb(jx),ny+1,ifld)) = nurlxn * 
     .             ( ave(ni(ixlb(jx),ny,ifld),ni(ixlb(jx)+1,ny+1,ifld))
     .                        - ni(ixlb(jx),ny+1,ifld) ) / n0(ifld)
              endif
            enddo
          endif
          if (ixmxbcl*iymxbcl.eq.1) then
            do jx = 1, nxpt
              if(isnionxy(ixrb(jx)+1,ny+1,ifld)==1) then
	         yldot(idxn(ixrb(jx)+1,ny+1,ifld)) = nurlxn * 
     .               ( ave(ni(ixrb(jx)+1,ny,ifld),ni(ixrb(jx),ny+1,ifld))
     .                         - ni(ixrb(jx)+1,ny+1,ifld) ) / n0(ifld)
              endif
            enddo
          endif
  280 continue             # end loop over ifld for ni BC

c...  Do the parallel velocity BC along iy = ny+1

      do ifld = 1, nusp
         do ix = i4+1-ixmnbcl, i8-1+ixmxbcl
            if (isuponxy(ix,ny+1,ifld)==1) then
               iv2 = idxu(ix,ny+1,ifld)
               if (isupwoix(ix,ifld)==1) then  #zero parallel momentum flux
                  yldot(iv2) = nurlxu * fmiy(ix,ny,ifld) / 
     .                                 (vpnorm*sy(ix,ny)*fnorm(ifld))
               elseif (isupwoix(ix,ifld)==2) then  #  d(up)/dy = 0
                  yldot(iv2) = nurlxu * nm(ix,ny,ifld) / fnorm(ifld) *
     .                           (up(ix,ny,ifld) - up(ix,ny+1,ifld))
               elseif (isupwoix(ix,ifld)==3) then  # d(up)/dy = up/lyup
                  yldot(iv2) = -nurlxu * nm(ix,ny,ifld) / fnorm(ifld) *
     .                             ( up(ix,ny+1,ifld) - up(ix,ny,ifld)*
     .                                (2*gyf(ix,ny)*lyup(2)-1)/
     .                                      (2*gyf(ix,ny)*lyup(2)+1) )
               else
                  yldot(iv2) = nurlxu * nm(ix,ny,ifld) / fnorm(ifld) *
     .                                      (0. - up(ix,ny+1,ifld))
               endif
            endif 
         enddo
      enddo

c ... Do boundary conditions for impurities along iy = ny+1.
cc   - - - - - - - - - - - - - - - - - -
c     Force corner cells (0,ny+1) and (nx+1,ny+1) to relax to an average
c     of the adjacent cells.
      do ifld = 1, nzspt
        iimp = nhsp + ifld
        if (isimpon .ge. 3 .and. isimpon .le. 7) then
          do ix = i2, i5  # i2 and i5 limits omit the ix=0 and nx+1 corners
            if (isnionxy(ix,ny+1,iimp)==1) then
               iv = idxn(ix,ny+1,iimp)  
               if (isnwconoix(ix,iimp) .eq. 0) then       #fix flux
                  yldot(iv) = nurlxn *
     .                        (fniy(ix,ny,iimp) + fnzyso(ix,ifld)) /
     .                            (sy(ix,ny) * n0(iimp) * vpnorm)
               elseif (isnwconoix(ix,iimp) .eq. 1) then   #set to nwallo
                  yldot(iv) =nurlxn*(nwallo(ix) - ni(ix,ny+1,iimp))/
     .                                                      n0(iimp)
               elseif (isnwconoix(ix,iimp) .eq. 2) then   #extrapolation
                  nbound = ni(ix,ny,iimp) + gyf(ix,ny-1)*
     .                      (ni(ix,ny,iimp)-ni(ix,ny-1,iimp))/gyf(ix,ny)
                  nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                         ni(ix,ny,iimp)-1) ) ) + 0.2*ni(ix,ny,iimp)
                  yldot(iv) = nurlxn*(nbound - ni(ix,ny+1,iimp))/n0(iimp)
               elseif (isnwconoix(ix,iimp) .eq. 3) then   #spec. gradient 
                  yldot(iv) = -nurlxn*( niy1(ix,ny,iimp) -
     .              niy0(ix,ny,iimp)*(2*gyf(ix,ny)*lyniix(2,ix,iimp)-1)/
     .                               (2*gyf(ix,ny)*lyniix(2,ix,iimp)+1) -
     .                         nwomin(iimp) ) / n0(iimp)
               endif
            endif       # end if for isnionxy
          enddo
c...   set the corner values here
          if (ixmnbcl*iymxbcl.eq.1) then
             do jx = 1, nxpt
               if (isnionxy(ixlb(jx),ny+1,iimp)==1) then
                  iv = idxn(ixlb(jx),ny+1,iimp)
                  yldot(iv) = nurlxn *
     .              ( ave(ni(ixlb(jx),ny,iimp),ni(ixlb(jx)+1,ny+1,iimp))
     .                           - ni(ixlb(jx),ny+1,iimp) ) / n0(iimp)
               endif
             enddo
          endif
          if (ixmxbcl*iymxbcl.eq.1) then
             do jx = 1, nxpt
               if (isnionxy(ixrb(jx)+1,ny+1,iimp)==1) then
                  iv = idxn(ixrb(jx)+1,ny+1,iimp)
                  yldot(iv) = nurlxn *
     .                 ( ave(ni(ixrb(jx)+1,ny,iimp),ni(ixrb(jx),ny+1,iimp))
     .                           - ni(ixrb(jx)+1,ny+1,iimp) ) / n0(iimp)
               endif
             enddo
          endif
        endif     # end for if (isimpon .ge. 3 ...
      enddo        # end of impurities along iy = ny+1

ccc  Now do Te, Ti
ccc  - - - - - - - - - - - - - -
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  # loops over te,ti,ng,phi
         ix1 = ixm1(ix,ny+1)
         if(isteonxy(ix,ny+1) .eq. 1) then
           iv1 = idxte(ix,ny+1)
           if (istewcix(ix) .eq. 0) then       # fix flux to zero
             yldot(iv1) = nurlxe*(feey(ix,ny)/(n0(1)*vpnorm*sy(ix,ny)))
     .                                                      / (temp0*ev)
           elseif (istewcix(ix) .eq. 1) then   # set to tewallo
             yldot(iv1) = nurlxe*(tewallo(ix)*ev-te(ix,ny+1))/(temp0*ev)
           elseif (istewcix(ix) .eq. 2) then   # extrapolation
             tbound = te(ix,ny) + gyf(ix,ny-1)*(te(ix,ny)-
     .                           te(ix,ny-1))/gyf(ix,ny)
             tbound = max(tbound, tbmin*ev)
             yldot(iv1) = nurlxe*(tbound - te(ix,ny+1)) / (temp0*ev) 
           elseif (istewcix(ix) .eq. 3) then   #spec. gradient 
             yldot(iv1) = nurlxe*( (te(ix,ny) - te(ix,ny+1)) -
     .                         0.5*(te(ix,ny) + te(ix,ny+1))/
     .                           (gyf(ix,ny)*lyteix(2,ix)) )/(temp0*ev)
           elseif (istewcix(ix) .eq. 4) then   #heat flux ~bcee*te*elec_flux
             yldot(iv1) = nurlxe*( feey(ix,ny) - bceew*ne(ix,ny)*
     .                             vey(ix,ny)*sy(ix,ny)*te(ix,ny+1) )/
     .                                      (vpnorm*ennorm*sy(ix,ny))
            endif
         endif

         if(istionxy(ix,ny+1) .eq. 1) then
           iv2 = idxti(ix,ny+1)
           if (istiwcix(ix) .eq. 0) then       # fix flux to zero
             yldot(iv2) = nurlxi*(feiy(ix,ny)/(n0(1)*vpnorm*sy(ix,ny)))
     .                                                      / (temp0*ev)
           elseif (istiwcix(ix) .eq. 1) then   # set to tiwallo
             yldot(iv2) = nurlxi*(tiwallo(ix)*ev-ti(ix,ny+1))/(temp0*ev)
           elseif (istiwcix(ix) .eq. 2) then   # extrapolation
             tbound = ti(ix,ny) + gyf(ix,ny-1)*(ti(ix,ny)-
     .                           ti(ix,ny-1))/gyf(ix,ny)
             tbound = max(tbound, tbmin*ev)
             yldot(iv2) = nurlxi*(tbound - ti(ix,ny+1)) / (temp0*ev)
           elseif (istiwcix(ix) .eq. 3) then   #spec. gradient 
             yldot(iv2) = nurlxi*( (ti(ix,ny) - ti(ix,ny+1)) -
     .                         0.5*(ti(ix,ny) + ti(ix,ny+1))/
     .                            (gyf(ix,ny)*lytiix(2,ix)) )/(temp0*ev)
           elseif (istiwcix(ix) .eq. 4) then   #heat flux ~bcei*ti*elec_flux
             t0 = max(tg(ix,ny+1,1), temin*ev)
             fngyw=0.25*sqrt(8*t0/(pi*mg(1)))*ng(ix,ny,1)*sy(ix,ny)
             yldot(iv2) = nurlxi*( feiy(ix,ny) - bceiw*ne(ix,ny)*
     .                             vey(ix,ny)*sy(ix,ny)*ti(ix,ny+1) -
     .              cftiexclg*( bcenw*fniy(ix,ny,iigsp)*tg(ix,ny+1,1) +
     .                              cgengw*2.*tg(ix,ny+1,1)*fngyw ) ) /
     .                                       (vpnorm*ennorm*sy(ix,ny))
           endif
         endif
       enddo   # ix-loop for Te,Ti

ccc  Now do the diffusive neutral density (ng) and Tg equations
ccc  - - - - - - - - - - - - - -          
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  # ix-loop for ng & Tg
        nzsp_rt = nhsp
        do igsp = 1, ngsp
          jz = max(igsp - nhgsp, 1)     #identify impurity index
	  if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fniy
          if(isngonxy(ix,ny+1,igsp) .eq. 1) then #skip do-loop if isngon=0
            if (iscplo(ix) .eq. 1) call wsmodo(igsp)
            iv = idxg(ix,ny+1,igsp)
	    t0 = max(cdifg(igsp)*tg(ix,ny+1,igsp), tgmin*ev)
            vyn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
            if (igsp .gt. nhgsp) then
              zflux = 0.
              do iimp = 1, nzsp(jz)
		zflux = zflux + max(fniy(ix,ny,nzsp_rt+iimp), 0.)
              enddo
            endif
c ... prepare chemical sputtering info
            fng_chem = 0.
            do igsp2 = 1, 1+ishymol  #include hydrogen neut sputt
	      t0 = max(tg(ix,ny+1,igsp2), tgmin*ev)
              flx_incid=ng(ix,ny+1,igsp2)*0.25*sqrt(8*t0/(pi*mg(igsp2)))
	      if (isch_sput(igsp).ne.0) then
                call sputchem (isch_sput(igsp), t0/ev, tvwallo(ix),
     .                                      flx_incid, yld_carbo(ix))
                fng_chem = fng_chem + fchemygwo(igsp2)*yld_carbo(ix)*
     .                                            flx_incid*sy(ix,ny)
              else
                fng_chem = fng_chem + chemsputo(igsp,igsp2)*
     .                                            flx_incid*sy(ix,ny) 
              endif
            enddo
c ... add ion sputtering to gas BC
	    sputflxw(ix,igsp) = 0.
            if (isi_sputw(igsp).ge.1 .and. igsp>1) then  # incl phys sput from ions
              do ifld = ipsputt_s, ipsputt_e 
                eng_sput = ( 2*ti(ix,ny+1) + zi(ifld)*
     .                                          3*te(ix,ny+1) )/ev
                if(zi(ifld) > 0) sputflxw(ix,igsp) = sputflxw(ix,igsp) +
     .                            max(fniy(ix,ny,ifld),0.)*fphysyiwo(igsp)*
     .                                    yld96(matp,matt,eng_sput)
              enddo
              if (isi_sputw(igsp) .ge. 2) then  # add chem sput from ions
                do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                  eng_sput = ( 2*ti(ix,ny+1) + zi(ifld)*
     .                                           3*te(ix,ny+1) )/ev
		  flx_incid = max(fniy(ix,ny,ifld)/sy(ix,ny), 0.)
                  call sputchem (isch_sput(igsp), eng_sput, tvwallo(ix),
     .                           flx_incid, yld_chm)
                  sputflxw(ix,igsp) = sputflxw(ix,igsp) + fchemyiwo(igsp)*
     .                                   max(fniy(ix,ny,ifld),0.)*yld_chm
                enddo
              endif
            endif
            nharmave = 2.*(ng(ix,ny,igsp)*ng(ix,ny+1,igsp)) /
     .                    (ng(ix,ny,igsp)+ng(ix,ny+1,igsp))
            fng_alb = (1-albedoo(ix,igsp))*nharmave*vyn*sy(ix,ny)
            yldot(iv) = nurlxg*( fngy(ix,ny,igsp) - fng_alb + fng_chem +
     .                                             sputflxw(ix,igsp) )
     .                                      /(vyn*sy(ix,ny)*n0g(igsp))
            if(matwallo(ix) .gt. 0) then
              if (recycwot(ix,igsp) .gt. 0.) then
ccc
ccc   MER 01 Apr 2002: need to correct fniy below when drifts are included
ccc
                fniy_recy = fac2sp*fniy(ix,ny,1)
                if(isrefluxclip==1) fniy_recy=max(fniy_recy,0.)
                if (igsp .gt. nhgsp) fniy_recy = zflux
                if (ishymol.eq.1 .and. igsp.eq.2) then # 2 atoms per molecule
                  if (isupgon(1) .eq. 1) then
                    fniy_recy = 0.5*( fniy(ix,ny,1) + fniy(ix,ny,2) )
                  else
                    fniy_recy = 0.5*( fniy(ix,ny,1) + fngy(ix,ny,1) )
                  endif
                  if(isrefluxclip==1) fniy_recy=max(fniy_recy,0.)
                endif
                yldot(iv) = nurlxg*( fngy(ix,ny,igsp) + fniy_recy*
     .                        recycwot(ix,igsp) + fngyso(ix,igsp) +
     .                        fngyo_use(ix,igsp) - fng_alb +
     .                        fng_chem + sputflxw(ix,igsp) ) / 
     .                             (vyn*n0g(igsp)*sy(ix,ny))
              elseif (recycwot(ix,igsp) < -1) then
                yldot(iv)=nurlxg*(ngbackg(igsp)-ng(ix,ny+1,igsp))/
     .                                                        n0g(igsp)
             elseif (recycwot(ix,igsp) .le. 0.) then # treat recycw as albedo
                nharmave = 2.*(ng(ix,ny,igsp)*ng(ix,ny+1,igsp)) /
     .                        (ng(ix,ny,igsp)+ng(ix,ny+1,igsp))              
                yldot(iv) = nurlxg*( fngy(ix,ny,igsp) -
     .                  (1+recycwot(ix,igsp))*nharmave*vyn*sy(ix,ny) )/
     .                                    (vyn*n0g(igsp)*sy(ix,ny))
              endif 
            endif
            if(fngyso(ix,igsp)+fngyo_use(ix,igsp).ne.0. 
     .                                          .and. matwallo(ix).eq.0) 
     .                   yldot(iv) = nurlxg*(fngy(ix,ny,igsp) +
     .                         fngyo_use(ix,igsp) + fngyso(ix,igsp) ) 
     .                                      /(vyn*sy(ix,ny)* n0g(igsp))
          endif     # end if-test on isngonxy
 
c... BC for neutral gas temperature/energy at iy=ny+1
          if (istgonxy(ix,ny+1,igsp) == 1) then
            iv = idxtg(ix,ny+1,igsp)
            if (istgwcix(ix,igsp) == 0) then    # fixed Tg
              yldot(iv) = nurlxg*(tgwall(igsp)*ev-tg(ix,ny+1,igsp))/
     .                                                    (temp0*ev)
            elseif (istgwcix(ix,igsp) == 1)    # extrapolation
              tbound = tg(ix,ny,igsp) + gyf(ix,ny)*
     .                     (tg(ix,ny,igsp)-tg(ix,ny-1,igsp))/gyf(ix,ny)
              tbound = max(tbound, 0.25*tbmin*ev)  #tbmin=.1 eV
              yldot(iv) = nurlxi *(tbound - tg(ix,ny+1,igsp))/(temp0*ev)
            elseif (istgwcix(ix,igsp) == 2)    # specified gradient
              yldot(iv) = nurlxi*( (tg(ix,ny,igsp) - tg(ix,ny+1,igsp)) -
     .                         0.5*(tg(ix,ny,igsp) + tg(ix,ny+1,igsp))/
     .                         (gyf(ix,ny)*lytg(2,igsp)) )/(temp0*ev)
            elseif (istgwcix(ix,igsp) == 3)  #Maxwell thermal flux to wall
              t0 = max(cdifg(igsp)*tg(ix,ny,igsp), temin*ev)
              vyn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
              yldot(iv) =  nurlxg*( fegy(ix,ny,igsp) - 2*cgengmw*
     .                              ng(ix,ny,igsp)*vyn*t0*sy(ix,ny) )/
     .                                      (sy(ix,ny)*vpnorm*ennorm)
            elseif (istgwcix(ix,igsp) == 4) 
	      t0 = max(tg(ix,ny+1,igsp),tgmin*ev)
              vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
              nharmave = 2.*(ng(ix,ny,igsp)*ng(ix,ny+1,igsp)) /
     ,                      (ng(ix,ny,igsp)+ng(ix,ny+1,igsp))
              fng_alb = (1-albedoo(ix,1))*nharmave*vyn*sy(ix,ny)
              fng_chem = 0.
              yldot(iv) = nurlxg*(fegy(ix,ny,igsp) - cfalbedo*fng_alb*t0
     .                                             + 2.*fng_chem*t0)
     .                                    /(vpnorm*ennorm*sy(ix,ny))
              fniy_recy = 0.
              if (matwallo(ix) .gt. 0) then
                if (recycwot(ix,igsp) .gt. 0.) then
                  fniy_recy = recycwot(ix,igsp)*fac2sp*fniy(ix,ny,1)
                  if (isrefluxclip==1) fniy_recy=max(fniy_recy,0.)
                  yldot(iv)=nurlxg*(fegy(ix,ny,igsp) - cfalbedo*fng_alb*t0
     .                                             + 2.*fng_chem*t0
     .                                             + fniy_recy*(1.-cfdiss)
     .                                              *cfalbedo*recycwe
     .                                              *ti(ix,ny) )
     .                                     /(vpnorm*ennorm*sy(ix,ny))
                endif
              endif
	          elseif (istgwc(igsp) == 5) then    # set tg=ti*cftgtiwc
              yldot(iv) = nurlxg*(ti(ix,ny+1)*cftgtiwc(igsp) -
     .                                   tg(ix,ny+1,igsp))/(temp0*ev)
            elseif (istgwcix(ix,igsp) > 5)
               call xerrab("***Input error: invalid istgwc ***")
            endif
          endif

        enddo  # igsp loop over gas species
      enddo  # ix-loop for ng and Tg


ccc  Now do the potential
ccc  - - - - - - - - - - - -
      do ix = i4+1-ixmnbcl, i8-1+ixmxbcl  # ix-loop for phi
         if(isphionxy(ix,ny+1) .eq. 1) then
            iv3 = idxphi(ix,ny+1)
	    iv2 = idxphi(ix,ny)
            if(iphibcwoix(ix) == 0) then
              yldot(iv3) = nurlxp*(phi(ix,ny) - phi(ix,ny+1))/temp0
            elseif(iphibcwoix(ix) == 1) then
              yldot(iv3) = nurlxp*
     .               (phintewo*te(ix,ny+1)/ev - phi(ix,ny+1))/temp0
	    elseif(iphibcwoix(ix) == 3) then
              yldot(iv3) = nurlxp*( (phi(ix,ny) - phi(ix,ny+1)) -
     .                         0.5*(phi(ix,ny) + phi(ix,ny+1))/
     .                            (gyf(ix,ny)*lyphiix(2,ix)) )/temp0
	    elseif(iphibcwoix(ix) == 4) then  #fixed prof = phiwo
              yldot(iv3) = nurlxp*(phiwo(ix) - phi(ix,ny+1))/temp0
ccc              yldot(iv2) = nurlxp*(phiwo(ix) - phi(ix,ny))/temp0
            endif
         endif
      enddo  # ix-loop for phi

c...  Now do the corners: the ion density corners are done elsewhere
       if ((xcnearlb.or.openbox) .and. ixmnbcl*iymxbcl.eq.1) then
         do ifld = 1, nusp
           do jx = 1, nxpt
             if(isuponxy(ixlb(jx),ny+1,ifld)==1) then
               yldot(idxu(ixlb(jx),ny+1,ifld)) = - nurlxu*
     .                      ( up(ixlb(jx),ny+1,ifld) -
     .             0.5*(up(ixlb(jx),ny,ifld)+up(ixlb(jx)+1,ny+1,ifld)) )
     .                                         /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
            if(isteonxy(ixlb(jx),ny+1)==1) then
               yldot(idxte(ixlb(jx),ny+1)) = nurlxe*
     .              ( 0.5*(te(ixlb(jx)+1,ny+1)+te(ixlb(jx),ny))
     .                               - te(ixlb(jx),ny+1) ) / (temp0*ev)
            endif
         enddo
         do jx = 1, nxpt
            if(istionxy(ixlb(jx),ny+1)==1) then
               yldot(idxti(ixlb(jx),ny+1)) = nurlxi*
     .              ( 0.5*(ti(ixlb(jx)+1,ny+1)+ti(ixlb(jx),ny))
     .                               - ti(ixlb(jx),ny+1) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
             if(isngonxy(ixlb(jx),ny+1,igsp)==1) then
               yldot(idxg(ixlb(jx),ny+1,igsp)) = nurlxg*
     .            ( ng(ixlb(jx)+1,ny+1,igsp) - ng(ixlb(jx),ny+1,igsp) )
     .                                                      / n0g(igsp)
             endif
             if(istgonxy(ixlb(jx),ny+1,igsp)==1) then 
               yldot(idxtg(ixlb(jx),ny+1,igsp)) = nurlxg*
     .            ( 0.5*(tg(ixlb(jx)+1,ny+1,igsp)+tg(ixlb(jx),ny,igsp))
     .                                       - tg(ixlb(jx),ny+1,igsp) )
     .                                                      / (temp0*ev)
             endif
           enddo
         enddo
       endif

       if ((xcnearrb.or.openbox) .and. ixmxbcl*iymxbcl.eq.1) then
         do ifld = 1, nusp
           do jx = 1, nxpt
             if(isuponxy(ixrb(jx),ny+1,ifld)==1) then
               yldot(idxu(ixrb(jx),ny+1,ifld)) = - nurlxu*
     .                      ( up(ixrb(jx),ny+1,ifld) -
     .             0.5*(up(ixrb(jx)-1,ny+1,ifld)+up(ixrb(jx),ny,ifld)) )
     .                                         /vpnorm
             endif
             if(isuponxy(ixrb(jx)+1,ny+1,ifld)==1) then
               yldot(idxu(ixrb(jx)+1,ny+1,ifld)) = - nurlxu*
     .                   ( up(ixrb(jx)+1,ny+1,ifld)-up(ixrb(jx),ny+1,ifld) )
     .                                            /vpnorm
             endif
           enddo
         enddo
         do jx = 1, nxpt
            if (isteonxy(ixrb(jx)+1,ny+1)==1) then
               yldot(idxte(ixrb(jx)+1,ny+1)) = nurlxe*
     .              ( 0.5*(te(ixrb(jx),ny+1)+te(ixrb(jx)+1,ny))
     .                               - te(ixrb(jx)+1,ny+1) ) / (temp0*ev)
            endif
         enddo
         do jx = 1, nxpt
            if (istionxy(ixrb(jx)+1,ny+1)==1) then
               yldot(idxti(ixrb(jx)+1,ny+1)) = nurlxi*
     .              ( 0.5*(ti(ixrb(jx),ny+1)+ti(ixrb(jx)+1,ny))
     .                               - ti(ixrb(jx)+1,ny+1) ) / (temp0*ev)
            endif
         enddo
         do igsp = 1, ngsp
           do jx = 1, nxpt
             if(isngonxy(ixrb(jx)+1,ny+1,igsp)==1) then
               yldot(idxg(ixrb(jx)+1,ny+1,igsp)) = nurlxg*
     .            ( ng(ixrb(jx),ny+1,igsp) - ng(ixrb(jx)+1,ny+1,igsp) )
     .                                                      / n0g(igsp)
             endif
             if(istgonxy(ixrb(jx)+1,ny+1,igsp)==1) then 
               yldot(idxtg(ixrb(jx)+1,ny+1,igsp)) = nurlxg*
     .            ( 0.5*(tg(ixrb(jx),ny+1,igsp)+tg(ixrb(jx)+1,ny,igsp)) 
     .                               - tg(ixrb(jx)+1,ny+1,igsp) )
     .                                                      / (temp0*ev)
             endif
           enddo
         enddo
       endif  # end of if-test on xcnearrb
      endif   # end of iy = ny+1 boundary conditions except if isnewpot=1

c...  If isnewpot=1, phi boundary involves two eqns at iy=ny and ny+1
CCC   NOTE: special coding for testing alternate phi B.C. (requires 
CCC   isnewpot*isphion=1000, so one can generally ignore this if section) 
      if (isnewpot*isphion.eq.1000 .and. j7.ge.ny-2) then 
        do ix = i4+1-ixmnbcl, i8-1+ixmxbcl
          if(isphionxy(ix,ny+1)*isphionxy(ix,ny)==1) then
            iv  = idxphi(ix,ny+1)
            iv1 = idxphi(ix,ny)
            ix3 = ixm1(ix,ny)
            ix4 = ixm1(ix,ny-1)
            yldot(iv) = nurlxp*( (phi(ix,ny) - phi(ix,ny+1)) -
     .                         0.5*(phi(ix,ny) + phi(ix,ny+1))/
     .                             (gyf(ix,ny)*lyphiix(2,ix)) )/temp0
            yldot(iv1) = -nurlxp*( phi(ix,ny) - phi(ix,ny+1) +
     .                           gyf(ix,ny-1)*(phi(ix,ny)-phi(ix,ny-1))/
     .                                                      gyf(ix,ny) )
          endif    # end of if-test ons isphionxy
        enddo 
      endif

 1200 continue   # jump to here if iymxbcl = 0; interior bdry

c ====================================================================
c ======================== The ix = 0 boundary =======================
c ====================================================================

      if (ixmnbcl .eq. 0) goto 1300   #skip setting Eqn because interior bdry

c********************************************************************
c...  First check if ix=0 has fixed boundary values, no potential
c...  isfixlb=1 sets all profiles; isfixlb=2 sets reflection boundary 
c...  conditions
c********************************************************************
      if (i3 .le. 0 .and. isfixlb(1) .ne. 0) then
      do 176 iy = j2, j5
         do 174 ifld = 1 , nisp
            if(isnionxy(0,iy,ifld) .eq. 1) then
               iv1 = idxn(0,iy,ifld)
               yldot(iv1) = nurlxn *
     .                (nib(ifld)*nibprof(iy)-ni(0,iy,ifld))/n0(ifld)
               if(isfixlb(1).eq.2) yldot(iv1) = nurlxn * (1/n0(ifld)) *
     .                               (ni(1,iy,ifld) - ni(0,iy,ifld))
            endif
 174     continue
         do ifld = 1, nusp
            if(isuponxy(0,iy,ifld) .eq. 1) then
               iv2 = idxu(0,iy,ifld)
               yldot(iv2) = nurlxu *
     .           (upb(ifld)*upbprof(iy) - up(0,iy,ifld))/vpnorm
               if(isfixlb(1).eq.2) yldot(iv2) = nurlxu *
     .                           (0. - up(0,iy,ifld))/vpnorm
               if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
                  cs = sqrt( (te(0,iy)+ti(0,iy))/mi(ifld) )
                  yldot(iv2) = nurlxu*
     .                          (-cs -up(0,iy,ifld))/vpnorm

               endif
            endif
         enddo

c...  now do the gas and temperatures
         if(isteonxy(0,iy) .eq. 1) then
           iv1 = idxte(0,iy)
           yldot(iv1) = nurlxe * ne(0,iy) *
     .                     (teb*ev*tebprof(iy) - te(0,iy))/ennorm
           if(isfixlb(1).eq.2) yldot(iv1) = nurlxe * ne(0,iy) *
     .                               (te(1,iy) - te(0,iy))/ennorm
           if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
              yldot(iv1) = - nurlxe*(feex(0,iy)/sx(0,iy) - bcee*
     .                               ne(0,iy)*vex(0,iy)*te(0,iy))/
     .                                 (vpnorm*ennorm)
           endif
         endif
         if(istionxy(0,iy) .eq. 1) then
           iv2 = idxti(0,iy)
            yldot(iv2) = nurlxi * ne(0,iy) *
     .                     (tib*ev*tibprof(iy) - ti(0,iy))/ennorm
            if(isfixlb(1).eq.2) yldot(iv2) = nurlxi * ne(0,iy) *
     .                               (ti(1,iy) - ti(0,iy))/ennorm
            if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
               yldot(iv2) = -nurlxi*
     .          ( feix(0,iy) - bcei*ti(0,iy)*fac2sp*fnix(0,iy,1) ) / 
     .                                       (vpnorm*ennorm*sx(0,iy))
            endif
         endif
         do igsp = 1, ngsp
            if(isngonxy(0,iy,igsp) .eq. 1) then
               iv = idxg(0,iy,igsp)
               yldot(iv) = nurlxg * (ngbackg(igsp) - 
     .                                         ng(0,iy,igsp))/n0g(igsp)
               if(isfixlb(1).eq.2) yldot(iv) = nurlxg * 
     .                        (ng(1,iy,igsp) - ng(0,iy,igsp))/n0g(igsp)
               if(isfixlb(1).eq.2 .and. yylb(iy,1).gt.rlimiter) then
                  t1 = engbsr * max(tg(1,iy,igsp),tgmin*ev)
                  vxn = 0.25 * sqrt( 8*t1/(pi*mg(igsp)) )
                  flux_inc = fac2sp*fnix(0,iy,1)
                  if (ishymol.eq.1 .and. igsp.eq.2) then
                    if (isupgon(1) .eq. 1) then  # two atoms per molecule
                      flux_inc = 0.5*( fnix(0,iy,1) + fnix(0,iy,2) )
                    else
                      flux_inc = 0.5*( fnix(0,iy,1) + fngx(0,iy,1) ) 
                    endif
                  endif
                  areapl = isoldalbarea*sx(0,iy) + (1-isoldalbarea)*sxnp(0,iy)
                  yldot(iv) = -nurlxg * ( fngx(0,iy,igsp) -
     .                                           fngxlb_use(iy,igsp,1) +
     .                 fngxslb(iy,igsp,1) + recylb(iy,igsp,1)*flux_inc +
     .                  (1-alblb(iy,igsp,1))*ng(1,iy,igsp)*vxn*areapl )
     .                                     / (vpnorm*n0g(igsp)*sx(0,iy))
               endif
               if (is1D_gbx.eq.1) yldot(iv) = nurlxg*(ng(1,iy,igsp) -
     .                                    ng(0,iy,igsp))/n0g(igsp)
            endif
         enddo
c ... Neutral temperature - test if tg eqn is on, then set BC
	 do igsp = 1, ngsp
           if (istgonxy(0,iy,igsp) == 1) then
             iv = idxtg(0,iy,igsp)
             yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(0,iy,igsp))/(temp0*ev)
             if(isfixlb(1)==2) then #just above applies if isfixlb=1
               yldot(iv)=nurlxg*(tg(1,iy,igsp)-tg(0,iy,igsp))/(temp0*ev)
             endif
             if(isfixlb(1)==2 .and. yylb(iy,1) > rlimiter) then
               yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(0,iy,igsp))/(temp0*ev)
             endif
           endif
         enddo

         if (isphionxy(0,iy) .eq. 1) then
            iv = idxphi(0,iy)
            yldot(iv) = nurlxp*(phi(1,iy) - phi(0,iy))/temp0
         endif

 176  continue

      endif         # end of ix = 0, isfixlb.ne.0 boundary conditions

c...  If isfixlb=2, check if i2,i5 range for yldot loop in pandf includes ixpt2(1)
c...  Then overwrite pandf value is up(ixpt2(1)) --> 0 eqn.
      if (isfixlb(1) .eq. 2) then
         if (i2.le.ixpt2(1) .and. i5.ge.ixpt2(1) .and. j2.le.iysptrx2(1)) then  
            do 179 ifld = 1, nusp
               do 178 iy = 0+1-iymnbcl, iysptrx2(1)
                 if(isuponxy(ixpt2(1),iy,ifld)==1) then
                     iv = idxu(ixpt2(1),iy,ifld)
                     yldot(iv) = nurlxu*(0.-up(ixpt2(1),iy,ifld))/vpnorm
                 endif
 178           continue
 179        continue
         endif 
      endif            # end of isfixlb=2, check for ix=ixpt2
c***********************************************************************

c************************************************************************
c************************************************************************
c     begin standard divertor plate conditions for left boundaries 
c************************************************************************
c************************************************************************
      do jx = 1, nxpt # loop over nxpt mesh regions
      if (xcnearlb .or. openbox) then
c     First, the density equations --
        do ifld = 1, nisp
          if ( i3.le.(ixlb(jx)+isextrnp) .and. isfixlb(jx)==0) then
          do iy = j2, j5
       
            ixt  = ixlb(jx)        # analog of ix=0
            ixt1 = ixp1(ixt,iy)    # analog of ix=1
            ixt2 = ixp1(ixt1,iy)   # analog of ix=2

            if(isnionxy(ixt,iy,ifld)==1) then
              iv1 = idxn(ixt,iy,ifld)
              if (isupgon(1)==1 .and. zi(ifld)==0.0) then   ## neutrals
                if (recylb(iy,1,jx) .gt. 0.) then           # recycling
                  t0 = max(tg(ixt1,iy,1),tgmin*ev) 
                  vxn = 0.25 * sqrt( 8*t0/(pi*mi(ifld)) )
                  areapl = isoldalbarea*sx(ixt,iy) + (1-isoldalbarea)*sxnp(ixt,iy)
                  yldot(iv1) = -nurlxg *
     .             (fnix(ixt,iy,ifld) + recylb(iy,1,jx)*fnix(ixt,iy,1) -
     .                                               fngxlb_use(iy,1,jx) +
     .              (1-alblb(iy,1,jx))*ni(ixt1,iy,ifld)*vxn*areapl -
     .                 fngxslb(iy,1,jx) ) / (vpnorm*n0(ifld)*sx(ixt,iy))
                elseif (recylb(iy,1,jx) <=  0. .and. 
     .                  recylb(iy,1,jx) >= -1.) then  # recylb is albedo
                  t0 = max(tg(ixt,iy,1),tgmin*ev) 
                  vyn = sqrt( 0.5*t0/(pi*mi(1)) )
                  yldot(iv1) = -nurlxg * ( fnix(ixt,iy,ifld) +
     .             (1+recylb(iy,1,jx))*ni(ixt,iy,ifld)*vyn*sx(ixt,iy) )/
     .               (vpnorm*n0(ifld)*sx(ixt,iy))
                elseif (recylb(iy,1,jx) < -1. .and.
     .                 recylb(iy,1,jx) > -2. ) then  #fix density nglfix
                  yldot(iv1)=nurlxg*(nglfix - ni(ixt,iy,ifld))/n0(ifld)
                elseif (recylb(iy,1,jx) <= -2.) then  #zero gradient
                  yldot(iv1) =nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                endif # end if-test on recylb
              else                                     ## ions
                if (isextrnp==0) then                  # zero x-gradient
                   yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                else                                   # extrapolation
                   nbound =  ni(ixt1,iy,ifld) - gxf(ixt1,iy)*
     .                   (ni(ixt2,iy,ifld)-ni(ixt1,iy,ifld))/gxf(ixt,iy)
                   nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                    ni(ixt,iy,ifld)-1) ) ) + 0.2*ni(ixt,iy,ifld)
                   yldot(iv1) = nurlxn*(nbound-ni(ixt,iy,ifld))/n0(ifld)
                endif # end if-test on isextrnp
              endif # end if-test on isupgon and zi
            endif   # end if-test on isnionxy
          enddo     # end do-loop on iy
          endif     # end if-test on i3 and isfixlb and isnion
        enddo       # end do-loop on ifld

c     Now do the parallel velocity and other variables --
      if ( (i3 .le. ixlb(jx)) .and. isfixlb(jx)==0 ) then
      do iy = j2, j5 # begin big do-loop on iy

        ixt  = ixlb(jx)        # analog of ix=0
        ixt1 = ixp1(ixt,iy)    # analog of ix=1
        ixt2 = ixp1(ixt1,iy)   # analog of ix=2

        kfeix = 0.
        sumb = 0.
        do ifld = 1, nfsp  # set up generalized Bohm condition
          ueb = cfueb*( cf2ef*v2ce(ixt,iy,ifld)*rbfbt(ixt,iy) -
     .            vytan(ixt,iy,ifld) ) / rrv(ixt,iy) 
          sumb = sumb + ni(ixt,iy,ifld)*zi(ifld)**2*te(ixt,iy) /
     .                (mi(ifld)*(upi(ixt,iy,ifld)+ueb)**2 - ti(ixt,iy))
        enddo # end do-loop on ifld for Bohm condition
        sumb = sqrt(abs(sumb/ne(ixt,iy)))

c       Next, the momentum equations --
        do ifld = 1, nusp
          if (isuponxy(ixt,iy,ifld)==1) then
            iv2 = idxu(ixt,iy,ifld)
            cs = csfaclb(ifld,jx)*sqrt( (te(ixt,iy)+
     .                                  csfacti*ti(ixt,iy))/mi(ifld) )
            if (isupgon(1)==1 .and. zi(ifld)==0.0) then  ## neutrals
              if (recycmlb(iy,1,jx) > -9.9) then  # backscatter with recycm
                yldot(iv2) = -nurlxu*(recycmlb(iy,1,jx)*up(ixt,iy,1) + 
     .                                       up(ixt,iy,ifld))/vpnorm
              elseif (recycmlb(iy,1,jx) <= -9.9 .and. 
     .                    recycmlb(iy,1,jx) > -10.1) then # zero x-gradient
                yldot(iv2) = nurlxu*(up(ixt1,iy,ifld) -
     .                                          up(ixt,iy,ifld))/vpnorm
	      else  # neutral thermal flux to wall if recycm < -10.1
                t0 = max(tg(ixt1,iy,1),tgmin*ev)
                vxn = cgmompl*0.25*sqrt( 8*t0/(pi*mi(ifld)) ) 
		vparn = up(ixt,iy,ifld)
                yldot(iv2) = -nurlxu*( fmix(ixt1,iy,ifld) + vparn*vxn*
     .                        0.5*(nm(ixt1,iy,ifld)+nm(ixt,iy,ifld))*
     .                                                 sx(ixt,iy) ) /
     .                               (vpnorm*fnorm(ifld)*sx(ixt,iy))
              endif
            else                                         ## ions
              ueb = cfueb*( cf2ef*v2ce(ixt,iy,ifld)*rbfbt(ixt,iy) -
     .                vytan(ixt,iy,ifld) ) / rrv(ixt,iy) 
              yldot(iv2) = -nurlxu*(sumb - 1.)  # multispecies Bohm
              if (isbohmms==0) then          # simple Bohm condition
                yldot(iv2) = nurlxu * (-cs-ueb-up(ixt,iy,ifld))/vpnorm
              endif
              if(isupss(ifld)==1 .and. up(ixt1,iy,ifld)+ueb .lt. -cs)
                                             # dup/dx=0 if supersonic
     .          yldot(iv2) = nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/
     .                                                           vpnorm
              if (isupss(ifld)==-1) then     # slip boundary conditions
                yldot(iv2) = nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/
     .                                                           vpnorm
              elseif (isupss(ifld)==-2) then # extrap. + no pos. uu
                vbound = up(ixt1,iy,ifld) - gx(ixt2,iy)*
     .                   (up(ixt2,iy,ifld)-up(ixt1,iy,ifld))/gx(ixt1,iy) 
                vbound = min(vbound, -ueb)   # forces uu & fnix =< 0
                yldot(iv2) = nurlxu*(vbound - up(ixt,iy,ifld))/vpnorm
              elseif (isupss(ifld)==-3) then # modified Bohm condition
                vbound = -ueb -2*cs*uu(ixt1,iy,ifld)/(uu(ixt1,iy,ifld)-cs)
                vbound = min(vbound, -ueb)   # forces uu & fnix =< 0
                yldot(iv2) = nurlxu * (vbound-up(ixt,iy,ifld))/vpnorm
              endif # end if-test on isupss
            endif # end if-test on isupgon
          endif # end if-test on isuponxy
          if (zi(ifld)==0.0) then 
            kfeix = kfeix - cftiexclg*cfvcsx(ifld)*0.5*sx(ixt,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt,iy,ifld)**2 )
          else
            kfeix = kfeix - cfvcsx(ifld)*0.5*sx(ixt,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt,iy,ifld)**2 ) 
          endif
        enddo # end do-loop on ifld

c       Next, the potential equation --
        if(isphion+isphiofft .eq. 1) then
cc           fqpsatlb(iy,jx) = - qe*isfdiax*( ne(ixt,iy)*v2ce(ixt,iy,1)*
cc     .                              rbfbt(ixt,iy)*sx(ixt,iy) + fdiaxlb(iy,jx) )
cc           do ifld = 1, nfsp   # note fqp,fqpsat are poloidal proj. of || curr
cc              fqpsatlb(iy,jx) = fqpsatlb(iy,jx) - qe*zi(ifld)*ni(ixt,iy,ifld)*
cc     .                                   upi(ixt,iy,ifld)*sx(ixt,iy)*rrv(ixt,iy)
cc           enddo
           lambdae = 2e16*(te(ixt,iy)/ev)**2/ne(ixt,iy)  #approx mfp
           kincorlb(iy,jx) = 1./(1 + cfkincor*(lambdae/lcone(ixt,iy))*
     .                                       abs(ev*phi(ixt,iy)/te(ixt,iy)))
           fqpsate = qe*ne(ixt,iy)*sqrt(te(ixt,iy)/(2*pi*me))*
     .                                    kincorlb(iy,jx)*sx(ixt,iy)*rrv(ixt,iy)
c          NOTE: by definition, fqpsate is always > 0
           if (ikapmod==0) then
cc              if (fqp(ixt,iy) < 0.) then #limit to saturation current
cc                fqp(ixt,iy)=-( abs((fqp(ixt,iy)*fqpsatlb(iy,jx)))**exjbdry/
cc     .                      (abs(fqp(ixt,iy))**exjbdry +
cc     .                       abs(fqpsatlb(iy,jx))**exjbdry) )**(1/exjbdry)
cc              endif
              if ( fqpsatlb(iy,jx)+(1.-gamsec)*fqp(ixt,iy) > 0) then  # force +ve log argument
                 arglgphi=( ((fqpsatlb(iy,jx)+(1.-gamsec)*fqp(ixt,iy))
     .                          /fqpsate)**2 + expkmx**2 )**0.5
              else
                 arglgphi = expkmx
              endif
              if (iskaplex .eq. 0) kappal(iy,jx) = - log(arglgphi)
              if (newbcl(jx).eq.0 .and. iskaplex.eq.0) kappal(iy,jx) = 3.0
           elseif (ikapmod==1) then
              kappal(iy,jx) = kappa(fqpsatlb(iy,jx), fqpsate, 
     .                                         -(1.-gamsec)*fqp(ixt,iy))
           endif

           if (isphionxy(ixt,iy) .eq. 1) then
              iv = idxphi(ixt,iy)
              if ((isnewpot==100) .and. ((iy==1) .or. (iy==ny))) then
                 continue
              elseif (isphilbc==1) then
                 yldot(iv) = - nurlxp * (phi(ixt,iy) - phi0l(iy,jx))/temp0 
              else
                 yldot(iv) = -nurlxp*(1.-bctype(iy)) *
     .            (phi(ixt,iy)-kappal(iy,jx)*te(ixt,iy)/ev-phi0l(iy,jx))/temp0
     .                    -nurlxp * bctype(iy) * (1.-gamsec)*fqp(ixt,iy) 
     .                                         / (fqpsatlb(iy,jx)+cutlo)
              endif
           endif
	else    # corresponding to (isphion+isphiofft .ne. 1)
           fqpsate = 0.
	   kappal(iy,jx) = 3.  # some default of e*phi/Te - used for sputtering
        endif
        isphion2 = isphion + isphiofft
        bcel(iy,jx) = (1-newbcl(jx)*isphion2) * bcee
     .              + newbcl(jx)*isphion2 * (2. + kappal(iy,jx))
        if (iskaplex .eq. 1) bcel(iy,jx) = (2. + kappal(iy,jx))
        bcil(iy,jx) = (1-newbcl(jx)*isphion2) * bcei
     .              + newbcl(jx)*isphion2 * (2.5)

        t0 = te(ixt,iy)/ev
        f_cgpld = .5*(1. - cos(pi*(t0 - temin)/(.3 - temin)))
                                # ramp cgpld effect down at low guard cell
                                # Te; rampdown might not be necessary if
                                # guard cell were allowed to access negative
                                # temperatures (ETM 7 Jan 2014)
        if(t0 < temin) f_cgpld = 0.
        if(t0 > 0.3) f_cgpld = 1.
        t0 = max(tg(ixt1,iy,1),tgmin*ev)
        vxn = f_cgpld * 0.25 * sqrt(8*t0/(pi*mg(1)))

c   Do the electron temp Eqn -----------------------------------
        if (isteonxy(ixt,iy) == 1) then
	  iv1 = idxte(ixt,iy)
          if (ibctepl == 1) then
            faceel =  bcel(iy,jx)*(fqpsate/qe)*exp(-kappal(iy,jx))
            faceel2 = bcel(iy,jx)*(fqpsate/qe)*exp(-kappamx+2)
            totfeexl(iy,jx) = feex(ixt,iy) + cfeexdbo*( 
     .             2.5*fdiaxlb(iy,jx) + floxebgt(ixt,iy) )*te(ixt,iy)
	    totfnex = ne(ixt,iy)*vex(ixt,iy)*sx(ixt,iy)
            if (isphion+isphiofft==1) then
              yldot(iv1) =-nurlxe*(totfeexl(iy,jx)
     .                        +faceel*te(ixt,iy)
     .                        +faceel2*(te(ixt,iy)-te(ixt1,iy))
     .                        -cmneut*fnix(ixt,iy,1)*recycp(1)*eedisspl*ev
     .                                     )/(sx(ixt,iy)*vpnorm*ennorm)
            else
              yldot(iv1) =-nurlxe*(totfeexl(iy,jx) 
     .                        -totfnex*te(ixt,iy)*bcel(iy,jx)
     .                        +cgpld*sx(ixt,iy)*0.5*ng(ixt1,iy,1)*vxn*ediss*ev
     .                        -cmneut*fnix(ixt,iy,1)*recycp(1)*eedisspl*ev
     .                                     )/(sx(ixt,iy)*vpnorm*ennorm)
            endif
          elseif (ibctepl .eq. 0) then
             yldot(iv1) = nurlxe*(tepltl*ev-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctepl .eq. 2) then
             yldot(iv1) = nurlxe*(te(ixt1,iy)-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          endif # end loop for ibctepl
        endif  # end loop for isteon=1

c  Do the ion temp Eqn ---------------------------
        if (istionxy(ixt,iy) == 1) then
          iv2 = idxti(ixt,iy)
          if (ibctipl == 1) then
c            totfeix is thermal + kinetic + viscous flux
            totfeixl(iy,jx) = feix(ixt,iy) + ckinfl*kfeix 
            totfnix = 0.
            do ifld = 1, nfsp
            # the condition on zi must include use of v2cd which is not properly set to zero in pandf1
            if(zi(ifld)>1e-10) then
              totfeixl(iy,jx) = totfeixl(iy,jx) + cfeixdbo*(
     .                      2.5* ni(ixt,iy,ifld)*v2cd(ixt,iy,ifld)*
     .                                sx(ixt,iy)*rbfbt(ixt,iy)+
     .                             floxibgt(ixt,iy,ifld) )*ti(ixt,iy)
              totfnix = totfnix + fnix(ixt,iy,ifld)
              endif
            enddo
            if (isupgon(1)==1) then
c              Different boundary conditions for neutral momentum equation
cc               if (recylb(iy,1,jx) .gt. 0.) then
cc                  bcen = recyce*bcil(iy,jx) - 0.5*mi(1)*(up(ixt,iy,2)**2 -
cc     .               recyce*upi(ixt,iy,1)**2)/(recylb(iy,1,jx)*ti(ixt,iy))
cc               endif
cc               if (recyce .le. 0) bcen = 0.  # gets back to old case

                 yldot(iv2) = -nurlxi*(totfeixl(iy,jx)
     .                        -totfnix*ti(ixt,iy)*bcil(iy,jx)+ 
     .            cftiexclg*( -cfneut*fnix(ixt,iy,iigsp)*tg(ixt,iy,1)*bcen
     .                        +(cgengpl*2.*tg(ixt,iy,1) - cgpld*eion*ev)*
     .                                ng(ixt1,iy,1)*vxn*sx(ixt,iy)
     .                        -cmneut*fnix(ixt,iy,1)*recycp(1)*
     .                                cmntgpl*(ti(ixt,iy)-eidisspl*ev)
     .                                  ) )/(vpnorm*ennorm*sx(ixt,iy))
            else
                 yldot(iv2) = -nurlxi*(totfeixl(iy,jx) 
     .                         -totfnix*bcil(iy,jx)*ti(ixt,iy)+
     .             cftiexclg*( -cmneut*fnix(ixt,iy,1)*recycp(1)*
     .                                cmntgpl*(ti(ixt,iy)-eidisspl*ev)
     .                                  ) )/(vpnorm*ennorm*sx(ixt,iy))
            endif #end loop for isupgon(1)

          elseif (ibctipl .eq. 0) then
             yldot(iv2) = nurlxi*(tipltl*ev - ti(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctipl .eq. 2) then
             yldot(iv2) = nurlxi*(ti(ixt1,iy) - ti(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          endif # end loop for ibctipl
        endif  # end loop for istion=1

c       Do hydrogenic gas equations --
        do igsp = 1, nhgsp  # imp gas below
           if (isngonxy(ixt,iy,igsp) .eq. 1) then
             iv = idxg(ixt,iy,igsp)
             if (recylb(iy,igsp,jx) .gt. 0.) then  # normal recycling
               flux_inc = fac2sp*fnix(ixt,iy,1)
               if (ishymol.eq.1 .and. igsp.eq.2) then
                 if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                   flux_inc = 0.5*( fnix(ixt,iy,1) + fnix(ixt,iy,2) ) 
                 else
                   flux_inc = 0.5*( fnix(ixt,iy,1) + fngx(ixt,iy,1) ) 
                 endif
               endif
               t0 = max(tg(ixt1,iy,igsp), tgmin*ev)
               vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
               areapl = isoldalbarea*sx(ixt,iy) + (1-isoldalbarea)*sxnp(ixt,iy)
               yldot(iv) = -nurlxg * ( fngx(ixt,iy,igsp) - 
     .                                           fngxlb_use(iy,igsp,jx) -
     .               fngxslb(iy,igsp,jx) + recylb(iy,igsp,jx)*flux_inc +
     .               (1-alblb(iy,igsp,jx))*ng(ixt1,iy,igsp)*vxn*areapl )
     .                                   / (vpnorm*n0g(igsp)*sx(ixt,iy))
             elseif (recylb(iy,igsp,jx) <=  0. .and.
     .               recylb(iy,igsp,jx) >= -1.) then # recylb is albedo
               t0 = max(tg(ixt,iy,igsp), tgmin*ev)
               vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
               yldot(iv) = -nurlxg*( fngx(ixt,iy,igsp) +
     .           (1+recylb(iy,igsp,jx))*ng(ixt,iy,igsp)*vxn*sx(ixt,iy) )
     .                                      / (vxn*sx(ixt,iy)*n0g(igsp))
             elseif (recylb(iy,igsp,jx) < -1. .and.
     .               recylb(iy,igsp,jx) > -2.) then  #fix density nglfix
               yldot(iv)=nurlxg*(nglfix - ng(ixt,iy,igsp))/n0g(igsp)
             elseif (recylb(iy,igsp,jx) <= -2.) then #zero gradient
               yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                         n0g(ifld)
             endif
             if (is1D_gbx.eq.1) yldot(iv) = nurlxg*(ng(ixt1,iy,igsp) -
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
c   Special coding for Maxim
             if (isnglf==1) yldot(iv) = nurlxg*(nglfix - 
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
           endif # end if-test on isngon
        enddo # end do-loop on igsp

c       Finally, do impurity gas equations --
c       Set neutral gas impurity flux to be the sum of the fluxes of the
c       sputtered impurities plus recycled impurities from all charge states.
        if (isimpon .ge. 4 .and. isimpon .le. 7 .and. nzspt.ge.1) then
           nzsp_rt = nhsp
           do igsp = nhgsp+1, ngsp
	      jz = max(igsp - nhgsp, 1)   # identify impurity index
	      if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fnix
              if (isngonxy(ixt,iy,igsp) .eq. 1) then
                 iv = idxg(ixt,iy,igsp)
                 hflux = 0.
                 do ihyd = 1, nhsp
                    if (zi(ihyd).gt.0.) hflux = hflux + fnix(ixt,iy,ihyd)
                 enddo
                 zflux = 0.
                 do iimp = 1, nzsp(jz) # loop limited by numb. species
                    zflux = zflux + fnix(ixt,iy,nzsp_rt+iimp)
                 enddo
                 sputflxlb(iy,igsp,jx) = 0.
                 zflux_chm = 0.
                 if (isph_sput(igsp) .ge. 1) then  # use fits for phys sput
                   do ifld = ipsputt_s, ipsputt_e 
                     eng_sput = ( 0.5*mi(ifld)*up(ixt,iy,ifld)**2 + 
     .                            ti(ixt,iy) +  zi(ifld)*
     .                            kappal(iy,jx)*te(ixt,iy) )/ev
                     if(zi(ifld)>0.) sputflxlb(iy,igsp,jx) = 
     .                                          sputflxlb(iy,igsp,jx) +
     .                                           fnix(ixt,iy,ifld)*
     .                                           fphysylb(igsp,jx)*
     .                                        yld96(matp,matt,eng_sput)
                   enddo
                   if (isph_sput(igsp) .ge. 2) then  # add chem sput from ions
                     do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                       eng_sput = ( 0.5*mi(ifld)*up(ixt,iy,ifld)**2 + 
     .                            ti(ixt,iy) +  zi(ifld)*
     .                            kappal(iy,jx)*te(ixt,iy) )/ev
                       flx_incid = abs(fnix(ixt,iy,ifld))/sx(ixt,iy)
                       call sputchem (isch_sput(igsp),eng_sput,tvplatlb(iy,jx),
     .                             flx_incid, yld_chm)
                       sputflxlb(iy,igsp,jx) = sputflxlb(iy,igsp,jx) +
     .                       fchemylb(igsp,jx)*fnix(ixt,iy,ifld)*yld_chm
                     enddo
                   endif
                   if (isph_sput(igsp) .eq. 3) then #add chem sput from h neuts
                     do igsp2 = 1, 1+ishymol  #hydrogen neut fluxes only
                       t0p = max(tg(ixt1,iy,igsp2),tgmin*ev)
                       flx_incid = ng(ixt,iy,igsp2)*.25*
     .                                        sqrt(8*t0p/(pi*mg(igsp2)))
                       call sputchem (isch_sput(igsp),t0p/ev,tvplatlb(iy,jx),
     .                                               flx_incid, yld_chm)
                       zflux_chm = zflux_chm - flx_incid*
     .                              fchemylb(igsp,jx)*yld_chm*sx(ixt,iy)
                     enddo
                   endif
                 endif
                 if (sputtlb(iy,igsp,jx) .ge. 0. .or. 
     .                              abs(sputflxlb(iy,igsp,jx)).gt. 0.) then
                    t0 = max(cdifg(igsp)*tg(ixt1,iy,igsp), tgmin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    areapl = isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy)
                    zflux = - sputtlb(iy,igsp,jx) * hflux - 
     .                        sputflxlb(iy,igsp,jx) -
     .                   recylb(iy,igsp,jx) * zflux -
     .                (1-alblb(iy,igsp,jx))*ng(ixt1,iy,igsp)*vxn*areapl-
     .                   zflux_chm + fngxslb(iy,igsp,jx)+fngxlb_use(iy,igsp,jx)
                    yldot(iv) = -nurlxg * (fngx(ixt,iy,igsp) - zflux) /
     .                         (n0(igsp) * vpnorm * sx(ixt,iy))
                 elseif (sputtlb(iy,igsp,jx).ge.-9.9) then # neg. sputtlb ==> albedo
                    t0 = max(cdifg(igsp)*tg(ixt1,iy,igsp), tgmin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    yldot(iv) = -nurlxg*( fngx(ixt,iy,igsp) -
     .                (1+sputtlb(iy,igsp,jx))*ng(ixt,iy,igsp)*vxn*sx(ixt,iy) )
     .                                       / (vxn*sx(ixt,iy)*n0g(igsp))
                 else                                # sputtlb < -9.9 ==> fix dens
                    yldot(iv) = -nurlxg*(ng(ixt,iy,igsp)-ngplatlb(igsp,jx))/
     .                                                        n0g(igsp)
                 endif  # end if-test on sputtlb
              endif  # end if-test on isngon
           enddo  # end do-loop in igsp
        endif  # end if-test on isimpon and nzspt

c ... Neutral temperature - test if tg eqn is on, then set BC
	do igsp = 1, ngsp
          if (istgonxy(ixt,iy,igsp) == 1) then
            iv = idxtg(ixt,iy,igsp)
            if (istglb(igsp) == 0) then  #set tg=tgwall
              yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istglb(igsp) == 1)  #extrap. tg from interior cells
              tbound = tg(ixt1,iy,igsp) - gyf(ixt1,iy)*
     .                   (tg(ixt2,iy,igsp)-tg(ixt1,iy,igsp))/gxf(ixt,iy)
              tbound = max(tbound,0.5*temin*ev)
              yldot(iv) = nurlxg*(tbound - tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istglb(igsp) == 2)  #placeholder for gradient BC
              call xerrab("**INPUT ERROR: istglb=2 grad opt not implemented")
            elseif (istglb(igsp) == 3)  #Maxwell thermal flux to wall
              t0 = max(cdifg(igsp)*tg(ixt1,iy,igsp), tgmin*ev)
              vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
              yldot(iv) =  -nurlxg*( fegx(ixt,iy,igsp) + 2*cgengmpl*
     .                            ng(ixt1,iy,igsp)*vxn*t0*sx(ixt,iy) )/
     .                                     (sx(ixt,iy)*vpnorm*ennorm)
            elseif (istglb(igsp) == 4) 
	      if (isupgon(igsp)==1) then
                if (recylb(iy,igsp,jx) .gt. 0.) then
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                  fng_alb=(1-alblb(iy,igsp,jx))*ng(ixt1,iy,igsp)
     .                                         *vxn*sx(ixt,iy)
                  yldot(iv) = -nurlxg*( fegx(ixt,iy,igsp)
     .                                 +cfalbedo*fng_alb*t0
     .                                 +recylb(iy,igsp,jx)*(1.-cfdiss)
     .                                 *fnix(ixt,iy,1)
     .                                 *recyce*cfalbedo
     .                                 *( kappal(iy,jx)*zi(1)*te(ixt,iy)
     .                                   +ti(ixt,iy) ) )/
     .                                 (vpnorm*ennorm*sx(ixt,iy))
                elseif (recylb(iy,igsp,jx) <=  0. .and.
     .                  recylb(iy,igsp,jx) >= -1.) then  # recylb is albedo
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
                  fng_alb = (1+recylb(iy,igsp,jx))*ng(ixt1,iy,igsp)
     .                                            *vyn*sx(ixt,iy)
                  yldot(iv) = -nurlxg*( fegx(ixt,iy,igsp)
     .                                 +cfalbedo*fng_alb*t0 )/
     .                                 (vpnorm*ennorm*sx(ixt,iy))
                elseif (recylb(iy,igsp,jx) < -1.) then  #..half Maxwellian
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
                  yldot(iv) = -nurlxg*( fegx(ixt,iy,igsp)
     .                                 +cfalbedo*fnix(ixt,iy,iigsp)*t0 )/
     .                                 (vpnorm*ennorm*sx(ixt,iy))
                endif
              endif
            elseif (istglb(igsp) == 5) then  #set tg=ti*cftgtipltl
              yldot(iv) = nurlxg*(ti(ixt,iy)*cftgtipltl(igsp) -
     .                                    tg(ixt,iy,igsp))/(temp0*ev)
            else
              call xerrab("**INPUT ERROR: istglb set to unknown option")
            endif
          endif
        enddo

      enddo # end big do-loop on iy

c ...   Make electric field at plate uniform over last 2 cells --
        do iy = j2, j5
           if (isphionxy(ixt1,iy)==1 .and. isexunif==1) then
              iv = idxphi(ixt1,iy)  # next formula assumes uniform mesh
              yldot(iv)=-nurlxp*(3*phi(ixt1,iy)-2*phi(ixt,iy)-
     .                                            phi(ixt2,iy))/temp0
           endif
        enddo

      endif # end if-test on i3 and isfixlb

      endif # end if-test on xcnearlb and openbox
      enddo # end do-loop over nxpt mesh regions
c************************************************************************
c     end standard divertor plate conditions for left boundaries
c************************************************************************

c ... Special case for no divertor leg (ixpt1(1)=0, iy.le.iysptrx1(1)), just
c ... use continuation boundary conditions for the guard cell (ix=0)
ccc      if (ixpt1(1).eq.0 .and. isfixlb(1).eq.0) then  # should test i3.le.0 (no extrap)
ccc
ccc      do iy = j2, min(j5, iysptrx1(1))
ccc        do ifld = 1, nisp
ccc          if (isnionxy(0,iy,ifld) .eq. 1) then
ccc            iv = idxn(0,iy,ifld)
ccc            yldot(iv) = -nurlxn*(ni(0,iy,ifld)-ni(1,iy,ifld))/n0(ifld)
ccc          endif
ccc        enddo
ccc        do ifld = 1, nusp
ccc          if (isuponxy(0,iy,ifld) .eq. 1) then
ccc            iv = idxu(0,iy,ifld)
ccc            yldot(iv) = -nurlxu*(up(0,iy,ifld)-up(1,iy,ifld))/vpnorm
ccc          endif
ccc        enddo
ccc        if (isteonxy(0,iy) .eq. 1) then
ccc          iv = idxte(0,iy)
ccc          yldot(iv) = -nurlxe*(te(0,iy)-te(1,iy))*ne(0,iy)/ennorm
ccc        endif
ccc        if (istionxy(0,iy) .eq. 1) then
ccc          iv = idxti(0,iy)
ccc          yldot(iv) = -nurlxi*(ti(0,iy)-ti(1,iy))*ne(0,iy)/ennorm
ccc        endif
ccc        do igsp = 1, ngsp
ccc          if (isngonxy(0,iy,igsp) .eq. 1) then
ccc            iv = idxg(0,iy,igsp)
ccc            yldot(iv) = -nurlxg*(ng(0,iy,igsp)-ng(1,iy,igsp))/n0g(igsp)
ccc          endif
ccc        enddo
ccc        if (isphionxy(0,iy)+isphiofft .eq. 1) then
ccc          iv = idxphi(0,iy)
ccc          yldot(iv) = -nurlxp*(phi(0,iy)-phi(1,iy))/temp0
ccc        endif
ccc      enddo   # large loop over iy
ccc
ccc      endif
c ... End special coding for no divertor leg at ix = 0

 1300 continue    #jump over entire ix=0 cases for ixmnbcl=0; interior bdry

c ====================================================================
c ======================== The ix=nx+1 boundary ======================
c ====================================================================

      if (ixmxbcl .eq. 0) goto 1400   #skip setting Eqn because interior bdry

c********************************************************************
c...  First, check if isfixrb=2 for using symmetry BC at ix = nx+1
c*******************************************************************
      if (i6 .ge. nx+1 .and. isfixrb(1) .gt. 0) then  #begin symmetry BC at nx+1

       do iy = j2, j5
c...  First do the ion density
         do ifld = 1, nisp
            if(isnionxy(nx+1,iy,ifld) .eq. 1) then
               iv1 = idxn(nx+1,iy,ifld)
               yldot(iv1) = nurlxn *
     .              (nib(ifld)*nibprof(iy)-ni(nx+1,iy,ifld))/n0(ifld)
               if(isfixrb(1).eq.2) yldot(iv1) = nurlxn * (1/n0(ifld)) *
     .              (ni(nx,iy,ifld) - ni(nx+1,iy,ifld))
            endif
         enddo

c...  Now do the parallel velocity
         do ifld = 1, nusp
            if(isuponxy(nx+1,iy,ifld) .eq. 1) then
               iv2 = idxu(nx+1,iy,ifld)
               yldot(iv2) = nurlxu*(up(nx,iy,ifld) - up(nx+1,iy,ifld))/
     .                                                         vpnorm
               if(isfixrb(1).eq.2) yldot(iv2) = nurlxu*
     .                            (0. - up(nx+1,iy,ifld))/vpnorm
               if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
                  cs = sqrt( (te(nx+1,iy)+ti(nx+1,iy))/mi(ifld) )
                  yldot(iv2) = nurlxu*(cs -up(nx+1,iy,ifld))/vpnorm
               endif
            endif
         enddo

c...  now do the gas and temperatures
         if(isteonxy(nx+1,iy) .eq. 1) then
           iv1 = idxte(nx+1,iy)
           yldot(iv1) = nurlxe * ne(nx+1,iy) *
     .                     (teb*ev*tebprof(iy) - te(nx+1,iy))/ennorm
           if(isfixrb(1).eq.2) yldot(iv1) = nurlxe * ne(nx+1,iy) *
     .                               (te(nx,iy) - te(nx+1,iy))/ennorm
           if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
              yldot(iv1) = - nurlxe*(feex(nx,iy)/sx(nx,iy) - bcee*
     .                           ne(nx+1,iy)*vex(nx,iy)*te(nx+1,iy))/
     .                                 (vpnorm*ennorm)
           endif
         endif
         if(istionxy(nx+1,iy) .eq. 1) then
            iv2 = idxti(nx+1,iy)
            yldot(iv2) = nurlxi * ne(nx+1,iy) *
     .                     (tib*ev*tibprof(iy) - ti(nx+1,iy))/ennorm
            if(isfixrb(1).eq.2) yldot(iv2) = nurlxi * ne(nx+1,iy) *
     .                               (ti(nx,iy) - ti(nx+1,iy))/ennorm
            if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
               yldot(iv2) = -nurlxi*
     .          ( feix(nx,iy) - bcei*ti(nx+1,iy)*fac2sp*fnix(nx,iy,1) ) /
     .                                         (vpnorm*ennorm*sx(nx,iy))
            endif
         endif
         do igsp = 1, nhgsp # not valid for ngsp > nhgsp; only on hydrog. gas
            if (isngonxy(nx+1,iy,igsp) .eq. 1) then
               iv = idxg(nx+1,iy,igsp)
               yldot(iv) = nurlxg * (ngbackg(igsp) - 
     .                                     ng(nx+1,iy,igsp)) / n0g(igsp)
               if(isfixrb(1).eq.2) yldot(iv) = nurlxg * 
     .                     (ng(nx,iy,igsp) - ng(nx+1,iy,igsp))/n0g(igsp)
               if(isfixrb(1).eq.2 .and. yyrb(iy,1).gt.rlimiter) then
                  t1 = engbsr * max(tg(nx,iy,1),tgmin*ev)
                  vxn = 0.25 * sqrt( 8*t1/(pi*mg(igsp)) )
                  areapl = isoldalbarea*sx(nx,iy) + (1-isoldalbarea)*sxnp(nx,iy)
                  flux_inc = fac2sp*fnix(nx,iy,1)
                  if (ishymol.eq.1 .and. igsp.eq.2) then
                    if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                      flux_inc = 0.5*( fnix(nx,iy,1) + fnix(nx,iy,2) ) 
                    else
                      flux_inc = 0.5*( fnix(nx,iy,1) + fngx(nx,iy,1) ) 
                    endif
                  endif
                  yldot(iv) = -nurlxg * ( fngx(nx,iy,igsp) +
     .                                            fngxrb_use(iy,igsp,1) -
     .                      fngxsrb(iy,igsp,1) + recyrb(iy,igsp,1)*flux_inc -
     .                  (1-albrb(iy,igsp,nxpt))*ng(nx,iy,igsp)*vxn*areapl ) 
     .                                     / (vpnorm*n0g(igsp)*sx(nx,iy))
               endif
            endif
         enddo      # end of igsp loop over gas
c ... Neutral temperature - test if tg eqn is on, then set BC
	 do igsp = 1, ngsp
           if (istgonxy(nx+1,iy,igsp) == 1) then
             iv = idxtg(nx+1,iy,igsp)
             yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(nx+1,iy,igsp))/(temp0*ev)
             if(isfixrb(1)==2) then #just above applies if isfixrb=1
               yldot(iv)=nurlxg*(tg(nx,iy,igsp)-tg(nx+1,iy,igsp))/
     .                                                      (temp0*ev)
             endif
             if(isfixrb(1)==2 .and. yyrb(iy,1) > rlimiter) then
               yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                    tg(nx+1,iy,igsp))/(temp0*ev)
             endif
           endif
         enddo

         if (isphionxy(nx+1,iy) .eq. 1) then
            iv = idxphi(nx+1,iy)
            yldot(iv) = nurlxp*(phi(nx,iy) - phi(nx+1,iy))/temp0
         endif

       enddo         # end of main loop starting with iy = j2, j5

      endif         # end of ix = nx+1, isfixrb.ne.0 boundary conditions

c...  If isfixrb=2, check if i2,i5 range for yldot in pandf includes ixpt1(1);
c...  the overwrite ix=ixpt1(1) up eqn with up(ixpt1(1)) --> 0 eqn.
      if (isfixrb(1) .eq. 2) then
         if (i2.le.ixpt1(1) .and. i5.ge.ixpt1(1) .and. j2.le.iysptrx1(1)) then  
           do ifld = 1, nusp
             if (isupcore(ifld).eq.0) then
               do iy = 0+1-iymnbcl, iysptrx1(1)
                 if(isuponxy(ixpt1(1),iy,ifld)==1) then
                   iv = idxu(ixpt1(1),iy,ifld)
                   yldot(iv) = nurlxu*(0.-up(ixpt1(1),iy,ifld))/vpnorm
                 endif
               enddo
             endif
           enddo
         endif 

      endif                             #end symm. BC at ix=nx+1 for isfixrb=2

c************************************************************************
c     begin standard divertor plate conditions for right boundaries
c************************************************************************
      do jx = 1, nxpt # loop over nxpt mesh regions
      if (xcnearrb .or. openbox) then
c     First, the density equations --
        do ifld = 1, nisp
          if ( i6.ge.(ixrb(jx)+1-isextrnp) .and. isfixrb(jx)==0) then
          do iy = j2, j5

            ixt  = ixrb(jx) + 1    # analog of ix=nx+1
            ixt1 = ixm1(ixt,iy)    # analog of ix=nx
            ixt2 = ixm1(ixt1,iy)   # analog of ix=nx-1
            ixt3 = ixm1(ixt2,iy)   # analog of ix=nx-2

	    if (isnionxy(ixt,iy,ifld)==1) then
              iv1 = idxn(ixt,iy,ifld)
              if (isupgon(1)==1 .and. zi(ifld)==0.0) then   ## neutrals
                if (recyrb(iy,1,jx) .gt. 0.) then           # recycling
                  t0 = max(tg(ixt1,iy,1),tgmin*ev) 
                  vxn = 0.25 * sqrt( 8*t0/(pi*mi(ifld)) )
                  areapl = isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy)
                  yldot(iv1) = nurlxg *
     .               (fnix(ixt1,iy,ifld) + recyrb(iy,1,jx)*fnix(ixt1,iy,1) +
     .                                            fngxrb_use(iy,1,jx) -
     .                (1-albrb(iy,1,jx))*ni(ixt1,iy,ifld)*vxn*areapl
     .                - fngxsrb(iy,1,jx) ) / (vpnorm*n0(ifld)*sx(ixt1,iy))
                elseif (recyrb(iy,1,jx) <=  0. .and. 
     .                  recyrb(iy,1,jx) >= -1.) then   # recyrb is albedo
                  t0 = max(tg(ixt1,iy,1),tgmin*ev) 
                  vyn = sqrt( 0.5*t0/(pi*mi(1)) )
                  yldot(iv1) = nurlxg * ( fnix(ixt1,iy,ifld) -
     .                 (1+recyrb(iy,1,jx))*ni(ixt,iy,ifld)*vyn*sx(ixt1,iy) )/
     .                 (vpnorm*n0(ifld)*sx(ixt1,iy))
                elseif (recyrb(iy,1,jx) < -1. .and.
     .                  recyrb(iy,1,jx) > -2.) then    #fix density ngrfix
                  yldot(iv1)=nurlxg*(ngrfix - ni(ixt,iy,ifld))/n0(ifld)
                elseif (recyrb(iy,1,jx) <= -2.) then   #zero gradient
                  yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                endif # end if-test on recyrb
              else                                     ## ions
                if (isextrnp==0) then                  # zero x-gradient
                  yldot(iv1) = nurlxn*(ni(ixt1,iy,ifld)-ni(ixt,iy,ifld))/
     .                                                          n0(ifld)
                else                                   # extrapolation
                  nbound =  ni(ixt1,iy,ifld) + gxf(ixt2,iy)*
     .                   (ni(ixt1,iy,ifld)-ni(ixt2,iy,ifld))/gxf(ixt1,iy)
                  nbound = 1.2*nbound/( 1+0.5*exp( -2*(nbound/
     .                    ni(ixt,iy,ifld)-1) ) ) + 0.2*ni(ixt,iy,ifld)
                  yldot(iv1) = nurlxn*(nbound-ni(ixt,iy,ifld))/n0(ifld)
                endif # end if-test on isextrnp
              endif   # end if-test on isupgon and zi
            endif     # end if-test on isnionxy  
          enddo       # end do-loop in iy
          endif       # end if-test on i6 and isfixrb and isnion
        enddo         # end do-loop on ifld

c     Now do the parallel velocity and other variables --
      if ( (i6 .ge. (ixrb(jx)+1)) .and. isfixrb(jx)==0 ) then
      do iy = j2, j5 # begin big do-loop on iy

        ixt  = ixrb(jx) + 1    # analog of ix=nx+1
        ixt1 = ixm1(ixt,iy)    # analog of ix=nx
        ixt2 = ixm1(ixt1,iy)   # analog of ix=nx-1
        ixt3 = ixm1(ixt2,iy)   # analog of ix=nx-2

        kfeix = 0.
        sumb = 0.
        do ifld = 1, nfsp  # set up generalized Bohm condition
          upi(ixt,iy,ifld) = upi(ixt1,iy,ifld) # not set before; upi not var
          upi(ixt,iy,1   ) =  up(ixt ,iy,1   ) # need to keep up(,,1) as var
          ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .            vytan(ixt1,iy,ifld) ) / rrv(ixt1,iy) 
          sumb = sumb + ni(ixt,iy,ifld)*zi(ifld)**2*te(ixt,iy) /
     .                (mi(ifld)*(upi(ixt,iy,ifld)+ueb)**2 - ti(ixt,iy))
        enddo # end do-loop on ifld for Bohm condition
        sumb = sqrt(abs(sumb/ne(ixt,iy)))

c       Next, the momentum equations --
        do ifld = 1, nusp
          if (isuponxy(ixt,iy,ifld)==1) then
            iv2 = idxu(ixt1,iy,ifld)  #ixt1 ~ nx
	    iv = idxu(ixt,iy,ifld)    #ixt ~ nx+1
            cs = csfacrb(ifld,jx)*sqrt( (te(ixt,iy)+
     .                                  csfacti*ti(ixt,iy))/mi(ifld) )
            if (isupgon(1)==1 .and. zi(ifld)==0.0) then  ## neutrals
              if (recycmrb(iy,1,jx) > -9.9) then  # backscatter with recycm
                yldot(iv2) = -nurlxu*(recycmrb(iy,1,jx)*up(ixt1,iy,1) + 
     .                                       up(ixt1,iy,ifld))/vpnorm
              elseif (recycmrb(iy,1,jx) <= -9.9 .and. 
     .                      recycmrb(iy,1,jx) > -10.1) then # zero x-gradient
                yldot(iv2) = nurlxu*(up(ixt2,iy,ifld) -
     .                                          up(ixt1,iy,ifld))/vpnorm
              else  #neutral thermal flux to wall if recycm < -10.1
                t0 = max(tg(ixt,iy,1),tgmin*ev)
                vxn = cgmompl*0.25*sqrt( 8*t0/(pi*mi(ifld))) 
c...              if up > 0, leave unchanged; if up<0, big reduction
cc                vparn = up(ixt,iy,ifld)*( 1./(1. + 
cc     .                             exp(-up(ixt,iy,ifld)/vgmomp)) )
                vparn = up(ixt,iy,ifld)
                yldot(iv2) = -nurlxu*( fmix(ixt1,iy,ifld) - vparn*vxn*
     .                          0.5*(nm(ixt1,iy,ifld)+nm(ixt,iy,ifld))*
     .                                              sx(ixt1,iy) ) /
     .                               (vpnorm*fnorm(ifld)*sx(ixt1,iy))
              endif
              yldot(iv) =nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/vpnorm
            else                                         ## ions
              ueb = cfueb*( cf2ef*v2ce(ixt1,iy,ifld)*rbfbt(ixt,iy) -
     .                vytan(ixt1,iy,ifld) ) / rrv(ixt1,iy) 
              yldot(iv2) = -nurlxu*(sumb - 1.)  # multispecies Bohm
              if (isbohmms==0) then          # simple Bohm condition
                yldot(iv2) = nurlxu * (cs-ueb-up(ixt1,iy,ifld))/vpnorm
              endif
              if(isupss(ifld)==1 .and. up(ixt2,iy,ifld)+ueb .gt. cs)
                                             # dup/dx=0 if supersonic
     .          yldot(iv2) = nurlxu*(up(ixt2,iy,ifld)-up(ixt1,iy,ifld))/
     .                                                           vpnorm
              if (isupss(ifld)==-1) then     # slip boundary conditions
                yldot(iv2) = nurlxu*(up(ixt2,iy,ifld)-up(ixt1,iy,ifld))/
     .                                                           vpnorm
              elseif (isupss(ifld)==-2) then # extrap. + no pos. uu
                vbound = up(ixt2,iy,ifld) - gx(ixt2,iy)*
     .                   (up(ixt3,iy,ifld)-up(ixt2,iy,ifld))/gx(ixt1,iy) 
                vbound = max(vbound, -ueb)   # forces uu & fnix >= 0
                yldot(iv2) = nurlxu*(vbound - up(ixt1,iy,ifld))/vpnorm
              elseif (isupss(ifld)==-3) then # modified Bohm condition
                vbound = -ueb +2*cs*uu(ixt2,iy,ifld)/
     .                             (uu(ixt2,iy,ifld)+rrv(ixt,iy)*cs)
                vbound = max(vbound, -ueb)   # forces uu & fnix >= 0
                yldot(iv2) = nurlxu * (vbound-up(ixt1,iy,ifld))/vpnorm
              endif # end if-test on isupss

c Finally set unused up(ixt,,) = up(ixt1,,); note ixt~nx+1
              yldot(iv) =nurlxu*(up(ixt1,iy,ifld)-up(ixt,iy,ifld))/vpnorm
            endif # end if-test on isupgon
          endif # end if-test on isupon
          if (zi(ifld)==0.0) then 
            kfeix = kfeix - cftiexclg*cfvcsx(ifld)*0.5*sx(ixt1,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt2,iy,ifld)**2 ) 
          else
            kfeix = kfeix - cfvcsx(ifld)*0.5*sx(ixt1,iy)
     .                     *visx(ixt1,iy,ifld)*gx(ixt1,iy)
     .           *( up(ixt1,iy,ifld)**2 - up(ixt2,iy,ifld)**2 )
          endif
        enddo # end do-loop on ifld

c       Next, the potential equation --
        if(isphion+isphiofft .eq. 1) then
cc           fqpsatrb(iy,jx) = qe*isfdiax*( ne(ixt,iy)*v2ce(ixt1,iy,1)*
cc     .                           rbfbt(ixt,iy)*sx(ixt1,iy) + fdiaxrb(iy,jx) )
cc           do ifld = 1, nfsp   # note fqp,fqpsat are poloidal proj. of || curr
cc            fqpsatrb(iy,jx) = fqpsatrb(iy,jx) + qe*zi(ifld)*ni(ixt,iy,ifld)*
cc     .                        upi(ixt1,iy,ifld)*sx(ixt1,iy)*rrv(ixt1,iy)  
cc           enddo
           lambdae = 2e16*(te(ixt,iy)/ev)**2/ne(ixt,iy)  #approx mfp
           kincorrb(iy,jx) = 1./(1 + cfkincor*(lambdae/lcone(ixt,iy))*
     .                                       abs(ev*phi(ixt,iy)/te(ixt,iy)))
           fqpsate = qe*ne(ixt,iy)*sqrt(te(ixt,iy)/(2*pi*me))*
     .                               kincorrb(iy,jx)*sx(ixt1,iy)*rrv(ixt1,iy)
c          NOTE: by definition, fqpsate is always > 0
           if (ikapmod==0) then
cc           if (fqp(ixt1,iy) > 0.) then #limit to saturation current
cc             fqp(ixt1,iy)= ( abs((fqp(ixt1,iy)*fqpsatrb(iy,jx)))**exjbdry/
cc     .                  (abs(fqp(ixt1,iy))**exjbdry +
cc     .                   abs(fqpsatrb(iy,jx))**exjbdry) )**(1/exjbdry)
cc           endif
              if ( fqpsatrb(iy,jx)-(1.-gamsec)*fqp(ixt1,iy) > 0) then  # force +ve log argument
                arglgphi=(((fqpsatrb(iy,jx)-(1.-gamsec)*fqp(ixt1,iy))
     .                             /fqpsate)**2 + expkmx**2)**(0.5)
              else
                arglgphi = expkmx
              endif
              if (iskaprex.eq.0) kappar(iy,jx) = - log(arglgphi)
cccTDR            if ((isnewpot==1) .and. ((iy==1) .or. (iy==ny))) then
           elseif (ikapmod==1) then
              kappar(iy,jx) = kappa( fqpsatrb(iy,jx), fqpsate, 
     .                                       (1.-gamsec)*fqp(ixt1,iy) )
           endif
           if (isphionxy(ixt,iy) .eq. 1) then
              iv = idxphi(ixt,iy)
              if ((isnewpot==100) .and. ((iy==1) .or. (iy==ny))) then
                 continue
              elseif (isphirbc==1)
                 yldot(iv) = - nurlxp * (phi(ixt,iy) - phi0r(iy,jx))/temp0 
              else
                 yldot(iv) = -nurlxp*(1.-bctype(iy))*
     .              ( phi(ixt,iy)-kappar(iy,jx)*te(ixt,iy)/ev
     .                                     -phi0r(iy,jx) )/temp0
     .                       -nurlxp*bctype(iy)*(1.-gamsec)*fqp(ixt1,iy)
     .                                 / (fqpsatrb(iy,jx)+cutlo)
              endif
           endif
	else    # corresponding to (isphion+isphiofft .ne. 1)
           fqpsate = 0.
	   kappar(iy,jx) = 3.  # estimate for e*phi/Te; may be used for sputtering
        endif
        isphion2 = isphion + isphiofft
        bcer(iy,jx) = (1-newbcr(jx)*isphion2) * bcee
     .              + newbcr(jx)*isphion2 * (2. + kappar(iy,jx))
        if (iskaprex.eq.1) bcer(iy,jx) = 2. + kappar(iy,jx)
        bcir(iy,jx) = (1-newbcr(jx)*isphion2) * bcei
     .              + newbcr(jx)*isphion2 * (2.5)

        t0 = te(ixt,iy)/ev
        f_cgpld = .5*(1. - cos(pi*(t0 - temin)/(.3 - temin)))
                                # ramp cgpld effect down at low guard cell
                                # Te; rampdown might not be necesasry if
                                # guard cell were allowed to access negative
                                # temperatures (ETM 7 Jan 2014)
        if(t0 < temin) f_cgpld = 0.
        if(t0 > 0.3) f_cgpld = 1.
        t0 = max(tg(ixt1,iy,1),tgmin*ev)
        vxn = f_cgpld * 0.25 * sqrt(8*t0/(pi*mg(1)))

c   Do the electron temp Eqn -----------------------------------
        if (isteonxy(ixt,iy) == 1) then
	  iv1 = idxte(ixt,iy)
          if (ibctepr == 1) then
             faceel =  bcer(iy,jx)*(fqpsate/qe)*exp(-kappar(iy,jx))
             faceel2 = bcer(iy,jx)*(fqpsate/qe)*exp(-kappamx+2) 
             totfeexr(iy,jx) = feex(ixt1,iy) + cfeexdbo*( 
     .             2.5*fdiaxrb(iy,jx) + floxebgt(ixt1,iy) )*te(ixt,iy)
cc           if (feex(ixt1,iy) < 0.) then #if totfeex<0;force te(ixt)=.5*te(ixt1)
cc             totfeex = 0.1*faceel*te(ixt1,iy)
cc           else
cc             totfeex =(feex(ixt1,iy)**2+(0.1*faceel*te(ixt1,iy))**2)**0.5
cc           endif
	     totfnex = ne(ixt,iy)*vex(ixt1,iy)*sx(ixt1 ,iy)
             if (isphion+isphiofft==1) then
              yldot(iv1) = nurlxe*(totfeexr(iy,jx) 
     .                       -faceel*te(ixt,iy)
     .                       -faceel2*(te(ixt,iy)-te(ixt1,iy))
     .                       -cmneut*fnix(ixt1,iy,1)*recycp(1)*eedisspr*ev
     .                                  )/(sx(ixt1,iy)*vpnorm*ennorm)
              else
               yldot(iv1) = nurlxe*(totfeexr(iy,jx)
     .                        -totfnex*te(ixt,iy)*bcer(iy,jx)
     .                        -cgpld*sx(ixt1,iy)*0.5*ng(ixt1,iy,1)*vxn*ediss*ev
     .                        -cmneut*fnix(ixt1,iy,1)*recycp(1)*eedisspr*ev
     .                                   )/(sx(ixt1,iy)*vpnorm*ennorm)
             endif
          elseif (ibctepr .eq. 0) then
             yldot(iv1) = nurlxe*(tepltr*ev-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctepr .eq. 2) then
             yldot(iv1) = nurlxe*(te(ixt1,iy)-te(ixt,iy))
     .                                    * 1.5 * ne(ixt,iy)/ennorm
          endif 
        endif  # end loop for isteon=1

c  Do the ion temp Eqn ---------------------------
        if (istionxy(ixt,iy) == 1) then
          iv2 = idxti(ixt,iy)
          if (ibctipr == 1) then
c            totfeix is thermal + kinetic + viscous flux
            totfeixr(iy,jx) = feix(ixt1,iy) + ckinfl*kfeix 
            totfnix = 0.
            do ifld = 1, nfsp
            # the condition on zi must include use of v2cd which is not properly set to zero in pandf1
            if(zi(ifld)>1e-10) then
              totfeixr(iy,jx) = totfeixr(iy,jx) + cfeixdbo*(
     .                    2.5* ni(ixt,iy,ifld)*v2cd(ixt1,iy,ifld)*
     .                          sx(ixt1,iy)*rbfbt(ixt,iy)+
     .                           floxibgt(ixt1,iy,ifld) )*ti(ixt,iy)
               totfnix = totfnix + fnix(ixt1,iy,ifld)
               endif
            enddo
            if (isupgon(1)==1) then
c              Different boundary conditions for neutral momentum equation
cc               if (recyrb(iy,1,jx) .gt. 0.) then
cc                  bcen = recyce*bcir(iy,jx) - 0.5*mi(1)*(up(ixt1,iy,2)**2 -
cc     .              recyce*upi(ixt1,iy,1)**2)/(recyrb(iy,1,jx)*ti(ixt,iy))
cc               else
cc                  bcen = 0.
cc               endif
cc               if (recyce .le. 0) bcen = 0.  # gets back to old case

                 yldot(iv2) = nurlxi*(totfeixr(iy,jx)
     .                        -totfnix*bcir(iy,jx)*ti(ixt,iy)+
     .            cftiexclg*( -cfneut*fnix(ixt1,iy,iigsp)*bcen*tg(ixt,iy,1)
     .                        -(cgengpl*2.*tg(ixt,iy,1) - cgpld*eion*ev)*
     .                                 ng(ixt1,iy,1)*vxn*sx(ixt1,iy) 
     .                        -cmneut*fnix(ixt1,iy,1)*recycp(1)*
     .                                 cmntgpr*(ti(ixt,iy)-eidisspr*ev)
     .                                  ) )/(vpnorm*ennorm*sx(ixt1,iy))
            else
                 yldot(iv2) = nurlxi*(totfeixr(iy,jx) 
     .                        -totfnix*bcir(iy,jx)*ti(ixt,iy)+
     .            cftiexclg*( -cmneut*fnix(ixt1,iy,1)*recycp(1)*
     .                                 cmntgpr*(ti(ixt,iy)-eidisspr*ev)    
     .                                   ) )/(vpnorm*ennorm*sx(ixt1,iy))
            endif # end test for isupgon(1)==1

          elseif (ibctipr .eq. 0) then
             yldot(iv2) = nurlxi*(tipltr*ev - ti(ixt,iy))
     .                                  * 1.5 * ne(ixt,iy)/ennorm
          elseif (ibctipr .eq. 2) then
             yldot(iv2) = nurlxi*(ti(ixt1,iy) - ti(ixt,iy))
     .                                  * 1.5 * ne(ixt,iy)/ennorm
          endif # end loop for bctipr
        endif  # end loop for istion=1

c       Next, the hydrogenic gas equations --
        do igsp = 1, nhgsp  # imp gas below
           if (isngonxy(ixt,iy,igsp) .eq. 1) then
             iv = idxg(ixt,iy,igsp)
             if (recyrb(iy,igsp,jx) .gt. 0.) then  # normal recycling
               flux_inc = fac2sp*fnix(ixt1,iy,1)
               if (ishymol.eq.1 .and. igsp.eq.2) then
                 if (isupgon(1) .eq. 1) then  # two atoms for one molecule
                   flux_inc = 0.5*( fnix(ixt1,iy,1) + fnix(ixt1,iy,2) ) 
                 else
                   flux_inc = 0.5*( fnix(ixt1,iy,1) + fngx(ixt1,iy,1) ) 
                 endif
               endif
               t0 = max(tg(ixt1,iy,igsp), tgmin*ev)
               vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
               areapl = isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy) 
               yldot(iv) = nurlxg *  ( fngx(ixt1,iy,igsp) +
     .                                          fngxrb_use(iy,igsp,jx) -
     .               fngxsrb(iy,igsp,jx) + recyrb(iy,igsp,jx)*flux_inc -
     .               (1-albrb(iy,igsp,jx))*ng(ixt1,iy,igsp)*vxn*areapl )
     .                                  / (vpnorm*n0g(igsp)*sx(ixt1,iy))
             elseif (recyrb(iy,igsp,jx) <=  0. .and.
     .               recyrb(iy,igsp,jx) >= -1.) then  # recyrb is albedo
               t0 = max(tg(ixt,iy,igsp), tgmin*ev)
               vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
               yldot(iv) =  nurlxg*( fngx(ixt1,iy,igsp) -
     .          (1+recyrb(iy,igsp,jx))*ng(ixt,iy,igsp)*vxn*sx(ixt1,iy) )
     .                                     / (vxn*sx(ixt1,iy)*n0g(igsp))
             elseif (recyrb(iy,igsp,jx) < -1. .and.
     .               recyrb(iy,igsp,jx) > -2.)  then  #fix density ngrfix
               yldot(iv)=nurlxg*(ngrfix - ng(ixt,iy,igsp))/n0g(igsp)
             elseif (recyrb(iy,igsp,jx) <= -2.) then   #zero gradient
               yldot(iv) = nurlxg*(ng(ixt1,iy,igsp) -
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
             endif
             if (is1D_gbx.eq.1) yldot(iv) = nurlxg*(ng(ixt1,iy,igsp) -
     .                                      ng(ixt,iy,igsp))/n0g(igsp)
c   Special coding for Maxim
             if (isngrf==1) yldot(iv) = nurlxg*(ngrfix - 
     .                                    ng(ixt,iy,igsp))/n0g(igsp)
           endif # end if-test on isngon
        enddo # end do-loop on igsp

c       Finally, do impurity gas equations --
c       Set neutral gas impurity flux to be the sum of the fluxes of the
c       sputtered impurities plus recycled impurities from all charge states.
        if (isimpon .ge. 4 .and. isimpon .le. 7 .and. nzspt.ge.1) then
           nzsp_rt = nhsp
           do igsp = nhgsp+1, ngsp
	      jz = max(igsp - nhgsp, 1)   # identify impurity index
	      if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fnix
              if (isngonxy(ixt,iy,igsp) .eq. 1) then
                 iv = idxg(ixt,iy,igsp)
                 hflux = 0.
                 do ihyd = 1, nhsp
                    if (zi(ihyd).gt.0.) hflux = hflux + fnix(ixt1,iy,ihyd)
                 enddo
                 zflux = 0.
                 do iimp = 1, nzsp(jz) # need to limit this loop per species
                   zflux = zflux + fnix(ixt1,iy,nzsp_rt+iimp)
                 enddo
                 sputflxrb(iy,igsp,jx) = 0.
                 zflux_chm = 0.
                 if (isph_sput(igsp) .ge. 1) then  # use fits for phys sput
                   do ifld = ipsputt_s, ipsputt_e 
                     eng_sput = ( 0.5*mi(ifld)*up(ixt1,iy,ifld)**2 + 
     .                             ti(ixt,iy) + zi(ifld)*
     .                             kappar(iy,jx)*te(ixt,iy) )/ev
                     if(zi(ifld)>0.) sputflxrb(iy,igsp,jx) = sputflxrb(iy,igsp,jx) +
     .                                             fnix(ixt1,iy,ifld)*
     .                            fphysyrb(igsp,jx)*yld96(matp,matt,eng_sput)
                   enddo
                   if (isph_sput(igsp) .ge. 2) then  # add chem sput from ions
                      do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                        eng_sput = ( 0.5*mi(ifld)*up(ixt1,iy,ifld)**2 + 
     .                               ti(ixt,iy) + zi(ifld)*
     .                               kappar(iy,jx)*te(ixt,iy) )/ev
                        flx_incid = abs(fnix(ixt1,iy,ifld))/sx(ixt1,iy)
                        call sputchem (isch_sput(igsp),eng_sput,tvplatrb(iy,jx),
     .                             flx_incid, yld_chm)
                        sputflxrb(iy,igsp,jx) = sputflxrb(iy,igsp,jx) + 
     .                      fchemyrb(igsp,jx)*fnix(ixt1,iy,ifld)*yld_chm
                      enddo
                   endif   #test on isph_sput >=2
                   if (isph_sput(igsp) .eq. 3) then # add chem sput from h. neut
                     do igsp2 = 1, 1+ishymol
                       t0p = max(tg(ixt1,iy,igsp2),tgmin*ev)
                       flx_incid=ng(ixt,iy,igsp2)*.25*
     .                                        sqrt(8*t0p/(pi*mg(igsp2)))
                       call sputchem (isch_sput(igsp),t0p/ev,tvplatrb(iy,jx),
     .                                               flx_incid, yld_chm)
                       zflux_chm = zflux_chm + flx_incid*
     .                             fchemyrb(igsp,jx)*yld_chm*sx(ixt1,iy)
                     enddo
                   endif    #test on isph_sput = 3
                 endif   #test on isph_sput > 1              
                 if (sputtrb(iy,igsp,jx) .ge. 0. .or.
     .                              abs(sputflxrb(iy,igsp,jx)).gt.0.) then
                    t0 = max(cdifg(igsp)*tg(ixt1,iy,igsp), tgmin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    areapl = isoldalbarea*sx(ixt1,iy) + (1-isoldalbarea)*sxnp(ixt1,iy)
                    zflux = - sputtrb(iy,igsp,jx) * hflux - 
     .                        sputflxrb(iy,igsp,jx) -
     .                   recyrb(iy,igsp,jx) * zflux +
     .                (1-albrb(iy,igsp,jx))*ng(ixt1,iy,igsp)*vxn*areapl-
     .                   zflux_chm + fngxsrb(iy,igsp,jx)-fngxrb_use(iy,igsp,jx)
                    yldot(iv) = nurlxg * (fngx(ixt1,iy,igsp) - zflux) /
     .                         (n0(igsp) * vpnorm * sx(ixt1,iy))
                 elseif (sputtrb(iy,igsp,jx).ge.-9.9) then # neg. sputtrb ==> albedo
                    t0 = max(cdifg(igsp)*tg(ixt1,iy,igsp), tgmin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    yldot(iv) =  nurlxg*( fngx(ixt1,iy,igsp) -
     .                (1+sputtrb(iy,igsp,jx))*ng(ixt,iy,igsp)*vxn*sx(ixt1,iy) )
     .                                       / (vxn*sx(ixt1,iy)*n0g(igsp))
                 else                                # sputtrb < -9.9 ==> fix dens
                    yldot(iv) = -nurlxg*(ng(ixt,iy,igsp)-ngplatrb(igsp,jx))/
     .                                                        n0g(igsp)    
                 endif  # end if-test on sputtrb
              endif  # end if-test on isngon
           enddo  # end do-loop in igsp
        endif  # end if-test on isimpon and nzspt

c ... Neutral temperature - test if tg eqn is on, then set BC
	do igsp = 1, ngsp
          if (istgonxy(ixt,iy,igsp) == 1) then
            iv = idxtg(ixt,iy,igsp)
            if (istgrb(igsp) == 0) then  #set tg=tgwall
              yldot(iv) = nurlxg*(tgwall(igsp)*ev -
     .                                 tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istgrb(igsp) == 1) then  #extrapolate tg from interior
              tbound = tg(ixt1,iy,igsp) + gxf(ixt2,iy)*
     .               (tg(ixt1,iy,igsp)-tg(ixt2,iy,igsp))/gxf(ixt1,iy)
              tbound = max(tbound,0.5*temin*ev)
              yldot(iv) = nurlxg*(tbound - tg(ixt,iy,igsp))/(temp0*ev)
            elseif (istgrb(igsp) == 2)  #placeholder for gradient BC
              call xerrab("**INPUT ERROR: istgrb=2 grad opt not implemented")
            elseif (istgrb(igsp) == 3)  #Maxwell thermal flux to wall
              t0 = max(cdifg(igsp)*tg(ixt1,iy,igsp), temin*ev)
              vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
              yldot(iv) =  nurlxg*( fegx(ixt1,iy,igsp) - 2*cgengmpl*
     .                            ng(ixt1,iy,igsp)*vxn*t0*sx(ixt1,iy) )/
     .                                     (sx(ixt1,iy)*vpnorm*ennorm)
            elseif (istgrb(igsp) == 4) 
	      if (isupgon(igsp)==1) then
                if (recyrb(iy,igsp,jx) .gt. 0.) then
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                  fng_alb = (1-albrb(iy,igsp,jx))*ng(ixt1,iy,igsp)
     .                                           *vxn*sx(ixt1,iy)
                  yldot(iv) = nurlxg * ( fegx(ixt1,iy,igsp)
     .                                  -cfalbedo*fng_alb*t0
     .                                  +recyrb(iy,igsp,jx)*(1.-cfdiss)
     .                                  *fnix(ixt1,iy,1)
     .                                  *recyce*cfalbedo
     .                                  *( kappar(iy,jx)*zi(1)*te(ixt,iy)
     .                                    +ti(ixt,iy) ) )/
     .                                  (vpnorm*ennorm*sx(ixt1,iy))
                elseif (recyrb(iy,igsp,jx) <=  0. .and.
     .                  recyrb(iy,igsp,jx) >= -1.) then   # recyrb is albedo
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
                  fng_alb = (1+recyrb(iy,igsp,jx))*ng(ixt1,iy,igsp)
     .                                            *vyn*sx(ixt1,iy)
                  yldot(iv) = nurlxg * ( fegx(ixt1,iy,igsp)
     .                                  -cfalbedo*fng_alb*t0 )/
     .                                   (vpnorm*ennorm*sx(ixt1,iy))
                elseif (recyrb(iy,igsp,jx) < -1.) then  #..half Maxwellian
                  t0 = max(tg(ixt1,iy,igsp),tgmin*ev)
                  vyn = sqrt( 0.5*t0/(pi*mg(igsp)) )
                  yldot(iv) = nurlxg * ( fegx(ixt1,iy,igsp)
     .                                  -cfalbedo*fnix(ixt1,iy,iigsp)*t0 )/
     .                                   (vpnorm*ennorm*sx(ixt1,iy))
                endif
              endif
            elseif (istgrb(igsp) == 5) then  #set tg=ti*cftgtipltr
              yldot(iv) = nurlxg*(ti(ixt,iy)*cftgtipltr(igsp) -
     .                                 tg(ixt,iy,igsp))/(temp0*ev)
            else
              call xerrab("**INPUT ERROR: istgrb set to unknown option")
            endif
          endif
        enddo

      enddo # end big do-loop on iy

c ...   Make electric field at plate uniform over last 2 cells --
        if (isexunif.eq.1) then
           do iy = j2, j5
	     if(isphionxy(ixt1,iy)==1) then
               iv = idxphi(ixt1,iy)  # next formula assumes uniform mesh
               yldot(iv)=-nurlxp*(3*phi(ixt1,iy)-2*phi(ixt,iy)-
     .                                            phi(ixt2,iy))/temp0
             endif
           enddo
        endif

      endif # end if-test on i6 and isfixrb

      endif # end if-test on xcnearrb and openbox
      enddo # end do-loop over nxpt mesh regions
c************************************************************************
c     end standard divertor plate conditions for right boundaries
c************************************************************************

c ... Special case for no divertor leg (ixpt2(1)=nx, iy.le.iysptrx2(1)), just
c ... use continuation boundary conditions for the guard cell (ix=nx+1)
ccc      if (ixpt2(1).eq.nx .and. isfixrb(1).eq.0) then  # should test i6.ge.nx+1
ccc
ccc      do iy = j2, min(j5, iysptrx2(1))
ccc        do ifld = 1, nisp
ccc          if (isnion(ifld) .eq. 1) then
ccc            iv = idxn(nx+1,iy,ifld)
ccc            yldot(iv) = -nurlxn*(ni(nx+1,iy,ifld)-ni(nx,iy,ifld))/n0(ifld)
ccc          endif
ccc        enddo
ccc        do ifld = 1, nusp
ccc          if (isupon(ifld) .eq. 1) then
ccc            iv = idxu(nx+1,iy,ifld)
ccc            yldot(iv) = -nurlxu*(up(nx+1,iy,ifld)-up(nx,iy,ifld))/vpnorm
ccc          endif
ccc        enddo
ccc        if (isteon .eq. 1) then
ccc          iv = idxte(nx+1,iy)
ccc          yldot(iv) = -nurlxe*(te(nx+1,iy)-te(nx,iy))*ne(nx+1,iy)/ennorm
ccc        endif
ccc        if (istion .eq. 1) then
ccc          iv = idxti(nx+1,iy)
ccc          yldot(iv) = -nurlxi*(ti(nx+1,iy)-ti(nx,iy))*ne(nx+1,iy)/ennorm
ccc        endif
ccc        do igsp = 1, ngsp
ccc          if (isngon(igsp) .eq. 1) then
ccc            iv = idxg(nx+1,iy,ifld)
ccc            yldot(iv) = -nurlxg*(ng(nx+1,iy,igsp)-ng(nx,iy,igsp))/n0g(igsp)
ccc          endif
ccc        enddo
ccc        if (isphion+isphiofft .eq. 1) then
ccc          iv = idxphi(nx+1,iy)
ccc          yldot(iv) = -nurlxp*(phi(nx+1,iy)-phi(nx,iy))/temp0
ccc        endif
ccc      enddo   # large loop over iy
ccc
ccc      endif
c ... End special coding for no divertor leg at ix = nx

 1400 continue    #jump over all ix=nx+1 if ixmxbcl=0 since interior bdry

c************************************************************************
c...  do ix = nxc (normally nx/2) boundary if geometry=dnbot (double null)
c************************************************************************

      if((isudsym==1.or.(geometry .eq. "dnXtarget"))
     &                                     .and. isfixlb(1).eq.0) then
         if (i2.le.nxc+1 .and. i5.ge.nxc-1) then
c...  this "if" test assumes both xlinc and xrinc are at least 1
            do 194 iy = j1, j6  # was j2, j5, but this missed corners
c...  First do the ion density
               do 192 ifld = 1, nisp
                 if(isnionxy(nxc,iy,ifld)*isnionxy(nxc+1,iy,ifld)==1) then
                   if (.not. ((isnicore(ifld)==3).and.(iy==0))) then
                     # do not over-write corner cell b.c. from iy=0
                     # because it connects inboard and outboard core te
		     # force ni(nxc+1,,) = ni(nxc,,) if isnicore(ifld)=3
                     iv = idxn(nxc,iy,ifld)
                     iv2 = idxn(nxc+1,iy,ifld)
                     if(isnicore(ifld)==3 .and. iy==0) then #ni=ni(nxc+1
                       yldot(iv) = nurlxn*( ni(nxc+1,iy,ifld) -
     .                                        ni(nxc,iy,ifld) )/n0(ifld) 
                     else                                   #ni=ni(nxc-1
                       yldot(iv) = nurlxn*( ni(nxc-1,iy,ifld) -
     .                                        ni(nxc,iy,ifld) )/n0(ifld)
                     endif
                     yldot(iv2) = nurlxn *
     .                    (ni(nxc+2,iy,ifld)-ni(nxc+1,iy,ifld))/n0(ifld)
                   endif
                 endif
 192           continue
c...  Now do the parallel velocity
               do ifld = 1, nusp
                  if(isuponxy(nxc,iy,ifld)*isuponxy(nxc-1,iy,ifld)*
     .                           isuponxy(nxc+1,iy,ifld) .eq. 1) then 
                     iv3 = idxu(nxc,iy,ifld)
                     yldot(iv3) = nurlxu*(0.-up(nxc,iy,ifld))/vpnorm
                     iv3 = idxu(nxc-1,iy,ifld)
                     yldot(iv3) = nurlxu*(0.-up(nxc-1,iy,ifld))/vpnorm
                     iv3 = idxu(nxc+1,iy,ifld)
                     yldot(iv3) = nurlxu*(0.-up(nxc+1,iy,ifld))/vpnorm
                  endif
               enddo

               if(isteonxy(nxc,iy)*isteonxy(nxc+1,iy) .eq. 1) then
                  if ( .not. ((iflcore==1).and.(iy==0)) ) then
                     # do not over-write corner cell b.c. from iy=0
                     # because it connects inboard and outboard core te
                     iv = idxte(nxc,iy)
                     iv2 = idxte(nxc+1,iy)
                     yldot(iv) = nurlxe*(te(nxc-1,iy) - te(nxc,iy))
     .                            *1.5*ne(nxc,iy)/ennorm
                     yldot(iv2) = nurlxe*(te(nxc+2,iy) - te(nxc+1,iy))
     .                            *1.5*ne(nxc+1,iy)/ennorm
                  endif
               endif
               if(istionxy(nxc,iy)*istionxy(nxc+1,iy) .eq. 1) then
                  if ( .not. ((iflcore==1).and.(iy==0)) ) then
                     # do not over-write corner cell b.c. from iy=0
                     # because it connects inboard and outboard core ti
                     iv = idxti(nxc,iy)
                     iv2 = idxti(nxc+1,iy)
                     yldot(iv) = nurlxi*(ti(nxc-1,iy) - ti(nxc,iy))
     .                            *1.5*ne(nxc,iy)/ennorm
                     yldot(iv2) = nurlxi*(ti(nxc+2,iy) - ti(nxc+1,iy))
     .                            *1.5*ne(nxc+1,iy)/ennorm
                  endif
               endif
   
               do igsp = 1, ngsp
                  if(isngonxy(nxc,iy,igsp)*isngonxy(nxc+1,iy,igsp)==1) then
                     iv =  idxg(nxc  ,iy,igsp)
                     iv2 = idxg(nxc+1,iy,igsp)
                     yldot(iv ) = nurlxg*(ng(nxc-1,iy,igsp) - 
     .                                    ng(nxc  ,iy,igsp)) / n0g(igsp)
                     yldot(iv2) = nurlxg*(ng(nxc+2,iy,igsp) - 
     .                                    ng(nxc+1,iy,igsp)) / n0g(igsp)
                  endif
                  if(istgonxy(nxc,iy,igsp)*istgonxy(nxc+1,iy,igsp)==1) then
                     iv =  idxtg(nxc  ,iy,igsp)
                     iv2 = idxtg(nxc+1,iy,igsp)
                     yldot(iv ) = nurlxg*(tg(nxc-1,iy,igsp) - 
     .                                    tg(nxc  ,iy,igsp))*
     .                                1.5*ng(nxc,iy,igsp) / ennorm
                     yldot(iv2) = nurlxg*(tg(nxc+2,iy,igsp) - 
     .                                    tg(nxc+1,iy,igsp))*
     .                                1.5*ng(nxc+1,iy,igsp)/ ennorm
                  endif
               enddo

c...  Do boundary condition for potential along ix=nxc and ix=nxc+1
               if (isphionxy(nxc,iy)*isphionxy(nxc+1,iy)==1) then
                  if ( .not. (isnewpot==1 .and. (iy==0 .or. iy==1)) ) then
                     # do not over-write BC from iy=0 & iy=1 because
                     # these connect inboard and outboard phi values
                     iv = idxphi(nxc,iy)
                     iv2 = idxphi(nxc+1,iy)
                     yldot(iv) = nurlxi*(phi(nxc-1,iy) - phi(nxc,iy))
     .                            *1.5*ne(nxc,iy)/ennorm
                     yldot(iv2) = nurlxi*(phi(nxc+2,iy) - phi(nxc+1,iy))
     .                            *1.5*ne(nxc+1,iy)/ennorm
                  endif
               endif
   
c...  Do boundary condition for impurities along ix=nxc
               do ifld = 1, nzspt
                  if (isimpon .ge. 3 .and. isimpon .le. 7 .and.
     .                      isnionxy(nxc,iy,nhsp+ifld)*
     .                      isnionxy(nxc+1,iy,nhsp+ifld) .eq. 1 ) then
                     iv  = idxn(nxc  ,iy,nhsp+ifld)
                     iv2 = idxn(nxc+1,iy,nhsp+ifld)
                     if(isnicore(nhsp+ifld)==3 .and. iy==0) then #ni=ni(nxc+1
                       yldot(iv) = nurlxn*( ni(nxc+1,iy,nhsp+ifld) -
     .                                         ni(nxc,iy,nhsp+ifld) ) 
     .                                                 /n0(nhsp+ifld)
                     else                                   #ni=ni(nxc-1
                       yldot(iv) = nurlxn*( ni(nxc-1,iy,nhsp+ifld) -
     .                                         ni(nxc,iy,nhsp+ifld) )
     .                                                 /n0(nhsp+ifld)
                     endif
                     yldot(iv2) = nurlxn*( ni(nxc+2,iy,nhsp+ifld) -
     .                                     ni(nxc+1,iy,nhsp+ifld) )
     .                                                 /n0(nhsp+ifld)
                  endif
               enddo            # end of impurities along ix=nxc

 194        continue
         endif
         endif                  # end of ix = nxc b.c. for double-null
         
cc    if (islimon .ne. 0) call subroutine limterbc  #to-be-done subroutine
cc                                                  #to replace next
c...############## Include bdry cond for limiter if present ##########
c...##################################################################
      if( islimon .ne. 0) then  #if-test 8888
        if (i2.le.ix_lim+1 .and. i5.ge.ix_lim-1) then #if-test 8889
c...     last "if" test assumes both xlinc and xrinc are at least 1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...  For flux tubes that do not intersect the limiter:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do 198 iy = j2, min(j5,iy_lims-1)
c...  First do the ion density
            do ifld = 1, nisp
              if(isnionxy(ix_lim,iy,ifld)*
     .                     isnionxy(ix_lim+1,iy,ifld) .eq. 1) then
                iv = idxn(ix_lim,iy,ifld)
                iv2 = idxn(ix_lim+1,iy,ifld)
                if (islbcn .eq. 0) then
ccc islbcn=0:
                  yldot(iv) = nurlxn * (-ni(ix_lim,iy,ifld)
     .              +0.5*(ni(ix_lim-1,iy,ifld)+ni(ix_lim+2,iy,ifld)))/
     .                                                      n0(ifld)
                  yldot(iv2) = nurlxn * (-ni(ix_lim+1,iy,ifld)
     .              +0.5*(ni(ix_lim-1,iy,ifld)+ni(ix_lim+2,iy,ifld)))/
     .                                                      n0(ifld)
                elseif (islbcn .eq. 1) then
ccc islbcn=1:
                  yldot(iv) = nurlxn * (-ni(ix_lim,iy,ifld)
     .              +0.5*(ni(ix_lim-1,iy,ifld)+ni(ix_lim+2,iy,ifld)))/
     .                                                      n0(ifld)
                   yldot(iv2) = nurlxn * (-ni(ix_lim+1,iy,ifld)
     .                                    +ni(ix_lim,iy,ifld))/
     .                                                      n0(ifld)
                elseif (islbcn .eq. 2) then
ccc islbcn=2:
                  yldot(iv) = nurlxn * (-fnix(ix_lim-1,iy,ifld)
     .                                   +fnix(ix_lim+1,iy,ifld))/
     .                             (sx(ix_lim-1,iy)*n0(ifld)*vpnorm)
                  yldot(iv2) = nurlxn * (-ni(ix_lim+1,iy,ifld)
     .                                    +ni(ix_lim,iy,ifld))/
     .                                                      n0(ifld)
                endif
              endif
            enddo
c...........................
c...  Now do the parallel velocity for iy beyond limiter
            do ifld = 1, nusp
              if(isuponxy(ix_lim-1,iy,ifld)*isuponxy(ix_lim,iy,ifld)*
     .           isuponxy(ix_lim+1,iy,ifld) .eq. 1) then
                   iv1 = idxu(ix_lim-1,iy,ifld)
                   iv2 = idxu(ix_lim  ,iy,ifld)
                   iv3 = idxu(ix_lim+1,iy,ifld)
                if (islbcu .eq. 0) then
ccc islbcu=0:
                   yldot(iv1) = nurlxu * (-up(ix_lim-1,iy,ifld)
     .               +0.5*(up(ix_lim-2,iy,ifld)+up(ix_lim+2,iy,ifld)))/
     .                                                         vpnorm
                   yldot(iv2) = nurlxu * (-up(ix_lim,iy,ifld)
     .               +0.5*(up(ix_lim-2,iy,ifld)+up(ix_lim+2,iy,ifld)))/
     .                                                         vpnorm
                   yldot(iv3) = nurlxu * (-up(ix_lim+1,iy,ifld)
     .               +0.5*(up(ix_lim-2,iy,ifld)+up(ix_lim+2,iy,ifld)))/
     .                                                         vpnorm
                elseif (islbcu .eq. 1) then
ccc islbcu=1:
                   yldot(iv1) = nurlxu * (-up(ix_lim-1,iy,ifld)
     .               +0.5*(up(ix_lim-2,iy,ifld)+up(ix_lim+2,iy,ifld)))/
     .                                                         vpnorm
                   yldot(iv2) = nurlxu * (-up(ix_lim,iy,ifld)
     .                               +up(ix_lim-1,iy,ifld))/vpnorm
                   yldot(iv3) = nurlxu * (-up(ix_lim+1,iy,ifld)
     .                               +up(ix_lim,iy,ifld))/vpnorm
                elseif (islbcu .eq. 2) then
ccc islbcu=2:
                   yldot(iv1) = nurlxu * (-fmix(ix_lim-1,iy,ifld)
     .                                      +fmix(ix_lim,iy,ifld))/
     .                                    (sx(ix_lim-1,iy)*fnorm(ifld))
                   yldot(iv2) = nurlxu * (-up(ix_lim,iy,ifld)
     .                               +up(ix_lim-1,iy,ifld))/vpnorm
                   yldot(iv3) = nurlxu * (-up(ix_lim+1,iy,ifld)
     .                               +up(ix_lim,iy,ifld))/vpnorm
                elseif (islbcu.eq. 3) then
ccc islbcu=3:
                   yldot(iv2) = nurlxu * (-fmix(ix_lim,iy,ifld)
     .                                      +fmix(ix_lim-1,iy,ifld))/
     .                                    (sx(ix_lim,iy)*fnorm(ifld))
                elseif (islbcu.eq. 4) then
ccc islbcu=4:
                   yldot(iv2) = nurlxu * ( fmix(ix_lim,iy,ifld)
     .                                      -fmix(ix_lim-1,iy,ifld))/
     .                                    (sx(ix_lim,iy)*fnorm(ifld))
                elseif (islbcu.eq. 5) then
ccc islbcu=5:
                   yldot(iv2) = nurlxu * (-fmix(ix_lim,iy,ifld)
     .                +fmixy(ix_lim,iy,ifld)-fmixy(ix_lim+1,iy,ifld)
     .                                      +fmix(ix_lim+1,iy,ifld))/
     .                                    (sx(ix_lim,iy)*fnorm(ifld))
                elseif (islbcu.eq. 6) then
ccc islbcu=6:
                   yldot(iv2) = nurlxu * ( fmix(ix_lim,iy,ifld)
     .                                      -fmix(ix_lim+1,iy,ifld))/
     .                                    (sx(ix_lim,iy)*fnorm(ifld))
                endif
              endif
            enddo
c................................
c...  Now do electron temperature for iy beyond limiter
            if(isteonxy(ix_lim,iy)*isteonxy(ix_lim+1,iy)==1) then
                iv = idxte(ix_lim,iy)
                iv2 = idxte(ix_lim+1,iy)
              if (islbce .eq. 0) then
ccc islbce=0:
                yldot(iv) = nurlxe * (-te(ix_lim,iy)
     .                     +0.5*(te(ix_lim-1,iy)+te(ix_lim+2,iy)))
     .                                   *1.5*ne(ix_lim,iy)/ennorm
                yldot(iv2) = nurlxe * (-te(ix_lim+1,iy)
     .                     +0.5*(te(ix_lim-1,iy)+te(ix_lim+2,iy)))
     .                                 *1.5*ne(ix_lim+1,iy)/ennorm
              elseif (islbce .eq. 1) then
ccc islbce=1:
                yldot(iv) = nurlxe * (-te(ix_lim,iy)
     .                     +0.5*(te(ix_lim-1,iy)+te(ix_lim+2,iy)))
     .                                   *1.5*ne(ix_lim,iy)/ennorm
                yldot(iv2) = nurlxe * (-te(ix_lim+1,iy)
     .                                   +te(ix_lim,iy))
     .                                 *1.5*ne(ix_lim+1,iy)/ennorm
              elseif (islbce .eq. 2) then
ccc islbce=2:
                yldot(iv) = nurlxe * (-feex(ix_lim-1,iy)
     .                                  +feex(ix_lim+1,iy))/
     .                              (sx(ix_lim-1,iy)*vpnorm*ennorm)
                yldot(iv2) = nurlxe * (-te(ix_lim+1,iy)
     .                                   +te(ix_lim,iy))
     .                                 *1.5*ne(ix_lim+1,iy)/ennorm
              endif
            endif
c...............................
c...  Now do ion temperature for iy beyond limiter
            if(istionxy(ix_lim,iy)*istionxy(ix_lim+1,iy)==1) then
              iv = idxti(ix_lim,iy)
              iv2 = idxti(ix_lim+1,iy)
              if (islbci .eq. 0) then
ccc islbci=0:
                yldot(iv) = nurlxi * (-ti(ix_lim,iy)
     .                     +0.5*(ti(ix_lim-1,iy)+ti(ix_lim+2,iy)))
     .                                   *1.5*ne(ix_lim,iy)/ennorm
                yldot(iv2) = nurlxi * (-ti(ix_lim+1,iy)
     .                     +0.5*(ti(ix_lim-1,iy)+ti(ix_lim+2,iy)))
     .                                 *1.5*ne(ix_lim+1,iy)/ennorm
              elseif (islbci .eq. 1) then
ccc islbci=1:
                yldot(iv) = nurlxi * (-ti(ix_lim,iy)
     .                     +0.5*(ti(ix_lim-1,iy)+ti(ix_lim+2,iy)))
     .                                   *1.5*ne(ix_lim,iy)/ennorm
                yldot(iv2) = nurlxi * (-ti(ix_lim+1,iy)
     .                                   +ti(ix_lim,iy))
     .                                 *1.5*ne(ix_lim+1,iy)/ennorm
              elseif (islbci .eq. 2) then
ccc islbci=2:
                yldot(iv) = nurlxi * (-feix(ix_lim-1,iy)
     .                                  +feix(ix_lim+1,iy))/
     .                              (sx(ix_lim-1,iy)*vpnorm*ennorm)
                yldot(iv2) = nurlxi * (-ti(ix_lim+1,iy)
     .                                   +ti(ix_lim,iy))
     .                                 *1.5*ne(ix_lim+1,iy)/ennorm
              endif
              endif
              
c.....................................
c...  Do neutral gas for iy beyond limiter: => continuity conds.
            do igsp = 1, ngsp
              if(isngonxy(ix_lim,iy,igsp)*
     .                 isngonxy(ix_lim+1,iy,igsp) .eq. 1) then
                iv =  idxg(ix_lim  ,iy,igsp)
                iv2 = idxg(ix_lim+1,iy,igsp)
                if (islbcg .eq. 0) then
ccc islbcg=0:
                   yldot(iv) = nurlxg * (-ng(ix_lim,iy,igsp)
     .               +0.5*(ng(ix_lim-1,iy,igsp)+ng(ix_lim+2,iy,igsp)))/
     .                                                        n0g(igsp)
                   yldot(iv2) = nurlxg * (-ng(ix_lim+1,iy,igsp)
     .               +0.5*(ng(ix_lim-1,iy,igsp)+ng(ix_lim+2,iy,igsp)))/
     .                                                        n0g(igsp)
                elseif (islbcg .eq. 1) then
ccc islbcg=1:
                  yldot(iv) = nurlxg * (-ng(ix_lim,iy,igsp)
     .               +0.5*(ng(ix_lim-1,iy,igsp)+ng(ix_lim+2,iy,igsp)))/
     .                                                        n0g(igsp)
                  yldot(iv2) = nurlxg * (-ng(ix_lim+1,iy,igsp)
     .                                      +ng(ix_lim,iy,igsp))/
     .                                                        n0g(igsp)
                elseif (islbcg .eq. 2) then
ccc islbcg=2:
                  yldot(iv) = nurlxg * (-fngx(ix_lim-1,iy,igsp)
     .                                     +fngx(ix_lim+1,iy,igsp))/
     .                               (sx(ix_lim-1,iy)*n0g(igsp)*vpnorm)
                  yldot(iv2) = nurlxg * (-ng(ix_lim+1,iy,igsp)
     .                                      +ng(ix_lim,iy,igsp))/
     .                                                        n0g(igsp)
                endif
              endif
            enddo
 198      continue   # end of loop over iy

c..............................................
c...  Now do potential; start at iy=2 since 0,1 already set as core BCs
          do iy = max(j2p,2), min(j5p,iy_lims-1)
            if(isphionxy(ix_lim,iy)*isphionxy(ix_lim+1,iy)==1) then
              iv = idxphi(ix_lim,iy)
              iv2 = idxphi(ix_lim+1,iy)
              if (islbcp .eq. 0) then
                yldot(iv) = nurlxp * (-phi(ix_lim,iy)
     .                  +0.5*(phi(ix_lim-1,iy)+phi(ix_lim+2,iy)))
     .                                                   /temp0
                yldot(iv2) = nurlxp * (-phi(ix_lim+1,iy)
     .                  +0.5*(phi(ix_lim-1,iy)+phi(ix_lim+2,iy)))
     .                                                   /temp0
              elseif (islbcp .eq. 1) then
                yldot(iv) = nurlxp * (-phi(ix_lim,iy)
     .                  +0.5*(phi(ix_lim-1,iy)+phi(ix_lim+2,iy)))
     .                                                   /temp0
                yldot(iv2) = nurlxp * (-phi(ix_lim+1,iy)
     .                                +phi(ix_lim,iy)) / temp0
              elseif (islbcp .eq. 2) then
                yldot(iv) = nurlxp * (-fqx(ix_lim-1,iy)
     .                               +fqx(ix_lim+1,iy))/
     .                              (sx(ix_lim-1,iy)*0.1*vpnorm)
                yldot(iv2) = nurlxp * (-phi(ix_lim+1,iy)
     .                                +phi(ix_lim,iy)) / temp0
              endif
            endif
            enddo
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...  For flux tubes that do intersect physical limiter; iy>1 required
c...  Limit iy< ny; avoids non-physical flux from double-guard-cell corner
c...  Eq for iy=ny at ix=ix_lim,ix_lim+1 is at end;set Grad_y(Pg)=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do 197 iy = max(j2,iy_lims), min(j5,ny-1)  #iy loop around limiter
c...  Do the ion density on limiter
            do ifld = 1, nisp
              if(isnionxy(ix_lim,iy,ifld)*
     .                isnionxy(ix_lim+1,iy,ifld) .eq. 1) then
                iv = idxn(ix_lim,iy,ifld)
                iv2 = idxn(ix_lim+1,iy,ifld)
                if (isupgon(1).eq.1 .and. zi(ifld).eq.0.) then
                  yldot(iv) = nurlxg* (fnix(ix_lim-1,iy,ifld) +
     .                     recyllim(iy,1)*fnix(ix_lim-1,iy,1))/
     .                     (vpnorm*n0(ifld)*sx(ix_lim-1,iy))
                  yldot(iv2) =-nurlxg* (fnix(ix_lim+1,iy,ifld) +
     .                     recyrlim(iy,1)*fnix(ix_lim+1,iy,1))/
     .                     (vpnorm*n0(ifld)*sx(ix_lim+1,iy))
                else
                  yldot(iv) = nurlxn *
     .                     (ni(ix_lim-1,iy,ifld)-ni(ix_lim,iy,ifld))/
     .                                                        n0(ifld)
                  yldot(iv2) = nurlxn *
     .                     (ni(ix_lim+2,iy,ifld)-ni(ix_lim+1,iy,ifld))/
     .                                                        n0(ifld)
                endif
              endif
            enddo
              
c...  Do the parallel velocity on limiter
            do ifld = 1, nusp
              if(isuponxy(ix_lim,iy,ifld)*isuponxy(ix_lim-1,iy,ifld)*
     .                          isuponxy(ix_lim+1,iy,ifld) .eq. 1) then 
                iv3 = idxu(ix_lim,iy,ifld)
                yldot(iv3) = nurlxu*(0.-up(ix_lim,iy,ifld))/vpnorm
                iv3 = idxu(ix_lim-1,iy,ifld)
ccc   Apply sonic flow condition in "smooth" manner (MER 08 Apr 2002)
                csfac = cslim*(1.-exp(-(iy-iy_lims+1)/(cutlo+dcslim)))
                cs = csfac*sqrt( (te(ix_lim-1,iy)+
     .                            csfacti*ti(ix_lim-1,iy))/ mi(ifld) )
                if (isupgon(1).eq.1 .and. zi(ifld).eq.0) then
                  if (recycml.gt.-9.9) then
                    yldot(iv3) = nurlxu*(-recycml*cs-
     .                                   up(ix_lim-1,iy,ifld))/vpnorm
                  else
                    yldot(iv3) = nurlxu*(up(ix_lim-2,iy,ifld) -
     .                                    up(ix_lim-1,iy,ifld))/vpnorm
                  endif
                else
                  yldot(iv3) = nurlxu*(cs-up(ix_lim-1,iy,ifld))/
     .                                                         vpnorm
                endif
                iv3 = idxu(ix_lim+1,iy,ifld)
                cs = csfac*sqrt( (te(ix_lim+1,iy)+
     .                            csfacti*ti(ix_lim+1,iy))/ mi(ifld) )
                if (isupgon(1).eq.1 .and. zi(ifld).eq.0) then
                  if (recycml.gt.-9.9) then
                    yldot(iv3) = nurlxu*(recycml*cs-
     .                                  up(ix_lim+1,iy,ifld))/vpnorm
                  else
                    yldot(iv3) = nurlxu*(up(ix_lim+2,iy,ifld)-
     .                                  up(ix_lim+1,iy,ifld))/vpnorm
                  endif
                else
                  yldot(iv3) = nurlxu*(-cs-up(ix_lim+1,iy,ifld))/
     .                                                        vpnorm
                endif
              endif
            enddo   #loop over ifld for up
 
c...  Do Te on physical limter
            if(isteonxy(ix_lim,iy)*isteonxy(ix_lim+1,iy) .eq. 1) then
              iv = idxte(ix_lim,iy)
              iv2 = idxte(ix_lim+1,iy)
              totfeex = feex(ix_lim-1,iy)
              yldot(iv) = nurlxe*(totfeex/sx(ix_lim-1,iy) -bcee*
     .                           ne(ix_lim,iy)*vex(ix_lim-1,iy)*
     .                           te(ix_lim,iy))/(vpnorm*ennorm)
              totfeex = feex(ix_lim+1,iy)
              yldot(iv2) =-nurlxe*(totfeex/sx(ix_lim+1,iy) -bcee*
     .                         ne(ix_lim+1,iy)*vex(ix_lim+1,iy)*
     .                           te(ix_lim+1,iy))/(vpnorm*ennorm)
            endif
     
c... Do Ti on physical limiter
            if(istionxy(ix_lim,iy)*istionxy(ix_lim+1,iy) .eq. 1) then
              iv = idxti(ix_lim,iy)
              iv2 = idxti(ix_lim+1,iy)
              totfeix = feix(ix_lim-1,iy)  # reduced model/no neutrals
              yldot(iv) = nurlxi*(totfeix -bcei*
     .                         fnix(ix_lim-1,iy,1)*ti(ix_lim,iy))/
     .                         (vpnorm*ennorm*sx(ix_lim-1,iy))
              totfeix = feix(ix_lim+1,iy)  # reduced model/no neutrals
              yldot(iv2) = -nurlxi*(totfeix -bcei*
     .                         fnix(ix_lim+1,iy,1)*ti(ix_lim+1,iy))/
     .                         (vpnorm*ennorm*sx(ix_lim+1,iy))
            endif
   
c... Do hydrogenic ng gas on physical limiter
            do igsp = 1, nhgsp  # impurities done next
              if(isngonxy(ix_lim,iy,igsp)*isngonxy(ix_lim+1,iy,igsp)
     .                                                     .eq. 1) then
                iv =  idxg(ix_lim  ,iy,igsp)
                iv2 = idxg(ix_lim+1,iy,igsp)
                yldot(iv)  =  nurlxg*(fngx(ix_lim-1,iy,igsp) +
     .                         recyllim(iy,igsp)*fnix(ix_lim-1,iy,1) ) /
     .                              (vpnorm*n0g(igsp)*sx(ix_lim-1,iy))
                yldot(iv2) = -nurlxg*( fngx(ix_lim+1,iy,igsp) +
     .                         recyrlim(iy,igsp)*fnix(ix_lim+1,iy,1) ) /
     .                               (vpnorm*n0g(igsp)*sx(ix_lim+1,iy))
              endif
            enddo

c...############################################################            
c...  Do impurity ng gas BC touching physical limiter with sputtering
c...  First do left limiter surface at ix=ixlim, similar to right plate
c...  with positive uu into the divertor surface
         ixtl = ix_lim     #analog of ixrb+1 for right divertor plate
         ixtl1 = ix_lim-1  #analog of ixrb for right divertor plate            
         if (isimpon .ge. 4 .and. isimpon .le. 7 .and. nzspt.ge.1) then
            nzsp_rt = nhsp
            do igsp = nhgsp+1, ngsp
              jz = max(igsp - nhgsp, 1)   # identify impurity index
	      if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fnix
              if (isngonxy(ixtl,iy,igsp) .eq. 1) then
                 iv = idxg(ixtl,iy,igsp)
                 hflux = 0.
                 do ihyd = 1, nhsp
                    if (zi(ihyd).gt.0.) hflux=hflux+fnix(ixtl1,iy,ihyd)
                 enddo
                 zflux = 0.
                 do iimp = 1, nzsp(jz) # loop limited by numb. species
                    zflux = zflux + fnix(ixtl1,iy,nzsp_rt+iimp)
                 enddo
                 sputflxllim(iy,igsp) = 0.
                 zflux_chm = 0.
                 if (islim_sput(igsp) .ge. 1) then  # use fits for phys sput
                   do ifld = ipsputt_s, ipsputt_e 
                     eng_sput = ( 0.5*mi(ifld)*up(ixtl1,iy,ifld)**2 + 
     .                            ti(ixtl,iy) +  zi(ifld)*
     .                            kappallim(iy)*te(ixtl,iy) )/ev
                     if(zi(ifld)>0.) sputflxllim(iy,igsp) =
     .                                sputflxllim(iy,igsp) +
     .                                             fnix(ixtl1,iy,ifld)*
     .                      fphysyllim(igsp)*yld96(matp,matt,eng_sput)
                   enddo
                   if (islim_sput(igsp) .ge. 2) then  # add chem sput from ions
                     do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                       eng_sput = ( 0.5*mi(ifld)*up(ixtl1,iy,ifld)**2 + 
     .                            ti(ixtl,iy) +  zi(ifld)*
     .                            kappallim(iy)*te(ixtl,iy) )/ev
                       flx_incid = abs(fnix(ixtl1,iy,ifld))/sx(ixtl1,iy)
                       call sputchem (isch_sput(igsp),eng_sput,
     .                                tvliml(iy),flx_incid, yld_chm)
                       sputflxllim(iy,igsp) = sputflxllim(iy,igsp) +
     .                     fchemyllim(igsp)*fnix(ixtl1,iy,ifld)*yld_chm
                     enddo
                   endif  #test on islim_sput >= 2
c... Add chem sputtering from neutral H                     
                   if (islim_sput(igsp) .eq. 3) then #add chem sput from h neuts
                     do igsp2 = 1, 1+ishymol  #hydrogen neut fluxes only
                       t0p = max(tg(ixtl1,iy,igsp2),temin*ev)
                       flx_incid = ng(ixtl,iy,igsp2)*.25*
     .                                        sqrt(8*t0p/(pi*mg(igsp2)))
                       call sputchem (isch_sput(igsp),t0p/ev,tvliml(iy),
     .                                               flx_incid, yld_chm)
                       zflux_chm = zflux_chm - flx_incid*
     .                             fchemyllim(igsp)*yld_chm*sx(ixtl1,iy)
                     enddo
                   endif  #test on islim_sput = 3
                 endif   #test on islim_sput > 1
                 if (sputllim_use(iy,igsp) .ge. 0. .or. 
     .                           abs(sputflxllim(iy,igsp)).gt. 0.) then
                    t0 = max(cdifg(igsp)*tg(ixtl1,iy,igsp), temin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    zflux = - sputllim_use(iy,igsp) * hflux - 
     .                        sputflxllim(iy,igsp) -
     .                   recyllim(iy,igsp) * zflux -
     .                   zflux_chm + fngxsllim(iy,igsp) +
     .                                       fngxllim_use(iy,igsp)
                    yldot(iv) = -nurlxg*(fngx(ixtl1,iy,igsp) - zflux) /
     .                         (n0(igsp) * vpnorm * sx(ixtl1,iy))
                 elseif (sputllim_use(iy,igsp).lt.0. .and.
     .                   sputllim_use(iy,igsp).ge.-1.) then # sputllim_use; ==> -albedo
                    t0 = max(cdifg(igsp)*tg(ixtl1,iy,igsp), temin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    yldot(iv) = -nurlxg*( fngx(ixtl1,iy,igsp) -
     .                                   (1+sputllim_use(iy,igsp))*
     .                              ng(ixtl,iy,igsp)*vxn*sx(ixtl1,iy) )
     .                                   / (vxn*sx(ixtl1,iy)*n0g(igsp))
                 else                                # sputllim_use < -1 ==> fix dens
                    yldot(iv) = -nurlxg*(ng(ixtl,iy,igsp)-ngllim(igsp))/
     .                                                        n0g(igsp)
                 endif  # end if-test on sputllim_use
              endif    # if-test on isngon = 1
            enddo  #loop over igsp
         endif  #test on isimpon and nzspt
c...##################################################################
c...  Now do sputter to impurity gas on right side of limiter surface
c...  where -uu gives flow onto the limiter as for the left-side plate

         ixtr = ix_lim+1   #anolog of ixlb for left divertor plate
         ixtr1 = ix_lim+2  #anolog of ixlb+1 for left divertor plate
         if (isimpon .ge. 4 .and. isimpon .le. 7 .and. nzspt.ge.1) then
           nzsp_rt = nhsp
            do igsp = nhgsp+1, ngsp
              jz = max(igsp - nhgsp, 1)   # identify impurity index
	      if (jz > 1) nzsp_rt = nzsp_rt + nzsp(jz-1) #prev index for fnix
              if (isngonxy(ixtr,iy,igsp) .eq. 1) then
                 iv = idxg(ixtr,iy,igsp)
                 hflux = 0.
                 do ihyd = 1, nhsp
                    if (zi(ihyd).gt.0.) hflux = hflux+fnix(ixtr,iy,ihyd)
                 enddo
                 zflux = 0.
                 do iimp = 1, nzsp(jz) # loop limited by numb. species
                    zflux = zflux + fnix(ixtr,iy,nzsp_rt+iimp)
                 enddo
                 sputflxrlim(iy,igsp) = 0.
                 zflux_chm = 0.
                 if (islim_sput(igsp) .ge. 1) then  # use phys sput fits
                   do ifld = ipsputt_s, ipsputt_e 
                     eng_sput = ( 0.5*mi(ifld)*up(ixtr,iy,ifld)**2 + 
     .                            ti(ixtr1,iy) +  zi(ifld)*
     .                            kapparlim(iy)*te(ixtr1,iy) )/ev
                     if(zi(ifld)>0.) sputflxrlim(iy,igsp) = 
     .                                sputflxrlim(iy,igsp) +
     .                                              fnix(ixtr,iy,ifld)*
     .                      fphysyrlim(igsp)*yld96(matp,matt,eng_sput)
                   enddo
                   if (islim_sput(igsp) .ge. 2) then  # add chem sput from ions
                     do ifld = 1,1  #(ipsputt_s, ipsputt_e) place holder for imp/imp sputt
                       eng_sput = ( 0.5*mi(ifld)*up(ixtr,iy,ifld)**2 + 
     .                            ti(ixtr1,iy) +  zi(ifld)*
     .                            kapparlim(iy)*te(ixtr1,iy) )/ev
                       flx_incid = abs(fnix(ixtr,iy,ifld))/sx(ixtr,iy)
                       call sputchem (isch_sput(igsp),eng_sput,
     .                                tvliml(iy),flx_incid, yld_chm)
                       sputflxrlim(iy,igsp) = sputflxrlim(iy,igsp) +
     .                       fchemyrlim(igsp)*fnix(ixtr,iy,ifld)*yld_chm
                     enddo
                   endif   #test on islim_sput >=2
                   if (islim_sput(igsp) .eq. 3) then #add chem sput from h neuts
                     do igsp2 = 1, 1+ishymol  #hydrogen neut fluxes only
                       t0p = max(tg(ixtr1,iy,igsp2),temin*ev)
                       flx_incid = ng(ixtr,iy,igsp2)*.25*
     .                                        sqrt(8*t0p/(pi*mg(igsp2)))
                       call sputchem (isch_sput(igsp),t0p/ev,tvlimr(iy),
     .                                               flx_incid, yld_chm)
                       zflux_chm = zflux_chm - flx_incid*
     .                              fchemyrlim(igsp)*yld_chm*sx(ixtr,iy)
                     enddo
                   endif   #test on islim_sput = 3
                 endif    #test on islim_sput >= 1
                 if (sputrlim_use(iy,igsp) .ge. 0. .or. 
     .                           abs(sputflxrlim(iy,igsp)).gt. 0.) then
                    t0 = max(cdifg(igsp)*tg(ixtr1,iy,igsp), temin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    zflux = - sputrlim_use(iy,igsp) * hflux - 
     .                        sputflxrlim(iy,igsp) -
     .                   recyrlim(iy,igsp) * zflux -
     .                   zflux_chm + fngxsrlim(iy,igsp) +
     .                                       fngxrlim_use(iy,igsp)
                    yldot(iv) = nurlxg * (fngx(ixtr,iy,igsp) - zflux) /
     .                         (n0(igsp) * vpnorm * sx(ixtr,iy))
                 elseif (sputrlim_use(iy,igsp).lt.0. .and. 
     .                   sputrlim_use(iy,igsp).ge.-1) then # sputrlim_use; ==> -albedo
                    t0 = max(cdifg(igsp)*tg(ixtr1,iy,igsp), temin*ev)
                    vxn = 0.25 * sqrt( 8*t0/(pi*mg(igsp)) )
                    yldot(iv) = nurlxg*( fngx(ixtr,iy,igsp) -
     .                                   (1+sputrlim_use(iy,igsp))*
     .                               ng(ixtr1,iy,igsp)*vxn*sx(ixtr,iy) )
     .                                     / (vxn*sx(ixtr,iy)*n0g(igsp))
                 else                                # sputrlim_use < -1 ==> fix dens
                    yldot(iv) = nurlxg*(ng(ixtr,iy,igsp)-ngrlim(igsp))/
     .                                                        n0g(igsp)
                 endif  # end if-test on sputrlim_use
           
             endif  # end if-test on isngon
           enddo  # end do-loop in igsp
         endif  # if-test on isimpon and nzspt for impurity gas
cc ##########end of most limiter BC on right-side of limiter ###

 197        continue  #large loop for all iy on physical limiter
     
c ... Do potential BC for the limiter surfaces
            do iy = max(j2p,iy_lims), min(j5p,ny-1)
	       if(isphionxy(ix_lim,iy)*isphionxy(ix_lim+1,iy)==1) then  
                  iv = idxphi(ix_lim,iy)
                  iv2 = idxphi(ix_lim+1,iy)
                  yldot(iv) = -nurlxp*( phi(ix_lim,iy) - kappallim(iy)*
     .                                   te(ix_lim,iy)/ev ) / temp0  
                  yldot(iv2) =-nurlxp*(phi(ix_lim+1,iy) - kapparlim(iy)*
     .                                   te(ix_lim+1,iy)/ev ) / temp0  
               endif
            enddo 

c ... Corner guard-cell reset: Grad(Pg)=0 at iy=ny face for fngy~0 there
c ... Here normalize pressure by ev*n0g, implying Tg ~ 1 eV
	    do igsp = 1, ngsp
              if(isngonxy(ix_lim,ny+1,igsp)*
     .                          isngonxy(ix_lim+1,ny+1,igsp)==1) then
	        iv = idxg(ix_lim,ny+1,igsp)
	        iv2 = idxg(ix_lim+1,ny+1,igsp)
	        yldot(iv) = -nurlxg*(pg(ix_lim,ny+1,igsp)-
     .                         pg(ix_lim,ny,igsp))/(ev*n0g(igsp))
	        yldot(iv2) =-nurlxg*(pg(ix_lim+1,ny+1,igsp)-
     .                         pg(ix_lim+1,ny,igsp))/(ev*n0g(igsp))
              endif
            enddo
              
        endif         #test on ix2 & i5 near limiter; if-loop 8889
      endif           # end of limiter case; if-loop 8888
cc ####################################################################
cc ######################### End of limiter BCs #######################

      ncrhs = ncrhs + 1

      if (igas .eq. 1) then
         do 715 iy = j1+1-iymnbcl, j6-1+iymxbcl
           do 714 ix = i1+1-ixmnbcl, i6-1+ixmxbcl
             if(isngonxy(ix,iy,1)==1) then
               iv = idxg(ix,iy,1)
               yldot(iv) = nurlxg * (nginit(ix,iy) - ng(ix,iy,1))/n0g(1)
             endif
  714      continue
  715    continue
      endif

      return
      end
c***** end of subroutine bouncon *****
c-----------------------------------------------------------------------
      subroutine idalg

c************************************************************************
c     This subroutine sets components of the array iseqalg(neq) to be unity 
c     corresponding to the boundary and potential equation, i.e., those
c     without time derivatives. 
c************************************************************************

      implicit none
      Use(Dim)      # nx,ny,nhsp,nisp,ngsp,nusp
      Use(Math_problem_size)   # neqmx
      Use(Lsode)
      Use(UEpar)    # isnion,isupon,isteon,istion,isngon,isnionxy,isuponxy,
                    # isteonxy,istionxy,isngonxy,isphionxy 
      Use(Indexes)  # iseqalg
      Use(Compla)
      Use(Ynorm)
      Use(Selec) 
      Use(Bcond)    # isfixlb, isfixrb
      Use(Parallv)  # nxg,nyg
      Use(Xpoint_indices)  # ixpt1,ixpt2, iysptrx1,iysptrx2
      Use(Indices_domain_dcl)   #ixmxbcl,ixmnbcl,iymxbcl,iymnbcl
      Use(Share)   # islimon, ix_lim, iy_lims, geometry, nxc

      integer ifld, ii
      #Former Aux module variables
      integer ix,iy,igsp
c ... Initialize the iseqalg array to zero (==>differential equation)
      do ii = 1, neqmx
         iseqalg(ii) = 0
      enddo

c ... Now set iseqalg=1 for boundary equations and potential equation

        if (iymnbcl .eq. 1) then    # true iy=0 boundary for domain dc
         do 11 ix = 1-ixmnbcl, nx+ixmxbcl
           do ifld = 1, nisp
	    if (isnionxy(ix,0,ifld).eq.1) then
               iseqalg(idxn(ix,0,ifld)) = 1
            endif
           enddo
           do ifld = 1, nusp
	    if (isuponxy(ix,0,ifld).eq.1) then
	       iseqalg(idxu(ix,0,ifld)) = 1
	    endif
           enddo

	    if (isteonxy(ix,0).eq.1) then
               iseqalg(idxte(ix,0)) = 1
            endif
	    if (istionxy(ix,0).eq.1) then
               iseqalg(idxti(ix,0)) = 1
	    endif

           do 9 igsp = 1, ngsp
	    if (isngonxy(ix,0,igsp).eq.1) then
               iseqalg(idxg(ix,0,igsp)) = 1
	    endif 
	    if (istgonxy(ix,0,igsp).eq.1) then 
               iseqalg(idxtg(ix,0,igsp)) = 1
	    endif
   9       continue

	    if (isphionxy(ix,0).eq.1) then
               iseqalg(idxphi(ix,0)) = 1
	    endif 
   11    continue
        endif    # for iymnbcl test

        if (ixmnbcl .eq. 1) then    # true ix=0 boundary for domain dc
         do 18 iy = 1, ny
           do ifld = 1, nisp
	    if (isnionxy(0,iy,ifld).eq.1) then
               iseqalg(idxn(0,iy,ifld)) = 1
	    endif
           enddo
           do ifld = 1, nusp
	    if (isuponxy(0,iy,ifld).eq.1) then
               iseqalg(idxu(0,iy,ifld)) = 1
	    endif
           enddo

	    if (isteonxy(0,iy).eq.1) then
               iseqalg(idxte(0,iy)) = 1
	    endif
	    if (istionxy(0,iy).eq.1) then
               iseqalg(idxti(0,iy)) = 1
	    endif
             
           do 13 igsp = 1, ngsp
	    if (isngonxy(0,iy,igsp).eq.1) then
	       iseqalg(idxg(0,iy,igsp)) = 1
	    endif
	    if (istgonxy(0,iy,igsp).eq.1) then 
               iseqalg(idxtg(0,iy,igsp)) = 1
	    endif
   13      continue
	    do 15 ix = 1-ixmnbcl, nx
               if(isphionxy(ix,iy)==1) iseqalg(idxphi(ix,iy)) = 1
   15       continue
   18    continue
        endif    # for ixmnbcl test

        if (ixmxbcl .eq. 1) then    # true ix=nx+1 boundary for domain dc
         do 19 iy = 1, ny
           do ifld = 1, nisp
	    if (isnionxy(nx+1,iy,ifld).eq.1) then
               iseqalg(idxn(nx+1,iy,ifld)) = 1
	    endif
           enddo
           do ifld = 1, nusp
	    if (isuponxy(nx+1,iy,ifld)*isuponxy(nx,iy,ifld).eq.1) then
               iseqalg(idxu(nx+1,iy,ifld)) = 1
	       iseqalg(idxu(nx,iy,ifld)) = 1
	    endif
           enddo
	    if (isteonxy(nx+1,iy).eq.1) then
               iseqalg(idxte(nx+1,iy)) = 1
	    endif
	    if (istionxy(nx+1,iy).eq.1) then
               iseqalg(idxti(nx+1,iy)) = 1
	    endif
           do 17 igsp = 1, ngsp
	    if (isngonxy(nx+1,iy,igsp).eq.1) then
               iseqalg(idxg(nx+1,iy,igsp)) = 1
	    endif
	    if (istgonxy(nx+1,iy,igsp).eq.1) then 
               iseqalg(idxtg(nx+1,iy,igsp)) = 1
	    endif
   17      continue
	    if (isphionxy(nx+1,iy).eq.1) then
               iseqalg(idxphi(nx+1,iy)) = 1
	    endif
   19    continue
        endif    # for ixmxbcl test

        if (iymxbcl .eq. 1) then    # true iy=ny+1 boundary for domain dc
	 do 24 ix = 1-ixmnbcl, nx+ixmxbcl
           do ifld = 1, nisp
	    if (isnionxy(ix,ny+1,ifld).eq.1) then
               iseqalg(idxn(ix,ny+1,ifld)) = 1
	    endif
           enddo
           do ifld = 1, nusp
	    if (isuponxy(ix,ny+1,ifld).eq.1) then
               iseqalg(idxu(ix,ny+1,ifld)) = 1
	    endif
           enddo
	    if (isteonxy(ix,ny+1).eq.1) then
               iseqalg(idxte(ix,ny+1)) = 1
	    endif
	    if (istionxy(ix,ny+1).eq.1) then
               iseqalg(idxti(ix,ny+1)) = 1
	    endif
           do 22 igsp = 1, ngsp
            if (isngonxy(ix,ny+1,igsp).eq.1) then
               iseqalg(idxg(ix,ny+1,igsp)) = 1
	    endif
	    if (istgonxy(ix,ny+1,igsp).eq.1) then 
               iseqalg(idxtg(ix,ny+1,igsp)) = 1
	    endif
   22      continue
	    if (isphionxy(ix,ny+1).eq.1) then
               iseqalg(idxphi(ix,ny+1)) = 1
	    endif
   24    continue
        endif    # for iymxbcl test

c ... Check for interior boundaries with set velocities; used in rscalf
      if (isfixlb(1).eq.2 .or. isfixrb(1).eq.2) then
         if (isfixlb(1).eq.2) then
            ix = ixpt2(1)
         else
            ix = ixpt1(1)
         endif
         if (ix.ge.0) then
            do ifld = 1, nusp
               do iy = 0+1-iymnbcl, iysptrx1(1)
                  if(isuponxy(ix,iy,ifld)==1) then
                     iseqalg(idxu(ix,iy,ifld)) = 1
                  endif
               enddo
            enddo
         endif
      endif

c ... Check for double-null symmetry boundaries
      if (isudsym==1 .and. isfixlb(1) == 0) then
        do iy=1,ny   # whole ix=nxc and nxc+1 surfaces use algebraic eqns.
          do ifld = 1, nisp
            if (isnionxy(nxc,iy,ifld)*isnionxy(nxc+1,iy,ifld)==1) then
              iseqalg(idxn(nxc,iy,ifld)) = 1
              iseqalg(idxn(nxc+1,iy,ifld)) = 1
            endif
          enddo
          do ifld = 1, nusp
            if (isuponxy(nxc,iy,ifld)==1) then
              iseqalg(idxu(nxc-1,iy,ifld)) = 1
              iseqalg(idxu(nxc,iy,ifld)) = 1
              iseqalg(idxu(nxc+1,iy,ifld)) = 1
            endif
          enddo
          if (isteonxy(nxc,iy)*isteonxy(nxc+1,iy)==1) then
            iseqalg(idxte(nxc,iy)) = 1
            iseqalg(idxte(nxc+1,iy)) = 1
          endif
          if (istionxy(nxc,iy)*istionxy(nxc+1,iy)==1) then
            iseqalg(idxti(nxc,iy)) = 1
            iseqalg(idxti(nxc+1,iy)) = 1
          endif
          do igsp = 1, ngsp
            if (isngonxy(nxc,iy,igsp)*isngonxy(nxc+1,iy,igsp)==1) then
              iseqalg(idxg(nxc,iy,igsp)) = 1
              iseqalg(idxg(nxc+1,iy,igsp)) = 1
            endif
            if (istgonxy(nxc,iy,igsp)*istgonxy(nxc+1,iy,igsp)==1) then 
              iseqalg(idxtg(nxc,iy,igsp)) = 1
              iseqalg(idxtg(nxc+1,iy,igsp)) = 1
            endif
          enddo
          if (isphionxy(nxc,iy)*isphionxy(nxc+1,iy)==1) then
            iseqalg(idxphi(nxc,iy)) = 1
            iseqalg(idxphi(nxc+1,iy)) = 1
          endif
        enddo
      endif

c ... Check for an interior limiter surface
      if (islimon .ne. 0) then
        do iy = 1, ny   # whole ix=ix_lim surface uses algebraic eqns.
          do ifld = 1, nisp
            if (isnionxy(ix_lim,iy,ifld)*isnionxy(ix_lim+1,iy,ifld)==1) then
              iseqalg(idxn(ix_lim,iy,ifld)) = 1
              iseqalg(idxn(ix_lim+1,iy,ifld)) = 1
            endif
          enddo
          do ifld = 1, nusp
            if (isuponxy(ix_lim,iy,ifld)==1) then
              iseqalg(idxu(ix_lim,iy,ifld)) = 1
              if (iy .ge. iy_lims .and. isuponxy(ix_lim-1,iy,ifld)*
     .                                   isuponxy(ix_lim+1,iy,ifld)==1) then
                 iseqalg(idxu(ix_lim-1,iy,ifld)) = 1
                 iseqalg(idxu(ix_lim+1,iy,ifld)) = 1
              endif
            endif
          enddo
          if (isteonxy(ix_lim,iy)*isteonxy(ix_lim+1,iy)==1) then
            iseqalg(idxte(ix_lim,iy)) = 1
            iseqalg(idxte(ix_lim+1,iy)) = 1
          endif
          if (istionxy(ix_lim,iy)*istionxy(ix_lim+1,iy)==1) then
            iseqalg(idxti(ix_lim,iy)) = 1
            iseqalg(idxti(ix_lim+1,iy)) = 1
          endif
          do igsp = 1, ngsp
            if (isngonxy(ix_lim,iy,igsp)*isngonxy(ix_lim+1,iy,igsp)==1) then
              iseqalg(idxg(ix_lim,iy,igsp)) = 1
              iseqalg(idxg(ix_lim+1,iy,igsp)) = 1
            endif
            if (istgonxy(ix_lim,iy,igsp)*istgonxy(ix_lim+1,iy,igsp)==1) then 
              iseqalg(idxtg(ix_lim,iy,igsp)) = 1
              iseqalg(idxtg(ix_lim+1,iy,igsp)) = 1
            endif
          enddo
          if (isphionxy(ix_lim,iy)*isphionxy(ix_lim+1,iy)==1) then
            iseqalg(idxphi(ix_lim,iy)) = 1
            iseqalg(idxphi(ix_lim+1,iy)) = 1
          endif
        enddo
      endif

c ... Check for upper target plates
      if (geometry == "dnull" .or. geometry(1:9)=="snowflake" .or.
     .    geometry=="dnXtarget" .or. geometry=="isoleg") then
        do iy = 1, ny
          do ifld = 1, nisp
            if (isnionxy(ixrb(1),iy,ifld)==1) 
     .                           iseqalg(idxn(ixrb(1)+1,iy,ifld)) = 1
            if (isnionxy(ixlb(2),iy,ifld)==1) 
     .                           iseqalg(idxn(ixlb(2),iy,ifld)) = 1
          enddo
          do ifld = 1, nusp
            if (isuponxy(ixrb(1),iy,ifld)==1) 
     .                           iseqalg(idxu(ixrb(1),iy,ifld)) = 1
            if (isuponxy(ixrb(1)+1,iy,ifld)==1) 
     .                           iseqalg(idxu(ixrb(1)+1,iy,ifld)) = 1
            if (isuponxy(ixlb(2),iy,ifld)==1) 
     .                           iseqalg(idxu(ixlb(2),iy,ifld)) = 1
          enddo
          if(isteonxy(ixrb(1)+1,iy)==1) iseqalg(idxte(ixrb(1)+1,iy))=1
          if(isteonxy(ixlb(2),iy)==1)   iseqalg(idxte(ixlb(2),iy)) = 1
          if(istionxy(ixrb(1)+1,iy)==1) iseqalg(idxti(ixrb(1)+1,iy))=1
          if(istionxy(ixlb(2),iy)==1)   iseqalg(idxti(ixlb(2),iy)) = 1
          do igsp = 1, ngsp
            if (isngonxy(ixrb(1)+1,iy,igsp)==1) 
     .                            iseqalg(idxg(ixrb(1)+1,iy,igsp)) = 1
            if (isngonxy(ixlb(2),iy,igsp)==1) 
     .                            iseqalg(idxg(ixlb(2),iy,igsp)) = 1
            if (istgonxy(ixrb(1)+1,iy,igsp)==1) 
     .                            iseqalg(idxtg(ixrb(1)+1,iy,igsp)) = 1
            if (istgonxy(ixlb(2),iy,igsp)==1) 
     .                            iseqalg(idxtg(ixlb(2),iy,igsp)) = 1
          enddo
          if (isphionxy(ixrb(1)+1,iy)==1) 
     .                            iseqalg(idxphi(ixrb(1)+1,iy)) = 1
          if (isphionxy(ixlb(2),iy)==1) 
     .                            iseqalg(idxphi(ixlb(2),iy)) = 1
        enddo
      endif

      return
      end
c****** end of subroutine idalg *********
c-----------------------------------------------------------------------
      subroutine walsor

*     WALSOR defines the profiles for wall sources corresponding to gas
*     puffing or pumping; zero albedo implies a recycling wall in bouncon

      implicit none

      Use(Dim)                 # nx,ny,ngsp,nxpt
      Use(Share)               # nyomitmx
      Use(Xpoint_indices)      # ixlb,ixrb,ixpt1,ixpt2,iysptrx1,iysptrx2
      Use(Math_problem_size)   # neqmx(for arrays not used here) 
      Use(Selec)    # ixp1
      Use(Phyvar)   # pi,qe
      Use(Comgeo)   # gy,sy,xcwi,xcwo,xcpf,yylb,yyrb
      Use(Noggeo)   # vtag
      Use(Bcond)    # nwsor,igasi,igaso,issorlb,xgasi,xgaso,wgasi,wgaso,
                    # jxsori,jxsoro
                    # albdsi,albdso,matwsi,matwso,
                    # issori,iesori,issoro,iesoro,ncpli,ncplo,cplsori,cplsoro,
                    # iscpli,iscplo,fwsori,fwsoro,fngysi,fngyso,
                    # albedoo,albedoi,matwallo,matwalli,
                    # ywnii,ywnio,ywupi,ywupo,ywtei,ywteo,ywtii,ywtio,
                    # nibprof,upbprof,tebprof,tibprof
      Use(Parallv)  # nxg,nyg
      Use(Npes_mpi) # npes,mype
      Use(Rccoef)   # albedo_by_user
      Use(Indices_domain_dcg) #ndomain


*  -- local scalars --
      integer isor, jx, jxlbi, jxrbi, jxlbo, jxrbo, ixbegi, ixendi, ixbego, ixendo
      #Former Aux module variables
      integer ix,iy,igsp
      real argi, argo, xnoti, xnoto
*  -- local arrays --
      real sycosi(10), sycoso(10)

c     Calculate distances along left and right poloidal boundaries --
      do jx = 1, nxpt
         yylb(iysptrx1(jx),jx) = - 0.5/
     .        ( gy(ixlb(jx),iysptrx1(jx)) * cos(vtag(ixlb(jx),iysptrx1(jx))) )
         do iy = iysptrx1(jx)+1, ny+1
            yylb(iy,jx) = yylb(iy-1,jx) + 0.5 *
     .                       ( 1/gy(ixlb(jx),iy-1) + 1/gy(ixlb(jx),iy) ) /
     .                                       cos(vtag(ixlb(jx),iy-1))
         enddo
         do iy = iysptrx1(jx)-1, 0, -1
            yylb(iy,jx) = yylb(iy+1,jx) - 0.5 *
     .                       ( 1/gy(ixlb(jx),iy+1) + 1/gy(ixlb(jx),iy) ) /
     .                                       cos(vtag(ixlb(jx),iy))
         enddo
         yyrb(iysptrx2(jx),jx) = - 0.5/
     .    ( gy(ixrb(jx)+1,iysptrx2(jx)) * cos(vtag(ixrb(jx),iysptrx2(jx))) )
         do iy = iysptrx2(jx)+1, ny+1
            yyrb(iy,jx) = yyrb(iy-1,jx) + 0.5 *
     .                       ( 1/gy(ixrb(jx)+1,iy-1) + 1/gy(ixrb(jx)+1,iy) ) /
     .                                       cos(vtag(ixrb(jx),iy-1))
         enddo
         do iy = iysptrx2(jx)-1, 0, -1
            yyrb(iy,jx) = yyrb(iy+1,jx) - 0.5 *
     .                       ( 1/gy(ixrb(jx)+1,iy+1) + 1/gy(ixrb(jx)+1,iy) ) /
     .                                       cos(vtag(ixrb(jx),iy))
         enddo
      enddo

c ... Return if this is a parallel local domain; arrays filled from mype=0
      if (ndomain > 1 .and. nxg .ne. nx) return

c...  calculate profiles for fixed left-hand boundary conditions
      do 27 iy = iysptrx1(1)+1, ny+1
         nibprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywnio)**2 )
         upbprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywupo)**2 )
         tebprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywteo)**2 )
         tibprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywtio)**2 )
 27   continue
      do 28 iy = iysptrx1(1), 0, -1
         nibprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywnii)**2 )
         upbprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywupi)**2 )
         tebprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywtei)**2 )
         tibprof(iy) = exp( -((yylb(iy,1)-yylb0)/ywtii)**2 )
 28    continue

c...  Calculate geometry factors for boundary sources, and then boundary
c...  sources and albedos; presently, only fgnysi,o(,1) is non-zero (F-C)

      do 29 ix = 0, nx+1
         matwalli(ix) = 0
         matwallo(ix) = 0
         iscpli(ix) = 0
         iscplo(ix) = 0
         if (tvapi(ix) == 0) tvapi(ix)=1.
         if (tvapo(ix) == 0) tvapo(ix)=1.
         do 283 igsp = 1, ngsp
	    if(albedo_by_user == 0) then
              albedoi(ix,igsp) = 1.
              albedoo(ix,igsp) = 1.
            endif
            fngysi(ix,igsp) = 0.
            fngyso(ix,igsp) = 0.
 283     continue    
         do 286 isor = 1, nwsor
            fwsori(ix,isor) = 0.
            fwsoro(ix,isor) = 0.
 286     continue
 29   continue
       
c...  We allow for 10 separate sources each on the inner and outer wall
c...  If nwsor > 10, arrays must be enlarged
      if(nwsor .gt. 10) then
         call xerrab ('nwsor > 10, must increase wall source arrays')
      endif

c ... Set possible offset wall positions if parallel domain decomposed
c      xsorlbpfshift = 0.
c      xsorrbpfshift = 0.
c      xsorlbpltshift = 0.
c      xsorrbpltshift = 0.
      do 33 isor = 1, nwsor
         sycosi(isor) = 0.
         sycoso(isor) = 0.
         issori(isor) = -1
         iesori(isor) = -1
         issoro(isor) = -1
         iesoro(isor) = -1
c ... jxlb and jxrb specify the mesh regions where the wall begins and ends:
c ... private flux walls may extend over two regions for double-nulls
         jxlbi = jxsori(isor)
         if (jxlbi==1) then
            jxrbi = nxpt
         else
            jxrbi = jxlbi - 1
         endif
         ixbegi = ixlb(jxlbi)+1     # omit guard cell
         ixendi = ixrb(jxrbi)       # omit guard cell
         xnoti = issorlb(isor)*xgasi(isor)
     .             + (1-issorlb(isor))*(xcpf(ixendi)-xgasi(isor))
c ... outer walls begin and end in the same mesh region
         jxlbo = jxsoro(isor)
         jxrbo = jxlbo
         ixbego = ixlb(jxlbo)+1     # omit guard cell
         ixendo = ixrb(jxrbo)       # omit guard cell
         xnoto = issorlb(isor)*xgaso(isor)
     &             + (1-issorlb(isor))*(xcwo(ixendo)-xgaso(isor))
c...  calculate the normalization factors
         ix=ixbegi
         do
            argi = (xcpf(ix) - xnoti) * pi / (wgasi(isor)+1.e-20)
            if(ix <= ixpt1(1) .or. ix > ixpt2(1)) then
              if( abs(argi)<pi/2. ) sycosi(isor) = sycosi(isor) +
     .                                            cos(argi) * sy(ix,0)
            endif
            if (ix==ixendi) break
            ix=ixp1(ix,0)
         enddo
         ix=ixbego
         do
            argo = (xcwo(ix) - xnoto) * pi / (wgaso(isor)+1.e-20)
            if( abs(argo) .lt. pi/2. ) sycoso(isor) = sycoso(isor) +
     .                                            cos(argo) * sy(ix,ny)
            if (ix==ixendo) break
            ix=ixp1(ix,ny)
         enddo
c...  The gas is injected via fluxes fngyso,i
         ix=ixbego
         do
            argo = (xcwo(ix) - xnoto) * pi / (wgaso(isor)+1.e-20)
            if( abs(argo) .lt. pi/2. .and. sycoso(isor) .gt. 0.) then
               if(issoro(isor) .lt. 0) issoro(isor) = ix
               iesoro(isor) = ix
               if(ncplo(isor).gt.0 .or. cplsoro(isor).gt. 0.)
     .                                                    iscplo(ix) = 1
               fwsoro(ix,isor) = cos(argo)*sy(ix,ny)/(sycoso(isor)*qe)
               fngyso(ix,igspsoro(isor)) = fngyso(ix,igspsoro(isor)) + 
     .                                       igaso(isor)*fwsoro(ix,isor) +
     .                                  fvapo(isor)*sy(ix,ny)*avapo(isor)*
     .                        exp(-bvapo(isor)/tvapo(ix))/sqrt(tvapo(ix))
               if(albedo_by_user==0) 
     .                        albedoo(ix,igspsoro(isor)) = albdso(isor)
               matwallo(ix) = matwallo(ix) + matwso(isor)
            endif
            if (ix==ixendo) break
            ix=ixp1(ix,ny)
         enddo
         ix=ixbegi
         do
            argi = (xcpf(ix) - xnoti) * pi / (wgasi(isor)+1.e-20)
            if( abs(argi) .lt. pi/2. .and. sycosi(isor) .gt. 0.) then
               if(issori(isor) .lt. 0) issori(isor) = ix
               iesori(isor) = ix
               if(ncpli(isor).gt.0 .or. cplsori(isor).gt. 0.)
     .                                                   iscpli(ix) = 1
               if(ix <= ixpt1(1) .or. ix > ixpt2(1)) then
                 fwsori(ix,isor) = cos(argi)*sy(ix,0)/(sycosi(isor)*qe)
                 fngysi(ix,igspsori(isor)) = fngysi(ix,igspsori(isor)) +
     .                                     igasi(isor)*fwsori(ix,isor) +
     .                                fvapi(isor)*sy(ix,ny)*avapi(isor)*
     .                       exp(-bvapi(isor)/tvapi(ix))/sqrt(tvapi(ix))
               endif
               if(albedo_by_user==0) 
     .                     albedoi(ix,igspsori(isor)) = albdsi(isor)
               matwalli(ix) = matwalli(ix) + matwsi(isor)
            endif
            if (ix==ixendi) break
            ix=ixp1(ix,0)
         enddo
 33   continue

      return
      end
c ***** End of subroutine walsor ******
c-----------------------------------------------------------------------
      subroutine wsmodi(ig)

*     wsmodi calculates the current from source regions on the inner wall
*     boundary and resets the wall source if the coupling switch is on

      implicit none

      integer ig

      Use(Dim)      # nx,ny,ngsp
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(Selec)    # ixp1
      Use(Phyvar)   # qe,ev,pi
      Use(Compla)   # tg,mi
      Use(Comgeo)   # sy
      Use(Indexes)  # idxg
      Use(Xpoint_indices)  # ixpt1,ixpt2
###      Use(Ynorm)    # nnorm
      Use(UEpar)    # nurlxg
      Use(Comflo)   # fngy
      Use(Bcond)    # ncpli,iwalli,issori,iesori,cplsori,fngysi
      Use(Parallv)  # nxg,nyg

*  -- local scalars --
      integer isor, jsor, ixt
      real t0, vyn

      do isor = 1, nwsor
       if (igspsori(isor)==ig) then

c...  Set inner source index for coupling
        jsor = ncpli(isor)  # jsor > 0 for inner bdry, jsor < 0 for outer bdry

c...  calculate the current out at source isor
        iwalli(isor) = 0.
        if (jsor > 0) then
          do ixt = issori(isor), iesori(isor)
            if(ixt <= ixpt1(1) .or. ixt > ixpt2(1)) then
	      t0 = max(tg(ixt,0,ig),tgmin*ev)
              vyn = 0.25 * sqrt( 8*t0/(pi*mg(ig)) ) # mass only for scaling here
              iwalli(isor) = iwalli(isor) + qe*(1 - albedoi(ixt,ig))*
     .                                  ng(ixt,1,ig)*vyn*sy(ixt,0)
            endif
          enddo
	elseif (jsor < 0) then # signal for outer bdry;now change to jsor > 0
          jsor = -jsor
          do ixt = issoro(isor), iesoro(isor)
	    t0 = max(tg(ixt,ny+1,ig),tgmin*ev)
            vyn = 0.25 * sqrt( 8*t0/(pi*mg(ig)) ) # mass only for scaling here
            iwalli(isor) = iwalli(isor) + qe*(1 - albedoo(ixt,ig))*
     .                              ng(ixt,ny,ig)*vyn*sy(ixt,ny)
          enddo
        endif
c...  reinject current at coupled source jsor
        if (abs(jsor) > 0) then
          do ixt = issori(jsor), iesori(jsor)
            if(ixt <= ixpt1(1) .or. ixt > ixpt2(1)) then
	      t0 = max(tg(ixt,0,ig),tgmin*ev)
              vyn = 0.25 * sqrt( 8*t0/(pi*mg(ig)) ) # mass only for scaling here
              fngysi(ixt,ig)= iwalli(isor)*cplsori(jsor)*fwsori(ixt,jsor)
ccc     .                                   - (1 - albedoi(ixt,ig))*
ccc     .                                     ng(ixt,1,ig)*vyn*sy(ixt,0)
            endif
          enddo
        endif

       endif     # end if-test on igspori(isor)
      enddo

      return
      end
c **** End of Subroutine wsmodi *************
c-----------------------------------------------------------------------
      subroutine wsmodo(ig)

*     wsmodo calculates the current from source regions on the outer wall
*     boundary and resets the wall source if the coupling switch is on

      implicit none

      integer ig

      Use(Dim)      # nx,ny,ngsp
      Use(Math_problem_size)   # neqmx(for arrays not used here)
      Use(Selec)    # ixp1
      Use(Phyvar)   # qe
      Use(Compla)   # tg,mi
      Use(Comgeo)   # sy
      Use(Indexes)  # idxg
###      Use(Ynorm)    # nnorm
      Use(UEpar)    # nurlxg
      Use(Comflo)   # fngy
      Use(Bcond)    # ncplio,iwallio,issoro,iesoro,cplsoro,fwsoro
      Use(Parallv)  # nxg,nyg


*  -- local scalars --
      integer isor, jsor, ixt
      real t0, vyn

      do isor = 1, nwsor
       if (igspsoro(isor)==ig) then

c...  Set outer source index for coupling
        jsor = ncplo(isor)  # jsor > 0 for outer bdry, jsor < 0 for inner bdry

c...  calculate the current out at source isor
        iwallo(isor) = 0.
        if (jsor > 0) then
          do ixt = issoro(isor), iesoro(isor)
	    t0 = max(tg(ixt,ny+1,ig),tgmin*ev)
            vyn = 0.25 * sqrt( 8*t0/(pi*mg(ig)) ) # mass only for scaling here
            iwallo(isor) = iwallo(isor) + qe*(1 - albedoo(ixt,ig))*
     .                               ng(ixt,ny,ig)*vyn*sy(ixt,ny)
          enddo
	elseif (jsor < 0) then # signal for inner bdry;now change to jsor > 0 
          jsor = -jsor
          do ixt = issori(isor), iesori(isor)
	    t0 = max(tg(ixt,0,ig),tgmin*ev)
            vyn = 0.25 * sqrt( 8*t0/(pi*mg(ig)) ) # mass only for scaling here
            iwallo(isor) = iwallo(isor) + qe*(1 - albedoi(ixt,ig))*
     .                                    ng(ixt,1,ig)*vyn*sy(ixt,0)
          enddo
        endif
c...  reinject current at coupled source jsor
        if (abs(jsor) > 0) then
          do ixt = issoro(jsor), iesoro(jsor)
	    t0 = max(tg(ixt,ny+1,ig),tgmin*ev)
            vyn = 0.25 * sqrt( 8*t0/(pi*mg(ig)) ) # mass only for scaling here
            fngyso(ixt,ig)=iwallo(isor)*cplsoro(jsor)*fwsoro(ixt,jsor)
ccc     .                                 - (1 - albedoo(ixt,ig))*
ccc     .                                  ng(ixt,ny,ig)*vyn*sy(ixt,ny)  
          enddo
        endif

       endif     # end if-test on igsporo(isor)
      enddo

      return
      end
c **** End of Subroutine wsmodo *************
c-----------------------------------------------------------------------
      subroutine pltsor

*     PLTSOR defines the profiles for plate sources corresponding to gas
*     puffing or pumping; zero albedo implies a recycling wall in bouncon

      implicit none

      Use(Dim)                 # nx,ny,ngsp,nxpt
      Use(Share)               # nyomitmx,ix_lim,iy_lims
      Use(Xpoint_indices)      # ixlb,ixrb
##      Use(Math_problem_size)   # neqmx(for arrays not used here) 
##      Use(Selec)    # ixp1
      Use(Rccoef)   # igasl,rb; igspsorl,rb; ygasl,rb; wgasl,rb;fvaplb,avaplb,
                    # tvaplb,tvaprb,
      Use(Phyvar)   # pi,qe
      Use(Comgeo)   # gy,sy,xcwi,xcwo,xcpf,yylb,yyrb
      Use(Noggeo)   # vtag
      Use(Bcond)    # nwsor,igaslb,igasrb,xgaslb,xgasrb,wgaslb,wgasrb,
                    # fwsori,fsorrb,fngxslb,fngxsrb
      Use(Parallv)  # nxg,nyg

*  -- local scalars --
      integer isor,jx,igw
      #Former Aux module variables
      integer igsp,iy
      real arglb, argrb
*  -- local arrays --
      real sxcoslb(10), sxcosrb(10)
      real fsorlb(0:ny+1,npltsor), fsorrb(0:ny+1,npltsor)

c     Calculate distances along left and right boundaries --
      do jx = 1, nxpt

         yylb(iysptrx1(jx),jx) = - 0.5/
     .        ( gy(ixlb(jx),iysptrx1(jx)) * cos(vtag(ixlb(jx),iysptrx1(jx))) )
         do iy = iysptrx1(jx)+1, ny+1
            yylb(iy,jx) = yylb(iy-1,jx) + 0.5 *
     .                       ( 1/gy(ixlb(jx),iy-1) + 1/gy(ixlb(jx),iy) ) /
     .                                       cos(vtag(ixlb(jx),iy-1))
         enddo
         do iy = iysptrx1(jx)-1, 0, -1
            yylb(iy,jx) = yylb(iy+1,jx) - 0.5 *
     .                       ( 1/gy(ixlb(jx),iy+1) + 1/gy(ixlb(jx),iy) ) /
     .                                       cos(vtag(ixlb(jx),iy))
         enddo
         yyrb(iysptrx2(jx),jx) = - 0.5/
     .    ( gy(ixrb(jx)+1,iysptrx2(jx)) * cos(vtag(ixrb(jx),iysptrx2(jx))) )
         do iy = iysptrx2(jx)+1, ny+1
            yyrb(iy,jx) = yyrb(iy-1,jx) + 0.5 *
     .                       ( 1/gy(ixrb(jx)+1,iy-1) + 1/gy(ixrb(jx)+1,iy) ) /
     .                                       cos(vtag(ixrb(jx),iy-1))
         enddo
         do iy = iysptrx2(jx)-1, 0, -1
            yyrb(iy,jx) = yyrb(iy+1,jx) - 0.5 *
     .                       ( 1/gy(ixrb(jx)+1,iy+1) + 1/gy(ixrb(jx)+1,iy) ) /
     .                                       cos(vtag(ixrb(jx),iy))
         enddo
      enddo


c...  Calculate geometry factors for boundary sources, and then boundary
c...  sources

      do jx = 1, nxpt
        do iy = 0, ny+1
          if(tvaplb(iy,jx) == 0.) tvaplb(iy,jx) = 1.  #avoid divide-by-zero
          if(tvaprb(iy,jx) == 0.) tvaprb(iy,jx) = 1.  #avoid divide-by-zero
          do igsp = 1, ngsp
            fngxslb(iy,igsp,jx) = 0.
            fngxsrb(iy,igsp,jx) = 0.
          enddo  
        enddo
      enddo

c...  Allow for 10 separate sources each on the inner and outer plates
c...  If npltsor > 10, arrays must be enlarged
      if(npltsor .gt. 10) then
         call xerrab ('npltsor > 10, increase max wall source arrays')
      endif

      do jx = 1, nxpt
       do isor = 1, npltsor
         sxcoslb(isor) = 0.
         sxcosrb(isor) = 0.

c...  calculate the normalization factors
         do iy = 1, ny
            arglb = (yylb(iy,jx) - ygaslb(isor,jx))*pi / 
     .                                              (wgaslb(isor,jx)+1.e-20)
            if( abs(arglb) .lt. pi/2. ) sxcoslb(isor) = sxcoslb(isor) +
     .                                            cos(arglb) * sx(ixlb(jx),iy)
            argrb = (yyrb(iy,jx) - ygasrb(isor,jx))*pi / 
     .                                              (wgasrb(isor,jx)+1.e-20)
            if( abs(argrb) .lt. pi/2. ) sxcosrb(isor) = sxcosrb(isor) +
     .                                            cos(argrb) * sx(ixrb(jx),iy)
         enddo

c...  The gas is injected at plates via fluxes fngxslb,rb
         do iy = 1, ny
            arglb = (yylb(iy,jx) - ygaslb(isor,jx))*pi / 
     .                                      (wgaslb(isor,jx)+1.e-20)
            igw = igspsorlb(isor,jx)
            if( abs(arglb) .lt. pi/2. .and. sxcoslb(isor) .gt. 0.) then
               fsorlb(iy,isor) = cos(arglb)*sx(ixlb(jx),iy)/(sxcoslb(isor)*qe)
               fngxslb(iy,igw,jx) =  fngxslb(iy,igw,jx) + 
     .                                       igaslb(isor,jx)*fsorlb(iy,isor)
            endif
            argrb = (yyrb(iy,jx) - ygasrb(isor,jx))*pi / 
     .                                             (wgasrb(isor,jx)+1.e-20)
            igw = igspsorrb(isor,jx)
            if( abs(argrb) .lt. pi/2. .and. sxcosrb(isor) .gt. 0.) then
               fsorrb(iy,isor) = cos(argrb)*sx(ixrb(jx),iy)/(sxcosrb(isor)*qe)
               fngxsrb(iy,igw,jx) =  fngxsrb(iy,igw,jx) - 
     .                                       igasrb(isor,jx)*fsorrb(iy,isor)
            endif
         enddo
       enddo     #loop over isor for source number
      enddo      #loop over jx for X-pt index

c ... Compute and add evaporative gas flux from liquid plates
      do jx = 1, nxpt
        do igsp = 1, ngsp
         if (fvaplb(igsp,jx) + fvaprb(igsp,jx) > 1.e-20) then
           do iy = 0, ny+1
             if (tvaplb(iy,jx)<=0. .or. tvaprb(iy,jx)<=0.) then
               call remark('**ERR: tvaplb  or tvaprb = 0; must set positive')
             endif
             fngxslb(iy,igsp,jx) =  fngxslb(iy,igsp,jx) + fvaplb(igsp,jx)*
     .             sxnp(ixlb(jx),iy)*avaplb(igsp,jx)*exp(-bvaplb(igsp,jx)/
     .                                 tvaplb(iy,jx))/sqrt(tvaplb(iy,jx))
             fngxsrb(iy,igsp,jx) = fngxsrb(iy,igsp,jx) - fvaprb(igsp,jx)*
     .             sxnp(ixrb(jx),iy)*avaprb(igsp,jx)*exp(-bvaprb(igsp,jx)/
     .                                 tvaprb(iy,jx))/sqrt(tvaprb(iy,jx))
           enddo
         endif
        enddo
        enddo
        
c...  Initial limiter evaporation temp & flux arrays; single-limiter only
      do iy = 0, ny+1
        if(tvaprlim(iy) == 0) tvaprlim(iy)=1.  #prevent rate eval outside range
        if(tvapllim(iy) == 0) tvapllim(iy)=1.
        do igsp = 1, ngsp
          fngxsllim(iy,igsp) = 0.
          fngxsrlim(iy,igsp) = 0.
        enddo
      enddo      
        
c ... Compute and add evaporative gas flux from liquid-wall limiter
c ... Note below fvprlim,llim is not added recursively as no limiter
c ... gas sources unlike fngxslb,rb.

      do igsp = 1, ngsp
       if (fvapllim(igsp) + fvaprlim(igsp) > 1.e-20) then
         do iy = iy_lims, ny
           fngxsllim(iy,igsp) =  -fvapllim(igsp)*
     .           sxnp(ix_lim-1,iy)*avapllim(igsp)*exp(-bvapllim(igsp)/
     .                             tvapllim(iy))/sqrt(tvapllim(iy))
           fngxsrlim(iy,igsp) = fvaprlim(igsp)*
     .           sxnp(ix_lim+1,iy)*avaprlim(igsp)*exp(-bvaprlim(igsp)/
     .                             tvaprlim(iy))/sqrt(tvaprlim(iy))
         enddo
       endif
      enddo

      return
      end
c ***** End of subroutine pltsor ******
c------------------------------------------------------------------------
      subroutine recyprof

***   Recyprof calculates the profiles of the recycling and albedo 
***   coefficients on the inner and outer divertor plate surfaces.
***   Also computes evaporative source-fluxes fngxslb,rb from liquids.

      implicit none

      Use(Dim)      # nx,ny,nhsp,ngsp,nxpt
      Use(Comgeo)   # yylb,yyrb
      Use(Rccoef)   # recycp,recycflb,recycfrb,recylb,recyrb
                    # ndatlb,ndatrb,ydatlb,ydatrb,rdatlb,rdatrb
                    # alblb,albrb,adatlb,adatrb,sptttlb,sputtrb,
                    # albedo_by_user
      Use(Xpoint_indices)  # ixlb,ixrb
      Use(Selec)    # ixm1

*  -- local scalars --
      integer iy, ix, ifld, igsp, idat, jx

c...  Initialize plt recycling, albedo; may replace if ndatrlb,rb>0
      do jx = 1, nxpt  # loop over number of X-points
        do igsp = 1, ngsp
          do iy = 0, ny+1
            recylb(iy,igsp,jx) = recycp(igsp)*recycflb(igsp,jx) +
     .                                    recylb_use(iy,igsp,jx)
            recyrb(iy,igsp,jx) = recycp(igsp)*recycfrb(igsp,jx) + 
     .                                    recyrb_use(iy,igsp,jx)
            recycmlb(iy,igsp,jx) = recycm + recycmlb_use(iy,igsp,jx)
            recycmrb(iy,igsp,jx) = recycm + recycmrb_use(iy,igsp,jx)
          enddo
          if (abs(sputtr) .gt. 0.) then #use sputtr; neg => pump or const
            call sfill (ny+2, sputtr, sputtlb(0,igsp,jx), 1)
            call sfill (ny+2, sputtr, sputtrb(0,igsp,jx), 1)
          endif
          if(albedo_by_user==0) then
           call sfill (ny+2, albedolb(igsp,jx), alblb(0,igsp,jx), 1)
           call sfill (ny+2, albedorb(igsp,jx), albrb(0,igsp,jx), 1)
          endif
        enddo
      enddo

c...  Instead use data arrays on inner plt recycling coeff., albedo if ndatlb>0
      do jx = 1, nxpt
        do igsp = 1, ngsp
          if (ndatlb(igsp,jx).gt.50 .or. ndatrb(igsp,jx).gt.50) then
            call xerrab ('*** Exceeding storage for ndatlb,rb (>50)*')
          endif
          do idat = 1, ndatlb(igsp,jx)
            do iy = 0, ny+1
              if ( yylb(iy,jx) .gt. ydatlb(igsp,idat,jx) .and.
     .             yylb(iy,jx) .le. ydatlb(igsp,idat+1,jx) ) then
                recylb(iy,igsp,jx) = rdatlb(igsp,idat,jx) + (
     .                (rdatlb(igsp,idat+1,jx) - rdatlb(igsp,idat,jx))/
     .                (ydatlb(igsp,idat+1,jx) - ydatlb(igsp,idat,jx)) ) 
     .                    * (yylb(iy,jx) - ydatlb(igsp,idat,jx)) 
                if (adatlb(igsp,idat  ,jx) .lt. 0.99999999 .or.
     .              adatlb(igsp,idat+1,jx) .lt. 0.99999999) then
                  if (albedo_by_user == 0) then
                     alblb(iy,igsp,jx) = adatlb(igsp,idat,jx) + (
     .                 (adatlb(igsp,idat+1,jx) - adatlb(igsp,idat,jx))/
     .                 (ydatlb(igsp,idat+1,jx) - ydatlb(igsp,idat,jx)) ) 
     .                    * (yylb(iy,jx) - ydatlb(igsp,idat,jx)) 
                  endif
                endif                   
              endif
            enddo
          enddo
        enddo
      enddo  # loop over jx X-pts

c...  Instead use data arrays on outer plt recycling coeff., albedo if ndatrb>0
      do jx = 1, nxpt  # loop over x-points
        do igsp = 1, ngsp
          do idat = 1, ndatrb(igsp,jx)
            do iy = 0, ny+1
              if ( yyrb(iy,jx) .gt. ydatrb(igsp,idat,jx) .and.
     .             yyrb(iy,jx) .le. ydatrb(igsp,idat+1,jx) ) then
                recyrb(iy,igsp,jx) = rdatrb(igsp,idat,jx) + (
     .                (rdatrb(igsp,idat+1,jx) - rdatrb(igsp,idat,jx))/
     .                (ydatrb(igsp,idat+1,jx) - ydatrb(igsp,idat,jx)) ) 
     .                    * (yyrb(iy,jx) - ydatrb(igsp,idat,jx)) 
                if (adatrb(igsp,idat  ,jx) .lt. 0.99999999 .or.
     .              adatrb(igsp,idat+1,jx) .lt. 0.99999999) then
                  if (albedo_by_user == 0) then
                    albrb(iy,igsp,jx) = adatrb(igsp,idat,jx) + (
     .                (adatrb(igsp,idat+1,jx) - adatrb(igsp,idat,jx))/
     .                (ydatrb(igsp,idat+1,jx) - ydatrb(igsp,idat,jx)) ) 
     .                    * (yyrb(iy,jx) - ydatrb(igsp,idat,jx)) 
                  endif
                endif                   
              endif
            enddo
          enddo
        enddo     
      enddo  # loop over jx X-pts

#...  Set tot recycling arrays for private flux (wi) and outer wall (wo)
      do igsp = 1, ngsp
        do ix = 0, nx+1
          recycwot(ix,igsp) = recycw(igsp) + recywall_use(ix,igsp)
          do jx = 1, nxpt
            recycwit(ix,igsp,jx) = recycw(igsp) + recypf_use(ix,igsp,jx)
          enddo
        enddo
      enddo

      return
      end
c **** End of Subroutine recyprof *************
c-----------------------------------------------------------------------
      real function kappa (jsati, jsate, j)

c     Calculates the sheath drop (in units of electron temperature) with 
c     modified form that allows j > jsati and also jsati < 0.
c     Inputs are:
c        jsati = ion saturation current
c        jsate = electron saturation current
c        j     = net electrical current
c     Sign convention is current  > 0 for flow from the plasma into the plate.

      implicit none
c     input arguments
      real jsati, jsate, j
c     common blocks
      Use(Dim)
      Use(Bcond)    # kappa0, kappamx
      Use(Parallv)  # nxg,nyg

c     local variables
      real d, d0, fac, facmin

      d0 = exp(-kappa0)          # defines transition between forms
      facmin = exp(-kappamx)     # imposes an upper limit on kappa

      d = (jsati-j)/jsate
      if (d .ge. d0) then
c        standard algorithm for j < jsati - jsate*exp(-kappa0)
         fac = d + facmin
      else
c        modified form allows j > jsati - jsate*exp(-kappa0)
c        with smooth transition to standard algorithm
         fac = d0 * exp((d-d0)/d0) / 
     .              (1.+.5*((d-d0)/d0)**2+2.*((d-d0)/d0)**4) + facmin
      endif
      kappa = - log(fac)

      return
      end
c **** End of Function kappa  *************

c -----------------------------------------------------------------
      real function fluxsurfav2(inarray)
c    returns the flux surface average of 2D inarray (dimensioned nx by ny)
c       at the core boundary, averaged over
c     the core flux surface

      implicit none

c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
      Use(Dim)                  # nx,ny
      Use(Xpoint_indices)       # ixpt1, ixpt2
      Use(Indices_domain_dcl)   # iymnbcl, nx_loc,ny_loc
      Use(Npes_mpi)             # npes
      Use(Comgeo)               # vol, volv

      integer ifake  #forces Forthon scripts to put implicit none above here
CC c_mpi      include 'mpif.h'

c_mpi      integer typeind, ierr, status(MPI_STATUS_SIZE)
c_mpi      data typeind/1/
c    Assuming nx, ny are resized to be the same as nx_loc,ny_loc, per TDR.
      real inarray(0:nx,0:ny)
      real sumarr(2),gsumarr(2)
      integer ibeg,iend,i
c    sumarr(1) will be numerator, sumarr(2) will be denominator
      sumarr(1) = 0.
      sumarr(2) = 0.
      fluxsurfav2 = 0.
c    only calculate something if this domain touches the inner boundary.
c    Following coding makes use of min(ixpt1) = 0.
      if (iymnbcl == 1) then
           ibeg = ixpt1(1)+1
           iend = min(ixpt2(1),nx)
c         only calculate something if this domain is in the core, in which case ibeg <= iend
           if (iend >= ibeg) then
              do i = ibeg,iend
                  sumarr(1) = sumarr(1)+volv(i,0)*inarray(i,0)
                  sumarr(2) = sumarr(2) + vol(i,0)
              enddo
           endif
       endif
c    copy so we have a result if we are not parallel
      gsumarr(1) = sumarr(1)
      gsumarr(2) = sumarr(2)
c    Sum results from all processors if in parallel
c_mpi      if(npes > 1) then
c_mpi        call mpi_allreduce(sumarr,gsumarr,2,MPI_DOUBLE_PRECISION,MPI_SUM,
c_mpi     .     uedgeComm,ierr)
c  JRC: 11nov09: Do NOT USE MPI_COMM_WORLD
cc_mpi     .     MPI_COMM_WORLD,ierr)
c_mpi      endif
      if (gsumarr(2) == 0) return
      fluxsurfav2 = gsumarr(1)/gsumarr(2)
      return
      end

c------End of function fluxsurfav2

c----------------------------------------------------------------------------
      real function fluxsurfav1(inarray)
c    returns the flux surface average of 1D inarray which runs from global ixpt1(1)+1 to ixpt2(1)
c       at the core boundary, averaged over the core flux surface
c    For running in parallel, the array runs from ixpt1(1)+1 to iend = min(ixpt2(1),nx).
c    Coding below uses current convention that ixpt1 = 0 if x point is to left of this domain

      implicit none

c_mpi      Use(MpiVars)  #module defined in com/mpivarsmod.F.in
      Use(Dim)                  # nx,ny
      Use(Indices_domain_dcl)   # iymnbcl, nx_loc,ny_loc
      Use(Xpoint_indices)       # ixpt1, ixpt2
      Use(Npes_mpi)             # npes
      Use(Comgeo)               # vol, volv

      integer ifake  #forces Forthon scripts to put implicit none above here
CC c_mpi      include 'mpif.h'
c_mpi      integer typeind, ierr, status(MPI_STATUS_SIZE)
c_mpi      data typeind/1/
c    Assuming nx, ny are resized to be the same as nx_loc,ny_loc, per TDR
cc      real inarray(ixpt1(1)+1:ixpt2(1))
      real inarray(*)
      real sumarr(2),gsumarr(2)
      integer ibeg,iend,i
c    sumarr(1) will be numerator, sumarr(2) will be denominator
      sumarr(1) = 0.
      sumarr(2) = 0.
      fluxsurfav1 = 0.
c    only calculate something if this domain touches the inner boundary.
      if (iymnbcl == 1) then
           ibeg = ixpt1(1)+1
           iend = min(ixpt2(1),nx)
c         only calculate something if this domain is in the core, in which case ibeg <= iend
           if (iend >= ibeg) then
              do i = ibeg,iend
                  sumarr(1) = sumarr(1)+volv(i,0)*inarray(i)
                  sumarr(2) = sumarr(2) + vol(i,0)
              enddo
           endif
       endif
c    copy so we have a result if we are not parallel
      gsumarr(1) = sumarr(1)
      gsumarr(2) = sumarr(2)
c    Sum results from all processors if in parallel
c_mpi      if(npes > 1) then
c_mpi        call mpi_allreduce(sumarr,gsumarr,2,MPI_DOUBLE_PRECISION,MPI_SUM,
c_mpi     .     uedgeComm,ierr)
c  JRC: 11nov09: Do NOT USE MPI_COMM_WORLD
cc_mpi     .     MPI_COMM_WORLD,ierr)
c_mpi      endif
      if (gsumarr(2) == 0) return
      fluxsurfav1 = gsumarr(1)/gsumarr(2)
      return
      end

c----------------------------------------------------------------------c
c------------------------------------------------------------------------
      subroutine outwallflux

***   Computes particle and energy fluxes that can be passed to a
***   wall-response code like WallPSI

      implicit none

      Use(Dim)      # nx,ny,nhsp,ngsp,nxpt
      Use(Comgeo)   # sx,sy
      Use(Noggeo)   # angfx
      Use(Phyvar)   # qe,pi
      Use(UEpar)    # ebind
      Use(Comflo)   # fqx,fqy,fnix,fniy,feex,feey,feix,feiy,fngx,fngy
      Use(Xpoint_indices)  # ixlb,ixrb,ixpt1,ixpt2
      Use(Selec)    # ixm1
      Use(Outpwall) # ue_part_fluxelb, etc.
      Use(Compla)   # mg,te,ng
      Use(Poten)    # kappalb,kapparb

*  -- local scalars --
      integer ix, iy, ifld, igsp, ipf, jx, ixs, ixe, nxt, nxt1
      real sxcos, syu, ue_ion_flux_sum

c   Compute fluxes along left (inner divertor - need to change sign)
      do jx = 1, nxpt  #number of x-points; =1 for single null
        do iy = 1, ny
          sxcos = sx(ixlb(jx),iy)*cos(angfx(ixlb(jx),iy))
          ue_ion_flux_sum = 0.
          do ifld = 1, 1  #placeholder; needs extension for impurities
            ue_part_fluxh2p1lb(iy,jx) = -fnix(ixlb(jx),iy,1)/sxcos
            ue_ion_flux_sum = ue_ion_flux_sum-fnix(ixlb(jx),iy,ifld)/sxcos
          enddo
          ue_part_fluxelb(iy,jx) = ue_ion_flux_sum - fqx(ixlb(jx),iy)/
     .                                                        qe/sxcos
          ue_part_fluxh2lb(iy,jx) = 0.25*sqrt(8*tg(ixlb(jx),iy,1)/
     .                                        (pi*mg(1)))*ng(nxt1,iy,1)
          ue_heat_fluxh2p1lb(iy,jx) = -(feix(ixlb(jx),iy) +0.5*mi(1)*
     .               up(ixlb(jx),iy,1)**2*fnix(ixlb(jx),iy,1))/sxcos +
     .               ue_part_fluxelb(iy,jx)*kappal(iy,jx)*te(ixlb(jx),iy)
          ue_heat_fluxelb(iy,jx) = -feex(ixlb(jx),iy)/sxcos -
     .               ue_part_fluxelb(iy,jx)*kappal(iy,jx)*te(ixlb(jx),iy)
          ue_heat_fluxh2lb(iy,jx)=tg(ixlb(jx),iy,1)*ue_part_fluxh2lb(iy,jx)
          ue_mean_engh2p1lb(iy,jx) =  ue_heat_fluxh2p1lb(iy,jx)/
     .                                         ue_part_fluxh2p1lb(iy,jx)
          ue_mean_engelb(iy,jx) =  ue_heat_fluxelb(iy,jx)/
     .                                          ue_part_fluxelb(iy,jx)
          ue_mean_engh2lb(iy,jx) =  ue_heat_fluxh2lb(iy,jx)/
     .                                         ue_part_fluxh2lb(iy,jx)
          ue_pot_engh2p1lb(iy,jx) = ebind*ue_part_fluxh2p1lb(iy,jx)
        enddo
      enddo

c   Compute fluxes along right (outer divertor)
      do jx = 1, nxpt  #number of x-points; =1 for single null
        nxt = ixrb(jx)-1   #effec nx
        nxt1 = ixrb(jx)    #effec nx+1
        do iy = 1, ny
          sxcos = sx(nxt,iy)*cos(angfx(nxt,iy))
          ue_ion_flux_sum = 0.
          do ifld = 1, 1  #placeholder; needs extension for impurities
            ue_part_fluxh2p1rb(iy,jx) = fnix(nxt,iy,1)/sxcos
            ue_ion_flux_sum = ue_ion_flux_sum + fnix(nxt,iy,ifld)/sxcos
          enddo
          ue_part_fluxerb(iy,jx) = (ue_ion_flux_sum - fqx(nxt,iy)/qe)/sxcos
          ue_part_fluxh2rb(iy,jx) = 0.25*sqrt(8*tg(nxt1,iy,1)/(pi*mg(1)))*
     .                                                       ng(nxt1,iy,1)
          ue_heat_fluxh2p1rb(iy,jx) = (feix(nxt,iy) + 0.5*mi(1)*
     .               up(nxt,iy,1)**2*fnix(nxt,iy,1))/sxcos +
     .               ue_part_fluxerb(iy,jx)*kappar(iy,jx)*te(ixrb(jx),iy)
          ue_heat_fluxerb(iy,jx) = feex(nxt,iy)/sxcos -
     .               ue_part_fluxerb(iy,jx)*kappar(iy,jx)*te(ixrb(jx),iy)
          ue_heat_fluxh2rb(iy,jx) = tg(nxt1,iy,1)*ue_part_fluxh2rb(iy,jx)
          ue_mean_engh2p1rb(iy,jx) =  ue_heat_fluxh2p1rb(iy,jx)/
     .                                          ue_part_fluxh2p1rb(iy,jx)
          ue_mean_engerb(iy,jx) =  ue_heat_fluxerb(iy,jx)/
     .                                          ue_part_fluxerb(iy,jx)
          ue_mean_engh2rb(iy,jx) = ue_heat_fluxh2rb(iy,jx)/
     .                                          ue_part_fluxh2rb(iy,jx)
          ue_pot_engh2p1rb(iy,jx) = ebind*ue_part_fluxh2p1rb(iy,jx)
        enddo
      enddo

c   Compute fluxes along outer wall 
      do ix = 1, nx
        syu = sy(ix,ny)
        ue_ion_flux_sum = 0.
        do ifld = 1, 1  #needs extension to ifld -> nfsp for impurities
          ue_part_fluxh2p1yo(ix) = fniy(ix,ny,ifld)/syu
          ue_ion_flux_sum = ue_ion_flux_sum + fniy(ix,ny,ifld)/syu
        enddo
        ue_part_fluxeyo(ix) = ue_ion_flux_sum - (fqy(ix,ny)/qe)/syu
        ue_part_fluxh2yo(ix) = 0.25*sqrt(8*tg(ix,ny+1,1)/(pi*mg(1)))*
     .                                               ng(ix,ny+1,1)
        ue_heat_fluxh2p1yo(ix) = feiy(ix,ny)/syu +
     .                          ue_part_fluxeyo(ix)*3.0*te(ix,ny+1)
        ue_heat_fluxeyo(ix) = feey(ix,ny)/syu -
     .                          ue_part_fluxeyo(ix)*3.0*te(ix,ny+1)
        ue_heat_fluxh2yo(ix) = tg(ix,ny+1,1)*ue_part_fluxh2yo(ix)
        ue_mean_engh2p1yo(ix) =  ue_heat_fluxh2p1yo(ix)/
     .                                         ue_part_fluxh2p1yo(ix)
        ue_mean_engeyo(ix) =  ue_heat_fluxeyo(ix)/
     .                                        ue_part_fluxeyo(ix)
        ue_mean_engh2yo(ix) =  ue_heat_fluxh2yo(ix)/
     .                                         ue_part_fluxh2yo(ix)
        ue_pot_engh2p1yo(ix) = ebind*ue_part_fluxh2p1yo(ix)
      enddo

c   Compute fluxes along inner wall (need to change sign; some segms are core)
      do jx = 1, nxpt
        do ipf = 1, 2
        if (ipf == 1) then
          ixs = ixlb(jx)+1
          ixe = ixpt1(jx)
        elseif (ipf == 2) then
          ixs = ixpt2(jx)+1
          ixe = ixrb(jx)
        endif
        do ix = ixs, ixe
          syu = sy(ix,0)
          ue_ion_flux_sum = 0.  
          do ifld = 1, 1  #needs extension to ifld -> nfsp for impurities
            ue_part_fluxh2p1yi(ix) = -fniy(ix,0,ifld)/syu
            ue_ion_flux_sum = ue_ion_flux_sum + fniy(ix,0,ifld)/syu
          enddo
          ue_part_fluxeyi(ix) = -ue_ion_flux_sum + (fqy(ix,0)/qe)/syu
          ue_part_fluxh2yi(ix) = 0.25*sqrt(8*tg(ix,0,1)/(pi*mg(1)))*
     .                                                    ng(ix,0,1)
          ue_heat_fluxh2p1yi(ix) = -feiy(ix,0)/syu +
     .                          ue_part_fluxeyi(ix)*3.0*te(ix,0)
          ue_heat_fluxeyi(ix) = -feey(ix,0)/syu -
     .                          ue_part_fluxeyi(ix)*3.0*te(ix,0)
          ue_heat_fluxh2yi(ix) = tg(ix,0,1)*ue_part_fluxh2yi(ix)
          ue_mean_engh2p1yi(ix) =  ue_heat_fluxh2p1yi(ix)/
     .                                         ue_part_fluxh2p1yi(ix)
          ue_mean_engeyi(ix) =  ue_heat_fluxeyi(ix)/
     .                                        ue_part_fluxeyi(ix)
          ue_mean_engh2yi(ix) =  ue_heat_fluxh2yi(ix)/
     .                                         ue_part_fluxh2yi(ix)
          ue_pot_engh2p1yi(ix) = -ebind*ue_part_fluxh2p1yi(ix)
        enddo   # assocated with poloidal index ix
        enddo   # associated with private flux index ipf
      enddo

      return
      end

***** End of subroutine outwallflux ***********
c----------------------------------------------------------------------c
