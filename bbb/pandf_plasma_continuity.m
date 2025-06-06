c!include "bbb.h"
c!include "../com/com.h"
c!include "../mppl.h"
c!include "../sptodp.h"


      SUBROUTINE initialize_ismcnon(ylneq)
      IMPLICIT NONE
      integer ifld
      real ylneq
      Use(MCN_dim)
      Use(MCN_sources)
      Use(UEpar)
      Use(Comflo)
      Use(Comtra)
      Use(Compla)
      Use(Coefeq)
      Use(PNC_params)
      Use(Time_dep_nwt)
      Use(Ext_neutrals)
      Use(Dim)
      

c     Set switches for neutrals-related source terms in plasma equations
c     (MER 1996/10/28)
c     (IJ  2015/04/06) add ismcnon>=3 for external call to run_neutrals 
      if (ismcnon .eq. 1) then        # use MC sources only:
         cfneut=0.
         if (isupgon(1) .eq. 1) then
            cfvgpx(iigsp)=0.
            cfvgpy(iigsp)=0.
            cfvcsx(iigsp)=0.
            cfvcsy(iigsp)=0.
         endif
         cmneut=1.
      else if (ismcnon .eq. 2) then    # switch between two models:
         if (ylneq .gt. 0) then   # use fluid model for Jacobian
            cfneut=1.
            if (isupgon(1) .eq. 1) then
               cfvgpx(iigsp)=1.
               cfvgpy(iigsp)=1.
               cfvcsx(iigsp)=1.
               cfvcsy(iigsp)=1.
            endif
            cmneut=0.
         elseif (ylneq .lt. 0) then     # use MC model for evaluation
            cfneut=0.
            if (isupgon(1) .eq. 1) then
               cfvgpx(iigsp)=0.
               cfvgpy(iigsp)=0.
               cfvcsx(iigsp)=0.
               cfvcsy(iigsp)=0.
            endif
            cmneut=1.
         else
            call xerrab('*** PANDF: ismcnon=2 & yl(neq+1)=0 ???')
         endif
      else if (ismcnon .eq. 3) then         # switch between two neutral models internally
         if (ylneq .gt. 0) then         # use fluid model for preconditioner
            if (extneutmeth .eq. 1) then				#fluid source & implicit MC flux
               cfneut=1.     #turn on  fluid sources
               cfneutdiv=0.  #turn off fluid div fluxes  
               cmneut=0.     #turn off MC sources
               cmneutdiv=1.  #turn on  MC div fluxes    
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif         
            else         								#fluid source & fluid flux
               cfneut=1.     #turn on fluid sources
               cfneutdiv=1.  #turn on fluid div fluxes  
               cmneut=0.     #turn off MC sources
               cmneutdiv=0.  #turn off MC div fluxes    
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            endif 
         elseif (ylneq .lt. 0) then     # use MC model for RHS evaluation
            if (extneutmeth .eq. 1) then				#fluid source & implicit MC flux
               cfneut=1.     #turn on  fluid sources
               cfneutdiv=0.  #turn off fluid div fluxes  
               cmneut=0.     #turn off MC sources
               cmneutdiv=1.  #turn on  MC div fluxes             
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=1.
                  cfvgpy(iigsp)=1.
                  cfvcsx(iigsp)=1.
                  cfvcsy(iigsp)=1.
               endif
            else         								#MC source & fluid flux
               cfneut=0.     #turn off fluid sources
               cfneutdiv=1.  #turn on  fluid div fluxes  
               cmneut=1.     #turn on  MC sources
               cmneutdiv=0.  #turn off MC div fluxes    
               if (isupgon(1) .eq. 1) then
                  cfvgpx(iigsp)=0.
                  cfvgpy(iigsp)=0.
                  cfvcsx(iigsp)=0.
                  cfvcsy(iigsp)=0.
               endif
            endif 

         else
            call xerrab('*** PANDF: ismcnon=3 & yl(neq+1)=0 ???')
         endif
         if (ylneq .lt. 0) then       # RHS eval (not precon eval)
			if (isextneuton .ne. 0) then  # implicit use of external neutrals inside exmain
                #Neutral step
                dtold=dtreal
                dtreal=dtneut
                call store_neutrals
                call run_neutrals		  # extneutopt sets choice of model
                call update_neutrals
                dtreal=dtold
            endif
        endif
      else if (ismcnon .eq. 4) then # test a different fluid model in the preconditioner
         if (ylneq .gt. 0) then   # Precon eval
            parvis=parvis*pnc_cfparvis
            travis=travis*pnc_cftravis
            do ifld=1,nisp
              ni(:,:,ifld)=ni(:,:,ifld)*pnc_cfni(ifld)
              up(:,:,ifld)=up(:,:,ifld)*pnc_cfup(ifld)
            enddo
c            write(*,*) 'ismcnon=4'
c            write(*,*) parvis
         endif
      end if #ismcnon


      END SUBROUTINE initialize_ismcnon



      SUBROUTINE calc_volumetric_sources(xc, yc)
      IMPLICIT NONE
      Use(Selec)
      Use(Dim)
      Use(Rhsides)
      Use(Share)
      Use(Compla)
      Use(Comgeo)
      Use(Conduc)
      Use(Fixsrc)
      Use(UEpar)
      Use(Coefeq)
      Use(Comtra)
      Use(Phyvar)
      Use(Imprad)
      Use(Timing)
      Use(Xpoint_indices)
      Use(PandfTiming)
      Use(Lsode)
      Use(Bcond)
      Use(Gradients)
      Use(Cfric)
      Use(Comflo)
      Use(Jacobian_restore)
      Use(ParallelEval)
      integer xc, yc
      integer iy, ix, ifld, igsp, j2pwr, j5pwr, i2pwr, i5pwr, 
     .  ix1, ix2, jg, ifld_lcs, jz, ifld_fcs, izch,  z1fac, 
     .  iyp1, iym1, iy1
      real rdumx, dr1, dr2, rdumy, ne_sgvi, t0, t1, tsimp, nevol, ngvol, 
     .  krecz, kcxrz, kionz0, kionz, kcxrzig, niz_floor, pscx0, massfac, 
     .  kionm, krecm, kcxrm, pxri, tsnpg, t1old, t2old, t1new, t2new, 
     .  vyiy0, vyiym1, v2ix0, v2ixm1, t2, nexface, nizm_floor, tv
      real rsa, rra, rcx, tick, svdiss, sv_crumpet, tock, ave
      external rsa, rra, rcx, svdiss, sv_crumpet
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)
************************************************************************
*   We Calculate the source terms now.
************************************************************************
*  ---------------------------------------------------------------------
*  Coefficients for the source terms.
*  ---------------------------------------------------------------------

       snic = 0.0
       sniv = 0.0
       psori = 0.0
       smoc = 0.0
       smov = 0.0
       seec(ix,iy) = 0.0
       seev(ix,iy) = 0.0
       seic(ix,iy) = 0.0
       seiv(ix,iy) = 0.0
       seik(ix,iy) = 0.0
       seid(ix,iy) = 0.0
       seidh(ix,iy) = 0.0
       psorbgz(ix,iy) = 0.    # diagnostic only

************************************************************************
*  -- steady sources
************************************************************************
*  ---------------------------------------------------------------------
*  volume sources. (old PHYSRC)
*  ---------------------------------------------------------------------

c ... Calculate effective Lyman-alpha optical depth factors used if
c ... istabon=14 or 15 for hydr. rate look-ups set rtauxfac<=0 to bypass
c ----------------------------------------------------------------------
      if (rtauxfac .gt. 0.) then

CC .. FIRST GET THE POLOIDAL OPTICAL DEPTH FACTOR
c ------------------------------------------------
         do iy = 0, ny+1

c ... get optical-depth to left (ix=0) boundary
            rdumx = 0.
            do ix = 0, nx+1
               dr1  = 0.5*dx(ix,iy)
               dr2  = 0.5*dx(ix,iy)
               if(ix==(nx+1)/2) rdumx=0.  #make two plate regions independent
               rdumx = rdumx + rtauxfac*dr1*ng(ix,iy,1)
               rtaux(ix,iy) = rdumx*rt_scal
               rdumx = rdumx + rtauxfac*dr2*ng(ix,iy,1)
            enddo

c ... get optical-depth to right (ix=nx+1) bdry; initial selection of min rtaux
            rdumx = 0.
            do ix = nx+1, 0, -1
               dr1  = 0.5*dx(ix,iy)
               dr2  = 0.5*dx(ix,iy)
               if(ix==(nx+1)/2) rdumx=0.  #make two plate regions independent
               rdumx = rdumx + rtauxfac*dr1*ng(ix,iy,1)
               rtaux(ix,iy) = min(rdumx*rt_scal, rtaux(ix,iy))
               rdumx = rdumx + rtauxfac*dr2*ng(ix,iy,1)
            enddo
         enddo   # iy loop

CC .. NOW GET THE RADIAL OPTICAL DEPTH FACTOR
c --------------------------------------------
         do ix = 0, nx+1

c ... get optical-depth to inner (iy=0) bdry
            rdumy = 0.
            do iy = 0, ny+1
               dr1  = 0.5*dy(ix,iy)
               dr2  = 0.5*dy(ix,iy)
               rdumy = rdumy + rtauyfac*dr1*ng(ix,iy,1)
               rtauy(ix,iy) = rdumy*rt_scal
               rdumy = rdumy + rtauyfac*dr2*ng(ix,iy,1)
            enddo

c ... get optical-depth to outer (iy=ny+1) bdry; selection of min rtau
            rdumy = 0.
            do iy = ny+1, 0, -1
               dr1  = 0.5*dy(ix,iy)
               dr2  = 0.5*dy(ix,iy)
               rdumy = rdumy + rtauyfac*dr1*ng(ix,iy,1)
               rtauy(ix,iy) = min(rdumy*rt_scal, rtauy(ix,iy))
               rtau(ix,iy) = min(rtaux(ix,iy), rtauy(ix,iy))
               rdumy = rdumy + rtauyfac*dr2*ng(ix,iy,1)
            enddo
         enddo   # ix loop

      endif     # test on rtauxfac, skip if rtauxfac is negative

*     The following is a temporary recycling model.

*  -- recalculate particle source psor if ifixpsor=0 --
c...  The particle source can be frozen if ifixpsor.ne.0
      if(ifixpsor .eq. 0) then
            
        igsp = 0
        do ifld = 1, nhsp  # Hydrogen-only loop
         if (zi(ifld) > 0.) then  #calc only for hydrogen ions
          igsp = igsp + 1
          do iy = iys1, iyf6
            do ix = ixs1, ixf6

c     Ionization of neutral hydrogen by electrons and recombination--
               if (icnuiz .eq. 0) then
                  ne_sgvi = ne(ix,iy)
                  if (ifxnsgi.eq.1) ne_sgvi = cne_sgvi  # fix density dependence
                  nuiz(ix,iy,igsp) = chioniz *  ne(ix,iy) * (
     .                           rsa(te(ix,iy),ne_sgvi,rtau(ix,iy),0)
     .                         + sigvi_floor )
               elseif (icnuiz .eq. 1) then
                  nuiz(ix,iy,igsp) = cnuiz
               endif
               if (isrecmon == 1) then
                  nurc(ix,iy,igsp) = cfrecom * ne(ix,iy) 
     .                         * rra(te(ix,iy),ne(ix,iy),rtau(ix,iy),1)
               else
                   nurc(ix,iy,igsp) = 0.
               endif
               psorbgg(ix,iy,igsp) = ngbackg(igsp)*( (0.9 + 0.1*
     .                            (ngbackg(igsp)/ng(ix,iy,igsp))**ingb) ) * 
     .                             nuiz(ix,iy,igsp) * vol(ix,iy)
               psorgc(ix,iy,igsp) = -ng(ix,iy,igsp)*nuiz(ix,iy,igsp)*vol(ix,iy) +
     .                              psorbgg(ix,iy,igsp)
               psorc(ix,iy,ifld) = - psorgc(ix,iy,igsp)
               psordis(ix,iy,2) = cfdiss*psorc(ix,iy,1)  # changed below if ishymol=1
               psorxrc(ix,iy,ifld) = -ni(ix,iy,ifld)*nurc(ix,iy,igsp)*vol(ix,iy)
               psorrgc(ix,iy,igsp) = -psorxrc(ix,iy,ifld)
               msor(ix,iy,ifld) = 0.
               msorxr(ix,iy,ifld) = 0.


c     Charge exchange on neutral hydrogen --
              if (icnucx .eq. 0) then
	         t0 = max(ti(ix,iy),temin*ev)
ccc   we omit the weak velocity dependence as it brings in ni(ix+1) in Jac
                 t1 = t0/(mi(ifld)/mp)
                 nucx(ix,iy,igsp) = ni(ix,iy,ifld) * rcx(t1,ni(ix,iy,ifld),1)
              elseif (icnucx .eq. 1) then
                 nucx(ix,iy,igsp) = cnucx
              elseif (icnucx == 2) then
	         t0 = max(ti(ix,iy),temin*ev)
                 nucx(ix,iy,igsp) = sqrt(t0/mi(ifld))*
     .                         sigcx*(ni(ix,iy,ifld)+rnn2cx*ng(ix,iy,igsp))
              endif
              nuix(ix,iy,igsp) = fnuizx*nuiz(ix,iy,igsp) + 
     .                           fnucxx*nucx(ix,iy,igsp)
                              #dont use neutral-neutral collisions here
c
c   neutral particle source/sink for isupgon=1 (reset below if multispecies
c   models are on [isimpon = 5 or 6 or 7])
              if(isupgon(igsp) .eq. 1)then # inertia gas species is ifld+1
                 psorc(ix,iy,ifld+1)= -psorc(ix,iy,ifld)
                 psorxrc(ix,iy,ifld+1)= -psorxrc(ix,iy,ifld)
                 msor(ix,iy,ifld+1)= 0.
                 msorxr(ix,iy,ifld+1)= 0.
              endif
c
            enddo   #end loop over ix
          enddo     #end loop over iy
         endif      #omit whole loop if zi(ifld) = 0. (neutrals)
        enddo       #end loop over hydrogen species (ifld)

c*****************************************************************
c ... Average psorgc and psorc over cell vol with simple 5pt ave
c*****************************************************************
        igsp = 0
        do ifld = 1, nhsp  # Hydrogen-only loop
         if (zi(ifld) > 0.) then  #calc only for hydr ions, not neuts
           igsp = igsp + 1
           if (ispsorave.eq.0.) then  #use only single-cell value
             do iy = iys1, iyf6
               do ix = ixs1, ixf6
                 psorg(ix,iy,igsp) = psorgc(ix,iy,igsp)
                 psor(ix,iy,ifld) =  psorc(ix,iy,ifld)
                 psorxr(ix,iy,ifld) = psorxrc(ix,iy,ifld)
                 psorrg(ix,iy,igsp) = psorrgc(ix,iy,igsp)
                 if(isupgon(igsp) .eq. 1) then # inertia gas is ifld+1
                   psor(ix,iy,ifld+1)= -psor(ix,iy,ifld)
                   psorxr(ix,iy,ifld+1)= -psorxr(ix,iy,ifld)
                 endif
               enddo
             enddo

           elseif (ispsorave > 0.) # use 5pt ave; first divide by vol

             if (xc < 0) then  #full RHS eval
               j2pwr = j2
               j5pwr = j5
             else  # Jacobian eval
               j2pwr = max(1, yc-1)
               j5pwr = min(ny, yc+1)
             endif 
             if (ParallelPandfCall.gt.0) then
                j2pwr = j2omp
                j5pwr = j5omp
             endif
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
               endif
               do ix = i2pwr, i5pwr
                 ix1 = ixm1(ix,iy)
                 ix2 = ixp1(ix,iy)
                 psorg(ix,iy,igsp) = (1.-ispsorave*0.5)*
     .                                  psorgc(ix,iy,igsp)+ 
     .                               0.125*ispsorave*vol(ix,iy)*
     .                          ( psorgc(ix,iy-1,igsp)/vol(ix,iy-1) + 
     .                            psorgc(ix,iy+1,igsp)/vol(ix,iy+1) +
     .                            psorgc(ix1,iy,igsp)/vol(ix1,iy)   + 
     .                            psorgc(ix2,iy,igsp)/vol(ix2,iy) )
                 psorxr(ix,iy,ifld) = (1.-ispsorave*0.5)*
     .                                 psorxrc(ix,iy,ifld) + 
     .                                  0.125*ispsorave*vol(ix,iy)*
     .                           ( psorxrc(ix,iy-1,ifld)/vol(ix,iy-1) + 
     .                             psorxrc(ix,iy+1,ifld)/vol(ix,iy+1) +
     .                             psorxrc(ix1,iy,ifld)/vol(ix1,iy)   + 
     .                             psorxrc(ix2,iy,ifld)/vol(ix2,iy) )
                 psor(ix,iy,ifld) = -psorg(ix,iy,igsp)
                 psorrg(ix,iy,igsp) = -psorxr(ix,iy,ifld)
                 if(isupgon(igsp) .eq. 1) then # inertia gas is ifld+1
                   psor(ix,iy,ifld+1)= -psor(ix,iy,ifld)
                   psorxr(ix,iy,ifld+1)= -psorxr(ix,iy,ifld)
                 endif

               enddo   #end loop over ix
             enddo     #end loop over iy
           endif       #if-loop on ipsorave
         endif         #omit whole loop if zi(ifld) = 0. (neutrals)
        enddo          #end loop over hydrogen species (ifld)

c ... Can now calc current from nucx since it is updated
      if (cfqyn .gt. 0.) call calc_curr_cx

c ... Ionization and recombination of impurities.
c     (Note that charge-exchange recombination is implemented for
c     impurities but the corresponding terms need to be added to
c     equations for hydrogenic-ion and gas continuity and electron
c     and ion energy.) ???
c     The total source is written as psor+psorxr, where
c     psor = n_(z-1) ne K^i_(z-1) - n_z ne K^i_z     # ionization gain/loss
c     psorxr = -n_z[ne K^r_z+ng K^cx_z]              # cx/r loss from z to z-1
c              +n_(z+1)[ne K^r_(z+1)+ng K^cx_(z+1)]  # cx/r gain to z from z+1

        if (isimpon .ge. 5 .and. nzspt .gt. 0) then
          do iy = iys1, iyf6
             do ix = ixs1, ixf6

                  if (istimingon .eq. 1) tsimp = tick()
                  nevol = ne(ix,iy) * vol(ix,iy)
                  ngvol = ng(ix,iy,1) * vol(ix,iy)

                  jg = nhgsp
                  ifld_lcs = nhsp
                  do jz = 1, ngspmx-1           # for all impurity isotopes
                     if (nzsp(jz)==0) break
                     ifld_fcs = ifld_lcs + 1
                     ifld_lcs = ifld_fcs + nzsp(jz) - 1
                     if (ngsp .gt. nhgsp) then  # impurity gas is present
                         jg = jg + 1
                         if (ismctab .eq. 1) then  # rates for Z=0 gas
                            call imprates(te(ix,iy), 0, nzsp(jz), kionz0,
     .                                    krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy),te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   0, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz0, krecz, kcxrz)
                         endif
                         kionz0 = kionz0 + sigvi_floor
			 psorbgg(ix,iy,jg)= ngbackg(jg)*
     .                     (0.9+0.1*(ngbackg(jg)/ng(ix,iy,jg))**ingb) * 
     .                                                      nevol*kionz0
                         psorg(ix,iy,jg) = -ng(ix,iy,jg)*nevol*kionz0 +
     .                                      psorbgg(ix,iy,jg)
                         psor(ix,iy,ifld_fcs) = - psorg(ix,iy,jg)
                         msor(ix,iy,ifld_fcs)= 0.  # zero gas mom. assumed
                         if (ismctab .eq. 1) then  # rates for Z=1 ions
                            call imprates(te(ix,iy), 1, nzsp(jz), kionz,
     .                                    krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy), te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   1, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz, krecz, kcxrz)
                            kcxrz = cfcximp1*kcxrz  # rescale if desired
                         endif
                         kionz = kionz + sigvi_floor # only to set kionm below
                         kcxrzig = rcxighg(jg)*kcxrz  # K_cx of ng(jg)+ni(1)->
                         niz_floor = nzbackg(ifld_fcs) * (0.9 + 0.1*
     .                          (nzbackg(ifld_fcs)/ni(ix,iy,ifld_fcs))**inzb)
                         pscx0 = ngvol*(ni(ix,iy,ifld_fcs)-niz_floor)*kcxrz - 
     .                           ng(ix,iy,jg)*ni(ix,iy,1)*vol(ix,iy)*
     .                                                        kcxrzig
                         psorcxg(ix,iy,jg) = pscx0
                         psorcxg(ix,iy,1) = -pscx0
                         psorrg(ix,iy,jg) = nevol*(ni(ix,iy,ifld_fcs)-
     .                                                   niz_floor)*krecz
                         psorxr(ix,iy,ifld_fcs)= -psorrg(ix,iy,jg) - pscx0 
                         psorxr(ix,iy,1) = psorxr(ix,iy,1) + pscx0
cc                    Note: summed over ion/neutrals here backgrd source=0
                         massfac = cfmassfac*16*mi(1)/(3*(mg(jg)+mi(1)))
                         nuiz(ix,iy,jg) = kionz0*ne(ix,iy)
                         nuix(ix,iy,jg) = fnuizx*nuiz(ix,iy,jg) + 
     .                                   kcxrzig*ni(ix,iy,1) +
     .                         massfac*( kelighi(jg)*ni(ix,iy,1) +
     .                                   kelighg(jg)*ng(ix,iy,1) )
                         nucxi(ix,iy,ifld_fcs) = sigcxms(ifld_fcs,jg)*
     .                          sqrt(ti(ix,iy)/mi(ifld_fcs))*ng(ix,iy,jg)
                         nucx(ix,iy,jg) = sigcxms(ifld_fcs,jg)*
     .                        sqrt(ti(ix,iy)/mi(ifld_fcs))*ni(ix,iy,ifld_fcs)
                         massfac = cfmassfac*16*mi(ifld_fcs)/
     .                                               (3*(mg(jg)+mi(ifld_fcs)))
                         nueli(ix,iy,ifld_fcs) = massfac*( keligii(jg)*
     .                                                    ng(ix,iy,jg) )
                         nuelg(ix,iy,jg) = massfac*( keligii(jg)*
     .                                                    ni(ix,iy,ifld_fcs) )

                     else                       # no impurity gas present
                         izch = nint(zi(ifld_fcs))
                         if (ismctab .eq. 1) then
                            call imprates(te(ix,iy), izch, nzsp(jz), 
     .                                   kionz, krecz, kcxrz)
                         elseif (ismctab .eq. 2) then
                            call mcrates(ne(ix,iy), te(ix,iy),
     .                                   ti(ix,iy)*mp/mi(1),
     .                                   izch, nzsp(jz), znucl(ifld_fcs),
     .                                   kionz, krecz, kcxrz)
                         endif
                         kionz = kionz + sigvi_floor
                         psor(ix,iy,ifld_fcs) = 0.   # in case nzspt=1
                         psorxr(ix,iy,ifld_fcs) = 0.
                         msor(ix,iy,ifld_fcs) = 0.   # in case nzspt=1
                         msorxr(ix,iy,ifld_fcs) = 0.
                         krecz = 0.
                         kcxrz = 0.
                     endif   # end if-branches for Z=1 with/wo impurity gas
                     if (znucl(ifld_fcs)==1) then
                         # hydrogenic impurity: include cx on ifld_fcs in nuix
                         nuix(ix,iy,jg) = nuix(ix,iy,jg) + 
     .                                   kcxrz*ni(ix,iy,ifld_fcs)
                         nuix(ix,iy,1) = nuix(ix,iy,1) + 
     .                                   kcxrzig*ni(ix,iy,ifld_fcs)
                     endif

                     kionm = kionz     # set values as previous charge state
                     krecm = krecz
                     kcxrm = kcxrz

                     do ifld = ifld_fcs + 1, ifld_lcs  # for charge states Z > 1
                        izch = nint(zi(ifld))
                        if (ismctab .eq. 1) then
                           call imprates(te(ix,iy), izch, nzsp(jz), 
     .                                  kionz, krecz, kcxrz)
                        elseif (ismctab .eq. 2) then
                           call mcrates(ne(ix,iy), te(ix,iy),
     .                                  ti(ix,iy)*mp/mi(1),
     .                                  izch, nzsp(jz), znucl(ifld),
     .                                  kionz, krecz, kcxrz)
                            kcxrz = cfcximp2*kcxrz   #rescale if desired
                        endif
                        kionz = kionz + sigvi_floor
			if (ifld==ifld_lcs) kionz = 0. #ensure no lcs ioniz
                        pxri = 0.    # gets reset if ifld.eq.ifld_fcs+1
                        z1fac = 1.   # gets reset = 0 if ifld.eq.ifld_fcs+1

                        if (ifld .eq. ifld_fcs+1) then #for 2nd charge-state
                           nizm_floor = nzbackg(ifld_fcs) * (0.9 + 0.1*
     .                         (nzbackg(ifld_fcs)/ni(ix,iy,ifld_fcs))**inzb)
                           psor(ix,iy,ifld_fcs) = psor(ix,iy,ifld_fcs)-
     .                            nevol*(ni(ix,iy,ifld_fcs)-nizm_floor)*
     .                                                            kionm
			   psorbgz(ix,iy) = psorbgz(ix,iy)-nevol*nizm_floor*
     .                                                            kionm
	                   msor(ix,iy,ifld_fcs) = msor(ix,iy,ifld_fcs)-
     .                            nevol*(ni(ix,iy,ifld_fcs))*
     .                            kionm*mi(ifld_fcs)*up(ix,iy,ifld_fcs)    
                           pxri = psorxr(ix,iy,ifld_fcs) #set in Z=1 loop
                           z1fac = 0.
                        endif

                        niz_floor = nzbackg(ifld) * (0.9 + 0.1*
     .                            (nzbackg(ifld)/ni(ix,iy,ifld))**inzb)
                        psor(ix,iy,ifld) = nevol *
     .                                     ( ni(ix,iy,ifld-1) * kionm -
     .                             (ni(ix,iy,ifld)-niz_floor) * kionz )
			psorbgz(ix,iy) = psorbgz(ix,iy)-nevol*niz_floor*kionz
                        msor(ix,iy,ifld) = nevol *
     .                                  ( ni(ix,iy,ifld-1) * kionm *
     .                                   mi(ifld-1) * up(ix,iy,ifld-1) -
     .                            (ni(ix,iy,ifld)) * kionz *
     .                                     mi(ifld) * up(ix,iy,ifld) )
                        psorxr(ix,iy,ifld-1) = pxri-z1fac*(nevol*krecm +
     .                                                     ngvol*kcxrm) *
     .                                    (ni(ix,iy,ifld-1)-nizm_floor) +
     .                                   (nevol*krecz + ngvol*kcxrz) *
     .                                                   ni(ix,iy,ifld)
			psorbgz(ix,iy) = psorbgz(ix,iy) + z1fac*
     .                                   (nevol*krecm + ngvol*kcxrm) *
     .                                   nizm_floor
                        msorxr(ix,iy,ifld-1) = 0. - (nevol*krecm +
     .                                               ngvol*kcxrm) *
     .                                    (ni(ix,iy,ifld-1)) *
     .                                      mi(ifld-1)*up(ix,iy,ifld-1) +
     .                        (nevol*krecz + ngvol*kcxrz)*ni(ix,iy,ifld)*
     .                                            mi(ifld)*up(ix,iy,ifld)
                        psorxr(ix,iy,1) = psorxr(ix,iy,1) + ngvol*
     .                                               ni(ix,iy,ifld)*kcxrz
                        psorcxg(ix,iy,1) = psorcxg(ix,iy,1) - ngvol*
     .                                               ni(ix,iy,ifld)*kcxrz
                        nucxi(ix,iy,ifld) = sigcxms(ifld,jg)*
     .                              sqrt(ti(ix,iy)/mi(ifld))*ng(ix,iy,jg)
                        nucx(ix,iy,jg) = nucx(ix,iy,jg) + sigcxms(ifld,jg)*
     .                          sqrt(ti(ix,iy)/mi(ifld))*ni(ix,iy,ifld)
                        massfac = cfmassfac*16*mi(ifld)/
     .                                              (3*(mg(jg)+mi(ifld)))
                        nueli(ix,iy,ifld) = massfac*( keligii(jg)*
     .                                                   ng(ix,iy,jg) )
                        nuelg(ix,iy,jg) = nuelg(ix,iy,jg) + massfac*
     .                                  ( keligii(jg)*ni(ix,iy,ifld) )

                        kionm = kionz
                        krecm = krecz
                        kcxrm = kcxrz
                        nizm_floor = niz_floor
                        if (ifld .eq. ifld_lcs) then  # last charge-state
                          psorxr(ix,iy,ifld) = -(nevol * krecz +
     .                                           ngvol * kcxrz) *
     .                                          (ni(ix,iy,ifld)-niz_floor) 
			  psorbgz(ix,iy) = psorbgz(ix,iy) + niz_floor *
     .                                      (nevol*krecz + ngvol*kcxrz)
                          msorxr(ix,iy,ifld) = -(nevol * krecz +
     .                                           ngvol * kcxrz) *
     .                                          (ni(ix,iy,ifld))*
     .                                             mi(ifld)*up(ix,iy,ifld)
                          nuix(ix,iy,jg) = nuix(ix,iy,jg) + nucx(ix,iy,jg) +
     .                                     nuelg(ix,iy,jg)
                        endif


                     enddo   # end do-loop on charge states for isotope jz
                  enddo  # end do-loop on impurity isotopes

                  if (istimingon .eq. 1) call timimpfj (tsimp, xc)
c
c   neutral particle source/sink for isupgon=1; make consistent with impurity
c   contributions just calculated for multispecies
              if (isupgon(1) .eq. 1) then #should be generalized to D & T
                 psor(ix,iy,iigsp)= -psor(ix,iy,1)
                 psorxr(ix,iy,iigsp)= -psorxr(ix,iy,1)
              endif
c
            end do
           end do
        endif            # end of if (isimpon .ge. 5 .and. nzspt .gt. 0)

c ... Add volume particle loss terms for quasi 1-D radial model
      if (l_parloss .le. 1e9) then
        do iy = iys1, iyf6  # core region has no loss
          do ix = ixs1, ixf6
            do ifld = 1, nfsp
              if (iy .le. iysptrx) then # inside the LCFS when nxpt>1
                                        # (see definition of iysptrx in nphygeo)
                nuvl(ix,iy,ifld) = 0.
              else
                nuvl(ix,iy,ifld) = ( cfvlh*
     .                             sqrt((te(ix,iy)+ti(ix,iy))/mi(1))+
     .                             cfvli(ifld)*
     .              sqrt((zi(ifld)*te(ix,iy)+ti(ix,iy))/mi(ifld)) ) /
     .                                                      l_parloss
              endif
            enddo
          enddo
        enddo
      endif

c ... Set up nuiz & sources for hydrogen molecular gas
      if (ishymol .eq. 1) then
         if (nhgsp .eq. 1) then
          call xerrab('*** nhgsp must exceed 1 for ishymol=1 ***')
        endif
        do iy = iys1, iyf6
         do ix = ixs1, ixf6
           nuiz(ix,iy,2) = 
     .          ne(ix,iy) * (  
     .              (1-ishymol*ismolcrm)*(svdiss( te(ix,iy) )
     .              + cfizmol*rsa(te(ix,iy),ne_sgvi,rtau(ix,iy),0)
     .              + sigvi_floor ) 
     .          ) - ishymol*ismolcrm*sv_crumpet( te(ix,iy), ne(ix,iy),10)
           massfac = 16*mi(1)/(3*(mg(2)+mi(1)))
           nuix(ix,iy,2)= fnuizx*nuiz(ix,iy,2) + 
     .                           massfac*( kelighi(2)*ni(ix,iy,1)+
     .                                     kelighg(2)*ng(ix,iy,1) )
c ...  molecule-molecule collisions would enter viscosity, not nuix
           psorbgg(ix,iy,2) = ngbackg(2)* 
     .                     (0.9+0.1*(ngbackg(2)/ng(ix,iy,2))**ingb ) * 
     .                                        nuiz(ix,iy,2) * vol(ix,iy)
           psorgc(ix,iy,2) = - ng(ix,iy,2)*nuiz(ix,iy,2)*vol(ix,iy) +
     .                        psorbgg(ix,iy,2)
                psorg(ix,iy,2) = psorgc(ix,iy,2)  # no mol sor averaging
                psordisg(ix,iy,2) = -ng(ix,iy,2)*nuiz(ix,iy,2)*vol(ix,iy)
                psordis(ix,iy,2) = 
     .                  - (1-ishymol*ismolcrm)*2*psordisg(ix,iy,2) 
     .                  + ishymol*ismolcrm*cfcrma*ng(ix,iy,2)*vol(ix,iy)
     .                  * sv_crumpet(te(ix,iy),ne(ix,iy),11)
                # 2 atoms per molecule in old model, rates from CRM for new
                psordisg(ix,iy,1)=psordis(ix,iy,2)
                psordis(ix,iy,1) = -cfcrmi*(2*psordisg(ix,iy,2)+
     .                              psordis(ix,iy,2))
                psor(ix,iy,1) = psor(ix,iy,1) + psordis(ix,iy,1)
c ... TODO: How to deal with diffusive atom model - is it maintained?
                if(isupgon(1) .eq. 1) then
                  psor(ix,iy,iigsp) = psor(ix,iy,iigsp) + psordis(ix,iy,iigsp)
                endif
         enddo
        enddo 
      endif  # end of loop for ishymol=1 (hydrogen molecules on)


c  *** Now integrate sources over cell volume if ishosor=1 & yl(neq+1)=-1,
c  *** where the last condition means this is only a full RHS eval, not
c  *** a Jacobian calculation

         if (ishosor.eq.1) then  #full RHS eval

          if (yl(neq+1).lt.0) then  #full RHS eval

c ...    integ. sources over cells (but not for Jac) for higher-order accuracy

             do ifld = 1, nfsp  # loop over ions
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                 ixm1(0:nx+1,0:ny+1),fsprd, vol(0:nx+1,0:ny+1), 
     .                 psor_tmpov(0:nx+1,0:ny+1), psor(0:nx+1,0:ny+1,ifld))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                 ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                 psor_tmpov(0:nx+1,0:ny+1), psorxr(0:nx+1,0:ny+1,ifld))
             enddo

c *** Now do the gas
             do igsp = 1, ngsp  # now loop over gas species
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                psor_tmpov(0:nx+1,0:ny+1), psorg(0:nx+1,0:ny+1,igsp))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                psor_tmpov(0:nx+1,0:ny+1), psorrg(0:nx+1,0:ny+1,igsp))
                call volave(nx, ny, j2, j5, i2, i5, ixp1(0:nx+1,0:ny+1), 
     .                ixm1(0:nx+1,0:ny+1), fsprd, vol(0:nx+1,0:ny+1), 
     .                psor_tmpov(0:nx+1,0:ny+1), psorcxg(0:nx+1,0:ny+1,igsp))
             enddo
       
          endif   # end of if (yl(neq+1).lt.0) test
         endif    # end of integrating over sources and ishosor test
      
      endif              # end of big loop starting if (ifixpsor .eq. 0) 

*-----------------------------------------------------------------------
*  -- Calculates the fixed source if it is on
*-----------------------------------------------------------------------

      if (ifixsrc .ne. 0) then
         do iy = j2omp, j5omp
            do ix = i2omp, i5omp
               snic(ix,iy,1) = snic(ix,iy,1) + vol(ix,iy) * a1n *
     .                          exp(-b1n*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1n*(yyc(iy)-yysrc)**2)
               seic(ix,iy) = seic(ix,iy) + vol(ix,iy) * a1i * ev *
     .                          exp(-b1i*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1i*(yyc(iy)-yysrc)**2)
               seec(ix,iy) = seec(ix,iy) + vol(ix,iy) * a1e * ev *
     .                          exp(-b1e*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1e*(yyc(iy)-yysrc)**2)
            end do
        end do
      endif



*****************************************************************
*  Other volume sources calculated in old SRCMOD
*****************************************************************
*  ---------------------------------------------------------------------
*  electron-ion transfer terms and an
*  approximation to the ion-ion thermal force.
*  cfw: For the neutral momentum eq there is also a v_gas*grad(p_gas)
*       term which is evaluated using ng, ti and gpiy
*  ---------------------------------------------------------------------

c...  Force fluxes and gradients on cuts to be zero for half-space problems
      if (isfixlb(1).eq.2.or. isfixrb(1).eq.2) then
         if (isfixlb(1).eq.2) then
            ix = ixpt2(1)
         else
            ix = ixpt1(1)
         endif
         if (ix.ge.i2 .and. ix.le.i5+1 .and. iysptrx1(1) > 0) then
            do iy = 0, iysptrx1(1)
               gpex(ix,iy) = 0.
               frice(ix,iy) = 0.
               ex(ix,iy) = 0.
               upe(ix,iy) = 0.
               do ifld = 1, nfsp
                  gpix(ix,iy,ifld) = 0.
                  frici(ix,iy,ifld) = 0.
                  uu(ix,iy,ifld) = 0.
                  upi(ix,iy,ifld) = 0.
               enddo
            enddo
         endif
      endif

*  -- Set up electron parallel contribution to seec & smoc
      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            nexface = 0.5*(ne(ix2,iy)+ne(ix,iy))
            t1old =.5*cvgp*(upe(ix,iy)*rrv(ix,iy)*
     .          ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .          upe(ix1,iy)*rrv(ix1,iy)*
     .          ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2old = 1.e-20* 0.25*(fqp(ix,iy)+fqp(ix1,iy))*
     .          (ex(ix,iy)+ex(ix1,iy))/gx(ix,iy)
            iyp1 = min(iy+1,ny+1)
            iym1 = max(iy-1,0)
            t1new = .5*cvgp*( vex(ix,iy)*
     .          ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .          vex(ix1,iy)*
     .          ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2new = .5*cvgp*( vey(ix,iy)*
     .          ave(gy(ix,iy),gy(ix,iyp1))*gpey(ix,iy)/gyf(ix,iy) +
     .          vey(ix,iy)*
     .          ave(gy(ix,iy),gy(ix,iym1))*gpey(ix,iym1)/gyf(ix,iym1)
     .                                                            )
            seec(ix,iy) = seec(ix,iy)
     .          + (t1old*vol(ix,iy) - t2old)*oldseec
     .          + ((t1new+t2new)*vol(ix,iy))*(1-oldseec)
            if (nusp-isupgon(1).eq.1) smoc(ix,iy,1)=(( -cpgx*gpex(ix,iy)-
     .                   qe*nexface*gpondpotx(ix,iy) )*rrv(ix,iy)  +
     .                     pondomfpare_use(ix,iy) )*sx(ix,iy)/gxf(ix,iy)
         enddo 
      enddo

*  -- Now loop over all ion species for seec, seic, and smoc --

      do ifld = 1, nusp  #not nfsp; up only for ifld<=nusp
* ------ *
        if(zi(ifld) > 1.e-20) then  #only ions here; atoms follow
* ------ *
*     -- coupling in the x-direction --
*     -- (note sign change in pondomfpari_use term starting 031722)
           do iy = j2omp, j5omp
             do ix = i2omp, i5omp
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               tv = gpix(ix ,iy,ifld)/gxf(ix,iy)
               t1 = gpix(ix1,iy,ifld)/gxf(ix1,iy)
               t1 = .5*cvgp*( up(ix,iy,ifld)*rrv(ix,iy)*
     .                                 ave(gx(ix2,iy),gx(ix,iy))*tv
     .                      + up(ix1,iy,ifld)*rrv(ix1,iy)*
     .                                 ave(gx(ix,iy),gx(ix1,iy))*t1 )
               seic(ix,iy) = seic(ix,iy) + cfvgpx(ifld)*t1*vol(ix,iy)
               t0 = - cpiup(ifld)*( gpix(ix,iy,ifld)*rrv(ix,iy) -
     .                                  pondomfpari_use(ix,iy,ifld) )*
     .                                            sx(ix,iy)/gxf(ix,iy)
               if (nusp-isupgon(1) .eq. 1) then  # single ion mom. eq.
                  smoc(ix,iy,1) = smoc(ix,iy,1) + cpgx*t0
               else                # multiple mom. eq., so nusp=nisp
                  t0 = t0 +( qe*zi(ifld)*0.5*( ni(ix2,iy,ifld)+
     .                       ni(ix,iy,ifld) )*ex(ix,iy)*rrv(ix,iy) +
     .                       frici(ix,iy,ifld) )* sx(ix,iy)/gxf(ix,iy)
                  if (ifld <= nusp) smoc(ix,iy,ifld) = 
     .                                      smoc(ix,iy,ifld) + cpgx*t0 
               endif
c...  Add friction part of Q_e here
               tv = 0.25*(frice(ix,iy)+frice(ix1,iy))*
     .              ( upe(ix,iy)     + upe(ix1,iy) -
     .              upi(ix,iy,ifld) - upi(ix1,iy,ifld) )
               seec(ix,iy) = seec(ix,iy) - zi(ifld)**2*ni(ix,iy,ifld)*
     .                                        tv*vol(ix,iy)/nz2(ix,iy)
               
            end do
        end do

*     -- coupling in the x & y-directions --
           do iy = j2omp, j5omp
            do ix = i2omp, i5omp
             if (isgpye == 0) then
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               vyiy0 = fracvgpgp*vygp(ix,iy,ifld) 
     .                    +(1.-fracvgpgp)*vycb(ix,iy,ifld)
               vyiym1 = fracvgpgp*vygp(ix,iy-1,ifld) 
     .                    +(1.-fracvgpgp)*vycb(ix,iy-1,ifld)
               v2ix0 = fracvgpgp*v2xgp(ix,iy,ifld) 
     .                    +(1.-fracvgpgp)*v2cb(ix,iy,ifld)
               v2ixm1 = fracvgpgp*v2xgp(ix1,iy,ifld) 
     .                    +(1.-fracvgpgp)*v2cb(ix1,iy,ifld)
               t1 =.5*cvgp*( vygp(ix,iy  ,ifld)*gpiy(ix,iy  ,ifld) + 
     .                       vygp(ix,iy-1,ifld)*gpiy(ix,iy-1,ifld) +
     .                    v2xgp(ix ,iy,ifld)*ave(gx(ix,iy),gx(ix2,iy))*
     .                               gpix(ix ,iy,ifld)/gxf(ix ,iy) +    
     .                    v2xgp(ix1,iy,ifld)*ave(gx(ix,iy),gx(ix1,iy))*
     .                               gpix(ix1,iy,ifld)/gxf(ix1,iy) )     
               t2 = t1
             elseif (isgpye == 1) then    # Old B2 model with Jperp=0
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
               t2 = t1
             elseif (isgpye == 2) then    # Knoll expression
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
             endif
             seec(ix,iy) = seec(ix,iy) - fluxfacy*t1 * vol(ix,iy)
             seic(ix,iy) = seic(ix,iy) + fluxfacy*cfvgpy(ifld)*t2*
     .                                                     vol(ix,iy)
            end do
        end do

        endif  #test on zi(ifld) > 0, so only ion terms
        end do  #ifld loop over ion species

c ... Now include seic contribution from hydrogen atoms if isupgon=1
c ... Then "ion" species 2 (redundant as gas species 1) is hydr atom

      if(isupgon(1)==1 .and. zi(2)<1.e-20) then # .and. istgon(1)==0) then 
        if(cfvgpx(2) > 0.) then
          do iy = j2omp, j5omp
            do ix = i2omp, i5omp
              ix1 = ixm1(ix,iy)
              ix2 = ixp1(ix,iy)
              iy1 = max(0,iy-1)
              seic(ix,iy) = seic(ix,iy) + cftiexclg
     .                                   *0.5*cfvgpx(2)*( 
     .                       uuxg(ix, iy,1)*gpix(ix,iy,2) +
     .                       uuxg(ix1,iy,1)*gpix(ix1,iy,2) )*vol(ix,iy)
              seic(ix,iy) = seic(ix,iy) + cftiexclg
     .                                   *0.5*cfvgpy(2)*( 
     .                        vyg(ix, iy,1)*gpiy(ix,iy,2) +
     .                        vyg(ix,iy1,1)*gpiy(ix,iy1,2) )*vol(ix,iy)
            enddo
          enddo
        else  # Here if cfvgpx(2)=0, old vpar_g*grad_Pg only => ifld=2
          do iy = j2omp, j5omp
            do ix = i2omp, i5omp
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               tv = gpix(ix ,iy,2)/gxf(ix,iy)
               t1 = gpix(ix1,iy,2)/gxf(ix1,iy)
               t1 = .5*cvgp*( up(ix,iy,2)*rrv(ix,iy)*
     .                               ave(gx(ix2,iy),gx(ix,iy))*tv
     .                    + up(ix1,iy,2)*rrv(ix1,iy)*
     .                               ave(gx(ix,iy),gx(ix1,iy))*t1 )
               seic(ix,iy) = seic(ix,iy) + cftiexclg*t1*vol(ix,iy)
            enddo
          enddo
        endif  #test on cfvgpx(2) > 0 or = 0
      endif   #test for inertial neutrals


      END SUBROUTINE calc_volumetric_sources



      SUBROUTINE calc_srcmod
      IMPLICIT NONE
      Use(Selec)
      Use(Dim)
      Use(Rhsides)
      Use(Share)
      Use(Compla)
      Use(Comgeo)
      Use(Conduc)
      Use(Fixsrc)
      Use(UEpar)
      Use(Coefeq)
      Use(Comtra)
      Use(Phyvar)
      Use(Imprad)
      Use(Timing)
      Use(Xpoint_indices)
      Use(PandfTiming)
      Use(Lsode)
      Use(Bcond)
      Use(Gradients)
      Use(Cfric)
      Use(Comflo)
      Use(Jacobian_restore)
      integer iy, ix, ifld, igsp, j2pwr, j5pwr, i2pwr, i5pwr, 
     .  ix1, ix2, jg, ifld_lcs, jz, ifld_fcs, izch,  z1fac, 
     .  iyp1, iym1, iy1
      real rdumx, dr1, dr2, rdumy, ne_sgvi, t0, t1, tsimp, nevol, ngvol, 
     .  krecz, kcxrz, kionz0, kionz, kcxrzig, niz_floor, pscx0, massfac, 
     .  kionm, krecm, kcxrm, pxri, tsnpg, t1old, t2old, t1new, t2new, 
     .  vyiy0, vyiym1, v2ix0, v2ixm1, t2, nexface, nizm_floor, tv
      real rsa, rra, rcx, tick, svdiss, sv_crumpet, tock, ave
      external rsa, rra, rcx, svdiss, sv_crumpet
      ave(t0,t1) = 2*t0*t1 / (cutlo+t0+t1)
************************************************************************
*   We Calculate the source terms now.
************************************************************************
*  ---------------------------------------------------------------------
*  Coefficients for the source terms.
*  ---------------------------------------------------------------------


      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            do ifld = 1, nfsp
               snic(ix,iy,ifld) = 0.0
            enddo
            do ifld = 1, nusp
               smoc(ix,iy,ifld) = 0.0
            enddo
            seec(ix,iy) = 0.0
            seic(ix,iy) = 0.0
         end do
        end do

*****************************************************************
*  Other volume sources calculated in old SRCMOD
*****************************************************************
*  ---------------------------------------------------------------------
*  electron-ion transfer terms and an
*  approximation to the ion-ion thermal force.
*  cfw: For the neutral momentum eq there is also a v_gas*grad(p_gas)
*       term which is evaluated using ng, ti and gpiy
*  ---------------------------------------------------------------------
*-----------------------------------------------------------------------
*  -- Calculates the fixed source if it is on
*-----------------------------------------------------------------------

      if (ifixsrc .ne. 0) then
         do iy = j2omp, j5omp
            do ix = i2omp, i5omp
               snic(ix,iy,1) = snic(ix,iy,1) + vol(ix,iy) * a1n *
     .                          exp(-b1n*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1n*(yyc(iy)-yysrc)**2)
               seic(ix,iy) = seic(ix,iy) + vol(ix,iy) * a1i * ev *
     .                          exp(-b1i*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1i*(yyc(iy)-yysrc)**2)
               seec(ix,iy) = seec(ix,iy) + vol(ix,iy) * a1e * ev *
     .                          exp(-b1e*(xcs(ix)-xxsrc)**2) *
     .                          exp(-c1e*(yyc(iy)-yysrc)**2)
            end do
        end do
      endif




*  -- Set up electron parallel contribution to seec & smoc
      do iy = j2omp, j5omp
         do ix = i2omp, i5omp
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            nexface = 0.5*(ne(ix2,iy)+ne(ix,iy))
            t1old =.5*cvgp*(upe(ix,iy)*rrv(ix,iy)*
     .          ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .          upe(ix1,iy)*rrv(ix1,iy)*
     .          ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2old = 1.e-20* 0.25*(fqp(ix,iy)+fqp(ix1,iy))*
     .          (ex(ix,iy)+ex(ix1,iy))/gx(ix,iy)
            iyp1 = min(iy+1,ny+1)
            iym1 = max(iy-1,0)
            t1new = .5*cvgp*( vex(ix,iy)*
     .          ave(gx(ix,iy),gx(ix2,iy))*gpex(ix,iy)/gxf(ix,iy) +
     .          vex(ix1,iy)*
     .          ave(gx(ix,iy),gx(ix1,iy))*gpex(ix1,iy)/gxf(ix1,iy) )
            t2new = .5*cvgp*( vey(ix,iy)*
     .          ave(gy(ix,iy),gy(ix,iyp1))*gpey(ix,iy)/gyf(ix,iy) +
     .          vey(ix,iy)*
     .          ave(gy(ix,iy),gy(ix,iym1))*gpey(ix,iym1)/gyf(ix,iym1)
     .                                                            )
            seec(ix,iy) = seec(ix,iy)
     .          + (t1old*vol(ix,iy) - t2old)*oldseec
     .          + ((t1new+t2new)*vol(ix,iy))*(1-oldseec)
            if (nusp-isupgon(1).eq.1) smoc(ix,iy,1)=(( -cpgx*gpex(ix,iy)-
     .                   qe*nexface*gpondpotx(ix,iy) )*rrv(ix,iy)  +
     .                     pondomfpare_use(ix,iy) )*sx(ix,iy)/gxf(ix,iy)
         enddo 
      enddo

*  -- Now loop over all ion species for seec, seic, and smoc --

      do ifld = 1, nusp  #not nfsp; up only for ifld<=nusp
* ------ *
        if(zi(ifld) > 1.e-20) then  #only ions here; atoms follow
* ------ *
*     -- coupling in the x-direction --
*     -- (note sign change in pondomfpari_use term starting 031722)
           do iy = j2omp, j5omp
             do ix = i2omp, i5omp
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               tv = gpix(ix ,iy,ifld)/gxf(ix,iy)
               t1 = gpix(ix1,iy,ifld)/gxf(ix1,iy)
               t1 = .5*cvgp*( up(ix,iy,ifld)*rrv(ix,iy)*
     .                                 ave(gx(ix2,iy),gx(ix,iy))*tv
     .                      + up(ix1,iy,ifld)*rrv(ix1,iy)*
     .                                 ave(gx(ix,iy),gx(ix1,iy))*t1 )
               seic(ix,iy) = seic(ix,iy) + cfvgpx(ifld)*t1*vol(ix,iy)
               t0 = - cpiup(ifld)*( gpix(ix,iy,ifld)*rrv(ix,iy) -
     .                                  pondomfpari_use(ix,iy,ifld) )*
     .                                            sx(ix,iy)/gxf(ix,iy)
               if (nusp-isupgon(1) .eq. 1) then  # single ion mom. eq.
                  smoc(ix,iy,1) = smoc(ix,iy,1) + cpgx*t0
               else                # multiple mom. eq., so nusp=nisp
                  t0 = t0 +( qe*zi(ifld)*0.5*( ni(ix2,iy,ifld)+
     .                       ni(ix,iy,ifld) )*ex(ix,iy)*rrv(ix,iy) +
     .                       frici(ix,iy,ifld) )* sx(ix,iy)/gxf(ix,iy)
                  if (ifld <= nusp) smoc(ix,iy,ifld) = 
     .                                      smoc(ix,iy,ifld) + cpgx*t0 
               endif
c...  Add friction part of Q_e here
               tv = 0.25*(frice(ix,iy)+frice(ix1,iy))*
     .              ( upe(ix,iy)     + upe(ix1,iy) -
     .              upi(ix,iy,ifld) - upi(ix1,iy,ifld) )
               seec(ix,iy) = seec(ix,iy) - zi(ifld)**2*ni(ix,iy,ifld)*
     .                                        tv*vol(ix,iy)/nz2(ix,iy)
               
            end do
        end do

*     -- coupling in the x & y-directions --
           do iy = j2omp, j5omp
            do ix = i2, i5
             if (isgpye == 0) then
               ix1 = ixm1(ix,iy)
               ix2 = ixp1(ix,iy)
               vyiy0 = fracvgpgp*vygp(ix,iy,ifld) 
     .                    +(1.-fracvgpgp)*vycb(ix,iy,ifld)
               vyiym1 = fracvgpgp*vygp(ix,iy-1,ifld) 
     .                    +(1.-fracvgpgp)*vycb(ix,iy-1,ifld)
               v2ix0 = fracvgpgp*v2xgp(ix,iy,ifld) 
     .                    +(1.-fracvgpgp)*v2cb(ix,iy,ifld)
               v2ixm1 = fracvgpgp*v2xgp(ix1,iy,ifld) 
     .                    +(1.-fracvgpgp)*v2cb(ix1,iy,ifld)
               t1 =.5*cvgp*( vygp(ix,iy  ,ifld)*gpiy(ix,iy  ,ifld) + 
     .                       vygp(ix,iy-1,ifld)*gpiy(ix,iy-1,ifld) +
     .                    v2xgp(ix ,iy,ifld)*ave(gx(ix,iy),gx(ix2,iy))*
     .                               gpix(ix ,iy,ifld)/gxf(ix ,iy) +    
     .                    v2xgp(ix1,iy,ifld)*ave(gx(ix,iy),gx(ix1,iy))*
     .                               gpix(ix1,iy,ifld)/gxf(ix1,iy) )     
               t2 = t1
             elseif (isgpye == 1) then    # Old B2 model with Jperp=0
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
               t2 = t1
             elseif (isgpye == 2) then    # Knoll expression
               t1 = -0.5*( vy(ix,iy  ,ifld)*gpey(ix,iy  ) +
     .                     vy(ix,iy-1,ifld)*gpey(ix,iy-1) )
             endif
             seec(ix,iy) = seec(ix,iy) - fluxfacy*t1 * vol(ix,iy)
             seic(ix,iy) = seic(ix,iy) + fluxfacy*cfvgpy(ifld)*t2*
     .                                                     vol(ix,iy)
            end do
        end do

        endif  #test on zi(ifld) > 0, so only ion terms
        end do  #ifld loop over ion species

      END SUBROUTINE calc_srcmod


      SUBROUTINE calc_plasma_particle_residuals
      IMPLICIT NONE
      Use(Dim)
      Use(Selec)
      Use(UEpar)
      Use(Rhsides)
      Use(Compla)
      Use(Volsrc)
      Use(Coefeq)
      Use(MCN_sources)
      Use(Conduc)
      Use(Comgeo)
      Use(Phyvar)
      Use(Comflo)
      Use(Comtra)
      Use(Ext_neutrals)
      integer ifld, iy, ix, ix1, jfld

*  -- compute the residual if isnion = 1 --

      do ifld = 1, nfsp
       do iy = j2omp, j5omp
         do ix = i2omp, i5pomp
	   if(isnionxy(ix,iy,ifld) == 1) then
              resco(ix,iy,ifld) = 
     .           snic(ix,iy,ifld)+sniv(ix,iy,ifld)*ni(ix,iy,ifld) +
     .           volpsor(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psor(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psorxr(ix,iy,ifld) +
     .           cfneut * cfneutsor_ni * cnsor * psori(ix,iy,ifld) -
     .           nuvl(ix,iy,ifld)*vol(ix,iy)*ni(ix,iy,ifld) +
     .           voljcsor(ix,iy)/qe
           endif
c           if (ifld .ne. iigsp) then
	       if(zi(ifld) .ne. 0) then # IJ 2016 skip if neutral zi(ifld)=0
              resco(ix,iy,ifld) = resco(ix,iy,ifld) + cmneut * uesor_ni(ix,iy,ifld)
           else # IJ 2016 zi==0, assume neutral->ifld and ion->ifld-1
              resco(ix,iy,ifld) = resco(ix,iy,ifld) - cmneut * uesor_ni(ix,iy,ifld-1)
           endif
          end do
        end do

       do iy = j2omp, j5omp
         do ix = i2omp, i5pomp
	       if(isnionxy(ix,iy,ifld) == 1) then
              ix1 = ixm1(ix,iy)
	        if(zi(ifld) .ne. 0) then # IJ 2016 skip if neutral  zi(ifld)=0
                 resco(ix,iy,ifld) = resco(ix,iy,ifld)
     .                                   - ( (fnix(ix,iy,ifld) - fnix(ix1,iy, ifld))
     .                            + fluxfacy*(fniy(ix,iy,ifld) - fniy(ix,iy-1,ifld)) )
              else ## zi==0
                 resco(ix,iy,ifld) = resco(ix,iy,ifld)
     .              - cfneutdiv*cfneutdiv_fng*( (fnix(ix,iy,ifld)-fnix(ix1,iy, ifld))
     .                               + fluxfacy*(fniy(ix,iy,ifld)-fniy(ix,iy-1,ifld)) )
           
c ... IJ 2016/10/19 add MC neutral flux if flags set
                 if (get_neutral_moments .and. cmneutdiv_fng .ne. 0.0) then 
                    jfld=1  ## assume main ions in ifld=1
                    sng_ue(ix,iy,jfld) = - ( (fngx_ue(ix,iy,jfld) - fngx_ue(ix1,iy, jfld))
     .                        +   fluxfacy*(fngy_ue(ix,iy,jfld) - fngy_ue(ix,iy-1,jfld)) )
     .                        *( (ng(ix,iy,jfld)*ti(ix,iy))/(ng(ix,iy,jfld)*ti(ix,iy)) )
c                   if (ix .eq. 1 .and. iy .eq. 1) write(*,*) 'sng_ue', ifld, jfld
                    resco(ix,iy,ifld) = resco(ix,iy,ifld) + 
     .                                  cmneutdiv*cmneutdiv_fng*sng_ue(ix,iy,jfld)
                 endif
              endif
           endif
          end do
        end do
       enddo       # end of ifld loop
      END SUBROUTINE calc_plasma_particle_residuals


      SUBROUTINE calc_plasma_transport
      IMPLICIT NONE
      Use(Dim)
      Use(UEpar)
      Use(Selec)
      Use(Compla)
      Use(Comgeo)
      Use(Comflo)
      Use(Share)
      Use(Coefeq)
      Use(Bfield)
      Use(Xpoint_indices)
      Use(Comtra)
      Use(Bcond)
      integer ifld, methnx, methny, iy, ix, ix2, iym1, iyp1, iyp2
      real dndym1, dndy0, dndyp1, d2ndy20, d2ndy2p1, d3ndy3, 
     .  fniycboave(1:nx), sycore, corecells, t0, t1, t2
   
*****************************************************************
*****************************************************************
*  Here starts the old routine PARBAL
*****************************************************************
      do ifld = 1, nfsp
*  ---------------------------------------------------------------------
*     compute flux, residual
*     The residual is: res := snic + sniv * ni - outflow(ni).
*  ---------------------------------------------------------------------

*  -- compute fnix --

         methnx = mod(methn, 10)
         methny = methn/10
         do iy = j4omp, j8omp
            do ix = i1momp, i5pomp
              if ( (zi(ifld).eq.0.) .and. (1.-rrv(ix,iy) > 1.e-4) ) then
                 fnix(ix,iy,ifld) = fngx(ix,iy,1)
              else
               ix2 = ixp1(ix,iy)

               if (methnx .eq. 2) then   # central differencing
                  t2 = ( ni(ix, iy,ifld) + ni(ix2,iy,ifld) ) / 2

               elseif (methnx .eq. 3) then   # upwind differencing

                  if( uu(ix,iy,ifld) .ge. 0.) then
                     t2 = ni(ix,iy,ifld)
                  else
                     t2 = ni(ix2,iy,ifld)
                  endif

               elseif (methnx .eq. 6) then   # log central differencing
                  t2 = exp(0.5*
     .                ( log(ni(ix,iy,ifld)) + log(ni(ix2,iy,ifld)) ))

               else   # interp. ave or harmonic ave depending on wind*grad

                  t0 = ( ni(ix, iy,ifld)*gx(ix, iy) +
     .                   ni(ix2,iy,ifld)*gx(ix2,iy) ) / 
     .                                      ( gx(ix,iy)+gx(ix2,iy) )
                  t1 = ( gx(ix,iy)+gx(ix2,iy) ) * ni(ix,iy,ifld) *
     .                   ni(ix2,iy,ifld) / ( cutlo + ni(ix,iy,ifld)*
     .                   gx(ix2,iy) + ni(ix2,iy,ifld)*gx(ix,iy) )
                  if( uu(ix,iy,ifld)*(ni(ix,iy,ifld)-ni(ix2,iy,ifld)) 
     .                                                     .ge. 0.) then
                     t2 = t0
                  else
                     t2 = t1
                  endif

               endif

               fnix(ix,iy,ifld) = cnfx*uu(ix,iy,ifld) * sx(ix,iy) * t2
               fnixcb(ix,iy,ifld)=cnfx*sx(ix,iy) * t2 * 0.5*
     .                 (rbfbt(ix,iy) + rbfbt(ix2,iy))*v2cb(ix,iy,ifld)
                  fnix(ix,iy,ifld) = fnix(ix,iy,ifld)/sqrt( 1 +
     .              (nlimix(ifld)*ni(ix ,iy,ifld)/ni(ix2,iy,ifld))**2 +
     .              (nlimix(ifld)*ni(ix2,iy,ifld)/ni(ix ,iy,ifld))**2 )
              endif
            end do
           if ((isudsym==1.or.geometry.eq.'dnXtarget') .and. nxc > 1) then
              fnix(nxc-1,iy,ifld)=0.
              fnix(nxc,  iy,ifld)=0.
              fnix(nxc+1,iy,ifld)=0.
              uu(nxc-1,iy,ifld) = 0.
              uu(nxc  ,iy,ifld) = 0.
              uu(nxc+1,iy,ifld) = 0.
              vytan(nxc-1,iy,ifld) = 0.
              vytan(nxc  ,iy,ifld) = 0.
              vytan(nxc+1,iy,ifld) = 0.
 
           endif
           if (islimon.ne.0 .and. iy.ge.iy_lims) fnix(ix_lim,iy,ifld)=0.
           if (nxpt==2) fnix(ixrb(1)+1,iy,ifld)=0.
        end do

*  -- compute fniy  --

         do iy = j1omp1, j5omp
            do ix = i4omp, i8pomp
               if (zi(ifld).eq.0.) then #inertial gas must follow ion index
                  fniy(ix,iy,ifld) = fngy(ix,iy,ifld-1)
               else
                  if (methny .eq. 2) then   # central differencing
                     t2 = ( niy0(ix,iy,ifld) + niy1(ix,iy,ifld) ) / 2

                  elseif (methny .eq. 3) then   # upwind differencing

                     if( vy(ix,iy,ifld) .ge. 0.) then
                        t2 = niy0(ix,iy,ifld)
                     else
                        t2 = niy1(ix,iy,ifld)
                     endif

                  elseif (methny .eq. 6) then   # log central differencing
                     t2 = exp( 0.5*
     .                   (log(niy0(ix,iy,ifld))+log(niy1(ix,iy,ifld))) )

                  else    # interp. ave or harmonic ave depending on wind*grad

                     t0 = ( niy0(ix,iy,ifld)*gy(ix,iy  ) +
     .                    niy1(ix,iy,ifld)*gy(ix,iy+1) ) / 
     .                    ( gy(ix,iy)+gy(ix,iy+1) )
                     t1 = ( gy(ix,iy)+gy(ix,iy+1) ) * niy0(ix,iy,ifld)*
     .                    niy1(ix,iy,ifld) / ( cutlo + niy0(ix,iy,ifld)*
     .                    gy(ix,iy+1) + niy1(ix,iy,ifld)*gy(ix,iy) )
                     if( (niy0(ix,iy,ifld)-niy1(ix,iy,ifld))*
     .                    vy(ix,iy,ifld) .ge. 0.) then
                        t2 = t0
                     else
                        t2 = t1
                     endif
                  
                  endif
               
                  fniy(ix,iy,ifld) = cnfy*vy(ix,iy,ifld)*sy(ix,iy)*t2
                  fniycb(ix,iy,ifld) = cnfy*vycb(ix,iy,ifld)*sy(ix,iy)*t2
                  if (vy(ix,iy,ifld)*(ni(ix,iy,ifld)-ni(ix,iy+1,ifld))
     .                                                      .lt. 0.) then
                     fniy(ix,iy,ifld) = fniy(ix,iy,ifld)/( 1 +
     .                               (nlimiy(ifld)/ni(ix,iy+1,ifld))**2 +
     .                               (nlimiy(ifld)/ni(ix,iy  ,ifld))**2 )
                  endif
c...  Note: nonorthogonality comes in through calc. of vy
               endif
            end do
        end do

c ... cosmetic setting of fniy - not used         
         do ix = i4omp, i8omp
            fniy(ix,ny+1,ifld) = 0.0e0
         enddo

        end do

c ... Add rad flux of 4th order diff operator; damp grid-scale oscillations
      do ifld = 1, nfsp
        if (abs(dif4order(ifld)) > 1.e-50) then
          do iy = j2p, j5m   #limits to range iy=1:ny-1 for fniy4ord
            iym1 = max(iy-1,0)
            iyp1 = min(iy+1,ny+1)
            iyp2 = min(iy+2,ny+1)
            do ix = i4omp, i8omp
              dndym1 = (ni(ix,iy,ifld)-ni(ix,iym1,ifld))*gyf(ix,iym1)
              dndy0 = (ni(ix,iyp1,ifld)-ni(ix,iy,ifld))*gyf(ix,iy)
              dndyp1 = (ni(ix,iyp2,ifld)-ni(ix,iyp1,ifld))*gyf(ix,iyp1)
              d2ndy20 = (dndy0 - dndym1)*gy(ix,iy)
              d2ndy2p1 = (dndyp1 - dndy0)*gy(ix,iyp1)
              d3ndy3 = (d2ndy2p1 - d2ndy20)*gyf(ix,iy)
              fniy4ord(ix,iy,ifld) = dif4order(ifld)*d3ndy3*sy(ix,iy)/
     .                                                  gyf(ix,iy)**2
              fniy(ix,iy,ifld) = fniy(ix,iy,ifld) + fniy4ord(ix,iy,ifld)
            enddo
          enddo
        endif
      enddo

 
      END SUBROUTINE calc_plasma_transport


      SUBROUTINE calc_fniycbo
      IMPLICIT NONE
      Use(Dim)
      Use(UEpar)
      Use(Selec)
      Use(Compla)
      Use(Comgeo)
      Use(Comflo)
      Use(Share)
      Use(Coefeq)
      Use(Bfield)
      Use(Xpoint_indices)
      Use(Comtra)
      Use(Bcond)
      integer ifld, iy, ix, jx
      real sycore, corecells, fniycboave
c ... Setup a correction to surface-flux for grad_B and grad_P effects at iy=0
      fniycbo = 0.0
      do ifld = 1, nfsp
        fniycbo(:,ifld) = ( ni(:,0,ifld)*sy(:,0) ) *
     .                         ( (1-cfniybbo)*cfybf*vycb(:,0,ifld) -
     .                            cfniydbo*(1-cfydd)*vycp(:,0,ifld) )
      enddo

c ... Normalize core flux to zero to avoid introducing artifical core source/sink
      do ifld = 1, nfsp
          if (isfniycbozero(ifld) .gt. 0) then 
            fniycboave = 0
            corecells = 0
            sycore = 0
            do jx = 1, nxpt
                fniycboave = fniycboave + SUM(fniycbo(ixpt1(jx)+1:ixpt2(jx), ifld))
                corecells =  corecells + (ixpt2(jx) - ixpt1(jx))
            end do
c ... TODO: Add double-null fix here (now only does one half-mesh...)            
            fniycboave = fniycboave / corecells
            do jx = 1, nxpt
                do ix = ixpt1(jx)+1, ixpt2(jx)
                  fniycbo(ix, ifld) = fniycbo(ix, ifld) - isfniycbozero(ifld)*fniycboave
                end do
            end do
          else if (isfniycbozero(ifld) .lt. 0) then 
            do jx = 1, nxpt
              fniycbo(ixpt1(jx)+1:ixpt2(jx), ifld) = 0

            end do
          end if
      end do
             
 
      END SUBROUTINE calc_fniycbo


