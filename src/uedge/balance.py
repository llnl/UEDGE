



class UeBalance():
    def __init__(self):
        from numpy import zeros
        from uedge import com
        for var in [
            'fetx', 'fety', 'engerr', 'pmloss', 'pmrada', 'pmradm', 'pmpot',
            'peirad', 'pmomv', 'engerr', 'pradrc', 'pradiz', 'pradht', 'prdiss',
            'pibirth', 'pbinde', 'pbindrc', 'pradzbind', 'pradff'
        ]:
            self.__dict__[var] = zeros((com.nx+2, com.ny+2))
        for var in ['icxgas', 'iion', 'irecomb']:
            self.__dict__[var] = zeros((com.nx+2, com.ny+2, com.ngsp))

        self.pradimp = zeros((com.nx+2, com.ny+2, sum(com.nzsp), sum(com.nzsp)))
            
        return

    def ave(x, y, cutlo=1e-300):
        return (x * y) / (x + y + cutlo)


    def engbal(self):
        """ Calculates various components of the 2-D energy flow and the 
        ionization and radiation for use in the postprocessing file
        balancee to determine energy balance; these 2-D loops become 
        expensive when done from the parser.
        """
        from uedge import bbb, com
        from numpy import zeros, cos

        '''
 c-----------------------------------------------------------------------
      subroutine engbal(pwrin)
        '''

        
        if bbb.ishosor != 0:
            raise NotImplementedError("Option ishosor>0 not implemented")
            """  
            call volavenv(nx, ny, 1, ny, 1, nx, ixp1(0:nx+1,0:ny+1), 
                ixm1(0:nx+1,0:ny+1), fsprd, psor_tmpov(0:nx+1,0:ny+1), prad)
            do igsp = nhgsp+1, ngsp
                jz = igsp - nhgsp
                do iimp = 0, nzsp(jz)       
                    call volavenv(nx, ny, 1, ny, 1, nx, ixp1(0:nx+1,0:ny+1), 
                        ixm1(0:nx+1,0:ny+1), fsprd, psor_tmpov(0:nx+1,0:ny+1), 
                        pradz(0:nx+1,0:ny+1,iimp,jz))
                enddo
            enddo
            """
        upi = bbb.upi

        # Set arrays to check energy conserv; add ion parallel drift and visc heat
        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx], com.ixrb[jx]+1):
                for iy in range(0, com.ny+1): 
                    ix1 = bbb.ixm1[ix, iy]
                    ix2 = bbb.ixp1[ix, iy]
                    for ii in range(0, com.nusp):
                        thetaix = 0.5 * (com.angfx[ix1,iy] + com.angfx[ix,iy])
                        thetaix2 = 0.5 * (com.angfx[ix,iy] + com.angfx[ix2,iy])
                        eta_dup2dy = 0.25*(
                            bbb.visy[ix, iy+1, ii]*(
                                upi[ix1, iy+1, ii] + upi[ix, iy+1, ii]
                            )**2 - bbb.visy[ix,iy,ii]*(
                                upi[ix1,iy,ii] + upi[ix,iy,ii]
                            )**2 )

                        self.fety[ix,iy] += (bbb.mi[ii]/32)*(
                                upi[ix1,iy,ii] + upi[ix,iy,ii] \
                                + upi[ix1,iy+1,ii] + upi[ix,iy+1,ii]\
                            )**2 * bbb.fnix[ix,iy,ii] \
                        -bbb.cfvisy*0.5*com.sy[ix,iy]*com.gyf[ix,iy]\
                            *eta_dup2dy
                    
                        self.fetx[ix,iy] += 0.5*bbb.mi[ii]*upi[ix,iy,ii]**2 \
                            *bbb.fnix[ix,iy,ii] - bbb.cfvisx*0.25*com.sx[ix,iy]*( \
                                bbb.visx[ix,iy,ii]*com.gx[ix,iy]*cos(thetaix)* \
                                (upi[ix,iy,ii]**2 - upi[ix1,iy,ii]**2) +
                                bbb.visx[ix2,iy,ii]*com.gx[ix2,iy]*cos(thetaix2)* \
                                (upi[ix2,iy,ii]**2 - upi[ix,iy,ii]**2)
                            ) - upi[ix,iy,ii] * bbb.fmixy[ix,iy,ii]
                    self.fety[ix,iy] += bbb.feey[ix,iy] + bbb.feiy[ix,iy]
                    self.fetx[ix,iy] += bbb.feex[ix,iy] + bbb.feix[ix,iy]

        '''
      do jx = 1, nxpt
        do ix=ixlb(jx),ixrb(jx)
          do iy=0,ny
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
            fety(ix,iy) = 0.
            fetx(ix,iy) = 0.
            do id = 1, nusp
	       thetaix =  0.5*(angfx(ix1,iy) + angfx(ix,iy))
               thetaix2 = 0.5*(angfx(ix,iy) + angfx(ix2,iy))
               eta_dup2dy = 0.25*( visy(ix,iy+1,id)*
     .                       (upi(ix1,iy+1,id)+upi(ix,iy+1,id))**2 -
     .                             visy(ix,iy  ,id)*
     .                       (upi(ix1,iy  ,id)+upi(ix,iy  ,id))**2 )
               fety(ix,iy) = fety(ix,iy) + (mi(id)/32)*
                            ( upi(ix1,iy,id)+
     .                         upi(ix,iy,id)+upi(ix1,iy+1,id)+
     .                         upi(ix,iy+1,id) )**2*fniy(ix,iy,id) -
     .                         cfvisy*0.5*sy(ix,iy)*gyf(ix,iy)*eta_dup2dy



               fetx(ix,iy) = fetx(ix,iy) + 0.5*mi(id)*upi(ix,iy,id)**2*
     .                          fnix(ix,iy,id) - cfvisx*0.25*sx(ix,iy)*(
     .                         visx(ix ,iy,id)*gx(ix ,iy)*cos(thetaix)* 
     .                          ( upi(ix,iy,id)**2 - upi(ix1,iy,id)**2 ) +

     .                        visx(ix2,iy,id)*gx(ix2,iy)*cos(thetaix2)* 
     .                          ( upi(ix2,iy,id)**2 - upi(ix,iy,id)**2 ) )
               fetx(ix,iy) = fetx(ix,iy) - upi(ix,iy,id)*fmixy(ix,iy,id)
            enddo
            fety(ix,iy) = fety(ix,iy) + feey(ix,iy) + feiy(ix,iy)
            fetx(ix,iy) = fetx(ix,iy) + feex(ix,iy) + feix(ix,iy)
          enddo
        enddo
      enddo
    '''
    # Correct the boundary x-fluxes if non-unity ckinfl
        up = bbb.up
        if (abs(bbb.ckinfl - 1) > 1e-10):
            for jx in range(com.nxpt):
                ixt = com.ixlb[jx]
                ixt1 = ixt + 1
                ixr = com.ixrb[jx]
                ixr1 = ixr - 1
                for iy in range(com.ny+1):
                    self.fetx[ixt,iy] = 0
                    self.fetx[ixr,iy] = 0
                    for ii in range(com.nfsp):
                        self.fetx[ixt,iy] += 0.5*bbb.mi[ii]*up[ixt,iy,ii]**2*bbb.fnix[ixt,iy,ii] \
                        - bbb.ckinfl*0.5*com.sx[ixt,iy]*bbb.visx[ixt1,iy,ii]*com.gx[ixt1,iy]* \
                        (up[ixt1, iy, ii]**2 - up[ixt, iy, ii]**2)

                        self.fety[ixr,iy] += 0.5*bbb.mi[ii]*up[ixr,iy,ii]**2*bbb.fnix[ixr,iy,ii] \
                        - bbb.ckinfl*0.5*com.sx[ixr,iy]*bbb.visx[ixr,iy,ii]*com.gx[ixr,iy] * \
                        (up[ixr,iy,ii]**2 - up[ixr1,iy,ii]**2)
                self.fetx[ixt,iy] += bbb.feex[ixt,iy] + bbb.feix[ixt,iy]
                self.fetx[ixr,iy] += bbb.feex[ixr,iy] + bbb.feix[ixr,iy]
        '''

# Now correct the boundary x-fluxes if non-unity ckinfl
      if (abs(ckinfl-1.) > 1.e-10) then
       do jx = 1, nxpt
         ixt  = ixlb(jx)
         ixt1 = ixt + 1
         ixr  = ixrb(jx)
         ixr1 = ixr - 1
         do 15 iy = 0, ny
            fetx(ixt,iy) = 0.
            fetx(ixr,iy) = 0.
            do id = 1, nfsp
               fetx(ixt,iy) = fetx(ixt,iy) +
     .                        0.5*mi(id)*up(ixt,iy,id)**2*fnix(ixt,iy,id) -
     .                   ckinfl*0.5*sx(ixt,iy)*visx(ixt1,iy,id)*gx(ixt1,iy)*
     .                          ( up(ixt1,iy,id)**2 - up(ixt,iy,id)**2 )



               fetx(ixr,iy) = fetx(ixr,iy) +
     .                        0.5*mi(id)*up(ixr,iy,id)**2*fnix(ixr,iy,id) - 
     .                   ckinfl*0.5*sx(ixr,iy)*visx(ixr,iy,id)*gx(ixr,iy)*
     .                          ( up(ixr,iy,id)**2 - up(ixr1,iy,id)**2 )
            enddo
            fetx(ixt,iy) = fetx(ixt,iy) + feex(ixt,iy) + feix(ixt,iy)
            fetx(ixr,iy) = fetx(ixr,iy) + feex(ixr,iy) + feix(ixr,iy)
 15      continue
       enddo  # end do-loop over nxpt mesh regions
      endif   # test on ckinfl-1

    
        '''
        pwrin = 1

        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx], com.ixrb[jx]+1):
                for iy in range(1,com.ny+1):
                    if bbb.ishymol:
                        self.pmloss[ix,iy] = (1-bbb.ismolcrm)*bbb.cnsor*( \
                                bbb.ediss*bbb.ev*(0.5*bbb.psordis[ix,iy,1]) \
                                + bbb.ceisor*bbb.eion*bbb.ev*(bbb.psordis[ix,iy,1])
                            ) + bbb.ismolcrm*bbb.cnsor*( \
                                bbb.cmesori*(bbb.emolia[ix,iy,0] + bbb.emolia[ix,iy,1]) \
                                + bbb.cmsore*bbb.edisse[ix,iy]
                            )
                        self.pmpot[ix,iy] = bbb.ismolcrm*bbb.ng[ix,iy,1]*com.vol[ix,iy] \
                            *sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], 22)
                        self.pmrada[ix,iy] = bbb.ismolcrm*bbb.ng[ix,iy,1]*com.vol[ix,iy] \
                            *sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], 23)
                        self.pmradm[ix,iy] = bbb.ismolcrm*bbb.ng[ix,iy,1]*com.vol[ix,iy] \
                            *sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], 24)
# Here peirad includes sum of electron and ion energy losses; note that binding
# energy is included in eeli term, but it is carried by the ions.
# Note also that eion and ediss generally balance in the next line
# because should have ediss=2*eion - transfer from electron to ion energy
                    self.peirad[ix,iy] = bbb.cnsor*(
                            bbb.ebind*bbb.ev*bbb.psor[ix,iy,0] \
                            - bbb.ebind*bbb.ev*bbb.psorrg[ix,iy,0] \
                            + self.pmloss[ix,iy]
                        )

                    if (bbb.isupgon[0] == 0):
                        self.pmomv[ix,iy] = bbb.cngmom[0]*up[ix,iy,0]*com.sx[ix,iy] \
                            *com.rrv[ix,iy]*(
                                bbb.ng[ix2,iy,0]*bbb.tg[ix2,iy,0] \
                                - bbb.ng[ix,iy,0]*bbb.tg[ix,iy,0]
                            ) + bbb.cmwall[0]*0.125*bbb.mi[0]*(
                                up[ix,iy,0] + up[ix1,iy,1] 
                            )**2*bbb.ng[ix,iy,0]*bbb.nucx[ix,iy,0]*com.vol[ix,iy]
                    self.engerr[ix,iy] = (
                            self.fetx[ix1,iy] - self.fetx[ix,iy] + self.fety[ix,iy-1] \
                            - self.fety[ix,iy] - self.peirad[ix,iy] - bbb.png2ni[ix,iy]
                        ) / abs(pwrin)
                    if bbb.isimpon != 0:
                        self.engerr[ix,iy] -= bbb.prad[ix,iy]*com.vol[ix,iy]/abs(pwrin)                        
        '''
      pvmomcx = 0.e0
      ptjdote = 0.e0 
      do jx = 1, nxpt
        do ix=ixlb(jx)+1,ixrb(jx)
          do iy=1,ny
            ix1 = ixm1(ix,iy)
            ix2 = ixp1(ix,iy)
                pmloss(ix,iy) =(1-ismolcrm)*cnsor*(ediss*ev*(0.5*psordis(ix,iy,2))+
     .                      ceisor*eion*ev*(psordis(ix,iy,2)) ) + 
     .                      ismolcrm*cnsor*(
     .                            cmesori*(emolia(ix,iy,1) + emolia(ix,iy,2))
     .                          + cmesore*edisse(ix,iy)
     .                      )
                pmpot(ix,iy) = ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                          sv_crumpet(te(ix,iy), ne(ix,iy), 22)
                pmrada(ix,iy) = ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                          sv_crumpet(te(ix,iy), ne(ix,iy), 23)
                pmradm(ix,iy) = ismolcrm*ng(ix,iy,2)*vol(ix,iy)*
     .                          sv_crumpet(te(ix,iy), ne(ix,iy), 24)
            peirad(ix,iy) = cnsor*( erliz(ix,iy) + erlrc(ix,iy) +
     .                              ebind*ev*psor(ix,iy,1) -
     .                              ebind*ev*psorrg(ix,iy,1) +
     .                              pmloss(ix,iy))
# other energy diagnostics are given below
cc            jdote(ix,iy) = -   # this energy is included in resee, not lost
cc     .                  0.5 * fqx(ix ,iy)*(phi(ix2,iy  )+phi(ix ,iy)) +
cc     .                  0.5 * fqx(ix1,iy)*(phi(ix ,iy  )+phi(ix1,iy)) -
cc     .                  0.5 * fqy(ix ,iy)*(phi(ix ,iy+1)+phi(ix ,iy)) +
cc     .                  0.5 * fqy(ix,iy-1)*(phi(ix,iy)+phi(ix,iy-1))
            ptjdote = ptjdote + wjdote(ix,iy)

            if (isupgon(1) .eq. 0) then
               pmomv(ix,iy)=cngmom(1)*up(ix,iy,1)*sx(ix,iy)*rrv(ix,iy)*
     .                           ( ng(ix2,iy,1)*tg(ix2,iy,1)- 
     .                             ng(ix ,iy,1)*tg(ix ,iy,1) ) +
     .             cmwall(1)*0.125*mi(1)*(up(ix,iy,1)+up(ix1,iy,1))**2*
     .                ng(ix,iy,1)*nucx(ix,iy,1)*vol(ix,iy)
            elseif (isupgon(1) .eq. 1) then    # inertial neutrals
               pmomv(ix,iy) = 0.  # coupled back to therm eng for inertial neut
            endif
            pvmomcx = pvmomcx + pmomv(ix,iy)

            engerr(ix,iy) = ( fetx(ix1,iy  )-fetx(ix,iy)+
     .                        fety(ix ,iy-1)-fety(ix,iy)-
     .                        peirad(ix,iy)-png2ni(ix,iy) ) /
     .                         abs(pwrin)
            if (isimpon.ne.0) then  # prad allocated only if isimpon.ne.0
               engerr(ix,iy) = engerr(ix,iy) - prad(ix,iy)*vol(ix,iy)/
     .                                                     abs(pwrin)
            endif
          enddo
        enddo
      enddo

        '''

        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx]+1, com.ixrb[jx]+1):
                for iy in range(1, com.ny+1):
                    for ig in range(1, com.ngsp):
                        if ((bbb.ishymol == 0) or (ig != 1)):
                            self.iion[ix,iy,ig] -= bbb.cnsor*bbb.qe*bbb.psorg[ix,iy,ig]
                            self.irecomb[ix,iy,ig] -= bbb.cnsor*bbb.qe*bbb.psorrg[ix,iy,ig]
                            self.icxgas[ix,iy,ig] -= bbb.qe*bbb.psorcxg[ix,iy,ig]
        
                    '''
# ionization and background sources

      iion_tot = 0.e0
      irecomb_tot = 0.e0
      icxgas_tot = 0.e0
      pradrc = 0.e0
      pradiz = 0.e0
      pradht = 0.e0
      prdiss = 0.e0
      pibirth = 0.e0
      pbinde = 0.e0
      pbindrc = 0.e0
      pradzbind = 0.e0
      do igsp = 1, ngsp
         iion(igsp) = 0.e0
         irecomb(igsp) = 0.e0
         icxgas(igsp) = 0.e0
      enddo
      do igsp = 1, max(1, ngsp-nhgsp)
         pradimpt(igsp) = 0.
         if (nzsp(igsp) .ne. 0) then
            do iimp = 0, nzsp(igsp)
               pradimp(iimp,igsp) = 0.
            enddo
         endif
      enddo    
      pradfft = 0.
      
      do jx = 1, nxpt
        do ix = ixlb(jx)+1,ixrb(jx)
          do iy = 1, ny
            do igsp = 1, ngsp
              if (ishymol.eq.0 .or. igsp.ne.2) then
               iion(igsp) = iion(igsp) - cnsor*qe*psorg(ix,iy,igsp)
               irecomb(igsp) = irecomb(igsp) -cnsor*qe*psorrg(ix,iy,igsp)
               icxgas(igsp) = icxgas(igsp) - qe*psorcxg(ix,iy,igsp)
              endif
            enddo
          enddo
        enddo
      enddo

      do igsp = 1, ngsp
         iion_tot = iion_tot + iion(igsp)
         irecomb_tot = irecomb_tot + irecomb(igsp)
         icxgas_tot = icxgas_tot + icxgas(igsp)
      enddo
                    '''
                    self.pradrc[ix,iy] += bbb.cnsor*bbb.erlrc[ix,iy]
                    self.pradiz[ix,iy] += (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev) * bbb.psor[ix,iy,0]
                    self.pbinde[ix,iy] += bbb.ebind*bbb.ev*bbb.psor[ix,iy,0]
                    self.pbindrc[ix,iy] += bbb.ebind*bbb.ev*bbb.psorrg[ix,iy,0]
                    self.prdiss[ix,iy] += (1-bbb.ismolcrm)*(bbb.ediss*bbb.ev*(0.5*bbb.psordis[ix,iy,1])) \
                            + bbb.ismolcrm*bbb.cmesore*bbb.edisse[ix,iy]
                    self.pibirth[ix, iy] += (1-bbb.ismolcrm) * (bbb.ceisor*bbb.eion*bbb.ev \
                            *bbb.psordis[ix,iy,1] - bbb.ccoldsor*bbb.ng[ix,iy,0]*(
                                1.5*bbb.ti[ix,iy] - bbb.eion*bbb.ev)*bbb.nucx[ix,iy,0]*com.vol[ix,iy]) \
                            + bbb.ismolcrm*(
                                bbb.ceisor*bbb.cmesore*(bbb.emolia[ix,iy,0] + bbb.emolia[ix,iy,1]) \
                                - bbb.ccoldsor*bbb.ng[ix,iy,0]*(
                                    1.5*bbb.ti[ix,iy] - bbb.eion*bbb.ev
                                )*bbb.nucx[ix,iy,0]*com.vol[ix,iy]
                            )
                                
            
        '''

      do iy = 1, ny
        do jx = 1, nxpt
          do ix = ixlb(jx)+1, ixrb(jx)
            ix1 = ixm1(ix,iy)
            pradrc = pradrc + cnsor*erlrc(ix,iy)
            pradiz = pradiz + (eeli(ix,iy)-ebind*ev) * psor(ix,iy,1) 
            pbinde = pbinde + ebind*ev * psor(ix,iy,1)
            pbindrc = pbindrc + ebind*ev*psorrg(ix,iy,1)
            prdiss = prdiss+(1-ismolcrm)*(ediss*ev*(0.5*psordis(ix,iy,2)))+
     .                              ismolcrm*cmesore*edisse(ix,iy)
            pibirth = pibirth+(1-ismolcrm)*(ceisor* eion*ev * (psordis(ix,iy,2)) -
     .                ccoldsor*ng(ix,iy,1)*(1.5*ti(ix,iy)-eion*ev)*
     .                                          nucx(ix,iy,1)*vol(ix,iy) )+
     .                  ismolcrm*( ceisor*cmesore*(emolia(ix,iy,1) + emolia(ix,iy,2))
     .                   - ccoldsor*ng(ix,iy,1)*(1.5*ti(ix,iy)-eion*ev)*
     .                                          nucx(ix,iy,1)*vol(ix,iy) )
         enddo
        enddo
      enddo
      pradht = pradiz + pradrc
        '''
        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx], com.ixrb[jx]+1):
                for iy in range(1,com.ny+1):
                    if (bbb.isimpon in [2,7]):
                        self.pradff[ix,iy] += bbb.pracff[ix,iy]*com.vol[ix,iy]
                    if bbb.isimpon > 2:
                        for ig in range(com.nhgsp+1, com.ngsp+1):
                            for jz in range(com.ngsp - com.nhgsp):
                                for iimp in range(com.nzsp[jz]):
                                    ii = iimp + com.nhgsp + sum(com.nzsp[:jz]) - 1
                                    self.pradimp[ix,iy, ii, jz] += bbb.pradz[ix,iy,ii,jz]*com.vol[ix,iy]
                            
                        

        '''
   
      if (isimpon .eq. 2 .or. isimpon .eq. 7) then #fixed fraction model
         do iy = 1, ny
           do jx = 1, nxpt
             do ix = ixlb(jx)+1, ixrb(jx)
               pradfft = pradfft + pradcff(ix,iy)*vol(ix,iy)
             enddo
           enddo
         enddo
      endif

      if (isimpon .gt. 2) then  #separate impurity species
         do iy = 1, ny
           do jx = 1, nxpt
             do ix = ixlb(jx)+1, ixrb(jx)
               do igsp = nhgsp+1, ngsp
                  jz = igsp - nhgsp
                  do iimp = 0, nzsp(jz)
                     pradimp(iimp,jz) = pradimp(iimp,jz) +
     .                             pradz(ix,iy,iimp,jz)*vol(ix,iy)
                  enddo
               enddo
               pradzbind = pradzbind + (pwrze(ix,iy)-prad(ix,iy))*
     .                                vol(ix,iy) # only total pradzbind calc
             enddo
           enddo
         enddo
         do igsp = nhgsp+1, ngsp
            jz = igsp - nhgsp
            do iimp = 0, nzsp(jz)
               pradimpt(jz) = pradimpt(jz) + pradimp(iimp,jz)
            enddo
         enddo
      endif
  
      return
      end          

******* end of subroutine engbal *******

        '''

    def pradpltwl(self):
        """ Calc radiation power to divertor and outer wall surfaces """
        from uedge import bbb, com
        from numpy import zeros, arctan2, pi, cos
        from copy import deepcopy
        nz = com.nxomit 
        prdu = deepcopy(bbb.prad)
        '''
c----------------------------------------------------------------------c
      subroutine pradpltwl
c ... Local variables
      real prdu(0:nx+1,0:ny+1)   
      real theta_ray1, theta_ray2, dthgy, dthgx, sxo, frth
      integer ixv, iyv, nj, ix, iy, ip, iodd

# Initialize arrays
      call sfill ((ny+2)*2*nxpt, 0., pwr_pltz, 1)   
      call sfill ((ny+2)*2*nxpt, 0., pwr_plth, 1)   
      call sfill ((nx+2), 0., pwr_wallz, 1)   
      call sfill ((nx+2), 0., pwr_wallh, 1) 
      call sfill ((nx+2)*nxpt, 0., pwr_pfwallz, 1)   
      call sfill ((nx+2)*nxpt, 0., pwr_pfwallh, 1) 
  
      if (isimpon > 0) then   # use prdu since prad might not be dimensioned
         call s2copy (nx+2,ny+2, prad,1,nx+2, prdu,1,nx+2) #prad --> prdu
      else
         call s2fill (nx+2, ny+2, 0., prdu, 1, nx+2)
      endif

      nj = nxomit
        '''
        """ Divertor plates """
        for ip in range(2*com.nxpt):
            iodd = int((ip % 2)==1)
            if iodd:
                ixv = com.ixlb[ip/2 + 1]
            else:
                ixv = ixrb[ip/2] + 1
            for iyv in range(1, com.ny+1):
                for iy in range(1, com.ny):
                    for ix in range(1,com.nx+1):
                        thetaray1 = arctan2(
                            com.zm[ixv+nj, iyv, 1] - com.zm[ix+nj, iy, 0],
                            com.rm[ixv+nf, iyv, 1] - com.rm[ix+nj, iy, 0]
                        )
                        thetaray2 = arctan2(
                            com.zm[ixv+nj, iyv, 3] - com.zm[ix+nj, iy, 0],
                            com.rm[ixv+nf, iyv, 3] - com.rm[ix+nj, iy, 0]
                        )
                        dthgy = abs(thetaray1 - thetaray2)
                        frth = min(dthgy, 2*pi - dthgy)/2/pi
                        sxo = com.sx[ixv, iyv]/cos(com.angfx[ixv,iyv])
                        pwr_pltz[iyv,ip] += prdu[ix,iy]*com.vol[ix,iy]*fth/sxo
                        pwr_plth[iyv,ip] += ( (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev) \
                            *bbb.psor[ix,iy,0]) + bbb.erlrc[ix,iy]*frth/sxo
            # Set corner values
            pwr_pltz[0,ip] = pwr_pltz[1,ip]
            pwr_pltz[com.ny+1,ip] = pwr_pltz[com.ny,ip]
            pwr_plth[0,ip] = pwr_plth[1,ip]
            pwr_plth[com.ny+1,ip] = pwr_plth[com.ny,ip]
                            
            


        '''

c ... First do the divertor plate surfaces
      do ip = 1, 2*nxpt
        iodd = (ip+1)/2 - ip/2  # =1 if ip odd, =0 if ip even
	if (iodd==1) then
          ixv = ixlb(ip/2+1)    # viewing ix position
	else
	  ixv = ixrb(ip/2) + 1
	endif
        do iyv = 1, ny		# loop over viewing ix
          do iy = 1, ny	        # loop over source iy
            do ix = 1, nx	# loop over source ix
              theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0), 
     .                            rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
              theta_ray2 = atan2( zm(ixv+nj,iyv,3)-zm(ix+nj,iy,0), 
     .                            rm(ixv+nj,iyv,3)-rm(ix+nj,iy,0) )
              dthgy = abs(theta_ray1-theta_ray2)
              frth = min(dthgy, 2*pi-dthgy)/(2*pi)  # frac.; need angle < pi
              sxo = sx(ixv,iyv)/(cos(angfx(ixv,iyv)))
              pwr_pltz(iyv,ip) = pwr_pltz(iyv,ip) + 
     .                                 prdu(ix,iy)*vol(ix,iy)*frth/sxo
              pwr_plth(iyv,ip) = pwr_plth(iyv,ip) + (
     .                           (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1)
     .                                        + erlrc(ix,iy))*frth/sxo
            enddo
          enddo
        enddo
c ... Set corner values
      pwr_pltz(0,ip)    = pwr_pltz(1,ip)	
      pwr_pltz(ny+1,ip) = pwr_pltz(ny,ip)
      pwr_plth(0,ip)    = pwr_plth(1,ip)	
      pwr_plth(ny+1,ip) = pwr_plth(ny,ip)

      enddo             # end of ip loop over divertor plates

        '''
        """ OUTER WALL """
        iyv = com.ny + 1
        for ixv in range(1, nx+1):
            for iy in range(1, com.ny+1):
                for ix in range(1, com.nx+1):
                    thetaray1 = arctan2(
                        com.zm[ixv+nj, iyv, 1] - com.zm[ix+nj, iy, 0],
                        com.rm[ixv+nj, iyv, 1] - com.rm[ix+nj, iy, 0],
                    )
                    thetaray2 = arctan2(
                        com.zm[ixv+nj, iyv, 2] - com.zm[ix+nj, iy, 0],
                        com.rm[ixv+nj, iyv, 2] - com.rm[ix+nj, iy, 0],
                    )
                    dthgz = abs(thetaray1 - thetaray2)
                    frth = min(dthgz, 2*pi - dthgx)/2/pi
                    pwr_wallz[ixv] += prdu[ix,iy]*com.vol[ix,iy]*frth/com.sy[ixv,iyv]
                    pwr_wallh[ixv] += ( (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev)*bbb.psor[ix,iy,0] \
                        +bbb.erlrc[ix,iy])*frth/com.sy[ixv,iyv]
        pwr_wallz[0] = pwr_wallz[1]	# Because prad(0,) is not calculated
        pwr_wallz[com.nx+1] = pwr_wallz[com.nx]
        pwr_wallh[0] = pwr_wallh[1]
        pwr_wallh[com.nx+1] = pwr_wallh[com.nx]
    

        '''

c ... Now do the "outer" wall surface, i.e., iy=ny+1
      iyv = ny+1                # viewing iy position
      do ixv = 1, nx		# loop over viewing ix
        do iy = 1, ny	        # loop over source iy
          do ix = 1, nx	        # loop over source ix
            theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
            theta_ray2 = atan2( zm(ixv+nj,iyv,2)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,2)-rm(ix+nj,iy,0) )
            dthgx = abs(theta_ray1-theta_ray2)             
            frth = min(dthgx, 2*pi-dthgx)/(2*pi)  # frac; need angle < pi
            pwr_wallz(ixv) = pwr_wallz(ixv) + 
     .                           prdu(ix,iy)*vol(ix,iy)*frth/sy(ixv,iyv) 
            pwr_wallh(ixv) = pwr_wallh(ixv) + (
     .                       (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1) +
     .                          erlrc(ix,iy) )*frth/sy(ixv,iyv) 
          enddo
        enddo
      enddo
      pwr_wallz(0) = pwr_wallz(1)	# Because prad(0,) is not calculated
      pwr_wallz(nx+1) = pwr_wallz(nx)
      pwr_wallh(0) = pwr_wallh(1)
      pwr_wallh(nx+1) = pwr_wallh(nx)

        '''
        """ PFR WALL """
        iyv = 0
        for ip in range(com.nxpt):
            for ixv in range(1, com.nx+1):
                for iy in range(1, com.ny+1):
                    for ix in range(1, com.nx+1):
                        thetaray1 = arctan2(
                            com.zm[ixv+nj, iyv, 1] - com.zm[ix+nj,iy,0],
                            com.rm[ixv+nj, iyv, 1] - com.rm[ix+nj,iy,0],
                        )
                        thetaray2 = arctan2(
                            com.zm[ixv+nj, iyv, 2] - com.zm[ix+nj,iy,0],
                            com.rm[ixv+nj, iyv, 2] - com.rm[ix+nj,iy,0],
                        )
                        dthgx = abs(thetaray1 - thetaray2)
                        frth = min(dthgx, 2*pi - dthgx)/2/pi
                        pwr_pfwallz[ixv,ip] += prdu[ix,iy]*com.vol[ix,iy]*frth/com.sy[ixv,iyv]
                        pwr_pfwallh[ixv,ip] += (
                            (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev)*bbb.psor[ix,iy,0] \
                            + bbb.erlrc[ix,iy])*frth/com.sy[ixv,iyv]
                if ((ixv>com.ixpt1[ip]) and (ixv<com.ixpt2[ip]+1)):
                    pwr_pfwallh[ixv,ip] = 0
                    pwr_pfwallz[ixv,ip] = 0
            pwr_pfwallz[0,ip] = pwr_pfwallz[1,ip]
            pwr_pfwallz[com.nx+1,ip] = pwr_pfwallz[com.nx,ip]
            pwr_pfwallh[0,ip] = pwr_pfwallh[1,ip]
            pwr_pfwallh[com.nx+1,ip] = pwr_pfwallh[com.nx,ip]

        '''
c ... Finally do the private-flux wall surfaces, i.e., iy=0
      iyv = 0                        # viewing iy position
      do ip = 1, nxpt                # loop over number of x-points
        do ixv = 1, nx   	     # loop over viewing ix
          do iy = 1, ny	             # loop over source iy
            do ix = 1, nx	     # loop over source ix
              theta_ray1 = atan2( zm(ixv+nj,iyv,1)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,1)-rm(ix+nj,iy,0) )
              theta_ray2 = atan2( zm(ixv+nj,iyv,2)-zm(ix+nj,iy,0),
     .                          rm(ixv+nj,iyv,2)-rm(ix+nj,iy,0) )
              dthgx = abs(theta_ray1-theta_ray2)             
              frth = min(dthgx, 2*pi-dthgx)/(2*pi)  # frac; need angle < pi
              pwr_pfwallz(ixv,ip) = pwr_pfwallz(ixv,ip) + 
     .                           prdu(ix,iy)*vol(ix,iy)*frth/sy(ixv,iyv) 
              pwr_pfwallh(ixv,ip) = pwr_pfwallh(ixv,ip) + (
     .                       (eeli(ix,iy)-ebind*ev)*psor(ix,iy,1) +
     .                          erlrc(ix,iy) )*frth/sy(ixv,iyv) 
            enddo
          enddo
          if(ixv>ixpt1(ip) .and. ixv <ixpt2(ip)+1) then  # 0 in core region
            pwr_pfwallh(ixv,ip) = 0.
            pwr_pfwallh(ixv,ip) = 0.
            pwr_pfwallz(ixv,ip) = 0.
            pwr_pfwallz(ixv,ip) = 0.
          endif
        enddo
        pwr_pfwallz(0,ip) = pwr_pfwallz(1,ip)  # prad(0,) is not calculated
        pwr_pfwallz(nx+1,ip) = pwr_pfwallz(nx,ip)
        pwr_pfwallh(0,ip) = pwr_pfwallh(1,ip)
        pwr_pfwallh(nx+1,ip) = pwr_pfwallh(nx,ip)
      enddo

      return
      end
******* end of subroutine pradpltwl *******
        '''




    def balancee(self):
        # balancee as of 22Dec21, now includes molecular energy fluxes to plate/walls
        # sdmin, sdmout, pwallm, and ppfm, as well as impurity ion binding energy
        # heating on plates (sbindzin, sbindzout) (GDP 11 June 2018).
        # As before, also ncludes neutral energy fluxes to plates and walls.
        #
        # added sbindzin, sbindzout, pbindzin, pbindzout	11 June 2018 (GDP)
        # TOTAL POWER AND PARTICLE FLOWS
        from uedge import bbb, com
        from numpy import cos, zeros, sqrt, pi

        def have(x1, x2):
           # The function ave gives the harmonic average of two numbers
           return (2*x1*x2)/(x1+x2)

        ix, iy, id, igsp = 0, 0, 0, 0
        ixdivin=0
        iycore=0
        ixineut=1
        iybegin=1

        ixdivout=com.nx
        ixoneut=com.nx
        iysep=com.iysptrx
        gy = com.gy
        nx = com.nx
        ny = com.ny
        angfx = com.angfx
        visy = bbb.visy
        ixpt1 = com.ixpt1[0]
        ixpt2 = com.ixpt2[0]
        ixm1 = bbb.ixm1
        mi = bbb.mi
        znucl = bbb.znucl
        pwr_wallh = bbb.pwr_wallh
        pwr_wallz = bbb.pwr_wallz
        fngy = bbb.fngy
        pradht = bbb.pradht
        pradrc = bbb.pradrc
        pbinde = bbb.pbinde
        pbindrc = bbb.pbindrc
        znucl = bbb.znucl
        pradfft = bbb.pradfft
        pradimpt = bbb.pradimpt
        pradzbind = bbb.pradzbind
        prdiss = bbb.prdiss
        pibirth = bbb.pibirth
        pvmomcx = bbb.pvmomcx
        ptjdote = bbb.ptjdote
        te = bbb.te
        iycore = 1*(bbb.isguardc == 0)
        cfvisy = bbb.cfvisy
        sy = com.sy
        fniy = bbb.fniy
        feey = bbb.feey
        feiy = bbb.feiy
        up = bbb.up
        fqx = bbb.fqx
        phi0r = bbb.phi0r
        phi0l = bbb.phi0l
        fnix = bbb.fnix
        pwrsore = bbb.pwrsore
        pwrsori = bbb.pwrsori
        volpsor = bbb.volpsor
        nuvl = bbb.nuvl
        vol = com.vol
        ne = bbb.ne
        ni = bbb.ni
        ti = bbb.ti
        upi = bbb.upi
        tg = bbb.tg
        visx = bbb.visx
        ebindz = bbb.ebindz
        ng = bbb.ng
        tg = bbb.tg
        phi = bbb.phi
        fngx = bbb.fngx
        hcxn = bbb.hcxn
        sx = com.sx
        mg = bbb.mg
        zi = bbb.zi
        rrv = com.rrv
        gx = com.gx
        feix = bbb.feix
        feex = bbb.feex
        gxf = com.gxf
        ebind = bbb.ebind
        ev = bbb.ev
        qe = bbb.qe
        l_parloss = bbb.l_parloss
        engbal = bbb.engbal
        nfsp = com.nfsp
        pradpltwl = bbb.pradpltwl
        ishymol = bbb.ishymol
        nusp = com.nusp
        ckinfl = bbb.ckinfl
        nzdf = com.nzdf
        nhsp = com.nhsp
        nzsp = com.nzsp
        isupon = bbb.isupon
        isupgon = bbb.isupgon
        ngsp = com.ngsp
        isimpon = bbb.isimpon
        iion = bbb.iion
        irecomb = bbb.irecomb
        icxgas = bbb.icxgas
        iion_tot = bbb.iion_tot
        irecomb_tot = bbb.irecomb_tot
        icxgas_tot = bbb.icxgas_tot
        fegy = bbb.fegy
        pwr_plth = bbb.pwr_plth
        pwr_pltz = bbb.pwr_pltz


        try:
           ishymoleng = bbb.ishymoleng
        except:
           ishymoleng = 0
        try:
           cenggpl = bbb.cenggpl
        except:
           cenggpl = 0
        try:
           cenggw = bbb.cenggw
        except:
           cenggw = 0


        #two ifs to be able to used old executables
        if(com.islimon == 1) and (com.nyomitmx != 0):
           ixdivout = com.ix_lim-1
           ixdivin = com.ix_lim+1
           iybegin = com.iy_lims
          

        # Determine if molecular hydrogen energy fluxes are present
        if(com.ngsp >= 2) and (bbb.ishymol*bbb.istgon[1] == 1):
           ishymoleng = 1
        else: 
           ishymoleng = 0

        # here we calculate the distance along the inner and outer divertor plates,
        # ydpin and ydpout in meters
        #
        dysepi = 1/(gy[0,0]*cos(angfx[0,0]))
        dysepo = 1/(gy[nx+1,0]*cos(angfx[nx+1,0]))
        for iy in range(1, iysep+1):
           dysepi += 1/(gy[0,iy]*cos(angfx[0,iy]))
           dysepo += 1/(gy[nx,iy]*cos(angfx[nx,iy]))


        ydpin, ydpout = zeros((ny+2,)), zeros((ny+2,))
        ydpin[0]  = -dysepi
        ydpout[0] = -dysepo
        for iy in range(1,ny+2):
           ydpin[iy]  = ydpin[iy-1] + ( 1/gy[0,iy-1] + 1/gy[0   ,iy]) / \
                                       ( cos(angfx[0,iy-1])+cos(angfx[0,iy]) )
           ydpout[iy] = ydpout[iy-1] + ( 1/gy[nx,iy-1] + 1/gy[nx,iy] ) / \
                                       ( cos(angfx[nx,iy-1])+cos(angfx[nx,iy]) )


        # power outflow from separatrix
        psepi, psepe, psepb, pbcorei, pbcoree, pcorebd = 0, 0, 0, 0, 0, 0
        try:
           fluxfacy = fluxfacy
        except:
           fluxfacy=1.




        for ix in range(max(ixpt1+1,0), ixpt2+1):
           ix1 = ixm1[ix,iysep]
           psepi += fluxfacy*( feiy[ix,iysep] + (mi[0]/32)*
                       (up[ix1,iysep  ,0]+up[ix,iysep  ,0]+
                        up[ix1,iysep+1,0]+up[ix,iysep+1,0])**2*fniy[ix,iysep,0] -
               cfvisy*0.125*sy[ix,iysep]*have( visy[ix,iysep,0]*gy[ix,iysep],
                                   visy[ix,iysep+1,0]*gy[ix,iysep+1]) *
                    ( (up[ix1,iysep+1,0]+up[ix,iysep+1,0])**2 -
                      (up[ix1,iysep  ,0]+up[ix,iysep  ,0])**2 ) )
           ix1 = ixm1[ix,iycore]
           psepe += fluxfacy*feey[ix,iysep]
           pbcorei += fluxfacy*( feiy[ix,iycore] + (mi[0]/32)*
                       (up[ix1,iycore  ,0]+up[ix,iycore  ,0]+
                        up[ix1,iycore+1,0]+up[ix,iycore+1,0])**2*fniy[ix,iycore,0] -
                 cfvisy*0.125*sy[ix,iycore]*have( visy[ix,iycore  ,0]*gy[ix,iycore],
                                             visy[ix,iycore+1,0]*gy[ix,iycore+1] ) *
                    ( (up[ix1,iycore+1,0]+up[ix,iycore+1,0])**2 -
                      (up[ix1,iycore  ,0]+up[ix,iycore  ,0])**2 ) )
           pbcoree += fluxfacy*feey[ix,iycore]
           pcorebd += fluxfacy*fniy[ix,iycore,0]*ebind*ev
           psepb += fluxfacy*fniy[ix,iysep,0]*ebind*ev

        #
        # particle outflow from separatrix
        isephyd = qe*sum(fniy[max(ixpt1+1,0):ixpt2+1,iysep,0])


        #       POWER INPUT FROM BIAS SUPPLY
        p_bias = -sum(fqx[nx,1:ny+1]*(phi0r[1:ny+1]-phi0l[1:ny+1]))

        #       TOTAL BIAS CURRENT
        i_bias = sum(fqx[nx,1:ny+1])

        #       ION SATURATION CURRENT
        i_sat_outer = qe*zi[0]*sum(fnix[nx,1:ny+1,0])
        i_sat_inner = qe*zi[0]*sum(fnix[0,1:ny+1,0])

        #       CURRENT AND POWER FROM FIXED VOLUME SOURCES

        p_e_vol = sum(sum(pwrsore))
        p_i_vol = sum(sum(pwrsori))
        i_vol = qe*sum(sum(volpsor[:,:,0]))

        if l_parloss <= 1e9:
          p_e_vol -= sum(nuvl[:,:,0]*vol*ne*bcee*te)
          p_i_vol -= sum(nuvl[:,:,0]*vol*ni[:,:,0]*bcei*ti)
          i_vol -= sum(nuvl[:,:,0]*vol*ni[:,:,0])*qe


        #########################################################################
        # 2-D arrays for checking energy conservation, and for calculation of
        # ionization and radiation sources 

        ptotin = pbcoree+pbcorei ## +p_i_vol+p_e_vol
        ptotin = pbcoree+pbcorei+p_i_vol+p_e_vol

        engbal(ptotin)

        #########################################################################

        # power incident on divertor plate

        pdiviin=0
        pdiviout=0
        pdiveout=0
        pdivein=0.
        pdivmin=0.
        pdivmout=0.
        pbindout=0.
        pbindin=0.
        pdivnout=0.
        pdivnin=0.
        pbindzout=0.
        pbindzin=0.
        pneutout=0.
        pneutin=0.
        pradhout=0.
        pradzout=0.
        pradhin=0.
        pradzin=0.
        sdrout=zeros((ny+2,))
        sdrin=zeros((ny+2,))
        ixi=ixdivin
        ixo=ixdivout
        # note: ixi=ixdivin has been set to 1 to allow velocity derivatives to be calc.
        sdiout=zeros((ny+2,))
        sdeout=zeros((ny+2,))
        sdiin=zeros((ny+2,))
        sdmout=zeros((ny+2,))
        sdmin=zeros((ny+2,))
        sdein=zeros((ny+2,))
        sdtout=zeros((ny+2,))
        sdtin=zeros((ny+2,))
        sbindout=zeros((ny+2,))
        sbindin=zeros((ny+2,))
        sbindzout=zeros((ny+2,))
        sbindzin=zeros((ny+2,))
        sneutout=zeros((ny+2,))
        sneutin=zeros((ny+2,))
        sdioutd=zeros((ny+2,nfsp))
        sdiind=zeros((ny+2,nfsp))
        sdnout=zeros((ny+2,))
        sdnin=zeros((ny+2,))

        # Allow use of "old" or "new" switches for neutral energy loss


        # First get radiation power in pwr_plth and pwr_pltz
        pradpltwl()
        sdrin = pwr_plth[:,0]+pwr_pltz[:,0]
        sdrout = pwr_plth[:,1]+pwr_pltz[:,1]


        # here the sds are ion and electron poloidal power fluxes in W/m**2
        for iy in range(iybegin, ny+1):
           sxo = sx[ixo,iy]/(cos(angfx[ixo,iy]))
           sxi = sx[ixi,iy]/(cos(angfx[ixi,iy]))
           vxno =  0.25*sqrt(8*tg[ixoneut,iy,0]/(pi*mg[0]))
           vxni =  0.25*sqrt(8*tg[ixineut,iy,0]/(pi*mg[0]))

           for id in range(nfsp): # note: upi=0 for the netural species
              if zi[id] > 0:
                 sdioutd[iy,id] = ( 0.5*mi[id]*upi[ixo,iy,id]**2*fnix[ixo,iy,id] )/sxo
                 sdiind[iy,id] = ( -0.5*mi[id]*upi[ixi,iy,id]**2*fnix[ixi,iy,id] )/sxi
              else:  #zero-out neutrals as fnix will be into plasma 5/27/08 - TDR
                 if ishymol == 0: # recompute parallel fnix
                    sdioutd[iy,id] = 0.5*mi[id]*abs(up[ixo,iy,id])**3*ni[ixo,iy,id]*rrv[ixo,iy]
                    sdiind[iy,id] = 0.5*mi[id]*abs(up[ixi,iy,id])**3*ni[ixi,iy,id]*rrv[ixi,iy]
                 else:   # for molecules, fnix should be ok
                    sdiout[iy,id]=0#1e-20*( -0.5*mi[id]*up[ixo,iy,id]**2*fnix[ixo,iy,id] )/sxo
                    sdiind[iy,id]=0#1e-20*( -0.5*mi[id]*up[ixi,iy,id]**2*fnix[ixi,iy,id] )/sxi
              sdiout[iy] += sdioutd[iy,id]
              sdiin[iy]  += sdiind[iy,id]

           # AH: ckinfl treated as array, although it is double. Legacy switch?
           for id in range(nusp):      # note: up for the netural species in nonzero
              tempvo =  - ckinfl*0.5*sx[ixo,iy]*visx[ixo,iy,id]*gx[ixo,iy]*\
                            ( up[ixo,iy,id]**2 - up[ixo-1,iy,id]**2 ) /sxo
              tempvi =  + ckinfl*0.5*sx[ixi,iy]*visx[ixi+1,iy,id]*gx[ixi+1,iy]*\
                           ( up[ixi+1,iy,id]**2 - up[ixi,iy,id]**2 ) / sxi
              sdioutd[iy,id] += tempvo
              sdiout[iy] += tempvo
              sdiind[iy,id] += tempvi
              sdiin[iy] += tempvi

           sdiout[iy] += feix[ixo,iy]/sxo
           sbindout[iy] = fnix[ixo,iy,0] * ebind*ev / sxo
           if ishymoleng==1:  #mol heat flux; drift eng small,<0
             sdmout[iy] += fegx[ixo,iy,1]/sxo

        #  Compute binding-energy energy fluxes for impurities
           for id in range(nzdf):
              if (id == 0):
                 id2min = nhsp
                 id2max = id2min +nzsp[id]-1
              else:
                 id2min = nhsp+sum(nzsp[1:id-1])
                 id2max = id2min + nzsp[id]-1
              for id2 in range(id2min, id2max+1):
                 for id3 in range(znucl[id]):
                    sbindzout[iy] += fnix[ixo,iy,id2]*ebindz(id3,znucl[id2])*ev/sxo
                    sbindzin[iy] -= fnix[ixi,iy,id2]*ebindz(id3,znucl[id2])*ev/sxo

           sneutout[iy] = cenggpl*2.*vxno*ng[ixoneut,iy,0]*tg[ixoneut,iy,0]
           sdeout[iy] = ( feex[ixo,iy]+fqx[ixo,iy]*(phi[ixo,iy]-phi0r[iy]) )/sxo
           sdtout[iy] = sdeout[iy] + sdiout[iy] + sbindout[iy] + sdmout[iy] + \
                           sbindzout[iy] + pwr_plth[iy,1] + pwr_pltz[iy,1]
           sdiin[iy]  -= feix[ixi,iy]/sxi 
           if ishymoleng==1: #mol heat flux; drift eng small,<0
             sdmin[iy] -= fegx[ixi,iy,1]/sxo
           
           sbindin[iy] = - fnix[ixi,iy,1] * ebind*ev / sxi
           sneutin[iy] = cenggpl*2.*vxni*ng[ixineut,iy,0]*tg[ixineut,iy,0]
           sdein[iy]  = -( feex[ixi,iy] + .001*fqx[ixi,iy]*(phi[ixi,iy]-phi0l[iy]) )/sxi
           sdtin[iy] = sdein[iy] + sdiin[iy] + sdmin[iy] + sbindin[iy] + \
                          sbindzin[iy] + pwr_plth[iy,0] + pwr_pltz[iy,0]
           pdiviout += sdiout[iy]*sxo ## + sdioutd[iy,1]*sxo
           pdiveout += sdeout[iy]*sxo
           pdivmout += sdmout[iy]*sxo
           pdiviin  += sdiin[iy]*sxi ## + sdiind[iy,1]*sxi
           pdivein  += sdein[iy]*sxi
           pdivmin  += sdmin[iy]*sxi
           pbindout += sbindout[iy]*sxo
           pbindzout += sbindzout[iy]*sxo
           pbindin += sbindin[iy]*sxi
           pbindzin += sbindzin[iy]*sxi
           pneutout += 0*sneutout[iy]*sxo  # included in pdiviout
           pneutin += 0.*sneutin[iy]*sxo     # included in pdiviin
           pradhout += pwr_plth[iy,1]*sxo
           pradzout += pwr_pltz[iy,1]*sxo
           pradhin += pwr_plth[iy,0]*sxi
           pradzin += pwr_pltz[iy,0]*sxi
           if isupgon[0] == 1: # Approx. neutral energy flux
              sdnout[iy]=sdioutd[iy,1] + ( sx[ixo,iy]*hcxn[ixo,iy]*gxf[ixo,iy]* \
                                                    (ti[ixo,iy]-ti[ixo+1,iy]) + \
                                          2.5*fnix[ixo,iy,1]*ti[ixo+1,iy] ) / sxo
              sdnin[iy]=sdiind[iy,1] + ( sx[ixi,iy]*hcxn[ixi,iy]*gxf[ixi,iy]* \
                                                  (ti[ixi,iy]-ti[ixi+1,iy]) + \
                                          2.5*fnix[ixi,iy,1]*ti[ixi,iy] ) / sxi
              pdivnout += sdnout[iy]*sxo
              pdivnin += sdnin[iy]*sxi

        # 
        # added, 11Jun2018 GDP - only diagnostic; not used below
        ptot=pdiviout+pdiveout+pdiviin+pdivein+pbindout+pbindin+\
                   pbindzin+pbindzout+pdivmout+pdivmin+pneutout+pneutin
        #

        # Fix up boundary values
        sdtin[0] = sdtin[1]
        sdtin[ny+1] = sdtin[ny]
        sdein[0] = sdein[1]
        sdein[ny+1] = sdein[ny]
        sdiin[0] = sdiin[1]
        sdiin[ny+1] = sdiin[ny]
        sdtout[0] = sdtout[1]
        sdtout[ny+1] = sdtout[ny]
        sdeout[0] = sdeout[1]
        sdeout[ny+1] = sdeout[ny]
        sdiout[0] = sdiout[1]
        sdiout[ny+1] = sdiout[ny]


        #
        ptotpartin = pdiviin+pdivein+pbindin+pbindzin+pdivmin  ##pneutin
        ptotpartout = pdiviout+pdiveout+pbindout+pbindzout+pdivmout  ##pneutout
        ptotpart = ptotpartin+ptotpartout  ##pneutout+pneutin
        ptotrad = pradhout+pradzout+pradhin+pradzin
        #
        # ion current to divertor plate
        idivout = zeros((nfsp,))
        idivin= zeros((nfsp,))
        igasout= zeros((ngsp,))
        igasin= zeros((ngsp,))

        for iy in range(iybegin,ny+1):
           for id in range(nfsp):
              idivout[id] += fnix[ixdivout,iy,id]
              idivin[id] += fnix[ixdivin,iy,id]
           for igsp in range(ngsp):
              if isupgon[igsp] == 0:
                 igasout[igsp] += fngx[ixdivout,iy,igsp]
                 igasin[igsp] += fngx[ixdivin,iy,igsp]
              else:
                 igasout[igsp] += fnix[ixdivout,iy,1]
                 igasin[igsp] += fnix[ixdivin,iy,1]

        idivout *= qe
        idivin *= qe
        igasout *= qe
        igasin *= qe
        #
        # ion current to the core
        icore = zeros((nfsp,))
        for ix in range(max(0, ixpt1+1), ixpt2+1):
           for id in range(nfsp):
              icore[id] += fluxfacy*fniy[ix,0,id]
           
        icore *= qe
        #
        iywall=ny       # DEFINITION
        iypf=0
        #
        # power flow to vessel and private flux wall
        pwalli = 0.
        ppfi = 0.
        pwalle = 0.
        ppfe = 0.
        pwallm = 0.
        ppfm = 0.
        pwallbd = 0.
        ppfbd = 0.
        pradhwall = 0.
        pradzwall = 0.
        pradhpf = 0.
        pradzpf = 0.
        sneutpf = zeros((nx+2))
        sneutw = zeros((nx+2))
        pneutpf = 0.
        pneutw = 0.

        for ix in range(nx+2):
           ix1 = ixm1[ix,iywall]
           vynw = 0.25*sqrt(8*tg[ix,ny+1,0]/(pi*mg[0]))
           pwalli += fluxfacy*( feiy[ix,iywall] +
                    0.125*mi[0]*(up[ix1,iywall,0]+up[ix,iywall,0])**2*fniy[ix,iywall,0] -
                    cfvisy*0.125*sy[ix,iywall]*have( visy[ix,iywall,0]*gy[ix,iywall],
                                             visy[ix,iywall+1,0]*gy[ix,iywall+1] ) *
                    ( (up[ix1,iywall+1,0]+up[ix,iywall+1,0])**2 -
                      (up[ix1,iywall  ,0]+up[ix,iywall  ,0])**2 ) )
           pwalle += fluxfacy*feey[ix,iywall]
           pwallbd += fluxfacy*fniy[ix,iywall,0]*ebind*ev
           pradhwall += fluxfacy*pwr_wallh[ix]*sy[ix,iywall]
           pradzwall += fluxfacy*pwr_wallz[ix]*sy[ix,iywall]
           sneutw[ix] = cenggw*2.*vynw*tg[ix,ny+1,0]*sx[ix,ny]
           pneutw += sneutw[ix]
           if ishymoleng == 1:  #molec temp eqn active
              pwallm += fegy[ix,iywall,1]


        for ix in range(0, ixpt1+1):
           ix1 = ixm1[ix,iypf]
           vynpf = 0.25*sqrt(8*tg[ix,0,1]/(pi*mg[0]))
           ppfi -= fluxfacy*(feiy[ix,iypf] -
                    0.125*mi[0]*(up[ix1,iypf,0]+up[ix,iypf,0])**2*fniy[ix,iypf,0] +
                    cfvisy*0.125*sy[ix,iypf]*have( visy[ix,iypf,0]*gy[ix,iypf ],
                                           visy[ix,iypf+1,0]*gy[ix,iypf+1] ) *
                    ( (up[ix1,iypf+1,0]+up[ix,iypf+1,0])**2 -
                      (up[ix1,iypf  ,0]+up[ix,iypf  ,0])**2 ) )
           ppfe -= fluxfacy*feey[ix,iypf]
           ppfbd -= fluxfacy*fniy[ix,iypf,0]*ebind*ev
        ##   pradhpf += fluxfacy*pwr_pfh[ix]*sy[ix,iywall]
        ##   pradzpf += fluxfacy*pwr_pfz[ix]*sy[ix,iywall]
           sneutpf[ix] = cenggw*2.*vynw*tg[ix,0,1]*sx[ix,0]
           pneutpf += sneutpf[ix]
           if ishymoleng==1:  #molec temp eqn active
             ppfm -= fegy[ix,iypf,1]
           

        for ix in range(ixpt2+1, nx+2):
           ix1 = ixm1[ix,iypf]
           vynpf = 0.25*sqrt(8*tg[0,iy,0]/(pi*mg[0]))
           ppfi -= fluxfacy*( feiy[ix,iypf] -
                    0.125*mi[0]*(up[ix1,iypf,0]+up[ix,iypf, 0])**2*fniy[ix,iypf,0] +
                    cfvisy*0.125*sy[ix,iypf]*have( visy[ix,iypf,0]*gy[ix,iypf],
                                           visy[ix,iypf+1,0]*gy[ix,iypf+1] ) *
                    ( (up[ix1,iypf+1,0]+up[ix,iypf+1,0])**2 -
                      (up[ix1,iypf  ,0]+up[ix,iypf  ,0])**2 ) )
           ppfe -= fluxfacy*feey[ix,iypf]
           ppfbd -= fluxfacy*fniy[ix,iypf,0] * ebind*ev
           if ishymoleng==1: #molec temp eqn active
             ppfm -= fegy[ix,iypf,1]
           
        #
        # ion current to vessel wallst.
        iwall = zeros((nfsp,))
        ipf = zeros((nfsp,))
        igaswall = zeros((nfsp,))
        igaspf = zeros((nfsp,))
        igascr = zeros((nfsp,))

        for ix in range(0, nx+2):
           for id in range(nfsp):
              iwall[id] += fluxfacy*fniy[ix,iywall,id]
           for igsp in range(ngsp):
              igaswall[igsp] += fluxfacy*fngy[ix,iywall,igsp]

        for ix in range(0, ixpt1+1):
           for id in range(0, nfsp):
              ipf[id] += fluxfacy*fniy[ix,iypf,id]
           for igsp in range(ngsp):
              igaspf[igsp] += fluxfacy*fngy[ix,iypf,igsp]

        for ix in range(ixpt2+1, nx+2):
           for id in range(0, nfsp):
              ipf[id] += fluxfacy*fniy[ix,iypf,id]
           for igsp in range(ngsp):
              igaspf[igsp] += fluxfacy*fngy[ix,iypf,igsp]

        for ix in range(max(ixpt1+1,0), ixpt2+1):
           for igsp in range(ngsp):
              igascr[igsp] =  fluxfacy*fngy[ix,iypf,igsp]
           

        iwall *= qe
        igaswall *= qe
        ipf *= qe
        igaspf *= qe
        igascr *= qe
        #
        print(" ")
        print("Power Flow [Watts] from Core to Scrape-off Layer is:")
        print("   Total: {:.2e}".format(pbcorei+pbcoree+pcorebd))
        print("   Ion contribution: {:.2e}".format(pbcorei))
        print("   Electron contribution: {:.2e}".format(pbcorei))
        print("   Binding energy contribution: {:.2e}".format(pbcorei))
        print()
        print("Power Flow [Watts] over the separatrix to Scrape-off Layer is:")
        print("Power Flow [Watts] over the separatrix to Scrape-off Layer is:")
        print("   Total: {:.2e}".format(psepi+psepe+psepb))
        print("   Ion contribution: {:.2e}".format(psepi))
        print("   Electron contribution: {:.2e}".format(psepe))
        print("   Binding energy contribution: {:.2e}".format(psepb))
        #
        print(" ")
        print("Power Input [Watts] from Fixed Volume Sources (Sinks):")
        print(p_i_vol)
        print(p_e_vol)
        #
        print(" ")
        print("Power Flows [Watts] incident on Divertor Plates are:")
        print("Inner plate:")
        print("   Total: {:.2e}".format(pdiviin+pdivein+pdivmin,pbindin+pneutin+pradhin+pradzin))
        print("   Ion contribution: {:.2e}".format(pdiviin))
        print("   Electron contribution: {:.2e}".format(pdivein))
        print("   Molecular contribution: {:.2e}".format(pdivmin))
        print("   Neutral atom contribution: {:.2e}".format(pneutin))
        print("   Binding energy contribution: {:.2e}".format(pbindin))
        print("   Hydrogenic radiation contribution: {:.2e}".format(pradhin))
        print("   Impurity radiation contribution: {:.2e}".format(pradzin))

        print("Outer plate;")
        print("   Total: {:.2e}".format(pdiviout+pdiveout+pdivmout,pbindout+pneutout+pradhout+pradzout))
        print("   Ion contribution: {:.2e}".format(pdiviout))
        print("   Electron contribution: {:.2e}".format(pdiveout))
        print("   Molecular contribution: {:.2e}".format(pdivmout))
        print("   Neutral atom contribution: {:.2e}".format(pneutout))
        print("   Binding energy contribution: {:.2e}".format(pbindout))
        print("   Hydrogenic radiation contribution: {:.2e}".format(pradhout))
        print("   Impurity radiation contribution: {:.2e}".format(pradzout))

        print("Totals for both plates:")
        print("   Total: {:.2e}".format(ptotpart+ptotrad))
        print("   Particle fluxes: {:.2e}".format(ptotpart))
        print("   Radiation: {:.2e}".format(ptotrad))

        #
        if isupgon[0] == 1:
           print(" ")
           print("Est. of pdiviout and pdivin from atoms (backscatter) [Watts]:")
           print("   Outer: {:.2e}".format(pdivnout))
           print("   Inner: {:.2e}".format(pdivnin))

        #
        print(" ")
        print("Power [W] lost via ionization & recombination radiation is:")
        print("{:.2e}".format(pradht))
        #
        print(" ")
        print("Power [W] lost via recombination only (included in pradht above):")
        print("{:.2e}".format(pradrc))
        #
        print(" ")
        print("Power [W] lost at ionization but carried as ion binding-energy:")
        print("{:.2e}".format(pbinde))
        #
        print(" ")
        print("Power [W] gained by electrons in 3-body recombination (via binding eng):")
        print("{:.2e}".format(pbindrc))
        #
        if isimpon != 0:
           print(" ")
           print("Power [W] lost via impurity radiation is:")
           id2=nhsp-1
           try:
              pradfft = pradfft
           except:
              pradfft = 0
           for id in range(nzdf):
              id2 += nzsp[id]
              print("  for nuclear charge = {};  Power = {}".format(znucl[id2], pradimpt[id]))
              print("  for fixed-fraction species;  Power = ", pradfft)
           #
           print(" ")
           print("Electron Power [W] lost at impur. ioniz. but carried as binding-energy:")
           print("{:.2e}".format(pradzbind))
        #
        print(" ")
        print("Power [W] lost via dissociation of molecules is:")
        print("{:.2e}".format(prdiss))
        #
        print(" ")
        print("Power [W] gained by ions from initial Franck-Condon Energy:")
        print("{:.2e}".format(pibirth))
        #

        print(" ")
        print("Power [W] lost in parallel momentum exhange via charge exchange:")
        print("{:.2e}".format(pvmomcx))
        #
        print(" ")
        print("Power [W] from J.E Joule heating - goes to electrons:")
        print("{:.2e}".format(ptjdote))
        #
        print(" ")
        print("Power Flow [Watts] incident on Vessel Wall is:")
        ###print(pwalli,pwalle,pwallm,pwallbd,pradhwall,pradzwall,pneutw)
        print("   Outerwall_sum: {:.2e}".format(pwalli+pwalle+pwallm+pwallbd+pradhwall+pradzwall+pneutw))
        print("      pwalli: {:.2e}".format(pwalli))
        print("      pwalle: {:.2e}".format(pwalle))
        print("      pwallm: {:.2e}".format(pwallm))
        print("      pwallbd: {:.2e}".format(pwallbd))
        print("      pradhwall: {:.2e}".format(pradhwall))
        print("      pradzwall: {:.2e}".format(pradzwall))
        print("      pneutw: {:.2e}".format(pneutw))


        print(" ")
        print("Power Flow [Watts] incident on Private Flux Wall is:")


        ###print("Total power: {:.2e}" .format(ppfi+ppfe+ppfm+ppfbd+pneutpf+ppfi+ppfe)
        print("   PFwall_sum: {:.2e}".format(ppfi+ppfe+ppfm+ppfbd+pneutpf))
        print("      ppfi: {:.2e}".format(ppfi))
        print("      ppfe: {:.2e}".format(ppfe))
        print("      ppfm: {:.2e}".format(ppfm))
        print("      ppfbd: {:.2e}".format(ppfbd))
        print("      pneutpf: {:.2e}".format(pneutpf))

        # Calculate power into plasma volume for normalizing factor in power balance
        pnormpb = 0.
        if ptotpartin > 0.:
            pnormpb = pnormpb + ptotpartin
        if ptotpartout < 0.:
            pnormpb = pnormpb - ptotpartout
        if pbcoree+pbcorei+pcorebd+p_i_vol > 0.:
            pnormpb = pnormpb+pbcoree+pbcorei+pcorebd+p_i_vol
        if p_i_vol+p_e_vol > 0.:
            pnormpb = pnormpb+p_i_vol+p_e_vol

        # Here particle power and radiation power are separate terms
        powbal = (pbcoree+pbcorei+pcorebd+p_i_vol+p_e_vol+ptjdote - ptotpart -
                  pwalli-pwalle-pwallm-pwallbd-ppfi-ppfe-ppfm-ppfbd-sum(pradimpt)-
                  pradfft-pradht-prdiss-pvmomcx+pibirth)/ pnormpb
        if abs(pdivein+pdiviin) > 10.*abs(ptotin+ptjdote):   #for cases with no radial pwr input
           powbal = powbal*abs(ptotin+ptjdote)/abs(pdivein+pdiviin)
           
        ##print(" ")
        ##print("Total power out")
        ##print( (pbcoree+pbcorei+pcorebd+p_i_vol+p_e_vol+ptjdote - ptotpart -
        ##          pwalli-pwalle-pwallm-pwallbd-ppfi-ppfe-ppfm-ppfbd-sum(pradimpt)-
        ##          pradfft-pradht-prdiss-pvmomcx+pibirth))

        print(" ")
        print("Power Balance: (Pin-Pout)/Pin")
        print("    {:.2e}".format(powbal))
        #
        print("============================================================== ")
        print(" ")
        print("Particle Flow [Amps] from Core to Scrape-off Layer is:")
        for id in range(nfsp):
           if zi[id] > 0:
              print("  for nuclear charge = {}; ion charge = {}".format(znucl[id],zi[id]))
              print("  icore = ", icore[id])
              print(" ")

        for id in range(ngsp):
           print("  for gas species = {}; igascr = {:.2e}".format(id, igascr[id]))
        print("  hydrogen ion current at separatrix, isephyd = ", isephyd)
        #
        print(" ")
        print("Current from Fixed Volume Source (Sinks) [Amps]:")
        print("  ion current, i_vol = ", i_vol)
        #
        print(" ")
        print("Ionization current [Amps] created by re-ionization of gas is:")

        for id in range(ngsp):
           print("  for gas species = {}; iion = {:.2e}".format(id, iion[id]))
        print(" ")
        print("Recomb. current [Amps] from ions --> gas by electron-ion recombination:")

        for id in range(ngsp):
           print("  for gas species = {}; irecomb = {:.2e}".format(id, irecomb[id]))

        print(" ")
        print("Charge exchange current [Amps] from ions --> gas by iter-species cx:")

        for id in range(ngsp):
           print("  for gas species = {}; icxgas = {:.2e}".format(id, icxgas[id]))

        print(" ")
        print("Particle Flow [Amps] incident on Divertor Plate is:")

        for id in range(nfsp):
           if zi[id] > 0:
              print("  for nuclear charge = {}; ion charge = {}".format(znucl[id], zi[id]))
              print("  idivin = {:.2e}   idivout = {:.2e}".format(idivin[id],idivout[id]))
              print(" ")
         

        print("Neutral Flow [Amps] away from Divertor Plate is:")

        for id in range(ngsp):
           print("  for gas species = {}; igasin = {:.2e},  igasout = {:.2e}".format(id,  igasin[id], igasout[id]))


        #
        print(" ")
        print("Particle Flow [Amps] incident on Vessel Wall is:")
        for id in range(nfsp):
           print("  for ion species = {}; iwall = {:.2e}".format(id, iwall[id]))

        for id in range(ngsp):
           print("  for gas species = {}; igaswall = {:.2e}".format(id, igaswall[id]))


        print(" ")
        print("Particle Flow [Amps] incident on Private Flux Wall is:")
        for id in range(nfsp):
           print("  for ion species = {}; ipf = {:.2e}".format(id, ipf[id]))

        for id in range(ngsp):
           print("  for gas species = {}; igaspf = {:.2e}".format(id, igaspf[id]))


        #

        igastot = 0.
        igasdenom = 0.
        for id in range(ngsp):
           igastot += igasin[id]-igasout[id]-igaswall[id]+igaspf[id]\
                             + igascr[id] + i_vol
           igasdenom += igasin[id] -igasout[id]

        if ishymol == 1:  # add id=2 case again because of 2 atoms/molecule
           igastot += igasin[1]-igasout[1]-igaswall[1]+igaspf[1]\
                             + igascr[1]
           igasdenom += (igasin[1] -igasout[1])


        deligas = (igastot - iion_tot - irecomb_tot - icxgas_tot) / igasdenom
        print(" ")
        print("Particle Balance for Neutrals: Itotal/Iinput")
        print("    {:.2e}".format(deligas))
        #
        print(" ")
        print("Peak electron temperatures [eV] on divertor plates")
        print("    Outer: {:.2e}".format(max(te[nx+1,1:ny+1])/ev))
        print("    Inner: {:.2e}".format(max(te[0,1:ny+1])/ev))
        #
        print(" ")
        print("Peak ion temperatures [eV] on divertor plates")
        print("    Outer: {:.2e}".format(max(ti[nx+1,1:ny+1])/ev))
        print("    Inner: {:.2e}".format(max(ti[0,1:ny+1])/ev))
        #
        print(" ")
        print("Peak power flux [MW/m**2] on divertor plates")
        print("    Outer: {:.2e}".format(1.e-6*max(sdtout[1:ny+1])))
        print("    Inner: {:.2e}".format(1.e-6*max(sdtin[1:ny+1])))
        #
        print(" ")
        print("Peak ion densities [m**(-3)] on divertor plates")
        print("    Outer: {:.2e}".format(max(ni[nx+1,1:ny+1,0])))
        print("    Inner: {:.2e}".format(max(ni[0,1:ny+1,0])))
        #
        print(" ")
        print("Peak gas densities [m**(-3)] on divertor plates")
          # note: guard cell gas density often meaningless, so use one cell in
        print("    Outer: {:.2e}".format(max(ng[nx,1:ny+1,0])))
        print("    Inner: {:.2e}".format(max(ng[1,1:ny+1,0])))




        '''

c-------------------------------------------------------------------------c

      subroutine plateflux

      Implicit none

c ... Calc particle & heat fluxes to divertor plates
c ... Use as diagnostic called from BASIS/Python parser (alt. to balancee)

      Use(Dim)            # nx,ny
      Use(Comgeo)         # sx,vol,gx,gxf
      Use(Noggeo)         # angfx
      Use(Phyvar)         # pi, ev
      Use(Postproc)       # pwr_plth,Pwr_pltz,sdel,rb;sdil,rb; sdbindl,rb,
                          # sdtl,rb;gdil,rd
      Use(Xpoint_indices) # ixlb, ixrb
      Use(Compla)         # mi,ni,up,te,ti
      Use(Comflo)         # fnix,feex,feix
      Use(Conduc)         # visx,hcxn
      Use(UEpar)          # ebind
      Use(Bcond)          # ckinfl
      Use(Parallv)        # nxg,nyg
      Use(Poten)          # phi0l,phi0r

c  Local variables
      integer iu,jx,id,ixi,ixo 
      #Former Aux module variables
      integer ix,iy
      real tempvo,tempvi
      real pdivilb(1:nxpt),pdivirb(1:nxpt),pdivelb(1:nxpt),
     .     pdiverb(1:nxpt)
      real pbindlb(1:nxpt),pbindrb(1:nxpt),pdivnlb(1:nxpt),
     .     pdivnrb(1:nxpt)
      real sdilbd(0:ny+1,1:nfsp,1:nxpt),sdirbd(0:ny+1,1:nfsp,1:nxpt)
      real sdnlb(0:ny+1,1:nxpt),sdnrb(0:ny+1,1:nxpt)
      real sxi(0:ny+1,1:nxpt),sxo(0:ny+1,1:nxpt)

########################################
# First do the particle flux
##############################################################
      do jx=1,nxpt
        ixi=ixlb(jx)	# ixi=0
        ixo=ixrb(jx)	# ixo=nx
        do iy=1,ny+1
          sxo(iy,jx) = sx(ixo,iy)/(cos(angfx(ixo,iy)))
          sxi(iy,jx) = sx(ixi,iy)/(cos(angfx(ixi,iy)))
          do id = 1, nfsp
	    gdilb(iy,id,jx) = -fnix(ixi,iy,id)/sxi(iy,jx)
	    gdirb(iy,id,jx) =  fnix(ixo,iy,id)/sxo(iy,jx)
            engilb(iy,id,jx) = (2.*ti(ixi,iy)/ev + zi(id)*
     .                          (phi(ixi,iy)-phi0l(iy,jx)) )
            engirb(iy,id,jx) = 2.*ti(ixo+1,iy)/ev + zi(id)*
     .                          (phi(ixo+1,iy)-phi0r(iy,jx))
          enddo
        enddo
      enddo

# Fix corner boundary values
      do jx = 1, nxpt
        do id = 1, nfsp
          gdilb(0,id,jx) =     gdilb(1,id,jx)
          gdilb(ny+1,id,jx) =  gdilb(ny,id,jx)
          gdirb(0,id,jx) =     gdirb(1,id,jx)
          gdirb(ny+1,id,jx) =  gdirb(ny,id,jx)
          engilb(0,id,jx) =    engilb(1,id,jx)
          engilb(ny+1,id,jx) = engilb(ny,id,jx)
          engirb(0,id,jx) =    engirb(1,id,jx)
          engirb(ny+1,id,jx) = engirb(ny,id,jx)
        enddo
      enddo

#######################################
# Now do the heat flux
##############################################################
#

# Get radiation power in pwr_plth and pwr_plt
      call pradpltwl

# Fill radiation heat flux arrays and initialize tot powers
      do jx = 1, nxpt
        pdivirb(jx) = 0.
        pdiverb(jx) = 0. 
        pdivilb(jx) = 0. 
        pdivelb(jx) = 0. 
        pbindrb(jx) = 0. 
        pbindlb(jx) = 0. 
        do iy = 0, ny+1
          iu = 2*(jx/2)+1  # iu/iu+1 gives "odd/even" plates (in/out)
          sdrlb(iy,jx) = pwr_plth(iy,iu)+pwr_pltz(iy,iu)
          sdrrb(iy,jx) = pwr_plth(iy,iu+1)+pwr_pltz(iy,iu+1)
        enddo
      enddo

# here the sds are ion and electron poloidal power fluxes in W/m**2
      do jx=1,nxpt
        ixi=ixlb(jx)	# ixi=0
        ixo=ixrb(jx)	# ixo=nx
        do iy=1,ny+1
          sdirb(iy,jx) = 0.
          sdilb(iy,jx) = 0.
          do id = 1, nfsp
           if (zi(id) .gt. 0) then
	    sdirbd(iy,id,jx) = ( 0.5*mi(id)*upi(ixo,iy,id)**2*
     .                           fnix(ixo,iy,id) )/sxo(iy,jx)
	    sdilbd(iy,id,jx) = (-0.5*mi(id)*upi(ixi,iy,id)**2*
     .                           fnix(ixi,iy,id) )/sxi(iy,jx)
           else    # note: upi=0 for neutral species; use up instead
	    sdirbd(iy,id,jx) = ( 0.5*mi(id)*up(ixo,iy,id)**2*
     .                           fnix(ixo,iy,id) )/sxo(iy,jx)
	    sdilbd(iy,id,jx) = (-0.5*mi(id)*up(ixi,iy,id)**2*
     .                           fnix(ixi,iy,id) )/sxi(iy,jx)
           endif
           sdirb(iy,jx) = sdirb(iy,jx) + sdirbd(iy,id,jx)
           sdilb(iy,jx) = sdilb(iy,jx) + sdilbd(iy,id,jx)
         enddo
         do id = 1, nusp      # note: up neutral species in nonzero
           tempvo =  - ckinfl*0.5*sx(ixo,iy)*visx(ixo,iy,id)*gx(ixo,iy)*
     .               ( up(ixo,iy,id)**2 - up(ixo-1,iy,id)**2 ) /sxo(iy,jx)
	   tempvi =  + ckinfl*0.5*sx(ixi,iy)*visx(ixi+1,iy,id)*gx(ixi+1,iy)*
     .                ( up(ixi+1,iy,id)**2 - up(ixi,iy,id)**2 ) / sxi(iy,jx)
	   sdirbd(iy,id,jx) = sdirbd(iy,id,jx) + tempvo
	   sdirb(iy,jx) = sdirb(iy,jx) + tempvo
	   sdilbd(iy,id,jx) = sdilbd(iy,id,jx) + tempvi
	   sdilb(iy,jx) = sdilb(iy,jx) + tempvi
         enddo
         sdirb(iy,jx) = sdirb(iy,jx) + feix(ixo,iy)/sxo(iy,jx)
         sbindrb(iy,jx) =  fnix(ixo,iy,1) * ebind*ev / sxo(iy,jx)
         sderb(iy,jx) =  ( feex(ixo,iy)+fqx(ixo,iy)*
     .                     (phi(ixo+1,iy)-phi0r(iy,jx)) )/sxo(iy,jx)
         sdtrb(iy,jx) = sderb(iy,jx) + sdirb(iy,jx) + sbindrb(iy,jx)
         sdilb(iy,jx) = sdilb(iy,jx) - feix(ixi,iy)/sxi(iy,jx) 
         sbindlb(iy,jx) = -fnix(ixi,iy,1) * ebind*ev / sxi(iy,jx)
         sdelb(iy,jx) = -( feex(ixi,iy)+fqx(ixi,iy)*
     .                     (phi(ixi  ,iy)-phi0l(iy,jx)) )/sxi(iy,jx)
         sdtlb(iy,jx) = sdelb(iy,jx) + sdilb(iy,jx) + sbindlb(iy,jx)
         pdivirb(jx) = pdivirb(jx) + sdirb(iy,jx)*sxo(iy,jx)
         pdiverb(jx) = pdiverb(jx) + sderb(iy,jx)*sxo(iy,jx)
         pdivilb(jx) = pdivilb(jx) + sdilb(iy,jx)*sxi(iy,jx)
         pdivelb(jx) = pdivelb(jx) + sdelb(iy,jx)*sxi(iy,jx)
         pbindrb(jx) = pbindrb(jx) + sbindrb(iy,jx)*sxo(iy,jx)
         pbindlb(jx) = pbindlb(jx) + sbindlb(iy,jx)*sxi(iy,jx)
         if (isupgon(1).eq.1) then    # Approx. neutral energy flux
	   sdnrb(iy,jx)=sdirbd(iy,jx,2) + ( sx(ixo,iy)*hcxn(ixo,iy)*
     .                    gxf(ixo,iy)*(ti(ixo,iy)-ti(ixo+1,iy)) + 
     .   		  2.5*fnix(ixo,iy,2)*ti(ixo+1,iy) ) / sxo(iy,jx)
	   sdnlb(iy,jx)=sdilbd(iy,jx,2) - ( sx(ixi,iy)*hcxn(ixi,iy)*
     .                    gxf(ixi,iy)*(ti(ixi,iy)-ti(ixi+1,iy)) + 
     .                    2.5*fnix(ixi,iy,2)*ti(ixi  ,iy) ) / sxi(iy,jx)
	   pdivnrb(jx) = pdivnrb(jx) + sdnrb(iy,jx)*sxo(iy,jx)
	   pdivnlb(jx) = pdivnlb(jx) + sdnlb(iy,jx)*sxi(iy,jx)
        endif
        enddo  # end do-loop for iy=1,ny+1
      enddo  # end do-loop for jx=1,nxpt

# Fix corner boundary values
      do jx = 1, nxpt
         sdtlb(0,jx) = sdtlb(1,jx)
         sdtlb(ny+1,jx) = sdtlb(ny,jx)
         sdelb(0,jx) = sdelb(1,jx)
         sdelb(ny+1,jx) = sdelb(ny,jx)
         sdilb(0,jx) = sdilb(1,jx)
         sdilb(ny+1,jx) = sdilb(ny,jx)
         sdrlb(ny+1,jx) = sdrlb(ny,jx)
         sdtrb(0,jx) = sdtrb(1,jx)
         sdtrb(ny+1,jx) = sdtrb(ny,jx)
         sderb(0,jx) = sderb(1,jx)
         sderb(ny+1,jx) = sderb(ny,jx)
         sdirb(0,jx) = sdirb(1,jx)
         sdirb(ny+1,jx) = sdirb(ny,jx)
         sdrrb(ny+1,jx) = sdrrb(ny,jx)
      enddo

      return
      end
******* end of subroutine plateflux *******
c-------------------------------------------------------------------------c

      subroutine wallflux

      Implicit none

c ... Calc particle & heat fluxes to outer wall surfaces; only PF rad flux?
c ... Use as diagnostic called from BASIS/Python parser (alt. to balancee)


      Use(Dim)            # nx,ny,nxpt
      Use(Comgeo)         # sy
      Use(Phyvar)         # pi, ev
      Use(Postproc)       # pwr_wallh, Pwr_wallz
                          # swallr, swalli, swalle, swbind, swallt
                          # spfwallr
      Use(Compla)         # mi,ni,up,te,ti
      Use(Comflo)         # fniy,feey,feiy
      Use(UEpar)          # ebind

c  Local variables
      integer jx,id,ixi,ixo,ip
      #Former Aux module variables
      integer ix,iy
      real pwallr,pwalli,pwalle,pwbind
      real swallid(0:nx+1,nfsp)

########################################
# First do the wall particle flux
##############################################################
      do ix = 1, nx
        do id = 1, nfsp
	  gwalli(ix,id) = fniy(ix,ny,id)/sy(ix,ny)
          engwalli(ix,id) = 2.*ti(ix,ny+1)/ev + zi(id)*phi(ix,ny+1)
        enddo
      enddo

# Fix corner boundary values
      do id = 1, nfsp
        gwalli(0,id) = gwalli(1,id)
        gwalli(nx+1,id) = gwalli(nx,id)
        engwalli(0,id) = engwalli(1,id)
        engwalli(nx+1,id) = engwalli(nx,id)
      enddo

#######################################
# Now do the wall heat flux
##############################################################
#
# Initialize total powers (local use only)
      pwallr = 0.
      pwalli = 0.
      pwalle = 0.
      pwbind = 0.

# Get radiation power in pwr_plth and pwr_plt
      call pradpltwl

# Fill radiation heat flux arrays
      do ix = 0, nx+1
	swallr(ix) = pwr_wallh(ix) + pwr_wallz(ix)
        do ip = 1, nxpt
          spfwallr(ix,ip) = pwr_pfwallz(ix,ip)+pwr_pfwallh(ix,ip)
        enddo
      enddo

# swalls are ion and electron radial power fluxes in W/m**2
      do ix = 1, nx
        swalli(ix) = 0.
        do id = 1, nfsp
          if (zi(id) .gt. 0) then
            swallid(ix,id) = ( 0.5*mi(id)*upi(ix,ny,id)**2*
     .                        fniy(ix,ny,id) )/sy(ix,ny)
          else    # note: upi=0 for neutral species; use up instead
	    swallid(ix,id) = ( 0.5*mi(id)*up(ix,ny,id)**2*
     .                         fniy(ix,ny,id) )/sy(ix,ny)
          endif
          swalli(ix) = swalli(ix) + swallid(ix,id)
        enddo
        swalli(ix) = swalli(ix) + feiy(ix,ny)/sy(ix,ny)
        swbind(ix) = fniy(ix,ny,1) * ebind*ev / sy(ix,ny)
        swalle(ix) = (feey(ix,ny)+fqy(ix,ny)*phi(ix,ny+1) )/sy(ix,ny)
        swallt(ix) = swalle(ix) + swalli(ix) + swbind(ix) + swallr(ix)
        pwallr = pwallr + swallr(ix)*sy(ix,ny)
        pwalle = pwalle + swalle(ix)*sy(ix,ny)
        pwalli = pwalli + swalli(ix)*sy(ix,ny)
        pwbind = pwbind + swbind(ix)*sy(ix,ny)
      enddo  # end do-loop for ix=1,nx

# Fix corner boundary values
      swallr(0) = swallr(1)
      swalli(0) = swalli(1)
      swalle(0) = swalle(1)
      swbind(0) = swbind(1)
      swallt(0) = swallt(1)
      swallr(nx+1) = swallr(nx)
      swalli(nx+1) = swalli(nx)
      swalle(nx+1) = swalle(nx)
      swbind(nx+1) = swbind(nx)
      swallt(nx+1) = swallt(nx)

      return
      end
******* end of subroutine wallflux *******
c-------------------------------------------------------------------------c

    '''


def checkbal():
    from uedge import bbb
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    plt.ion()
    
    new = UeBalance()
    new.engbal()
    bbb.engbal(1)
    for var in ['fetx', 'fety', 'engerr']:
        dif = abs( (new.__dict__[var] - bbb.__getattribute__(var))/ bbb.__getattribute__(var)**0)
        f, ax = plt.subplots()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        im = ax.imshow(dif, cmap='viridis')
        f.colorbar(im, cax=cax, orientation='vertical')
        ax.set_title(var)
        print(var,  dif.max())
