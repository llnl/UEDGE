



class UeBalance():
    def __init__(self):
        from numpy import zeros
        from uedge import com
        for var in [
            'fetx', 'fety', 'engerr', 'pmloss', 'pmrada', 'pmradm', 'pmpot',
            'peirad', 'pmomv', 'pradrc', 'pradiz', 'prdiss',
            'pibirth', 'pbinde', 'pbindrc', 'pradzbind', 'pradff'
        ]:
            self.__dict__[var] = zeros((com.nx+2, com.ny+2))
        for var in ['icxgas', 'iion', 'irecomb']:
            self.__dict__[var] = zeros((com.nx+2, com.ny+2, com.ngsp))
        self.pradimp = zeros((com.nx+2, com.ny+2, sum(com.nzsp), sum(com.nzsp)))
        self.pwr_plth = zeros((com.ny+2, 2*com.nxpt))
        self.pwr_pltz = zeros((com.ny+2, 2*com.nxpt))
        self.pwr_wallh = zeros((com.nx+2))
        self.pwr_wallz = zeros((com.nx+2))
        self.pwr_pfwallh = zeros((com.nx+2, com.nxpt))
        self.pwr_pfwallz = zeros((com.nx+2, com.nxpt))
        return

    def engbal(self):
        """ Calculates various components of the 2-D energy flow and the 
        ionization and radiation for use in the postprocessing file
        balancee to determine energy balance; these 2-D loops become 
        expensive when done from the parser.
        """
        from uedge import bbb, com, aph
        from numpy import zeros, cos

        for var in ['icxgas', 'iion', 'irecomb', 'pradimp',
            'fetx', 'fety', 'engerr', 'pmloss', 'pmrada', 'pmradm', 'pmpot',
            'peirad', 'pmomv', 'pradrc', 'pradiz', 'prdiss', 'pradzbind',
            'pibirth', 'pbinde', 'pbindrc', 'pradzbind', 'pradff']:
            self.__dict__[var] *= 0


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
                            )**2 * bbb.fniy[ix,iy,ii] \
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

                        self.fetx[ixr,iy] += 0.5*bbb.mi[ii]*up[ixr,iy,ii]**2*bbb.fnix[ixr,iy,ii] \
                        - bbb.ckinfl*0.5*com.sx[ixr,iy]*bbb.visx[ixr,iy,ii]*com.gx[ixr,iy] * \
                        (up[ixr,iy,ii]**2 - up[ixr1,iy,ii]**2)
                self.fetx[ixt,iy] += bbb.feex[ixt,iy] + bbb.feix[ixt,iy]
                self.fetx[ixr,iy] += bbb.feex[ixr,iy] + bbb.feix[ixr,iy]

        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx]+1, com.ixrb[jx]+1):
                for iy in range(1,com.ny+1):
                    ix1 = bbb.ixm1[ix,iy]
                    ix2 = bbb.ixp1[ix,iy]
                    self.pmloss[ix,iy] = (1-bbb.ismolcrm)*bbb.cnsor*( \
                            bbb.ediss*bbb.ev*(0.5*bbb.psordis[ix,iy,1]) \
                            + bbb.ceisor*bbb.eion*bbb.ev*(bbb.psordis[ix,iy,1])
                        ) + bbb.ismolcrm*bbb.cnsor*( \
                            bbb.cmesori*(bbb.emolia[ix,iy,0] + bbb.emolia[ix,iy,1]) \
                            + bbb.cmesore*bbb.edisse[ix,iy]
                        )
                    if bbb.ishymol:
                        self.pmpot[ix,iy] = bbb.ismolcrm*bbb.ng[ix,iy,1]*com.vol[ix,iy] \
                            *aph.sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], 22)
                        self.pmrada[ix,iy] = bbb.ismolcrm*bbb.ng[ix,iy,1]*com.vol[ix,iy] \
                            *aph.sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], 23)
                        self.pmradm[ix,iy] = bbb.ismolcrm*bbb.ng[ix,iy,1]*com.vol[ix,iy] \
                            *aph.sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], 24)
# Here peirad includes sum of electron and ion energy losses; note that binding
# energy is included in eeli term, but it is carried by the ions.
# Note also that eion and ediss generally balance in the next line
# because should have ediss=2*eion - transfer from electron to ion energy
                    self.peirad[ix,iy] = bbb.cnsor*(
                            bbb.erliz[ix,iy] + bbb.erlrc[ix,iy]
                            + bbb.ebind*bbb.ev*bbb.psor[ix,iy,0] \
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

        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx]+1, com.ixrb[jx]+1):
                for iy in range(1, com.ny+1):
                    for ig in range(com.ngsp):
                        if ((bbb.ishymol == 0) or (ig != 1)):
                            self.iion[ix,iy,ig] -= bbb.cnsor*bbb.qe*bbb.psorg[ix,iy,ig]
                            self.irecomb[ix,iy,ig] -= bbb.cnsor*bbb.qe*bbb.psorrg[ix,iy,ig]
                            self.icxgas[ix,iy,ig] -= bbb.qe*bbb.psorcxg[ix,iy,ig]
        

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
                                

        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx], com.ixrb[jx]+1):
                for iy in range(1,com.ny+1):
                    if (bbb.isimpon in [2,7]):
                        self.pradff[ix,iy] += bbb.pradcff[ix,iy]*com.vol[ix,iy]
                    if bbb.isimpon > 2:
                        for ig in range(com.nhgsp+1, com.ngsp+1):
                            for jz in range(com.ngsp - com.nhgsp):
                                for iimp in range(com.nzsp[jz]):
                                    ii = iimp + com.nhgsp + sum(com.nzsp[:jz]) - 1
                                    self.pradimp[ix,iy, ii, jz] += bbb.pradz[ix,iy,ii,jz]*com.vol[ix,iy]
                self.pradzbind[ix,iy] = (bbb.pwrze - bbb.prad)[ix,iy]*com.vol[ix,iy]
                            
                        

    def calc_engerr(self, pwrin=1, redo=True):
        from uedge import bbb, com

        self.engerr *= 0
        if redo:
            self.engbal()

        for jx in range(com.nxpt):
            for ix in range(com.ixlb[jx]+1, com.ixrb[jx]+1):
                for iy in range(1,com.ny+1):
                    ix1 = bbb.ixm1[ix,iy]
                    self.engerr[ix,iy] = (
                            self.fetx[ix1,iy] - self.fetx[ix,iy] + self.fety[ix,iy-1] \
                            - self.fety[ix,iy] - self.peirad[ix,iy] - bbb.png2ni[ix,iy]
                        ) / abs(pwrin)
                    if bbb.isimpon != 0:
                        self.engerr[ix,iy] -= bbb.prad[ix,iy]*com.vol[ix,iy]/abs(pwrin)                        

    def pradpltwl(self):
        """ Calc radiation power to divertor and outer wall surfaces """
        from uedge import bbb, com
        from numpy import zeros, arctan2, pi, cos, sum
        from copy import deepcopy
        nj = com.nxomit 
        prdu = deepcopy(bbb.prad)

        # Initialize arrays for subsequent runs
        for var in ['pwr_pltz', 'pwr_plth', 'pwr_wallh', 'pwr_wallz', 
            'pwr_pfwallh', 'pwr_pfwallz']:
            self.__dict__[var] *= 0

        for ip in range(2*com.nxpt):
            if (ip % 2) == 0: # Even numbers
                ixv = com.ixlb[int( ip > 1) + int(ip > 3)]
            else:
                ixv = com.ixrb[int( ip > 1) + int(ip > 3)]+1
            for iyv in range(1, com.ny+1):
                for iy in range(1, com.ny+1):
                    for ix in range(1,com.nx+1):
                        thetaray1 = arctan2(
                            com.zm[ixv+nj, iyv, 1] - com.zm[ix+nj, iy, 0],
                            com.rm[ixv+nj, iyv, 1] - com.rm[ix+nj, iy, 0]
                        )
                        thetaray2 = arctan2(
                            com.zm[ixv+nj, iyv, 3] - com.zm[ix+nj, iy, 0],
                            com.rm[ixv+nj, iyv, 3] - com.rm[ix+nj, iy, 0]
                        )
                        dthgy = abs(thetaray1 - thetaray2)
                        frth = min(dthgy, 2*pi - dthgy)/2/pi
                        sxo = com.sx[ixv, iyv]/cos(com.angfx[ixv,iyv])
                        self.pwr_pltz[iyv,ip] += prdu[ix,iy]*com.vol[ix,iy]*frth/sxo
                        self.pwr_plth[iyv,ip] += ( (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev) \
                            *bbb.psor[ix,iy,0] + bbb.erlrc[ix,iy])*frth/sxo
            # Set corner values
            self.pwr_pltz[0,ip] = self.pwr_pltz[1,ip]
            self.pwr_pltz[com.ny+1,ip] = self.pwr_pltz[com.ny,ip]
            self.pwr_plth[0,ip] = self.pwr_plth[1,ip]
            self.pwr_plth[com.ny+1,ip] = self.pwr_plth[com.ny,ip]


        """ OUTER WALL """
        iyv = com.ny + 1
        for ixv in range(1, com.nx+1):
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
                    frth = min(dthgz, 2*pi - dthgz)/2/pi
                    self.pwr_wallz[ixv] += prdu[ix,iy]*com.vol[ix,iy]*frth/com.sy[ixv,iyv]
                    self.pwr_wallh[ixv] += ( (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev)*bbb.psor[ix,iy,0] \
                        +bbb.erlrc[ix,iy])*frth/com.sy[ixv,iyv]
        self.pwr_wallz[0] = self.pwr_wallz[1]	# Because prad(0,) is not calculated
        self.pwr_wallz[com.nx+1] = self.pwr_wallz[com.nx]
        self.pwr_wallh[0] = self.pwr_wallh[1]
        self.pwr_wallh[com.nx+1] = self.pwr_wallh[com.nx]
    

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
                        self.pwr_pfwallz[ixv,ip] += prdu[ix,iy]*com.vol[ix,iy]*frth/com.sy[ixv,iyv]
                        self.pwr_pfwallh[ixv,ip] += (
                            (bbb.eeli[ix,iy] - bbb.ebind*bbb.ev)*bbb.psor[ix,iy,0] \
                            + bbb.erlrc[ix,iy])*frth/com.sy[ixv,iyv]
                if ((ixv>com.ixpt1[ip]) and (ixv<com.ixpt2[ip]+1)):
                    self.pwr_pfwallh[ixv,ip] = 0
                    self.pwr_pfwallz[ixv,ip] = 0
            self.pwr_pfwallz[0,ip] = self.pwr_pfwallz[1,ip]
            self.pwr_pfwallz[com.nx+1,ip] = self.pwr_pfwallz[com.nx,ip]
            self.pwr_pfwallh[0,ip] = self.pwr_pfwallh[1,ip]
            self.pwr_pfwallh[com.nx+1,ip] = self.pwr_pfwallh[com.nx,ip]

    def platedist(self):
        from numpy import sum, cos, zeros, cumsum, concatenate
        from uedge import com

        self.ydpin, self.ydpout = zeros((com.ny+2,)), zeros((com.ny+2,))
        
        self.ydpin = concatenate(([0],cumsum( ( 1/com.gy[0,:com.ny+1] + 1/com.gy[0,1:com.ny+2]) / \
                    ( cos(com.angfx[0,:com.ny+1])+cos(com.angfx[0,1:com.ny+2]) )))) \
                - sum( 1/(com.gy[0,:com.iysptrx+1]*cos(com.angfx[0,:com.iysptrx+1])))
        self.ydpout = concatenate(([0],cumsum( ( 1/com.gy[com.nx,:com.ny+1] + 1/com.gy[com.nx,1:com.ny+2]) / \
                    ( cos(com.angfx[com.nx,:com.ny+1])+cos(com.angfx[com.nx,1:com.ny+2]) )))) \
                    - sum( 1/(com.gy[com.nx,:com.iysptrx+1]*cos(com.angfx[com.nx,:com.iysptrx+1])))


        return

    def balance_power(self, redo=True):
        # converted into UeBalance by AH on March 25, 2025
        # balancee as of 22Dec21, now includes molecular energy fluxes to plate/walls
        # sdmin, sdmout, pwallm, and ppfm, as well as impurity ion binding energy
        # heating on plates (sbindzin, sbindzout) (GDP 11 June 2018).
        # As before, also ncludes neutral energy fluxes to plates and walls.
        #
        # added sbindzin, sbindzout, pbindzin, pbindzout	11 June 2018 (GDP)
        # TOTAL POWER AND PARTICLE FLOWS
        from uedge import bbb, com
        from numpy import cos, zeros, sqrt, pi, sum

        def have(x1, x2):
           # The function ave gives the harmonic average of two numbers
           return (2*x1*x2)/(x1+x2)

        ixdivin=0
        iycore=0
        ixineut=1
        iybegin=1

        """ INITIALIZE VARIABLES """
        for var in [
            'ppfi', 'pwalle', 'ppfe', 'pwallm', 'ppfm', 'pwallbd', 'ppfbd',
            'pradhwall', 'pradzwall', 'pradhpf', 'pradzpf', 'pneutpf', 'pneutw',
            'pdiviin', 'pdiviout', 'pdivein', 'pdiveout', 'pdivmin', 'pdivmout',
            'pbindout', 'pbindin', 'pdivnout', 'pdivnin', 'pbindzout', 
            'pbindzin', 'pneutout', 'pneutin', 'pradhout', 'pradhin',
            'pradzin', 'pradzout', 'pwalli', 'psepi', 'psepe', 'psepbd',
            'pbcorei', 'pbcoree', 'pbcorebd'
        ]:
            self.__dict__[var] = 0

        for var in ['sneutpf', 'sneutw']:
            self.__dict__[var] = zeros((com.nx+2,))

        for var in [
            'sdrout', 'sdrin', 'sdiout', 'sdeout', 'sdiin', 'sdmout',
            'sdmin', 'sdein', 'sdtout', 'sdtin', 'sbindout', 'sbindin',
            'sbindzout', 'sbindzin', 'sneutout', 'sneutin',
            'sdioutd', 'sdiind', 'sdnout', 'sdnin'
        ]:
            self.__dict__[var] = zeros((com.ny+2,))

        for var in ['sdioutd', 'sdiind']:
            self.__dict__[var] = zeros((com.ny+2,com.nfsp))

        #two ifs to be able to used old executables
        if(com.islimon == 1) and (com.nyomitmx != 0):
           com.nx = com.ix_lim-1
           ixdivin = com.ix_lim+1
           iybegin = com.iy_lims
          
        # note: ixi=ixdivin has been set to 1 to allow velocity derivatives to be calc.
        ixi=ixdivin
        ixo=com.nx

        iycore = int(bbb.isguardc == 0)
        # Determine if molecular hydrogen energy fluxes are present
        if(com.ngsp >= 2) and (bbb.ishymol*bbb.istgon[1] == 1):
           ishymoleng = 1
        else: 
           ishymoleng = 0

        # power outflow from separatrix
        try:
           fluxfacy = bbb.fluxfacy
        except:
           fluxfacy=1.
        try:
           cenggpl = bbb.cenggpl
        except:
           cenggpl = 0
        try:
           cenggw = bbb.cenggw
        except:
           cenggw = 0

        for ix in range(max(com.ixpt1[0]+1,0), com.ixpt2[0]+1):
           ix1 = bbb.ixm1[ix,com.iysptrx]
           self.psepi += fluxfacy*( bbb.feiy[ix,com.iysptrx] + (bbb.mi[0]/32)*
                       (bbb.up[ix1,com.iysptrx  ,0]+bbb.up[ix,com.iysptrx  ,0]+
                        bbb.up[ix1,com.iysptrx+1,0]+bbb.up[ix,com.iysptrx+1,0])**2*bbb.fniy[ix,com.iysptrx,0] -
               bbb.cfvisy*0.125*com.sy[ix,com.iysptrx]*have( bbb.visy[ix,com.iysptrx,0]*com.gy[ix,com.iysptrx],
                                   bbb.visy[ix,com.iysptrx+1,0]*com.gy[ix,com.iysptrx+1]) *
                    ( (bbb.up[ix1,com.iysptrx+1,0]+bbb.up[ix,com.iysptrx+1,0])**2 -
                      (bbb.up[ix1,com.iysptrx  ,0]+bbb.up[ix,com.iysptrx  ,0])**2 ) )
           ix1 = bbb.ixm1[ix,iycore]
           self.psepe += fluxfacy*bbb.feey[ix,com.iysptrx]
           self.pbcorei += fluxfacy*( bbb.feiy[ix,iycore] + (bbb.mi[0]/32)*
                       (bbb.up[ix1,iycore  ,0]+bbb.up[ix,iycore  ,0]+
                        bbb.up[ix1,iycore+1,0]+bbb.up[ix,iycore+1,0])**2*bbb.fniy[ix,iycore,0] -
                 bbb.cfvisy*0.125*com.sy[ix,iycore]*have( bbb.visy[ix,iycore  ,0]*com.gy[ix,iycore],
                                             bbb.visy[ix,iycore+1,0]*com.gy[ix,iycore+1] ) *
                    ( (bbb.up[ix1,iycore+1,0]+bbb.up[ix,iycore+1,0])**2 -
                      (bbb.up[ix1,iycore  ,0]+bbb.up[ix,iycore  ,0])**2 ) )
           self.pbcoree += fluxfacy*bbb.feey[ix,iycore]
           self.pbcorebd += fluxfacy*bbb.fniy[ix,iycore,0]*bbb.ebind*bbb.ev
           self.psepbd += fluxfacy*bbb.fniy[ix,com.iysptrx,0]*bbb.ebind*bbb.ev

        #

        #       POWER INPUT FROM BIAS SUPPLY
        self.p_bias = -sum(bbb.fqx[com.nx,1:com.ny+1]*(bbb.phi0r[1:com.ny+1]-bbb.phi0l[1:com.ny+1]))

        #       TOTAL BIAS CURRENT
        self.i_bias = sum(bbb.fqx[com.nx,1:com.ny+1])

        #       ION SATURATION CURRENT
        self.i_sat_outer = bbb.qe*bbb.zi[0]*sum(bbb.fnix[com.nx,1:com.ny+1,0])
        self.i_sat_inner = bbb.qe*bbb.zi[0]*sum(bbb.fnix[0,1:com.ny+1,0])

        #       CURRENT AND POWER FROM FIXED VOLUME SOURCES
        self.p_e_vol = sum(bbb.pwrsore)
        self.p_i_vol = sum(bbb.pwrsori)

        if bbb.l_parloss <= 1e9:
          self.p_e_vol -= sum(bbb.nuvl[:,:,0]*com.vol*bbb.ne*bbb.bcee*bbb.te)
          self.p_i_vol -= sum(bbb.nuvl[:,:,0]*com.vol*bbb.ni[:,:,0]*bbb.bcei*bbb.ti)


        #########################################################################
        # 2-D arrays for checking energy conservation, and for calculation of
        # ionization and radiation sources 

        self.ptotin = self.pbcoree + self.pbcorei + self.p_i_vol + self.p_e_vol

        if redo:
            self.engbal(self.ptotin)

        #########################################################################

        # power incident on divertor plate
        # Allow use of "old" or "new" switches for neutral energy loss


        # First get radiation power in pwr_plth and pwr_pltz
        if redo:
            self.pradpltwl()
        self.sdrin = self.pwr_plth[:,0]+self.pwr_pltz[:,0]
        self.sdrout = self.pwr_plth[:,1]+self.pwr_pltz[:,1]


        # here the sds are ion and electron poloidal power fluxes in W/m**2
        for iy in range(iybegin, com.ny+1):
           sxo = com.sx[ixo,iy]/(cos(com.angfx[ixo,iy]))
           sxi = com.sx[ixi,iy]/(cos(com.angfx[ixi,iy]))
           vxno =  0.25*sqrt(8*bbb.tg[com.nx,iy,0]/(pi*bbb.mg[0]))
           vxni =  0.25*sqrt(8*bbb.tg[ixineut,iy,0]/(pi*bbb.mg[0]))

           for id in range(com.nfsp): # note: upi=0 for the netural species
              if bbb.zi[id] > 0:
                 self.sdioutd[iy,id] = ( 0.5*bbb.mi[id]*bbb.upi[ixo,iy,id]**2*bbb.fnix[ixo,iy,id] )/sxo
                 self.sdiind[iy,id] = ( -0.5*bbb.mi[id]*bbb.upi[ixi,iy,id]**2*bbb.fnix[ixi,iy,id] )/sxi
              else:  #zero-out neutrals as fnix will be into plasma 5/27/08 - TDR
                 if bbb.ishymol == 0: # recompute parallel fnix
                    self.sdioutd[iy,id] = 0.5*bbb.mi[id]*abs(bbb.up[ixo,iy,id])**3*bbb.ni[ixo,iy,id]*com.rrv[ixo,iy]
                    self.sdiind[iy,id] = 0.5*bbb.mi[id]*abs(bbb.up[ixi,iy,id])**3*bbb.ni[ixi,iy,id]*com.rrv[ixi,iy]
                 else:   # for molecules, fnix should be ok
                    self.sdiout[iy,id]=0#1e-20*( -0.5*bbb.mi[id]*bbb.up[ixo,iy,id]**2*bbb.fnix[ixo,iy,id] )/sxo
                    self.sdiind[iy,id]=0#1e-20*( -0.5*bbb.mi[id]*bbb.up[ixi,iy,id]**2*bbb.fnix[ixi,iy,id] )/sxi
              self.sdiout[iy] += self.sdioutd[iy,id]
              self.sdiin[iy]  += self.sdiind[iy,id]

           # AH: bbb.ckinfl treated as array, although it is double. Legacy switch?
           for id in range(com.nusp):      # note: up for the netural species in nonzero
              tempvo =  - bbb.ckinfl*0.5*com.sx[ixo,iy]*bbb.visx[ixo,iy,id]*com.gx[ixo,iy]*\
                            ( bbb.up[ixo,iy,id]**2 - bbb.up[ixo-1,iy,id]**2 ) /sxo
              tempvi =  + bbb.ckinfl*0.5*com.sx[ixi,iy]*bbb.visx[ixi+1,iy,id]*com.gx[ixi+1,iy]*\
                           ( bbb.up[ixi+1,iy,id]**2 - bbb.up[ixi,iy,id]**2 ) / sxi
              self.sdioutd[iy,id] += tempvo
              self.sdiout[iy] += tempvo
              self.sdiind[iy,id] += tempvi
              self.sdiin[iy] += tempvi

           self.sdiout[iy] += bbb.feix[ixo,iy]/sxo
           self.sbindout[iy] = bbb.fnix[ixo,iy,0] * bbb.ebind*bbb.ev / sxo
           if ishymoleng==1:  #mol heat flux; drift eng small,<0
             self.sdmout[iy] += bbb.fegx[ixo,iy,1]/sxo

        #  Compute binding-energy energy fluxes for impurities
           for id in range(com.nzdf):
              if (id == 0):
                 id2min = com.nhsp
                 id2max = id2min +com.nzsp[id]-1
              else:
                 id2min = com.nhsp+sum(com.nzsp[1:id-1])
                 id2max = id2min + com.nzsp[id]-1
              for id2 in range(id2min, id2max+1):
                 for id3 in range(bbb.znucl[id]):
                    self.sbindzout[iy] += bbb.fnix[ixo,iy,id2]*bbb.ebindz(id3,bbb.znucl[id2])*bbb.ev/sxo
                    self.sbindzin[iy] -= bbb.fnix[ixi,iy,id2]*bbb.ebindz(id3,bbb.znucl[id2])*bbb.ev/sxo

           self.sneutout[iy] = cenggpl*2.*vxno*bbb.ng[com.nx,iy,0]*bbb.tg[com.nx,iy,0]
           self.sdeout[iy] = ( bbb.feex[ixo,iy]+bbb.fqx[ixo,iy]*(bbb.phi[ixo,iy]-bbb.phi0r[iy]) )/sxo
           self.sdtout[iy] = self.sdeout[iy] + self.sdiout[iy] + self.sbindout[iy] + self.sdmout[iy] + \
                           self.sbindzout[iy] + self.pwr_plth[iy,1] + self.pwr_pltz[iy,1]
           self.sdiin[iy]  -= bbb.feix[ixi,iy]/sxi 
           if ishymoleng==1: #mol heat flux; drift eng small,<0
             self.sdmin[iy] -= bbb.fegx[ixi,iy,1]/sxo
           
           self.sbindin[iy] = - bbb.fnix[ixi,iy,1] * bbb.ebind*bbb.ev / sxi
           self.sneutin[iy] = cenggpl*2.*vxni*bbb.ng[ixineut,iy,0]*bbb.tg[ixineut,iy,0]
           self.sdein[iy]  = -( bbb.feex[ixi,iy] + .001*bbb.fqx[ixi,iy]*(bbb.phi[ixi,iy]-bbb.phi0l[iy]) )/sxi
           self.sdtin[iy] = self.sdein[iy] + self.sdiin[iy] + self.sdmin[iy] + self.sbindin[iy] + \
                          self.sbindzin[iy] + self.pwr_plth[iy,0] + self.pwr_pltz[iy,0]
           self.pdiviout += self.sdiout[iy]*sxo ## + self.sdioutd[iy,1]*sxo
           self.pdiveout += self.sdeout[iy]*sxo
           self.pdivmout += self.sdmout[iy]*sxo
           self.pdiviin  += self.sdiin[iy]*sxi ## + self.sdiind[iy,1]*sxi
           self.pdivein  += self.sdein[iy]*sxi
           self.pdivmin  += self.sdmin[iy]*sxi
           self.pbindout += self.sbindout[iy]*sxo
           self.pbindzout += self.sbindzout[iy]*sxo
           self.pbindin += self.sbindin[iy]*sxi
           self.pbindzin += self.sbindzin[iy]*sxi
           self.pneutout += 0*self.sneutout[iy]*sxo  # included in self.pdiviout
           self.pneutin += 0.*self.sneutin[iy]*sxo     # included in self.pdiviin
           self.pradhout += self.pwr_plth[iy,1]*sxo
           self.pradzout += self.pwr_pltz[iy,1]*sxo
           self.pradhin += self.pwr_plth[iy,0]*sxi
           self.pradzin += self.pwr_pltz[iy,0]*sxi
           if bbb.isupgon[0] == 1: # Approx. neutral energy flux
              self.sdnout[iy]=self.sdioutd[iy,1] + ( com.sx[ixo,iy]*bbb.hcxn[ixo,iy]*com.gxf[ixo,iy]* \
                                                    (bbb.ti[ixo,iy]-bbb.ti[ixo+1,iy]) + \
                                          2.5*bbb.fnix[ixo,iy,1]*bbb.ti[ixo+1,iy] ) / sxo
              self.sdnin[iy]=self.sdiind[iy,1] + ( com.sx[ixi,iy]*bbb.hcxn[ixi,iy]*com.gxf[ixi,iy]* \
                                                  (bbb.ti[ixi,iy]-bbb.ti[ixi+1,iy]) + \
                                          2.5*bbb.fnix[ixi,iy,1]*bbb.ti[ixi,iy] ) / sxi
              self.pdivnout += self.sdnout[iy]*sxo
              self.pdivnin += self.sdnin[iy]*sxi

        # 
        # added, 11Jun2018 GDP - only diagnostic; not used below
        ptot=self.pdiviout+self.pdiveout+self.pdiviin+self.pdivein+self.pbindout+self.pbindin+\
                   self.pbindzin+self.pbindzout+self.pdivmout+self.pdivmin+self.pneutout+self.pneutin
        #

        # Fix up boundary values
        for var in ['sdtin', 'sdein', 'sdiin', 'sdtout', 'sdeout', 'sdiout']:
            self.__dict__[var][0] = self.__dict__[var][1]
            self.__dict__[var][-1] = self.__dict__[var][-2]


        #
        self.ptotpartin = self.pdiviin + self.pdivein + self.pbindin + self.pbindzin + self.pdivmin  ##self.pneutin
        self.ptotpartout = self.pdiviout + self.pdiveout + self.pbindout + self.pbindzout + self.pdivmout  ##self.pneutout
        self.ptotpart = self.ptotpartin + self.ptotpartout  ##self.pneutout+self.pneutin
        self.ptotrad = self.pradhout + self.pradzout + self.pradhin + self.pradzin
        #
       

        #
        iywall=com.ny       # DEFINITION
        iypf=0
        #
        # power flow to vessel and private flux wall


        for ix in range(com.nx+2):
           ix1 = bbb.ixm1[ix,iywall]
           vynw = 0.25*sqrt(8*bbb.tg[ix,com.ny+1,0]/(pi*bbb.mg[0]))
           self.pwalli += fluxfacy*( bbb.feiy[ix,iywall] +
                    0.125*bbb.mi[0]*(bbb.up[ix1,iywall,0]+bbb.up[ix,iywall,0])**2*bbb.fniy[ix,iywall,0] -
                    bbb.cfvisy*0.125*com.sy[ix,iywall]*have( bbb.visy[ix,iywall,0]*com.gy[ix,iywall],
                                             bbb.visy[ix,iywall+1,0]*com.gy[ix,iywall+1] ) *
                    ( (bbb.up[ix1,iywall+1,0]+bbb.up[ix,iywall+1,0])**2 -
                      (bbb.up[ix1,iywall  ,0]+bbb.up[ix,iywall  ,0])**2 ) )
           self.pwalle += fluxfacy*bbb.feey[ix,iywall]
           self.pwallbd += fluxfacy*bbb.fniy[ix,iywall,0]*bbb.ebind*bbb.ev
           self.pradhwall += fluxfacy*self.pwr_wallh[ix]*com.sy[ix,iywall]
           self.pradzwall += fluxfacy*self.pwr_wallz[ix]*com.sy[ix,iywall]
           self.sneutw[ix] = cenggw*2.*vynw*bbb.tg[ix,com.ny+1,0]*com.sx[ix,com.ny]
           self.pneutw += self.sneutw[ix]
           if ishymoleng == 1:  #molec temp eqn active
              self.pwallm += bbb.fegy[ix,iywall,1]


        for ix in range(0, com.ixpt1[0]+1):
           ix1 = bbb.ixm1[ix,iypf]
           vynpf = 0.25*sqrt(8*bbb.tg[ix,0,1]/(pi*bbb.mg[0]))
           self.ppfi -= fluxfacy*(bbb.feiy[ix,iypf] -
                    0.125*bbb.mi[0]*(bbb.up[ix1,iypf,0]+bbb.up[ix,iypf,0])**2*bbb.fniy[ix,iypf,0] +
                    bbb.cfvisy*0.125*com.sy[ix,iypf]*have( bbb.visy[ix,iypf,0]*com.gy[ix,iypf ],
                                           bbb.visy[ix,iypf+1,0]*com.gy[ix,iypf+1] ) *
                    ( (bbb.up[ix1,iypf+1,0]+bbb.up[ix,iypf+1,0])**2 -
                      (bbb.up[ix1,iypf  ,0]+bbb.up[ix,iypf  ,0])**2 ) )
           self.ppfe -= fluxfacy*bbb.feey[ix,iypf]
           self.ppfbd -= fluxfacy*bbb.fniy[ix,iypf,0]*bbb.ebind*bbb.ev
        ##   self.pradhpf += fluxfacy*pwr_pfh[ix]*com.sy[ix,iywall]
        ##   self.pradzpf += fluxfacy*pwr_pfz[ix]*com.sy[ix,iywall]
           self.sneutpf[ix] = cenggw*2.*vynw*bbb.tg[ix,0,1]*com.sx[ix,0]
           self.pneutpf += self.sneutpf[ix]
           if ishymoleng==1:  #molec temp eqn active
             self.ppfm -= bbb.fegy[ix,iypf,1]
           

        for ix in range(com.ixpt2[0]+1, com.nx+2):
           ix1 = bbb.ixm1[ix,iypf]
           vynpf = 0.25*sqrt(8*bbb.tg[0,iy,0]/(pi*bbb.mg[0]))
           self.ppfi -= fluxfacy*( bbb.feiy[ix,iypf] -
                    0.125*bbb.mi[0]*(bbb.up[ix1,iypf,0]+bbb.up[ix,iypf, 0])**2*bbb.fniy[ix,iypf,0] +
                    bbb.cfvisy*0.125*com.sy[ix,iypf]*have( bbb.visy[ix,iypf,0]*com.gy[ix,iypf],
                                           bbb.visy[ix,iypf+1,0]*com.gy[ix,iypf+1] ) *
                    ( (bbb.up[ix1,iypf+1,0]+bbb.up[ix,iypf+1,0])**2 -
                      (bbb.up[ix1,iypf  ,0]+bbb.up[ix,iypf  ,0])**2 ) )
           self.ppfe -= fluxfacy*bbb.feey[ix,iypf]
           self.ppfbd -= fluxfacy*bbb.fniy[ix,iypf,0] * bbb.ebind*bbb.ev
           if ishymoleng==1: #molec temp eqn active
             self.ppfm -= bbb.fegy[ix,iypf,1]
           
        #
        # ion current to vessel wallst.

        # Calculate power into plasma volume for normalizing factor in power balance
        self.pnormpb = 0.
        if self.ptotpartin > 0.:
            self.pnormpb += self.ptotpartin
        if self.ptotpartout < 0.:
            self.pnormpb += self.ptotpartout
        if self.pbcoree + self.pbcorei + self.pbcorebd + self.p_i_vol > 0.:
            self.pnormpb += self.pbcoree + self.pbcorei + self.pbcorebd + self.p_i_vol
        if self.p_i_vol+self.p_e_vol > 0.:
            self.pnormpb += self.p_i_vol + self.p_e_vol

        self.newpnormpb = 0
        for jx in range(com.nxpt):
            0

        # Here particle power and radiation power are separate terms
        self.powbal = (self.pbcoree + self.pbcorei + self.pbcorebd + self.p_i_vol + self.p_e_vol + sum(bbb.wjdote[1:-1,1:-1]) - self.ptotpart -
                  self.pwalli-self.pwalle-self.pwallm-self.pwallbd-self.ppfi-self.ppfe-self.ppfm-self.ppfbd-sum(self.pradimp) +
                  sum( - self.pradff - self.pradiz - self.pradrc - self.prdiss - self.pmomv + self.pibirth))/ self.pnormpb
        if abs(self.pdivein+self.pdiviin) > 10.*abs(self.ptotin + sum(bbb.wjdote[1:-1,1:-1])):   #for cases with no radial pwr input
           self.powbal *= abs(self.ptotin + sum(bbb.wjdote[1:-1,1:-1]))/abs(self.pdivein+self.pdiviin)




    def balance_particles(self, redo=True):
        # converted into UeBalance by AH on March 25, 2025
        # balancee as of 22Dec21, now includes molecular energy fluxes to plate/walls
        # sdmin, sdmout, pwallm, and ppfm, as well as impurity ion binding energy
        # heating on plates (sbindzin, sbindzout) (GDP 11 June 2018).
        # As before, also ncludes neutral energy fluxes to plates and walls.
        #
        # added sbindzin, sbindzout, pbindzin, pbindzout	11 June 2018 (GDP)
        # TOTAL POWER AND PARTICLE FLOWS

        from uedge import bbb, com
        from numpy import sum, zeros
        
        if redo:
            self.engbal()

        ixdivin=0
        iycore=0
        ixineut=1
        iybegin=1
        #two ifs to be able to used old executables
        if(com.islimon == 1) and (com.nyomitmx != 0):
           com.nx = com.ix_lim-1
           ixdivin = com.ix_lim+1
           iybegin = com.iy_lims

        for var in ['iwall', 'ipf', 'igaswall', 'igaspf', 'igascr',
            'idivout', 'idivin', 'icore'
        ]:
            self.__dict__[var] = zeros((com.nfsp,))
           
        for var in ['igasout', 'igasin']:
            self.__dict__[var] = zeros((com.ngsp,))
  
        for var in ['igastot', 'igasdenom']:
            self.__dict__[var] = 0

        try:
           fluxfacy = bbb.fluxfacy
        except:
           fluxfacy=1.

        iywall=com.ny       # DEFINITION
        iypf=0
        for ix in range(0, com.nx+2):
           for id in range(com.nfsp):
              self.iwall[id] += fluxfacy*bbb.fniy[ix,iywall,id]
           for igsp in range(com.ngsp):
              self.igaswall[igsp] += fluxfacy*bbb.fngy[ix,iywall,igsp]

        for ix in range(0, com.ixpt1[0]+1):
           for id in range(0, com.nfsp):
              self.ipf[id] += fluxfacy*bbb.fniy[ix,iypf,id]
           for igsp in range(com.ngsp):
              self.igaspf[igsp] += fluxfacy*bbb.fngy[ix,iypf,igsp]

        for ix in range(com.ixpt2[0]+1, com.nx+2):
           for id in range(0, com.nfsp):
              self.ipf[id] += fluxfacy*bbb.fniy[ix,iypf,id]
           for igsp in range(com.ngsp):
              self.igaspf[igsp] += fluxfacy*bbb.fngy[ix,iypf,igsp]

        for ix in range(max(com.ixpt1[0]+1,0), com.ixpt2[0]+1):
           for igsp in range(com.ngsp):
              self.igascr[igsp] =  fluxfacy*bbb.fngy[ix,iypf,igsp]

        # ion current to divertor plate

        for iy in range(iybegin,com.ny+1):
           for id in range(com.nfsp):
              self.idivout[id] += bbb.fnix[com.nx,iy,id]
              self.idivin[id] += bbb.fnix[ixdivin,iy,id]
           for igsp in range(com.ngsp):
              if bbb.isupgon[igsp] == 0:
                 self.igasout[igsp] += bbb.fngx[com.nx,iy,igsp]
                 self.igasin[igsp] += bbb.fngx[ixdivin,iy,igsp]
              else:
                 self.igasout[igsp] += bbb.fnix[com.nx,iy,1]
                 self.igasin[igsp] += bbb.fnix[ixdivin,iy,1]

        # ion current to the core
        for ix in range(max(0, com.ixpt1[0]+1), com.ixpt2[0]+1):
           for id in range(com.nfsp):
              self.icore[id] += fluxfacy*bbb.fniy[ix,0,id]
 
           
        for var in [
            'idivout', 'idivin', 'igasout', 'igasin', 'icore', 'iwall',
            'igaswall', 'ipf', 'igaspf', 'igascr'
        ]:
            self.__dict__[var] *= bbb.qe

        self.i_vol = bbb.qe*sum(bbb.volpsor[:,:,0])

        if bbb.l_parloss <= 1e9:
          self.i_vol -= sum(bbb.nuvl[:,:,0]*com.vol*bbb.ni[:,:,0])*bbb.qe


        for id in range(com.ngsp):
           self.igastot += self.igasin[id]-self.igasout[id]-self.igaswall[id]+self.igaspf[id]\
                             + self.igascr[id] + self.i_vol
           self.igasdenom += self.igasin[id] -self.igasout[id]

        if bbb.ishymol == 1:  # add id=2 case again because of 2 atoms/molecule
           self.igastot += self.igasin[1]-self.igasout[1]-self.igaswall[1]+self.igaspf[1]\
                             + self.igascr[1]
           self.igasdenom += (self.igasin[1] - self.igasout[1])

        self.deligas = (self.igastot - sum(self.iion) - sum(self.irecomb) -sum(self.icxgas)) / self.igasdenom



    def print_balance_power(self, redo=True):
        from numpy import sum
        from uedge import bbb, com
        if redo:
            self.balance_power()
        print("\nPower Flow [Watts] from Core to Scrape-off Layer is:")
        print("   Total: {:.2e}".format(self.pbcorei + self.pbcoree + self.pbcorebd))
        print("   Ion contribution: {:.2e}".format(self.pbcorei))
        print("   Electron contribution: {:.2e}".format(self.pbcorei))
        print("   Binding energy contribution: {:.2e}".format(self.pbcorei))
        print("\nPower Flow [Watts] over the separatrix to Scrape-off Layer is:")
        print("Power Flow [Watts] over the separatrix to Scrape-off Layer is:")
        print("   Total: {:.2e}".format(self.psepi + self.psepe + self.psepbd))
        print("   Ion contribution: {:.2e}".format(self.psepi))
        print("   Electron contribution: {:.2e}".format(self.psepe))
        print("   Binding energy contribution: {:.2e}".format(self.psepbd))
        #
        print("\nPower Input [Watts] from Fixed Volume Sources (Sinks):")
        print(f"   {self.p_i_vol}")
        print(f"   {self.p_e_vol}")
        #
        print("\nPower Flows [Watts] incident on Divertor Plates are:")
        print(" Inner plate:")
        print("   Total: {:.2e}".format(self.pdiviin+self.pdivein+self.pdivmin,self.pbindin+self.pneutin+self.pradhin+self.pradzin))
        print("   Ion contribution: {:.2e}".format(self.pdiviin))
        print("   Electron contribution: {:.2e}".format(self.pdivein))
        print("   Molecular contribution: {:.2e}".format(self.pdivmin))
        print("   Neutral atom contribution: {:.2e}".format(self.pneutin))
        print("   Binding energy contribution: {:.2e}".format(self.pbindin))
        print("   Hydrogenic radiation contribution: {:.2e}".format(self.pradhin))
        print("   Impurity radiation contribution: {:.2e}".format(self.pradzin))

        print(" Outer plate;")
        print("   Total: {:.2e}".format(self.pdiviout+self.pdiveout+self.pdivmout,self.pbindout+self.pneutout+self.pradhout+self.pradzout))
        print("   Ion contribution: {:.2e}".format(self.pdiviout))
        print("   Electron contribution: {:.2e}".format(self.pdiveout))
        print("   Molecular contribution: {:.2e}".format(self.pdivmout))
        print("   Neutral atom contribution: {:.2e}".format(self.pneutout))
        print("   Binding energy contribution: {:.2e}".format(self.pbindout))
        print("   Hydrogenic radiation contribution: {:.2e}".format(self.pradhout))
        print("   Impurity radiation contribution: {:.2e}".format(self.pradzout))

        print(" Totals for both plates:")
        print("   Total: {:.2e}".format(self.ptotpart + self.ptotrad))
        print("   Particle fluxes: {:.2e}".format(self.ptotpart))
        print("   Radiation: {:.2e}".format(self.ptotrad))

        if bbb.isupgon[0] == 1:
           print("\nEst. of pdiviout and pdivin from atoms (backscatter) [Watts]:")
           print("   Outer: {:.2e}".format(self.pdivnout))
           print("   Inner: {:.2e}".format(self.pdivnin))

        print("\nPower [W] lost via ionization & recombination radiation is:")
        print("   {:.2e}".format(sum(self.pradiz + self.pradrc)))
        
        print("\nPower [W] lost via recombination only (included in pradht above):")
        print("   {:.2e}".format(sum(self.pradrc)))

        print("\nPower [W] lost at ionization but carried as ion binding-energy:")
        print("   {:.2e}".format(sum(self.pbinde)))
        
        print("\nPower [W] gained by electrons in 3-body recombination (via binding eng):")
        print("   {:.2e}".format(sum(self.pbindrc)))
        
        if bbb.isimpon != 0:
           print("\nPower [W] lost via impurity radiation is:")
           id2=com.nhsp-1
           for id in range(com.nzdf):
              id2 += com.nzsp[id]
              print("  for nuclear charge = {};  Power = {}".format(bbb.znucl[id2], sum(self.pradimp, axis=(0,1,2))[id]))
              print("  for fixed-fraction species;  Power = ", sum(self.pradff))
           
           print("\nElectron Power [W] lost at impur. ioniz. but carried as binding-energy:")
           print("   {:.2e}".format(sum(self.pradzbind)))
        
        print("\nPower [W] lost via dissociation of molecules is:")
        print("   {:.2e}".format(sum(self.prdiss)))
        
        print("\nPower [W] gained by ions from initial Franck-Condon Energy:")
        print("   {:.2e}".format(sum(self.pibirth)))

        print("\nPower [W] lost in parallel momentum exhange via charge exchange:")
        print("   {:.2e}".format(sum(self.pmomv)))
        
        print("\nPower [W] from J.E Joule heating - goes to electrons:")
        print("   {:.2e}".format(sum(bbb.wjdote[1:-1,1:-1])))

        print("\nPower Flow [Watts] incident on Vessel Wall is:")
        print("   Outerwall_sum: {:.2e}".format(self.pwalli + self.pwalle + self.pwallm + self.pwallbd + self.pradhwall + self.pradzwall + self.pneutw))
        print("      pwalli: {:.2e}".format(self.pwalli))
        print("      pwalle: {:.2e}".format(self.pwalle))
        print("      pwallm: {:.2e}".format(self.pwallm))
        print("      pwallbd: {:.2e}".format(self.pwallbd))
        print("      pradhwall: {:.2e}".format(self.pradhwall))
        print("      pradzwall: {:.2e}".format(self.pradzwall))
        print("      pneutw: {:.2e}".format(self.pneutw))

        print("\nPower Flow [Watts] incident on Private Flux Wall is:")
        print("   PFwall_sum: {:.2e}".format(self.ppfi + self.ppfe + self.ppfm + self.ppfbd + self.pneutpf))
        print("      ppfi: {:.2e}".format(self.ppfi))
        print("      ppfe: {:.2e}".format(self.ppfe))
        print("      ppfm: {:.2e}".format(self.ppfm))
        print("      ppfbd: {:.2e}".format(self.ppfbd))
        print("      pneutpf: {:.2e}".format(self.pneutpf))

        print("\nPower Balance: (Pin-Pout)/Pin")
        print("    {:.2e}".format(self.powbal))
        print(50*"=")


    def print_balance_particles(self, redo=True):
        from numpy import sum
        from uedge import bbb, com
        if redo:
            self.balance_particles()
        # particle outflow from separatrix
        self.isephyd = bbb.qe*sum(bbb.fniy[max(com.ixpt1[0]+1,0):com.ixpt2[0]+1,com.iysptrx,0])


        print("\nParticle Flow [Amps] from Core to Scrape-off Layer is:")
        for id in range(com.nfsp):
           if bbb.zi[id] > 0:
              print("  for nuclear charge = {}; ion charge = {}".format(bbb.znucl[id],bbb.zi[id]))
              print("  icore = {}\n".format(self.icore[id]))

        for id in range(com.ngsp):
           print("  for gas species = {}; igascr = {:.2e}".format(id, self.igascr[id]))
        print("  hydrogen ion current at separatrix, isephyd = ", self.isephyd)
        #
        print("\nCurrent from Fixed Volume Source (Sinks) [Amps]:")
        print("  ion current, i_vol = ", self.i_vol)
        #
        print("\nIonization current [Amps] created by re-ionization of gas is:")

        for id in range(com.ngsp):
           print("  for gas species = {}; iion = {:.2e}".format(id, sum(self.iion[id])))
        print("\nRecomb. current [Amps] from ions --> gas by electron-ion recombination:")

        for id in range(com.ngsp):
           print("  for gas species = {}; irecomb = {:.2e}".format(id, sum(self.irecomb[id])))

        print("\nCharge exchange current [Amps] from ions --> gas by iter-species cx:")

        for id in range(com.ngsp):
           print("  for gas species = {}; icxgas = {:.2e}".format(id, sum(self.icxgas[id])))

        print("\nParticle Flow [Amps] incident on Divertor Plate is:")

        for id in range(com.nfsp):
           if bbb.zi[id] > 0:
              print("  for nuclear charge = {}; ion charge = {}".format(bbb.znucl[id], bbb.zi[id]))
              print("  idivin = {:.2e}   idivout = {:.2e}\n".format(self.idivin[id],self.idivout[id]))
         
        print("Neutral Flow [Amps] away from Divertor Plate is:")

        for id in range(com.ngsp):
           print("  for gas species = {}; igasin = {:.2e},  igasout = {:.2e}".format(id,  self.igasin[id], self.igasout[id]))

        print("\nParticle Flow [Amps] incident on Vessel Wall is:")
        for id in range(com.nfsp):
           print("  for ion species = {}; iwall = {:.2e}".format(id, self.iwall[id]))

        for id in range(com.ngsp):
           print("  for gas species = {}; igaswall = {:.2e}".format(id, self.igaswall[id]))


        print("\nParticle Flow [Amps] incident on Private Flux Wall is:")
        for id in range(com.nfsp):
           print("  for ion species = {}; ipf = {:.2e}".format(id, self.ipf[id]))

        for id in range(com.ngsp):
           print("  for gas species = {}; igaspf = {:.2e}".format(id, self.igaspf[id]))


        print("\nParticle Balance for Neutrals: (Iinput-Itotal)/Iinput")
        print("    {:.2e}".format(self.deligas))
        print(50*"=")
        #

    def print_peak_values(self, redo = True):
        from uedge import bbb, com
        if redo:
            self.balance_power()
            self.balance_particles()
        print("\nPeak electron temperatures [eV] on divertor plates")
        print("    Outer: {:.2e}".format(max(bbb.te[com.nx,1:com.ny+1])/bbb.ev))
        print("    Inner: {:.2e}".format(max(bbb.te[1,1:com.ny+1])/bbb.ev))
        #
        print("\nPeak ion temperatures [eV] on divertor plates")
        print("    Outer: {:.2e}".format(max(bbb.ti[com.nx,1:com.ny+1])/bbb.ev))
        print("    Inner: {:.2e}".format(max(bbb.ti[1,1:com.ny+1])/bbb.ev))
        #
        print("\nPeak power flux [MW/m**2] on divertor plates")
        print("    Outer: {:.2e}".format(1.e-6*max(self.sdtout[1:com.ny+1])))
        print("    Inner: {:.2e}".format(1.e-6*max(self.sdtin[1:com.ny+1])))
        #
        print("\nPeak ion densities [m**(-3)] on divertor plates")
        print("    Outer: {:.2e}".format(max(bbb.ni[com.nx,1:com.ny+1,0])))
        print("    Inner: {:.2e}".format(max(bbb.ni[1,1:com.ny+1,0])))
        #
        print("\nPeak gas densities [m**(-3)] on divertor plates")
        print("    Outer: {:.2e}".format(max(bbb.ng[com.nx,1:com.ny+1,0])))
        print("    Inner: {:.2e}".format(max(bbb.ng[1,1:com.ny+1,0])))
        print(50*"=")

    
    def balancee_new(self):
        self.engbal()
        self.pradpltwl()
        self.balance_power(redo=False)
        self.balance_particles(redo=False)
        self.print_balance_power(redo=False)
        self.print_balance_particles(redo=False)
        self.print_peak_values(redo=False)

 
    def balancee(self):
        self.engbal()
        self.pradpltwl()
        self.balance_power(redo=False)
        self.balance_particles(redo=False)
        self.print_balance_power(redo=False)
        self.print_balance_particles(redo=False)
        self.print_peak_values(redo=False)

        
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
    from numpy import sum
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import contextlib

    plt.ion()
    

    new = UeBalance()


    with open('new_balancee.txt', 'w') as f:
        with contextlib.redirect_stdout(f):
            new.balancee_new()

    with open('old_balancee.txt', 'w') as f:
        with contextlib.redirect_stdout(f):
            new.balancee()
    '''
    new.balancee_new()
    '''

    bbb.engbal(1)
    bbb.pradpltwl()
    new.engbal()
    new.calc_engerr(1, False)
    new.pradpltwl()

    for var in ['fetx', 'fety', 'engerr', 'peirad',
        'pmloss', 'pmpot', 'pmrada', 'pmradm', 'pmomv'
    ]:
        denom = bbb.__getattribute__(var)
        denom[denom==0] = 1e-100
        dif = abs( (new.__dict__[var] - bbb.__getattribute__(var))/ denom)
        dif[denom == 1e-100] = 0
        if dif.max() > 1e-6:
            print(var,  dif.max())
            f, ax = plt.subplots()
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            im = ax.imshow(dif, cmap='viridis')
            f.colorbar(im, cax=cax, orientation='vertical')
            ax.set_title(var)

    varlist = [
        'iion', 'irecomb', 'icxgas', 'pradrc', 'pradiz', 'pbinde', 'pbindrc', 'prdiss', 'pibirth', 'pradff'
    ]
    if (bbb.isimpon > 2):
        varlist += ['pradimp']
    for var in varlist:
        if var == 'pradff':
            dif = abs(sum(new.__dict__[var]) - bbb.pradfft)/max(bbb.pradfft,1e-100)
        else:
            denom = sum(bbb.__getattribute__(var))
            if len(denom.shape)>1:
                denom[denom==0] = 1e-100
            else:
                denom += 1e-100
            dif = abs(sum(new.__dict__[var]) - sum(bbb.__getattribute__(var)))/ denom
            if len(denom.shape)>1:
                dif[denom==0] = 0
            else:
                dif -= 1e100
        if dif > 1e-6:
            print(var, dif)
    
    for var in ['pwr_plth', 'pwr_pltz', 'pwr_wallh', 'pwr_wallz', 'pwr_pfwallh', 'pwr_pfwallz' 
    ]:
        denom = bbb.__getattribute__(var)
        denom[denom==0] = 1e-100
        dif = abs( (new.__dict__[var] - bbb.__getattribute__(var))/ denom)
        dif[denom == 1e-100] = 0
        if dif.max() > 1e-6:
            print(var,  dif.max())
            f, ax = plt.subplots()
            colors = ['k', 'r', 'b', 'c', 'm', 'g']
            for ix in range(bbb.__getattribute__(var).shape[-1]):
                ax.plot(dif[:, ix], color=colors[ix])
#            divider = make_axes_locatable(ax)
#            cax = divider.append_axes('right', size='5%', pad=0.05)
#            im = ax.imshow(dif, cmap='viridis')
#            f.colorbar(im, cax=cax, orientation='vertical')
#            ax.set_title(var)





