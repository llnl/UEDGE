



class UeBalance():
    def __init__(self):
        from numpy import zeros, array
        from uedge import com

        yy = []
        y = list(range(com.ny+2))
        for _ in range(com.nx+2):
            yy.append(y)
        self.yy=array(yy)
        return


    def engbal(self):
        """ Calculates various components of the 2-D energy flow and the 
        ionization and radiation for use in the postprocessing file
        balancee to determine energy balance; these 2-D loops become 
        expensive when done from the parser.
        """
        from uedge import bbb, com, aph
        from copy import deepcopy
        from numpy import zeros, cos

        for var in [
            'fetx', 'fety', 'engerr', 'pmloss', 'pmrada', 'pmradm', 'pmpot',
            'peirad', 'pmomv', 'pradrc', 'pradiz', 'prdiss',
            'pibirth', 'pbinde', 'pbindrc', 'pradzbind', 'pradff'
        ]:
            self.__dict__[var] = zeros((com.nx+2, com.ny+2))
        for var in ['icxgas', 'iion', 'irecomb']:
            self.__dict__[var] = zeros((com.nx+2, com.ny+2, com.ngsp))
        self.pradimp = zeros((com.nx+2, com.ny+2, sum(com.nzsp)+1, max(sum(com.nzsp),1)))

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
        up = bbb.up
        upip1 = bbb.upi[bbb.ixp1, self.yy]
        upim1 = bbb.upi[bbb.ixm1, self.yy]
        if bbb.isimpon > 2:
            prad = deepcopy(bbb.prad)
            pwrze = deepcopy(bbb.pwrze)
        else:
            prad = zeros(bbb.ne.shape)
            pwrze = zeros(bbb.ne.shape)

        # TODO: remove nested loops and replace with indirect referencing
        # Set arrays to check energy conserv; add ion parallel drift and visc heat
        thetaix = 0.5*(com.angfx[bbb.ixm1, self.yy] + com.angfx)
        thetaix2 = 0.5*(com.angfx[bbb.ixp1, self.yy] + com.angfx)
        eta_dup2dy = 0.25*(
            (bbb.visy[:,1:]*(upi[bbb.ixm1[:,:-1], self.yy[:,1:]] + upi[:,1:])**2) \
#           NOTE: unsure whether the iy-indices offset causes an X-point glitch in original 
#           implementation or not...
#            (bbb.visy*(upim1 + upi)**2)[:,1:] \
            - (bbb.visy*(upim1 + upi)**2)[:,:-1])
        for ii in range(com.nusp):
            self.fety[:,:-1] += (bbb.mi[ii]/32)*(
                (upim1 + upi)[:,:-1,ii] + \
                    (upi[bbb.ixm1[:,:-1], self.yy[:,1:],ii] + upi[:,1:,ii])\
#                    (upim1 + upi)[:,1:,ii]\
                    )**2*bbb.fniy[:,:-1,ii] \
                - (bbb.cfvisy*0.5*com.sy[:,:-1]*com.gyf[:,:-1]*eta_dup2dy[:,:,ii])

            self.fetx[:,:-1] += 0.5*bbb.mi[ii]*(upi**2*bbb.fnix)[:,:-1,ii] \
                - bbb.cfvisx*0.25*com.sx[:,:-1]*(
                    bbb.visx[:,:-1,ii]*com.gx[:,:-1]*cos(thetaix[:,:-1])*(\
                        upi**2 - upim1**2)[:,:-1,ii] \
                    + bbb.visx[bbb.ixp1,self.yy][:,:-1,ii]\
                        *com.gx[bbb.ixp1,self.yy][:,:-1]*cos(thetaix2[:,:-1])*( \
                        upip1**2 - upi**2)[:,:-1,ii]) \
                - (upi[:,:-1,ii]*bbb.fmixy[:,:-1,ii])
        self.fety[:,:-1] += (bbb.feey + bbb.feiy)[:,:-1]
        self.fetx[:,:-1] += (bbb.feex + bbb.feix)[:,:-1]


        # Correct the boundary x-fluxes if non-unity ckinfl
        if (abs(bbb.ckinfl - 1) > 1e-10):
            for jx in range(com.nxpt):
                ixt = com.ixlb[jx]
                ixr = com.ixrb[jx]
                self.fetx[ixt] = 0
                self.fetx[ixr] = 0
                for ii in range(com.nfsp):
                    self.fetx[ixt] +=  \
                        0.5*bbb.mi[ii]*(bbb.fnix*up**2)[ixt,:,ii] \
                        - bbb.ckinfl*0.5*com.sx[ixt]*bbb.visx[ixt+1,:,ii] \
                        *com.gx[ixt+1]*( up[ixt,:,ii]**2 - up[ixt+1,:,ii]**2)

                    self.fetx[ixr] +=  \
                        0.5*bbb.mi[ii]*(bbb.fnix*up**2)[ixr,:,ii] \
                        - bbb.ckinfl*0.5*com.sx[ixr]*bbb.visx[ixr,:,ii] \
                        *com.gx[ixr]*( up[ixr,:,ii]**2 - up[ixr-1,:,ii]**2)
                self.fetx[ixt] += (bbb.feex + bbb.feix)[ixt]
                self.fetx[ixr] += (bbb.feex + bbb.feix)[ixr]


        self.pmloss = (1-bbb.ismolcrm)*bbb.cnsor*( \
                bbb.ediss*bbb.ev*(0.5*bbb.psordis[:,:,1]) \
                + bbb.ceisor*bbb.eion*bbb.ev*(bbb.psordis[:,:,1])
            ) + bbb.ismolcrm*bbb.cnsor*( \
                bbb.cmesori*(bbb.emolia[:,:,0] + bbb.emolia[:,:,1]) \
                + bbb.cmesore*bbb.edisse
            )
    
        sv_crumpet = {
            22: zeros(bbb.ne.shape),
            23: zeros(bbb.ne.shape),
            24: zeros(bbb.ne.shape),
        }
        for ix in range(com.nx+2):
            for iy in range(com.ny+2):
                for ii in [22, 23, 24]:
                    sv_crumpet[ii][ix,iy] = \
                        aph.sv_crumpet(bbb.te[ix,iy], bbb.ne[ix,iy], ii)

        if bbb.ishymol:
            self.pmpot = bbb.ismolcrm*bbb.ng[:,:,1]*com.vol*sv_crumpet[22]
            self.pmrada = bbb.ismolcrm*bbb.ng[:,:,1]*com.vol*sv_crumpet[23]
            self.pmradm = bbb.ismolcrm*bbb.ng[:,:,1]*com.vol*sv_crumpet[24]
            # Here peirad includes sum of electron and ion energy losses; note that binding
            # energy is included in eeli term, but it is carried by the ions.
            # Note also that eion and ediss generally balance in the next line
            # because should have ediss=2*eion - transfer from electron to ion energy

        self.peirad = bbb.cnsor*(
            bbb.erliz + bbb.erlrc + self.pmloss \
            + bbb.ebind*bbb.ev*(bbb.psor[:,:,0] - bbb.psorrg[:,:,0]) 
        )

        if not bbb.isupgon[0]:
            self.pmomv = bbb.cngmom[0]*bbb.up[:,:,0]*com.sx*com.rrv*( \
                    (bbb.ng*bbb.tg)[bbb.ixp1,self.yy,0] - (bbb.ng*bbb.tg)[:,:,0]
                ) + bbb.cmwall[0]*0.125*bbb.mi[0]*(
                    up[:,:,0] + up[bbb.ixm1,self.yy,1]
                )**2*bbb.ng[:,:,0]*bbb.nucx[:,:,0]*com.vol
                    
        if ((bbb.ishymol == 0) or (ig != 1)):
            self.iion -= bbb.cnsor*bbb.qe*bbb.psorg
            self.irecomb -= bbb.cnsor*bbb.qe*bbb.psorrg
            self.icxgas -= bbb.qe*bbb.psorcxg
        self.pradrc += bbb.cnsor*bbb.erlrc
        self.pradiz += (bbb.eeli - bbb.ebind*bbb.ev) * bbb.psor[:,:,0]
        self.pbinde += bbb.ebind*bbb.ev*bbb.psor[:,:,0]
        self.pbindrc += bbb.ebind*bbb.ev*bbb.psorrg[:,:,0]
        self.prdiss += (1-bbb.ismolcrm)*(bbb.ediss*bbb.ev*(0.5*bbb.psordis[:,:,1])) \
                + bbb.ismolcrm*bbb.cmesore*bbb.edisse
        self.pibirth += (1-bbb.ismolcrm) * ( \
            bbb.ceisor*bbb.eion*bbb.ev*bbb.psordis[:,:,1] \
            - bbb.ccoldsor*bbb.ng[:,:,0]*( \
                1.5*bbb.ti - bbb.eion*bbb.ev \
            )*bbb.nucx[ix,iy,0]*com.vol) + bbb.ismolcrm*( \
                bbb.ceisor*bbb.cmesore*(  \
                    bbb.emolia[:,:,0] + bbb.emolia[:,:,1]
                ) - bbb.ccoldsor*bbb.ng[:,:,0]*(
                    1.5*bbb.ti - bbb.eion*bbb.ev
                )*bbb.nucx[:,:,0]*com.vol
            )
        if bbb.isimpon in [2,7]:
            self.pradff += bbb.pradcff*com.vol
        if bbb.isimpon > 2:
            for jz in range(com.ngsp - com.nhgsp):
                for iimp in range(com.nzsp[jz]):
                    ii = iimp + com.nhgsp + sum(com.nzsp[:jz]) - 1
                    self.pradimp[:,:, ii+1, jz] += bbb.pradz[:,:,ii,jz]*com.vol
                self.pradimp[:,:,0,jz] += bbb.pradz[:,:, com.nzsp[jz], jz]*com.vol
        self.pradzbind = (pwrze - prad)*com.vol

                           
    def calc_engerr(self, pwrin=1, redo=True):
        from uedge import bbb, com
        from numpy import zeros

        self.engerr = zeros((com.nx+2, com.ny+2))
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
        from numpy import zeros, arctan2, pi, cos, sum, minimum
        from copy import deepcopy
        nj = com.nxomit 
#        if bbb.isimpon > 2:
        prdu = deepcopy(bbb.prad)
#        else:
#            prdu = zeros(bbb.ne.shape)
        

        # Initialize arrays for subsequent runs
        self.pwr_plth = zeros((com.ny+2, 2*com.nxpt))
        self.pwr_pltz = zeros((com.ny+2, 2*com.nxpt))
        self.pwr_wallh = zeros((com.nx+2))
        self.pwr_wallz = zeros((com.nx+2))
        self.pwr_pfwallh = zeros((com.nx+2))
        self.pwr_pfwallz = zeros((com.nx+2))

        thetapl1_arr = zeros((com.ny+2, com.nx+2, com.ny+2))
        thetapl2_arr = zeros((com.ny+2, com.nx+2, com.ny+2))
        # TODO: remove nested loops and replace with indirect referencing
        for ip in range(2*com.nxpt):
            if (ip % 2) == 0: # Even numbers
                ixv = com.ixlb[int( ip > 1) + int(ip > 3)]
            else:
                ixv = com.ixrb[int( ip > 1) + int(ip > 3)]+1
            for iyv in range(1, com.ny+1):
                thetapl1_arr[iyv, 1:-1, 1:-1] = arctan2( \
                    com.zm[ixv+nj, iyv, 1] - com.zm[1:-1, 1:-1, 0],
                    com.rm[ixv+nj, iyv, 1] - com.rm[1:-1, 1:-1, 0],
                )
                thetapl2_arr[iyv, 1:-1, 1:-1] = arctan2( \
                    com.zm[ixv+nj, iyv, 3] - com.zm[1:-1, 1:-1, 0],
                    com.rm[ixv+nj, iyv, 3] - com.rm[1:-1, 1:-1, 0],
                )
            dthgy_arr = abs(thetapl1_arr - thetapl2_arr)
            frth = minimum(dthgy_arr, 2*pi - dthgy_arr)/2/pi
            sxo = com.sx[ixv]/cos(com.angfx[ixv])
            for iyv in range(1,com.ny+1):
                self.pwr_pltz[iyv, ip] += sum(
                    (prdu*com.vol*frth[iyv]/sxo[iyv])[nj:]
                )
                self.pwr_plth[iyv, ip] += sum(  
                    ((
                        (bbb.eeli - bbb.ebind*bbb.ev)*bbb.psor[:,:,0] \
                        + bbb.erlrc
                    )*frth[iyv]/sxo[iyv])[nj:]
                )
            # Set corner values
            self.pwr_pltz[0,ip] = self.pwr_pltz[1,ip]
            self.pwr_pltz[com.ny+1,ip] = self.pwr_pltz[com.ny,ip]
            self.pwr_plth[0,ip] = self.pwr_plth[1,ip]
            self.pwr_plth[com.ny+1,ip] = self.pwr_plth[com.ny,ip]
    


        thetapf1_arr = zeros((com.nx+2, com.nx+2, com.ny+2))
        thetapf2_arr = zeros((com.nx+2, com.nx+2, com.ny+2))
        thetaw1_arr = zeros((com.nx+2, com.nx+2, com.ny+2))
        thetaw2_arr = zeros((com.nx+2, com.nx+2, com.ny+2))
        for ixv in range(1, com.nx+1):
            thetapf1_arr[ixv, 1:-1, 1:-1] = arctan2( \
                com.zm[ixv+nj, 0, 1] - com.zm[1:-1, 1:-1, 0],
                com.rm[ixv+nj, 0, 1] - com.rm[1:-1, 1:-1, 0],
            )
            thetapf2_arr[ixv, 1:-1, 1:-1] = arctan2( \
                com.zm[ixv+nj, 0, 2] - com.zm[1:-1, 1:-1, 0],
                com.rm[ixv+nj, 0, 2] - com.rm[1:-1, 1:-1, 0],
            )
            thetaw1_arr[ixv, 1:-1, 1:-1] = arctan2( \
                com.zm[ixv+nj, -1, 1] - com.zm[1:-1, 1:-1, 0],
                com.rm[ixv+nj, -1, 1] - com.rm[1:-1, 1:-1, 0],
            )
            thetaw2_arr[ixv, 1:-1, 1:-1] = arctan2( \
                com.zm[ixv+nj, -1, 2] - com.zm[1:-1, 1:-1, 0],
                com.rm[ixv+nj, -1, 2] - com.rm[1:-1, 1:-1, 0],
            )
        dthgz = abs(thetaw1_arr - thetaw2_arr)
        dthgx = abs(thetapf1_arr - thetapf2_arr)
        frthw = minimum(dthgz, 2*pi - dthgz)/2/pi
        frthpf = minimum(dthgx, 2*pi - dthgx)/2/pi
        for ixv in range(1, com.nx+2):
            self.pwr_wallz[ixv] += sum(
                (prdu*com.vol)*frthw[ixv]/com.sy[ixv,-1]
            )
            self.pwr_wallh[ixv] += sum( ( \
                (bbb.eeli - bbb.ebind*bbb.ev)*bbb.psor[:,:,0] \
                +bbb.erlrc)*frthw[ixv]/com.sy[ixv,-1])
            self.pwr_pfwallz[ixv] += sum(
                (prdu*com.vol)*frthpf[ixv]/com.sy[ixv,0]
            )
            self.pwr_pfwallh[ixv] += sum( ( \
                (bbb.eeli - bbb.ebind*bbb.ev)*bbb.psor[:,:,0] \
                +bbb.erlrc)*frthpf[ixv]/com.sy[ixv,0])

        self.pwr_wallz[0] = self.pwr_wallz[1]	# Because prad(0,) is not calculated
        self.pwr_wallz[com.nx+1] = self.pwr_wallz[com.nx]
        self.pwr_wallh[0] = self.pwr_wallh[1]
        self.pwr_wallh[com.nx+1] = self.pwr_wallh[com.nx]
        for ip in range(com.nxpt):
            ixc = slice(com.ixpt1[ip]+1, com.ixpt2[ip]+1)
            self.pwr_pfwallh[ixc] = 0
            self.pwr_pfwallz[ixc] = 0
            self.pwr_pfwallh[com.ixlb[ip]] = self.pwr_pfwallh[com.ixlb[ip]+1]
            self.pwr_pfwallh[com.ixrb[ip]+1] = self.pwr_pfwallh[com.ixrb[ip]]
            self.pwr_pfwallz[com.ixlb[ip]] = self.pwr_pfwallz[com.ixlb[ip]+1]
            self.pwr_pfwallz[com.ixrb[ip]+1] = self.pwr_pfwallz[com.ixrb[ip]]



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


        if redo:
            self.engbal()
            self.pradpltwl()
        """ INITIALIZE VARIABLES """
        for var in [
            'pradhpf', 'pradzpf', 
            'pdiviin', 'pdiviout', 'pdivein', 'pdiveout', 'pdivmin', 'pdivmout',
            'pbindout', 'pbindin', 'pdivnout', 'pdivnin', 'pbindzout', 
            'pbindzin', 'pneutout', 'pneutin', 'pradhout', 'pradhin',
            'pradzin', 'pradzout', 'pwalli', 'psepi', 'psepe', 'psepbd',
            'pbcorei', 'pbcoree', 'pbcorebd'
        ]:
            self.__dict__[var] = 0

        for var in ['sneutpf', 'sneutw', 'ppfi', 'ppfe', 'pwalle', 'pwallm',
            'ppfm', 'pwallbd', 'ppfbd', 'pneutpf'
            'pwalli', 'pwalle', 'pwallm', 'pwallbd', 'pradhwall', 'pradzwall',
            'pneutw', 'pwalli', 'ppfi', 'ppfe', 'ppfm', 'ppfbd', 'pneutpf'
        ]:
            self.__dict__[var] = zeros((com.nx+2,))

        for var in [
            'sdrout', 'sdrin', 'sdeout', 'sdmout',
            'sdmin', 'sdein', 'sdtout', 'sdtin', 'sbindout', 'sbindin',
            'sbindzout', 'sbindzin', 'sneutout', 'sneutin',
            'sdnout', 'sdnin'
        ]:
            self.__dict__[var] = zeros((com.ny+2,))

        for var in ['sdiout', 'sdiin']:
            self.__dict__[var] = zeros((com.ny+2,com.nfsp))

          
        # note: ixi=ixdivin has been set to 1 to allow velocity derivatives to be calc.

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

        """ POWER ACROSS SEPARATRIX """
        iy = com.iysptrx
        self.psepi += fluxfacy*(
            bbb.feiy[:,iy] + (bbb.mi[0]/32)*(\
                (bbb.up[bbb.ixm1, self.yy] + bbb.up)[:,iy,0] \
                + (bbb.up[bbb.ixm1, self.yy] + bbb.up)[:,iy+1,0] \
            )**2*bbb.fniy[:,iy,0] - bbb.cfvisy*0.125*com.sy[:,iy] \
            * have( 
                (bbb.visy[:,:,0]*com.gy)[:,iy], 
                (bbb.visy[:,:,0]*com.gy)[:,iy+1]
            ) * ( \
                (bbb.up[bbb.ixm1,self.yy]+bbb.up)[:,iy+1,0]**2 \
                - (bbb.up[bbb.ixm1, self.yy]+bbb.up)[:,iy,0]**2
        ))
        self.psepe += fluxfacy*bbb.feey[:,iy]
        self.psepbd += fluxfacy*bbb.fniy[:,iy,0]*bbb.ebind*bbb.ev

        """ POWER ACROSS CORE BOUNDARY """
        iyc = iycore
        self.pbcorei += fluxfacy*(
            bbb.feiy[:,iyc] + (bbb.mi[0]/32)*(
                (bbb.up[bbb.ixm1, self.yy] + bbb.up)[:,iyc,0] \
                + (bbb.up[bbb.ixm1,self.yy] + bbb.up)[:, iyc+1,0] \
            )**2 * bbb.fniy[:,iyc,0] - bbb.cfvisy*0.125*com.sy[:,iyc]\
            * have(
                (bbb.visy[:,:,0]*com.gy)[:,iyc],
                (bbb.visy[:,:,0]*com.gy)[:,iyc+1]
            ) * ( \
                (bbb.up[bbb.ixm1,self.yy]+bbb.up)[:,iyc,0]**2 \
                - (bbb.up[bbb.ixm1,self.yy]+bbb.up)[:,iy,0]**2
        ))
        self.pbcoree += fluxfacy*bbb.feey[:,iyc]
        self.pbcorebd += fluxfacy*bbb.fniy[:,iyc,0]*bbb.ebind*bbb.ev

        for var in ['psepi', 'pbcoree', 'pbcorebd', 'psepbd', 'psepe',
            'pbcorei']:
            self.__dict__[var] = sum(self.__dict__[var]) 

        """ POWER INPUT FROM BIAS SUPPLY """
        # NOTE: Generalization to dnull by AH 25/04/01
        self.p_bias = 0
        for nx in range(com.nxpt):
            self.p_bias -= sum(
                bbb.fqx[com.ixrb[nx],1:com.ny+1]*(\
                        bbb.phi0r[1:com.ny+1, nx] - bbb.phi0l[1:com.ny+1, nx]
            ))

        """ TOTAL BIAS CURRENT """
        self.i_bias = sum(bbb.fqx[com.nx,1:com.ny+1])

        """ POWER FROM FIXED VOLUME SOURCES """
        self.p_e_vol = sum(bbb.pwrsore)
        self.p_i_vol = sum(bbb.pwrsori)

        if bbb.l_parloss <= 1e9:
            self.p_e_vol -= sum(\
                bbb.nuvl[:,:,0]*com.vol*bbb.ne*bbb.bcee*bbb.te
            )
            self.p_i_vol -= sum(
                bbb.nuvl[:,:,0]*com.vol*bbb.ni[:,:,0]*bbb.bcei*bbb.ti
            )

    
        """ POWER INCIDENT ON THE DIVERTOR PLATES """
        # Allow use of "old" or "new" switches for neutral energy loss


        # First get radiation power in pwr_plth and pwr_pltz
        self.sdrin = self.pwr_plth[:,0]+self.pwr_pltz[:,0]
        self.sdrout = self.pwr_plth[:,1]+self.pwr_pltz[:,1]


        #two ifs to be able to used old executables
        ixdivin=0
        iybegin=1
        if(com.islimon == 1) and (com.nyomitmx != 0):
           com.nx = com.ix_lim-1
           ixdivin = com.ix_lim+1
           iybegin = com.iy_lims
        ixi=ixdivin
        ixo=com.nx
        ixineut=1
        divsl = slice(iybegin, com.ny+1)
        
        sxo = (com.sx/cos(com.angfx))[com.nx, divsl]
        sxi = (com.sx/cos(com.angfx))[ixdivin, divsl]
        vxno = 0.25*sqrt(8*bbb.tg[com.nx,divsl,0]/(pi*bbb.mg[0]))
        vxni = 0.25*sqrt(8*bbb.tg[1,divsl,0]/(pi*bbb.mg[0]))
        # KINETIC POWER
        for id in range(com.nfsp):
            if bbb.zi[id] > 0:
                self.sdiout[divsl,id] += 0.5*bbb.mi[id]\
                    *bbb.upi[com.nx,divsl,id]**2*bbb.fnix[com.nx,divsl,id]/sxo
                self.sdiin[divsl,id] -= 0.5*bbb.mi[id]\
                    *bbb.upi[ixdivin,divsl,id]**2*bbb.fnix[ixdivin,divsl,id]/sxi
            else:
                if bbb.ishymol:
                    self.sdiout[divsl,id] = 0.5*bbb.mi[id]*abs(\
                        bbb.up[com.nx,divsl,id]**3*bbb.ni[com.nx,divsl,id])/\
                        com.rrv[com.nx,divsl]
                    self.sdiin[divsl,id] = 0.5*bbb.mi[id]*abs(\
                        bbb.up[ixdivin,divsl,id]**3*bbb.ni[ixdivin,divsl,id])/\
                        com.rrv[ixdivin,divsl]
                else:
                    self.sdiout[divsl,id] -= 0*(0.5*bbb.mi[id] \
                        *bbb.up[com.nx,divsl,id]**2*bbb.fnix[com.nx, divsl,id] \
                        )/sxo
                    self.sdiin[divsl,id] -= 0*(0.5*bbb.mi[id] \
                        *bbb.up[ixdivin,divsl,id]**2*bbb.fnix[ixdivin, divsl,id] \
                        )/sxi
        # VISCOUS POWER
        # AH: bbb.ckinfl treated as array, although it is double. Legacy switch?
        for id in range(com.nusp):
            self.sdiout[divsl,id] -= bbb.ckinfl*0.5*(
                    com.sx*bbb.visx[:,:,id]*com.gx
                )[com.nx, divsl] * (
                    bbb.up[com.nx,divsl,id]**2 - bbb.up[com.nx-1, divsl,id]**2
                ) / sxo
            self.sdiin[divsl,id] += bbb.ckinfl*0.5 \
                    * com.sx[ixdivin,divsl]*bbb.visx[ixdivin+1,divsl,id]\
                    * com.gx[ixdivin+1, divsl] * (
                    bbb.up[ixdivin+1,divsl,id]**2 - bbb.up[ixdivin, divsl,id]**2
                ) / sxi
        # CONVECTIVE-CONDUCTIVE POWER
        self.sdiout[divsl,0] += bbb.feix[com.nx,divsl]/sxo
        self.sdiin[divsl,0]  -= bbb.feix[ixdivin,divsl]/sxi 
        # BINDING POWER
        self.sbindout[divsl] += bbb.fnix[com.nx,divsl,0]*bbb.ebind*bbb.ev / sxo
        self.sbindin[divsl] -= bbb.fnix[ixdivin,divsl,1] * bbb.ebind*bbb.ev / sxi
        # NOTE: Should probably be index 0
        # MOLECULAR POWER
        if ishymoleng==1:  #mol heat flux; drift eng small,<0
             self.sdmout[divsl] += bbb.fegx[com.nx,divsl,1]/sxo
             self.sdmin[divsl] -= bbb.fegx[ixdivin,divsl,1]/sxi
        # IMPURITY BINDING ENERGY
        for id in range(com.nzdf):
            if (id == 0):
                id2min = com.nhsp
                id2max = id2min + com.nzsp[id] - 1
            else:
                id2min = com.nhsp+sum(com.nzsp[1:id-1])
                id2max = id2min + com.nzsp[id]-1
            for id2 in range(id2min, id2max+1):
                for id3 in range(bbb.znucl[id]):
                    ebindz = bbb.ebindz(id3, bbb.znucl[id2])
                    self.sbindzout[divsl] += bbb.fnix[com.nx,divsl,id2] * \
                        ebindz*bbb.ev/sxo
                    self.sbindzin[divsl] -= bbb.fnix[ixdivin,divsl,id2] * \
                        ebindz*bbb.ev/sxi
        # NEUTRAL CONVECTION
        # TODO: replcae this with convective-conductive contribution?
        self.sneutout[divsl] = cenggpl*2.*vxno\
                *bbb.ng[com.nx,divsl,0]*bbb.tg[com.nx,divsl,0]
        self.sneutin[divsl] = cenggpl*2.*vxni\
                *bbb.ng[ixineut,divsl,0]*bbb.tg[ixineut,divsl,0]
        # ELECTRON POWER
        self.sdeout[divsl] += ( 
                bbb.feex[com.nx,divsl] \
                +bbb.fqx[com.nx,divsl]*(bbb.phi[com.nx,divsl]-bbb.phi0r[divsl,0]) \
            ) / sxo 
        self.sdein[divsl]  -= ( 
                bbb.feex[ixdivin,divsl] \
                + .001*bbb.fqx[ixdivin,divsl]*(bbb.phi[ixdivin,divsl]-bbb.phi0l[divsl,0])\
            )/sxi
        # NOTE: No idea why there is a factor of 1e-3 here...
        # TOTAL DIVERTOR POWER
        self.sdtout += self.sdeout + sum(self.sdiout, axis=(1)) \
                + self.sbindout + self.sdmout + self.sbindzout[iy] \
                + self.pwr_plth[:,1] + self.pwr_pltz[:,1]
        self.sdtin += self.sdein + sum(self.sdiin, axis=(1)) \
                + self.sdmin + self.sbindin + self.sbindzin \
                + self.pwr_plth[:,0] + self.pwr_pltz[:,0]

        # INTEGRATE OVER PLATES
        self.pdiviout += sum(sum(self.sdiout[divsl], axis=(1))*sxo) 
        self.pdiveout += sum(self.sdeout[divsl]*sxo)
        self.pdivmout += sum(self.sdmout[divsl]*sxo)
        self.pdiviin  += sum(sum(self.sdiin[divsl], axis=(1))*sxi) 
        self.pdivein  += sum(self.sdein[divsl]*sxi)
        self.pdivmin  += sum(self.sdmin[divsl]*sxi)
        self.pbindout += sum(self.sbindout[divsl]*sxo)
        self.pbindzout += sum(self.sbindzout[divsl]*sxo)
        self.pbindin += sum(self.sbindin[divsl]*sxi)
        self.pbindzin += sum(self.sbindzin[divsl]*sxi)
        self.pneutout += 0*sum(self.sneutout[divsl]*sxo)  # included in self.pdiviout
        self.pneutin += 0.*sum(self.sneutin[divsl]*sxo)     # included in self.pdiviin
        self.pradhout += sum(self.pwr_plth[divsl,1]*sxo)
        self.pradzout += sum(self.pwr_pltz[divsl,1]*sxo)
        self.pradhin += sum(self.pwr_plth[divsl,0]*sxi)
        self.pradzin += sum(self.pwr_pltz[divsl,0]*sxi)
        # APPROXIMATE NEUTRAL THERMAL FLUXES  
        # NOTE: add convective-conductive flows here?
        if bbb.isupgon[0] == 1: # Approx. neutral energy flux
            self.sdnout[divsl] += self.sdiout[divsl,1] + ( \
                    com.sx[com.nx,divsl]*bbb.hcxn[com.nx,divsl]\
                    * com.gxf[com.nx,divsl]*(\
                        bbb.ti[com.nx,divsl]-bbb.ti[com.nx+1,divsl]\
                    ) + 2.5*bbb.fnix[com.nx,divsl,1]*bbb.ti[com.nx+1,divsl] \
            ) / sxo
            self.sdnin[divsl] += self.sdiin[divsl,1] + ( \
                    com.sx[ixdivin,divsl]*bbb.hcxn[ixdivin,divsl]\
                    * com.gxf[ixdivin,divsl]*(\
                        bbb.ti[ixdivin,divsl]-bbb.ti[ixdivin+1,divsl]\
                    ) + 2.5*bbb.fnix[ixdivin,divsl,1]*bbb.ti[ixdivin,divsl]\
            ) / sxi
            self.pdivnout += sum(self.sdnout[divsl]*sxo)
            self.pdivnin += sum(self.sdnin[divsl]*sxi)

        # Fix up boundary values
        for var in ['sdtin', 'sdein', 'sdiin', 'sdtout', 'sdeout', 'sdiout']:
            self.__dict__[var][0] = self.__dict__[var][1]
            self.__dict__[var][-1] = self.__dict__[var][-2]

        #
        iywall=com.ny       # DEFINITION
        iypf=0
        #
        # power flow to vessel and private flux wall

        vynw = 0.25*sqrt(8*bbb.tg[:,iywall+1,0]/(pi*bbb.mg[0]))
        self.pwalli += fluxfacy*( \
                bbb.feiy[:,iywall] + 0.125*bbb.mi[0]*(\
                    bbb.up[bbb.ixm1,self.yy] + bbb.up  \
                )[:,iywall,0]**2 * bbb.fniy[:,iywall,0] \
                - bbb.cfvisy*0.125*com.sy[:,iywall]* have( \
                        (bbb.visy[:,:,0]*com.gy)[:,iywall],\
                        (bbb.visy[:,:,0]*com.gy)[:,iywall+1] \
                ) * ( \
                    (bbb.up[bbb.ixm1, self.yy] + bbb.up)[:,iywall+1, 0]**2 \
                    - (bbb.up[bbb.ixm1, self.yy] + bbb.up)[:,iywall, 0]**2 \
        ) )
        self.pwalle += fluxfacy*bbb.feey[:,iywall]
        self.pwallbd += fluxfacy*bbb.fniy[:,iywall,0]*bbb.ebind*bbb.ev
        self.pradhwall += fluxfacy*self.pwr_wallh*com.sy[:,iywall]
        self.pradzwall += fluxfacy*self.pwr_wallz*com.sy[:,iywall]
        self.sneutw += cenggw*2.*vynw*bbb.tg[:,com.ny+1,0]*com.sx[:,com.ny]
        self.pneutw += self.sneutw
        if ishymoleng == 1:  #molec temp eqn active
            self.pwallm += bbb.fegy[:,iywall,1]



        for xsl in [slice(0,com.ixpt1[0]+1), slice(com.ixpt2[0]+1, com.nx+1)]:
#           NOTE: tg species index was 1/0, assuming it should be 0
            vynpf = 0.25*sqrt(8*bbb.tg[xsl,0,0]/(pi*bbb.mg[0]))
            self.ppfi[xsl] -= fluxfacy*( \
                    bbb.feiy[xsl,iypf] - 0.125*bbb.mi[0]*( \
                        (bbb.up[bbb.ixm1,self.yy] + bbb.up)[xsl,iypf,0] \
                    )**2*bbb.fniy[xsl,iypf,0] + bbb.cfvisy*0.125\
                    * com.sy[xsl,iypf]*have( \
                        (bbb.visy[:,:,0]*com.gy)[xsl,iypf ],\
                        (bbb.visy[:,:,0]*com.gy)[xsl,iypf+1 ] \
                    ) * ( 
                        (bbb.up[bbb.ixm1,self.yy]+bbb.up)[xsl,iypf+1,0]**2 -
                        (bbb.up[bbb.ixm1,self.yy]+bbb.up)[xsl,iypf,0]**2 
            ) )
            self.ppfe[xsl] -= fluxfacy*bbb.feey[xsl,iypf]
            self.ppfbd[xsl] -= fluxfacy*bbb.fniy[xsl,iypf,0]*bbb.ebind*bbb.ev
#            self.pradhpf += fluxfacy*pwr_pfh[xsl]*com.sy[xsl,iywall]
#            self.pradzpf += fluxfacy*pwr_pfz[xsl]*com.sy[xsl,iywall]
            self.sneutpf[xsl] = cenggw*2.*vynpf*bbb.tg[xsl,0,0]*com.sx[xsl,0]
            self.pneutpf[xsl] += self.sneutpf[xsl]
            if ishymoleng==1:  #molec temp eqn active
                 self.ppfm[xsl] -= bbb.fegy[xsl,iypf,1]
            
            self.__dict__[var] = sum(self.__dict__[var])

        """ CALCULATE TOTAL POWER BALANCE """
        # TODO: is there any way to tidy this up?
        # E.g. calculate divergences for outflows only?
        #
        self.ptotpartin = \
                self.pdiviin + self.pdivein + self.pbindin \
                + self.pbindzin + self.pdivmin  ##self.pneutin
        self.ptotpartout = self.pdiviout + self.pdiveout \
                + self.pbindout + self.pbindzout + self.pdivmout  ##self.pneutout
        self.ptotpart = self.ptotpartin + self.ptotpartout  ##self.pneutout+self.pneutin
        self.ptotrad = self.pradhout + self.pradzout + self.pradhin \
                + self.pradzin
        #
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

        #########################################################################
        # 2-D arrays for checking energy conservation, and for calculation of
        # ionization and radiation sources 

        self.ptotin = self.pbcoree + self.pbcorei \
                + self.p_i_vol + self.p_e_vol


        #########################################################################


        # Here particle power and radiation power are separate terms
        self.powbal = ( \
                self.pbcoree + self.pbcorei + self.pbcorebd \
                + self.p_i_vol + self.p_e_vol + sum(bbb.wjdote[1:-1,1:-1]) \
                - self.ptotpart - sum(self.pwalli) -sum(self.pwalle) -sum(self.pwallm) \
                - sum(self.pwallbd) - sum(self.ppfi) - sum(self.ppfe) - sum(self.ppfm) \
                - sum(self.ppfbd) - sum(self.pradimp) + sum( \
                    - self.pradff - self.pradiz - self.pradrc - self.prdiss \
                    - self.pmomv + self.pibirth \
                ))/ self.pnormpb
        #for cases with no radial pwr input
        if (abs(self.pdivein+self.pdiviin) \
            > 10.*abs(self.ptotin + sum(bbb.wjdote[1:-1,1:-1]))
        ):   
            self.powbal *= abs(self.ptotin + sum(bbb.wjdote[1:-1,1:-1]))\
                /abs(self.pdivein+self.pdiviin)



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

        #two ifs to be able to used old executables
#        if(com.islimon == 1) and (com.nyomitmx != 0):
#           com.nx = com.ix_lim-1
#           ixdivin = com.ix_lim+1
#           iybegin = com.iy_lims

        for var in ['igaswall', 'igaspf', 'igascr']:
            self.__dict__[var] = zeros((com.nx+2,com.ngsp))

        for var in ['iwall', 'ipf', 'icore']:
            self.__dict__[var] = zeros((com.nx+2,com.nfsp))

        for var in['idivout', 'idivin']:
            self.__dict__[var] = zeros((com.ny+2,com.nfsp, com.nxpt))
           
        for var in ['igasout', 'igasin']:
            self.__dict__[var] = zeros((com.ny+2,com.ngsp, com.nxpt))
  
        for var in ['igastot', 'igasdenom', 'isephyd']:
            self.__dict__[var] = 0

        for var in ['i_sat_outer', 'i_sat_inner']:
            self.__dict__[var] = zeros((com.nxpt))

        try:
           fluxfacy = bbb.fluxfacy
        except:
           fluxfacy=1.

        iywall=com.ny       # DEFINITION
        iypf=0
        self.iwall += fluxfacy*bbb.fniy[:,com.ny]
        self.igaswall += fluxfacy*bbb.fngy[:,com.ny]

        for ixp in range(com.nxpt):
            self.ipf[:com.ixpt1[ixp]+1] += \
                fluxfacy*bbb.fniy[:com.ixpt1[ixp]+1,0]
            self.igaspf[:com.ixpt1[ixp]+1] += \
                fluxfacy*bbb.fngy[:com.ixpt1[ixp]+1,0]
            self.ipf[com.ixpt2[ixp]+1:] += \
                fluxfacy*bbb.fniy[com.ixpt2[ixp]+1:,0]
            self.igaspf[com.ixpt2[ixp]+1:] += \
                fluxfacy*bbb.fngy[com.ixpt2[ixp]+1:,0]

            self.igascr[max(com.ixpt1[ixp]+1,0):com.ixpt2[ixp]+1] = \
                fluxfacy*bbb.fngy[max(com.ixpt1[ixp]+1,0):com.ixpt2[ixp]+1,0]
            self.icore[max(0,com.ixpt1[ixp]+1):com.ixpt2[ixp]+1] += \
                fluxfacy*bbb.fniy[max(0,com.ixpt1[ixp]+1):com.ixpt2[ixp]+1,0]
            self.isephyd = bbb.qe*sum(\
                bbb.fniy[max(com.ixpt1[ixp]+1,0):com.ixpt2[ixp]+1,com.iysptrx,0])

            """ ION SATURATION CURRENT """
            for ixp in range(com.nxpt): 
                self.i_sat_outer[ixp] += bbb.qe*bbb.zi[0]\
                                    *sum(bbb.fnix[com.ixrb[ixp],1:com.ny+1,0])
                self.i_sat_inner[ixp] += bbb.qe*bbb.zi[0]\
                                    *sum(bbb.fnix[com.ixlb[ixp],1:com.ny+1,0])
            # ion current to divertor plate
            self.idivout[1:-1,:,ixp] += bbb.fnix[com.ixrb[ixp],1:com.ny+1]
            self.idivin[1:-1,:,ixp] += bbb.fnix[com.ixlb[ixp],1:com.ny+1]
            self.igasout[1:-1,:,ixp] += bbb.fngx[com.ixrb[ixp],1:com.ny+1]
            self.igasin[1:-1,:,ixp] += bbb.fngx[com.ixlb[ixp],1:com.ny+1]
            if bbb.isupgon[1] == 0:
                 self.igasout[1:-1,0,ixp] = bbb.fnix[com.ixrb[ixp],1:com.ny+1,1]
                 self.igasin[1:-1,0,ixp] = bbb.fnix[com.ixlb[ixp],1:com.ny+1,1]

        # ion current to the core

        for var in [
            'idivout', 'idivin', 'igasout', 'igasin', 'icore', 'iwall',
            'igaswall', 'ipf', 'igaspf', 'igascr'
        ]:
            self.__dict__[var] *= bbb.qe

        self.i_vol = bbb.qe*sum(bbb.volpsor[:,:,0])
        if bbb.l_parloss <= 1e9:
          self.i_vol -= sum(bbb.nuvl[:,:,0]*com.vol*bbb.ni[:,:,0])*bbb.qe

        self.igastot += (sum(self.igasin - self.igasout) \
                + sum(-self.igaswall + self.igaspf + self.igascr) + self.i_vol)
        self.igasdenom += sum((self.igasin - self.igasout))

        if bbb.ishymol == 1:  # add id=2 case again because of 2 atoms/molecule
           self.igastot += self.igasin[:,1]-self.igasout[:,1]-self.igaswall[:,1] \
                        +self.igaspf[:,1]+ self.igascr[:,1]
           self.igasdenom += (self.igasin[:,1] - self.igasout[:,1])

        self.deligas = (self.igastot - sum(self.iion) - sum(self.irecomb) \
                -sum(self.icxgas)) / self.igasdenom

        # particle outflow from separatrix




    def print_balance_power(self, redo=True):
        from numpy import sum
        from uedge import bbb, com
        if redo:
            self.balance_power()
        print("\nPower Flow [Watts] from Core to Scrape-off Layer is:")
        print("   Total: {:.2e}".format(self.pbcorei + self.pbcoree + self.pbcorebd))
        print("   Ion contribution: {:.2e}".format(self.pbcorei))
        print("   Electron contribution: {:.2e}".format(self.pbcoree))
        print("   Binding energy contribution: {:.2e}".format(self.pbcorebd))
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
        
        print(self.pradimp.shape)
        if bbb.isimpon != 0:
           print("\nPower [W] lost via impurity radiation is:")
           id2=com.nhsp-1
           for id in range(com.nzdf):
              id2 += com.nzsp[id]
              print("  for nuclear charge = {};  Power = {}".format(bbb.znucl[id2], sum(self.pradimp, axis=(0,1,2))[id]))
              print("  for fixed-fraction species;  Power = {:.2e}".format( sum(self.pradff)))
           
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
        print("   Outerwall_sum: {:.2e}".format(sum(self.pwalli + self.pwalle \
                        + self.pwallm + self.pwallbd + self.pradhwall \
                        + self.pradzwall + self.pneutw)))
        print("      pwalli: {:.2e}".format(sum(self.pwalli)))
        print("      pwalle: {:.2e}".format(sum(self.pwalle)))
        print("      pwallm: {:.2e}".format(sum(self.pwallm)))
        print("      pwallbd: {:.2e}".format(sum(self.pwallbd)))
        print("      pradhwall: {:.2e}".format(sum(self.pradhwall)))
        print("      pradzwall: {:.2e}".format(sum(self.pradzwall)))
        print("      pneutw: {:.2e}".format(sum(self.pneutw)))

        print("\nPower Flow [Watts] incident on Private Flux Wall is:")
        print("   PFwall_sum: {:.2e}".format(sum(self.ppfi + self.ppfe \
                        + self.ppfm + self.ppfbd + self.pneutpf)))
        print("      ppfi: {:.2e}".format(sum(self.ppfi)))
        print("      ppfe: {:.2e}".format(sum(self.ppfe)))
        print("      ppfm: {:.2e}".format(sum(self.ppfm)))
        print("      ppfbd: {:.2e}".format(sum(self.ppfbd)))
        print("      pneutpf: {:.2e}".format(sum(self.pneutpf)))

        print("\nPower Balance: (Pin-Pout)/Pin")
        print("    {:.2e}".format(self.powbal))
        print(50*"=")


    def print_balance_particles(self, redo=True):
        from numpy import sum
        from uedge import bbb, com
        if redo:
            self.balance_particles()

        icore = sum(self.icore, axis=(0))
        print("\nParticle Flow [Amps] from Core to Scrape-off Layer is:")
        for id in range(com.nfsp):
           if bbb.zi[id] > 0:
              print("  for nuclear charge = {}; ion charge = {}".format(bbb.znucl[id],bbb.zi[id]))
              print("  icore = {:.2e}\n".format(icore[id]))

        igascr = sum(self.igascr, axis=(0))
        for id in range(com.ngsp):
           print("  for gas species = {}; igascr = {:.2e}".format(id, igascr[id]))
        print("  hydrogen ion current at separatrix, isephyd = {:.2e}".format(self.isephyd))
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

        idivin = sum(self.idivin, axis=(0,2))
        idivout = sum(self.idivout, axis=(0,2))
        for id in range(com.nfsp):
           if bbb.zi[id] > 0:
              print("  for nuclear charge = {}; ion charge = {}".format(bbb.znucl[id], bbb.zi[id]))
              print("  idivin = {:.2e}   idivout = {:.2e}\n".format(idivin[id], idivout[id]))
         
        print("Neutral Flow [Amps] away from Divertor Plate is:")

        igasin = sum(self.igasin, axis=(0,2))
        igasout = sum(self.igasout, axis=(0,2))
        for id in range(com.ngsp):
           print("  for gas species = {}; igasin = {:.2e},  igasout = {:.2e}".format(id,  igasin[id], igasout[id]))

        print("\nParticle Flow [Amps] incident on Vessel Wall is:")
        for id in range(com.nfsp):
           print("  for ion species = {}; iwall = {:.2e}".format(id, sum(self.iwall[id],axis=(0))))

        for id in range(com.ngsp):
           print("  for gas species = {}; igaswall = {:.2e}".format(id, sum(self.igaswall[id], axis=(0))))


        print("\nParticle Flow [Amps] incident on Private Flux Wall is:")
        for id in range(com.nfsp):
           print("  for ion species = {}; ipf = {:.2e}".format(id, sum(self.ipf[id],axis=(0))))

        for id in range(com.ngsp):
           print("  for gas species = {}; igaspf = {:.2e}".format(id, sum(self.igaspf[id],axis=0)))


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
    from uedge import bbb, com
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from numpy import sum
    import contextlib

    plt.ion()
    

    new = UeBalance()


    '''
    with open('new_balancee.txt', 'w') as f:
        with contextlib.redirect_stdout(f):
            new.balancee_new()

    with open('old_balancee.txt', 'w') as f:
        with contextlib.redirect_stdout(f):
            new.balancee()
    '''
    new.balancee_new()

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
            dif = abs(sum(new.__dict__[var]) - sum(bbb.__getattribute__(var)))/ abs(denom)
#            if len(denom.shape)>1:
#                dif[denom==0] = 0
#            else:
#                dif -= 1e100
#            print(var, dif)
        if dif > 1e-6:
            print(var, dif)
    
    colors = ['k', 'r', 'b', 'c', 'm', 'g']
    for var in ['pwr_plth', 'pwr_pltz', 'pwr_wallz', 'pwr_wallh', 'pwr_pfwallh', 'pwr_pfwallz' 
    ]:
        vardim = bbb.__getattribute__(var).shape
        if len(vardim)>1:
            for ix in range(vardim[-1]):
                denom = bbb.__getattribute__(var)[:,ix]
                denom[denom==0] = 1e-100
                try:
                    dif = abs( (new.__dict__[var][:,ix] - bbb.__getattribute__(var)[:,ix])/ denom)
                except:
                    dif = abs( (new.__dict__[var] - bbb.__getattribute__(var)[:,ix])/ denom)
                if '_pf' in var:
                    for nx in range(com.nxpt):
                        nxsl = slice(com.ixpt1[nx]+1, com.ixpt2[nx]+1)
                        dif[com.ixlb[nx]] = 0
                        dif[com.ixrb[nx]+1]=0
                        dif[nxsl] = 0
                dif[denom == 1e-100] = 0
                if dif.max() > 1e-6:
                    f, ax = plt.subplots()
                    ax.plot(dif, color=colors[ix])
                    ax.set_title(f"{var}: {ix}")
        else:
            denom = bbb.__getattribute__(var)
            denom[denom==0] = 1e-100
            dif = abs( (new.__dict__[var] - bbb.__getattribute__(var))/ denom)
            dif[denom == 1e-100] = 0
            if dif.max() > 1e-6:
                print(var,  dif.max())
                f, ax = plt.subplots()
                ax.plot(dif, color='k')
                ax.set_title(var)
    return new.pwr_pfwallz





