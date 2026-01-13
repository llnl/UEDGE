###character*6 probname="itfa40"
#Case of full  toroidal equilibrium for ITER-FEAT with multi-species Carbon
#

####package flx;package grd;package bbb
# Initial pyuedge
from uedge import *
from pdb_restore import *


bbb.mhdgeo=1
bbb.isfixlb=0
bbb.isfixrb=0
os.system('rm -f aeqdsk neqdsk')
os.system('cp aeq_iter-feat aeqdsk')
os.system('cp geq_iter-feat neqdsk')

###character*9 machine="bbb.iter-feat"

# Set the geometry
bbb.ngrid = 1
grd.kxmesh = 1		#=4 for exponential grid in leg regions
grd.dxgas[0] = 1.2e-03
grd.dxgas[1] = 1.2e-03
grd.nxgas[0] = 11
grd.nxgas[1] = 11
grd.alfx[0] = .64
grd.alfx[1] = .64
com.nxleg[0,0]=17
com.nxleg[0,1]=17
com.nxcore[0,0]=14
com.nxcore[0,1]=14
com.nycore[0]=10
com.nysol[0]=16
flx.psi0min1 = 0.95	#core minimum psi
flx.psi0min2 = 0.992	#private flux minimum psi
flx.psi0max = 1.035	#maximum flux at wall
flx.alfcy = 2.0		#nonuniformity factor for radial mesh
grd.slpxt=1.2

# Mesh construction--non orthogonal mesh
com.ismmon=3
grd.istream=0
grd.iplate=1
grd.nsmooth=3
grd.wtmesh1=0.75
grd.dmix0=1.0
execfile('plate.iter-feat.py')
com.isnonog=1	# non orthogonal differencing

# Boundary conditions
bbb.isnicore[0] = 1		#=1uses ncore for density BC
bbb.isngcore[0]=2		# use ionization scale length for gas
bbb.ncore[0] = 6.0e19	        #value of core density if isnicore=1
bbb.curcore = 0.		#core particle current if isnicore=0
bbb.iflcore = 1		        #specify core power 
bbb.pcoree = 5.0e7		#electron power across core
bbb.pcorei = 5.0e7		#ion power across core
bbb.recycp=1.0
bbb.istepfc=3;bbb.istipfc=3	#priv. flux has fixed temperature scale length.
bbb.istewc=3;bbb.istiwc=3	#wall has fixed temperature scale length.
bbb.isnwcono=3;bbb.isnwconi=3	#walls have fixed density scale length
bbb.lyni=0.05;bbb.lyte=0.05;bbb.lyti=0.05
# set walls into 7 zones
bbb.nwsor=7
bbb.xgaso[0:7]=[1.533E-01, 4.599E-01, 3.506E+00, 9.291E+00, 1.508E+01, 1.817E+01, 1.858E+01]
bbb.wgaso[0:7]=[3.066E-01, 3.066E-01, 5.790E+00, 5.790E+00, 5.790E+00, 5.000E-01, 4.089E-01]
bbb.albdso[0:7]=1.
bbb.matwso[0:7]=1.
bbb.xgasi[0:7]=[1.213E-01, 3.638E-01, 6.063E-01, 9.024E-01, 1.246E+00, 1.541E+00, 1.837E+00]
bbb.wgasi[0:7]=[2.425E-01, 2.425E-01, 2.425E-01, 3.705E-01, 2.955E-01, 2.955E-01, 2.955E-01]
bbb.albdsi[0:7]=[0.98,0.98,1.0,1.0,1.0,0.98,0.98]
bbb.matwsi[0:7]=1
bbb.recycw[0]=1.0

bbb.bcee = 5.; bbb.bcei = 3.5	#energy transmission coeffs.
bbb.isupss = 1		        #parallel vel can be supersonic

# Transport coefficients
bbb.difni[0] = 0.3
bbb.kye = 1.; bbb.kyi = 1.
bbb.travis[0]=1.;bbb.parvis[0]=1.

# Flux limits
bbb.flalfe=0.21;bbb.flalfi=0.21;bbb.flalfgx=1.;bbb.flalfgy=1.;bbb.flalfgxy=1.;bbb.flalfv=0.5
bbb.lgmax=0.05
bbb.isplflxl=0

# Finite difference algorithms
bbb.methe=33;bbb.methu=33;bbb.methg=66
bbb.methn=33;bbb.methi=33

# Solver package
bbb.svrpkg="nksol"
bbb.mfnksol = 3
bbb.epscon1 = .005
bbb.ftol = 1.e-8
bbb.iscolnorm = 3	# set to 3 for nksol
bbb.premeth="ilut"
bbb.lfililut = 100
bbb.lenpfac=75
bbb.runtim=1.e-7
bbb.rlx=0.9
###bbb.del=1.e-8

# Neutral gas properties
bbb.tfcx=5.;bbb.tfcy=5.		#Franck-Condon temperatures
bbb.eion = 5.		        #F-C energy to each born ion
bbb.ediss = 10.		        #diss. energy from elec. (ediss=2*eion)
bbb.isrecmon = 1		#e-i recombination (=1 is on)
bbb.ngbackg=1.e12		# minimum floor neutral density
bbb.ingb=4			# parameter used to force floor density

# Inertial neutral model
bbb.isupgon[0]=1;bbb.isngon[0]=0;com.ngsp=1;com.nhsp=2;bbb.ziin[com.nhsp-1]=0
bbb.cngmom=0;bbb.cmwall=0;bbb.cngtgx=0;bbb.cngtgy=0;bbb.cngflox=0;bbb.cngfloy=0;bbb.cfbgt=0
bbb.kxn=0;bbb.kyn=0
bbb.flalftgx=10.0;bbb.flalftgy=10.0

# Currents and potential parameters
bbb.isphion=0
bbb.rsigpl=1.e-8	#anomalous cross-field conductivity
bbb.cfjhf=0.		#turn-on heat flow from current (fqp)
bbb.jhswitch=0		#Joule Heating switch

# Hydrogenic ions
bbb.minu[0:2] = 2.5

# Atomic physics packages
com.istabon=10		#Stotler's rates for istabon=10
aph.isaphdir = 0        #=0 specifies atomic data file is in run directory

## Impurity gas
com.ngsp = 2                #total number of gas species
bbb.isngon[1] = 1           #turns on impurity gas
bbb.ngbackg[1] = 1.e9
bbb.ingb = 2
bbb.istgcon[1] = 1	    #=1 for constant tg(2) at tgas(2)
bbb.tgas[1] = 1.	    #value for tg when istgcon=1
bbb.rcxighg = 0.            # best value; ratio of imp cx to hyd cx
bbb.kelighi[1] = 5.e-16     #elastic sig_v for imp_gas/h_ion
bbb.kelighg[1] = 5.e-16     #elastic sig_v for imp_gas/h_gas
bbb.n0g[1] = 1.e16          #imp. gas density normalization
# Impurity gas boundary conditions
bbb.recycp[1] = 0.01        #plate recycling of impurities
bbb.recycw[1] = 1e-4        #wall recycling; matwsi,o set above for hyd
bbb.isch_sput[1]=7;bbb.isph_sput[1]=3 # Haasz/Davis sputtering model
bbb.t_wall=300;bbb.t_plat=500   #wall and plate temperatures
bbb.crmb=bbb.minu[0]            #mass of bombarding particles
bbb.allocate()		        #allocate chemywi,o etc
bbb.fchemywi=1.;bbb.fchemywo=1.	#scaling factor for chem sputt walls
bbb.fchemylb=1.;bbb.fchemyrb=1.	#scaling factor for chem sputt plates

## Impurity ions
bbb.isimpon = 6             #Use force-balance only
com.nzsp[0] = 6             #number chrg states impurity isotope #1
bbb.csfaclb[2:8] = 2.191    #at plate csout=sqrt(mi_imp/m_dt) for up=up_imp
bbb.csfacrb[2:8] = 2.191    #at plate csout=sqrt(mi_imp/m_dt) for up=up_imp
bbb.minu[2:8] = 12.         #mass in AMU; python doesn't fill last place of [2:8]
bbb.ziin[2:8] = [1,2,3,4,5,6] #charge per ion
bbb.znuclin[0:2] = 1        #nuclear charge for hydrogen species (python fills 2-1)
bbb.znuclin[2:8] = 6        #nuclear charge for impurity (python fills 8-1)
bbb.n0[2:8] = 1.e17         #global density normalization
bbb.nzbackg = 1.e9          #background density for impurities
bbb.inzb = 2		    #exponent for switching on nzbackg
bbb.ismctab = 2             # use Braams' rate tables
com.mcfilename[0] = "C"     # Imp rate file name
com.mcfilename[1] = "_"
com.mcfilename[2] = "r"
com.mcfilename[3] = "a"
com.mcfilename[4] = "t"
com.mcfilename[5] = "e"
com.mcfilename[6] = "s"
com.mcfilename[7] = "."
com.mcfilename[8] = "a"
com.mcfilename[9] = "d"
com.mcfilename[10] = "a"
com.mcfilename[11] = "s"
com.mcfilename = 'C_rates_adas\0'
# Impurity ion boundary conditions
###bbb.isnicore(com.nhsp+com.nzsp[0])= 3  #=3 for flux=curcore & constant ni
###bbb.curcore(com.nhsp+com.nzsp[0]) = 0. #Rcurrent for isnicore=3
bbb.isnicore[com.nhsp+com.nzsp[0]-1]= 3  #=3 for flux=curcore & constant ni
bbb.curcore[com.nhsp+com.nzsp[0]-1] = 0. #Rcurrent for isnicore=3
bbb.isnwcono[2:8] = 3       #use lyni scale-length (set above); outer wall
bbb.isnwconi[2:8] = 3       #use lyni scale-length (set above); inner wall
bbb.nwomin[2:8] = 1e7       #minimum ni at outer sidewall
bbb.nwimin[2:8] = 1e7       #minimum ni at inner sidewall

# Restart from a save file
bbb.restart = 1
bbb.allocate()
pdb_restore ('pfiter_msC.15')

# Execute UEDGE for this case
bbb.exmain()

# Print out a few variables across outer midplane
print''
print'*** Printing variables versus radial index at outer midplane'
print''
print '************\nradial position relative to separatrix [m]'
print(com.yyc)
print '************\n ion densities, ni [m**-3] = '
print(bbb.ni[bbb.ixmp,])
print '************\n parallel ion velocity, up [m/s] = '
print(bbb.up[bbb.ixmp,])
print '************\n electron temp, te [eV] = '
print(bbb.te[bbb.ixmp,]/bbb.ev)
print '************\n ion temp, ti [eV] = '
print(bbb.ti[bbb.ixmp,]/bbb.ev)
print '************\n gas densities, ng [m**-3] = '
print(bbb.ng[bbb.ixmp,])
