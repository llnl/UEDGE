#Typical mesh (nx=64, ny=32) for DIII-D MHD equilibrium
#Uses inertial neutrals, so six variables (ni,upi,Te,Ti,ng,upg)
#
##package flx;package grd;package bbb

# Initialize pyuedge
from uedge import *

# Set the geometry
bbb.mhdgeo = 1 		#use MHD equilibrium
os.system('rm -f aeqdsk neqdsk')  #change names of MHD eqil. files 
os.system('cp aeqdskd3d aeqdsk')  # (Cannot tab or indent these 3 lines)
os.system('cp neqdskd3d neqdsk')

flx.psi0min1 = 0.98		#normalized flux on core bndry
flx.psi0min2 = 0.98		#normalized flux on pf bndry
flx.psi0sep = 1.00001	        #normalized flux at separatrix
flx.psi0max = 1.07		#normalized flux on outer wall bndry
flx.alfcy = 3.
bbb.ngrid = 1		        #number of meshes (always set to 1)
com.nxleg[0,0] = 16		#pol. mesh pts from inner plate to x-point
com.nxcore[0,0] = 16		#pol. mesh pts from x-point to top on inside
com.nxcore[0,1] = 16		#pol. mesh pts from top to x-point on outside
com.nxleg[0,1] = 16		#pol. mesh pts from x-point to outer plate
com.nysol[0] = 24		#rad. mesh pts in SOL
com.nycore[0] =8		#rad. mesh pts in core

# Finite-difference algorithms (upwind, central diff, etc.)
bbb.methn = 33		        #ion continuty eqn
bbb.methu = 33		        #ion parallel momentum eqn
bbb.methe = 33		        #electron energy eqn
bbb.methi = 33		        #ion energy eqn
bbb.methg = 33		        #neutral gas continuity eqn

# Boundary conditions
bbb.ncore[0] = 4.0e19	        #hydrogen ion density on core
bbb.isnwcono = 3		#=3; use lyne
bbb.lyni = 0.05
##	iflcore = 0	        #flag; =0, fixed Te,i; =1, fixed power on core
bbb.tcoree = 300.	        #core Te 
bbb.tcorei = 300.	        #core Ti
bbb.istewc = 3		        #=3 for gradient-length = lyte
bbb.istiwc = 3		        #=3 for gradient-length = lyti
bbb.lyte = 0.05		
bbb.lyti = 0.05
bbb.recycp[0] = 1.0		#hydrogen recycling coeff at plates
bbb.recycw[0] = 1.0
bbb.nwsor = 1
bbb.matwso[0] = 1

# Transport coefficients
bbb.difni[0] = 0.3		#D for radial hydrogen diffusion
bbb.kye = 0.5		        #chi_e for radial elec energy diffusion
bbb.kyi = 0.5		        #chi_i for radial ion energy diffusion
bbb.travis[0] = 1.		#eta_a for radial ion momentum diffusion

# Flux limits
bbb.flalfe = 0.21		#electron parallel thermal conduct. coeff
bbb.flalfi = 0.21		#ion parallel thermal conduct. coeff
bbb.flalfv = 0.5		#ion parallel viscosity coeff
bbb.flalfgx = 1.		#neut. gas in poloidal direction
bbb.flalfgy = 1.		#neut. gas in radial direction
bbb.lgmax = 0.1
bbb.lgtmax = 0.1
bbb.lgvmax = 0.1

# Solver package
bbb.svrpkg = "nksol"	        #Newton solver using Krylov method
bbb.premeth = "ilut"	        #Solution method for precond. Jacobian matrix

# Neutral parallel momentum eqn
bbb.isupgon[0] = 1
bbb.isngon[0] = 0
com.ngsp = 1
com.nhsp = 2
bbb.ziin[com.nhsp-1] = 0
bbb.cfnidh = 0.                 #coeff. heating from charge-exchange

# Restart from a pfb savefile
bbb.restart = 1		        #Begin from savefile, not estimated profiles
bbb.allocate()		        #allocate space for savevariables
restore('pfd3d_ex_upg.64x32m')	#read in the solution from pfb file

###os.system('ln -s ~/Uedge/uedge/in/aph aph6')
aph.isaphdir = 0                #=0 if atomic data file in run directory
com.istabon = 10

# Execute uedge
bbb.exmain()

# Print out a few variables across outer midplane
print''
print'*** Printing variables versus radial index at outer midplane'
print''
print '************\nradial position relative to separatrix [m]'
print(com.yyc)
print '************\n ion density, ni [m**-3] = '
print(bbb.ni[bbb.ixmp,])
print '************\n parallel ion velocity, up [m/s] = '
print(bbb.up[bbb.ixmp,])
print '************\n electron temp, te [eV] = '
print(bbb.te[bbb.ixmp,]/bbb.ev)
print '************\n ion temp, ti [eV] = '
print(bbb.ti[bbb.ixmp,]/bbb.ev)
print '************\n gas density, ng [m**-3] = '
print(bbb.ng[bbb.ixmp,])
