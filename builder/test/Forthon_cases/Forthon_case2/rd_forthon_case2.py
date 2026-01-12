#Coarse mesh (nx=16, ny=8) for DIII-D MHD equilibrium
#Uses diffusive neutrals, so five variables (ni,upi,Te,Ti,ng)
#
##package flx;package grd;package bbb # Initialize pyuedge
from uedge import *
#from uefacets import *
from uedge.pdb_restore import *
from uedge.hdf5 import *

#bbb.uedge_petscInit()

# Set the geometry
bbb.mhdgeo = 1 		#use MHD equilibrium
os.system('rm -f aeqdsk neqdsk')  #change names of MHD eqil. files 
os.system('cp aeqdskd3d aeqdsk')  # (Cannot tab or indent these 3 lines)
os.system('cp neqdskd3d neqdsk')

flx.psi0min1 = 0.98		#normalized flux on core bndry
flx.psi0min2 = 0.98		#normalized flux on pf bndry
flx.psi0sep = 1.00001	#normalized flux at separatrix
flx.psi0max = 1.07		#normalized flux on outer wall bndry
bbb.ngrid = 1		#number of meshes (always set to 1)
com.nxleg[0,0] = 4		#pol. mesh pts from inner plate to x-point
com.nxcore[0,0] = 4		#pol. mesh pts from x-point to top on inside
com.nxcore[0,1] = 4		#pol. mesh pts from top to x-point on outside
com.nxleg[0,1] = 4		#pol. mesh pts from x-point to outer plate
com.nysol[0] = 6		#rad. mesh pts in SOL
com.nycore[0] =2		#rad. mesh pts in core

# Finite-difference algorithms (upwind, central diff, etc.)
bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 33		#neutral gas continuity eqn

# Boundary conditions
bbb.ncore[0] = 2.5e19	#hydrogen ion density on core
##	iflcore = 0		#flag; =0, fixed Te,i; =1, fixed power on core
bbb.tcoree = 100.		#core Te 
bbb.tcorei = 100.		#core Ti
bbb.tedge = 2.		#fixed wall,pf Te,i if istewcon=1, etc
bbb.recycp[0] = 0.8		#hydrogen recycling coeff at plates

# Transport coefficients
bbb.difni[0] = 1.		#D for radial hydrogen diffusion
bbb.kye = 1.		#chi_e for radial elec energy diffusion
bbb.kyi = 1.		#chi_i for radial ion energy diffusion
bbb.travis[0] = 1.		#eta_a for radial ion momentum diffusion

# Flux limits
bbb.flalfe = 0.21		#electron parallel thermal conduct. coeff
bbb.flalfi = 0.21		#ion parallel thermal conduct. coeff
bbb.flalfv = 1.		#ion parallel viscosity coeff
bbb.flalfgx = 1.e20		#neut. gas in poloidal direction
bbb.flalfgy = 1.e20		#neut. gas in radial direction

# Solver package
bbb.svrpkg = "nksol"	#Newton solver using Krylov method
#bbb.svrpkg = "petsc"	#Newton solver using Krylov method
bbb.premeth = "banded"	#Solution method for precond. Jacobian matrix

# Restart from a pfb savefile
bbb.restart = 1		#Begin from savefile, not estimated profiles
bbb.allocate()		#allocate space for savevariables
bbb.nis.shape
#pdb_restore('./pfd3d_ex.16x8')	#read in the solution from pfb file
hdf5_restore('./h5d3d_ex.16x8')
#hdf5_restore('./testhdf5')
#hdf5_save('testhdf5')
#restore('pfd3d_ex.16x8')	#read in the solution from pfb file

os.system('ln -s ~/Uedge/uedge/in/aph aph6')
com.istabon = 10
aph.isaphdir = 0        #=0 specifies atomic data file is in run directory
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

