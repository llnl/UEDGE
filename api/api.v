api   # Atomic physics (impurities)
{
NTEMP   = 48
NDEN    = 11
NCASET  = 40
NCASENO = 40
NCASENT = 40
NZSORMX = 10   # maximum number of impurity sources
MXMISO=5       # maximum number of charged isotopes; also must set in com.v
MXNZCH=26      # maximum number of charge states for any isotope
MXMINZ=MXMISO*MXNZCH
KXA=3
MXMISO1=MXMISO+1
KNX=3*KXA*KXA*MXNZCH
KMXZ=KXA*MXMINZ
NBA=5		# used in fmombal
}

***** Com_Dim_Vars hidden:    
dim_vars_hidden     integer    # Do not edit this group. It is used to build
                               # the Basis version of the code. 

***** Physical_constants2:
# Add the 2 to distinquish from same constants in package bbb
ev2	real	/1.6022e-19/   # 1 electron volt in Joules
qe2	real	/1.6022e-19/   # elementary charge in Coulombs

***** Normalization_constants:
# Constants set within subroutine inelinput and used in some impurity
# routines to convert to MKS units
crni	real	# multiplicative constant to convert densities to MKS
ctemp	real	# multiplicative constant to convert temperatures to MKS

***** Impfcn:
# routines for computing quantities relevant to impurity radiation
getatau(nx,ny,uu,gx,ixpt1,ixpt2,iysptrx,atau,tau1,tau2)	subroutine
  	    # computes ion/impurity lifetime for pol. flow along flux surfaces,
	    # evaluates arrays atau, tau1, and tau2
getprad(nx,ny,ngsp,te,ne:real,ng:real,afrac,atau,
        prad,na:real,ntau:real,nratio:real)	subroutine
	    # computes the impurity radiation loss for electrons,
 	    # evaluates arrays prad, na, ntau, and nratio

***** Impurity_transport:
# variable specifying radial transport rate of impurities
dnimp		real	/1./	[m**2/s]
methimp         integer /33/   +input # specifies interp. for finite diff. for (y,x)
                               # 66 is log interp., 77 inverse interp.,
                               # otherwise linear interp.
csexpn		real    /0./   # exponent for reducing impurity || vel from
			       # being too supersonic (1+x**(4-csexpn))/(1+x**4)
                               # with x=v||/cs; range should be 0->1, with 0.3 

***** Impurity_source:
# array for impurity source
simpfix(nx,ny)	_real	# vol. source (avg-ion or Z=1) input by user

***** Sources_at_walls:
# arrays for impurity sources on inner and outer walls
nzsor			integer	/0/	# number of impurity sources
iszsorlb(nzspt,nzsor)	_integer /1/	# =1 if origin of ximp* is at left bndry
jxzsori(nzspt,nzsor)    _integer /1/	# ximpi=0. is located at left boundary
                                        # of mesh region jxzsori for iszsorlb=1
jxzsoro(nzspt,nzsor)    _integer /1/ 	# ximpo=0. is located at left boundary
                                        # of mesh region jxzsori for iszsorlb=1
ixzbegi(nzspt,nzsor)    _integer        # ximpi=0. is located at cell index 
                                        # ix=ixzbegi for iszsorlb=1
ixzendi(nzspt,nzsor)    _integer        # ximpi=0. is located at cell index 
                                        # ix=ixzendi for iszsorlb=0
ixzbego(nzspt,nzsor)    _integer        # ximpo=0. is located at cell index 
                                        # ix=ixzbego for iszsorlb=1
ixzendo(nzspt,nzsor)    _integer        # ximpo=0. is located at cell index 
                                        # ix=ixzendo for iszsorlb=0
ximpi(nzspt,nzsor)	_real	[m]	# center of source profile (inner wall)
ximpo(nzspt,nzsor)	_real	[m]	# center of source profile (outer wall)
wimpi(nzspt,nzsor)	_real	[m]	# width of source profile (inner wall)
wimpo(nzspt,nzsor)	_real	[m]	# width of source profile (outer wall)
impsori(nzspt,nzsor)	_real	[Amp]	# impurity source strength (inner wall)
impsoro(nzspt,nzsor)	_real	[Amp]	# impurity source strength (outer wall)

***** Input:
# variables for input file names
inelrates(1)		Filename	/'edgrat.dat'/  # (1) needed for pymac
inelrad(1)		Filename	/'carbavg.dat'/
inelmc(1)		Filename	/'carbmc.dat'/

***** Radiation:
# dimension variables and arrays for radiation
# read by namelist so cannot be allocated dynamically
ncaset		integer	/40/
ncaseno		integer	/40/
ncasent		integer	/40/
terad(NCASET)		real
xno(NCASENO)		real
rntau(NCASENT)		real
radrate(NCASET,NCASENO,NCASENT)	real
avgz(NCASET,NCASENO,NCASENT)	real
avgz2(NCASET,NCASENO,NCASENT)	real

***** Impdata:
# name of directory containing impurity atomic rate data
apidir character*120 /"."/  # name of directory containing impurity data files
coronalimpfname  character*120 /"mist.dat"/   +input # name coronal impurity rate file
                                                     # File looked for relative to coronalimpfname
coronalimppath  character*120 /"."/   +input    # path to folders of coronal imp files
                                                # file name defined by coronalimpfname
  

***** MC_subs:
readmc(nzdf:integer,mcfilename:string)		subroutine
	# reads formatted data tables for one or more impurities from files
	# specified by mcfilename(1:nzdf) (default filename is 'b2frates')
mcrates(ne:real,te,ti,\
        za:integer,zamax:integer,zn:integer,rion,rrec,rcxr)  subroutine
	# computes rate parameters for impurity transitions that originate
	# from atomic charge state za of the impurity with nuclear charge zn,
	# where the maximum atomic charge state is zamax.
radmc(nz:integer, znuc:integer, te, dene, denz, radz)	function
	# computes the radiation rates radz(0:nz) in [W/m**3] for various
	# charge states of an impurity with nuclear charge znuc;
	# this function returns the total electron energy loss rate in [W/m**3]
rcxr_zn6(tmp, za:integer)	function
	# compute cx rate for carbon ions on neutral hydrogen using
	# polynomial fit to curves in thesis by C.F. Maggi (1997).
	# input neutral hydrogen temperature, tmp, is in Joules/AMU;
	# za is initial carbon ion charge;
	# output is <sigma*v> in m**3/sec.

***** Impurity_charge:
nnz     integer  /2/
zq(nnz)	_real

***** P93dat:
# variables associated with impurity radiation from Post '93 tables
atn	integer			# atomic number
atw	integer			# atomic weight
nt	integer 		# number of temperature values
nr	integer			# number of density ratio values
nn	integer			# number of n*tau values
tdatm(nt,nr,nn)		_real	[J]		# temperature
rdatm(nt,nr,nn)		_real			# density ratio
ndatm(nt,nr,nn)		_real	[sec/m**3]	# n*tau
emdatm(nt,nr,nn)	_real	[Watts-m**3]	# emissivity
z1datm(nt,nr,nn)	_real			# average Z
z2datm(nt,nr,nn)	_real			# average Z**2

***** Imslwrk:
# working arrays for 3-d spline interpolation
nxdata_api	integer
nydata_api	integer
nzdata	integer
xdata_api(1:nxdata_api)	_real
ydata_api(1:nydata_api)	_real
zdata(1:nzdata)	_real
fdata_api(1:nxdata_api,1:nydata_api,1:nzdata)	_real
ldf_api	integer			# first dimension of 3-d data array
mdf	integer			# second dimension of 3-d data array
iflagi	integer			# input/output flag for 3-d spline routines
nwork2	integer			# size of array work2
work2(nwork2)	_real		# work array for B3VAL
nwork3	integer			# size of array work3
work3(nwork3)	_real		# work array for B3INT
iworki(10)	integer		# work array for B3VAL
icont		integer	/0/	# input flag for B3VAL
kxords_api	integer	/4/	# order of spline fit versus x
#	kxords_api=4 (default) is cubic interpolation
kyords_api	integer	/4/	# order of spline fit versus y
#	kyords_api=4 (default) is cubic interpolation
kzords	integer	/4/	# order of spline fit versus z
#	kzords=4 (default) is cubic interpolation
xknots_api(1:nxdata_api+kxords_api)	_real
yknots_api(1:nydata_api+kyords_api)	_real
zknots(1:nzdata+kzords)	_real
emcoef(1:nxdata_api,1:nydata_api,1:nzdata)	_real	# spline coeff's for emissivity 
z1coef(1:nxdata_api,1:nydata_api,1:nzdata)	_real	# spline coeff's for average-Z
z2coef(1:nxdata_api,1:nydata_api,1:nzdata)	_real	# spline coeff's for average-Z**2

***** P93fcn:
# setup and evaluation routines for spline representation of data on impurity
# radiation and charge state in POST '93 tables
readpost(fname:string)	subroutine
	# read formatted data table for one impurity
splinem		subroutine
	# construct 3-d spline representation for impurity radiation and charge
emissbs(te:real,nratio:real,ntau:real)	real function
	# radiated power per impurity atom per electron [Watts-m**3]
	# te = e-temperature		[J]
	# nratio = ng / ne		[none]
	# ntau = ne * tau-impurity	[sec/m**3]
z1avgbs(te:real,nratio:real,ntau:real)	real function
	# average Z of impurity
	# te = e-temperature		[J]
	# nratio = ng / ne		[none]
	# ntau = ne * tau-impurity	[sec/m**3]
z2avgbs(te:real,nratio:real,ntau:real)	real function
	# average Z**2 of impurity
	# te = e-temperature		[J]
	# nratio = ng / ne		[none]
	# ntau = ne * tau-impurity	[sec/m**3]

***** Reduced_ion_constants:
# Constants used within reduced-ion routines
 coulom		real	[coulomb]	# Elementary charge
 epsilo		real	[farad/m]	# Permittivity of free space
 promas		real	[kg]	# Proton mass
 xj7kv		real	[J/keV]	# Joules/keV units conversion
 one		real	[none]	# Unity
 pi0		real	[none]	# Pi
 zero		real	[none]	# Zero
 sumforce	real	# 
 totmass	real	# 
 anorm		real	# 
 acci		real	# 
 acci0		real	# 
 al32(3)	real	# 
 miso		integer	# 
 nzch		integer	# 
 mise		integer	# "isotope number" for electrons (=1)
 ilam1		integer	# =1
 ilam2		integer	# =2
 ilam3		integer	# =3
 iacci		integer	# =4
 iforc		integer	# =5
 natom(MXMISO)	integer	# 

***** Reduced_ion_variables:
# Variables used within reduced-ion routines
 capm(KXA*MXMISO*KXA*MXMISO)	real	# 
 capn(KXA*MXMISO*KXA*MXMISO)	real	# 
 caplam(KXA*MXMISO)		real	# 
 fmomenta(KMXZ)			real	# 
 denz(MXMISO*MXNZCH)		real	# 
 denmass(MXMISO*(MXNZCH+1))	real	# 
 ela(KXA*KXA*MXMISO)		real	# 
 elab(KXA*MXMISO*KXA*MXMISO)	real	# 
 mntau(MXMISO*MXMISO)		real	# 
 usol(KXA*MXNZCH*MXMISO)	real	# 
 sbar(KXA*MXMISO1)		real	# 
 zi_api(MXMISO*MXNZCH)		real	# 

***** Cyield:
# Variables used for DIVIMP physical sputtering models
ceth(7,12)	real [eV]	# threshold energy for physical sputtering
cetf(7,12)	real [eV]	# Thomas-Fermi potential
cq(7,12)	real 		# yield factor for phys. sputt. (atoms/ion)
ntars		integer		# number of target materials
cidata(7,12)	logical		# flag indicating avaiability of data
redf_haas	real /0.2/	# low-energy reduction of Haas97 chem sput 
			        # for IOPTCHEM=7

***** Sputt_subs:
# Subroutines and functions used for sputtering yields
syld96(matt:integer,matp:integer,cion:integer,\
       cizb:integer,crmb:real)				subroutine
yld96(matt:integer,matp:integer,energy:real)		real function
sputchem(ioptchem:integer,ee0:real,temp:real,\
         flux:real,ychem:real)				subroutine

***** Emissivities:
ntemp		integer	# number of temperature values
nlam		integer	# number of wavelengths
nden		integer	# number of density values
etemp(ntemp)	  	_real	[eV]		# electron temperature
lamb(nlam)	  	_real	[Angstrom]	# wavelength of line
eden(nden)		_real	[m^-3]		# electron density
rate(nlam,ntemp,nden)	_real	[ph/s]		# emission rate coeff
emiss(nlam,0:nx+1,0:ny+1)	_real	[ph/m^3/s]	# emission rate
readrates(apidir:string,impfname:string)        subroutine
#read Isler data files
calcrates(ne:real,te:real,density:real)         subroutine
#calculate emiss(1:nlam,0:nx+1,0:ny+1) via (ne,te) interpolation on Isler data

***** Pixels:
nrpix   integer   # number of pixels in horizontal direction
nzpix   integer   # number of pixels in vertical direction
npd(nrpix,nzpix) _integer /0/ # number of path data points per pixel
rp1(nrpix,nzpix) _real /0./   # R of 1st path-end-point for this pixel
zp1(nrpix,nzpix) _real /0./   # Z of 1st path-end-point for this pixel
rp2(nrpix,nzpix) _real /0./   # R of 2nd path-end-point for this pixel
zp2(nrpix,nzpix) _real /0./   # Z of 2nd path-end-point for this pixel
wt(nrpix,nzpix)  _real /0./   # path length through this pixel
rminpix   real   [m]   # left boundary (R) of pixel domain
rmaxpix   real   [m]   # right boundary (R) of pixel domain
zminpix   real   [m]   # lower boundary (Z) of pixel domain
zmaxpix   real   [m]   # upper boundary (Z) of pixel domain
drpix     real   [m]   # pixel size in horizontal direction
dzpix     real   [m]   # pixel size in vertical direction
lineintegral(arg:real,rvertex:real,zvertex:real)  real function
# computes the line integral of arg(nrpix,nzpix) along the path
# (rvertex(1:2),zvertex(1:2)) where arg(ii,jj) is the pixel
# representation of any UEDGE array arg_ue(ix,iy), obtained by
# mapping via DCE subroutine rzxform, e.g., arg =
# rzxform(arg_ue,rm,zm,nrpix,nzpix,rminpix,rmaxpix,zminpix,zmaxpix)
