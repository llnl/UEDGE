com   # Data and routines used in multiple packages
{
MXMISO=5        # maximum number of charged isotopes; also must set in api.v
ngspmx = 6	# maximum number of gas species; also must set in bbb.v
}

***** OMFIT:
iomfit  integer /1/     # switch used by OMFIT


***** COMroutines:
glbwrlog(ioun)          function
   # writes cpu,io,sys,mem to i/o unit number ioun

***** Dim:
# Variables that contain widely-used dimensions
nx	integer	/8/ +griddata	# number of cells in x (poloidal) direc. (see also nxm)
ny	integer	/4/ +griddata	# number of cells in y (radial) direction
nxm	integer	/8/	# total number of cells in x direc.; nxm = nx+|nxomit|
nym     integer /4/     # total number of cells in y direc.; nym = ny+nyomitmx
nxpt	integer /1/	# number of x-points in (R,Z) simulation domain
nhsp	integer	/1/	+input # number of hydrogenic species
nzsp(1:ngspmx-1) integer /5*0/ +input # number of impurity species per gas species
nzspt   integer /1/     # total number of impurity species
nzspmx  integer	   /10/	# maximum of nzsp(igsp) used for storage allocation
nisp	integer	/1/	+input # number of ion species
nusp	integer		+input # number of parallel momentum equations
nfsp	integer		+input # number of cont. eqns or flux species (calc internal)
ngsp	integer	/1/ +input +regrid	# number of gas species
nhgsp   integer /1/     +input # number of hydrogen gas species (prepare for tritium)
imx     integer /50/    # size in x of Zagorski arrays
imy     integer /40/    # size in y of Zagorski arrays
lnst	integer /43/    # size of impurity species(?) for Zagorski
num_elem integer /0/    # size of expt profile vs psi arrays

***** Dimflxgrd:
# dimensions used in both flx and grd packages
jdim   integer	  # (calculated) number of flux surfaces in horiz direction
npts   integer	  # (calculated) maximum number of points on a flux surface
noregs integer	  /2/	    # number of mesh regions in horizontal direction
nxefit integer	  # number of EFIT grid surfaces in x (horizontal) direction
nyefit integer	  # number of EFIT grid surfaces in y (vertical) direction
nlim   integer	  # number of EFIT points on limiter/vessel boundary
nbdry  integer    # number of EFIT points on last closed flux surface
nwork  integer	  # length of 2-d spline workspace array work

***** Comflxgrd:
# variables used in both flx and grd packages
isfw	    integer /0/
     # transfer data via flx-grd file if isfw=1
runid        character*60
     # output label from EFIT in neqdsk file
xold(nxefit)    _real    [m]
     # horizontal position (R) of grid surfaces in EFIT
yold(nyefit)    _real    [m]
     # vertical position (Z) of grid surfaces in EFIT
fold(nxefit,nyefit)     _real    [volt-sec/radian]
     # poloidal magnetic flux at EFIT grid points
bcentr  real    [Tesla]
     # reference toroidal field
rcentr  real    [m]
     # major radius where vacuum toroidal field is bcentr
rmagx   real    [m] +griddata
     # major radius of magnetic axis
zmagx   real    [m] +griddata
     # vertical position of magnetic axis
simagx  real    [volt-sec/radian] +griddata
     # poloidal magnetic flux at the magnetic axis
sibdry  real    [volt-sec/radian] +griddata
     # poloidal magnetic flux at primary x-point
sibdry1  real    [volt-sec/radian] +griddata
     # poloidal magnetic flux at the lower x-point
sibdry2  real    [volt-sec/radian] +griddata
     # poloidal magnetic flux at the upper x-point
xdim    real    [m]
     # horizontal (radial) extent of EFIT computational grid
zdim    real    [m]
     # vertical extent of EFIT computational grid
zmid	real	[m] +griddata
     # vertical height of midplane above bottom of EFIT mesh
zshift	real	[m] +griddata
     # vertical shift of EFIT data to put z=0 at the bottom
workk(nxefit)   _real
     # dummy array for neqdsk output
fpol(nxefit)    _real	[m-T]
     # poloidal current F(psi) = R*B_toroidal from geqdsk
pres(nxefit)    _real
     # dummy array for neqdsk output
qpsi(nxefit)    _real
     # dummy array for neqdsk output
rgrid1          real    [m]
     # major radius at inner edge of EFIT grid
cpasma          real    [Amps]
     # toroidal current in the core plasma
rbdry(nbdry)    _real   [m]
     # radial coordinates of last closed flux surface from EFIT
zbdry(nbdry)    _real   [m]
     # vertical coordinates of last closed flux surface from EFIT
xlim(nlim)      _real   [m] +griddata
     # radial coordinates of limiter/vessel boundary from EFIT
ylim(nlim)      _real   [m] +griddata
     # vertical coordinates of limiter/vessel boundary from EFIT
bscoef(nxefit,nyefit)    _real
     # 2-d spline coefficients 
kxord	integer	/4/
     # order of the 2-d spline in the x-direction
kyord	integer	/4/
     # order of the 2-d spline in the y-direction
xknot(nxefit+kxord)		_real
     # knot sequence of the 2-d spline in the x-direction
yknot(nyefit+kyord)		_real
     # knot sequence of the 2-d spline in the y-direction
ldf		integer
     # leading dimension of 2-d spline coefficient array
iflag		integer
     # input/output flag for 2-d spline subroutines
work(nwork)	_real
     # work array for 2-d spline subroutines
jmin(2)		integer
     # index of first flux surface for inboard and outboard half of mesh
jmax(2)		integer
     # index of last flux surface for inboard and outboard half of mesh
jsptrx(2)	integer
     # separatrix flux surface index for inboard and outboard half of mesh
jaxis           integer
     # magnetic axis flux surface index that separates two halves of the mesh
xlbnd   real    [m]
     # radial cutoff for contours at inner edge of efit domain
xubnd   real    [m]
     # radial cutoff for contours at outer edge of efit domain
ylbnd   real    [m]
     # vertical cutoff for contours at lower edge of efit domain
yubnd   real    [m]
     # vertical cutoff for contours at upper edge of efit domain
xcurve(npts,jdim)       _real    [m]
     # radial position of nth data point on jth contour segment
ycurve(npts,jdim)       _real    [m]
     # vertical position of nth data point on jth contour segment
npoint(jdim)    _integer
     # number of data points on jth contour segment
geqdskfname  character*128 /'neqdsk'/ 
     #File name for geqdsk/neqdsk EFIT file
  
***** Aeqflxgrd:
     # information read from the A-file produced by the EFIT code
vmonth	integer		# EFIT version month
vday	integer		# EFIT version day
vyear	integer		# EFIT version year
eshot   integer
     # shot number from EFIT
etime   real    [msec]
     # time slice from EFIT
rseps   real    [m] +griddata
     # major radius of lower x-point
zseps   real    [m] +griddata
     # vertical position of lower x-point
rseps1   real    [m] +griddata
     # major radius of lower x-point
zseps1  real    [m] +griddata
     # vertical position of lower x-point
rseps2   real    [m] +griddata
     # major radius of upper x-point
zseps2  real    [m] +griddata
     # vertical position of upper x-point
rvsin   real    [m] +griddata
     # major radius of inboard strike point
zvsin   real    [m] +griddata
     # vertical position of inboard strike point
rvsout  real    [m] +griddata
     # major radius of outboard strike point
zvsout  real    [m] +griddata
     # vertical position of outboard strike point
mco2v	integer	/3/
     # number of vertical co2 chords
mco2r	integer	/1/
     # number of radial co2 chords
rco2v(mco2v)	_real	
     # radial positions of vertical co2 chords
dco2v(mco2v)	_real
     # densities for vertical co2 chords
rco2r(mco2r)	_real	
     # vertical positions of radial co2 chords
dco2r(mco2r)	_real
     # densities for radial co2 chords
nsilop	integer	/41/
     # number of flux loops
csilop(nsilop)	_real
     # data from flux loops
magpri	integer	/60/
     # number of magnetic pickup loops
cmpr2(magpri)	_real
     # data from magnetic pickup loops
nfcoil	integer	/18/
     # number of f-coils
ccbrsp(nfcoil)	_real
     # data from f-coils
nesum	integer	/2/
     # number of e-coils
eccurt(nesum)	_real
     # data from e-coils
aeqdskfname  character*128 /'aeqdsk'/ 
     #File name for aeqdsk EFIT file

***** RZ_grid_info:
# RZ grid location data and magnetic-field info at those locations
rm(0:nxm+1,0:nym+1,0:4)     _real [m]    +griddata # radial cell position, 0 is center
zm(0:nxm+1,0:nym+1,0:4)     _real [m]    +griddata # vertical cell position, 0 is center
rmt(0:nxm+1,0:nym+1,0:4)    _real [m]    +griddata # temp rad cell position, 0 is center
zmt(0:nxm+1,0:nym+1,0:4)    _real [m]    +griddata # temp vert cell position, 0 is center
rv(0:nxm+2,-1:nym+1)        _real [m]    # rad position velocity cell corner
zv(0:nxm+2,-1:nym+1)        _real [m]    # vert position velocity cell corner
psi(0:nxm+1,0:nym+1,0:4)    _real [Tm^2] +griddata # polodial magnetic flux
br(0:nxm+1,0:nym+1,0:4)     _real [T]    +griddata # radial magnetic field, 0 is center
bz(0:nxm+1,0:nym+1,0:4)     _real [T]    +griddata # vert magnetic field, 0 is center
bpol(0:nxm+1,0:nym+1,0:4)   _real [T]    +griddata # pol. magnetic field, 0 is center
bphi(0:nxm+1,0:nym+1,0:4)   _real [T]    +griddata # tor.l magnetic field, 0 is center
b(0:nxm+1,0:nym+1,0:4)      _real [T]    +griddata # total magnetic field, 0 is center
bsqr(0:nxm+1,0:nym+1)       _real [T]    # B**2 at density-cell center
b12(0:nxm+1,0:nym+1)        _real [T]    # B**0.5 at vel-cell center/dens face
b12ctr(0:nxm+1,0:nym+1)     _real [T]    # B**0.5 at density-cell center
b32(0:nxm+1,0:nym+1)        _real [T]    # B**1.5 at velocity-cell center

***** RZ_grid_global:
# Global RZ grid location data and magnetic-field info for domain decomposition
rmg(0:nxm+1,0:nym+1,0:4)    _real [m]    # global rad cell pos., 0 is center
zmg(0:nxm+1,0:nym+1,0:4)    _real [m]    # global vert. cell pos., 0 is center
psig(0:nxm+1,0:nym+1,0:4)   _real [Tm^2] # global polodial magnetic flux
brg(0:nxm+1,0:nym+1,0:4)    _real [T]    # global rad B field, 0 is center
bzg(0:nxm+1,0:nym+1,0:4)    _real [T]    # global vert. B field, 0 is center
bpolg(0:nxm+1,0:nym+1,0:4)  _real [T]    # global pol. B field, 0 is center
bphig(0:nxm+1,0:nym+1,0:4)  _real [T]    # global tor. B field, 0 is center
bg(0:nxm+1,0:nym+1,0:4)     _real [T]    # global total B field, 0 is center

***** Share:
# Variables used by more than one package
nycore(30)	integer	/30*0/	+regrid +gridgen
   # number of radial zones in the core region of the edge plasma
nysol(30)	integer	/30*2/	+regrid +gridgen
   # number of radial zones in the open flux surface region of the edge plasma
nyout(30)	integer	/30*0/	+regrid +gridgen
   # number of radial zones beyond second separatrix of full double-null
nxleg(30,2)	integer	/0, 59*2/	+regrid +gridgen
   # number of cells along divertor legs:  (,1) inside, (,2) outside
nxcore(30,2)	integer	/0, 59*4/	+regrid +gridgen
   # number of cells along core boundary:  (,1) inside, (,2) outside
nxomit		integer /0/ +gridgen
   # number of x-cells to omit from ix=0; if <0, omit from ix=nx for fluid eqns
   # done typically for isfixlb=1,2 cases with reflection b.c. at left boundary
nxxpt +gridgen
   # number of extra poloidal cells at x-point (per quadrant)
nyomitmx	integer /0/ +gridgen
   # number of y-cells to omit from iy=ny; used to do core only for fluid eqns
igrid		integer	/1/	# loop index for which grid; e.g., nysol(igrid)
isgriduehdf5   integer /0/  
geometry	character*16 /"snull"/ +input
				# specifies magnetic configuration, e.g.,
				# ='snull' for lower single null
				# ='uppersn' for upper single null
                                # ='dnbot' for bottom half of double null
				# ='dnull' for full double null
				# ='snowflake15' for Ryutov's theta~15 deg
				# ='snowflake45' for Ryutov's theta~45 deg
				# ='snowflake75' for Ryutov's theta~75 deg
				# ='snowflake105' for Ryutov's theta~105 deg
				# ='snowflake135' for Ryutov's theta~135 deg
				# ='snowflake165' for Ryutov's theta~165 deg
				# ='dnXtarget' for dnbot with Xtarget
nxc		integer	/4/	# center index of x-grid, normally nx/2,
                                # OR
                                # guard cell index for 'dnbot' symmetry b.c.
simagxs         real            +griddata # "shared" value of simagx from flx package
sibdrys         real            +griddata # "shared" value of sibdry from flx package
nxpt2msh        integer /0/     # number of 1D cells above ixpt2 reset/no-cross
nxpt2psh        integer /0/     # number of 1D cells below ixpt2 reset/no-cross
zxpt2msh        real /0.01/ [m] # spacing iy=1 1D vertices reset above ixpt2
rxpt2msh        real /0.0/  [m] # spacing iy=1 1D vertices reset above ixpt2
zxpt2psh        real /0.01/ [m] # spacing iy=1 1D vertices reset below ixpt2
rxpt2psh        real /0.0/  [m] # spacing iy=1 1D vertices reset below ixpt2
nxpt1msh        integer /0/     # number of 1D cells above ixpt1 reset/no-cross
nxpt1psh        integer /0/     # number of 1D cells below ixpt1 reset/no-cross
zxpt1msh        real /0.01/ [m] # spacing iy=1 1D vertices reset above ixpt1
rxpt1msh        real /0.0/  [m] # spacing iy=1 1D vertices reset above ixpt1
zxpt1psh        real /0.01/ [m] # spacing iy=1 1D vertices reset below ixpt1
rxpt1psh        real /0.0/  [m] # spacing iy=1 1D vertices reset below ixpt1
ismpsym         integer /0/     +gridgen # =1 re-constructs "guard" cells at midplane
                                # of "dnbot" via up/down symmetry
isudsym         integer /0/     +input #=1 up-down symmetric setup (only down part is modeled)
islimon		integer	/0/ +input
	# option switch to apply b.c.'s at ix_lim and ix_lim+1
	# =1 turns on limiter/continuity b.c.'s
ix_lim          integer
	# poloidal index of 1st limiter guard cell
iy_lims         integer	/9999/
	# radial index of 1st flux surface where limiter b.c.'s apply
theta_split	real	/1.570796326794896/
	# (computed) poloidal angle where in/outboard regions of the mesh meet
isnonog		integer /0/	+input # is 9-point (non-orthog) differencing on
ismmon          integer /0/ +gridgen
	# flag that controls mesh modification
	#      = 0    strictly orthogonal mesh and divertor plates
	#      = 1    non-orthogonal mesh, compressed distrib'n on each surface
	#      = 2    non-orthogonal mesh, standard distrib'n on all surfaces
	#      = 3    combination of options 1 and 2 using weight factor
        #             wtmesh1
isoldgrid       integer /0/     # =1 returns to pre-1/25/96 nonorthog grid
isgrdsym        integer /0/	+gridgen # =1 makes symmetric grid for mhdgeo.ne.1,
                                # requires nxm,ixpt2 even
cutlo           real  /1.e-300/ #rough range of a small number(avoids R1MACH(2))
                                #CAUTION; this could be machine dependent
                                #1.e-300 is generally conservative for 64-bit
                                #words; on SUN double-p (e-330); C-90 (e-2500)
epslon		real /1e-6/     #scale factor for relative guard-cell size
spheromak	integer /0/	+input # =0 for tokamak, use vacuum toroidal field
				# =1 for spheromak, use bphi=fpoloidal(psi)/rm
isfrc	integer	/0/	+gridgen # =1 for field-reversed-configuration mesh options
ishalfm	integer	/0/	+gridgen # =1 for doing outboard half-mesh only
isbphicon integer /0/   +gridgen # =0 gives bphi=bcentrg*rcentrg/rm; =1 gives bphi=bcentrg
nhdf      integer /1/   # number of hydrogenic b2frates-format data files
hdfilename(1:12)   Filename  /12*'b2frates'/  # names of hydrogenic data files
nzdf	  integer /1/	# number of impurity   b2frates-format data files
mcfilename(1:12)   Filename  /12*'b2frates'/  +input # names of impurity   data files
istabon	integer	/7/ +input
                       #turns on look-up table for hydrogenic rate coefficients
#	= 0  simple analytic rates and constant energy loss per ionization
#	= 1  table look-up from ADPAK; rates.adpak
#	= 2  table look-up from STRAHL; rates.strahl
#	= 3  table look-up old DEGAS (created 84/09/14); eh.dat & atmc.dat
#	= 4,5,6  table look-up from new DEGAS (created 93/05/06): nwfits
#	=   4  linear interpolation for  rsa  vs log(te) and log10(ne)
#	=   5  spline fit for log10(rsa) vs log(te) and log10(ne)
#  temp.   = 6  Hindmarsh spline    log10(rsa) vs log(te) and log10(ne)
# disabled = 6  spline fit for       rsa  vs log(te) and log10(ne)
#	= 7  Campbell's poly. fit for rsa, etc., vs log10(te) and log10(ne)
#	= 8 tab look-up from latest DEGAS (created 93/04/08); ehr1.dat
#	= 9 tab look-up; Stotler PPPL (~95/07/10) with log(Te)-sigv; ehr2.dat
#	= 10 tab look-up;Stotler PPPL (~95/07/10);log(Te)-log(sigv); ehr2.dat
#	= 11 same as istabon=10 with data for n=5-9 excited states; thin.dat
#	= 12 as istabon=11, ex. Lyman-alpha local absb (thick); thickLyA.dat
#	= 13 as istabon=11, ex. all Lyman lines loc absorbed; thickAllLy.dat
#	= 14 H.Scott data; add rtau=Ly-a opacity, lin. interp; ehrtau.dat
#	= 15 H.Scott data; add rtau=Ly-a opacity, log interp; ehrtau.dat
#	= 16 table look-up using 'b2frates_hyd' file data, log-log interp
#	= 17 tab look-up; latest DEGAS2/PPPL rates: ehr5.dat (23 Apr 2015)

reset_core_og			integer /0/
        # flag to force diff stencil orthog in core; use if isupstreamx=1

***** Xpoint_indices:
# target plate, x-point and separatrix indices
ixlb(1:nxpt)	_integer +griddata	# ix index of left boundary
ixpt1(1:nxpt)	_integer +griddata	# ix index of first x-point
ixmdp(1:nxpt)	_integer +griddata	# ix index of midplane
ixpt2(1:nxpt)	_integer +griddata	# ix index of second x-point
ixrb(1:nxpt)	_integer +griddata	# ix index of right boundary
iysptrx1(1:nxpt)	_integer +griddata	# iy index of first x-point
iysptrx2(1:nxpt)	_integer +griddata	# iy index of second x-point
iysptrx			integer	 +griddata	# iy index of last closed flux surface

***** Cut_indices:
# Indices for cut location (in x) and length (in y)
ixcut1	integer		# ix index of 1st  cut
ixcut2	integer		# ix index of 2nd  cut
ixcut3	integer		# ix index of 3rd  cut
ixcut4	integer		# ix index of 4th  cut
iycut1	integer		# iy length of 1st  cut
iycut2	integer		# iy length of 2nd  cut
iycut3	integer		# iy length of 3rd  cut
iycut4	integer		# iy length of 4th  cut

***** Comgeo:
#Variables in common -- geometry
vol(0:nx+1,0:ny+1)  _real   [m^3]   #Volume of the (ix,iy) primary cell
gx(0:nx+1,0:ny+1)   _real   [m^-1]  #1/(x-diameter) of the (ix,iy) primary cell
                                     #and x-distance between velocity points
gy(0:nx+1,0:ny+1)   _real   [m^-1]  #1/(y-diameter) of the (ix,iy) primary cell
dx(0:nx+1,0:ny+1)   _real   [m]     # x-distance of (ix,iy) primary cell;=1/gx
dxvf(0:nx+1,0:ny+1) _real   [m]     # x-distance of between vel cell ctrs
dy(0:nx+1,0:ny+1)   _real   [m]     # y-distance of (ix,iy) primary cell;=1/gy
dynog(0:nx+1,0:ny+1) _real  [m]     #sep btwn nonog interp pts normal to y-face
gxf(0:nx+1,0:ny+1)  _real   [m^-1]  #1/(x-distance) between density centers
##gxfn(0:nx+1,0:ny+1) _real   [m^-1]  #1/dist. btwn interp. pts, normal to x-face
dxnog(0:nx+1,0:ny+1) _real  [m]     #sep btwn nonog interp pts normal to x-face
gyf(0:nx+1,0:ny+1)  _real   [m^-1]  #1/(y-distance) perp to y-face between
                                    #niy0 and niy1 variables, etc.
gxc(0:nx+1,0:ny+1)  _real   [m^-1]  #1/(x-distance) between density vertices
gyc(0:nx+1,0:ny+1)  _real   [m^-1]  #1/(y-distance) between density vertices
xnrm(0:nx+1,0:ny+1) _real   [ ]     #x-grid, norm. in three separate regions
xvnrm(0:nx+1,0:ny+1) _real  [ ]     #x-vel-grd, norm. in three separate regions
ynrm(0:nx+1,0:ny+1) _real   [ ]     #y-grid, norm. in two separate regions
yvnrm(0:nx+1,0:ny+1) _real  [ ]     #y-vel-grid, norm. in two separate regions
sx(0:nx+1,0:ny+1)   _real   [m^2]   #area normal to 'east' flux-surface of
		  	            #(ix,iy) primary cell;vtag=0 for orthog.
sxnp(0:nx+1,0:ny+1)  _real   [m^2]  #area of tilted 'east' pol face
sy(0:nx+1,0:ny+1)   _real   [m^2]   #area of the 'north' surface of the
                                     #(ix,iy) primary cell
rr(0:nx+1,0:ny+1)   _real           #pitch of the field line at the center of
				     #the (ix,iy) primary cell (Bpol/B)
xcs(0:nx+1)         _real   [m]     #x-coordinate of the center of the (ix,iy)
				     #primary cell just outside separatrix
xfs(0:nx+1)         _real   [m]     #coordinate of the x-face of the (ix,iy)
 				     #primary cell just outside separatrix
xcwi(0:nx+1)        _real   [m]     #x-coord. cell center inner (PF,core) wall (iy=0)
xfwi(0:nx+1)        _real   [m]     #x-coord. cell face inner (PF,core) wall (iy=0)
xfpf(0:nx+1)        _real   [m]     #x-coord. of cell face on p.f. wall (iy=0)
xcpf(0:nx+1)        _real   [m]     #x-coord. of cell center on p.f. wall (iy=0)
xcwo(0:nx+1)        _real   [m]     #x-coord. of cell center on outer wall (ny+1)
xfwo(0:nx+1)        _real   [m]     #x-coord. of cell face on outer wall (ny+1)

yyc(0:ny+1)         _real   [m]     #y-coordinate of the center of the (ix,iy)
                                     #primary cell at midplane
yyf(0:ny+1)         _real   [m]     #coordinate of the y-face of the (ix,iy)
                                     #primary cell at midplane
yylb(0:ny+1,1:nxpt) _real   [m]     #radial dist. on left boundary from septrx.
yyrb(0:ny+1,1:nxpt) _real   [m]     #radial dist. on right boundary from septrx.
xcv(0:nx+1)         _real
xfv(0:nx+1)         _real
psinormc(0:ny+1)    _real           #psi at cell center; norm to unity on septrx
psinormf(0:ny+1)    _real           #psi at cell face; norm to unity on septrx
rrv(0:nx+1,0:ny+1)  _real
volv(0:nx+1,0:ny+1) _real
hxv(0:nx+1,0:ny+1)  _real   [1/m]   #harmonic average of neighboring gx's
                                     # grid (not vel. points); not used
syv(0:nx+1,0:ny+1)  _real
sygytotc             real  [1/m**3] #total volume for iy=1 cell of core region
area_core            real  [1/m**3] #area of core boundary at iy=0
ghxpt                real   [m^-1]  #1/(horiz dist.) between x-pt velocity centers
gvxpt                real   [m^-1]  #1/(vert. dist.) between x-pt velocity centers
sxyxpt               real   [m^2]   #average velocity cell area touching x-point
ghxpt_lower          real   [m^-1]  #1/(horiz dist.) between x-pt velocity centers
gvxpt_lower          real   [m^-1]  #1/(vert. dist.) between x-pt velocity centers
sxyxpt_lower         real   [m^2]   #average velocity cell area touching x-point
ghxpt_upper          real   [m^-1]  #1/(horiz dist.) between x-pt velocity centers
gvxpt_upper          real   [m^-1]  #1/(vert. dist.) between x-pt velocity centers
sxyxpt_upper         real   [m^2]   #average velocity cell area touching x-point
linelen              real   [m]     #field-line length on iy = (iysptrx + ny) / 2
isxptx(0:nx+1,0:ny+1) _integer      #=0 if x-face of this cell touches an x-point
isxpty(0:nx+1,0:ny+1) _integer      #=0 if y-face touches Xpt;=-1 at iysptrx+1 &
                                    #ixpt1+1 or ixpt2+1 for dupdy in Ti vis heat
lcon(0:nx+1,0:ny+1)   _real  [m]    #par dist per pol circ/bwt plates (conn leng)
lconi(0:nx+1,0:ny+1)  _real  [m]    #banana-width corrected ion conn length
lcone(0:nx+1,0:ny+1)  _real  [m]    #banana-width corrected elec conn length
lconneo(0:nx+1,0:ny+1) _real [m]    #neoclassical conn lngth; Int(dx/[2pi*rr])
epsneo(0:nx+1,0:ny+1) _real  [ ]    #inverse aspect ratio
isixcore(0:nx+1)      _integer      #=1 if iy=0 bdry is core; =0 otherwise
dxmin         real /1.e-20/  [m]    #minimum dx & maximum gx=1/dx

***** Comgeo_g:
#Global variables in common -- geometry; used for parallel domain decomp
lcong(0:nx+1,0:ny+1)  _real [m]   #glob par dist per pol circ/bwt plt(conn leng)
lconig(0:nx+1,0:ny+1) _real [m]   #glob banana-wid corrected ion conn leng
lconeg(0:nx+1,0:ny+1) _real [m]   #glob banana-wid corrected elec conn leng

***** Noggeo:
#Geometrical variables used for the nonorthogonal grid stencil
vtag(0:nx+1,0:ny+1)    _real  [ ]  #angle at upper-right vertex; 0 rad. is orthg
                                   #clockwise rot. of face => vtag > 0
angfx(0:nx+1,0:ny+1)   _real  [ ]  #angle on x-face;=.5*(vtag(,iy)+vtag(,iy-1))
fxm(0:nx+1,0:ny+1,0:1) _real  [ ]  #frac. of (ix-1,iy+k) pt used for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fx0(0:nx+1,0:ny+1,0:1) _real  [ ]  #frac. of (ix,iy+k) pt used for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fxp(0:nx+1,0:ny+1,0:1) _real  [ ]  #frac. of (ix+1,iy+k) pt for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fym(0:nx+1,0:ny+1,0:1) _real  [ ]  #frac. of (ix+k,iy-1) pt used for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fy0(0:nx+1,0:ny+1,0:1) _real  [ ]  #frac. of (ix+k,iy) pt used for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fyp(0:nx+1,0:ny+1,0:1) _real  [ ]  #frac. of (ix+k,iy+1) pt for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fxmy(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix-1,iy+1-k) pt for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fxpy(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1,iy+1-k) pt for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fymx(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1-k,iy-1) pt for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fypx(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1-k,iy+1) pt for y deriv,ave
                                   # for (ix,iy) cell where k is third index
fxmv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix-1,iy+k) pt used for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fx0v(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix,iy+k) pt used for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fxpv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1,iy+k) pt for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fymv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+k,iy-1) pt used for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fy0v(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+k,iy) pt used for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fypv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+k,iy+1) pt for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fxmyv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix-1,iy+1-k) pt for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fxpyv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1,iy+1-k) pt for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fymxv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1-k,iy-1) pt for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
fypxv(0:nx+1,0:ny+1,0:1) _real [ ]  #frac. of (ix+1-k,iy+1) pt for y deriv,ave
                                   # for (ix,iy) vel. cell where k is third index
redopltvtag  integer       /0/    #if=1, redo plate vtag if numerically inaccur

***** Timing:
istimingon integer         /1/    # =1 calcs timing data; call wtottim to output
iprinttim  integer         /0/    # =1 to write timing report to terminal
ttotfe     real            /0./   # time spent in full f evaluation (pandf)
ttimpfe    real            /0./   # time spent in full f eval. for impurities
ttotjf     real            /0./   # time spent in Jacobian f evaluation (pandf)
ttimpjf    real            /0./   # time spent in Jac. eval. for impurities
ttmatfac   real            /0./   # time spent in factoring Jacobian
ttmatsol   real            /0./   # time spent in backsolve for using Jacobian
ttjstor    real            /0./   # time spent in storing Jacobian
ttjrnorm   real            /0./   # time spent in row normalization
ttjreorder real            /0./   # time spent in row and column reordering
ttimpc     real	           /0./	  # time spent in impurity calculations
tstart     real            /0./   # initial time from function second
tend       real            /0./   # final time from function second
ttnpg      real            /0./   # time spent in neudifpg
ttngxlog   real            /0./   # time spent for fngx in neudifpg
ttngylog   real            /0./   # time spent for fngy in neudifpg
ttngfd2    real            /0./   # time spent in fd2tra for neudifpg
ttngfxy    real            /0./   # time spent for fngxy in neudifpg

***** Linkbbb:
# information shared by bbb and wdf packages
nxbbb	integer
# number of cells in x
nybbb	integer
# number of cells in y
nycorebbb		integer
# number of core cells in y
nysolbbb		integer
# number of SOL cells in y
nxleg1bbb		integer
# number of cells in x for inboard divertor leg
nxcore1bbb	integer
# number of cells in x for inboard half of core region
nxleg2bbb		integer
# number of cells in x for outboard divertor leg
nxcore2bbb	integer
# number of cells in x for outboard half of core region
geometrybbb	character*8
# same as geometry variable
nibbb(0:nx+1,0:ny+1)	_real
# hydrogen ion density
tibbb(0:nx+1,0:ny+1)	_real
# ion temperature
nebbb(0:nx+1,0:ny+1)	_real
# electron density
tebbb(0:nx+1,0:ny+1)	_real
# electron temperature
vflowxbbb(0:nx+1,0:ny+1)	_real
# ion flow velocity in horizontal (R) direction
vflowybbb(0:nx+1,0:ny+1)	_real
# ion flow velocity in toroidal direction
vflowzbbb(0:nx+1,0:ny+1)	_real
# ion flow velocity in vertical (Z) direction
fnixbbb(0:nx+1,0:ny+1)	_real
# total ion particle current through the east face of each cell
fngysibbb(0:nx+1)	_real
# gas puffing source on innermost flux surface
fngysobbb(0:nx+1)	_real
# gas puffing source on outermost flux surface

***** Timespl:
# Timing data for splines
totb2val   real            /0./   # time spent in spline b2val routine
totintrv   real            /0./   # time spent in spline intrv routine

***** Limiter:
# variables that define a limiter surface
nlimu	integer
# number of user-specified data points on limiter surface
rlimu(1:nlimu)	_real	[m]
# radial position of user-specified data point on limiter surface
zlimu(1:nlimu)	_real	[m]
# vertical position of user-specified data point on limiter surface
nptnma		integer
# computed index of point on (rlimu,zlimu) nearest the magnetic axis
rptnma		real	[m]
# computed radial position of point on limiter nearest the magnetic axis
zptnma		real	[m]
# computed vertical position of point on limiter nearest the magnetic axis
dslims		real	[m]	/.0001/
# minimum 'gap' in poloidal mesh that defines start point for limiter b.c.'s
nsplit1			integer
# dimension of arrays r,zsplit1
rsplit1(nsplit1)	_real	[m]
# radial position of data points on the inboard side of the
# interface between inboard and outboard regions of the mesh
zsplit1(nsplit1)	_real	[m]
# vertical position of data points on the inboard side of the
# interface between inboard and outboard regions of the mesh
nsplit2			integer
# dimension of arrays r,zsplit2
rsplit2(nsplit2)	_real	[m]
# radial position of data points on the outboard side of the
# interface between inboard and outboard regions of the mesh
zsplit2(nsplit2)	_real	[m]
# vertical position of data points on the outboard side of the
# interface between inboard and outboard regions of the mesh

***** Multicharge:
# parameters and arrays for multi-charge-model rate tables
ntev		integer	/101/	# number of temperature values
nz		integer	/2/	# number of charge states
tevb(ntev)		_real   # temperatures for rate arrays
rsi(ntev,0:nz-1)	_real	# ionization rate
rre(ntev,1:nz)		_real	# recombination rate
rpwr(ntev,0:nz)		_real	# radiative power rate
rrcx(ntev,1:nz)		_real	# CX recombination rate
labelrt(1:12)	character*120	# header information from 'b2frates' data file
rtnt	integer		# number of intervals in 'b2frates' temperature data
rtnn	integer		# number of intervals in 'b2frates' density data
rtns	integer		# number of species in a 'b2frates' file
rtnsd	integer		# total number of species from all 'b2frates' files
rtza(0:rtnsd-1)	_real	# atomic charge state
rtzn(0:rtnsd-1)	_real	# nuclear charge state
rtza2(0:rtnsd-1)	_real	# atomic charge state squared
rtt(0:rtnt)	_real	[eV]	# temperature data in 'b2frates' table
rtn(0:rtnn)	_real	[/m**3]	# density data in 'b2frates' table
rtlt(0:rtnt)	_real	# ln(rtt) where rtt[eV] is 'b2frates' temperature data
rtln(0:rtnn)	_real	# ln(rtn) where rtn[/m**3] is 'b2frates' density data
rtlsa(0:rtnt,0:rtnn,0:rtnsd-1)	_real
     # ln(rtsa) where rtsa[m**3/s] is 'b2frates' data for ionization
rtlra(0:rtnt,0:rtnn,0:rtnsd-1)	_real
     # ln(rtra) where rtra[m**3/s] is 'b2frates' data for recombination
rtlqa(0:rtnt,0:rtnn,0:rtnsd-1)	_real
     # ln(rtqa) where rtqa[eV*m**3/s] is 'b2frates' data for elec energy loss
rtlcx(0:rtnt,0:rtnn,0:rtnsd-1)	_real
     # ln(rtcx); rtcx[m**3/s] is 'b2frates' data for c-x on neutral hydrogen
iscxfit	 real /1./ +input	# flag for C ion c-x on neutral H
     #=0 use analytic forms in Braams' rate package
     #=1 use polynomial fit to C.F. Maggi curves (1997)
     #=2 same as =1, except Z=1 has lower rate from Pigarov
isrtndep  integer /1/ +input   # flag for turning on ne-dep atomic rates when
     # ismctab=2 (impurities) and/or istabon=16 (hydrogen)
mcfformat(1:12) integer /12*0/ # flag read from mcfilename rate files to distingu
     # pre-2012 format (0) from 2012 format (1) for rates;
     # assoc. with how ebindz is used.
chgstate_format(0:rtnsd-1) _integer +input # flag denoting data format used for
     # each charge state; see mcfformat
ispradextrap integer /0/ +input #=1 extrapolated neutral prad below table min

***** Fitdata:
# Arrays and variables for the tanh and b-spline fits to DIII-D radial profiles
isprof_coef       integer /0/   #if=0, use interp of data; if=1, use coeffs for fits
ncoefne_tanh	  integer /1/   #Number of tanh coeffs for ne
ncoefte_tanh	  integer /1/   #Number of tanh coeffs for te
numt_bs		  integer /1/   #Size of b-spline knot vector t
numc_bs		  integer /1/   #Number of b-spline coeffs
numk_bs      	  integer /1/   #Order of b-spline
fcoefne_tanh(ncoefne_tanh) _real #Coeffs for tanh fit of ne
fcoefte_tanh(ncoefte_tanh) _real #Coeffs for tanh fit of te
fit_paramne_tanh  character*8   #Special ne fit parameter; rarely used
fit_paramte_tanh  character*8   #Special te fit parameter; rarely used
fit_t_bs(numt_bs)     _real 	#B-spline knot vector t
fcoef_bs(numc_bs)     _real   	#B-spline coeffs
nefit(0:ny+1,2)       _real [m**-3]/0/ #ne from DIII-D tanh fit
tefit(0:ny+1,2)       _real [keV]  /0/ #Te from DIII-D tanh fit
tifit(0:ny+1,2)       _real [keV]  /0/ #Ti from DIII-D B-spline fit
nefituse(0:ny+1)      _real /0/ #ne fit used combin sets 1&2 with fitfrac1
tefituse(0:ny+1)      _real /0/ #te fit used combin sets 1&2 with fitfrac1
tifituse(0:ny+1)      _real /0/ #ti fit used combin sets 1&2 with fitfrac1
fitfrac1              real /1./ #fraction of data set 1 used for interp fits
dumfit(0:ny+1)        _real /0/ #1D array to transf data from tanh routines
isdatfmtnew       integer /1/   #d3d data fit opts tanh te, ne;=0 before 2009;
                                #=1 up to Oct09;=2 from Oct09c
psishift	  real   /0./   #shift in psi=1 used for tanh, Bspline fits
isdndtfitdat     integer /0/    #=0 if one data set; =1 if two for sets for dt
ifitset          integer /1/    #2nd index for nefit etc; 1st or 2nd data set
tim_interval_fit  real [s] /3.e-2/ #time between profile data fits
tim_chng_max      real [s] /1./ #max time assigned to d/dt for data fits
taudndt(0:ny+1)  _real [s] /1.e50/ #calc dens rise time for data fits
taudeedt(0:ny+1) _real [s] /1.e50/ #calc elec energy rise time for data fits
taudeidt(0:ny+1) _real [s] /1.e50/ #calc ion energy rise time for data fits
isprofvspsi      integer   /1/    #=1 if expt fit vs psi, =0 if fit vs r
epsi_fit(num_elem) _real  /0./    #expt pol flux data values
psi_s(num_elem)    _real  /0./    #value of epsi_fit shifted by psishift
yyc_fit(num_elem)  _real  /0./    #radial distance from outermp for expt prof
eprofile_fit(num_elem) _real /0./ #expt profile values at epsi_fit

***** Subs2:
# Subroutines that can be called from parser
xerrab(msg:string)                  subroutine
        # interface to remark and kaboom
readne_dat(fname:string)			subroutine
        # reads tanh coeffs for radial elec density profile fits 
        # in fname   the filename
readte_dat(fname:string)			subroutine
        # reads tanh coeffs for radial elec temp profile fits 
        # in fname   the filename
readti_dat(fname:string)			subroutine
        # reads spline coeffs for radial ion temp profile fits 
        # in fname   the filename
fit_neteti()                                    subroutine
        # fits ne and te to radial tanh profile; ti to radial spline
        # in fname   the filename
tanh_multi(i:integer,a:real,j:integer,b:real,fname:string,d:real) subroutine
        # provides general tanh fuction given input coeff
        # in i       number of coeffs
        # in a       coeff values
        # in j       number of eval pts
        # in b       psi values at eval pts
        # in fname   the filename
        # out d      values of fit

***** Flags:
# Common flags used by UEDGE
iprint  integer     /1/     # Flag controlling whether to be verbose or not
