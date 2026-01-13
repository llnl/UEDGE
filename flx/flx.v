flx
     # This is the flux contouring package that converts the EQDSK output
     # from the EFIT code into flux surface data that forms the basic structure
     # of the plasma zone grid used by DEGAS and B2/UEDGE.
{
MFLX = 100       # maximum number of flux surfaces
NOREGS = 2	 # number of horizontal regions
JDIM = 2 * MFLX
NPTS = 2000      # maximum number of data points on a contour
NSEARCH = 4      # maximum number of search regions
}

***** Com_Dim_Vars hidden:    
dim_vars_hidden     integer    # Do not edit this group. It is used to build
                               # the Basis version of the code. 
***** Grd_Dim_Vars hidden:    
grd_vars_hidden     integer    # Do not edit this group. It is used to build
                               # the Basis version of the code. 

***** Dimflx:
     # dimensioning parameter used in flx package only
nsearch integer /NSEARCH/
     # maximum number search regions for contouring package

***** Flxin:
istchkon	integer	/0/ +gridgen
     # switch for imposing limits on the polar angle about (rmagx,zmagx)
     # of (r,z) points on flux contours.
     #	=0 no limits
     #	=1 limits are defined by theta_split, thetax, dtheta_exclude,
     #             dtheta_overlap_sol and dtheta_overlap_pf
isthmmxn	integer	/1/ +gridgen
     # switch that controls expressions for thetamin and thetamax when
     # istchkon=1 for 'dnbot' configurations:
     # =0     uses same theta_split for in/outboard halves
     # =1     uses pi/twopi in place of theta_split for in/outboard halves
dtheta_exclude(1:2)	real	[radians]	/2*1.5/ +gridgen
     # angular width of region where SOL flux contours are excluded
     # index 1 refers to inboard flux contours; 2 refers to outboard contours
dtheta_overlap_sol(1:2)	real	[radians]	/2*0.5/ +gridgen
     # angular width over which SOL flux contours can overlap
     # with flux contours in the adjacent region.
     # index 1 refers to inboard flux contours; 2 refers to outboard contours
dtheta_overlap_pf(1:2)	real	[radians]	/2*0.25/ +gridgen
     # angular width over which p.f. flux contours can overlap
     # with flux contours in the adjacent region.
thetax			real	[radians]
     # computed poloidal angle of xpoint relative to magnetic axis
thetamin(1:2)		real	[radians]
     # computed minimum poloidal angle of flux contour data
     # index 1 refers to inboard flux contours; 2 refers to outboard contours
thetamax(1:2)		real	[radians]
     # computed maximum poloidal angle of flux contour data
     # index 1 refers to inboard flux contours; 2 refers to outboard contours
theta1fac	real	/1.0/ +gridgen
     # multiplicative factor for theta1 boundary of inner-half pf region
theta2fac	real	/0.0/ +gridgen
     # multiplicative factor for theta2 boundary of outer-half pf region
ymax1fac        real    /1.0/ +gridgen
     # multiplies zmagx to set upper bound on search region 1
ymax2fac        real    /3.0/ +gridgen
     # multiplies zseps to set upper bound on search region 2
ymax3fac        real    /1.0/ +gridgen
     # multiplies zseps to set upper bound on search region 3
ymax4fac        real    /1.0/ +gridgen
     # multiplies zseps to set upper bound on search region 4
imagx		integer
     # horizontal (R) index of the refined-EFIT cell containing the magnetic axis
jmagx		integer
     # vertical (Z) index of the refined-EFIT cell containing the magnetic axis
iseps		integer
     # horizontal (R) index of the refined-EFIT cell containing the x-point
jseps		integer
     # vertical (Z) index of the refined-EFIT cell containing the x-point
icutoff1	integer
     # horizontal (R) index of the refined-EFIT cell containing xcutoff1
jcutoff1	integer
     # vertical (Z) index of the refined-EFIT cell containing ycutoff1
slpyt		real	/1.0/ +gridgen
# scaling factor for radial mesh near double-null separatrices -
#    - mesh size is scaled by slpyt near inner separatrix
#    - mesh size is scaled by 1/slpyt near outer separatrix
slp2fac		real	/1.0/ +gridgen
# scale factor for radial mesh in asymmetric double-null configurations
# (see function rho3dn for details)
# = 1. for uniform mesh between separatrices
# < 1. for finer mesh at innermost separatrix
# > 1. for coarser mesh at innermost separatrix
slp3fac		real	/1.0/ +gridgen
# scale factor for radial mesh in asymmetric double-null configurations
# (see function rho3dn for details)
# = 1. for uniform mesh between separatrices
# < 1. for finer mesh at outermost separatrix
# > 1. for coarser mesh at outermost separatrix
psifac		real	/1.0005/ +gridgen
# multiplicative factor in normalized separatrix flux
# to ensure that contour is just outside the x-point.
psi0sep1	real
# normalized flux value at lower x-point
psi0sep2	real
# normalized flux value at upper x-point
psi0max_inner		real	/1.07/ +regrid +gridgen
# normalized flux value at wall on inboard side of magnetic axis
psi0max_outer		real	/1.07/ +regrid +gridgen
# normalized flux value at wall on outboard side of magnetic axis
psi0min2_upper		real	/0.98/ +regrid +gridgen
# normalized flux value at private flux boundary near upper x-point
psi0min2_lower		real	/0.98/ +regrid +gridgen
# normalized flux value at private flux boundary near lower x-point
psi0min1	real	/0.98/ +regrid +gridgen
     # normalized flux value at innermost core flux surface
psi0min2	real	/0.98/ +regrid +gridgen
     # normalized flux value at innermost private flux surface
psi0sep		real	/1.00001/ +regrid +gridgen
     # normalized flux value at separatrix flux surface (just slightly outside)
psi0max			real	/1.07/ +regrid +gridgen
     # normalized flux value at wall on outboard side of magnetic axis
psi0lim		real +gridgen
     # normalized flux value where core radial mesh can be concentrated
sfaclim		real	/1.0/ +gridgen
     # multiplicative factor for slope of psi0(iy) distribution at psi0lim
     #	sfaclim = 1	-->	nearly uniform
     #	sfaclim < 1	-->	more concentrated near psi0lim
alfcy_inner	real	/.0001/ +regrid +gridgen
     # exponential factor for  distribution of SOL flux contours (inboard)
     #	alfcy << 1	-->	uniform
     #	alfcy >> 1	-->	concentrated near separatrix
alfcy		real	/.0001/ +regrid +gridgen
     # exponential factor for  distribution of SOL flux contours (outboard)
     #	alfcy << 1	-->	uniform
     #	alfcy >> 1	-->	concentrated near separatrix
xoverlap(2)	real	/5.0,4.0/ +gridgen +gridgen
     # overlap parameters for inboard and outboard flux contours
     # (in units of dxefit)
rho(0:nym)	_real
     # temporary array for normalized flux values on B2/UEDGE flux surfaces
tflx(0:nym)	_real
     # temporary array for indicies of B2/UEDGE flux surfaces
psitop(1:jdim)	_real +gridgen
     # array of un-normalized pol. mag. flux values across the midplane 
psibot(1:jdim)	_real +gridgen
     # array of un-normalized pol. mag. flux values across both divertor plates
iseqdskr	integer	/0/ +gridgen
     # switch (=1) for reflecting eqdsk data about the midplane;
     # used with geometry='uppersn' and for upper half-mesh of geometry='dnull'
kymesh		integer	/1/ +gridgen
     # option flag for setting flux contour values that define radial mesh
     # =0 user sets arrays psitop (midplane) and psibot (div plates)
     # =1 psitop and psibot are analytic forms (see rho1); adjust via alfcy
     # =2 psitop and psibot are analytic forms (see rho5); uniform core
xcutoff1	real	/0./	[m] +gridgen
     # inboard cutoff radius for flux contours
ycutoff1	real	/0./	[m] +gridgen
     # lower vertical cutoff for flux contours
mdsefit		integer	/0/ +gridgen
     # Set this flag to 1 if the data is read in from Mdsplus
     # Prevents inflx from calling readefit.

***** Workdn:
     # some working arrays for full double-null configurations
psi0_mp_inner(0:nym)	_real +gridgen
# normalized flux values at radial mesh points of inner midplane
psi0_mp_outer(0:nym)	_real +gridgen
# normalized flux values at radial mesh points of outer midplane
psi0_dp_lower_inner(0:nym)	_real +gridgen
# normalized flux values at radial mesh points of lower inner divertor plate
psi0_dp_lower_outer(0:nym)	_real +gridgen
# normalized flux values at radial mesh points of lower outer divertor plate
psi0_dp_upper_inner(0:nym)	_real +gridgen
# normalized flux values at radial mesh points of upper inner divertor plate
psi0_dp_upper_outer(0:nym)	_real +gridgen
# normalized flux values at radial mesh points of upper outer divertor plate

***** Inpf0:
# flux contour control variables for each sub-region, nsearch:
plflux0(jdim,nsearch)	_real
     # poloidal flux contour values
ncmin0(nsearch)	_integer
     # flux index of first poloidal flux curve
ncmax0(nsearch)	_integer
     # flux index of last poloidal flux curve
iserch0(nsearch) _integer
# horizontal (x) index of starting point for flux contour search
# on refined EFIT grid
jserch0(nsearch) _integer
# vertical (y) index of starting point for flux contour search
# on refined EFIT grid
istepf0(nsearch) _integer
# horizontal step size (on refined EFIT mesh) for contour search path
jstepf0(nsearch) _integer
# vertical step size (on refined EFIT mesh) for contour search path
xminf0(nsearch)      _real   [m]
     # radial cutoff for contours at inner edge of efit domain
xmaxf0(nsearch)      _real   [m]
     # radial cutoff for contours at outer edge of efit domain
yminf0(nsearch)      _real   [m]
     # vertical cutoff for contours at upper edge of efit domain
ymaxf0(nsearch)      _real   [m]
     # vertical cutoff for contours at lower edge of efit domain
istcvon		integer	/0/ +gridgen
     # obsolete option istcvon=1; use altsearch=2 instead.
altsearch	integer /0/ +gridgen
     # flag to change search path for private flux surfaces
     # altsearch=0	--> search vertically up toward x-point
     # altsearch=1	--> search vertically down from x-point
     # altsearch=2	--> search diagonally down and in from x-point
isetpath     integer /0/ +gridgen
     # set path for finding flux surfaces in midplane
     # isetpath=0    --> search radially outward from left boundary
     # isetpath=1    --> search radially inward from magnetic axis

***** Inpf:
# local variables for flux contouring; see inpf0 description
ncmin   integer
     # flux index of first poloidal flux curve
ncmax   integer
     # flux index of last poloidal flux curve
iserch(jdim)  _integer
# horizontal (x) index of starting point for flux contour search
# on refined EFIT grid
jserch(jdim)  _integer
# vertical (y) index of starting point for flux contour search
# on refined EFIT grid
istepf	integer
# horizontal step size (on refined EFIT mesh) for contour search path
jstepf	integer
# vertical step size (on refined EFIT mesh) for contour search path
leadir(jdim)  _integer
# specifies initial direction for flux contour search, subsequent movement :
     # leadir = 1    upward vertical search - contouring right, then left
     # leadir = 2    rightward horizontal search - contouring down, then up
     # leadir = 3    upward vertical search - contouring left, then right
     # leadir = 4    rightward horizontal search - contouring up, then down
     # leadir = 5    downward vertical search - contouring left, then right
     # leadir = 6    downward vertical search - contouring right, then left
ncmin1  integer
     # flux index of first poloidal flux curve
ncmax1  integer
     # flux index of last poloidal flux curve

***** Polflx:
     # flux contouring variables on the refined grid
plflux(jdim)     _real
     # temporary array of poloidal flux values in a search region
mrfac	integer	/4/ +gridgen
     # mesh refinement factor relative to EFIT
nx4     integer
     # number of grid surfaces in x-direction on refined grid
ny4     integer
     # number of grid surfaces in y-direction on refined grid
x(nx4)  _real
     # radial position of grid surfaces in x-direction on refined grid
y(ny4)  _real
     # vertical position of grid surfaces in y-direction on refined grid
f(nx4,ny4)      _real
     # poloidal magnetic flux values on refined grid
ijumpf(jdim)	_integer
     # core/p.f. discontinuity on flux contour j occurs between
     # data points ijumpf(j) and ijumpf(j)+1
dsjumpf		real	/0.1/	[m] +gridgen
     # criterion for "jump" discontinuity in core/p.f. contours
ilast(jdim)     _integer
     # comment needed
xcn(npts)        _real   [m]
     # comment needed
ycn(npts)        _real   [m]
     # comment needed
imin    integer
     # comment needed
imax    integer
     # comment needed
istart  integer
     # comment needed
jmins	integer
     # comment needed
jmaxs	integer
     # comment needed
jstart  integer
     # comment needed
ncontr  integer
     # comment needed
pfr     real
     # comment needed
twopie  real
     # comment needed
xminf   real   [m]
     # comment needed
xmaxf   real   [m]
     # comment needed
yminf   real   [m]
     # comment needed
ymaxf   real   [m]
     # comment needed
rs_com	real [m]  #rm output from subr. findstrike via common for python
zs_com	real [m]  #zm output from subr. findstrike via common for python 

***** Rho:
     # various analytic forms for choosing rho versus t
rho1(t,rho,nt,t1,t2,t3,r1,r2,r3,alf)        subroutine
     # Defines rho(t) as a rational function of the form (f+gt)/(p+qt) on
     # the interval t1 < t < t2 (the plasma core region) and
     # an exponential function of the form sinh(at+b) on
     # the interval t2 < t < t3 (the scrape-off layer region) with
     # continuous first derivative at t2.
     # input array is t(0:nt) and output array is rho(0:nt).
     # alf sets the strength of the exponential behavior which
     # becomes linear in the limit alf ---> 0.
rho1dn(t,rho,nt,t1,t2,t3,t4,r1,r2,r3,r4,alf)        subroutine
     # Defines rho(t) for unbalanced double-null configurations:
     # rho(t) is a rational function of the form (f+gt)/(p+qt) on
     # the interval t1 < t < t2 (the inner core region) and also on the
     # interval t2 < t < t3 (between separatrices), and has an exponential
     # form sinh(at+b) on the interval t3 < t < t4 (the scrape-off layer
     # region) with continuous first derivatives at t2 and t3.
     # t(0:nt) is an input array and rho(0:nt) is the output array.
     # alf sets the strength of the exponential behavior which
     # becomes linear in the limit alf ---> 0.
rho2dn(t,rho,nt,t1,t2,t3,t4,r1,r2,r3,r4,fac)        subroutine
     # Defines rho(t) for unbalanced double-null configurations:
     # rho(t) is a piece-wise rational function of the form (f+gt)/(p+qt)
     # on the intervals t1 < t < t2 (the inner core region), t2 < t < t3
     # (between separatrices), and  t3 < t < t4 (the outer scrape-off layer
     # region) with continuous first derivatives at t2 and t3.
     # t(0:nt) is an input array and rho(0:nt) is the output array.
     # fac is a scaling factor for the slopes at t2 and t3:
     #      rho'(t2)=((r3-r2)/(t3-t2))*fac
     #      rho'(t3)=((r3-r2)/(t3-t2))/fac
rho3dn(t,rho,nt,t1,t2,t3,t4,r1,r2,r3,r4,slp2fac,slp3fac,r2p,r3p)  subroutine
     # Defines rho(t) for unbalanced double-null configurations:
     # rho(t) is a cubic function on the interval t2 < t < t3 (between
     # separatrices) and a piece-wise rational function of the form
     # (f+gt)/(p+qt) on the intervals t1 < t < t2 (the inner core region)
     # and  t3 < t < t4 (the outer scrape-off layer region) with
     # continuous first derivatives at t2 and t3.
     # t(0:nt) is an input array and rho(0:nt) is the output array.
     # slp2fac & slp3fac are scaling factors for the slopes at t2 and t3:
     #      rho'(t2)=((r3-r2)/(t3-t2))*slp2fac
     #      rho'(t3)=((r3-r2)/(t3-t2))*slp3fac
rho1l(t,rho,nt,t1,t2,r1,r2,r1p)		subroutine
     # rho(t) is a rational function of the form
     # (f+gt)/(p+qt) on the interval t1 < t < t2
     # with specified derivative at the left boundary.
rho1r(t,rho,nt,t1,t2,r1,r2,r2p)		subroutine
     # rho(t) is a rational function of the form
     # (f+gt)/(p+qt) on the interval t1 < t < t2
     # with specified derivative at the right boundary.
rho2(t,rho,nt,t1,t2,t3,r1,r2,r3)        subroutine
     # Defines rho(t) as quadratic on the interval t1 < t < t2
     # and linear on the interval t2 < t < t3 with continuous
     # first derivative at t2.
     # nt specifies the length of input array t and output array rho.
rho3(t,rho,nt,t1,t2,t3,r1,r2,r3)        subroutine
     # Defines rho(t) as cubic on the interval t1 < t < t2
     # and linear on the interval t2 < t < t3 with continuous
     # first and second derivatives at t2.
     # nt specifies the length of input array t and output array rho.
rho4(t,rho,nt,t1,t2,t3,r1,r2,r3,s2)        subroutine
     # Defines rho(t) as a cubic on the interval t1 < t < t3 
     # with slope at t2 modified relative to quadratic 
rho5(t,rho,nt,t1,t2,t3,r1,r2,r3,r2p)        subroutine
     # Defines rho(t) as a piece-wise rational function of the form
     # (f+gt)/(p+qt) on each of the intervals t1 < t < t2 (core region)
     # and t2 < t < t3 (SOL region) with specified derivative rho'(t2).
***** Efit:
     # subroutines that the user can call directly
aeqdsk          subroutine
     # This subroutine reads Efit data from file aeqdsk.
neqdsk          subroutine
     # This subroutine reads Efit data from file neqdsk.
readefit        subroutine
     # This subroutine calls neqdsk and aeqdsk and procefit.
procefit        subroutine
     # This subroutine assumes the equivalent of neqdsk and aeqdsk have
     # already been called. It then defines the vertical 
     # coordinate with the z=0 surface located at the bottom of the mesh
     # rather than at the midplane. 
refine		subroutine
     # This subroutine subdivides the original EFIT grid and evaluates the
     # original spline fit on the refined grid (see mrfac)
contours(ns:integer)	subroutine
     # This subroutine computes the data for flux contours in search region ns
flxrun		subroutine
     # This is the main driver routine for the flx package.
flxfin		subroutine
     # This subroutine writes the flx-grd file.
inflx		subroutine
     # This routine selects the psi values for the flux surfaces
theta_ok(r:real, z:real, n:integer) logical function
     # test for theta [polar angle about (rmagx,zmagx)] within acceptable range
     # in contouring region n where:
     #	n=1	inboard SOL and core
     #	n=2	inboard private flux
     #	n=3	outboard SOL and core
     #	n=4	outboard private flux
efitvers(vmonth:integer, vday:integer, vyear:integer)	integer function
     # test for updated version of EFIT (on or after 24 May 1997)
     # =1 means new version
     # =0 means old version
findstrike(js:integer, rs:real, zs:real)	subroutine
     # Finds the first intersection of a poloidal flux surface with
     # the limiter surface from EFIT; input is flux surface index js;
     # outputs are (rs,zs); assumes x,ycurve and x,ylim data exist.

