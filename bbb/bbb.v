bbb
{
nispmx = 31	# maximum number of ion species
ngspmx = 6	# maximum number of gas species; also must set in com.v
nmcmx = 12	# maximum number of EIRENE test species in data file 'fort.44'
ndomainmx = 32	# maximum number of domains for domain decomposition
nxptmx = 2      # maximum number of x-points in R-Z domain
ndcsmx = nispmx*nxptmx	# data dimension for csfaclb and csfacrb
nispmxngspmx = nispmx*ngspmx # tot numb ion*gas species
nstramx = 10 # maximum number of strata for MC neutrals code
}

***** Com_Dim_Vars hidden:    
dim_vars_hidden     integer    # Do not edit this group. It is used to build
                               # the Basis version of the code. 

***** Math_problem_size:
neqmx		integer		# number of math. eqns to be solved/integrated
numvar		integer		# number of physical variables per cell
numvarbwpad     integer    /1/  # add to numvar for bandwidth calc;safety param

***** UEpar:
#  Parameters and switches for the simulation
csfaclb(nispmx,nxptmx)  real /ndcsmx*1./  +input #frac of cs used for Bohm sheath b.c.
csfacrb(nispmx,nxptmx)  real /ndcsmx*1./  +input #frac of cs used for Bohm sheath b.c.
csfacti   real            /1./  +input #Bohm speed = sqrt((te+csfacti*ti)/mi)
cslim     real            /1./  +input #frac of cs used for limiter Bohm sheath b.c.
dcslim    real            /0./  +input 
                                #reduce sonic flow at limiter by the factor
                                #cslim*[1-exp(-(iy-iy_lims+1)/dcslim)]
islnlamcon integer        /0/   +input #=0, loglambda=Braginskii;if=1,loglambda=lnlam
lnlam     real            /12./ +input #Coulomb log;shouldn't be constant
methe     integer         /33/  +input #elec. eng. eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
methu     integer         /33/  +input #ion mom. eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
methn     integer         /33/  +input #ion cont. eqn: 22-harmonic average, 33-uw
methi     integer         /33/  +input #ion eng. eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
methg     integer         /33/  +input 
                                #neut. gas eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
                                #66 nonorth. log intrp, 77 nonorth. 1/ng intrp
methp     integer         /33/  +input #potential eqn: 22-cd, 33-uw, 44-hyb, 55-p-law
isgxvon   integer         /0/   +input #=0 uses gx in fmix; =1 for harmonic ave of gxf
ishavisy  integer         /1/   +input #=1 uses harmonic ave for conxi up
cfaccony  real            /1./  +input #scales conxi for up
isvylog   real            /0/   +input #=0 for vy~(1/n)dn/dy; =1 for vy~d(log(n))/dy
isintlog  integer         /0/   +input #nonog logrithm interp for remaining terms
concap    integer         /0/   +input #=1 fixes Te and Ti to afix for thermal cond.
convis    integer         /0/   +input #=1 fixes Te to afix for ion viscosity
cniatol   real            /1./  +input #multiplier for atol for ni
cngatol   real            /1./  +input #multiplier for atol for ng
cupatol   real            /1./  +input #multiplier for atol for up
cteatol   real            /1./  +input #multiplier for atol for te
ctiatol   real            /1./  +input #multiplier for atol for ti
cphiatol  real            /1./  +input #multiplier for atol for phi
tolbf     real            /1./  +input #multiplier for atol&rtol for the boundary eqns
tadj      real            /10./ +input #reduces time step by 1/tadj if iopts=1
icnuiz    integer         /0/   +input #=1 constant ioniz. freq., cnuiz; =2 freezes
icnucx    integer         /0/   +input 
                                #=0, var nucx;=1 const. nucx=cnucx;
                                # =2, use sigcx, so nucx~(Tg)**.5
cnuiz     real [1/s]      /5.e+4/ +input #constant ioniz. freq. for icnuiz=1
cnucx     real [1/s]      /1.e+0/ +input #constant charge exhange freq. for icnucx=1
isrecmon  integer         /0/  +input 
				 #flag to turn-on recombination (yes=1); use
				 #cfrecom to turn-off recomb after isrecmon was on
cfrecom   real            /1./ +input 
                                  #scale factor multiplying recombination freq.
igas      integer         /0/     +input #=1 invokes local rate eqn. for ng
ngbackg(ngspmx) real [1/m**3] /ngspmx*1.e14/ +input 
                                  #background gas density
ingb      integer         /2/     +input 
                                  #background gas source=nuiz*ngbackg*
                                  #                  (.9+.1*(ngbackg/ng)**ingb)
inflbg    integer         /4/     +input 
                                  #expon to force flalfg large near ng~ngback
                				  #ex:flalfgx,y*(1.+(cflgb*ngbackg/ng)**inflbg)
cflbg     real            /10./   +input #scaling fac for flalfgx,y using inflbg
facngbackg2ngs(ngspmx) real /ngspmx*1.e-8/ +input 
                                  #fraction of ngbackg add to initial ngs
nzbackg(nispmx) real [1/m**3] /nispmx*1.e9/ +input #background impurity density
inzb      integer         /2/     +input
                                  #background impurity source=nuiz*nzbackg*
                                  #                 (.9+.1*(nzbackg/nzi)**ingb)
facnzbackg2nis(nispmx) real /nispmx*1.e-8/ +input 
                                  #fraction of nzbackg add to initial nis
upclng(nispmx) real [m/s] /nispmx*1.e8/ +input 
                                  #max ion vel at beginning of iteration
facupclng2ups(nispmx) real /nispmx*1.e-8/ +input 
                                  #fraction of upclng subtract from initial ups
tebg      real [eV]    /1.e-20/   +input #backgrd elec eng sor to limit te~tebg
tibg      real [eV]    /1.e-20/   +input #backgrd ion eng sor to limit te~tebg
iteb      integer         /2/     +input #exponent of (tebg*ev/te)**iteb for bkg sor
temin     real [eV]      /0.03/   +input #min value of te allow; if less, reset to
temin2    real [eV]      /0.03/   +input #soft floor with te=sqrt[te**2+(temin2*ev)**2]
tgmin     real [eV]      /0.03/   # min value of tg allowed
pwrbkg_c  real [W/m**3]  /1.e3/   +input #const background factor in pwrebkg express
pwribkg_c real [W/m**3]  /1.e3/   +input #const background factor in pwribkg express
cfwjdotelim real         /1./     +input #factor scaling reduction of wjdote if te<tebg
nlimix(nispmx) real  /nispmx*0./  +input #factor to prevent ion density pump out in x
nlimiy(nispmx) real  /nispmx*0./  +input #factor to prevent ion density pump out in y
nlimgx	  real            /0./    +input #factor to prevent ion density pump out in x
nlimgy	  real            /0./    +input #factor to prevent ion density pump out in y
is1D_gbx  integer         /0/     +input #=1 turns on 1-D gas-box model
xgbx      real [m]        /.25/   +input #poloidal location of 1-D gas box
ixgb      integer         /0/     +input #poloidal index of xgbx 1-D gas box (calc)
agdc      real            /1./    +input #exp. decay factor ng from gas-box edge
pcolwid   real [m]        /0./    +input #width of plasma column for 1-D gas-box model
eion      real [eV]       /5./    +input #energy that ionized ion is born with
ediss     real [eV]       /10./   +input
                                  #elec eng lost by mol. dissoc; should = 2*eion
cfdiss    real            /1./    # fraction of neutrals that are from dissociations.
ebind     real [eV]     /13.6/    +input
                                  #binding energy carried by hydrogen ion
afix      real [e]        /50./   +input #Te,i for fixed cond.(concap), visc.(convis)
coef      real            /0.96/  +input #factor for ion viscosity: was 1.92 ???
ce        real            /3.16/  +input #factor for electron thermal conductivity
ci        real            /3.9/   +input #factor for ion thermal conductivity
                #The zeff dependence of ce has been explicitly added in zcoef,
                #thus ce should always be left as 3.16 even if zeff is not 1,
                #provided zeff is less than or equal to 4.
ncrhs     integer
istep     integer
iter      integer
dp1       real
qfl       real
csh       real
qsh       real
mfl       real
msh       real
ro        real
cs        real
ctaue(0:nx+1,0:ny+1,nisp)  _real  #calc factor for elec Coulomb coll
ctaui(0:nx+1,0:ny+1,nisp)  _real  #calc factor for ion Coulomb coll
fxe       real
fxi       real
zcoef     real          #factor (calc) give zeff dependence of elec thermal c.
coef1     real          #factor (calc) for energy equipartion rate
cnurn     real    /1./  +input #scales nurlx rate for ion continuity eqn.
cnuru     real    /1./  +input #scales nurlx rate for ion mom. eqn.
cnure     real    /1./  +input #scales nurlx rate for elec. eng. eqn.
cnuri     real    /1./  +input #scales nurlx rate for ion eng. eqn.
cnurg     real    /1./  +input #scales nurlx rate for gas eqn.
cnurp     real    /1./  +input #scales nurlx rate for potential eqn.
nurlxn    real
nurlxu    real
nurlxe    real
nurlxi    real
nurlxg    real
nurlxp    real
label     character*72 #code name and run time, date, and machine
rnewpot   real    /0./      +input #mixture of fqy=(1-rnewpot)*fqy_old+rnewpot*fqy_new
r0slab	  real [m]/1e-20/   +input #effect. major radius for isnewpot j_r calc in slab
ishymol   integer  /0/      +input #=1 turns on hydr. mol; requires nhgsp=2
te_s_dis  real     /5./     +input #Te shift of ioniz curve to approx dissociation curve
isfqpave  integer  /0/      +input #=0 for lin interp for fqp terms; =1 for simple ave.
svrpkg    character*8 /"nksol"/ +input #use solver pkg daspk,vodpk,nksol,newton
                                #reset to newton if inewton=1
petscoptfile  character*80 /""/     +input #specify options file for petsc code
petscopts     character*2048 /""/   +input #specify options for petsc code
isgpye    integer  /0/ +input #change -vy*dP/dy eng. terms; =1 for old B2; =2 for Knoll
iigsp     integer      # Ion index for neutrals when isupgon=1 (zi(iigsp)=0)
fnnuiz    real    /1./  +solver #fraction of new nuiz used for Jacobian
itrap_negni integer /1/ +solver #flag to trap negative ni condition
itrap_negt integer /1/  +solver #flag to trap negative Te,i condition
itrap_negng integer /1/ +solver #flag to trap negative ng condition
isybdryog integer  /0/  +input #=1 sets fx0, fmx stencil to orthog values at iy=0 & ny
isybdrywd integer  /0/  +input #=1 vy diffusion-only for iy=0 & ny if matwalli,o=1
isxmpog   integer  /0/  +input 
                        #=1 sets fy0, fmy stencil to orthog values at ix=nxc-1
                        # and ix=nxc+1 for geometry='dnbot'
iexclnxc1 integer  /0/  +input 
                        #if=0; include nxc+1 for fee,iytotc if geometry=dnbot;
                        #if=1; exclude nxc+1 for fee,iytotc
ineudif   integer  /2/  +input 
                #=1 gas sub. neudif uses ng, tg for gas vel & fngx->fnix
		        #=2 gas sub. neudifgp uses pg for gas vel & fngx->fnix
		        #=3 gas sub. neudifl use log_ng, tg for gas vel
		       #otherwise, old case has ug=ui (strong cx coupling)
thetar    real    /0./  +input #rotate (R,Z) coordinates by angle theta (degrees)
isbcwdt   integer /0/   +solver #include dtreal in B.C. if isbcwdt=1
ishosor   integer /0/   +input #if=1, integrate hydr. sources over cell; full RHS only
iseesorave real   /0./  +input #cell ave factor; 0 ctr only; 1 5 pt ave elec eng sors
ispsorave real    /0./  +input #cell ave factor; 0 ctr only; 1 5 pt ave of psorg,psor,etc.
fsprd     real /.0625 / +input #fraction of eng. sor. spread to each of 4 neighbors
issyvxpt0 integer /0/   +input #if=1, set syv=0 around x-point; ambig. rad. mom. flux
isrrvave  integer /0/   +input 
                        #if=0, rrv from vertex B's; 
                        #if=1, rrv=0.5*(rr_1+rr_2);
                        #if=2, average of cases 0 and 1
rr_fac    real    /1./ +input #scale factor to multiple rr and rrv
rrmin     real    /0./ +input #min rr used in calc of u_tor & fqy for potential calc.
isdtsfscal integer /0/ +input #if=1, dt is included in sfscal Jac scaling factor
frfqpn    real    /1./ +input #frac. of new fqp at ix=0,nx using grad at ix=1,nx-1
cffqpsat  real    /1./ +input #factor by which fqp can exceed fqpsatlb,rb (sat. cur)
isplflxl  integer /0/  +input #=0, flalfe,i not active at ix=0 & nx;=1 active all ix
isplflxlv integer /0/  +input #=0, flalfv not active at ix=0 & nx;=1 active all ix
isplflxlgx integer /0/ +input #=0, flalfgx not active at ix=0 & nx;=1 active all ix
isplflxlgxy integer /0/ +input #=0, flalfgxy not active at ix=0 & nx;=1 active all ix
islimflxlgx integer /0/ +input #=0, flalfgx,xy not active at ix=ix_lim
iswflxlgy   integer /0/ +input #=0, flalfgy not active at iy=0 & ny;=1 active all iy
isplflxlvgx integer /0/ +input #=0, flalfvgx not active at ix=0 & nx;=1 active all ix
isplflxlvgxy integer /0/ +input #=0, flalfvgxy not active at ix=0 & nx;=1 active all ix
iswflxlvgy  integer /1/ +input #=0, flalfvgy not active at iy=0 & ny;=1 active all iy
isplflxltgx integer /0/ +input #=0, flalfvgx not active at ix=0 & nx;=1 active all ix
isplflxltgxy integer /0/ +input #=0, flalfvgxy not active at ix=0 & nx;=1 active all ix
iswflxltgy  integer /1/ +input #=0, flalfvgy not active at iy=0 & ny;=1 active all iy
flalfipl  real /1.e20/ +input #ion therm flux lim factor on plates when isplflxl=0
flalfepl  real /1.e20/ +input #elec therm flux lim factor on plates when isplflxl=0
isfeexpl0 integer  /0/ +input #if=1, feex cannot be out of inner/outer plates
isfeixpl0 integer  /0/ +input #if=1, feix cannot be out of inner/outer plates
isofric   integer /0/  +input #If =1, use old (B2) interspecies up drag expression
del_te_ro real /1e10/  +input #te width in eV of tanh which turns off pwrze below 1 eV
iskaboom  integer /0/  #=1 turns on ijmgetmr "k" or "kaboom" stopping option
isnglf	  integer /0/  +input #=1 gives ng=nglfix at ix=0
nglfix    real  /1e15/ +input #value of ng at ix=0 if isnglf=1
isngrf	  integer /0/  +input #=1 gives ng=nglfix at ix=nx+1
ngrfix    real  /1e19/ +input #value of ng at ix=nx+1 if isnglf=1
isnion(nispmx) integer  /nispmx*1/   +input #user:turns on (=1) ion continuity eqn.
isupon(nispmx) integer  /nispmx*1/   +input #user:turns on (=1) parallel vel. eqn.
isup1up2       integer   /0/         +input #=1 sets up2=rup21*up1
rup21          real      /1./        +input #rup21=up2/up1 if isup1up2
isupgon(ngspmx) integer /ngspmx*0/   +input #user:=1 for par neutral vel. eqn.; index igsp
isteon      integer  /1/             +input #user:turns on (=1) electron energy eqn.
istion      integer  /1/             +input #user:turns on (=1) ion enegy eqn.
isngon(ngspmx) integer  /6*1/        +input #user:turns on (=1) neutral eqn.; index igsp
istgon(ngspmx) integer  /6*0/        +input #user:turns on (=1) gas enegy eqn.
isphion     integer  /0/             +input #user:turns on (=1) potential eqn.
isphiofft integer  /0/ +input   #user: set =1 only when continuing a simulation that 
                                #previously had isphion=1, but now the user wants to
                                #freeze the potential by setting isphion=0, while 
                                #continuing to use the poloidal electric field 
                                #ex=-d(phis)/dx and the previous parallel current

#user:=1 leaves old cur. on & ex=-d(phis)/dx; must be used
                       #with isphion=0
isnewpot  integer  /0/ +input #user:turns on (=1) new potential; J_r from tor. mom. bal.
                       #=-2 sets phi constant on core boundary with
                       #total core current = icoreelec
isnionxy(0:nx+1,0:ny+1,nisp)  _integer   #calc:=1 for ni eqn on; =0 for eqn off
isuponxy(0:nx+1,0:ny+1,nisp)  _integer   #calc:=1 for up eqn on; =0 for eqn off
isteonxy(0:nx+1,0:ny+1)	      _integer   #calc:=1 for te eqn on; =0 for eqn off
istionxy(0:nx+1,0:ny+1)	      _integer   #calc:=1 for ti eqn on; =0 for eqn off
isngonxy(0:nx+1,0:ny+1,ngsp)  _integer   #calc:=1 for ng eqn on; =0 for eqn off
istgonxy(0:nx+1,0:ny+1,ngsp)  _integer   #calc:=1 for tg eqn on; =0 for eqn off
isphionxy(0:nx+1,0:ny+1)      _integer   #calc:=1 for phi eqn on; =0 for eqn off
isnioffxy(0:nx+1,0:ny+1,nisp) _integer /0/ #user:=1, ni eqn off; =0 for eqn on
isupoffxy(0:nx+1,0:ny+1,nisp) _integer /0/ #user:=1, up eqn off; =0 for eqn on
isteoffxy(0:nx+1,0:ny+1)      _integer /0/ #user:=1, te eqn off; =0 for eqn on
istioffxy(0:nx+1,0:ny+1)      _integer /0/ #user:=1, ti eqn off; =0 for eqn on
isngoffxy(0:nx+1,0:ny+1,ngsp) _integer /0/ #user:=1, ng eqn off; =0 for eqn on
istgoffxy(0:nx+1,0:ny+1,ngsp) _integer /0/ #user:=1, tg eqn off; =0 for eqn on
isphioffxy(0:nx+1,0:ny+1)     _integer /0/ #user:=1, phi eqn off;=0 for eqn on
fdtnixy(0:nx+1,0:ny+1,nisp) _real /0./ #user:=1 for ni eqn off; =0 for eqn on
fdtupxy(0:nx+1,0:ny+1,nisp) _real /0./ #user:=1 for up eqn off; =0 for eqn on
fdttexy(0:nx+1,0:ny+1)      _real /0./ #user:=1 for te eqn off; =0 for eqn on
fdttixy(0:nx+1,0:ny+1)      _real /0./ #user:=1 for ti eqn off; =0 for eqn on
fdtngxy(0:nx+1,0:ny+1,ngsp) _real /0./ #user:=1 for ng eqn off; =0 for eqn on
fdttgxy(0:nx+1,0:ny+1,ngsp) _real /0./ #user:=1 for tg eqn off; =0 for eqn on
fdtphixy(0:nx+1,0:ny+1)     _real /0./ #user:=1 for phi eqn off; =0 for eqn on
isugfm1side               integer /0/   +input 
                                        #=0, use pol ave gas vels in par up eqn
                                        #=1, use 1-sided vals for domain decomp
isnupdot1sd               integer /0/   +input 
                                        #=0, use 2-pt ndot for (n*up)_dot;
                                        #=1, use 1-sided n_dot for (n*up)_dot
isphicore0		  integer /0/  +input #=1 sets phi=phi_mp in core if isphion=1
is_z0_imp_const           integer /0/  +input #=0 use hydr Keilhacker;=1 z0_imp_const
z0_imp_const              real    /1./ +input #z0 in therm force if is_z0_imp_const=1
					
***** Model_choice restart:
#Flags for choosing one or another calculation of a part of the model
iondenseqn	character*8	/"llnl"/	# ion continuity equation

***** Coefeq restart:
#Coefficients for the different terms in the fluid equations.
cnfx      real      /1./    +input #X-flux coef for conv. in n-eq.
cnfy      real      /1./    +input #Y-flux coef for conv. in n-eq.
cnsor     real      /1./    +input #Coef for particle src. in n-eq.
cfneut    real      /1./    +input #Coef for fluid neutrals contrib's to resid's
cfnidh    real      /1./    +input #Coef for neutral-ion drift heating
cfnidh2   real      /0./    +input #the above coef (cfnidh=1.0) is not exactly the real coef for neutral-ion drift heating term. That's why we introduce cfnidh2 but only for testing. Default =0.0: nothing; =1.0 (only for testing), remove the drift heating term.
cfnidhgy  real      /0./    +input # =1, consider vgy(,,1)**2 for n-i drift heating, assuming vy(,,0) negligible
cfnidhg2  real      /0./    +input # =1, consider vg2(,,1)**2 for n-i drift heating, assuming v2(,,0) negligible
cfnidhdis real      /0./    +input # =1, consider drift heating term (vg(,,1)-vg(,,2))^2 for molecular dissociation.
cfnidhmol real      /0./    +input # =1, consider vg2(,,2)**2 for n-i drift heating, e.g. H2,D2 dissociation.
cfupcx    real      /1./    +input #Coef for nucx*(up_ion - up_gas) momentum coupling
cfticx    real      /1./    +input #Coef for nucx*(up_ion-up_gas)**2 heating in Ti Eq
cfupimpg  real      /0./    +input #Coef for impur up Cx/elast drag on up=0 imp gas
cftiimpg  real      /0./    +input #Coef for Ti cooling CX/elast loss to cold imp gas
cmneut    real      /0./    +input #Coef for Monte Carlo neutrals contrib's to resid's
cnflux(ngspmx) real /ngspmx*1./ +input #coef for particle flux in n-eq. (resco)
chradi    real      /1./    +input #Coef for hyd. ioniz. rad. loss in elec. eng. eq.
chradr    real      /1./    +input #Coef for hyd. recomb. rad. loss in elec. eng. eq.
chioniz   real      /1./    +input #Coef for hydrogen ionization in elec. eng. eq.
cfizmol   real      /0./    #..Tom: Coef adding hyd ioniz rate to molec dissociation
                            #       rate to mimic ioniz of mols not in svdiss.
		            #       Tom added it for me, however surprised it is not
			    #       present in V8.0.0
ifxnsgi   integer   /0/	    +input #=1 sets ne for <sig*v>_i to cne_sgvi
cne_sgvi  real [1/m**3] /1.e18/ +input #ne for <sig*v>_i if ifxnsgi=1
ctsor     real      /1./    +input #Coef for eng. src. in Ti eq. 0.5*mi*up**2*psor
ceisor    real      /1./    +input #scale fac for ion energy source term (nu_i & eion)
ccoldsor  real      /0./    +input #scale fac for ion eng loss from cold cx
cngfx(ngspmx)  real /ngspmx*1./+input #scale fac for flux from grad_x T_g in gas eqn.
cngfy(ngspmx)  real /ngspmx*1./+input #scale fac for flux from grad_y T_g in gas eqn.
cngmom(nispmx) real /nispmx*0./+input #mom. cx-loss coeff for diffusve-neut hydr only
cmwall(nispmx) real /nispmx*0./+input #mom. wall-loss coeff for diff-neut hydr only
cngtgx(ngspmx) real /ngspmx*0./+input #X-flux coef for gas comp. of Ti eqn.
cngtgy(ngspmx) real /ngspmx*0./+input #Y-flux coef for gas comp. of Ti eqn.
sxgsol    real      /1./    +input #stretches x-coord. for gas in sol & core regions
sxgpr     real      /1./    +input #stretches x-coord. for gas in private flux region
xstscal   real  [m] /1./    +input #scale-length with stretch-coord decays from plates
rld2dxg(ngspmx) real /ngspmx*0./ +input #ratio of gas decay-length to dx via artificial diff.
rld2dyg(ngspmx) real /ngspmx*0./ +input #ratio of gas decay-length to dy via artificial diff.
cngflox(ngspmx) real /ngspmx*0./ +input #fac for x-flux from convection in ng-eqn.
cngfloy(ngspmx) real /ngspmx*1./ +input #fac for y-flux from convection in ng-eqn.
cngniflox(nispmx,ngspmx) real /nispmxngspmx*0./ +input #fac for rel ion-neut x-vel in ng-eqn
cngnifloy(nispmx,ngspmx) real /nispmxngspmx*0./ +input #fac for rel ion-neut y-vel in ng-eqn
cngsor          real    /1./     +input #Coef for part. src. in ng-eq.
cdifg(ngspmx)   real /ngspmx*1./ +input #scale factor for gas diffusion coeff.
lgmax(ngspmx)   real /ngspmx*1.e20/ +input #max gas scale length for calc particle D_g
lgtmax(ngspmx)  real /ngspmx*1.e20/ +input #max gas scale length for calc. thermal D_g
lgvmax          real   /1.e20/   +input #max gas scale length for calc. viscous D_g
rtg2ti(ngspmx)  real /ngspmx*1./ +input #ratio of gas temp to ion temp
tgas(ngspmx) real [eV] /ngspmx*1./ +input #value of tg if istgcon=1
cmfx      real      /1./    +input #X-flux coef for conv. in up-eq.
cmfy      real      /1./    +input #Y-flux coef for conv. in up-eq.
cpgx      real      /1./    +input #Coef for Grad(p) in up-eq.
cfvisx    real      /1./    +input #Coef. for x-visc. in ti-eq.
cfvisy    real      /1./    +input #Coef. for y-visc. in ti-eq.
cfanomvisxg real    /0./    +input #Coef. for neut x-visc ~travis(2)
cfanomvisyg real    /0./    +input #Coef. for neut y-visc ~travis(2)
cfvisxn   real      /1./    +input #Coef. for neutral x-visc. in up(,,iispg) eqn
cfvisxy(1:nispmx) real /nispmx*1./    +input #Coef. mult fmixy(ifld)
isvisxn_old integer /0/     +input #=1 uses sigcx,rrfac=1; =0 uses kelhihg, rrfac=rr**2
cfvxnrr   real      /1./    +input #=1 gives rr**2 in visx gas; =0 gives old 1 factor
cfvisyn   real      /1./    +input #Coef. for neutral y-visc. in up(,,iispg) eqn
cfvcsx(1:nispmx) real /nispmx*1./  +input #Coefs for x-visc. in ti-eq. with ismcnon>0
cfvcsy(1:nispmx) real /nispmx*1./  +input #Coefs for y-visc. in ti-eq. with ismcnon>0
isvhyha   integer   /0/     +input #switch (=1) for harmonic y-ave of up in visc heat
upvhflr   real      /1.e2/  +input #min denom for up harmc ave (isvhyha=1); visc heat
vboost    real      /1./    +input #previously scaled eqp; no longer in use
cvgp      real      /1./    +input #Coef for v.Grad(p) ion/elec eng. terms
cvgpg     real      /1./    +input #Coef for v.Grad(pg) gas eng. terms
cfvgpx(1:nispmx) real /nispmx*1./ +input #Coefs for x components of v.grad(p) in ti-eq
cftiexclg real      /1./    +input #Coef =1.0 for including atom gas contrib. in Ti eq. Make it 0.0 when turn on atom temperature.
cfvgpy(1:nispmx) real /nispmx*1./ +input #Coefs for y components of v.grad(p) in ti-eq
cfbgt     real      /0./    +input #Coef for the B x Grad(T) terms.
cfjhf     real      /1./    +input #Coef for convective cur (fqp) heat flow
jhswitch  integer   /0/     +input #Coef for the Joule-heating terms
oldseec   real      /0./    +input #Switch for Joule-heating bugfix
cf2ef     real      /0./    +input #Coef for ExB drift in 2-direction
cfyef     real      /0./    +input #Coef for ExB drift in y-direction
cftef     real      /0./    +input #Coef for ExB drift in toroidal direction
cf2bf     real      /0./    +input #Coef for Grad B drift in 2-direction
cfybf     real      /0./    +input #Coef for Grad B drift in y-direction
cfcbti    real      /0./    +input #Coef for adding fnixcb & fniycb to Ti eqn.
cfcurv	  real      /1./    +input #Coef for curvature part of Grad_B drift
cfgradb   real      /1./    +input #Coef for p_perp part of Grad_B drift	
cfq2bf    real      /0./    +input #Coef for Grad_B current in 2-direction
cfqybf    real      /0./    +input #Coef for Grad_B current in y-direction
cfqyn     real      /0./    +input #Coef for cx coll. rad current in y-direction
cfqym     real      /1./    +input #Coef for spatial inertial rad current in y-dir.
cfqydt    real      /0./    +input #Coef for time-dep inertial rad current in y-dir.
cf2dd     real      /0./    +input #Coef for diamagnetic drift in 2-direction
cfydd     real      /0./    +input #Coef for diamagnetic drift in y-direction
cftdd     real      /0./    +input #Coef for diamagnetic drift in toroidal direction
cfrd      real      /0./    +input #Coef for resistive cross-field drift
cfvisxneov real     /0./    +input #Coef for v-driven parallel viscosity
cfvisxneoq real     /0./    +input #Coef for q-driven parallel viscosity
cfvycr    real      /0./    +input #Coef for thermal force class. vel. vycr
cfvycf    real      /0./    +input #Coef for visc. force class. vel. vycf
cfvyavis  real      /0./    +input #Coef for vy from anom perp viscosity 
cfjve     real      /0./    +input #Coef for J-contribution to ve.
cfjp2     real      /0./    +input #Coef for B x gradP terms in div(J) eqn
cfjpy     real      /0./    +input #Coef for B x gradP terms in div(J) eqn
isnfmiy   integer   /0/     +input #diff fmiy for symmetry for vel. cells touching x-pt.
cfnfmiy   real      /1./    +input #Coef for new fmiy for vel. cells touching x-pt.
cnimp	  real	    /1./    +input #Coef for impurity radiation loss
fac2sp    real  [ ] /1./    +input #factor to test 2-species model; for
                            #equal densities, set fac2sp=2
$$$cfw(1:10) real   /0.1,9*1./ +input #Coeff's for the parallel neutral momentum eq.
$$$                            #cfw(1)*(the local sound speed) is the parallel
$$$                            #neutral velocity out of the plate (bound. con.)
cftnm     real  [ ] /0./    +input #Coef for neutral cx in toroidal mom. eq for fqy
cfupjr    real  [ ] /0./    +input #coef to include u_par in Jr calc.
cfcximp1  real  [ ] /1./    +input #coef multi. kcxrz for imp(+1)+D(0)->imp(0)+D(+1)
cfcximp2  real  [ ] /1./    +input #coef mult. kcxrz;imp(+p)+D(0)->imp(p-1)+D(+1),p>1
cfnetap   real  [ ] /1./    +input #coef mult. netap*fqp term in frice express.
fcdif     real  [ ] /1./    +input #coef mult all constant anomal diff coef
cfmsor    real  [ ] /1./    +input #coef mult msor and msorxr in up eqn.
cpiup(nispmx) real /nispmx*1./ +input #mult. press. grad term in up eqn
cfloyi    real  [ ] /2.5/   +input #coef mult ion radial convective energy flow
cfloye    real  [ ] /2.5/   +input #coef mult elec radial convective energy flow
cfcvte	  real  [ ] /1.0/   +input #coef mult elec poloidal convect(~5/2) energy flow
cfcvti	  real  [ ] /1.0/   +input #coef mult ion & neut pol convect(~5/2) energy flow
cfcvtg	  real  [ ] /1.0/   +input #coef mult gas pol convect(~5/2) energy flow
cfloxiplt real  [ ] /.0/    +input #coef mult neutral convect engy from plates
cfloygwall real [ ] /.0/    +input #coef mult neutral convect engy from walls
cftgdiss(ngspmx)  real  [ ] /ngspmx*1./ +input # coef mult tg*nu_diss eng loss
exjbdry   real  [ ] /10./   +input #exponent pwr to limit fqp < fqpsat at plates
cfgpijr   real  [ ] /1./    +input #scalar factor for grad_Pi term in fqya
iszeffcon integer   /0/     +input #if =1, zeff=zeffcon
zeffcon   real      /1./    +input #value of zeff if iszeffcon=1
alftng    real      /0./    +input #neutral thermal force coeff; careful of sign

cfqya     real      /1./    +input #Coef for anomalous current in y-direction (new model)
cfqyao    real      /1./    +input #Coef for anomalous current in y-direction (old model)
cfqyae    real      /1./    +input #Coef for anomalous electron current in y-direction
cfqyai    real      /0./    +input #Coef for anomalous ion current in y-direction
cfhcxgc(ngspmx) real /ngspmx*0./ +input # Coef constant pol heat conduct (chixg_use)
cfhcygc(ngspmx) real /ngspmx*0./ +input # Coef constant rad heat conduct (chiyg_use)
cftgcond  real      /1./    +input #Coef for gas thermal cond (usually molecules)
cftgeqp  real       /1.5/   +input #Coef for gas thermal equipartion (usually molecules)
 

***** Bcond restart:
#Variables for setting the boundary conditions.
ispfbcvsix integer      /0/ +input #=0 for uniform BC on PF wall;
                                   #=1 for user variable BC vs pol ix index
iswobcvsix integer      /0/ +input #=0 for uniform BC on outer wall;
                                   #=1 for user variable BC vs pol ix index
ibctepl   integer	/1/	+input 
                #Switch for ix=0 energy flux bc's
				#=0, fixed te (see tepltl)
				#=1, standard sheath transmission b.c.
				#=2, zero poloidal gradients for te
ibctipl	  integer       /1/     +input # Same as ibctepl, with te --> ti
ibctepr   integer	/1/	+input 
                #Switch for ix=nx+1 energy flux bc's
				#=0, fixed te (see tepltr)
				#=1, standard sheath transmission b.c.
				#=2, zero poloidal gradients for te
ibctipr   integer       /1/	+input # Same as ibctepr, with te --> ti
isphilbc	integer	/0/     +input 
                #Switch for ix=0 b.c. on phi
				#=0, phi = phi0l + kappal * te
				#=1, phi = phi0l
isphirbc	integer	/0/     +input 
                #Switch for ix=nx+1 b.c. on phi
				#=0, phi = phi0r + kappar * te
				#=1, phi = phi0r
iphibcc 	integer /3/	+input #core BC at iy=1 when isnewpot=1;iy=0 
                                #=1, d^2(ey)/dy^2=0
                                #=2, te=constant & ey(ixmp,0)=eycore
                                #=3, phi=constant & ey(ixmp,0)=eycore
                                #>3 or < 1 now unavailable, previously
				#dphi(ix,1)=dphi_iy1,isutcore ctrls ix=ixmp
iphibcwi integer /0/   +input #=0, d(ey)/dy=0
				#=1, phi(ix,0) = phintewi*te(ix,0)/ev
				#=3, d(phi)/dy/phi = 1/lyphi(1)
				#=4, phi(ix,0)=phiwi(ix) in PF region
iphibcwiix(0:nx+1) _integer +input # pol depend iphibcwi if ispfbcvsix=1
iphibcwo integer /0/    +input  #=0, d(ey)/dy=0
				#=1, phi(ix,ny+1) = phintewi*te(ix,ny+1)/ev
				#=3, d(phi)/dy/phi = 1/lyphi(2)
				#=4, phi(ix,ny+1)=phiwo(ix)
iphibcwoix(0:nx+1) _integer +input # pol depend iphibcwo if iswobcvsix=1
phiwi(0:nx+1) _real [eV] +input #/(nx+2)*0./
                                #PF wall phi profile if iphicwi=4;user set
phiwo(0:nx+1) _real [eV] +input #/(nx+2)*0./
                                #outer wall phi profile if iphicwo=4;user set
phintewi        real /3./	+input #phi/te on inner wall if iphibcwi=1
phintewo        real /3./	+input #phi/te on outer wall if iphibcwo=1
ncore(nispmx)  real [m**-3]/nispmx*7.e19/ +input #core ion dens if isnicore=1
upcore(nispmx) real [m/s]  /nispmx*0./    +input #core ion parall vel. if isupcore=0
ngcore(ngspmx) real [m**-3]/ngspmx*1.e15/ +input #core gas dens if isngcore=1
tgcore(ngspmx) real [eV]  /ngspmx*100./   +input #core gas temp if istgcore=1
isnicore(nispmx)  integer   /1,30*0/      +input #switch for ion-density core B.C.
				    #=1, set uniform, fixed density, ncore
				    #=0, set flux to curcore/sy locally in ix
				    #=2, set flux & ni over range
				    #=3, set icur=curcore-recycc*fngy, const ni
				    #=4, use impur. source terms (impur only)
				    #=5, set d(ni)/dy=-ni/lynicore at midp &
                                    #    ni constant poloidally
isfniycbozero(nispmx)   real /nispmx*0./ +input # Switch for divergence-free fluxes on core boundary
                    #=0, allows divergence-free fluxes to modify net core flux
                    #=1, redistributes fluxes due to divergence-free term 
                    #    without affecting the net core boundary flux
                    #=-1,assumes no divergence-free fluxes on the core boundary
isupcore(nispmx) integer /nispmx*0/ +input #=0 sets up=upcore on core bdry
				    #=1 sets d(up)/dy=0 on the core bdry
				    #=2 sets d^2(up)/dy^2 = 0
 				    #=3 sets fmiy = 0
			   	    #=4 sets tor. ang mom flux = lzflux & n*up/R=const
			   	    #=5 sets ave tor vel = utorave & n*up/R=const
isngcore(ngspmx) integer /ngspmx*0/ +input #switch for neutral-density core B.C.
				    #=0, set loc flux= -(1-albedoc)*ng*vtg/4
				    #=1, set uniform, fixed density, ngcore
				    #=2, not available
				    #=3, extrapolation, but limited
				    #=anything else, set zero deriv which was
				    #prev default inert hy
				    # anything else same as =0
istgcore(ngspmx) integer /ngspmx*1/ +input #switch for neutral-density core B.C.
                                    #=0, set tg(ixcore,0,igsp)=ti(ixcore,0)*cftgticore
				    #=1, set fixed temp tgcore(igsp)
				    #=2, set energy flux = cfalbedo(=2.0)*fng_alb*tg
				    #if > 2, set zero grad; tg(,0,)=tg(,1,)
cftgticore(ngspmx) real /ngspmx*1./ +input #set tg(ixcore,0,igsp)=ti(ixcore,0)*cftgticore(igsp) when istgcore(igsp) = 0.
curcore(1:nispmx) real [A] /0.,30*0./ +input #value of current from core if isnicore=0
lzcore(1:nispmx)  real [kg/ms] /nispmx*0./ +input #tor. ang. mom dens core bdry; phi eqn
lzflux(1:nispmx)  real [kg/s**2]/nispmx*0./ +input #tor. ang. mom dens flux core bdry; up eqn
utorave(1:nispmx) real [m/s] /nispmx*0./ +input #ave tor ave vel = utorave; up eqn
icoreelec real [A]   /0./    +input # electrical current from core
eycore	  real [V/m] /0./    +input # rad E-field core BC is iphibcc=3
cfniybbo  real       /0./    +input # factor to includ. vycb in fniy,feiy at iy=0 only
cfniydbo  real       /0./    +input # factor to includ. vycp in fniy,feiy at iy=0 only
cfeeybbo  real       /0./    +input # factor to includ. vycb in feey at iy=0 only
cfeeydbo  real       /0./    +input # factor to includ. vycp in feey at iy=0 only
cfeexdbo  real       /0./    +input # factor includ v2cde & BxgradTe in BC at ix=0,nx
cfeixdbo  real       /0./    +input # factor includ v2cdi & BxgradTi in BC at ix=0,nx
cfqybbo   real       /0./    +input # factor to includ. fqyb in core current B.C. only
cfqydbo   real       /0./    +input # factor to includ. fqyd in core current B.C. only
nfqya0core integer   /0/     +input # num iy core cells beyond iy=0 where force fqya=0
nfqya0pf  integer    /1/     +input # num. iy pf cells beyond iy=0 where force fqya=0
nfqya0ow  integer    /1/     +input # num iy outer wall cell below iy=ny+1 with fqya=0
ixfixnc   integer    /1/     +input # ix where ni=ncore if isnicore=2 begins
incixc	  integer    /0/     +input # ix range for ni=ncore from ixfixnc if isnicore=2
tcoree    real [eV] /50./    +input #core elecron temp if iflcore=0
tcorei    real [eV] /50./    +input #core ion temp if iflcore=0
tedge     real [eV] /2./     +input #edge ion and elec. temp; used for te,iwalli,o
                             #arrays if last element zero (as in interpolation)
tepltl    real [eV]  /2./    +input #left plate Te B.C. if ibctepl=0
tipltl    real [eV]  /2./    +input #left plate Ti B.C. if ibctipl=0
tepltr    real [eV]  /2./    +input #right plate Te B.C. if ibcteplr=0
tipltr    real [eV]  /2./    +input #right plate Ti B.C. if ibctiplr=0
tbmin     real [eV] /.1/     +input #min. wall & pf temp for extrap. b.c.(isextrt..)
nbmin     real [m**-3] /1.e17/ +input #min. wall & pf den for extrap. b.c.(isextrn..)
ngbmin    real [m**-3] /1.e10/ +input #min. core gas den for extrap. b.c.(isextrngc)
istewc integer /1/ +input # switch for outer-wall BC on Te
			     # =0, set zero energy flux
	 		     # =1, set fixed temp to tedge or tewallo
			     # =2, use extrapolation BC
			     # =3, set Te scale length to lyte
                             # =4, set feey = bceew*fniy*te
istewcix(0:nx+1) _integer +input # pol depend BC for istewc if iswobcvsix=1
istiwc integer /1/ +input #switch for outer-wall BC on Ti, see istewc detail
istiwcix(0:nx+1) _integer +input # pol depend BC for istiwc if iswobcvsix=1
istgwc(ngspmx) integer /ngspmx*0/ +input #switch for outer-wall BC on Tg(,0,igsp)
	 		     # =0, set fixed temp to tgwall
			     # =1, use extrapolation BC
			     # =2, set Tg scale length to lytg(2,
			     # =3, eng flux = 2Tg*Maxw-flux
			     # =4, eng flux = sum of all parts e.g. recycled,spttered,pumped, assuming a half-Maxw for each
			     # =5, tg = ti*cftgtiwc
           # >5, report error in input

istgwcix(0:nx+1,ngspmx) _integer +input #pol depend BC for chamber wall tg if iswobcvsix=1
istepfc integer  /0/ +input  # switch for priv.-flux BC on Te

			     # =0, set zero energy flux
	 		     # =1, set fixed temp to tedge or tewalli
			     # =2, use extrapolation BC
			     # =3, set Te scale length to lyte
           # =4, set feey = bceew*fniy*te
istepfcix(0:nx+1) _integer #/(nx+2)*0/ +input  #pol depend BC for pf wall te if ispfbcvsix=1 
istipfc integer  /0/ +input  #switch for priv.-flux To BC, see istepfc detail
istipfcix(0:nx+1) _integer +input  #pol depend BC for pf wall ti if ispfbcvsix=1
istgpfc(ngspmx) integer /ngspmx*0/   +input #switch for PF BC on Tg(,0,igsp)
	 		     # =0, set fixed temp to tgwall
			     # =1, use extrapolation BC
			     # =2, set Tg scale length to lytg(1,
			     # =3, eng flux = 2Tg*Maxw-flux
           # =4, eng flux = sum of all parts e.g. recycled,spttered,pumped, 
           # ...  assuming a half-Maxw for each
           # =5, tg = ti*cftgtipfc
           # >5, report error in input           
istgpfcix(0:nx+1,ngspmx) _integer +input # pol depend BC for pf wall tg if ispfbcvsix=1
cftgtiwc(ngspmx)   real   /ngspmx*1./    +input #wall Tg B.C.: tg = ti*cftgtiwc if istgwc=5
cftgtipfc(ngspmx)  real   /ngspmx*1./    +input #pfc Tg B.C.: tg = ti*cftgtipfc if istgpfc=5
tewalli(0:nx+1) _real [eV] +input #/(nx+2)*0./
                             #inner wall Te for istepfc=1.; = tedge if not set
tiwalli(0:nx+1) _real [eV] +input #/(nx+2)*0./
                             #inner wall Ti for istipfc=1.; = tedge if not set
tewallo(0:nx+1) _real [eV] +input #/(nx+2)*0./
                             #outer wall Te for istewc=1.; = tedge if not set
tiwallo(0:nx+1) _real [eV] +input #/(nx+2)*0./
                             #outer wall Ti for istiwc=1.; = tedge if not set
tgwall(ngspmx)   real [eV] /ngspmx*0.025/ +input #Wall gas temp BC
lyte(1:2)  real /2*1e20/ [m] +input #decaying rad Te grad leng;(1,2) istepfc,wc=3
lyteix(2,0:nx+1)    _real [m] +input # pol dep radial te grad length if set < 1e5
                             # istepfc,wc=3: 1:2=i:o, 2nd dim ix
lyti(1:2)  real /2*1e20/ [m] +input #decaying rad Ti grad leng;(1,2) istipfc,wc=3
lytg(1:2,ngspmx) real /12*1e20/ +input #rad tg scale length: PF (1,; Outer(2,
lytiix(2,0:nx+1)    _real [m] +input # pol dep radial ti grad length if set < 1e5
                             # istipfc,wc=3: 1:2=i:o, 2nd dim ix
lyphi(1:2) real /2*1e20/ [m] +input #decaying rad phi grad leng;(1,2) iphibcwi,o=3
lyphiix(2,0:nx+1)   _real [m] +input # pol dep radial phi grad length if set < 1e5
                             # isphipfc,wc=3: 1:2=i:o, 2nd dim ix
isextrnp  integer   /0/      +input #=1 sets extrap. b.c. at div. plate bound'y for ni
isextrnpf integer   /0/      +input #=1 sets extrap. b.c. at p.f. bound'y for ni
isextrtpf integer   /0/      +input #=1 sets extrap. b.c. at p.f. bound'y for Te & Ti
isextrngc integer   /0/      +input #=1 sets extrap. b.c. on core bdry for ng
isextrnw  integer   /0/      +input #=1 sets extrap. b.c. at outer wall for ni
isextrtw  integer   /0/      +input #=1 sets extrap. b.c. at outer wall for Te & Ti
iflcore   integer   /0/      +input #=0, core Te,i=tcoree,i; =1 core power=pcoree,i;
			     #=-1, core d(Te,i)/dy=0
pcoree    real [W]  /4e4/    +input #electron power from core if iflcore=1
pcorei    real [W]  /4e4/    +input #ion power from core if iflcore=1
ifluxni   integer  /1/       +input #flag for setting iy=0,ny+1 dens flux to 0 (=1,yes)
ckinfl    real    /1./       +input #includes kinetic viscosity in energy bound. cond.
isupss(nispmx) integer /nispmx*0/ +input
                             #=0, up=cs; =1, up>=1; =-1, dup/dx=0 at plates
isbohmms  integer /0/       +input  #=0 for single-species Bohm; =1 for multispecies B
isnwconi(1:nispmx) integer /nispmx*0/ +input 
		#switch for private-flux wall (iy=0) density B.C.
		#=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0
		#=1, fixed density to nwalli(ix) array
		#=2, extrapolation B.C.
		#=3, approx grad-length lyni, but limited by nwimin
isnwconiix(0:nx+1,1:nispmx) _integer +input  #pol depen isnwconi if ispfvsix=1
isnwcono(1:nispmx) integer /nispmx*0/ +input
		#switch for outer wall (iy=ny+1) density B.C.
		#=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0
		#=1, fixed density to nwallo(ix) array
		#=2, extrapolation B.C.
		#=3, approx grad-length lyni, but limited by nwomin
isnwconoix(0:nx+1,1:nispmx) _integer +input  #pol depen isnwcono if iswobcvsix=1
nwalli(0:nx+1) _real [m**-3] +input #inner wall dens set by isnwconi
nwallo(0:nx+1) _real [m**-3] +input #outer wall dens set by isnwcono
nwimin(nispmx)  real [m**-3] /nispmx*1e16/  +input # min inner wall dens if isnwconi=3
nwomin(nispmx)  real [m**-3] /nispmx*1e16/  +input # min outer wall dens if isnwcono=3
ncoremin(nispmx) real [m**-3] /nispmx*1e10/ +input # min ncore for isnicore=5
lyni(2)         real [m]     /2*1e20/       +input #rad dens grad length -isnwconi,o=3
lyniix(2,0:nx+1,nisp) _real [m] +input # pol dep radial dens grad length if set < 1e5
			         # isnwconi,o=3: 1:2=i:o, 2nd dim ix, 3rd spec
lynicore(nispmx) real [m] /nispmx*1e20/     +input # ni core BC rad scale-length if
					    # isnicore=5
lyup(2) real    /2*1e20/     +input #radial up grad length if isupwi,o=3: 1:2=i:o
isulyupx        integer      /0/            #if=0, lyupx filled with lyup
					    #if=1, user values of lynup used
lyupix(2,0:nx+1,nisp) _real [m] +input # pol dep radial up grad length if set < 1e5
			       # isupwi,o=3: indices,1:2=i:o, 2nd dim ix, 3rd spec
nwsor     integer    /1/     +input #number of sources on wall; must be < 10
igasi(10) real [Amp] /10*0./ +input #Gas currents from inner wall (iy=0)
igaso(10) real [Amp] /10*0./ +input #Gas currents from outer wall (iy=ny+1)
igspsori(10) integer /10*1/  +input #index of gas species for inner wall sources
igspsoro(10) integer /10*1/  +input #index of gas species for outer wall sources
issorlb(10) integer  /10*1/  +input
                             #flag for coord. origin of source.;=1, left plate;
                             #=0, right plate
jxsori(10)  integer  /10*1/  +input # xgasi=0. is located at left boundary
                             # of mesh region jxsori
jxsoro(10)  integer  /10*1/  +input # xgaso=0. is located at left boundary
                             # of mesh region jxsoro
xgasi(10) real [m]   /10*0./ +input
                             #location of inner wall sources; if issorlb(i)=1,0
                             #measured from left plate, right plate
xgaso(10) real [m]   /10*0./ +input 
                             #location of outer wall sources; if issorlb(i)=1,0,
                             #measured from left plate, right plate
wgasi(10)  real [m] /10*100./+input 
                             #total cosine widths of inner wall gas sources
wgaso(10)  real [m] /10*100./+input 
                             #total cosine widths of outer wall gas sources
albdsi(10) real [m] /10*1.0/ +input #albedos at inner gas source locations
albdso(10) real [m] /10*1.0/ +input #albedos at outer gas source locations
chemsputi(10,10)   /100*0/   +input #chem sputt coeff, priv flux surface, flux(i)=
                             # sum(chemsputi(i,j)*ng(j)*vt*sy)
chemsputo(10,10)   /100*0/   +input #chem sputt coeff, outer wall - see chemsputi def.
matwsi(10) integer  /10*0/   +input #material wall at inner gas source locations
matwso(10) integer  /10*0/   +input #material wall at outer gas source locations
issori(10) integer  /10*0/   +input #starting ix cell index for inner source
iesori(10) integer  /10*0/   +input #ending ix cell index for inner source
issoro(10) integer  /10*0/   +input #starting ix cell index for outer source
iesoro(10) integer  /10*0/   +input #ending ix cell index for outer source
iwalli(10) real     /10*0./  +input #current from inner source region isor for coupling
iwallo(10) real     /10*0./  +input #current from outer source region isor for coupling
ncpli(10)  integer  /10*0/   +input #flag for coupling between inner srce isor & ncpli
ncplo(10)  integer  /10*0/   +input #flag for coupling between outer srce isor & ncpli
cplsori(10)  real   /10*0./  +input #coeff. giving coupling from inner isor to ncpli
cplsoro(10)  real   /10*0./  +input #coeff. giving coupling from outer isor to ncpli
iscpli(0:nx+1)      _integer +maybeinput #(=1) => ix pt involved in inner bndry coupling
iscplo(0:nx+1)      _integer +maybeinput #(=1) => ix pt involved in outer bndry coupling
fwsori(0:nx+1,10)    _real   +maybeinput #profile of inner wall source isor (missing igasi)
fwsoro(0:nx+1,10)    _real   +maybeinput #profile of outer wall source isor (missing igasi)
fngysi(0:nx+1,ngsp)  _real   +maybeinput
                             #gas input flux from igasi on inner wall (calc)
fngyi_use(0:nx+1,ngsp) _real [1/m**3s] +input #user supplied gas input flux*area
fngysig(0:nxg+1,ngsp) _real  +maybeinput #global value of fngysi if domain decomp (parll)
fngyso(0:nx+1,ngsp)  _real   +maybeinput 
                             #gas input flux from igaso on outer wall (calc)
fngyo_use(0:nx+1,ngsp) _real [1/m**3s] +input #user supplied gas input flux*area
fngysog(0:nxg+1,ngsp) _real  +maybeinput #global value of fngyso if domain-decomp (parll)
albedoo(0:nx+1,ngsp) _real   +input
                             #albedo outer iy=ny+1 surface for neutrals (calc)
albedoog(0:nxg+1,ngsp) _real +maybeinput #global val albedoo if domain-decomp (parll)
albedoi(0:nx+1,ngsp) _real   +input 
                             #albedo of inner iy=0 surface for neutrals (calc)
albedoig(0:nxg+1,ngsp) _real +maybeinput #global val albedoi if domain-decomp (parll)
matwallo(0:nx+1)     _integer +maybeinput 
                             #flag (=1) denoting outer material side wall
matwallog(0:nxg+1)  _integer +maybeinput #global val matwallo if domain-decomp (parll)
matwalli(0:nx+1)    _integer +maybeinput 
                             #flag (=1) denoting inner material side wall (pf)
matwallig(0:nxg+1)  _integer +maybeinput #global val matwalli if domain-decomp (parll)
sinphi    real         /.1/     +input #sine of angle between side wall and flux surf.
isfixlb(nxptmx) integer /nxptmx*0/+input #=1 fixes values left bndry;=2 for symm. pt.
isfixrb(nxptmx) integer /nxptmx*0/+input #=2 for symmetry pt. at ix=nx+1
nib(20)   real [m**-3]/20*5.e19/+input #left plate density for isfixlb=1
upb(20)   real [m/s]  /20*1e2/  +input #left plate parallel velocity for isfixlb=1
teb       real [eV]    /50./    +input #left plate electron temp for isfixlb=1
tib       real [eV]    /50./    +input #left plate ion temp for isfixlb=1
yylb0     real [m]     /0./     +input #radial shift in LHB profiles for isfixlb=1
ywnii     real [m]     /.01/    +input #inner Gaussian radial width of nib
ywnio     real [m]     /.04/    +input #outer Gaussian radial width of nib
ywupi     real [m]     /.01/    +input #inner Gaussian radial width of upb
ywupo     real [m]     /.04/    +input #outer Gaussian radial width of upb
ywtei     real [m]     /.01/    +input #inner Gaussian radial width of teb
ywteo     real [m]     /.04/    +input #outer Gaussian radial width of teb
ywtii     real [m]     /.01/    +input #inner Gaussian radial width of tib
ywtio     real [m]     /.04/    +input #outer Gaussian radial width of tib
nibprof(0:ny+1) _real           +input #radial profile of nib
upbprof(0:ny+1) _real           +input #radial profile of upb
tebprof(0:ny+1) _real           +input #radial profile of teb
tibprof(0:ny+1) _real           +input #radial profile of tib
engbsr    real [ ]     /1./     +input #energy factor for backscattered neutrals to Ti
epsbs     real [ ]     /1.e-5/  +input #small fac added (substracted) from Rbs (Rfc)
isguardc  integer      /1/      +input #using guard cells? (=1 yes, =0 no)
rlimiter  real [m]     /1.e20/  +input #position of limiter at ix=0 for isfixlb=2
islimsor  integer      /0/      +input #=1 extends sources into limiter region
isutcore  integer      /0/      +input #Used for ix=ixcore phi BC ONLY IF iphibcc > 3
				#=0, tor mom=lzcore on core;
                #=1, d<uz>/dy=0;
				#>1, d^2(Ey)/dy^2=0 at outer midplane
isupwi(nispmx) integer /nispmx*2/ +input #options for up BC on inner (PF) wall
                        #=0 sets up=0 on inner wall
                        #=1 sets fmiy=0 (parallel mom-dens y-flux)
              			#=2 sets dup/dy=0 on inner wall
              			#=3 sets (1/up)dup/dy=1/lyup(1) scale length
isupwiix(0:nx+1,nispmx) _integer +input #pol depend BC for isupwi if ispfbcvsix=1
isupwo(nispmx) integer /nispmx*2/ +input #options for up BC on outer (main-chamber) wall
                        #=0 sets up=0 on outer wall
                        #=1 sets fmiy=0 (parallel mom-dens y-flux)
              			#=2 sets dup/dy=0 on outer wall
              			#=3 sets (1/up)dup/dy=1/lyup(2) scale length
isupwoix(0:nx+1,nispmx) _integer +input #pol depend BC for isupwo if iswobcvsix=1
islbcn    integer      /2/	+input 
                # b.c. for ni at limiter guard cells;
				# =0,1 set ni in 2 cells
				# =2 set ni in 1 cell, fnix at interface
islbcu    integer      /5/	+input 
                # b.c. for up at limiter guard cells;
				# =0,1 set up in 3 cells
				# =2 set up in 2 cells, fmix at interface
				# =3,4,6 set fmix at interface
				# =5 set fmix-fmixy at interface
islbce    integer      /2/	+input 
                # b.c. for te at limiter guard cells;
				# =0,1 set te in 2 cells
				# =2 set te in 1 cell, feex at interface
islbci    integer      /2/	+input 
                # b.c. for ti at limiter guard cells;
				# =0,1 set ti in 2 cells
				# =2 set ti in 1 cell, feix at interface
islbcg    integer      /2/	+input 
                # b.c. for ng at limiter guard cells;
				# =0,1 set ng in 2 cells
				# =2 set ng in 1 cell, fngx at interface
islbcp    integer      /2/	+input 
                # b.c. for phi at limiter guard cells;
				# =0,1 set phi in 2 cells
				# =2 set phi in 1 cell, fqx at interface
isph_sput(ngspmx) integer /ngspmx*0/  +input 
                #flag for plate sputtering;
				#0=old fixed case; 1=DIVIMP/JET phys sputt fits
				#=2 adds h-ion chem sputt;=3 adds h-neut c_sput
isi_sputw(ngspmx) integer /ngspmx*0/ +input 
                #flag for outer wall ion-based sputter;
				#=0, no ion sputtering
				#=1 adds phys ion sputt; =2 adds chem ion sputt
isi_sputpf(ngspmx) integer /ngspmx*0/ +input 
                #flag for priv flux ion-based sputter;
				#=0, no ion sputtering
				#=1 adds phys ion sputt; =2 adds chem ion sputt
islim_sput(ngspmx) integer /ngspmx*0/ +input
                #flag for limiter flux ion-based sputter;
				#=0, no ion sputtering
				#=1 adds phys ion sputt; =2 adds chem ion sputt

matt	  integer               #output flag from syld96 for sputt. target mat.
matp	  integer               #output flag from syld96 for sputt. plasma
cion      integer     /6/       +input #input to syld96; atom num. of sputt. target
cizb      integer     /1/       +input #input to syld96; max charge state of plasma
crmb	  real        /2./[AMU] +input #input to syld96; mass of plasma ions
isch_sput(ngspmx) integer /ngspmx*0/ +input 
                    #chem sputt. opt; 0=old;
			        #5=Roth,G-R; 
                    #6=Haasz97; 
                    #7=Haasz97+Davis at low E
eincid    real         [eV]     +input #incident energy of ion or neut. for chem sputt
t_wall    real       /300./ [K] +input #temp. of side wall; now use tvwallo,i
t_plat    real       /300./ [K] +input #temp. of divertor plate; now use tvplatlb,rb
tvwallo(0:nx+1) _real /300./ [K] +input #user outer wall temp if iswalltempc=0
tvwalli(0:nx+1) _real /300./ [K] +input #user inner wall temp if iswalltempc=0
tvplatlb(0:ny+1,nxptmx) _real /300./ [K] +input #user left plate temp if isplttempc=0
tvplatrb(0:ny+1,nxptmx) _real /300./ [K] +input #user left plate temp if isplttempc=0
flux_in   real        [1/m**2s] #incident ion or neutral flux for chem sputt
ychem     real                  #chem sputt. yield output from sputchem
yld_carbi(0:nx+1) _real         #chem sputt. yield, inner wall if isch_sput=5,6
yld_carbo(0:nx+1) _real         #chem sputt. yield, outer wall if isch_sput=5,6
fchemygwi(ngspmx) _real   /1./ +input #fac mult pf wall gas chem yield if isch_sput>0
fchemygwo(ngspmx) _real   /1./ +input #fac mult outer wall gas chem yield; isch_sput>0
fchemyiwi(ngspmx) _real   /1./ +input #fac mult pf wall ion chem yield if isch_sput>0
fchemyiwo(ngspmx) _real   /1./ +input #fac mult outer wall ion chem yield; isch_sput>0
fphysyiwi(ngspmx) _real   /1./ +input #fac mult pf wall ion phys yield if isch_sput>0
fphysyiwo(ngspmx) _real   /1./ +input #fac mult outer wall ion phys yield; isch_sput>0
fchemylb(ngspmx,nxptmx) _real /1./ +input #fac*inner plt gas chem yield; isch_sput>0
fchemyrb(ngspmx,nxptmx) _real /1./ +input #fac*outer plt gas chem yield; isch_sput>0
fphysylb(ngspmx,nxptmx) _real /1./ +input #fac*inner plt ion phys sp yield;isch_sput>0
fphysyrb(ngspmx,nxptmx) _real /1./ +input #fac*outer plt ion phys sp yield;isch_sput>0
fchemyllim(ngspmx)      _real /1./ +input #fac*left-limiter gas chem yield; isch_sput>0
fchemyrlim(ngspmx)      _real /1./ +input #fac*right-limiter gas chem yield; isch_sput>0
fphysyllim(ngspmx)      _real /1./ +input #fac*left limiter ion phys sp yield;isch_sput>0
fphysyrlim(ngspmx)      _real /1./ +input #fac*right limiter plt ion phys sp yield;isch_sput>0
isexunif  integer      /0/      +maybeinput #=1 forces ex ~ uniform at div. plates
xcnearlb  logical    /FALSE/    #=TRUE if Jac'n "box" overlaps a left boundary
xcnearrb  logical    /FALSE/    #=TRUE if Jac'n "box" overlaps a right boundary
openbox   logical    /FALSE/    #=TRUE if Jac'n "box" is wide open
kappa0	  real       /3.0/	+maybeinput #modified sheath drop (allows j>jsat) for kappa > kappa0
kappamx   real       /10.0/	+maybeinput #maximum kappa value
fqpsatlb(0:ny+1,nxptmx)	_real	#ion saturation current at left boundary
fqpsatrb(0:ny+1,nxptmx)	_real	#ion saturation current at right boundary
cfueb		real	/1./	+input #scale factor for ueb in plate b.c.'s
ikapmod   integer	/0/	    +input #=1 for new kappa model; =0 for qpfac model
fvapi(10) real		/10*0./	+input #scale factor for inner evap vapor source
avapi(10) real		/10*1./	+input #linear coeff. for inner evap vapor source
bvapi(10) real		/10*1./ +input #exponent coeff. for inner evap vapor source
fvapo(10) real		/10*0./	+input #scale factor for outer evap vapor source
avapo(10) real		/10*1./	+input #linear coeff. for outer evap vapor source
bvapo(10) real		/10*1./ +input #exponent coeff. for outer evap vapor source
tvapi(0:nx+1) _real [K]   	+input #inner wall temp for evap; input after alloc
tvapo(0:nx+1) _real [K]   	+input #outer wall temp for evap; input after alloc
cfvytanbc real          /1./    +input #factor for adding vytan to plate B.C.
totfeexl(0:ny+1,nxpt) _real [W] +maybeinput #elec polod energy flux*area on "left" plate
totfeexr(0:ny+1,nxpt) _real [W] +maybeinput #elec polod energy flux*area on "right" plate
totfeixl(0:ny+1,nxpt) _real [W] +maybeinput #elec polod energy flux*area on "left" plate
totfeixr(0:ny+1,nxpt) _real [W] +maybeinput #elec polod energy flux*area on "right" plate
cgpl      real           /0./   +input #scale fac atom eng plate loss; experim.
cgpld     real           /0./   +input #scale fac disso eng loss; experim.
cgengpl   real		 /0./       +input #new scale fac atom eng plate loss; old cgpl
cgengw    real           /0./   +input #new scale fac atom eng wall loss
cgmompl   real           /1./   +input #scale fac atom par mom plate loss
vgmomp    real     [m/s] /2.e3/ +input #vel used in exp factor of atom mom loss
istglb(ngspmx) _integer  /0/    +input #=0 for tg=tgwall; 
                                       #=1 for extrap;
			               #=3, Maxw flux;
			               #=4 for assuming half-Maxwellian for all kinds of neutral sources e.g. recycled, sputtered, pumped etc.;
				       #=5 for tg = ti*cftgtipltl.
istgrb(ngspmx) _integer  /0/    +input #=0 for tg=tgwall;
                                       #=1 for extrap;
			               #=3, Maxw flux;
			               #=4 for assuming half-Maxwellian for all kinds of neutral sources e.g. recycled, sputtered, pumped etc.;
				       #=5 for tg = ti*cftgtipltr.
cftgtipltl(ngspmx)  real   /ngspmx*1./    +input #left plate Tg B.C.: tg = cftgtipltl*ti if istglb=5
cftgtipltr(ngspmx)  real   /ngspmx*1./    +input #right plate Tg B.C.: tg = cftgtipltr*ti if istgrb=5
cgengmpl  real		 /1./   +input #scale fac mol plate eng loss for Maxw
cgengmw  real		 /1./   +input #scale fac mol wall eng loss for Maxw

***** Outpwall:
# Arrays used to communicate wall fluxes to a wall simulation code
ue_part_fluxelb(0:ny+1,nxptmx) _real [1/m**2s] #inner plt elec flux
ue_part_fluxerb(0:ny+1,nxptmx) _real [1/m**2s] #outer plt elec flux
ue_part_fluxeyi(0:nx+1)        _real [1/m**2s] #inner (PF) wall elec part flux
ue_part_fluxeyo(0:nx+1)        _real [1/m**2s] #outer (PF) wall elec part flux
ue_part_fluxh2p1lb(0:ny+1,nxptmx) _real [1/m**2s] #inner plt deut ion part flux
ue_part_fluxh2p1rb(0:ny+1,nxptmx) _real [1/m**2s] #outer plt deut ion part flux
ue_part_fluxh2p1yi(0:nx+1)        _real [1/m**2s] #inner (PF) wall deut ion part flux
ue_part_fluxh2p1yo(0:nx+1)        _real [1/m**2s] #outer (PF) wall deut ion part flux
ue_part_fluxh2lb(0:ny+1,nxptmx) _real [1/m**2s] #inner plt deut neut part flux
ue_part_fluxh2rb(0:ny+1,nxptmx) _real [1/m**2s] #outer plt deut neut part flux
ue_part_fluxh2yi(0:nx+1)        _real [1/m**2s] #inner (PF) wall deut neut part flux
ue_part_fluxh2yo(0:nx+1)        _real [1/m**2s] #outer (PF) wall deut neut part flux
ue_heat_fluxelb(0:ny+1,nxptmx) _real [W/m**2] #inner plt elec heat flux
ue_heat_fluxerb(0:ny+1,nxptmx) _real [W/m**2] #outer plt elec heat flux
ue_heat_fluxeyi(0:nx+1)        _real [W/m**2] #inner (PF) wall elec heat flux
ue_heat_fluxeyo(0:nx+1)        _real [W/m**2] #outer (PF) wall elec heat flux
ue_heat_fluxh2p1lb(0:ny+1,nxptmx) _real [1/m**2s] #inner plt deut ion heat flux
ue_heat_fluxh2p1rb(0:ny+1,nxptmx) _real [1/m**2s] #outer plt deut ion heat flux
ue_heat_fluxh2p1yi(0:nx+1)        _real [1/m**2s] #inner (PF) wall deut ion heat flux
ue_heat_fluxh2p1yo(0:nx+1)        _real [1/m**2s] #outer (PF) wall deut ion heat flux
ue_heat_fluxh2lb(0:ny+1,nxptmx) _real [1/m**2s] #inner plt deut neut heat flux
ue_heat_fluxh2rb(0:ny+1,nxptmx) _real [1/m**2s] #outer plt deut neut heat flux
ue_heat_fluxh2yi(0:nx+1)        _real [1/m**2s] #inner (PF) wall deut neut heat flux
ue_heat_fluxh2yo(0:nx+1)        _real [1/m**2s] #outer (PF) wall deut neut heat flux
ue_mean_engelb(0:ny+1,nxptmx) _real [J] #inner plt elec mean energy
ue_mean_engerb(0:ny+1,nxptmx) _real [J] #outer plt elec mean energy
ue_mean_engeyi(0:nx+1)        _real [J] #inner (PF) wall elec mean energy
ue_mean_engeyo(0:nx+1)        _real [J] #outer (PF) wall elec mean energy
ue_mean_engh2p1lb(0:ny+1,nxptmx) _real [J] #inner plt deut ion mean energy
ue_mean_engh2p1rb(0:ny+1,nxptmx) _real [J] #outer plt deut ion mean energy
ue_mean_engh2p1yi(0:nx+1)        _real [J] #inner (PF) wall deut ion mean energy
ue_mean_engh2p1yo(0:nx+1)        _real [J] #outer (PF) wall deut ion mean energy
ue_mean_engh2lb(0:ny+1,nxptmx) _real [J] #inner plt deut neut mean energy
ue_mean_engh2rb(0:ny+1,nxptmx) _real [J] #outer plt deut neut mean energy
ue_mean_engh2yi(0:nx+1)        _real [J] #inner (PF) wall deut neut mean energy
ue_mean_engh2yo(0:nx+1)        _real [J] #outer (PF) wall deut neut mean energy
ue_pot_engh2p1lb(0:ny+1,nxptmx)  _real [J] #inner plt deut ion pot energy
ue_pot_engh2p1rb(0:ny+1,nxptmx)  _real [J] #outer plt deut ion pot energy
ue_pot_engh2p1yi(0:nx+1)         _real [J] #inner (PF) wall deut ion pot energy
ue_pot_engh2p1yo(0:nx+1)         _real [J] #outer (PF) wall deut ion pot energy
  
***** Rccoef:
#Variables for recycling coeff. profiles on divertor plates
#Set for ngspmx gas species
recylb(0:ny+1,ngspmx,nxptmx) _real     +input
                       #tot inner plate recycling coeff. (calc)
				       #if recylb > 0, recycling coeff
				       #if in range [-1,0], acts as albedo
				       #if in range (-2,-1), gives ng=nglfix
				       #if recylb <= -2, gives ng(1)=ng(0)
recyrb(0:ny+1,ngspmx,nxptmx) _real     +input 
                       #tot outer plate recycling coeff. (calc)
				       #if recyrb > 0, recycling coeff
				       #if in range [-1,0], acts as albedo
				       #if in range (-2,-1), gives ng=ngrfix
				       #if recyrb <= -2, gives ng(nx+1)=ng(nx)
recylb_use(0:ny+1,ngspmx,nxptmx) _real +input #inner plate recycling coeff. user input
recyrb_use(0:ny+1,ngspmx,nxptmx) _real +input #outer plate recycling coeff. user input
recycp(ngspmx)   real       /.9,5*0./  +input #recycling coef at plates if ndatlb,rb=0
recyllim(0:ny+1,ngspmx) _real /1./     +input #recyling coeff on left-side of limiter
recyrlim(0:ny+1,ngspmx) _real /1./     +input #recyling coeff on right-side of limiter		 
recycflb(ngspmx,nxptmx) _real /1./     +maybeinput #extra factor for recycling at ix=0
recycfrb(ngspmx,nxptmx) _real /1./     +maybeinput #extra factor for recycling at ix=nx+1
recycm           real       /-0.9/     +input #mom recycling inertial gas; at plates
                                       # =up(,,2)/up(,,1)=-recycm.
				                               #if recycm betwn -9.9 & -10.1 d(up)/dx=0
                                       #if recycm < -10.1, therm mom flux used
recycmlb_use(0:ny+1,ngspmx,nxptmx) _real +maybeinput #inner plt mom-recycl coeff user input
recycmrb_use(0:ny+1,ngspmx,nxptmx) _real +maybeinput #outer plt mom-recycl coeff user input
recycmlb(0:ny+1,ngspmx,nxptmx) _real   +maybeinput #total inner plt mom recycling coeff
recycmrb(0:ny+1,ngspmx,nxptmx) _real   +maybeinput #total outer plt mom recycling coeff
recyce           real       /0./       +input #energy recycling/Rp for inertial gas
recycwe          real       /0./       +input #energy recycling/Rp for inertial gas on walls and prfs
recycml          real       /0.1/      +input #momentum recycling/Rp for gas at limtr
recycc(ngspmx)   real       /6*1./     +input #core recycling coeff. if isnicore=3
albedoc(ngspmx)  real       /6*1./     +input #core neut albedo for isngcore=0
albedolb(ngspmx,nxptmx) _real /1./     +input #albedo at inner plate if ndatlb=0
albedorb(ngspmx,nxptmx) _real /1./     +input #albedo at outer plate if ndatrb=0
cfalbedo         real       /2./       +input #coef. for Tg Eq. due to albedo
ndatlb(ngspmx,nxptmx)    _integer   /0/  +maybeinput #number of recycp data pts on inner plt
ndatrb(ngspmx,nxptmx)    _integer   /0/  +maybeinput #number of recycp data pts on outer plt
ydatlb(ngspmx,50,nxptmx) _real [m]  /0./ +maybeinput #inner data pt location from sep.
ydatrb(ngspmx,50,nxptmx) _real [m]  /0./ +maybeinput #outer data pt location from sep.
rdatlb(ngspmx,50,nxptmx) _real      /0./ +maybeinput #inner recycp data for each ydatlb
rdatrb(ngspmx,50,nxptmx) _real      /0./ +maybeinput #outer recycp data for each ydatrb
alblb(0:ny+1,ngspmx,nxptmx) _real        +input #inner plate albedo; used if <1 (calc)
albrb(0:ny+1,ngspmx,nxptmx) _real        +input #outer plate albedo; used if <1 (calc)
areapl          real   /0./ +work   # Work variable for plate projection area for 
                                    # albedo-like recycling
isoldalbarea    real   /0./ +input # Switch whether to use old (wrong) albedo
                                    # area which is perpendicular-to-poloidal flux
                                    # tube area (=1) or the correct area projected
                                    # onto the target plate (=0)
albedo_by_user   integer            /0/  +input #if=1, user fills albedoo,i & albdlb,rb
fngxslb(0:ny+1,ngspmx,nxptmx) _real [1/s]+input #inner plt liq vapor gas sour. if sputtlb>0
fngxsrb(0:ny+1,ngspmx,nxptmx) _real [1/s]+input #outer plt liq vapor gas sour. if sputtlb>0
fngxsllim(0:ny+1,ngspmx)      _real [1/s]+input #left limiter liq vapor gas sour. if sputtlb>0
fngxsrlim(0:ny+1,ngspmx)      _real [1/s]+input #right limiter liq vapor gas sour. if sputtlb>0
fngxlb_use(0:ny+1,ngspmx,nxptmx) _real [1/s] +input #user external left plate source
fngxrb_use(0:ny+1,ngspmx,nxptmx) _real [1/s] +input #user external left plate source
fngxllim_use(0:ny+1,ngspmx)      _real [1/s] +input #user external left limiter source
fngxrlim_use(0:ny+1,ngspmx)      _real [1/s] +input #user external right limiter source
adatlb(ngspmx,50,nxptmx) _real    /1./   +maybeinput #inner albdedo data for each ydati
adatrb(ngspmx,50,nxptmx) _real    /1./   +maybeinput #outer albdedo data for each ydati
recycw(ngspmx)   real     /ngspmx*1e-10/ +input #recycling coef. at side walls
recypf_use(0:nx+1,ngspmx,nxptmx) _real   +maybeinput #priv flux recycling coef; user input
recywall_use(0:nx+1,ngspmx)      _real   +maybeinput #outer wall recycling coef; user input
recycwit(0:nx+1,ngspmx,nxptmx)   _real   +maybeinput  #tot recyc coeff on PF wall
recycwot(0:nx+1,ngspmx)         _real    +maybeinput #tot recyc coeff on outer wall
isrefluxclip     integer       /1/       +input #=1 prohib outward gas for inward ion
gamsec           real       /0./         +input #secondary elec emiss coeff on plates
sputtr           real       /0./         +input #sputtering coef. at plates
sputtlb(0:ny+1,ngspmx,nxptmx)    _real   +input #set sputt coef. inner plate (iy,igsp)
sputtrb(0:ny+1,ngspmx,nxptmx)    _real   +input #set sputt coef. outer plate (iy,igsp)
sputllim_use(0:ny+1,ngspmx)      _real   +input #user-set sputt coef. left 
					        #limiter if>1; albedo if <0, >-1;
                                                #ng=ngllim on limiter if < -1 
sputrlim_use(0:ny+1,ngspmx)      _real   +input #user-set sputt coef. right 
					        #limiter if>1; albedo if <0, >-1;
                                                #ng=ngrlim on limiter if < -1		       	  
sputflxlb(0:ny+1,ngspmx,nxptmx)  _real   +maybeinput #calc sput flux inner plate (iy,igsp)
sputflxrb(0:ny+1,ngspmx,nxptmx)  _real   +maybeinput #calc sput flux outer plate (iy,igsp)
sputflxllim(0:ny+1,ngspmx)       _real   #calc sput flux inner plate (iy,igsp)
sputflxrlim(0:ny+1,ngspmx)       _real   #calc sput flux outer plate (iy,igsp)
sputflxw(0:nx+1,ngspmx)          _real   +maybeinput #calc sput flux outer wall (ix,igsp)
sputflxpf(0:nx+1,ngspmx)         _real   +maybeinput #calc sput flux PF wall (ix,igsp)
ngplatlb(ngspmx,nxptmx)          _real   +input #ng on inner plate if sputti < -1
ngplatrb(ngspmx,nxptmx)          _real   +input #ng on outer plate if sputto < -1
ngllim(ngspmx)                   _real   +input #ng on left limiter if sputllim_use < -1
ngrlim(ngspmx)                   _real   +input #ng on right limiter if sputrlim_use < -1
ipsputt_s		integer /1/      +input #start dens-index phys sputt species
ipsputt_e		integer /1/      +input #end dens-index of phys sputt species
npltsor     integer    /1/               +input #number sources on plates; must be <= 10
igaslb(10,nxptmx) _real [Amp] /0./     +input #Gas cur from left-hand plate(s) (ix=0)
igasrb(10,nxptmx) _real [Amp] /0./     +input #Gas cur from right-hand plate(s) (ix=nx)
igspsorlb(10,nxptmx) _integer /1/      +input #gas species index, left-hand plate sources
igspsorrb(10,nxptmx) _integer /1/      +input #gas species index, right-hand plate sources
ygaslb(10,nxptmx)  _real [m]  /0./     +input #loc of left-plate sources wrt strike pt
ygasrb(10,nxptmx)  _real [m]  /0./     +input #loc of right-plate sources wrt strike pt
wgaslb(10,nxptmx)  _real [m] /100./    +input #total cos width of left-plate gas sources
wgasrb(10,nxptmx)  _real [m] /100./    +input #total cos  width of right-plate gas sources
fvaplb(ngspmx,nxptmx) _real     /0./   +input #scale factor left-plate evap vapor source
avaplb(ngspmx,nxptmx) _real [k**.5/(m**2s)] /1./ +input #lin coeff left-plate evapor sor
bvaplb(ngspmx,nxptmx) _real [K] /1./   +input #expon. coeff. left-plate evap vapor source
fvaprb(ngspmx,nxptmx) _real     /0./   +input #scale factor right-plate evap vapor source
avaprb(ngspmx,nxptmx) _real [k**.5/(m**2s)] /1./ +input #lin coeff right-plate evapor sor
bvaprb(ngspmx,nxptmx) _real [K] /1./   +input #expon coeff. right-plate evap vapor source
tvaplb(0:ny+1,nxptmx)   _real [K]      +input #left-plate temp for evap; input after alloc
tvaprb(0:ny+1,nxptmx)   _real [K]      +input #right-plate temp for evap; input after alloc
fvapllim(ngspmx)        _real     /0./ +input #scale factor left-plate evap vapor source
avapllim(ngspmx) _real [k**.5/(m**2s)] /1./ +input #lin coeff left-plate evapor sor
bvapllim(ngspmx) _real [K] /1./        +input #expon. coeff. left-plate evap vapor source
fvaprlim(ngspmx) _real     /0./        +input #scale factor right-plate evap vapor source
avaprlim(ngspmx) _real [k**.5/(m**2s)] /1./ +input #lin coeff right-plate evapor sor
bvaprlim(ngspmx) _real [K] /1./        +input #expon coeff. right-plate evap vapor source
tvapllim(0:ny+1) _real [K]             +input #left-lim-face temp for evap; input after alloc
tvaprlim(0:ny+1) _real [K]             +input #right-lim-fac temp for evap; input after alloc
tvliml(0:ny+1)   _real /300./ [K]      +input #user left limiter face temp
tvlimr(0:ny+1)   _real /300./ [K]     +input #user right limiter face temp
isextpltmod            integer  /0/   +input #=1 use ext gas plate fluxes fngxextlb,rb
                                      # and feixextlb,rb
isextwallmod           integer  /0/   +input #=1 use ext gas wall fluxes fngyexti,o
                                      # and feiyexti,o
isoutwall              integer  /0/   +input #=1 call outwallflux to export wall fluxes
fngxextlb(0:ny+1,ngspmx,nxptmx) _real [1/s]+maybeinput #inner plt external particle flux*A
fngxextrb(0:ny+1,ngspmx,nxptmx) _real [1/s]+maybeinput #outer plt external particle flux*A
fngyexti(0:nx+1,ngspmx)         _real [1/s]+maybeinput #inner wall external particle flux*A
fngyexto(0:nx+1,ngspmx)         _real [1/s]+maybeinput #outer wall external particle flux*A
feixextlb(0:ny+1,nxptmx)        _real [J/s]+maybeinput #inner plt external energy flux*A
feixextrb(0:ny+1,nxptmx)        _real [J/s]+maybeinput #outer plt external energy flux*A
feiyexti(0:nx+1)                _real [J/s]+maybeinput #inner wall external energy flux*A
feiyexto(0:nx+1)                _real [J/s]+maybeinput #outer wall external energy flux*A

***** Fixsrc:
#Variables for including a fix-source at (xxsrc,yysrc) with int.=a*
ifixsrc     integer         /0/        #=1 turns on fixed Gaussian source
ifixpsor    integer         /0/        #=1 freezes part. source to initial val.
xxsrc       real            /3./
yysrc       real            /0.2/
c1n         real            /0./
c1e         real            /0./
c1i         real            /0./
a1n         real            /.5e+23/
a1e         real            /150.e+23/
a1i         real            /45.e+23/
b1n         real            /.1111/
b1e         real            /.1111/
b1i         real            /.1111/

***** Selec:
#Variables for the calculation of the Jacobian locally.
i1          integer
i2          integer
i2p         integer    #used for 4th-order diffusion in x
i3          integer
i4          integer
i5          integer
i5m         integer    #same as i5, except restricted to ix<nx
i6          integer
i7          integer
i8          integer
j1          integer
j2          integer
j3          integer
j4          integer
j5          integer
j5m         integer	#same as j5, except restricted to iy<ny
j6          integer
j7          integer
j8          integer
j1p         integer	#y-index lower range for potential eqn
j2p         integer	#y-index lower range for potential eqn
j5p         integer	#y-index upper range for potential eqn
j6p         integer	#y-index upper range for potential eqn
ixs1        integer
ixf6        integer
iys1        integer
iyf6        integer
xlinc       integer         /2/
xrinc       integer         /1/
yinc        integer         /2/
isjaccorall integer         /1/ #if=1 uses all ix cells for iy=0 Jac
ixm1(0:nx+1,0:ny+1)       _integer
ixp1(0:nx+1,0:ny+1)       _integer
iym1a(0:nx+1,0:ny+1)      _integer  # for mdsplus use only
iyp1a(0:nx+1,0:ny+1)      _integer
stretcx(0:nx+1,0:ny+1)    _real      #array for stretching gas x-coord.

***** Aux:
#Variables for auxiliary storage.
ixmp 	    integer		# poloidal index of outer midplane; for yyc,f


***** Err_msg_out restart:
#Controls for output of error messages
errmsgflag	integer	/1/	# =0 turns off error messages, =1 turns them on
errunit		integer	/6/	# output unit for error messages
				# (nksol ignores errunit, sending everything to 6)

***** Opt_input:
#provides optional input for vodpk and daspk through non-eraseable arrays
inopt          integer   /0/   #resets iopts for solvers (vodpk, daspk)
iworkin(1:25)  integer /25*0/  #optional integer array input for vodpk
rworkin(1:25)  real    /25*0./ #optional real array input for vodpk

***** Lsode:
#Parameters and variables for the solver routines (originally lsodpk)
runtim       real /1.e-7/     #time of first output; total time=runtim*trange
trange       real /1.e7/      #factor multiplying runtim to give total sim. time
neq          integer          #total number of equations over whole domain
jacflg       integer /1/      #flag for computing Jacobian in vodpk
jpre         integer /1/      #flag for using the preconditioning step in vodpk
rtolv(30)    real  /30*1.e-4/ #relative tol. vector used in convert.m
itol         integer /4/
itask        integer /1/
istate       integer /1/
iopts        integer /0/   #internally set to inopt, an input variable
mf           integer /24/  #vodpk flag: mf=21, full user J; mf=22,full lsode J
                           #mf=24, banded user J; mf=25, banded lsode J
idid         integer /0/
ires         integer /0/
ts           real    /0./  #start time for ODE solvers
tout         real	   #output times for ODE solvers
dtmax        real /1.e20/  #maximum allowed dt for daspk if info(7)=1
dtinit       real /1.e-10/ #starting dt for daspk if info(8)=1
maxpoly      integer /5/   #maximum polynomial power used in daspk timestepping
ipar(3)      integer  #integer parameter communication
rpar(1)      real     #real parameter communication
rtol(neqmx)  _real    #relative tol. used by solvers; rtol=rtolv(igrid)*tolbf
atol(neqmx)  _real    #absolute tolerance for lsode-like routines
yl(neqmx)    _real    #primary variables for ODE's
yldot(neqmx) _real    #time-derivatives of yl's; RHS of ODE's
delta(neqmx) _real
iextra(5)    integer
info(20)     integer /0,1,0,0,1,1,0,0,0,1, 1,1,0,0,1,0,0,1,0,0/ #set for daspk
		# info(1) =0, indicates start of new problem to initializd code
		# info(2) =0, rtol,atol scalar; = 1, rtol,atol vectors
		# info(3) =0, output only at tout; =1, output at intermed. t
		# info(4) =0, no restrict on t; =1, stop at t=tstop=rwork(1)
		# info(5) =0, daspk gen Jac, only if info(12)=0; =1, ext Jac
		# info(6) =0, full mat sol, only if info(12)=0; =1, band sol
		#             iwork(1)=lbw, iwork(2)=ubw
		# info(7) =0, code sets dt_max; =1, rwork(2)=dt_max
		# info(8) =0, code sets dt_init; =1, rwork(3)=dt_init
		# info(9) =0, maxord=5; =1 iwork(3)=maxord (=<5)
		# info(10)=0, no constraints; =1, constr initial cond.;
		#             =2, constr Y>0; =3, options 1 & 2 added
		# info(11)=0, consist init cond; =1, calc init cond(alg + der)
		#             =2, calc init cond (use Y', calc Y)
		# info(12)=0, direct mat sol; =1, Krylov method
		# info(13)=0, default Krylov param; =1, iwork(24)=maxl,
		#             iwork(25)=kmp, iwork(26)=nrmax; rwork(10)=eplidpk
                #             default Krylov and direct pram: rwork(16)=epnldpk
		# info(14)=0, proc after init; =1, stop after init, then
		#             reset info(1)=0 to avoid recalc init
		#             (or info(14)=0 & info(11)=0 to avoid recalc init?)
		# info(15)=0, no Jac routine; =1 Jac provided
		# info(16)=0, err chk on all Y; =1, no err on alg Y
		# info(17)=0, init calc control default; =1,input init control
		# info(18)=0. no init print; =1, min print; =2, full print
		# info(19)    not specified
		# info(20)    not specified

eplidpk      real /5.e-2/  # optional input for daspk when info(13) = 1
                           # tolerance for linear Krylov iteration.
epnldpk      real /1.e0/   # optional input for daspk when info(13) = 1
                           # tolerance for Newton iteration convergence.

srtolpk      real /1.e-4/  #del=srtolpk*rtol for num. diff. (daspk,vodpk)
                           #Now set internally as srtolpk=del/rtolv
efacn        real /1.e0/   #scaling factor for Newton error test in vodpk
ftol         real /1.e-8/  #stop tolerance of su*f for nksol (=epsmch**(1/3))
                           #  ( maxnorm(f) .le. ftol to stop. )
stptol       real /0.e0/   #stop tolerance of yl(k)-yl(k-1) in nksol
                           #  ( maxnorm(yl(k)-yl(k-1)) .le. stptol to stop. )
epscon1      real /1.e-1/  #linear solve tolerance in nksol, epsfac =
                           # epscon1*min(epscon2,frnm)
epscon2      real /1.e-2/  #linear solve tolerance in nksol, epsfac =
                           # epscon1*min(epscon2,frnm)
iterm        integer     #output flag for nksol
mdif         integer /0/ #nksol flag for user-supplied j*v product (0=internal)
ipflag       integer /1/ #nksol flag to precondition (1=yes)
mfnksol      integer /-3/#nksol method flag; =1 means dogleg strategy,
                         #=2 means linesearch with Arnoldi method,
                         #=3 means linesearch with GMRES method.
			 #=4 full direct solve by RSmirnov;set premeth=banded
                         #negative mfnksol ignores global constaints
xiprint       integer /1/ #nksol optional statistics flag.
                         #=0 means no optional statistics are printed.
                         #=1 means iteration count, norm of F(u) and
                         #   no. of F evaluations are printed.
                         #=2 means irpint=1 statistics are printed, and
                         #   statistics regarding the convergence of the
                         #   Krylov iteration, dogleg strategy.  See
                         #   nksol documentation for more details.
itermx       integer /30/# maximum number of nonlinear iterations for nksol.
stepmx       real /1.e9/ # maximum length of a Newton step for nksol.
del2nksol    real /1.e-14/ # if nonzero, size of del**2 for diff. quot. Jac
taunksol     real /1.e9/ # initial size of trust region for dogleg strategy
                         # (mfnksol = 1) in nksol.
incpset      integer /5/ # maximum number of nonlinear iterations before
                         # the preconditioner is reevaluated within nksol.
ismmaxuc     integer /1/ #=1 for calc. mmaxu internally from nx and ny
mmaxu        integer /25/# maximum Krylov subspace dimension.
                         # currently, only used in nksol			                         # If ismmaxuc=1, calc. internally; ismmaxuc=0 use input
icntnunk     integer /0/ #nksol continuation call flag.
                         #=1 tells nksol not to call the preconditioner routine
                         #   pset on the current call.  In this case, nksol
                         #   assumes that the preconditioner was evaulated
                         #   on an earlier call, and is to be used for as
                         #   many steps as it is successful on this call.
                         #=0 tells nksol that this is not a continuation call.
                         #   The preconditioner routine pset is called to
                         #   evaluate and factor the Jacobian matrix.

***** Parallv:
#Variables used for the parallel version utilizing pvode or kinsol
nlocal       integer       #number of equations on given processor
neqg         integer       #total number of equations over all processors
nxg	     integer       #number of global poloidal mesh points = nxg+2
nyg	     integer       #number of global radial mesh points = nyg+2g
meth         integer /2/   #input for fpvmalloc; spec. method (lmm)
itmeth       integer /2/   #input for fpvmalloc; spec. interation method (iter)
iatol        integer /1/   #input for fpvmalloc; spec. error array type
igs          integer /1/   #input for fcvspgrm2; Gram-Schmidt process
maxkd        integer /50/  #maximum Krylov dimension for kinsol
maxlrst      integer /2/   # for kinsol
msbpre       integer /0/   #preconditioner flag for kinsol
globalstrat  integer /0/   #global strategy flag for kinsol
iopt(40)     integer /40*0/#opt. input/output array
ropt(40)     real  /40*0./ #opt. input/output array
rtol_pv      real /1e-4/   #relative tol. for parallel pvode
atol_pv      real /1e-6/   #relative tol. for parallel pvode
delt_pv      real /5e-4/   #linear converg. error-test param. for pvode

***** Constraints:
#Variables for checking constraints, i.e., negativity
icflag          integer /1/ #flag to use constraint that ni, etc. not < 0
                            #=1 turns on for nksol(with rlx) and vodpk(no rlx)
                            #=2 adds rlx constraint to vodpk
rlx             real  /.4/  #fractional change allowed per iteration
rlxv            real  /.4/  #fractional change in up allowed for svrpkg=newton
icnstr(neqmx)  _integer     #nksol constraint array; =1 means must be > 0
constr(neqmx)  _real        #kinsol constraint array; =1 means must be > 0
ylprevc(neqmx) _real        #yl vector from previous call from vodpk
ylchng(neqmx)  _real        #change in yl from prev. step, yl-ylprevc
adjf1           real  /1.2/ #if mfnksol=3 glob strat, frnm_new/adjf1>=fnrm_old


***** Ynorm:
#Variables for the normalization of the y's:
iscolnorm	integer /3/	# =0 for no implicit scaling (suscal=1)
				# =1 for scaling by normalization constants
				# =2 for scaling by max(abs(yl),floors)
				# =3 combination of global scaling with nnorm,
                                #    etc, followed by local scaling by each yl
norm_cons(numvar)	_real	# normalization constants (calculated)
floor_cons(numvar)	_real	# floor constants (calculated)
var_scale_floor	  real	/1.e-7/ # factor from norm_c to floor_c except for up
	# factor multiplied by normalization constants to get floors for scaling
vsf_up            real  /1.e0/  #var_scale_floor factor for up eqns
vsf_phi           real  /1.e0/  #var_scale_floor factor for phi eqns
n0(nispmx)        real  [m**-3] /nispmx*1.e20/ #normalization ion density
n0g(ngspmx)       real  [m**-3] /ngspmx*1.e20/ #normalization gas density
temp0           real  [eV]     /100./    #normalization temperature
isflxvar        integer        /0/       #sets variables for ODE, Jacobian
                                         #=1 for yl=n,nv,nT; =0 for yl=n,v,T
                                         #=2 for yl=n,v,nT
isrscalf        integer        /1/       #rescales ODE rhs if isflxvar.ne.1
dx0             real  [m]      /30./     #norm. grid spacing factor for phi eqn
nnorm           real  [m**-3]            #normalization density(calc)
ennorm          real  [J/m**3]           #normalization energy density(calc)
sigbar0         real  [Mho/m]            #normalization parallel cond. (calc)
vpnorm          real  [m/s]              #normalization ion paral.velocity(calc)
fnorm(1:nusp)  _real  [kg/m**2 s]        #normalization momentum flux(calc)
suscal(neqmx)  _real      #scale factors for yl's in nksol routine (nominally=1)
sfscal(neqmx)  _real      #scale factors for f's in nksol routine (nominally=1)
rdoff           real /0./    #ranf-induced roundoff error compared to unity
yloext(neqmx)  _real         #last var for d(yl)/dt set externally
isyloext        integer /0/  #=1 allows d(yl)/dt using ext. yloext

***** Phyvar:
#Values of some physical constants
pi              real            /3.14159265358979323e0/ #Pi
me              real            /9.1096e-31/            #Electron mass
mp              real            /1.6726e-27/            #Proton mass
ev              real            /1.6022e-19/            #1 electron volt
qe              real            /1.6022e-19/            #Elementary charge
mu0             real            /1.25663706e-6/         #Vac. magnetic perm.
eps0            real            /8.8542e-12/            #Vac. dielectric perm.
rt8opi          real            /1.595769121606e0/      #sqrt(8/pi)

***** Comtra restart:
#Variables that contain the transport parameters.
parvis(1:nispmx)  real          /nispmx*1./     +input #factor times parallel visc.
travis(1:nispmx)  real [m**2/s] /nispmx*1./     +input #value of perp. visc.
difni(1:nispmx)   real [m**2/s] /nispmx*0.5/    +input #value of density radial diff. coef.
dif4order(1:nispmx) real [m**4/s] /nispmx*0./   +input #4th ord ion density radial diff. coef.
difgy4order(1:ngspmx) real [m**4/s] /ngspmx*0./ +input #4th ord gas density radial diff. coef.
difgx4order(1:ngspmx) real [m**4/s] /ngspmx*0./ +input #4th ord gas density poloid diff. coef.
difax(1:nispmx)   real [m*2/s] /nispmx*0./      +input #poloid. diff coeff.,scaled with dn/dx
difnit(1:nispmx)  real [none]   /nispmx*0./     +input #turb. radial diff. multiplier
cdifnit           real [none]   /1./            +input #=1 for all turb., =0 all fixed, D coef
difpr(1:nispmx)   real [m**2/s] /nispmx*0./     +input #value of pressure radial diff. coef.
difni2(1:nispmx)  real [m**2/s] /nispmx*0./     +input #value of density e_ll x e_r diff. coef.
difpr2(1:nispmx)  real [m**2/s] /nispmx*0./     +input #value of pressure e_ll x e_r diff. coef
difutm(1:nispmx)  real [m**2/s] /nispmx*1./     +input #value of toroidal mom. diff. coef
difniv(0:ny+1,1:nisp)  _real [m**2/s]   /0./	+input #dens diff. if isbohmcalc=3, varys w/B
difprv(0:ny+1,1:nisp)  _real [m**2/s]   /0./	+input #press diff if isbohmcalc=3, varys w/B
difniv2(0:ny+1,1:nisp) _real [m**2/s]   /0./	+input #dens2 diff. if isbohmcalc=3, varys w/B
travisv(0:ny+1,1:nisp) _real [m**2/s]   /0./	+input #viscosity. if isbohmcalc=3, varys w/B
kyev(0:ny+1)	_real [m**2/s]   /0./	+input #elec eng chi if isbohmcalc=3, vary(0:ny+1)s w/B
kyiv(0:ny+1)	_real [m**2/s]   /0./	+input #ion eng chi if isbohmcalc=3, varys w/B
difutmv(0:ny+1,1:nisp)	_real [m**2/s]   /0./	+input #dens diff. if isbohmcalc=3, varys w/B
vconyv(0:ny+1,1:nisp) _real [m/s] /0./   +input #convec radial vel if isbohmcalc=3
inbtdif		real		 /2/	+input #if isbohmcalc=3, D,chi ~1/Bt**inbtdif
inbpdif		integer		 /1/	+input #if isbohmcalc=3, D,chi ~1/Bp**inbpdif
ixbpmin         integer          /0/    +input #isbohmcalc=3, min bpol(ixpt2-ixbpmin,
isbohmcalc      integer          /1/	+input 
                    #if=1, calc Bohm diff if facb... > 0
					#if=2, harmonic ave of Bohm, difni, etc.
					#if=3, D=difniv*(B0/B)**inbdif, etc
facbni        real [ ]      /0./	+input #factor for Bohm density y-diff. coeff.
facbup        real [ ]      /0./	+input #factor for Bohm parll v y-diff. coeff.
facbni2       real [ ]      /0./	+input #factor for Bohm density 2-diff. coeff.
facbee        real [ ]      /0./	+input #factor for Bohm Te diff. coeff.
facbei	      real [ ]      /0./	+input #factor for Bohm Ti diff. coeff.
vcony(1:nispmx) real [m/s] /nispmx*0./  +input #value of constant radial velocity
difcng	      real [m**2/s] /50./       +input #constant gas diff. coeff if isgasdc=1
isgasdc       integer       /0/         +input #switch to turn on constant gas dif coef
flalfe        real          /0.21/      +input #|| heat flux limit factor for elec.
flalfi        real          /0.21/    +input #|| heat flux limit factor for ions
lxtemax       real [m]      /1.e10/	+input #max pol. scale len of elec heat-flux lim
lxtimax       real [m]      /1.e10/	+input #max pol. scale len of ion heat-flux lim
lxtgmax       real [m]      /1.e10/	+input #max pol. scale len of gas heat-flux lim
flalftf	      real          /1.e20/     +input #elec. thermal force flux-lim factor
flgam         real          /1./        +input #exponent for ion flux-limit expression
flgamv        real          /2./        +input #exponent for vel flux-limit expression
flgamg        real          /2./        +input #exponent for gas dens flux-limit
flgamvg       real          /2./        +input #exponent for gas visc flux-limit
flgamtg       real          /1./        +input #exponent for gas temp flux-limit
fricflf       real          /1./        +input #flux-limiting factor for inputs to
                                        #multispecies friction (and upi) calc
isflxlde      integer       /0/         +input #=1,elec flux limit diff;=0, conv/diff
isflxldi      integer       /2/         +input #=1,ion flux limit diff;=0, conv/diff
                                        #=2, diff on individ hxcij
kxe           real          /1./        +input #pol Braginsk elec heat conduc factor;
                                        #prev 1.35->Balescu explain by M.Zhao 
alfkxi        real          /0./        +input #reduces ion thermal conduc, K_||, if
                                        #|ti(ix+1)-ti(ix)|<alfkxi*ti(ix)
alfkxe        real          /0./        +input #reduces elec thermal conduc, K_||, if
                                        #|te(ix+1)-te(ix)|<alfkxe*te(ix)
rkxecore      real          /1./        +input #pol elec heat diff. reduc fac in core
inkxc         integer       /2/         +input #expon on yyf/yyf(0) fac for core kxe
kye           real [m**2/s] /0.7/       +input #radial electron heat diffusivity
kye4order     real [m**2/s] /0./        +input #4th order Te radial diff. coef.
kyet          real [none]   /0./        +input #turb. radial elec. heat diff. multiplier
ckyet         real [none]   /1./        +input #=1 for all turb., =0 all fixed, chi_e
kxi           real          /1./        +input #poloidal ion heat diff. multi. fac
kxicore       real          /1./        +input #poloidal ion heat diff. factor in core
kxn           real          /0./        +input #poloidal cx-neutral heat diff. factor
kyi           real [m**2/s] /.2/        +input #radial ion heat diffusivity
kyi4order     real [m**2/s] /0./        +input #4th order Ti radial diff. coef.
kyit          real [none]   /0./        +input #turb. radial ion heat diff. multiplier
ckyit         real [none]   /1./        +input #=1 for all turb., =0 all fixed, chi_i
kyn           real          /0./        +input #radial cx-neutral heat diff. factor
feqp          real          /1./        +input #(Te-Ti) equipartition multiplier
alfeqp        real          /0./        +input #reduces equipart. term if te~ti
flalfgx(10)   real          /10*1./     +input #poloidal gas diff. flux limit
flalfgy(10)   real          /10*1./     +input #radial gas diff. flux limit
flalfgxy(10)  real          /10*1./     +input #nonorthog pol-face gas flux limit
flalfgnx      real          /1.e20/	+input #flux-limit on total fngx;for safety
flalfgny      real          /1.e20/	+input #flux-limit on total fngy;for safety
gcfacgx	      real	    /1./	+input #mult tot conv gas x-flux at ix=0 & nx
gcfacgy	      real	    /1./	+input #mult tot conv gas y-flux at iy=0 & ny
gcfacgtx      real	    /0./	+input #mult grad Ti conv gas x-flux ix=0 & nx
gcfacgty      real	    /0./	+input #mult grad Ti conv gas y-flux iy=0 & ny
isdifxg_aug   integer       /0/         +input #=1 enhances D_xgas with flx-lim factor
isdifyg_aug   integer       /0/         +input #=1 enhances D_ygas with flx-lim factor
flalfv        real          /0.5/    +input #parallel velocity flux limit factor
isupdrag      integer       /0/         +input #=1 adds nonunif B-field drag on v_||
con_leng      real          /1e20/ [m]  +input #connect length used for coll trans fac
frac_pt       real          /0.3/       +input #invers aspect ratio for colliless drag
flalftgx      real          /1./        +input #poloidal atom temp. diff. flux limit
flalftgy      real          /1./        +input #radial atom temp. diff. flux limit
flalftmx      real          /1./        +input +input #poloidal mol temp. diff. flux limit
flalftmy      real          /1./        +input #radial mol temp. diff. flux limit
flalfvgx      real          /1./        +input #poloidal gas parall viscosity flux lim
flalfvgy      real          /1./        +input #radial gas parall viscosity flux limit
flalfvgxy     real          /1./        +input #FL for neutral fmixy
flalftxy      real          /1.e20/	+input #FL for feixy (addition to hcy FL)
flalftgxy(ngspmx) real   /ngspmx*1./    +input #FL for fegxy
rnn2cx	      real          /0.2/       +input #ratio of neut-neut coll. to cx coll.
rscat2cx      real          /1./        +input #fraction of cx counted as visx scatt
sigcx         real   [m**2] /1.e-18/    +input #cx cross-sect if icnucx=2
sigcxms(nisp,ngsp)   _real [m**2] /1e-19/ +input #cx x-sect for (ifld,igsp) coll
rcxighg(ngspmx) real      /ngspmx*0./   +input #ratio of charge-exchange rate for
                                        +input #ng_imp+ni_hydrn -> ng_hydrn+ni_Z=1_imp
kelhihg         real [m**3/s] /5e-16/   +input #elastic coll. coeff:hyd_ion+hyd_atom
kelhghg         real [m**3/s] /5e-16/   +input #elastic coll. coeff:hyd_atm+hyd_atom
kelhmhg         real [m**3/s] /5e-16/   +input #elastic coll. coeff:hyd_mol+hyd_atom
kelhmhm         real [m**3/s] /5e-16/   +input #elastic coll. coeff:hyd_mol+hyd_mol
kelighi(ngspmx) real [m**3/s]/0.,5*5e-16/+input #elastic coll. coeff:imp_gas+hyd_ion
kelighg(ngspmx) real [m**3/s]/0.,5*5e-16/+input #elastic coll. coeff:imp_gas+hyd_gas
keligii(ngspmx) real [m**3/s]/0.,5*5e-16/+input #elastic coll. coeff:imp_gas+imp_ion
keligig(ngspmx) real [m**3/s]/6*5e-16/  +input #elastic coll. coeff:imp_gas+imp_gas
cfmassfac       real         /1./       +input #scales elas scat factor 16mi/(3mg+mi)
sigvi_floor   real [m**3/s]  /0./       +input #minimum of ioniz. rates allowed(1e-18)
fupe_cur      real []        /1./       +input #=1 fixes cur err to upe for isimpon=6
diffusivity(0:nx+1,0:ny+1)  _real 
   # anomalous (turbulent) diffusivity (calculated during rhs eval)
diffusivwrk(0:nx+1,0:ny+1)  _real -restart
   # anomalous (turbulent) diffusivity (mixed w/difni using cdifnit)
diffusivloc(0:nx+1,0:ny+1)  _real -restart
   # anomalous (turbulent) diffusivity (local values for isturbcons=2)
cfnus_e         real      /1.e20/ +input # factor mult nu_star_e for elec coll_fe
cfnus_i         real      /1.e20/ +input # factor mult nu_star_i for ion coll_fi
coll_fe(0:nx+1,0:ny+1) _real      +input # nu_star_e/(1+nu_star_e) for elec CF drifts
coll_fi(0:nx+1,0:ny+1) _real      +input # nu_star_i/(1+nu_star_i) for ion CF drifts
tibsep	 	[eV] real /100./  +input # Ion temp on sep for banana width in lconi
tebsep	 	[eV] real /100./  +input # Elec temp on sep for banana width in lcone
cfelecbwd	     real /10./	  +input # Factor for elec banana width in lcone; makes
				  # elec banana width not to small for mesh
fluxfacy        real  /1./        +input # multiples y-fluxes & v(dP/dy) for 1D sims
isdifbetap      integer /0/       +input #=1 turns on betap-dependent & difniv diffusion
iexpbp          real   /1./       +input #exponent for diff ~ betap**iexpbp
dfacbp          real [m**2/s] /0./+input #diff. coeff dens *betap**iexpbp;dif_use,dif_use
trfacbp         real [m**2/s] /0./+input #diff. coeff up *betap**iexpbp;tray_use,trax_use
kefacbp         real [m**2/s] /0./+input #diff. coeff Te *betap**iexpbp;kye_use,kxe_use
kifacbp         real [m**2/s] /0./+input #diff. coeff Ti *betap**iexpbp;kyi_use,kxi_use
flalfea(0:nx+1)      _real /0./+input #calc:elec thermal flux-limit array (see flalfe)
flalfia(0:nx+1)      _real /0./+input #calc:ion thermal flux-limit array (see flalfi)
flalfva(0:nx+1)      _real /0./+input #calc:ion visc flux-limit array (see flalfv)
flalfgxa(0:nx+1,10)  _real /0./+input #calc:neut pol flux-limit array (see flalfgx)
flalfgxya(0:nx+1,10) _real /0./+input #calc:neut xy flux-limit array (see flalfgxy)
flalfgya(0:ny+1,10)  _real /0./+input #calc:neut rad flux-limit array (see flalfgy)
flalfvgxa(0:nx+1)  _real /0./+input #calc:neut part pol flux-limit array (see flalfgx)
flalfvgxya(0:nx+1) _real /0./+input #calc:neut part xy flux-limit array (see flalfgxy)
flalfvgya(0:ny+1)  _real /0./+input #calc:neut part rad flux-limit array (see flalfgy)
flalftgxa(0:nx+1)  _real /0./+input #calc:neut part pol flux-limit array (see flalfgx)
flalftgxya(0:nx+1) _real /0./+input #calc:neut part xy flux-limit array (see flalfgxy)
flalftgya(0:ny+1)  _real /0./+input #calc:neut part rad flux-limit array (see flalfgy)
cfmolcool       real     /0./+input  #scale factor for molec cooling if ishymol=1


***** Interprettrans:
# Variables used for interpretive analysis of expt profiles
isdifuseinterp  integer         /1/    # =1 use dif_int in dif_use, etc, only
isadjsolprof    integer         /1/   # adjs SOL nis,tes,tis to smooth connect
denrdrop        real            /.9/  # rel drop in nis rad prof for adjsolprof
terdrop         real            /.9/  # rel drop in tes rad prof for adjsolprof
tirdrop         real           /.95/  # rel drop in tis rad prof for adjsolprof
maxchgdiff      real            /.5/   # max change allowed in dif_int per cell
dif_use_max real       [m**2/s] /1./   # max used values of diff. coeff
dif_use_min real       [m**2/s] /.02/  # max used values of diff. coeff
kye_use_max real       [m**2/s] /3./   # max used values of kye coeff
kye_use_min real       [m**2/s] /.1/   # max used values of kye coeff
kyi_use_max real       [m**2/s] /3./   # max used values of kyeicoeff
kyi_use_min real       [m**2/s] /.1/   # max used values of kye coeff
difni_sol   real       [m**2/s] /0./   # added to SOL dif_use if interp mode
difni_pf    real       [m**2/s] /0./   # added to PF dif_use if interp mode
kye_sol     real       [m**2/s] /0./   # added to SOL kye_use if interp mode
kye_pf      real       [m**2/s] /0./   # added to PF kye_use if interp mode
kyi_sol     real       [m**2/s] /0./   # added to SOL kyi_use if interp mode
kyi_pf      real       [m**2/s] /0./   # added to PF kye_use if interp mode

##taudndt  real           [s]    /1.e10/ # global density rise-time
del_sp(0:ny+1)    _real [1/m**3s] /0./ # dN/dt for ion source term
del_witot(0:ny+1) _real [W/m**3]  /0./ # Total dWi/dt; not used
del_wetot(0:ny+1) _real [W/m**3]  /0./ # Total dWe/dt; not used
del_dndt(0:ny+1)  _real [1/m**3s] /0./ # dN/dt deduced from data (input)
del_deedt(0:ny+1) _real [W/m**3s] /0./ # dNTe/dt deduced from data (input)
del_deidt(0:ny+1) _real [W/m**3s] /0./ # dNTi/dt deduced from data (input)
del_wicd(0:ny+1)  _real [W/m**3]  /0./ # dWi/dt for ions rad diff heat flux
del_wicv(0:ny+1)  _real [W/m**3]  /0./ # dWi/dt for ions rad conv heat flux
del_wecd(0:ny+1)  _real [W/m**3]  /0./ # dWe/dt for elec rad diff heat flux
del_wecv(0:ny+1)  _real [W/m**3]  /0./ # dWe/dt for elec rad conv heat flux
del_wicdd(0:ny+1) _real [W/m**3]  /0./ # diag: dWi/dt ion rad diff heat flux
del_wecdd(0:ny+1) _real [W/m**3]  /0./ # diag: dWe/dt elec rad diff heat flux
del_cei(0:ny+1)   _real [W/m**3]  /0./ # dWe/dt=-dWi/dt elec-ion coll exchange
del_wivdp(0:ny+1) _real [W/m**3]  /0./ # dWi/dt p-v type pressure-work terms
del_wevdp(0:ny+1) _real [W/m**3]  /0./ # dWe/dt p-v type pressure-work terms
dif_int(0:ny+1)   _real	[m**2/s]  /0./ # interp particle diff, D
kyi_int(0:ny+1)   _real	[m**2/s]  /0./ # interp ion radial conduc., Chi_i
kye_int(0:ny+1)   _real	[m**2/s]  /0./ # interp elec radial conduc., Chi_e
vyn_int(0:ny+1)   _real	[m/s]     /0./ # interp ion radial particle drift vel
vyei_int(0:ny+1)  _real	[m/s]     /0./ # interp ion radial energy drift vel
vyee_int(0:ny+1)  _real	[m/s]     /0./ # interp elec radial energy drift vel
gamp(0:ny+1)      _real [1/m**2s] /0./ # radial ion particle flux
gamei(0:ny+1)     _real [W/m**2s] /0./ # ion radial energy flux
gamee(0:ny+1)     _real [W/m**2s] /0./ # elec radial energy flux
pfmpg(0:ny+1)     _real [1/m**2s] /0./ # radial neutral particle flux
facgam(0:nx+1,0:ny+1) _real       /0./ # geom factor for flux-surf averaging
floyd(0:nx+1,0:ny+1)  _real [1/s] /0./ # equiv radial ion particle current
							
***** Turbulence:
# Variables used in calculating anomalous diffusivities
kappabar      real [1/m]  /0.003/ +maybeinput # field-line avg'd curvature
lambdan       real [none]   /4./  +maybeinput # dens(divertor) / dens(midplane)
lambdat       real [none]   /4./  +maybeinput # temp(midplane) / temp(divertor)
gammasi       real [none]   /0./  +maybeinput # secondary emis. coef. from ion bombard.
lambdap       real [none]   /3.5/ +maybeinput # e (dPhi0 / dr0) / (dTed / dr0)
suppress_lmode   integer    /0/   +maybeinput # =1 to suppress L-mode turbulence in SOL
maxmag_lmode  real          /2./  +maybeinput # max magn. of step in bracketing kymax
nky           integer       /30/  +maybeinput # number of ky values in search for kymax
kybeg         real [none]  /0.05/ +maybeinput # lower limit of acceptable kymax
kyend         real [none]   /3./  +maybeinput # upper limit of acceptable kymax
kya           real [none]   /1.0/ +maybeinput # one initial point in search for kymax
kyb           real [none]   /1.1/ +maybeinput # other initial point in search for kymax
iprint_lmode  integer  /0/  +maybeinput # =1 for diagnostic output, =2 for more output
tol_lmode     real [none] /1.e-6/ +maybeinput # abs & rel tolerance in search for kymax
isturbnloc    integer  /1/  +maybeinput # =1 to turn on nonlocal dependence of D,chi
isturbcons    integer  /1/  +maybeinput # =1 to make turbulent D,chi const within SOL
                            # =2 to apply radial digital filter to D,chi
diffusrange   real [m]  /0.01/  +maybeinput # radial range of digital filter
diffuslimit   integer       +maybeinput # no. of surfaces in half-range of digital filter
diffuswgts(-9:9)  real [none]  +maybeinput # weights for radial digital filter
islmodebeta   integer  /1/  +maybeinput # =1 to turn on finite-beta correction
gradvconst    real [none] /0.005/ +maybeinput # factor involving rad. grad. of v(parallel)

***** Turbulence_comm:
# Communication of variables to minimization routine for turbulent growth rate
epsilon       real [none]         # rhos / lte
turbdelta     real [none]         # param. depending on temp. & length ratios
ssqthsqavg    real [none]   /1./  # field-line avg of s**2 * theta**2
kxconst       real [none]  /2.47/ # constant appearing in calc. of kx
cubrtnu       real [none]         # cube root of parameter nu
bcoef0        complex [none]
ccoef1        complex [none]
ccoef2        real [none]
ccoef3        real [none]

***** Turbulence_diagnostics:
chinorml(0:nx+1,0:ny+1)  _real # norm. anom. diffusivity (L-mode turbulence)
chinormh(0:nx+1,0:ny+1)  _real # norm. anom. diffusivity (H-mode turbulence)

***** Timary:
#Vars stored at output times for time-dependent calculations; also control vars
nsteps     integer         /100/  #number of logarithmically spaced output times
n_stor	   integer         /1/    #number of storage pts for solution
dt_init_rundt  real      /1.e30/  #init dtreal for first call to rundt
nist1(nsteps,0:nx+1,0:ny+1,nisp) _real [1/m**3] #density for ODE output times
upst1(nsteps,0:nx+1,0:ny+1,nisp) _real [m/s]    #parallel vel for ODE output times
test1(nsteps,0:nx+1,0:ny+1)      _real [eV]     #Te at for ODE output times
tist1(nsteps,0:nx+1,0:ny+1)      _real [eV]     #Ti at for ODE output times
ngst1(nsteps,0:nx+1,0:ny+1,ngsp) _real [1/m**3] #ng at for ODE output times
phist1(nsteps,0:nx+1,0:ny+1)  _real [V]      #phi at various times
toutlsod(nsteps)              _real [s]      #time which nist1, etc, are filled
yldnmx(nsteps)                _real [1/s]    #max. rate-of-change, yldot/yl
iyldnmx(nsteps)               _integer       #index of vector yldnmx
istep_nk                       integer  /0/  #array index for time-dep. nksol
nsteps_nk                      integer  /1/  #number of nksol time-steps
ni_stor(n_stor,0:nx+1,0:ny+1,nisp) _real [1/m**3] #density for rundt output times
up_stor(n_stor,0:nx+1,0:ny+1,nisp) _real [m/s]    #par vel for rundt output times
te_stor(n_stor,0:nx+1,0:ny+1)      _real [eV]     #Te for rundt output times
ti_stor(n_stor,0:nx+1,0:ny+1)      _real [eV]     #Ti for rundt output times
ng_stor(n_stor,0:nx+1,0:ny+1,ngsp) _real [1/m**3] #ng for rundt output times
phi_stor(n_stor,0:nx+1,0:ny+1)     _real [V]      #phi for rundt times
tim_stor(n_stor)		   _real [t]      #output times for rundt
nfe_stor(n_stor)		   _real          #num func evals in rundt interval
dtreal_stor(n_stor)		   _real [t]	  #dtreal for rundt output times
rdtphidtr	real  /1.e20/  	#ratio dtphi/dtreal
ismfnkauto	integer /1/	# if =1, mfnksol=3 for dtreal<dtmfnk3, otherwise=-3
dtmfnk3		real /1.e-3/ 	# dtreal for mfnksol sign change if ismfnkauto=1
mult_dt         real /3.4/	# factor expanding dtreal after ii2max steps
ii1max 		integer /100/	# number of changes to dtreal
ii2max  	integer /5/	# number of timesteps at current dtreal
itermxrdc       integer /7/     # value of itermx used by rdcontdt
itermx_dt       integer /10/    # sets itermx for rundt; max iters per dt try
ftol_dt 	real /1.e-5/	# fnrm tolerance for the time-dependent steps
ftol_min 	real /1.e-9/	# value of fnrm where time advance will stop
dt_tot 		real /0./	# total time accumulated for run (output, not input)
t_stop 		real /100./	# value of dt_tot (sec) where calculation will stop
dt_max 		real /100./	# max time step for dtreal
dt_kill 	real /1.e-14/	# min allowed time step; rdcontdt stops if reached
deldt_min 	real /0.04/	# min relative change allowed for model_dt > 0
initjac 	integer /0/	# if=1, calc initial Jac upon reading rdcontdt
irev  		integer /-1/	# flag allows reduced dt advance after cutback
numrevjmax      integer /2/     # num dt reducts before Jac recalculated
numfwdjmax      integer /1/	# num dt increases before Jac recalculated
numrev          integer /0/     # count dt reducts in rdcontdt
numfwd          integer /0/     # count dt increases in rdcontdt
numrfcum        integer /0/     # number of cumulative reducts/increase in dt
tstor_s 	real /1.e-3/	# beginning time for storing solution
tstor_e  	real /4.e-2/	# ending time for storing solution
ipt 		integer /1/	# index of variable printed at each out timestep
savefname       character*5 /"it333"/ # name of pfb save file pfdt_"savefname"
iprtrundt       integer /1/     # =1, then print rundt diag; .ne.1, no printing

***** Compla:
#Variables in common -- plasma parameters
mi(1:nisp)         _real [kg] /1.67e-27/ #ion mass in kg, calculated from minu
zi(1:nisp)         _real [ ]  /1./       #ion charge number, calc. from ziin
mg(1:ngsp)         _real [kg] /1.67e-27/ #gas species mass, calc. fr minu
facmg(1:nispmx)        real /nispmx*1./  #scale factor for mg to recov old case
znucl(1:nisp)              _integer [ ]   #tot. nucl. charge, calc. from znuclin
ni(0:nx+1,0:ny+1,1:nisp)   _real  [m^-3]  #ion density in primary cell (ix,iy)
lni(0:nx+1,0:ny+1,1:nisp)  _real  [m^-3]  #log(ion dens) in prim. cell (ix,iy)
nm(0:nx+1,0:ny+1,1:nisp)   _real [kg*m^-3]#mass density [nm(,,1) is sum, exclud.
                                          #gas, if nusp=1, isimpon=5] in cell
nz2(0:nx+1,0:ny+1)         _real  [m^-3]  #sum of ni*zi**2 over all ion species
uu(0:nx+1,0:ny+1,1:nisp)   _real  [m/s]   #ratio ion-flux/density at x-face;
                                          #if orthog mesh, poloidal ion velocity
uup(0:nx+1,0:ny+1,1:nisp)  _real  [m/s]   #poloidal ion vel (|| flow contrib)
up(0:nx+1,0:ny+1,1:nisp)   _real  [m/s]   #par ion vel if full mom eqn on
                                          # (mass-dens. avg if isimpon = 5)
upi(0:nx+1,0:ny+1,1:nisp)  _real  [m/s]   #inter. par ion vel even if force bal
upifmb(0:nx+1,0:ny+1,1:nisp) _real [m/s]  #par ion vel fmombal if isimpon=5
uz(0:nx+1,0:ny+1,1:nisp)   _real  [m/s]   #toroidal ion vel in pol X rad direct
v2(0:nx+1,0:ny+1,1:nisp)   _real  [m/s]   #vel normal to parallel & rad. direc.
v2xgp(0:nx+1,0:ny+1,1:nisp) _real [m/s]   #v2 ion vel for v2x_gradx_P eng terms
v2ce(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of v2 from ExB
v2cb(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of v2 from grad_B
ve2cb(0:nx+1,0:ny+1)       _real  [m/s]   #electron v2 from grad_B
v2cd(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of ion v2 from grad_PxB
ve2cd(0:nx+1,0:ny+1,1:nisp) _real [m/s]   #portion of elec v2 from grad_PxB
q2cd(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #ion heat flux from grad_PxB
v2rd(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of v2 from resistive drift
v2dd(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of v2 from anomalous drift
vy(0:nx+1,0:ny+1,1:nisp)   _real  [m/s]   #radial ion velocity
vygp(0:nx+1,0:ny+1,1:nisp) _real  [m/s]    #radial ion vel for vy_grady_P eng terms
vytan(0:nx+1,0:ny+1,1:nisp)_real  [m/s]   #radial ion vel.*tan(vtag) on x-face
vygtan(0:nx+1,0:ny+1,1:ngsp)_real [m/s]   #radial gas grad-T vel.*tan(vtag) on
					  #x-face
vyce(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of vy from ExB
vycb(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of vy from grad_B
veycb(0:nx+1,0:ny+1)       _real  [m/s]   #electron vy from grad_B
vycp(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #ion vy from grad_PixB
veycp(0:nx+1,0:ny+1)       _real  [m/s]   #electron vy from grad_PeXB
vyrd(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of vy from resistive drift
vydd(0:nx+1,0:ny+1,1:nisp) _real  [m/s]   #portion of vy from anomalous drift
vyavis(0:nx+1,0:ny+1,1:nisp) _real [m/s]  #rad vel from anom perp vis (ExB,P)
vex(0:nx+1,0:ny+1)         _real  [m/s]   #Poloidal electron velocity
upe(0:nx+1,0:ny+1)         _real  [m/s]   #parallel electron velocity
vep(0:nx+1,0:ny+1)         _real  [m/s]   #old parallel electron velocity-remove
ve2(0:nx+1,0:ny+1)         _real  [m/s]   #old "2" electron velocity-remove
vey(0:nx+1,0:ny+1)         _real  [m/s]   #Radial electron velocity
vycf(0:nx+1,0:ny+1)	   _real  [m/s]   #radial vel from class. viscosity
vycr(0:nx+1,0:ny+1)	   _real  [m/s]   #radial vel from class. thermal force
te(0:nx+1,0:ny+1)          _real  [J]	  #electron temperature in primary cell
ti(0:nx+1,0:ny+1)          _real  [J]	  #ion temperature in primary cell
ng(0:nx+1,0:ny+1,1:ngsp)   _real  [m^-3]  #gas density in primary cell (ix,iy)
lng(0:nx+1,0:ny+1,1:ngsp)  _real  [m^-3]  #log(gas dens) in prim. cell (ix,iy)
uug(0:nx+1,0:ny+1,1:ngsp)  _real  [m/s]   #ratio gas-flux/density at x-face;
                                          #if orthog mesh, poloidal gas velocity
uuxg(0:nx+1,0:ny+1,1:ngsp) _real  [m/s]   #poloidal-only component of uug;
                                          #for uuxg*gradx_Pg coomp of seic
vyg(0:nx+1,0:ny+1,1:ngsp)  _real  [m/s]   #radial gas velocity
tg(0:nx+1,0:ny+1,1:ngsp)   _real  [J]	  #gas temperature in primary cell
istgcon(ngspmx)       real /ngspmx*0/ #=0, set tg(,,i)=rtg2ti*ti; if >0, set
                                      #tg=(1-istgcon)*rtg2ti*ti+istgcon*tgas*ev
tev(0:nx+1,0:ny+1)         _real  [J]	  #ion temperature at vertex of cell
niv(0:nx+1,0:ny+1,1:nisp)  _real  [m^-3]  #ion dens up-right vert[rm,zm(,,4)]
upv(0:nx+1,0:ny+1,1:nisp)  _real  [m/s]   #ion par vel up-right vert[rm,zm(,,4)]
ngv(0:nx+1,0:ny+1,1:ngsp)  _real  [m^-3]  #gas dens up-right vert[rm,zm(,,4)]
tiv(0:nx+1,0:ny+1)         _real  [J]	  #ion temperature at vertex of cell
niy0(0:nx+1,0:ny+1,1:nisp) _real  [m^-3]  #ion density below y-face center
niy1(0:nx+1,0:ny+1,1:nisp) _real  [m^-3]  #ion density above y-face center
niy0s(0:nx+1,0:ny+1,1:nisp) _real [m^-3]  #old ion density below y-face center
niy1s(0:nx+1,0:ny+1,1:nisp) _real [m^-3]  #old ion density above y-face center
ney0(0:nx+1,0:ny+1)        _real  [m^-3]  #elec density below y-face center
ney1(0:nx+1,0:ny+1)        _real  [m^-3]  #elec density above y-face center
nity0(0:nx+1,0:ny+1)       _real  [m^-3]  #total ion density below y-face center
nity1(0:nx+1,0:ny+1)       _real  [m^-3]  #total ion density above y-face center
tey0(0:nx+1,0:ny+1)        _real  [eV]    #elec temp below y-face center
tey1(0:nx+1,0:ny+1)        _real  [eV]    #elec temp above y-face center
tiy0(0:nx+1,0:ny+1)        _real  [eV]    #ion temp below y-face center
tiy1(0:nx+1,0:ny+1)        _real  [eV]    #ion temp above y-face center
tiy0s(0:nx+1,0:ny+1)       _real  [eV]    #old ion temp below y-face center
tiy1s(0:nx+1,0:ny+1)       _real  [eV]    #old ion temp above y-face center
tgy0(0:nx+1,0:ny+1,1:ngsp) _real  [eV]    #atom temp below y-face center
tgy1(0:nx+1,0:ny+1,1:ngsp) _real  [eV]    #atom temp above y-face center
ngy0(0:nx+1,0:ny+1,1:ngsp) _real  [m^-3]  #gas density below y-face center
ngy1(0:nx+1,0:ny+1,1:ngsp) _real  [m^-3]  #gas density above y-face center
pgy0(0:nx+1,0:ny+1,1:ngsp) _real  [J/m^3] #gas pressure below y-face center
pgy1(0:nx+1,0:ny+1,1:ngsp) _real  [J/m^3] #gas pressure above y-face center
pg(0:nx+1,0:ny+1,1:ngsp)   _real  [J/m^3] #gas pressure at cell center
phiy0(0:nx+1,0:ny+1)       _real  [V]     #potential below y-face center
phiy1(0:nx+1,0:ny+1)       _real  [V]     #potential above y-face center
phiy0s(0:nx+1,0:ny+1)      _real  [V]     #old potential below y-face center
phiy1s(0:nx+1,0:ny+1)      _real  [V]     #old potential above y-face center
pr(0:nx+1,0:ny+1)          _real  [J/m^3] #total pressure at center of cell
prev(0:nx+1,0:ny+1)        _real  [J/m^3] #elec pressure at vertex of cell
prtv(0:nx+1,0:ny+1)        _real  [J/m^3] #total pressure at vertex of cell
pri(0:nx+1,0:ny+1,1:nisp)  _real  [J/m^3] #ion plasma pressure
priv(0:nx+1,0:ny+1,1:nisp) _real  [J/m^3] #ion pressure at vertex of cells
priy0(0:nx+1,0:ny+1,1:nisp) _real [J/m^3] #ion pressure below y-face center
priy1(0:nx+1,0:ny+1,1:nisp) _real [J/m^3] #ion pressure above y-face center
pre(0:nx+1,0:ny+1)         _real  [J/m^3] #el. plasma pressure
ne(0:nx+1,0:ny+1)          _real  [m^-3]  #electron dens in primary cell (ix,iy)
nit(0:nx+1,0:ny+1)         _real  [m^-3]  #tot ion dens in primary cell (ix,iy)
nginit(0:nx+1,0:ny+1)      _real  [m^-3]  #init gas dens in primary cell (ix,iy)
phi(0:nx+1,0:ny+1)         _real  [V]     #potential in primary cell (ix,iy)
phiv(0:nx+1,0:ny+1)        _real  [V]     #potential at vertex of cell
zeff(0:nx+1,0:ny+1)        _real  [ ]     #Z_effective charge in cell (ix,iy)
loglambda(0:nx+1,0:ny+1)   _real  [ ]     #Coulomb logarithm on "east" x-face
netap(0:nx+1,0:ny+1)       _real  [ ]     #ne*parallel resistivity
znot(0:nx+1,0:ny+1)        _real  [ ]     #=Sum(n_z * Z^2)/n_i in cell
zimpc(0:nx+1,0:ny+1)       _real  [ ]     #Zimp (avg-ion model) in cell (ix,iy)
nil(0:nx+1,0:ny+1,1:nisp)  _real  [m^-3]  #ion density at last output
upl(0:nx+1,0:ny+1,1:nisp)  _real  [m/s]   #parallel ion velocity at last output
tel(0:nx+1,0:ny+1)         _real  [J]	  #electron temperature at last output
til(0:nx+1,0:ny+1)         _real  [J]	  #ion temperature at last output
ngl(0:nx+1,0:ny+1,1:ngsp)  _real  [m^-3]  #gas density at last output
phil(0:nx+1,0:ny+1)        _real  [V]     #potential at last output
upxpt(1:nusp,1:nxpt)       _real  [m/s]   #parallel velocity at x-point
nixpt(1:nusp,1:nxpt)       _real  [m^-3]  #ion density at x-point
visyxpt(1:nusp,1:nxpt)     _real          #ion viscosity at x-point
vyhxpt(1:nusp,1:nxpt)      _real  [m/s]   #horiz. ion drift vel. at x-point
vyvxpt(1:nusp,1:nxpt)      _real  [m/s]   #vert. ion drift vel. at x-point
fmihxpt(1:nusp,1:nxpt)     _real  [Nwt]   #horiz. mom. flux at x-point
fmivxpt(1:nusp,1:nxpt)     _real  [Nwt]   #vert. mom. flux at x-point
rtauxfac                    real  /0./    #fac*rtaux, Ly-a optic depth to plate
					  #=1 standard; <=0 skips rtau calc.
rtauyfac                    real  /1./    #fac*rtauy, Ly-a optic depth to wall
rt_scal			    real  /1.e-16/#factor to scale rtaux,y & thus rtau
rtaux(0:nx+1,0:ny+1) _real [1e-16 m^-2]/0./ #Norm. poloidal neutral line-dens.,
					    #Ly-a opacity to plates
rtauy(0:nx+1,0:ny+1) _real [1e-16 m^-2]/0./ #Norm. radial neutral line-dens.,
					    #norm. Ly-a opacity to radial wall
rtau(0:nx+1,0:ny+1)  _real [1e-16 m^-2]/0./ #Min. norm neutral line-dens.,
					    #min. Ly-a  opacity; min(rtaux,rtauy)
betap(0:nx+1,0:ny+1) _real        /0./      #poloidal plasma beta
fracvgpgp             real        /1./      #frac of vgp in vgradp eng terms

***** Postproc:
#Variables used in the postprocessing of data to check energy and particle bal.
fetx(0:nx+1,0:ny+1)   _real [W]  #total energy flow through a poloidal cell face
fety(0:nx+1,0:ny+1)   _real [W]  #total energy flow through a radial cell face
pdrift(0:nx+1,0:ny+1) _real [W/m^3] #power in bringing new ion to flow velocity
peirad(0:nx+1,0:ny+1) _real [W]     #tot. power lost by electrons and ions in
                                    #rad., ion. and dissoc.
png2ni(0:nx+1,0:ny+1) _real [W]     #power exchange bwt. neutral and ion
pmomv(0:nx+1,0:ny+1)  _real [W/m^3] #power exchange bwt. neutal & ion from flow
jdote(0:nx+1,0:ny+1)  _real [W]     #power from J.E heating
engerr(0:nx+1,0:ny+1) _real         #local error in power balance
ptjdote                real [W]     #sum of J.E heating
ptigas                 real [W/m^3] #sum of png2ni heating
pvmomcx                real [W/m^3] #sum of pmomv heating
iion(ngsp)            _real [A]     #net ionization current per gas isotope
irecomb(ngsp)         _real [A]     #net recombination current of gas isotope
icxgas(ngsp)          _real [A]     #net charge exchange current of gas isotope
iion_tot               real [A]     #total ioniz. current over all isotopes
irecomb_tot            real [A]     #total recomb. current over all isotopes
icxgas_tot             real [A]     #total cx current over all gas isotopes
pradht                 real [W]     #total H photon rad. loss (- binding eng.)
pradiz                 real [W]     #ionization radiation energy loss
pradrc                 real [W]     #recombination radiation energy loss
pradimpt(ngsp)        _real [W]     #net impurity photon rad. loss
pradfft                real [W]     #net radiation loss via fixed-fraction impurity
pradzbind              real [W]     #elec loss at imp ioniz carried as bind eng
pradimp(0:nzspmx,ngsp-1) _real [W]  #rad. loss for each impurity charge state
pbinde                 real [W]     #pwr stored in ion binding pot. eng. from ioniz.
pbindrc                real [W]     #binding eng. pwr released to elec. from recomb.
prdiss                 real [W]     #net photon pwr lost in dissoc.
pibirth                real [W]     #net ion energy gain from dissoc.
pwr_plth(0:ny+1,2*nxpt) _real [W/m**2]#hydrog rad pwr flux on divertor plate
pwr_pltz(0:ny+1,2*nxpt) _real [W/m**2]#impur rad pwr flux on divertor plate
pwr_wallh(0:nx+1)       _real [W/m**2]#hydrog rad pwr flux on outer wall
pwr_wallz(0:nx+1)       _real [W/m**2]#impur rad pwr flux on outer wall
pwr_pfwallh(0:nx+1,nxpt) _real [W/m**2]#hydrog rad pwr flux on PF wall
pwr_pfwallz(0:nx+1,nxpt) _real [W/m**2]#impur rad pwr flux on PF wall
sdelb(0:ny+1,1:nxpt)	_real [W/m**2]#elec pwr flux to left div
sderb(0:ny+1,1:nxpt)	_real [W/m**2]#elec ion pwr flux to right div
sdilb(0:ny+1,1:nxpt)	_real [W/m**2]#tot ion pwr flux to left div
sdirb(0:ny+1,1:nxpt)	_real [W/m**2]#tot ion pwr flux to right div
sbindlb(0:ny+1,1:nxpt)	_real [W/m**2]#tot bind eng pwr flux to left div
sbindrb(0:ny+1,1:nxpt)	_real [W/m**2]#tot bind eng pwr flux to right div
sdrlb(0:ny+1,1:nxpt)	_real [W/m**2]#tot rad pwr flux to left div
sdrrb(0:ny+1,1:nxpt)	_real [W/m**2]#tot rad flux to right div
sdtlb(0:ny+1,1:nxpt)	_real [W/m**2]#tot pwr flux to left div
sdtrb(0:ny+1,1:nxpt)	_real [W/m**2]#tot pwr flux to right div
gdilb(0:ny+1,nisp,nxpt) _real [1/m**2s]#particle flux to left div
gdirb(0:ny+1,nisp,nxpt) _real [1/m**2s]#particle flux to right div
engilb(0:ny+1,nisp,nxpt) _real [Volts] #ave ion energy to left div
engirb(0:ny+1,nisp,nxpt) _real [Volts] #ave ion energy to right div
gwalli(0:nx+1,nisp)     _real [1/m**2s]#particle flux to left div
engwalli(0:nx+1,nisp)   _real [Volts] #ave ion energy to outer wall
swallr(0:nx+1)		_real [W/m**2]#radiation pwr flux to outer wall
swalli(0:nx+1)		_real [W/m**2]#ion pwr flux to outer wall
swalle(0:nx+1)		_real [W/m**2]#elec pwr flux to outer wall
swbind(0:nx+1)		_real [W/m**2]#binding energy flux to outer wall
swallt(0:nx+1)		_real [W/m**2]#total pwr flux to outer wall
spfwallr(0:nx+1,nxpt)   _real [W/m**2]#radiation pwr flux to PF wall

***** Volsrc:
#Parameters for volume particle, mom. & power sources for electrons and ions
pwrsore(0:nx+1,0:ny+1) _real [W]       +input #power src into electrons in cell ix,iy
pwrsori(0:nx+1,0:ny+1) _real [W]       +input #power src into ions in cell ix,iy
volpsor(0:nx+1,0:ny+1,1:nisp) _real [1/s]+input #current src into ions in cell ix,iy
volmsor(0:nx+1,0:ny+1,1:nisp) _real [kg m/s**2] +input #up mom src in cell ix,iy
voljcsor(0:nx+1,0:ny+1) _real [A]      +input #uniform core-region curr sor. in ix,iy
volpsorg(0:nx+1,0:ny+1,1:ngsp) _real [1/s]+input #curr source for gas in cell ix,iy
pondpot(0:nx+1,0:ny+1)  _real [V] /0./ +input #elec ponderomotive potential
psgov_use(0:nx+1,0:ny+1,1:ngsp) _real [1/m**3 s]+input #user-specified gas source
jcvsor                  real [A]  /0./ +input #total core-region current for voljcsor
ix_sjcsor	        integer   /0/  +input #if nonzero, beginning ix for voljcsor
ix_ejcsor	        integer   /0/  +input #if nonzero, ending ix for voljcsor
iy_sjcsor	        integer   /0/  +input #if nonzero, beginning iy for voljcsor
iy_ejcsor	        integer   /0/  +input #if nonzero, ending iy for voljcsor
pvole                   real [W]  /0./ +input #total power into electrons
pvoli                   real [W]  /0./ +input #total power into ions
z0pe                    real [m]  /0./ +input #axial or x loc. of elec. power profile
z0pi                    real [m]  /0./ +input #axial or x loc. of ion power profile
r0pe                    real [m]  /0./ +input #radial or y loc. of elec. power profile
r0pi                    real [m]  /0./ +input #radial or y loc. of ion power profile
zwpe                    real [m] /3./  +input #axial or x Gauss. 1/2 width of e-power
zwpi                    real [m] /3./  +input #axial or y Gaussian 1/2 width ion power
rwpe                    real [m] /.05/ +input #rad. or x Gaussian 1/2 width e-power
rwpi                    real [m] /.05/ +input #rad. or y Gaussian 1/2 width ion power
ivolcur(1:nisp)        _real [A]  /0./ +input #total volume current
mvolcur(1:nisp)        _real [kgA m/s] /0./ +input #total volume parallel mom. curr.
z0ni                    real [m]  /0./ +input #axial or x loc. of ion particle profile
r0ni                    real [m]  /0./ +input #rad. or y loc. of ion particle profile
zwni                    real [m] /3./  +input #axial or y Gaussian 1/2 width ion prtcl
rwni                    real [m] /.05/ +input #rad. or y Gaussian 1/2 width ion prtcl
z0up                    real [m]  /0./ +input #axial or x loc. of ion mom. profile
r0up                    real [m]  /0./ +input #rad. or y loc. of ion mom. profile
zwup                    real [m] /3./  +input #axial or y Gaussian 1/2 width ion mom.
rwup                    real [m] /.05/ +input #rad. or y Gaussian 1/2 width ion mom.
ponderompot             real [V] /0./  +input #peak elec ponderomotive potential
z0pondp                 real [m]  /0./ +input #axial or x loc. of ion mom. profile
r0pondp                 real [m]  /0./ +input #rad. or y loc. of ion mom. profile
zwpondp                 real [m] /3./  +input #axial or y Gaussian 1/2 width ion mom.
rwpondp                 real [m] /.05/ +input #rad. or y Gaussian 1/2 width ion mom.
thetarot		real [rad]/0./ +input #rotation angle for R,Z with effec. R,Z
		                       # R_e= R0+(R-R0)cos(th)+(Z-Z0)sin(th),
		                       # Z_e= Z0-(R-R0)sin(th)+(Z-Z0)cos(th),
rcutmin			real [m] /0./  +input #source zero if R<rcutmin
zcutmin			real [m] /0./ +input  #source zero if Z<zcutmin
effvng(1:ngsp)	       _real [m] /0./  # normalizing factor of gas source; calc
ivolcurg(1:ngsp)       _real [A] /0./  +input #tot. volumn gas source strength
z0ng(1:ngsp)	       _real [m] /0./  +input #axial or x loc. of gas particle profile
r0ng(1:ngsp)           _real [m] /0./  +input #rad. or y loc. of gas particle profile
zwng(1:ngsp)           _real [m] /3./  +input #axial or y Gaussian 1/2 width gas prtcl
rwng(1:ngsp)           _real [m] /.05/ +input #rad. or y Gaussian 1/2 width gas prtcl

***** Bfield:
#Variables for the B-field and grad_B drift geo/B-field factors
b0                     real /1./ +input
                                 #scale factor for magnetic field
b02d(0:nx+1,0:ny+1)   _real      #net B-field scale fac. =b0+b0_use present iter.
b0old(0:nx+1,0:ny+1)  _real      #net B-field scale factor at last iter.
b0_use(0:nx+1,0:ny+1) _real      #spatial B-field scale factor; user input
rbpol(0:nx+1,0:ny+1)  _real [mT]    #major radius*poloidal magnetic field
btot(0:nx+1,0:ny+1)   _real [T]     #total magnetic field strength
rbfbt(0:nx+1,0:ny+1)  _real [ ]     #ratio bphi/btot at density cell center
rbfbt2(0:nx+1,0:ny+1) _real [1/T]   #ratio bphi/btot**2 at density cell center
curvrby(0:nx+1,0:ny+1) _real [1/mT] #curvature drift factor on y-face
curvrb2(0:nx+1,0:ny+1) _real [1/mT] #curvature drift factor on x-face
gradby(0:nx+1,0:ny+1)  _real [1/mT] #grad_B drift (p_perp) drift factor, y-face
gradb2(0:nx+1,0:ny+1)  _real [1/mT] #grad_B drift (p_perp) drift factor, x-face
dbm2dx(0:nx+1,0:ny+1) _real [1/T**2m] #pol deriv of 1/B**2 on x-face
dbm2dy(0:nx+1,0:ny+1) _real [1/T**2m] #rad deriv of 1/B**2 on y-face
isrozhfac            integer  /0/   +maybeinput #=0 sets bfacx,yrozh=1; =1 computes ave
bfacxrozh(0:nx+1,0:ny+1) _real [ ]  #[1-B**2/B**2_ave], x-face; Rozhansky
bfacyrozh(0:nx+1,0:ny+1) _real [ ]  #[1-B**2/B**2_ave], y-face; Rozhansky

***** Oldpla:
#Old value of some of the plasma quantities.
# ISSUE: use nx, ny or nx,ny; latter works but ???
ni0(0:nx+1,0:ny+1,1:nisp)  _real [m^-3]  #old ion density
ng0(0:nx+1,0:ny+1,1:ngsp)  _real [m^-3]  #old neutral density
te0(0:nx+1,0:ny+1)         _real [J]     #old electron temperature
ti0(0:nx+1,0:ny+1)         _real [J]     #old ion temperature
phi0(0:nx+1,0:ny+1)        _real [V]     #old electrostatic potential
up0(0:nx+1,0:ny+1,1:nisp)  _real [m/s]   #old parallel velocity
vy0(0:nx+1,0:ny+1,1:nisp)  _real [m/s]   #old radial velocity

***** Comflo:
#Variables in common -- flows
fqp(0:nx+1,0:ny+1)         _real [Amp]  #pol proj of par cur, east face
cfparcur                    real /0./   +input #scale fac fqp=cfparcur*parcurrent if	
                                        #isimpon=5 (fmombal from Hirshman)
fq2(0:nx+1,0:ny+1)         _real [Amp]  #pol proj of 2 cur, east face
fqx(0:nx+1,0:ny+1)         _real [Amp]  #net poloidal current, east face
fqxb(0:nx+1,0:ny+1)        _real [Amp]  #poloidal cur from grad_B, east face
fdiaxlb(0:ny+1,1:nxpt)     _real [Amp]  #left boundary Dia current for bc
fdiaxrb(0:ny+1,1:nxpt)     _real [Amp]  #right boundary Dia current for bc
floxebgt(0:nx+1,0:ny+1)    _real [W]    #BxgradTe diamag part floxe (-> feex)
floxibgt(0:nx+1,0:ny+1,1:nisp) _real [W]#BxgradTi diamag part floxi (-> feex)
fqy(0:nx+1,0:ny+1)         _real [Amp]  #net radial current, north face
fqyb(0:nx+1,0:ny+1)        _real [Amp]  #radial current from grad_B, north face
fqyn(0:nx+1,0:ny+1)        _real [Amp]  #radial cur from cx coll, north face
fqym(0:nx+1,0:ny+1)        _real [Amp]  #radial cur from inertia, north face
fqymi(0:nx+1,0:ny+1,1:nisp) _real [Amp] #spec rad cur from inertia, north face
fqya(0:nx+1,0:ny+1)        _real [Amp]  #anomalous visc rad cur, north face
fqydt(0:nx+1,0:ny+1)       _real [Amp]  #time-dep inertial rad cur, north face
fqydti(0:nx+1,0:ny+1,1:nisp) _real [Amp]#spec time-dep inert rad cur, north face
fqyao(0:nx+1,0:ny+1)       _real [Amp]  #old anom mobil rad current, north face
fqyae(0:nx+1,0:ny+1)       _real [Amp]  #anom mobil rad current for electrons, north face
fqyai(0:nx+1,0:ny+1)       _real [Amp]  #anom mobil rad current for ions, north face
fqyd(0:nx+1,0:ny+1)        _real [Amp]  #diamag radial current; north face
fqygp(0:nx+1,0:ny+1)       _real [Amp]  #net radial curr. uses grad_P, north face
fq2d(0:nx+1,0:ny+1)        _real [Amp]  #diamag 2-current; east face
fqypneo(0:nx+1,0:ny+1)     _real [Amp]  #rad-cur from neo particle flux
fq2pneo(0:nx+1,0:ny+1)     _real [Amp]  #2-cur from neo particle flux
fqyqneo(0:nx+1,0:ny+1)     _real [Amp]  #rad-cur from neo heat flux
fq2qneo(0:nx+1,0:ny+1)     _real [Amp]  #2-cur from neo heat flux
fnix(0:nx+1,0:ny+1,1:nisp) _real [1/s]  #ion poloidal current, east face
fnixcb(0:nx+1,0:ny+1,1:nisp) _real [1/s]  #ion grad-B pol. current, east face
fniy(0:nx+1,0:ny+1,1:nisp) _real [1/s]  #ion radial current, north face
fniy4ord(0:nx+1,0:ny+1,1:nisp) _real [1/s] #4th ord ion radial current, north face
fniycb(0:nx+1,0:ny+1,1:nisp) _real [1/s]  #ion  grad-B rad. current, north face
flnix(0:nx+1,0:ny+1,1:nisp) _real [1/s] #ion poloidal log-current, east face
flniy(0:nx+1,0:ny+1,1:nisp) _real [1/s] #ion radial log-current, north face
fmix(0:nx+1,0:ny+1,1:nusp) _real [Nwt]  #ion poloidal momentum current,east face
fmiy(0:nx+1,0:ny+1,1:nusp) _real [Nwt]  #ion radial momentum current, north face
fmixy(0:nx+1,0:ny+1,1:nusp) _real [Nwt] #nonorthog ion pol. mom. curr., east f.
fmity(0:nx+1,0:ny+1,1:nisp) _real [ ]   #rad flux of cross-field tor. mom*R/Bp;
                                        #nisp dimen, not nusp as for pot eqn
fmgx(0:nx+1,0:ny+1,ngsp)   _real [Nwt]  #pol. neutral mom. current, east face ### IJ 2016/10/11
fmgy(0:nx+1,0:ny+1,ngsp)   _real [Nwt]  #rad. neutral mom. current, north face  ### IJ 2016/10/11
feex(0:nx+1,0:ny+1)        _real [J/s]  #poloidal electron thermal current,
                                         #east face
feey(0:nx+1,0:ny+1)        _real [J/s]  #radial electron thermal current,
                                         #north face
feexy(0:nx+1,0:ny+1)       _real [J/s]  #nonorthog elec. pol. therm cur, east f.
feey4ord(0:nx+1,0:ny+1)    _real [J/s]  #elec. pol. kye4order therm cur, east f.
feix(0:nx+1,0:ny+1)        _real [J/s]  #poloidal ion thermal current, east face
feiy(0:nx+1,0:ny+1)        _real [J/s]  #radial ion thermal current, north face
fegx(0:nx+1,0:ny+1,ngsp)   _real [J/s]  #poloidal neut thermal curr, east face ### IJ 2016/09/2
fegy(0:nx+1,0:ny+1,ngsp)   _real [J/s]  #radial neut thermal curr, north face  ### IJ 2016/09/22
fegxy(0:nx+1,0:ny+1,ngsp)  _real [J/s]  #pol. nonog neut thermal curr, north face 
isfegxyqflave            integer  /0/   +input #=0fegxy T*vt,ng ave;=1, use harm aves
cfegxy                      real  /1./  +input #coeff multiple fegxy
qipar(0:nx+1,0:ny+1,nisp)  _real [J/m**2s] #parallel conductive ion heat flux
qgpar(0:nx+1,0:ny+1,ngsp)  _real [J/m**2s] #parallel conductive gas heat flux
fniycbo(0:nx+1,1:nisp)     _real [1/s]  #fniy cor. iy=0 bdry for grad_B, grad_P
feiycbo(0:nx+1)            _real [J/s]  #feiy cor. iy=0 bdry for grad_B, grad_P
feeycbo(0:nx+1)            _real [J/s]  #feey cor. iy=0 bdry for grad_B, grad_P
feixy(0:nx+1,0:ny+1)       _real [J/s]  #nonorthog ion pol. thermal cur, east f.
feiy4ord(0:nx+1,0:ny+1)    _real [J/s]  #ion pol. kyi4order therm cur, east f.
fngx(0:nx+1,0:ny+1,1:ngsp)  _real [1/s] #neutral polodial current, east face
fngx4ord(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #4th ord gas radial current, north face
flngx(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #neutral pol. log-current, east face
fngxs(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #neutral pol cur w/o fngxy, east face
fngy(0:nx+1,0:ny+1,1:ngsp)  _real [1/s] #neutral radial current, north face
fngy4ord(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #4th ord gas radial current, north face
flngy(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #neutral radial log-current, north face
fngxy(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #nonorthog gas pol. cur., east face
flngxy(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #nonorthog gas pol.log-cur., east face
fngyx(0:nx+1,0:ny+1,1:ngsp) _real [1/s] #nonorthog gas rad. cur., north face
fnixtot(0:nx+1,0:ny+1)      _real [1/s] #total poloidal ion cur.
fniytot(0:nx+1,0:ny+1)      _real [1/s] #total radial ion cur.

***** Indexes:
#Indices that help the calculation
idxn(0:nx+1,0:ny+1,1:nisp)	_integer  # index of yl vector for ni(ix,iy)
idxg(0:nx+1,0:ny+1,1:ngsp)	_integer  # index of yl vector for ng(ix,iy,ig)
idxtg(0:nx+1,0:ny+1,1:ngsp)	_integer  # index of yl vector for tg(ix,iy,ig)
idxu(0:nx+1,0:ny+1,1:nusp)	_integer  # index of yl vector for up(ix,iy)
idxti(0:nx+1,0:ny+1)		_integer  # index of yl vector for ti(ix,iy)
idxte(0:nx+1,0:ny+1)		_integer  # index of yl vector for te(ix,iy)
idxphi(0:nx+1,0:ny+1)		_integer  # index of yl vector for phi(ix,iy)
ivfirst(0:nx+1,0:ny+1)		_integer  # first eqn number iv at (ix,iy)
igyl(neqmx,2)			_integer  # ix,iy indices for vector yl(iv)
iseqalg(neqmx)                  _integer  # flag(=1) for eqn being algebraic
isvarup(numvar)			_integer  # flag(=1) for variable being up
isvarphi(numvar)		_integer  # flag(=1) for variable being phi

***** Stat:
#Variables for statistics of solver performance
hu(nsteps,ngrid)   _real    [1/s]  #present timestep
gpe(nsteps,ngrid)  _real           #ratio, linear to nonlinear iter, nli/nni
npe(nsteps,ngrid)  _integer        #cumulative number of precond. eval.
nps(nsteps,ngrid)  _integer        #cumulative number of precond. solves
nfe(nsteps,ngrid)  _integer        #cumulative number of RHS eval.
nst(nsteps,ngrid)  _integer        #cum. number of steps taken
nni(nsteps,ngrid)  _integer        #cumulative number of nonlinear iter.
nli(nsteps,ngrid)  _integer        #cumulative number of linear iter.
nje(nsteps,ngrid)  _integer        #cumulative number of Jacobian eval.
ncfn(nsteps,ngrid) _integer        #number of nonlinear converg. failures
ncfl(nsteps,ngrid) _integer        #number of linear converg. failures
nqu(nsteps,ngrid)  _integer        #method order last used
iddas(nsteps,ngrid)  _integer      #idid for daspk
eqmxer(nsteps,ngrid) _integer      #eqn. number giving maximum error in lsode
lacor               integer        #location in rwork where error vector begins
lewt                integer        #location in rwork where ewt**-1 begins
npsn(ngrid)        _integer        #cum. number of Jacobian solves for Newton
njen(ngrid)        _integer        #cum. Jacobian evals. for Newton iter.

***** Poten:
#Variables required for the calculation of the potential.
newbcl(nxptmx)     integer /nxptmx*0/ +input #switch on new sheath model
newbcr(nxptmx)     integer /nxptmx*0/ +input #switch on new sheath model
iskaplex           integer /0/  +input #=1 if kappal is set externally (from parser)
iskaprex           integer /0/  +input #=1 if kappar is set externally (from parser)
bcee               real  /4./ +input
                                #electron sheath energy trans. factor(newbc=0)
bcei               real  /2.5/ +input
                                #ion sheath energy trans. factor(newbc=0)
bceew              real  /4./   +input #elec wall energy trans factor
bceiw              real  /2.5/  +input #ion wall energy trans factor
bcen               real  /0./   +input #neut energy trans. factor on plates 
                                #For combined neutral+ion energy equation
bcenw              real  /0./   +input #neut eng trans fac on walls
isfdiax            real  /0./   +input #switch to turn on diamagnetic drift for sheath
				#potential calculation
cthe               real  /0.71/ +input #electron thermal force coeff.
cthi               real  /4.05/ +input #?ion thermal force coeff.?; enters if zeff.ne.1
sigma1             real [1/(eV**1.5_Ohm_m)] /1490./+input #parallel conductivity coeff.
                        #sigma1=1/(5.19e-5*Z*ln_lambda) & Z=1, ln_lambda=12.9
cfsigm		   real   /1./  +input #scale factor for parallel cond. sigma1
rsigpl             real  /0./   +input #ad hoc radial electrical conductivity - global
rsigplcore         real  /0./   +input #ad hoc radial electrical conduct - core only
                                #ratio of perp to parallel conductivity
bcel(0:ny+1,nxpt)      _real  [ ]   +maybeinput #electron sheath energy transmission factor
                                    #on the left boundary
bcer(0:ny+1,nxpt)      _real  [ ]   +maybeinput #electron sheath energy transmission factor
                                    #on the right boundary
bcil(0:ny+1,nxpt)      _real  [ ]   +maybeinput #ion sheath energy transmission factor
                                    #on the left boundary
bcir(0:ny+1,nxpt)      _real  [ ]   +maybeinput #ion sheath energy transmission factor
                                    #on the right boundary
kappal(0:ny+1,nxpt)  _real  [ ]	+maybeinput #sheath pot'l drop on left  boundary, phi/Te
kappar(0:ny+1,nxpt)  _real  [ ]	+maybeinput #sheath pot'l drop on right boundary, phi/Te
kappallim(0:ny+1)    _real  [ ]	#sheath pot'l drop on left-side limiter, phi/Te
kapparlim(0:ny+1)    _real  [ ]	#sheath pot'l drop on right-side limiter, phi/Te
bctype(0:ny+1)    _integer #/0,ny*0,0/+maybeinput 
phi0r(0:ny+1,nxpt)	_real [V] /0./ +maybeinput #plate pot'l at right poloidal boundary
phi0l(0:ny+1,nxpt)	_real [V] /0./ +maybeinput #plate pot'l at left  poloidal boundary
capx(1:ny)        _real    #/ny*0.0/+maybeinput 
dphi_iy1(0:nx+1)  _real [V] #/(nx+2)*0./  +maybeinput #incremental phi at iy=1 to have
                                          #Te=constant for second phi BC
kincorlb(0:ny+1,nxpt)   _real [ ]   +maybeinput # kinetic corr. factor for elec part. loss, left b
kincorrb(0:ny+1,nxpt)   _real [ ]   +maybeinput # kinetic corr. factor for elec part. loss, right b
cfkincor            real     [ ] /0.5/ +input # factor for kincorlb,rb denom. factor
#Variables for the grid-sequencing.
#yet to be defined?

***** Gradients:
#Gradients of the different physical quantities.
ex(0:nx+1,0:ny+1)           _real  [V/m]  #poloidal electric field
ey(0:nx+1,0:ny+1)           _real  [V/m]  #radial electric field
eymask1d(0:nx+1,0:ny+1)     _real  [V/m]  #set ey=0 in core+sep if isphicore0=1
einduc			     real  [V/m]  +input #inductive tor. E-field - input
gpix(0:nx+1,0:ny+1,1:nisp)  _real  [Pa/m] #X-gradient of ion pressure
gpiy(0:nx+1,0:ny+1,1:nisp)  _real  [Pa/m] #Y-gradient of ion pressure
gpex(0:nx+1,0:ny+1)         _real  [Pa/m] #X-gradient of el. pressure
gpey(0:nx+1,0:ny+1)         _real  [Pa/m] #Y-gradient of el. pressure
gprx(0:nx+1,0:ny+1)         _real  [Pa/m] #X-gradient of total pressure
gpry(0:nx+1,0:ny+1)         _real  [Pa/m] #Y-gradient of total pressure
gtex(0:nx+1,0:ny+1)         _real  [J/m]  #X-gradient of el. temperature
gtey(0:nx+1,0:ny+1)         _real  [J/m]  #Y-gradient of el. temperature
gtix(0:nx+1,0:ny+1)         _real  [J/m]  #X-gradient of ion temperature
gtiy(0:nx+1,0:ny+1)         _real  [J/m]  #Y-gradient of ion temperature
gpondpotx(0:nx+1,0:ny+1)    _real  [V/m]  #X-gradient of elec pondom pot

***** Cfric:
#Coulomb friction terms for parallel transport
frice(0:nx+1,0:ny+1)      _real [J/m**4]  +maybeinput #Electron parallel Coulomb friction
frici(0:nx+1,0:ny+1,nisp) _real [J/m**4]  +maybeinput #Ion parallel Coulomb friction
fricnrl(0:nx+1,0:ny+1,nusp) _real [J/m**4] +maybeinput #NRL ion par fric ni*mi*nu*(up1-up2)
cfgti             /1./     real           +input #scale factor for ion thermal force
cfgte             /1./     real           +input #scale factor for elec. thermal force
cftaud            /1./     real           +input #scale factor for ion-ion drag time
isalfecalc(1:nisp) /1/    _integer        +input #=1 for internal calc of alfe
isbetaicalc(1:nisp)/1/    _integer        +input #=1 for internal calc of betai
alfe(1:nisp)       /1./   _real         +input #grad_Te thm force coeff isalfecalc=0
betai(1:nisp)      /1./   _real         +input #grad_Ti thm force coeff isbetaicalc=0

***** Grid:
ngrid          /1/	 integer +regrid
ig             /1/       integer  #counter for mesh-seq number
inewton(30)   /30*0/     integer  #=1 for Newton iter., =0 for time-dependent
                                  #reset=1 internally if svrpkg=nksol or newton
imeth            /0/     integer  #imeth=inewton(igrid)
nurlx           /1.e8/   real    [1/s] #rate coeff. to relax to boundary cond.
ijac(ngrid)             _integer
ijactot          /0/     integer  # tot Jac calcs, used as check when icntnunk=1

***** Wkspace:
#Workspace arrays
w(0:nx+1,0:ny+1)        _real
w0(0:nx+1,0:ny+1)       _real
w1(0:nx+1,0:ny+1)       _real
w2(0:nx+1,0:ny+1)       _real
w3(0:nx+1,0:ny+1)       _real

***** Locflux:
#Local arrays for the calculation of the fluxes and other quantities.
flox(0:nx+1,0:ny+1)     _real
floy(0:nx+1,0:ny+1)     _real
conx(0:nx+1,0:ny+1)     _real
cony(0:nx+1,0:ny+1)     _real
floxe(0:nx+1,0:ny+1)    _real
floye(0:nx+1,0:ny+1)    _real
floxi(0:nx+1,0:ny+1)    _real
floyi(0:nx+1,0:ny+1)    _real
floxg(0:nx+1,0:ny+1)    _real
floyg(0:nx+1,0:ny+1)    _real
fgtdx(0:nx+1)		_real	#scale factor for gas grad-x T vel
fgtdy(0:ny+1)		_real	#scale factor for gas grad-x T vel
conxe(0:nx+1,0:ny+1)    _real
conye(0:nx+1,0:ny+1)    _real
conxi(0:nx+1,0:ny+1)    _real
conyi(0:nx+1,0:ny+1)    _real
conxg(0:nx+1,0:ny+1)    _real
conyg(0:nx+1,0:ny+1)    _real
floxge(0:nx+1,0:ny+1,1:ngsp) _real
floyge(0:nx+1,0:ny+1,1:ngsp) _real
conxge(0:nx+1,0:ny+1,1:ngsp) _real
conyge(0:nx+1,0:ny+1,1:ngsp) _real

***** Conduc:
#Variables for the common -- conduc
visx(0:nx+1,0:ny+1,1:nisp)  _real [kg/m s]#poloidal viscosity coeff.
visy(0:nx+1,0:ny+1,1:nisp)  _real [kg/m s]#radial viscosity coeff.
hcxe(0:nx+1,0:ny+1)         _real [1/m s] #poloidal elec. therm. conduct.
hcye(0:nx+1,0:ny+1)         _real [1/m s] #radial elec. therm. conduct.
hcxij(0:nx+1,0:ny+1,1:nisp) _real [1/m s] #j-species pol. ion therm. conduct.
hcyij(0:nx+1,0:ny+1,1:nisp) _real [1/m s] #j-species rad. ion therm. conduct.
hcxg(0:nx+1,0:ny+1,1:ngsp)  _real [1/m s] #j-species pol. gas therm. conduct.
hcyg(0:nx+1,0:ny+1,1:ngsp)  _real [1/m s] #j-species rad. gas therm. conduct.
hcxi(0:nx+1,0:ny+1)         _real [1/m s] #summed pol. ion+neut therm. conduct.
hcxineo(0:nx+1,0:ny+1)      _real [1/m s] #neocl. pol. ion+neut therm. conduct.
hcyi(0:nx+1,0:ny+1)         _real [1/m s] #summed rad. ion+neut therm. conduct.
hcxn(0:nx+1,0:ny+1)         _real [1/m s] #poloidal neutral therm. conduct.
hcyn(0:nx+1,0:ny+1)         _real [1/m s] #radial neutral therm. conduct.
kxbohm(0:nx+1,0:ny+1)	    _real [m**2/s] +input #spatially depend. diff. on x-face
					  #set by user; Bohm if isbohmcalc=1
kybohm(0:nx+1,0:ny+1)	    _real [m**2/s] +input #spatially depend. diff. on y-face
					  #set by user; Bohm if isbohmcalc=1
vybohm(0:nx+1,0:ny+1)       _real [m/s]   +input #spatially depend. convect. y-vel
					  #set user if isbohmcalc=0; else =0
dif_use(0:nx+1,0:ny+1,1:nisp) _real [m**2/s] +input #spatially depend. diff; if
					  #isbohmcalc=1, user input if all
					  #facbni+facbup+facbee+facbei =0,
					  # or kybohm if facbni, etc. > 0;
					  # if isbohmcalc=2, then
					  # D = difni*kybohm/(difni+kybohm)
difp_use(0:nx+1,0:ny+1,1:nisp) _real [m**2/s]   +input #for gen pr diff; see dif_use comment
dif2_use(0:nx+1,0:ny+1,1:nisp) _real [m**2/s]   +input #for dif2; see dif_use comment
tray_use(0:nx+1,0:ny+1,1:nisp) _real [m**2/s]   +input #for travis; see dif_use comment
trax_use(0:nx+1,0:ny+1,1:nisp) _real [m**2/s]   +input #pol. analog to tra_use
kye_use(0:nx+1,0:ny+1)        _real [m**2/s]    +input #for kye; see dif_use comment
kyi_use(0:nx+1,0:ny+1)        _real [m**2/s]    +input #for kyi; see dif_use comment
kxe_use(0:nx+1,0:ny+1)        _real [m**2/s]    +input #user elec pol. heat cond
kxi_use(0:nx+1,0:ny+1)        _real [m**2/s]    +input #user ion pol. heat cond.
kxg_use(0:nx+1,0:ny+1,1:ngsp) _real [m**2/s]    +input #user gas pol. heat cond.
kyg_use(0:nx+1,0:ny+1,1:ngsp) _real [m**2/s]    +input #user gas rad. heat cond.
dutm_use(0:nx+1,0:ny+1,1:nisp) _real [m**2/s]   +input #for difutm; see dif_use comment
vy_use(0:nx+1,0:ny+1,1:nisp) _real [m/s]  +input #user-set rad vel;for isbohmcalc=0
vyup_use(0:nx+1,0:ny+1)     _real  [m/s]  +input #user-set conv vel of ion || vel, up
vyte_use(0:nx+1,0:ny+1)     _real  [m/s]  +input #user-set rad elec eng vel
vyti_use(0:nx+1,0:ny+1)     _real  [m/s]  +input #user-set rad ion eng vel
pondomfpare_use(0:nx+1,0:ny+1) _real [N/m**3] +input #user-in elec parallel pondero force
pondomfpari_use(0:nx+1,0:ny+1,1:nisp) _real [N/m**3] +input #user-in parallel ion ponder. force
fniyos_use(0:nx+1,0:ny+1,1:nisp) _real [1/s m**2] +input #user-set particle flux
feeyosn_use(0:nx+1,0:ny+1)  _real  [J/s m**2]    +input #user-set Te energy flux
feiyosn_use(0:nx+1,0:ny+1)  _real  [J/s m**2]    +input #user-set Ti energy flux
vy_cft(0:nx+1,0:ny+1,1:nisp) _real [m/s]  #calc vy from fniyos_use (fix flux)
vyte_cft(0:nx+1,0:ny+1)     _real  [m/s]  #calc vyte from feeyos_use (fix flux)
vyti_cft(0:nx+1,0:ny+1)     _real  [m/s]  #calc vyte from feiyos_use (fix flux)
nuiz(0:nx+1,0:ny+1,ngsp)    _real  [1/s]  #ionization rate (=ne*sigma*v)
nucx(0:nx+1,0:ny+1,ngsp)    _real  [1/s]  #charge-exchg rate for neut(sigv*ni)
nucxi(0:nx+1,0:ny+1,nisp)   _real  [1/s]  #charge-exchg rate for ion (sigv*ng)
nueli(0:nx+1,0:ny+1,nisp)   _real  [1/s]  #elast scatt rate for ion (sigv*ng)
nuelg(0:nx+1,0:ny+1,ngsp)   _real  [1/s]  #elast scatt rate for gas (sigv*nimp)
nuix(0:nx+1,0:ny+1,ngsp)    _real  [1/s]  #fnuizx*nuiz+fnucxx*nucx
fnuizx                       real    /0./ +input #fraction of nuiz in nuix (see nuix)
fnucxx                       real    /1./ +input #fraction of nucx in nuix (see nuix)
nurc(0:nx+1,0:ny+1,ngsp)    _real  [1/s]  #recombination rate
nuvl(0:nx+1,0:ny+1,nisp)    _real  [1/s]  #vol loss rate, ~cs/l_parloss for 1-D
cfvlh			     real         +input #scal fac for hyd rate in nuvl
cfvli(nisp)		    _real #/nisp*0./+input #scal fac for individ ion rate nuvl
l_parloss		     real [m] /1.e20/ +input #parall length for nuvl loss rate
eqp(0:nx+1,0:ny+1)          _real [1/m**3]#Te,i equipart. fact; needs *(Te-Ti)*vol
eqpg(0:nx+1,0:ny+1,ngsp)    _real [1/m**3]#Tg,i equipart. fact; needs *(Tg-Ti)*vol
                                          #..: modified to incorporate the separation of Tg and Ti.
engcoolm(0:nx+1,0:ny+1)     _real [J/s]   #cool rate ion/atoms by mols if ishymol=1
eeli(0:nx+1,0:ny+1)         _real  [J]    #electron energy loss per ionization
pradhyd(0:nx+1,0:ny+1)      _real [W/m**3] /0./#power radiated by hydrogen
tdiflim                      real [s] /0./ +input #lim on hcxe/ne; reduces hcxe if >0
lmfplim			     real [m] /1.e20/+input #hcxe,i -> hcxe,i/(1+lmfp/lmfelim)
eta1(0:nx+1,0:ny+1)	    _real [J-s/m**3] +maybeinput #Braginskii ion visc coeff eta_1
cfeta1                       real   /0./  +input # scale factor for eta1
rtaue(0:nx+1,0:ny+1)        _real [s/kg]  +maybeinput #Brag. R coeff (t_e/me)/(w_ce*t_e)**2
cfrtaue			     real   /0./  +input # scale factor for cfrtaue
dclass_e(0:nx+1,0:ny+1)     _real [m**2/s]#classical elec perp heat conduc.
dclass_i(0:nx+1,0:ny+1)     _real [m**2/s]#classical ion perp heat conduc.
cfcl_e	                     real  /0./   +input #scale fac for dclass_e
cfcl_i	                     real  /0./   +input #scale fac for dclass_i
omgci_taui		     real  /10./  +input #ion gy_freq*coll_rate for cl_model
omgce_taue		     real  /10./  +input #elec gy_freq*coll_rate for cl_model
nuneo			     real  /0./   +input #neoclass pol. damping rate; for fqyn
visxneo(0:nx+1,0:ny+1,1:nisp) _real [kg/m s] #Braginskii eta_0 neo-modified
visvol_v(0:nx+1,0:ny+1,1:nisp) _real #vel-based viscosity in (n*m*up)^dot eqn
visvol_q(0:nx+1,0:ny+1,1:nisp) _real #heat-flux-based viscosity (n*m*up)^dot eqn
nuii(0:nx+1,0:ny+1,1:nisp)    _real  #Braginski nuii coll freq.
nuiistar(0:nx+1,0:ny+1,1:nisp)_real  #neoclassical nuii coll freq.
alfneo(0:nx+1,0:ny+1,1:nisp)  _real  #neoclassical factor for q-based visc.
k2neo(0:nx+1,0:ny+1,1:nisp)   _real  #neoclassical coeff reducing therm cond
ktneo(0:nx+1,0:ny+1,1:nisp)   _real  #neoclassical coeff of grad Ti

***** Rhsides:
#Variables to evaluate the sources and RHS's.
snic(0:nx+1,0:ny+1,1:nisp)    _real
sniv(0:nx+1,0:ny+1,1:nisp)    _real
psorc(0:nx+1,0:ny+1,1:nisp)   _real  [part/s]  # cell ctr ioniz. sor plasma (>0)
psor(0:nx+1,0:ny+1,1:nisp)    _real  [part/s]  # cell ave ioniz. sor plasma (>0)
psort(0:nx+1,0:ny+1,1:nisp)   _real  [part/s]  # ioniz. source for plasma (>0)
psorxrc(0:nx+1,0:ny+1,1:nisp) _real  [part/s]  # cell ctr cx &recomb. for ions (<0)
psorxr(0:nx+1,0:ny+1,1:nisp)  _real  [part/s]  # cell ave cx &recomb. for ions (<0)
psor_tmpov(0:nx+1,0:ny+1)     _real  [part/s]  # work array for psor,etc for ave
psorgc(0:nx+1,0:ny+1,1:ngsp)  _real  [part/s]  # cell ctr ioniz. sor neutral (<0)
psorg(0:nx+1,0:ny+1,1:ngsp)   _real  [part/s]  # cell ave ioniz. sor neutral (<0)
psorrgc(0:nx+1,0:ny+1,1:ngsp) _real  [part/s]  # cell ctr recomb. source for neutrals
psorrg(0:nx+1,0:ny+1,1:ngsp)  _real  [part/s]  # cell ave recomb. source for neutrals
psorcxgc(0:nx+1,0:ny+1,1:ngsp) _real [part/s]  # cell ctr cx source for neutrals
psorcxg(0:nx+1,0:ny+1,1:ngsp) _real  [part/s]  # cell ave cx source for neutrals
psori(0:nx+1,0:ny+1,1:nisp)   _real  [part/s]  # impurity gas source
psordis(0:nx+1,0:ny+1)        _real  [part/s]  # diss. source of hydrogen
psorbgg(0:nx+1,0:ny+1,1:ngsp) _real  [part/s]  # diag artific neut backg source
psorbgz(0:nx+1,0:ny+1)        _real  [part/s]  # diag artific impur backg source
erliz(0:nx+1,0:ny+1)          _real  [J/s]     # H rad'n loss for ioniz'n
erlrc(0:nx+1,0:ny+1)          _real  [J/s]     # H rad'n loss for recom'n
vsoreec(0:nx+1,0:ny+1)	      _real  [J/s]     # cell ctr tot elec vol eng source
vsoree(0:nx+1,0:ny+1)	      _real  [J/s]     # cell ave tot elec vol eng source
pwrebkg(0:nx+1,0:ny+1)	      _real  [W/m**3] 
                               # elec energy backgrd source; limits te~tebg
pwribkg(0:nx+1,0:ny+1)	      _real  [W/m**3] 
                               # ion energy backgrd source; limits ti~tibg
wjdote(0:nx+1,0:ny+1)         _real  [J/s]     # Joule heating rate
wvh(0:nx+1,0:ny+1,1:nusp)     _real  [kg/m-s**3]  #ion viscous heating
smoc(0:nx+1,0:ny+1,1:nusp)    _real
smov(0:nx+1,0:ny+1,1:nusp)    _real
msor(0:nx+1,0:ny+1,1:nisp)    _real [kg-m/s**2]# ioniz. mom. source for ions
msorxr(0:nx+1,0:ny+1,1:nisp)  _real [kg-m/s**2]# cx&recomb. mom. sink for ions
seec(0:nx+1,0:ny+1)           _real
seev(0:nx+1,0:ny+1)           _real
seic(0:nx+1,0:ny+1)           _real
seiv(0:nx+1,0:ny+1)           _real
segc(0:nx+1,0:ny+1,1:ngsp)    _real [J/(sm**3)]#v_grad_P for neutral eng. eqn
resco(0:nx+1,0:ny+1,1:nisp)   _real
resng(0:nx+1,0:ny+1,1:ngsp)   _real
reseg(0:nx+1,0:ny+1,1:ngsp)   _real
resmo(0:nx+1,0:ny+1,1:nusp)   _real
resee(0:nx+1,0:ny+1)          _real
resei(0:nx+1,0:ny+1)          _real
resphi(0:nx+1,0:ny+1)         _real

***** MCN_dim:
# array bounds used in connection with Monte Carlo Neutrals
nstra		integer		/2/
# number of 'strata' or 'source groups' in Monte-Carlo-Neutrals model;
# i.e., a surface or volume element where neutrals originate;
# for multi-species neutrals, each is a separate source group.
nfl		integer		/1/
# number of plasma fluids recognized by Monte-Carlo-Neutrals model
natmi	integer		/1/
# number of atomic neutral species in EIRENE code
nmoli	integer		/1/
# number of molecular neutral species in EIRENE code
nioni	integer		/1/
# number of molecular ion species in EIRENE code
nxf	integer
# ix dimension from EIRENE file fort.44 or DEGAS2 file testdata.out
nyf	integer
# iy dimension from EIRENE file fort.44 or DEGAS2 file testdata.out
nmcsp	integer		/1/
# number of Monte Carlo species

***** MCN_sources:
# plasma source terms from Monte-Carlo-Neutrals model
ismcnon		integer		/0/
# flag for turning on plasma source terms from Monte-Carlo-Neutrals
# ismcnon=0  --> MCN plasma source terms are OFF (default)
# ismcnon=1  --> MCN-only is used for both Jac'n and RHS in pandf
# ismcnon=2  --> MCN-only is used for RHS, fluid-only is used for Jac'n
ismcnvar	integer		/0/
# flag for scaling plasma source terms from Monte-Carlo-Neutrals
# ismcnvar=0  --> MCN plasma source terms are constant (default)
# ismcnvar=1  --> MCN plasma source terms scale with plate currents

#Special case of neutral atoms emitted with finite energy from walls
eedisspl 	real /0/ [eV] +mcinput #energy loss for prompt dissociation at left plate
eedisspr 	real /0/ [eV] +mcinput #energy loss for prompt dissociation at right plate
eidisspl 	real /0/ [eV] +mcinput #energy gain for prompt dissociation at left plate
eidisspr 	real /0/ [eV] +mcinput #energy gain for prompt dissociation at right plate
cmntgpl 	real /0/ +mcinput #coeff. for neutral energy at left plate:  cmntipl*ti
cmntgpr 	real /0/ +mcinput #coeff. for neutral energy at right plate: cmntipr*ti
edisswo		real /0/ [eV] +mcinput #energy for prompt dissociation loss at outer wall
edisswi		real /0/ [eV] +mcinput #energy for prompt dissociation loss at private flux wall
cmntgwo 	real /0/ +mcinput #coeff. for neutral energy at outer plate:  cmntiwo*ti
cmntgwi 	real /0/ +mcinput #coeff. for neutral energy at private flux plate: cmntiwi*ti

#cfneut    		real /1./   +mcinput  #Coef to turn on all fluid neutrals contrib's to resid's
cfneutsor_ni	real /1/	+mcinput #coeff. for fluid neutral particle source in resco
cfneutsor_mi	real /1/	+mcinput #coeff. for fluid neutral momentum source in resmo
cfneutsor_ei	real /1/	+mcinput #coeff. for fluid neutral energy source in resei
cfneutsor_ee	real /1/	+mcinput #coeff. for fluid neutral energy source in resee

#cmneut    		real /0./    +mcinput #Coef to turn on all Monte Carlo neutral sources
cmneutsor_ni	real /1/	+mcinput #coeff. for MC neutral particle source in resco
cmneutsor_mi	real /1/	+mcinput #coeff. for MC neutral momentum source in resmo
cmneutsor_ei	real /1/	+mcinput #coeff. for MC neutral energy source in resei
cmneutsor_ee	real /1/	+mcinput #coeff. for MC neutral energy source in resee

cfneutdiv		real /1/	+mcinput #coeff. to turn on divergence of all fluid neutral fluxes
cfneutdiv_fng	real /1/	+mcinput #coeff. for div. fluid neutral particle flux in resng
cfneutdiv_fmg	real /1/	+mcinput #coeff. for div. fluid neutral momentum flux in resmo
cfneutdiv_feg	real /1/	+mcinput #coeff. for div. fluid neutral energy flux in resei

cmneutdiv		real /0/	+mcinput #coeff. to turn on divergence of all MC neutral fluxes
cmneutdiv_fng	real /1/	+mcinput #coeff. for div. fluid neutral particle flux in resng
cmneutdiv_fmg	real /1/	+mcinput #coeff. for div. fluid neutral momentum flux in resmo
cmneutdiv_feg	real /1/	+mcinput #coeff. for div. fluid neutral energy flux in resei


mcalpha_ng  real /2/  +mcinput #coeff. for blending kinetic and fluid ng
mcalpha_pg  real /2/  +mcinput #coeff. for blending kinetic and fluid pg
mcalpha_fng real /2/  +mcinput #coeff. for blending kinetic and fluid fng
mcalpha_fmg real /2/  +mcinput #coeff. for blending kinetic and fluid fmg
mcalpha_feg real /2/  +mcinput #coeff. for blending kinetic and fluid feg

### Scalars ###
ng_mc(0:nx+1,0:ny+1,nfl)		_real	[part/m**3]
# neutral gas density from Monte-Carlo-Neutrals model
ng_mc_rsd(0:nx+1,0:ny+1,nfl)	_real	[#]
# neutral gas density rsd from Monte-Carlo-Neutrals model
pg_mc(0:nx+1,0:ny+1,nfl)		_real	[Pa]
# neutral gas pressure from Monte-Carlo-Neutrals model
pg_mc_rsd(0:nx+1,0:ny+1,nfl)	_real	[#]
# neutral gas pressure rsd from Monte-Carlo-Neutrals model


ng_ue(0:nx+1,0:ny+1,nfl)		_real	[part/m**3]
# neutral gas density from Monte-Carlo-Neutrals model, blended with fluid result
ng_ue_rsd(0:nx+1,0:ny+1,nfl)	_real	[#]
# neutral gas density rsd from Monte-Carlo-Neutrals model, blended with fluid result

upg_ue(0:nx+1,0:ny+1,nfl)		_real	[m/s]
# neutral gas parallel velocity from Monte-Carlo-Neutrals model, blended with fluid result
upg_ue_rsd(0:nx+1,0:ny+1,nfl)	_real	[#]
# neutral gas parallel velocity rsd from Monte-Carlo-Neutrals model, blended with fluid result

pg_ue(0:nx+1,0:ny+1,nfl)		_real	[Pa]
# neutral gas pressure from Monte-Carlo-Neutrals model, blended with fluid result
pg_ue_rsd(0:nx+1,0:ny+1,nfl)	_real	[#]
# neutral gas pressure rsd from Monte-Carlo-Neutrals model, blended with fluid result

tg_ue(0:nx+1,0:ny+1,nfl)		_real	[J]
# neutral gas temperature from Monte-Carlo-Neutrals model
tg_ue_rsd(0:nx+1,0:ny+1,nfl)	_real	[#]
# neutral gas temperature rsd from Monte-Carlo-Neutrals model

sng_ue(0:nx+1,0:ny+1,1:nfl)		_real	[part/m**3-s]	#neutral particle source density (convective only)
smg_ue(0:nx+1,0:ny+1,1:nfl)		_real	[N/m**3]		#neutral parallel momentum source density 
seg_ue(0:nx+1,0:ny+1,1:nfl)		_real	[W/m**3]		#neutral energy source density (convective only)


### Vectors ###

jng_mc(0:nx+1,0:ny+1,nfl,3)			_real	[part/s*m**2]
# neutral gas particle flux density from Monte-Carlo-Neutrals model
jng_mc_rsd(0:nx+1,0:ny+1,nfl,3)		_real	[part/s*m**2]
# relative standard deviation
jng_ue(0:nx+1,0:ny+1,nfl,3)			_real	[part/s*m**2]
# neutral gas particle flux density interpolated to UEDGE grid
jng_ue_rsd(0:nx+1,0:ny+1,nfl,3)		_real	[#]
# relative standard deviation

vg_mc(0:nx+1,0:ny+1,nfl,3)			_real	[m/s]
# neutral gas velocity  from Monte-Carlo-Neutrals model
vg_mc_rsd(0:nx+1,0:ny+1,nfl,3)		_real	[#]
# relative standard deviation
vg_ue(0:nx+1,0:ny+1,nfl,3)			_real	[m/s]
# neutral gas velocity interpolated to UEDGE grid
vg_ue_rsd(0:nx+1,0:ny+1,nfl,3)		_real	[#]
# relative standard deviation

fngx_mc(0:nx+1,0:ny+1,nfl)	 	 	_real	[part/s]
# blended poloidal neutral gas particle flux from Monte-Carlo-Neutrals model
fngx_mc_rsd(0:nx+1,0:ny+1,nfl)		_real	[#]
# relative standard deviation
fngy_mc(0:nx+1,0:ny+1,nfl)	 	 	_real	[part/s]
# blended poloidal neutral gas particle flux from Monte-Carlo-Neutrals model
fngy_mc_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation

fngx_ue(0:nx+1,0:ny+1,nfl)	 	 	_real	[part/s]
# blended poloidal neutral gas particle flux on UEDGE grid
fngx_ue_rsd(0:nx+1,0:ny+1,nfl)		_real	[#]
# relative standard deviation
fngy_ue(0:nx+1,0:ny+1,nfl)	 	 	_real	[part/s]
# blended poloidal neutral gas particle flux on UEDGE grid
fngy_ue_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation

fmgx_mc(0:nx+1,0:ny+1,nfl)	 		_real	[N]
# blended poloidal neutral gas momentum flux from Monte-Carlo-Neutrals model
fmgx_mc_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation
fmgy_mc(0:nx+1,0:ny+1,nfl)	 		_real	[N]
# blended radial neutral gas momentum flux from Monte-Carlo-Neutrals model
fmgy_mc_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation

fmgx_ue(0:nx+1,0:ny+1,nfl)	 		_real	[N]
# blended poloidal neutral gas momentum flux on UEDGE grid
fmgx_ue_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation
fmgy_ue(0:nx+1,0:ny+1,nfl)	 		_real	[N]
# blended radial neutral gas momentum flux on UEDGE grid
fmgy_ue_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation
fmgxy_ue(0:nx+1,0:ny+1,nfl)	 		_real	[N]
# blended poloidal neutral gas momentum flux on nonorthog. UEDGE grid
fmgxy_ue_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation

jeg_mc(0:nx+1,0:ny+1,nfl,3)			_real	[W/m**2]
# neutral gas energy flux density vector from Monte-Carlo-Neutrals model
jeg_mc_rsd(0:nx+1,0:ny+1,nfl,3)		_real	[#]
# relative standard deviation
jeg_ue(0:nx+1,0:ny+1,nfl,3)			_real	[W/m**2]
# neutral gas energy flux density vector interpolated to UEDGE grid
jeg_ue_rsd(0:nx+1,0:ny+1,nfl,3)		_real	[#]
# relative standard deviation

fegx_mc(0:nx+1,0:ny+1,nfl)	 		_real	[W]
# blended poloidal neutral gas heat flux from Monte-Carlo-Neutrals model
fegx_mc_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation
fegy_mc(0:nx+1,0:ny+1,nfl)	 		_real	[W]
# blended radial neutral gas heat flux from Monte-Carlo-Neutrals model
fegy_mc_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation

fegx_ue(0:nx+1,0:ny+1,nfl)	 		_real	[W]
# blended poloidal neutral gas heat flux on UEDGE grid
fegx_ue_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation
fegy_ue(0:nx+1,0:ny+1,nfl)	 		_real	[W]
# blended radial neutral gas heat flux on UEDGE grid
fegy_ue_rsd(0:nx+1,0:ny+1,nfl)	 	_real	[#]
# relative standard deviation



 

### Tensors ###

stressg_mc(0:nx+1,0:ny+1,nfl,1:3,1:3)		_real	[Pa]
# neutral gas stress tensor from Monte-Carlo-Neutrals model
stressg_mc_rsd(0:nx+1,0:ny+1,nfl,1:3,1:3)	_real	[Pa]
# neutral gas stress tensor from Monte-Carlo-Neutrals model

stressg_ue(0:nx+1,0:ny+1,nfl,1:3,1:3)		_real	[Pa]
# neutral gas stress tensor interpolated to UEDGE grid
stressg_ue_rsd(0:nx+1,0:ny+1,nfl,1:3,1:3)	_real	[Pa]
# neutral gas stress tensor interpolated to UEDGE grid

pxz_mc(0:nx+1,0:ny+1,nfl)					_real	[Pa]
# neutral gas pressure from Monte-Carlo-Neutrals model
pxz_mc_rsd(0:nx+1,0:ny+1,nfl)				_real	[Pa]
# neutral gas pressure rsd from Monte-Carlo-Neutrals model

mcnsor_ni(0:nx+1,0:ny+1,1:nisp,1:nstra)	_real	[part/s]
# ion particle source from Monte-Carlo-Neutrals model
mcnsor_up(0:nx+1,0:ny+1,1:nisp,1:nstra)	_real	[kg-m/s**2]
# ion parallel momentum source from Monte-Carlo-Neutrals model
mcnsor_te(0:nx+1,0:ny+1,1:nstra)	_real	[J/s]
# electron thermal energy source from Monte-Carlo-Neutrals model
mcnsor_ti(0:nx+1,0:ny+1,1:nstra)	_real	[J/s]
# ion thermal energy source from Monte-Carlo-Neutrals model
mcncurr(1:nstra)	_real	[part/s]
# neutral source current from each strata in Monte-Carlo-Neutrals model
uesor_ni(0:nx+1,0:ny+1,1:nisp)	_real	[part/s]
# scaled ion particle source from Monte-Carlo-Neutrals model
uesor_up(0:nx+1,0:ny+1,1:nisp)	_real	[kg-m/s**2]
# scaled ion parallel momentum source from Monte-Carlo-Neutrals model
uesor_te(0:nx+1,0:ny+1)		_real	[J/s]
# scaled electron thermal energy source from Monte-Carlo-Neutrals model
uesor_ti(0:nx+1,0:ny+1)		_real	[J/s]
# scaled ion thermal energy source from Monte-Carlo-Neutrals model
uecurr(1:nstra)		_real	[part/s]
# neutral source current from each strata according to UEDGE plasma model
olduecurr(1:nstra)		_real	[part/s]
# neutral source current from each strata according to UEDGE plasma model
strascal(1:nstra)	_real
# scaling factor for plasma source terms due to each strata
wsor(1:nstra)	_real	[part/s]
# normalization constant for plasma source terms from EIRENE file 'fort.32'
esor(1:nstra)	_real
# unused constant from EIRENE file 'fort.32'
sni(0:nx+1,0:ny+1,1:nfl,1:nstra)	_real	[part/s]
# normalized ion particle sources from EIRENE file 'fort.32'
# or absolute ion particle source from DEGAS2
smo(0:nx+1,0:ny+1,1:nfl,1:nstra)	_real	[kg-m/s**2]
# normalized ion parallel momentum sources from EIRENE file 'fort.32',
# or absolute ion parallel momentum source from DEGAS2
smor(0:nx+1,0:ny+1,1:nfl,1:nstra)	_real	[kg-m/s**2]
# "radial" component of ion momentum source from DEGAS2
smophi(0:nx+1,0:ny+1,1:nfl,1:nstra)	_real	[kg-m/s**2]
# "toroidal" component of ion momentum source from DEGAS2
smoz(0:nx+1,0:ny+1,1:nfl,1:nstra)	_real	[kg-m/s**2]
# "vertical" component of ion momentum source from DEGAS2
see(0:nx+1,0:ny+1,1:nstra)	_real		[J/s]
# normalized electron energy source from EIRENE file 'fort.32'
# or absolute electron energy source from DEGAS2
sei(0:nx+1,0:ny+1,1:nstra)	_real		[J/s]
# normalized ion energy source from EIRENE file 'fort.32'
# or absolute ion energy source from DEGAS2
labela(1:nmcmx)		character*8
# data from Monte Carlo neutrals code:
$ C90 compiler does not allow dynamic character arrays, e.g., labela(1:natmi)
labelm(1:nmcmx)		character*8
# data from Monte Carlo neutrals code:
$ C90 compiler does not allow dynamic character arrays, e.g., labelm(1:nmoli)
labeli(1:nmcmx)		character*8
# data from Monte Carlo neutrals code:
$ C90 compiler does not allow dynamic character arrays, e.g., labeli(1:nioni)
naf(1:nxf,1:nyf,1:natmi)		_real
# data from Monte Carlo neutrals code:
# atomic neutral density
taf(1:nxf,1:nyf,1:natmi)		_real
# data from Monte Carlo neutrals code:
# atomic neutral temperature
nmf(1:nxf,1:nyf,1:nmoli)		_real
# data from Monte Carlo neutrals code:
# molecular neutral density
tmf(1:nxf,1:nyf,1:nmoli)		_real
# data from Monte Carlo neutrals code:
# molecular neutral temperature
ntf(1:nxf,1:nyf,1:nioni)		_real
# data from Monte Carlo neutrals code:
# molecular ion density
ttf(1:nxf,1:nyf,1:nioni)		_real
# data from Monte Carlo neutrals code:
# molecular ion temperature
fnax(1:nxf,1:nyf,1:natmi)	_real
# data from Monte Carlo neutrals code:
# x-particle flux of atomic neutrals
fnmx(1:nxf,1:nyf,1:nmoli)	_real
# data from Monte Carlo neutrals code:
# x-particle flux of molecular neutrals
fntx(1:nxf,1:nyf,1:nioni)	_real
# data from Monte Carlo neutrals code:
# x-particle flux of molecular ions
fnay(1:nxf,1:nyf,1:natmi)	_real
# data from Monte Carlo neutrals code:
# y-particle flux of atomic neutrals
fnmy(1:nxf,1:nyf,1:nmoli)	_real
# data from Monte Carlo neutrals code:
# y-particle flux of molecular neutrals
fnty(1:nxf,1:nyf,1:nioni)	_real
# data from Monte Carlo neutrals code:
# y-particle flux of molecular ions
fnaz(1:nxf,1:nyf,1:natmi)	_real
# data from Monte Carlo neutrals code:
# z-particle flux of atomic neutrals
fnmz(1:nxf,1:nyf,1:nmoli)	_real
# data from Monte Carlo neutrals code:
# z-particle flux of molecular neutrals
fntz(1:nxf,1:nyf,1:nioni)	_real
# data from Monte Carlo neutrals code:
# z-particle flux of molecular ions
feay(1:nxf,1:nyf,1:natmi)	_real
# data from EIRENE file fort.44:  y-energy flux of atomic neutrals
femy(1:nxf,1:nyf,1:nmoli)	_real
# data from EIRENE file fort.44:  y-energy flux of molecular neutrals
feax(1:nxf,1:nyf,1:natmi)	_real
# data from EIRENE file fort.44:  x-energy flux of atomic neutrals
femx(1:nxf,1:nyf,1:nmoli)	_real
# data from EIRENE file fort.44:  x-energy flux of molecular neutrals
hatm(1:nxf,1:nyf)		_real
# data from EIRENE file fort.44:  h-alpha radiation from atomic neutrals
hmol(1:nxf,1:nyf)		_real
# data from EIRENE file fort.44:  h-alpha radiation from molecular neutrals

***** MCN_test:
# densities, temperatures and fluxes from DEGAS2 code
labelmc(1:nmcmx)		character*8
# labels for Monte Carlo species
$ C90 compiler does not allow dynamic character arrays, e.g., labelmc(1:nmcsp)
nmc(1:nxf,1:nyf,1:nmcsp)		_real
# density from Monte Carlo neutrals code
tmc(1:nxf,1:nyf,1:nmcsp)		_real
# temperature from Monte Carlo neutrals code
fnmcx(1:nxf,1:nyf,1:nmcsp)	_real
# x-component of particle flux from Monte Carlo neutrals code
fnmcy(1:nxf,1:nyf,1:nmcsp)	_real
# y-component of particle flux from Monte Carlo neutrals code
fnmcz(1:nxf,1:nyf,1:nmcsp)	_real
# z-component of particle flux from Monte Carlo neutrals code
femcx(1:nxf,1:nyf,1:nmcsp)	_real
# x-component of energy flux from Monte Carlo neutrals code
femcy(1:nxf,1:nyf,1:nmcsp)	_real
# y-component of energy flux from Monte Carlo neutrals code
femcz(1:nxf,1:nyf,1:nmcsp)	_real
# z-component of energy flux from Monte Carlo neutrals code

***** MCN_bkgd:
# plasma background velocities for DEGAS2 Monte-Carlo-Neutrals model
v2c(1:nx,1:ny,1:nisp)	_real	[m/s]
# v2 velocity component at cell centers
vyc(1:nx,1:ny,1:nisp)	_real	[m/s]
# vy velocity component at cell centers
upc(1:nx,1:ny,1:nisp)	_real	[m/s]
# up velocity component at cell centers
uuc(1:nx,1:ny,1:nisp)	_real	[m/s]
# uu velocity component at cell centers
utc(1:nx,1:ny,1:nisp)	_real	[m/s]
# ut velocity component at cell centers
vr(1:nx,1:ny,1:nisp)	_real	[m/s]
# vr velocity component at cell centers
vphi(1:nx,1:ny,1:nisp)	_real	[m/s]
# vphi velocity component at cell centers
vz(1:nx,1:ny,1:nisp)	_real	[m/s]
# vz velocity component at cell centers
v2tg1(1:ny,1:nisp)	_real	[m/s]
# v2 velocity component at target plate number 1 (ix=0)
vytg1(1:ny,1:nisp)	_real	[m/s]
# vy velocity component at target plate number 1 (ix=0)
uptg1(1:ny,1:nisp)	_real	[m/s]
# up velocity component at target plate number 1 (ix=0)
uutg1(1:ny,1:nisp)	_real	[m/s]
# uu velocity component at target plate number 1 (ix=0)
uttg1(1:ny,1:nisp)	_real	[m/s]
# ut velocity component at target plate number 1 (ix=0)
vrtg1(1:ny,1:nisp)	_real	[m/s]
# vr velocity component at target plate number 1 (ix=0)
vphitg1(1:ny,1:nisp)	_real	[m/s]
# vphi velocity component at target plate number 1 (ix=0)
vztg1(1:ny,1:nisp)	_real	[m/s]
# vz velocity component at target plate number 1 (ix=0)
v2tg2(1:ny,1:nisp)	_real	[m/s]
# v2 velocity component at target plate number 2 (ix=nx)
vytg2(1:ny,1:nisp)	_real	[m/s]
# vy velocity component at target plate number 2 (ix=nx)
uptg2(1:ny,1:nisp)	_real	[m/s]
# up velocity component at target plate number 2 (ix=nx)
uutg2(1:ny,1:nisp)	_real	[m/s]
# uu velocity component at target plate number 2 (ix=nx)
uttg2(1:ny,1:nisp)	_real	[m/s]
# ut velocity component at target plate number 2 (ix=nx)
vrtg2(1:ny,1:nisp)	_real	[m/s]
# vr velocity component at target plate number 2 (ix=nx)
vphitg2(1:ny,1:nisp)	_real	[m/s]
# vphi velocity component at target plate number 2 (ix=nx)
vztg2(1:ny,1:nisp)	_real	[m/s]
# vz velocity component at target plate number 2 (ix=nx)

***** Ext_neutrals:
# External Neutrals API
extneutopt  integer   		/0/						#specifies which external neutral program to use
isextneuton integer   		/0/						#whether to use external neutrals implicitly within exmain
extneutmeth integer			/0/						#method for external neutrals: default=sources, 1=div. fluxes

uedgecmd  	character*16 	/"xuedge6dev"/			#uedge command
uedgescript character*16 	/"setup_neutrals.bas"/	#uedge script to run
uedgefile   character*16 	/"uedge_out.pdb"/		#uedge output file
uedgesave   character*16 	/"uedge_in.pdb"/		#uedge save file

eirenecmd  	character*16 	/"./eirene"/			#eirene command
eirenefile  character*16 	/"fort.40"/				#eirene output file

degas2cmd  	character*32 	/"./flighttest"/		#degas2 MC executable
degas2mpi  	character*32 	/"./flighttest_mpi"/	#degas2 MC executable for use with MPI
degas2file 	character*32    /"sources.out"/			#degas2 output file
gecmd		character*32    /"./readgeometry"/      #degas2 readgeometry executable
geufile		character*32 	/"readgeometry_input"/	#geometry input file for degas2
gedfile	 	character*32  	/"ge_uers.nc"/ 			#degas2 readgeometry output file
bkcmd  		character*32 	/"./readbackground"/ 	#degas2 readbackground executable
bkufile		character*32 	/"uedata.u"/		 	#uedge output file for degas2 readbackground
	#NOTE: the same file(grid) must also be used for readgeometry, as specified in the geufile
bkdfile 	character*32  	/"bk_uers.nc"/ 			#degas2 readbackground output file
degas2outcmd character*32   /"./outputbrowser"/     #degas2 outputbrowser executable 
degas2outscript character*32   /"output.input"/     #degas2 outputbrowser input file 
degas2outsh  character*32   /"seddata.sh *.dat"/          #sed script to clean up output files

mcnflights(1:nstramx) integer /nstramx*500/  		#number of mc pseudo-particle trajectories
ncsetcmd  	character*32 	/"ncset"/				#netcdf file editor command
ncsetvar 	character*32 	/"source_num_flights"/	#variable to edit in bkdfile netcdf file

ext_verbose logical			/.true./				# whether to print system call commands
istimecmdon logical			/.true./				# whether to time system call commands
ismpicmdon  logical 	 	/.false./				# whether to use MPI for external system call
mpicmd     	character*64 	/"mpirun"/ 				# MPI command
npopt       character*8     /"-np"/                 # option to specify # procs
npext      	integer 		/2/						# number of procs for external system call

runid_ext   character*80    /" "/					# description of run

get_neutral_sources		logical		/.true./			#whether to use neutral source data
get_neutral_moments		logical		/.false./			#whether to use neutral moment data

neut_output_dir	character*64	/"output"/					#output directory
neut_ng_file		character*64	/"neutral_density.dat"/		#neutral density file
neut_pg_file		character*64	/"neutral_pressure.dat"/	#neutral pressure file
neut_jng1_file	character*64	/"neutral_particle_flux_1.dat"/		#neutral particle flux: R
neut_jng2_file	character*64	/"neutral_particle_flux_2.dat"/		#neutral particle flux: T
neut_jng3_file	character*64	/"neutral_particle_flux_3.dat"/		#neutral particle flux: Z
neut_pg11_file	character*64	/"neutral_stress_11.dat"/	#neutral stress: RR
neut_pg22_file	character*64	/"neutral_stress_22.dat"/	#neutral stress: TT
neut_pg33_file	character*64	/"neutral_stress_33.dat"/	#neutral stress: ZZ
neut_pg23_file	character*64	/"neutral_stress_23.dat"/	#neutral stress: TZ
neut_pg31_file	character*64	/"neutral_stress_31.dat"/	#neutral stress: ZR
neut_pg12_file	character*64	/"neutral_stress_12.dat"/	#neutral stress: RT
neut_jeg1_file	character*64	/"neutral_heat_flux_1.dat"/		#neutral heat flux: R
neut_jeg2_file	character*64	/"neutral_heat_flux_2.dat"/		#neutral heat flux: T
neut_jeg3_file	character*64	/"neutral_heat_flux_3.dat"/		#neutral heat flux: Z

***** PNC_params:
#Plasma-Neutral Coupling Algorithm
pnc_verbose logical		/.true./					#print diagnostic info
pnc_opt		integer			/0/						# specifies choice of plasma-neutral coupling
pnc_step	integer			/0/						# step count for plasma-neutral coupling
pnc_maxstep	integer			/10/					# maximum number of coupled plasma+neutral steps
pnc_time	real			/0/						# time since beginning of coupled run
#pnc_ftol	real			/1e-4/					# ftol for PNC 
dtneut	    real	[s] 	/1.e20/	 				# time step for neutrals
dtplasma	real	[s]		/1e-6/					# time step for plasma-neutral coupling
dtold	    real	[s] 	/1.e20/	 				# old time step
relax_p     real            /1./					# relaxation parameter for plasma
relax_g     real            /1./					# relaxation parameter for neutral gas

pnc_ngs_mc   logical          /.false./             # replace fluid density with MC value
pnc_upgs_mc  logical          /.false./				# replace fluid parallel velocity with MC value
pnc_tgs_mc   logical          /.false./				# replace fluid temperature with MC value

pnc_histfile character*64 	/"pnc_hist.dat"/		# file name
pnc_fp		integer			/1001/					# file unit

pnc_print_norm	integer		/2/						# choice of normalization: 0=absolute, 1=relative to max

pnc_nsave    integer        /50/					# number of steps before saving data
pnc_savefile character*64 	/"pnc_step"/			# default pdb filename for saving pnc data
pnc_dobalance logical	    /.true./				# run dobalance function MUST BE READ FIRST!!!
pnc_balancefile character*64 /"pnc_balance.dat"/         # file to store diagnositic info for each step

## test preconditioner with alternate neutrals model
pnc_cfparvis(1:nispmx)  real  /nispmx*1./     		#factor for parallel visc. in preconditioner
pnc_cftravis(1:nispmx)  real  /nispmx*1./     		#factor for perp. visc. in preconditioner
pnc_cfni(1:nispmx)		real  /nispmx*1./			#factor for ni in preconditioner
pnc_cfup(1:nispmx)		real  /nispmx*1./			#factor for up in preconditioner



***** PNC_data:
#Plasma-Neutral Coupling Data Storage
ni_pnc(0:nx+1,0:ny+1,1:nisp)   _real  [1/m**3]     	#ion density in primary cell (ix,iy) at last pnc step
up_pnc(0:nx+1,0:ny+1,1:nisp)   _real  [m/s]     	#parallel velocity in primary cell (ix,iy) at last pnc step
ti_pnc(0:nx+1,0:ny+1)          _real  [J]     		#ion temperature in primary cell (ix,iy) at last pnc step
te_pnc(0:nx+1,0:ny+1)          _real  [J]     		#electron temperature in primary cell (ix,iy) at last pnc step
phi_pnc(0:nx+1,0:ny+1)         _real  [V]     		#potential in primary cell (ix,iy) at last pnc step

ng_pnc(0:nx+1,0:ny+1,1:nfl)   	_real  	[1/m**3]   	#neutral density in primary cell (ix,iy) at last pnc step
upg_pnc(0:nx+1,0:ny+1,1:nfl)  	_real  	[m/s]     	#parallel neutral velocity in primary cell (ix,iy) at last pnc step
tg_pnc(0:nx+1,0:ny+1,1:nfl)     _real  	[J]     	#neutral temperature in primary cell (ix,iy) at last pnc step
sng_pnc(0:nx+1,0:ny+1,1:nfl)	_real	[1/s]		#neutral particle source at last pnc step
smg_pnc(0:nx+1,0:ny+1,1:nfl)	_real	[N]			#neutral momentum source at last pnc step
seg_pnc(0:nx+1,0:ny+1,1:nfl)	_real	[W]			#neutral energy source at last pnc step

sni_pnc(0:nx+1,0:ny+1,1:nfl,1:nstra)     _real 	[1/s]	#density source in primary cell (ix,iy) at last pnc step
smor_pnc(0:nx+1,0:ny+1,1:nfl,1:nstra)    _real 	[N]  	#radial momentum source in primary cell (ix,iy) at last pnc step
smophi_pnc(0:nx+1,0:ny+1,1:nfl,1:nstra)  _real 	[N]  	#toroidal momentum source in primary cell (ix,iy) at last pnc step
smoz_pnc(0:nx+1,0:ny+1,1:nfl,1:nstra)    _real 	[N]  	#vertical momentum source in primary cell (ix,iy) at last pnc step
sei_pnc(0:nx+1,0:ny+1,1:nstra)           _real 	[W]    	#ion energy source in primary cell (ix,iy) at last pnc step
see_pnc(0:nx+1,0:ny+1,1:nstra)           _real 	[W]    	#electron energy source in primary cell (ix,iy) at last pnc step


#fngx_pnc(0:nx+1,0:ny+1,1:nfl)			_real	[1/s]		#neutral particle flux
#fngy_pnc(0:nx+1,0:ny+1,1:nfl)			_real	[1/s]		#neutral particle flux
#feg_pnc(0:nx+1,0:ny+1,1:3,1:nfl)		_real	[W]			#total neutral heat flux
#vg_pnc(0:nx+1,0:ny+1,1:3,1:nfl)		_real	[m/s]		#neutral velocity
#qg_pnc(0:nx+1,0:ny+1,1:3,1:nfl)		_real	[W/m**2]	#neutral heat flux


res_ni  	real  		#standard deviation of relative change in density
res_up  	real  		#standard deviation of relative change in parallel velocity
res_ti 	 	real  		#standard deviation of relative change in ion temperature
res_te  	real 		#standard deviation of relative change in electron temperature
res_phi 	real  		#standard deviation of relative change in electric potential

res_ng  	real  		#standard deviation of relative change in neutral density
res_upg  	real  		#standard deviation of relative change in neutral parallel velocity
res_tg  	real  		#standard deviation of relative change in neutral temperature
res_sng  	real  		#standard deviation of relative change in neutral density source
res_smg  	real  		#standard deviation of relative change in neutral parallel momentum source
res_seg  	real  		#standard deviation of relative change in neutral energy source

res_sni    	real		#standard deviation of relative change in ion particle source
res_smor   	real 		#standard deviation of relative change in radial ion momentum source
res_smophi 	real 		#standard deviation of relative change in toroidal ion momentum source
res_smoz   	real 		#standard deviation of relative change in vertical ion momentum source
res_sei    	real 		#standard deviation of relative change in ion energy source
res_see   	real	 	#standard deviation of relative change in electron energy source

del_ni 		real  		#maximum absolute change in density
del_up  	real  		#maximum absolute change in parallel velocity
del_ti  	real  		#maximum absolute change in ion temperature
del_te  	real  		#maximum absolute change in electron temperature
del_phi 	real  		#maximum absolute change in electric potential

del_ng  	real  		#maximum absolute change in neutral density
del_upg  	real  		#maximum absolute change in neutral parallel velocity
del_tg  	real  		#maximum absolute change in neutral temperature
del_sng  	real  		#maximum absolute change in neutral density source
del_smg  	real  		#maximum absolute change in neutral parallel momentum source
del_seg  	real  		#maximum absolute change in neutral energy source

del_sni    	real		#maximum absolute change in ion particle source
del_smor   	real		#maximum absolute change in radial ion momentum source
del_smophi 	real		#maximum absolute change in toroidal ion momentum source
del_smoz  	real		#maximum absolute change in vertical ion momentum source
del_sei  	real		#maximum absolute change in ion energy source 
del_see    	real		#maximum absolute change in electron energy source

***** Save_terms:
#Arrays to hold unperturbed values of particle-source terms
psorold(1:nisp)		_real	[part/s]  # unpert. ioniz. sources
psorxrold(1:nisp)	_real	[part/s]  # unpert. recom. & cx sources
msorold(1:nisp)		_real	[kg-m/s**2]  # unpert. ioniz. mom. sources
msorxrold(1:nisp)	_real	[kg-m/s**2]  # unpert. recom. & cx mom. sources

***** Time_dep_nwt:
#Old variables and time step for Newton iteration
ylodt(neqmx)  _real              #primary variables for ODE's at last output
nufak          real [1/s] /0./ 	 #pseudo freq. on precond.-Jac diag for nksol
nufak0         real [1/s]   	 #initial value of nufak0 saved (calc)
inufaknk       integer    /1/    #flag for using nufak in Krylov step of nksol
dtreal	       real [s] /1.e20/  #real timestep (both Jac and RHS) for nksol
				 #Do not use large nufak and small dtreal simult
dtdamp         real [s]   /0./   #mix old/new as frac=1/(1+(dtdamp/dtreal)**itdamp)
itdamp         real       /1./   #exponent for mix of old/new dt solutions
dtreal_old     real [s] /1.e20/  #previous value of dtreal
dtphi	       real [s] /1.e20/  #additional dt to relax phi equation
ydt_max        real              #maximum of yldot*sfscal
ydt_max0       real              #old value of ydt_max
alfnuf         real     /1./     #dtnewt->dtnewt*alfdtn*exp(ydt_max0/ydt_max)
                                 # **expdtn
expnuf         real     /0./     # see alfdtn
deldt          real     /0.3/    # frac. of var. change per cell for var. dt
dtoptx(0:nx+1,0:ny+1) _real      # spatial-depend. time step, min. in a cell
dtoptv(neqmx) _real 	 	 # variable-dependent time step; each var. diff
dtuse(neqmx)  _real 		 # time step used based on model_dt value
model_dt       integer  /0/      # determines dtuse for nksol time-step:
				 #=0, use dtreal
				 #=1, use dtreal*dtoptv/(dtreal+dtoptv)
				 #=2, use dtoptv
				 #=3, use sqrt(dtreal*dtoptv);
				 #=4, use dtreal*dtoptx/(dtreal+dtoptx)
				 #=5, use dtoptx
				 #=6, use sqrt(dtreal*dtoptx)

***** Decomp:
#Arrays required for the L-U decomposition.
ubw                  integer
lbw                  integer

***** Condition_number:
rconds(300,ngrid)	_real	# condition numbers of Jacobians

***** Jacobian:
#Jacobian matrix in compressed sparse row format
neqp1		integer		# Dimension (=neq+1) of jaci
nnzmx		integer		# Maximum no. of nonzeros in Jacobian matrix.
jac(nnzmx)	_real		# Nonzero entries of the Jacobian matrix.
				# This array, together with jacj and jaci,
				# contain the Jacobian in compressed sparse
				# row format.
jaci(neqp1)	_integer	# Nonzero structure of Jacobian matrix jac.
				# jaci(i+1) - jaci(i) = no. of nonzeros
				# in row i of jac.
jacj(nnzmx)	_integer	# Column indices of nonzero entries in jac.
isjacstnlon     integer    /0/  # Compute 9-pt stencil in ivl2gstnl - serial

***** Jacobian_csc:
#Jacobian matrix in compressed sparse column format
rcsc(nnzmx)	_real		# Nonzero entries of the Jacobian matrix.
				# This array, together with jcsc and icsc,
				# contain the Jacobian in compressed sparse
				# column format.
jcsc(neq+1)	_integer	# Nonzero structure of Jacobian matrix rcsc.
				# jcsc(j+1) - jcsc(j) = no. of nonzeros
				# in column j of rcsc.
icsc(nnzmx)	_integer	# Row indices of nonzero entries in rcsc.
yldot_pert(neqmx) _real         # Perturbed yldot within Jac_calc (diagnostic)
yldot_unpt(neqmx) _real		# Initial yldot with Jac_calc (diagnostic)


***** Jacobian_part:
#Part of Jacobian matrix arising from a particular equation
nnz1mx		integer		# Length of arrays in Jacobian_part
jac1(nnz1mx)	_real		# Nonzero elements of Jacobian
ia1(nnz1mx)	_integer	# Row indices of elements in jac1, or
				# nonzero structure of jac1 in csr format
ja1(nnz1mx)	_integer	# Column indices of elements in jac1

***** Jacreorder:
#Arrays used in performing a reordering of the rows and columns in the
#Jacobian matrix
perm(neq)	_integer	# Integer array containing the permutation
				# used in reordering the rows and columns of
				# the Jacobian matrix.
qperm(neq)	_integer	# Integer array holding the inverse of the
				# permutation in array perm.
levels(neq)	_integer	# Work array used by the bfs reordering
				# subroutine.   See subroutine bfs for
				# more details.
nlev		integer		# Number of levels in levels array.
				# See subroutine bfs for more details.
mask(neq)	_integer	# Work array used by the bfs reordering
				# subroutine.  See bfs subroutine.
maskval		integer		# Scalar used with mask.
ireorder	integer	    /1/ # Flag used to determine if a reordering
				# of the Jacobian matrix is desired.
				# = 1 means a reverse Cuthill-McKee
				#     reordering of the rows and columns
				#     of the Jacobian is done.
				# = 0 means no reordering.

***** Jacobian_full:
#Jacobian matrix in full storage format
jacfull(neq,neq) _real

***** Preconditioning:
#Parameters for type of preconditioning and sizes of matrices
premeth character*8 /"ilut"/  # type of preconditioning used in the
                              # linear iteration:
                              # ="banded" means use full banded jacobian as
                              #  preconditioner. Also used with mfnksol=4
                              # ="ilut" means use ilut preconditioning.
                              # ="inel" means use INEL ILU preconditioning
lenpfac        integer   /60/ # fudge factor to multiply neq by to get an
                              # estimate for the number of nonzeros in the
                              # preconditioner matrix.
lenplufac      integer  /100/ # fudge factor to multiply neq by to get an
                              # estimate for the number of nonzeros in the
                              # factored preconditioner matrix.
lenplumx       integer        # maximum number of nonzeros in the
                              # factored preconditioner matrix
                              # lenplumx = nnzmx + lenplufac*neq.
***** Ilutv:
#Control parameters for ILUT preconditioner
tolilut        real   /1.e-3/ # threshold tolerance for ILUT.
lfililut       integer   /50/ # fill-in parameter used in ILUT.  ILUT
                              # will allow up to lfililut additional nonzeros
                              # in each row of L and U.

***** Nonzero_diagonals:
ndiagmx        integer  /100/ # maximum number of nonzero diagonals in the
                              # Jacobian matrix
ndiag          integer        # actual number of nonzero diagonals in the
                              # Jacobian matrix
lfilinel       integer   /0/  # fill-in parameter used in INEL preconditioner
                              # lfilinel= number of additional diagonals
                              # used in the INEL ILU preconditioner
                              # lfilinel+ndiag .le. ndiagmx.
ndiagm         integer        # number of nonzero diagonals stored in the
                              # INEL ILU preconditioner
                              # = min(lfilinel+ndiag,ndiagmx)
adiag(neq,ndiagmx) _real      # diagonals of the Jacobian matrix
siginel(neq)   _real          # work array used by INEL precond5
fmuinel(neq)   _real          # work array used by INEL precond5
rwkd(ndiagmx)  _real          # work array used by cdiagsrt
iwkd1(2*neq-1) _integer       # number of nonzeros in each diagonal
iwkd2(ndiagmx) _integer       # work array used by cdiagsrt

***** UEint:
#Auxiliary variables for Ueinit.
GridFileName   character*200 /"gridue"/ +input
                              # name of Grid file to be read
newgeo         integer   /1/  +setup #flag to calculate new grid (1=yes)
mhdgeo         integer  /-1/  +input #flag for grid geometry
                              #mhdgeo =  2 ==> toroidal circular limiter
                              #mhdgeo =  1 ==> toroidal MHD equilibrium
                              #mhdgeo =  0 ==> cylindrical geometry
                              #mhdgeo = -1 ==> cartesian geometry
                              #mhdgeo = -2 ==> mag mirror (FRC-annulus)
gengrid        integer   /1/  +input #flag to generate grid, else read from file GridFileName
manualgrid     integer   /0/  +setup #flag whether to read grid values from gridue or memory
isgindx        integer   /1/  #=1 for interpolating grid based on indices
nfmax          integer   /10/
restart        integer   /0/  +input #flag for restart from previous case(yes=1)
initsol        integer   /0/  #flag to initially solve algebraic eqns for
                              #DASPK (yes=1)
ttbeg          real                  #initial Te in Joules = tinit/ev (calc)
tinit          real      /40./       #initial electron temperature Te in eV
tscal          real      /.5/        #ratio of initial Ti & Tg to Te
ngscal(ngspmx) real   /ngspmx*.1/    #ratio of initial gas density to ion dens
xgscal         real      /1./        #exponential scale of initial gas (m)
nibeg(1:nispmx) real  /nispmx*2.e19/ #initial ion density
minu(1:nispmx)  real  /nispmx*2./ +input
                                     #ion mass in units of proton mass (AMU)
ziin(1:nispmx)  real  /nispmx*1./ +input
                                     #ion charge read in, used to reset zi in
                                     #group Compla which gets erased on gallot
znuclin(1:nispmx) integer /nispmx*1./ +input #total nuclear charge of ion (i.d. isotope)
isallloc		integer   /0/        #=1 for local process. allocation with mpi
newaph			integer  /1/ +input #=1 calls aphread for hyd. atomic data;=0 not
newapi		integer /1/	     +input #=1, call readmc for new imp. data;=0, no						
pyrestart_file    character*80 /""/ #Python file that can also be used to restart
read_diffs		integer /0/	     +maybeinput #=0,a flag to signal whether to read diffusivities
dif_io		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write dif_use
tra_io		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write tra_use
dutm_io		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write dutm_use
kye_io		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kye_use
kyi_io		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
vy_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
vyup_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
vyte_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
vyti_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
fniyos_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
feeyosn_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
feiyosn_io 		integer /0/	     +maybeinput #=0,a flag to signal whether to read/write kyi_use
isvolsorext             integer /0/   +maybeinput #volsor sources if =0; or user sors if =1
isimpwallsor            integer /0/   +maybeinput #user imp ion wall fluxes; may have
					          # prob with ixendo for nyomitmx>0

***** Interp:
#Variables for the interpolation

uedge_savefile character*64 /"uedge_save.pdb"/ #default pdb filename for saving uedge data
isnintp                      integer /1/ +restart
                                       #switch to turn on new interpol. (=1)
                                       #also check isgindx switch in UEint
isimesh                      integer /1/ #flag for initial mesh => must copy
                                         #save variables and not interpolate
isumesh2		     integer /0/ #for parallel vers;=1, interp new mesh
nxold                        integer
nyold                        integer
nxoldg                       integer
nyoldg                       integer
ixlbo(1:nxpt)                _integer  #prev. grid value for ixlb
ixpt1o(1:nxpt)               _integer  #prev. grid value for ixpt1
ixpt2o(1:nxpt)               _integer  #prev. grid value for ixpt2
ixrbo(1:nxpt)                _integer  #prev. grid value for ixrb
iysptrxo                     integer   #prev. grid value for iysptrx
ixst(1:6)                    integer   #starting ix for 6 poloid region interp
ixsto(1:6)		     integer   #value of ixst on previous grid
ixend(1:6)                   integer   #end ix for 6 poloid region interp
ixendo(1:6)		     integer   #value of ixend on previous grid
xnrmo(0:nxold+1,0:nyold+1)   _real     #norm. x-grd; old x-grid, old y-grid
xvnrmo(0:nxold+1,0:nyold+1)  _real     #norm. xv-grd; old x-grid, old y-grid
xnrmox(0:nxold+1,0:ny+1)     _real     #norm. x-grd;nxold grd interp. to new ny
xvnrmox(0:nxold+1,0:ny+1)    _real     #norm. xv-grd;nxold grd interp.to new ny
xnrmnx(0:nx+1,0:ny+1)        _real     #norm. x-grd; second intermed. grid
xvnrmnx(0:nx+1,0:ny+1)       _real     #norm. xv-grd; second intermed. grid
ynrmo(0:nxold+1,0:nyold+1)   _real     #norm. y-grd; old x-grid, old y-grid
yvnrmo(0:nxold+1,0:nyold+1)  _real     #norm. yv-grd; old x-grid, old y-grid
ynrmox(0:nxold+1,0:ny+1)     _real     #norm. y-grd; old x-grid, new y-grid
yvnrmox(0:nxold+1,0:ny+1)    _real     #norm. yv-grd; old xv-grid, new y-grid
ynrmnx(0:nx+1,0:ny+1)        _real     #norm. y-grd; second intermed. grid
yvnrmnx(0:nx+1,0:ny+1)       _real     #norm. yv-grd; second intermed. grid
wrkint(0:nxold+1,0:ny+1)     _real     #wrk array; vars on old x, new y grid
wrkint2(0:nx+1,0:ny+1)       _real     #wrk array; vars on second interm. grid
ixmg(0:nxold+1,0:ny+1)       _integer  #ix index used for (ixo,iy) pt.
iyomg(0:nxold+1,0:ny+1)      _integer  #iyo index used for (ixo,iy) pt.
ixvmg(0:nxold+1,0:ny+1)      _integer  #ixv index used for (ixvo,iy) pt.
iyvomg(0:nxold+1,0:ny+1)     _integer  #iyvo index used for (ixvo,iy) pt.
ix2g(0:nx+1,0:ny+1)          _integer  #ix index for sec. interm. (ix,iy) pt.
iy2g(0:nx+1,0:ny+1)          _integer  #iy index for sec. interm. (ixo,iy) pt.
ixv2g(0:nx+1,0:ny+1)         _integer  #ixv index for sec. interm.(ixvo,iy) pt.
iyv2g(0:nx+1,0:ny+1)         _integer  #iyv index for sec.interm.(ixvo,iy) pt.
nis(0:nxold+1,0:nyold+1,1:nisp) _real [m^-3] +state
                                             #ion dens at last success. calc
tes(0:nxold+1,0:nyold+1)        _real [J]    #elec. temp at last success. calc
tis(0:nxold+1,0:nyold+1)        _real [J]    #ion temp at last success. calc
tgs(0:nxold+1,0:nyold+1,1:ngsp) _real [J]    #gas temp at last success. calc
phis(0:nxold+1,0:nyold+1)       _real [V]    #potential at last success. calc
ups(0:nxold+1,0:nyold+1,1:nisp) _real [m/s]  #parall. vel at last success. calc
ngs(0:nxold+1,0:nyold+1,1:ngsp) _real [m^-3] #gas dens at last success. calc.
afracs(0:nxold+1,0:nyold+1)     _real [ ]    +input #rel. imp. frac at last succ. calc

***** Global_vars:
# Arrays for primary variables over full mesh for domain decomposition
nisg(0:nxoldg+1,0:nyoldg+1,1:nisp) _real [m^-3] #global array for nis
tesg(0:nxoldg+1,0:nyoldg+1)        _real [J]    #global array for tes
tisg(0:nxoldg+1,0:nyoldg+1)        _real [J]    #global array for tis
tgsg(0:nxoldg+1,0:nyoldg+1)        _real [J]    #global array for tgs
phisg(0:nxoldg+1,0:nyoldg+1)       _real [V]    #global array for phis
upsg(0:nxoldg+1,0:nyoldg+1,1:nisp) _real [m/s]  #global array for ups
ngsg(0:nxoldg+1,0:nyoldg+1,1:ngsp) _real [m^-3] #global array for ngs
afracsg(0:nxoldg+1,0:nyoldg+1)     _real [ ]    #global array for afracs

***** Global_input:
# Arrays for real and integer input variables to be passed for domain decomp.
ipassin(1:100)		integer 	#integer input variables to be passed
rpassin(1:100)		real 		#real input variables to be passed
cpassin(1:30)		character*8 	#character input variables to be passed

***** Npes_mpi:
# Processor numbers for parallel version with mpi
npes		integer	/0/	#total number of processors
mype		integer	/-1/	#processor number of local processor (domain)
ismpion	        integer /0/     #flag to indicate using MPI (if=1)
hascomm         integer /0/     #flag indicates communicator has been set (if=1)
isparmultdt     integer /0/     #=1 for multistep parallel beyond 1st step

***** Indices_domain_dcg:
# Indices used for domain decomposition on the global mesh
isddcon		  integer   /0/	   #switch to turn on domain decomposition
ndleg(1:10,1:2)   integer   /20*1/ #number of x-domains in nxleg regions
ndxcore(1:10)     integer   /10*1/ #number of x-domains in nxcore(,1:2) regions
ndycore(1:10)     integer   /10*1/ #number of y-domains in core
ndysol(1:10)      integer   /10*1/ #number of y-domains in sol
idxpt(1:2)        integer   /2*0/  #PF/core domains with up touching X-point
ndxt              integer          #total number of x-domains
ndyt              integer          #total number of y-domains
ndomain           integer   /1/    #total number of domains
ndomain_orig      integer   /1/    #tot num orig domains before par_data gather
nvrsend           integer /10000/  #size of global real send/recv array for MPI
nvisend 	  integer /10000/  #size of global integer send/recv array for MPI
ixmin(ndomainmx)	integer    #min global ix for given domain
ixmax(ndomainmx)	integer    #max global ix for given domain
iymin(ndomainmx)	integer    #min global iy for given domain
iymax(ndomainmx)	integer    #max global iy for given domain
ixmnbcg(ndomainmx)	integer /ndomainmx*1/
                               	#B.C. type at ix=ixmin bdry;=0 inter.,=1 ex.
ixmxbcg(ndomainmx)	integer /ndomainmx*1/
			       	#B.C. type at ix=ixmax bdry;=0 inter.,=1 ex.
iymnbcg(ndomainmx)	integer /ndomainmx*1/
				#B.C. type at iy=iymin bdry;=0 inter.,=1 ex.
iymxbcg(ndomainmx)	integer /ndomainmx*1/
				#B.C. type at iy=iymax bdry;=0 inter.,=1 ex.
ncell(ndomainmx)	integer    #number of cells for given domain
idxp1g(ndomainmx)	integer    #domain to the right of given domain (ix+1)
idxm1g(ndomainmx)	integer    #domain to the left of given domain (ix-1)
idyp1g(ndomainmx)	integer    #domain above given domain (iy+1)
idym1g(ndomainmx)	integer    #domain below given domain (iy-1)
idcorng(ndomainmx,1:4)	integer    #domains touching corners; from lower left,
				   #numbering as in rm,zm: (l,r bot=1,2; top=3,4)
ixpt1g(ndomainmx)	integer    #ixpt1 for a given domain
ixpt2g(ndomainmx)	integer    #ixpt2 for a given domain
iysptrxg(ndomainmx)	integer    #iysptrx for a given domain
vrsend(nvrsend)        _real	   #real array used for passing global data via MPI
visend(nvisend)        _integer	   #int array used for passing global data via MPI
neq_locg(ndomainmx)    _integer    #number of vars per domain
neq_locgmx              integer    #maximum of neq_locg
ispwrbc(ndomainmx)      integer    #=1 for core pwr flux BC if corresp to ixpt2g

***** Indices_loc_glob_map:
# Indices that provide maps from loc-vars to glob-var and Jac entries
ivcum(ndomainmx)       _integer     #counter to build yl-local to yl-global map
ivloc2sdg(neqmx,ndomainmx) _integer #map loc-var to glob-var, single domain
ivloc2mdg(neqmx,ndomainmx) _integer #map loc-var to glob-var, mult domain
ivl2gstnl(neq_locgmx,9*numvar,ndomainmx) _integer /0/ # 1st arg loc-eqn number;

                            # 2nd arg poss Jac vars - global-mp; 3rd arg domain
iellast(neqmx,ndomainmx)   _integer #last meaningful entry into ivl2gstnl

***** Indices_domain_dcl:
# Indices used to connect domain with "neighbors"; known on local processor
nx_loc		integer		#number of ix cells for given processor
ny_loc		integer		#number of iy cells for given processor
nvrsendl	integer	/10000/	#size of local real send/recv array for MPI
nvisendl	integer	/10000/	#size of local integer send/recv array for MPI
ixmnbcl		integer   /1/   #B.C. type at ix=ixmin bdry;=0 intern,=1 extern
ixmxbcl		integer   /1/   #B.C. type at ix=ixmax bdry;=0 intern,=1 extern
iymnbcl		integer   /1/   #B.C. type at ix=iymin bdry;=0 intern,=1 extern
iymxbcl		integer   /1/   #B.C. type at iy=iymax bdry;=0 intern,=1 extern
idxp1		integer   	#domain to the right of given domain (ix+1)
idxm1		integer 	#domain to the left of given domain (ix-1)
idyp1		integer 	#domain to the above given domain (iy+1)
idym1		integer 	#domain to the below given domain (iy-1)
idcorn(1:4)	integer		#domains touching corners; from lower left,
				#numbering as in rm,zm: (l,r bot=1,2; top=3,4)
iv_totbdy(1:8)  integer  /8*0/  #number of elems. in bdry messages vrsendl
typebdyi(1:4)   integer /21,22,23,24/ #mpi tags for bdry iv_totbdy along edges
typecni(1:4)    integer /25,26,27,28/ #mpi tags for bdry iv_totbdy at corners
typebdy(1:4)    integer /11,12,13,14/ #mpi tags for bdry vrsendl along edges
typecn(1:4)     integer /15,16,17,18/ #mpi tags for bdry vrsendl at corners
vrsendl(nvrsendl) _real	        #real array used for passing local data via MPI
visendl(nvisendl) _integer      #int array used for passing local data via MPI
neq_locl        integer    /1/  #number of variables on local processor
numvarl         integer    /1/  #=numvar global via MPI_BCAST for parallel
ivloc2sdgl(nvisendl) _integer   #maps loc-var to glob-var, single domain
ivloc2mdgl(nvisendl) _integer   #maps loc-var to glob-var, mult domain
ivl2gstnll(neq_locl,9*numvarl) _integer /0/ # 1st arg loc-eqn number;
                                            # 2nd arg poss Jac vars-global-mp
ispwrbcl        integer   /1/   #=1 if domain has cell for core power BC
ixpt1l          integer   /0/   #local ixpt1 before par_data gather
ixpt2l          integer   /1/   #local ixpt2 before par_data gather
iysptrx1l       integer   /1/   #local iysptrx1 before par_data gather
ixlbl           integer   /0/   #local ixlb before par_data gather
ixrbl           integer   /1/   #local ixrb before par_data gather

***** Jacaux:
#Internal variables of jacnw
scrit             real            /1.e-4/
ylold(neqmx)     _real
yldot1(neqmx)    _real
yldot0(neqmx)    _real
normtype          integer /0/  #0,1,2 for max-norm, 1-norm, or 2-norm row scaling
fnormnw(neqmx)   _real
issfon            integer /1/  #=1 calc sfscal for row scaling (norml.) by nksol
isrnorm           integer /1/  #=1 causes row normaliza. of Jac. (see normtype)
jscalcol          integer /0/  #=1 causes column scaling for daspk
ExtendedJacPhi  integer /1/    #=1 extend bandwith for jacobian calculation when isphion=1

***** Variable_perturbation restart:
del             real	/1.e-8/		# fractional change for finite diffs
delpy           real	/-1.e-8/	# Forthon del; used to set del if > 0
dylconst        real    /1./            # factor in floor term in dyl
isjacreset      integer  /1/            # if=1, pandf1 reset for last variable

***** Jacobian_clipping restart:
jaccliplim	real	  /0./	# rel. value of elements to be retained
istopjac        integer   /0/   # flag to stop if non-zero elem at irstop,icstop
irstop          integer   /0/   # row (or eqn) index of non-zero stopping test
icstop          integer   /0/   # column (or var.) index of n-z stopping test

***** Newtaux:
icsum              integer
rwmin              real    /1.e-7/ #value of sumnew1 to stop Newton iter.
saux               real
saux1              real
sumnew             real
sumrdy             real
nmaxnewt           integer /15/    #max number of Newton iterations
ysave(2,neqmx)    _real            #last two yl's in Newton (1,) most recent
ycor(neqmx)       _real
saux2(nmaxnewt)    _real           #maximum update allowed in Newton
sumf(0:nmaxnewt)   _real           #ave value of right-hand-sides after Newton
irwd(nmaxnewt,2)   _integer
rwdmax(nmaxnewt,2) _real
rwdmin(nmaxnewt,2) _real
sumnew1(nmaxnewt)  _real           #average change in variables for Newton iter.
sumr1dy(nmaxnewt)  _real

***** Cdv:
#Commonly defined variables
ifexmain           integer /0/  #scalar to indicate if subroutine allocate
                                #is called by exmain.
                                #=1 means allocate is called by exmain,
                                #=0 means it is not.
exmain_aborted logical /.false./ # Set to .true. in Python version on control-C abort
iallcall	   integer /0/  #flag to signal first call to allocate

***** RZ_cell_info:
# RZ grid-cell center and face locations
rcn(0:nxm+1,0:nym+1)        _real [m]    # radial position of density cell
zcn(0:nxm+1,0:nym+1)        _real [m]    # vertical position of density cell
rfn(-1:nxm+1,-1:nym+1)      _real [m]    # radial position of density face
zfn(-1:nxm+1,-1:nym+1)      _real [m]    # vertical position of density face
rcv(0:nxm+1,0:nym+1)        _real [m]    # radial position of velocity cell
zcv(0:nxm+1,0:nym+1)        _real [m]    # vertical position of velocity cell
rfv(0:nxm+2,0:nym+2)        _real [m]    # radial position of velocity face
zfv(0:nxm+2,0:nym+2)        _real [m]    # vertical position of velocity face

***** Subs:
# Subroutines that can be called from the BASIS parser
#
###outputstats                                      subroutine
exmain                                           subroutine
exmain_prelims                                   subroutine
uedriv()                                         subroutine
convert()                                        subroutine
guardc()                                         subroutine
convsr_vo(i,j,yl:real)                           subroutine
	# in i
	# in j
	# in yl(*)  variables
convsr_aux(i,j)                          subroutine
	# in i
	# in j
pandf(i,j,neq:integer,t:real,yl:real,yldot:real) subroutine
	# in i   polodial index for perturbed variable
	# in j   radial index for perturbed vars for Jacobian calc (-1 for full RHS)
	# in  neq      total number of variables
	# in  t        present physical time
	# in  yl(*)    vector of unknowns
	# out yldot(*) RHS of ODE solver or RHS=0 for Newtown solver (NKSOL)
pandf1(i,j,ieq,neq:integer,t:real,yl:real,yldot:real) 	subroutine
	# in i   polodial index for perturbed variable
	# in j   radial index for perturbed vars for Jacobian calc (-1 for full RHS)
	# in  ieq        equation index for Jacobian calc
	# in  neq        total number of variables
	# in  t          present physical time
	# in  yl(*)      vector of unknowns
	# inout yldot(*) RHS of ODE solver or RHS=0 for Newtown solver (NKSOL)
bouncon(neq:integer,yl:real,yldot:real)          subroutine
	# Updates RHS (yldot) for boundary cells
	# in    neq        number of equations
	# in    yl(neq)    variables
	# inout yldot(neq) RHS values
poteneq(i,a:real,b:real)                         subroutine
	# Calculates RHS for potential equations
	# in    i          number of equations(neq)
	# in    a          yl variablevalues
	# inout b          yldot (RHS) values for potential eqn
potent_1dsol()					 subroutine
        # User diagnostic routine to calc potential from ex-field obtained
        # from parallel Ohms law when isphion+isphiofft=0; assume zero par current
ffun(neq:integer,t:real,yl:real,yldot:real)      subroutine
	# Used to calculate RHS for (old) direct Newton solve
	# in    neq        number of equations
	# in    t          physical time
	# in    yl(*)      most recent iterate of solution vector
	# in    yldot(neq) right hand sides
resid(t:real,y:real,yp:real,delta:real,ires,rp:real,ip) subroutine
	# Used to calculate RHS for DASPK solver
	# in    t          physical time
	# in    y(*)	   variables
	# out   yp(*)	   RHS
	# in	delta	   physical time
	# in 	ires	   error flag
	# in 	rp	   real parameters
	# in    ip	   integer parameters
ueinit()                                         subroutine
set_var_norm (job, neq, nvars, yl:real, norm_cons:real, \
              floor_cons:real, su:real)          subroutine
  #  Set su column scale factor
  #   in    job          methods 0: su=1, 1: global, 2 and 3 (default):local
  #   in    neq          total number of equations
  #   in    nvars        total number of variables
  #   in    yl(*)        most recent iterate of solution vector
  #   in    norm_cons(nvars)    global normalization constants
  #   in    floor_cons(nvars)   minimum normalization constants
  #   out   su(*)        column scale factors
gridseq()                                        subroutine
nphygeo()                                        subroutine
###init_par_meshg()                                 subroutine
  #  builds initial global mesh & Bcasts area_core
jacnw(neq,yl:real,f0:real,dt:real,wk:real,\
      wp:real,iwp)                   		 subroutine
  # calc LU of Jacobian at yl
  #   in    neq     total number of equations
  #   in    yl(*)   most recent iterate of solution vector
  #   in    f0(neq) function values f(yl)
  #   in    dt      false timestep to improve condition number
  #   inout wk(neq) work space
  #   inout wp(*)   matrix element of LU
  #   inout iwp(*)  array indices for elements of LU
psolnw(neq,yl:real,wk:real,wp:real,iwp,\
       bl:real,ierr)                		 subroutine
	# in    neq       total number of equations
	# in    yl(neq)   most recent iterate of solution vector
	# inout wk(neq)   work space
	# in    wp(*)     matrix elements of LU
	# in    iwp(*)    dimensions and array indices for elements of LU
	# inout bl(neq)   on input c of P*x=c, on output x
	# out   ierr      error flag
psolbody(neq,using:logical,suscal:real,wk:real,\
       wp:real,iwp,bl:real,ierr) 		 subroutine
	# in    neq          total number of equations
	# in    using        if su if used
	# in    suscal       scale factors for yl
	# inout wk(neq)      work array
	# in    wp(*)        matrix elements of LU
	# in    iwp(*)       dimensions and array indices for elements of LU
	# inout bl(neq)      on input, c of P*x=c; on output, x
	# out   ierr         error flag
csrcsc(neq,job,ipos,rcsc:real,icsc,jcsc,\
       jac:real,jacj,jaci)                       subroutine
  	# Note that csrcsc resides in uedge/svr/svrut4.f
  	# in    neq     dimension of matrix
  	# in    job     integer to indicate whether or not to fill
  	# in    ijob    starting position in ao, jao of transposed matrix
  	# in    rcsc(*) the matrix coefficients
  	# in    icsc(*) outgoing matrix column (nnz)
  	# in    jcsc(neq+1)
  	# inout jac(*)  outgoing matrix elements
  	# inout jacj(*) outgoing matrix column (nnz)
  	# inout jaci(neq+1)
allocate()                                	 subroutine
walsor()                                  	 subroutine
volsor()                                  	 subroutine
write_profs()                                  	 subroutine
read_profs()                                  	 subroutine
write_profs_boris(fname:string)		      	 subroutine
  	# in fname   the filename
read_profs_boris(fname:string,ierr)		       	 subroutine
  	# in fname   the filename
  	# in ierr    version flag
quadsvr(neq,a:real,b:real,c:real,d:real,yl:real,\
        yldot:real,ylprev:real,ylchng:real,\
        sfscal:real)         			 subroutine
  	# actual args are (neq,xs,xs1,xs2,xs3,ys1,ys2,ys3,yq1,yq2)
  	# in neq          number of equations
  	# in a            new xs
 	# in b            previous value of xs
  	# in c            previous value of xs
  	# in d            previous value of xs
  	# in yl(neq)      previous solution
  	# in yldot(neq)   previous solution
  	# in ylprev(neq)  previous solution
  	# inout y1chng(neq) one solution of quadratic (closest to solution ys1)
  	# inout sfscal(neq) other solution of quadratic
aplsb(nrow,ncol,a:real,ja,ia,s:real,b:real,\
      jb,ib,c:real,jc,ic,nzmax,iw,ierr)   	 subroutine
  	# Calculate C = A+s*B, where A, B, and C are matrices in
  	# compressed sparse row format (e.g., A is in a,ja,ia).
  	# Dimensions for output and work arrays:  c(nnzmx),jc(nnzmx),
  	# ic(nrow+1),iw(ncol).  ierr=0 if nnzmx was large enough.
  	#
  	# in nrow          row dimension of a
  	# in ncol          column dimension of a
  	# in a(*)          matrix a in compressed sparse row format
  	# in ja(*)         column number of each element of a (sparse row format)
  	# in ia(nrow+1)    ia(j) the starting index of elements of row j in a, ja
  	# in s             scalar applied to matrix b
  	# in b(*)          matrix b in compressed sparse row format
  	# in jb(*)         column number of each element of b (sparse row format)
  	# in ib(nrow+1)    ib(j) the starting index of elements of row j in b, jb
  	# inout c(nzmax)   resulting of a+s+b is stored here
  	# inout jc(nzmax)  column indices of c matrix (sparse row format)
  	# inout ic(nrow+1) starting index of each row of c
  	# inout iw(ncol)   integer workspace
  	# in nzmax         maximum number of nonzeros for c
  	# out ierr         error flag
jacmap()                                  	 subroutine
  	# output Jacobian map to file
map_var_jac1d()                                	 subroutine
  	# compute Jacobian stencil ivl2gstnl with 1 where elements
jacstnlout()                                  	 subroutine
  	# output Jacobian stencil with 4 where elements
jacout()                                  	 subroutine
      	# output Jacobian matrix in Boeing-Harwell format to a file
radintp(i,j,k,l,m,n,ii,jj,kk,ll,a:real,b:real,\
        c:real,d:real) 				 subroutine
      	# does interpolation in radial direction
      	# in i
      	# in j
      	# in k
      	# in l
      	# in m
      	# in n
      	# in ii
      	# in jj
     	# in kk
      	# in ll
      	# in a(kk+2,ll+2)
      	# in b(kk+2,ll+2)
      	# in c(kk+2,ll+2)
      	# inout d(kk+2,jj+2)
polintp(i,j,k,l,m,n,ii,jj,kk,ll,a:real,b:real,\
        c:real,d:real) 				 subroutine
      	# does interpolation in poloidal direction
      	# in i
      	# in j
      	# in k
      	# in l
      	# in m
      	# in n
      	# in ii
      	# in jj
      	# in kk
      	# in ll
      	# in a(ii+2,jj+2)
      	# in b(kk+2,jj+2)
      	# in c(kk+2,jj+2)
      	# inout d(ii+2,jj+2)
intpvar(a:real,b:real,i,j,k)	 		 subroutine
      	# does combined radial and poloidal interpolation
      	# in a(j+2, k+2)
      	# inout b(*,*)
      	# in i
      	# in j
      	# in k
engbal(a:real)                                 	 subroutine
      	# calc. arrays for postproc. energy balance
      	# in a  total input power for normalization
pradpltwl()                                      subroutine
      	# calc. radiation flux on plates from impurities and hydrogen
ebindz(za:integer, zn:integer)                   real function
      	# calculates the ionization energy for atomic charge state za-->za+1
      	# of impurity with total nuclear charge zn.
      	# in za   atomic charge
      	# in zn   nuclear charge
wtottim()					 subroutine
      	# writes out timing information
domain_dc()					 subroutine
      	# calculates indices of domains for domain decomposition
map_var_jac()					 subroutine
        # calculates indices of Jacobian; array ivl2gstnl is main output
bbb2wdf()                                        subroutine
      	# write file containing plasma information for DEGAS namelist
write30 (fname:string, runid:string)		subroutine
      	# write geometry data file 'fname' for EIRENE code
      	# in fname
      	# in runid
write31 (fname:string, runid:string)		subroutine
      	# write plasma data file 'fname' for EIRENE code
      	# in fname
      	# in runid
write_eirene					subroutine
      	# write geometry and plasma data files for EIRENE code
read32(fname:string)		subroutine
      	# read source term data file fname='fort.32' from EIRENE code
	# in fname
read44(fname:string)		subroutine
      	# read diagnostic data file fname='fort.44' from EIRENE code
	# in fname
writemcnfile(fname:string, runid:string)		subroutine
      	# write geometry and plasma background data for DEGAS2 code
      	# in fname			filename
      	# in runid			case id
readmcntest(fname:string)							subroutine
      	# read neutral density, temperature and flux data from DEGAS2 code
      	# in fname			filename
readmcnsor(fname:string)							subroutine
      	# read plasma source data from DEGAS2 code
      	# in fname			filename
readmcndens(fname:string)							subroutine
      	# read plasma density and pressure data from DEGAS2 code
      	# in fname			filename
readmcnoutput(fname:string,skip:integer,var:real,rsd:real)		subroutine
      	# read data from DEGAS2 outputbrowser data file
      	# in fname			filename
		# out var			data to be read
		# out rsd			rel. std. dev. of data
readmcnmoments(dname:string)									subroutine
      	# read plasma data from DEGAS2 outputbrowser files
      	# in fname			filename
lmode_roots(bcoef:complex, ccoef:complex, omega:complex)   subroutine
      	# in bcoef
      	# in ccoef
      	# inout omega(2)
lmode_chi_norm(kappabar:real, lte:real, rhos:real,
   cubrtnu:real, ti0:real, ted:real, zavg:real,
   lpi:real, lambdap:real,
   maxmag:real, nky:integer, kybeg:real, kyend:real, kya:real, kyb:real,
   tol:real, iprint:integer, islmodebeta:integer, kt:real,
   lmodechin:real, gammamax:real, kymax:real)   subroutine
	# in kappabar           field-line-averaged curvature [1/m]
	# in lte                L_Te = Ted / (dTed / dr0) [m]
	# in rhos               ion gyroradius at Ted [m]
	# in cubrtnu            cube root of collisionarility nu
	# in ti0                ti at "mid-plain" [eV]
	# in ted                Te at divertor plate [eV]
	# in zavg               average Z
	# in lpi                Pi / (dPi / dr) at "mid-plane" [m]
	# in lambdap            e (dPhi0 / dr0) / (dTed / dr0)
	# in maxmag             max magnitude of parab. step in bracketing ky
	# in nky                number of ky's used in maximizing growth rate
	# in kybeg              lower limit of acceptable kymax [none]
	# in kyend              upper limit of acceptable kymax [none]
	# in kya                one initial point in search for kymax [none]
	# in kyb                other initial point in search for kymax [none]
	# in tol                abs & rel tolerance in search for kymax
	# in iprint
	# in islmodebeta        =1 to turn on finite-beta correction
	# in kt                 parameter in finite-beta correction [none]
	# out lmodechin
	# out gammamax          maximum growth rate
	# out kymax             ky at maximum growth rate
hmode_chi_norm (gradvconst:real, cubrtnu:real, epsilon:real, lambdap:real,
   hmodechin:real)   subroutine
	# in gradvconst           factor involing rad. grad of v (parallel)
	# in cubrtnu              cube root of collisionality nu
	# in epsilon              rhos / L_Te
	# in lambdap              e (dPhi0 / dr0) / (dTed / dr0)
	# out hmodechin           normalized chi for H-mode turbulence
turb_chi (lmodechin:real, hmodechin:real, rhos:real, csed:real, lte:real,
   lambdap:real, cubrtnu:real, chi:real)   subroutine
	# in lmodechin             normalized chi for L-mode turbulence
	# in hmodechin             normalzied chi for H-mode turbulence
	# in rhos                  ion gyroradius at Ted [m]
	# in csed                  sound speed cs at divertor plate [m/s]
	# in lte                   Ted / (dTed / dr0) [m]
	# in lambdap               e (dPhi0 / dr0) / (dTed / dr0)
	# in cubrtnu               cube root of collisionarity nu
	# out chi                  turbulent transport coefficient chi (SI units)
read_zag()				   subroutine
     	# reads output for Zagorski's edge code
kappa (fqpsati:real, fqpsate:real, fqp:real)	real function
     # Calculates sheath drop (in units of Te) with modified form
     # that allows fqp > fqpsati and fqpsati < 0.
     # Adjustable form factors are kappamx and kappa0
mombal(ix:integer, ix1:integer, iy:integer)     subroutine
     # calculates impurity upi and frici from Hirschman's fmombal
     # ix1=ixp1(ix,iy)
     # in ix          index
     # in iy          index
     # out ix1      result
mombalni(ix:integer, ix1:integer, iy:integer)   subroutine
     # calculates impurity upi and frici from force balance
     # ix1=ixp1(ix,iy)
     # in ix          index
     # in iy          index
     # out ix1      result
###fluxsurfav1(inarray:real)    real function
     #  calcs flux surface average of 1D array inarray(ixpt1+1:ixpt2)
     #  Inarray dimensions nx*ny, distributed to processors if parallel
     # in inarray(*)       1-D flux array
###fluxsurfav2(inarray:real)   real function
     # calcs flux surface average of 2D array inarray.
     #  Inarray dimensions nx*ny, distributed to processors if parallel
     # in inarray(*,*)     2-D flux array
###interptrans()                                   subroutine
     # Deduces(interprets) radial transport coeff. from present profiles
     # Model is only applicable in core region; assumes plasma flux functions
fitdata2svar()                                  subroutine
     # evals expt profile fits; fills nis,tes,tis;
     # fitfrac1 cntrls mix of profiles at times 1 and 2
onedconteq()                                    subroutine
     # solves 1D convection/diffusion eqn; see Convdiffeqn group for vars

fit_neteti()                                    subroutine
        # interpolates expt profile data to UEDGE mesh for ne, Te, Ti
        # in fname   the filename
###build_global_soln()                             subroutine
     # evals expt profile fits; fills nis,tes,tis;
###reset_pe0_vars()                                subroutine
     # resets PE0's variables to local after call to build_global_soln
###getixiyloc(ixg:integer,iyg:integer,ixl:integer,iyl:integer,iownit:integer)   subroutine
     # Determines if the local processor owns the global index (ixg,iyg); if yes,
     # set iownit to one and return local index pair (ixl,iyl)
     # in ixg           global index
     # in iyg           global index
     # out iownit       1 if ixg,iyg are in the local processor's space, else 0
     # out ixl          local index
     # out iyl          local index
###getbdyindexlims()                               subroutine
     # calculates values of running index corresponding to start and end of various
     # portions of the edge boundary
###getixiybdy(lindex:integer,ix:integer,iy:integer,surfacename:character)   subroutine
     # returns poloidal and radial indices corresponding to single running
     # index that moves around the UEDGE boundary.  Also returns the name
     # of the bounding surface, for single null or  upper single null
     # in lindex      1-d index that runs around periphery of domain starting at junction
     #    of inner divertor and private flux wall, across inner divertor, around outer wall,
     #    across outer divertor, and around private flux wall to starting point
     # out ix               global poloidal index corresponding to 1D running lindex
     # out iy               global radial index corresponding to 1D running lindex
     # out surfacename      name of bounding surface corresponding to lindex
###set2dat2dpoint(darray:real,ix:integer,iy:integer,val:real)   subroutine 
     # Sets value of 2D array "darray" at global index point ix,iy to value "val".
     # Assumes that darray is dimensioned (0:nx+1,0:ny+1).
     # out darray(*,*)  global array being set
     # in ix       poloidal index at which darray is set
     # in iy       radial index at which darray is set
     # in val      value to which darray(ix,iy) is set
###set1dat1dpoint(darray:real,lindex:integer,val:real)   subroutine 
     # Sets value of 1D array "darray" at global index point ix,iy to value "val".
     # out darray(*)       global array being set
     # in index           poloidal index at which darray is set
     # in val             value to which darray(lindex) is set
###getat2dpoint(darray:real,ix:integer,iy:integer)   real function
     # Returns value of 2D array "darray" at global index point ix,iy
     # Assumes that darray is dimensioned (0:nx+1,0:ny+1).
     # in darray(*,*)      global array being queried
     # in ix              poloidal index at which darray is set
     # in iy              radial index at which darray is set
###getat1dpoint(darray:real,lindex:integer)   real function
     # Returns value of 2D array "darray" at index "lindex" that runs around periphery
     # in darray(*,*)     global darray being set
     # in lindex          poloidal index at which darray is set
ru_active(amumass:integer,znucleus:integer,charge:integer)  integer function
     # tests if given mass, charge, znucleus ion is active                  
     # in amumass is particle mass in AMU
     # in znucleus is the total charge of the atomic nucleus
     # in charge is particle charge in abs value of fundamental charge

upvisneo()                                      subroutine
     # computes ion neoclassical viscosity terms

jvisneo()                                       subroutine
     # computes neoclassical current terms

# Ext_neutrals ###
init_neutrals                					subroutine
     # call external neutral program (typically with system call)
     # variable extneutopt determines which neutral program option to use
init_degas2										subroutine
	 # initialize degas2 via system call: prepare options & setup files needed for readbackground + flighttest
init_eirene										subroutine
	 # initialize eirene via system call -- not yet working! -- use for testing purposes

run_neutrals                					subroutine
     # call external neutral program (typically with system call)
     # variable extneutopt determines which neutral program option to use
run_uedge										subroutine
	 # run a new uedge session via system call
run_degas2										subroutine
	 # run degas2 via system call
run_eirene										subroutine
	 # run eirene via system call -- not yet working! -- use for testing purposes

convertmcnsor									subroutine
	 # convert degas2 sources into uedge sources

# PNC ### Plasma-Neutral Coupling

uedge_plasma											subroutine
     # run uedge plasma fluid model without neutrals -- NOTE: uses timestep dtplasma
uedge_neutrals											subroutine
     # run uedge neutral fluid model without plasma  -- NOTE: uses timestep dtneut
uedge_uedge										subroutine
     # run coupled uedge plasma fluid + uedge neutral fluid model
uedge_degas2									subroutine
     # run coupled uedge plasma fluid + degas2 monte carlo neutral model
run_pnc											subroutine
     # run coupled uedge plasma fluid + external neutral physics model


uedge_save_pdb(filename:string)				integer function
	 # save standard uedge data: nis, ups, tis, tes, phis, ngs to pdb file
uedge_save 									integer function
	 # save standard uedge data: nis, ups, tis, tes, phis, ngs to pdb file
	 # filename set in uedge_savefile
uedge_read_pdb(filename:string)				integer function
	 # read standard uedge data: nis, ups, tis, tes, phis, ngs from pdb file
uedge_read 									integer function
	 # read standard uedge data: nis, ups, tis, tes, phis, ngs from pdb file
	 # filename set in uedge_savefile
mcnsor_save_pdb(filename:string)			integer function
	 # save neutral source data: sni, smor, smophi, smoz, sei, see
mcnsor_append_pdb(filename:string)			integer function
	 # save neutral source data: sni, smor, smophi, smoz, sei, see
pnc_save_pdb(filename:string)				integer function
	 # save uedge plasma and neutral source data

get_fnrm(dtreal_try:real)					real function
	# get fnrm given present value of yl and yldot for time step dtreal_try

test_opt(;optarg:string)					subroutine
	# test optional argument

test_parser(;optarg:string)					subroutine
	# test optional argument

interpmcnvec(mcvar:real,uevar:real,mcvar_rsd:real,uevar_rsd:real)	subroutine
	# interpolate vector components from MC grid to uedge grid

convertmcnvec(mcvar:real,uevar:real,mcvar_rsd:real,uevar_rsd:real,sgn:integer)	subroutine
	# convert vector components from RTZ to XYP

convertmcnvector(mcvar:real,uevar:real,mcvar_rsd:real,uevar_rsd:real)	subroutine
	# convert vector components from MC grid to uedge grid

convertmcntensor(mcvar:real,uevar:real,mcvar_rsd:real,uevar_rsd:real)	subroutine
	# convert vector components from MC grid to uedge grid

convertmcnmoments							subroutine
	# convert moments from MC grid to uedge grid

mcndivide(out:real,var:real,dens:real,out_rsd:real,var_rsd:real,dens_rsd:real)  subroutine
	# out = MC kinetic variable divided by density, also calculates rel. std. dev.

mcuedivide(out:real,var:real,dens:real,out_rsd:real,var_rsd:real,dens_rsd:real)	subroutine
	# out = MC kinetic variable divided by density, also calculates rel. std. dev.

mcnrsdfix(mcrsd:real)						subroutine
	# fix range of rel. std. dev. to to lie within (0,1]

mcnblend(out:real,uevar:real,mcvar:real,out_rsd:real,mcrsd:real,alpha:real)	subroutine
	# interpolation between fluid and MC kinetic results based on rel. std. dev.
	# out = mcvar*(1-mcrsd)**alpha + uevar*(1-(1-mcrsd)**alpha)

mult23(var2:real,var3:real,n3:integer)			function	
	# component-wise multiplication of 2d*3d variable along x and y directions

mult24(var2:real,var4:real,n3:integer,n4:integer)	function
	# component-wise multiplication of 2d*4d variable along x and y directions

mult34(var2:real,var4:real,n3:integer,n4:integer)	function
	# component-wise multiplication of 3d*4d variable along x, y, and third directions

wallflux					subroutine
	# diagnostic calc particle/heat flux to walls; alt. to balancee

plateflux					subroutine
	# diagnostic calc particle/heat flux to plates; alt. to balancee

***** Imprad:
# variables associated with impurity radiation
isimpon                 integer   /0/ +input
                               # switch for impurity model:
                               # 0 for no impurities
                               # 2 for fixed-fraction model
                               # 3 for average-impurity-ion model(disabled)
                               # 4 for INEL multi-charge-state model(disabled)
                               # 5 for Hirshman's reduced-ion model
                               # 6 for force-balance model or nusp_imp > 0;
                               #   see also isofric for full-Z drag term
                               # 7 for simultaneous fixed-fraction and
                               #       multi-charge-state (isimpon=6) models
nusp_imp        integer  /0/   +input # fixes nusp for total num. of par. mom. eqns.
isupimpap       integer  /1/   +input # =1 includes imp atm phys in up eqn; =0, omits
ismctab		integer		/1/	+input
# Determines which data is used for multi-charge-state rates.
#	=1  tables originally generated by R. Campbell for D. Knoll,
#	    data file name is specified by inelmc=....
#	    corresponding rate evaluation routines are imprates and radimpmc.
#	=2  tables generated by code from B. Braams,
#	    data file name is specified by mcfilename=...,
#	    corresponding rate evaluation routines are mcrates and radmc.
nzloc(0:nzspmx)		_real	[/m**3]
                               # imp. dens. for each Z at one grid cell
impradloc(0:nzspmx)	_real	[Watts/m**3]
                  # rad. power loss density for each Z at one grid cell
pwrzec(0:nx+1,0:ny+1)	_real	[Watts/m**3]
                               # elec energy loss via impurities at cell-cntr
pwrze(0:nx+1,0:ny+1)	_real	[Watts/m**3]
                               # elec energy loss via impurities; cell-ave
pradc(0:nx+1,0:ny+1)	_real	[Watts/m**3]
                               # cell ctr total impurity radiation
pradcff(0:nx+1,0:ny+1)	_real	[Watts/m**3]
                               # cell ctr impurity radiation (fixed-fraction)
prad(0:nx+1,0:ny+1)	_real	[Watts/m**3]
                               # cell ave total impurity radiation
pradzc(0:nx+1,0:ny+1,0:nzspmx,1:ngsp-1)	_real	[Watts/m**3]
                               # cell ctr imp rad due to each imp. ch. state
pradz(0:nx+1,0:ny+1,0:nzspmx,1:ngsp-1)	_real	[Watts/m**3]
                               # cell ave imp rad due to each imp. ch. state
na(0:nx+1,0:ny+1)	_real	[/m**3]
                               # atomic density of impurity (=afrac*ne)
ntau(0:nx+1,0:ny+1)	_real	[sec/m**3]
                               # confinement parameter for impurity (=atau*ne)
nratio(0:nx+1,0:ny+1)	_real
                               # ratio of neutrals to electrons
afrac(0:nx+1,0:ny+1)	_real	/.00/ +maybeinput
                               # atomic impur conc; set internally to afracs
atau(0:nx+1,0:ny+1)	_real	/1./	[sec] +input
                               # lifetime of impurity
tau1(0:nx+1,0:ny+1)	_real	# time to escape to inboard divertor plate
tau2(0:nx+1,0:ny+1)	_real	# time to escape to outboard divertor plate

***** Impurity_source_flux:
# Arrays for impurity-source fluxes calculated from api variables
fnzysi(0:nx+1,nzspt)	_real	# profiles along inner wall
fnzyso(0:nx+1,nzspt)	_real	# profiles along outer wall

***** Reduced_ion_interface:
# Arrays for info transfer between subroutines mombal and fmombal
misotope		integer  # number of isotopes (including electrons)
nchstate		integer  # maximum charge state among all isotopes
natomic(1:MXMISO)	integer  # maximum charge state of each isotope
amu(1:misotope)		_real	[none]     # atomic mass, relative to proton
tempa(1:misotope)	_real	[J]        # temperature
qneut(1:misotope)	_real	[J/m**2-s] # parallel heat flux of neutral
uneut(1:misotope)	_real	[m/s]      # parallel flow speed of neutral
den(1:misotope,0:nchstate)	_real	[1/m**3] # density
gradp(1:misotope,1:nchstate)	_real	[J/m**4] # parallel pressure grad
gradt(1:misotope,1:nchstate)	_real	[J/m**4] # parallel temp gradient
friction(1:misotope,1:nchstate)	_real	[J/m**4] # parallel friction force
friccomp(1:misotope,1:nchstate,1:5) _real [J/m**4] # par friction components
						 #friccomp(,,1)~ upi-upj
						 #friccomp(,,2)~ qcond
						 #friccomp(,,3)~ h;higher mom
						 #friccomp(,,4)~ caplam;elec?
				  		 #friccomp(,,5)~ ioniz/recomb
nuion(1:misotope,0:nchstate)	_real	[1/s]    # ionization rate
nurec(1:misotope,1:nchstate)	_real	[1/s]    # recombination rate
qcond(1:misotope,1:nchstate)	_real	[J/m**2-s] # parallel heat flux
ucond(1:misotope,1:nchstate)	_real	[m/s]      # parallel flow speed
dztot(1:misotope)               _real   [1/m**3] # total local isotope density

***** Bdy_indexlims:
# Limits of running index that goes around boundary
ib_idiv			integer    /0/         # begin inner divertor
ie_idiv			integer    /0/         # end inner divertor
ib_comwall		integer    /0/         # begin common flux wall
ie_comwall		integer    /0/         # end common flux wall
ib_odiv			integer    /0/         # begin outer divertor
ie_odiv			integer    /0/         # end outer divertor
ib_opfwall			integer    /0/         # begin outer part of pf wall
ie_opfwall			integer    /0/         # end outer part of pf wall
ib_ipfwall			integer    /0/         # begin inner part of pf wall
ie_ipfwall			integer    /0/         # end inner part of pf wall

***** Solver_work_arrays:
# Work arrays used to communicate with solver routines
liw			integer	  /1/	# length of iwork
lrw			integer	  /1/	# length of rwork
iwork(liw)		_integer	# integer work array
rwork(lrw)		_real		# real work array

***** Jac_work_arrays:
# work arrays needed by Jac routines when using cvode (FPRECO)
liwp			integer   /1/	# length of iwwp
lwp			integer	  /1/	# length of wwp
iwwp(liwp)	       _integer		# integer work array
wwp(lwp)	       _real		# real work array

***** Temporary_work_arrays:
rwk1(neq+1)		_real
rwk2(neq+1)		_real
iwk1(neq+1)		_integer
iwk2(neq+1)		_integer
iwk3(neq+1)		_integer

***** Zag_output:
# Arrays used to store output from Zagorski's FPIT edge code
# PLASMA common block
nezag(imx+1,imy+1)	_real	# Zagorski's electron density
nizag(imx+1,imy+1)	_real	# Zagorski's hydrogen ion density
tezag(imx+1,imy+1)	_real	# Zagorski's electron temperature
tizag(imx+1,imy+1)	_real	# Zagorski's electron temperature
vizag(imx+1,imy+1)	_real	# Zagorski's velocity 1
uizag(imx+1,imy+1)	_real	# Zagorski's velocity 2
v0zag(imx+1,imy+1,5) 	_real	# Zagorski's velocity 3
u0zag(imx+1,imy+1,5) 	_real	# Zagorski's velocity 4
zeffzag(imx+1,imy+1)	_real	# Zagorski's Zeff
elfzag(imx+1,imy+1)	_real	# Zagorski's elf
vezag(imx+1,imy+1)	_real	# Zagorski's ve
uezag(imx+1,imy+1)	_real	# Zagorski's ue
curxzag(imx+1,imy+1)	_real	# Zagorski's curx
curyzag(imx+1,imy+1)	_real	# Zagorski's cury
t0zag(imx+1,imy+1,5) 	_real	# Zagorski's t0
n0zag(imx+1,imy+1,5) 	_real	# Zagorski's n0
tz0zag(imx+1,imy+1) 	_real	# Zagorski's tz0
# IMPUR common block
nzzag(imx+1,imy+1,lnst+1) _real	# Zagorski's nz
vzzag(imx+1,imy+1,lnst+1) _real	# Zagorski's vz
uzzag(imx+1,imy+1,lnst+1) _real	# Zagorski's uz
# GEOMETRY common block
rxzag(imx+1,imy+1)	_real	# Zagorski's rx
ryzag(imx+1,imy+1)	_real	# Zagorski's ry
ggzag(imx+1,imy+1)	_real	# Zagorski's gg
bratiozag(imx+1,imy+1)	_real	# Zagorski's bratio
viparzag(imx+1,imy+1)	_real	# Zagorski's vipar
veparzag(imx+1,imy+1)	_real	# Zagorski's vepar
vzparzag(imx+1,imy+1,lnst+1) _real # Zagorski's vzpar
# NET common block
xpzag		 	 real   # Zagorski's xp
xszag			 real	# Zagorski's xs
xkzag			 real	# Zagorski's xk
ypzag			 real	# Zagorski's yp
ykzag			 real	# Zagorski's yk
xzag(imx+2)		_real	# Zagorski's x
yzag(imy+2)		_real	# Zagorski's y
# WARIANT common block
inh			integer
inz		/0/	integer
ihf			integer
istype			integer
ibound			integer
iboundz			integer
zs			real
spuff			real
# PARAM common block
lst			integer
wx			real
wy			real
ht			real
rmach			real
recyc			real
zrecyc			real
imap			integer
snz			real
sn			real
sqe			real
sqi			real
fe(imx+1)	       _real
fi(imx+1)	       _real
fe0			real
fi0			real
flime			real
flimi			real
wkr			real
# OUT common block
xsi			real
zxsi			real
sdod(imx+1,imy+1,lnst+1) _real
yielh(imx+1)	       _real
yielz(imx+1,lnst+1)    _real

***** Ident_vars:
uedge_ver  character*80 /'$Name: 8.0.6-beta.0$'/
uedge_date character*80 /'Version date in README_Uedge_vers in dir uedge'/
session_id  integer /0/ # Identifier for use with uetools
max_session_id  integer /0/ # Identifier for max allocated runs, use with uetools
exmain_evals  integer /0/ # Number of successfull exmain evaluations

***** Last_group_ex_sav_var:
# Last group in bbb where new variables from read save files get stored

***** MpiComponent:
# Link to setting and getting communicators
###set_uedgeComm(comm:integer) integer function   # in comm the new communicator

###uedge_mpiInit	subroutine
     # Do mpiInit if it is not a component
###uedge_petscInit	subroutine
     # Enable the initialization of petsc from
###uedge_mpiFinal subroutine
     # call MPI_Finalize() to exit MPI
###uedge_petscFinal subroutine
     # calls PetscFinalizeWrap() to exit petsc
###uedge_petscInsertOpts subroutine
     # calls PetscInsertOpts() to set petscoptions
###uedge_reset	subroutine
     # Does the reset

***** Logging:
# Variables/Methods required for logging output to file instead of console
logfname            character*64         /'uedgelog'/
     # name of the log file to which to write

***** Convdiffeqn:

# Variables used in the standalone sub onedconteq solving cont eqn
nxx              integer  /50/      # number of mesh points
alfz             real     /10./
vrfac            real     /1./      # scales convective velocity
sp               real     /20./
courant          real     /0.9/
tendoned         real     /0.1/     # final time
ndtmax           integer  /10000/   # max number of timesteps allowed
ntim             integer  /50/      # number of output times
ito              integer   /1/
xcz(1:nxx)       _real  
xfz(1:nxx)       _real
vrz(1:nxx)       _real
drz(1:nxx)       _real
dens(1:nxx)      _real
vrhs(1:nxx)      _real
drhs(1:nxx)      _real
gampz(1:nxx)     _real
gampzt(1:nxx,1:ntim) _real
nnt(1:nxx,1:ntim)    _real     /0./     # solution at output times
timo(1:ntim)         _real     /0./     # output times

setLogFile(filename:string)   subroutine
     # sets the filename base to which to direct log output
writeToLog(message:string)    subroutine
     # writes a message to the uedge log file

***** ParallelEval: # added by J.Guterl
ParallelJac     integer /0/        # [0]: serial jacobian calc [1] parallel jacobian calc
ParallelPandf1     integer /0/        # [0]: serial pandf1 calc [1] parallel pandf1 calc

**** PandfTiming: # added by J.Guterl
TimingPandfOn integer /1/
TimingPandf integer /1/
TimePandf real /0.0/
TotTimePandf real /0.0/
TimeConvert0 real /0.0/
TotTimeConvert0 real /0.0/
TimeConvert1 real /0.0/
TotTimeConvert1 real /0.0/
TimeNeudif real /0.0/
TotTimeNeudif real /0.0/
Timefd2tra real /0.0/
TotTimefd2tra real /0.0/
PrintTimingPandf() subroutine


**** Uetools:
# Variables for explicit use by uetools
dummy           real    /0./     # Dummy variable for testing etc.
gridmorph       real    /0./     # Grid morphing factor to be used with UETOOLS 

