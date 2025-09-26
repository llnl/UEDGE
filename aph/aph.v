aph   # Atomic physics (hydrogenic)

***** Physical_constants:
ev_aph	real	/1.6022e-19/   # 1 electron volt in Joules
m_prot  real    /1.67e-27/     # proton mass

***** Data_input:
# type of data and directory where it is located
issgvcxc integer       /0/    #=1 fixes sig-v_cx; =2 fixes sig=sfvcxc
sgvcxc   real [m^2/s] /2.e-14/# const. value of sigv_cx for issgvcxc=1; 
                              # const. sig for issgvcxc=2, sig-v=sig*sqrt(T/mi)
			      # with units now [m^2]
isaphdir integer       /1/    #=1 uses aphdir; =0 uses explicit rate file names
crmnfname character*120 /"crumpet_nrates.dat"/ # Name of CRM rate data file
crmefname character*120 /"crumpet_Erates.dat"/ # Name of CRM E-sink data file
aphdir	character*120 # name of directory containing data files
data_directory	character*120 # another dirname containing data files. This is to be be passed in

***** Ionization_energy:
erad	real [eV] /25./ +input # tot elec engy loss/ioniz (rad+binding) if istabon=0

***** Rtdata:
# hydrogenic rate table data from ADPAK via Braams' rate code
labelht	character*120	# header information from ADPAK data file
htns	integer		# number of atomic species
htnn	integer		# number of intervals in ADPAK density data
htnt	integer		# number of intervals in ADPAK temperature data
htza(0:htns-1)	_real	# atomic charge state of species
htzn(0:htns-1)	_real	# nuclear charge state of species
htn(0:htnn)	_real	# density data in ADPAK table
htt(0:htnt)	_real	# temperature data in ADPAK table
htln(0:htnn)	_real	# log(htn) where htn is ADPAK density data
htlt(0:htnt)	_real	# log(htt) where htt is ADPAK temperature data
htlsa(0:htnt,0:htnn,0:htns-1)	_real	
	# log(htsa) where htsa is ADPAK rate parameter data for ionization
htlra(0:htnt,0:htnn,0:htns-1)	_real	
	# log(htra) where htra is ADPAK rate parameter data for recombination
htlcx(0:htnt,0:htnn,0:htns-1)	_real	
	# log(htcx) where htcx is ADPAK rate parameter data for charge exchange
htlqa(0:htnt,0:htnn,0:htns-1)	_real	
	# log(htqa) where htqa is ADPAK rate parameter data for energy loss


***** Rtcrumpet:
# data from CRUMPET rate tables
cmpe	integer	/60/		# number of energy points for crm params
cmpd	integer /15/		# number of density points for crm params
crmdiss(cmpe,cmpd) _real [m**3/s] # Mol diss rate vs. temp and dens
crmarate(cmpe,cmpd) _real [m**3/s] # a-create rate from mol vs tem, dens
crmselm(cmpe,cmpd) _real [W m**3] # el e-change due to mol diss vs temp, dens
crmsiam(cmpe,cmpd) _real [W m**3] # ia e-change due to mol diss vs temp, dens
crmspotm(cmpe,cmpd) _real [W m**3] # Epot change due to mol diss vs temp, dens
crmsrada(cmpe,cmpd) _real [W m**3] # Radiation due to mol diss vs temp, dens
crmsradm(cmpe,cmpd) _real [W m**3] # Radiation due to mol diss vs temp, dens
cekpt(cmpe)	_real		# natural log of temperature (eV)
crlemin		real		# minimum of ekpt
crlemax		real		# maximum of ekpt
cerefmin	real	[eV]	# minimum temperature in table data
cerefmax	real	[eV]	# maximum temperature in table data
cdelekpt	real		# interval size for ekpt in table data
cdkpt(cmpd)	_real		# log10 of density (/m**3)
crldmin		real		# minimum of dkpt
crldmax		real		# maximum of dkpt
cdrefmin	real	[/m**3]	# minimum density in table data
cdrefmax	real	[/m**3]	# maximum density in table data
cdeldkpt	real		# dkpt interval in atomic data tables
ctaumin		real		# minimum tau in rtau opt-dep table data
ctaumax		real		# maximum tau in trau opt-dep table data
cdeltau	    real		# log int. size of rtau in opt-dep table dat

***** Rtdegas:
# data from DEGAS rate tables eh.dat, ehr1.dat, nwfits and atmc.dat
mpe	integer	/48/		# number of energy points
mpd	integer /11/		# number of density points
mpr     integer /1/		# number of optical-depth points
wsveh(mpe,mpd,mpr)  _real [m**3/s] # ioniz rate param vs temp and dens
wsveh0(mpe,mpd,mpr) _real [m**3/s] # recomb rate param vs temp and dens
wlemiss(mpe,mpd) _real [W m**3] # hydro. line rad. rate param vs temp, dens
welms(mpe,mpd)	 _real  [eV]    # elec. eng loss per ioniz. vs temp and dens
welms1(mpe,mpd,mpr) _real [J/sec] # elec rad ioniz-loss rate vs temp & dens
				  # Stotler's "coupling to ground state"
welms2(mpe,mpd,mpr) _real [J/sec] # elec rad recomb-loss rate vs temp & dens
			 	  # Stotler's "coupling to continuum"
pne3(mpe,mpd)	_real		# n=3 excited frac vs temp and dens; old DEGAS
pne31(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=3" popu
pne32(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=3" popu.
pne2(mpe,mpd)	_real		# n=2 excited frac vs temp and dens; old DEGAS
pne21(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=2" popu
pne22(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=2" popu
pne41(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=4" popu
pne42(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=4" popu.
pne51(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=5" popu
pne52(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=5" popu
pne61(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=6" popu
pne62(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=6" popu.
pne71(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=7" popu
pne72(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=7" popu
pne81(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=8" popu
pne82(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=8" popu.
pne91(mpe,mpd)	_real		# ground-state-coupl. (ng) coeff for "n=9" popu
pne92(mpe,mpd)	_real		# continuum-coupl. (ni) coeff for "n=9" popu
svdum2(mpe,mpe)	_real		# dummy array for un-used data in atmc.dat
svphcx(mpe,mpe)	_real  [m**3/s] # charge exc rate parameter vs Ti and neut eng
ekpt(mpe)	_real		# natural log of temperature (eV)
rlemin		real		# minimum of ekpt
rlemax		real		# maximum of ekpt
erefmin		real	[eV]	# minimum temperature in table data
erefmax		real	[eV]	# maximum temperature in table data
delekpt		real		# interval size for ekpt in table data
dkpt(mpd)	_real		# log10 of density (/m**3)
rldmin		real		# minimum of dkpt
rldmax		real		# maximum of dkpt
drefmin		real	[/m**3]	# minimum density in table data
drefmax		real	[/m**3]	# maximum density in table data
deldkpt		real		# dkpt interval in atomic data tables
taumin		real		# minimum tau in rtau opt-dep table data
taumax		real		# maximum tau in trau opt-dep table data
deltau		real		# log int. size of rtau in opt-dep table dat

***** Rtfcn:
# evaluation routines for ionization, recombination and charge-exchange rates
rsa(te:real,ne:real,r:real,k:integer)	real function	[m**3/sec]
	# rate parameter <sigma*v> for k--->k+1 ionization by electrons
	# in    te = electron temperature [J]
	# in	ne = electron density [/m**3]
        # in    r  = effective min. norm. neut. line-dens (~opacity) [m**-2]
	# in    k  = initial charge state
rra(te:real,ne:real,r:real,k:integer)	real function	[m**3/sec]
	# rate parameter <sigma*v> for k--->k-1 recombination
	# in	te = electron temperature [J]
	# in	ne = electron density [/m**3]
        # in    r  = effective min. norm. neut. line-dens (~opacity) [m**-2]
	# in    k  = initial charge state
rcx(t0:real,n0:real,k:integer)	real function	[m**3/sec]
	# rate param <sigma*v> for k--->k-1 charge exchange on neutral hydrogen
	# in	t0 = effective temperature per atomic mass [J/AMU]
	# in	n0 = density [/m**3] 
	# in    k  = initial charge state
rqa(te:real,ne:real,k:integer)	real function	[J*m**3/sec]
	# electron energy loss rate parameter <sigma*v>*dE for processes 
	# starting from charge state k, where dE is electron energy loss
	# in	te = electron temperature [J]
	# in	ne = electron density [/m**3]
	# in    k  = initial charge state
erl1(te:real,ne:real,r:real)	real function	[J/sec]
	# electron radiation loss rate per neutral H atom for "ionization"
	#   (D. Stotler's "coupling to the ground state")
	# in	te = electron temperature [J]
	# in	ne = electron density [/m**3]
        # in    r  = effective min. norm. neut. line-dens (~opacity) [m**-2]
erl2(te:real,ne:real,r:real)	real function	[J/sec]
	# electron radiation loss rate per H ion for "recombination"
	#   (D. Stotler's "coupling to the continuum")
	# in	te = electron temperature [J]
	# in	ne = electron density [/m**3]
        # in    r  = effective min. norm. neut. line-dens (~opacity) [m**-2]
svdiss(te:real)		real function	[m**3/s]
	# molecular dissociation by electrons (Janev's reaction 2.2.5)
	# in	te = electron temperature [J]
readrt(fname:string)	subroutine
	# read ADPAK (or STRAHL) file rates.adpak (or rates.strahl)
        # in    fname  filename to read
readeh(fname:string)	subroutine
	# read DEGAS file eh.dat
        # in    fname  filename to read
readnw(fname:string)	subroutine
	# read POST93 file nwfits
	# in    fname  filename to read
readehr1(fname:string)	subroutine
	# read DEGAS93 file ehr1.dat or ehr2.dat with n=2-3
	# in    fname  filename to read
readehr2(fname:string)	subroutine
	# read 2004 data (n=2-9) files thin.dat, thickLyA.dat, thickAllLy.dat
	# in    fname  filename to read
readcrumpete(fname:string)  subroutine
    # read CRUMPET energy sink/source file
    # in    fname  filename to read
readcrumpetn(fname:string)  subroutine
    # read CRUMPET rate file
    # in    fname  filename to read
sv_crumpet(te:real,ne:real,rate:integer) real function [X m**3/s]
    # calculate CRUMPET atomic creation form molecules
    # in    te=electron temperature [J]
    # in    ne=electron density [/m**3]
    # in    rate=data table to be read
    #            10: mol diss rate, [X]=1
    #            11: mol->atom rate, [X]=1
    #            20: el E change due to mol interaction, [X]=J
    #            21: i/a E change due to mol interaction, [X]=J
    #            22: Epot change due to mol interaction, [X]=J
    #            23: radiation due to mol interaction, [X]=J

***** Aphwrk:
# working arrays for 2-d spline interpolation
nxdata_aph  integer
nydata_aph  integer
xdata_aph(1:nxdata_aph)             _real
ydata_aph(1:nydata_aph)             _real
fdata_aph(1:nxdata_aph,1:nydata_aph)    _real
ldf_aph                          integer
iflag_aph                        integer
kxords_aph  integer /4/               # order of spline fit versus log(te)
                                  # kxords_aph=4 (default) is cubic interpolation
kyords_aph  integer /4/               # order of spline fit versus log10(ne)
                                  # kyords_aph=4 (default) is cubic interpolation
xknots_aph(1:nxdata_aph+kxords_aph)     _real
yknots_aph(1:nydata_aph+kyords_aph)     _real
workh(1:nxdata_aph*nydata_aph+2*kxords_aph*(nxdata_aph+1)) _real # work array
rsacoef(1:nxdata_aph,1:nydata_aph)  _real # spline coeff's for ionization
rracoef(1:nxdata_aph,1:nydata_aph)  _real # spline coeff's for recombination
rqacoef(1:nxdata_aph,1:nydata_aph)  _real # spline coeff's for line emission

***** Subs:
# Subroutines that can be called from the parser
aphread                                     subroutine
crumpetread                                 subroutine


