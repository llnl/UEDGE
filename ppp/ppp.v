ppp
{
}

***** ParallelSettings:
OMPParallelPandf1 integer /0/   # [0]: serial pandf1 rhs calc [1] omp parallel pandf1 rhs calc
OMPParallelJac    integer /0/   # [0]: serial jacobian calc [1] omp parallel jacobian calc
ParallelWarning   integer /0/   # Warning for users who wish to use it
CheckJac          integer /0/   # [0/1]: Turn on on-the-fly comparison of parallel vs serial 
                                # evaluation of Jacobian. If differences between para and 
                                # serial Jacobians, dump both Jacs in serialjac.dat and 
                                # paralleljac.dat with routine jac_write in current working 
                                # folder.
Nthreads          integer /64/  # Number of threads to be used to calculate the Jacobian
CheckPandf1       integer /0/   # [0/1]: Turn on on-the-fly comparison of parallel vs 
                                # serial evaluation of pandf1.

***** ParallelDebug:
OMPJacDebug       integer   /0/     # Print debug info for omp constructs
iidebugprint      integer   /-1/    # Index ii of jacobian dyldot(ii)/yl(iv) at which 
                                    # threadprivate variables are printed after calculation of the 
                                    # jacobian element. iv is determined by ivdebugprint
ivdebugprint      integer   /-1/    # index iv of jacobian dyldot(ii)/yl(iv) at which threadprivate 
                                    # variables are printed after calculation of the jacobian element. 
                                    # ii is determined by iidebugprint
ForceSerialCheck  integer   /0/     # [0/1]: Force two sequential serial evaluations of the Jacobian to
                                    # verify that Jacobian evaluation is reproducible (fail when e.g.
                                    # variables are not properly initialized in pandf).
DebugJac          integer   /0/
DumpJac           integer   /0/     # [0/1]: Turn on dumping of data for the diverging element of serial
                                    # and parallel jacobian (only available when CheckJac is on). 
DumpFullJac       integer   /0/     # [0/1]: Turn on dumping of full serial jacobian for analysis of 
                                    #bandwidth (dumping in file ). 
WriteJacobian     integer   /0/     # Write jacobian in an ascii text file
OMPCopyArray      integer   /1/     # For Debug purpose: turn on/off(0/1) copy of threadprivate arrays 
                                    # before jacobian calculation (WARNING:could cause numerical inacurarry 
                                    # if turned on)
OMPCopyScalar     integer    /1/    # For Debug purpose: turn on/off copy(0/1) of threadprivate scalar 
                                    # before jacobian calculation (WARNING:could cause numerical 
                                    #inacurracy if turned on)

***** OMPJacSettings:
OMPJacNchunks   integer     /0/     # Number of chunks to be used to calculate the Jacobian. If 
                                    # Nchunks.lt.0, Nchunks=Nthreads elif Nchunks=0, Nchunks=neq
OMPJacVerbose   integer     /0/     # Print info for omp jacobian calculation
OMPlenpfac      integer     /5/     # Factor to increase nnzmxperchunk
OMPCheckNaN     integer     /0/     # Check whether jacobian terms are NaN after jacobian calculation
OMPLoadBalance  integer     /0/     # Enable user defined weights for each OMP tasks (overrided by 
                                    # MPIAutoBalance)
OMPAutoBalance  integer     /1/     # Automatic load balancing for OMP thread tasks (if OMPLoadWeight=)
OMPBalanceStrength  real    /1.0/   # Strenght s of the load balance 
                                    # (Loadweight=Loadweight*(t_thread/<t_thread>)**s)
OMPTimingJacRow     integer /0/     # Profile execution time of calculatation of each row of the jacobian
OMPLoopJacNchunk    integer /1/     
OMPJacStamp     character*20    /"*OMPJac* "/ 
                                    # Stamp for hybrid output (not an user input)

***** OMPJac:
NchunksJac      integer     /1/
nnzmxperchunk   integer             # Maximum number of jacobian elements which 
                                    # can be stored per thread. Can be increased 
                                    # with omplenpfac
OMPivmin(NchunksJac)    _integer    # jacobian rows with 
                                    # ivmin(ithread)<=iv<=ivmax(ithread) are calculated
                                    # on thread ithread
OMPivmax(NchunksJac)    _integer    # jacobian rows with ivmin(ithread)<=iv<=ivmax(ithread) 
                                    # are calculated on thread ithread
OMPLoadWeight(1:NchunksJac)     _real   # weight for load distribution of jacobian 
                                        # calculation among threads
OMPTimeLocalJac(1:NchunksJac)   _real   # runtime for jac calculation on each threads. 
                                        # Used to optimize load distribution of 
                                        # jacobian calculation among threads when 
                                        # AutoBalance=1
iJacRow(neq)    _integer            #
OMPTimeJacRow(neq)      _real       #
nnz(NchunksJac)         _integer    # 
nnzcum(NchunksJac)      _integer    #
iJacCol(nnzmxperchunk,NchunksJac)   _integer    #
rJacElem(nnzmxperchunk,NchunksJac)  _real       #

**** OMPPandf1Settings:
OMPPandf1Stamp      character*20    /"*OMPPandf1* "/ #
OMPPandf1Debug      integer         /0/
OMPPandf1Verbose    integer         /0/
OMPTimeParallelPandf1   real        /0.0/
OMPTimeSerialPandf1     real        /0.0/
OMPPandf1LoopNchunk     integer     /1/
OMPPandf1Nychunks   integer         /0/
OMPPandf1Nxchunks   integer         /1/
xpadding            integer         /2/
ypadding            integer         /2/

**** OMPPandf1:
isnionxy_old(0:nx+1,0:ny+1,nisp) _integer /0/
isngonxy_old(0:nx+1,0:ny+1,ngsp) _integer /0/
isuponxy_old(0:nx+1,0:ny+1,nisp) _integer /0/
istionxy_old(0:nx+1,0:ny+1) _integer /0/
isteonxy_old(0:nx+1,0:ny+1) _integer /0/
istgonxy_old(0:nx+1,0:ny+1,ngsp) _integer /0/
isphionxy_old(0:nx+1,0:ny+1) _integer /0/
nisp_old integer /0/
ngsp_old integer /0/
nx_old integer /0/
ny_old integer /0/
Nxchunks_old integer /0/
Nychunks_old integer /0/
Nxptchunks_old integer /0/
neq_old             integer     /0/


chunks(neq,3)       _integer
Nchunks             integer
Nychunks            integer     /0/
Nxchunks            integer     /1/
NchunksPandf1       integer     /1/
Nchunksmax          integer     /1/
Nixychunksmax       integer     /1/
yincchunk(NchunksPandf1)    _integer
xincchunk(NchunksPandf1)    _integer
ixchunk(NchunksPandf1)      _integer
iychunk(NchunksPandf1)      _integer
rangechunk(NchunksPandf1, 4)  _integer
Nivchunk(NchunksPandf1)     _integer
ivchunk(NchunksPandf1,Nchunksmax)  _integer
Nixychunk(NchunksPandf1)    _integer
ixychunk(NchunksPandf1,Nixychunksmax,2)  _integer
iymaxchunk(NchunksPandf1)   _integer
ixmaxchunk(NchunksPandf1)   _integer
iyminchunk(NchunksPandf1)   _integer
ixminchunk(NchunksPandf1)   _integer

**** OMPTiming:
DebugTime       integer /0/     # Display execution times of various subroutines
ShowTime        integer /1/     # Show execution time of routines
SerialDebug     integer /0/     # Show execution time of routines
ParaTime        real    /0./    
SerialTime      real    /0./

#OMPTotJacCalc      real  /0./ # time to calculate jacobian in jac_calc_omp
#OMPTotTimeCollect  real  /0./ # time to collect jacobian elements in jac_calc_omp
#OMPTotTimeBuild    real  /0./ # time to calculate elements of jacobian in jac_calc_omp

***** JacDebug:
EvalDumpJac(FileName:string) subroutine # Dumps jacobian to file for further processing


***** Subs:
Make2DChunks(i,j)                           subroutine

***** PandfCopies:
up_cp(0:nx+1,0:ny+1,1:nisp)                    _real
ne_cp(0:nx+1,0:ny+1)                           _real
phi_cp(0:nx+1,0:ny+1)                          _real
ti_cp(0:nx+1,0:ny+1)                           _real
tg_cp(0:nx+1,0:ny+1,1:ngsp)                    _real
te_cp(0:nx+1,0:ny+1)                           _real
nit_cp(0:nx+1,0:ny+1)                          _real
nz2_cp(0:nx+1,0:ny+1)                          _real
ng_cp(0:nx+1,0:ny+1,1:ngsp)                    _real
nm_cp(0:nx+1,0:ny+1,1:nisp)                    _real
lng_cp(0:nx+1,0:ny+1,1:ngsp)                   _real
ni_cp(0:nx+1,0:ny+1,1:nisp)                    _real
pg_cp(0:nx+1,0:ny+1,1:ngsp)                    _real
pr_cp(0:nx+1,0:ny+1)                           _real
pre_cp(0:nx+1,0:ny+1)                          _real
pri_cp(0:nx+1,0:ny+1,1:nisp)                   _real
pgy0_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
tgy0_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
gpry_cp(0:nx+1,0:ny+1)                         _real
niy1_cp(0:nx+1,0:ny+1,1:nisp)                  _real
nity1_cp(0:nx+1,0:ny+1)                        _real
gpix_cp(0:nx+1,0:ny+1,1:nisp)                  _real
ex_cp(0:nx+1,0:ny+1)                           _real
ey_cp(0:nx+1,0:ny+1)                           _real
tiy1_cp(0:nx+1,0:ny+1)                         _real
phiy0_cp(0:nx+1,0:ny+1)                        _real
tiy0s_cp(0:nx+1,0:ny+1)                        _real
gtex_cp(0:nx+1,0:ny+1)                         _real
ney0_cp(0:nx+1,0:ny+1)                         _real
gtiy_cp(0:nx+1,0:ny+1)                         _real
zeff_cp(0:nx+1,0:ny+1)                         _real
tiy0_cp(0:nx+1,0:ny+1)                         _real
pgy1_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
phiy0s_cp(0:nx+1,0:ny+1)                       _real
phiy1_cp(0:nx+1,0:ny+1)                        _real
ngy0_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
gpiy_cp(0:nx+1,0:ny+1,1:nisp)                  _real
priy0_cp(0:nx+1,0:ny+1,1:nisp)                 _real
gtey_cp(0:nx+1,0:ny+1)                         _real
gtix_cp(0:nx+1,0:ny+1)                         _real
phiv_cp(0:nx+1,0:ny+1)                         _real
niy0s_cp(0:nx+1,0:ny+1,1:nisp)                 _real
tey1_cp(0:nx+1,0:ny+1)                         _real
tiy1s_cp(0:nx+1,0:ny+1)                        _real
priy1_cp(0:nx+1,0:ny+1,1:nisp)                 _real
tgy1_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
gpex_cp(0:nx+1,0:ny+1)                         _real
niy0_cp(0:nx+1,0:ny+1,1:nisp)                  _real
tiv_cp(0:nx+1,0:ny+1)                          _real
priv_cp(0:nx+1,0:ny+1,1:nisp)                  _real
ney1_cp(0:nx+1,0:ny+1)                         _real
gpondpotx_cp(0:nx+1,0:ny+1)                    _real
tev_cp(0:nx+1,0:ny+1)                          _real
phiy1s_cp(0:nx+1,0:ny+1)                       _real
prtv_cp(0:nx+1,0:ny+1)                         _real
prev_cp(0:nx+1,0:ny+1)                         _real
tey0_cp(0:nx+1,0:ny+1)                         _real
gprx_cp(0:nx+1,0:ny+1)                         _real
gpey_cp(0:nx+1,0:ny+1)                         _real
znot_cp(0:nx+1,0:ny+1)                         _real
niy1s_cp(0:nx+1,0:ny+1,1:nisp)                 _real
nity0_cp(0:nx+1,0:ny+1)                        _real
ngy1_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
dutm_use_cp(0:nx+1,0:ny+1,1:nisp)              _real
difp_use_cp(0:nx+1,0:ny+1,1:nisp)              _real
dif_use_cp(0:nx+1,0:ny+1,1:nisp)               _real
kyi_use_cp(0:nx+1,0:ny+1)                      _real
trax_use_cp(0:nx+1,0:ny+1,1:nisp)              _real
kxbohm_cp(0:nx+1,0:ny+1)                       _real
kxe_use_cp(0:nx+1,0:ny+1)                      _real
kxi_use_cp(0:nx+1,0:ny+1)                      _real
kye_use_cp(0:nx+1,0:ny+1)                      _real
vy_use_cp(0:nx+1,0:ny+1,1:nisp)                _real
betap_cp(0:nx+1,0:ny+1)                        _real
kybohm_cp(0:nx+1,0:ny+1)                       _real
tray_use_cp(0:nx+1,0:ny+1,1:nisp)              _real
dif2_use_cp(0:nx+1,0:ny+1,1:nisp)              _real
eta1_cp(0:nx+1,0:ny+1)                         _real
ctaui_cp(0:nx+1,0:ny+1,nisp)                   _real
dclass_e_cp(0:nx+1,0:ny+1)                     _real
loglambda_cp(0:nx+1,0:ny+1)                    _real
rtaue_cp(0:nx+1,0:ny+1)                        _real
ctaue_cp(0:nx+1,0:ny+1,nisp)                   _real
dclass_i_cp(0:nx+1,0:ny+1)                     _real
vydd_cp(0:nx+1,0:ny+1,1:nisp)                  _real
veycb_cp(0:nx+1,0:ny+1)                        _real
coll_fe_cp(0:nx+1,0:ny+1)                      _real
vyte_cft_cp(0:nx+1,0:ny+1)                     _real
vyrd_cp(0:nx+1,0:ny+1,1:nisp)                  _real
vycf_cp(0:nx+1,0:ny+1)                         _real
vygp_cp(0:nx+1,0:ny+1,1:nisp)                  _real
diffusivwrk_cp(0:nx+1,0:ny+1)                  _real
vy_cp(0:nx+1,0:ny+1,1:nisp)                    _real
vycr_cp(0:nx+1,0:ny+1)                         _real
vyce_cp(0:nx+1,0:ny+1,1:nisp)                  _real
coll_fi_cp(0:nx+1,0:ny+1)                      _real
veycp_cp(0:nx+1,0:ny+1)                        _real
vyti_cft_cp(0:nx+1,0:ny+1)                     _real
vycb_cp(0:nx+1,0:ny+1,1:nisp)                  _real
vycp_cp(0:nx+1,0:ny+1,1:nisp)                  _real
vy_cft_cp(0:nx+1,0:ny+1,1:nisp)                _real
v2dd_cp(0:nx+1,0:ny+1,1:nisp)                  _real
ve2cb_cp(0:nx+1,0:ny+1)                        _real
vytan_cp(0:nx+1,0:ny+1,1:nisp)                 _real
v2rd_cp(0:nx+1,0:ny+1,1:nisp)                  _real
fdiaxlb_cp(0:ny+1,1:nxpt)                      _real
ve2cd_cp(0:nx+1,0:ny+1,1:nisp)                 _real
v2xgp_cp(0:nx+1,0:ny+1,1:nisp)                 _real
v2cd_cp(0:nx+1,0:ny+1,1:nisp)                  _real
vyavis_cp(0:nx+1,0:ny+1,1:nisp)                _real
v2ce_cp(0:nx+1,0:ny+1,1:nisp)                  _real
q2cd_cp(0:nx+1,0:ny+1,1:nisp)                  _real
v2cb_cp(0:nx+1,0:ny+1,1:nisp)                  _real
v2_cp(0:nx+1,0:ny+1,1:nisp)                    _real
fdiaxrb_cp(0:ny+1,1:nxpt)                      _real
fq2d_cp(0:nx+1,0:ny+1)                         _real
fmity_cp(0:nx+1,0:ny+1,1:nisp)                 _real
fqym_cp(0:nx+1,0:ny+1)                         _real
fqyai_cp(0:nx+1,0:ny+1)                        _real
fqyb_cp(0:nx+1,0:ny+1)                         _real
fqydti_cp(0:nx+1,0:ny+1,1:nisp)                _real
fqymi_cp(0:nx+1,0:ny+1,1:nisp)                 _real
fqydt_cp(0:nx+1,0:ny+1)                        _real
fqyao_cp(0:nx+1,0:ny+1)                        _real
fqya_cp(0:nx+1,0:ny+1)                         _real
fqyae_cp(0:nx+1,0:ny+1)                        _real
fqygp_cp(0:nx+1,0:ny+1)                        _real
fq2_cp(0:nx+1,0:ny+1)                          _real
fqyd_cp(0:nx+1,0:ny+1)                         _real
fqy_cp(0:nx+1,0:ny+1)                          _real
netap_cp(0:nx+1,0:ny+1)                        _real
fqp_cp(0:nx+1,0:ny+1)                          _real
dphi_iy1_cp(0:nx+1)                            _real
fqpsatrb_cp(0:ny+1,2)                          _real
fqx_cp(0:nx+1,0:ny+1)                          _real
fqpsatlb_cp(0:ny+1,2)                          _real
fqxb_cp(0:nx+1,0:ny+1)                         _real
upi_cp(0:nx+1,0:ny+1,1:nisp)                   _real
uz_cp(0:nx+1,0:ny+1,1:nisp)                    _real
uu_cp(0:nx+1,0:ny+1,1:nisp)                    _real
frici_cp(0:nx+1,0:ny+1,nisp)                   _real
uup_cp(0:nx+1,0:ny+1,1:nisp)                   _real
frice_cp(0:nx+1,0:ny+1)                        _real
vey_cp(0:nx+1,0:ny+1)                          _real
vex_cp(0:nx+1,0:ny+1)                          _real
upe_cp(0:nx+1,0:ny+1)                          _real
psorrgc_cp(0:nx+1,0:ny+1,1:ngsp)               _real
psorgc_cp(0:nx+1,0:ny+1,1:ngsp)                _real
nuelg_cp(0:nx+1,0:ny+1,ngsp)                   _real
nurc_cp(0:nx+1,0:ny+1,ngsp)                    _real
psordis_cp(0:nx+1,0:ny+1,1:nisp)               _real
psorbgg_cp(0:nx+1,0:ny+1,1:ngsp)               _real
psordisg_cp(0:nx+1,0:ny+1,1:ngsp)              _real
smoc_cp(0:nx+1,0:ny+1,1:nusp)                  _real
nuvl_cp(0:nx+1,0:ny+1,nisp)                    _real
psorc_cp(0:nx+1,0:ny+1,1:nisp)                 _real
nuiz_cp(0:nx+1,0:ny+1,ngsp)                    _real
nucxi_cp(0:nx+1,0:ny+1,nisp)                   _real
snic_cp(0:nx+1,0:ny+1,1:nisp)                  _real
seic_cp(0:nx+1,0:ny+1)                         _real
psor_cp(0:nx+1,0:ny+1,1:nisp)                  _real
rtauy_cp(0:nx+1,0:ny+1)                        _real
psorbgz_cp(0:nx+1,0:ny+1)                      _real
seec_cp(0:nx+1,0:ny+1)                         _real
msor_cp(0:nx+1,0:ny+1,1:nisp)                  _real
psorg_cp(0:nx+1,0:ny+1,1:ngsp)                 _real
rtaux_cp(0:nx+1,0:ny+1)                        _real
psorcxg_cp(0:nx+1,0:ny+1,1:ngsp)               _real
psorrg_cp(0:nx+1,0:ny+1,1:ngsp)                _real
msorxr_cp(0:nx+1,0:ny+1,1:nisp)                _real
rtau_cp(0:nx+1,0:ny+1)                         _real
nucx_cp(0:nx+1,0:ny+1,ngsp)                    _real
nueli_cp(0:nx+1,0:ny+1,nisp)                   _real
psorxr_cp(0:nx+1,0:ny+1,1:nisp)                _real
psorxrc_cp(0:nx+1,0:ny+1,1:nisp)               _real
nuix_cp(0:nx+1,0:ny+1,ngsp)                    _real
floxg_cp(0:nx+1,0:ny+1)                        _real
conyg_cp(0:nx+1,0:ny+1)                        _real
vyg_cp(0:nx+1,0:ny+1,1:ngsp)                   _real
fngy4ord_cp(0:nx+1,0:ny+1,1:ngsp)              _real
fngy_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
fngxy_cp(0:nx+1,0:ny+1,1:ngsp)                 _real
conxg_cp(0:nx+1,0:ny+1)                        _real
fngx_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
vygtan_cp(0:nx+1,0:ny+1,1:ngsp)                _real
floyg_cp(0:nx+1,0:ny+1)                        _real
uug_cp(0:nx+1,0:ny+1,1:ngsp)                   _real
fngx4ord_cp(0:nx+1,0:ny+1,1:ngsp)              _real
uuxg_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
alfneo_cp(0:nx+1,0:ny+1,1:nisp)                _real
visy_cp(0:nx+1,0:ny+1,1:nisp)                  _real
ktneo_cp(0:nx+1,0:ny+1,1:nisp)                 _real
visxneo_cp(0:nx+1,0:ny+1,1:nisp)               _real
nuiistar_cp(0:nx+1,0:ny+1,1:nisp)              _real
nuii_cp(0:nx+1,0:ny+1,1:nisp)                  _real
w_cp(0:nx+1,0:ny+1)                            _real
visx_cp(0:nx+1,0:ny+1,1:nisp)                  _real
k2neo_cp(0:nx+1,0:ny+1,1:nisp)                 _real
hcxineo_cp(0:nx+1,0:ny+1)                      _real
w2_cp(0:nx+1,0:ny+1)                           _real
hcyij_cp(0:nx+1,0:ny+1,1:nisp)                 _real
hcxij_cp(0:nx+1,0:ny+1,1:nisp)                 _real
hcyi_cp(0:nx+1,0:ny+1)                         _real
hcyn_cp(0:nx+1,0:ny+1)                         _real
hcye_cp(0:nx+1,0:ny+1)                         _real
qipar_cp(0:nx+1,0:ny+1,nisp)                   _real
hcxi_cp(0:nx+1,0:ny+1)                         _real
w1_cp(0:nx+1,0:ny+1)                           _real
hcxe_cp(0:nx+1,0:ny+1)                         _real
hcxn_cp(0:nx+1,0:ny+1)                         _real
eqp_cp(0:nx+1,0:ny+1)                          _real
hcxg_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
hcyg_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
floxge_cp(0:nx+1,0:ny+1,1:ngsp)                _real
segc_cp(0:nx+1,0:ny+1,1:ngsp)                  _real
floyge_cp(0:nx+1,0:ny+1,1:ngsp)                _real
conxge_cp(0:nx+1,0:ny+1,1:ngsp)                _real
conyge_cp(0:nx+1,0:ny+1,1:ngsp)                _real
fegx_cp(0:nx+1,0:ny+1,ngsp)                    _real
fegxy_cp(0:nx+1,0:ny+1,ngsp)                   _real
fegy_cp(0:nx+1,0:ny+1,ngsp)                    _real
fniy_cp(0:nx+1,0:ny+1,1:nisp)                  _real
fniy4ord_cp(0:nx+1,0:ny+1,1:nisp)              _real
fnixcb_cp(0:nx+1,0:ny+1,1:nisp)                _real
fniycbo_cp(0:nx+1,1:nisp)                      _real
fniycb_cp(0:nx+1,0:ny+1,1:nisp)                _real
fnix_cp(0:nx+1,0:ny+1,1:nisp)                  _real
conx_cp(0:nx+1,0:ny+1,nusp)                    _real
floy_cp(0:nx+1,0:ny+1,nusp)                    _real
flox_cp(0:nx+1,0:ny+1,nusp)                    _real
cony_cp(0:nx+1,0:ny+1,nusp)                    _real
fmivxpt_cp(1:nusp,1:nxpt)                      _real
fmihxpt_cp(1:nusp,1:nxpt)                      _real
fmixy_cp(0:nx+1,0:ny+1,1:nusp)                 _real
vyvxpt_cp(1:nusp,1:nxpt)                       _real
nixpt_cp(1:nusp,1:nxpt)                        _real
vyhxpt_cp(1:nusp,1:nxpt)                       _real
visyxpt_cp(1:nusp,1:nxpt)                      _real
upxpt_cp(1:nusp,1:nxpt)                        _real
fmix_cp(0:nx+1,0:ny+1,1:nusp)                  _real
fmiy_cp(0:nx+1,0:ny+1,1:nusp)                  _real
wvh_cp(0:nx+1,0:ny+1,1:nusp)                   _real
prad_cp(0:nx+1,0:ny+1)                         _real
pwrzec_cp(0:nx+1,0:ny+1)                       _real
eeli_cp(0:nx+1,0:ny+1)                         _real
floxibgt_cp(0:nx+1,0:ny+1,1:nisp)              _real
vsoree_cp(0:nx+1,0:ny+1)                       _real
pwrze_cp(0:nx+1,0:ny+1)                        _real
edisse_cp(0:nx+1,0:ny+1)                       _real
feexy_cp(0:nx+1,0:ny+1)                        _real
feeycbo_cp(0:nx+1)                             _real
feey4ord_cp(0:nx+1,0:ny+1)                     _real
seak_cp(0:nx+1,0:ny+1)                         _real
seik_cp(0:nx+1,0:ny+1)                         _real
emolia_cp(0:nx+1,0:ny+1,1:nisp)                _real
ntau_cp(0:nx+1,0:ny+1)                         _real
pradcff_cp(0:nx+1,0:ny+1)                      _real
pwribkg_cp(0:nx+1,0:ny+1)                      _real
psicx_cp(0:nx+1,0:ny+1)                        _real
floye_cp(0:nx+1,0:ny+1)                        _real
na_cp(0:nx+1,0:ny+1)                           _real
feixy_cp(0:nx+1,0:ny+1)                        _real
pradc_cp(0:nx+1,0:ny+1)                        _real
pradzc_cp(0:nx+1,0:ny+1,0:nzspmx,1:nzspmx+1)   _real
pradhyd_cp(0:nx+1,0:ny+1)                      _real
seid_cp(0:nx+1,0:ny+1)                         _real
conxi_cp(0:nx+1,0:ny+1)                        _real
floxe_cp(0:nx+1,0:ny+1)                        _real
sead_cp(0:nx+1,0:ny+1)                         _real
conxe_cp(0:nx+1,0:ny+1)                        _real
seit_cp(0:nx+1,0:ny+1)                         _real
pwrebkg_cp(0:nx+1,0:ny+1)                      _real
pradz_cp(0:nx+1,0:ny+1,0:nzspmx,1:nzspmx+1)    _real
feey_cp(0:nx+1,0:ny+1)                         _real
seg_ue_cp(0:nx+1,0:ny+1,1:nfl)                 _real
vsoreec_cp(0:nx+1,0:ny+1)                      _real
floyi_cp(0:nx+1,0:ny+1)                        _real
seidh_cp(0:nx+1,0:ny+1)                        _real
seadh_cp(0:nx+1,0:ny+1)                        _real
conye_cp(0:nx+1,0:ny+1)                        _real
conyi_cp(0:nx+1,0:ny+1)                        _real
feiy4ord_cp(0:nx+1,0:ny+1)                     _real
feiycbo_cp(0:nx+1)                             _real
floxebgt_cp(0:nx+1,0:ny+1)                     _real
feex_cp(0:nx+1,0:ny+1)                         _real
nratio_cp(0:nx+1,0:ny+1)                       _real
feiy_cp(0:nx+1,0:ny+1)                         _real
feix_cp(0:nx+1,0:ny+1)                         _real
floxi_cp(0:nx+1,0:ny+1)                        _real
erliz_cp(0:nx+1,0:ny+1)                        _real
w0_cp(0:nx+1,0:ny+1)                           _real
nzloc_cp(0:nzspmx)                             _real
eqpg_cp(0:nx+1,0:ny+1,ngsp)                    _real
eiamoldiss_cp(0:nx+1,0:ny+1,1:nisp)            _real
sng_ue_cp(0:nx+1,0:ny+1,1:nfl)                 _real
resco_cp(0:nx+1,0:ny+1,1:nisp)                 _real
resng_cp(0:nx+1,0:ny+1,1:ngsp)                 _real
fricnrl_cp(0:nx+1,0:ny+1,nusp)                 _real
resmo_cp(0:nx+1,0:ny+1,1:nusp)                 _real
reseg_cp(0:nx+1,0:ny+1,1:ngsp)                 _real
resei_cp(0:nx+1,0:ny+1)                        _real
resee_cp(0:nx+1,0:ny+1)                        _real
resphi_cp(0:nx+1,0:ny+1)                       _real
