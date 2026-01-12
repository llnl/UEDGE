# This input file sets the parameters for the parallel test case
#
from uedge import *

# flags also set in the workstation version.
bbb.mhdgeo=-1
bbb.isnion=1
bbb.isupon=1
bbb.isteon=1
bbb.istion=1
bbb.isngon=0
bbb.svrpkg="vodpk"
bbb.premeth="ilut"
bbb.epscon1=1.e-2
bbb.ireorder=0
bbb.ncore=2.e19
bbb.tcoree=100
bbb.tcorei=100
bbb.isfixlb=2
bbb.recycp=.9
bbb.runtim=1.e-7
##bbb.trange=4.e6
bbb.trange=4.e3
bbb.nsteps=30
bbb.n0g=1.e16
bbb.difni=1.
bbb.kye=1.
bbb.flalfe=0.21
bbb.flalfi=0.21
bbb.flalfgx=1.e10
bbb.flalfgy=1.e10
com.nycore=0
com.nysol=10
com.nxleg[0,0]=0
com.nxleg[0,1]=2
com.nxcore[0,0]=0
com.nxcore[0,1]=4
grd.zax=1.
grd.zaxpt=.75
grd.alfyt=-1.e-5
print "Finished setting variables"

print "Allocate Storage."
bbb.allocate ()
bbb.restart=0

bbb.exmain()

