# This is a simple example exercising the methods in uefacets.py
from uefacets import *
from Numeric import *

print "Creating UEDGE instance"
ue = Uedge()
print "Reading input file 'testin.py' which is in the working directory"
print "(Note: alternatively give the full path)"
ue.ReadParams("testin.py")
print "Set up MPI. This does nothing for now since UEDGE is running serial"
ue.SetComm()
print "Set some boundary conditions different from those in testin.py"
ni = array([3.e19])
ng = array([1.e15])
ti = array([100.])
ui = array([0.])
print "Specifically, ni = ",ni[0], ", ng = ",ng[0],","
print "  Te=Ti=", ti[0], ", ui = ",ui[0]
ue.SetData(ni,ui,ng,ti,ti,"val","val","val","val")
print "Take a step to time t=0.001"
niflux,uiflux,ngflux,tiflux,teflux = ue.Advance(.001)
print "The output fluxes are:"
print "niflux = ",niflux
print "uiflux = ",uiflux
print "ngflux = ",ngflux
print "tiflux = ",tiflux
print "teflux = ",teflux

print "Let's not accept this result and reset"
ue.Reset()
print "Change n_i to 2e19 and retake step"
ni = array([2.e19])
ue.SetData(ni,ui,ng,ti,ti,"val","val","val","val")

niflux,uiflux,ngflux,tiflux,teflux = ue.Advance(.001)
print "The output fluxes are:"
print "niflux = ",niflux
print "uiflux = ",uiflux
print "ngflux = ",ngflux
print "tiflux = ",tiflux
print "teflux = ",teflux

print "create dumpfiles dump.pdb and dump.hdf"
print "Note this requires that data to be dumped has been flagged in dot v files"
ue.Dump("dump")

print "Illustrating use of 'DoCommand' to print the sum of tiflux and teflux"
ue.DoCommand("print 'tiflux+teflux = ',tiflux+teflux")
