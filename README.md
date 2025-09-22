# UEDGE
A 2D time-dependent fluid simulation code of plasma and neutrals in magnetic fusion devices.
## Brief description 
UEDGE is an interactive suite of physics packages using the Python or BASIS scripting systems. 
The original (circa 2007) Python version was further developed under the FACETS SciDAC project 
[Cary 2008, McCourt 2012]. The plasma is described by time-dependent 2D plasma fluid equations 
that include equations for density, velocity, ion temperature, electron temperature, electrostatic 
potential, and gas density in the edge region of a magnetic fusion energy confinement device. slab, 
cylindrical, and toroidal geometries are allowed, and closed and open magnetic field-line regions 
are included. Classical transport is assumed along magnetic field lines, and anomalous transport 
is assumed across field lines.  Multi-charge state impurities can be included with the corresponding 
line-radiation energy loss. 

## Method of solution
A fully implicit numerical algorithm is used that allows both Newton-like iterations to steady state 
and time-dependent solutions with large time-steps.  A preconditioning matrix is obtained by approximate 
(ILUT) inversion of a numerical finite-difference Jacobian, which is then used in a Newton-Krylov 
solution algorithm. A finite-volume differencing algorithm is used. Over 95% of the coding is in 
Fortran with the remainder being C.

## Related and auxiliary software

Although UEDGE is written in Fortran, for efficient execution and analysis of results, it utilizes 
either Python or BASIS scripting shells. Python is easily available for many platforms 
(http://www.Python.org/). The features and availability of BASIS are described in "Basis Manual Set" 
by P.F. Dubois, Z.C. Motteler, et al., Lawrence Livermore National Laboratory report UCRL-MA-118541, 
June, 2002 and http://basis.llnl.gov/), however, BASIS is deprecated. Contact one of the UEDGE developers
if you insist on running it within that environment. 
The Python version of UEDGE uses the same source files but utilizes Forthon to produce a Python-compatible 
source.  Forthon has been developed by D.P. Grote (see http://hifweb.lbl.gov/Forthon/ and Grote et al. 
in the references below), and it is freely available. The graphics can be performed by any package importable 
to Python, such as PYGIST.  The parallel version of UEDGE available through Python also uses the PETSc linear 
algebra solver whose development has been led by ANL (https://www.mcs.anl.gov/petsc/).

UEDGE can also be coupled to other codes.  An excellent example is couplling to the DUSTT code from
UCSD (contact rsmirnov@eng.ucsd.edu) that follows the trajectories and ablation of dust particles in
the background UEDGE plasma and provides impurity sources to UEDGE.  For an example, see R. Smirnov
et al., Phys. Plasmas 22 (2015) 012506.

## Getting started 
### Installing UEDGE using PIP
The Python version of UEDGE can be installed using pip:

    pip install uedge

To install a specific UEDGE version using pip:

    pip install uedge==X.X.X

where X.X.X is the version number. For a list of the available UEDGE versions:

    pip install uedge==

### Installing UEDGE from the source
PIP is also used to install UEDGE from the source files:

    pip install [-v] /path/to/uedge

where the optional -v makes the build verbose.

Building UEDGE from the source requires fopenmp to be available to the
compiler, unless a serial build is specified. To ensure the right compiler 
is linked, set the CC environment variable:

    export CC=path/to/clang

For most OSX builds, the llvm clang installation is recommended. 

There are a number of build options available for building the code:

  serial - builds the code without OMP parallel capabilities
  noclean - does not clean the build directory between builds
  checkbounds - turns on boundary checks for non-debug builds 
  debug - turns on compiler debugging flags

The options are passed to the build using the following syntax

    pip install -v . -C "--global-option=--arg1 --arg2"

where arg1 and/or arg2 are one of the options from the above list.

## Using the OMP parallel version of UEDGE
The pip-distributed UEDGE packages starting from V8.3 have OMP capabilities.

The following runtime environment variables _must_ be defined before using
the OMP capabilities to avoid running out of memory on the threads

    ulimit -s unlimited
    export OMP_STACKSIZE=2G

In case of unexpected crashes, increase OMP_STACKSIZE further. To control the 
number of processes, set

    export OMP_NUM_THREADS=32


The following variables are used to turn on the OMP capabilities for UEDGE:
  ppp.OMPParallelJac - =1 activates OMP paralellization of the Jacobian
  ppp.OMPParallelPandf1 - =1 activates OMP paralellization of the pandf1 routine

Additional information on the OMP evaluation can be displayed by setting
  ppp.ParallelWarning - =1 outputs OMP parallelization information
  ppp.OMPPandf1Verbose - =[1,2] outputs diagnostics on parallel pandf1 evaluation
  ppp.OMPPandf1Verbose - =[1,2] outputs diagnostics on parallel pandf1 evaluation
  ppp.OMPJacVerbose - =[1,2] outputs diagnostics on parallel pandf1 evaluation

The OMP settings are controlled using the variables:
  ppp.OMPlenpfac - size allocated to each process
  ppp.OMPJacNchunks - number of chunks Jacobian is split into (defaults to neq)
  ppp.OMPPandf1Nychunks - number of chunks to decompose domain into (defaults to ny)

Checks on the OMP parallelization for code development testing can be turned
on using:
  bbb.isjacreset - checks on subsequent parallel Jacobian (computes Jacobian twice)
  ppp.CheckJac - =1 compares parallel Jacobian to serial (runs both)
  ppp.CheckPandf1 - =1 compares parallel pandf1 to serial (runs both)

## How to get involved and contribute
Sent email to one of the developers listed below expressing your interest in modifying or developing packages 
for UEDGE.  Either new or improved physics models or numerical algorithms are most welcome.

## Authors contributing to V7 release
T.D. Rognlien, I. Joseph, W.H. Meyer, M.E. Rensink, and M.V. Umansky, LLNL  
(trognlien@llnl.gov, joseph5@llnl.gov, meyer8@llnl.gov, rensink1@llnl.gov, umansky1@llnl.gov)

## Acknowledgements to previous contributors
P.N. Brown, R.H. Cohen, D.P. Grote, A.C. Hindmarsh, L.L. LoDestro, J.L. Milovich, 
A. Pankin, G.D. Porter, and G.R. Smith, all presently or formerly at LLNL; M. McCourt, 
L.C. McInnes, and H. Zhang, ANL; J.R. Cary, A.H. Hakim, S.E. Kruger, and A. Pankin, Tech-X; 
D.A. Knoll, INEEL; D.P. Stotler, PPPL; B.J. Braams, retired, IAEA; A.Yu. Pigarov and 
R. Smirnov, UCSD; J.D. Elder, U. Toronto; M. Groth, Aalto Univ.; and R.B. Campbell, Sandia.

## References
**_UEDGE development_**   
T.D. Rognlien, J.L. Milovich, M.E. Rensink, and G.D. Porter, J. Nucl. Mat. 196-198 (1992) 347-351.  
G.R. Smith, P.N. Brown, R.B. Campbell, D.A. Knoll, P.R. McHugh, M.E. Rensink, and T.D. Rognlien, J. Nucl. Mater. 220-222 (1995) 1024.  
M.E. Rensink and T.D. Rognlien, J. Nucl. Mater. 266-269 (1999) 1180.  
T.D. Rognlien, D.D. Ryutov, N. Mattor, and G.D. Porter, Phys. Plasmas 6, (1999) 1851.  
T.D. Rognlien, M.E. Rensink, and G.R. Smith, "User manual for the UEDGE edge-plasma transport code," January 2000, LLNL Rpt. UCRL-ID-137121, lastest revision May 1, 2013.  

**_Forthon development_** 
D. P. Grote, A. Friedman, I. Haber, ``Methods used in WARP3d, a Three-Dimensional PIC/Accelerator Code'', Proceedings of the 1996 Computational Accelerator Physics Conference, AIP Conference Proceedings 391, p. 51.  
See also: http://hifweb.lbl.gov/Forthon/  .

**_FACETS project_**    
J.R. Cary, J. Candy, R.H. Cohen et al., J. Phys.: Conf. Ser. 125 (2008) 012040.  
A.H. Hakim, T.D. Rognlien, R.J. Groebner et al., Phys. Plasmas 19 (2012) 032505.  
M. McCourt, T.D. Rognlien, L.C. McInnes, and H. Zhang, Computational Science & Discovery 5 (2012) 014012.  

## Release 

UEDGE is released under an LGPL license.  For more details see the
NOTICE and LICENSE files.

``LLNL-CODE-845914``
------
--------
