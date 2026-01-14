UEDGE

NAME
     UEDGE - Unified Tokamak Edge Physics Modeling Code
             A collaborative project involving LLNL, INEL, and 
             Sandia - New Mexico


SYNOPSIS
  uedge [-vers] [<initial command>]

  where
     vers = Version to be run.  Possible values:
               old    (run the version in /usr/local/old)
       	       pub    (run the version in /usr/local/pub) (default)
               new    (run the version in /usr/local/new)

     <initial command> = the first command for UEDGE to execute
                         after it completes initialization; this
                         command often has the form read <filename>,
                         where <filename> contains more commands


DESCRIPTION
     UEDGE solves fluid equations describing transport and sources and
     sinks of charged and neutral species in the edge region of an
     axisymmetric tokamak plasma.  A fully implicit numerical scheme is
     used for either a time-dependent simulation or a Newton interation 
     to steady state. The fluid equations are based on those of 
     Braginski, but have added to them anomalous transport across 
     magnetic field lines.  The region is two-dimensional (radial and 
     poloidal), generally contains a separatrix, and has boundaries 
     of various kinds.

FILES
     UEDGE resides in public in the directory /usr/local/pub
     For easy access modify your search path:
     C-shell users:
	  set path = (oldpath /usr/local/pub)
     Bourne and K-shell users:
	  PATH=oldpath:/usr/local/pub; export PATH

EXAMPLES
     Test cases for UEDGE are available in cfs in tar files that can
     be retrieved by typing the following commands:
          % cfs
          ? get dir=/003030/uedge tc_*.tar
          ? end

LIMITATIONS
     UEDGE must be run from your work directory.  If you try to run it
     in your home directory it will run out of disk space.  Unless you
     have run cpulim to set a longer limit, interactive jobs must
     complete in less than 20 minutes.  Longer jobs need to be
     submitted to NQS.  Read the man pages on qsub for an explanation
     of the memory and time command line options, which must be
     specified before submittal. Set the environmental varible NCPUS=1.

BUGS
     No known bugs, but convergence is difficult in some cases.

SEE ALSO
     A manual for UEDGE users and developers is available in file
     uedge.man in the cfs directory /003030/uedge.

     In case of problems contact a UEDGE developer by e-mail:
          rognlien@wobegon.llnl.gov or gsmith@llnl.gov

Last Modified 3/4/93
