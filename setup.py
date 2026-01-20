#!/usr/bin/env python
# To use:
#       python setup.py install
#
import os
import glob
import os.path
from os import remove, environ
import string
import site
from Forthon.compilers import FCompiler
import getopt
import logging
from subprocess import call
from sys import hexversion, argv, platform
import sys, subprocess
import numpy
try:
    from setuptools import Extension, setup, Distribution, find_packages
    from setuptools.command.build import build
    from setuptools.command.bdist_wheel import bdist_wheel
except:
    from distutils.core import Extension, setup
    from distutils.dist import Distribution
    from distutils.command.build import build
    from distutils.command.bdist_wheel import bdist_wheel

try:
    # setuptools vendored distutils (what you have in the traceback)
    from setuptools._distutils.ccompiler import CCompiler
except Exception:
    # fallback (older environments)
    from distutils.ccompiler import CCompiler
# Check Python version
if hexversion < 0x03000000:
    raise SystemExit("Python versions < 3 not supported")

# *.o-files to be included in the build (in addition to pkg.o and pkg_p.o)
uedgeobjects = []


os.environ['PATH'] += os.pathsep + site.USER_BASE + '/bin'


arglist = {
        'debug': False,
        'serial': False,
        'checkbounds': False,
        'noclean': False,
#        'petsc': False,
}


optlist, args = getopt.getopt(argv[1:], 'gt:F:', list(arglist.keys()))
machine = platform
fcomp = None
parallel = 0

for o in optlist:
    if o[0] == '-g':
        arglist['debug'] = True
    elif o[0] == '-t':
        machine = o[1]
    elif o[0] == '-F':
        fcomp = o[1]
    elif o[0] == '--parallel':
        parallel = 1
    elif o[0].replace('--','') in arglist.keys():
        arglist[o[0].replace('--','')] = True

# Capture and extract PIP options passed using -C
for arg, _ in arglist.items():
    if "--{}".format(arg) in argv:
        argv.remove("--{}".format(arg))
        arglist[arg] = True

        
# Flags for makefile. Flags are easier to handle from setup.py and it prevents dealing with the makefile.)

FLAGS  = {
    # Specify FORTRAN compiler arguments
    'FARGS': [ 
        '-g -fmax-errors=15', 
        '-DFORTHON',
        '-cpp',
        '-Wconversion',
        '-fimplicit-none',
        '-fopenmp'
    ],
    # Specify C compiler arguments
    'CARGS': [],
    # Specify debugging flags
    'FARGSDEBUG': [
        '-fbacktrace',
        '-ffree-line-length-0', 
        '-fcheck=all',
        '-fbounds-check',
        '-ffpe-trap=invalid,overflow,underflow -finit-real=snan',
        '-Og'
    ],
    # Specify optional flags used if not debugging
    'FARGSOPT': ['-O3'],
    # Specify OMP args
    'OMPargs': ['--omp'],
}

fcompiler = FCompiler(machine=machine,
                      debug=int(arglist['debug']),
                      fcompname=fcomp)


dummydist = Distribution()
dummydist.parse_command_line()
dummybuild = dummydist.get_command_obj('build')
dummybuild.finalize_options()
builddir = dummybuild.build_temp

uedgeobjects = map(lambda p: os.path.join(builddir, p), uedgeobjects)

if os.getenv('PACT_DIR') != None:
    library_dirs = fcompiler.libdirs + [
        os.path.join(os.getenv('PACT_DIR'), 'lib')]
    libraries = ['pdb', 'pml', 'score', 'blas', 'm'] + fcompiler.libs
else:
    library_dirs = fcompiler.libdirs
    libraries = fcompiler.libs

 
if parallel:
    library_dirs = fcompiler.libdirs + ['/usr/lpp/ppe.poe/lib']
    libraries = fcompiler.libs + ['mpi']
    # uedgeobjects = uedgeobjects + ['/usr/local/mpi/ifc_farg.o']

define_macros=[("WITH_NUMERIC", "0"),
               ("FORTHON_PKGNAME", '\"uedgeC\"'),
               ("FORTHON","1")]

# check for readline
rlncom = "echo \"int main(){}\" | gcc -x c -lreadline - "
rln = os.system(rlncom)
if rln == 0: 
   define_macros = define_macros + [("HAS_READLINE","1")]
   os.environ["READLINE"] = "-l readline"
   libraries = ['readline'] + libraries




""" ==========================================================================
                    CONVERSION OF MPPL FILES TO FORTRAN
========================================================================== """
""" Begin old filelist """
def filesublist(path,suffix):
# Methods to create a list of files in a directory, and a sub-list containing sspecified suffixes.
    from os import listdir
    sublist = []
    biglist = listdir(path)
    for filename in biglist:
        filesplit = filename.split(".")
        if filesplit[-1] == suffix:
            sublist.append(filename)
    return sublist
""" End old filelist """

""" Begin old localrules """

# rules for converting mppl to f90
# This is generic down to the UEDGE section
class Use2use:
    def parse(self, s):
        #Return s if this is a comment
        if (not s[0].isspace()) and s[0] != "U":
            return s
        #Do a substitution if line contains "Use" and it is not
        #part of a comment
        if (s.find("Use")+1):
            if (s.find("!")==-1) or (s.find("Use")<s.find("!")):
                sout=s.replace("("," ")
                sout=sout.replace(")"," ")
                sout=sout.replace("Use","      use")
                return sout
        return s

class Allot:
    def parse(self, s):
        # converts MPPL allot to F90 allocate
        if (s.find("allot")+1):
            if (s.find("call allot")+1):
                s=s.replace("call allot","allocate")
                s=s.replace('"','')
                s=s.replace("'","")
                s=s.replace(",","(")
                s=s.replace(")","))")
                return s
            s=s.replace(", allot","")
            s=s.replace(",allot","")
        return s

class Nopdb:
    def parse(self, s):
        from os import getenv
        if getenv('PACT_DIR') == None:
            if s.startswith("c!nopdb"):
                s=s.replace("c!nopdb","       ")
        return s

class Petsc:
    def parse(self, s):
        from os import getenv
        if getenv('PETSC_DIR') != None:
            if getenv('PARALLEL') == None:
               if s.startswith("cunipetsc"):
                   s=s.replace("cunipetsc","")
            if s.startswith("cpetsc"):
                s=s.replace("cpetsc","")
        return s

class Omp:
    def parse(self, s):
        from os import getenv
#        if getenv("OMP", 'False').lower() in ('true', '1', 't'):
        if True:
            if s.startswith("c!omp"):
                s=s.replace("c!omp","     ")
        return s

class MoveDecs:
    def __init__(self):
        self.saved_dec = 0
        self.in_uses = 0
        self.savedlines = []

    def parse(self, s):
        from copy import copy
        # Return if this is a comment
        if (not s[0].isspace()) and (not self.saved_dec) and (not self.in_uses):
            return s
        # collect lines for declarations
        # if we find an "implicit none" statement, store it and remove the line
        sls=s.lstrip().lower()
        indfunc=sls.find("function")
        indcom=sls.find("!")
        functest= (indfunc == -1 or -1<indcom<indfunc)
        # tests to exclude "real function" but allow "real ! function as
        # part of declaration block
        if (sls[0:8]=="implicit") or (sls[0:4]=="real" and functest) \
           or (sls[0:7]=="integer" and functest) \
           or (sls[0:9]=="character") or (sls[0:9]=="parameter") or \
           (sls[0:8]=="external") or (sls[0:9] == "intrinsic") or \
           (sls[0:7]=="logical" and functest) or (sls[0:9]=="dimension") or \
           (sls[0:4] == "data"):
            self.savedlines.append(s)
            self.saved_dec=1
            self.in_uses=0
            return None
        # if we are in the midst of declarations, save also comments (except for
        # "Common block") and continuations and blank lines as part of
        #  what is moved.)
        if (self.saved_dec==1) and (sls == "" or s[0].lower() == "c" or sls[0]=="!" \
          or s[0]=="*") and (self.in_uses == 0) \
          and (s.find("Common block")==-1):
            self.savedlines.append(s)
            return None
        # Check for continuation line in midst of declarations
        if (self.saved_dec==1) and (len(s)>6):
            if (s[5].isspace() == 0):
                self.savedlines.append(s)
                return None

        if (sls[0:3] == "use"):
            self.in_uses=1

        if (self.saved_dec==1) and (sls != "") and s[0] != "c" and sls[0] != "!" and \
           (sls[0:3] != "use"):
            #This is our first executable statement.  Add it to our saved
            # declarations lines and return them now
            templines = copy(self.savedlines)
            templines.append(s)
            #empty out savedlines
            del self.savedlines[0:]
            self.saved_dec = 0
            self.in_uses=0
            return templines
        return s
        
class Elseifthen:
    def __init__(self):
        self.inelseif = 0
        self.savedelselines=""

    def parse(self, s):
        # put a "then" at the end of an elseif if it isn't already there
        # need to check to see if next line is a continue
        # return s if this is a comment
        if (not s[0].isspace()):
            return s
        if s.find("elseif")+1:
            if s.find("then")+1:
                return s
            # set a flag that we are in an "in-else-if" block that needs work
            self.inelseif = 1
            # If there is no "then" we need to save it to test if the next
            # line is a continuation
            self.savedelselines = s
            return(None)
        if (self.inelseif and len(s)>6 and not s[5].isspace()):
            # This is a continue line, add it to savedelselines
            self.savedelselines += s 
            return(None)
        if (self.inelseif and (len(s)<6 or s[5].isspace())):
            # No longer in a continue, so process lines
            if self.savedelselines.find("then")+1:
                self.savedelselines +=  s
                self.inelseif = 0
                return self.savedelselines
            if self.savedelselines.split("\n")[-2].find("!")+1:
                # if last line in saved lines has a comment,
                #  find index of last comment sign
                last = self.savedelselines.rfind("!")
                self.savedelselines=self.savedelselines[0:last]+ \
                    " then "+self.savedelselines[last:] + s
                self.inelseif=0
                return self.savedelselines
            #Otherwise the last line has no comment so insert "then" at end
            self.savedelselines = self.savedelselines[0:-1]+" then\n" + s
            self.inelseif=0
            return self.savedelselines
        return s


#-------------------------------------
# Special for UEDGE
class grdproc:
    def __init__(self):
        self.wordsizectr=0
    def parse(self, s):
        # process to eliminate ifelse write construction
        if (s.find("ifelse([WORDSIZE]")+1):
            s="#if WORDSIZE == 64\n 2001 format(1p3e23.15)\n#else\n 2001 format(1p3d23.15)\n#endif\n"
    #        s="#ifndef WORDSIZE\n 2001 format(1p3e23.15)\n#else\n 2001 format(1p3d23.15)\n#endif\n"
            self.wordsizectr=4
            return s
        elif self.wordsizectr > 0:
            self.wordsizectr -= 1
            return None
        else:
            self.wordsizectr = 0
            return s



""" End old localrules """

""" Begin old convert """
# A Python translator that will convert files following a set
# of file-suffix-specific substitution rules
# Written by Ron Cohen, February 2007.  Last revision March 15, 2007.

# import utilities to create filtered list of files 

"""
create a list of classes to enstantiate which will have conversion rules
for different language files.

Usage:

  Create a file "localrules.py" that defines globalsubrules to apply to
  all file types and any of cppsubrules, pythonsubrules,fortransubrules,
  mpplsubrules,F90subrules that you want to apply to specific file
  types.  These rules are of the form:
  globalsubrules = [("in1","out1"),("in2","out2"),
    ...] where string in1 will be converted to string out1, etc.  It is
    an ordered list which will be executed first to last.
    Another option is that an entry in globalsubrules can be the name
    of a function of one argument (the string being processed).
  Note currently, if the ordered list returns "None" (the python None,
    not the string "None")for a line,
    no further proceessing of that line is done and no line is written.
    This can be used to re-order lines, since a function entry used in
    globalsubrules can be used to return a list of lines, or "None". But
    no further processing is done of the returned list of lines, so any
    such re-ordering must be done as the last of the rules..
  Note in additon to duples, one can also include in the subrules methods
    which take as its sole argument the string the method operates on;
    use this for more complicated operations than simple replacements
  Also create a local script in the directory where you want to
  convert files. Begin that script with "from convert import *",
  making sure that convert.py and localrules.py are in your path.
  Warning: don't end the local script with the suffix .py or the convertor
  will appy your rules to `your local script (assuming you use the default
  suffix as the suffxin for class Py). 

  Further notes

  1.If desired, the default suffixes searched to designate files of
    a specific type can be edited in your local script.  to do so
    import convert, edit the suffix for the sub-class (for example,
    convert.py.suffix = "newsuffix"), then type "from convert import *".

  2.Create sub-classes that inherit from class generic, for any additional
    file suffix that will be converted.  In each subclass, set suffixin
    to a string that gives the suffix of files to be converted, e.g. ".cc"
    or ".py".  If the processed files are to have a different suffix,
    set suffixout to the desired suffix.

  3.Append processall.classlist, a list of sub-class names that will
    be processed to do the conversion, to add the names of any classes
    added in step 2.

  Then, to proces files:
  1.Execute processall(indir,outdir) to process all the specified file
    types.  indir and outdir are optional arguments to specify paths
    for the directory to be processed and the directory where processed
    files will be written.  By default indir = ".", the current directory,
    and outdir = "./converted".
  2.You can also create instances of specific classes created in step 2, and
    just process those files.  For example a = Py(indir,outdir).  See
    documentation for class generic to see available methods.
    WARNING: If you've run the script before and already created the outdir,
    the script will overwrite files in outdir with the same name
  3.The script will only process files if the source file is newer than
    the target (or the target doesn't exist).
  4.You can also use processdirs(dirlist) to run processall in a list of
    directory names specified in dirlist (names in quotes).
"""

# substitution rules to be appended to globalsubrules for each language
# subrules of form [("in1","out1"),("in2","out2"), ...]
# e.g. cppsubrules = [ ("hier::Box<NDIM>",       "Box"),
#    ("hier::IntVector<NDIM>", "IntVect"),
#    ("hier::Index<NDIM>",     "IntVect") ]
subrules = {
    'global': [],
    'cpp': [],
    'python': [],
    'fortran': [],
    'F90': [],
    'F90_90': [],
    'MPPL': [],
}

def fnconvert(name,suffixout):
    # Converts a file name "name" by substituting suffixout for the
    # existing suffix, or if there is no suffix, appending suffixout
    suffixin = name.split(".")[-1]
    if (suffixin != name):
        nameout = name.rstrip(suffixin)+suffixout
        # this coding makes sure only last occurence of suffix is repalced
    else:
        nameout = name+"."+suffixout
    return nameout

def newer(file0,file1):
    from stat import ST_MTIME
    from os import stat
    from numpy import greater
    # returns 1 if file0 is newer than file1, 0 otherwise
    # if file0 does not exist, raise a standard python error
    # if file1 does not exist, return 2
    time0 = stat(file0)[ST_MTIME]
    try:
        time1 = stat(file1)[ST_MTIME]
    except:
        return 2
    return greater(time0,time1)

class processdirs:
    # Process a list of directories to do conversions with processall
    # The default list of directories is ".", the current directory
    # Establishes a set of subdirectories with processed files
    # in a root directory whose default is called "converted" and
    # is parallel to ".".
    def __init__(self,indirs = ["."],outroot = "../converted",clean=None):
        from os import getenv, mkdir, chdir
        # create the output root directory if it needs to be
        curpath = getenv("PWD")
        try:
            mkdir(outroot)
            print("Creating directory "+outroot)
        except:
            try:
                chdir(outroot)
                print("Output directory already exists; proceeding")
                chdir(curpath)
            except:
                raise "Can't create output directory"
        for direc in indirs:
            print("Entering directory "+direc)
            if direc == ".":
                curdir = curpath.split("/")[-1]
                outdir = outroot+"/"+curdir
            else:
                outdir = outroot + "/"+direc
            processall(direc,outdir,clean)

class processall:
    # Process all files in a directory with rules according to file type
    # as designated by suffix.
    # List of classes of distict file types
    classlist = ["Py","Cpp","Cpp_h","Cpp_hh","Fortran","F90","F90_90","MPPL"]
    def __init__(self,indir=".",outdir="./converted",clean=None):
        for entry in self.classlist:
            a=eval(entry+"("+repr(indir)+","+repr(outdir)+",clean="+repr(clean)+")")
            # Only process file types with non-empty substition rules
            if (a.subrules != []):
                print("processing for file type "+entry)
                a.process()

class generic:
    """
    generic class for replacing strings in a series of files.
    Methods:
       process(): processes the list of files in the directory indir with
         suffix suffixin, and for each one creates a new file in directory
         outdir with the same root name appended by
         the string suffixout (which by default is suffixin).
       processfile(filename): processes a specific file
    Notes:
       indir is string with path.  By default it is ".", the current directory
       outdir is by default "./converted".   If the output directory does not
          exist it will be created.  Can be overwritten after instance created
    """

    suffixin = ""
    suffixout = ""
    subrules = subrules['global']
    def __init__(self,indir=".",outdir="./converted",clean=None):
        from os import mkdir, remove
        self.indir = indir
        self.outdir = outdir
        self.clean=clean
        self.doclean = 0   # default is no removal of duplicate files
        if (self.suffixout == ""):
            self.suffixout = self.suffixin
        # create the output directory if it hasn't been created
        try:
            mkdir(self.outdir)
        except:
            pass
    def process(self):
        # get the list of files to process
        self.filelist = filesublist(self.indir,self.suffixin)
        print( "processing directory '" + self.indir+ "' to directory '" + self.outdir +"'")
        # exclude convert.py, localrules.py and filelists.py from filelist
        if (self.suffixin == "py"):
            try:
                self.filelist.remove("convert.py")
            except:
                pass
            try:
                self.filelist.remove("localrules.py")
            except:
                pass
            try:
                self.filelist.remove("filelists.py")
            except:
                pass
        # Set a flag to remove unchanged files if input and output directories are different
        if self.clean and (self.outdir != self.indir):
            # if the input and output directories are different, eliminate files
            # from the output directory that duplicate those in the input directory.
            self.doclean = 1
        else:
            self.doclean = 0
        # Now process files remaining
        for file in self.filelist:
            self.processfile(file)
    def processfile(self,filename):
        from filecmp import cmp
        # Convert an individual file
        self.infile=self.indir+"/"+filename
        if (self.suffixout == self.suffixin):
            outfilename=filename
        else:
            outfilename=fnconvert(filename,self.suffixout)
        self.outfile = self.outdir+"/"+outfilename
        if (newer(self.infile,self.outfile)):
            # Only process infile if it is newer than outfile
            print("converting file "+self.infile+" to "+self.outfile)
            f=open(self.infile,"r")
            g=open(self.outfile,"w")
            lines=f.readlines()
            iline = 0   # just a diagnostic, counts lines
            for line in lines:
              iline = iline + 1
              for rule in self.subrules:
                if (line != None):
                  if type(rule)==type((0,0)) or type(rule)==type([0,0]):
                    line = line.replace(rule[0],rule[1])
                  else:
                    # if it's not a tuple or list, assume it is a method
                    # that executes something
                    # more complicated than a simple replace
                    line = rule.parse(line)
              if (line != None):
                  if (type(line) == type("fubar")):
                      g.write(line)
                  if (type(line) == type([0,0])):
                      g.writelines(line)
                    
            g.close()
            f.close()
            if self.doclean == 1:
              if cmp(self.outfile,self.infile) and self.outfile != self.infile:
                print(self.outfile + " is the same as " + self.infile + ", removing")
                remove(self.outfile)

class Py(generic):
    suffixin="py"
    subrules = subrules['global'] + subrules['python']

class Cpp(generic):
    suffixin = "C"
    subrules = subrules['global'] + subrules['cpp']

class Cpp_h(Cpp):
    suffixin = "h"

class Cpp_hh(Cpp):
    suffixin = "hh"

class Fortran(generic):
    suffixin = "f"
    subrules = subrules['global'] + subrules['fortran']

class F90(generic):
    suffixin = "F"
    subrules = subrules['global'] + subrules['F90']

class F90_90(F90):
    suffixin = "F90"

class MPPL(generic):
    suffixin = "m"
    subrules = subrules['global'] + subrules['MPPL']

""" End old convert """


""" Begin old convertor """
class M2F(generic):
   suffixin = "m"
   suffixout = "F"
   subrules = subrules['global'] +  [ 
            ("#","!"),
            Use2use(),
            ("c!ifdef","#ifdef"),
            ("do i1=","do_i1: do i1="),
            ("break (2) ! exit do_i1", "exit do_i1"),
            ("enddo ! do_i1", "enddo do_i1"),
            ("dfloat","real"),
            ("float","real"),
            ("c!else","#else"),
            ("c!endif","#endif"),
            ("(Size4)","(kind=4)::"),
            grdproc(),
            (":: function"," function"),
            (" break "," exit "),
            (" break\n"," exit\n"),
            ("while (","do while ("),
            ("endwhile","end do"),
            ("      call ruthere","c     call ruthere"),
            ("c!include ","#include "),
            Nopdb(),
            Petsc(),
            Omp(),
            Allot(),
            Elseifthen(),
            MoveDecs(),
        ]

def mppl2f90(
    dirlist = ["com", "bbb", "aph", "api", "flx", "grd", "svr", "wdf", "ncl"],
    targetdir = "."
    ):
    """ Script to convert MPPL files to .F for uedge. """
    #create the directory list
    print("mppl2f90 directory list ="+repr(dirlist))
    # define the mppl to f90 class


    processall.classlist.append("M2F")
    pd = processdirs(dirlist,targetdir)

""" End old convertor """



""" ==========================================================================
                    END CONVERSION OF MPPL FILES TO FORTRAN
========================================================================== """

def _libgomp_postargs():
    fc = os.environ.get("FC", "gfortran")
    if sys.platform == "darwin":
        lib = subprocess.check_output([fc, "-print-file-name=libgomp.1.dylib"], text=True).strip()
    else:
        lib = subprocess.check_output([fc, "-print-file-name=libgomp.so.1"], text=True).strip()

    if os.sep not in lib:
        return ["-lgomp"]

    libdir = os.path.dirname(lib)
    return [f"-L{libdir}", "-lgomp", f"-Wl,-rpath,{libdir}"]


def enable_libgomp_for_fortran_link():
    try:
        from numpy.distutils.fcompiler import FCompiler
    except Exception:
        return  # numpy not present / not using numpy.distutils

    orig = FCompiler.link_shared_object
    extra = _libgomp_postargs()

    def patched(self, objects, output_libname, output_dir=None,
                libraries=None, library_dirs=None, runtime_library_dirs=None,
                export_symbols=None, debug=0, extra_preargs=None,
                extra_postargs=None, build_temp=None, target_lang=None):
        if extra_postargs is None:
            extra_postargs = []
        # Append if missing
        for e in extra:
            if e not in extra_postargs:
                extra_postargs.append(e)

        return orig(self, objects, output_libname, output_dir=output_dir,
                    libraries=libraries, library_dirs=library_dirs,
                    runtime_library_dirs=runtime_library_dirs,
                    export_symbols=export_symbols, debug=debug,
                    extra_preargs=extra_preargs, extra_postargs=extra_postargs,
                    build_temp=build_temp, target_lang=target_lang)

    FCompiler.link_shared_object = patched



class uedgeWheel(bdist_wheel):
    from os import environ
    if arglist['noclean']:
        environ['NOCLEAN'] = "TRUE"
    else:
        environ['NOCLEAN'] = "FALSE"



# Modify Build to inject call to building UEDGE source
class uedgeBuild(build):
    from os import environ
    if arglist['noclean']:
        environ['NOCLEAN'] = "TRUE"
    else:
        environ['NOCLEAN'] = "FALSE"

    def run(self):
        # Convert MPPL files to F90 files
        '''
        # Perform PETSc build [deprecated]
        if arglist['petsc']:
            if arglist['noclean'] is False:
                if call(['make', '-f', 'Makefile.PETSc', 'clean']) != 0:
                        raise SystemExit("Pre-PETSc UEDGE build clean failure")
            if os.getenv('PETSC_DIR') == None:
                raise SystemExit("PETSc requested but PETSC_DIR not set")
            if os.getenv('PETSC_ARCH') == None:
                raise SystemExit("PETSc requested but PETSC_ARCH not set")
            PETSC_DIR = os.getenv('PETSC_DIR')
            PETSC_ARCH = os.getenv('PETSC_ARCH')
            library_dirs = fcompiler.libdirs + [
                    os.path.join(PETSC_DIR, PETSC_ARCH, 'lib')
            ]
            libraries = [
                    'petscts', 'petscsnes', 'petscksp', 'petscdm', 'petscmat',
                    'petscvec', 'petsc', 'HYPRE', 'mpich', 'lapack', 'blas', 
                    'X11', 'pthread', 'rt', 'stdc++', 'm'] + fcompiler.libs
            libraries = ['petsc'] + fcompiler.libs
            if call(['make', '-f', 'Makefile.PETSc']) !=0:
                raise SystemExit("PETSc build failure")
            build.run(self)
            return
        '''
        # If serial buid requested, remove parallel flags
        if arglist['serial']:
            FLAGS['OMPargs'].remove('--omp')
            FLAGS['FARGS'].remove('-fopenmp')
            # Flag for localrules
            os.environ["OMP"] = "FALSE"
        else:
            # Flag for localrules
            os.environ["OMP"] = "TRUE"
        # Add boundary checking to non-debug build
        if arglist['checkbounds']:
            FLAGS['FARGSOPT'] = FLAGS['FARGSOPT']  + [
                '-fbounds-check',
                '-fcheck=all',
            ]
        # Check whether to add debug flags
        if arglist['debug']:
            FLAGS['FARGS'] = FLAGS['FARGS'] + FLAGS['FARGSDEBUG']
        else:
            FLAGS['FARGS'] = FLAGS['FARGS'] + FLAGS['FARGSOPT']
        # Clean build directory unless flagged not to
        if arglist['noclean'] is False:
            print("Cleaning UEDGE build directory")
            status = call(['make', '-f', 'Makefile.Forthon', 'clean'])
            if status != 0: 
                raise SystemExit("Pre-UEDGE build clean failure")
            print("Build directory cleaned")

    
        if os.environ.get("UEDGE_CLEAN_BUILD", "") == "1":
            print("UEDGE_CLEAN_BUILD=1: cleaning build directory (env)")
            shutil.rmtree("build", ignore_errors=True)


        print("Converting MPPL files to Fortran in build")
        mppl2f90()
        print("MPPL files converted to Fortran")

        # Compile arguments into compiler flags    
        OMPFLAGS='OMPFLAGS = {}'.format(' '.join(FLAGS['OMPargs']))
        COMPILEFLAGS ='DEBUG = -v --fargs "{}"'.format(' '.join(FLAGS['FARGS']))
        if FLAGS['CARGS']!=[]:
            COMPILEFLAGS = COMPILEFLAGS+' --cargs="{}"'.format(' '.join(FLAGS['CARGS']))
        # Fix libraries
        if not arglist['serial']:
            enable_libgomp_for_fortran_link()

        py_tag = sys.implementation.cache_tag  # cpython-39, pypy39, etc.
        #NOTE: removing flags causes crash: unclear why
#        # Call compiler
        status = call(['make',COMPILEFLAGS,OMPFLAGS,'-f','Makefile.Forthon', f"PYTAG={py_tag}", f"MYPYTHON={sys.executable}"]) 
        if status != 0: 
            raise SystemExit("Build failure")
        # Run the normal setuptools build (creates build_lib, builds ext_modules, etc.)
        build.run(self)


        # Sanity check: ensure generated modules are where we expect for packaging
        pkg_out = os.path.join(self.build_lib, "uedge")
        need = [f"{pkg}*.so" for pkg in ['uedgeC', 'compy', 'grdpy', 'flxpy', 'bbbpy', 'ppppy', 'svrpy', 'wdfpy', 'aphpy', 'apipy', 'nclpy']]
        missing = []
        for pat in need:
            if not glob.glob(os.path.join(pkg_out, pat)):
                missing.append(pat)

        if missing:
            raise SystemExit(
                f"Expected generated modules missing from {pkg_out}: {missing}\n"
                f"Found: {sorted(os.listdir(pkg_out))}"
            )


# Modify clean to use make clean
class uedgeClean(build):
    def run(self):
#        if call(['make', '-f', 'Makefile.'+ "PETSc"*(arglist['petsc']) + "Forthon"*(not arglist['petsc']), 'clean']) != 0:
        if call(['make', '-f', 'Makefile.Forthon', 'cleanbuild']) != 0:
            raise SystemExit("Clean failure")

setup(  name="uedge",
        ext_modules= [
            Extension(
                'uedge.uedgeC',
                [   'uedgeC_Forthon.c',
                    os.path.join(builddir, 'Forthon.c'),
                    'com/handlers.c', 
                    'com/vector.c',
                    'bbb/exmain.c'
                ],
                include_dirs = [
                    builddir, 
                    numpy.get_include()
                ],
                library_dirs=library_dirs,
                libraries=libraries,
                define_macros=define_macros,
                extra_objects=uedgeobjects,
                extra_link_args=FLAGS['CARGS']+['-g','-DFORTHON'] +
                fcompiler.extra_link_args,
                extra_compile_args=fcompiler.extra_compile_args
            )
        ],
        cmdclass={
            'build': uedgeBuild,
            'clean': uedgeClean,
            'bdist_wheel': uedgeWheel,
        },
)



