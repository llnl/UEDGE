#!/usr/bin/env python
# To use:
#       python setup.py install
#
import sys
import os
import os.path
import string
import site
from Forthon.compilers import FCompiler
import getopt
import glob
import shutil
import logging
from subprocess import call
import numpy
try:
    from setuptools import Extension, setup, Distribution, find_packages
    from setuptools.command.build import build
#    from setuptools.command.develop import develop
#    from setuptools.command.egg_info import egg_info
except:
    from distutils.core import Extension, setup
    from distutils.dist import Distribution
    from distutils.command.build import build
os.environ['PATH'] += os.pathsep + site.USER_BASE + '/bin'


optlist, args = getopt.getopt(sys.argv[1:], 'gt:F:', ['parallel', 'omp'])
machine = sys.platform
debug = 0
fcomp = None
parallel = 0
#petsc = 0

for o in optlist:
    if o[0] == '-g':
        debug = 1
    elif o[0] == '-t':
        machine = o[1]
    elif o[0] == '-F':
        fcomp = o[1]
    elif o[0] == '--parallel':
        parallel = 1
#    elif o[0] == '--petsc':
#        petsc = 1
    elif o[0] == '--omp':
        os.putenv("OMP","1")
        


#if petsc == 1 and os.getenv('PETSC_DIR') == None:
#    raise SystemExit("PETSc requested but PETSC_DIR not set")
#if os.getenv('PETSC_DIR') != None:
#    petsc = 1
#if petsc == 1 and os.getenv('PETSC_ARCH') == None:
#    raise SystemExit("PETSc requested but PETSC_ARCH not set")


sys.argv = ['setup2.py']+args
fcompiler = FCompiler(machine=machine,
                      debug=debug,
                      fcompname=fcomp)


#class uedgeInstall(build):
#    def run(self):
#       install.run(self)
#       logging.basicConfig(stream=sys.stderr,level=logging.INFO)
#       log = logging.getLogger()
#       log.info("test")

class uedgeBuild(build):
    def run(self):
        print("Building UEDGE...")
        # with python2 everything is put into a single uedgeC.so file
        if sys.hexversion < 0x03000000:
            raise SystemExit("Python versions < 3 not supported")

        if os.environ.get("UEDGE_CLEAN_BUILD", "") == "1":
            print("UEDGE_CLEAN_BUILD=1: cleaning build directory")
            shutil.rmtree(builddir, ignore_errors=True)

        # ABI-specific build dir so stamps don't collide between cp38/cp39/etc
        py_tag = sys.implementation.cache_tag  # e.g. "cpython-39", "pypy39"
        builddir = f"build-{py_tag}"

        status = call([
            "make", "-f", "Makefile.Forthon",
            f"BUILDDIR={builddir}",
            f"MYPYTHON={sys.executable}",
        ])
        if status != 0: 
            raise SystemExit("Build failure")
        # Run the normal setuptools build (creates build_lib, builds ext_modules, etc.)
        build.run(self)



        # Sanity check: ensure generated modules are where we expect for packaging
        pkg_out = os.path.join(self.build_lib, "uedge")
        need = ["compy*.so", "uedgeC*.so"]
        missing = []
        for pat in need:
            if not glob.glob(os.path.join(pkg_out, pat)):
                missing.append(pat)

        if missing:
            raise SystemExit(
                f"Expected generated modules missing from {pkg_out}: {missing}\n"
                f"Found: {sorted(os.listdir(pkg_out))}"
            )




class uedgeClean(build):
    def run(self):
        print("Cleaning UEDGE")
        if sys.hexversion < 0x03000000:
            raise SystemExit("Python versions < 3 not supported")
        else:
#            if petsc == 0:
            status = call(['make', '-f', 'Makefile.Forthon', 'clean'])
#            else:
#                status = call(['make', '-f', 'Makefile.PETSc', 'clean'])
            if status != 0: 
                raise SystemExit("Clean failure")


#class CustomDevelopCommand(develop):
#    def run(self):
#        print("Develop called")
#        develop.run(self)


#class CustomEggInfoCommand(egg_info):
#    def run(self):
#        print("EggInfo called")
#        egg_info.run(self)

uedgepkgs = ['aph', 'api', 'bbb', 'com', 'flx', 'grd', 'svr', 'wdf', 'ncl']


def makeobjects(pkg):
    return [pkg+'_p.o', pkg+'pymodule.o']


uedgeobjects = []

# add here any extra dot o files other than pkg.o, pkg_p.o


if sys.hexversion < 0x03000000:
    raise SystemExit("Python versions < 3 not supported")
else:
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

#if petsc:
    # PETSC_DIR = '/homes/mccomic/petsc-uedge'
    # PETSC_ARCH = 'linux-uedge'
#    PETSC_DIR = os.getenv('PETSC_DIR')
#    PETSC_ARCH = os.getenv('PETSC_ARCH')
#    library_dirs = fcompiler.libdirs + \
#        [os.path.join(PETSC_DIR, PETSC_ARCH, 'lib')]
#    libraries = ['petscts', 'petscsnes', 'petscksp', 'petscdm', 'petscmat',
#                 'petscvec', 'petsc', 'HYPRE', 'mpich', 'lapack', 'blas', 'X11',
#                 'pthread', 'rt', 'stdc++', 'm'] + fcompiler.libs
#    libraries = ['petsc'] + fcompiler.libs

if parallel:
    library_dirs = fcompiler.libdirs + ['/usr/lpp/ppe.poe/lib']
    libraries = fcompiler.libs + ['mpi']
    # uedgeobjects = uedgeobjects + ['/usr/local/mpi/ifc_farg.o']

#with open('pyscripts/__version__.py','w') as ff:
#    ff.write("__version__ = '%s'\n"%version)
#with open('pyscripts/__src__.py','w') as ff:
#    ff.write("__src__ = '%s'\n"%os.getcwd())

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


setup(
        packages=find_packages(),
        include_package_data=True,
        package_data={
             "uedge": ["*.so", "*.dylib", "*.pyd", "VERSION", "uedgeC.py"],
        },
        ext_modules=[Extension('uedge.uedgeC',
                            ['uedgeC_Forthon.c',
                            os.path.join(builddir, 'Forthon.c'),
                            'com/handlers.c', 'com/vector.c','bbb/exmain.c'],
                            include_dirs=[builddir, numpy.get_include()],
                            library_dirs=library_dirs,
                            libraries=libraries,
                            define_macros=define_macros,
                            extra_objects=uedgeobjects,
                            extra_link_args=['-g','-DFORTHON'] +
                            fcompiler.extra_link_args,
                            extra_compile_args=fcompiler.extra_compile_args
                    )],
        cmdclass={
            'build': uedgeBuild, 
            'clean': uedgeClean,
#            'develop': CustomDevelopCommand,
#            'egg_info': CustomEggInfoCommand,
        },
#      test_suite="pytests",
)
