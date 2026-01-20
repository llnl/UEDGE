# import all of the neccesary packages
from os import path, RTLD_LAZY, RTLD_GLOBAL
from uedge import __path__
from resource import setrlimit
from sys import setdlopenflags
try:
    from json import loads
    from  urllib import request
except:
    pass
try:
    from  importlib import metadata
except:
    import pkg_resources

# --- With Python3, the so files of each Fortran package are imported
# --- separately. The dlopen flag needs to be set so that cross references
# --- among the packages can be satisfied.
setdlopenflags(RTLD_LAZY | RTLD_GLOBAL)
from . import uedgeC
from .uedgeC import *
from .compy import com
from .grdpy import grd
from .flxpy import flx
from .bbbpy import bbb
from .svrpy import svr
from .wdfpy import wdf
from .apipy import api
from .aphpy import aph
from .ppppy import ppp
from .nclpy import ncl


# --- The UEDGE modules must be imported in the order below because of
# --- linking dependencies.


# --- Check if the compiler was ifort - if so, set the stacksize unlimited
# --- The fcompname is not yet be available yet if Forthon is not up to date
try:
    if fcompname == 'ifort':
        setrlimit(resource.RLIMIT_STACK, (-1, -1))
except:
    pass

pkg = "uedge"
try:
    thisver = metadata.version(pkg)
except:
    thisver = pkg_resources.get_distribution(pkg).version
    
# Hard-write version to equal that set in VERSION
with open(path.join(__path__[0],"VERSION")) as f:
    __version__ = f.read().replace('\n', '').strip()

# Load the startup file .uedgerc.py from cwd or home.
_homefile = path.join(path.expanduser('~'), '.uedgerc.py')
_localfile = path.join(path.expanduser('.'), '.uedgerc.py')

if path.exists(_localfile):
   with open(_localfile) as f:
      exec(open(_localfile).read())
elif path.exists(_homefile):
   with open(_homefile) as f:
      exec(open(_homefile).read())
      
# Ensure that the UEDGE internal version matches the Python release
bbb.uedge_ver = __version__

# Set up crmpath to point to the included 
aph.crmdir = f"{__path__[0]}/rates/aph"
#
# Check version: included from checkver
# 
def check_uedge_ver():
    try:
        contents = request.urlopen("https://pypi.org/pypi/" + pkg + "/json").read()
        data = loads(contents.decode())
        thatver = data["info"]["version"]

        if thisver < thatver:
            print()
            print("Uedge version " + thisver + ", an update is available to " + thatver)
            print()

    except Exception as err:
        print()
        print("Error checking pypi version: {}".format(err))
        print()


def set_prompt(prompt='UEDGE>>> '):
    import sys
    try:
       import IPython
       from IPython.terminal.prompts import Prompts,Token
       from IPython.terminal.embed import InteractiveShellEmbed
    except:
       pass
    try:

       class MyPrompt(Prompts):
         def in_prompt_tokens(self, cli=None):
             return [(Token.Prompt, prompt)]
         def out_prompt_tokens(self, cli=None):
             return [(Token.Prompt, prompt)]

       get_ipython
    except:
       sys.ps1='UEDGE>>> '
    else:
       ip = get_ipython()
       ip.prompts = MyPrompt(ip)

#set_prompt()
check_uedge_ver()

