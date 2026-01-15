from __future__ import annotations
from .uedge import *
from os import path
from uedge import __path__
from ._banner import maybe_print_banner, set_banner_config
from ._provenance import get_uedge_provenance
from ._checkver import _check_newer_uedge_ver

with open(path.join(__path__[0],"VERSION")) as f:
    __version__ = f.read().replace('\n', '').strip()

# Load the startup file .uedgerc.py from cwd or home.
#
_homefile = path.join(path.expanduser('~'), '.uedgerc.py')
_localfile = path.join(path.expanduser('.'), '.uedgerc.py')

if path.exists(_localfile):
   with open(_localfile) as f:
      exec(open(_localfile).read())
elif path.exists(_homefile):
   with open(_homefile) as f:
      exec(open(_homefile).read())
      
#
# Ensure that the UEDGE internal version matches the Python release
#

uedge.bbb.uedge_ver = __version__

maybe_print_banner()


def cite(
    fmt: str = "version",
    *,
    repo_url: str | None = "https://github.com/LLNL/UEDGE",
    key: str = "uedge",
    csl_id: str = "uedge",
    print_to: object | None = sys.stderr,
) -> str:
    """
    Return a copy-pasteable reference for this installed UEDGE.

    Parameters
    ----------
    fmt:
        One of: "version", "bibtex", "apa", "csl"
        - "version": provenance string (prov.as_text())
        - "bibtex": BibTeX @software entry
        - "apa": APA-ish software reference
        - "csl": CSL-JSON record
    repo_url:
        Canonical repository URL to embed in formatted citations.
    key:
        BibTeX entry key (only used for fmt="bibtex").
    csl_id:
        CSL-JSON record id (only used for fmt="csl").
    print_to:
        If provided, also write the returned string to this stream (e.g. sys.stderr).

    Returns
    -------
    str
        The requested citation/provenance text.
    """
    fmt_norm = fmt.strip().lower()
    prov = get_uedge_provenance(repo_url=repo_url)

    if fmt_norm in {"version", "text", "provenance"}:
        out = prov.as_text() + "\n"
    elif fmt_norm in {"bibtex", "bib"}:
        out = prov.bibtex(key=key)
    elif fmt_norm in {"apa"}:
        out = prov.apa()
    elif fmt_norm in {"csl", "csl-json", "csljson", "json"}:
        out = prov.csl_json(id=csl_id)
    else:
        raise ValueError(
            f"Unknown fmt={fmt!r}. Expected one of: 'version', 'bibtex', 'apa', 'csl'."
        )

    if print_to is not None:
        try:
            print_to.write(out)
        except Exception:
            pass

    return 


