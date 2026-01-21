# Citing UEDGE

If you use **UEDGE** in a scientific publication, please cite it appropriately so results can be traced to a specific code version.

## Recommended citation text

At minimum, include:

- Code name: **UEDGE**
- Version number
- Distribution/source (PyPI wheel, custom build, repository)
- Git commit hash (when available)

### Example (free-form)

```
UEDGE version X.Y.Z, Lawrence Livermore National Laboratory (LLNL), Git commit <SHA> (YEAR).
```

## Programmatic provenance (recommended)

UEDGE provides a built-in mechanism to retrieve complete provenance information in preferred citation form
for the version you are using for the following kwargs:

- 'version' - Lists build info available
- 'APA' - APA citation
- 'bibtex' - Bibtex citation
- 'CSL' - CSL-JSON citation

To invoke the automatic citation modules:

```python
import uedge
uedge.cite('version'|'apa'|'bibtex'|'csl')
```

This reports (when available):
- Version number
- Build type (wheel / editable / source)
- Git branch and commit
- Whether local modifications were present

You can copy the output verbatim into manuscripts, reports, or supplementary material. 

> **Note"** This is an automated best-effort to produce a citation for UEDGE. The author is responsible for amending and
> correcting the automated citation as needded.

## BibTeX

> **Note:** Please update the `year`, `url`, and (if you have it) `doi` fields to match your release / archival record.
> If you archive releases on Zenodo, replace `url` and add `doi` accordingly.

```bibtex
@software{uedge,
  title        = {UEDGE},
  author       = {{UEDGE Development Team}},
  organization = {Lawrence Livermore National Laboratory},
  year         = {YEAR},
  version      = {X.Y.Z},
  url          = {https://github.com/LLNL/UEDGE},
  note         = {Git commit: <SHA>; build type: wheel/editable/source}
}
```

## CSL-JSON (for Zotero, Pandoc, etc.)

> **Note:** Please update the `issued` year, `URL`, `version`, and add a `DOI` if available.

```json
{
  "id": "uedge",
  "type": "software",
  "title": "UEDGE",
  "author": [
    { "literal": "UEDGE Development Team" }
  ],
  "publisher": "Lawrence Livermore National Laboratory",
  "version": "X.Y.Z",
  "URL": "https://github.com/LLNL/UEDGE",
  "issued": { "date-parts": [[YEAR]] },
  "note": "Git commit: <SHA>; build type: wheel/editable/source"
}
```

## Prebuilt wheels vs custom builds

- **Prebuilt wheels** usually report a version number but may not include Git metadata.
- **Editable or source builds** from the repository can include branch and commit information and are preferred for reproducibility.

## Acknowledgements

See [ACKNOWLEDGEMENTS.md](ACKNOWLEDGEMENTS.md).
