# src/uedge/_provenance.py
from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


def _read_build_info(pkg_root: Path) -> dict | None:
    """
    Read metadata baked into the wheel/sdist at build time.
    Expected path: src/uedge/_build_info.json (installed into package).
    """
    p = pkg_root / "_build_info.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None

@dataclass(frozen=True)
class UEDGEProvenance:
    name: str
    version: str
    build_type: str                 # "wheel", "editable", "source"
    git_branch: str | None
    git_commit: str | None
    git_commit_date_iso: str | None  # e.g. "2026-01-14T03:21:10-08:00"
    git_dirty: bool | None
    repo_url: str | None            # e.g. "https://github.com/LLNL/UEDGE"

    # ---------- human-readable ----------
    def as_text(self) -> str:
        parts: list[str] = [f"{self.name} version {self.version}"]

        if self.build_type:
            parts.append(f"Build type: {self.build_type}")
        if self.git_commit:
            parts.append(f"Git commit: {self.git_commit}")
        if self.git_branch and self.git_branch != "HEAD":
            parts.append(f"Git branch: {self.git_branch}")
        if self.git_dirty:
            parts.append("Working tree contained local modifications")

        return "; ".join(parts)

    # ---------- BibTeX ----------
    def bibtex(self, *, key: str = "uedge") -> str:
        """
        Returns a BibTeX @software entry suitable for copy/paste.
        """
        year = self._year_or_none()
        url = self.repo_url or "https://github.com/LLNL/UEDGE"

        note_bits: list[str] = []
        if self.git_commit:
            note_bits.append(f"Git commit: {self.git_commit}")
        if self.git_branch and self.git_branch != "HEAD":
            note_bits.append(f"branch: {self.git_branch}")
        if self.build_type:
            note_bits.append(f"build: {self.build_type}")
        if self.git_dirty:
            note_bits.append("local modifications present")

        note = "; ".join(note_bits) if note_bits else None

        def f(field: str, value: str | None) -> str:
            if not value:
                return ""
            return f"  {field:<12}= {{{_escape_bibtex(value)}}},\n"

        out = "@software{" + key + ",\n"
        out += f("title", self.name)
        out += f("author", "UEDGE Development Team")
        out += f("organization", "Lawrence Livermore National Laboratory")
        out += f("year", year)
        out += f("version", self.version if self.version != "unknown" else None)
        out += f("url", url)
        out += f("note", note)
        out += "}\n"
        return out

    def write_bibtex(self, path: str | Path, *, key: str = "uedge") -> Path:
        p = Path(path)
        p.write_text(self.bibtex(key=key), encoding="utf-8")
        return p

    # ---------- CSL-JSON ----------
    def csl_json(self, *, id: str = "uedge") -> str:
        """
        Returns a CSL-JSON record as a JSON string for Zotero/Pandoc etc.
        """
        url = self.repo_url or "https://github.com/LLNL/UEDGE"

        note_bits: list[str] = []
        if self.git_commit:
            note_bits.append(f"Git commit: {self.git_commit}")
        if self.git_branch and self.git_branch != "HEAD":
            note_bits.append(f"branch: {self.git_branch}")
        if self.build_type:
            note_bits.append(f"build: {self.build_type}")
        if self.git_dirty:
            note_bits.append("local modifications present")
        note = "; ".join(note_bits) if note_bits else None

        year_int = self._year_int_or_none()

        record: dict = {
            "id": id,
            "type": "software",
            "title": self.name,
            "author": [{"literal": "UEDGE Development Team"}],
            "publisher": "Lawrence Livermore National Laboratory",
            "version": self.version if self.version != "unknown" else None,
            "URL": url,
            "note": note,
        }
        # Remove Nones
        record = {k: v for k, v in record.items() if v is not None}

        if year_int is not None:
            record["issued"] = {"date-parts": [[year_int]]}

        return json.dumps(record, indent=2, sort_keys=False) + "\n"

    def write_csl_json(self, path: str | Path, *, id: str = "uedge") -> Path:
        p = Path(path)
        p.write_text(self.csl_json(id=id), encoding="utf-8")
        return p

    # ---------- APA-style ----------
    def apa(self) -> str:
        """
        Returns an APA-ish reference string for software, suitable for copy/paste.

        Example:
          UEDGE Development Team. (2026). UEDGE (Version 7.2.1) [Computer software].
          Lawrence Livermore National Laboratory. https://github.com/LLNL/UEDGE
          (Git commit: abc123; build: editable).
        """
        year = self._year_or_nd()
        url = self.repo_url or "https://github.com/LLNL/UEDGE"

        version_part = ""
        if self.version and self.version != "unknown":
            version_part = f" (Version {self.version})"

        trailing_bits: list[str] = []
        if self.git_commit:
            trailing_bits.append(f"Git commit: {self.git_commit}")
        if self.git_branch and self.git_branch != "HEAD":
            trailing_bits.append(f"branch: {self.git_branch}")
        if self.build_type:
            trailing_bits.append(f"build: {self.build_type}")
        if self.git_dirty:
            trailing_bits.append("local modifications present")
        trailing = f" ({'; '.join(trailing_bits)})." if trailing_bits else "."

        return (
            f"UEDGE Development Team. ({year}). {self.name}{version_part} "
            f"[Computer software]. Lawrence Livermore National Laboratory. {url}{trailing}\n"
        )

    def write_apa(self, path: str | Path) -> Path:
        p = Path(path)
        p.write_text(self.apa(), encoding="utf-8")
        return p

    # ---------- helpers ----------
    def _year_int_or_none(self) -> int | None:
        y = self._year_or_none()
        if not y:
            return None
        try:
            return int(y)
        except Exception:
            return None

    def _year_or_none(self) -> str | None:
        # Prefer git commit date (most defensible); otherwise none.
        if self.git_commit_date_iso:
            try:
                # Accept "YYYY-..." formats; parse first 4 digits.
                return str(int(self.git_commit_date_iso[:4]))
            except Exception:
                return None
        return None

    def _year_or_nd(self) -> str:
        return self._year_or_none() or "n.d."


def _escape_bibtex(s: str) -> str:
    # Minimal BibTeX escaping (good enough for these fields).
    return (
        s.replace("\\", "\\\\")
        .replace("{", "\\{")
        .replace("}", "\\}")
    )


def _run_git(cmd: list[str], cwd: Path) -> str | None:
    try:
        out = subprocess.check_output(
            ["git"] + cmd,
            cwd=cwd,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return out or None
    except Exception:
        return None


def _detect_git_info(pkg_root: Path) -> tuple[str | None, str | None, str | None, bool | None]:
    """
    Returns (branch, commit, commit_date_iso, dirty) if .git exists; else (None, None, None, None).
    """
    if not (pkg_root / ".git").exists():
        return None, None, None, None

    commit = _run_git(["rev-parse", "HEAD"], pkg_root)
    branch = _run_git(["rev-parse", "--abbrev-ref", "HEAD"], pkg_root)

    # ISO 8601 commit date of HEAD (committer date)
    commit_date_iso = _run_git(["show", "-s", "--format=%cI", "HEAD"], pkg_root)

    dirty_txt = _run_git(["status", "--porcelain"], pkg_root)
    dirty = bool(dirty_txt) if dirty_txt is not None else None

    return branch, commit, commit_date_iso, dirty


def _detect_build_type(pkg_root: Path) -> str:
    # Heuristic: if .git exists, likely a checkout / editable use.
    if (pkg_root / ".git").exists():
        return "editable"
    return "wheel"


def _read_version(pkg_root: Path) -> str:
    try:
        return (pkg_root / "VERSION").read_text(encoding="utf-8").strip()
    except Exception:
        return "unknown"

def get_uedge_provenance(*, repo_url: str | None = "https://github.com/LLNL/UEDGE") -> UEDGEProvenance:
    pkg_root = Path(__file__).resolve().parent

    version = _read_version(pkg_root)

    build_info = _read_build_info(pkg_root)
    if build_info:
        branch = build_info.get("git_branch")
        commit = build_info.get("git_sha")
        commit_date_iso = build_info.get("git_commit_date_iso") or None
        dirty = build_info.get("git_dirty")  # optional; usually omitted for wheels
        build_type = build_info.get("source") or "wheel"
        repo_url = build_info.get("repo_url") or repo_url
    else:
        build_type = _detect_build_type(pkg_root)
        branch, commit, commit_date_iso, dirty = _detect_git_info(pkg_root)

    return UEDGEProvenance(
        name="UEDGE",
        version=version,
        build_type=build_type,
        git_branch=branch,
        git_commit=commit,
        git_commit_date_iso=commit_date_iso,
        git_dirty=dirty,
        repo_url=repo_url,
    )

