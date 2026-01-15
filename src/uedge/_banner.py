from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from uedge._checkver import _check_newer_uedge_ver

# Sentinel: ensures banner prints only once per interpreter session
_ALREADY_SHOWN = False


ASCII_LOGO = r"""
    __  ____________  ____________
   / / / / ____/ __ \/ ____/ ____/
  / / / / __/ / / / / / __/ __/
 / /_/ / /___/ /_/ / /_/ / /___
 \____/_____/_____/\____/_____/
""".strip("\n")


def _truthy_env(name: str, default: str = "1") -> bool:
    raw = os.getenv(name, default)
    val = raw.strip().lower()
    if val in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if val in {"0", "false", "f", "no", "n", "off"}:
        return False
    return default.strip().lower() in {"1", "true", "t", "yes", "y", "on"}


@dataclass(frozen=True)
class BannerConfig:
    enabled: bool = True
    show_logo: bool = True
    show_version: bool = True
    show_messages: bool = True
    stream: object = sys.stderr


DEFAULT_CONFIG = BannerConfig(
    enabled=_truthy_env("UEDGE_BANNER", "1"),
    show_logo=_truthy_env("UEDGE_LOGO", "1"),
    show_version=_truthy_env("UEDGE_VERSION", "1"),
    show_messages=_truthy_env("UEDGE_MESSAGES", "1"),
    stream=sys.stderr,
)


def set_banner_config(**kwargs) -> None:
    """Runtime override of banner behavior."""
    global DEFAULT_CONFIG
    DEFAULT_CONFIG = BannerConfig(**{**DEFAULT_CONFIG.__dict__, **kwargs})


def read_version_from_file() -> str:
    """
    Reads src/uedge/VERSION (matches setuptools.dynamic.version).
    """
    version_path = Path(__file__).resolve().parent / "VERSION"
    try:
        return version_path.read_text(encoding="utf-8").strip()
    except Exception:
        return "unknown"


def maybe_print_banner(*, config: BannerConfig | None = None) -> None:
    """
    Print banner once per interpreter session (not on reloads).
    """
    global _ALREADY_SHOWN
    if _ALREADY_SHOWN:
        return
    _ALREADY_SHOWN = True

    cfg = config or DEFAULT_CONFIG
    if not cfg.enabled:
        return

    parts: list[str] = []

    if cfg.show_logo:
        parts.append(ASCII_LOGO)

    if cfg.show_version:
        parts.append(f"v{read_version_from_file()}".rjust(31))
        parts.append(f"{_check_newer_uedge_ver()}")


    if cfg.show_messages:
        parts.append(
            "\nFor citing UEDGE in publications, see CITATION.md "
            "(https://github.com/LLNL/UEDGE/blob/main/CITATION.md)."
        )
        parts.append(
            "Tip: run uedge.cite('version'|'bibtex'|'apa') for a copy-pasteable reference."
        )

    if not parts:
        return

    try:
        cfg.stream.write("\n".join(parts) + "\n")
    except Exception:
        # Never break import due to printing issues
        pass

