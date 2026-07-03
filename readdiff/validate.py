"""Up-front validation helpers run before any `run X` subprocess.

Each helper either passes silently (returning None) or raises ValidationError
with a user-facing message. The CLI catches ValidationError once, prints it as
a Rich panel, and exits with code 2 (distinct from a Snakemake failure code).
"""

from __future__ import annotations

import shutil
from pathlib import Path


class ValidationError(Exception):
    """User-facing validation error. The message is shown verbatim."""


def check_writable(workdir: Path) -> None:
    try:
        workdir.mkdir(parents=True, exist_ok=True)
        probe = workdir / ".poc_cli_write_probe"
        probe.write_text("ok")
        probe.unlink(missing_ok=True)
    except OSError as e:
        raise ValidationError(f"workdir is not writable: {workdir} ({e})") from None


def check_disk_space(workdir: Path, min_mb: int = 100) -> None:
    target = workdir if workdir.exists() else workdir.parent
    free_mb = shutil.disk_usage(target).free // (1024 * 1024)
    if free_mb < min_mb:
        raise ValidationError(
            f"insufficient disk space at {target}: {free_mb} MB free, "
            f"{min_mb} MB required"
        )


def check_binary(name: str) -> None:
    if shutil.which(name) is None:
        raise ValidationError(
            f"required binary not found in PATH: {name!r}. "
        )


def preflight(workdir: Path) -> None:
    """Run all checks. Raises ValidationError on the first failure."""
    check_writable(workdir)
    check_disk_space(workdir, min_mb=100)
    check_binary("snakemake")
