"""YAML configuration with fallback resolution.

Resolution order (first match wins, no merging):
  1. Project : <WORKDIR>/config.yaml
  2. Global  : ~/.config/readdiff/config.yaml
  3. Default : readdiff/default_config.yaml (embedded resource)

Write commands target one scope at a time:
  - default scope = project (workdir)
  - --global scope = ~/.config/readdiff/

Files are auto-created on the first write/edit; there is no `init` step.
The default file is read-only (embedded). `validate` checks the active
config against the default's keys and types.
"""

from __future__ import annotations

import os
import subprocess
from importlib.resources import files as pkg_files
from pathlib import Path

import yaml
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

GLOBAL_CONFIG_DIR = Path.home() / ".config" / "readdiff"
GLOBAL_CONFIG_PATH = GLOBAL_CONFIG_DIR / "config.yaml"
PROJECT_CONFIG_NAME = "config.yaml"


def workdir_config_path(workdir: Path) -> Path:
    return workdir / PROJECT_CONFIG_NAME


def load_default() -> dict:
    text = pkg_files("readdiff").joinpath("default_config.yaml").read_text()
    data = yaml.safe_load(text)
    return data if isinstance(data, dict) else {}


def _load_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text())
    return data if isinstance(data, dict) else {}


def resolve_config(workdir: Path) -> tuple[list[str], Path | None, dict]:
    """Return (scope, path, config_dict). scope ∈ {project, global, default}.

    `path` is None when the embedded default resource is used.
    """
    config = load_default()

    if GLOBAL_CONFIG_PATH.exists():
        config.update(_load_yaml(GLOBAL_CONFIG_PATH))

    wpath = workdir_config_path(workdir)
    if wpath.exists():
        config.update(_load_yaml(wpath))

    return config


def get_value(workdir: Path, key: str):
    cfg = resolve_config(workdir)
    return cfg.get(key)


def _target_path(workdir: Path, is_global: bool) -> Path:
    return GLOBAL_CONFIG_PATH if is_global else workdir_config_path(workdir)


def _cast_scalar(value: str):
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    if value.lower() in ("true", "false"):
        return value.lower() == "true"
    return value


def set_value(workdir: Path, key: str, value: str, is_global: bool) -> Path:
    """Write KEY=VALUE to the project (default) or global config. Auto-creates."""
    path = _target_path(workdir, is_global)
    cfg = _load_yaml(path)
    cfg[key] = _cast_scalar(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f, default_flow_style=False, sort_keys=True)
    return path


def unset_value(workdir: Path, key: str, is_global: bool) -> tuple[Path, bool]:
    """Remove KEY from the chosen scope. Returns (path, was_present).

    No-op if the file doesn't exist or the key wasn't there — does not
    create empty config files.
    """
    path = _target_path(workdir, is_global)
    if not path.exists():
        return path, False
    cfg = _load_yaml(path)
    if key not in cfg:
        return path, False
    cfg.pop(key)
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f, default_flow_style=False, sort_keys=True)
    return path, True


def edit_config(workdir: Path, is_global: bool) -> Path:
    """Open the chosen scope's config in $EDITOR. Auto-creates with a stub."""
    path = _target_path(workdir, is_global)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# readdiff configuration\n")
    editor = os.environ.get("EDITOR", "vi")
    subprocess.call([editor, str(path)])
    return path


def validate_config(workdir: Path) -> tuple[bool, str, list[str]]:
    """Return (ok, messages). Validates the *active* config."""
    cfg = resolve_config(workdir)
    defaults = load_default()
    known = set(defaults.keys())
    msgs: list[str] = []

    for k in cfg:
        if k not in known:
            msgs.append(f"unknown key: {k!r}")

    for k, v in cfg.items():
        if k not in defaults:
            continue
        expected = type(defaults[k])
        if isinstance(defaults[k], (int, float)) and isinstance(v, (int, float)) and not isinstance(v, bool):
            continue
        if not isinstance(v, expected):
            msgs.append(f"type mismatch for {k!r}: expected {expected.__name__}, got {type(v).__name__}")

    return (len(msgs) == 0), msgs


# ── Display ──────────────────────────────────────────────────────────


def show_config(workdir: Path) -> None:
    console = Console()
    cfg = resolve_config(workdir)

    header = ["[dim]Config files:[/]",]
    wpath = workdir_config_path(workdir)
    for label, p in [("project", wpath), ("global", GLOBAL_CONFIG_PATH), ("default", "readdiff/default_config.yaml")]:
        marker = "[green]✓ present[/]" if os.path.exists(p) else "[dim]— absent[/]"
        header.append(f"  {label:7s} : {p}  {marker}")


    console.print(Panel("\n".join(header), title="readdiff config", border_style="cyan"))

    if not cfg:
        console.print("[dim](no values; falling back to next scope)[/]")
        return
    
    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("key", style="cyan")
    table.add_column("value", style="green")
    for k in sorted(cfg):
        table.add_row(k, str(cfg[k]))
    console.print(table)
