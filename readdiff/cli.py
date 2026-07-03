"""readdiff — Typer CLI entry point.

All commands take a positional WORKDIR (default `.`). The workdir is the only
location concept: outputs, logs, and the local config file all live inside it.

Commands:
  run index|features|estimate|diff|annotate|explorer [WORKDIR] [-v N] [--dry-run] [--unlock]
  status [WORKDIR]
  clean  [WORKDIR] [--pipeline X] [--all] [--yes]
  config show|get|set|unset|edit|validate [--global]
  version
"""

from __future__ import annotations

import shutil
from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel

from readdiff.__init__ import __version__
from readdiff.config import (
    edit_config,
    get_value,
    resolve_config,
    set_value,
    show_config,
    unset_value,
    validate_config,
)
from readdiff.status_render import render as render_status
from readdiff.pipelines import PIPELINES, pipeline_names
from readdiff.runner import run as run_pipeline, run_explorer as run_explorer_app
from readdiff.status import all_done, status_all
from readdiff.validate import ValidationError, preflight

app = typer.Typer(
    name="readdiff",
    help="Minimizer-based, alignment-last differential expression testing for sequencing reads.",
    add_completion=False,
    no_args_is_help=True,
)
run_app = typer.Typer(help="Run a pipeline.", no_args_is_help=True)
config_app = typer.Typer(help="Configuration management.", no_args_is_help=True)
app.add_typer(run_app, name="run")
app.add_typer(config_app, name="config")


WORKDIR_OPT = typer.Argument(Path("."), help="Workdir (default: current directory).")
VERBOSE_OPT = typer.Option(1, "-v", "--verbose", min=0, max=2,
                           help="0=silent, 1=progress bar, 2=progress+logs")
DRY_RUN_OPT = typer.Option(False, "--dry-run", help="Show DAG plan without executing.")
UNLOCK_OPT = typer.Option(False, "--unlock", help="Unlock working directory.")
GLOBAL_OPT = typer.Option(False, "--global",
                          help="Target ~/.config/readdiff/config.yaml instead of the workdir.")


# ── run ─────────────────────────────────────────────────────────────────


def _do_run(pipeline: str, workdir: Path, verbose: int, dry_run: bool, unlock: bool) -> None:
    console = Console()
    workdir = workdir.resolve()
    try:
        preflight(workdir)
    except ValidationError as e:
        console.print(Panel(str(e), title="[bold red]Validation failed[/]", border_style="red"))
        raise typer.Exit(code=2)

    cfg = resolve_config(workdir)
    rc = run_pipeline(
        pipeline=pipeline,
        workdir=workdir,
        config=dict(cfg),
        verbose=verbose,
        dry_run=dry_run,
        unlock=unlock,
    )
    raise typer.Exit(code=rc)

def _do_run_explorer(workdir: Path) -> None:
    console = Console()
    workdir = workdir.resolve()
    try:
        preflight(workdir)
    except ValidationError as e:
        console.print(Panel(str(e), title="[bold red]Validation failed[/]", border_style="red"))
        raise typer.Exit(code=2)

    cfg = resolve_config(workdir)
    rc = run_explorer_app(
        workdir=workdir,
        config=dict(cfg),

    )
    raise typer.Exit(code=rc)


@run_app.command("index")
def run_index(
    workdir: Path = WORKDIR_OPT,
    verbose: int = VERBOSE_OPT,
    dry_run: bool = DRY_RUN_OPT,
    unlock: bool = UNLOCK_OPT,
):
    """Build Needle index"""
    _do_run("index", workdir, verbose, dry_run, unlock)

@run_app.command("features")
def run_features(
    workdir: Path = WORKDIR_OPT,
    verbose: int = VERBOSE_OPT,
    dry_run: bool = DRY_RUN_OPT,
    unlock: bool = UNLOCK_OPT,
):
    """Perform random feature sampling"""
    _do_run("features", workdir, verbose, dry_run, unlock)

@run_app.command("estimate")
def run_estimate(
    workdir: Path = WORKDIR_OPT,
    verbose: int = VERBOSE_OPT,
    dry_run: bool = DRY_RUN_OPT,
    unlock: bool = UNLOCK_OPT,
):
    """Estimate query sequences"""
    _do_run("estimate", workdir, verbose, dry_run, unlock)

@run_app.command("diff")
def run_diff(
    workdir: Path = WORKDIR_OPT,
    verbose: int = VERBOSE_OPT,
    dry_run: bool = DRY_RUN_OPT,
    unlock: bool = UNLOCK_OPT,
):
    """Test for differences between conditions"""
    _do_run("diff", workdir, verbose, dry_run, unlock)

@run_app.command("annotate")
def run_annotate(
    workdir: Path = WORKDIR_OPT,
    verbose: int = VERBOSE_OPT,
    dry_run: bool = DRY_RUN_OPT,
    unlock: bool = UNLOCK_OPT,
):
    """Annotate differentially expressed features"""
    _do_run("annotate", workdir, verbose, dry_run, unlock)

@run_app.command("explorer")
def run_explorer(
    workdir: Path = WORKDIR_OPT,
):
    """Run volcano explorer app"""
    _do_run_explorer(workdir)


# ── status ──────────────────────────────────────────────────────────────


@app.command("status")
def cmd_status(workdir: Path = WORKDIR_OPT):
    """Show pipeline state for WORKDIR. Exit 0 iff every pipeline is done."""
    workdir = workdir.resolve()
    console = Console()
    states = status_all(workdir)
    render_status(console, states, workdir)
    raise typer.Exit(code=0 if all_done(states) else 1)


# ── clean ───────────────────────────────────────────────────────────────


@app.command("clean")
def cmd_clean(
    workdir: Path = WORKDIR_OPT,
    pipeline: str = typer.Option(None, "--pipeline", "-p", help="Clean only this pipeline."),
    all_: bool = typer.Option(False, "--all", help="Also remove .readdiff/ (logs, snakemake cache)."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation and actually delete."),
):
    """Remove pipeline outputs. Dry-run by default; pass --yes to execute."""
    workdir = workdir.resolve()
    console = Console()

    if pipeline is not None and pipeline not in PIPELINES:
        raise typer.BadParameter(f"unknown pipeline {pipeline!r}; choose from {pipeline_names()}")

    pipelines = [pipeline] if pipeline else pipeline_names()
    targets: list[Path] = [workdir / p for p in pipelines]
    if all_:
        targets.append(workdir / ".readdiff")
        targets.append(workdir / ".snakemake")

    existing = [t for t in targets if t.exists()]
    if not existing:
        console.print("[dim]Nothing to clean.[/]")
        raise typer.Exit(0)

    console.print("[bold]Targets:[/]")
    for t in existing:
        console.print(f"  • {t}")

    if not yes:
        console.print("\n[yellow]Dry-run.[/] Re-run with [bold]--yes[/] to actually delete.")
        raise typer.Exit(0)

    for t in existing:
        shutil.rmtree(t)
        console.print(f"[red]removed[/] {t}")


# ── version ─────────────────────────────────────────────────────────────


@app.command("version")
def cmd_version():
    """Show CLI version and the snakemake binary in PATH."""
    console = Console()
    sm = shutil.which("snakemake")
    console.print(f"[bold]readdiff[/] {__version__}")
    console.print(f"snakemake: {sm if sm else '[red]not found[/]'}")


# ── config sub-app ──────────────────────────────────────────────────────


@config_app.command("show")
def cmd_config_show(workdir: Path = WORKDIR_OPT):
    """Display the active config (after fallback) plus other candidates."""
    show_config(workdir.resolve())


@config_app.command("get")
def cmd_config_get(
    key: str = typer.Argument(..., help="Configuration key"),
    workdir: Path = WORKDIR_OPT,
):
    """Read a value from the active config."""
    val = get_value(workdir.resolve(), key)
    if val is None:
        typer.echo(f"key {key!r} not found", err=True)
        raise typer.Exit(1)
    typer.echo(val)


@config_app.command("set")
def cmd_config_set(
    key: str = typer.Argument(..., help="Configuration key"),
    value: str = typer.Argument(..., help="Value to set (auto-cast int/float/bool)"),
    global_: bool = GLOBAL_OPT,
    workdir: Path = WORKDIR_OPT,
):
    """Set KEY to VALUE. Targets the workdir config by default, --global for ~/.config."""
    path = set_value(workdir.resolve(), key, value, global_)
    typer.echo(f"Set {key} = {value} in {path}")


@config_app.command("unset")
def cmd_config_unset(
    key: str = typer.Argument(..., help="Configuration key to remove"),
    global_: bool = GLOBAL_OPT,
    workdir: Path = WORKDIR_OPT,
):
    """Remove a key from the workdir or global config."""
    path, was_present = unset_value(workdir.resolve(), key, global_)
    if was_present:
        typer.echo(f"Removed {key} from {path}")
    else:
        typer.echo(f"Key {key!r} was not present in {path}")


@config_app.command("edit")
def cmd_config_edit(global_: bool = GLOBAL_OPT, workdir: Path = WORKDIR_OPT):
    """Open the workdir or global config in $EDITOR (defaults to vi). Auto-creates."""
    path = edit_config(workdir.resolve(), global_)
    typer.echo(f"Edited {path}")


@config_app.command("validate")
def cmd_config_validate(workdir: Path = WORKDIR_OPT):
    """Check the active config for unknown keys or type mismatches."""
    console = Console()
    ok, msgs = validate_config(workdir.resolve())
    if ok:
        console.print(f"[green]✓[/] config is valid")
        raise typer.Exit(0)
    console.print(f"[red]✗[/] config has issues:")
    for m in msgs:
        console.print(f"  • {m}")
    raise typer.Exit(1)
