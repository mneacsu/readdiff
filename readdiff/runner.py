"""Snakemake subprocess launcher with banner, progress, recap and persistent logs.

A `run X` invocation goes through:
  1. Banner (panel with workdir, cores, dependency states).
  2. Snakemake subprocess (silent / streaming with progress / dry-run echo).
  3. Logs always persisted at WORKDIR/.readdiff/runs/run-<pipeline>-<ts>.log.
  4. On failure: red panel with the `Error in rule …:` block and log path.
  5. On success: recap panel listing produced outputs and their sizes.
"""

from __future__ import annotations

import datetime as dt
import re
import subprocess
import time
from pathlib import Path
import os
import json

from rich.console import Console
from rich.panel import Panel

from readdiff.pipelines import PIPELINES
from readdiff.progress import MultiPipelineTracker
from readdiff.status import STATE_GLYPH, STATE_STYLE, quick_status

WORKFLOWS_DIR = Path(__file__).resolve().parent / "workflows"
SNAKEFILE = WORKFLOWS_DIR / "Snakefile"


def run(
    pipeline: str,
    workdir: Path,
    config: dict,
    verbose: int = 0,
    dry_run: bool = False,
    unlock: bool = False,
) -> int:
    """Execute `pipeline` for `workdir`. Return the snakemake exit code."""
    console = Console()
    workdir = workdir.resolve()

    profile = config.get("snakemake_profile")
    cores = int(config.get("snakemake_cores"))

    _show_banner(console, pipeline, workdir, profile, cores, dry_run)

    runs_dir = workdir / ".readdiff" / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = runs_dir / f"run-{pipeline}-{timestamp}.log"

    cmd = _build_cmd(pipeline, workdir, config, profile, cores, dry_run, unlock)

    start = time.time()
    if dry_run:
        rc, full_log = _run_dry(cmd, console)
    elif verbose == 0:
        rc, full_log = _run_silent(cmd)
    else:
        rc, full_log = _run_streaming(cmd, pipeline, console, verbose)
    duration = time.time() - start

    log_path.write_text(full_log)

    if rc != 0:
        _show_error(console, full_log, log_path)
        return rc

    if not dry_run:
        _show_recap(console, pipeline, workdir, log_path, duration)
    return rc


# ── Banner ──────────────────────────────────────────────────────────────


def _show_banner(
    console: Console,
    pipeline: str,
    workdir: Path,
    profile: Path | None,
    cores: int,
    dry_run: bool,
) -> None:
    states = quick_status(workdir)
    deps = PIPELINES[pipeline]["deps"]

    if profile:
        rows = [
            f"[bold]workdir[/] : {workdir}",
            #f"[bold]config[/]  : {config_path or '[dim]embedded default[/]'}",
            f"[bold]profile[/]   : {profile}",
        ]
    else:
        rows = [
            f"[bold]workdir[/] : {workdir}",
            #f"[bold]config[/]  : {config_path or '[dim]embedded default[/]'}",
            f"[bold]cores[/]   : {cores}",
        ]
    if deps:
        chunks = []
        for d in deps:
            s = states[d]
            chunks.append(f"{d} [{STATE_STYLE[s.state]}]{STATE_GLYPH[s.state]}[/]")
        suffix = ""
        if any(states[d].state != "done" for d in deps):
            suffix = "  [dim]→ snakemake will rebuild what's needed[/]"
        rows.append(f"[bold]deps[/]    : {'  '.join(chunks)}{suffix}")
    if dry_run:
        rows.append("[bold yellow]DRY RUN[/] — no files will be written")

    console.print(
        Panel(
            "\n".join(rows),
            title=f"[bold]pipeline: {pipeline}[/]",
            border_style="cyan",
            expand=False,
        )
    )


# ── Snakemake invocation ────────────────────────────────────────────────


def _build_cmd(
    pipeline: str, workdir: Path, config: dict, profile: Path | None, cores: int, dry_run: bool, unlock: bool
) -> list[str]:
    cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--directory", str(workdir),
        "--profile" if profile else "--cores", str(profile) if profile else str(cores),
        pipeline,
    ]
    if dry_run:
        cmd.append("-n")
    if unlock:
        cmd.append("--unlock")

    pairs = [f"workdir={workdir}"]
    for k, v in config.items():
        if k.startswith("snakemake_"):
            continue
        if v in ["None", "none", "null", "", None]:
            continue
        if isinstance(v, (dict, list)):
            v = json.dumps(v, separators=(",", ":"))
        else:
            v = str(v)
        pairs.append(f"{k}={v}")

    cmd += ["--config", *pairs]
    return cmd


def _run_silent(cmd: list[str]) -> tuple[int, str]:
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, (result.stdout or "") + (result.stderr or "")


def _run_streaming(
    cmd: list[str], pipeline: str, console: Console, verbose: int
) -> tuple[int, str]:
    tracker = MultiPipelineTracker(console=console)
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    captured: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        stripped = line.rstrip()
        captured.append(stripped)
        tracker.parse_line(stripped)
        if verbose >= 2:
            tracker.print_with_bar(stripped)
        else:
            tracker.update_bar()
    proc.wait()
    tracker.finish()
    return proc.returncode, "\n".join(captured) + "\n"


def _run_dry(cmd: list[str], console: Console) -> tuple[int, str]:
    result = subprocess.run(cmd, capture_output=True, text=True)
    body = ((result.stdout or "") + (result.stderr or "")).rstrip() or "(empty)"
    console.print(
        Panel(
            body,
            title="[bold yellow]dry-run plan[/]",
            border_style="yellow",
            expand=False,
        )
    )
    return result.returncode, body + "\n"


# ── Error rendering ─────────────────────────────────────────────────────


def _show_error(console: Console, full_log: str, log_path: Path) -> None:
    blocks = _extract_error_blocks(full_log)
    body: list[str] = []
    if blocks:
        for block in blocks:
            body.extend(block)
            body.append("")
    else:
        body.append("(no `Error in rule` block found — showing last 30 lines)")
        body.extend(full_log.splitlines()[-30:])
    body.append("")
    body.append(f"[dim]Full log:[/] {log_path}")
    console.print(
        Panel(
            "\n".join(body),
            title="[bold red]Pipeline failed[/]",
            border_style="red",
            expand=False,
        )
    )


_ERROR_HEADER_RE = re.compile(r"^\s*Error in rule \w+:")


def _extract_error_blocks(log: str) -> list[list[str]]:
    """Extract `Error in rule …:` blocks; snakemake repeats each block in its
    final summary so we de-duplicate by content."""
    blocks: list[list[str]] = []
    seen: set[tuple[str, ...]] = set()
    lines = log.splitlines()
    i = 0
    while i < len(lines):
        if _ERROR_HEADER_RE.match(lines[i]):
            block = [lines[i]]
            i += 1
            saw_content = False
            while i < len(lines):
                if not lines[i].strip():
                    if saw_content:
                        break
                else:
                    saw_content = True
                block.append(lines[i])
                i += 1
            key = tuple(block)
            if key not in seen:
                seen.add(key)
                blocks.append(block)
        else:
            i += 1
    return blocks


# ── Recap ───────────────────────────────────────────────────────────────


def _show_recap(
    console: Console, pipeline: str, workdir: Path, log_path: Path, duration: float
) -> None:
    outputs = PIPELINES[pipeline]["outputs"]

    rows = []
    for rel in outputs:
        path = workdir / rel
        if path.exists():
            size = _fmt_size(path.stat().st_size)
            rows.append(f"[green]●[/] {rel:32s} [dim]{size:>10}[/]")
        else:
            rows.append(f"[red]○ missing[/] {rel}")
    rows.append("")
    try:
        log_rel = log_path.relative_to(workdir)
        rows.append(f"[dim]log:[/] {log_rel}")
    except ValueError:
        rows.append(f"[dim]log:[/] {log_path}")

    console.print(
        Panel(
            "\n".join(rows),
            title=f"[bold]{pipeline} — done in {duration:.1f}s[/]",
            border_style="green",
            expand=False,
        )
    )


def _fmt_size(n_bytes: int) -> str:
    size = float(n_bytes)
    for u in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.1f} {u}" if u != "B" else f"{int(size)} B"
        size /= 1024
    return f"{size:.1f} TB"

def run_explorer(
    workdir: Path,
    config: dict,
) -> int:
    """Run volcano explorer app for `workdir`. Return the exit code."""
    console = Console()
    workdir = workdir.resolve()

    runs_dir = workdir / ".readdiff" / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = runs_dir / f"run-explorer-{timestamp}.log"

    query = (
        config["query_file"] 
        if "query_file" in config and config["query_file"]
        else f"{workdir}/features/random_{config["n_features"]}features_seed{config["random_seed"]}.fasta"
    )
    query_name = os.path.basename(query).removesuffix(".fasta")

    if not os.path.exists(f"{workdir}/differential_testing/{query_name}/deseq2_merged_results.tsv"):
        raise ValueError("Error: DESeq2 results not available. Run the differential testing pipeline first.")
    if not os.path.exists(f"{workdir}/differential_testing/{query_name}/deseq2_norm_counts.h5ad"):
        raise ValueError("Error: DESeq2 normalized counts not available. Run the differential testing pipeline first.")
    if not os.path.exists(f"{workdir}/annotations/{query_name}.analysis_table.tsv"):
        print("Warning: Sequence annotations not available. Proceeding without annotations. If required, run the annotation pipeline to generate them.")
        
    volcano_config = (
        config["volcano_config"] 
        if "volcano_config" in config and config["volcano_config"]  
        else f"{Path(__file__).resolve().parent}/volcano_explorer/volcano_config.yaml"
    )
    if not os.path.isabs(volcano_config):
        volcano_config = os.path.join(workdir, volcano_config)

    print(f"Configuration: {volcano_config}")

    cmd = [
        "python3", f"{Path(__file__).resolve().parent}/volcano_explorer/explore_volcano.py",
        "--config", volcano_config,
        "--counts", f"{workdir}/differential_testing/{query_name}/deseq2_norm_counts.h5ad",
        "--deseq2-results", f"{workdir}/differential_testing/{query_name}/deseq2_merged_results.tsv",
    ]

    if os.path.exists(f"{workdir}/annotations/{query_name}.analysis_table.tsv"):
        cmd.extend(["--annotations", f"{workdir}/annotations/{query_name}.analysis_table.tsv"])

    start = time.time()
    rc, full_log = _run_streaming(cmd, "explorer", console, 2)
    duration = time.time() - start

    log_path.write_text(full_log)

    if rc != 0:
        _show_error(console, full_log, log_path)
        return rc

    if not dry_run:
        _show_recap(console, pipeline, workdir, log_path, duration)
    return rc