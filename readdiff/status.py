"""Compute per-pipeline state for `status` and the run banner.

States and their visual encoding:

  done    — ● green  : all outputs present, dependencies done, no rerun planned
  stale   — ● orange : outputs present but Snakemake would rerun something here
                       (input newer than output, or upstream pipeline stale)
  missing — ○ default: at least one output absent; dependencies are done so
                       this pipeline can run now
  pending — ○ dim    : at least one dependency not done; can't run yet

`done`/`stale` share the filled glyph (the work was once finished); `missing`
and `pending` share the empty glyph (the work isn't finished). Color carries
the secondary signal (OK vs needs redo, your turn vs waiting).

Stale detection is delegated to Snakemake itself: we run `snakemake -n delta`
once per `status` invocation and parse which rules would rerun. This matches
exactly what running the pipeline would do.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from readdiff.pipelines import PIPELINES

WORKFLOWS_DIR = Path(__file__).resolve().parent.parent / "workflows"
SNAKEFILE = WORKFLOWS_DIR / "Snakefile"

STATE_GLYPH = {
    "done":    "●",
    "stale":   "●",
    "missing": "○",
    "pending": "○",
}

STATE_STYLE = {
    "done":    "bold green",
    "stale":   "bold orange3",
    "missing": "default",
    "pending": "grey50",
}


@dataclass
class PipelineState:
    name: str
    state: str              # done | stale | missing | pending
    detail: str             # one-line why
    last_run: float | None  # mtime of most recent output, or None
    outputs_present: int = 0
    outputs_total: int = 0
    rerun_count: int = 0    # rules that snakemake would (re)run for this pipeline


_RULE_LINE_RE = re.compile(r"^\s*(?:local)?rule (\w+):")


def _planned_rules(workdir: Path, config: dict) -> set[str] | None:
    """Run `snakemake -n delta` and return the set of rule names that would
    execute. Returns None if Snakemake is unavailable or the call fails.

    Targeting `delta` plans the entire DAG (cumulative target rules ensure
    every upstream pipeline is also considered).
    """
    if not SNAKEFILE.exists():
        return None
    cmd = [
        "snakemake",
        "--snakefile", str(SNAKEFILE),
        "--directory", str(workdir),
        "-n",
        "delta",
    ]
    config_pairs = [f"workdir={workdir}"]
    for k, v in config.items():
        if not k.startswith("snakemake_"):
            config_pairs.append(f"{k}={v}")
    cmd += ["--config", *config_pairs]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None
    if result.returncode != 0:
        return None
    rules: set[str] = set()
    for line in (result.stdout or "").splitlines():
        m = _RULE_LINE_RE.match(line)
        if m:
            rules.add(m.group(1))
    return rules


def compute_one(
    workdir: Path,
    name: str,
    upstream: dict[str, PipelineState],
    planned: set[str] | None,
) -> PipelineState:
    info = PIPELINES[name]
    outputs = [workdir / p for p in info["outputs"]]
    missing = [p for p in outputs if not p.exists()]

    deps = info["deps"]
    deps_done = all(upstream[d].state == "done" for d in deps)

    n_total = len(outputs)
    n_present = n_total - len(missing)
    rerun = len(set(info["rules"]) & planned) if planned else 0

    common = dict(
        outputs_present=n_present,
        outputs_total=n_total,
        rerun_count=rerun,
    )

    if missing:
        if not deps_done:
            return PipelineState(name, "pending", f"{len(deps)} dep(s) not ready", None, **common)
        return PipelineState(
            name, "missing",
            f"{len(missing)}/{n_total} output(s) missing",
            None,
            **common,
        )

    own_max = max(p.stat().st_mtime for p in outputs)

    dep_stale = any(upstream[d].state == "stale" for d in deps)
    self_planned = planned is not None and bool(set(info["rules"]) & planned)

    if dep_stale or self_planned:
        detail = "upstream stale" if dep_stale else "snakemake would rerun"
        return PipelineState(name, "stale", detail, own_max, **common)

    return PipelineState(name, "done", "ready", own_max, **common)


def compute_states(
    workdir: Path, planned: set[str] | None = None
) -> dict[str, PipelineState]:
    """Compute state for every pipeline. PIPELINES dict order respects deps."""
    out: dict[str, PipelineState] = {}
    for name in PIPELINES:
        out[name] = compute_one(workdir, name, out, planned)
    return out


def status_all(workdir: Path) -> dict[str, PipelineState]:
    """Full status with Snakemake-based stale detection. Used by `status`."""
    from readdiff.config import resolve_config

    cfg = resolve_config(workdir)
    planned = _planned_rules(workdir, cfg)
    return compute_states(workdir, planned)


def quick_status(workdir: Path) -> dict[str, PipelineState]:
    """Fast file-existence-only status. Used by the run banner."""
    return compute_states(workdir, planned=None)


def all_done(states: dict[str, PipelineState]) -> bool:
    return all(s.state == "done" for s in states.values())
