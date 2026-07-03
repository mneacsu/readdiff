"""Render pipeline status as a Rich table.

Topology is communicated via the `deps` column rather than ASCII DAG art:
the table format is what most production CLIs (kubectl, docker, terraform)
use for similar multi-row state displays — compact, scannable, easy to
extend with new columns.
"""

from __future__ import annotations

import time
from pathlib import Path

from rich.console import Console
from rich.table import Table
from rich.text import Text

from readdiff.pipelines import PIPELINES
from readdiff.status import STATE_GLYPH, STATE_STYLE, PipelineState


def _fmt_age(mtime: float | None) -> str:
    if mtime is None:
        return "—"
    elapsed = time.time() - mtime
    if elapsed < 60:
        return f"{int(elapsed)}s ago"
    if elapsed < 3600:
        return f"{int(elapsed / 60)}m ago"
    if elapsed < 86400:
        return f"{int(elapsed / 3600)}h ago"
    return f"{int(elapsed / 86400)}d ago"


def _rules_cell(s: PipelineState) -> Text:
    n = len(PIPELINES[s.name]["rules"])
    if s.state == "done":
        return Text(f"{n}/{n} ok", style="dim")
    if s.state == "stale":
        return Text(f"{s.rerun_count}/{n} rerun", style="orange3")
    if s.state == "missing":
        return Text(f"{s.rerun_count}/{n} to run", style="default")
    return Text(f"{s.rerun_count}/{n} queued", style="grey50")


def _outputs_cell(s: PipelineState) -> Text:
    txt = f"{s.outputs_present}/{s.outputs_total}"
    if s.outputs_present == s.outputs_total and s.outputs_total > 0:
        return Text(txt, style="green")
    if s.outputs_present == 0:
        return Text(txt, style="red")
    return Text(txt, style="orange3")


def render(console: Console, states: dict[str, PipelineState], workdir: Path) -> None:
    console.print(f"[dim]workdir:[/] {workdir}")
    console.print()

    table = Table(show_header=True, header_style="bold dim", box=None, padding=(0, 2))
    table.add_column("pipeline")
    table.add_column("state")
    table.add_column("rules", justify="right")
    table.add_column("outputs", justify="right")
    table.add_column("deps")
    table.add_column("last run", justify="right")

    for name, s in states.items():
        deps = PIPELINES[name]["deps"]
        deps_cell = Text(", ".join(deps) if deps else "—", style="dim")
        table.add_row(
            Text(name, style="bold"),
            Text(f"{STATE_GLYPH[s.state]} {s.state}", style=STATE_STYLE[s.state]),
            _rules_cell(s),
            _outputs_cell(s),
            deps_cell,
            Text(_fmt_age(s.last_run), style="dim"),
        )

    console.print(table)
