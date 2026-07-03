"""Per-pipeline progress bars driven by Snakemake's stdout.

Snakemake reports its plan upfront in a `Job stats:` block listing how many
times each rule will run. We parse that to size one bar per pipeline (rules
are mapped back to pipelines via PIPELINES). Bars are stacked vertically
and advance independently as Snakemake reports `Finished jobid: N (Rule: X)`.

Pipelines that don't appear in the plan (already up to date) don't get a
bar at all — keeps the display tied to the actual work.

Each active pipeline gets a randomly-picked spinner from a curated list,
displayed as the leftmost column. When a pipeline is finished it switches
to a green ●; when it hasn't started yet it shows a grey ○.
"""

from __future__ import annotations

import random
import re

from rich.console import Console
from rich.live import Live
from rich.progress import (
    BarColumn,
    Progress,
    ProgressColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
)
from rich.spinner import SPINNERS, Spinner
from rich.text import Text

from readdiff.pipelines import PIPELINES


# ── Custom spinners (injected before any Spinner(name) call) ─────────


SPINNERS["pulse_dot"] = {
    "interval": 120,
    "frames": ["·", "∙", "•", "●", "•", "∙"],
}
SPINNERS["pulse_circle"] = {
    "interval": 120,
    "frames": ["○", "◌", "⊙", "●", "⊙", "◌"],
}
SPINNERS["braille_pulse"] = {
    "interval": 100,
    "frames": ["⠁", "⠉", "⠋", "⠛", "⠟", "⠿", "⠟", "⠛", "⠋", "⠉"],
}
SPINNERS["sparkle"] = {
    "interval": 100,
    "frames": ["·", "•", "✦", "✶", "✦", "•"],
}
SPINNERS["dot_grow"] = {
    "interval": 100,
    "frames": [".", ":", "•", "●", "◉", "●", "•", ":"],
}
SPINNERS["asterisk"] = {
    "interval": 100,
    "frames": ["·", "+", "*", "✱", "*", "+"],
}

CURATED_SPINNERS = [
    "balloon", "balloon2",
    "pulse_dot", "pulse_circle", "braille_pulse",
    "sparkle", "dot_grow", "asterisk",
]


# ── Status column (spinner / done / idle) ────────────────────────────


class PipelineSpinnerColumn(ProgressColumn):
    """Per-task indicator:
      - the task's own randomly-picked spinner while at least one rule of
        the pipeline is running;
      - [green]●[/] once every rule has finished;
      - [grey]○[/] before the first rule of the pipeline starts.
    """

    def render(self, task) -> Text:
        if task.finished:
            return Text("●", style="green")
        if task.fields.get("active", False):
            spinner: Spinner | None = task.fields.get("spinner")
            if spinner is not None:
                return spinner.render(task.get_time())
        return Text("○", style="grey50")


# ── Snakemake log patterns ───────────────────────────────────────────


RULE_START_RE = re.compile(r"^\s*(?:local)?rule (\w+):")
JOB_FINISH_RE = re.compile(r"Finished (?:job |jobid: )(\d+)(?: \(Rule: (\w+)\))?\.?")
JOB_STATS_LINE_RE = re.compile(r"^\s*([A-Za-z_][\w]*)\s+(\d+)\s*$")


def _build_rule_index() -> dict[str, str]:
    return {rule: pname for pname, info in PIPELINES.items() for rule in info["rules"]}


class MultiPipelineTracker:
    """One progress bar per pipeline involved in the current run."""

    def __init__(self, console: Console | None = None):
        self.console = console or Console()
        self.progress = Progress(
            PipelineSpinnerColumn(),
            TextColumn("[bold]{task.description}[/]"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TextColumn("{task.completed}/{task.total} jobs"),
            TextColumn("[cyan]{task.fields[rule]}[/]"),
            TimeElapsedColumn(),
            console=self.console,
            transient=False,
        )
        self.tasks: dict[str, TaskID] = {}
        self.live: Live | None = None
        self._rule_to_pipeline = _build_rule_index()
        self._pipeline_totals: dict[str, int] = {}
        self._in_job_stats = False
        self._initialized = False

    # ── lifecycle ────────────────────────────────────────────────────

    def _ensure_live(self) -> None:
        if self.live is None:
            self.live = Live(self.progress, console=self.console, refresh_per_second=12)
            self.live.start()

    def _initialize_tasks(self) -> None:
        if self._initialized:
            return
        self._initialized = True
        active_pipelines = [p for p in PIPELINES if self._pipeline_totals.get(p, 0) > 0]
        if not active_pipelines:
            return

        # Random spinner per pipeline, sampled without replacement when possible.
        n = len(active_pipelines)
        if n <= len(CURATED_SPINNERS):
            spinner_names = random.sample(CURATED_SPINNERS, k=n)
        else:
            spinner_names = (
                random.sample(CURATED_SPINNERS, k=len(CURATED_SPINNERS))
                + [random.choice(CURATED_SPINNERS) for _ in range(n - len(CURATED_SPINNERS))]
            )

        self._ensure_live()
        for pname, spinner_name in zip(active_pipelines, spinner_names):
            self.tasks[pname] = self.progress.add_task(
                pname.ljust(5),
                total=self._pipeline_totals[pname],
                rule="—",
                spinner=Spinner(spinner_name, style="cyan"),
                active=False,
            )

    def finish(self) -> None:
        # Force every bar to its declared total (in case parsing missed a line).
        for task in list(self.progress.tasks):
            self.progress.update(
                task.id, completed=task.total or 0, rule="done", active=False
            )
        if self.live is not None:
            self.live.stop()

    # ── parsing ──────────────────────────────────────────────────────

    def parse_line(self, line: str) -> None:
        stripped = line.strip()

        if stripped == "Job stats:":
            self._in_job_stats = True
            return

        if self._in_job_stats:
            if stripped.startswith("total"):
                self._in_job_stats = False
                self._initialize_tasks()
                return
            if (
                not stripped
                or stripped.startswith("---")
                or stripped.startswith("job ")
            ):
                return
            m = JOB_STATS_LINE_RE.match(line)
            if m:
                rule, count = m.group(1), int(m.group(2))
                pipeline = self._rule_to_pipeline.get(rule)
                if pipeline:
                    self._pipeline_totals[pipeline] = (
                        self._pipeline_totals.get(pipeline, 0) + count
                    )
            return

        m = RULE_START_RE.match(line)
        if m:
            rule = m.group(1)
            pipeline = self._rule_to_pipeline.get(rule)
            if pipeline and pipeline in self.tasks:
                # Once a pipeline has started, it stays "active" until the bar
                # reaches 100% (then `task.finished` takes over the rendering).
                self.progress.update(self.tasks[pipeline], rule=rule, active=True)
            return

        m = JOB_FINISH_RE.search(line)
        if m:
            rule = m.group(2)
            if rule:
                pipeline = self._rule_to_pipeline.get(rule)
                if pipeline and pipeline in self.tasks:
                    self.progress.advance(self.tasks[pipeline], 1)
            return

    # ── refresh helpers (verbose 1 vs 2) ─────────────────────────────

    def update_bar(self) -> None:
        if not self.tasks:
            return
        self._ensure_live()
        if self.live is not None:
            self.live.refresh()

    def print_with_bar(self, line: str) -> None:
        self.console.print(line, highlight=False)
