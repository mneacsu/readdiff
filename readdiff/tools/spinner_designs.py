"""Spinner showcase — pick the one to use in `progress.py`.

Renders many spinners side-by-side, all animating, so you can compare and
choose. Run from main/script_poc_cli/:

    $EXEC python tools/spinner_designs.py [duration_seconds]

Custom spinners (those marked `(custom)`) are appended to Rich's `SPINNERS`
dict at import time so they can be referenced by name like the built-ins.
To adopt one in production, copy its `frames` and `interval` into a similar
SPINNERS injection at the top of `poc_cli/progress.py`, then pass
`spinner_name="<name>"` to the `SpinnerColumn(...)` constructor.

The candidates here are all "single character morphing shape" style — no
multi-dot rotation, in line with the brief.
"""

from __future__ import annotations

import sys
import time

from rich.console import Console
from rich.live import Live
from rich.spinner import SPINNERS, Spinner
from rich.table import Table


# ── Custom spinners (mutate Rich's global SPINNERS dict) ─────────────


SPINNERS["pulse_dot"] = {
    "interval": 120,
    "frames": ["·", "∙", "•", "●", "•", "∙"],
}
SPINNERS["pulse_circle"] = {
    "interval": 120,
    "frames": ["○", "◌", "⊙", "●", "⊙", "◌"],
}
SPINNERS["block_density"] = {
    "interval": 120,
    "frames": ["░", "▒", "▓", "█", "▓", "▒"],
}
SPINNERS["vertical_full"] = {
    "interval": 80,
    "frames": ["▁", "▂", "▃", "▄", "▅", "▆", "▇", "█", "▇", "▆", "▅", "▄", "▃", "▂"],
}
SPINNERS["sparkle"] = {
    "interval": 100,
    "frames": ["·", "•", "✦", "✶", "✦", "•"],
}
SPINNERS["star_twinkle"] = {
    "interval": 100,
    "frames": ["·", "✦", "✶", "✷", "✶", "✦"],
}
SPINNERS["braille_pulse"] = {
    "interval": 100,
    "frames": ["⠁", "⠉", "⠋", "⠛", "⠟", "⠿", "⠟", "⠛", "⠋", "⠉"],
}
SPINNERS["dot_grow"] = {
    "interval": 100,
    "frames": [".", ":", "•", "●", "◉", "●", "•", ":"],
}
SPINNERS["asterisk"] = {
    "interval": 100,
    "frames": ["·", "+", "*", "✱", "*", "+"],
}


# ── Candidate list (label → description) ─────────────────────────────


CANDIDATES = [
    # Built-in
    ("dots",           "rotating braille dots ⠋⠙⠹  (current default)"),
    ("arc",            "rotating arc ◜◠◝◞◡◟"),
    ("circleHalves",   "rotating filled half ◐◓◑◒"),
    ("triangle",       "rotating triangle ◢◣◤◥"),
    ("balloon",        "growing dot . o O @ *"),
    ("balloon2",       "balloon with shrink phase"),
    ("growVertical",   "vertical bar grow ▁▃▄▅▆▇"),
    ("squareCorners",  "rotating filled square ◰◳◲◱"),
    ("circleQuarters", "rotating quarter ◴◷◶◵"),
    # Custom
    ("pulse_dot",      "(custom) dot pulses size · ∙ • ● • ∙"),
    ("pulse_circle",   "(custom) circle fills/empties ○◌⊙●"),
    ("block_density",  "(custom) block density ░▒▓█"),
    ("vertical_full",  "(custom) full-height vertical bar ▁→█→▁"),
    ("sparkle",        "(custom) dot morphs into star · • ✦ ✶"),
    ("star_twinkle",   "(custom) star twinkle ✦ ✶ ✷"),
    ("braille_pulse",  "(custom) braille column pulse ⠁→⠿"),
    ("dot_grow",       "(custom) ascii-friendly dot grow . : • ● ◉"),
    ("asterisk",       "(custom) asterisk pulse · + * ✱"),
]


def main() -> None:
    duration = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0
    console = Console()

    table = Table.grid(padding=(0, 4))
    table.add_column(justify="center", min_width=2)
    table.add_column()
    table.add_column(style="dim")

    for name, desc in CANDIDATES:
        table.add_row(
            Spinner(name, style="cyan"),
            f"[bold]{name}[/]",
            desc,
        )

    console.print(
        f"[dim]Showing {len(CANDIDATES)} spinners for {duration:.0f}s. "
        f"Pick one and tell me which to use in `progress.py`.[/]"
    )
    console.print()

    with Live(table, console=console, refresh_per_second=12, transient=False):
        time.sleep(duration)


if __name__ == "__main__":
    main()
