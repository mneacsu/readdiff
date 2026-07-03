#!/usr/bin/env python3
"""Interactive volcano plot explorer dashboard — Clean Scientific theme.

Launches a local web server that serves an interactive volcano plot dashboard
for exploring DESeq2 results stored in annotated h5ad files.

Usage:
    python explore_volcano.py --config config.yaml
    # Then open http://localhost:8050 in your browser
"""

import argparse
import atexit
import io
import os

import anndata
import dash
import dash_mantine_components as dmc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import yaml
from dash import (ALL, ClientsideFunction, Dash, Input, Output, Patch, State,
                  callback_context, dcc, html)


# ─── Constants ────────────────────────────────────────────────────────────
DE_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "factor", "test"]

NUCLEOTIDE_COLORS_DEFAULT = {
    "A": "#86efac",  # pastel green  (light theme)
    "T": "#fca5a5",  # pastel red
    "U": "#fca5a5",
    "C": "#93c5fd",  # pastel blue
    "G": "#fde68a",  # pastel yellow
    "default": "#e5e7eb",
}
NUCLEOTIDE_COLORS_DEFAULT_DARK = {
    "A": "#16a34a",  # deeper green   (dark theme)
    "T": "#dc2626",  # deeper red
    "U": "#dc2626",
    "C": "#2563eb",  # deeper blue
    "G": "#ca8a04",  # deeper amber
    "default": "#52525b",
}

# Hand-curated discrete colormaps for the volcano plot.
# Plotly's qualitative module only provides Dark24/Light24 (24) and Alphabet
# (26) beyond ~12 colors, so we add Polychrome36 (Trubetskoy / Glasbey-style),
# Kelly20 (Kelly 1965 max-contrast set) and a 48-color combo for the largest
# class counts.
_POLYCHROME36 = [
    "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
    "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
    "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
    "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0",
    "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5",
    "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB",
]
_KELLY20 = [
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
    "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97",
    "#F6A600", "#B3446C", "#DCD300", "#882D17", "#8DB600", "#654522",
    "#E25822", "#2B3D26",
]

DISCRETE_PALETTES = {
    "Dark24":     list(px.colors.qualitative.Dark24),     # 24, deep   — default
    "Light24":    list(px.colors.qualitative.Light24),    # 24, light
    "Alphabet":   list(px.colors.qualitative.Alphabet),   # 26, max distinct
    "Polychrome": _POLYCHROME36,                            # 36
    "Kelly":      _KELLY20,                                 # 20
    "Dark+Light": (list(px.colors.qualitative.Dark24)
                    + list(px.colors.qualitative.Light24)),  # 48
}

# Feather-style icons (24×24 viewBox), rendered via CSS mask in _svg_icon
ICON_SEARCH    = "M21 21l-6-6m2-5a7 7 0 1 1-14 0 7 7 0 0 1 14 0Z"
ICON_BOOKMARK  = "M19 21l-7-5-7 5V5a2 2 0 0 1 2-2h10a2 2 0 0 1 2 2v16Z"
ICON_PROJECT   = "M3 9l9-7 9 7v11a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V9Z"
ICON_SETTINGS  = "M12 15a3 3 0 1 0 0-6 3 3 0 0 0 0 6Z"
ICON_HELP      = "M9.1 9a3 3 0 1 1 5.8 1c0 2-3 3-3 3M12 17h0M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Z"
ICON_SUN       = "M12 3v2m0 14v2m9-9h-2M5 12H3m15.4-6.4l-1.4 1.4M7 17l-1.4 1.4m0-12.8L7 7m10 10l1.4 1.4M16 12a4 4 0 1 1-8 0 4 4 0 0 1 8 0Z"
ICON_TRASH     = "M3 6h18M8 6V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v2M6 6l1 14a2 2 0 0 0 2 2h6a2 2 0 0 0 2-2l1-14"
ICON_CARET     = "M6 9l6 6 6-6"


def _svg_icon(d, size=17):
    """Inline SVG icon via CSS mask — picks up `currentColor`."""
    url = (
        "url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' "
        "viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='1.6' "
        f"stroke-linecap='round' stroke-linejoin='round'><path d='{d}'/></svg>\")"
    )
    return html.Span(style={
        "width": f"{size}px", "height": f"{size}px",
        "display": "inline-block", "backgroundColor": "currentColor",
        "WebkitMaskImage": url, "maskImage": url,
        "WebkitMaskRepeat": "no-repeat", "maskRepeat": "no-repeat",
        "WebkitMaskSize": "contain", "maskSize": "contain",
        "WebkitMaskPosition": "center", "maskPosition": "center",
    })


# ─── Helpers ──────────────────────────────────────────────────────────────
def _render_filter_table(history, active_index=None):
    if not history:
        return dmc.Text("No active filter.", size="xs", c="dimmed")
    rows = []
    for i in reversed(range(len(history))):
        entry = history[i]
        is_active = (i == active_index) if active_index is not None else (i == len(history) - 1)
        rows.append(html.Div(
            className="cs-history-row" + (" active" if is_active else ""),
            children=[
                html.Button(
                    [
                        html.Span(f"{i + 1}", className="num"),
                        html.Span(entry["desc"], className="desc",
                                  title=entry["desc"]),
                        html.Span(f'{entry["n_points"]:,}', className="count"),
                    ],
                    id={"type": "filter-row", "index": i},
                    className="cs-history-restore",
                    title="Click to restore this state",
                ),
                html.Button(
                    _svg_icon(ICON_TRASH, 11),
                    id={"type": "filter-row-delete", "index": i},
                    className="cs-btn-danger cs-history-delete",
                    title="Remove this filter from history",
                ),
            ],
        ))
    return html.Div(rows, className="cs-history")


def _apply_text_filter(df, annotation_cols, value, relation, column):
    if not value:
        return df["feature"].tolist()
    target_cols = annotation_cols + ["feature"] if column == "any" else [column]
    target_cols = [c for c in target_cols if c in df.columns]
    try:
        num_val = float(value)
        is_numeric = True
    except ValueError:
        num_val = None
        is_numeric = False
    negate = relation in ("not contains", "!=")
    if negate and column == "any":
        mask = pd.Series(True, index=df.index)
        for col in target_cols:
            s = df[col]
            if relation == "not contains":
                mask &= ~s.astype(str).str.contains(value, case=False, na=False)
            elif relation == "!=":
                mask &= s.astype(str) != value
    else:
        mask = pd.Series(False, index=df.index)
        for col in target_cols:
            s = df[col]
            if relation == "contains":
                mask |= s.astype(str).str.contains(value, case=False, na=False)
            elif relation == "not contains":
                mask |= ~s.astype(str).str.contains(value, case=False, na=False)
            elif relation == "=":
                mask |= s.astype(str) == value
            elif relation == "!=":
                mask |= s.astype(str) != value
            elif relation in (">=", "<=") and is_numeric:
                nc = pd.to_numeric(s, errors="coerce")
                mask |= (nc >= num_val) if relation == ">=" else (nc <= num_val)
    return df.loc[mask, "feature"].tolist()


def _render_subfilters(subfilters, annotation_cols):
    """Render the dynamic list of subfilter groups for the current filter.

    Each subfilter has its own column / relation / value inputs (pattern-
    matched IDs ``sf-col``, ``sf-rel``, ``sf-val``) plus a trash button
    (only when there are 2+ subfilters) and a ``+`` button to add another
    OR-clause.
    """
    if not subfilters:
        subfilters = [{"column": "any", "relation": "contains", "value": ""}]
    n = len(subfilters)
    rows = []
    for i, sf in enumerate(subfilters):
        if i > 0:
            rows.append(html.Div("or", className="cs-or-divider"))
        line2 = [
            dmc.TextInput(
                id={"type": "sf-val", "i": i},
                placeholder="e.g. EGFR…",
                value=sf.get("value", "") or "",
                size="xs",
                style={"flex": 1, "minWidth": 0},
            ),
        ]
        if n > 1:
            line2.append(html.Button(
                _svg_icon(ICON_TRASH, 11),
                id={"type": "sf-del", "i": i},
                className="cs-btn-danger cs-sf-del",
                title="Remove this subfilter",
            ))
        line2.append(html.Button(
            "+",
            id={"type": "sf-add", "i": i},
            className="cs-action-btn cs-sf-add",
            title="Add another subfilter (OR with existing ones)",
        ))
        rows.append(html.Div(className="cs-subfilter", children=[
            html.Div(className="cs-sf-line", children=[
                dmc.Select(
                    id={"type": "sf-col", "i": i},
                    data=[{"label": "any", "value": "any"},
                          {"label": "feature", "value": "feature"}]
                          + [{"label": c, "value": c} for c in annotation_cols],
                    value=sf.get("column") or "any",
                    searchable=True, size="xs",
                    style={"flex": 1, "minWidth": 0},
                ),
                dmc.Select(
                    id={"type": "sf-rel", "i": i},
                    data=["contains", "not contains", "=", "!=", ">=", "<="],
                    value=sf.get("relation") or "contains",
                    size="xs", w=90,
                ),
            ]),
            html.Div(className="cs-sf-line", children=line2),
        ]))
    return rows


def _filtered_volcano_df(condition, selected_reads, show_all,
                          var_df, annotation_cols, config):
    """Build the volcano DF for the active condition then apply the
    history scope (``selected_reads``) and the show_all toggle.
    Returned DF is the dataset on which any subfilter matching runs.
    """
    df = build_volcano_df(var_df, annotation_cols, condition, config)
    if selected_reads is not None:
        df = df[df["feature"].isin(selected_reads)]
    if not show_all:
        df = df[df["regulation"] != "Not significant"]
    return df


def _compute_subfilter_matches(df, sfs, annotation_cols):
    """Apply OR across subfilters on a pre-filtered df.

    Each subfilter is a (column, relation, value) tuple; empty values are
    skipped. Returns the list of read names that match *any* of the
    non-empty subfilters, preserving DF order.
    """
    if not sfs:
        return []
    matched = set()
    for sf in sfs:
        v = (sf.get("value") or "").strip()
        if not v:
            continue
        m = _apply_text_filter(df, annotation_cols, v,
                                sf.get("relation"), sf.get("column"))
        matched.update(m)
    return [r for r in df["feature"].tolist() if r in matched]


def _build_subfilters(cols, rels, vals):
    """Zip pattern-matched State arrays into a list of subfilter dicts."""
    n = max(len(cols or []), len(rels or []), len(vals or []))
    out = []
    for i in range(n):
        out.append({
            "column": cols[i] if i < len(cols or []) else "any",
            "relation": rels[i] if i < len(rels or []) else "contains",
            "value": vals[i] if i < len(vals or []) else "",
        })
    return out


def _has_nonempty_subfilter(sfs):
    return any((sf.get("value") or "").strip() for sf in sfs or [])


def _normalize_text_subfilters(params):
    """Return the list of subfilters for a ``text`` row, supporting both the
    new ``subfilters`` schema and the legacy single (column/relation/value).
    """
    if params and "subfilters" in params:
        return params["subfilters"] or []
    return [{
        "column": (params or {}).get("column"),
        "relation": (params or {}).get("relation"),
        "value": (params or {}).get("value"),
    }]


def _compute_filter_matches(sfs, condition, selected_reads, show_all,
                             var_df, annotation_cols, config):
    """High-level: build the filtered df + apply OR across subfilters."""
    df = _filtered_volcano_df(condition, selected_reads, show_all,
                               var_df, annotation_cols, config)
    return _compute_subfilter_matches(df, sfs, annotation_cols)


def _compute_legend_matches(color_col, visible_classes, max_legend,
                              condition, selected_reads, show_all,
                              var_df, annotation_cols, config):
    """Re-apply a legend filter (set of visible class names on a color column).

    Replays the same max_legend_classes binning as ``build_volcano_figure``
    so a legend filter row can be replayed deterministically when an
    earlier row is removed from history.
    """
    df = build_volcano_df(var_df, annotation_cols, condition, config)
    if selected_reads is not None:
        df = df[df["feature"].isin(selected_reads)]
    if not show_all:
        df = df[df["regulation"] != "Not significant"]
    if color_col not in df.columns:
        return df["feature"].tolist()
    if color_col != "regulation" and max_legend and max_legend > 0:
        top = df[color_col].value_counts().head(max_legend).index
        binned = df[color_col].where(df[color_col].isin(top), other="Other")
    else:
        binned = df[color_col]
    visible_str = [str(v) for v in (visible_classes or [])]
    mask = binned.astype(str).isin(visible_str)
    return df.loc[mask, "feature"].tolist()


def _format_filter_desc(params):
    """Recreate the human-readable desc string from a row's params dict.

    Used by ``_recompute_history_chain`` so deleted rows recompute desc
    consistently with how new filters render them. Text filters can have
    multiple subfilters joined by OR.
    """
    if not params:
        return "?"
    t = params.get("type")
    if t == "text":
        sfs = _normalize_text_subfilters(params)
        parts = []
        for sf in sfs:
            v = (sf.get("value") or "").strip()
            if not v:
                continue
            col = sf.get("column") or "any"
            rel = sf.get("relation") or "contains"
            parts.append(f'{col} {rel} "{v}"')
        return " OR ".join(parts) if parts else "?"
    if t == "legend":
        col = params.get("color_column", "?")
        visible = params.get("visible_classes") or []
        ml = params.get("max_legend") or 20
        half = max(1, (ml if ml > 0 else 20) // 2)
        if len(visible) == 0:
            return f"{col}: none"
        if len(visible) == 1:
            return f"{col} = {visible[0]}"
        if len(visible) <= half:
            return f"{col} IN ({', '.join(map(str, visible))})"
        return f"{col} IN ({len(visible)} classes)"
    return "?"


def _recompute_history_chain(rows, condition, show_all,
                              var_df, annotation_cols, config):
    """Re-apply each row's filter from scratch, chaining through the list.

    Starts from the full dataset (selected_reads=None) and replays each
    filter using its stored ``params``. Rows missing ``params`` (legacy
    saved-views, etc.) keep their stored ``reads`` snapshot. Returns a
    new list of rows with refreshed ``reads`` and ``n_points``.
    """
    selected = None
    out = []
    for r in rows or []:
        params = r.get("params") if isinstance(r, dict) else None
        if not params:
            row = dict(r) if isinstance(r, dict) else {}
            selected = row.get("reads", selected)
            out.append(row)
            continue
        t = params.get("type")
        if t == "text":
            sfs = _normalize_text_subfilters(params)
            matches = _compute_filter_matches(
                sfs, condition, selected, show_all,
                var_df, annotation_cols, config,
            )
        elif t == "legend":
            matches = _compute_legend_matches(
                params.get("color_column", "regulation"),
                params.get("visible_classes") or [],
                params.get("max_legend", -1),
                condition, selected, show_all,
                var_df, annotation_cols, config,
            )
        else:
            matches = list(selected) if selected is not None else []
        row = {
            "desc": _format_filter_desc(params),
            "n_points": len(matches),
            "features": matches,
            "params": params,
        }
        out.append(row)
        selected = matches
    return out


def _persist_saved_views(config_path, saved_views):
    """Re-write the YAML config with an updated saved_views list (atomic)."""
    if not config_path:
        return False
    try:
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        cfg["saved_views"] = saved_views
        tmp = config_path + ".tmp"
        with open(tmp, "w") as f:
            yaml.safe_dump(cfg, f, default_flow_style=False, sort_keys=False)
        os.replace(tmp, config_path)
        return True
    except Exception as e:
        print(f"  Warning: failed to persist saved_views: {e}")
        return False


# ─── Data loading ─────────────────────────────────────────────────────────
def identify_columns(var_columns):
    de_cols = set(var_columns).intersection(DE_COLUMNS)

    return [c for c in var_columns if c not in de_cols], de_cols


def strip_digit_suffix(s):
    if len(s) >= 2 and s[-2] == "_" and s[-1].isdigit():
        return s[:-2]
    return s


def load_data(config):
    deseq2_h5ad = config["deseq2_h5ad"]
    deseq2_results = config["deseq2_results"]
    analysis_table_path = config["analysis_table"]
    blacklist = set(config.get("blacklist_columns", []))

    adata = anndata.read_h5ad(deseq2_h5ad, backed="r")

    var_df = pd.read_csv(deseq2_results, sep='\t')

    if analysis_table_path:
        at = pd.read_csv(analysis_table_path, sep="\t", index_col=0)
        var_df = var_df.merge(at, left_on="feature", right_index=True, how="left")
    
    var_df.index = var_df["feature"]

    n_before = len(var_df["feature"].unique())
    for col, vals in (config.get("filter_var_keep") or {}).items():
        if col not in var_df.columns:
            print(f"  Warning: filter_var_keep column '{col}' not found, skipping")
            continue
        sv = [str(v) for v in vals]
        mask = var_df[col].astype(str).isin(sv)
        if "nan" in sv:
            mask = mask | pd.isna(var_df[col])
        var_df = var_df[mask]
    for col, vals in (config.get("filter_var_exclude") or {}).items():
        if col not in var_df.columns:
            print(f"  Warning: filter_var_exclude column '{col}' not found, skipping")
            continue
        sv = [str(v) for v in vals]
        mask = ~var_df[col].astype(str).isin(sv)
        if "nan" in sv:
            mask = mask & ~pd.isna(var_df[col])
        var_df = var_df[mask]
    if len(var_df["feature"].unique()) < n_before:
        print(f"  Var filtering: {n_before} → {len(var_df["feature"].unique())} features")

    obs_columns = [c for c in adata.obs.columns if c not in blacklist]
    layer_names = list(adata.layers.keys()) if adata.layers else []
    all_idx = {n: i for i, n in enumerate(adata.var_names)}
    var_name_to_idx = {n: all_idx[n] for n in var_df["feature"].unique() if n in all_idx}

    annotation_cols_raw, _ = identify_columns(var_df.columns)
    annotation_cols = [c for c in annotation_cols_raw if c not in blacklist]
    if blacklist:
        print(f"  Blacklisted {len(blacklist)} columns")

    available = [f"{factor}: {test}" for factor, test in var_df[["factor", "test"]].drop_duplicates().values]

    default_condition = available[0]

    return {
        "adata": adata, "var_df": var_df,
        "annotation_cols": annotation_cols, "blacklist": blacklist,
        "available_conditions": available, "default_condition": default_condition,
        "obs_columns": obs_columns, "layer_names": layer_names,
        "var_name_to_idx": var_name_to_idx,
    }


def build_volcano_df(var_df, annotation_cols, condition, config):

    factor = condition.split(': ')[0]
    test = condition.split(': ')[1]

    padj_t = config.get("padj_threshold", 0.05)
    lfc_t = config.get("log2fc_threshold", 1.0)

    df = var_df[(var_df["factor"] == factor) & (var_df["test"] == test)]

    keep_cols = [c for c in df.columns if c in (DE_COLUMNS + annotation_cols)]

    df = df[keep_cols]

    pos = df.loc[df["padj"] > 0, "padj"]
    if len(pos) > 0:
        df["padj"] = df["padj"].replace(0, pos.min() / 10)
    df["-log10(padj)"] = -np.log10(df["padj"].clip(lower=1e-300))
    sig = (df["padj"] < padj_t) & (np.abs(df["log2FoldChange"]) > lfc_t)
    df["regulation"] = "Not significant"
    df.loc[sig & (df["log2FoldChange"] > 0), "regulation"] = "Upregulated"
    df.loc[sig & (df["log2FoldChange"] < 0), "regulation"] = "Downregulated"
    for c in df.columns:
        if df[c].dtype.name == "category":
            df[c] = df[c].astype(str)

    return df


# ─── Volcano figure ───────────────────────────────────────────────────────
def build_volcano_figure(df, color_col, point_size, opacity, tooltip_cols, config,
                         show_all, max_nonsig, max_legend_classes, theme="light",
                         selected_reads=None, filter_revision=0,
                         palette_name="Dark24"):
    padj_t = config.get("padj_threshold", 0.05)
    lfc_t = config.get("log2fc_threshold", 1.0)
    x_min = config.get("x_min", -20)
    x_max = config.get("x_max", 20)

    full = df.copy()
    if selected_reads is not None:
        full = full[full["feature"].isin(selected_reads)]
    sig_mask = full["regulation"] != "Not significant"
    sig_df = full[sig_mask].copy()
    ns_df = full[~sig_mask].copy()
    if show_all and len(ns_df) > max_nonsig:
        ns_df = ns_df.sample(n=max_nonsig, random_state=42)

    plot_df = sig_df
    if len(plot_df) == 0 and (not show_all or len(ns_df) == 0):
        return go.Figure().update_layout(title="No data to display",
                                          paper_bgcolor="rgba(0,0,0,0)",
                                          plot_bgcolor="rgba(0,0,0,0)")

    actual_color = color_col if color_col in plot_df.columns else "regulation"
    is_num = (actual_color != "regulation"
              and pd.api.types.is_numeric_dtype(plot_df[actual_color]))
    if not is_num and max_legend_classes > 0 and actual_color != "regulation":
        top = plot_df[actual_color].value_counts().head(max_legend_classes).index
        plot_df[actual_color] = plot_df[actual_color].where(
            plot_df[actual_color].isin(top), other="Other")

    # Pre-format tooltip text per row: keys ljust, values right-aligned to
    # the longest value of that row. We render in monospace so the spaces
    # align in the SVG hover label.
    visible_tt_cols = [c for c in (tooltip_cols or []) if c in plot_df.columns]
    if visible_tt_cols:
        key_w = max(len(c) for c in visible_tt_cols) + 3
        val_series = [plot_df[c].astype(str).fillna("").str.slice(0, 60)
                      for c in visible_tt_cols]
        max_val_len = pd.concat([s.str.len() for s in val_series], axis=1).max(axis=1)
        lines_acc = None
        for c, vs in zip(visible_tt_cols, val_series):
            pad = (max_val_len - vs.str.len()).clip(lower=0)
            spaces = pad.map(lambda n: " " * int(n))
            line = c.ljust(key_w) + spaces + vs
            if lines_acc is None:
                lines_acc = line
            else:
                lines_acc = lines_acc + "<br>" + line
        plot_df = plot_df.assign(_tt=lines_acc.values)

    color_map = None
    if actual_color == "regulation":
        color_map = {
            "Upregulated":   "#f87171" if theme == "dark" else "#dc2626",
            "Downregulated": "#60a5fa" if theme == "dark" else "#2563eb",
            "Not significant": "#a8a29e" if theme == "light" else "#52525b",
        }

    custom_data_cols = ["feature"] + (["_tt"] if visible_tt_cols else [])
    # Disable auto-generated hovers (we set hovertemplate ourselves below)
    hover_data = {"log2FoldChange": False, "-log10(padj)": False}
    kwargs = dict(
        x="log2FoldChange", y="-log10(padj)", color=actual_color,
        hover_data=hover_data, custom_data=custom_data_cols,
        labels={"log2FoldChange": "log₂(fold change)",
                "-log10(padj)": "−log₁₀(adjusted p)"},
        render_mode="webgl",
    )
    if color_map:
        kwargs["color_discrete_map"] = color_map
    elif not is_num:
        # Discrete categorical (other than "regulation") → use the chosen
        # palette so >15-class annotations get distinguishable colors.
        kwargs["color_discrete_sequence"] = DISCRETE_PALETTES.get(
            palette_name, DISCRETE_PALETTES["Dark24"])
    if is_num:
        kwargs["color_continuous_scale"] = "Viridis"

    fig = px.scatter(plot_df, **kwargs)
    if visible_tt_cols:
        fig.update_traces(hovertemplate="%{customdata[1]}<extra></extra>")
    else:
        fig.update_traces(hovertemplate="<extra></extra>")

    if show_all and len(ns_df) > 0:
        fig.add_trace(go.Scattergl(
            x=ns_df["log2FoldChange"], y=ns_df["-log10(padj)"],
            mode="markers",
            marker=dict(size=point_size,
                        color="#52525b" if theme == "dark" else "#d4d4d8",
                        opacity=opacity * 0.5),
            customdata=ns_df[["feature"]].values,
            name="Not significant", showlegend=True, hoverinfo="skip",
        ))
        fig.data = (fig.data[-1],) + fig.data[:-1]

    fig.update_traces(marker=dict(size=point_size, opacity=opacity))

    threshold_color = "#52525b" if theme == "dark" else "#a8a29e"
    fig.add_hline(y=-np.log10(padj_t), line_dash="dash", line_color=threshold_color, opacity=0.6)
    fig.add_vline(x=lfc_t, line_dash="dash", line_color=threshold_color, opacity=0.6)
    fig.add_vline(x=-lfc_t, line_dash="dash", line_color=threshold_color, opacity=0.6)

    text_color = "#f5f5f4" if theme == "dark" else "#1c1917"
    grid_color = "rgba(255,255,255,0.04)" if theme == "dark" else "rgba(0,0,0,0.04)"
    axis_color = "#52525b" if theme == "dark" else "#a8a29e"

    fig.update_xaxes(range=[x_min, x_max], gridcolor=grid_color, zerolinecolor=grid_color,
                     linecolor=axis_color, tickcolor=axis_color,
                     tickfont=dict(color=text_color, size=11),
                     title_font=dict(color=text_color, size=12))
    fig.update_yaxes(gridcolor=grid_color, zerolinecolor=grid_color,
                     linecolor=axis_color, tickcolor=axis_color,
                     tickfont=dict(color=text_color, size=11),
                     title_font=dict(color=text_color, size=12))
    surface_color = "#141416" if theme == "dark" else "#ffffff"
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Inter, system-ui, sans-serif", color=text_color, size=12),
        margin=dict(l=56, r=20, t=24, b=44),
        legend=dict(
            title=dict(
                # Column name (mono, theme color) + a small italic hint
                # describing the legend interaction, placed on a second line.
                text=(
                    f"{actual_color}"
                    f"<br><span style=\"font-size:11px;color:{axis_color}\">"
                    f"<i>click ⟶ toggle · 2× ⟶ isolate</i></span>"
                ),
                font=dict(
                    family="JetBrains Mono, ui-monospace, SFMono-Regular, monospace",
                    color=text_color, size=13,
                ),
            ),
            uirevision=f"{actual_color}-{filter_revision}",
            font=dict(color=text_color, size=12),
            bgcolor="rgba(0,0,0,0)",
            # Push the legend down ~30px so it doesn't crowd the Plotly
            # modebar that floats at the top-right corner of the plot.
            yanchor="top",
            y=0.95,
        ),
        hovermode="closest", dragmode="pan",
        uirevision=f"volcano-{filter_revision}",
        hoverlabel=dict(
            bgcolor=surface_color,
            bordercolor=axis_color,
            font=dict(family="JetBrains Mono, ui-monospace, SFMono-Regular, monospace",
                      size=11, color=text_color),
            align="left",
            namelength=-1,
        ),
    )

    # Per-trace bordercolor = trace's marker color (so the surrounding box
    # shows which class the point belongs to, instead of the default opaque
    # background which is unreadable for light categories).
    for trace in fig.data:
        marker = getattr(trace, "marker", None)
        if marker is None:
            continue
        mc = getattr(marker, "color", None)
        if isinstance(mc, str) and mc.startswith(("#", "rgb")):
            trace.update(hoverlabel=dict(bordercolor=mc, bgcolor=surface_color,
                                          font=dict(color=text_color)))

    # Theme-aware ring colors:
    # - prefilter: muted gray that reads on either bg (lighter on dark, darker on light)
    # - highlight: a cyan accent distinct from the regulation red/blue palette
    prefilter_ring_color = "#a8a29e" if theme == "dark" else "#78716c"
    highlight_ring_color = "#22d3ee" if theme == "dark" else "#0e7490"

    # pre-filter rings (data[-2])
    fig.add_trace(go.Scattergl(
        x=[], y=[], mode="markers",
        marker=dict(size=point_size * 2, color="rgba(0,0,0,0)",
                    line=dict(color=prefilter_ring_color,
                               width=max(2, point_size * 0.4))),
        showlegend=False, hoverinfo="skip", name="_prefilter",
    ))
    # selected-feature highlight (data[-1])
    fig.add_trace(go.Scattergl(
        x=[], y=[], mode="markers",
        marker=dict(size=point_size * 2, color="rgba(0,0,0,0)",
                    line=dict(color=highlight_ring_color,
                               width=max(2, point_size * 0.6))),
        showlegend=False, hoverinfo="skip", name="_highlight",
    ))
    return fig


# ─── Detail panel rendering helpers ───────────────────────────────────────
def _render_sequence(seq, colors):
    """Render a nucleotide sequence with per-base background colors."""
    if seq is None or (isinstance(seq, float) and np.isnan(seq)):
        return html.Div("(no sequence)", className="cs-sequence",
                        style={"color": "var(--muted)", "fontStyle": "italic"})
    s = str(seq)
    if not s.strip() or s.lower() == "nan":
        return html.Div("(no sequence)", className="cs-sequence",
                        style={"color": "var(--muted)", "fontStyle": "italic"})
    fallback = colors.get("default", NUCLEOTIDE_COLORS_DEFAULT["default"])
    spans = []
    for nt in s:
        bg = colors.get(nt.upper(), fallback)
        spans.append(html.Span(nt, className="nt",
                                style={"backgroundColor": bg}))
    return html.Div(spans, className="cs-sequence")


def _resolve_nt_colors(config, theme):
    cfg = config or {}
    if (theme or "light") == "dark":
        return {**NUCLEOTIDE_COLORS_DEFAULT_DARK,
                **(cfg.get("nucleotide_colors_dark") or {})}
    return {**NUCLEOTIDE_COLORS_DEFAULT,
            **(cfg.get("nucleotide_colors") or {})}


def _render_detail_head(feat, var_df=None, config=None, theme="light"):
    if not feat:
        return [
            html.Div(className="cs-detail-head-row", children=[
                html.Div("Selected sequence", className="cs-eyebrow"),
            ]),
            html.Div("Click a point on the plot",
                     className="cs-detail-name",
                     style={"color": "var(--muted)", "fontStyle": "italic"}),
        ]
    seq_col = (config or {}).get("sequence_column", "seq")
    colors = _resolve_nt_colors(config, theme)
    seq = None
    if var_df is not None and seq_col in var_df.columns and feat in var_df["feature"]:
        seq = var_df.loc[var_df["feature"]==feat, seq_col].iloc[0]
    return [
        html.Div(className="cs-detail-head-row", children=[
            html.Div("Selected sequence", className="cs-eyebrow"),
            html.Button("Unselect", id="unselect-btn",
                        className="cs-action-btn",
                        title="Clear the current point selection"),
        ]),
        _render_sequence(seq, colors),
    ]


def _render_stat_tiles(feat, var_df, condition, config):
    if not feat or feat not in var_df["feature"]:
        return []
    factor = condition.split(": ")[0]
    test = condition.split(": ")[1]
    fd = var_df.loc[(var_df["feature"]==feat) & (var_df["factor"]==factor) & (var_df["test"]==test), :].iloc[0]
    lfc = fd.get("log2FoldChange", np.nan)
    padj = fd.get("padj", np.nan)
    bm = fd.get(f"baseMean", np.nan)

    try:
        lfc_v = float(lfc)
    except (TypeError, ValueError):
        lfc_v = np.nan
    try:
        padj_v = float(padj)
    except (TypeError, ValueError):
        padj_v = np.nan
    try:
        bm_v = float(bm)
    except (TypeError, ValueError):
        bm_v = np.nan

    padj_t = config.get("padj_threshold", 0.05)
    lfc_t = config.get("log2fc_threshold", 1.0)
    is_sig = (not np.isnan(padj_v) and padj_v < padj_t
              and not np.isnan(lfc_v) and abs(lfc_v) > lfc_t)
    if is_sig and lfc_v > 0:
        lfc_color = "var(--accent-up)"
    elif is_sig and lfc_v < 0:
        lfc_color = "var(--accent-dn)"
    else:
        lfc_color = "var(--text)"

    if not np.isnan(padj_v):
        neglogp = -np.log10(max(padj_v, 1e-300))
        padj_str = f"{padj_v:.2e}" if padj_v < 0.001 else f"{padj_v:.3f}"
        neglogp_str = f"{neglogp:.1f}"
    else:
        padj_str = "—"
        neglogp_str = "—"
    lfc_str = f"{lfc_v:+.2f}" if not np.isnan(lfc_v) else "—"
    bm_str = f"{bm_v:.0f}" if not np.isnan(bm_v) else "—"

    return [
        html.Div([
            html.Div("log₂ FC", className="cs-stat-label"),
            html.Div(lfc_str, className="cs-stat-value", style={"color": lfc_color}),
        ]),
        html.Div([
            html.Div("p adj", className="cs-stat-label"),
            html.Div(padj_str, className="cs-stat-value"),
        ]),
        html.Div([
            html.Div("−log₁₀(p)", className="cs-stat-label"),
            html.Div(neglogp_str, className="cs-stat-value"),
        ]),
        html.Div([
            html.Div("base mean", className="cs-stat-label"),
            html.Div(bm_str, className="cs-stat-value"),
        ]),
    ]


def _plot_title_children(condition, available_conditions):
    """Build the volcano plot's title children.
    """
    
    return [condition]

    


def _render_metadata_rows(feat, var_df, blacklist, config=None):
    if not feat or feat not in var_df["feature"]:
        return []
    fd = var_df[var_df["feature"]==feat].iloc[0]
    config = config or {}
    annotation_cols = [
        c for c in fd.index
        if c not in (list(blacklist) + DE_COLUMNS + ["feature",]) 
    ]
    order = [c for c in (config.get("metadata_column_order") or [])
             if c in annotation_cols]
    rest = sorted(c for c in annotation_cols if c not in order)
    rows = []
    for col in order + rest:
        rows.append(html.Div(className="cs-meta-row", children=[
            html.Span(str(col), className="k"),
            html.Span(str(fd[col]), className="v"),
        ]))
    # Read name (var_df index) goes last
    rows.append(html.Div(className="cs-meta-row", children=[
        html.Span("feature", className="k"),
        html.Span(str(feat), className="v mono"),
    ]))
    return rows


# ─── Layout (Clean Scientific) ────────────────────────────────────────────
def build_layout(data, config):
    annotation_cols = sorted(data["annotation_cols"])
    obs_columns = data["obs_columns"]
    layer_names = data["layer_names"]
    available = data["available_conditions"]
    default_condition = data["default_condition"]
    saved_views_init = config.get("saved_views") or []

    color_options = [{"label": "regulation", "value": "regulation"}] + [
        {"label": c, "value": c} for c in annotation_cols]
    layer_options = [{"label": "DESeq2 normalized counts", "value": "X"}] + [
        {"label": k, "value": k} for k in layer_names]
    obs_options = [{"label": c, "value": c} for c in obs_columns]
    configured_tt = config.get("default_tooltip_columns", [])
    default_tooltip = ([c for c in configured_tt if c in annotation_cols] if configured_tt
                       else annotation_cols[:5])
    quick_cols = [c for c in config.get("quick_color_columns", []) if c in annotation_cols]

    # Header
    project_pill_target = html.Button(
        id="project-pill-target", className="cs-project-pill",
        children=[
            html.Span(className="dot"),
            html.Span(default_condition, id="project-pill-condition",
                      style={"fontWeight": 500}),
            html.Span("·", className="sub"),
            html.Span(f"{len(data['var_df']["feature"].unique()):,} features", className="sub"),
            _svg_icon(ICON_CARET, 11),
        ],
    )
    condition_menu = dmc.Menu(
        children=[
            dmc.MenuTarget(project_pill_target),
            dmc.MenuDropdown([
                dmc.MenuLabel("Switch test"),
                dmc.TextInput(
                    id="project-pill-search",
                    placeholder="Search dataset…",
                    size="xs",
                    style={"padding": "4px 8px"},
                ),
                html.Div(
                    id="project-pill-items",
                    style={"maxHeight": "320px", "overflowY": "auto"},
                    children=[
                        dmc.MenuItem(c, id={"type": "condition-item", "value": c})
                        for c in available
                    ],
                ),
            ]),
        ],
        position="bottom-start",
        shadow="md",
        width=280,
        closeOnItemClick=True,
        closeOnClickOutside=True,
    )

    header = html.Div(
        className="cs-header",
        children=[
            html.Div("Volcano Explorer", className="cs-brand"),
            html.Div(className="cs-divider"),
            condition_menu,
            html.Button(
                id="header-saved-btn",
                className="cs-saved-trigger",
                children=[
                    _svg_icon(ICON_BOOKMARK, 13),
                    html.Span("Saved views"),
                ],
                title="Open the Saved views drawer",
            ),
            html.Div(id="recent-views-pills", className="cs-recent-pills"),
            html.Div(style={"flex": 1}),
            html.Button(
                id="cs-search-trigger", className="cs-search-trigger",
                children=[
                    _svg_icon(ICON_SEARCH, 13),
                    html.Span("Quick filter"),
                    html.Span("⌘K", className="kbd"),
                ],
            ),
            html.Button(id="cs-theme-toggle-btn", className="cs-icon-btn",
                        children=_svg_icon(ICON_SUN, 15),
                        title="Toggle light/dark"),
            html.Button(id="cs-help-btn", className="cs-icon-btn",
                        children=_svg_icon(ICON_HELP, 15),
                        title="Keyboard shortcuts"),
        ],
    )

    # ── Filter tab content (CNF builder + history) ────────────────────
    # Filters in the history are conjunctive (AND). Within a single filter
    # we let the user OR several subfilters together — a subfilter is a
    # (column, relation, value) tuple. The form below is the *current*
    # filter being composed. Submitted via Highlight / Filter.
    filter_tab = html.Div(id="left-filter-tab", children=[
        html.Div("new filter", className="cs-eyebrow"),
        html.Div(id="subfilters-container"),
        html.Div(style={"height": 8}),
        html.Div(style={"display": "flex", "gap": 6}, children=[
            html.Button("Highlight", id="prefilter-btn",
                        className="cs-action-btn", style={"flex": 1}),
            html.Button("Filter", id="update-selection-btn",
                        className="cs-btn-primary", style={"flex": 1}),
        ]),
        html.Div(style={"height": 6}),
        html.Div(style={"display": "flex", "gap": 6}, children=[
            html.Button("Unhighlight", id="undo-prefilter-btn",
                        className="cs-action-btn", style={"flex": 1}),
            html.Button("Reset", id="reset-selection-btn",
                        className="cs-action-btn", style={"flex": 1}),
        ]),
        html.Div(id="prefilter-loading-target", style={"display": "none"}),
        html.Div(style={"height": 8}),
        html.Div(id="selection-status"),
        html.Div(style={"height": 18}),
        html.Div("Filter history", className="cs-eyebrow"),
        html.Div(id="filter-history-table", children=_render_filter_table([])),
    ])

    # ── Display tab content (display options + tooltip cols) ──────────
    display_tab = html.Div(id="left-display-tab",
                            style={"display": "none"},
                            children=[
        html.Div("Display options", className="cs-eyebrow"),
        dmc.Stack([
            dmc.Checkbox(id="show-all-toggle",
                         label="Show non-significant points",
                         checked=config.get("show_all_points", False)),
            dmc.Text("Point size", size="xs", c="dimmed"),
            dmc.Slider(id="size-slider", min=1, max=15, step=1,
                       value=config.get("default_point_size", 5)),
            dmc.Text("Opacity", size="xs", c="dimmed"),
            dmc.Slider(id="opacity-slider", min=0.1, max=1.0, step=0.05,
                       value=config.get("default_opacity", 0.7)),
            dmc.NumberInput(id="max-legend-input",
                            label="Max legend classes",
                            value=config.get("max_legend_classes", 20),
                            min=-1, step=1, size="xs"),
            dmc.Select(
                id="palette-dropdown",
                label="Color palette",
                data=[
                    {"label": f"{name} ({len(colors)})", "value": name}
                    for name, colors in DISCRETE_PALETTES.items()
                ],
                value=config.get("default_palette", "Dark24"),
                size="xs",
                searchable=False,
            ),
        ], gap="xs"),
        html.Div(style={"height": 14}),
        html.Div("Tooltip columns", className="cs-eyebrow"),
        html.Div(style={"display": "flex", "gap": 6, "marginBottom": 8},
                 children=[
            html.Button("All", id="select-all-btn",
                        className="cs-action-btn", style={"flex": 1}),
            html.Button("None", id="deselect-all-btn",
                        className="cs-action-btn", style={"flex": 1}),
        ]),
        dmc.CheckboxGroup(
            id="tooltip-checklist", value=default_tooltip,
            children=dmc.Stack([
                dmc.Checkbox(label=c, value=c) for c in annotation_cols
            ], gap=4),
        ),
    ])

    controls = html.Div(
        id="cs-controls",
        className="cs-controls",
        children=[
            html.Div(className="cs-section", children=[
                html.Div("Color points by", className="cs-eyebrow"),
                dmc.Select(id="color-dropdown", data=color_options,
                           value=config.get("default_color_column", "regulation"),
                           searchable=True, size="xs"),
                html.Div(style={"height": 8}),
                html.Div(id="quick-chip-group", className="cs-chip-group", children=[
                    html.Button(c, id={"type": "quick-chip", "value": c},
                                className="cs-chip")
                    for c in quick_cols
                ]) if quick_cols else html.Div(),
            ]),
            html.Div(className="cs-section cs-section-tabs", children=[
                html.Div(className="cs-tab-bar", children=[
                    html.Button("Filter", id="left-tab-filter-btn",
                                className="cs-tab active"),
                    html.Button("Display", id="left-tab-display-btn",
                                className="cs-tab"),
                ]),
                html.Div(className="cs-tab-body", children=[
                    filter_tab,
                    display_tab,
                ]),
            ]),
        ],
    )

    # Center plot
    main = html.Div(
        className="cs-main",
        children=[
            html.Div(className="cs-plot-header", children=[
                html.Div([
                    html.H1(id="cs-plot-title",
                            children=_plot_title_children(
                                default_condition, available),
                            className="cs-plot-title"),
                    html.Div(id="cs-plot-meta", className="cs-plot-meta"),
                ]),
                html.Div(style={"display": "flex", "gap": 4}, children=[
                    html.Button("Export", id="export-btn",
                                className="cs-action-btn"),
                    html.Button("Reset", id="reset-view-btn",
                                className="cs-action-btn"),
                ]),
            ]),
            html.Div(className="cs-plot-area", children=[
                dcc.Graph(id="volcano-plot",
                          config={"scrollZoom": True, "displaylogo": False,
                                  "showTips": False},
                          style={"height": "100%"}),
                # Filter-visible button — absolutely positioned over the plot
                # and snapped to the bottom of the Plotly legend by a
                # clientside observer (see legend-btn-positioner callback).
                # Commits the currently-visible legend classes as a filter.
                html.Button(
                    "Filter visible classes",
                    id="filter-legend-btn",
                    className="cs-action-btn cs-legend-action",
                    title=("Commit the currently-visible legend classes as "
                            "a filter on the volcano plot"),
                    style={"display": "none"},
                ),
                html.Div(id="legend-btn-positioner",
                         style={"display": "none"}),
            ]),
            html.Div(className="cs-stats-strip", id="cs-stats-strip"),
        ],
    )

    # Right detail panel — both tab contents always present (stable IDs)
    boxplot_tab_content = html.Div(
        id="boxplot-tab-content",
        style={"display": "none"},
        children=[
            dmc.Stack([
                dmc.Select(id="violin-obs-dropdown",
                           label="Group by (obs column)",
                           data=obs_options,
                           value=obs_columns[0] if obs_columns else None,
                           searchable=True, size="xs"),
                dmc.Select(id="violin-layer-dropdown",
                           label="Expression layer",
                           data=layer_options, value="X", size="xs"),
                dmc.Checkbox(id="boxplot-log-toggle",
                             label="Log transform (log1p)", checked=False),
                dmc.Select(id="boxplot-sort-dropdown",
                           label="Order groups by",
                           data=[
                               {"value": "median_desc", "label": "Median (desc)"},
                               {"value": "median_asc", "label": "Median (asc)"},
                               {"value": "mean_desc", "label": "Mean (desc)"},
                               {"value": "mean_asc", "label": "Mean (asc)"},
                               {"value": "alpha_asc", "label": "Alphabetic (asc)"},
                               {"value": "alpha_desc", "label": "Alphabetic (desc)"},
                               {"value": "max_desc", "label": "Max (desc)"},
                               {"value": "max_asc", "label": "Max (asc)"},
                               {"value": "min_desc", "label": "Min (desc)"},
                               {"value": "min_asc", "label": "Min (asc)"},
                           ],
                           value="median_desc", size="xs"),
                html.Button("Generate Box Plot", id="violin-button",
                            className="cs-btn-primary"),
                html.Div(id="violin-status"),
                dcc.Loading(
                    dcc.Graph(id="violin-plot",
                              config={"displayModeBar": False}),
                    type="default",
                ),
            ], gap="xs"),
        ],
    )

    detail = html.Div(
        id="cs-detail",
        className="cs-detail cs-no-feature",
        children=[
            html.Div(id="cs-detail-head", className="cs-detail-head",
                     children=_render_detail_head(None)),
            html.Div(id="cs-detail-stats", className="cs-detail-stats", children=[]),
            html.Div(className="cs-tab-bar", children=[
                html.Button("Metadata", id="tab-metadata-btn",
                            className="cs-tab active"),
                html.Button("Box plot", id="tab-boxplot-btn",
                            className="cs-tab"),
            ]),
            html.Div(className="cs-detail-body", children=[
                html.Div(id="metadata-tab-content",
                         children=[html.Div(id="metadata-rows", children=[])]),
                boxplot_tab_content,
            ]),
        ],
    )

    # Quick-filter modal (⌘K) — applies as pre-filter "any contains <value>"
    search_modal = dmc.Modal(
        id="search-modal", opened=False,
        title="Quick filter — any contains",
        size="md",
        children=[
            dmc.Text(
                "Type a value; pressing Enter or clicking Highlight is "
                "equivalent to setting the left panel's highlight to "
                "'any contains <value>'.",
                size="xs", c="dimmed",
            ),
            html.Div(style={"height": 8}),
            dmc.TextInput(id="search-input",
                          placeholder="e.g. EGFR, LINE, chr4 …",
                          size="sm"),
            html.Div(style={"height": 10}),
            html.Div(style={"display": "flex", "gap": 6,
                            "justifyContent": "flex-end"}, children=[
                html.Button("Cancel", id="search-cancel-btn",
                            className="cs-action-btn"),
                html.Button("Highlight", id="search-apply-btn",
                            className="cs-btn-primary",
                            style={"width": "auto", "padding": "7px 14px"}),
            ]),
        ],
    )

    # Help modal — keyboard shortcuts
    help_modal = dmc.Modal(
        id="help-modal", opened=False,
        title="Keyboard shortcuts",
        size="md",
        children=[
            html.Div(className="cs-shortcut-row", children=[
                html.Span("Open quick filter (any contains)", className="k"),
                html.Span("⌘K · Ctrl+K", className="kbd"),
            ]),
            html.Div(className="cs-shortcut-row", children=[
                html.Span("Submit highlight from input", className="k"),
                html.Span("Enter", className="kbd"),
            ]),
            html.Div(className="cs-shortcut-row", children=[
                html.Span("Close modal / drawer", className="k"),
                html.Span("Esc", className="kbd"),
            ]),
            html.Div(className="cs-shortcut-row", children=[
                html.Span("Pan / zoom plot", className="k"),
                html.Span("Drag · Scroll", className="kbd"),
            ]),
            html.Div(className="cs-shortcut-row", children=[
                html.Span("Toggle a legend class on the plot", className="k"),
                html.Span("Click", className="kbd"),
            ]),
            html.Div(className="cs-shortcut-row", children=[
                html.Span("Isolate a single legend class", className="k"),
                html.Span("Double-click", className="kbd"),
            ]),
        ],
    )

    # Saved views drawer
    saved_drawer = dmc.Drawer(
        id="saved-drawer", opened=False,
        title="Saved views", position="left", size=380,
        children=[
            dmc.Text(
                "Saved views capture the active condition, color, filter history "
                "and display options. They persist in the YAML config file.",
                size="xs", c="dimmed",
            ),
            html.Div(style={"height": 12}),
            dmc.TextInput(id="save-view-name", placeholder="View name…", size="xs"),
            html.Div(style={"height": 6}),
            html.Button("+ Save current view", id="save-view-btn",
                        className="cs-btn-primary"),
            html.Div(id="save-view-status",
                     style={"fontSize": 11, "color": "var(--muted)",
                            "marginTop": 6, "minHeight": 14}),
            html.Div(id="saved-views-list", className="cs-saved-list"),
        ],
    )

    return dmc.MantineProvider(
        html.Div([
            # Stores
            dcc.Store(id="theme-store", data="light", storage_type="local"),
            dcc.Store(id="selected-condition", data=default_condition),
            dcc.Store(id="selected-reads", data=None),
            dcc.Store(id="filter-history", data=[]),
            dcc.Store(id="prefilter-reads", data=None),
            dcc.Store(id="active-detail-tab", data="metadata"),
            dcc.Store(id="active-controls-tab", data="filter"),
            dcc.Store(id="active-feature", data=None),
            # current-subfilters drives the form rendered in the "new filter"
            # section. Initialised with one empty subfilter so the form is
            # never blank on first load.
            dcc.Store(
                id="current-subfilters",
                data=[{"column": "any", "relation": "contains", "value": ""}],
            ),
            dcc.Store(id="saved-views", data=saved_views_init),
            # recent-views persists across page reloads (per browser);
            # render_recent_pills falls back to the tail of saved-views when
            # this store is empty (e.g. fresh browser, cleared cache, or a
            # different project whose names don't match any persisted entry).
            dcc.Store(id="recent-views", data=[], storage_type="local"),
            # Hidden div for clientside keyboard init
            html.Div(id="kbd-init", style={"display": "none"}),
            # Download for export
            dcc.Download(id="export-download"),
            # Top
            header,
            html.Div(className="cs-body", children=[
                controls,
                html.Div(id="cs-resize-left", className="cs-resize-handle"),
                main,
                html.Div(id="cs-resize-right", className="cs-resize-handle"),
                detail,
            ]),
            search_modal,
            help_modal,
            saved_drawer,
            # Hidden div used by the resize-init clientside callback
            html.Div(id="resize-init", style={"display": "none"}),
        ]),
    )


# ─── App factory ──────────────────────────────────────────────────────────
def create_app(config, config_path=None):
    data = load_data(config)
    var_df = data["var_df"]
    adata = data["adata"]
    annotation_cols = data["annotation_cols"]
    blacklist = data["blacklist"]
    var_name_to_idx = data["var_name_to_idx"]
    available_conditions = data["available_conditions"]
    quick_cols = [c for c in config.get("quick_color_columns", [])
                  if c in annotation_cols]

    atexit.register(lambda: adata.file.close())

    app = Dash(__name__, assets_folder="assets", suppress_callback_exceptions=True)
    app.title = "Volcano Plot Explorer"
    app.layout = build_layout(data, config)

    # ── Clientside: theme toggle (flip data-theme on <html>) ────────────
    app.clientside_callback(
        """
        function(n, current) {
            if (!n) {
                if (current === 'dark')
                    document.documentElement.setAttribute('data-theme', 'dark');
                else
                    document.documentElement.removeAttribute('data-theme');
                return current || 'light';
            }
            const next = current === 'dark' ? 'light' : 'dark';
            if (next === 'dark')
                document.documentElement.setAttribute('data-theme', 'dark');
            else
                document.documentElement.removeAttribute('data-theme');
            return next;
        }
        """,
        Output("theme-store", "data"),
        Input("cs-theme-toggle-btn", "n_clicks"),
        State("theme-store", "data"),
    )

    # ── Clientside: drag-to-resize for left/right sidebars ─────────────
    app.clientside_callback(
        """
        function(_) {
            if (window._volcanoResize) return '';
            window._volcanoResize = true;
            const init = (handleId, sideId, fromRight) => {
                const handle = document.getElementById(handleId);
                const side = document.getElementById(sideId);
                if (!handle || !side) return;
                handle.addEventListener('mousedown', function(e) {
                    e.preventDefault();
                    handle.classList.add('dragging');
                    document.body.classList.add('cs-dragging');
                    const startX = e.clientX;
                    const startW = side.getBoundingClientRect().width;
                    const total = side.parentElement.getBoundingClientRect().width;
                    const minW = total * 0.10, maxW = total * 0.50;
                    const move = ev => {
                        let dx = ev.clientX - startX;
                        if (fromRight) dx = -dx;
                        const w = Math.max(minW, Math.min(maxW, startW + dx));
                        side.style.width = w + 'px';
                        side.style.flexBasis = w + 'px';
                        // Plotly listens for window resize to relayout the figure
                        // when its container changes width.
                        window.dispatchEvent(new Event('resize'));
                    };
                    const up = () => {
                        document.removeEventListener('mousemove', move);
                        document.removeEventListener('mouseup', up);
                        handle.classList.remove('dragging');
                        document.body.classList.remove('cs-dragging');
                        window.dispatchEvent(new Event('resize'));
                    };
                    document.addEventListener('mousemove', move);
                    document.addEventListener('mouseup', up);
                });
            };
            init('cs-resize-left', 'cs-controls', false);
            init('cs-resize-right', 'cs-detail', true);
            return '';
        }
        """,
        Output("resize-init", "children"),
        Input("resize-init", "id"),
    )

    # ── Clientside: ⌘K / Ctrl+K opens search ───────────────────────────
    app.clientside_callback(
        """
        function(_) {
            if (window._volcanoKbd) return '';
            window._volcanoKbd = true;
            document.addEventListener('keydown', function(e) {
                const k = (e.key || '').toLowerCase();
                if ((e.metaKey || e.ctrlKey) && k === 'k') {
                    e.preventDefault();
                    const btn = document.getElementById('cs-search-trigger');
                    if (btn) btn.click();
                }
            });
            return '';
        }
        """,
        Output("kbd-init", "children"),
        Input("kbd-init", "id"),
    )

    # ── Clientside: snap "Filter visible classes" button under the legend ─
    # The Plotly legend's vertical extent depends on the number of classes
    # currently displayed, so we observe its rendered SVG bounding box and
    # position the button at its bottom-right edge. Re-runs whenever the
    # figure is rebuilt and on window resize.
    app.clientside_callback(
        """
        function(_fig) {
            const place = () => {
                const btn = document.getElementById('filter-legend-btn');
                const area = document.querySelector('.cs-plot-area');
                if (!btn || !area) return;
                // Plotly renders the legend as <g class="legend"> inside
                // .infolayer. Take the first one we find.
                const legend = area.querySelector('g.legend');
                if (!legend) {
                    btn.style.display = 'none';
                    return;
                }
                const lb = legend.getBoundingClientRect();
                const ab = area.getBoundingClientRect();
                if (lb.width === 0 || lb.height === 0) {
                    btn.style.display = 'none';
                    return;
                }
                btn.style.top   = (lb.bottom - ab.top + 8) + 'px';
                btn.style.right = (ab.right - lb.right) + 'px';
                btn.style.display = 'inline-block';
            };
            // Wait one paint frame so Plotly has finished rendering.
            requestAnimationFrame(() => requestAnimationFrame(place));
            if (!window._volcanoLegendBtnObserver) {
                window._volcanoLegendBtnObserver = true;
                window.addEventListener('resize', place);
            }
            return '';
        }
        """,
        Output("legend-btn-positioner", "children"),
        Input("volcano-plot", "figure"),
    )

    # ── Tooltip select-all/none ─────────────────────────────────────────
    @app.callback(
        Output("tooltip-checklist", "value"),
        [Input("select-all-btn", "n_clicks"),
         Input("deselect-all-btn", "n_clicks")],
        prevent_initial_call=True,
    )
    def toggle_tt(_s, _d):
        if callback_context.triggered_id == "select-all-btn":
            return annotation_cols
        return []

    # ── Quick chip → color dropdown ────────────────────────────────────
    @app.callback(
        Output("color-dropdown", "value", allow_duplicate=True),
        Input({"type": "quick-chip", "value": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def chip_pick(_clicks):
        t = callback_context.triggered_id
        if not t or not isinstance(t, dict):
            return dash.no_update
        if not any(_clicks or []):
            return dash.no_update
        return t["value"]

    # ── Quick chip active class sync ───────────────────────────────────
    if quick_cols:
        @app.callback(
            Output({"type": "quick-chip", "value": ALL}, "className"),
            Input("color-dropdown", "value"),
            State({"type": "quick-chip", "value": ALL}, "id"),
        )
        def sync_chip_active(color_val, ids):
            return ["cs-chip active" if i["value"] == color_val else "cs-chip"
                    for i in ids]

    # ── Condition switching from header menu ───────────────────────────
    @app.callback(
        [Output("selected-condition", "data"),
         Output("project-pill-condition", "children"),
         Output("cs-plot-title", "children")],
        Input({"type": "condition-item", "value": ALL}, "n_clicks"),
        State("selected-condition", "data"),
        prevent_initial_call=True,
    )
    def switch_condition(_clicks, current):
        t = callback_context.triggered_id
        if not t or not isinstance(t, dict) or not any(_clicks or []):
            return dash.no_update, dash.no_update, dash.no_update
        new_cond = t["value"]
        if new_cond == current:
            return dash.no_update, dash.no_update, dash.no_update
        return (
            new_cond, new_cond,
            _plot_title_children(new_cond, available_conditions),
        )

    # ── Project pill search filter ─────────────────────────────────────
    @app.callback(
        Output("project-pill-items", "children"),
        Input("project-pill-search", "value"),
    )
    def filter_project_pill(search):
        s = (search or "").strip().lower()
        items = ([c for c in available_conditions if s in c.lower()] if s
                 else list(available_conditions))
        if not items:
            return [dmc.MenuLabel("No match")]
        return [dmc.MenuItem(c, id={"type": "condition-item", "value": c})
                for c in items]

    # ── Left sidebar tab switching (Filter / Display) ──────────────────
    @app.callback(
        [Output("active-controls-tab", "data"),
         Output("left-filter-tab", "style"),
         Output("left-display-tab", "style"),
         Output("left-tab-filter-btn", "className"),
         Output("left-tab-display-btn", "className")],
        [Input("left-tab-filter-btn", "n_clicks"),
         Input("left-tab-display-btn", "n_clicks")],
        prevent_initial_call=True,
    )
    def switch_controls_tab(_f, _d):
        t = callback_context.triggered_id
        if t == "left-tab-display-btn":
            return ("display",
                    {"display": "none"}, {"display": "block"},
                    "cs-tab", "cs-tab active")
        return ("filter",
                {"display": "block"}, {"display": "none"},
                "cs-tab active", "cs-tab")

    # ── Render subfilter form from the current-subfilters store ─────────
    @app.callback(
        Output("subfilters-container", "children"),
        Input("current-subfilters", "data"),
    )
    def render_subfilters_cb(subfilters):
        return _render_subfilters(subfilters, annotation_cols)

    # ── Add / remove subfilter groups (CNF builder) ─────────────────────
    # Captures the current values from the form before mutating the store
    # so the user doesn't lose what they had already typed.
    @app.callback(
        Output("current-subfilters", "data"),
        [Input({"type": "sf-add", "i": ALL}, "n_clicks"),
         Input({"type": "sf-del", "i": ALL}, "n_clicks")],
        [State({"type": "sf-col", "i": ALL}, "value"),
         State({"type": "sf-rel", "i": ALL}, "value"),
         State({"type": "sf-val", "i": ALL}, "value")],
        prevent_initial_call=True,
    )
    def update_subfilters_store(add_clicks, del_clicks, cols, rels, vals):
        triggered = callback_context.triggered_id
        if not isinstance(triggered, dict):
            return dash.no_update
        # Snapshot current values from the rendered form
        current = _build_subfilters(cols, rels, vals)
        t = triggered.get("type")
        if t == "sf-add":
            if not any(add_clicks or []):
                return dash.no_update
            current.append({"column": "any",
                             "relation": "contains",
                             "value": ""})
            return current
        if t == "sf-del":
            if not any(del_clicks or []):
                return dash.no_update
            idx = triggered.get("i", -1)
            if 0 <= idx < len(current):
                current.pop(idx)
            if not current:  # never let it go empty
                current = [{"column": "any",
                            "relation": "contains",
                            "value": ""}]
            return current
        return dash.no_update

    # ── Highlight (rings on plot) / Unhighlight ─────────────────────────
    @app.callback(
        [Output("prefilter-reads", "data"),
         Output("volcano-plot", "figure", allow_duplicate=True),
         Output("prefilter-loading-target", "children")],
        [Input("prefilter-btn", "n_clicks"),
         Input("undo-prefilter-btn", "n_clicks")],
        [State({"type": "sf-col", "i": ALL}, "value"),
         State({"type": "sf-rel", "i": ALL}, "value"),
         State({"type": "sf-val", "i": ALL}, "value"),
         State("selected-condition", "data"),
         State("selected-reads", "data"),
         State("show-all-toggle", "checked"),
         State("size-slider", "value")],
        prevent_initial_call=True,
    )
    def apply_prefilter(_a, _b, cols, rels, vals,
                        condition, selected_reads, show_all, point_size):
        triggered = callback_context.triggered_id
        if triggered == "undo-prefilter-btn":
            p = Patch()
            p["data"][-2]["x"] = []
            p["data"][-2]["y"] = []
            return None, p, ""
        sfs = _build_subfilters(cols, rels, vals)
        if not _has_nonempty_subfilter(sfs):
            return dash.no_update, dash.no_update, ""
        condition = condition or data["default_condition"]
        df = _filtered_volcano_df(condition, selected_reads, bool(show_all),
                                    var_df, annotation_cols, config)
        matching = _compute_subfilter_matches(df, sfs, annotation_cols)
        match_df = df[df["feature"].isin(matching)]
        p = Patch()
        p["data"][-2]["x"] = match_df["log2FoldChange"].tolist()
        p["data"][-2]["y"] = match_df["-log10(padj)"].tolist()
        p["data"][-2]["marker"]["size"] = (point_size or 5) * 2
        return matching, p, ""

    # ── Quick-filter (⌘K) modal: replaces the form with one subfilter ──
    @app.callback(
        [Output("current-subfilters", "data", allow_duplicate=True),
         Output("search-modal", "opened", allow_duplicate=True),
         Output("volcano-plot", "figure", allow_duplicate=True),
         Output("prefilter-reads", "data", allow_duplicate=True)],
        [Input("search-apply-btn", "n_clicks"),
         Input("search-input", "n_submit")],
        [State("search-input", "value"),
         State("selected-condition", "data"),
         State("selected-reads", "data"),
         State("show-all-toggle", "checked"),
         State("size-slider", "value")],
        prevent_initial_call=True,
    )
    def quick_filter_apply(_a, _s, search_value, condition, selected_reads,
                            show_all, point_size):
        val = (search_value or "").strip()
        if not val:
            return (dash.no_update,) * 4
        sfs = [{"column": "any", "relation": "contains", "value": val}]
        condition = condition or data["default_condition"]
        df = _filtered_volcano_df(condition, selected_reads, bool(show_all),
                                    var_df, annotation_cols, config)
        matching = _compute_subfilter_matches(df, sfs, annotation_cols)
        match_df = df[df["feature"].isin(matching)]
        p = Patch()
        p["data"][-2]["x"] = match_df["log2FoldChange"].tolist()
        p["data"][-2]["y"] = match_df["-log10(padj)"].tolist()
        p["data"][-2]["marker"]["size"] = (point_size or 5) * 2
        return sfs, False, p, matching

    # ── Filter (commit current subfilters to history) / Reset ──────────
    # The Filter button always recomputes from the form's subfilters
    # (joined by OR within this filter; AND between filters in history).
    @app.callback(
        [Output("selected-reads", "data"),
         Output("filter-history", "data"),
         Output("selection-status", "children"),
         Output("filter-history-table", "children"),
         Output("prefilter-reads", "data", allow_duplicate=True),
         Output("current-subfilters", "data", allow_duplicate=True)],
        [Input("update-selection-btn", "n_clicks"),
         Input("reset-selection-btn", "n_clicks")],
        [State("filter-history", "data"),
         State({"type": "sf-col", "i": ALL}, "value"),
         State({"type": "sf-rel", "i": ALL}, "value"),
         State({"type": "sf-val", "i": ALL}, "value"),
         State("selected-condition", "data"),
         State("selected-reads", "data"),
         State("show-all-toggle", "checked")],
        prevent_initial_call=True,
    )
    def update_selection(_u, _r, history, cols, rels, vals,
                         condition, selected_reads_state, show_all):
        if callback_context.triggered_id == "reset-selection-btn":
            return (None, [],
                    dmc.Text("Selection cleared.", size="xs", c="dimmed"),
                    _render_filter_table([]), None,
                    [{"column": "any", "relation": "contains", "value": ""}])
        sfs = _build_subfilters(cols, rels, vals)
        # Drop empty subfilters from the committed row (kept in form for UX)
        committed = [sf for sf in sfs if (sf.get("value") or "").strip()]
        if not committed:
            return (dash.no_update, dash.no_update,
                    dmc.Text("Enter at least one filter value.",
                             size="xs", c="orange"),
                    dash.no_update, dash.no_update, dash.no_update)
        condition = condition or data["default_condition"]
        matches = _compute_filter_matches(
            committed, condition, selected_reads_state, bool(show_all),
            var_df, annotation_cols, config,
        )
        n = len(matches)
        params = {"type": "text", "subfilters": committed}
        desc = _format_filter_desc(params)
        history = (history or []) + [
            {"desc": desc, "n_points": n, "reads": matches, "params": params}
        ]
        # Reset the form to one empty subfilter — the row was just
        # committed, the user starts the next filter from scratch.
        fresh_form = [{"column": "any", "relation": "contains", "value": ""}]
        return (matches, history,
                dmc.Text(f"{n:,} features kept.", size="xs", c="dimmed"),
                _render_filter_table(history), None, fresh_form)

    # ── Filter visible classes (legend-based, dedicated button) ────────
    @app.callback(
        [Output("selected-reads", "data", allow_duplicate=True),
         Output("filter-history", "data", allow_duplicate=True),
         Output("selection-status", "children", allow_duplicate=True),
         Output("filter-history-table", "children", allow_duplicate=True)],
        Input("filter-legend-btn", "n_clicks"),
        [State("volcano-plot", "figure"),
         State("color-dropdown", "value"),
         State("max-legend-input", "value"),
         State("filter-history", "data")],
        prevent_initial_call=True,
    )
    def filter_legend_visible(_n, fig, color_col, max_legend, history):
        if not _n or not fig or "data" not in fig:
            return (dash.no_update,) * 4
        visible_reads = []
        visible_names = []
        hidden_names = []
        for tr in fig["data"]:
            if "customdata" not in tr:
                continue
            vis = tr.get("visible", True)
            name = tr.get("name", "")
            if vis == "legendonly" or vis is False:
                hidden_names.append(name)
            else:
                visible_names.append(name)
                for pt in tr["customdata"]:
                    visible_reads.append(pt[0])
        n = len(visible_reads)
        col = color_col or "?"
        half = max(1, (max_legend if max_legend and max_legend > 0 else 20) // 2)
        if not hidden_names:
            desc = f"{col}: all classes"
        elif len(hidden_names) == 1:
            desc = f"{col} != {hidden_names[0]}"
        elif len(visible_names) == 1:
            desc = f"{col} = {visible_names[0]}"
        elif len(visible_names) <= half:
            desc = f"{col} IN ({', '.join(visible_names)})"
        else:
            desc = f"{col} NOT IN ({', '.join(hidden_names)})"
        params = {
            "type": "legend",
            "color_column": col,
            "visible_classes": visible_names,
            "max_legend": max_legend if max_legend is not None else -1,
        }
        history = (history or []) + [
            {"desc": desc, "n_points": n, "reads": visible_reads,
             "params": params}
        ]
        return (visible_reads, history,
                dmc.Text(f"{n:,} features kept (legend).",
                         size="xs", c="dimmed"),
                _render_filter_table(history))

    # ── Delete a single row from filter history (trash icon) ───────────
    # The remaining rows are *re-applied* in order from the initial state
    # (no filter), so each row's reads/n_points are recomputed against the
    # current condition + show_all toggle. Avoids the inconsistencies you
    # would see if we just spliced and kept the cached snapshots.
    @app.callback(
        [Output("selected-reads", "data", allow_duplicate=True),
         Output("filter-history", "data", allow_duplicate=True),
         Output("selection-status", "children", allow_duplicate=True),
         Output("filter-history-table", "children", allow_duplicate=True)],
        Input({"type": "filter-row-delete", "index": ALL}, "n_clicks"),
        [State("filter-history", "data"),
         State("selected-condition", "data"),
         State("show-all-toggle", "checked")],
        prevent_initial_call=True,
    )
    def delete_filter_row(n_clicks, history, condition, show_all):
        if not any(n_clicks or []):
            return (dash.no_update,) * 4
        idx = callback_context.triggered_id["index"]
        history = list(history or [])
        if idx >= len(history):
            return (dash.no_update,) * 4
        history.pop(idx)
        recomputed = _recompute_history_chain(
            history,
            condition or data["default_condition"],
            bool(show_all),
            var_df, annotation_cols, config,
        )
        new_selected = recomputed[-1]["reads"] if recomputed else None
        if recomputed:
            n = recomputed[-1]["n_points"]
            status = dmc.Text(f"{n:,} features kept (recomputed).",
                               size="xs", c="dimmed")
        else:
            status = dmc.Text("All filters removed.", size="xs", c="dimmed")
        return (new_selected, recomputed, status,
                _render_filter_table(recomputed))

    # ── Unselect the current point (clears highlight ring + active feat) ─
    @app.callback(
        [Output("active-feature", "data", allow_duplicate=True),
         Output("volcano-plot", "figure", allow_duplicate=True)],
        Input("unselect-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def unselect_feature(_n):
        if not _n:
            return dash.no_update, dash.no_update
        p = Patch()
        p["data"][-1]["x"] = []
        p["data"][-1]["y"] = []
        return None, p

    # ── Restore filter from history click ──────────────────────────────
    @app.callback(
        [Output("selected-reads", "data", allow_duplicate=True),
         Output("filter-history", "data", allow_duplicate=True),
         Output("selection-status", "children", allow_duplicate=True),
         Output("filter-history-table", "children", allow_duplicate=True)],
        Input({"type": "filter-row", "index": ALL}, "n_clicks"),
        State("filter-history", "data"),
        prevent_initial_call=True,
    )
    def restore_filter(n_clicks, history):
        if not any(n_clicks or []):
            return (dash.no_update, dash.no_update,
                    dash.no_update, dash.no_update)
        idx = callback_context.triggered_id["index"]
        history = (history or [])[:idx + 1]
        n = history[idx]["n_points"]
        return (history[idx]["reads"], history,
                dmc.Text(f"Restored: {n:,} features.", size="xs", c="dimmed"),
                _render_filter_table(history))

    # ── Reset view (clear filters + zoom) ───────────────────────────────
    @app.callback(
        [Output("selected-reads", "data", allow_duplicate=True),
         Output("filter-history", "data", allow_duplicate=True),
         Output("prefilter-reads", "data", allow_duplicate=True),
         Output("selection-status", "children", allow_duplicate=True),
         Output("filter-history-table", "children", allow_duplicate=True),
         Output("volcano-plot", "figure", allow_duplicate=True)],
        Input("reset-view-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_view(n):
        if not n:
            return (dash.no_update,) * 6
        x_min = config.get("x_min", -20)
        x_max = config.get("x_max", 20)
        p = Patch()
        p["layout"]["xaxis"]["range"] = [x_min, x_max]
        p["layout"]["yaxis"]["autorange"] = True
        return (None, [], None,
                dmc.Text("View reset.", size="xs", c="dimmed"),
                _render_filter_table([]),
                p)

    # ── Volcano figure update ───────────────────────────────────────────
    @app.callback(
        [Output("volcano-plot", "figure"),
         Output("cs-plot-meta", "children"),
         Output("cs-stats-strip", "children")],
        [Input("color-dropdown", "value"), Input("size-slider", "value"),
         Input("opacity-slider", "value"), Input("tooltip-checklist", "value"),
         Input("show-all-toggle", "checked"), Input("max-legend-input", "value"),
         Input("selected-reads", "data"), Input("theme-store", "data"),
         Input("selected-condition", "data"),
         Input("palette-dropdown", "value")],
        State("filter-history", "data"),
    )
    def update_volcano(color_col, ps, op, tt, sa, ml, sel, theme, condition,
                       palette, history):
        condition = condition or data["default_condition"]
        df = build_volcano_df(var_df, annotation_cols, condition, config)
        max_ns = config.get("max_nonsignificant_points", 10000)
        ml = ml if ml is not None else -1
        show_all = bool(sa)
        fig = build_volcano_figure(
            df, color_col, ps or 5, op or 0.7, tt or [], config,
            show_all, max_ns, ml, theme=theme or "light",
            selected_reads=sel,
            filter_revision=len(history) if history else 0,
            palette_name=palette or "Dark24",
        )
        # Stats reflect what is *displayed* on the plot:
        # apply selected_reads filter, then the show_all toggle.
        visible = df.copy()
        if sel is not None:
            visible = visible[visible["feature"].isin(sel)]
        n_up = int((visible["regulation"] == "Upregulated").sum())
        n_dn = int((visible["regulation"] == "Downregulated").sum())
        if show_all:
            n_ns_visible = int((visible["regulation"] == "Not significant").sum())
            n_ns_displayed = min(n_ns_visible, max_ns)
        else:
            n_ns_displayed = 0
        meta = [
            "p", html.Sub("adj"), f" < {config.get('padj_threshold', 0.05)}",
            f" · |log₂FC| > {config.get('log2fc_threshold', 1.0)}",
            f" · {n_up + n_dn:,} significant of {len(df):,}",
        ]
        stats = [
            html.Span("Stats", className="cs-stats-eyebrow"),
            *[html.Div(className="cs-stats-item", children=[
                html.Span(label, className="label"),
                html.Span(f"{n:,}", className="count"),
            ]) for label, n in [
                ("Upregulated", n_up),
                ("Downregulated", n_dn),
                ("Not significant", n_ns_displayed),
            ]],
        ]
        return fig, meta, stats

    # ── Active feature: from plot click ─────────────────────────────────
    @app.callback(
        Output("active-feature", "data"),
        Input("volcano-plot", "clickData"),
        State("active-feature", "data"),
        prevent_initial_call=True,
    )
    def set_active_feature(click_data, current):
        if click_data and "points" in click_data and click_data["points"]:
            return click_data["points"][0]["customdata"][0]
        return current

    # ── Detail panel (head, stats, metadata rows) ──────────────────────
    @app.callback(
        [Output("cs-detail", "className"),
         Output("cs-detail-head", "children"),
         Output("cs-detail-stats", "children"),
         Output("metadata-rows", "children")],
        [Input("active-feature", "data"),
         Input("selected-condition", "data"),
         Input("theme-store", "data")],
    )
    def render_detail(feat, condition, theme):
        condition = condition or data["default_condition"]
        theme = theme or "light"
        if not feat:
            return ("cs-detail cs-no-feature",
                    _render_detail_head(None, var_df, config, theme), [], [])
        return (
            "cs-detail",
            _render_detail_head(feat, var_df, config, theme),
            _render_stat_tiles(feat, var_df, condition, config),
            _render_metadata_rows(feat, var_df, blacklist, config),
        )

    # ── Highlight ring on click / search ───────────────────────────────
    @app.callback(
        Output("volcano-plot", "figure", allow_duplicate=True),
        [Input("active-feature", "data")],
        [State("size-slider", "value"),
         State("selected-condition", "data"),
         State("volcano-plot", "figure")],
        prevent_initial_call=True,
    )
    def highlight_active_feature(feat, point_size, condition, fig):
        if not feat:
            return dash.no_update
        if feat not in var_df["feature"]:
            return dash.no_update
        condition = condition or data["default_condition"]
        try:
            factor = condition.split(": ")[0]
            test = condition.split(": ")[1]
            row = var_df.loc[(var_df["feature"]==feat) & (var_df["factor"]==factor) & (var_df["test"]==test), :].iloc[0]
            lfc = float(row["log2FoldChange"])
            padj = float(row["padj"])
            y = -np.log10(max(padj, 1e-300))
        except (KeyError, TypeError, ValueError):
            return dash.no_update
        p = Patch()
        p["data"][-1]["x"] = [lfc]
        p["data"][-1]["y"] = [y]
        p["data"][-1]["marker"]["size"] = (point_size or 5) * 2
        p["data"][-1]["marker"]["line"]["width"] = max(2, (point_size or 5) * 0.6)
        return p

    # ── Tab switching ──────────────────────────────────────────────────
    @app.callback(
        [Output("active-detail-tab", "data"),
         Output("metadata-tab-content", "style"),
         Output("boxplot-tab-content", "style"),
         Output("tab-metadata-btn", "className"),
         Output("tab-boxplot-btn", "className")],
        [Input("tab-metadata-btn", "n_clicks"),
         Input("tab-boxplot-btn", "n_clicks")],
        prevent_initial_call=True,
    )
    def switch_tab(_m, _b):
        t = callback_context.triggered_id
        if t == "tab-boxplot-btn":
            return ("boxplot",
                    {"display": "none"}, {"display": "block"},
                    "cs-tab", "cs-tab active")
        return ("metadata",
                {"display": "block"}, {"display": "none"},
                "cs-tab active", "cs-tab")

    # ── Box plot generation ────────────────────────────────────────────
    @app.callback(
        [Output("violin-plot", "figure"),
         Output("violin-status", "children")],
        Input("violin-button", "n_clicks"),
        [State("active-feature", "data"),
         State("violin-obs-dropdown", "value"),
         State("violin-layer-dropdown", "value"),
         State("boxplot-log-toggle", "checked"),
         State("boxplot-sort-dropdown", "value"),
         State("theme-store", "data")],
        prevent_initial_call=True,
    )
    def generate_boxplot(_n, feat, obs_col, layer, log_t, sort_by, theme):
        if not feat:
            return go.Figure(), dmc.Text("Select a point first.",
                                          c="dimmed", size="xs")
        if not obs_col:
            return go.Figure(), dmc.Text("Select a grouping column.",
                                          c="dimmed", size="xs")
        if feat not in var_name_to_idx:
            return go.Figure(), dmc.Text(f"'{feat}' not found.",
                                          c="red", size="xs")
        idx = var_name_to_idx[feat]

        try:
            col = adata.X[:, idx] if layer == "X" else adata.layers[layer][:, idx]
            vals = col.toarray().flatten() if hasattr(col, "toarray") else \
                np.asarray(col).flatten()
        except Exception as e:
            return go.Figure(), dmc.Text(f"Error: {e}", c="red", size="xs")

        if log_t:
            vals = np.log1p(vals)
        groups = adata.obs[obs_col].values
        if hasattr(groups, "categories"):
            groups = groups.astype(str)

        box_df = pd.DataFrame({"expression": vals, obs_col: groups})
        sort_by = sort_by or "median_desc"
        gs = box_df.groupby(obs_col)["expression"]
        if sort_by.startswith("median"):
            order_vals = gs.median()
        elif sort_by.startswith("mean"):
            order_vals = gs.mean()
        elif sort_by.startswith("max"):
            order_vals = gs.max()
        elif sort_by.startswith("min"):
            order_vals = gs.min()
        else:
            order_vals = None
        ascending = sort_by.endswith("_asc")
        if order_vals is not None:
            sorted_groups = order_vals.sort_values(ascending=ascending).index.tolist()
        else:
            sorted_groups = sorted(box_df[obs_col].unique(), reverse=not ascending)
        sorted_groups = list(reversed(sorted_groups))

        text_color = "#f5f5f4" if theme == "dark" else "#1c1917"
        axis_color = "#52525b" if theme == "dark" else "#a8a29e"
        grid_color = "rgba(255,255,255,0.04)" if theme == "dark" \
            else "rgba(0,0,0,0.04)"

        fig = go.Figure()
        for g in sorted_groups:
            gd = box_df.loc[box_df[obs_col] == g, "expression"]
            fig.add_trace(go.Box(x=gd, name=str(g), boxpoints="outliers",
                                  orientation="h"))
        n_groups = len(sorted_groups)
        plot_height = max(220, 32 * n_groups + 80)
        x_label = "Expression (log1p)" if log_t else "Expression"
        fig.update_layout(
            xaxis=dict(title=x_label, gridcolor=grid_color,
                       tickfont=dict(color=text_color),
                       title_font=dict(color=text_color)),
            yaxis=dict(tickfont=dict(color=text_color), gridcolor=grid_color),
            showlegend=False,
            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            font=dict(family="Inter, sans-serif", color=text_color),
            margin=dict(l=10, r=10, t=10, b=40), height=plot_height,
        )
        return fig, ""

    # ── Quick-filter modal: open via header trigger / cancel ───────────
    @app.callback(
        Output("search-modal", "opened"),
        [Input("cs-search-trigger", "n_clicks"),
         Input("search-cancel-btn", "n_clicks")],
        State("search-modal", "opened"),
        prevent_initial_call=True,
    )
    def open_search_modal(_open, _cancel, opened):
        if callback_context.triggered_id == "search-cancel-btn":
            return False
        return not bool(opened)

    # ── Help modal: open / close ──────────────────────────────────────
    @app.callback(
        Output("help-modal", "opened"),
        Input("cs-help-btn", "n_clicks"),
        State("help-modal", "opened"),
        prevent_initial_call=True,
    )
    def toggle_help(_n, opened):
        return not bool(opened)

    # ── Saved views: open drawer (from header bookmark) ────────────────
    @app.callback(
        Output("saved-drawer", "opened"),
        Input("header-saved-btn", "n_clicks"),
        State("saved-drawer", "opened"),
        prevent_initial_call=True,
    )
    def toggle_drawer(_n, opened):
        return not bool(opened)

    # ── Recent-views pills (header) ────────────────────────────────────
    @app.callback(
        Output("recent-views-pills", "children"),
        [Input("recent-views", "data"),
         Input("saved-views", "data")],
    )
    def render_recent_pills(recent, views):
        views = views or []
        valid = {v.get("name") for v in views}
        # 1) start with explicit recent activity (load OR save), most recent
        #    first, filtered to names that still exist in saved-views.
        ordered = []
        for n in (recent or []):
            if n in valid and n not in ordered:
                ordered.append(n)
        # 2) pad with the rest of saved-views in reverse-append order
        #    (= most recently saved first) until we have 5 pills.
        for v in reversed(views):
            if len(ordered) >= 5:
                break
            n = v.get("name")
            if n and n not in ordered:
                ordered.append(n)
        if not ordered:
            return []
        return [
            html.Button(
                n,
                id={"type": "recent-pill", "name": n},
                className="cs-recent-pill",
                title=f"Load saved view '{n}'",
            )
            for n in ordered[:5]
        ]

    # ── Saved views: render list ───────────────────────────────────────
    @app.callback(
        Output("saved-views-list", "children"),
        Input("saved-views", "data"),
    )
    def render_saved_views(views):
        views = views or []
        if not views:
            return dmc.Text("No saved views yet.", size="xs", c="dimmed",
                             style={"marginTop": 12})
        return [
            html.Div(className="cs-saved-card", children=[
                html.Div(className="info", children=[
                    html.Div(v.get("name", f"View {i + 1}"), className="name"),
                    html.Div(
                        f"{v.get('condition', '?')} · {len(v.get('filter_history', []))} filter(s)",
                        className="sub",
                    ),
                ]),
                html.Div(className="actions", children=[
                    html.Button("Load",
                                id={"type": "view-load", "index": i},
                                className="cs-action-btn"),
                    html.Button(
                        _svg_icon(ICON_TRASH, 13),
                        id={"type": "view-delete", "index": i},
                        className="cs-btn-danger",
                        title="Delete",
                    ),
                ]),
            ])
            for i, v in enumerate(views)
        ]

    # ── Saved views: save current ──────────────────────────────────────
    @app.callback(
        [Output("saved-views", "data"),
         Output("save-view-status", "children"),
         Output("save-view-name", "value"),
         Output("recent-views", "data")],
        Input("save-view-btn", "n_clicks"),
        [State("save-view-name", "value"),
         State("saved-views", "data"),
         State("recent-views", "data"),
         State("selected-condition", "data"),
         State("color-dropdown", "value"),
         State("filter-history", "data"),
         State("show-all-toggle", "checked"),
         State("size-slider", "value"),
         State("opacity-slider", "value"),
         State("max-legend-input", "value"),
         State("tooltip-checklist", "value")],
        prevent_initial_call=True,
    )
    def save_view(_n, name, views, recent, condition, color_col, history,
                  show_all, ps, op, ml, tt):
        if not _n:
            return dash.no_update, dash.no_update, dash.no_update, dash.no_update
        name = (name or "").strip()
        if not name:
            return (dash.no_update,
                    dmc.Text("Please enter a name.", size="xs", c="red"),
                    dash.no_update, dash.no_update)
        new_view = {
            "name": name,
            "condition": condition,
            "color_column": color_col,
            "filter_history": history or [],
            "show_all": bool(show_all),
            "point_size": ps,
            "opacity": op,
            "max_legend_classes": ml,
            "tooltip_columns": tt or [],
        }
        views = (views or []) + [new_view]
        ok = _persist_saved_views(config_path, views)
        msg = (dmc.Text(f'Saved "{name}".', size="xs", c="dimmed") if ok else
               dmc.Text("Saved (in memory only — config write failed)",
                        size="xs", c="orange"))
        new_recent = [n for n in (recent or []) if n != name]
        new_recent.insert(0, name)
        return views, msg, "", new_recent[:5]

    # ── Saved views: delete ────────────────────────────────────────────
    @app.callback(
        [Output("saved-views", "data", allow_duplicate=True),
         Output("save-view-status", "children", allow_duplicate=True),
         Output("recent-views", "data", allow_duplicate=True)],
        Input({"type": "view-delete", "index": ALL}, "n_clicks"),
        [State("saved-views", "data"),
         State("recent-views", "data")],
        prevent_initial_call=True,
    )
    def delete_view(_clicks, views, recent):
        if not any(_clicks or []):
            return dash.no_update, dash.no_update, dash.no_update
        idx = callback_context.triggered_id["index"]
        views = list(views or [])
        if idx >= len(views):
            return dash.no_update, dash.no_update, dash.no_update
        removed = views.pop(idx)
        removed_name = removed.get("name", "")
        ok = _persist_saved_views(config_path, views)
        msg = (dmc.Text(f'Deleted "{removed_name}".',
                         size="xs", c="dimmed") if ok else
               dmc.Text("Deleted (in memory only — config write failed)",
                        size="xs", c="orange"))
        new_recent = [n for n in (recent or []) if n != removed_name]
        return views, msg, new_recent

    # ── Saved views: load (from drawer OR from a header recent pill) ──
    @app.callback(
        [Output("selected-condition", "data", allow_duplicate=True),
         Output("project-pill-condition", "children", allow_duplicate=True),
         Output("cs-plot-title", "children", allow_duplicate=True),
         Output("color-dropdown", "value", allow_duplicate=True),
         Output("filter-history", "data", allow_duplicate=True),
         Output("filter-history-table", "children", allow_duplicate=True),
         Output("selected-reads", "data", allow_duplicate=True),
         Output("show-all-toggle", "checked"),
         Output("size-slider", "value"),
         Output("opacity-slider", "value"),
         Output("max-legend-input", "value"),
         Output("tooltip-checklist", "value", allow_duplicate=True),
         Output("saved-drawer", "opened", allow_duplicate=True),
         Output("save-view-status", "children", allow_duplicate=True),
         Output("recent-views", "data", allow_duplicate=True)],
        [Input({"type": "view-load", "index": ALL}, "n_clicks"),
         Input({"type": "recent-pill", "name": ALL}, "n_clicks")],
        [State("saved-views", "data"),
         State("recent-views", "data")],
        prevent_initial_call=True,
    )
    def load_view(_drawer_clicks, _pill_clicks, views, recent):
        t = callback_context.triggered_id
        if not isinstance(t, dict):
            return (dash.no_update,) * 15
        v = None
        if t.get("type") == "view-load":
            if not any(_drawer_clicks or []):
                return (dash.no_update,) * 15
            idx = t.get("index", -1)
            if 0 <= idx < len(views or []):
                v = views[idx]
        elif t.get("type") == "recent-pill":
            if not any(_pill_clicks or []):
                return (dash.no_update,) * 15
            v = next((vv for vv in (views or [])
                       if vv.get("name") == t.get("name")), None)
        if v is None:
            return (dash.no_update,) * 15
        cond = v.get("condition") or data["default_condition"]
        if cond not in available_conditions:
            cond = data["default_condition"]
        history = v.get("filter_history", []) or []
        last_reads = history[-1]["reads"] if history else None
        msg = dmc.Text(f'Loaded "{v.get("name", "")}".', size="xs", c="dimmed")
        loaded_name = v.get("name", "")
        new_recent = [n for n in (recent or []) if n != loaded_name]
        if loaded_name:
            new_recent.insert(0, loaded_name)
        return (
            cond, cond,
            _plot_title_children(cond, available_conditions),
            v.get("color_column") or "regulation",
            history,
            _render_filter_table(history),
            last_reads,
            bool(v.get("show_all", False)),
            v.get("point_size") or 5,
            v.get("opacity") or 0.7,
            v.get("max_legend_classes") if v.get("max_legend_classes") is not None else 20,
            v.get("tooltip_columns") or [],
            False,
            msg,
            new_recent[:5],
        )

    # ── Export CSV ─────────────────────────────────────────────────────
    @app.callback(
        Output("export-download", "data"),
        Input("export-btn", "n_clicks"),
        [State("selected-condition", "data"),
         State("selected-reads", "data"),
         State("show-all-toggle", "checked")],
        prevent_initial_call=True,
    )
    def export_csv(_n, condition, selected_reads, show_all):
        if not _n:
            return dash.no_update
        condition = condition or data["default_condition"]
        df = build_volcano_df(var_df, annotation_cols, condition, config)
        if selected_reads is not None:
            df = df[df["feature"].isin(selected_reads)]
        if not show_all:
            df = df[df["regulation"] != "Not significant"]
        df = df.copy()
        df.insert(1, "condition", condition)
        buf = io.StringIO()
        df.to_csv(buf, sep="\t", index=False)
        fname = f"volcano_{condition}_{len(df)}features.tsv"
        return dict(content=buf.getvalue(), filename=fname)

    return app


# ─── Main ─────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Interactive volcano plot explorer (Clean Scientific)")
    parser.add_argument("--counts", required=True, help="Path to .h5ad file containing counts")
    parser.add_argument("--deseq2-results", required=True, help="Path to DESeq2 results")
    parser.add_argument("--annotations",default=None, help="Path to feature annotation")
    parser.add_argument("--config", required=True, help="Path to YAML configuration file")
    parser.add_argument("-p", "--port", type=int, default=None,
                        help="Port to serve on (overrides config file)")
    args = parser.parse_args()

    config_path = os.path.realpath(args.config)
    with open(config_path) as f:
        config = yaml.safe_load(f) or {}

    config["deseq2_h5ad"] = args.counts
    config["deseq2_results"] = args.deseq2_results
    config["analysis_table"] = args.annotations


    app = create_app(config, config_path=config_path)
    app.run(
        host=config.get("host", "127.0.0.1"),
        port=args.port if args.port is not None else config.get("port", 8050),
        debug=config.get("debug", False),
    )


if __name__ == "__main__":
    main()
