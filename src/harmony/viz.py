"""Linear HA protein map visualization utilities."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SITE_RANGES = {
    "Site A": [(122, 122), (124, 124), (126, 126), (128, 146), (150, 154)],
    "Site B": [(155, 160), (163, 163), (164, 164), (186, 198)],
    "Site C": [(50, 54), (273, 280)],
    "Site D": [(121, 121), (167, 167), (171, 174), (201, 248)],
    "Site E": [(62, 94), (260, 262)],
}

SITE_COLORS = {
    "Site A": "#3b82f6",
    "Site B": "#dc2626",
    "Site C": "#16a34a",
    "Site D": "#eab308",
    "Site E": "#c026d3",
}

STRUCTURAL_REGION_COLORS = {
    "Signal Peptide": "#d1d5db",
    "HA1": "#dbeafe",
    "Cleavage Site": "#ef4444",
    "Fusion Peptide": "#f97316",
    "HA2 Stalk": "#bbf7d0",
    "Transmembrane Domain": "#4b5563",
    "Cytoplasmic Tail": "#a855f7",
}

H3_REGION_RANGES = {
    "Signal Peptide": [(0, 0)],
    "HA1": [(1, 326)],
    "Cleavage Site": [(327, 330)],
    "Fusion Peptide": [(330, 352)],
    "HA2 Stalk": [(353, 513)],
    "Transmembrane Domain": [(514, 540)],
    "Cytoplasmic Tail": [(541, 800)],
}


def _coord_to_int(coord: str) -> int | None:
    match = re.match(r"(\d+)", str(coord))
    return int(match.group(1)) if match else None


def _merge_ranges(coords: list[int]) -> list[tuple[int, int]]:
    if not coords:
        return []
    out: list[tuple[int, int]] = []
    start = prev = coords[0]
    for c in coords[1:]:
        if c <= prev + 1:
            prev = c
            continue
        out.append((start, prev))
        start = prev = c
    out.append((start, prev))
    return out


def _structural_ranges(rows: list[dict], coord_col: str, ref_subtype: str) -> dict[str, list[tuple[int, int]]]:
    by_region: dict[str, list[int]] = {}
    for row in rows:
        region = str(row.get("Structural_Region", "")).strip()
        coord = _coord_to_int(row.get(coord_col, ""))
        if not region or coord is None:
            continue
        by_region.setdefault(region, []).append(coord)
    if not by_region and ref_subtype.upper() == "H3":
        return H3_REGION_RANGES
    out: dict[str, list[tuple[int, int]]] = {}
    for region, coords in by_region.items():
        out[region] = _merge_ranges(sorted(set(coords)))
    return out


def generate_protein_map(
    mapping_df: Iterable[dict],
    output_path: str | Path,
    sequence_id: str = "Unknown",
    ref_subtype: str = "H3",
    pathogenicity_score: str = "",
) -> None:
    """Render a linear protein schematic from mapped residue rows."""
    rows = list(mapping_df)
    coord_col = f"{ref_subtype.upper()}_Coord"
    coords = [_coord_to_int(r.get(coord_col, "")) for r in rows]
    coords = [c for c in coords if c is not None]
    if not coords:
        return

    cmin, cmax = min(coords), max(coords)
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    fig, ax = plt.subplots(figsize=(14, 3.8))

    # Background protein bar
    main_lw = 12
    site_lw = main_lw * 1.2  # 20% taller than backbone bar

    # Structural region shading track.
    region_ranges = _structural_ranges(rows, coord_col, ref_subtype)
    for region, ranges in region_ranges.items():
        color = STRUCTURAL_REGION_COLORS.get(region)
        if not color:
            continue
        for start, end in ranges:
            if end < cmin or start > cmax:
                continue
            s = max(start, cmin)
            e = min(end, cmax)
            ax.axvspan(s, e, color=color, alpha=0.18, zorder=0)

    ax.hlines(y=0, xmin=cmin, xmax=cmax, color="#d1d5db", linewidth=main_lw, zorder=1)

    # Site overlays use Burke & Smith H3 coordinates.
    for site_name, ranges in SITE_RANGES.items():
        color = SITE_COLORS[site_name]
        for start, end in ranges:
            if end < cmin or start > cmax:
                continue
            s = max(start, cmin)
            e = min(end, cmax)
            ax.hlines(y=0, xmin=s, xmax=e, color=color, linewidth=site_lw, zorder=2)
            ax.hlines(y=0, xmin=s, xmax=e, color="#111827", linewidth=0.7, zorder=3)

    # HA1/HA2 cleavage marker at H3 329/330.
    if cmin <= 330 <= cmax:
        hpai = pathogenicity_score.startswith("HPAI")
        marker_color = "#ef4444" if hpai else "#6b7280"
        label = "HA1/HA2"
        ax.axvline(330, color=marker_color, linestyle="--", linewidth=2.0 if hpai else 1.5, zorder=3)
        ax.text(330 + 1.5, 0.26, label, fontsize=9, color=marker_color if hpai else "#4b5563")
        if hpai:
            ax.text(
                0.99,
                0.95,
                "WARNING: HPAI Detected",
                transform=ax.transAxes,
                ha="right",
                va="top",
                color="#991b1b",
                fontweight="bold",
                fontsize=9,
                bbox={
                    "boxstyle": "round,pad=0.25",
                    "facecolor": "#fee2e2",
                    "edgecolor": "#991b1b",
                    "linewidth": 1.0,
                },
            )

    # Glycan "lollipop" markers.
    glycan_points: list[tuple[int, str]] = []
    for row in rows:
        gly = str(row.get("Glycan", "")).strip()
        if not gly:
            continue
        c = _coord_to_int(row.get(coord_col, ""))
        if c is not None:
            glycan_points.append((c, gly))

    drawn = set()
    used_label_x: list[int] = []
    for x, label in glycan_points:
        if x in drawn:
            continue
        drawn.add(x)
        ax.vlines(x, 0.02, 0.22, color="#111827", linewidth=1.0, alpha=0.7, zorder=4)
        if "Conserved" in label:
            t = "Conserved"
            marker = "o"  # circle
        elif "Gained" in label:
            t = "Gained"
            marker = "D"  # diamond
        elif "Lost" in label:
            t = "Lost"
            marker = "D"
        else:
            t = "Detected"
            marker = "o"

        ax.scatter([x], [0.24], s=65, c="#f97316", edgecolors="#111827", marker=marker, zorder=5)

        close_count = sum(1 for prev in used_label_x if abs(x - prev) <= 10)
        y_offset = 0.34 + (0.06 * (close_count % 3))
        used_label_x.append(x)
        ax.text(x, y_offset, t, fontsize=8, rotation=45, ha="left", va="bottom")

    tick_step = 50 if (cmax - cmin) <= 800 else 100
    start_tick = (cmin // tick_step) * tick_step
    ax.set_xticks(list(range(start_tick, cmax + tick_step, tick_step)))
    ax.set_xlim(cmin - 5, cmax + 5)
    ax.set_ylim(-0.35, 0.58)
    ax.set_yticks([])
    ax.set_xlabel(f"{ref_subtype.upper()} Coordinates")
    ax.set_title(f"HA Structural Map | Sequence: {sequence_id} | Reference: {ref_subtype.upper()}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.3)

    site_legend = [Patch(color=SITE_COLORS[s], label=s) for s in ("Site A", "Site B", "Site C", "Site D", "Site E")]
    region_legend = [
        Patch(color=STRUCTURAL_REGION_COLORS[r], label=r)
        for r in (
            "Signal Peptide",
            "HA1",
            "Cleavage Site",
            "Fusion Peptide",
            "HA2 Stalk",
            "Transmembrane Domain",
            "Cytoplasmic Tail",
        )
    ]
    leg1 = ax.legend(handles=site_legend, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.16))
    ax.add_artist(leg1)
    ax.legend(handles=region_legend, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.29))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_comparison(
    df1: Iterable[dict],
    df2: Iterable[dict],
    output_path: str | Path,
    seq1_id: str = "Sequence_1",
    seq2_id: str = "Sequence_2",
    ref_subtype: str = "H3",
) -> None:
    """Plot stacked side-by-side structural maps with mutation connectors."""
    rows1 = list(df1)
    rows2 = list(df2)
    coord_col = f"{ref_subtype.upper()}_Coord"

    def to_map(rows: list[dict]) -> dict[int, dict]:
        out: dict[int, dict] = {}
        for r in rows:
            c = _coord_to_int(r.get(coord_col, ""))
            if c is not None:
                out[c] = r
        return out

    m1 = to_map(rows1)
    m2 = to_map(rows2)
    coords = sorted(set(m1) | set(m2))
    if not coords:
        return

    cmin, cmax = min(coords), max(coords)
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    fig, ax = plt.subplots(figsize=(14, 5.6))

    y1, y2 = 0.65, 0.05
    main_lw = 11
    site_lw = main_lw * 1.2

    region_ranges = _structural_ranges(rows1 + rows2, coord_col, ref_subtype)
    for region, ranges in region_ranges.items():
        color = STRUCTURAL_REGION_COLORS.get(region)
        if not color:
            continue
        for start, end in ranges:
            if end < cmin or start > cmax:
                continue
            s = max(start, cmin)
            e = min(end, cmax)
            ax.axvspan(s, e, color=color, alpha=0.14, zorder=0)

    ax.hlines(y1, cmin, cmax, color="#d1d5db", linewidth=main_lw, zorder=1)
    ax.hlines(y2, cmin, cmax, color="#d1d5db", linewidth=main_lw, zorder=1)

    for site_name, ranges in SITE_RANGES.items():
        color = SITE_COLORS[site_name]
        for start, end in ranges:
            if end < cmin or start > cmax:
                continue
            s = max(start, cmin)
            e = min(end, cmax)
            for y in (y1, y2):
                ax.hlines(y, s, e, color=color, linewidth=site_lw, zorder=2)
                ax.hlines(y, s, e, color="#111827", linewidth=0.7, zorder=3)

    for c in coords:
        r1 = m1.get(c, {})
        r2 = m2.get(c, {})
        if str(r1.get("Residue", "-")) != str(r2.get("Residue", "-")):
            ax.vlines(c, y2 + 0.06, y1 - 0.06, color="#dc2626", alpha=0.35, linewidth=1.0, zorder=0)

    def draw_glycans(rows: list[dict], ybase: float) -> None:
        used_x: list[int] = []
        seen: set[int] = set()
        for row in rows:
            gly = str(row.get("Glycan", "")).strip()
            if not gly:
                continue
            x = _coord_to_int(row.get(coord_col, ""))
            if x is None or x in seen:
                continue
            seen.add(x)
            marker = "D" if ("Gained" in gly or "Lost" in gly) else "o"
            ax.vlines(x, ybase + 0.02, ybase + 0.14, color="#111827", linewidth=0.9, alpha=0.7, zorder=4)
            ax.scatter([x], [ybase + 0.16], marker=marker, s=58, c="#f97316", edgecolors="#111827", zorder=5)
            label = "Gained" if "Gained" in gly else "Conserved" if "Conserved" in gly else "Detected"
            close = sum(1 for px in used_x if abs(x - px) <= 10)
            used_x.append(x)
            ax.text(x, ybase + 0.22 + 0.05 * (close % 3), label, fontsize=7, rotation=45, ha="left", va="bottom")

    draw_glycans(rows1, y1)
    draw_glycans(rows2, y2)

    if cmin <= 330 <= cmax:
        ax.axvline(330, color="#6b7280", linestyle="--", linewidth=1.4, zorder=3)

    tick_step = 50 if (cmax - cmin) <= 800 else 100
    start_tick = (cmin // tick_step) * tick_step
    ax.set_xticks(list(range(start_tick, cmax + tick_step, tick_step)))
    ax.set_xlim(cmin - 5, cmax + 5)
    ax.set_ylim(-0.2, 1.1)
    ax.set_yticks([y1, y2])
    ax.set_yticklabels([seq1_id, seq2_id])
    ax.set_xlabel(f"{ref_subtype.upper()} Coordinates")
    ax.set_title(f"HA Structural Comparison | {seq1_id} vs {seq2_id} | Reference: {ref_subtype.upper()}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.3)

    site_legend = [Patch(color=SITE_COLORS[s], label=s) for s in ("Site A", "Site B", "Site C", "Site D", "Site E")]
    region_legend = [
        Patch(color=STRUCTURAL_REGION_COLORS[r], label=r)
        for r in (
            "Signal Peptide",
            "HA1",
            "Cleavage Site",
            "Fusion Peptide",
            "HA2 Stalk",
            "Transmembrane Domain",
            "Cytoplasmic Tail",
        )
    ]
    leg1 = ax.legend(handles=site_legend, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.10))
    ax.add_artist(leg1)
    ax.legend(handles=region_legend, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.24))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
