"""3D structural script exporters for HArmony tracking outputs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import pandas as pd


def _coord_to_pdb_resi(coord: str) -> str | None:
    """Convert H3 coordinate labels to PDB residue identifiers.

    Supports:
      - '158' -> '158'
      - '329a' -> '329A'

    Ignores:
      - insertion labels like '158+ins' (no stable PDB numbering target)
    """

    s = str(coord).strip()
    if not s:
        return None
    if "+ins" in s:
        return None
    m = re.match(r"^(\d+)([a-zA-Z])?$", s)
    if not m:
        return None
    num = m.group(1)
    ins = m.group(2)
    return f"{num}{ins.upper()}" if ins else num


def _mutation_coords_from_matrix(df: pd.DataFrame, coord_col: str, descendant_cols: Sequence[str]) -> Set[str]:
    coords: Set[str] = set()
    if df.empty:
        return coords
    for _, row in df.iterrows():
        anc = str(row.get("Ancestor_AA", "-"))
        for col in descendant_cols:
            if str(row.get(col, "-")) != anc:
                coords.add(str(row.get(coord_col, "")).strip())
                break
    return {c for c in coords if c}


def _find_descendant_cols(df: pd.DataFrame) -> List[str]:
    core = {"Structural_Region", "Antigenic_Site", "Ancestor_AA"}
    coord_candidates = {c for c in df.columns if c.endswith("_Coord")}
    core |= coord_candidates
    return [c for c in df.columns if c not in core]


def _coord_col(df: pd.DataFrame, ref: str) -> str:
    ref = ref.upper()
    col = f"{ref}_Coord"
    if col in df.columns:
        return col
    # fallback: first *_Coord column
    for c in df.columns:
        if c.endswith("_Coord"):
            return c
    raise ValueError("Mutation matrix is missing a coordinate column (e.g., H3_Coord).")


def _site_ranges_h3() -> Dict[str, List[Tuple[int, int]]]:
    # Keep aligned with `harmony.viz.SITE_RANGES` (Burke & Smith).
    return {
        "Site A": [(122, 122), (124, 124), (126, 126), (128, 146), (150, 154)],
        "Site B": [(155, 160), (163, 163), (164, 164), (186, 198)],
        "Site C": [(50, 54), (273, 280)],
        "Site D": [(121, 121), (167, 167), (171, 174), (201, 248)],
        "Site E": [(62, 94), (260, 262)],
    }


def _canonical_structural_region_ranges_h3() -> Dict[str, List[Tuple[int, int]]]:
    return {
        "Signal Peptide": [(0, 0)],
        "HA1": [(1, 326)],
        "Cleavage Site": [(327, 330)],
        "Fusion Peptide": [(330, 352)],
        "HA2 Stalk": [(353, 513)],
        "Transmembrane Domain": [(514, 540)],
        "Cytoplasmic Tail": [(541, 900)],
    }


def _merge_int_ranges(coords: List[int]) -> List[Tuple[int, int]]:
    if not coords:
        return []
    coords = sorted(set(coords))
    out: List[Tuple[int, int]] = []
    start = prev = coords[0]
    for c in coords[1:]:
        if c <= prev + 1:
            prev = c
            continue
        out.append((start, prev))
        start = prev = c
    out.append((start, prev))
    return out


def _structural_region_ranges_from_matrix(df: pd.DataFrame, coord_col: str, ref: str) -> Dict[str, List[Tuple[int, int]]]:
    by_region: Dict[str, List[int]] = {}
    for _, row in df.iterrows():
        region = str(row.get("Structural_Region", "")).strip()
        coord = str(row.get(coord_col, "")).strip()
        if not region:
            continue
        m = re.match(r"^(\d+)", coord)
        if not m:
            continue
        by_region.setdefault(region, []).append(int(m.group(1)))

    if by_region:
        return {region: _merge_int_ranges(vals) for region, vals in by_region.items()}
    if ref.upper() == "H3":
        return _canonical_structural_region_ranges_h3()
    return {}


def _glycan_windows_by_coord(sequence: str, qpos_to_coord: Dict[int, str]) -> Set[str]:
    """Return coordinate labels for residues participating in any NX[S/T] motif window."""
    seq = sequence.strip().upper()
    motif_starts: List[int] = []
    for i in range(len(seq) - 2):
        tri = seq[i : i + 3]
        if tri[0] == "N" and tri[1] != "P" and tri[2] in {"S", "T"}:
            motif_starts.append(i + 1)  # 1-indexed

    coords: Set[str] = set()
    for start in motif_starts:
        for qpos in (start, start + 1, start + 2):
            coord = qpos_to_coord.get(qpos)
            if coord:
                coords.add(coord)
    return coords


def _qpos_to_coord_from_mapped_rows(rows: Iterable[dict], coord_col: str) -> Dict[int, str]:
    out: Dict[int, str] = {}
    for r in rows:
        qpos = r.get("Query_Pos")
        coord = str(r.get(coord_col, "")).strip()
        if isinstance(qpos, int) and coord:
            out[qpos] = coord
    return out


@dataclass(frozen=True)
class StructureScriptPaths:
    pymol_pml: Path
    chimerax_cxc: Path


def write_pymol_script(
    df: pd.DataFrame,
    output_path: Path,
    ref: str = "H3",
    reference_object: str = "ha",
    mutated_representation: str = "sticks",
) -> Path:
    """Write a PyMOL `.pml` that highlights all mutated residues in red."""
    coord_col = _coord_col(df, ref)
    descendant_cols = _find_descendant_cols(df)
    mut_coords = sorted(_mutation_coords_from_matrix(df, coord_col, descendant_cols))

    resi_tokens = [t for t in (_coord_to_pdb_resi(c) for c in mut_coords) if t]
    selection = "+".join(resi_tokens)

    region_ranges = _structural_region_ranges_from_matrix(df, coord_col, ref)
    region_colors = {
        "Signal Peptide": "grey80",
        "HA1": "palecyan",
        "Cleavage Site": "firebrick",
        "Fusion Peptide": "orange",
        "HA2 Stalk": "palegreen",
        "Transmembrane Domain": "grey30",
        "Cytoplasmic Tail": "purple",
    }

    lines: List[str] = []
    lines.append("# HArmony tracking: PyMOL macro")
    lines.append("# Load a reference HA first (e.g., 4FNK), then run: @mutation_tracking.pml")
    lines.append(f"# Expected reference object name: {reference_object}")
    lines.append("")
    lines.append("bg_color white")
    lines.append(f"color grey80, {reference_object}")
    lines.append(f"show cartoon, {reference_object}")
    lines.append("")
    if region_ranges:
        lines.append("# Structural region coloring")
        for region, ranges in region_ranges.items():
            color = region_colors.get(region)
            if not color:
                continue
            for start, end in ranges:
                lines.append(f"select harmony_{region.lower().replace(' ', '_')}_{start}_{end}, {reference_object} and resi {start}-{end}")
                lines.append(f"color {color}, harmony_{region.lower().replace(' ', '_')}_{start}_{end}")
        lines.append("")
    if selection:
        lines.append(f"select harmony_mutations, {reference_object} and resi {selection}")
        lines.append("color red, harmony_mutations")
        if mutated_representation == "spheres":
            lines.append("show spheres, harmony_mutations")
            lines.append("set sphere_scale, 0.35, harmony_mutations")
        else:
            lines.append("show sticks, harmony_mutations")
            lines.append("set stick_radius, 0.18, harmony_mutations")
    else:
        lines.append("# No mutations found in matrix (delta-filter produced empty set).")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def write_chimerax_script(
    df: pd.DataFrame,
    output_path: Path,
    ref: str = "H3",
    model: str = "#1",
    chain: str = "A",
    gained_glycan_coords: Iterable[str] = (),
) -> Path:
    """Write a ChimeraX `.cxc` that highlights mutations, antigenic sites, and gained glycans."""
    coord_col = _coord_col(df, ref)
    descendant_cols = _find_descendant_cols(df)
    mut_coords = sorted(_mutation_coords_from_matrix(df, coord_col, descendant_cols))

    mut_resi = [t for t in (_coord_to_pdb_resi(c) for c in mut_coords) if t]
    gained = [t for t in (_coord_to_pdb_resi(c) for c in gained_glycan_coords) if t]

    region_ranges = _structural_region_ranges_from_matrix(df, coord_col, ref)
    region_colors = {
        "Signal Peptide": "lightgray",
        "HA1": "lightblue",
        "Cleavage Site": "firebrick",
        "Fusion Peptide": "orange",
        "HA2 Stalk": "lightgreen",
        "Transmembrane Domain": "darkgray",
        "Cytoplasmic Tail": "purple",
    }

    # Pastel site colors (low-opacity). Mutations will be red, gained glycans blue.
    site_colors = {
        "Site A": "#a5b4fc",  # indigo-300
        "Site B": "#fca5a5",  # red-300
        "Site C": "#86efac",  # green-300
        "Site D": "#fde68a",  # amber-200
        "Site E": "#f0abfc",  # fuchsia-300
    }

    lines: List[str] = []
    lines.append("# HArmony tracking: ChimeraX macro")
    lines.append("# Open your reference structure first (e.g., 4FNK), then run: open mutation_tracking.cxc")
    lines.append("")
    lines.append(f"color {model} lightgray")
    lines.append("style %s cartoon" % model)
    lines.append("")
    if region_ranges:
        lines.append("# Structural region coloring")
        for region, ranges in region_ranges.items():
            color = region_colors.get(region)
            if not color:
                continue
            for start, end in ranges:
                sel = f"{model}/{chain}:{start}-{end}"
                lines.append(f"color {sel} {color}")
        lines.append("")
    lines.append("# Burke & Smith antigenic sites (H3 numbering), shown as low-opacity pastels")
    for site, ranges in _site_ranges_h3().items():
        color = site_colors.get(site, "lightgray")
        for start, end in ranges:
            sel = f"{model}/{chain}:{start}-{end}"
            lines.append(f"color {sel} {color}")
            # ChimeraX transparency is in percent (0 opaque, 100 fully transparent) for atoms.
            lines.append(f"transparency 70 {sel}")
    lines.append("")
    if mut_resi:
        lines.append("# Mutations vs ancestor")
        for r in mut_resi:
            sel = f"{model}/{chain}:{r}"
            lines.append(f"color {sel} red")
            lines.append(f"style {sel} stick")
    else:
        lines.append("# No mutations found in matrix (delta-filter produced empty set).")
    lines.append("")
    if gained:
        lines.append("# Gained glycans (NX[S/T] motif residues gained vs ancestor)")
        for r in gained:
            sel = f"{model}/{chain}:{r}"
            lines.append(f"color {sel} dodgerblue")
            lines.append(f"style {sel} stick")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def compute_gained_glycans(
    ancestor_sequence: str,
    ancestor_mapped_rows: Iterable[dict],
    descendant_sequences: Sequence[str],
    descendant_mapped_rows: Sequence[Iterable[dict]],
    coord_col: str,
) -> Set[str]:
    """Return coordinates that are in a glycan window in any descendant but not in ancestor."""
    anc_map = _qpos_to_coord_from_mapped_rows(ancestor_mapped_rows, coord_col)
    anc_win = _glycan_windows_by_coord(ancestor_sequence, anc_map)

    gained: Set[str] = set()
    for seq, rows in zip(descendant_sequences, descendant_mapped_rows):
        d_map = _qpos_to_coord_from_mapped_rows(rows, coord_col)
        d_win = _glycan_windows_by_coord(seq, d_map)
        gained |= (d_win - anc_win)
    return gained


def write_tracking_scripts(
    matrix_df: pd.DataFrame,
    output_dir: Path,
    ref: str = "H3",
    ancestor_sequence: str | None = None,
    ancestor_mapped_rows: Iterable[dict] | None = None,
    descendant_sequences: Sequence[str] | None = None,
    descendant_mapped_rows: Sequence[Iterable[dict]] | None = None,
) -> StructureScriptPaths:
    """Convenience wrapper to write both `.pml` and `.cxc` scripts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    coord_col = _coord_col(matrix_df, ref)

    gained: Set[str] = set()
    if (
        ancestor_sequence is not None
        and ancestor_mapped_rows is not None
        and descendant_sequences is not None
        and descendant_mapped_rows is not None
    ):
        gained = compute_gained_glycans(
            ancestor_sequence=ancestor_sequence,
            ancestor_mapped_rows=ancestor_mapped_rows,
            descendant_sequences=descendant_sequences,
            descendant_mapped_rows=descendant_mapped_rows,
            coord_col=coord_col,
        )

    pml = write_pymol_script(matrix_df, output_dir / "mutation_tracking.pml", ref=ref)
    cxc = write_chimerax_script(
        matrix_df,
        output_dir / "mutation_tracking.cxc",
        ref=ref,
        gained_glycan_coords=sorted(gained),
    )
    return StructureScriptPaths(pymol_pml=pml, chimerax_cxc=cxc)

