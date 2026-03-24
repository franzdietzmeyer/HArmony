"""Pairwise comparison logic for harmony."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO

from .engine import HarmonyMapper
from .formatter import events_to_rows


def _load_single_fasta(path: Path) -> Tuple[str, str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA sequence found in {path}")
    rec = records[0]
    return rec.id, str(rec.seq)


def _row_by_coord(rows: List[dict], ref_subtype: str) -> Dict[str, dict]:
    coord_col = f"{ref_subtype.upper()}_Coord"
    out: Dict[str, dict] = {}
    for r in rows:
        coord = str(r.get(coord_col, ""))
        if coord and coord != "-":
            out[coord] = r
    return out


def _glycan_window_by_coord(
    sequence: str, events: List, ref_subtype: str
) -> tuple[Dict[str, bool], Dict[str, str]]:
    """Map reference coordinates to glycan-window membership and motif context.

    For each NX[S/T] motif in sequence space, mark all three residues (N, X, S/T)
    as part of the glycan window after projection to reference coordinates.
    """
    sequence = sequence.upper()
    motif_start_to_context: Dict[int, str] = {}
    for i in range(len(sequence) - 2):
        tri = sequence[i : i + 3]
        if tri[0] == "N" and tri[1] != "P" and tri[2] in {"S", "T"}:
            motif_start_to_context[i + 1] = tri

    coord_col = f"{ref_subtype.upper()}_Coord"
    qpos_to_coord: Dict[int, str] = {}
    for row in events_to_rows(events, ref_subtype=ref_subtype):
        qpos = row.get("Query_Pos")
        coord = str(row.get(coord_col, ""))
        if isinstance(qpos, int) and coord:
            qpos_to_coord[qpos] = coord

    coord_is_glycan_window: Dict[str, bool] = {}
    coord_context: Dict[str, str] = {}
    for start, ctx in motif_start_to_context.items():
        for qpos in (start, start + 1, start + 2):
            coord = qpos_to_coord.get(qpos)
            if coord:
                coord_is_glycan_window[coord] = True
                coord_context[coord] = ctx

    return coord_is_glycan_window, coord_context


def compare_sequences(seq1: Path, seq2: Path, ref: str = "H3", only_deltas: bool = False) -> Dict:
    """Run pairwise map+compare between two sequences on a shared reference."""
    ref = ref.upper()
    mapper = HarmonyMapper(ref_subtype=ref)

    id1, s1 = _load_single_fasta(seq1)
    id2, s2 = _load_single_fasta(seq2)

    e1 = mapper.map_sequence(s1)
    g1 = mapper.detect_glycans(s1, e1)
    c1 = mapper.analyze_cleavage_site(e1)
    rows1 = events_to_rows(e1, ref_subtype=ref)

    e2 = mapper.map_sequence(s2)
    g2 = mapper.detect_glycans(s2, e2)
    c2 = mapper.analyze_cleavage_site(e2)
    rows2 = events_to_rows(e2, ref_subtype=ref)
    gly_win1, motif_ctx1 = _glycan_window_by_coord(s1, e1, ref)
    gly_win2, motif_ctx2 = _glycan_window_by_coord(s2, e2, ref)

    by1 = _row_by_coord(rows1, ref)
    by2 = _row_by_coord(rows2, ref)
    all_coords = sorted(set(by1) | set(by2), key=lambda x: (int("".join(ch for ch in x if ch.isdigit()) or 0), x))

    coord_col = f"{ref}_Coord"
    diffs: List[dict] = []
    site_counts = {"Site A": 0, "Site B": 0, "Site C": 0, "Site D": 0, "Site E": 0}
    glycan_change_coords: List[str] = []

    for coord in all_coords:
        r1 = by1.get(coord, {})
        r2 = by2.get(coord, {})
        res1 = str(r1.get("Residue", "-"))
        res2 = str(r2.get("Residue", "-"))
        has_gly1 = gly_win1.get(coord, False)
        has_gly2 = gly_win2.get(coord, False)
        seq1_glyco_status = "Present" if has_gly1 else "None"
        seq2_glyco_status = "Present" if has_gly2 else "None"
        if has_gly1 and has_gly2:
            relative = "Conserved"
        elif has_gly1 and not has_gly2:
            relative = "Lost in Seq 2"
        elif not has_gly1 and has_gly2:
            relative = "Gained in Seq 2"
        else:
            relative = "None"
        site = str(r1.get("Site", r2.get("Site", "None")))

        mutation = res1 != res2
        glycan_shift = relative in {"Lost in Seq 2", "Gained in Seq 2"}
        if only_deltas and not mutation and not glycan_shift:
            continue

        if mutation and site in site_counts:
            site_counts[site] += 1
        if glycan_shift:
            glycan_change_coords.append(coord)

        diffs.append(
            {
                coord_col: coord,
                "Seq1_ID": id1,
                "Seq2_ID": id2,
                "Seq1_Res": res1,
                "Seq2_Res": res2,
                "Site": site,
                "Seq1_Glyco_Status": seq1_glyco_status,
                "Seq2_Glyco_Status": seq2_glyco_status,
                "Relative_Glycan_Status": relative,
                "Is_Mutation": mutation,
                "Is_Glycan_Change": glycan_shift,
                "Motif_Context": motif_ctx1.get(coord) or motif_ctx2.get(coord, ""),
                "Pathogenicity_Seq1": c1["pathogenicity_score"],
                "Pathogenicity_Seq2": c2["pathogenicity_score"],
                "Cleavage_Motif_Seq1": c1["cleavage_motif"],
                "Cleavage_Motif_Seq2": c2["cleavage_motif"],
            }
        )

    pathogenicity_shift = c1["pathogenicity_score"] != c2["pathogenicity_score"]
    return {
        "seq1_id": id1,
        "seq2_id": id2,
        "ref_subtype": ref,
        "rows1": rows1,
        "rows2": rows2,
        "differences": diffs,
        "total_mutations": sum(1 for d in diffs if d["Is_Mutation"]),
        "antigenic_site_impact": site_counts,
        "glycan_change_coords": glycan_change_coords,
        "pathogenicity_shift": pathogenicity_shift,
        "pathogenicity_seq1": c1["pathogenicity_score"],
        "pathogenicity_seq2": c2["pathogenicity_score"],
        "cleavage_motif_seq1": c1["cleavage_motif"],
        "cleavage_motif_seq2": c2["cleavage_motif"],
        "glycan_summary_seq1": g1,
        "glycan_summary_seq2": g2,
    }


def write_comparison_report(path: Path, rows: List[dict], ref_subtype: str) -> None:
    """Write CSV report containing only differing coordinates."""
    coord_col = f"{ref_subtype.upper()}_Coord"
    fieldnames = [
        coord_col,
        "Seq1_ID",
        "Seq2_ID",
        "Seq1_Res",
        "Seq2_Res",
        "Site",
        "Seq1_Glyco_Status",
        "Seq2_Glyco_Status",
        "Relative_Glycan_Status",
        "Is_Mutation",
        "Is_Glycan_Change",
        "Motif_Context",
        "Pathogenicity_Seq1",
        "Pathogenicity_Seq2",
        "Cleavage_Motif_Seq1",
        "Cleavage_Motif_Seq2",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
