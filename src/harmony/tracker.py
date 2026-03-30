"""Multi-sequence evolutionary tracking for HArmony.

Builds a master mutation matrix comparing one ancestor against many descendants
on a shared reference coordinate grid (default: H3).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from Bio import SeqIO

from .engine import HarmonyMapper
from .formatter import events_to_rows


@dataclass(frozen=True)
class TrackedSequence:
    """A mapped sequence with cached coordinate lookup."""

    name: str
    sequence: str
    rows: List[dict]


def _load_single_fasta(path: Path) -> Tuple[str, str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA sequence found in {path}")
    rec = records[0]
    return rec.id, str(rec.seq)


def _load_multi_fasta(path: Path) -> List[Tuple[str, str]]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA sequences found in {path}")
    return [(rec.id, str(rec.seq)) for rec in records]


def _coord_sort_key(coord: str) -> tuple[int, str]:
    s = str(coord)
    m = re.match(r"^(\d+)(.*)$", s)
    if not m:
        return (10**9, s)
    return (int(m.group(1)), m.group(2))


def _rows_by_coord(rows: Iterable[dict], coord_col: str) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for r in rows:
        c = str(r.get(coord_col, "")).strip()
        if c:
            out[c] = r
    return out


def _map_rows_cached(
    mapper: HarmonyMapper, ref_subtype: str, seq: str, cache: Dict[str, List[dict]]
) -> List[dict]:
    key = seq.strip().upper()
    if key in cache:
        return cache[key]
    events = mapper.map_sequence(seq)
    mapper.detect_glycans(seq, events)
    mapper.analyze_cleavage_site(events)
    rows = events_to_rows(events, ref_subtype=ref_subtype)
    cache[key] = rows
    return rows


def track_ancestor_against_descendants(
    ancestor_fasta: Path,
    descendants_fasta: Path,
    ref: str = "H3",
) -> tuple[pd.DataFrame, str, List[dict], List[str], List[List[dict]], List[str]]:
    """Core tracking routine returning matrix + mapped rows for downstream exporters.

    Returns:
      (df, ancestor_sequence, ancestor_rows, descendant_names, descendant_sequences, descendant_rows)
    """

    ref = ref.upper()
    coord_col = f"{ref}_Coord"

    mapper = HarmonyMapper(ref_subtype=ref)
    cache: Dict[str, List[dict]] = {}

    _, anc_seq = _load_single_fasta(ancestor_fasta)
    anc_rows = _map_rows_cached(mapper, ref, anc_seq, cache)

    desc_records = _load_multi_fasta(descendants_fasta)
    if not desc_records:
        raise ValueError("No descendant records found.")

    descendant_names: List[str] = []
    descendant_sequences: List[str] = []
    descendant_rows: List[List[dict]] = []
    for name, seq in desc_records:
        descendant_names.append(name)
        descendant_sequences.append(seq)
        descendant_rows.append(_map_rows_cached(mapper, ref, seq, cache))

    anc_by = _rows_by_coord(anc_rows, coord_col)
    desc_by = {
        name: _rows_by_coord(rows, coord_col)
        for name, rows in zip(descendant_names, descendant_rows)
    }

    all_coords = set(anc_by.keys())
    for by in desc_by.values():
        all_coords |= set(by.keys())
    sorted_coords = sorted(all_coords, key=_coord_sort_key)

    records: List[dict] = []
    for c in sorted_coords:
        arow = anc_by.get(c, {})
        ancestor_aa = str(arow.get("Residue", "-"))
        domain = str(arow.get("Structural_Region", "Unknown"))
        site = str(arow.get("Site", "None"))

        row: dict = {
            coord_col: c,
            "Structural_Region": domain,
            "Antigenic_Site": site,
            "Ancestor_AA": ancestor_aa,
        }

        any_delta = False
        for name in descendant_names:
            drow = desc_by.get(name, {}).get(c, {})
            daa = str(drow.get("Residue", "-"))
            row[name] = daa
            if daa != ancestor_aa:
                any_delta = True

        if any_delta:
            records.append(row)

    df = pd.DataFrame.from_records(records)
    if not df.empty:
        ordered = [coord_col, "Structural_Region", "Antigenic_Site", "Ancestor_AA"] + descendant_names
        df = df.reindex(columns=ordered)

    return df, anc_seq, anc_rows, descendant_names, descendant_sequences, descendant_rows


def build_mutation_tracking_matrix(
    ancestor_fasta: Path,
    descendants_fasta: Path,
    ref: str = "H3",
) -> pd.DataFrame:
    """Compare ancestor against many descendants and return delta-filtered matrix.

    Output columns:
      - <ref>_Coord (H3_Coord by default)
      - Structural_Region
      - Antigenic_Site
      - Ancestor_AA
      - one column per descendant (using FASTA record id)
    """

    df, _, _, _, _, _ = track_ancestor_against_descendants(
        ancestor_fasta=ancestor_fasta,
        descendants_fasta=descendants_fasta,
        ref=ref,
    )
    return df


def write_mutation_tracking_matrix(df: pd.DataFrame, output_dir: Path, ref: str = "H3") -> Path:
    ref = ref.upper()
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "mutation_tracking_matrix.csv"
    df.to_csv(out_path, index=False)
    return out_path
