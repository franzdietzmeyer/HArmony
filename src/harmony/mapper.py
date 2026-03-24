"""Reference-map generation from 18-subtype aligned HA FASTA."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict

from Bio import SeqIO


SUBTYPE_PATTERN = re.compile(r"\b(H([1-9]|1[0-8]))N\d+\b", re.IGNORECASE)

SITE_A = {122, 124, 126, *range(128, 147), *range(150, 155)}
SITE_B = {*range(155, 161), 163, 164, *range(186, 199)}
SITE_C = {*range(50, 55), *range(273, 281)}
SITE_D = {121, 167, *range(171, 175), *range(201, 249)}
SITE_E = {*range(62, 95), *range(260, 263)}


def _extract_h_subtype(header: str) -> str:
    match = SUBTYPE_PATTERN.search(header)
    if not match:
        raise ValueError(f"Could not parse H-subtype from FASTA header: {header}")
    return match.group(1).upper().split("N", 1)[0]


def _site_for_coord(coord_numeric: int) -> str:
    if coord_numeric in SITE_A:
        return "Site A"
    if coord_numeric in SITE_B:
        return "Site B"
    if coord_numeric in SITE_C:
        return "Site C"
    if coord_numeric in SITE_D:
        return "Site D"
    if coord_numeric in SITE_E:
        return "Site E"
    return "None"


def _domain_for_coord(coord_numeric: int) -> str:
    return "HA1" if coord_numeric < 330 else "HA2"


def _build_reference_map(aligned_reference: str) -> Dict[str, Dict[str, str]]:
    mapping: Dict[str, Dict[str, str]] = {}
    base_idx = 0
    gap_suffix = 0

    for hmm_state, ref_char in enumerate(aligned_reference, start=1):
        if ref_char != "-":
            base_idx += 1
            gap_suffix = 0
            label = str(base_idx)
        else:
            suffix = chr(ord("a") + gap_suffix)
            label = f"{base_idx}{suffix}"
            gap_suffix += 1
        coord_numeric = max(0, base_idx)
        mapping[str(hmm_state)] = {
            "coord": label,
            "site": _site_for_coord(coord_numeric),
            "domain": _domain_for_coord(coord_numeric),
        }
    return mapping


def generate_maps(msa_fasta: str | Path, output_dir: str | Path) -> Dict[str, Path]:
    """Generate H1-H18 reference maps from aligned subtype FASTA.

    Output files are written as `h1_map.json` ... `h18_map.json`.
    """
    msa_fasta = Path(msa_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(str(msa_fasta), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {msa_fasta}")

    aln_len = len(records[0].seq)
    for rec in records:
        if len(rec.seq) != aln_len:
            raise ValueError("MSA sequences must have identical aligned length.")

    by_subtype: Dict[str, str] = {}
    for rec in records:
        subtype = _extract_h_subtype(f"{rec.id} {rec.description}")
        by_subtype[subtype] = str(rec.seq).upper()

    expected = {f"H{i}" for i in range(1, 19)}
    missing = sorted(expected - set(by_subtype))
    if missing:
        raise ValueError(f"Missing expected subtypes in MSA: {', '.join(missing)}")

    written: Dict[str, Path] = {}
    for subtype in sorted(expected, key=lambda x: int(x[1:])):
        mapping = _build_reference_map(by_subtype[subtype])
        out_path = output_dir / f"{subtype.lower()}_map.json"
        out_path.write_text(json.dumps(mapping, indent=2), encoding="utf-8")
        written[subtype] = out_path
    return written
