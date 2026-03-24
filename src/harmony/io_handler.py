"""Input and output path handling for harmony batch processing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO


@dataclass
class SequenceRecord:
    """Minimal input sequence record for mapping."""

    query_id: str
    sequence: str
    source_file: Path


FASTA_EXTENSIONS = {".fasta", ".fa", ".faa", ".fas"}


def collect_fasta_files(input_path: Path) -> List[Path]:
    """Collect FASTA files from a file path or directory."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    if input_path.is_file():
        if input_path.suffix.lower() not in FASTA_EXTENSIONS:
            raise ValueError(f"Input file must be FASTA ({FASTA_EXTENSIONS}): {input_path}")
        return [input_path]

    files = sorted(
        p for p in input_path.iterdir() if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    )
    if not files:
        raise ValueError(f"No FASTA files found in directory: {input_path}")
    return files


def load_sequences(input_path: Path) -> List[SequenceRecord]:
    """Load sequences from one FASTA or many FASTA files in directory."""
    records: List[SequenceRecord] = []
    for fasta_file in collect_fasta_files(input_path):
        parsed = list(SeqIO.parse(str(fasta_file), "fasta"))
        if not parsed:
            continue
        for rec in parsed:
            records.append(
                SequenceRecord(
                    query_id=rec.id,
                    sequence=str(rec.seq),
                    source_file=fasta_file,
                )
            )
    if not records:
        raise ValueError(f"No FASTA sequences found in input path: {input_path}")
    return records


def ensure_output_dir(output_dir: Path) -> Path:
    """Create output directory if needed and return resolved path."""
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def safe_name(name: str) -> str:
    """Normalize sequence ID to filesystem-safe name."""
    keep = []
    for ch in name:
        if ch.isalnum() or ch in ("-", "_", "."):
            keep.append(ch)
        else:
            keep.append("_")
    normalized = "".join(keep).strip("._")
    return normalized or "sequence"


def unique_names(records: Iterable[SequenceRecord]) -> dict[str, str]:
    """Build stable unique file stems for each query id."""
    counts: dict[str, int] = {}
    out: dict[str, str] = {}
    for record in records:
        base = safe_name(record.query_id)
        n = counts.get(base, 0)
        counts[base] = n + 1
        out[record.query_id] = base if n == 0 else f"{base}_{n+1}"
    return out
