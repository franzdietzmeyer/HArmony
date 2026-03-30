"""Alignment engine and dynamic reference-coordinate mapping logic."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import date
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Literal

import pyhmmer

MapEventType = Literal["Match", "Insertion", "Deletion"]


@dataclass
class MappingEvent:
    """Single mapped event between query and HMM."""

    query_pos: int | None
    residue: str
    ref_coord: str
    site: str
    domain: str
    event_type: MapEventType
    hmm_state: int | None
    confidence: float
    glycan: str = ""
    region_label: str = ""


class HarmonyMapper:
    """Map Influenza HA amino-acid sequences to chosen reference coordinates."""

    VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZJUO")
    H3_BASELINE_GLYCANS = {8, 22, 38, 45, 63, 122, 126, 133, 165, 246, 285}

    def __init__(
        self,
        hmm_path: str | Path | None = None,
        map_path: str | Path | None = None,
        ref_subtype: str = "H3",
        model_version: str = "harmony-1.1.0-alpha",
    ) -> None:
        self.hmm_path = self._resolve_data_path(hmm_path, "flumap.hmm")
        self.ref_subtype = ref_subtype.upper()
        self.map_path = self._resolve_ref_map_path(map_path, self.ref_subtype)
        self.model_version = model_version

        with pyhmmer.plan7.HMMFile(str(self.hmm_path)) as hmm_file:
            hmm = hmm_file.read()
        if hmm is None:
            raise ValueError(f"Could not load HMM from {self.hmm_path}")
        self.hmm = hmm
        self.alphabet = hmm.alphabet
        self.pipeline = pyhmmer.plan7.Pipeline(self.alphabet)

        self.ref_map: Dict[str, Dict[str, str]] = json.loads(
            self.map_path.read_text(encoding="utf-8")
        )

    @staticmethod
    def _resolve_data_path(candidate: str | Path | None, filename: str) -> Path:
        if candidate is not None:
            path = Path(candidate)
            if not path.exists():
                raise FileNotFoundError(f"Data file not found: {path}")
            return path

        package_path = Path(str(resources.files("harmony").joinpath("data", filename)))
        if package_path.exists():
            return package_path

        cwd_path = Path(filename)
        if cwd_path.exists():
            return cwd_path
        raise FileNotFoundError(
            f"Could not find {filename} in package data directory or current working directory."
        )

    @staticmethod
    def _resolve_ref_map_path(candidate: str | Path | None, ref_subtype: str) -> Path:
        if candidate is not None:
            path = Path(candidate)
            if not path.exists():
                raise FileNotFoundError(f"Map file not found: {path}")
            return path

        maps_dir = Path(__file__).resolve().parent / "data" / "maps"
        target = maps_dir / f"{ref_subtype.lower()}_map.json"
        if target.exists():
            return target

        available = []
        if maps_dir.exists():
            available = sorted(p.stem.replace("_map", "").upper() for p in maps_dir.glob("*_map.json"))
        avail_text = ", ".join(available) if available else "none found"
        raise FileNotFoundError(
            f"Reference map for {ref_subtype} not found at {target}. "
            f"Available subtypes: {avail_text}"
        )

    @classmethod
    def _clean_and_validate_sequence(cls, sequence: str) -> str:
        cleaned = sequence.strip().replace(" ", "").replace("\n", "").replace(".", "").upper()
        if not cleaned:
            raise ValueError("Sequence is empty after cleanup.")
        invalid = sorted({aa for aa in cleaned if aa not in cls.VALID_AA})
        if invalid:
            raise ValueError(f"Sequence contains invalid residues: {','.join(invalid)}")
        return cleaned

    def _annotation_for_state(self, state_idx: int) -> Dict[str, str]:
        default = {"coord": f"state_{state_idx}", "site": "None", "domain": "Unknown"}
        raw: Any = self.ref_map.get(str(state_idx))
        if raw is None:
            return default
        if isinstance(raw, dict):
            coord = str(raw.get("coord", default["coord"]))
            site = str(raw.get("site", "None"))
            structural_region = self._structural_region_for_coord(coord)
            return {"coord": coord, "site": site, "domain": structural_region}
        # Backward compatibility with legacy string-only maps.
        coord = str(raw)
        return {"coord": coord, "site": "None", "domain": self._structural_region_for_coord(coord)}

    def _structural_region_for_coord(self, coord: str) -> str:
        """Assign high-resolution structural regions from mapped H3 coordinates.

        H3 coordinate mapping:
          - Signal Peptide (SP): coords < 1
          - HA1: 1-329 (excluding cleavage site if explicitly modeled)
          - Cleavage Site: boundary bridge residues around HA1/HA2 (H3 327-330)
          - Fusion Peptide (FP): HA2 1-23  => H3 330-352
          - HA2 Stalk: HA2 24-184         => H3 353-513
          - Transmembrane Domain (TM): HA2 185-211 => H3 514-540
          - Cytoplasmic Tail (CT): HA2 212+ => H3 541+

        For non-H3 reference maps, falls back to coarse HA1/HA2.
        """

        numeric = self._coord_numeric(str(coord))
        if numeric is None:
            return "Unknown"

        if self.ref_subtype != "H3":
            return "HA1" if numeric < 330 else "HA2"

        if numeric < 1:
            return "Signal Peptide"

        # Cleavage loop spans end of HA1 and start of HA2 and may include insertions (e.g., 329a).
        if 327 <= numeric <= 330:
            return "Cleavage Site"

        if numeric <= 329:
            return "HA1"

        ha2_pos = numeric - 329  # H3 330 == HA2 1
        if 1 <= ha2_pos <= 23:
            return "Fusion Peptide"
        if 24 <= ha2_pos <= 184:
            return "HA2 Stalk"
        if 185 <= ha2_pos <= 211:
            return "Transmembrane Domain"
        return "Cytoplasmic Tail"

    @staticmethod
    def _confidence_value(conf_char: str) -> float:
        # HMMER's posterior string is symbolic confidence; map to a usable numeric scale.
        if conf_char in ("*", "9"):
            return 0.99
        if conf_char.isdigit():
            return max(0.0, min(1.0, int(conf_char) / 10.0))
        if conf_char in ("+",):
            return 0.85
        if conf_char in (".", "-", " "):
            return 0.0
        if conf_char.isalpha():
            return 0.75
        return 0.0

    def map_sequence(self, query_sequence: str) -> List[MappingEvent]:
        """Map one query sequence using alignment trace paths."""
        seq_text = self._clean_and_validate_sequence(query_sequence)
        digital = pyhmmer.easel.TextSequence(sequence=seq_text).digitize(self.alphabet)
        block = pyhmmer.easel.DigitalSequenceBlock(self.alphabet, [digital])
        hits = self.pipeline.search_hmm(self.hmm, block)

        best_hit = next(iter(hits), None)
        if best_hit is None or len(best_hit.domains) == 0:
            raise ValueError("No HA HMM hit found for query sequence.")

        domain = best_hit.domains[0]
        ali = domain.alignment
        hmm_path = ali.hmm_sequence
        target_path = ali.target_sequence
        conf_path = getattr(ali, "confidence", None)
        if conf_path is None:
            conf_path = getattr(ali, "posterior_probabilities", "")

        current_hmm_idx = ali.hmm_from
        current_target_idx = ali.target_from
        out: List[MappingEvent] = []
        prev_annotation = {"coord": "state_1", "site": "None", "domain": "Unknown"}

        for i, (h_char, t_char) in enumerate(zip(hmm_path, target_path)):
            conf = self._confidence_value(conf_path[i] if i < len(conf_path) else ".")
            if h_char != "-" and t_char != "-":
                annotation = self._annotation_for_state(current_hmm_idx)
                out.append(
                    MappingEvent(
                        query_pos=current_target_idx,
                        residue=t_char,
                        ref_coord=annotation["coord"],
                        site=annotation["site"],
                        domain=annotation["domain"],
                        event_type="Match",
                        hmm_state=current_hmm_idx,
                        confidence=conf,
                    )
                )
                prev_annotation = annotation
                current_hmm_idx += 1
                current_target_idx += 1
            elif h_char == "-" and t_char != "-":
                prev_state = max(1, current_hmm_idx - 1)
                inherited = prev_annotation
                if prev_state >= 1:
                    inherited = self._annotation_for_state(prev_state)
                    prev_annotation = inherited
                out.append(
                    MappingEvent(
                        query_pos=current_target_idx,
                        residue=t_char,
                        ref_coord=f"{inherited['coord']}+ins",
                        site=inherited["site"],
                        domain=inherited["domain"],
                        event_type="Insertion",
                        hmm_state=None,
                        confidence=conf,
                    )
                )
                current_target_idx += 1
            elif t_char == "-" and h_char != "-":
                annotation = self._annotation_for_state(current_hmm_idx)
                out.append(
                    MappingEvent(
                        query_pos=None,
                        residue="-",
                        ref_coord=annotation["coord"],
                        site=annotation["site"],
                        domain=annotation["domain"],
                        event_type="Deletion",
                        hmm_state=current_hmm_idx,
                        confidence=conf,
                    )
                )
                prev_annotation = annotation
                current_hmm_idx += 1

        if not out:
            raise ValueError("Alignment returned no trace events.")
        return out

    @staticmethod
    def _coord_numeric(coord: str) -> int | None:
        match = re.match(r"(\d+)", coord)
        if not match:
            return None
        return int(match.group(1))

    def detect_glycans(self, query_sequence: str, events: List[MappingEvent]) -> Dict[str, Any]:
        """Detect N-X-[S/T] motifs and annotate mapped N residues."""
        seq = self._clean_and_validate_sequence(query_sequence)
        motif_positions: List[int] = []
        for i in range(len(seq) - 2):
            if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in {"S", "T"}:
                motif_positions.append(i + 1)

        pos_to_event = {
            e.query_pos: e
            for e in events
            if e.query_pos is not None and e.residue.upper() == "N"
        }

        total_found = 0
        conserved = 0
        gained = 0
        detected_baseline: set[int] = set()

        for pos in motif_positions:
            event = pos_to_event.get(pos)
            if event is None:
                continue
            total_found += 1
            if self.ref_subtype == "H3":
                numeric = self._coord_numeric(event.ref_coord)
                if numeric is not None and numeric in self.H3_BASELINE_GLYCANS:
                    event.glycan = "N-Glyc (Conserved)"
                    conserved += 1
                    detected_baseline.add(numeric)
                else:
                    event.glycan = "N-Glyc (Gained)"
                    gained += 1
            else:
                event.glycan = "N-Glyc (Detected)"

        missing: List[int] = []
        if self.ref_subtype == "H3":
            missing = sorted(self.H3_BASELINE_GLYCANS - detected_baseline)

        return {
            "total_glycans_found": total_found,
            "conserved_glycans": conserved,
            "gained_glycans": gained,
            "missing_baseline_glycans": missing,
        }

    def analyze_cleavage_site(self, mapping_results: List[MappingEvent]) -> Dict[str, Any]:
        """Analyze HA1/HA2 cleavage loop for polybasic signatures."""
        loop_events: List[MappingEvent] = []
        for event in mapping_results:
            if event.event_type == "Deletion":
                continue
            numeric = self._coord_numeric(event.ref_coord)
            if numeric is None:
                continue
            # Capture residues around H3-327 to HA2 boundary including insertions.
            if 327 <= numeric <= 330:
                loop_events.append(event)

        cleavage_seq = "".join(e.residue for e in loop_events if e.residue != "-")
        for e in loop_events:
            e.region_label = "Cleavage_Site"

        basic_count = sum(1 for aa in cleavage_seq if aa in {"R", "K"})
        furin_found = re.search(r"R.[RK]R", cleavage_seq) is not None
        if basic_count >= 4 or furin_found:
            label = "HPAI (Polybasic)"
        else:
            label = "LPAI (Monobasic)"

        return {
            "cleavage_motif": cleavage_seq,
            "cleavage_basic_count": basic_count,
            "furin_motif_found": furin_found,
            "pathogenicity_score": label,
        }

    def build_metadata(self, query_id: str) -> Dict[str, str]:
        """Return shared metadata payload for JSON output."""
        return {
            "query_id": query_id,
            "date": date.today().isoformat(),
            "model_version": self.model_version,
            "reference_subtype": self.ref_subtype,
        }
