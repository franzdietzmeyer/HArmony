"""Formatting and output writers for harmony mapping results."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List

from rich.console import Console
from rich.table import Table

from .engine import MappingEvent


def _coord_col(ref_subtype: str) -> str:
    return f"{ref_subtype.upper()}_Coord"


def _site_style(site: str) -> str:
    mapping = {
        "Site A": "[blue]Site A[/blue]",
        "Site B": "[bold red]Site B[/bold red]",
        "Site C": "[green]Site C[/green]",
        "Site D": "[yellow]Site D[/yellow]",
        "Site E": "[magenta]Site E[/magenta]",
        "None": "[dim]None[/dim]",
    }
    return mapping.get(site, f"[dim]{site}[/dim]")


def events_to_rows(events: Iterable[MappingEvent], ref_subtype: str) -> List[dict]:
    """Convert mapping events to standard row dictionaries."""
    coord_col = _coord_col(ref_subtype)
    rows = []
    for e in events:
        rows.append(
            {
                "Query_Pos": e.query_pos if e.query_pos is not None else "-",
                "Residue": e.residue,
                coord_col: e.ref_coord,
                "Site": e.site,
                "Domain": e.domain,
                "Glycan": e.glycan,
                "Region": e.region_label,
                "Type": e.event_type,
                "Confidence": round(float(e.confidence), 4),
            }
        )
    return rows


def write_csv(path: Path, rows: List[dict], ref_subtype: str) -> None:
    """Write a row list to CSV."""
    fieldnames = [
        "Query_Pos",
        "Residue",
        _coord_col(ref_subtype),
        "Site",
        "Domain",
        "Glycan",
        "Region",
        "Type",
        "Confidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_combined_csv(path: Path, rows: List[dict], ref_subtype: str) -> None:
    """Write combined rows from multiple sequences."""
    fieldnames = [
        "Query_ID",
        "Query_Pos",
        "Residue",
        _coord_col(ref_subtype),
        "Site",
        "Domain",
        "Glycan",
        "Region",
        "Type",
        "Confidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Dict) -> None:
    """Write JSON payload with UTF-8 encoding."""
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_master_summary(path: Path, entries: List[dict]) -> None:
    """Write summary CSV spanning all processed sequences."""
    if not entries:
        fieldnames = [
            "Query_ID",
            "Source_File",
            "Total_Events",
            "Matches",
            "Insertions",
            "Deletions",
            "Mean_Confidence",
            "Site_A_Count",
            "Site_B_Count",
            "Site_C_Count",
            "Site_D_Count",
            "Site_E_Count",
            "Site_None_Count",
            "Total_Glycans_Found",
            "Pathogenicity_Score",
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
        return
    fieldnames = [
        "Query_ID",
        "Source_File",
        "Total_Events",
        "Matches",
        "Insertions",
        "Deletions",
        "Mean_Confidence",
        "Site_A_Count",
        "Site_B_Count",
        "Site_C_Count",
        "Site_D_Count",
        "Site_E_Count",
        "Site_None_Count",
        "Total_Glycans_Found",
        "Pathogenicity_Score",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(entries)


def print_sequence_table(
    query_id: str,
    rows: List[dict],
    console: Console,
    ref_subtype: str,
    pathogenicity_score: str = "",
) -> None:
    """Print first events table with insertion/deletion highlighting."""
    coord_col = _coord_col(ref_subtype)
    if pathogenicity_score.startswith("HPAI"):
        console.print(f"[bold red]WARNING: {pathogenicity_score} detected[/bold red]")
    table = Table(title=f"Mapping preview: {query_id}", show_lines=False)
    table.add_column("Query_Pos", style="cyan")
    table.add_column("Residue", style="white")
    table.add_column(coord_col, style="green")
    table.add_column("Site")
    table.add_column("Glycan")
    table.add_column("Type")
    table.add_column("Confidence", justify="right", style="magenta")

    for row in rows[:20]:
        event = str(row["Type"])
        event_style = (
            "[green]Match[/green]"
            if event == "Match"
            else "[yellow]Insertion[/yellow]"
            if event == "Insertion"
            else "[red]Deletion[/red]"
        )
        table.add_row(
            str(row["Query_Pos"]),
            str(row["Residue"]),
            str(row[coord_col]),
            _site_style(str(row["Site"])),
            str(row.get("Glycan", "")),
            event_style,
            f"{float(row['Confidence']):.2f}",
        )
    console.print(table)


def print_run_summary(summary_entries: List[dict], console: Console) -> None:
    """Print compact run-level summary table."""
    table = Table(title="harmony run summary")
    table.add_column("Query ID", style="bold cyan")
    table.add_column("Matches", justify="right")
    table.add_column("Insertions", justify="right")
    table.add_column("Deletions", justify="right")
    table.add_column("Glycans", justify="right", style="cyan")
    table.add_column("Mean confidence", justify="right", style="magenta")

    for entry in summary_entries:
        table.add_row(
            str(entry["Query_ID"]),
            str(entry["Matches"]),
            str(entry["Insertions"]),
            str(entry["Deletions"]),
            str(entry.get("Total_Glycans_Found", 0)),
            f"{float(entry['Mean_Confidence']):.2f}",
        )
    console.print(table)


def print_comparison_summary(report: Dict, console: Console) -> None:
    """Print concise high-level summary for sequence comparison."""
    console.print(f"[bold]Total Mutations:[/bold] {report.get('total_mutations', 0)}")
    impacts = report.get("antigenic_site_impact", {})
    site_text = ", ".join(
        f"{site}:{count}" for site, count in impacts.items() if count > 0
    ) or "None"
    console.print(f"[bold]Antigenic Site Impact:[/bold] {site_text}")
    glycan_coords = report.get("glycan_change_coords", [])
    glycan_text = ", ".join(str(c) for c in glycan_coords) or "None"
    console.print(f"[bold]Glycan Change:[/bold] {glycan_text}")
