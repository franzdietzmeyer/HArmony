"""Typer CLI entrypoints for harmony."""

from __future__ import annotations

import time
import shutil
from pathlib import Path
from typing import Literal

import typer
from rich import box
from rich.align import Align
from rich.cells import cell_len
from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
)
from rich.text import Text

from .engine import HarmonyMapper
from .compare import compare_sequences, write_comparison_report
from .formatter import (
    events_to_rows,
    print_comparison_summary,
    print_run_summary,
    print_sequence_table,
    write_combined_csv,
    write_csv,
    write_json,
    write_master_summary,
)
from .io_handler import ensure_output_dir, load_sequences, safe_name
from .mapper import generate_maps
from .viz import generate_protein_map, plot_comparison
from .tracker import track_ancestor_against_descendants
from .visualizer import write_tracking_scripts

app = typer.Typer(
    help=(
        "harmony: Influenza A HA renumbering.\n\n"
        "Use `harmony map ...` explicitly, or run `harmony --input ... --output-dir ...` "
        "to execute the default map workflow."
    ),
    no_args_is_help=True,
)
console = Console()
WAIT_DURATION_SECONDS = 15.0


def display_welcome() -> None:
    """Render HArmony logo and subtitle."""
    logo_lines = [
        "██   ██  █████  ██████  ███    ███  ██████  ███    ██ ██    ██ ",
        "██   ██ ██   ██ ██   ██ ████  ████ ██    ██ ████   ██  ██  ██  ",
        "███████ ███████ ██████  ██ ████ ██ ██    ██ ██ ██  ██   ████   ",
        "██   ██ ██   ██ ██   ██ ██  ██  ██ ██    ██ ██  ██ ██    ██    ",
        "██   ██ ██   ██ ██   ██ ██      ██  ██████  ██   ████    ██    ",
    ]
    colors = [
        "cyan",
        "bright_cyan",
        "deep_sky_blue1",
        "dodger_blue1",
        "blue",
    ]
    text = Text(justify="left")
    logo_width = max(cell_len(line) for line in logo_lines)
    for line, color in zip(logo_lines, colors):
        pad = max(0, logo_width - cell_len(line))
        text.append(line + (" " * pad) + "\n", style=f"bold {color}")
    console.print(Align.center(text))
    console.print(
        Align.center(
            "[bold white]v1.1.0-alpha[/]"
        )
    )


def ensure_hmmer_available() -> None:
    """Fail fast if required HMMER binary is not present."""
    if shutil.which("hmmsearch") is None:
        raise RuntimeError(
            "HMMER not found. Please ensure you installed the package using the provided "
            "Conda environment: conda env create -f environment.yml"
        )


def print_citation_panel() -> None:
    """Render closing citation panel."""
    body = (
        "Thanks for using HArmony! If this tool helped your research, please cite:\n"
        "[italic]Dietzmeyer, F. (2026). HArmony: Universal HA Numbering and Antigenic Site "
        "Mapping. [Link/DOI][/]"
    )
    panel = Panel(
        body,
        title="[bold green]Success![/bold green]",
        border_style="bright_blue",
        box=box.DOUBLE,
    )
    console.print(panel)


def _run_map(
    input_path: Path,
    output_dir: Path,
    format: Literal["csv", "json", "both"],
    ref: str,
    combined: bool,
    plot: bool,
    quiet: bool,
    wait: bool,
) -> None:
    """Core map workflow shared by command and default callback."""
    if not quiet:
        display_welcome()
    try:
        output_dir = ensure_output_dir(output_dir)
        ref = ref.upper()
        ensure_hmmer_available()
        mapper = HarmonyMapper(ref_subtype=ref)
        records = load_sequences(input_path)
    except Exception as exc:
        console.print(f"[red]Initialization error:[/red] {exc}")
        raise typer.Exit(code=1)

    summary_entries: list[dict] = []
    combined_csv_rows: list[dict] = []
    combined_json_payload: dict[str, dict] = {}
    wait_sleep_per_sequence = WAIT_DURATION_SECONDS / max(len(records), 1) if wait else 0.0

    progress_cm = (
        Progress(
            TextColumn("[progress.description]{task.description}"),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            BarColumn(
                style="bright_black",
                complete_style="bright_blue",
                finished_style="blue",
                pulse_style="cyan",
            ),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
            console=console,
        )
        if not quiet
        else None
    )

    if progress_cm is not None:
        progress_ctx = progress_cm
    else:
        class _NoopContext:
            def __enter__(self):  # noqa: D401
                return None

            def __exit__(self, exc_type, exc, tb):  # noqa: D401
                return False

        progress_ctx = _NoopContext()

    with progress_ctx as progress:
        task_id = None
        if progress is not None:
            task_id = progress.add_task("Processing sequences", total=len(records))

        for rec in records:
            if progress is not None and task_id is not None:
                progress.update(task_id, description=f"Processing {rec.query_id}")
            try:
                events = mapper.map_sequence(rec.sequence)
            except Exception as exc:
                console.print(f"[yellow]Skipped {rec.query_id}:[/yellow] {exc}")
                if progress is not None and task_id is not None:
                    progress.advance(task_id)
                continue
            glycan_summary = mapper.detect_glycans(rec.sequence, events)
            cleavage_summary = mapper.analyze_cleavage_site(events)

            rows = events_to_rows(events, ref_subtype=ref)
            matches = sum(1 for e in events if e.event_type == "Match")
            insertions = sum(1 for e in events if e.event_type == "Insertion")
            deletions = sum(1 for e in events if e.event_type == "Deletion")
            mean_conf = sum(e.confidence for e in events) / len(events)
            site_counts = {
                "Site_A_Count": sum(1 for e in events if e.site == "Site A"),
                "Site_B_Count": sum(1 for e in events if e.site == "Site B"),
                "Site_C_Count": sum(1 for e in events if e.site == "Site C"),
                "Site_D_Count": sum(1 for e in events if e.site == "Site D"),
                "Site_E_Count": sum(1 for e in events if e.site == "Site E"),
                "Site_None_Count": sum(1 for e in events if e.site == "None"),
            }
            summary_entries.append(
                {
                    "Query_ID": rec.query_id,
                    "Source_File": str(rec.source_file),
                    "Total_Events": len(events),
                    "Matches": matches,
                    "Insertions": insertions,
                    "Deletions": deletions,
                    "Mean_Confidence": round(mean_conf, 4),
                    "Total_Glycans_Found": glycan_summary["total_glycans_found"],
                    "Pathogenicity_Score": cleavage_summary["pathogenicity_score"],
                    **site_counts,
                }
            )

            payload = {
                "metadata": {
                    **mapper.build_metadata(rec.query_id),
                    **glycan_summary,
                    **cleavage_summary,
                },
                "rows": rows,
            }

            if combined:
                for row in rows:
                    combined_csv_rows.append({"Query_ID": rec.query_id, **row})
                combined_json_payload[rec.query_id] = payload
            else:
                stem = safe_name(rec.query_id)
                seq_dir = ensure_output_dir(output_dir / stem)
                if format in ("csv", "both"):
                    write_csv(seq_dir / f"{stem}.csv", rows, ref_subtype=ref)
                if format in ("json", "both"):
                    write_json(seq_dir / f"{stem}.json", payload)

            if plot:
                fig_path = output_dir / f"{safe_name(rec.query_id)}_map.png"
                generate_protein_map(
                    mapping_df=rows,
                    output_path=fig_path,
                    sequence_id=rec.query_id,
                    ref_subtype=ref,
                    pathogenicity_score=cleavage_summary["pathogenicity_score"],
                )

            print_sequence_table(
                rec.query_id,
                rows,
                console,
                ref_subtype=ref,
                pathogenicity_score=cleavage_summary["pathogenicity_score"],
            )
            if wait:
                time.sleep(wait_sleep_per_sequence)
            if progress is not None and task_id is not None:
                progress.advance(task_id)

    if not summary_entries:
        console.print("[red]No sequences were successfully mapped.[/red]")
        raise typer.Exit(code=2)

    write_master_summary(output_dir / "master_summary.csv", summary_entries)

    if combined:
        if format in ("csv", "both"):
            write_combined_csv(output_dir / "detailed_mappings.csv", combined_csv_rows, ref_subtype=ref)
        if format in ("json", "both"):
            write_json(output_dir / "detailed_mappings.json", combined_json_payload)

    print_run_summary(summary_entries, console)
    console.print(f"[green]Done.[/green] Wrote outputs to [bold]{output_dir}[/bold]")
    if not quiet:
        print_citation_panel()


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    input: Path | None = typer.Option(
        None,
        "--input",
        exists=True,
        help="FASTA file or directory of FASTA files (default map mode).",
    ),
    test: bool = typer.Option(
        False,
        "--test",
        hidden=True,
        help="Dev tool: render UI mock output and exit.",
    ),
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        help="Directory to write outputs (default map mode).",
    ),
    format: Literal["csv", "json", "both"] = typer.Option(
        "both",
        "--format",
        help="Output format for detailed mapping files (default map mode).",
    ),
    ref: str = typer.Option(
        "H3",
        "--ref",
        help="Reference subtype (H1-H18) for default map mode.",
    ),
    combined: bool = typer.Option(
        False,
        "--combined/--split",
        help="Default map mode: write one combined file or split per-sequence files.",
    ),
    plot: bool = typer.Option(
        False,
        "--plot",
        help="Generate per-sequence linear protein visualization PNG maps.",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        help="Suppress logo and progress bar output.",
    ),
    wait: bool = typer.Option(
        False,
        "--wait",
        hidden=True,
        help="Dev tool: Adds a delay to test UI",
    ),
) -> None:
    """Default CLI behavior: run `map` when no subcommand is provided."""
    if ctx.invoked_subcommand is not None:
        return
    if test:
        display_welcome()
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            BarColumn(
                style="bright_black",
                complete_style="bright_blue",
                finished_style="blue",
                pulse_style="cyan",
            ),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
            console=console,
        ) as progress:
            steps = int(WAIT_DURATION_SECONDS / 0.5) if wait else 12
            task = progress.add_task("UI test: processing sequences", total=steps)
            for _ in range(steps):
                time.sleep(0.5 if wait else 0.06)
                progress.advance(task)

        console.print("[bold]Mock outputs[/bold]")
        console.print(f"- mutation_tracking_matrix.csv → [dim]./evolution_results/[/dim]")
        console.print(f"- mutation_tracking.pml → [dim]./evolution_results/[/dim]")
        console.print(f"- mutation_tracking.cxc → [dim]./evolution_results/[/dim]")
        console.print("")
        print_run_summary(
            [
                {
                    "Query_ID": "A/Vietnam/1203/2004",
                    "Source_File": "data/Vietnam.fasta",
                    "Total_Events": 565,
                    "Matches": 542,
                    "Insertions": 12,
                    "Deletions": 11,
                    "Mean_Confidence": 0.93,
                    "Total_Glycans_Found": 9,
                    "Pathogenicity_Score": "HPAI (Polybasic)",
                    "Site_A_Count": 3,
                    "Site_B_Count": 5,
                    "Site_C_Count": 1,
                    "Site_D_Count": 2,
                    "Site_E_Count": 0,
                    "Site_None_Count": 531,
                },
                {
                    "Query_ID": "Descendant_1",
                    "Source_File": "data/H15.fasta",
                    "Total_Events": 560,
                    "Matches": 545,
                    "Insertions": 6,
                    "Deletions": 9,
                    "Mean_Confidence": 0.91,
                    "Total_Glycans_Found": 8,
                    "Pathogenicity_Score": "LPAI (Monobasic)",
                    "Site_A_Count": 1,
                    "Site_B_Count": 2,
                    "Site_C_Count": 0,
                    "Site_D_Count": 1,
                    "Site_E_Count": 1,
                    "Site_None_Count": 541,
                },
            ],
            console,
        )
        console.print("")
        print_citation_panel()
        raise typer.Exit(code=0)
    if input is None or output_dir is None:
        raise typer.BadParameter(
            "Default mode requires --input and --output-dir. "
            "Equivalent explicit command: harmony map --input PATH --output-dir PATH --split --ref H3"
        )
    _run_map(
        input_path=input,
        output_dir=output_dir,
        format=format,
        ref=ref,
        combined=combined,
        plot=plot,
        quiet=quiet,
        wait=wait,
    )


@app.command("generate-maps")
def generate_maps_cmd(
    msa_fasta: Path = typer.Argument(..., exists=True, help="Aligned FASTA containing H1-H18 references."),
    out_dir: Path = typer.Option(
        Path("src/harmony/data/maps"),
        "--out-dir",
        help="Directory to write generated reference maps.",
    ),
) -> None:
    """Generate h1_map.json ... h18_map.json from aligned subtype FASTA."""
    try:
        written = generate_maps(msa_fasta=msa_fasta, output_dir=out_dir)
    except Exception as exc:
        console.print(f"[red]Map generation failed:[/red] {exc}")
        raise typer.Exit(code=1)
    console.print(f"[green]Generated {len(written)} maps[/green] in {out_dir}")


@app.command("map")
def map_cmd(
    input: Path = typer.Option(..., "--input", exists=True, help="FASTA file or directory of FASTA files."),
    output_dir: Path = typer.Option(..., "--output-dir", help="Directory to write outputs."),
    format: Literal["csv", "json", "both"] = typer.Option(
        "both", "--format", help="Output format for detailed mapping files."
    ),
    ref: str = typer.Option("H3", "--ref", help="The reference subtype to use for numbering (H1-H18)."),
    combined: bool = typer.Option(
        True,
        "--combined/--split",
        help="Write one combined detailed file or split per-sequence files.",
    ),
    plot: bool = typer.Option(
        False,
        "--plot",
        help="Generate per-sequence linear protein visualization PNG maps.",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        help="Suppress logo and progress bar output.",
    ),
    wait: bool = typer.Option(
        False,
        "--wait",
        hidden=True,
        help="Dev tool: Adds a delay to test UI",
    ),
) -> None:
    """Map all input sequences and write CSV/JSON outputs."""
    _run_map(
        input_path=input,
        output_dir=output_dir,
        format=format,
        ref=ref,
        combined=combined,
        plot=plot,
        quiet=quiet,
        wait=wait,
    )


@app.command("compare")
def compare_cmd(
    seq1: Path = typer.Option(..., "--seq1", exists=True, help="First input FASTA file."),
    seq2: Path = typer.Option(..., "--seq2", exists=True, help="Second input FASTA file."),
    ref: str = typer.Option("H3", "--ref", help="Common reference subtype coordinate system."),
    output_dir: Path = typer.Option(..., "--output-dir", help="Directory to write comparison outputs."),
    plot: bool = typer.Option(True, "--plot/--no-plot", help="Generate stacked comparison plot."),
    only_deltas: bool = typer.Option(
        False,
        "--only-deltas",
        help="Show only coordinates with sequence or glycan differences.",
    ),
    quiet: bool = typer.Option(
        False,
        "--quiet",
        help="Suppress logo and progress bar output.",
    ),
    wait: bool = typer.Option(
        False,
        "--wait",
        hidden=True,
        help="Dev tool: Adds a delay to test UI",
    ),
) -> None:
    """Compares two sequences on a common reference grid.

    By default, it returns the full sequence map. Use --only-deltas to see only mutations.
    """
    if not quiet:
        display_welcome()
    output_dir = ensure_output_dir(output_dir)
    try:
        ensure_hmmer_available()
    except Exception as exc:
        console.print(f"[red]Initialization error:[/red] {exc}")
        raise typer.Exit(code=1)
    if not quiet:
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            BarColumn(
                style="bright_black",
                complete_style="bright_blue",
                finished_style="blue",
                pulse_style="cyan",
            ),
            MofNCompleteColumn(),
            TimeRemainingColumn(),
            console=console,
        ) as progress:
            steps = int(WAIT_DURATION_SECONDS / 0.5) if wait else 1
            task = progress.add_task(f"Comparing {seq1.stem} vs {seq2.stem}", total=steps)
            try:
                if wait:
                    for _ in range(steps - 1):
                        time.sleep(0.5)
                        progress.advance(task)
                result = compare_sequences(seq1=seq1, seq2=seq2, ref=ref, only_deltas=only_deltas)
            except Exception as exc:
                console.print(f"[red]Comparison failed:[/red] {exc}")
                raise typer.Exit(code=1)
            progress.advance(task)
    else:
        try:
            result = compare_sequences(seq1=seq1, seq2=seq2, ref=ref, only_deltas=only_deltas)
        except Exception as exc:
            console.print(f"[red]Comparison failed:[/red] {exc}")
            raise typer.Exit(code=1)

    report_path = output_dir / "comparison_report.csv"
    write_comparison_report(report_path, result["differences"], ref_subtype=result["ref_subtype"])

    if plot:
        plot_comparison(
            df1=result["rows1"],
            df2=result["rows2"],
            output_path=output_dir / "comparison_plot.png",
            seq1_id=result["seq1_id"],
            seq2_id=result["seq2_id"],
            ref_subtype=result["ref_subtype"],
        )

    print_comparison_summary(result, console)
    console.print(f"[bold]Pathogenicity Shift:[/bold] {result['pathogenicity_shift']}")
    console.print(
        f"[bold]Cleavage Motifs:[/bold] "
        f"{result['seq1_id']}={result['cleavage_motif_seq1']} | "
        f"{result['seq2_id']}={result['cleavage_motif_seq2']}"
    )
    console.print(f"[green]Done.[/green] Wrote comparison outputs to [bold]{output_dir}[/bold]")
    if not quiet:
        print_citation_panel()


@app.command("track", context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def track_cmd(
    ctx: typer.Context,
    ancestor: Path = typer.Option(..., "--ancestor", exists=True, help="FASTA with a single ancestral sequence."),
    descendants: Path = typer.Option(
        ...,
        "--descendants",
        exists=True,
        help=(
            "Descendants input: either (1) one multi-FASTA file, (2) a directory of FASTA files, "
            "or (3) a single FASTA file followed by extra FASTA paths (e.g. shell glob expansion)."
        ),
    ),
    ref: str = typer.Option("H3", "--ref", help="Reference subtype coordinate system (default: H3)."),
    output_dir: Path = typer.Option(..., "--output-dir", help="Directory to write outputs."),
    quiet: bool = typer.Option(False, "--quiet", help="Suppress logo output."),
) -> None:
    """Track one ancestor against multiple descendants.

    Writes:
      - mutation_tracking_matrix.csv
      - mutation_tracking.pml (PyMOL)
      - mutation_tracking.cxc (ChimeraX)
    """
    if not quiet:
        display_welcome()
    output_dir = ensure_output_dir(output_dir)
    try:
        ensure_hmmer_available()
    except Exception as exc:
        console.print(f"[red]Initialization error:[/red] {exc}")
        raise typer.Exit(code=1)

    try:
        # Build delta-filtered mutation matrix and reuse mapping outputs for visualization.
        desc_path = descendants
        extra_desc = [Path(a) for a in getattr(ctx, "args", [])]
        if extra_desc:
            # Shell glob like data/*.fasta expands into extra args; merge into one multi-FASTA.
            from Bio import SeqIO

            descendant_files: list[Path] = [desc_path, *extra_desc]
            records = []
            for p in descendant_files:
                if not p.exists():
                    raise FileNotFoundError(f"Descendant file not found: {p}")
                if p.is_dir():
                    raise ValueError(
                        f"Descendants path {p} is a directory; pass it directly as --descendants {p}"
                    )
                records.extend(list(SeqIO.parse(str(p), "fasta")))
            if not records:
                raise ValueError("No descendant FASTA records found.")
            merged = output_dir / "descendants_merged.fasta"
            SeqIO.write(records, str(merged), "fasta")
            desc_path = merged

        df, anc_seq, anc_rows, _desc_names, desc_seqs, desc_rows = track_ancestor_against_descendants(
            ancestor_fasta=ancestor,
            descendants_fasta=desc_path,
            ref=ref,
        )
        out_matrix = output_dir / "mutation_tracking_matrix.csv"
        df.to_csv(out_matrix, index=False)
        ref_u = ref.upper()

        write_tracking_scripts(
            matrix_df=df,
            output_dir=output_dir,
            ref=ref_u,
            ancestor_sequence=anc_seq,
            ancestor_mapped_rows=anc_rows,
            descendant_sequences=desc_seqs,
            descendant_mapped_rows=desc_rows,
        )
    except Exception as exc:
        console.print(f"[red]Tracking failed:[/red] {exc}")
        raise typer.Exit(code=1)

    console.print(f"[green]Done.[/green] Wrote tracking outputs to [bold]{output_dir}[/bold]")
    if not quiet:
        print_citation_panel()


if __name__ == "__main__":
    app()
