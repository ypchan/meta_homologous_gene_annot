#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
meta_homologous_gene_annot.py

Map reference proteins to one metagenomic contig assembly with miniprot,
filter and collapse redundant gene models, extract hit contigs, and export
contig-coordinate GFF3, standalone gene references, CDS, proteins, and
transcripts.

Required Python packages:
    rich
    rich-argparse

Required external programs:
    miniprot
    gffread
    pigz

Python >= 3.9 is recommended.
"""

from __future__ import annotations

import argparse
import csv
import errno
import gzip
import hashlib
import json
import os
import re
import resource
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import IO, Any, Iterable, Iterator, Sequence
from urllib.parse import quote, unquote

from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table
from rich.text import Text
from rich_argparse import RawDescriptionRichHelpFormatter


PROGRAM = "phi_contig_annotator"
PROGRAM_VERSION = "1.1.0"
GFF_SOURCE = "miniprot"
STAGE_ORDER = [
    "prepare_reference",
    "build_index",
    "map_proteins",
    "parse_hits",
    "extract_contigs",
    "export_sequences",
    "finalize",
]

INPUT_HELP = r"""
INPUT FORMAT
============

1. --proteins / -p
   Protein FASTA, optionally gzip-compressed (.gz/.bgz).

   Example:
       >PHI:1234 gene=MEP4 pathogen=Trichophyton_mentagrophytes
       MKFSLALALAVASASA...

   Rules:
   - The first non-whitespace token in each header is treated as original_id.
   - Duplicate original IDs are allowed; the script creates unique internal IDs.
   - Terminal '*' is removed.
   - Internal '*' and unsupported amino-acid symbols are replaced with X and reported.
   - Empty sequences are skipped and reported.

2. --contigs / -c
   One nucleotide FASTA assembly, optionally gzip-compressed (.gz/.bgz).

   Example:
       >contig_000001 length=18234 cov=13.5
       ACGTACGT...

   Rules:
   - The contig identifier is the first non-whitespace token after '>'.
   - Identifiers must be unique.
   - The FASTA may be very large. It is streamed and is never loaded wholly into RAM.
   - The same identifiers must appear in miniprot GFF3 and the FASTA.

3. Output directory
   All main outputs are written directly under --outdir. Only one hidden work
   directory is used temporarily. Existing compatible results can be resumed.

4. --organism_type
   Use "eukaryote" for splice-aware annotation (the default) or "prokaryote"
   to disable miniprot splicing. Eukaryotic GFF3 contains inferred introns;
   prokaryotic GFF3 contains exon/CDS features but no intron features.

Typical command
---------------
python3 meta_homologous_gene_annot.py \
    --proteins /share/cn1_fs/project/cyp_chenyanpeng/3domains/15_phibase/phi-base_current.fas \
    --contigs /share/data02/project/chenyanpeng/mangrove_2017_2025/01_contigs/201704_MF1.fasta.gz \
    --outdir 15_phibase/201704_MF1 \
    --sample 201704_MF1 \
    --organism_type eukaryote \
    --threads 24
"""

DEFAULT_HELP = r"""
DEFAULT PARAMETERS
==================

Alignment
---------
organism_type           eukaryote
threads                 available CPUs, capped at 32 when not specified
splice_model            1 for eukaryote; splicing disabled for prokaryote
max_intron               20000  maximum intron length in bp
index_subsample          1      miniprot -M; k-mer sampling rate is 1/2**M
max_hits                 50     retained/output alignments per query protein
secondary_ratio          0.50   retain secondary hits >= 50% of best score
prefilter_query_coverage 0.30   miniprot output prefilter
min_score_ratio          0.50   output hits >= 50% of best alignment score

Final filtering
---------------
min_identity             0.40   amino-acid identity fraction
min_query_coverage       0.60   aligned protein fraction
max_frameshift           1
max_stop_codon           0
locus_overlap            0.80   overlap/shorter interval used to merge references

High-confidence label
---------------------
high_identity            0.60
high_query_coverage      0.80
frameshift               must be 0
internal stop codon      must be 0

Execution
---------
resume                   enabled
keep_index               disabled after a successful run
keep_uncompressed        disabled
compression_threads      min(threads, 8)
log width                max(48, terminal_columns // 2)

Gene-reference coordinates
--------------------------
<sample>.gene.fasta      one genomic gene span per selected locus
<sample>.gene.gff3       the same models rebased to their locus sequence
                         (SeqID=locus_id; gene coordinates=1..gene_length)

Threshold interpretation
------------------------
These defaults are designed for sensitive homolog discovery, not species-level
identification. A PHI-base homolog does not by itself prove that the contig is
from the pathogen represented by the reference protein.
"""


class ShowPanelAction(argparse.Action):
    """Print a Rich panel and exit before required-argument validation."""

    def __init__(self, option_strings, dest, **kwargs):
        self.panel_text = kwargs.pop("panel_text")
        self.panel_title = kwargs.pop("panel_title")
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=0,
            default=argparse.SUPPRESS,
            **kwargs,
        )

    def __call__(self, parser, namespace, values, option_string=None):
        width = max(48, shutil.get_terminal_size(fallback=(120, 30)).columns // 2)
        console = Console(width=width)
        console.print(
            Panel.fit(
                self.panel_text.strip("\n"),
                title=self.panel_title,
                border_style="cyan",
                padding=(1, 2),
            )
        )
        parser.exit()


@dataclass
class StageMetric:
    stage: str
    status: str
    started_at: str
    finished_at: str
    wall_seconds: float
    cpu_seconds: float
    max_rss_mb: float
    detail: str


class RunLogger:
    """Rich console logger plus a plain-text persistent log."""

    def __init__(self, log_path: Path, width: int):
        self.log_path = log_path
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = open(log_path, "a", encoding="utf-8", buffering=1)
        self.console = Console(width=width, highlight=False)
        self.file_console = Console(
            file=self._handle,
            width=140,
            color_system=None,
            no_color=True,
            highlight=False,
        )

    def close(self) -> None:
        self._handle.close()

    @staticmethod
    def timestamp() -> str:
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def _write(self, level: str, message: str, style: str) -> None:
        prefix = f"[{self.timestamp()}] [{level}]"
        self.console.print(f"[{style}]{prefix}[/{style}] {message}")
        self.file_console.print(f"{prefix} {message}")

    def info(self, message: str) -> None:
        self._write("INFO", message, "bold cyan")

    def success(self, message: str) -> None:
        self._write("DONE", message, "bold green")

    def warning(self, message: str) -> None:
        self._write("WARN", message, "bold yellow")

    def error(self, message: str) -> None:
        self._write("ERROR", message, "bold red")

    def command(self, command: Sequence[str]) -> None:
        rendered = shlex.join([str(x) for x in command])
        self._write("CMD", rendered, "magenta")


class PipelineState:
    """JSON-backed stage state used for resumable execution."""

    def __init__(self, path: Path, signature: str):
        self.path = path
        self.signature = signature
        self.data: dict[str, Any] = {
            "program": PROGRAM,
            "program_version": PROGRAM_VERSION,
            "signature": signature,
            "created_at": iso_now(),
            "updated_at": iso_now(),
            "stages": {},
        }

        if path.exists():
            with open(path, "r", encoding="utf-8") as handle:
                loaded = json.load(handle)
            if loaded.get("signature") == signature:
                self.data = loaded

    def compatible(self) -> bool:
        return self.data.get("signature") == self.signature

    def is_done(self, stage: str, outputs: Sequence[Path]) -> bool:
        record = self.data.get("stages", {}).get(stage, {})
        if record.get("status") != "completed":
            return False
        saved = record.get("output_meta", {})
        for path in outputs:
            if not path.exists():
                return False
            expected = saved.get(str(path))
            if expected is not None and path.stat().st_size != expected.get("size"):
                return False
        return True

    def mark_completed(self, metric: StageMetric, outputs: Sequence[Path]) -> None:
        self.data.setdefault("stages", {})[metric.stage] = {
            **asdict(metric),
            "outputs": [str(path) for path in outputs],
            "output_meta": {
                str(path): {"size": path.stat().st_size}
                for path in outputs
                if path.exists()
            },
        }
        self.data["updated_at"] = iso_now()
        self.save()

    def reset(self) -> None:
        self.data = {
            "program": PROGRAM,
            "program_version": PROGRAM_VERSION,
            "signature": self.signature,
            "created_at": iso_now(),
            "updated_at": iso_now(),
            "stages": {},
        }
        self.save()

    def save(self) -> None:
        atomic_write_json(self.path, self.data)


class StageTimer:
    def __init__(self, stage: str):
        self.stage = stage
        self.started_at = iso_now()
        self.start_wall = time.perf_counter()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        self.start_cpu = usage.ru_utime + usage.ru_stime

    def finish(self, status: str, detail: str = "") -> StageMetric:
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_now = usage.ru_utime + usage.ru_stime
        max_rss_mb = usage.ru_maxrss / 1024.0
        return StageMetric(
            stage=self.stage,
            status=status,
            started_at=self.started_at,
            finished_at=iso_now(),
            wall_seconds=round(time.perf_counter() - self.start_wall, 3),
            cpu_seconds=round(max(0.0, cpu_now - self.start_cpu), 3),
            max_rss_mb=round(max_rss_mb, 3),
            detail=detail,
        )


def iso_now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with open(temp, "w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temp, path)


def move_replace(source: Path, destination: Path) -> None:
    """Replace destination safely, including across filesystems."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.replace(source, destination)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        incoming = destination.with_name(destination.name + ".incoming")
        incoming.unlink(missing_ok=True)
        shutil.copy2(source, incoming)
        os.replace(incoming, destination)
        source.unlink()


def open_text_auto(path: Path, mode: str = "rt") -> IO[str]:
    if path.suffix.lower() in {".gz", ".bgz", ".bgzf"}:
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with open_text_auto(path, "rt") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def write_wrapped(handle: IO[str], sequence: str, width: int = 80) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def file_signature(path: Path, content_hash: bool = False) -> dict[str, Any]:
    resolved = path.resolve()
    stat = resolved.stat()
    result: dict[str, Any] = {
        "path": str(resolved),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if content_hash:
        digest = hashlib.sha256()
        with open(resolved, "rb") as handle:
            while block := handle.read(8 * 1024 * 1024):
                digest.update(block)
        result["sha256"] = digest.hexdigest()
    return result


def stable_hash(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def infer_sample_name(contigs: Path) -> str:
    name = contigs.name
    for suffix in (".gz", ".bgz", ".bgzf", ".fasta", ".fna", ".fa", ".fas"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "dataset"


def human_bytes(value: int | float) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    size = float(value)
    for unit in units:
        if abs(size) < 1024.0 or unit == units[-1]:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} TiB"


def human_seconds(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m {sec:.1f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h {int(minutes)}m {sec:.1f}s"


def detect_threads() -> int:
    for key in ("SLURM_CPUS_PER_TASK", "NSLOTS", "OMP_NUM_THREADS"):
        value = os.environ.get(key)
        if value and value.isdigit() and int(value) > 0:
            return int(value)
    return min(os.cpu_count() or 1, 32)


def miniprot_annotation_options(args: argparse.Namespace) -> list[str]:
    if args.organism_type == "prokaryote":
        return ["-S"]
    return ["-j", str(args.splice_model), "-G", str(args.max_intron)]


def require_executable(name_or_path: str) -> str:
    if os.path.sep in name_or_path:
        path = Path(name_or_path).expanduser().resolve()
        if not path.exists() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f"Executable is missing or not executable: {path}")
        return str(path)
    found = shutil.which(name_or_path)
    if found is None:
        raise FileNotFoundError(f"Executable not found in PATH: {name_or_path}")
    return found


def capture_tool_version(command: str) -> str:
    probes = ([command, "--version"], [command, "-V"], [command, "-h"])
    for probe in probes:
        result = subprocess.run(
            probe,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        text = result.stdout.strip()
        if text:
            return text.splitlines()[0][:300]
    return "unknown"


def tail_text(path: Path, lines: int = 30) -> str:
    if not path.exists():
        return ""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        content = handle.readlines()
    return "".join(content[-lines:]).rstrip()


def run_external(
    command: Sequence[str],
    logger: RunLogger,
    stdout_path: Path | None,
    stderr_path: Path,
) -> None:
    logger.command(command)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    stdout_handle: IO[Any] | int
    if stdout_path is None:
        stdout_handle = subprocess.DEVNULL
    else:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        stdout_handle = open(stdout_path, "wb")

    with open(stderr_path, "wb") as stderr_handle:
        try:
            process = subprocess.run(
                [str(x) for x in command],
                stdout=stdout_handle,
                stderr=stderr_handle,
                check=False,
            )
        finally:
            if stdout_path is not None and hasattr(stdout_handle, "close"):
                stdout_handle.close()

    if process.returncode != 0:
        excerpt = tail_text(stderr_path)
        message = f"Command failed with exit code {process.returncode}: {shlex.join(command)}"
        if excerpt:
            message += f"\nLast stderr lines:\n{excerpt}"
        raise RuntimeError(message)


def prepare_reference(
    input_fasta: Path,
    clean_fasta: Path,
    map_tsv: Path,
) -> dict[str, Any]:
    clean_tmp = clean_fasta.with_name(clean_fasta.name + ".tmp")
    map_tmp = map_tsv.with_name(map_tsv.name + ".tmp")
    allowed = set("ABCDEFGHIKLMNPQRSTVWXYZUO")

    n_input = 0
    n_written = 0
    n_empty = 0
    n_replaced_records = 0
    n_replaced_residues = 0
    total_aa = 0

    with open(clean_tmp, "w", encoding="utf-8") as fasta_out, open(
        map_tmp, "w", encoding="utf-8", newline=""
    ) as map_out:
        writer = csv.writer(map_out, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "query_id",
                "original_id",
                "original_header",
                "protein_length",
                "replaced_residues",
            ]
        )

        for header, raw_sequence in read_fasta(input_fasta):
            n_input += 1
            sequence = re.sub(r"\s+", "", raw_sequence).upper().replace("-", "").replace(".", "")
            sequence = sequence.rstrip("*")
            if not sequence:
                n_empty += 1
                continue

            replaced = 0
            cleaned: list[str] = []
            for residue in sequence:
                if residue in allowed:
                    cleaned.append(residue)
                else:
                    cleaned.append("X")
                    replaced += 1
            clean_sequence = "".join(cleaned)

            n_written += 1
            query_id = f"PHIREF{n_written:09d}"
            original_id = header.split(maxsplit=1)[0] if header else query_id
            original_header = re.sub(r"[\t\r\n]+", " ", header).strip()

            fasta_out.write(f">{query_id}\n")
            write_wrapped(fasta_out, clean_sequence)
            writer.writerow(
                [query_id, original_id, original_header, len(clean_sequence), replaced]
            )

            total_aa += len(clean_sequence)
            n_replaced_residues += replaced
            if replaced:
                n_replaced_records += 1

    if n_written == 0:
        clean_tmp.unlink(missing_ok=True)
        map_tmp.unlink(missing_ok=True)
        raise ValueError("No usable protein sequences were found in the input FASTA")

    os.replace(clean_tmp, clean_fasta)
    os.replace(map_tmp, map_tsv)
    return {
        "input_records": n_input,
        "written_records": n_written,
        "empty_records": n_empty,
        "records_with_replacements": n_replaced_records,
        "replaced_residues": n_replaced_residues,
        "total_amino_acids": total_aa,
    }


def load_reference_map(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            records[row["query_id"]] = row
    return records


def parse_gff_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        result[key] = unquote(value)
    return result


def format_gff_attributes(attributes: dict[str, str]) -> str:
    return ";".join(
        f"{key}={quote(str(value), safe=':@|,._+-')}" for key, value in attributes.items()
    )


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def overlap_shorter(a: dict[str, Any], b: dict[str, Any]) -> float:
    overlap = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
    if overlap == 0:
        return 0.0
    shorter = min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)
    return overlap / shorter if shorter > 0 else 0.0


def best_hit_key(hit: dict[str, Any]) -> tuple[Any, ...]:
    confidence_rank = {"low": 0, "medium": 1, "high": 2}
    return (
        confidence_rank.get(hit["confidence"], 0),
        hit["identity"] * hit["query_coverage"],
        hit["query_coverage"],
        hit["identity"],
        hit["score"],
        -hit["frameshift"],
        -hit["stop_codon"],
        -hit["rank"],
    )


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], columns: Sequence[str]) -> None:
    temp = path.with_name(path.name + ".tmp")
    with open(temp, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(temp, path)


def load_miniprot_models(raw_gff: Path) -> dict[str, dict[str, Any]]:
    """Load miniprot feature blocks and their embedded detailed PAF lines."""
    models: dict[str, dict[str, Any]] = {}
    pending_alignment: list[str] = []

    with open(raw_gff, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("##PAF\t"):
                pending_alignment = [line]
                continue
            if line.startswith(("##ATN\t", "##ATA\t", "##AAS\t", "##AQA\t", "##STA\t")):
                pending_alignment.append(line)
                continue
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue
            attrs = parse_gff_attributes(fields[8])
            if fields[2] == "mRNA":
                model_id = attrs.get("ID", "")
                if not model_id:
                    pending_alignment = []
                    continue
                models[model_id] = {
                    "mrna": fields,
                    "children": [],
                    "alignment": pending_alignment,
                }
                pending_alignment = []
                continue

            for parent in (item for item in attrs.get("Parent", "").split(",") if item):
                if parent in models:
                    models[parent]["children"].append(fields)

    return models


def merge_feature_intervals(
    features: Sequence[list[str]], fallback: tuple[int, int]
) -> list[tuple[int, int]]:
    """Merge CDS/stop-codon intervals into exon intervals."""
    intervals = sorted(
        (safe_int(fields[3]), safe_int(fields[4]))
        for fields in features
        if fields[2] in {"CDS", "stop_codon"}
        and safe_int(fields[3]) > 0
        and safe_int(fields[4]) >= safe_int(fields[3])
    )
    if not intervals:
        return [fallback]

    merged: list[tuple[int, int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1] + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def write_selected_gff_files(
    raw_gff: Path,
    best_loci: Sequence[dict[str, Any]],
    best_gff: Path,
    gene_gff: Path,
    organism_type: str,
) -> dict[str, int]:
    """Write selected models in contig and standalone-gene coordinate systems."""
    models = load_miniprot_models(raw_gff)
    best_gff_tmp = best_gff.with_name(best_gff.name + ".tmp")
    gene_gff_tmp = gene_gff.with_name(gene_gff.name + ".tmp")
    exon_total = 0
    intron_total = 0

    with open(best_gff_tmp, "w", encoding="utf-8") as contig_out, open(
        gene_gff_tmp, "w", encoding="utf-8"
    ) as gene_out:
        contig_out.write("##gff-version 3\n")
        gene_out.write("##gff-version 3\n")

        for selected in best_loci:
            model = models.get(str(selected["model_id"]))
            if model is None:
                raise RuntimeError(
                    f"Selected miniprot model is absent from raw GFF3: {selected['model_id']}"
                )

            mrna = model["mrna"]
            children: list[list[str]] = model["children"]
            locus_id = str(selected["locus_id"])
            transcript_id = f"{locus_id}.t1"
            gene_start = safe_int(selected["start"])
            gene_end = safe_int(selected["end"])
            gene_length = gene_end - gene_start + 1
            strand = str(selected["strand"])
            reference_id = str(selected["original_id"] or selected["query_id"])
            source_attrs = {
                "gene_id": locus_id,
                "Name": locus_id,
                "SampleID": str(selected["sample_id"]),
                "ReferenceID": reference_id,
                "ReferenceHitCount": str(selected["n_reference_hits"]),
                "Confidence": str(selected["confidence"]),
                "QueryCoverage": f"{selected['query_coverage']:.6f}",
                "Identity": f"{selected['identity']:.6f}",
                "ModelID": str(selected["model_id"]),
                "SourceContig": str(selected["contig"]),
                "SourceStart": str(gene_start),
                "SourceEnd": str(gene_end),
                "SourceStrand": strand,
                "OrganismType": organism_type,
                "ReferenceAnnotation": str(selected["reference_annotation"])[:2000],
            }
            if selected["alternative_references"]:
                source_attrs["AlternativeReferences"] = str(
                    selected["alternative_references"]
                )
            exon_intervals = merge_feature_intervals(children, (gene_start, gene_end))
            exon_total += len(exon_intervals)

            exon_transcript_order = (
                exon_intervals if strand != "-" else list(reversed(exon_intervals))
            )
            exon_numbers = {
                interval: index for index, interval in enumerate(exon_transcript_order, 1)
            }

            introns: list[tuple[int, int, int]] = []
            if organism_type == "eukaryote":
                for left, right in zip(exon_intervals, exon_intervals[1:]):
                    intron_start = left[1] + 1
                    intron_end = right[0] - 1
                    if intron_start <= intron_end:
                        intron_number = min(exon_numbers[left], exon_numbers[right])
                        introns.append((intron_start, intron_end, intron_number))
            intron_total += len(introns)

            for coordinate_mode, output in (("contig", contig_out), ("gene", gene_out)):
                if coordinate_mode == "contig":
                    seqid = str(selected["contig"])
                    offset = 0
                    for alignment_line in model["alignment"]:
                        output.write(alignment_line + "\n")
                else:
                    seqid = locus_id
                    offset = gene_start - 1
                    output.write(f"##sequence-region {locus_id} 1 {gene_length}\n")

                def convert(start: int, end: int) -> tuple[int, int]:
                    return start - offset, end - offset

                local_gene_start, local_gene_end = convert(gene_start, gene_end)
                gene_fields = [
                    seqid,
                    PROGRAM,
                    "gene",
                    str(local_gene_start),
                    str(local_gene_end),
                    mrna[5],
                    strand,
                    ".",
                    format_gff_attributes({"ID": locus_id, **source_attrs}),
                ]
                output.write("\t".join(gene_fields) + "\n")

                raw_mrna_attrs = parse_gff_attributes(mrna[8])
                raw_mrna_attrs.pop("ID", None)
                transcript_attrs = {
                    "ID": transcript_id,
                    "Parent": locus_id,
                    "gene_id": locus_id,
                    "transcript_id": transcript_id,
                    **raw_mrna_attrs,
                    "ModelID": str(selected["model_id"]),
                    "ReferenceID": reference_id,
                    "Confidence": str(selected["confidence"]),
                    "QueryCoverage": f"{selected['query_coverage']:.6f}",
                }
                transcript_fields = [
                    seqid,
                    mrna[1],
                    "mRNA",
                    str(local_gene_start),
                    str(local_gene_end),
                    mrna[5],
                    strand,
                    ".",
                    format_gff_attributes(transcript_attrs),
                ]
                output.write("\t".join(transcript_fields) + "\n")

                derived_records: list[tuple[int, int, int, list[str]]] = []
                for exon_start, exon_end in exon_intervals:
                    exon_number = exon_numbers[(exon_start, exon_end)]
                    local_start, local_end = convert(exon_start, exon_end)
                    fields = [
                        seqid,
                        PROGRAM,
                        "exon",
                        str(local_start),
                        str(local_end),
                        ".",
                        strand,
                        ".",
                        format_gff_attributes(
                            {
                                "ID": f"{transcript_id}.exon{exon_number}",
                                "Parent": transcript_id,
                                "gene_id": locus_id,
                                "transcript_id": transcript_id,
                                "exon_number": str(exon_number),
                            }
                        ),
                    ]
                    derived_records.append((local_start, local_end, 0, fields))

                cds_children = [fields for fields in children if fields[2] == "CDS"]
                cds_order = sorted(
                    cds_children,
                    key=lambda fields: (safe_int(fields[3]), safe_int(fields[4])),
                    reverse=strand == "-",
                )
                cds_numbers = {id(fields): index for index, fields in enumerate(cds_order, 1)}
                for child in cds_children:
                    cds_number = cds_numbers[id(child)]
                    local_start, local_end = convert(safe_int(child[3]), safe_int(child[4]))
                    child_attrs = parse_gff_attributes(child[8])
                    child_attrs.pop("Parent", None)
                    child_attrs.pop("ID", None)
                    fields = [
                        seqid,
                        child[1],
                        "CDS",
                        str(local_start),
                        str(local_end),
                        child[5],
                        strand,
                        child[7],
                        format_gff_attributes(
                            {
                                "ID": f"{transcript_id}.cds{cds_number}",
                                "Parent": transcript_id,
                                "gene_id": locus_id,
                                "transcript_id": transcript_id,
                                **child_attrs,
                            }
                        ),
                    ]
                    derived_records.append((local_start, local_end, 1, fields))

                stop_children = [fields for fields in children if fields[2] == "stop_codon"]
                stop_order = sorted(
                    stop_children,
                    key=lambda fields: (safe_int(fields[3]), safe_int(fields[4])),
                    reverse=strand == "-",
                )
                stop_numbers = {id(fields): index for index, fields in enumerate(stop_order, 1)}
                for child in stop_children:
                    stop_number = stop_numbers[id(child)]
                    local_start, local_end = convert(safe_int(child[3]), safe_int(child[4]))
                    child_attrs = parse_gff_attributes(child[8])
                    child_attrs.pop("Parent", None)
                    child_attrs.pop("ID", None)
                    fields = [
                        seqid,
                        child[1],
                        "stop_codon",
                        str(local_start),
                        str(local_end),
                        child[5],
                        strand,
                        child[7],
                        format_gff_attributes(
                            {
                                "ID": f"{transcript_id}.stop{stop_number}",
                                "Parent": transcript_id,
                                "gene_id": locus_id,
                                "transcript_id": transcript_id,
                                **child_attrs,
                            }
                        ),
                    ]
                    derived_records.append((local_start, local_end, 2, fields))

                for intron_start, intron_end, intron_number in introns:
                    local_start, local_end = convert(intron_start, intron_end)
                    fields = [
                        seqid,
                        PROGRAM,
                        "intron",
                        str(local_start),
                        str(local_end),
                        ".",
                        strand,
                        ".",
                        format_gff_attributes(
                            {
                                "ID": f"{transcript_id}.intron{intron_number}",
                                "Parent": transcript_id,
                                "gene_id": locus_id,
                                "transcript_id": transcript_id,
                                "intron_number": str(intron_number),
                            }
                        ),
                    ]
                    derived_records.append((local_start, local_end, 3, fields))

                for _, _, _, fields in sorted(
                    derived_records, key=lambda record: (record[0], record[1], record[2])
                ):
                    output.write("\t".join(fields) + "\n")

    os.replace(best_gff_tmp, best_gff)
    os.replace(gene_gff_tmp, gene_gff)
    return {"exons": exon_total, "introns": intron_total}


def parse_and_filter_hits(
    raw_gff: Path,
    reference_map_tsv: Path,
    sample: str,
    args: argparse.Namespace,
    all_hits_tsv: Path,
    best_loci_tsv: Path,
    best_gff: Path,
    gene_gff: Path,
    hit_contig_ids: Path,
    query_summary_tsv: Path,
    contig_summary_tsv: Path,
) -> dict[str, Any]:
    reference = load_reference_map(reference_map_tsv)
    hits: list[dict[str, Any]] = []

    with open(raw_gff, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue

            attrs = parse_gff_attributes(fields[8])
            model_id = attrs.get("ID", "")
            if not model_id:
                continue

            target_fields = attrs.get("Target", "").split()
            query_id = target_fields[0] if target_fields else model_id.split("@", 1)[0]
            metadata = reference.get(query_id, {})
            query_length = safe_int(metadata.get("protein_length"), 0)
            query_start = safe_int(target_fields[1], 0) if len(target_fields) >= 3 else 0
            query_end = safe_int(target_fields[2], 0) if len(target_fields) >= 3 else 0
            aligned_query_length = (
                abs(query_end - query_start) + 1
                if query_start > 0 and query_end > 0
                else 0
            )
            query_coverage = aligned_query_length / query_length if query_length else 0.0
            identity = safe_float(attrs.get("Identity"), 0.0)
            positive = safe_float(attrs.get("Positive"), 0.0)
            frameshift = safe_int(attrs.get("Frameshift"), 0)
            stop_codon = safe_int(attrs.get("StopCodon"), 0)
            rank = safe_int(attrs.get("Rank"), 1)
            score = safe_float(fields[5], 0.0)

            fail_reasons: list[str] = []
            if identity < args.min_identity:
                fail_reasons.append("low_identity")
            if query_coverage < args.min_query_coverage:
                fail_reasons.append("low_query_coverage")
            if frameshift > args.max_frameshift:
                fail_reasons.append("too_many_frameshifts")
            if stop_codon > args.max_stop_codon:
                fail_reasons.append("internal_stop_codon")
            pass_filter = not fail_reasons

            if (
                identity >= args.high_identity
                and query_coverage >= args.high_query_coverage
                and frameshift == 0
                and stop_codon == 0
            ):
                confidence = "high"
            elif pass_filter:
                confidence = "medium"
            else:
                confidence = "low"

            hits.append(
                {
                    "sample_id": sample,
                    "organism_type": args.organism_type,
                    "model_id": model_id,
                    "query_id": query_id,
                    "original_id": metadata.get("original_id", ""),
                    "reference_annotation": metadata.get("original_header", ""),
                    "contig": fields[0],
                    "start": safe_int(fields[3]),
                    "end": safe_int(fields[4]),
                    "strand": fields[6],
                    "genomic_span": safe_int(fields[4]) - safe_int(fields[3]) + 1,
                    "score": score,
                    "rank": rank,
                    "identity": identity,
                    "positive": positive,
                    "query_start": query_start,
                    "query_end": query_end,
                    "aligned_query_length": aligned_query_length,
                    "query_length": query_length,
                    "query_coverage": query_coverage,
                    "frameshift": frameshift,
                    "stop_codon": stop_codon,
                    "confidence": confidence,
                    "pass_filter": int(pass_filter),
                    "fail_reason": ";".join(fail_reasons),
                }
            )

    all_columns = [
        "sample_id",
        "organism_type",
        "model_id",
        "query_id",
        "original_id",
        "reference_annotation",
        "contig",
        "start",
        "end",
        "strand",
        "genomic_span",
        "score",
        "rank",
        "identity",
        "positive",
        "query_start",
        "query_end",
        "aligned_query_length",
        "query_length",
        "query_coverage",
        "frameshift",
        "stop_codon",
        "confidence",
        "pass_filter",
        "fail_reason",
    ]
    write_tsv(all_hits_tsv, hits, all_columns)

    filtered = [hit for hit in hits if hit["pass_filter"] == 1]
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for hit in filtered:
        grouped[(hit["contig"], hit["strand"])].append(hit)

    clusters: list[list[dict[str, Any]]] = []
    for group_hits in grouped.values():
        group_hits.sort(key=lambda item: (item["start"], item["end"]))
        current: list[dict[str, Any]] = []
        current_end = -1
        for hit in group_hits:
            if not current:
                current = [hit]
                current_end = hit["end"]
                continue
            belongs = hit["start"] <= current_end and any(
                overlap_shorter(hit, member) >= args.locus_overlap for member in current
            )
            if belongs:
                current.append(hit)
                current_end = max(current_end, hit["end"])
            else:
                clusters.append(current)
                current = [hit]
                current_end = hit["end"]
        if current:
            clusters.append(current)

    clusters.sort(
        key=lambda cluster: (
            cluster[0]["contig"],
            min(item["start"] for item in cluster),
            max(item["end"] for item in cluster),
            cluster[0]["strand"],
        )
    )

    best_loci: list[dict[str, Any]] = []
    for index, cluster in enumerate(clusters, start=1):
        best = max(cluster, key=best_hit_key).copy()
        alternatives = sorted(
            {
                item["original_id"] or item["query_id"]
                for item in cluster
                if item["model_id"] != best["model_id"]
            }
        )
        best["locus_id"] = f"{sample}_PHILOCUS{index:07d}"
        best["locus_start"] = min(item["start"] for item in cluster)
        best["locus_end"] = max(item["end"] for item in cluster)
        best["n_reference_hits"] = len(cluster)
        best["alternative_references"] = ";".join(alternatives)
        best_loci.append(best)

    best_columns = [
        "sample_id",
        "organism_type",
        "locus_id",
        "model_id",
        "query_id",
        "original_id",
        "reference_annotation",
        "contig",
        "start",
        "end",
        "locus_start",
        "locus_end",
        "strand",
        "genomic_span",
        "score",
        "rank",
        "identity",
        "positive",
        "query_start",
        "query_end",
        "aligned_query_length",
        "query_length",
        "query_coverage",
        "frameshift",
        "stop_codon",
        "confidence",
        "n_reference_hits",
        "alternative_references",
    ]
    write_tsv(best_loci_tsv, best_loci, best_columns)

    feature_stats = write_selected_gff_files(
        raw_gff=raw_gff,
        best_loci=best_loci,
        best_gff=best_gff,
        gene_gff=gene_gff,
        organism_type=args.organism_type,
    )

    hit_contigs = sorted({item["contig"] for item in best_loci})
    temp_ids = hit_contig_ids.with_name(hit_contig_ids.name + ".tmp")
    with open(temp_ids, "w", encoding="utf-8") as handle:
        for contig in hit_contigs:
            handle.write(contig + "\n")
    os.replace(temp_ids, hit_contig_ids)

    raw_by_query = Counter(hit["query_id"] for hit in hits)
    pass_by_query = Counter(hit["query_id"] for hit in filtered)
    selected_by_query = Counter(hit["query_id"] for hit in best_loci)
    best_for_query: dict[str, dict[str, Any]] = {}
    for hit in hits:
        current = best_for_query.get(hit["query_id"])
        if current is None or best_hit_key(hit) > best_hit_key(current):
            best_for_query[hit["query_id"]] = hit

    query_rows: list[dict[str, Any]] = []
    for query_id, metadata in reference.items():
        best = best_for_query.get(query_id, {})
        query_rows.append(
            {
                "query_id": query_id,
                "original_id": metadata.get("original_id", ""),
                "reference_annotation": metadata.get("original_header", ""),
                "protein_length": metadata.get("protein_length", ""),
                "raw_hits": raw_by_query[query_id],
                "passing_hits": pass_by_query[query_id],
                "selected_loci": selected_by_query[query_id],
                "best_identity": best.get("identity", ""),
                "best_query_coverage": best.get("query_coverage", ""),
                "best_contig": best.get("contig", ""),
                "best_confidence": best.get("confidence", "unmapped"),
            }
        )
    query_columns = [
        "query_id",
        "original_id",
        "reference_annotation",
        "protein_length",
        "raw_hits",
        "passing_hits",
        "selected_loci",
        "best_identity",
        "best_query_coverage",
        "best_contig",
        "best_confidence",
    ]
    write_tsv(query_summary_tsv, query_rows, query_columns)

    contig_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for hit in best_loci:
        contig_groups[hit["contig"]].append(hit)
    contig_rows: list[dict[str, Any]] = []
    for contig, records in sorted(contig_groups.items()):
        contig_rows.append(
            {
                "sample_id": sample,
                "contig": contig,
                "n_loci": len(records),
                "n_high": sum(item["confidence"] == "high" for item in records),
                "n_medium": sum(item["confidence"] == "medium" for item in records),
                "best_identity": max(item["identity"] for item in records),
                "best_query_coverage": max(item["query_coverage"] for item in records),
                "reference_ids": ";".join(
                    sorted({item["original_id"] or item["query_id"] for item in records})
                ),
            }
        )
    contig_columns = [
        "sample_id",
        "contig",
        "n_loci",
        "n_high",
        "n_medium",
        "best_identity",
        "best_query_coverage",
        "reference_ids",
    ]
    write_tsv(contig_summary_tsv, contig_rows, contig_columns)

    return {
        "raw_hits": len(hits),
        "passing_hits": len(filtered),
        "selected_loci": len(best_loci),
        "high_confidence_loci": sum(item["confidence"] == "high" for item in best_loci),
        "medium_confidence_loci": sum(item["confidence"] == "medium" for item in best_loci),
        "hit_contigs": len(hit_contigs),
        "mapped_queries": len(raw_by_query),
        "passing_queries": len(pass_by_query),
        "total_queries": len(reference),
        **feature_stats,
    }


def extract_hit_contigs(
    contigs_fasta: Path,
    wanted_ids_path: Path,
    output_fasta: Path,
    best_loci_tsv: Path,
    gene_fasta: Path,
    organism_type: str,
) -> dict[str, Any]:
    with open(wanted_ids_path, "r", encoding="utf-8") as handle:
        wanted = {line.strip() for line in handle if line.strip()}

    loci_by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open(best_loci_tsv, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            loci_by_contig[row["contig"]].append(row)
    for records in loci_by_contig.values():
        records.sort(
            key=lambda row: (
                safe_int(row["start"]),
                safe_int(row["end"]),
                row["locus_id"],
            )
        )

    output_tmp = output_fasta.with_name(output_fasta.name + ".tmp")
    gene_tmp = gene_fasta.with_name(gene_fasta.name + ".tmp")
    found: set[str] = set()
    duplicate_ids: set[str] = set()
    scanned = 0
    written_bases = 0
    written_genes = 0
    written_gene_bases = 0

    with open(output_tmp, "w", encoding="utf-8") as output, open(
        gene_tmp, "w", encoding="utf-8"
    ) as gene_output:
        for header, sequence in read_fasta(contigs_fasta):
            scanned += 1
            seq_id = header.split(maxsplit=1)[0]
            if seq_id not in wanted:
                continue
            if seq_id in found:
                duplicate_ids.add(seq_id)
                continue
            found.add(seq_id)
            written_bases += len(sequence)
            output.write(f">{header}\n")
            write_wrapped(output, sequence)

            for locus in loci_by_contig.get(seq_id, []):
                start = safe_int(locus["start"])
                end = safe_int(locus["end"])
                if start < 1 or end < start or end > len(sequence):
                    output_tmp.unlink(missing_ok=True)
                    gene_tmp.unlink(missing_ok=True)
                    raise RuntimeError(
                        f"Gene coordinates are outside contig {seq_id} (length {len(sequence)}): "
                        f"{locus['locus_id']}={start}-{end}"
                    )
                gene_sequence = sequence[start - 1 : end]
                expected_length = end - start + 1
                if len(gene_sequence) != expected_length:
                    raise AssertionError(
                        f"Internal coordinate error for {locus['locus_id']}: "
                        f"expected {expected_length}, extracted {len(gene_sequence)}"
                    )
                gene_output.write(
                    f">{locus['locus_id']} source_contig={seq_id} source_start={start} "
                    f"source_end={end} source_strand={locus['strand']} "
                    f"organism_type={organism_type}\n"
                )
                write_wrapped(gene_output, gene_sequence)
                written_genes += 1
                written_gene_bases += len(gene_sequence)

    missing = wanted - found
    if missing:
        output_tmp.unlink(missing_ok=True)
        gene_tmp.unlink(missing_ok=True)
        example = ", ".join(sorted(missing)[:10])
        raise RuntimeError(
            f"{len(missing)} hit contig IDs were absent from the input FASTA. Examples: {example}"
        )
    if duplicate_ids:
        output_tmp.unlink(missing_ok=True)
        gene_tmp.unlink(missing_ok=True)
        example = ", ".join(sorted(duplicate_ids)[:10])
        raise RuntimeError(f"Duplicate contig identifiers detected. Examples: {example}")

    if written_genes != sum(len(records) for records in loci_by_contig.values()):
        output_tmp.unlink(missing_ok=True)
        gene_tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"Expected {sum(len(records) for records in loci_by_contig.values())} gene sequences "
            f"but extracted {written_genes}"
        )

    os.replace(output_tmp, output_fasta)
    os.replace(gene_tmp, gene_fasta)
    return {
        "requested_contigs": len(wanted),
        "found_contigs": len(found),
        "scanned_contigs": scanned,
        "written_bases": written_bases,
        "written_genes": written_genes,
        "written_gene_bases": written_gene_bases,
    }


def compress_with_pigz(
    source: Path,
    pigz: str,
    threads: int,
    logger: RunLogger,
    keep_source: bool,
) -> Path:
    destination = Path(str(source) + ".gz")
    temp_destination = Path(str(destination) + ".tmp")
    run_external(
        [pigz, "-p", str(threads), "-c", str(source)],
        logger=logger,
        stdout_path=temp_destination,
        stderr_path=Path(str(source) + ".pigz.log"),
    )
    os.replace(temp_destination, destination)
    if not keep_source:
        source.unlink(missing_ok=True)
    return destination


def write_summary_tsv(path: Path, summary: dict[str, Any]) -> None:
    temp = path.with_name(path.name + ".tmp")
    with open(temp, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(summary.keys())
        writer.writerow(summary.values())
    os.replace(temp, path)


def build_parser() -> argparse.ArgumentParser:
    RawDescriptionRichHelpFormatter.styles["argparse.args"] = "bold cyan"
    RawDescriptionRichHelpFormatter.styles["argparse.groups"] = "bold dark_orange"
    RawDescriptionRichHelpFormatter.styles["argparse.metavar"] = "bold green"
    RawDescriptionRichHelpFormatter.styles["argparse.prog"] = "bold magenta"

    program_name = Path(sys.argv[0]).name
    example_command = (
        f"python3 {program_name}" if program_name.endswith(".py") else program_name
    )

    parser = argparse.ArgumentParser(
        prog=program_name,
        description=(
            "Map reference proteins to one metagenomic assembly with miniprot, retain "
            "high-quality gene models, and export contig-coordinate plus standalone-gene "
            "GFF3/FASTA references."
        ),
        epilog=f"""
Example:
  {example_command} \\
      --proteins phi-base_current.fas \\
      --contigs 201704_MF1.fasta.gz \\
      --outdir 201704_MF1.phi_scan \\
      --sample 201704_MF1 \\
      --organism_type eukaryote \\
      --threads 24

Resume is automatic. Use --force to discard compatible checkpoints and rerun.
""",
        formatter_class=RawDescriptionRichHelpFormatter,
        add_help=True,
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-p",
        "--proteins",
        required=True,
        type=Path,
        metavar="FASTA",
        help="Reference protein FASTA, optionally .gz/.bgz.",
    )
    required.add_argument(
        "-c",
        "--contigs",
        required=True,
        type=Path,
        metavar="FASTA",
        help="One metagenomic contig FASTA, optionally .gz/.bgz.",
    )
    required.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        metavar="DIR",
        help="Output directory for this dataset.",
    )

    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "--sample",
        default=None,
        metavar="NAME",
        help="Sample name. Default: inferred from contig filename.",
    )
    optional.add_argument(
        "--organism_type",
        choices=("eukaryote", "prokaryote"),
        default="eukaryote",
        help=(
            "Annotation mode. Eukaryote enables splice-aware alignment and introns; "
            "prokaryote disables splicing and omits introns. Default: eukaryote."
        ),
    )
    optional.add_argument(
        "-t",
        "--threads",
        type=int,
        default=detect_threads(),
        metavar="INT",
        help="CPU threads. Default: scheduler allocation or min(system CPUs, 32).",
    )
    optional.add_argument(
        "--compression_threads",
        type=int,
        default=None,
        metavar="INT",
        help="pigz threads. Default: min(--threads, 8).",
    )
    optional.add_argument(
        "--tmpdir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Temporary parent directory. Default: $SLURM_TMPDIR or OUTDIR/.work.",
    )
    optional.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Resume completed compatible stages. Default: enabled.",
    )
    optional.add_argument(
        "--force",
        action="store_true",
        help="Rerun all stages and overwrite known outputs.",
    )
    optional.add_argument(
        "--keep_index",
        action="store_true",
        help="Keep the miniprot .mpi index after successful completion.",
    )
    optional.add_argument(
        "--keep_uncompressed",
        action="store_true",
        help="Keep uncompressed FASTA outputs in addition to .gz files.",
    )
    optional.add_argument(
        "--skip_sequence_export",
        action="store_true",
        help=(
            "Do not run gffread or export CDS/protein/transcript FASTA files. "
            "The standalone gene FASTA is still generated."
        ),
    )

    alignment = parser.add_argument_group("miniprot alignment")
    alignment.add_argument(
        "--splice_model",
        type=int,
        choices=(0, 1, 2),
        default=None,
        help="Eukaryotic splice model. Default: 1; omit or set 0 in prokaryote mode.",
    )
    alignment.add_argument("--max_intron", type=int, default=20_000, metavar="BP")
    alignment.add_argument("--index_subsample", type=int, default=1, metavar="INT")
    alignment.add_argument("--max_hits", type=int, default=50, metavar="INT")
    alignment.add_argument("--secondary_ratio", type=float, default=0.50, metavar="FLOAT")
    alignment.add_argument(
        "--prefilter_query_coverage", type=float, default=0.30, metavar="FLOAT"
    )
    alignment.add_argument("--min_score_ratio", type=float, default=0.50, metavar="FLOAT")
    alignment.add_argument(
        "--miniprot_extra",
        default="",
        metavar="'OPTIONS'",
        help="Additional miniprot mapping options parsed with shlex; use cautiously.",
    )

    filtering = parser.add_argument_group("final hit filtering")
    filtering.add_argument("--min_identity", type=float, default=0.40, metavar="FLOAT")
    filtering.add_argument(
        "--min_query_coverage", type=float, default=0.60, metavar="FLOAT"
    )
    filtering.add_argument("--max_frameshift", type=int, default=1, metavar="INT")
    filtering.add_argument("--max_stop_codon", type=int, default=0, metavar="INT")
    filtering.add_argument("--locus_overlap", type=float, default=0.80, metavar="FLOAT")
    filtering.add_argument("--high_identity", type=float, default=0.60, metavar="FLOAT")
    filtering.add_argument(
        "--high_query_coverage", type=float, default=0.80, metavar="FLOAT"
    )

    executables = parser.add_argument_group("external executables")
    executables.add_argument("--miniprot", default="miniprot", metavar="PATH")
    executables.add_argument("--gffread", default="gffread", metavar="PATH")
    executables.add_argument("--pigz", default="pigz", metavar="PATH")

    help_group = parser.add_argument_group("detailed help")
    help_group.add_argument(
        "--help_input",
        action=ShowPanelAction,
        panel_text=INPUT_HELP,
        panel_title="Input format help",
        help="Show detailed input-file requirements and exit.",
    )
    help_group.add_argument(
        "--help_default",
        action=ShowPanelAction,
        panel_text=DEFAULT_HELP,
        panel_title="Default parameter help",
        help="Show detailed default parameters and exit.",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    args.proteins = args.proteins.expanduser().resolve()
    args.contigs = args.contigs.expanduser().resolve()
    args.outdir = args.outdir.expanduser().resolve()
    if not args.proteins.is_file() or args.proteins.stat().st_size == 0:
        raise FileNotFoundError(f"Protein FASTA is missing or empty: {args.proteins}")
    if not args.contigs.is_file() or args.contigs.stat().st_size == 0:
        raise FileNotFoundError(f"Contig FASTA is missing or empty: {args.contigs}")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.compression_threads is None:
        args.compression_threads = min(args.threads, 8)
    if args.compression_threads < 1:
        raise ValueError("--compression_threads must be >= 1")
    if args.organism_type == "eukaryote":
        if args.splice_model is None:
            args.splice_model = 1
    else:
        if args.splice_model not in (None, 0):
            raise ValueError(
                "--splice_model is not used for prokaryotes; omit it or set it to 0"
            )
        args.splice_model = 0
    if args.max_intron < 1 or args.max_hits < 1 or args.index_subsample < 0:
        raise ValueError("Invalid miniprot integer parameter")
    for name in (
        "secondary_ratio",
        "prefilter_query_coverage",
        "min_score_ratio",
        "min_identity",
        "min_query_coverage",
        "locus_overlap",
        "high_identity",
        "high_query_coverage",
    ):
        value = getattr(args, name)
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"--{name} must be between 0 and 1")
    if args.high_identity < args.min_identity:
        raise ValueError("--high_identity must be >= --min_identity")
    if args.high_query_coverage < args.min_query_coverage:
        raise ValueError("--high_query_coverage must be >= --min_query_coverage")
    if args.max_frameshift < 0 or args.max_stop_codon < 0:
        raise ValueError("Frameshift and stop-codon limits must be >= 0")


def output_paths(outdir: Path, sample: str) -> dict[str, Path]:
    prefix = outdir / sample
    return {
        "log": Path(str(prefix) + ".run.log"),
        "state": Path(str(prefix) + ".state.json"),
        "metadata": Path(str(prefix) + ".run_metadata.json"),
        "metrics": Path(str(prefix) + ".stage_metrics.tsv"),
        "clean_proteins": Path(str(prefix) + ".reference.clean.faa"),
        "reference_map": Path(str(prefix) + ".reference_map.tsv"),
        "index": Path(str(prefix) + ".miniprot.mpi"),
        "index_log": Path(str(prefix) + ".miniprot.index.log"),
        "raw_gff": Path(str(prefix) + ".miniprot.raw.gff3"),
        "map_log": Path(str(prefix) + ".miniprot.map.log"),
        "all_hits": Path(str(prefix) + ".all_hits.tsv"),
        "best_loci": Path(str(prefix) + ".best_loci.tsv"),
        "best_gff": Path(str(prefix) + ".best_loci.gff3"),
        "gene_gff": Path(str(prefix) + ".gene.gff3"),
        "gene_fasta": Path(str(prefix) + ".gene.fasta"),
        "hit_ids": Path(str(prefix) + ".hit_contig_ids.txt"),
        "query_summary": Path(str(prefix) + ".query_summary.tsv"),
        "contig_summary": Path(str(prefix) + ".contig_summary.tsv"),
        "hit_contigs": Path(str(prefix) + ".hit_contigs.fasta"),
        "cds": Path(str(prefix) + ".genes.cds.fasta"),
        "proteins": Path(str(prefix) + ".genes.protein.fasta"),
        "transcripts": Path(str(prefix) + ".genes.transcript.fasta"),
        "gffread_log": Path(str(prefix) + ".gffread.log"),
        "summary": Path(str(prefix) + ".summary.tsv"),
        "done": Path(str(prefix) + ".done"),
    }


def remove_known_outputs(paths: dict[str, Path], keep_log: bool = True) -> None:
    for key, path in paths.items():
        if keep_log and key == "log":
            continue
        path.unlink(missing_ok=True)
        Path(str(path) + ".gz").unlink(missing_ok=True)
        Path(str(path) + ".tmp").unlink(missing_ok=True)


def print_parameter_table(
    logger: RunLogger,
    args: argparse.Namespace,
    sample: str,
    tool_versions: dict[str, str],
) -> None:
    table = Table(title="Run parameters", show_header=True, header_style="bold cyan")
    table.add_column("Parameter", style="bold")
    table.add_column("Value", overflow="fold")
    rows = [
        ("sample", sample),
        ("proteins", str(args.proteins)),
        ("contigs", str(args.contigs)),
        ("contig file size", human_bytes(args.contigs.stat().st_size)),
        ("outdir", str(args.outdir)),
        ("organism type", args.organism_type),
        ("threads", str(args.threads)),
        ("compression threads", str(args.compression_threads)),
        ("min identity", f"{args.min_identity:.3f}"),
        ("min query coverage", f"{args.min_query_coverage:.3f}"),
        ("high confidence", f"identity >= {args.high_identity:.3f}; qcov >= {args.high_query_coverage:.3f}; fs=0; stop=0"),
        ("max intron", f"{args.max_intron:,} bp"),
        ("resume", str(args.resume)),
        ("miniprot", tool_versions["miniprot"]),
        ("gffread", tool_versions["gffread"]),
        ("pigz", tool_versions["pigz"]),
    ]
    for key, value in rows:
        table.add_row(key, value)
    logger.console.print(table)
    logger.file_console.print(table)


def stage_outputs(paths: dict[str, Path], skip_sequence_export: bool) -> dict[str, list[Path]]:
    outputs = {
        "prepare_reference": [paths["clean_proteins"], paths["reference_map"]],
        "build_index": [paths["index"]],
        "map_proteins": [paths["raw_gff"]],
        "parse_hits": [
            paths["all_hits"],
            paths["best_loci"],
            paths["best_gff"],
            paths["gene_gff"],
            paths["hit_ids"],
            paths["query_summary"],
            paths["contig_summary"],
        ],
        "extract_contigs": [
            Path(str(paths["hit_contigs"]) + ".gz"),
            paths["gene_fasta"],
        ],
        "export_sequences": [],
        "finalize": [paths["summary"], paths["metadata"], paths["metrics"], paths["done"]],
    }
    if not skip_sequence_export:
        outputs["export_sequences"] = [
            Path(str(paths["cds"]) + ".gz"),
            Path(str(paths["proteins"]) + ".gz"),
            Path(str(paths["transcripts"]) + ".gz"),
        ]
    return outputs


def write_stage_metrics(path: Path, state: PipelineState) -> None:
    rows: list[dict[str, Any]] = []
    for stage in STAGE_ORDER:
        record = state.data.get("stages", {}).get(stage)
        if record:
            rows.append(
                {
                    "stage": stage,
                    "status": record.get("status", ""),
                    "started_at": record.get("started_at", ""),
                    "finished_at": record.get("finished_at", ""),
                    "wall_seconds": record.get("wall_seconds", ""),
                    "cpu_seconds": record.get("cpu_seconds", ""),
                    "max_rss_mb": record.get("max_rss_mb", ""),
                    "detail": record.get("detail", ""),
                }
            )
    write_tsv(
        path,
        rows,
        [
            "stage",
            "status",
            "started_at",
            "finished_at",
            "wall_seconds",
            "cpu_seconds",
            "max_rss_mb",
            "detail",
        ],
    )


def summarize_best_loci_tsv(path: Path) -> dict[str, int]:
    selected = 0
    high = 0
    medium = 0
    contigs: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            selected += 1
            high += row.get("confidence") == "high"
            medium += row.get("confidence") == "medium"
            if row.get("contig"):
                contigs.add(row["contig"])
    return {
        "selected_loci": selected,
        "high_confidence_loci": high,
        "medium_confidence_loci": medium,
        "hit_contigs": len(contigs),
    }


def summarize_gff_features(path: Path) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) == 9:
                counts[fields[2]] += 1
    return {"exons": counts["exon"], "introns": counts["intron"]}


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)

    sample = args.sample or infer_sample_name(args.contigs)
    sample = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample).strip("_")
    if not sample:
        raise ValueError("Sample name became empty after sanitization")

    args.outdir.mkdir(parents=True, exist_ok=True)
    terminal_width = shutil.get_terminal_size(fallback=(120, 30)).columns
    log_width = max(48, terminal_width // 2)
    paths = output_paths(args.outdir, sample)
    logger = RunLogger(paths["log"], width=log_width)

    try:
        miniprot = require_executable(args.miniprot)
        pigz = require_executable(args.pigz)
        gffread = require_executable(args.gffread) if not args.skip_sequence_export else args.gffread
        tool_versions = {
            "miniprot": capture_tool_version(miniprot),
            "gffread": capture_tool_version(gffread) if not args.skip_sequence_export else "skipped",
            "pigz": capture_tool_version(pigz),
        }

        signature_payload = {
            "program_version": PROGRAM_VERSION,
            "proteins": file_signature(args.proteins, content_hash=True),
            "contigs": file_signature(args.contigs, content_hash=False),
            "sample": sample,
            "parameters": {
                key: value
                for key, value in vars(args).items()
                if key not in {"force", "resume", "tmpdir", "outdir"}
            },
            "tools": tool_versions,
        }
        signature_payload["parameters"] = {
            key: str(value) if isinstance(value, Path) else value
            for key, value in signature_payload["parameters"].items()
        }
        signature = stable_hash(signature_payload)

        existing_signature = None
        if paths["state"].exists():
            with open(paths["state"], "r", encoding="utf-8") as handle:
                existing_signature = json.load(handle).get("signature")

        if args.force:
            logger.warning("--force enabled: removing known outputs and checkpoints")
            remove_known_outputs(paths, keep_log=True)
        elif existing_signature is not None and existing_signature != signature:
            logger.warning("Input, parameter, or tool signature changed; incompatible checkpoints will be rebuilt")
            remove_known_outputs(paths, keep_log=True)

        state = PipelineState(paths["state"], signature)
        if args.force or not args.resume or not state.compatible():
            state.reset()

        print_parameter_table(logger, args, sample, tool_versions)
        logger.info(f"Persistent log: {paths['log']}")
        logger.info(f"Checkpoint state: {paths['state']}")

        if args.tmpdir is not None:
            work_parent = args.tmpdir.expanduser().resolve()
        elif os.environ.get("SLURM_TMPDIR"):
            work_parent = Path(os.environ["SLURM_TMPDIR"]).resolve()
        else:
            work_parent = args.outdir / ".work"
        workdir = work_parent / f"{sample}.{os.getpid()}"
        workdir.mkdir(parents=True, exist_ok=True)

        outputs = stage_outputs(paths, args.skip_sequence_export)
        if args.resume and not args.force and state.is_done("finalize", outputs["finalize"]):
            logger.success(f"Pipeline already completed with matching inputs and parameters: {sample}")
            logger.info(f"Summary: {paths['summary']}")
            if workdir.exists():
                shutil.rmtree(workdir)
            return 0

        parse_stats: dict[str, Any] = {}
        reference_stats: dict[str, Any] = {}
        extraction_stats: dict[str, Any] = {}
        stage_count = len(STAGE_ORDER)

        progress = Progress(
            SpinnerColumn(style="cyan"),
            TextColumn("[bold]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            console=logger.console,
            transient=False,
        )

        with progress:
            task = progress.add_task("Initializing", total=stage_count)

            stage = "prepare_reference"
            progress.update(task, description="Preparing protein reference")
            if args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: protein reference already prepared")
            else:
                timer = StageTimer(stage)
                reference_stats = prepare_reference(
                    args.proteins, paths["clean_proteins"], paths["reference_map"]
                )
                metric = timer.finish("completed", json.dumps(reference_stats, sort_keys=True))
                state.mark_completed(metric, outputs[stage])
                logger.success(
                    f"Prepared {reference_stats['written_records']:,} proteins "
                    f"({reference_stats['total_amino_acids']:,} aa) in {human_seconds(metric.wall_seconds)}"
                )
            progress.advance(task)

            stage = "build_index"
            progress.update(task, description="Building miniprot index")
            if args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: miniprot index is complete")
            else:
                timer = StageTimer(stage)
                index_tmp = workdir / f"{sample}.miniprot.mpi"
                run_external(
                    [
                        miniprot,
                        "-t",
                        str(args.threads),
                        "-M",
                        str(args.index_subsample),
                        "-d",
                        str(index_tmp),
                        str(args.contigs),
                    ],
                    logger,
                    stdout_path=None,
                    stderr_path=paths["index_log"],
                )
                move_replace(index_tmp, paths["index"])
                metric = timer.finish(
                    "completed",
                    f"input_size={args.contigs.stat().st_size}; index_size={paths['index'].stat().st_size}",
                )
                state.mark_completed(metric, outputs[stage])
                rate = args.contigs.stat().st_size / max(metric.wall_seconds, 1e-9)
                logger.success(
                    f"Index built in {human_seconds(metric.wall_seconds)}; physical-input throughput {human_bytes(rate)}/s"
                )
            progress.advance(task)

            stage = "map_proteins"
            progress.update(task, description="Mapping proteins to contigs")
            if args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: miniprot mapping is complete")
            else:
                timer = StageTimer(stage)
                raw_tmp = workdir / f"{sample}.miniprot.raw.gff3"
                annotation_options = miniprot_annotation_options(args)
                command = [
                    miniprot,
                    "-t",
                    str(args.threads),
                    *annotation_options,
                    "-N",
                    str(args.max_hits),
                    "-p",
                    str(args.secondary_ratio),
                    "--outn",
                    str(args.max_hits),
                    "--outs",
                    str(args.min_score_ratio),
                    "--outc",
                    str(args.prefilter_query_coverage),
                    "--gff",
                    "--gff-delim",
                    "@",
                    *shlex.split(args.miniprot_extra),
                    str(paths["index"]),
                    str(paths["clean_proteins"]),
                ]
                run_external(command, logger, stdout_path=raw_tmp, stderr_path=paths["map_log"])
                if raw_tmp.stat().st_size == 0:
                    with open(raw_tmp, "w", encoding="utf-8") as handle:
                        handle.write("##gff-version 3\n")
                move_replace(raw_tmp, paths["raw_gff"])
                metric = timer.finish("completed", f"raw_gff_size={paths['raw_gff'].stat().st_size}")
                state.mark_completed(metric, outputs[stage])
                logger.success(f"Protein mapping completed in {human_seconds(metric.wall_seconds)}")
            progress.advance(task)

            stage = "parse_hits"
            progress.update(task, description="Filtering and collapsing loci")
            if args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: hit parsing and locus selection are complete")
            else:
                timer = StageTimer(stage)
                parse_stats = parse_and_filter_hits(
                    raw_gff=paths["raw_gff"],
                    reference_map_tsv=paths["reference_map"],
                    sample=sample,
                    args=args,
                    all_hits_tsv=paths["all_hits"],
                    best_loci_tsv=paths["best_loci"],
                    best_gff=paths["best_gff"],
                    gene_gff=paths["gene_gff"],
                    hit_contig_ids=paths["hit_ids"],
                    query_summary_tsv=paths["query_summary"],
                    contig_summary_tsv=paths["contig_summary"],
                )
                metric = timer.finish("completed", json.dumps(parse_stats, sort_keys=True))
                state.mark_completed(metric, outputs[stage])
                logger.success(
                    f"Selected {parse_stats['selected_loci']:,} loci on {parse_stats['hit_contigs']:,} contigs; "
                    f"high confidence: {parse_stats['high_confidence_loci']:,}"
                )
            progress.advance(task)

            stage = "extract_contigs"
            progress.update(task, description="Extracting hit contigs")
            if args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: hit contigs are already extracted")
            else:
                timer = StageTimer(stage)
                extraction_stats = extract_hit_contigs(
                    args.contigs,
                    paths["hit_ids"],
                    paths["hit_contigs"],
                    paths["best_loci"],
                    paths["gene_fasta"],
                    args.organism_type,
                )
                compressed = compress_with_pigz(
                    paths["hit_contigs"],
                    pigz=pigz,
                    threads=args.compression_threads,
                    logger=logger,
                    keep_source=args.keep_uncompressed,
                )
                metric = timer.finish("completed", json.dumps(extraction_stats, sort_keys=True))
                state.mark_completed(metric, [compressed, paths["gene_fasta"]])
                logger.success(
                    f"Extracted {extraction_stats['found_contigs']:,} contigs "
                    f"and {extraction_stats['written_genes']:,} standalone genes"
                )
            progress.advance(task)

            stage = "export_sequences"
            progress.update(task, description="Exporting CDS and proteins")
            if args.skip_sequence_export:
                metric = StageTimer(stage).finish("completed", "sequence export skipped by user")
                state.mark_completed(metric, [])
                logger.warning("Sequence export skipped by --skip_sequence_export")
            elif args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: CDS/protein/transcript FASTA files are complete")
            else:
                timer = StageTimer(stage)
                hit_contigs_uncompressed = paths["hit_contigs"]
                created_temporary_uncompressed = False
                if not hit_contigs_uncompressed.exists():
                    with gzip.open(Path(str(hit_contigs_uncompressed) + ".gz"), "rb") as source, open(
                        hit_contigs_uncompressed, "wb"
                    ) as target:
                        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
                    created_temporary_uncompressed = True

                has_loci = paths["hit_ids"].exists() and paths["hit_ids"].stat().st_size > 0
                if has_loci:
                    run_external(
                        [
                            gffread,
                            "-E",
                            str(paths["best_gff"]),
                            "-g",
                            str(hit_contigs_uncompressed),
                            "-x",
                            str(paths["cds"]),
                            "-y",
                            str(paths["proteins"]),
                            "-w",
                            str(paths["transcripts"]),
                        ],
                        logger,
                        stdout_path=None,
                        stderr_path=paths["gffread_log"],
                    )
                else:
                    logger.warning("No passing loci; creating empty sequence FASTA outputs")
                    for fasta_path in (paths["cds"], paths["proteins"], paths["transcripts"]):
                        fasta_path.touch()

                compressed_outputs: list[Path] = []
                for fasta_path in (paths["cds"], paths["proteins"], paths["transcripts"]):
                    if not fasta_path.exists():
                        fasta_path.touch()
                    compressed_outputs.append(
                        compress_with_pigz(
                            fasta_path,
                            pigz=pigz,
                            threads=args.compression_threads,
                            logger=logger,
                            keep_source=args.keep_uncompressed,
                        )
                    )
                if created_temporary_uncompressed and not args.keep_uncompressed:
                    hit_contigs_uncompressed.unlink(missing_ok=True)

                metric = timer.finish("completed", "gffread sequence export and pigz compression")
                state.mark_completed(metric, compressed_outputs)
                logger.success(f"Sequence export completed in {human_seconds(metric.wall_seconds)}")
            progress.advance(task)

            stage = "finalize"
            progress.update(task, description="Writing final summary")
            if args.resume and state.is_done(stage, outputs[stage]):
                logger.info("Resume: final summary already exists")
            else:
                timer = StageTimer(stage)

                if not parse_stats:
                    parse_stats = summarize_best_loci_tsv(paths["best_loci"])
                    parse_stats.update(summarize_gff_features(paths["best_gff"]))

                stage_wall_total = sum(
                    float(record.get("wall_seconds", 0.0))
                    for record in state.data.get("stages", {}).values()
                )
                summary = {
                    "sample_id": sample,
                    "organism_type": args.organism_type,
                    "proteins_input": str(args.proteins),
                    "contigs_input": str(args.contigs),
                    "contig_file_bytes": args.contigs.stat().st_size,
                    "selected_loci": parse_stats.get("selected_loci", 0),
                    "high_confidence_loci": parse_stats.get("high_confidence_loci", 0),
                    "medium_confidence_loci": parse_stats.get("medium_confidence_loci", 0),
                    "hit_contigs": parse_stats.get("hit_contigs", 0),
                    "exons": parse_stats.get("exons", 0),
                    "introns": parse_stats.get("introns", 0),
                    "threads": args.threads,
                    "completed_at": iso_now(),
                    "completed_stage_wall_seconds": round(stage_wall_total, 3),
                    "signature": signature,
                }
                write_summary_tsv(paths["summary"], summary)

                metadata = {
                    "program": PROGRAM,
                    "program_version": PROGRAM_VERSION,
                    "command": sys.argv,
                    "signature": signature,
                    "signature_payload": signature_payload,
                    "sample": sample,
                    "tool_versions": tool_versions,
                    "summary": summary,
                    "outputs": {
                        key: str(value)
                        for key, value in paths.items()
                        if value.exists() or Path(str(value) + ".gz").exists()
                    },
                }
                atomic_write_json(paths["metadata"], metadata)
                write_stage_metrics(paths["metrics"], state)
                with open(paths["done"], "w", encoding="utf-8") as handle:
                    handle.write(f"{sample}\t{iso_now()}\t{signature}\n")

                metric = timer.finish("completed", "summary, metadata, metrics, and completion marker")
                state.mark_completed(metric, outputs[stage])
                write_stage_metrics(paths["metrics"], state)
                state.mark_completed(metric, outputs[stage])
                logger.success("Final summary and metadata written")
            progress.advance(task)

        if not args.keep_index:
            paths["index"].unlink(missing_ok=True)
            logger.info("Removed miniprot index; use --keep_index to retain it")
        if workdir.exists():
            shutil.rmtree(workdir)
        if work_parent == args.outdir / ".work" and work_parent.exists() and not any(work_parent.iterdir()):
            work_parent.rmdir()

        logger.console.print(
            Panel.fit(
                Text.from_markup(
                    f"[bold green]Completed[/bold green]\n"
                    f"Sample: [cyan]{sample}[/cyan]\n"
                    f"Best loci: [bold]{parse_stats.get('selected_loci', 'see summary')}[/bold]\n"
                    f"Output: {args.outdir}"
                ),
                border_style="green",
            )
        )
        logger.success(f"Pipeline completed: {sample}")
        return 0

    except KeyboardInterrupt:
        logger.error("Interrupted by user; completed checkpoints were retained")
        return 130
    except Exception as error:
        logger.error(str(error))
        logger.error("Pipeline stopped; completed compatible checkpoints were retained")
        return 1
    finally:
        logger.close()


if __name__ == "__main__":
    raise SystemExit(main())
