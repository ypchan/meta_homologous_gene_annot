# Meta Homologous Gene Annotator

`meta_homologous_gene_annot.py` maps a collection of reference proteins to a metagenomic contig assembly, filters candidate protein-coding gene models, collapses redundant loci, and exports annotation tables, GFF3, hit contigs, CDS sequences, proteins, and transcripts.

The pipeline uses [miniprot](https://github.com/lh3/miniprot) for splice-aware protein-to-genome alignment, gffread for sequence extraction, and pigz for parallel compression. It supports stage-level checkpointing and records its inputs, parameters, tool versions, runtime, and resource usage.

> This tool identifies candidate genes that are homologous to the reference proteins. A match to a pathogen-associated protein does not, by itself, prove that a contig originated from the corresponding pathogen. It is not a substitute for taxonomic, pathogenicity, or experimental validation.

## Workflow

1. Read and clean the reference proteins, assigning a unique internal ID to every record.
2. Build a miniprot index for the contig FASTA.
3. Search for candidate gene models with miniprot.
4. Filter hits by amino acid identity, query coverage, frameshifts, and internal stop codons.
5. Cluster overlapping hits on the same contig and strand, retaining the best model at each locus.
6. Stream the assembly and extract all hit contigs without loading the entire FASTA into memory.
7. Export CDS, protein, and transcript FASTA files with gffread, then compress them with pigz.
8. Write summary tables, run metadata, stage metrics, and checkpoint state.

## Requirements

- Linux or another Unix-like environment capable of running the required command-line tools
- Python 3.9 or later
- Python packages: `rich` and `rich-argparse`
- External programs: `miniprot`, `gffread`, and `pigz`

gffread is optional when `--skip_sequence_export` is specified. miniprot and pigz are always required.

## Installation

### Option 1: Conda or Mamba (recommended)

Clone or download the repository, then run the following commands from its root directory:

```bash
mamba env create -f environment.yml
mamba activate meta-homologous-gene-annot
python3 meta_homologous_gene_annot.py --help
```

If Mamba is not available, replace `mamba` in the first command with `conda`.

### Option 2: Existing Python environment

Install miniprot, gffread, and pigz and ensure that they are available on `PATH`. Then install the Python dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

Verify the installation:

```bash
python3 --version
miniprot --version
gffread --version
pigz --version
python3 meta_homologous_gene_annot.py --help
```

If an external program is not on `PATH`, pass its executable path with `--miniprot`, `--gffread`, or `--pigz`. The script does not need to be modified.

## Input Files

### Reference proteins: `--proteins/-p`

The input must be a protein FASTA file. Uncompressed files and files ending in `.gz`, `.bgz`, or `.bgzf` are supported.

```text
>PHI:1234 gene=MEP4 pathogen=Trichophyton_mentagrophytes
MKFSLALALAVASASA...
```

The following cleaning rules apply:

- The first whitespace-delimited field in the header becomes `original_id`; the complete header is retained in the output.
- Duplicate original IDs are allowed. The program generates a unique `query_id`, such as `PHIREF000000001`, in input order.
- Sequences are converted to uppercase. Whitespace, `-`, `.`, and terminal `*` characters are removed.
- Supported amino acid characters are `ABCDEFGHIKLMNPQRSTVWXYZUO`. Other characters, including internal `*`, are replaced with `X`, and the number of replacements is recorded in the reference map.
- Empty sequences are skipped. The program stops if no usable protein sequences remain.

### Contig assembly: `--contigs/-c`

The input must be a nucleotide FASTA file. Uncompressed files and files ending in `.gz`, `.bgz`, or `.bgzf` are supported. The first whitespace-delimited field in each header becomes the contig ID. IDs must be unique and must match the sequence IDs reported in the miniprot GFF3 output.

The program streams the FASTA and retains only the set of requested hit IDs in memory. The extraction stage stops with an error if a hit ID from the GFF3 is absent from the FASTA or if a requested contig ID occurs more than once.

## Quick Start

```bash
python3 meta_homologous_gene_annot.py \
  --proteins data/phi-base_current.fas \
  --contigs data/201704_MF1.fasta.gz \
  --outdir results/201704_MF1 \
  --sample 201704_MF1 \
  --threads 24
```

`--sample` may be omitted. In that case, the sample name is inferred from the contig filename by removing `.gz/.bgz/.bgzf` and `.fasta/.fna/.fa/.fas` suffixes in sequence. Characters other than letters, digits, periods, underscores, and hyphens are replaced with `_`.

Display the command-line help and the built-in detailed guides:

```bash
python3 meta_homologous_gene_annot.py --help
python3 meta_homologous_gene_annot.py --help_input
python3 meta_homologous_gene_annot.py --help_default
```

## Parameter Reference

All proportion parameters use decimal values from `0` to `1`, not percentages. For example, `0.40` means 40%.

### Required parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `-p, --proteins FASTA` | None | Reference protein FASTA. gzip- and bgzip-compressed input is supported. The file must exist and must not be empty. |
| `-c, --contigs FASTA` | None | Contig FASTA for one sample. gzip- and bgzip-compressed input is supported. The file must exist and must not be empty. |
| `-o, --outdir DIR` | None | Output directory. It is created automatically when necessary, and the primary results are written directly into it. |

### Execution and output parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--sample NAME` | Inferred from the contig filename | Sample name used as the output prefix and in GFF3 locus IDs. |
| `-t, --threads INT` | Scheduler allocation or `min(CPU, 32)` | Number of miniprot threads. The program checks `SLURM_CPUS_PER_TASK`, `NSLOTS`, and `OMP_NUM_THREADS` in that order. If none is available, it uses the system CPU count, capped at 32. Must be at least 1. |
| `--compression_threads INT` | `min(--threads, 8)` | Number of pigz compression threads. Must be at least 1. |
| `--tmpdir DIR` | `$SLURM_TMPDIR` or `OUTDIR/.work` | Parent directory for temporary files. Each run creates `<sample>.<PID>` below it and removes that directory after successful completion. This parameter is excluded from the checkpoint signature. |
| `--resume / --no-resume` | `--resume` | Enable or disable reuse of compatible completed stages. Disabling resume resets the state and reruns every stage, replacing outputs stage by stage. |
| `--force` | Disabled | Delete known outputs and checkpoints for this sample before rerunning the complete pipeline. The primary run log is retained and appended to. |
| `--keep_index` | Disabled | Retain `<sample>.miniprot.mpi` after successful completion. The index is removed by default to conserve disk space. |
| `--keep_uncompressed` | Disabled | Retain uncompressed hit-contig, CDS, protein, and transcript FASTA files in addition to their `.gz` files. |
| `--skip_sequence_export` | Disabled | Do not run gffread or generate CDS, protein, and transcript FASTA files. Hit contigs are still exported and compressed. |

### Miniprot alignment parameters

| Parameter | Default | Mapping and effect |
| --- | --- | --- |
| `--splice_model {0,1,2}` | `1` | Passed to miniprot as `-j`. Model `1` is the general splice model and is the default for fungal homology searches. Consult the help for the installed miniprot version for the precise meaning of other values. |
| `--max_intron BP` | `20000` | Passed to `-G`. Sets the maximum allowed intron length and must be at least 1. Increase it when the target organisms are expected to have longer introns. |
| `--index_subsample INT` | `1` | Passed to `-M` during index construction. The k-mer sampling rate is `1/2**M`. Must be at least 0. Increasing it can reduce index size and runtime at the cost of sensitivity. |
| `--max_hits INT` | `50` | Passed to both `-N` and `--outn`. Limits the number of candidate alignments retained and reported per query. Must be at least 1. |
| `--secondary_ratio FLOAT` | `0.50` | Passed to `-p`. A secondary hit is retained only if its score reaches this fraction of the best score. |
| `--prefilter_query_coverage FLOAT` | `0.30` | Passed to `--outc`. Sets the query-coverage prefilter applied by miniprot before the final filtering performed by this program. |
| `--min_score_ratio FLOAT` | `0.50` | Passed to `--outs`. A hit is reported only if its score reaches this fraction of the best alignment score. |
| `--miniprot_extra 'OPTIONS'` | Empty | Additional raw miniprot mapping options, split according to shell quoting rules. Use this only for options not exposed directly by the script. Do not include the index or protein input paths. |

### Final hit filtering and locus parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--min_identity FLOAT` | `0.40` | Minimum amino acid identity. Records below this threshold are marked `low_identity`. |
| `--min_query_coverage FLOAT` | `0.60` | Minimum reference-protein coverage, calculated as the miniprot `Target` interval length divided by the cleaned reference-protein length. Records below this threshold are marked `low_query_coverage`. |
| `--max_frameshift INT` | `1` | Maximum permitted number of frameshifts. Must be at least 0. Records above this limit are marked `too_many_frameshifts`. |
| `--max_stop_codon INT` | `0` | Maximum permitted number of internal stop codons. Must be at least 0. Records above this limit are marked `internal_stop_codon`. |
| `--locus_overlap FLOAT` | `0.80` | Threshold used to merge redundant reference hits on the same contig and strand: `overlap length / shorter interval length`. |
| `--high_identity FLOAT` | `0.60` | Amino acid identity threshold for the high-confidence label. It cannot be lower than `--min_identity`. |
| `--high_query_coverage FLOAT` | `0.80` | Query-coverage threshold for the high-confidence label. It cannot be lower than `--min_query_coverage`. |

A hit passes the final filter only when it satisfies all four identity, coverage, frameshift, and stop-codon limits. A high-confidence hit must also satisfy `identity >= high_identity`, `query_coverage >= high_query_coverage`, `frameshift == 0`, and `stop_codon == 0`. Passing loci that do not meet the high-confidence criteria are labeled `medium`. Failed raw hits are labeled `low` and are retained only in `all_hits.tsv`.

When several reference hits belong to the same candidate locus, the representative model is selected by confidence level, `identity × query_coverage`, query coverage, identity, alignment score, fewer frameshifts, fewer stop codons, and lower numeric rank, in that order. The remaining reference IDs are reported in `alternative_references`.

### External executables and help parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--miniprot PATH` | `miniprot` | miniprot command name or executable path. |
| `--gffread PATH` | `gffread` | gffread command name or executable path. It is not checked when `--skip_sequence_export` is used. |
| `--pigz PATH` | `pigz` | pigz command name or executable path. |
| `-h, --help` | — | Display the complete command-line reference and exit. |
| `--help_input` | — | Display input formats, rules, and examples, then exit without requiring the mandatory arguments. |
| `--help_default` | — | Display the default thresholds and their interpretation, then exit without requiring the mandatory arguments. |

## Configuration Examples

Increase sensitivity for more distant homologs. This produces more candidates and requires more careful downstream review:

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/sensitive \
  --min_identity 0.30 \
  --min_query_coverage 0.50 \
  --secondary_ratio 0.30 \
  --min_score_ratio 0.30
```

Generate alignments, tables, GFF3, and hit contigs without exporting gene sequences:

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/tables_only \
  --skip_sequence_export
```

Use scheduler-local temporary storage and retain the miniprot index:

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/sample01 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --tmpdir "${SLURM_TMPDIR}" \
  --keep_index
```

## Output Files

The table below assumes a sample name of `sample01`. All files are written below `--outdir`.

| File | Contents |
| --- | --- |
| `sample01.summary.tsv` | One-row core summary containing the inputs, locus counts, confidence counts, hit-contig count, threads, completion time, and run signature. |
| `sample01.best_loci.tsv` | Representative model for each collapsed locus. This is the primary table for downstream analysis. |
| `sample01.best_loci.gff3` | Filtered, nonredundant mRNA and child features. Each mRNA receives `Locus`, `ReferenceID`, `Confidence`, `QueryCoverage`, and `ReferenceAnnotation` attributes. |
| `sample01.all_hits.tsv` | Every mRNA hit reported by miniprot, including records that failed the final filter, `pass_filter`, and `fail_reason`. |
| `sample01.query_summary.tsv` | Raw-hit, passing-hit, and selected-locus counts for each reference protein, plus its best hit. Unmapped references are labeled `unmapped`. |
| `sample01.contig_summary.tsv` | Locus count, high- and medium-confidence counts, best metrics, and reference ID set for each hit contig. |
| `sample01.hit_contig_ids.txt` | Sorted, deduplicated hit-contig IDs. |
| `sample01.hit_contigs.fasta.gz` | Complete sequences of all hit contigs. |
| `sample01.genes.cds.fasta.gz` | CDS sequences exported by gffread. Not generated with `--skip_sequence_export`. |
| `sample01.genes.protein.fasta.gz` | Translated proteins exported by gffread. Not generated with `--skip_sequence_export`. |
| `sample01.genes.transcript.fasta.gz` | Transcripts exported by gffread. Not generated with `--skip_sequence_export`. |
| `sample01.reference.clean.faa` | Cleaned and renumbered reference proteins used by miniprot. |
| `sample01.reference_map.tsv` | Mapping from internal `query_id` values to original IDs, complete headers, protein lengths, and replaced-residue counts. |
| `sample01.miniprot.raw.gff3` | Raw miniprot GFF3 before final filtering by this program. |
| `sample01.run_metadata.json` | Program version, complete command, input signature, parameters, tool versions, result paths, and summary. |
| `sample01.stage_metrics.tsv` | Start and finish times, wall-clock duration, child-process CPU time, maximum RSS, and details for each stage. |
| `sample01.state.json` | Resume state and output sizes for completed stages. Do not edit this file while the pipeline is running. |
| `sample01.run.log` | Primary run log and the external commands that were executed. |
| `sample01.miniprot.index.log` | miniprot index-building standard error log. |
| `sample01.miniprot.map.log` | miniprot mapping standard error log. |
| `sample01.gffread.log` | gffread log. Not generated when sequence export is skipped. |
| `*.pigz.log` | pigz log for each FASTA compression step. |
| `sample01.done` | Successful-completion marker containing the sample name, completion time, and run signature. |
| `sample01.miniprot.mpi` | miniprot index. Removed after successful completion unless `--keep_index` is specified. |

When `--keep_uncompressed` is specified, uncompressed versions of the four compressed FASTA output types are retained without the `.gz` suffix.

### Primary result fields

Common fields in `all_hits.tsv` and `best_loci.tsv` include:

- `query_id` / `original_id`: unique internal reference ID / original FASTA ID.
- `model_id`: model ID assigned by miniprot. `locus_id` is the stable, nonredundant output ID.
- `contig`, `start`, `end`, `strand`: genomic location of the gene model.
- `locus_start`, `locus_end`: complete interval spanned by the merged cluster.
- `score`, `rank`, `identity`, `positive`: miniprot alignment metrics.
- `query_start`, `query_end`, `aligned_query_length`, `query_length`, `query_coverage`: reference-protein coverage information.
- `frameshift`, `stop_codon`: counts of frameshifts and internal stop codons reported by miniprot.
- `confidence`: `high`, `medium`, or `low`.
- `pass_filter`, `fail_reason`: final-filter status and failure reasons. Present only in `all_hits.tsv`.
- `n_reference_hits`, `alternative_references`: number of reference hits at the merged locus and nonrepresentative reference IDs. Present only in `best_loci.tsv`.

## Checkpointing, Replacement, and Run Signatures

`--resume` is enabled by default. A stage is reused only when it is marked as completed, its output files still exist, and their sizes match the recorded values. If the compatible completion marker also exists, repeating the same command exits immediately.

The run signature includes the program version; the reference-protein path, size, modification time, and SHA-256 digest; the contig path, size, and modification time; the sample name; result-affecting parameters; and external tool versions. `--force`, `--resume`, `--tmpdir`, and `--outdir` are excluded from the parameter signature. If the signature is incompatible, the program removes known old outputs for that sample and rebuilds them.

Use `--force` to rerun unconditionally:

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/sample01 \
  --force
```

Do not run two processes with the same `--outdir` and `--sample` at the same time. They would write to the same state and result files.

## Updating

Retrieve code updates with Git and update the Conda environment:

```bash
git pull
mamba env update -f environment.yml --prune
mamba activate meta-homologous-gene-annot
python3 meta_homologous_gene_annot.py --help
```

For a `venv` and pip installation, update the Python dependencies as follows:

```bash
source .venv/bin/activate
python3 -m pip install --upgrade -r requirements.txt
```

After updating the script, miniprot, or gffread, use `--force` to recompute important results. The checkpoint signature records the declared program version and external tool versions, but it does not hash the script itself. If the code changes without a corresponding `PROGRAM_VERSION` change, old stages may still be considered compatible.

The script currently declares version `1.0.0`. Use `git log -1 --oneline` to inspect the local code revision and `git status` to identify uncommitted changes.

## Troubleshooting

### `Executable not found in PATH`

Confirm that the program is installed and can be executed directly, or pass an explicit path such as `--miniprot /absolute/path/to/miniprot`. `--skip_sequence_export` removes the gffread requirement only; miniprot and pigz are still required.

### The pipeline restarts after a parameter change

Changes to result-affecting inputs, parameters, or tool versions invalidate the run signature. The program then removes known old outputs and recomputes them to prevent incompatible stage results from being combined.

### No loci pass the filter

This is a valid result. The program still writes its summary and empty final sequence files. Inspect `fail_reason` in `all_hits.tsv` before deciding whether there is biological justification for changing a threshold. Do not lower thresholds solely to obtain a positive result.

### The gffread stage fails

Inspect `sample01.gffread.log`, `sample01.best_loci.gff3`, and the hit-contig FASTA. Common causes include inconsistent contig IDs, mismatched GFF3 and FASTA files, or behavior differences between gffread versions. After correcting the problem, repeat the same command to resume. Add `--force` only when every stage must be rebuilt.

## Changelog

### 1.0.0

- Added reference-protein cleaning, miniprot alignment, quality filtering, and redundant-locus collapsing.
- Added GFF3, result-table, hit-contig, CDS, protein, and transcript export.
- Added gzip/bgzip input, pigz output compression, stage-level resume support, and run metadata.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
