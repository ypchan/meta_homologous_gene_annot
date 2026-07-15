# Meta Homologous Gene Annotator

`meta_homologous_gene_annot.py` uses miniprot to align reference proteins to a single metagenomic assembly, filters the resulting candidate models, collapses redundant loci, and writes two coordinated GFF3 annotation sets:

- Contig coordinates for inspecting gene locations and protein alignments on the original assembly.
- Gene-reference coordinates in which each reference sequence is one gene. Gene references from multiple samples can be dereplicated so that reads only need to be mapped to a nonredundant gene set for metagenomic abundance or metatranscriptomic expression analysis.

The program supports eukaryotic and prokaryotic modes. Eukaryotic mode performs splice-aware alignment and writes `exon` and `intron` features. Prokaryotic mode disables miniprot splicing and omits `intron` features. These outputs are homology-supported candidate gene models; they are not taxonomic identification, evidence of pathogenicity, or a replacement for complete de novo gene annotation.

## Core Output Design

The table below assumes a sample name of `sample01`.

| File | Coordinate system and purpose |
| --- | --- |
| `sample01.miniprot.raw.gff3` | Raw miniprot GFF3 before final filtering. `##PAF` lines retain detailed `cg` and `cs` protein-to-genome alignments. |
| `sample01.best_loci.gff3` | Filtered and locus-collapsed GFF3 in contig coordinates. It contains the selected `##PAF` records and normalized `gene`, `mRNA`, `exon`, `CDS`, eukaryotic `intron`, and optional `stop_codon` features. |
| `sample01.gene.fasta` | One sequence per nonredundant within-sample locus, with IDs such as `sample01_PHILOCUS...`. Use this file for cross-sample dereplication and read mapping. It is always retained as uncompressed FASTA. |
| `sample01.gene.gff3` | Gene-coordinate GFF3 paired with `sample01.gene.fasta`. The SeqID equals the locus ID and every gene spans `1..gene_length`. |
| `sample01.best_loci.tsv` | Original contig and coordinates, reference protein, identity, query coverage, confidence, and alternative references for each selected locus. |
| `sample01.all_hits.tsv` | Every miniprot mRNA hit, including rejected records and their failure reasons. |
| `sample01.genes.cds.fasta.gz` | Spliced CDS sequences extracted by gffread. |
| `sample01.genes.transcript.fasta.gz` | Spliced transcripts extracted by gffread. Prefer this file when dereplicating eukaryotic genes across samples. |
| `sample01.genes.protein.fasta.gz` | Protein sequences translated by gffread. |
| `sample01.hit_contigs.fasta.gz` | Complete contigs containing at least one retained locus. |

The pipeline also writes query and contig summaries, cleaned reference proteins, reference-ID mappings, logs, tool versions, parameter signatures, stage metrics, checkpoint state, and a completion marker. See `sample01.run_metadata.json` for the complete result list.

### Relationship Between the Two Coordinate Systems

If a gene occupies `[S, E]` on a contig, GFF3 uses 1-based closed intervals:

```text
contig GFF3:  contig_7    gene    S            E
gene FASTA:   >sample01_PHILOCUS0000001, length E-S+1
gene GFF3:    locus_id     gene    1            E-S+1
conversion:   gene_pos = contig_pos - S + 1
```

`sample01.gene.fasta` is extracted with the exact Python slice `sequence[S-1:E]`. A negative-strand gene is not reverse-complemented; the sequence remains in contig orientation and both GFF3 files retain `strand=-`. This keeps FASTA, GFF3, and BAM coordinates consistent. minibwa and minimap2 search both strands during mapping.

Exons are produced by merging miniprot CDS and stop-codon intervals. Eukaryotic introns are the intervals between adjacent exons. Because miniprot is protein guided and generally does not predict UTRs, a gene boundary here means the boundary of the predicted protein-coding model, not necessarily the boundary of a complete transcription unit.

## Requirements and Installation

The annotation pipeline requires:

- Python 3.9 or later.
- The Python packages `rich` and `rich-argparse`.
- [miniprot](https://github.com/lh3/miniprot), [gffread](https://github.com/gpertea/gffread), and [pigz](https://github.com/madler/pigz).

`--skip_sequence_export` skips the gffread CDS, protein, and transcript outputs. `sample.gene.fasta` and both GFF3 files are still generated.

Install the Python dependencies with `pip` in a virtual environment:

```bash
gh repo clone ypchan/meta_homologous_gene_annot
cd meta_homologous_gene_annot

python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

For an existing checkout, start with `python3 -m venv .venv` in the repository root. Reactivate the environment in a new shell with `source .venv/bin/activate`.

The external bioinformatics programs are compiled executables and cannot be installed with `pip`. Install the GitHub CLI, a C/C++ build toolchain, GNU Make, zlib development headers, and the other system libraries required by the upstream projects. Rust and Cargo are required only for CoverM. The commands below clone the official repositories with `gh repo clone` and install executables under `$HOME/.local/bin`.

### Annotation and Mapping Tools

```bash
mkdir -p "$HOME/src" "$HOME/.local/bin"

gh repo clone lh3/miniprot "$HOME/src/miniprot"
make -C "$HOME/src/miniprot"
install -m 0755 "$HOME/src/miniprot/miniprot" "$HOME/.local/bin/miniprot"

gh repo clone gpertea/gffread "$HOME/src/gffread"
make -C "$HOME/src/gffread" release
install -m 0755 "$HOME/src/gffread/gffread" "$HOME/.local/bin/gffread"

gh repo clone madler/pigz "$HOME/src/pigz"
make -C "$HOME/src/pigz"
install -m 0755 "$HOME/src/pigz/pigz" "$HOME/.local/bin/pigz"

gh repo clone weizhongli/cdhit "$HOME/src/cdhit"
make -C "$HOME/src/cdhit" -j 8
install -m 0755 "$HOME/src/cdhit/cd-hit-est" "$HOME/.local/bin/cd-hit-est"

gh repo clone lh3/minibwa "$HOME/src/minibwa"
make -C "$HOME/src/minibwa"
install -m 0755 "$HOME/src/minibwa/minibwa" "$HOME/.local/bin/minibwa"

gh repo clone lh3/minimap2 "$HOME/src/minimap2"
make -C "$HOME/src/minimap2" -j 8
install -m 0755 "$HOME/src/minimap2/minimap2" "$HOME/.local/bin/minimap2"
install -m 0755 "$HOME/src/minimap2/misc/paftools.js" "$HOME/.local/bin/paftools.js"

export PATH="$HOME/.local/bin:$PATH"
```

The eukaryotic RNA example also calls `paftools.js`, which requires the `k8` JavaScript shell. Install the matching `k8` binary from the official minimap2 release assets, or replace that junction-conversion step with another GFF3-aware splice-junction converter.

### Quantification Tools

featureCounts is built as part of Subread:

```bash
gh repo clone ShiLab-Bioinformatics/subread "$HOME/src/subread"
make -C "$HOME/src/subread/src" -f Makefile.Linux -j 8
install -m 0755 "$HOME/src/subread/bin/featureCounts" "$HOME/.local/bin/featureCounts"
```

Build HTSlib and samtools from adjacent GitHub checkouts. A Git checkout of HTSlib requires its submodules plus Autoconf/Automake. HTSlib also requires zlib and normally uses bzip2, xz/liblzma, libcurl, and libdeflate development libraries.

```bash
gh repo clone samtools/htslib "$HOME/src/htslib" -- --recurse-submodules
cd "$HOME/src/htslib"
autoreconf -i
./configure --prefix="$HOME/.local"
make -j 8
make install

gh repo clone samtools/samtools "$HOME/src/samtools"
cd "$HOME/src/samtools"
autoheader
autoconf -Wno-syntax
./configure --prefix="$HOME/.local" --with-htslib="$HOME/src/htslib"
make -j 8
make install
```

Build CoverM from its cloned Rust source tree:

```bash
gh repo clone wwood/CoverM "$HOME/src/CoverM"
cargo install --locked --path "$HOME/src/CoverM" --root "$HOME/.local"
```

Ensure the user-local executable directory is available in every shell:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Check the installed versions and local help because available options can differ between releases:

```bash
python3 --version
python3 -m pip --version
gh --version
miniprot --version
gffread --version
pigz --version
cd-hit-est -h
minibwa --version
minimap2 --version
featureCounts -v
coverm --version
python3 meta_homologous_gene_annot.py --help
```

## Inputs

### Reference Proteins: `--proteins/-p`

Plain and `.gz/.bgz/.bgzf` FASTA files are supported. The first whitespace-delimited header field becomes `original_id`, while the complete header is retained in the results. Duplicate original IDs are allowed; the program assigns unique internal IDs such as `PHIREF000000001` in input order. Sequences are converted to uppercase, whitespace, `-`, `.`, and terminal `*` characters are removed, and unrecognized residues are converted to `X` and recorded.

### Single-Sample Contigs: `--contigs/-c`

Plain and compressed FASTA files are supported. The first header field is the contig ID and must be unique. Large FASTA files are streamed instead of loading the entire assembly into memory.

## Annotation Examples

For eukaryotic genes or other candidates that may contain introns:

```bash
python3 meta_homologous_gene_annot.py \
  --proteins references/target.faa \
  --contigs assemblies/sample01.fasta.gz \
  --outdir results/sample01 \
  --sample sample01 \
  --organism_type eukaryote \
  --splice_model 1 \
  --max_intron 20000 \
  --threads 24
```

For prokaryotic genes:

```bash
python3 meta_homologous_gene_annot.py \
  --proteins references/target.faa \
  --contigs assemblies/sample01.fasta.gz \
  --outdir results/sample01 \
  --sample sample01 \
  --organism_type prokaryote \
  --threads 24
```

`--organism_type prokaryote` passes `-S` to miniprot to disable splicing; do not set a nonzero `--splice_model` in this mode. Eukaryotic mode defaults to `--splice_model 1`, the general eukaryotic/fungal model. miniprot model `2` is intended for vertebrates and insects, while model `0` disables splice-signal scoring.

## Main Parameters

All proportions use decimal values from `0` to `1`, not percentages.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--organism_type` | `eukaryote` | `eukaryote` enables splice-aware models and introns; `prokaryote` passes miniprot `-S` and omits introns. |
| `--splice_model` | Eukaryotic `1` | miniprot `-j`; used only in eukaryotic mode. |
| `--max_intron` | `20000` | Maximum intron length passed to miniprot `-G` in eukaryotic mode. A value that is too small truncates models; a very large value increases runtime and spurious long gaps. |
| `--index_subsample` | `1` | miniprot `-M` during indexing; the sampling rate is `1/2**M`. |
| `--max_hits` | `50` | Passed to both miniprot `-N` and `--outn`. |
| `--secondary_ratio` | `0.50` | miniprot `-p`; the score threshold for secondary hits relative to the best hit. |
| `--prefilter_query_coverage` | `0.30` | miniprot `--outc` prefilter. The program subsequently applies a stricter final filter. |
| `--min_score_ratio` | `0.50` | miniprot `--outs`; the score threshold relative to the best alignment. |
| `--min_identity` | `0.40` | Final minimum amino-acid identity. |
| `--min_query_coverage` | `0.60` | Final minimum reference-protein coverage. |
| `--max_frameshift` | `1` | Maximum permitted number of frameshift events. |
| `--max_stop_codon` | `0` | Maximum permitted number of internal stop codons. |
| `--locus_overlap` | `0.80` | Same-contig, same-strand hits are placed in one locus when `overlap/shorter_interval` reaches this value. |
| `--high_identity` | `0.60` | Minimum identity for the high-confidence label. |
| `--high_query_coverage` | `0.80` | Minimum query coverage for high confidence; high confidence also requires zero frameshifts and zero internal stops. |
| `--miniprot_extra` | Empty | Additional miniprot arguments parsed with shell-quoting rules. Do not repeat the index or protein paths. |
| `--resume/--no-resume` | Resume | Reuse complete stages with matching input, parameter, and tool-version signatures. |
| `--force` | Disabled | Remove known results and checkpoints for this sample and rerun; the main log is appended rather than deleted. |
| `--keep_uncompressed` | Disabled | Also retain uncompressed hit-contig, CDS, protein, and transcript files. `sample.gene.fasta` is always uncompressed. |

Use the built-in detailed guides for additional input and default-value documentation:

```bash
python3 meta_homologous_gene_annot.py --help_input
python3 meta_homologous_gene_annot.py --help_default
```

## External Commands Used by the Pipeline

The key commands are expanded below so that parameters can be audited. Absolute executable paths and complete commands are also recorded in `sample01.run.log`.

Index construction:

```bash
miniprot -t 24 -M 1 \
  -d sample01.miniprot.mpi \
  assemblies/sample01.fasta.gz
```

Eukaryotic mapping:

```bash
miniprot -t 24 -j 1 -G 20000 \
  -N 50 -p 0.50 \
  --outn 50 --outs 0.50 --outc 0.30 \
  --gff --gff-delim '@' \
  sample01.miniprot.mpi sample01.reference.clean.faa \
  > sample01.miniprot.raw.gff3
```

Prokaryotic mapping replaces `-j/-G` with `-S`:

```bash
miniprot -t 24 -S \
  -N 50 -p 0.50 \
  --outn 50 --outs 0.50 --outc 0.30 \
  --gff --gff-delim '@' \
  sample01.miniprot.mpi sample01.reference.clean.faa \
  > sample01.miniprot.raw.gff3
```

The pipeline deliberately uses `--gff`, not `--gff-only`. The `##PAF` line before each model contains 0-based, right-open PAF coordinates, a `cg:Z:` protein CIGAR, and a `cs:Z:` difference string, whereas ordinary GFF3 feature lines use 1-based closed coordinates. CIGAR operations `N/U/V` represent introns in different phases, and `F/G` represent frameshifts. Do not interpret PAF start/end fields as GFF3 coordinates.

Sequence export:

```bash
gffread -E sample01.best_loci.gff3 \
  -g sample01.hit_contigs.fasta \
  -x sample01.genes.cds.fasta \
  -y sample01.genes.protein.fasta \
  -w sample01.genes.transcript.fasta
```

`-x`, `-y`, and `-w` write CDS, protein, and spliced transcript sequences, respectively. The program then compresses them with `pigz -p THREADS -c`.

## Cross-Sample Dereplication

Different samples from the same batch may assemble the same gene. Every locus ID includes the sample prefix, so files can be concatenated without collisions caused by generic contig names.

### 1. Combine Per-Sample Gene FASTA and GFF3 Files

```bash
find results -name '*.gene.fasta' -print0 \
  | sort -z \
  | xargs -0 cat \
  > batch.all.gene.fasta

printf '##gff-version 3\n' > batch.all.gene.gff3
find results -name '*.gene.gff3' -exec awk '!/^#/' {} + \
  >> batch.all.gene.gff3
```

### 2A. Prokaryotes: Dereplicate Genomic Genes Directly

```bash
cd-hit-est \
  -i batch.all.gene.fasta \
  -o batch.nr.gene.fasta \
  -c 0.95 -n 10 \
  -G 1 -aS 0.90 -s 0.80 \
  -g 1 -r 1 -d 0 \
  -T 24 -M 0
```

Parameter details:

- `-c 0.95` requires at least 95% nucleotide identity. Use `0.99` to remove only nearly identical sequences. Consider lower values only when the intended unit is a gene family rather than an individual gene or allele.
- `-n 10` is the nucleotide word length appropriate for identities from 0.95 to 1.0.
- `-G 1` uses global identity with the full shorter sequence as the denominator, preventing short local domains from clustering too easily.
- `-aS 0.90` requires the alignment to cover at least 90% of the shorter sequence.
- `-s 0.80` requires the shorter sequence to be at least 80% of the representative length. This permits partial assemblies to join a more complete representative while excluding very short fragments.
- `-g 1` assigns each sequence to the most similar qualifying cluster. It is slower but improves cluster assignment.
- `-r 1` checks both the forward sequence and reverse complement. Keep it enabled because this program preserves contig orientation.
- `-d 0` retains the complete first FASTA identifier so representatives can be linked back to GFF3.
- `-T 24` sets the thread count. `-M 0` removes the CD-HIT memory limit and should be adjusted to local scheduler policy when necessary.

CD-HIT processes sequences from longest to shortest. The first representative is therefore usually the longest and most complete candidate in that cluster. `-g 1` changes member assignment but not this representative-selection order. Preserve `batch.nr.gene.fasta.clstr`; it records every representative-to-member relationship.

### 2B. Eukaryotes: Prefer Dereplication of Spliced Transcripts

Long introns make full-length clustering of eukaryotic genomic genes difficult. The same coding gene in different assemblies can contain introns with very different sequences and lengths. Dereplicate intron-free transcripts first, then retrieve the genomic gene belonging to each representative transcript:

```bash
find results -name '*.genes.transcript.fasta.gz' -print0 \
  | sort -z \
  | xargs -0 zcat \
  > batch.all.transcript.fasta

cd-hit-est \
  -i batch.all.transcript.fasta \
  -o batch.nr.transcript.fasta \
  -c 0.95 -n 10 \
  -G 1 -aS 0.90 -s 0.80 \
  -g 1 -r 1 -d 0 \
  -T 24 -M 0

grep '^>' batch.nr.transcript.fasta \
  | sed 's/^>//; s/[[:space:]].*//; s/\.t1$//' \
  > batch.nr.gene.ids

awk 'NR==FNR {keep[$1]=1; next}
     /^>/ {id=substr($0,2); sub(/[[:space:]].*/, "", id); emit=(id in keep)}
     emit' \
  batch.nr.gene.ids batch.all.gene.fasta \
  > batch.nr.gene.fasta
```

This example assumes the single transcript produced for each locus has the ID `<locus_id>.t1`. If another tool rewrites the FASTA IDs, inspect the actual `grep '^>'` output before removing any suffix.

### 3. Generate the Nonredundant Gene GFF3

For the direct prokaryotic workflow, first extract representative gene IDs. The eukaryotic workflow above already creates this file.

```bash
grep '^>' batch.nr.gene.fasta \
  | sed 's/^>//; s/[[:space:]].*//' \
  > batch.nr.gene.ids

printf '##gff-version 3\n' > batch.nr.gene.gff3
awk 'NR==FNR {keep[$1]=1; next}
     !/^#/ && ($1 in keep)' \
  batch.nr.gene.ids batch.all.gene.gff3 \
  >> batch.nr.gene.gff3
```

Verify that FASTA and GFF3 contain exactly the same SeqIDs:

```bash
grep '^>' batch.nr.gene.fasta \
  | sed 's/^>//; s/[[:space:]].*//' \
  | sort -u > fasta.ids

awk '!/^#/ {print $1}' batch.nr.gene.gff3 \
  | sort -u > gff.ids

diff -u fasta.ids gff.ids
```

No output from `diff` means that the IDs match. The dereplication threshold defines the final reporting unit: a high threshold is closer to allele or strain level, while a lower threshold approaches gene-family level. Every sample being compared must use the same nonredundant reference and the same mapping, MAPQ, and multimapping rules.

## Read Mapping and Abundance Quantification

### A. minimap2 for Short Metagenomic Reads

```bash
minimap2 -x sr -d batch.nr.gene.sr.mmi batch.nr.gene.fasta

minimap2 -ax sr -t 24 --secondary=yes -N 20 \
  batch.nr.gene.sr.mmi \
  sample01_DNA_R1.fastq.gz sample01_DNA_R2.fastq.gz \
  | samtools sort -@ 8 -o sample01.DNA.nr_gene.bam

samtools index sample01.DNA.nr_gene.bam
samtools flagstat sample01.DNA.nr_gene.bam \
  > sample01.DNA.nr_gene.flagstat.txt
```

`-x sr` selects the short-read preset, `-a` writes SAM, and `--secondary=yes -N 20` retains candidate secondary locations so multimapping among homologous genes can be handled explicitly. Do not vary `-N` or MAPQ rules among samples.

### B. minibwa for Short Metagenomic Reads

```bash
minibwa index -t 24 \
  batch.nr.gene.fasta batch.nr.gene

minibwa map -x sr -t 24 \
  -R $'@RG\tID:sample01.DNA\tSM:sample01.DNA' \
  -N 20 --outn=20 \
  batch.nr.gene \
  sample01_DNA_R1.fastq.gz sample01_DNA_R2.fastq.gz \
  | samtools sort -@ 8 -o sample01.DNA.nr_gene.minibwa.bam

samtools index sample01.DNA.nr_gene.minibwa.bam
samtools flagstat sample01.DNA.nr_gene.minibwa.bam \
  > sample01.DNA.nr_gene.minibwa.flagstat.txt
```

minibwa parameter details:

- `index -t 24 REF PREFIX` builds `PREFIX.l2b` and `PREFIX.mbw` with 24 threads. The default libsais indexer is fast but requires approximately 18 times the reference length in memory. Use `minibwa index -l REF PREFIX` for lower memory usage at the cost of slower indexing.
- `map -x sr` explicitly selects the short-read preset. The default adaptive mode also adjusts parameters according to individual read lengths, but all samples in a comparison should use one fixed mode.
- `-t 24` sets the number of mapping worker threads; minibwa also uses a small number of I/O threads.
- `-R` writes a SAM read group. Bash `$'...'` quoting converts each `\t` into a literal tab.
- `-N 20` retains up to 20 candidate secondary alignments, and `--outn=20` writes those secondary records to SAM for inspecting multimapping among homologous genes. Omit both if only primary and supplementary records are needed; `--outn` defaults to zero.
- minibwa writes SAM by default. Do not add `-f`, which switches the output to PAF and cannot be piped directly to samtools or featureCounts.

minibwa is suitable for short metagenomic DNA reads and for prokaryotic metatranscriptomic reads when the reference genes contain no introns. It does not support spliced alignment and therefore must not be used to map eukaryotic RNA-seq directly to an intron-containing `batch.nr.gene.fasta`.

### C. Splice-Aware minimap2 for Eukaryotic Metatranscriptomes

Eukaryotic `batch.nr.gene.fasta` sequences retain introns, so RNA reads require a splice-aware aligner. Recent minimap2 releases provide the short-RNA preset `splice:sr`:

```bash
paftools.js gff2bed batch.nr.gene.gff3 \
  > batch.nr.gene.junctions.bed

minimap2 -x splice:sr -d batch.nr.gene.splice.mmi \
  batch.nr.gene.fasta

minimap2 -ax splice:sr -j batch.nr.gene.junctions.bed -t 24 \
  batch.nr.gene.splice.mmi \
  sample01_RNA_R1.fastq.gz sample01_RNA_R2.fastq.gz \
  | samtools sort -@ 8 -o sample01.RNA.nr_gene.bam

samtools index sample01.RNA.nr_gene.bam
```

`-j` supplies known junctions converted from GFF3. `splice:sr` writes long reference gaps as SAM CIGAR `N` operations. Run `minimap2 --help` to confirm that the installed release supports `splice:sr` and `-j`; otherwise upgrade or use a splice-aware tool such as STAR or HISAT2 with GFF3/GTF-derived junctions. minibwa and `minimap2 -x sr` cannot correctly align reads spanning exon junctions.

### D. featureCounts for Eukaryotic DNA Abundance and RNA Expression

For metagenomic DNA, count the complete genomic gene interval, including introns:

```bash
featureCounts \
  -T 16 -F GTF \
  -a batch.nr.gene.gff3 \
  -o sample01.DNA.gene_counts.tsv \
  -t gene -g gene_id \
  -p --countReadPairs \
  -s 0 -Q 10 \
  sample01.DNA.nr_gene.minibwa.bam
```

For metatranscriptomic RNA, aggregate exon counts by `gene_id`:

```bash
featureCounts \
  -T 16 -F GTF \
  -a batch.nr.gene.gff3 \
  -o sample01.RNA.gene_counts.tsv \
  -t exon -g gene_id \
  -p --countReadPairs \
  -s 0 -Q 10 \
  sample01.RNA.nr_gene.bam
```

Important parameters:

- `-t gene` counts the complete genomic gene for DNA. `-t exon -g gene_id` aggregates all exons of each gene for RNA.
- `-p --countReadPairs` counts a paired-end fragment once instead of counting its two reads independently.
- `-s 0/1/2` means unstranded, forward-stranded, or reverse-stranded. Select the value from the RNA library design rather than assuming every library is unstranded.
- `-Q 10` is an example minimum mapping quality, not a universal threshold. Evaluate its effect when many homologous genes remain.
- featureCounts excludes multimapping reads by default. Add `-M --fraction` only when fractional allocation is scientifically appropriate, and use and report the same rule for every sample.
- These commands omit `-B`. That option requires both mates to align and can discard additional fragments at the exact ends of gene-only references. Add it only when strict paired-end assignment is required.

featureCounts produces raw counts. Supply raw counts to count-based differential-expression models such as DESeq2 or edgeR; do not convert them to TPM first. For descriptive abundance, RPKM or TPM can be calculated from effective gene or exon length and library size.

### E. CoverM for Prokaryotic Gene-Level DNA or RNA Abundance

Every sequence in the nonredundant gene FASTA is one reporting unit. First use the minibwa workflow above to produce a reference-sorted and indexed BAM, then pass it to `coverm contig`:

```bash
coverm contig \
  --bam-files sample01.DNA.nr_gene.minibwa.bam \
  -t 24 \
  -m count mean covered_fraction rpkm tpm \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --contig-end-exclusion 0 \
  -o sample01.coverm.gene_abundance.tsv
```

This workflow can be used for metagenomic DNA and for intron-free prokaryotic metatranscriptomes. Run minibwa on the corresponding DNA or RNA FASTQ files and substitute the resulting BAM. CoverM cannot currently invoke minibwa as its internal `-p/--mapper`, so provide the minibwa BAM with `--bam-files`; do not use `-p minibwa`. `--min-read-percent-identity 95` and `--min-read-aligned-percent 75` require 95% identity and alignment of 75% of the read length. Validate these thresholds for the read quality, expected sequence divergence, and dereplication threshold of the project.

**Keep `--contig-end-exclusion 0`.** In `contig` mode, CoverM normally excludes bases at both reference ends from coverage calculations. When each reference sequence is one short gene, the default end exclusion can severely underestimate or eliminate coverage. `covered_fraction` can help reject hits restricted to one locally conserved domain. TPM and RPKM are relative to the current nonredundant reference and the alignments that pass filtering; they are not absolute proportions of the complete community, including unmapped organisms.

Recent CoverM releases also support per-feature reporting with `--gff`. This workflow already makes each gene an independent reference, so direct contig-level reporting is simpler and avoids reporting `gene`, `exon`, and `CDS` features as separate overlapping units.

## Boundary Alignment and Quantification Bias

1. **Exact gene boundaries lose fragments that extend beyond the model.** `sample.gene.fasta` contains no upstream or downstream flank. If one mate lies inside the gene and the other lies outside, the external mate cannot align completely. Strict proper-pair filtering or `featureCounts -B` discards additional fragments. The benefit is that flanking genes and noncoding regions cannot be assigned incorrectly to the target gene.
2. **Local alignment does not restore flanking sequence.** minibwa and minimap2 may soft-clip at a reference end, allowing part of a boundary-spanning read to remain aligned, but clipped bases do not contribute gene coverage. Keep the aligner, preset, identity, aligned-percent, MAPQ, and multimapping rules identical among samples.
3. **Eukaryotic RNA requires junction alignment.** Gene FASTA retains introns. Reads spanning exon junctions must be represented by `N` CIGAR operations from a splice-aware aligner. Ordinary short-DNA presets systematically underestimate multiexon genes.
4. **Protein-guided boundaries are not UTR boundaries.** miniprot primarily predicts coding models, so gene FASTA usually lacks complete 5-prime and 3-prime UTRs. Interpret RNA abundance over exon or CDS features rather than assuming that complete transcriptional regulatory regions were extracted.
5. **Homologous genes still cause multimapping.** Dereplication cannot eliminate every paralog. Ignoring, uniquely assigning, or fractionally assigning multimapping reads produces different abundance estimates. Preserve `.clstr` files and document the selected strategy.
6. **Assemblies can truncate genes.** Models at contig ends may be incomplete. CD-HIT's longest-first representative order helps retain complete candidates but cannot repair assembly errors. Consider query coverage, frameshifts, stop codons, and covered fraction during filtering.
7. **Do not mix coordinate systems.** A BAM mapped to `batch.nr.gene.fasta` must be quantified with `batch.nr.gene.gff3`. A BAM mapped to the original contigs must use `best_loci.gff3`. Their SeqIDs and coordinates are different.

## Checkpointing and Troubleshooting

`--resume` is enabled by default. A stage is reused only when it is marked complete, its output still exists, and file sizes match the recorded values. Changes to input path, size, or time; the reference-protein SHA-256 digest; result parameters; program version; or external-tool versions change the signature and rebuild incompatible results.

Force a complete rerun with:

```bash
python3 meta_homologous_gene_annot.py \
  -p references/target.faa \
  -c assemblies/sample01.fasta.gz \
  -o results/sample01 \
  --sample sample01 \
  --organism_type eukaryote \
  --force
```

Common checks:

```bash
# Summarize failure reasons
awk -F '\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="fail_reason") c=i; next}
              {print $c}' \
  results/sample01/sample01.all_hits.tsv | sort | uniq -c

# Count GFF3 feature types
awk '!/^#/ {n[$3]++} END {for (k in n) print k,n[k]}' \
  results/sample01/sample01.best_loci.gff3

# Check gene lengths in the gene-coordinate GFF3
awk '!/^#/ && $3=="gene" {print $1,$5-$4+1}' \
  results/sample01/sample01.gene.gff3 | sort -u

# Inspect external-tool errors
tail -n 50 results/sample01/sample01.miniprot.map.log
tail -n 50 results/sample01/sample01.gffread.log
```

Do not run two processes with the same `--outdir` and `--sample` simultaneously. If gffread fails, first confirm that GFF3 SeqIDs match the hit-contig FASTA IDs. No hits or no loci passing the filters is a valid result; the program writes empty gene/sequence files and complete summaries. Inspect `fail_reason` in `all_hits.tsv` before changing thresholds, and do not lower thresholds only to obtain a positive result.

## References

- [miniprot README and output documentation](https://github.com/lh3/miniprot)
- [miniprot manual](https://lh3.github.io/miniprot/miniprot.html)
- [CD-HIT User's Guide](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki)
- [minibwa README and command documentation](https://github.com/lh3/minibwa)
- [minimap2 README and manual](https://github.com/lh3/minimap2)
- [featureCounts examples](https://subread.sourceforge.net/featureCounts.html)
- [CoverM README](https://github.com/wwood/CoverM)
- [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

## Changelog

### 1.1.0

- Added `--organism_type {eukaryote,prokaryote}`; prokaryotic mode passes miniprot `-S` to disable splicing.
- Preserved selected miniprot `##PAF` alignment details in raw and filtered GFF3 output.
- Added a normalized `gene/mRNA/exon/CDS/intron` hierarchy and stable locus IDs to filtered GFF3.
- Added `sample.gene.fasta` and the independent gene-coordinate `sample.gene.gff3`.
- Added detailed CD-HIT cross-sample dereplication, minibwa/minimap2 mapping, featureCounts/CoverM quantification, and boundary-bias workflows to this README.

### 1.0.0

- Added reference cleaning, miniprot alignment, filtering, locus collapsing, GFF3/sequence export, and checkpointed resume support.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
