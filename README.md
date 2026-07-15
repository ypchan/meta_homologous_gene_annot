# Meta Homologous Gene Annotator

`meta_homologous_gene_annot.py` 使用 miniprot 将参考蛋白比对到单个宏基因组组装，过滤并合并冗余候选位点，同时输出两套互相对应的 GFF3：

- contig 坐标：用于检查原始组装上的基因位置和蛋白比对；
- gene-reference 坐标：每条参考序列就是一个基因，可在多个样本间去冗余后只把 reads 比对到非冗余基因集，计算宏基因组丰度或宏转录组表达量。

程序支持真核和原核两种模式。真核模式进行剪接感知比对并输出 `exon`/`intron`；原核模式关闭 miniprot 剪接，不输出 `intron`。这些结果是基于同源蛋白的候选基因模型，不等同于物种鉴定、致病性证据或完整的从头基因注释。

## 核心输出设计

假定样本名为 `sample01`：

| 文件 | 坐标和用途 |
| --- | --- |
| `sample01.miniprot.raw.gff3` | 未经过最终过滤的 miniprot GFF3。`##PAF` 行保留 `cg`/`cs` 详细蛋白—基因组比对。 |
| `sample01.best_loci.gff3` | 筛选、位点合并后的 contig 坐标 GFF3，包含所选 `##PAF`、`gene`、`mRNA`、`exon`、`CDS`、`intron`（真核）和可能的 `stop_codon`。 |
| `sample01.gene.fasta` | 每个非冗余样本内位点一条序列，序列 ID 为 `sample01_PHILOCUS...`。这是跨样本去冗余和 reads 比对的输入。始终保留为未压缩 FASTA。 |
| `sample01.gene.gff3` | 与 `sample01.gene.fasta` 配套的独立基因坐标 GFF3。SeqID 等于 locus ID，基因坐标固定为 `1..gene_length`。 |
| `sample01.best_loci.tsv` | 每个位点的来源 contig、原始坐标、参考蛋白、identity、query coverage、置信度和备选参考。 |
| `sample01.all_hits.tsv` | miniprot 的全部 mRNA 命中，包括未通过过滤的记录和失败原因。 |
| `sample01.genes.cds.fasta.gz` | gffread 提取的拼接 CDS。 |
| `sample01.genes.transcript.fasta.gz` | gffread 提取的拼接转录本；真核跨样本去冗余时优先使用它聚类。 |
| `sample01.genes.protein.fasta.gz` | gffread 翻译的蛋白。 |
| `sample01.hit_contigs.fasta.gz` | 至少包含一个保留位点的完整 contig。 |

另外还会输出 query/contig 汇总表、清洗后的参考蛋白、参考 ID 映射、运行日志、工具版本、参数签名、阶段耗时、断点状态和完成标记。完整列表可从 `sample01.run_metadata.json` 查看。

### 两套坐标如何对应

若 contig 上的一个基因位于 `[S, E]`，GFF3 使用 1-based 闭区间：

```text
contig GFF3:  contig_7    gene    S            E
gene FASTA:  >sample01_PHILOCUS0000001，长度 E-S+1
gene GFF3:   locus_id     gene    1            E-S+1
坐标换算:    gene_pos = contig_pos - S + 1
```

`sample01.gene.fasta` 精确使用 Python 切片 `sequence[S-1:E]`。负链基因不做反向互补，而是保留 contig 原方向，并在两套 GFF3 中保留 `strand=-`。这样 FASTA、GFF3 和 BAM 坐标始终一致；minibwa/minimap2 本身会搜索双链。

外显子由 miniprot 的 CDS/stop-codon 区间合并得到；真核内含子是相邻外显子之间的区间。miniprot 是蛋白引导的模型，通常不预测 UTR，因此这里的“gene 边界”是预测的蛋白编码模型边界，不应解释为完整转录单元边界。

## 依赖和安装

运行主流程需要：

- Python 3.9+；
- Python 包 `rich`、`rich-argparse`；
- [miniprot](https://github.com/lh3/miniprot)、[gffread](https://github.com/gpertea/gffread)、[pigz](https://github.com/madler/pigz)。

`--skip_sequence_export` 可跳过 gffread 生成的 CDS/蛋白/转录本，但 `sample.gene.fasta` 和两套 GFF3 仍会生成。

推荐使用 conda/mamba 安装完整的注释和定量环境：

```bash
mamba create -n meta_gene \
  -c conda-forge -c bioconda \
  python=3.11 rich rich-argparse \
  miniprot gffread pigz cd-hit minibwa minimap2 samtools subread coverm
conda activate meta_gene
```

也可以只安装 Python 依赖：

```bash
python3 -m pip install --user --upgrade -r requirements.txt
```

从源码安装核心工具的命令如下：

```bash
mkdir -p "$HOME/src" "$HOME/.local/bin"

git clone https://github.com/lh3/miniprot.git "$HOME/src/miniprot"
make -C "$HOME/src/miniprot"
install -m 0755 "$HOME/src/miniprot/miniprot" "$HOME/.local/bin/miniprot"

git clone https://github.com/gpertea/gffread.git "$HOME/src/gffread"
make -C "$HOME/src/gffread" release
install -m 0755 "$HOME/src/gffread/gffread" "$HOME/.local/bin/gffread"

git clone https://github.com/madler/pigz.git "$HOME/src/pigz"
make -C "$HOME/src/pigz"
install -m 0755 "$HOME/src/pigz/pigz" "$HOME/.local/bin/pigz"

git clone https://github.com/lh3/minibwa.git "$HOME/src/minibwa"
make -C "$HOME/src/minibwa"
install -m 0755 "$HOME/src/minibwa/minibwa" "$HOME/.local/bin/minibwa"

export PATH="$HOME/.local/bin:$PATH"
```

安装后检查实际版本和本机帮助；不同版本的可用参数可能不同：

```bash
python3 --version
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

## 输入

### 参考蛋白 `--proteins/-p`

支持普通或 `.gz/.bgz/.bgzf` FASTA。标题第一个空白分隔字段作为 `original_id`，完整标题保留到结果。重复 original ID 可以存在，程序会按输入顺序生成唯一的 `PHIREF000000001` 等内部 ID。序列转大写，去除空白、`-`、`.` 和末端 `*`；无法识别的残基转为 `X` 并记录。

### 单样本 contig `--contigs/-c`

支持普通或压缩 FASTA。标题第一个字段是 contig ID，必须唯一。程序流式读取大 FASTA，不会把整个组装同时载入内存。

## 注释运行示例

真核或含内含子的候选基因：

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

原核基因：

```bash
python3 meta_homologous_gene_annot.py \
  --proteins references/target.faa \
  --contigs assemblies/sample01.fasta.gz \
  --outdir results/sample01 \
  --sample sample01 \
  --organism_type prokaryote \
  --threads 24
```

`--organism_type prokaryote` 会给 miniprot 使用 `-S` 关闭剪接；此时不要设置非零 `--splice_model`。真核默认使用 `--splice_model 1`，适合一般真核/真菌；miniprot 的 `2` 是 vertebrate/insect 模型，`0` 不使用剪接信号模型。

## 主要参数

所有比例均写成 `0..1` 小数而不是百分数。

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--organism_type` | `eukaryote` | `eukaryote`：剪接感知并输出 intron；`prokaryote`：miniprot `-S`，不输出 intron。 |
| `--splice_model` | 真核 `1` | miniprot `-j`；只用于真核。 |
| `--max_intron` | `20000` | 真核 miniprot `-G` 最大内含子长度。过小会截断模型，过大会增加错误长间隔和计算量。 |
| `--index_subsample` | `1` | 建索引时 miniprot `-M`；采样率为 `1/2**M`。 |
| `--max_hits` | `50` | 同时传给 miniprot `-N` 和 `--outn`。 |
| `--secondary_ratio` | `0.50` | miniprot `-p`，次优命中相对最佳得分阈值。 |
| `--prefilter_query_coverage` | `0.30` | miniprot `--outc` 的预过滤；程序之后还会执行更严格的最终过滤。 |
| `--min_score_ratio` | `0.50` | miniprot `--outs`，相对最佳比对得分阈值。 |
| `--min_identity` | `0.40` | 最终氨基酸 identity 下限。 |
| `--min_query_coverage` | `0.60` | 最终参考蛋白覆盖度下限。 |
| `--max_frameshift` | `1` | 允许的移码事件上限。 |
| `--max_stop_codon` | `0` | 允许的内部终止密码子上限。 |
| `--locus_overlap` | `0.80` | 同 contig、同链命中的 `overlap/shorter_interval` 达到该值时归为同一位点。 |
| `--high_identity` | `0.60` | high-confidence identity 下限。 |
| `--high_query_coverage` | `0.80` | high-confidence query coverage 下限；high 还要求零移码和零内部终止。 |
| `--miniprot_extra` | 空 | 追加 miniprot 参数，使用 shell quoting 解析；不要重复加入索引或蛋白路径。 |
| `--resume/--no-resume` | resume | 按输入、参数、工具版本签名复用完整阶段。 |
| `--force` | 关闭 | 清除本样本已知结果和 checkpoint 后重跑，日志文件继续追加。 |
| `--keep_uncompressed` | 关闭 | 额外保留 hit-contig/CDS/蛋白/转录本的未压缩文件。`sample.gene.fasta` 不受此参数影响。 |

运行 `--help_input` 和 `--help_default` 可查看内置说明：

```bash
python3 meta_homologous_gene_annot.py --help_input
python3 meta_homologous_gene_annot.py --help_default
```

## 程序实际调用的外部命令

下面展开主流程的关键命令，便于审计参数。实际绝对路径和完整命令也记录在 `sample01.run.log`。

建索引：

```bash
miniprot -t 24 -M 1 \
  -d sample01.miniprot.mpi \
  assemblies/sample01.fasta.gz
```

真核比对：

```bash
miniprot -t 24 -j 1 -G 20000 \
  -N 50 -p 0.50 \
  --outn 50 --outs 0.50 --outc 0.30 \
  --gff --gff-delim '@' \
  sample01.miniprot.mpi sample01.reference.clean.faa \
  > sample01.miniprot.raw.gff3
```

原核比对把 `-j/-G` 换成 `-S`：

```bash
miniprot -t 24 -S \
  -N 50 -p 0.50 \
  --outn 50 --outs 0.50 --outc 0.30 \
  --gff --gff-delim '@' \
  sample01.miniprot.mpi sample01.reference.clean.faa \
  > sample01.miniprot.raw.gff3
```

这里特意使用 `--gff` 而不是 `--gff-only`：每个模型前的 `##PAF` 行含 0-based、右端不包含的 PAF 坐标，以及 `cg:Z:` 蛋白 CIGAR 和 `cs:Z:` 差异字符串；普通 GFF3 feature 行仍使用 1-based 闭区间。`N/U/V` CIGAR 操作表示不同 phase 的内含子，`F/G` 表示移码。不要直接把 PAF 的 start/end 当成 GFF3 坐标。

序列导出：

```bash
gffread -E sample01.best_loci.gff3 \
  -g sample01.hit_contigs.fasta \
  -x sample01.genes.cds.fasta \
  -y sample01.genes.protein.fasta \
  -w sample01.genes.transcript.fasta
```

`-x/-y/-w` 分别输出 CDS、蛋白和拼接转录本。随后程序用 `pigz -p THREADS -c` 压缩它们。

## 跨样本去冗余：推荐流程

同一批次的不同样本可能组装出相同基因。所有 locus ID 都带 sample 前缀，所以可以直接合并，不会因普通 contig 名称重复而冲突。

### 1. 合并单样本 gene FASTA/GFF3

```bash
find results -name '*.gene.fasta' -print0 \
  | sort -z \
  | xargs -0 cat \
  > batch.all.gene.fasta

printf '##gff-version 3\n' > batch.all.gene.gff3
find results -name '*.gene.gff3' -exec awk '!/^#/' {} + \
  >> batch.all.gene.gff3
```

### 2A. 原核：直接对 genomic gene 去冗余

```bash
cd-hit-est \
  -i batch.all.gene.fasta \
  -o batch.nr.gene.fasta \
  -c 0.95 -n 10 \
  -G 1 -aS 0.90 -s 0.80 \
  -g 1 -r 1 -d 0 \
  -T 24 -M 0
```

参数含义：

- `-c 0.95`：至少 95% nucleotide identity；若只想去掉近乎完全相同序列，可改为 `0.99`；若目标是基因家族丰度，才考虑更低阈值。
- `-n 10`：适合 0.95–1.0 identity 的 nucleotide word length。
- `-G 1`：global identity，以较短序列全长为分母，避免短局部结构域轻易聚类。
- `-aS 0.90`：比对至少覆盖较短序列的 90%。
- `-s 0.80`：较短序列长度至少为代表序列的 80%；允许部分组装基因归入更完整的代表，同时限制极短片段。
- `-g 1`：把序列分到达到阈值的最相似 cluster，速度较慢但 cluster 归属更准确。
- `-r 1`：检查正向和反向互补。虽然本程序保留 contig 方向，这一参数仍应开启。
- `-d 0`：保留完整的第一个 FASTA ID，便于回连 GFF3。
- `-T 24`：线程数；`-M 0` 表示不设置 CD-HIT 内存上限，应按集群策略调整。

CD-HIT 会先按序列从长到短处理，cluster 的第一个代表因此通常是该 cluster 最长、最完整的候选。`-g 1` 改变成员归属，但不改变这一代表选择规则。`batch.nr.gene.fasta.clstr` 必须保留，它记录代表和所有成员的对应关系。

### 2B. 真核：推荐先对拼接 transcript 去冗余

CD-HIT 官方说明也指出，长内含子会使真核 genomic gene 难以进行全长聚类。不同样本的同一编码基因可能有长度和序列差异很大的 intron，因此真核推荐对去除 intron 后的 transcript 聚类，再取代表 transcript 对应的 genomic gene：

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

这里假定本程序生成的单转录本 ID 为 `<locus_id>.t1`。如果后续工具改写了 FASTA ID，应先检查 `grep '^>'` 的实际内容，不要盲目删除后缀。

### 3. 生成非冗余 gene GFF3

若原核是直接对 gene FASTA 聚类，先提取代表 ID；真核上一步已经得到该文件：

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

必须验证 FASTA 和 GFF3 的 SeqID 完全对应：

```bash
grep '^>' batch.nr.gene.fasta \
  | sed 's/^>//; s/[[:space:]].*//' \
  | sort -u > fasta.ids

awk '!/^#/ {print $1}' batch.nr.gene.gff3 \
  | sort -u > gff.ids

diff -u fasta.ids gff.ids
```

`diff` 无输出才表示匹配。去冗余阈值决定最终统计单位：高阈值接近 allele/strain-level，较低阈值更接近 gene-family-level；所有样本必须使用同一个非冗余参考和同一套比对/多重比对规则，结果才可比较。

## reads 比对和丰度统计

### A. minimap2：宏基因组短 reads

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

`-x sr` 是短 DNA reads preset；`-a` 输出 SAM；`--secondary=yes -N 20` 保留候选次优位置，便于明确处理同源基因的 multi-mapping。不要在不同样本使用不同的 `-N` 或 MAPQ 规则。

### B. minibwa：宏基因组短 reads

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

minibwa 的参数含义：

- `index -t 24 REF PREFIX` 使用 24 线程建立 `PREFIX.l2b` 和 `PREFIX.mbw`。默认 libsais 建索引速度快，但大约需要参考序列长度 18 倍的内存；内存不足时使用 `minibwa index -l REF PREFIX`，代价是建索引更慢。
- `map -x sr` 明确使用 short-read preset；不指定 `-x` 时是 adaptive 模式，也能根据 read 长度调整参数，但同一批样本应固定一种模式。
- `-t 24` 是比对 worker 线程数；minibwa 还会使用少量 I/O 线程。
- `-R` 写入 SAM read group。这里使用 Bash 的 `$'...'`，把 `\t` 转成真正的制表符。
- `-N 20` 最多保留 20 个候选 secondary alignment；`--outn=20` 把这些 secondary records 实际写入 SAM，便于检查同源基因的 multi-mapping。若只要 primary/supplementary records，可省略二者，`--outn` 默认是 0。
- minibwa 默认输出 SAM；不要加入 `-f`，因为 `-f` 会改成 PAF，不能直接交给 samtools/featureCounts。

minibwa 可用于宏基因组短 DNA reads，也可用于没有内含子的原核宏转录组 reads。它不支持 spliced alignment，因此不能用来把真核 RNA-seq 比对到保留 intron 的 `batch.nr.gene.fasta`。

### C. 真核宏转录组：splice-aware minimap2

真核 `batch.nr.gene.fasta` 含 intron，RNA reads 必须进行剪接感知比对。较新的 minimap2 提供短 RNA preset `splice:sr`：

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

`-j` 使用 GFF3 转出的已知 junction；`splice:sr` 会把长参考缺口写成 SAM CIGAR `N`。先运行 `minimap2 --help` 确认安装版本支持 `splice:sr` 和 `-j`；旧版本请升级，或使用支持 GFF3/GTF junction 的 STAR/HISAT2。直接用 minibwa 或 `minimap2 -x sr` 会使跨 exon junction 的 reads 无法正确比对。

### D. 真核 featureCounts：DNA 丰度和 RNA 表达量

宏基因组 DNA 应按整个 gene 区间计数，包括 intron：

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

宏转录组 RNA 按 exon 汇总到 `gene_id`：

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

参数要点：

- `-t gene` 对 DNA 统计完整 genomic gene；`-t exon -g gene_id` 对 RNA 把多个 exon 合并为一个 gene。
- `-p --countReadPairs` 以 paired-end fragment 而不是两个 reads 计数。
- `-s 0/1/2` 分别表示非链特异、正向链特异、反向链特异。宏转录组必须根据建库方法选择，不能默认所有库都是 `0`。
- `-Q 10` 是最低 mapping quality 示例，不是通用真理；同源基因多时应比较不同阈值的影响。
- featureCounts 默认不计 multi-mapping reads。若研究问题需要分摊，可额外使用 `-M --fraction`，但必须对所有样本保持一致，并报告这一规则。
- 上述命令没有加 `-B`。`-B` 要求 paired fragment 两端都成功比对，在精确截断的 gene-only 参考上会额外丢弃跨基因边界的片段。只有需要严格双端比对时才加。

featureCounts 输出是 raw count。用于差异表达时应把 raw count 输入 DESeq2/edgeR 等模型；不要先把 count 转成 TPM 再做 count-based 差异检验。若只做描述性丰度，可按有效 gene/exon 长度和 library size 计算 RPKM/TPM。

### E. 原核 CoverM：DNA 或 RNA 的 gene-level 丰度

非冗余 gene FASTA 中每条序列就是一个统计单元。先用上面 minibwa 流程生成按参考坐标排序并建立索引的 BAM，再交给 `coverm contig`：

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

宏基因组和无内含子的原核宏转录组都可用这条命令：分别以对应 DNA/RNA FASTQ 运行 minibwa，并替换这里的 BAM。CoverM 当前不能把 minibwa 作为内部 `-p/--mapper` 直接调用，因此这里明确输入 minibwa 生成的 `--bam-files`，不能写成 `-p minibwa`。`--min-read-percent-identity 95` 和 `--min-read-aligned-percent 75` 分别要求 95% identity 和 75% read 长度完成比对，应根据 read 质量、物种差异和去冗余阈值验证。

**必须注意 `--contig-end-exclusion 0`。** CoverM 的 `contig` 模式默认会从每条参考序列两端排除一定碱基再计算覆盖度；当“每条 contig 就是一条短基因”时，默认端部排除会显著低估甚至清空短基因覆盖。`covered_fraction` 可以用于排除只有局部保守结构域有 reads 的假阳性。TPM/RPKM 是相对于当前非冗余 gene reference 和通过过滤的比对计算的，不代表未映射群落部分的绝对比例。

CoverM 新版本也支持 `--gff` 的 per-feature 统计；本流程已经把每条 gene 变成独立参考，因此直接按 contig 统计更简单，并避免一个 gene GFF3 中的 gene/exon/CDS 多种 feature 被重复报告。

## 边界比对和定量偏差：必须阅读

1. **精确 gene 边界会丢失跨边界片段。** `sample.gene.fasta` 不包含上下游 flank。一个 mate 在 gene 内、另一个 mate 在 gene 外时，外侧 mate 不能完整比对；严格 proper-pair/`featureCounts -B` 会进一步丢失该 fragment。优点是不会把相邻基因或非编码区覆盖错误归给目标基因。
2. **局部比对不是增加 flank。** minibwa/minimap2 可在参考末端 soft-clip，因此部分跨界 read 仍可保留，但 clipped 部分不提供 gene 内覆盖。要比较样本，必须固定 aligner、preset、identity、aligned-percent、MAPQ 和 multi-mapping 规则。
3. **真核 RNA junction。** gene FASTA 保留 intron；跨 exon junction 的 RNA read 必须由 splice-aware aligner 产生 `N` CIGAR。普通短 DNA preset 会系统性低估多 exon 基因。
4. **蛋白引导边界不是 UTR 边界。** miniprot 主要给出编码区模型，gene FASTA 通常不含完整 5'/3' UTR。RNA 表达量应按 exon/CDS解释，而不是假定提取了完整转录本调控区。
5. **同源基因 multi-mapping。** 去冗余不能消除所有 paralog。忽略、唯一分配或 fractional 分配会得到不同丰度；应保存 `.clstr`，并在项目方法中写明策略。
6. **组装截断。** 位于 contig 边缘的模型可能天然不完整。CD-HIT 从长到短选择代表有助于保留完整候选，但不能修复错误组装；应结合 query coverage、frameshift、stop codon 和 covered fraction 过滤。
7. **不要混用坐标。** reads 对 `batch.nr.gene.fasta` 的 BAM 必须使用 `batch.nr.gene.gff3`；reads 对原始 contig 的 BAM 必须使用 `best_loci.gff3`。两者 SeqID 和坐标不同。

## 断点续跑和故障排查

`--resume` 默认开启。只有阶段状态为 completed、输出仍存在且文件大小一致时才复用。输入路径/大小/时间、参考蛋白 SHA-256、结果参数、程序版本或外部工具版本变化会改变签名并重建不兼容结果。

强制重跑：

```bash
python3 meta_homologous_gene_annot.py \
  -p references/target.faa \
  -c assemblies/sample01.fasta.gz \
  -o results/sample01 \
  --sample sample01 \
  --organism_type eukaryote \
  --force
```

常见检查：

```bash
# 查看未通过原因
awk -F '\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="fail_reason") c=i; next}
              {print $c}' \
  results/sample01/sample01.all_hits.tsv | sort | uniq -c

# 检查 GFF3 feature 数量
awk '!/^#/ {n[$3]++} END {for (k in n) print k,n[k]}' \
  results/sample01/sample01.best_loci.gff3

# 检查 gene FASTA 长度是否等于 gene GFF3 长度
awk '!/^#/ && $3=="gene" {print $1,$5-$4+1}' \
  results/sample01/sample01.gene.gff3 | sort -u

# 查看外部工具错误
tail -n 50 results/sample01/sample01.miniprot.map.log
tail -n 50 results/sample01/sample01.gffread.log
```

不要在同一个 `--outdir/--sample` 上同时运行两个进程。若 gffread 失败，首先核对 GFF3 SeqID 与 hit-contig FASTA ID 是否一致。没有命中或没有位点通过过滤是合法结果，程序会生成空 gene/sequence 文件和完整汇总；应先检查 `all_hits.tsv` 的 `fail_reason`，不要只为得到阳性结果而降低阈值。

## 参考文档

- [miniprot README 和输出说明](https://github.com/lh3/miniprot)
- [miniprot man page](https://lh3.github.io/miniprot/miniprot.html)
- [CD-HIT User's Guide](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki)
- [minibwa README 和命令说明](https://github.com/lh3/minibwa)
- [minimap2 README/man page](https://github.com/lh3/minimap2)
- [featureCounts 官方示例](https://subread.sourceforge.net/featureCounts.html)
- [CoverM README](https://github.com/wwood/CoverM)
- [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

## Changelog

### 1.1.0

- 新增 `--organism_type {eukaryote,prokaryote}`；原核模式用 miniprot `-S` 关闭剪接。
- 原始和筛选 GFF3 保留所选 miniprot `##PAF` 详细比对。
- 筛选 GFF3 新增标准化 `gene/mRNA/exon/CDS/intron` 层级和稳定 locus ID。
- 新增 `sample.gene.fasta` 和独立 gene 坐标 `sample.gene.gff3`。
- README 新增 CD-HIT 跨样本去冗余、minibwa/minimap2、featureCounts、CoverM 和边界处理流程。

### 1.0.0

- 初始的参考清洗、miniprot 比对、过滤、位点合并、GFF3/序列导出和断点续跑。

## License

本项目采用 MIT License，见 [LICENSE](LICENSE)。
