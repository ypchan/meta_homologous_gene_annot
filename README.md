# Meta Homologous Gene Annotator

`meta_homologous_gene_annot.py` 用于将参考蛋白批量比对到一份宏基因组 contig 组装结果，筛选可信的蛋白编码基因模型，合并冗余位点，并导出注释表、GFF3、命中 contig、CDS、蛋白和转录本序列。

管线以 [miniprot](https://github.com/lh3/miniprot) 完成剪接感知的蛋白到基因组比对，以 gffread 提取序列，以 pigz 并行压缩结果。程序支持阶段级断点续跑，并记录输入、参数、工具版本、运行时间和资源占用。

> 本工具发现的是与参考蛋白同源的候选基因。命中某个病原相关蛋白并不能单独证明 contig 来自对应病原物，也不能替代分类学、致病性或实验验证。

## 工作流程

1. 读取并清洗参考蛋白，为每条记录分配唯一内部 ID。
2. 为 contig FASTA 构建 miniprot 索引。
3. 使用 miniprot 搜索参考蛋白的候选基因模型。
4. 按一致性、查询覆盖度、移码和内部终止密码子筛选命中。
5. 在相同 contig、相同链上按区间重叠聚类，每个位点保留最佳模型。
6. 流式提取所有命中 contig，避免把整份组装载入内存。
7. 使用 gffread 导出 CDS、蛋白和转录本 FASTA，并使用 pigz 压缩。
8. 写出汇总、运行元数据、阶段指标和断点状态。

## 环境要求

- Linux 或其他可运行下列命令行工具的类 Unix 环境
- Python 3.9 或更高版本
- Python 包：`rich`、`rich-argparse`
- 外部程序：`miniprot`、`gffread`、`pigz`

若指定 `--skip_sequence_export`，可以不安装 gffread；miniprot 和 pigz 始终需要。

## 安装

### 方法一：Conda/Mamba（推荐）

克隆或下载仓库后，在仓库根目录执行：

```bash
mamba env create -f environment.yml
mamba activate meta-homologous-gene-annot
python3 meta_homologous_gene_annot.py --help
```

没有 `mamba` 时，可将第一条命令中的 `mamba` 换成 `conda`。

### 方法二：已有 Python 环境

先安装 miniprot、gffread 和 pigz，并确保三者位于 `PATH`，然后安装 Python 依赖：

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -U pip
python3 -m pip install -r requirements.txt
```

检查环境：

```bash
python3 --version
miniprot --version
gffread --version
pigz --version
python3 meta_homologous_gene_annot.py --help
```

外部程序不在 `PATH` 时，无需修改脚本，可分别通过 `--miniprot`、`--gffread` 和 `--pigz` 传入可执行文件路径。

## 输入文件

### 参考蛋白：`--proteins/-p`

输入为蛋白 FASTA，可为未压缩文件或 `.gz`、`.bgz`、`.bgzf` 文件。

```text
>PHI:1234 gene=MEP4 pathogen=Trichophyton_mentagrophytes
MKFSLALALAVASASA...
```

处理规则如下：

- 标题中第一个空白分隔字段作为 `original_id`，完整标题保存在结果中。
- 允许原始 ID 重复；程序会按输入顺序生成 `PHIREF000000001` 形式的唯一 `query_id`。
- 序列统一转为大写，移除空白、`-`、`.` 和末端 `*`。
- 支持的氨基酸字符为 `ABCDEFGHIKLMNPQRSTVWXYZUO`；其他字符和内部 `*` 替换为 `X`，替换数量写入映射表。
- 空序列会跳过；若没有任何可用蛋白，程序终止。

### Contig 组装：`--contigs/-c`

输入为核酸 FASTA，同样支持未压缩文件或 `.gz`、`.bgz`、`.bgzf` 文件。标题中第一个空白分隔字段作为 contig ID。ID 必须唯一，并且必须与 miniprot GFF3 中的序列 ID 一致。

程序流式扫描 FASTA；只有命中 ID 集合保存在内存中。若 GFF3 中的命中 ID 在 FASTA 中缺失，或命中的 contig ID 在 FASTA 中重复，提取阶段会报错并停止。

## 快速开始

```bash
python3 meta_homologous_gene_annot.py \
  --proteins data/phi-base_current.fas \
  --contigs data/201704_MF1.fasta.gz \
  --outdir results/201704_MF1 \
  --sample 201704_MF1 \
  --threads 24
```

`--sample` 可省略；此时样本名从 contig 文件名推断，并依次移除 `.gz/.bgz/.bgzf` 和 `.fasta/.fna/.fa/.fas` 后缀。样本名中非字母、数字、点、下划线或连字符的字符会替换为 `_`。

查看命令行帮助和内置详细说明：

```bash
python3 meta_homologous_gene_annot.py --help
python3 meta_homologous_gene_annot.py --help_input
python3 meta_homologous_gene_annot.py --help_default
```

## 参数详解

比例参数均使用 `0` 到 `1` 的小数，而不是百分数。例如 `0.40` 表示 40%。

### 必需参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `-p, --proteins FASTA` | 无 | 参考蛋白 FASTA；支持 gzip/bgzip 压缩。文件必须存在且非空。 |
| `-c, --contigs FASTA` | 无 | 单个样本的 contig FASTA；支持 gzip/bgzip 压缩。文件必须存在且非空。 |
| `-o, --outdir DIR` | 无 | 输出目录。不存在时自动创建；主结果直接写在该目录中。 |

### 运行与输出参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--sample NAME` | 从 contig 文件名推断 | 输出文件前缀和 GFF3 locus ID 使用的样本名。 |
| `-t, --threads INT` | 调度器分配值或 `min(CPU, 32)` | miniprot 线程数。程序依次读取 `SLURM_CPUS_PER_TASK`、`NSLOTS`、`OMP_NUM_THREADS`；均不可用时取系统 CPU 数且最多 32。必须大于等于 1。 |
| `--compression_threads INT` | `min(--threads, 8)` | pigz 压缩线程数，必须大于等于 1。 |
| `--tmpdir DIR` | `$SLURM_TMPDIR` 或 `OUTDIR/.work` | 临时父目录。每次运行在其下建立 `<sample>.<PID>`，成功后删除。此参数不参与断点签名。 |
| `--resume / --no-resume` | `--resume` | 启用或禁用兼容阶段的断点复用。禁用时重置状态并重跑；已有结果按阶段覆盖。 |
| `--force` | 关闭 | 删除该样本已知的输出和检查点后完整重跑；保留并继续追加主运行日志。 |
| `--keep_index` | 关闭 | 成功完成后保留 `<sample>.miniprot.mpi`。默认删除索引以节省空间。 |
| `--keep_uncompressed` | 关闭 | 除 `.gz` 外，同时保留命中 contig、CDS、蛋白和转录本的未压缩 FASTA。 |
| `--skip_sequence_export` | 关闭 | 不运行 gffread，不生成 CDS、蛋白和转录本 FASTA；命中 contig 仍会导出并压缩。 |

### Miniprot 比对参数

| 参数 | 默认值 | 传递方式与作用 |
| --- | --- | --- |
| `--splice_model {0,1,2}` | `1` | 传给 miniprot `-j` 的剪接模型。`1` 为通用模型，适合默认的真菌同源搜索；其他取值的语义以所安装 miniprot 版本帮助为准。 |
| `--max_intron BP` | `20000` | 传给 `-G`，允许的最大内含子长度，必须大于等于 1。物种内含子较长时需相应调大。 |
| `--index_subsample INT` | `1` | 建索引时传给 `-M`；k-mer 采样率为 `1/2**M`。必须大于等于 0；增大可减小索引/加快搜索，但可能降低灵敏度。 |
| `--max_hits INT` | `50` | 同时传给 `-N` 和 `--outn`，限制每条查询保留/输出的候选比对数。必须大于等于 1。 |
| `--secondary_ratio FLOAT` | `0.50` | 传给 `-p`；相对最佳分数达到该比例的次级命中才保留。 |
| `--prefilter_query_coverage FLOAT` | `0.30` | 传给 `--outc` 的查询覆盖度预过滤阈值。它发生在本程序最终过滤之前。 |
| `--min_score_ratio FLOAT` | `0.50` | 传给 `--outs`；输出分数至少为最佳比对分数该比例的命中。 |
| `--miniprot_extra 'OPTIONS'` | 空 | 添加到 miniprot 映射命令中的原始选项，经 shell 风格拆词后传递。仅用于脚本尚未暴露的 miniprot 选项；不要在其中加入索引或蛋白输入路径。 |

### 最终命中过滤与位点参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--min_identity FLOAT` | `0.40` | 最低氨基酸一致性。低于阈值的记录标记为 `low_identity`。 |
| `--min_query_coverage FLOAT` | `0.60` | 最低查询蛋白覆盖度，计算为 miniprot `Target` 区间长度除以清洗后参考蛋白长度。低于阈值标记为 `low_query_coverage`。 |
| `--max_frameshift INT` | `1` | 允许的最大移码数；必须大于等于 0。超过阈值标记为 `too_many_frameshifts`。 |
| `--max_stop_codon INT` | `0` | 允许的最大内部终止密码子数；必须大于等于 0。超过阈值标记为 `internal_stop_codon`。 |
| `--locus_overlap FLOAT` | `0.80` | 在相同 contig 和相同链上合并冗余参考命中的阈值：`交叠长度 / 较短区间长度`。 |
| `--high_identity FLOAT` | `0.60` | 高置信度标签的一致性阈值，不能低于 `--min_identity`。 |
| `--high_query_coverage FLOAT` | `0.80` | 高置信度标签的查询覆盖度阈值，不能低于 `--min_query_coverage`。 |

最终通过条件为四项最低/最高阈值同时满足。高置信度还要求 `identity >= high_identity`、`query_coverage >= high_query_coverage`、`frameshift == 0` 且 `stop_codon == 0`。通过最终过滤但不满足高置信度规则的位点标为 `medium`；未通过的原始命中标为 `low`，只保留在 `all_hits.tsv` 中。

同一候选位点存在多个参考命中时，程序依次按置信等级、`identity × query_coverage`、覆盖度、一致性、分数、较少移码、较少终止密码子和较高排名选出代表模型。其余参考 ID 写入 `alternative_references`。

### 外部程序与帮助参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--miniprot PATH` | `miniprot` | miniprot 命令名或可执行文件路径。 |
| `--gffread PATH` | `gffread` | gffread 命令名或可执行文件路径。使用 `--skip_sequence_export` 时不检查。 |
| `--pigz PATH` | `pigz` | pigz 命令名或可执行文件路径。 |
| `-h, --help` | — | 显示完整命令行参数并退出。 |
| `--help_input` | — | 显示输入格式、规则和示例并退出，不要求提供必需参数。 |
| `--help_default` | — | 显示默认阈值及解释并退出，不要求提供必需参数。 |

## 常用配置示例

提高灵敏度、接受较远同源蛋白（结果也会更多，需加强后续核查）：

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/sensitive \
  --min_identity 0.30 \
  --min_query_coverage 0.50 \
  --secondary_ratio 0.30 \
  --min_score_ratio 0.30
```

仅生成比对、表格、GFF3 和命中 contig，不导出基因序列：

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/tables_only \
  --skip_sequence_export
```

在调度系统的本地临时盘上运行并保留索引：

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/sample01 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --tmpdir "${SLURM_TMPDIR}" \
  --keep_index
```

## 输出文件

假设样本名为 `sample01`，所有文件均位于 `--outdir`。

| 文件 | 内容 |
| --- | --- |
| `sample01.summary.tsv` | 单行核心汇总：输入、位点数、置信等级、命中 contig 数、线程、完成时间和签名。 |
| `sample01.best_loci.tsv` | 每个合并位点的代表模型，是下游分析的主要结果表。 |
| `sample01.best_loci.gff3` | 通过过滤并去冗余后的 mRNA/子特征；mRNA 增加 `Locus`、`ReferenceID`、`Confidence`、`QueryCoverage` 和 `ReferenceAnnotation`。 |
| `sample01.all_hits.tsv` | miniprot 输出的全部 mRNA 命中，包括未通过最终过滤的记录、`pass_filter` 和 `fail_reason`。 |
| `sample01.query_summary.tsv` | 每条参考蛋白的原始命中数、通过数、入选位点数和最佳命中；无命中时为 `unmapped`。 |
| `sample01.contig_summary.tsv` | 每条命中 contig 的位点数、高/中置信位点数、最佳指标和参考 ID 集合。 |
| `sample01.hit_contig_ids.txt` | 去重并排序后的命中 contig ID。 |
| `sample01.hit_contigs.fasta.gz` | 所有命中 contig 的完整序列。 |
| `sample01.genes.cds.fasta.gz` | gffread 导出的 CDS；使用 `--skip_sequence_export` 时不生成。 |
| `sample01.genes.protein.fasta.gz` | gffread 导出的翻译蛋白；使用 `--skip_sequence_export` 时不生成。 |
| `sample01.genes.transcript.fasta.gz` | gffread 导出的转录本；使用 `--skip_sequence_export` 时不生成。 |
| `sample01.reference.clean.faa` | 清洗并重新编号后的参考蛋白，供 miniprot 使用。 |
| `sample01.reference_map.tsv` | 内部 `query_id` 到原始 ID、完整标题、长度及替换残基数的映射。 |
| `sample01.miniprot.raw.gff3` | 未经本程序最终过滤的 miniprot GFF3。 |
| `sample01.run_metadata.json` | 程序版本、完整命令、输入签名、参数、工具版本、结果路径和汇总。 |
| `sample01.stage_metrics.tsv` | 各阶段开始/结束时间、墙钟时间、子进程 CPU 时间、最大 RSS 和阶段详情。 |
| `sample01.state.json` | 断点续跑状态及各阶段输出大小；不要在运行中手工修改。 |
| `sample01.run.log` | 主运行日志和实际执行的外部命令。 |
| `sample01.miniprot.index.log` | miniprot 建索引的标准错误日志。 |
| `sample01.miniprot.map.log` | miniprot 映射的标准错误日志。 |
| `sample01.gffread.log` | gffread 日志；跳过序列导出时不生成。 |
| `*.pigz.log` | 各 FASTA 压缩步骤的 pigz 日志。 |
| `sample01.done` | 成功完成标记，包含样本名、完成时间和运行签名。 |
| `sample01.miniprot.mpi` | miniprot 索引；默认在成功运行后删除，指定 `--keep_index` 才保留。 |

指定 `--keep_uncompressed` 后，还会保留上述四类压缩 FASTA 对应的不带 `.gz` 文件。

### 主要结果字段

`all_hits.tsv` 和 `best_loci.tsv` 的常用字段：

- `query_id` / `original_id`：内部唯一参考 ID / FASTA 中原始 ID。
- `model_id`：miniprot 生成的模型 ID；`locus_id` 是去冗余后的稳定输出 ID。
- `contig`, `start`, `end`, `strand`：基因模型的基因组位置。
- `locus_start`, `locus_end`：合并簇覆盖的整体范围。
- `score`, `rank`, `identity`, `positive`：miniprot 比对指标。
- `query_start`, `query_end`, `aligned_query_length`, `query_length`, `query_coverage`：查询蛋白覆盖信息。
- `frameshift`, `stop_codon`：miniprot 报告的移码和内部终止密码子计数。
- `confidence`：`high`、`medium` 或 `low`。
- `pass_filter`, `fail_reason`：是否通过最终过滤及失败原因；仅 `all_hits.tsv` 包含。
- `n_reference_hits`, `alternative_references`：同一合并位点中的参考命中数量和非代表参考 ID；仅 `best_loci.tsv` 包含。

## 断点续跑、覆盖与签名

默认启用 `--resume`。一个阶段只有在状态为完成、输出文件存在且文件大小与记录一致时才会复用。完成标记也兼容时，重复执行相同命令会立即退出。

运行签名包含：程序版本、参考蛋白路径/大小/修改时间/SHA-256、contig 路径/大小/修改时间、样本名、影响结果的参数以及外部工具版本。`--force`、`--resume`、`--tmpdir` 和 `--outdir` 不计入参数签名。签名不兼容时，程序会清理该样本的已知旧结果后重建。

需要无条件重跑时使用：

```bash
python3 meta_homologous_gene_annot.py \
  -p reference.faa -c assembly.fasta.gz -o results/sample01 \
  --force
```

不要让两个进程同时使用相同的 `--outdir` 和 `--sample`，否则它们会写入同一组状态和结果文件。

## 更新

使用 Git 获取代码更新：

```bash
git pull
mamba env update -f environment.yml --prune
mamba activate meta-homologous-gene-annot
python3 meta_homologous_gene_annot.py --help
```

若采用 `venv + pip` 安装，则更新 Python 依赖：

```bash
source .venv/bin/activate
python3 -m pip install -U -r requirements.txt
```

更新脚本或 miniprot/gffread 之后，建议对重要结果使用 `--force` 完整重跑。断点签名记录程序声明版本和外部工具版本，但不计算脚本文件本身的哈希；如果代码发生变化而 `PROGRAM_VERSION` 未变化，旧阶段仍可能被视为兼容。

当前脚本声明版本为 `1.0.0`。查看本地代码版本可使用 `git log -1 --oneline`；查看尚未提交的修改可使用 `git status`。

## 常见问题

### 提示 `Executable not found in PATH`

确认对应程序已安装且可直接执行，或者用 `--miniprot /绝对路径/miniprot` 等参数指定路径。使用 `--skip_sequence_export` 只能免除 gffread，不能免除 miniprot 或 pigz。

### 参数改变后为什么从头运行

影响结果的输入、参数或工具版本发生变化会使签名失效，程序会清理已知旧输出并重新计算。这是为了避免混用不兼容的阶段结果。

### 没有通过过滤的位点

这是合法结果。程序仍写出汇总和空的最终序列文件。先检查 `all_hits.tsv` 中的 `fail_reason`，再决定是否有生物学依据调整阈值；不要仅为了获得命中而盲目降低阈值。

### gffread 阶段失败

检查 `sample01.gffread.log`、`best_loci.gff3` 与命中 contig FASTA。常见原因是 contig ID 不一致、GFF3 与 FASTA 不匹配或 gffread 版本差异。修正后使用相同命令续跑；若需要重建所有阶段则加 `--force`。

## 版本记录

### 1.0.0

- 实现参考蛋白清洗、miniprot 比对、质量过滤和冗余位点合并。
- 导出 GFF3、结果表、命中 contig、CDS、蛋白和转录本。
- 支持 gzip/bgzip 输入、pigz 输出压缩、阶段级断点续跑和运行元数据记录。

## 许可证

本项目采用 MIT License，详见 [LICENSE](LICENSE)。
