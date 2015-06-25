---
layout: post
title: "估计RNA-seq转录本表达量和寻找差异表达基因"
date: 2015-06-03 16:58:01 +0800
comments: true
published: false
categories: Bioinfor
---

## 1. 标准化和计数 ##

### 1.1 HTseq ###

<!--more-->

**简介**：[HTseq](http://www-huber.embl.de/users/anders/HTSeq/doc/index.html) 是用于Python平台写成的处理高通量测序的平台。`htseq-count`{:.language-bash}可以用来对原始转录本计数，具体计数规则参考[Counting reads in features with htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)。

**平台**：Python跨平台使用。

**快速运行**：

{% codeblock lang:bash An Example of Running htseq-count%}
# 查看帮助
$ htseq-count --help

# 对bam文件按照reads名称排序
$ samtools sort -n accepted_filtered.bam accepted_sortname

# 注意“=”前后无空格
$ htseq-count --mode=union --stranded=no --type=exon --idattr=gene_id \
	          --format=bam accepted_sort.bam hg19USCS_ensembl.gtf > htseqcount_accepted.hsc
{% endcodeblock %}

**重要参数解释**：

* `--mode`{:.language-bash$}：统计落在某个基因上的reads数目的模型，默认值为“union”（[图示](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)）。作者认为“union”方法在绝大多数情况下都有很好的表现，建议使用。

* `--stranded`{:.language-bash$}：测序方法，默认为“yes”。

* `--type`{:.language-bash$}：计数单元类型，默认为[GTF](http://mblab.wustl.edu/GTF22.html)文件的`exon`{:.language-bash}。

* `--idattr`{:.language-bash$}：计数单元归类，默认为`gene_id`{:.language-bash}。比如把合并汇报多个exon对应的一个gene。

* `--format`{:.language-bash$}：可以输入bam或者sam文件，bam文件需要制定此参数。

**补充**：

* 如果输入的是bam文件，bam文件需要按照reads名称排序，操作方法为：

{% codeblock lang:bash bam to name-sorted bam%}
$ samtools sort -n accepted_filerted.bam acceted_sortname
{% endcodeblock %}

### 1.2  <span style="color:blue">GenomicAlignments</span> ###

[GenomicAlignments](http://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)是R/Bioconductor的一个包，其中`summarizeOverlaps`{:.language-R}函数用于对alignment数据进行计数。


### 1.3 Kallisto ###

**简介**：[kallisto](http://pachterlab.github.io/kallisto/)是由[Lior Pachter小组](http://pachterlab.github.io/)开发一款快速测量RNA-seq数据中转录本表达丰度的软件。因为使用了*pseudoalignment*的想法，可以不用alignment，直接测量原始测序数据，因此极大提高了运算速度。因为kallisto极高的速度，可以使用[bootstrap](http://en.wikipedia.org/wiki/Bootstrapping_(statistics))精确估计的“不确定性（uncertainty）”，可以配合下游软件sleuth确定差异表达基因。

**平台**：Linux和Mac

**快速运行**：

首先，建立索引。

{% codeblock lang:bash An Example of Building Index%}
# 查看帮助
$ kallisto index 

# 建立索引
$ kallisto index -i hg19.idx hg19UCSC_ensembl.fa
{% endcodeblock %}

**重要参数解释**：

* `-i/--index`{:.language-bash$}：kallisto输出索引文件。

* * `-k/--kmer-size`{:.language-bash$}：k-mer的长度（奇数），最大值为31，默认值为31。

之后，根据原始RNA-seq数据和转录组注释，测量表达量。

{% codeblock lang:bash An Example of Quantification%}
$ kallisto quant -i hg19.idx -o kallOut -b 1000 human1.fq.gz human2.fq.gz
{% endcodeblock %}

**重要参数解释**：

* `-i/--index`{:.language-bash$}：kallisto输出索引文件。

* `-o/--output-dir`{:.language-bash$}：输出目录。

* `-b/--bootstrap-samples`{:.language-bash$}：bootstrap的次数，默认为0。


## 2. 差异表达基因筛选 ##


### 2.1 <span style="color:blue">edgeR</span> ###

**简介**：[edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)是使用RNA-seq的reads的原始计数（count）和[负二项分布](http://en.wikipedia.org/wiki/Negative_binomial_distribution)模型，计算差异表达基因的R/Bioconductor包。edgeR能够处理RNA-seq、ChIP-seq和SAGE数据。根据edgeR的[说明文件](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)数据前处理有两点要求：

> 1. edgeR performs differential abundance analysis for pre-defined genomic features. Although not strictly necessary, it usually desirable that these genomic features are non-overlapping.（edgeR对预定义的基因组特征数据进行差异量分析。尽管不是严格必须，但是输入的基因特征最好能无相互重叠。）
>
> 2. The first step in an RNA-seq analysis is usually to align the raw sequence reads to a reference genome, although there are many variations on this process.（RNA-seq分析的第一步通常是把原始测序reads比对到参考基因组，当然，这一步有多种实现方法。）

以上两点表明输入的数据：1. 分析某个基因组特征，比如基因、外显子、[CDS](http://en.wikipedia.org/wiki/Coding_region)；2. 原始reads的计数，而不是一个估计量、FPKM或者RPKM，这些估计量会影响edgeR计算一些模型涉及的关键变量。对于第二点，edgeR指出也能处理类似[RSEM](http://deweylab.biostat.wisc.edu/rsem/)输出的计数估计值，但原始计数最好。


### 2.2 <span style="color:blue">DESeq2</span> ###

**简介**：[DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html)是R/Bioconductor的包，使用负二项分布模型和RNA-seq的原始count数值，寻找差异表达基因。[DESeq2说明文件](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf)指出：

> The count values must be raw counts of sequencing reads. This is important for DESeq2 ’s statistical model to hold, as only the actual counts allow assessing the measurement precision correctly. Hence, please do not supply other quantities, such as (rounded) normalized counts, or counts of covered base pairs – this will only lead to nonsensical results.（count数值必须为测序结果的原始count。只有使用真实的count才能正取评估模型参数，这对于DESeq2所使用的统计模型至关重要。因此，请不要使用诸如估计的（取整）count或者覆盖碱基对的count，否则只能得到没有意义的结果。）

**平台**：R跨平台使用。

**运行流程**：

{% codeblock lang:R Workflow of DESeq2 %}
# 载入包并且注册多核
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(4))

# “glioPR”是一个整理好的数据框，行表示count数，列标识样本（1~10列）和注释（11~14列）。
# “targets”是一个数据框，表示样本的分类，注意使用因子型（factor）标记不同组
> targets
        Treatment  Sample Patient
human1    primary  human1       1
human2  recurrent  human2       1
human3    primary  human3       2
human4  recurrent  human4       2
human5    primary  human5       3
human6  recurrent  human6       3
human7  recurrent  human7       4
human8  recurrent  human8       5
human9    primary  human9       6
human10   primary human10       7

# 创建DESeqDataSet对象，使用Treatment作为比较
glioPR <- DESeqDataSetFromMatrix(countData = htCountSelect[, 1:10], colData = targets, design = ~ Treatment)
mcols(glioPR) <- htCountSelect[, 11:14]

# 创建DESeqDataSet对象，使用Treatment和Pateint两个因子作为比较
glioPR <- DESeqDataSetFromMatrix(countData = htCountSelect[, 1:10], colData = targets, design = ~ Patient + Treatment)
mcols(glioPR) <- htCountSelect[, 11:14]

# DEGs
glioPR <- DESeq(glioPR)
res <- results(glioPR)
## 根据FDR值排序
res[order(res$padj), ]
summary(res)

# 两种标准化count值方法， rlog（regularized log）和VST（variance stabilizing transformations）
rld <- rlog(glioPR)
vsd <- varianceStabilizingTransformation(glioPR)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
{% endcodeblock %}





























### 参考网址 ###

* Biconductor DESeq2 workflows: [RNA-Seq workflow at the gene level](http://www.bioconductor.org/help/workflows/rnaseqGene/) 



### 更新记录 ###

2015年6月5日









