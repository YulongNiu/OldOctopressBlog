---
layout: post
title: "二代测序中的常用工具介绍"
date: 2014-07-28 12:23:42 -0400
comments: true
published: true
categories: Bioinfor
---

## 1. SAMtools ##

**简介**：[SAMtools](http://www.htslib.org/)是用于处理SAM (Sequence Alignment/Map)格式文件的一系列工具，主要用来储存大容量的核酸测序结果。BAM是SAM文件的binary格式文件。SAMtools的主要作者是[Heng Li](http://lh3lh3.users.sourceforge.net/)，Heng Li在2012年因为对二代测序领域的贡献获得[Benjamin Franklin Award](http://en.wikipedia.org/wiki/Benjamin_Franklin_Award_(Bioinformatics))。

**平台**：Mac OS/Linux

**安装**: 

* [下载](http://www.htslib.org/download/)SAMtools

<!--more-->

* 安装依赖的Ncurse库

{% codeblock lang:bash %}
# yum install ncurses*
{% endcodeblock %}

* 添加SAMtools到PATH（[设置参考](http://yulongniu.bionutshell.org/blog/2010/11/08/linux-command/)）

**快速运行**

{% codeblock lang:bash %}
# sam格式文件转换为bam格式文件
# 新版本samtools不用使用-S
$ samtools view -b samFile > bamFile

# bam格式文件转换为sam格式文件
$ samtools view -h bamFile > samFile

# bam格式文件转换为sorted bam格式文件，用于长期储存和后续分析
# 后一个‘sortedBamFilePrefix’是指需要存储文件名前缀，比如想存储“human_1.bam”，则输入“human_1”
$ samtools sort bamFile sortedBamFilePrefix

# bam文件按照reads名称排序
$ samtools sort -n bamFile sortedBamFilePrefix

# 直接查看bam文件
$ samtools view bamFile | head -2

# 创建bam的index文件
$ samtools index bamFile

# 输出alignment数目，配合-f和-F过滤reads
$ samtools view -c bamFile

# 统计bam文件map数量，可以用于评估mapping的质量。需要输入indexed的bam文件
# 输出结果每一列分别为参考序列名称、参考序列长度、map上的reads数目、未map上的reads数目
$ samtools idxstats sortedIndexedBamFile

{% endcodeblock %}

**重要参数解释**：

* `view -h`{:.language-bash}：打印bam文件头部，文件头部信息用于转换sam文件至bam文件。`view -H`{:.language-bash}：只打印bam文件头部。

**补充**：

* 一个sam文件式例子

{% codeblock lang:bash A small sample of sam file%}
$ head samFile
HD      VN:1.0          SO:coordinate
@SQ     SN:chr1         LN:249250621
@SQ     SN:chr10        LN:135534747
@SQ     SN:chr11        LN:135006516

HISEQ2000-02:436:C2PG3ACXX:3:1316:4503:75198    99      chr1    11527   3       100M    =       11624   197     CTGGGTTTAAAAGTAAAAAATAAATATGTTTAATTTGTGAACTGATTACCATCAGAATTGTACTGTTCTGTATCCCACCAGCAATGTCTAGGAATGCCTG  ?@@FFFD?D>DDHFGGEGGIJDHJCEGCHE@FCHIJJJGG@GEGHGJGHHIIIJJJIJJECGG@FHIIJJJIIGEADHHEEDF;BB@CACCCEDDDDDAC    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100     YT:Z:UU  NH:i:2  CC:Z:chr15      CP:i:102519544  HI:i:0

HISEQ2000-06:325:C2RC0ACXX:5:1101:16521:147635  355     chr1    11554   0       100M    =       11681   227     GTTTAATTTGTTAACTGATTACCATCAGAATTGTACTGTTCTGTATCCCACCAGCAATGTCTAGGAATGCCTGTTTCTCCACAAAGTGTTTACTTTTGGA  @@@DBDDDHFDDDHGGBA@AEFBBCGIC?BCGIAHBEFFHIHIDHC@DEGIIIHEG;8BDBFGGICFI:FG7@FHEHIICIEA;CA;77;7;;@;@A@6(    AS:i:-5 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:11G88   YT:Z:UU  NH:i:5  CC:Z:chr12      CP:i:93956      HI:i:0

HISEQ2000-06:325:C2RC0ACXX:5:2304:4393:52082    163     chr1    11579   3       100M    =       11643   164     CAGAATTGTACTGTTCTGTATCCCACCAGCAATGTCTAGGAATACCTGTTTCTCCACAAAGTGTTTACTTTTGGATTTTTGCCAGTCTAACAGGTGAAGC  CCCFFFFFHHHHHIJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJJJJJGIJJJIJJJJJFHFHHHIJJJJJIEHIJHHHHFFFFFFEEEEEDDDDDDD    AS:i:-6 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:43G56   YT:Z:UU  NH:i:2  CC:Z:chr15      CP:i:102519492  HI:i:0
{% endcodeblock %}


* sam文件内容各行解释，参考[SAM/BAM and related specifications](http://samtools.github.io/hts-specs/)和[NGS分析入门：操作SAM/BAM文件](http://blog.qiuworld.com:8080/archives/3050)。

* bam文件过滤，参考[如何统计BAM文件中的reads数](http://blog.qiuworld.com:8080/archives/3419)

## 2. Bowtie2 ##

Bowtie使用介绍，详见[二代测序中的短序列比对](http://yulongniu.bionutshell.org/blog/2014/07/26/short-sequence-alignment/) 。


## 3. Trinity ##

[Trinity](http://trinityrnaseq.github.io/)是[Broad Institute](http://www.broadinstitute.org/)开发的根据RNA-seq数据从头（*de novo*）组装转录组的工具。


## 4. miRDeep2 ##

[miRDeep2](https://www.mdc-berlin.de/36105849/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/documentation#mapper) 是用于高通量测序中检测[miRNA](http://en.wikipedia.org/wiki/MicroRNA) 的工具。

## 5. MISO ##

**简介**：[MISO](http://genes.mit.edu/burgelab/miso/)（Mixture-of-Isoforms）是用来计算和探测RNA-seq数据中不同样本的基因选择剪切。

**平台**：Python跨平台使用。

**快速运行**：

{% codeblock lang:bash An Example of Running MISO%}

{% endcodeblock %}

## 6. 质量检测 ##

### 6.1 RNA-SeQC ######

[RNA-SeQC](http://www.broadinstitute.org/cancer/cga/rna-seqc)是用于检测RNA-seq测序质量。

### 6.2 RSeQC ###

[RSeQC](http://rseqc.sourceforge.net/)用于对RNA-seq数据质量控制。

## 7. 小工具集锦 ##

### 7.1 操作fastq文件 ###

在Linux机器上配以一系列命令，操作fastq或者fastq.gz文件，比如`mySeq_1_1.fastq`测序文件。

{% codeblock lang:bash Manipulating fastq file with bash %}
# Fastq文件前几行
$ head mySeq_1_1.fastq
@HISEQ2000-06:325:C2RC0ACXX:5:1101:1217:2215 1:N:0:TGCCGGCT
GGCACAGTCCATGCTTTTAACCAGATTTGAACAGAAGAATGGCCACTTGNNNNNNGTAGNNNNNNANNNNGNNNTNNNTNNCATGTGTCACATAACTACC
+
@C@FFFFDFFHHHGFFIIBEGHJIFHIIGCGIJJIJIJJIJIJJIIIGI######-.<B######-####,###,###(##,,5<?CDDDDCCDDDCDDD
@HISEQ2000-06:325:C2RC0ACXX:5:1101:1621:2040 1:N:0:TGCCGGCT
NGGTGACCTTCTCCCGCCAGAAGCCAAAAACTGCAGCCTACTTTTCTGAAGTGGTTATCTTGGGACTGAGGTATGGGCTATCTTGGGCTGTTTCCTATTT
+
#1==D;ADBHDFDHFHGGIIB;FFHIIGGGIGADEGEIIAHEEHIGHH==C8BGC;=FCDHICHCAAC>CE..;A>AC?9@CCCCCBBCB?CCCCCCCCC
@HISEQ2000-06:325:C2RC0ACXX:5:1101:1663:2146 1:N:0:TGCCGGCT
AGAGCCAAATATTTCAACAAAACTGCAGTTTAATTTCAGAAAATGTTAAAATATATATTTATACATCAATTTCTGACATACACTTAATGTGTTAGTATAC

# 统计fastq文件中的reads数量
$ awk 'END{print NR/4}' mySeq_1_1.fastq

# 查看fastq.gz文件
$ zcat mySeq_1_1.fastq.gz

# 统计fastq.gz文件reads数量
# 以下两种方法等价
$ zcat mySeq_1_1.fastq.gz | awk 'END{print NR/4}'
$ zcat mySeq_1_1.fastq.gz | echo $((`wc -l`/4))
# 同时，由于一些测序文件会对read做特殊标记，比如“@HISEQ2000”。所以，可以使用正则匹配。
# 注意：这种方法有潜在错误，慎重使用。
$ zgrep -c '@HISEQ200' mySeq_1_1.fastq

# 统计fastq文件reads长度分布
$ zcat mySeq_1_1.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c
{% endcodeblock %}

### 7.2 操作fasta文件

{% codeblock lang:bash Manipulating fasta file with bash %}
# 序列数目
$ grep -c '>@' mySeq.fa

# 提取特定序列
samtools faidx mySeq.fa
samtools faidx mySeq.fa chr1 chr2 chr3
{% endcodeblock %}

### 7.3 Picard ###

[Picard](http://broadinstitute.github.io/picard/)是一个Java平台的工具包，包括一系列处理高通量测序的命令行工具。

### 7.4 FASTX-Toolkit ###

[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)是一系列用于处理大量短序列FASTA和FASTQ文件的工具。

## 8. 各种数据类型 ##

为了下游分析，高通量测序结果往往使用多种数据格式储存，详细参考[File Formats](http://www.broadinstitute.org/igv/?q=book/export/html/16)。

* GTF格式文件，比如[Ensembl类型](http://mblab.wustl.edu/GTF22.html)的GTF文件，每一列[说明](http://www.gencodegenes.org/gencodeformat.html)、生物学序列分类[GENCODE解释](http://www.gencodegenes.org/gencode_biotypes.html)和[VEGA解释](http://vega.sanger.ac.uk/info/about/gene_and_transcript_types.html)。GTF文件的第二行中，有[CDS（Coding DNA Sequence）](http://en.wikipedia.org/wiki/Coding_region)也有[exon](http://en.wikipedia.org/wiki/Exon)，这两者概念不同：CDS只包括翻译成蛋白质的序列；exon也包括[UTR](http://en.wikipedia.org/wiki/Untranslated_region)区域和ployA区域；exon也可以用来指示非编码RNA。

### 8.1 GTF和GFF文件互转 ###

使用[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/)的`gffread`{:.language-bash}直接进行转换。


{% codeblock lang:bash Convert the File Format Between GTF and GFF%}
# GTF转换为GFF/GFF3
$ gffreads myGtfFile.gtf -o myGffFile.gff3

# GFF/GFF3转换为GTF
$ gffread myGffFile.gff3 -T -o myGtfFile.gtf
{% endcodeblock %}

### 8.2 GTF文件提取序列信息 ###

使用[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/)的`gffread`{:.language-bash}直接进行转换。

{% codeblock lang:bash Extract sequences from GTF files%}
# 从GTF文件中提取序列信息
$ gffread -w myFastaFile.fa -g myGenome.fa myGtfFile.gtf
{% endcodeblock %}


### 参考资料 ###

* [NGS分析入门：操作SAM/BAM文件](http://blog.qiuworld.com:8080/archives/3050)

* [Essential AWK Commands for Next Generation Sequence Analysis](http://bioinformatics.cvr.ac.uk/blog/essential-awk-commands-for-next-generation-sequence-analysis/) 


### 更新记录 ###

2018年4月3日
