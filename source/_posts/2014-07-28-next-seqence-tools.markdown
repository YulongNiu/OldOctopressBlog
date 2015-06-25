---
layout: post
title: "二代测序中的常用工具介绍"
date: 2014-07-28 12:23:42 -0400
comments: true
published: false
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


## 3. TopHat ##


## 4. Cufflinks ##


## 5. Trinity ##

[Trinity](http://trinityrnaseq.github.io/)是[Broad Institute](http://www.broadinstitute.org/)开发的根据RNA-seq数据从头（*de novo*）组装转录组的工具。


## 6. miRDeep2 ##

[miRDeep2](https://www.mdc-berlin.de/36105849/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/documentation#mapper) 是用于高通量测序中检测[miRNA](http://en.wikipedia.org/wiki/MicroRNA) 的工具。

## 7. HTseq ##


## 8. MISO ##

**简介**：[MISO](http://genes.mit.edu/burgelab/miso/)（Mixture-of-Isoforms）是用来计算和探测RNA-seq数据中不同样本的基因选择剪切。

**平台**：Python跨平台使用。

**快速运行**：

{% codeblock lang:bash An Example of Running MISO%}

{% endcodeblock %}



## 8. 质量检测 ##

### 8.1 FastQC ###

**简介**：[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)是用于对二代测序数据质量快速检验的工具，可以输入fastq（fastq.gz）、sam或者bam文件。查看[输出结果解释](http://www.bioinformatics.nl/courses/RNAseq/FastQC_Manual.pdf)。

**平台**：所有平台。

**安装**： 依赖Java，下载后直接安装使用。

**快速运行**：

{% codeblock lang:bash %}
# 输出分析结果至特定文档
$ fastqc seqFile1 --outdir setFolder1

# 支持批量处理测序数据
$ fastqc seqFile1 seqFile2 seqFileN

# 查看帮助信息
$ fastqc --help

# 查看一共分析了多少个reads，比如fastqc文件为“accepted_filtered_fastqc.zip”
$ unzip -p accepted_filtered_fastqc.zip accepted_filtered_fastqc/fastqc_data.txt | \
      sed -n '7 p' | \
      awk '{print $3}'
{% endcodeblock %}

### 8.2 RNA-SeQC ######

[RNA-SeQC](http://www.broadinstitute.org/cancer/cga/rna-seqc)是用于检测RNA-seq测序质量。


### 8.3 RSeQC ###

[RSeQC](http://rseqc.sourceforge.net/)用于对RNA-seq数据质量控制。





## 9. 去除adapter ##

### 9.1 Trim Galore! ###

**简介**：[Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)是对[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)和[Cutadapt](https://cutadapt.readthedocs.org/en/stable/installation.html)的包装。可以处理Illumina、Nextera 3和smallRNA测序平台的双端和单端数据，包括去除adapter和低质量reads。

**平台**：Linux

**安装**:

* 需要先分别安装FastQC和Cutadapt，其中Cutadapt安装使用

{% codeblock lang:bash Install Cutadapt%}
# pip install cutadapt
{% endcodeblock %}

**快速运行**：

{% codeblock lang:bash Example code of runing the Cutadapt%}
# 处理双端测序结果
$ trim_galore --quality 20 --phred33 --stringency 3 --length 20 --paired \
              --gzip --output_dir human_trimgalore \
              mySeq_1_1.fastq.gz mySeq_1_2.fastq.gz

{% endcodeblock %}

**重要参数解释**：

* `--quality`{:.language-bash$}：设定Phred quality score阈值，默认为20。

* `--phred33`{:.language-bash$}：：选择`-phred33`{:.language-bash}或者`-phred64`{:.language-bash}，表示测序平台使用的Phred quality score。

* `--adapter`{:.language-bash$}：输入adapter序列。也可以不输入，Trim Galore!会自动寻找可能性最高的平台对应的adapter。自动搜选的平台三个，也直接显式输入这三种平台，即`--illumina`{:.language-bash}、`--nextera`{:.language-bash}和`--small_rna`{:.language-bash}。

* `--stringency`{:.language-bash$}：设定可以忍受的前后adapter重叠的碱基数，默认为1（非常苛刻）。可以适度放宽，因为后一个adapter几乎不可能被测序仪读到。

* `--length`{:.language-bash$}：设定输出reads长度阈值，小于设定值会被抛弃。

* `--paired`{:.language-bash$}：对于双端测序结果，一对reads中，如果有一个被剔除，那么另一个会被同样抛弃，而不管是否达到标准。

* `--retain_unpaired`{:.language-bash$}：对于双端测序结果，一对reads中，如果一个read达到标准，但是对应的另一个要被抛弃，达到标准的read会被单独保存为一个文件。

* `--gzip`{:.language-bash$}和`--dont_gzip`{:.language-bash$}：清洗后的数据zip打包或者不打包。

* `--output_dir`{:.language-bash$}：输入目录。需要提前建立目录，否则运行会报错。


### 9.2 Trimmomatic ###

**简介**：[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)是针对Illumina高通量测序平台设计的接头去除和低质量reads清洗软件。软件中包括有Illumina平台常见接头序列，可以很方便处理单端和双端RNA-seq数据。Trimmomatic也支持自己设计要去除的接头序列文件。

**平台**：Java跨平台使用

**快速运行**：

{% codeblock lang:bash Example code of runing the Trimmomatic %}
# 处理双端测序结果
$ java -jar /path/trimmomatic-0.33.jar PE\
       -threads 12 -phred33 -trimlog mySeq_1-trim.log \
       mySeq_1_1.fastq.gz mySeq_1_2.fastq.gz \
       mySeq_1_1-trim.fastq.gz mySeq_1_1-unpair.fastq.gz \
       mySeq_1_2-trim.fastq.gz mySeq_1_2-unpair.fastq.gz \
       ILLUMINACLIP:/path/TruSeq3-PE.fa:2:30:10 \
       LEADING:3 \
       TRAILING:3 \
       SLIDINGWINDOW:4:15 \
       MINLEN:51
{% endcodeblock %}

**重要参数解释**：

* `-threads`{:.language-bash}：设置线程数目。

* `-phred33`{:.language-bash}：选择`-phred33`{:.language-bash}或者`-phred64`{:.language-bash}，表示测序平台使用的Phred quality score。查询方法：首先，运行[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)，在结果报告第一项会“猜出”测序平台。之后，查询平台对应Phred[列表](http://en.wikipedia.org/wiki/FASTQ_format)。

* `-trimlog`{:.language-bash}：输出运行日志。日志中包括对每一个read具体选择数据，所以文件会比较大。

* `ILLUMINACLIP`{:.language-bash}：跟随四个参数，分别是`:<fastaWithAdaptersEtc>`{:.language-bash}为adaptesr文件完整路径（在Trimmomatic的默认安装目录下的 `adapter`{:.language-bash}，有整理好的）；`<fastaWithAdaptersEtc>`{:.language-bash}为seed matches（16bases）在匹配时的最大错配数目；`<palindrome clip threshold>`{:.language-bash}对于一对reads当得分超过30（约50 bases），seeds会被延伸和固定；`<simple clip threshold>`{:.language-bash}，对于单端reads当得分超过10（约17 bases），seeds会被延伸和固定。

* `LEADING`{:.language-bash}和`TRAILING`{:.language-bash}：分别为去除read头部和尾部的低质量（低于quality3）碱基数目。

* `SLIDINGWINDOW`{:.language-bash}：跟随两个参数，分别是 `<windowSize>`{:.language-bash}为扫描“窗口”长度；`<requiredQuality>`{:.language-bash}为窗口碱基质量的平均阈值，低于此会被删除。

* `MINLEN`{:.language-bash}：设置最短reads数目。需要根据下游alignment软件设定，比如[Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)适用于短序列，比如50bp以下；而[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)适用于50bp以上。[TopHat](http://ccb.jhu.edu/software/tophat/index.shtml) 则根据实际使用Bowtie或者Bowtie2选择。



## 10. 小工具集锦 ##

### 10.1 操作fastq文件 ###

在Linux机器上配以一系列命令，操作fastq或者fastq.gz文件，比如`mySeq_1_1.fastq`测序文件。

{% codeblock lang:bash Manipulating Fastq file with bash %}
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

# 统计Fastq文件中的reads数量
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


### 10.2 Picard ###

[Picard](http://broadinstitute.github.io/picard/)是一个Java平台的工具包，包括一系列处理高通量测序的命令行工具。

### 10.3 FASTX-Toolkit ###

[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)是一系列用于处理大量短序列FASTA和FASTQ文件的工具。


## 11. 各种数据类型 ##

为了下游分析，高通量测序结果往往使用多种数据格式储存，详细参考[File Formats](http://www.broadinstitute.org/igv/?q=book/export/html/16)。

* GTF格式文件，比如[Ensembl类型](http://mblab.wustl.edu/GTF22.html)的GTF文件，每一列[说明](http://www.gencodegenes.org/gencodeformat.html)、生物学序列分类[GENCODE解释](http://www.gencodegenes.org/gencode_biotypes.html)和[VEGA解释](http://vega.sanger.ac.uk/info/about/gene_and_transcript_types.html)。GTF文件的第二行中，有[CDS（Coding DNA Sequence）](http://en.wikipedia.org/wiki/Coding_region)也有[exon](http://en.wikipedia.org/wiki/Exon)，这两者概念不同：CDS只包括翻译成蛋白质的序列；exon也包括[UTR](http://en.wikipedia.org/wiki/Untranslated_region)区域和ployA区域；exon也可以用来指示非编码RNA。

### 11.1 GTF和GFF文件互转 ###

使用[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/)的`gffread`{:.language-bash}直接进行转换。


{% codeblock lang:bash Convert the File Format Between GTF and GFF%}
# GTF转换为GFF/GFF3
$ gffreads myGtfFile.gtf -o myGffFile.gff3

# GFF/GFF3转换为GTF
$ gffread myGffFile.gff3 -T -o myGtfFile.gtf
{% endcodeblock %}

### 11.2 GTF文件提取序列信息 ###

使用[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/file_formats/)的`gffread`{:.language-bash}直接进行转换。

{% codeblock lang:bash Extract sequences from GTF files%}
# 从GTF文件中提取序列信息
$ gffread -w myFastaFile.fa -g myGenome.fa myGtfFile.gtf
{% endcodeblock %}

### 11.3 GTF文件提取序列长度信息 ###

exon












### 参考资料 ###

* [NGS分析入门：操作SAM/BAM文件](http://blog.qiuworld.com:8080/archives/3050)

* [Essential AWK Commands for Next Generation Sequence Analysis](http://bioinformatics.cvr.ac.uk/blog/essential-awk-commands-for-next-generation-sequence-analysis/) 




### 更新记录 ###

2015年5月23日
