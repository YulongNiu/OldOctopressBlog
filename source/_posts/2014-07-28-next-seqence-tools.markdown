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

|MAPQ (tophat)   |Tag   |描述   |
|---+---+---|
|50   |NH:i:1   |   |
|30   |NH:i:2   |   |
|1   |NH:i:3 + NH:i:4   |   |
|0   |$\sum_{n>4}{NH:i:n}$   |   |


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

# bam格式文件转换未为sam格式文件
$ samtools view -h bamFile > samFile

# bam格式文件转换为sorted bam格式文件，用于长期储存和后续分析
# 后一个‘sortedBamFilePrefix’是指需要存储文件名前缀，比如想存储“human_1.bam”，则输入“human_1”
$ samtools sort bamFile sortedBamFilePrefix

# 直接查看bam文件
$ samtools view bamFile | head 2

# 输出alignment数目，配合-f和-F过滤reads
$ samtools view -c bamFile

# 统计bam文件map数量，可以用于评估mapping的质量。需要输入indexed的bam文件
# 输出结果每一列分别为参考序列名称、参考序列长度、map上的reads数目、未map上的reads数目
$ samtools idxstats sortedIndexedBamFile

{% endcodeblock %}

**重要参数解释**：

* `view -h`{:.language-bash}：打印sam文件头部，文件头部信息用于转换sam文件至bam文件。



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

**简介**：[TopHat](http://ccb.jhu.edu/software/tophat/index.shtml)是快速将RNA测序片段“对应（mapping）”到基因组上的工具，优势在于处理外显子间的剪切。内部首先使用bowtie或bowtie2把RNA测序片段“比对（alignment）”到基因组，之后再分析和鉴定剪切连接区域。

**平台**：Mac OS/Linux

**快速运行**

{% codeblock lang:bash %}
# 双端测序
$ tophat2 -p 8 -o human_1 hg19 human_1_1.fastq.gz human_1_2.fastq.gz

# 单端测序
$ tophat2 -p 8 -o human_1 hg19 human_1.fastq.gz
{% endcodeblock %}


**重要参数解释**：

* `-p`{:.language-bash}：设置线程数，用于多核计算。

* `-o`{:.language-bash}：执行输入文件夹。

* `-r/--mate-inner-dist`{:.language-bash}：mate-paired reads的间隔长度的期望平均值，默认值为50bp。具体解释，参考 [RNA-seq差异表达分析工作流程](http://blog.qiuworld.com:8080/archives/3007)和[Tophat中-r/–mate-inner-dist参数](http://www.plob.org/2012/12/04/4988.html)。

* `--mate-std-dev`{:.language-bash}：mate-paired reads的间隔长度分布的标准差，默认值为20bp。

* `--library-type`{:.language-bash}：测序仪器和方法，默认为标准Illumina平台的`fr-unstranded`{:.language-bash}。其他平台设置，详见[TopHat说明文档](http://ccb.jhu.edu/software/tophat/manual.shtml)和 [链特异性转录组原理](http://www.plob.org/2013/12/03/6731.html) 
。 

* `--no-discordant`{:.language-bash}：只对于paired reads，只报告concordant mappings。设置了也没用？？？？？

* `--no-mixed`{:.language-bash}：只对于paired reads，只报告paired reads都成功map。

**后续操作**：

* TopHat2运行后查看`align_summary.txt`{:.language-bash}获得比对结果。

* TopHat2会输出`accepted_hits.bam`{:.language-bash}（接受map的reads文件）和`unmapped.bam`{:.language-bash}（没有map上的reads文件）。对于后者，使用基因组浏览器，如[IGV](http://www.broadinstitute.org/igv/)或者[UCSC Genome Browser](http://genome.ucsc.edu/)大致看下是有无map，之后可以直接丢弃。

* 过滤双端测序的TopHat结果，参考：？？？？？？？？？？？？？？？？


**补充**：

* TopHat2运行时，可以将基因组的fasta格式文件（比如`hg.fa`{:.language-bash}），一起放在index的文件夹中（如果没有，TopHat2先生成一个），可以节省运行时间。具体生成方法，参考:????????????????????

* [不加思考地使用默认参数的下场](http://www.acgt.me/blog/2015/4/27/the-dangers-of-default-parameters-in-bioinformatics-lessons-from-bowtie-and-tophat?utm_content=bufferb2c35&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer)





## 4. Cufflinks ##


## 5. Trinity ##

[Trinity](http://trinityrnaseq.github.io/)是[Broad Institute](http://www.broadinstitute.org/)开发的根据RNA-seq数据从头（*de novo*）组装转录组的工具。


## 6. miRDeep2 ##

[miRDeep2](https://www.mdc-berlin.de/36105849/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/documentation#mapper) 是用于高通量测序中检测[miRNA](http://en.wikipedia.org/wiki/MicroRNA) 的工具。

## 7. HTseq ##

[HTseq](http://www-huber.embl.de/users/anders/HTSeq/doc/index.html) 是用于Python平台写成的处理高通量测序的平台。






## 7. 小工具集锦 ##

### 7.1 查看fastq文件 ###

查看fastq文件，比如`mySeq_1_1.fastq`测序文件，在Linux机器上可以使用`head`命令输出文件前几行，便于快速查看测序文件。

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
{% endcodeblock %}

### 7.2 FastQC ###

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)是用于对二代测序数据质量快速检验的工具，可以输入fastq、sam或者bam文件。

**平台**：所有平台。

**安装**： 依赖Java，下载后直接安装使用。

**快速运行**：

{% codeblock lang:bash %}
# 支持批量处理测序数据
$ fastqc seqFile1 seqFile2 seqFileN

# 查看帮助信息
$ fastqc --help
{% endcodeblock %}

### 7.3 RNA-SeQC ###


### 7.4 Picard ###

[Picard](http://broadinstitute.github.io/picard/)是一个Java平台的工具包，包括一系列处理高通量测序的命令行工具。

### 7.5 FASTX-Toolkit ###

[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)是一系列用于处理大量短序列FASTA和FASTQ文件的工具。
















### 参考网址 ###

* [一篇介绍cufflinks的中文博客](http://www.chenlianfu.com/?p=623)

* [cufflinks安装介绍](http://cole-trapnell-lab.github.io/cufflinks/getting_started/#common-uses-of-the-cufflinks-package)

* [NGS分析入门：操作SAM/BAM文件](http://blog.qiuworld.com:8080/archives/3050)

* [Essential AWK Commands for Next Generation Sequence Analysis](http://bioinformatics.cvr.ac.uk/blog/essential-awk-commands-for-next-generation-sequence-analysis/) 



### 更新记录 ###

2015年5月15日
