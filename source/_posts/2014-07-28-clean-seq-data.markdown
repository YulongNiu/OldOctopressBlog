---
layout: post
title: "清洗二代测序数据"
date: 2014-07-28 20:11:18 +0800
comments: true
categories: bioinfor
---

## 1. FastQC ##

**简介**：[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)是用于对二代测序数据质量快速检验的工具，可以输入fastq（fastq.gz）、sam或者bam文件。查看[输出结果解释](http://www.bioinformatics.nl/courses/RNAseq/FastQC_Manual.pdf)。

<!--more-->

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

## 2. Trim Galore! ##

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


## 3. Trimmomatic ##

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


### 更新记录 ###

2016年9月10日
