---
layout: post
title: "二代测序中的常用工具介绍"
date: 2014-07-28 12:23:42 -0400
comments: true
categories: Bioinfor
---

## 1. SAMtools ##

[SAMtools](http://samtools.sourceforge.net/)是用于处理SAM (Sequence Alignment/Map)格式文件的一系列工具，主要用来储存大容量的核酸测序结果。

**简介**：BAM是SAM文件的binary格式文件。SAMtools的主要作者是[Heng Li](http://lh3lh3.users.sourceforge.net/)，Heng Li在2012年因为对二代测序领域的贡献获得[Benjamin Franklin Award](http://en.wikipedia.org/wiki/Benjamin_Franklin_Award_(Bioinformatics))。

**平台**：Mac OS/Linux

**安装**: 

* [下载](http://sourceforge.net/projects/samtools/files/)SAMtools

<!--more-->

* 安装依赖的Ncurse库

{% codeblock lang:bash %}
# yum install ncurses*
{% endcodeblock %}

* 添加SAMtools到PATH（[设置参考](http://yulongniu.bionutshell.org/blog/2010/11/08/linux-command/)）

**快速运行**

{% codeblock lang:bash %}
# 将sam格式文件转换为bam格式文件
$ samtools view -bS samFile > bamFile

# 将bam格式文件转换为sorted bam格式文件，用于长期储存和后续分析 
$ samtools sort bamFile sortedBamFile

# 统计bam文件map数量，可以用于评估mapping的质量。需要输入sorted和indexed的bam文件
$ samtools idxstats sortedIndexedBamFile
{% endcodeblock %}


## 2.  ##



## 3. Tophat ##




## 4. Cufflinks ##




## 5. 小工具集锦 ##

### 5.1 查看fastq文件 ###

查看fastq文件，比如`mySeq_1_1.fastq`测序文件，在Linux机器上可以使用`head`命令输出文件前几行，便于快速查看测序文件。

{% codeblock lang:bash %}
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
{% endcodeblock %}

### 5.2 FastQC ###

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



### 参考网址 ###

* [一篇介绍cufflinks的中文博客](http://www.chenlianfu.com/?p=623)

* [cufflinks安装介绍](http://cole-trapnell-lab.github.io/cufflinks/getting_started/#common-uses-of-the-cufflinks-package)



### 更新记录 ###

2015年3月13日
