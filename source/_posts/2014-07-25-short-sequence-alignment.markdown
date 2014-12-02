---
layout: post
title: "二代测序中的短序列比对"
date: 2014-07-25 14:04:25 -0400
comments: true
categories: Bioinfor
---

在二代测序数据分析中，非常重要的一步是将测得的短序列“对应”到基因组上。所使用的工具被称为“短序列比对工具（short sequence aligners）”。以下是一些常用工具的介绍。

## 1. Bowtie ##

**简介**：[Bowtie2](http://bowtie-bio.sourceforge.net/)是现在广泛使用的序列比对工具。

**运行方式**：所有平台

**特点**：

* 相比较[Bowtie1](http://bowtie-bio.sourceforge.net/index.shtml)，处理大于50bp的短序列，速度更快、也更敏感。Bowtie1在处理小于50bp的短序列上，效果更好。

* [iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.ilmn)提供一些事先编排（index）的基因组。
<!--more-->

**快速运行**：

{% codeblock lang:bash %}
# 建立一系列FASTA文件目录
$ bowtie2-build /filePath/fastaFile indexName

# unpaired序列比对
$ bowtie2 -p 4 -x indexName -U readFiles -S samFileName

# paired序列比对
$ bowtie2 -p 4 -x indexName  -1 readFiles1 -2 readFiles2 -S eg2.sam
{% endcodeblock %}

* `-p`：多线程

* `-x`：之后跟index名称

* `-U`：测序文件（比如Fasta，Fastq文件）

* `-1`/`-2`：标识paired文件

* `-S`：SAM格式输出文件




### 参考网址 ###

* [Heng Li总结的工具列表](http://lh3lh3.users.sourceforge.net/NGSalign.shtml)

