---
layout: post
published: true
title: "RNA-seq数据分析的基本概念和流程"
date: 2014-07-22 17:03:15 -0400
comments: true
published: false
categories: bioinfor
---

RNA-Seq相比较基因芯片，价格虽然昂贵一些，但是效果可能会更好一些（在同样样本和重复基础上）。

## 1. 测序原理 ##



## 2. 基本概念 ##

* **Read**： 

* **Count**：在某次测序中，对于某个指标（比如某个基因），得到的reads数的总和。

* **RPKM**（Reads Per Kilobase per Million）和**FPKM**（Fragments Per Kilobase per Million）：

> 首先需要解释FPKM和RPKM的原理是相似的，区别在于FPKM对应的是DNA片段，比如在一个Illumina的pair-end（双尾）RNA-seq中，一对（两个）reads对应是一个DNA片段。有了FPKM（RPKM）概念，我们就能比较：同一个样本中基因A和基因B的相对表达量；或者不同样本中，同一个基因的相对表达量。
>
> 具体的原因是：引入“每一千碱基（per kilobase）”的原因在于，不同的RNA可能有不同长度，长度越长，对应的reads就越多。当每个RNA都除以自身长度（以1000碱基为单位）时，就可以比较同一个样本中不同基因的相对表达量了。相似地，引入“每一百万reads”的原因是，不同的样本可能测序的深度不一样，深度越深，当然对应的reads就越多了。如果结果除以各自库的数量（以一百万reads为单位），那么我们就能很好地衡量两个不同样本中同一个基因的相对表达量。

* **Alignment**：确定测得的reads在基因组上的位置的过程。

* **Mapping**：确定aligned reads对应的转录本。

## 3. 分析流程 ##

<!--more-->


### 参考网址 ###

* [Bioconductor详细流程](http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_12_16_2013/Rrnaseq/Rrna)

* RPKM和FPKM：[1](http://jefworks.com/rpkm-and-fpkm-explained/), [2](http://www.cureffi.org/2013/09/12/counts-vs-fpkms-in-rna-seq/)，[3 一个计算RPKM的例子](http://www.partek.com/Tutorials/microarray/User_Guides/UnderstandingReads.pdf)

* ENCODE推荐的RNA-seq数据分析指导 [The ENCODE Consortium: Standards, Guidelines and Best Practices for RNA-Seq](http://genome.ucsc.edu/ENCODE/protocols/dataStandards/ENCODE_RNAseq_Standards_V1.0.pdf)

* RNA-seq差异表达分析工作流程：[中文](http://216.49.144.90:8080/archives/3007/comment-page-2)，[英文](http://vallandingham.me/RNA_seq_differential_expression.html)

* [多少个read才够](http://www.rna-seqblog.com/how-many-reads-are-enough/)

* [高通量测序常用名词汇总](http://www.macrogencn.com/_d275872179.htm) 


### 更新记录 ###

2015年3月12日
