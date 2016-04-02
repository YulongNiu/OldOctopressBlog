---
layout: post
published: true
title: "RNA-seq基本概念和数据分析流程"
date: 2014-07-22 17:03:15 -0400
comments: true
styles: [data-table]
published: true
categories: bioinfor
---

RNA-seq相比较基因芯片，价格虽然昂贵一些，精度和灵敏度更高。同时，在测序深度足够时，也可以检测mRNA选择性剪切类型。

## 1. 样品制备 ##



<!--more-->


## 2. 测序 ##

Illumina双端测序[动画](http://www.illumina.com/technology/next-generation-sequencing/paired-end-sequencing_assay.html)和[图示](http://onetipperday.blogspot.sg/2013/12/illumina-hiseq2000-adapter-and.html)。

## 3. 分析流程 ##

RNA-seq[主要分析流程](http://blog.qiuworld.com:8080/archives/2343)：


{% codeblock lang:bash Workflow of RNA-seq data analysis %}
原始数据质量评估 -->
数据清洗（去除接头和低质量read） -->
清洗数据质量评估 -->
map测序结果至基因组（转录组） -->
map数据质量评估 -->
差异表达基因/选择性剪切/新基因/融合基因选择 -->
GO和pathway分析 -->
共表达网络分析
{% endcodeblock %}



### 3.1 序列清洗 ###

序列清洗主要是去除测序结果中的[adapter](http://onetipperday.blogspot.sg/2013/06/illumina-hiseq2000-adapter.html)，

以下使用Illumina HiSeq2000平台，对一个人类样本的RNA测序。统计各种序列清洗方法和选择后的reads数目。原始数据两端测序reads分别为54,492,228和54,492,228。

-----------------------

| Method                         | #Trimmed    | #Mapped*    | #Filtered  |
|--------------------------------+-------------+-------------+------------|
| r50-notrim                     | 108,984,456 | 109,278,388 | 79,143,942 |
| r50-nomixed-notrim             | 108,984,456 | 103,548,800 | 79,143,942 |
| r50-nomixed-trimmomatic-min20  | 104,164,622 | 116,315,394 | 80,337,256 |
| r50-nomixed-trimmomatic-min36  | 101,548,172 | 110,778,108 | 79,248,896 |
| r50-nomixed-trimmomatic-min50  | 98,525,312  | 106,659,988 | 77,779,424 |
| r50-nomixed-galore-min20       | 107,097,862 | 114,943,386 | 83,039,928 |
| r100-nomixed-galore-min20      | 107,097,862 | 114,944,672 | 87,899,316 |
| r165sd45-nomixed-galore-min20  | 107,097,862 | 114,201,366 | 90,750,208 |
| r165sd45G-nomixed-galore-min20 | 107,097,862 | 109,901,742 | 93,477,122 |
| r165sd45-nomixed-galore-min50  | 104,869,208 | 109,258,544 | 89,329,672 |

**\***：使用TopHat2把序列mapped到hs19基因组。TopHat2默认设置为，如果一个reads能mapped到多个位点，则都会报道。因此数目可能比原始数据多。

-----------------------

## 4. 分析过程注意事项 ##

1. 参考转录组注释文件（GFF/GTF）的染色体编号，与基因组信息一致。比如都用`chr1`{:.language-bash}，或者都用`1`{:.language-bash}标识1号染色体。可能出现问题地方：

> * [TopHat](http://ccb.jhu.edu/software/tophat/manual.shtml)的`-G/--GTF`{:.language-bash}参数。
>
> * [MISO](http://miso.readthedocs.org/en/fastmiso/#human-mouse-gene-models-for-isoform-centric-analyses)程序`pe_utils`{:.language-bash}的传入转录组注释文件。
>







## 5. 基本概念 ##

* **Read**：是组成测序的结果的基本单位。

* **Count**：在某次测序中，对于某个指标（比如某个基因），得到的reads数的总和。

* **RPKM**（Reads Per Kilobase per Million）和**FPKM**（Fragments Per Kilobase per Million）：

> 首先需要解释FPKM和RPKM的原理是相似的，区别在于FPKM对应的是DNA片段，比如在一个Illumina的pair-end（双尾）RNA-seq中，一对（两个）reads对应是一个DNA片段。有了FPKM（RPKM）概念，我们就能比较：同一个样本中基因A和基因B的相对表达量；或者不同样本中，同一个基因的相对表达量。
>
> 具体的原因是：引入“每一千碱基（per kilobase）”的原因在于，不同的RNA可能有不同长度，长度越长，对应的reads就越多。当每个RNA都除以自身长度（以1000碱基为单位）时，就可以比较同一个样本中不同基因的相对表达量了。相似地，引入“每一百万reads”的原因是，不同的样本可能测序的深度不一样，深度越深，当然对应的reads就越多了。如果结果除以各自库的数量（以一百万reads为单位），那么我们就能很好地衡量两个不同样本中同一个基因的相对表达量。

* **Alignment**：确定测得的reads在基因组上的位置的过程。

* **Mapping**：确定aligned reads对应的转录本。

* **Pair end (PE)**和**Mate-Pair (MP)**：

> 两种双端测序的方法，主要区别在样品库制备和测序上。比如PE制备的库是adapter在目标序列两端，而MP库中adapter在目标序列中间。因此，在数据分析时，MP类型测序必须注意剔除adapter。具体参考[论坛讨论](http://seqanswers.com/forums/showthread.php?t=503)和[Difference Between Paired-End and Mate-Pair Reads](http://scottmyourstone.blogspot.sg/2013/11/difference-between-paired-end-and-mate.html)。
>
> 在序列软件软件中，有时也称呼一对测序reads（对同一个目标片段分别测得的正义链和反义链），它们互为mate，分别存在两个对应的测序结果文件中。


* **Adapter（接头）**、**Barcode（标签）**和**Insert（插入片段）**：

> adapter是一段短的序列已知的核酸链，用于链接序列未知的目标测序片段。
> 
> barcode，也称为index，是一段很短的寡居核酸链，用于在多个样品混合测序时，标记不同的样品。
> insert是用于测序的目标片段，因为是包括在两个adapter之间，所以被称为“插入”片段。
>
> 一个常见测序片段类似与`adapter--barcode--insert--adapter`{:.language-bash}。测序开始时前几个碱基无法测得，第一个adapter在数据输出时被去除；由于测序仪读长限制，第二个adapter通常无法测得。所以，经常得到类似 `barcode--部分insert`{:.language-bash}的read。最后，把barcode去除，只保留测度insert的片段，这个操作的术语是demultiplexing。
>
> 需要注意的是，[Illumina TruSeq](http://www.illumina.com/documents/products/datasheets/datasheet_truseq_sample_prep_kits.pdf)样品库制备方法中，barcode是在adapter中部，而且是与insert分开测序。而[Illumina Nextera Mate Pair](http://res.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf)样品库制备中，adapter在目标序列中部。

* **Concordant Pairs**和**Discordant Pairs**：根据[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-inputs)的解释，concordant pairs表示一对reads在alignment时，既方向匹配又有合适距离（Bowtie2中是200bp～500bp）。如果上述方向和距离，任意一个条件不满足，则称为discordant pairs。

* **chrN_random**和**chrUn**：基因组文件中通常含有类似`chr9_gl000198_random`{:.language-bash}和`chrUn_gl000211`{:.language-bash}的基因组。根据[UCSC](https://genome.ucsc.edu/FAQ/FAQdownloads.html#download11)解释，`chrN_random`{:.language-bash}包括基因组已知但具体位置未知序列，或者位置已知但具体内容未完成序列。`ChrUn`{:.language-bash}中包括一些具体位置未知的序列。

* **chrN_xxx_hap1**：转录组注释文件中会出现类似`chr6_apd_hap1`{:.language-bash}和`chr6_dbb_hap3`{:.language-bash}的基因组注释。根据[UCSC解释](http://genome.ucsc.edu/cgi-bin/hgGateway?org=Human&db=hg19)这些基因组是[单倍型(haplotype)]基因组(http://hapmap.ncbi.nlm.nih.gov/originhaplotype.html.en)。



### 参考网址 ###

* 推荐幻灯片[RNA-seq quality control and pre-processing](http://www.slideshare.net/mikaelhuss/all-bio-rnaseqqc?from_action=save)

* [Bioconductor详细流程](http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_12_16_2013/Rrnaseq/Rrna)

* RPKM和FPKM：[1](http://jefworks.com/rpkm-and-fpkm-explained/)、 [2](http://www.cureffi.org/2013/09/12/counts-vs-fpkms-in-rna-seq/)和[3一个计算RPKM的例子](http://www.partek.com/Tutorials/microarray/User_Guides/UnderstandingReads.pdf)

* ENCODE推荐的RNA-seq数据分析指导 [The ENCODE Consortium: Standards, Guidelines and Best Practices for RNA-seq](http://genome.ucsc.edu/ENCODE/protocols/dataStandards/ENCODE_RNAseq_Standards_V1.0.pdf)

* RNA-seq差异表达分析工作流程：[中文](http://216.49.144.90:8080/archives/3007/comment-page-2)，[英文](http://vallandingham.me/RNA_seq_differential_expression.html)

* [多少个read才够](http://www.rna-seqblog.com/how-many-reads-are-enough/)

* [高通量测序常用名词汇总](http://www.macrogencn.com/_d275872179.htm)

* [Bowtie2对一些常用名词解释](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-inputs) 


* [二代测序中barcodes index的介绍](http://www.plob.org/2014/11/09/8672.html)

* Illumina样品制备参考：[Illumina TruSeq DNA Adapters De-Mystified](https://www.med.unc.edu/pharm/calabreselab/files/tufts-sequencing-primer)和[Illumina adapter and primer sequences](http://bioinformatics.cvr.ac.uk/blog/illumina-adapter-and-primer-sequences/)

* [Quality Control and Trimming of RNA-seq data](http://www.research.janahang.com/quality-control-and-trimming-of-rna-seq-data/)

* [Best Practices and Suggestions for Read Alignment](http://cgrlucb.wikispaces.com/file/view/readMapping.pdf)

* [Using RNA-seq to quantify gene levels and assay for differential expression](http://barcwiki.wi.mit.edu/wiki/SOPs/rna-seq-diff-expressions) 




### 更新记录 ###

2015年5月19日
