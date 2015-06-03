---
layout: post
title: "UCSC Table下载注释文件"
date: 2015-06-02 16:49:41 +0800
comments: true
categories: Bioinfor
---

在进行RNA-seq数据分析时，需要从[UCSC Table](http://genome.ucsc.edu/cgi-bin/hgTables)下载基因组注释信息。

## 1. 基因注释信息 ##

下载转录组Ensembl注释文件：

<img src="/images/hg19_ensembl.png" title="image" alt="UCSC下载Ensembl注释">

<!--more-->

在“table”选择中，`ensemblSource`{:.language-bash}表示Ensembl类型注释，`ensemblToGeneName`{:.language-bash}表示Ensembl到基因名注释对应。

## 2. rRNA注释信息 ##

分为两步：

第一步， “table”选择`wgEncodeGencodeBasicV19`{:.language-bash}。

<img src="/images/hg19_rRNA.png" title="image" alt="UCSC下载rRNA注释">

第二步，按照下图编辑“filter”。

<img src="/images/hg19_rRNA_maskTable.png" title="image" alt="UCSC下载rRNA注释Table">

<img src="/images/hg19_rRNA_mask.png" title="image" alt="UCSC下载rRNA注释筛选">


## 3. tRNA注释信息 ##

分为两步：

第一步，“track”选择 `tRNA Genes`{:.language-bash}。

<img src="/images/hg19_tRNA.png" title="image" alt="UCSC下载tRNA注释">

第二步，保留pseudo tRNA注释。

<img src="/images/hg19_tRNA_mask.png" title="image" alt="UCSC下载tRNA注释">


## 4. 线粒体基因组注释 ##

分为两步：

第一步， “table”选择`wgEncodeGencodeBasicV19`{:.language-bash}。

<img src="/images/hg19_chrM.png" title="image" alt="UCSC下载chrM注释">

第二步，按照下图编辑“filter”。

<img src="/images/hg19_chrM_mask.png" title="image" alt="UCSC下载chrM注释筛选">






### 参考网址 ###

* USCS Genome Browser的Google论坛：[1](https://groups.google.com/a/soe.ucsc.edu/forum/#!topic/genome/IL_aeOuPYU0)、[2](https://groups.google.com/a/soe.ucsc.edu/forum/#!msg/genome/jSAY8w1JVVo/P6lk4OJzDNEJ) 

* 另一种选择rRNA、tRNA和线粒体组注释的方法[How to get tRNA/rRNA/mitochondrial gene GTF file](http://onetipperday.blogspot.tw/2012/08/how-to-get-trnarrnamitochondrial-gene.html)

* [Extract rRNA and tRNA features from UCSC Browser](http://webappl.blogspot.tw/2015/02/extract-rrna-and-trna-features-from.html) 


### 更新记录 ###

2015年6月1日








