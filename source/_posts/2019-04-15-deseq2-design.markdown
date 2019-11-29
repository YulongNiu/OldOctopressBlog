---
layout: post
title: "DESeq2包的对比组设计"
date: 2019-04-15 14:58:58 +0200
comments: true
categories: bioinfor
published: true
---

## 1. 两两比对 ##

组A *vs.* 组B。

{% raw %}
```
DataFrame with 12 rows and 2 columns
         condition genotype
          <factor> <factor>
sample1          A        I
sample2          A        I
sample3          A        I
sample4          A       II
sample5          A       II
...            ...      ...
sample8          B        I
sample9          B        I
sample10         B       II
sample11         B       II
sample12         B       II
```
{% endraw %}

<!--more-->

{% codeblock pair-wise comparison lang:r %}
dds <- makeExampleDESeqDataSet(n = 100, m = 12)
dds$genotype <- factor(rep(rep(c('I', 'II'), each=3), 2))

## condition: A vs. B
design(dds) <- ~ condition
ddres <- DESeq(dds)
res <- results(ddres, contrast = c('condition', 'B', 'A'))

## genotype I vs. II
design(dds) <- ~ genotype
ddres <- DESeq(dds)
res <- results(ddres, contrast = c('genotype', 'I', 'II'))

## A vs. B at genotype II
dds$group <- factor(paste0(dds$genotype, dds$condition))
design(dds) <- ~ group
ddres <- DESeq(dds)
results(ddres, contrast = c('group', 'IIB', 'IIA'))
{% endcodeblock %}

## 2. 交叉项 ##

{% codeblock pair-wise comparison lang:r %}
dds <- makeExampleDESeqDataSet(n = 100, m = 12)
dds$genotype <- factor(rep(rep(c('I', 'II'), each=3), 2))
design(dds) <- ~ genotype + condition + genotype:condition
ddres <- DESeq(dds) 

## A vs. B at genotype I
res <- results(ddres, contrast = c('condition', 'B', 'A'))

## A vs. B at genotype II
res <- results(ddres, list(c('condition_B_vs_A', 'genotypeII.conditionB')))

## condition effect *different* across genotypes
res <- results(ddres, name = 'genotypeII.conditionB') 
{% endcodeblock %}

其中，第二例子中的`A vs. B at genotype II`与第一个的区别是，考虑了交叉项的影响。如果只是为了两两比对，可以考虑使用第一个例子的处理方法。

### <a id="Ref">参考网址</a> ###

* [Analyzing RNA-seq data with DESeq2](http://master.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#nested-indiv)

* [Question: DESeq2 factor design vs pair-wise comparison](https://support.bioconductor.org/p/64352/)

### 更新记录 ###

2019年04月15日

