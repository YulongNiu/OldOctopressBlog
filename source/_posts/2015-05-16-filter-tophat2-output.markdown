---
layout: post
title: "过滤双端测序的TopHat结果"
date: 2015-05-16 04:33:10 +0800
comments: true
published: false
categories: Bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

## 0. 结论 ##



## 1. TopHat输出sam文件的第五列 ##

[TopHat文档](http://ccb.jhu.edu/software/tophat/manual.shtml)没有解释其输出bam文件（比如`accepted_hits.bam`{:.language-bash}）的第五列的意义。按照[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)输出结果来看，是表示映射质量指标Mapping Quality scores（MAPQ），具体计算参考公式$\eqref{eq:1}$。

$$
\begin{align}
MAPQ = -10 \times log_{10}(pvalue)
\label{eq:1}
\end{align}
$$

<!--more-->

MAPQ值越大，表示对应的read的alignment质量越高。然而，在TopHat输出结果中，MAPQ所代表意义略有不同。

{% codeblock lang:bash "MAPQ" distribution from TopHat2 accepted mapping reads%}
# 查看accepted_hits.bam文件的MAPQ数值，并统计出现频数
$ samtools view accepted_hits.bam | awk '{print $5}' | sort --parallel=4 -n | uniq -c
5057430 0
3117731 1
8058500 3
93044727 50

# 查看前100位MAPQ数值和NH:i:n分布
$ samtools view accepted_hits.bam | awk '{print $5}' | head -100 | sort --parallel=4 -n | uniq -c
35 0
42 1
22 3
1 50

$ samtools view accepted_hits.bam | awk '{print $20}' | head -100 | sort --parallel=4 -n | uniq -c
1 NH:i:1
22 NH:i:2
3 NH:i:20
13 NH:i:3
29 NH:i:4
5 NH:i:5
21 NH:i:6
1 NH:i:7
2 NH:i:8
3 NH:i:9

# 查看前200位MAPQ数值和NH:i:n分布
$ samtools view accepted_hits.bam | awk '{print $5}' | head -200 | sort --parallel=4 -n | uniq -c
60 0
43 1
67 3
30 50

$ samtools view accepted_hits.bam | awk '{print $20}' | head -200 | sort --parallel=4 -n | uniq -c
30 NH:i:1
8 NH:i:12
7 NH:i:14
67 NH:i:2
3 NH:i:20
13 NH:i:3
30 NH:i:4
5 NH:i:5
31 NH:i:6
1 NH:i:7
2 NH:i:8
3 NH:i:9
{% endcodeblock %}

首先，我们可以看到TopHat输出的MAPQ只有四个数值，分别为`50`{:.language-bash}、`3`{:.language-bash}、`1`{:.language-bash}和`0`{:.language-bash}。根据[sam文件标准](http://samtools.github.io/hts-specs/SAMv1.pdf)，`NH:i:n`{:.language-bash}表示含有查询序列的alignment的数量。因此，通过上述前100位和前200位分析可以发现，MAPQ并不是按照公式$\eqref{eq:1}$计算，而有可能是以下关系


-----------------

|---------------+---------------+--------------+---------|
|**X/Y**        |**1(Presence)**|**0(Absence)**|**Sum**  |
|:--------------|:-------------:|:------------:|--------:|
|**1(Presence)**|a              |b             |a+b      |
|---------------|---------------|--------------|---------|
|**0(Absence)** |c              |d             |c+d      |
|---------------|---------------|--------------|---------|
|**Sum**        |a+c            |b+d           |n=a+b+c+d|
|---------------|---------------|--------------|---------|

-----------------



## 2. sam文件的第二列 ##

[快速解释第二列的bitwise FLAG](http://broadinstitute.github.io/picard/explain-flags.html)


## 3. 操作TopHat输出的bam文件命令集锦 ##


{% codeblock lang:bash Useful bash for bam files from TopHat%}
# samtools的view -c命令，其实就是输出sam文件有多少行
$ samtools view accepted_hits.bam | wc -l
109278388
$ samtools view -c accepted_hits.bam
109278388

{% endcodeblock %}









### 参考网址 ###

* [关于map当中的unique mapped reads问题](http://blog.qiuworld.com:8080/archives/3321)

* TopHat的bam输出文件第五列（类似MAPQ）的讨论： [tophat mapping quality](https://groups.google.com/forum/#!topic/tuxedo-tools-users/m0p1qXDEqKA)和[More madness with MAPQ scores](http://www.acgt.me/blog/2015/3/17/more-madness-with-mapq-scores-aka-why-bioinformaticians-hate-poor-and-incomplete-software-documentation)








### 更新记录 ###

2015年5月25日
