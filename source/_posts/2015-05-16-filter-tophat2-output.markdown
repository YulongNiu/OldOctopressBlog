---
layout: post
title: "过滤TopHat分析双端测序的输出"
date: 2015-05-16 04:33:10 +0800
comments: true
styles: [data-table]
published: true
categories: Bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

## 0. 结论 ##

在使用TopHat2匹配双端测序结果后，建议根据成对reads的map基因组位置唯一、方向正确和距离合适的标准，筛选得到的匹配结果。比如，TopHat2可能生成`accepted_hits.bam`{:.language-bash}文件，处理方法如下：


{% codeblock lang:bash Filter TopHat Outputs%}

# 首先查看bam文件头部有多少行
$ samtools view -H accepted_hits.bam | wc -l
86

# 筛选成对且成功map到基因组唯一位置的reads，按照上一条输出结果，调整“NR <= 86”
$ samtools view -h accepted_hits.bam | \
    awk '{if (NR <= 86) print $0}; {if($5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99)) print $0}' | \
    samtools view -b - > accepted_filtered.bam

$ samtools view accepted_filtered.bam | wc -l
79143942
{% endcodeblock %}


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

首先，我们可以看到TopHat输出的MAPQ只有四个数值，分别为`50`{:.language-bash}、`3`{:.language-bash}、`1`{:.language-bash}和`0`{:.language-bash}。根据[sam文件标准](http://samtools.github.io/hts-specs/SAMv1.pdf)，`NH:i:n`{:.language-bash}表示含有查询序列的alignment的数量。因此，通过上述前100位和前200位分析可以发现，MAPQ并不是按照公式$\eqref{eq:1}$计算，而有可能是以下关系：


-----------------

| MAPQ (tophat) | Tag            | 描述              |
|:---------------+:----------------:+-------------------:|
|            50 | NH:i:1         | map至唯一位置     |
|            3 | NH:i:2         | map至2个位置      |
|             1 | NH:i:3/NH:i:4  | map至3个或4个位置 |
|             0 | NH:i:n (n > 4) | map到多余4个位置  |

-----------------

展示一个`NH:i:1`{:.language-bash}的例子，注意Illumina双端测序平台`fr-unstranded`{:.language-bash}：


{% codeblock lang:bash Two pairs of unique mapped reads from Illumina HiSeq2000 %}
HISEQ2000-02:436:C2PG3ACXX:3:2313:10972:95322   163     chr1    637224  50      100M    =       637339  215     AAATGATCTGTACAATAACCCCCTGTGACACCAGTCTACCTATGTAACAAATGCCCCTAAACTTAAAATAAAAGTTAAAAAAAAAGAAAATTAAAATCTC  <@@BDABBDFBCDGEGHIGIIGIABAFHBGDFGGGHIIIFIEIGGGGIIIFFDHIGIIIIIIIICEHIIIIIIHEHHDHFFCBCCB@BCAACCCCCECCC    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100     YT:Z:UU  NH:i:1
HISEQ2000-02:436:C2PG3ACXX:3:2313:10972:95322   83      chr1    637339  50      100M    =       637224  -215    GTAATATGAAAAACACAAATCTTTCATTCATTCCTTTCAACTGATGAGGAAAATGAGGCATCGGGAGTTAGTAAAAGTCCACATTGAGATATGAGACCCA  CCADDDCCCCCCCDEEEECAEHEEGGIIHGFAAGGGHEF=IGGGIIGGHGCGIEIIIGIIIIIIIFIHIIIIIIGIIHFEAIGGIIGFFDFHDDDAD@@@    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100     YT:Z:UU  NH:i:1

HISEQ2000-06:325:C2RC0ACXX:5:2205:6961:88285    99      chr1    643662  50      100M    =       643707  145     CCTATCAAAATCTTAGCATTCCTCTTAGCCCTCAACAAAGCATTTCTAAAATGTGTATAGAAGACCAAAGGGCCAAAAGAGTCAACTTCTGAAGAAGCGC  CCCFFFFFHHHHHJJJJIJJHJJJJJJJIIJJJJJJJJJJJJIJJJJJIIIJJIIIJJJJJJJJJIJJJJIJJJIHHHHFFEFFFEEEEEEEDDCDDCDD    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100     YT:Z:UU  NH:i:1
HISEQ2000-06:325:C2RC0ACXX:5:2205:6961:88285    147     chr1    643707  50      100M    =       643662  -145    CTAAAATGTGTATAGAAGACCAAAGGGCCAAAAGAGTCAACTTCTGAAGAAGCGCAAAAAGAAAGTTGAGGAAATCTTAAAACATGTTATTGAGCTTAAA  CEEEDDDDEDFEEEEEEEBFFFFFHHHHEJJJJJJJJIJJJJGJIIIIHFJJJIIJJJJJIJJJJJJJJJJJJGHJJJJIJJJJJJJGHHHHFFFFFCCC    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100     YT:Z:UU  NH:i:1
{% endcodeblock %}


## 2. sam文件的第二列 ##

sam文件中的第二列提供了具体的map情况，下列表格摘自[sam文件标准](http://samtools.github.io/hts-specs/SAMv1.pdf)，sam/bam文件中第二列各种条件**求和**的**十进制**标识，[快速解释第二列的bitwise FLAG](http://broadinstitute.github.io/picard/explain-flags.html)：

----------------

| Decimal | Hexadecimal | Description                                                        |
|:---------:+:-------------:+:--------------------------------------------------------------------|
|       1 |         0x1 | template having multiple segments in sequencing                    |
|       2 |         0x2 | each segment properly aligned according to the aligner             |
|       4 |         0x4 | segment unmapped                                                   |
|       8 |         0x8 | next segment in the template unmapped                              |
|      16 |        0x10 | SEQ being reverse complemented                                     |
|      32 |        0x20 | SEQ of the next segment in the template being reverse complemented |
|      64 |        0x40 | the first segment in the template                                  |
|     128 |        0x80 | the last segment in the template                                   |
|     256 |       0x100 | secondary alignment                                                |
|     512 |       0x200 | not passing quality controls                                       |
|    1024 |       0x400 | PCR or optical duplicate                                           |
|    2048 |       0x800 | supplementary alignment                                            |

----------------

在map完成双端测序序列中，我们感兴趣的是一对reads都正确align到基因组上，而且方相匹配又距离合适。符合这样条件的reads，对应的第二列数值为99、147、83和163，具体图示参考[Directional RNA-seq— Part 1: SAM file flags](https://biobeat.wordpress.com/2013/04/29/directional-rna-seq-part-1-extract-strand-information-from-sam-file/)。下面表格解释四个数值的具体意义，其中`1`{:.language-bash}标识双端测序，`2`{:.language-bash}表示一对reads正确地map到基因组合适位置，表格中着重陈述`64`{:.language-bash}、`32`{:.language-bash}、`128`{:.language-bash}和`16`{:.language-bash}。

-------------------

| Flag | Composition | Explanation                                                      |
|------+-------------+------------------------------------------------------------------|
|   99 |   64+32+2+1 | 一对引物中第一个map到基因组正义链；第二个反方向map到基因组正义链 |
|  147 |  128+16+2+1 | 一对引物中第二个反方向map基因组正义链；第一个map到基因组正义链   |
|   83 |   64+16+2+1 | 一对引物中第一个map到基因组反义链；第二个反方向map到基因组反义链 |
|  163 |  128+32+2+1 | 一对引物中第二个反方向map基因组正义链；第一个map到基因组正义链   |

-------------------

之后，我们需要筛选含有这些flags的reads。由于我们通常需要操作bam文件，也希望输出是bam文件，中间过程不希望再重新生成sam文件。那么，就需要结合使用`awk`{:.language-bash}进行筛选，具体方法见本篇文章开头所示。当然，如果是只是查看，可以使用下面例子中的 `samtools view -f 0x2`{:.language-bash}。

{% codeblock lang:bash Count number of unique pair-mapped alignments%}
# map到基因组上唯一位置的reads数目
$ samtools view -q 50 accepted_hits.bam | wc -l
93044727
# 成对reads都map到基因组对应位置的reads数目
$ samtools view -f 0x2 accepted_hits.bam | wc -l
88793640
# 成对且唯一mapped的reads数目
$ samtools view -q 50 -f 0x2 accepted_hits.bam | wc -l
79143942
{% endcodeblock %}


## 3. 操作TopHat输出的文件命令集锦 ##

{% codeblock lang:bash Useful bash for bam files from TopHat%}
# samtools的view -c命令，其实就是输出sam文件有多少行
$ samtools view accepted_hits.bam | wc -l
109278388
$ samtools view -c accepted_hits.bam
109278388

# 查看bam文件中mapped的reads长度分布
# 第二种方法是运行FastQC，输出结果中也有显示
samtools view accepted_hits.bam | awk '{print length($10)}' | sort -n | uniq -c

# 查看bed文件前几行
$ head junction.bed
# 统计bed文件有多少行，需要去除第一行注释
# 以下两种方式相同，但不够完美
$ wc -l junction.bed
220648 junctions.bed
$ awk 'END {print NR}' junctions.bed
220648
{% endcodeblock %}


## 4. DNA链和mRNA链的称呼总结 ##

双链DNA和单链mRNA，对每条链都有特定的称呼。总结如下：

{% codeblock lang:bash DNA/RNA strands %}
3'~~~~~UCUGAU~~~~~ 5' mRNA的对应基因信息在reverse strand

5'-----AGACTA----------ATTGTT----- 3'
3'-----TCTGAT----------TAACAA----- 5'

                5'~~~~~AUUGUU~~~~~ 3' mRNA的对应基因信息在forward strand
{% endcodeblock %}

对于一条双链DNA，称呼列表如下：

| 方向<sup>a</sup> | 名称1   | 名称2 |
|------------------+---------+-------|
| 从左至右         | forward | plus  |
| 从右至左         | reverse | minus |

<sup>a</sup>：方向是指5'至3'的阅读方向，用于区分两条DNA链条

----------------

对于一条RNA链，其对应的双链DNA称呼如下：

| mRNA方向<sup>a</sup> | 名称1    | 名称2       | 名称3     | 名称4    | 名称5 |
|----------------------+----------+-------------+-----------+----------+-------|
| 同向                 | coding   | nontemplate | sense     | positive | +     |
| 反向                 | template | noncoding   | antisense | negative | -     |

<sup>a</sup>：mRNA方向是指5'至3'。

----------------

### 参考网址 ###

* [关于map当中的unique mapped reads问题](http://blog.qiuworld.com:8080/archives/3321)

* TopHat的bam输出文件第五列（类似MAPQ）的讨论： [tophat mapping quality](https://groups.google.com/forum/#!topic/tuxedo-tools-users/m0p1qXDEqKA)和[More madness with MAPQ scores](http://www.acgt.me/blog/2015/3/17/more-madness-with-mapq-scores-aka-why-bioinformaticians-hate-poor-and-incomplete-software-documentation)

* [Directional RNA-seq— Part 2: Explore SAM flags using samtools](https://biobeat.wordpress.com/2013/05/01/directional-rna-seq-part-2-using-samtools/)

* [Filtering a SAM file generated by TopHat to find uniquely mapped, concordant read pairs: AWK vs SAMtools](http://www.acgt.me/blog/2015/4/15/filtering-a-sam-file-generated-by-tophat-to-find-uniquely-mapped-concordant-read-pairs-awk-vs-samtools)

* [Question: Forward And Reverse Strand Conventions](https://www.biostars.org/p/3423/)

* [wiki](http://en.wikipedia.org/wiki/Sense_(molecular_biology))



### 更新记录 ###

2015年5月25日
