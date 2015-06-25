---
layout: post
title: "TopHat/Cufflinks/CummeRbund使用介绍"
date: 2014-06-26 14:51:56 +0800
comments: true
published: false
categories: Bioinfor
---

[TopHat](http://ccb.jhu.edu/software/tophat/index.shtml)、[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)和cummeRbund，被称为处理RNA-seq数据的“燕尾服（tuxedo）”。TopHat负责RNA-seq的reads映射比对到基因组，并且自动识别mRNA“内含子-外显子”剪切；Cufflinks擅长组装转录组和寻找差异表达基因（或转录起始位点TSS等）；cummeRbund主攻数据可视化。

<!--more-->

## 1. Tophat ##

**简介**：[TopHat](http://ccb.jhu.edu/software/tophat/index.shtml)是快速将RNA测序片段“对应（mapping）”到基因组上的工具，优势在于处理外显子间的剪切。内部首先使用bowtie或bowtie2把RNA测序片段“比对（alignment）”到基因组，之后再分析和鉴定剪切连接区域。

**平台**：Mac OS/Linux

**快速运行**：

{% codeblock lang:bash %}
# 双端测序
$ tophat2 -p 8 -o human_1 \
          --mate-inner-dist 165 --mate-std-dev 45 --no-mixed \
          hg19 human_1.fastq.gz human_2.fastq.gz 

# 单端测序
$ tophat2 -p 8 -o human_1 hg19 human_1.fastq.gz
{% endcodeblock %}


**重要参数解释**：

* `-p`{:.language-bash}：设置线程数，用于多核计算。

* `-o`{:.language-bash}：执行输入文件夹。

* `-r/--mate-inner-dist`{:.language-bash}：一对reads的间隔长度的期望平均值，默认值为50bp，**建议添加**。具体解释，参考[RNA-seq差异表达分析工作流程](http://blog.qiuworld.com:8080/archives/3007)和[Tophat中-r/–mate-inner-dist参数](http://www.plob.org/2012/12/04/4988.html)。

* `--mate-std-dev`{:.language-bash}：一对reads的间隔长度分布的标准差，默认值为20bp，**建议添加**。`-r/--mate-inner-dist`{:.language-bash}和`--mate-std-dev`{:.language-bash}的估计方法参考？？？？？？？

* `--library-type`{:.language-bash}：

> 测序仪器和方法，默认为标准Illumina平台的`fr-unstranded`{:.language-bash}。其他平台设置，详见[TopHat说明文档](http://ccb.jhu.edu/software/tophat/manual.shtml)、[How to tell which library type to use (fr-firststrand or fr-secondstrand)?](http://onetipperday.blogspot.sg/2012/07/how-to-tell-which-library-type-to-use.html)和[链特异性转录组原理](http://www.plob.org/2013/12/03/6731.html)。
>
> 如果分不清楚`fr-firststrand`{:.language-bash}和`fr-secondstrand`{:.language-bash}，推荐两种方法：[第一种](http://ccb.jhu.edu/software/tophat/faq.shtml)用两个参数试运行一个有1M reads的小样本，之后比较`junction.bed`{:.language-bash}大小；[第二种](http://onetipperday.blogspot.sg/2012/07/how-to-tell-which-library-type-to-use.html)在两个双端测序文件（`fastq.gz`{:.language-bash}）中抽取一些reads，之后[Blat](http://genome.ucsc.edu/cgi-bin/hgBlat?org=human)到USCS genomes上观察。

* `--no-discordant`{:.language-bash}：只对于paired reads，只报告concordant mappings。加入这个参数，tophat2在最后一步失败。也可以不加入这个参数，通过sam/bam文件第二列过滤discordant reads，方法参考[过滤TopHat分析双端测序的输出](http://yulongniu.bionutshell.org/blog/2015/05/16/filter-tophat2-output/)。

* `--no-mixed`{:.language-bash}：只对于paired reads，只报告paired reads都成功map。TopHat默认不加这个参数，即如果对于一个read，如果没有找到alignment的concordant或者discordant mate，那么这一对read将分别寻找和报道各自的alignment。这个参数与`--no-discordant`{:.language-bash}不同，因为加上`--no-mixed`{:.language-bash}也可能报道discordant pairs（例如一对reads都成功alignment，但是方向或者之间距离不对）。

* `-g/--max-multihits`{:.language-bash$}：对于多个map的reads，设定报告数目，默认数值为20。需要注意，尽管这个参数可以设定为1，也不能用于设定唯一map的read。因为某个read有可能map到基因组多个位点，当设定为1时，只会返回得分最高的情况。


**后续操作**：

* TopHat2运行后查看`align_summary.txt`{:.language-bash}获得比对结果。

* TopHat2会输出`accepted_hits.bam`{:.language-bash}（接受map的reads文件）和`unmapped.bam`{:.language-bash}（没有map上的reads文件）。对于后者，使用基因组浏览器，如[IGV](http://www.broadinstitute.org/igv/)或者[UCSC Genome Browser](http://genome.ucsc.edu/)大致看下是有无map，之后可以直接丢弃。

* 过滤双端测序的TopHat结果，参考[过滤TopHat分析双端测序的输出](http://yulongniu.bionutshell.org/blog/2015/05/16/filter-tophat2-output/)。


**补充**：

* TopHat2运行时，可以将基因组的fasta格式文件（比如`hg19.fa`{:.language-bash}），一起放在index的文件夹中（如果没有，Bowtie2先生成一个），可以节省运行时间。具体生成方法，参考：[二代测序中的短序列比对](http://yulongniu.bionutshell.org/blog/2014/07/26/short-sequence-alignment/)。


* [不加思考地使用默认参数的下场](http://www.acgt.me/blog/2015/4/27/the-dangers-of-default-parameters-in-bioinformatics-lessons-from-bowtie-and-tophat?utm_content=bufferb2c35&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer)



## 2. Cufflinks ##

**简介**：[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)是TopHat的下游工具，用于分析差异表达基因、差异转录起始位点、新基因和选择性剪切。一般可以分为三步：1. `cufflinks`{:.language-bash}对每个bam文件生成转录组；2. `cuffmerge`{:.language-bash}结合真实转录组和bam生成转录组，构建一个整合转录组；3. `cuffdiff`{:.language-bash}比较不同生物学样本，寻找差异表达基因。

**平台**：Mac OS/Linux

### 2.1 `cufflinks`{:.language-bash} 快速运行 ###

{% codeblock lang:bash Example code for cufflinks %}
$ cufflinks -p 8 -g hg19_ensembl.gtf -b hg19.fa -M hg19_rRNAtRNAchrM -u \
            -o outPutLinks accepted_filtered.bam
{% endcodeblock %}

**重要参数解释**：

* `-p/--num-threads`{:.language-bash$}：设置线程数，用于多核计算。

* `-o/--output-dir`{:.language-bash$}：输出文件夹。

* `-g/--GTF-guide`{:.language-bash$}：使用参考转录组注释文件指导组装过程，**建议添加**。加入这个参数和对应的注释文件，可以有效提高cufflinks转录本组装。输出的转录组是传入注释文件和cufflinks自行组装的联合。具体解释如下图[[Roberts et.al., 2011]](#Roberts et.al., 2011)：

<img src="/images/faux_reads_cufflinks.png" title="image" alt="cufflinks原理图">

* `-G/--GTF`{:.language-bash$}：只对传入的转录组注释文件进行表达量估计，即生成的结果如`genes.fpkm_tracking`{:.language-bash}只包括转录组注释文件，而不包括cufflinks自己可能的组装。如果RNA-seq分析的目的是为了获得差异表达基因，而不是寻找新的转录本，建议添加，以加快运行速度。

* `-b/--frag-bias-correct`{:.language-bash$}：提供一个multifasta文件，可以提高转录丰度的准确度，**建议添加**。

* `-M/–mask-file`{:.language-bash$}：需要去除的转录本，比如tRNA、rRNA和线粒体转录本，**建议添加**。在RNA-seq数据中，这些转录本量非常高。[cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html)提示剔除这些转录本，有助于提高转录量估计的准确程度。tRNA、rRNA和线粒体基因组注释信息提取，参考[UCSC Table下载注释文件](http://yulongniu.bionutshell.org/blog/2015/06/02/download-ucsc-gene-mask/)。

* `-u/--multi-read-correct`{:.language-bash$}：处理map到基因组多个位置的reads，**建议添加**。对于TopHat，建议过滤掉这部分reads，过滤方法见[过滤TopHat分析双端测序的输出](http://yulongniu.bionutshell.org/blog/2015/05/16/filter-tophat2-output/)。按照上述方法过滤这部分reads后，则不需要加这个参数。

* `--library-type`{:.language-bash$}：设定测序方法，默认为`fr-unstranded`{:.language-bash}，可以用于典型的Illumina双端（paired-end）测序。

### 2.2 `cuffmerge`{:.language-bash} 快速运行 ###

{% codeblock lang:bash Example code for cuffmerge %}
$ cuffmerge -p 8 -g hg19USCS_ensembl.gtf -s hg19.fa -o mergeFile assemblies.txt
{% endcodeblock %}

**重要参数解释**：

* `-p/--num-threads `{:.language-bash$}：设置线程数，用于多核计算。

* `-o/--output-dir`{:.language-bash$}：输出文件夹。

* `g/--ref-gtf`{:.language-bash$}：参考转录组，**建议添加**。

* `g/--ref-gtf`{:.language-bash$}：参考基因组序列，**建议添加**。

* `assemblies.txt`{:.language-bash$}：是一个文件，写入需要组装转录本的路径。


### 2.3 `cuffdiff`{:.language-bash} 快速运行 ###

{% codeblock lang:bash Example code for cuffd %}
$ cuffquant -p 14 -o quantOut -b hg19.fa \
            -M maskfile.gtf -u \
            merged.gtf h1.bam
{% endcodeblock %}

**重要参数解释**：

* `-p/--num-threads `{:.language-bash$}：设置线程数，用于多核计算。

* `-o/--output-dir`{:.language-bash$}：输出文件夹。

* `-b/--frag-bias-correct `{:.language-bash$}：参考基因组序列，**建议添加**。

* `-M/--mask-file `{:.language-bash$}：过滤序列GTF文件，放入这个文件中的序列都不会被计算，**建议添加**。

* `-u/--multi-read-correct`{:.language-bash$}：处理map到基因组多个位置的reads，**建议添加**。

### 2.4 `cuffdiff`{:.language-bash} 快速运行 ###

{% codeblock lang:bash Example code for cuffdiff %}
# cuffdiff可以输入cuffquant生成的cxb文件或者原始的bam文件
$ cuffdiff -p 14 -o diffOut -b hg19.fa \
           -M maskfile.gtf -u \
           -L P,R merged.gtf h1.cxb,h2.cxb,h5.cxb h2.cxb,h4.cxb,h6.cxb
{% endcodeblock %}

**重要参数解释**：

* `-p/--num-threads `{:.language-bash$}：设置线程数，用于多核计算。

* `-o/--output-dir`{:.language-bash$}：输出文件夹。

* `-b/--frag-bias-correct `{:.language-bash$}：参考基因组序列，**建议添加**。

* `-M/--mask-file `{:.language-bash$}：过滤序列GTF文件，放入这个文件中的序列都不会被计算，**建议添加**。

* `-u/--multi-read-correct`{:.language-bash$}：处理map到基因组多个位置的reads，**建议添加**。

* `-L/--labels `{:.language-bash$}：指定样品组名












### 参考资料 ###

* [一篇介绍cufflinks的中文博客](http://www.chenlianfu.com/?p=623)

* [cufflinks安装介绍](http://cole-trapnell-lab.github.io/cufflinks/getting_started/#common-uses-of-the-cufflinks-package)

* <a id="Roberts et.al., 2011">Roberts A</a>, Pimentel H, Trapnell C, Pachter L: **Identification of novel transcripts in annotated genomes using RNA-Seq.** *Bioinformatics*. 2011, 27(17):2325-9. [pdf](http://bioinformatics.oxfordjournals.org/content/early/2011/06/21/bioinformatics.btr355.full.pdf) 


### 更新记录 ###

2015年6月10日
