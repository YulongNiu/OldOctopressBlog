---
layout: post
title: "确定TopHat中mate-inner-dist参数"
date: 2015-05-23 19:08:19 +0800
comments: true
published: true
categories: Bioinfor
---

对于双端测序RNA-seq数据，TopHat在运行时候，有两个参数`-r/--mate-inner-dist`{:.language-bash}和`--mate-std-dev`{:.language-bash}分别标识一对reads的间隔长度的期望平均值和标准差，其默认值分别为50bp和20bp。这两个参数本身是个估计值，用于TopHat在map过程中确定一对reads是否匹配到基因组正确位置。如果能够准确设定这两个数值，将会提升TopHat结果的准确性和完整性。

有两种方法获得这对参数的准确值：

第一种：获取RNA-seq实验建库方法，之后按照以下网址说明计算，[RNA-seq差异表达分析工作流程](http://blog.qiuworld.com:8080/archives/3007)。

<!--more-->

第二种：根据RNA-seq数据进行估算，具体步骤为：


1. 使用TopHat默认参数先跑一遍。

2. 使用[MISO](http://miso.readthedocs.org/en/fastmiso/#computing-the-insert-length-distribution-and-its-statistics)的`pe_utils`{:.language-bash}工具估算。


以下详细介绍`pe_utils`{:.language-bash}使用方法。


第一步， 下载对应物种的基因注释文件GTF或者GFF，比如[USCS Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)（`output format`{:.language-bash}选择`GTF`{:.language-bash}）或者使用MISO提供的[Ensembl版本](http://miso.readthedocs.org/en/fastmiso/#human-mouse-gene-models-for-isoform-centric-analyses)。如果GTF文件，使用[Cufflinks](https://cole-trapnell-lab.github.io/cufflinks/file_formats/)的`gffread`{:.language-bash}工具进行转换。

第二步，确定TopHat运行结果的bam文件与基因注释GFF文件，两者基因组命名方法一致。有的使用类似`chr1`{:.language-bash}命名，而另外一些使用`1`{:.language-bash}。如果不一致，建议修改GFF文件。

{% codeblock lang:bash Summary Chromosomes' Names%}
# 查看GFF文件中基因组命名
$ awk '{print $1}' hg19USCS.gff | sort -n | uniq -c

# 查看bam文件中基因组命名
samtools view accepted_hits.bam | head -1000 | awk '{print $3}' | sort -n | uniq -c
{% endcodeblock %}

第三步，筛选较长外显子，比如长度大于1000bp。MISO提供了`exon_utils`{:.language-bash}工具用于提取长外显子，但是我们没有能够成功运行过。因此这里提供一个R版本的脚本，比如基因注释文件名为`hg19USCS.gff`{:.language-bash}，输出筛选的文件名为`hg19USCS_selected.gff`{:.language-bash}。

{% codeblock lang:R R Script of Select exons with length > 1000bp%}
gffFile <- read.table('hg19USCS.gff', stringsAsFactors = FALSE)

gffExon <- gffFile[gffFile[, 3] == 'exon', ]

exonLen <- abs(gffExon[, 5] - gffExon[, 4])

gffExonSelect <- gffExon[exonLen >= 1000, ]

write.table(gffExonSelect, 'hg19USCS_selected.gff', 
            row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = '\t')
{% endcodeblock %}


第四步，使用`pe_utils`{:.language-bash}，实例如下：

{% codeblock lang:bash An example of how to use pe_utils%}
# 输入bam文件和GFF文件
$ pe_utils --compute-insert-len accepted_hits.bam hg19USCS_selected.gff --output-dir insert-dist/

# 在insert-dist会出现类似accepted_filtered.bam.insert_len文件
# -r/--mate-inner-dist估计值为mean
# --mate-std-dev估计值为sdev
$ head -1 accepted_filtered.bam.insert_len
#mean=165.3,sdev=45.2,dispersion=3.5,num_pairs=5622239
{% endcodeblock %}


### 更新记录 ###

2015年5月23日
