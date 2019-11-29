---
layout: post
title: "Agilent的microRNA芯片分析"
date: 2014-12-31 16:51:28 -0500
comments: true
published: false
categories: bioinfor
---

## 1. 原理 ##

Agilent的[microRNA芯片](http://www.genomics.agilent.com/zh/miRNA-Microarrays/?pgid=AG-PG-5)采用的是茎环探针设计，具体原理是：成熟的miRNA（microRNA）在3'端进行荧光标记，标记位点是胞嘧啶（Cytosine，`C`）；成熟的miRNA会与芯片探针配对，形成独特的颈环结构，如下图（图片取自[博奥生物](http://bioservices.capitalbio.com/fwpt/Agilentpt/3984.shtml)）。

<img src="/images/agilent_microrna_how.jpg" width="500" height="500" title="image" alt="images">

<!--more-->

## 2. <span style="color: blue">AgiMicroRna</span>包 ##

<span style="color: blue">AgiMicroRna</span>是[R/Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/AgiMicroRna.html)中一个用于Agilent的MiRNA芯片数据前处理和差异表达miRNA筛选的包。

### 2.1 数据类型  ###

* Target File

Target File（索引文件）是一个需要用户自己建立的表格文件，其中包括了芯片原始数据（比如分组信息、分组编号、文件名称等）和生物学特征（比如性别、年龄、疾病等）。第一列是

以下是一个索引文件的示例：

{% codeblock lang:r %}
# tf是一个索引文件，以下文件展示了4个对照组和4个处理组
> tf[c(1:4, 11:14), ]
          FileName Treatment GErep Labnum Dye
Control_1   36_1_1   control     1     N1 Cy3
Control_2   36_2_1   control     1     N3 Cy3
Control_3   36_1_2   control     1     N4 Cy3
Control_4   36_2_2   control     1     N7 Cy3
G1_1        37_1_1        G1     2     T1 Cy3
G1_2        37_2_1        G1     2     T2 Cy3
G1_3        37_1_2        G1     2     T3 Cy3
G1_4        37_2_2        G1     2     T4 Cy3
{% endcodeblock %}

### 2.2 原始数据载入 ###

Agilent会使用自己开发的[AFE(Agilent Feature Extraction)](http://www.chem.agilent.com/Library/usermanuals/Public/G4460-90026_FE_Reference.pdf)软件，读取和前处理扫描的芯片图像数据。

读取后的芯片原始数据，被储存在`uRNAList`对象中，对象名称和解释如下（图片取自<span style="color: blue">AgiMicroRna</span>包[帮助文档](http://www.bioconductor.org/packages/release/bioc/vignettes/AgiMicroRna/inst/doc/AgiMicroRna.pdf)）：

<img src="/images/agilent_microrna_uRNAList_names.png" width="500" height="500" title="image" alt="images">

* `gMeanSignal`：原始探针数值。

* `gProcessedSignal`：在经过AFE芯片背景矫正后的数值。

* `gTotalProbeSignal`：所有重复探针的`gProcessedSignal`求一个robust average，详细可以参考[AFE帮助文档](http://www.chem.agilent.com/Library/usermanuals/Public/G4460-90026_FE_Reference.pdf)。

* `gTotalGeneSignal`：对同一个基因的所有不同的`gTotalProbeSignal`求和，即将不同的`gTotalProbeSignal`值加起来。这个数值在经过一些标准化处理之后（非必须），可以用于差异表达miRNA的筛选。一个展示`gTotalProbeSignal`和`gTotalGeneSignal`例子如下：

{% codeblock lang:r %}
# miRNA“hsa-miR-1”有两类探针组成，名称分别是“A_25_P00012149”和“A_25_P00012150”。
> tmp1 <- which(brainMiRNA$gene$GeneName == 'hsa-miR-1')
> tmp2 <- which(brainMiRNA$gene$ProbeName == 'A_25_P00012149')
> tmp3 <- which(brainMiRNA$gene$ProbeName == 'A_25_P00012150')

# A_25_P00012149的gProcessedSignal、gTotalProbeSignal和gTotalGeneSignal数值
> brainMiRNA$procS[tmp2, 1]
 [1] 23.29449 16.74556 22.56444 19.90511 29.57939 21.02361 23.14942 23.78994
 [9] 28.95353 25.48756 24.45641 26.21296 23.19507 24.81994 24.98999
> brainMiRNA$TPS[tmp2, 1]
 [1] 6.97601 6.97601 6.97601 6.97601 6.97601 6.97601 6.97601 6.97601 6.97601
[10] 6.97601 6.97601 6.97601 6.97601 6.97601 6.97601
> brainMiRNA$TGS[tmp2, 1]
 [1] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772
[10] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772

# A_25_P00012150的gProcessedSignal、gTotalProbeSignal和gTotalGeneSignal数值
> brainMiRNA$procS[tmp3, 1]
 [1] 102.72520 117.78560  95.13415 114.89230 109.46910  99.82609 102.22190
 [8] 114.36390 110.06570 113.38240 114.26590 109.54290 123.76020 102.68770
[15] 112.90630
> brainMiRNA$TPS[tmp3, 1]
 [1] 32.0012 32.0012 32.0012 32.0012 32.0012 32.0012 32.0012 32.0012 32.0012
[10] 32.0012 32.0012 32.0012 32.0012 32.0012 32.0012
> brainMiRNA$TGS[tmp3, 1]
 [1] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772
[10] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772

# hsa-miR-1的gTotalGeneSignal数值
> brainMiRNA$TGS[tmp1, 1]
 [1] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772
[10] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772
[19] 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772 38.9772
[28] 38.9772 38.9772 38.9772
{% endcodeblock %}




### 2.3 归一化处理 ###

根据[[López-Romero et.al., 2010]](#López-Romero et.al., 2010)比较了各种不同的归一化方法，最后有两种归一化方法有较好的表现。第一种是使用AFE算法得到的`gTotalGeneSignal`，简称为TGS;第二种方法是使用“无背景矫正的RMA（Robust Multiarray Average）算法”，简称RMA。TGS和RMA算法在实际操作中都有出色的表现，其中RMA算法对于探测低信号探针更加敏感。<span style="color: blue">AgiMicroRna</span>包提供了这两种方法。

* TGS

由于TGS可能负值，可以使用：1. 将所有TGS都加上一个小的正数（offset）；2. 将所有小于0.5的TGS都设置位0.5（half）。使用函数：

{% codeblock lang:r %}
# `dd`：uRNAList类，可以从readMicroRnaAFE()函数读取原始芯片数据得到。
# `offset`：使用加小正数方法。如果使用该方法，设置`half = FALSE`。
# `half`：如果设置为`TRUE`，则使用0.5为闸值。
# `makePLOT`：否绘制TGS的QC图。
# `verbose`：否打印出一些程序运行信息。
# 注意，此时得到的TGS数据是原始数值，没有经过log2转换。
tgsMicroRna(dd, offset, half, makePLOT=FALSE, verbose=FALSE)
{% endcodeblock %}

之后，对TGS进行组间归一化。组间归一化的目的是为了去除“芯片-芯片”之间的非生物学差异，比如芯片上样、杂交等。这些非生物学差异，可能会掩盖芯片之间真正的生物学差异。TGS组间归一化使用函数：

{% codeblock lang:r %}
# `ddTGS`：通过tgsMicroRna()函数提取的包含TGS的对象。
# `NORMmethod`：包括了三种方法，“quantile”、“scale”和“none”，之后将归一化数值做log2处理。
# `makePLOTpre`和`makePLOTpost`：否绘制处理前和处理后的密度图、盒箱图、MA图、RLE图和生物学样本聚类图。
# `targets`：索引文件。
# `verbose`：否打印出一些程序运行信息。
tgsNormalization(ddTGS, NORMmethod = "quantile", makePLOTpre = FALSE, makePLOTpost = FALSE, targets, verbose=FALSE)
{% endcodeblock %}

* RMA

MiRNA芯片也可以使用[[Gautier et.al., 2004]](#Gautier et.al., 2004)提出的RMA算法，包括三个步骤：1. 背景去除（可选，MiRNA芯片分析不推荐）；2. quantile组间归一化（可选）；3. 使用线性模型将探针信号转换成log2形式。使用函数：

{% codeblock lang:r %}
# `dd`：uRNAList类，可以从readMicroRnaAFE()函数读取原始芯片数据得到。
# `normalize`：否进行quantile归一化处理。
# `background`：否进行背景去除，MiRNA芯片分析建议将其设置为FALSE。
rmaMicroRna(dd, normalize, background)
{% endcodeblock %}

在进行TGS或RMA标准化处理后，可以进行“探针过滤”。探针过滤分为两类：1. 探针信号是否被探测到，使用的是由AFE得到的`gIsGeneDetected`，该值为1或者0，分别表示是否被探测到；2. 探针信号值是否达到最低限度，最低限度的设定使用“负对照”探针信号。使用函数：

{% codeblock lang:r %}
# `ddNORM`：TGS或者RMA标准化处理后的探针密度。
# `dd`：uRNAList类，可以从readMicroRnaAFE()函数读取原始芯片数据得到
# `control`：是否删除对照探针。
# `IsGeneDetected`：是否根据gIsGeneDetected数值过滤探针。
# `wellaboveNEG`：是否过滤未达到最低信号标准的探针。
# `limIsGeneDetected`：比例数值，表示所有样品中，至少有多少比例探针的gIsGeneDetected是1。
# `limNEG`：比例数值，表示所有样品中，至少有多少比例探针需要高于最低信号标准。
# `makePLOT`：是否绘图。
# `targets`：索引文件。
# `verbose`：否打印出一些程序运行信息。
# `writeout`：是否生成过滤探针相关文件。
filterMicroRna(ddNORM, 
               dd,
               control,
               IsGeneDetected,
               wellaboveNEG,
               limIsGeneDetected,
               limNEG,
               makePLOT,
               targets,
               verbose,
               writeout) 
{% endcodeblock %}

最后，将预处理后的探针转换为<span style="color: blue">Biobase</span>包的`ExpressionSet`对象，使用函数：

{% codeblock lang:r %}
# `uRNAList`：预处理后的uRNAList对象。
# `targets`：索引文件。
# `makePLOT`：是否生成前100个差异最大探针的heatmap、生物学样本PCA和聚类图。
esetMicroRna(uRNAList, targets, makePLOT=FALSE,verbose=FALSE)
{% endcodeblock %}

### 2.4 可视化展示 ###










### 参考资料 ###

* miRNA芯片不同归一化方法比较 <a id="López-Romero et.al., 2010"> López-Romero P</a>, González MA, Callejas S, Dopazo A, Irizarry RA: **Processing of Agilent microRNA array data.** *BMC Res Notes.* 2010, 3:18.

* AgiMicroRna包 López-Romero P: **Pre-processing and differential expression analysis of Agilent microRNA arrays using the AgiMicroRna Bioconductor library.** *BMC Genomics.* 2011, 12:64.

* <a id="Gautier et.al., 2004"> Gautier L</a>, Cope L, Bolstad BM, Irizarry RA: **affy--analysis of Affymetrix GeneChip data at the probe level.** *Bioinformatics.*  2004, 20(3):307-15.
