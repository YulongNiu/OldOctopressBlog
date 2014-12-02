---
layout: post
title: "构建进化树的算法"
date: 2014-09-17 14:14:59 -0400
comments: true
published: false
categories: bioinfor
---

## 1. 基础知识介绍 ##

* **生物界分类**

> 界（Kingdom）
> 
> 门（Phylum）
> 
> 纲（Class）
> 
> 目（Order）
> 
> 科（Family）
> 
> 属（Genus）
> 
> 种（Species）

<!--more-->

* **进化树（Phylogenetic Tree）**：进化树是一个包括了节点和分支的树图，其中每一个分支上链接了两个节点（A phylogenetic tree is a graph composed of nodes and branches, in which only one branch connects any two adjacent nodes.）。“有根树（rooted tree）”不仅包括了物种（DNA序列或者氨基酸序列）的亲缘关系，同时也显示了具体的进化时间。“无根树（unrooted tree）”只能显示物种亲缘关系。

假设有$N$个物种，所有可能的有根树和无根树的数目分别是：

$$
\begin{align}
N_{root} &= \frac{(2n-3)!}{2^{n-2}(n-2)!} \\
N_{unroot} &= \frac{(2n-5)!}{2^{n-3}(n-3!)}
\end{align}
$$

* **分子时钟理论（Molecular Clock Theory）**：

* OTU：Operational Taxonomic Unit.





### <a id="Ref">参考网址</a> ###

* 幻灯片：[1](https://www.site.uottawa.ca/~lucia/courses/5126-10/lecturenotes/14InferingPhylogeny_Turcotte.pdf)，[2](http://people.cs.missouri.edu/~chengji/slides9.pdf)，[3](http://www.comp.nus.edu.sg/~ksung/algo_in_bioinfo/slides/Ch7_phylogeny.pdf)，[4](http://www.comp.nus.edu.sg/~ksung/algo_in_bioinfo/slides/Ch8_phylogeny_comparison.pdf)







### 更新记录 ###

2014年9月19日
