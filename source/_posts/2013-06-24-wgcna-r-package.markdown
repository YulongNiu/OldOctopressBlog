---
layout: post
title: "使用WGCNA构建权重共表达网络"
date: 2013-06-24 16:34:12 -0500
comments: true
published: false
categories: r
---

[WGCNA](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)（Weighted correlation network analysis）提供了一个R的包[[Langfelder and Horvath, 2008]](#Langfelder and Horvath, 2008)，用来构建权重共表达网络。权重共表达网络（co-expression weighted network）是一个无向有度（undirected, weighed network）的网络。一般分为两种网络：

* unsigned网络，网络节点的相关度为

$$
a_{ij} = \lvert cor(x_i, x_j) \rvert^\beta
$$

* signed网络，网络节点的相关度为

$$
a_{ij} = \left\lvert
\frac{1 + cor(x_i, x_j)}{2}
\right\rvert^\beta
$$

与权重网络相对应的是“无权重网络（unweighted network）”，节点与节点之间的相关度只能是`0`或者`1`。`0`表示两个节点没有连接，而`1`表示有。可以看到，权重网络对“节点-节点”之间的判断比较严苛。比如，如果阀值选择为0.8，那么两个相关度0.79节点会被判断为`0`。





















### 参考资料 ###

* WGCNA的R包参考文献 <a id="Langfelder and Horvath, 2008">Langfelder P and Horvath S</a>: **WGCNA: an R package for weighted correlation network analysis.** *BMC Bioinformatics.* 2008, 9:559.
