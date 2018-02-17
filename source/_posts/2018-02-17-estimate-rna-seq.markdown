---
layout: post
title: "通过RNA-Seq评估基因表达量的模型"
date: 2018-02-17 18:23:29 +0800
comments: true
categories: bioinfor 
published: false
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

$$
\newcommand{\tildel}{\widetilde{l_t}}
$$

$$
\newcommand{\P}{\mathrm{P}}
$$

本文基于[参考资料1](#Ref)，展示RNA-Seq在评估基因表达量模型的细节。

## 1. 符号表示 ##

$K$个长度为$l_i$的转录序列$t_i$，构成转录本的集合$T=\\{t_1, t_2, \dots, t_K\\}$。单个转录组中，每个转录序列$t_i$有$c_i$个拷贝数，全部转录序列的总拷贝数为$M$。单个转录序列的相对丰度为$\rho_i=\frac{c_i}{\sum\limits_{t \in T}c_t} = \frac{c_i}{M}$，易得$\sum\limits_{i=1}^K \rho_i= 1$。

单个转录组中，全部转录片段构成集合$F$，总转录片段数目为$N=\|F\|$。比对到的转录序列$t_i$上的转录片段，构成集合$F_t \in F$，对应的转录片段数目为$X_t=\|F_t\|$。
<!--more-->

## 2. 简单模型 ##

简单模型为：单端RNA-Seq，每一个read只比对到一个转录序列上，且每个read的长度都为定值$m$。对于转录序列$t_i$，从`5'`到`3'`一共可能比对上的read数目为$\tildel = \widetilde{l_i} - m + 1$。建立模型的思路是：当给定一个read，它会被比对到某个转录序列的某个位置是一个随机事件。通过实际观测（即将read比对到转录序列），进而估计未知参数$\rho = \\{\rho_1, \rho_2, \dots, \rho_K\\}$。

通过read序列比对，可得观测数据类似如下矩阵。每一行表示某个read是否比对到某个转录序列的某个位置，行和为1。

$$
\left[
\begin{matrix}
0 & 0 & \cdots & 1 \\
0 & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1 \\
\end{matrix}
\right]
$$

对于某个read $f$，来自于转录序列$t$的概率为：

$$
\begin{align}
\begin{split}
\P(f \in t) &= \frac{\rho_t M \tildel}{\sum\limits_{s \in T} \rho_s M \widetilde{l_s}} \\
&= \frac{\rho_t \tildel}{\sum\limits_{s \in T} \rho_s \widetilde{l_s}} \\
&= \alpha_t
\end{split}
\label{eq:1}
\end{align}
$$

当$f$来自于转录序列$t$时，$f$比对该转录序列某个位置的概率为：

$$
\begin{align}
\begin{split}
\P(\mathrm{pos}|f \in t) = \frac{1}{\tildel}
\end{split}
\label{eq:2}
\end{align}
$$

联合$\eqref{eq:1}$和从$\eqref{eq:2}$，对于$f$来自于转录序列$t$的某个位置概率为：

$$
\begin{align}
\begin{split}
\P(\mathrm{pos}, f \in t) = \frac{\alpha_t}{\tildel}
\end{split}
\label{eq:3}
\end{align}
$$

因此，极大似然函数为：

$$
\begin{align}
\begin{split}
L &= \prod_{t \in T} \prod_{f \in F_t} \frac{\alpha_t}{\tildel} \\
&= \prod_{t \in T} \left( \frac{\alpha_t}{\tildel} \right)^{X_t} \\
&\propto \prod_{t \in T} \alpha_t^{X_t}
\end{split}
\label{eq:4}
\end{align}
$$

在约束条件$\sum\limits_{t \in T} \alpha_t= 1$，求得极大似然估计为$\alpha_t = \frac{X_t}{N}$。有趣的是，在简单模型条件下，该极大似然估计可以来源于multinoulli分布。




### <a id="Ref">参考资料</a> ###

* Lior Pachter: Models for transcript quantification from RNA-Seq. [arXiv:1104.3889v2](https://arxiv.org/abs/1104.3889) [q-bio.GN], 2011.

* Wing-Kin Sung: Algorithms for Next-Generation Sequencing. [Chapman and Hall/CRC](https://www.crcpress.com/Algorithms-for-Next-Generation-Sequencing/Sung/p/book/9781466565500), 2017. 

### 更新记录 ###

2018年2月17日
