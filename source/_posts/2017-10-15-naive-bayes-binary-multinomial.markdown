---
layout: post
title: "朴素贝叶斯分类器应用于二元数据类型的尝试"
date: 2017-10-15 17:38:16 +0800
comments: true
categories: bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

### 简介 ###

朴素贝叶斯分类器是一种简单、有效的分类器，其难点在于估算条件概率。比如，一个数据集拥有$N$个相互独立的特征，$C$个分组，对于$C_j$条件概率模型为：

<!--more-->

$$
\begin{align}
\begin{split}
p(C_j|F_1,\cdots,F_n) &= \frac{p(F_1,\cdots,F_n|C_j)p(C_j)}{p(F_1,\cdots,F_n)} \\
&= p(F_1|C_j) \cdots p(F_n|C_j)p(C_j)(1/p(F_1,\cdots,F_n))
\end{split}
\label{eq:1}
\end{align}
$$

由于$1/p(F_1,\cdots,F_n)$在不同分组中为定值，因此：

$$
\begin{align}
\begin{split}
p(C_j|F_1,\cdots,F_n) &\propto p(C_j)\prod_{i=1}^{N}p(F_i|C_j)
\end{split}
\label{eq:2}
\end{align}
$$

其中，$p(C_j)$通常容易求得，即$C_j$分组在测试数据集中出现的频率。而$p(F_i\ \vert C_j)$则根据不同的测试数据类型，有不同的估计值。

以下讨论两种二元数据类型，例如某个数据集有三种特征量：

$$
F = 
\left[
\begin{array}{f}
F_1\\
F_2\\
F_3
\end{array}
\right]
$$

### 伯努利分布 ###

每一个特征量的取值都为$0$或$1$。分组$C_j$含有两个已知样本为：

$$
C_{j1} = 
\left[
\begin{array}{cj1}
0\\
1\\
0
\end{array}
\right]
$$

$$
C_{j2} = 
\left[
\begin{array}{cj2}
0\\
0\\
1
\end{array}
\right]
$$

某个预测样本为：

$$
C_{jp1} = 
\left[
\begin{array}{cjp1}
1\\
0\\
1
\end{array}
\right]
$$

由于$p(F_i \vert C_j)$不能为0，根据[Rule of succession](https://en.wikipedia.org/wiki/Rule_of_succession)得各个特征的条件概率为：

$$
\begin{align*}
\begin{split}
p(F1|C_j) &= \frac{0+1}{2+2} &= 1/4 \\
p(F2|C_j) &= \frac{1+1}{2+2} &= 1/2 \\
p(F3|C_j) &= \frac{1+1}{2+2} &= 1/2
\end{split}
\end{align*}
$$

### 二项分布 ###

每一个特征量的取值都一个元素为$0$或$1$的向量（长度可不等）。分组$C_j$含有两个已知样本为：

$$
C_{j1} = 
\left[
\begin{array}{cj1}
0 & 1 & 0 & 1\\
1 & 0 & 1\\
0 & 0
\end{array}
\right]
$$

$$
C_{j2} = 
\left[
\begin{array}{cj2}
0 & 1 & 1 & 1\\
1 & 1 & 1\\
0 & 0
\end{array}
\right]
$$

某个预测样本为：

$$
C_{jp1} = 
\left[
\begin{array}{cjp1}
0 & 0 & 1 & 1\\
1 & 0 & 1\\
0 & 0
\end{array}
\right]
$$

各个特征的条件概率为：

$$
\begin{align*}
\begin{split}
p(F1|C_j) &= \left(\frac{3+1}{8+2}\right)^2 \times \left(\frac{5+1}{8+2}\right)^2 \\
p(F2|C_j) &= \left(\frac{5+1}{6+2}\right)^2 \times \left(\frac{1+1}{6+2}\right) \\
p(F3|C_j) &= \left(\frac{4+1}{4+2}\right)^2
\end{split}
\end{align*}
$$

### 优化 ###

当特征较多时，会面临多个小数（$p$值）相乘。可以取对数后再相加，即$\sum\log{p}$。

### 参考资料 ###

* [Text Classification using Naive Bayes](https://www.inf.ed.ac.uk/teaching/courses/inf2b/learnnotes/inf2b-learn-note07-2up.pdf) 

* [Naive Bayes text classification](https://nlp.stanford.edu/IR-book/html/htmledition/naive-bayes-text-classification-1.html) 


