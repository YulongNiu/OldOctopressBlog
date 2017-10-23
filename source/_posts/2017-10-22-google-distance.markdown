---
layout: post
title: "Normalized Google distance解释"
date: 2017-10-22 21:35:30 +0800
comments: true
categories: bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

本文尝试探索Normalized Google distance（简称NGD）的定义和拓展应用。

### 1. NGD原始定义 ###

[维基百科](https://en.wikipedia.org/wiki/Normalized_Google_distance)的定义为：

<!--more-->

$$
\begin{align}
\begin{split}
NGD(x, y) = \frac{\max\{\log f(x), \log f(y)\} - \log f(x, y)}{\log N - \min\{\log f(x), \log f(y)\}}
\end{split}
\label{eq:1}
\end{align}
$$

其中，$f(x)$和$f(y)$分别为关键词$x$和$y$出现的次数，$f(x,y)$为$x$和$y$同时出现的次数，$N$为全部搜索单词数目。根据$\eqref{eq:1}$可以得出：如果$x$和$y$几乎总是同时出现时，$NGD$趋近于$0$；如果$x$和$y$出现的次数很少，即$\log f(x,y)$趋近于负无穷，则$NGD$可能大于$1$。

### 2. NGD定义延伸 ###

Choi and Rashid在2008年的文章（参考资料1）提出一种针对向量的$NGD$定义：

$$
\begin{align}
\begin{split}
NGD(X, Y) &= \frac{\max\{\sum X, \sum Y\} - \sum \min(X, Y)}{\sum X + \sum Y - \min\{\sum X, \sum Y\}} \\
&= \frac{\max\{\sum X, \sum Y\} - \sum \min(X, Y)}{\max\{\sum X, \sum Y\}}
\end{split}
\label{eq:2}
\end{align}
$$

其中，$X$和$Y$分别为两个等长的数值向量，$\min\{(X, Y)\}$为$X$和$Y$中各个元素最小值所组成的数值向量。根据$\eqref{eq:2}$可以得出：$NGD$的取值范围为$[0, 1]$；当$X$和$Y$完全相同时，$NGD$为0；反之，$NGD$为1。

由此，可以得到normalized Google similarity（NGS）为：

$$
\begin{align}
\begin{split}
NGS(X, Y) &= 1 - NGD(X, Y) \\
&= \frac{\sum \min(X, Y)}{\max\{\sum X, \sum Y\}}
\end{split}
\label{eq:3}
\end{align}
$$

一个例子：$X = (1, 2, 0, 3)$和$Y = (0, 2, 1, 1)$，则$NGD = 0.5$和$NGS = 0.5$。

### 参考资料 ###

* [Adapting Normalized Google Similarity in Protein Sequence Comparison](https://www.cs.cmu.edu/afs/cs/user/aberger/www/html/tutorial/tutorial.html)

### 更新记录 ###

2017年10月21日
