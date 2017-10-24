---
layout: post
title: "Bray-Curtis distance解释"
date: 2017-10-24 12:26:00 +0800
comments: true
categories: bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

$$
\newcommand{\sumup}[1] {\sum\limits_{i=1}^{n} #1}
$$

Bray-Curtis distance（BCD）的定义为：

<!--more-->

$$
\begin{align}
\begin{split}
BCD(X, Y) = \frac{\sumup{|x_i - y_j|}}{\sumup{x_i} + \sumup{y_i}}
\end{split}
\label{eq:1}
\end{align}
$$

其中，$X$和$Y$分别为长度为$n$的数值向量。根据$\eqref{eq:1}$可以得出：$BCD$的取值范围为$[0, 1]$；当$X$和$Y$完全相同时，$BCD$为0；反之，$BCD$为1。

同样，Bray-Curtis similarity（BCS）或Bray-Curtis index为：

$$
\begin{align}
\begin{split}
BCS(X, Y) = 1 - BCD(X, Y)
\end{split}
\label{eq:2}
\end{align}
$$

### 参考资料 ###

* [Chapter 6 Measures of distance and correlation between variables](http://84.89.132.1/~michael/stanford/maeb6.pdf)

### 更新记录 ###

2017年10月22日
