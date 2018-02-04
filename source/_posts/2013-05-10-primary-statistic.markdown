---
layout: post
title: "统计学基本知识汇总"
date: 2013-05-10 00:29:47 +0800
comments: true
categories: bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

$$
\newcommand{\P}{\mathrm{P}}
$$


## 1. 完备事件 ##

对于完备事件组$X = \\{x_1, x_2, \dots, x_n\\}$：

$$
\begin{align}
\begin{split}
\P(x_1) + \P(x_2) + \dots + \P(x_n) = 1
\end{split}
\label{eq:1}
\end{align}
$$

常用的技巧构造乘法系数，例如$\P(Y) = \sum\limits_{i=1}^{n}P(x_i\|\theta)P(Y)$

## 2. 全概率公式 ##

对于完备事件组$X = \\{x_1, x_2, \dots, x_n\\}$，事件$Y$的全概率公式：

$$
\begin{align}
\begin{split}
\P(Y) &= \P(Y, x_1) + \P(Y, x_2) + \dots + \P(Y, x_n) \\
&= \P(x_1)\P(Y|x_1) + \P(x_2)\P(Y|x_2) + \dots + \P(x_n)\P(Y|x_n)
\end{split}
\label{eq:2}
\end{align}
$$

使用概率密度函数表示为：

$$
\begin{align}
\begin{split}
f(Y) &= \int_{-\infty}^{+\infty}f(Y, x) \mathrm{d}x \\
&= \int_{-\infty}^{+\infty}f(Y|x)f(x) \mathrm{d}x
\end{split}
\label{eq:3}
\end{align}
$$

<!--more-->

## 3. 独立事件 ##

$X$和$Y$为两事件独立，则：

$$
\begin{align}
\begin{split}
\P(X, Y) &= \P(X)\P(Y|X) \\
&= \P(Y)\P(X|Y) \\
&= \P(X)\P(Y)
\end{split}
\label{eq:4}
\end{align}
$$

在朴素贝叶斯分类的推到中，可以用到诸如等式$P(X\|Y, \theta) = P(X\|\theta)$。

## 4. 贝叶斯公式 ##

根据$\eqref{eq:3}$和$\eqref{eq:4}$得：

$$
\begin{align}
\begin{split}
\P(x_1|Y) &= \frac{\P(x_1, Y)}{\P(Y)} \\
&= \frac{\P(x_1)\P(Y|x_1)}{\sum\limits_{i=1}^{n} \P(x_i)\P(Y|x_i)}
\end{split}
\label{eq:5}
\end{align}
$$


### 更新记录 ###

2018年2月3日
