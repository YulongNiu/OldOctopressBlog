---
layout: post
title: "最大熵模型"
date: 2017-10-16 13:16:12 +0800
comments: true
categories: bioinfor
published: true
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

### 1. 熵和条件熵 ###

对于随机变量$X$，熵为：

$$
\begin{align}
\begin{split}
H(X) = -\sum_{x \in X}p(x)\log{p(x)}
\end{split}
\label{eq:1}
\end{align}
$$

<!--more-->

其中：

$$
\begin{align}
\begin{split}
\sum_{x \in X}p(x) = 1
\end{split}
\label{eq:2}
\end{align}
$$

同样道理，对于任意随机变量$X$和$Y$，联合熵为：

$$
\begin{align}
\begin{split}
H(X,Y) = -\sum_{x \in X, y \in Y}p(x,y)\log{p(x,y)}
\end{split}
\label{eq:3}
\end{align}
$$

基于$X$的$Y$的熵为条件熵：

$$
\begin{align}
\begin{split}
H(Y|X) &= H(X, Y) - H(X) \\
&= -\sum_{x \in X, y \in Y}p(x,y)\log{p(x,y)} + \sum_{x \in X}p(x)\log{p(x)} \\
&= -\sum_{x \in X, y \in Y}p(x,y)\log{p(x,y)} + \sum_{x \in X, y \in Y}p(x, y)\log{p(x)} \\
&= -\sum_{x \in X, y \in Y}p(x, y)\log \frac{p(x, y)}{p(x)}
\end{split}
\label{eq:4}
\end{align}
$$

### 2. 最大熵原理简介 ###

最大熵原理可以表述为，在满足$k+1$个约束条件的模型集合中，选取熵$H(p)$最大的模型。约束条件为：

$$
\begin{align}
\begin{split}
\sum_{x}p(x) &= 1 \\
\sum_{x}p(x)f_1(x) &= \tau_1 \\
\vdots \\
\sum_{x}p(x)f_k(x) &= \tau_k
\end{split}
\label{eq:5}
\end{align}
$$

使用拉格朗日乘子法求解带上述有约束的极值，即：

$$
\begin{align*}
\begin{split}
L(p) = -\sum_{x}p(x)\log{p(x)} &+ \\
\lambda_0(\sum_{x}p(x) - 1) &+ \\
\lambda_1(\sum_{x}p(x)f_1(x) - \tau_1) &+ \\
\cdots &+\\
\lambda_k(\sum_{x}p(x)f_1(x) - \tau_k)
\end{split}
\end{align*}
$$

$L(p)$对每一个$p(x)$偏导数$\frac{\partial L(p)}{\partial p(x)}$为0，即：

$$
-\log{p(x)} - 1 + \lambda_0 + \sum_{j=1}^{k}\lambda_j f_j(x) = 0
$$

解得

$$
p(x) = \frac{\exp\left(\sum\limits_{j=1}^{k}\lambda_j f_j(x)\right)}{\exp(1 - \lambda_0)}
$$

由约束条件$\sum\limits_x p(x)=1$得：

$$
\begin{align}
p(x) = \frac{1}{Z}\exp\left(\sum\limits_{j=1}^{k}\lambda_j f_j(x)\right)
\label{eq:6}
\end{align}
$$

其中

$$
\begin{align}
Z = \sum\limits_x \exp\left(\sum\limits_{j=1}^{k}\lambda_j f_j(x)\right)
\label{eq:7}
\end{align}
$$

将$\eqref{eq:3}$带入约束条件$\eqref{eq:2}$中，即可解得$\lambda_j$。

### 参考资料 ###

* [“熵”不起：从熵、最大熵原理到最大熵模型](http://spaces.ac.cn/archives/3552/)

### 更新记录 ###

2017年7月15日

