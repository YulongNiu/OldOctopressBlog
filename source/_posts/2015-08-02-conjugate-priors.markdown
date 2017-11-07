---
layout: post
title: "一些共轭先验"
date: 2015-08-02 15:53:54 +0800
comments: true
categories: bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

$$
\newcommand{\md}{\mathrm{d}}
$$

[共轭先验（Conjugate prior）](https://en.wikipedia.org/wiki/Conjugate_prior#cite_note-beta_rate-6)在贝叶斯估计中被广泛应用，本文尝试详细推理一些常见分布的共轭先验证。

<!--more-->

贝叶斯公式：

$$
\begin{align}
\begin{split}
f(\theta | x) &= \frac{f(\theta, x)}{f(x)} \\
&= \frac{f(x | \theta)f(\theta)}{\int f(x | \theta)f(\theta) \md \theta}
\end{split}
\label{eq:1}
\end{align}
$$


## 1. 离散分布 ##

### 1.1 伯努利分布 ###

伯努利分布（Bernoulli distribution）的概率质量函数为：

$$
\begin{align}
\begin{split}
f(k;p) = p^k (1-p)^{1-k} \quad \mathrm{for} \quad k \in \{0, 1\}
\end{split}
\label{eq:2}
\end{align}
$$

对于随机变量$X_i \in \\{X_1, X_2, \dots, X_m\\}$易得，$p$的极大似然估计（Maximum Likelihood Estimator, MLE）为$\hat{p}=\frac{\sum_{i=1}^{m}k_i}{m}$。

该分布的共轭先验为Beta分布$\mathrm{Beta}(\alpha, \beta)$，即对于随机变量$X_i$：

$$
\begin{align}
\begin{split}
f(p|X_i) &= \frac{p^{k_i} (1-p)^{1-k_i} \frac{1}{\mathrm{B}(\alpha, \beta)} p^{\alpha - 1} (1-p)^{\beta -1}}{f(X_i)} \\
&=\frac{\frac{1}{\mathrm{B}(\alpha, \beta)} p^{k_i+\alpha-1}(1-p)^{\beta - k_i}}{\int_0^1 \frac{1}{\mathrm{B}(\alpha, \beta)} p^{k_i+\alpha-1}(1-p)^{\beta - k_i} \md p} \\
&= \frac{p^{k_i+\alpha-1}(1-p)^{\beta - k_i}}{\mathrm{B}(k_i+\alpha, \beta+1 -k_i)} \\
&= \mathrm{Beta}(k_i+\alpha, \beta+1-k_i)
\end{split}
\label{eq:3}
\end{align}
$$

根据$\eqref{eq:3}$易得，$f(p \vert X_1, X_2, \dots, X_m) = \mathrm{Beta}(\sum_{i=1}^{m}k_i+\alpha, \beta+m-\sum_{i=1}^{m}k_i)$，期望$\hat{p}=\frac{\sum_{i=1}^{m}k_i+\alpha}{m+\alpha+\beta}$。特别，当$\alpha=1$和$\beta=1$时，即共轭先验为$0-1$之间的均匀分布，$\hat{p}=\frac{\sum_{i=1}^{m}k_i+1}{m+2}$。

### 1.2 二项分布 ###

二项分布（binomial distribution）的概率质量函数为：

$$
\begin{align}
\begin{split}
f(k;n,p) = \binom{n}{k} p^k (1-p)^{n-k}
\end{split}
\label{eq:4}
\end{align}
$$

对于随机变量$X_i \in \\{X_1, X_2, \dots, X_m\\}$易得，MLE为$\hat{p}=\frac{\sum_{i=1}^{m}k_i}{nm}$。

该分布的共轭先验为Beta分布$\mathrm{Beta}(\alpha, \beta)$。因此，类似于$\eqref{eq:3}$，$f(p \vert X_1, X_2, \dots, X_m) = \mathrm{Beta}(\sum_{i=1}^{m}k_i+\alpha, \beta+nm-\sum_{i=1}^{m}k_i)$，期望$\hat{p}=\frac{\sum_{i=1}^{m}k_i+\alpha}{nm+\alpha+\beta}$。

### 1.3 负二项分布 ###

负二项（negative binomial distribution）的概率质量函数为：

$$
\begin{align}
\begin{split}
f(k;r,p) = \binom{r+k-1}{k} p^k (1-p)^{r}
\end{split}
\label{eq:5}
\end{align}
$$

对于随机变量$X_i \in \\{X_1, X_2, \dots, X_m\\}$易得，MLE为$\hat{p}=\frac{\sum_{i=1}^{m}k_i}{mr+\sum_{i=1}^{m}k_i}$。

该分布的共轭先验为Beta分布$\mathrm{Beta}(\alpha, \beta)$。因此，类似于$\eqref{eq:3}$，$f(p \vert X_1, X_2, \dots, X_m) = \mathrm{Beta}(\sum_{i=1}^{m}k_i+\alpha, \beta+mr)$，期望$\hat{p}=\frac{\sum_{i=1}^{m}k_i+\alpha}{mr+\sum_{i=1}^{m}k_i+\alpha+\beta}$。

### 1.4 多项分布 ###

多项二项（multinomial distribution）的概率质量函数为：

$$
\begin{align}
\begin{split}
f(x_1, x_2,\dots,x_k;n,p_1,p_2,\dots,p_k) = \frac{n!}{x_1!x_2!\dots x_k!}p_1^{x_1}p_2^{x_2}\dots p_k^{x_k} \quad \mathrm{for} \quad \sum_{i=1}^k p_i=1 \quad \sum_{i=1}^k x_i=n
\end{split}
\label{eq:6}
\end{align}
$$

对于随机变量$X_j \in \\{X_1, X_2, \dots, X_m\\}$易得，MLE为$\hat{p_i}=\frac{\sum_{j=1}^{m}x_{ij}}{\sum_{i=1}^{k}\sum_{j=1}^{m}x_{ij}}$。

该分布的共轭先验为Dirichlet分布$\mathrm{Dir}(\alpha_1,\alpha_2,\dots,\alpha_k)$，即对于随机变量$X_i$：

$$
\begin{align}
\begin{split}
f(p_1,p_2,\dots,p_k|X_1) &= \frac{\frac{n!}{x_{11}!x_{21}!\dots x_{k1}!}p_1^{x_{11}}p_2^{x_{21}}\dots p_k^{x_{k1}} \frac{1}{\mathrm{B}(\boldsymbol{\alpha})} p_1^{\alpha_1-1} p_2^{\alpha_2-1} \dots p_k^{\alpha_k-1}}{f(X_i)} \\
&= \frac{p_1^{\alpha_1+x_{11}-1} p_2^{\alpha_2+x_{21}-1} \dots p_k^{\alpha_k+x_{k1}-1}}{\int_0^1 \int_0^{1-p_1} \cdots \int_0^{1-\sum_{i=1}^{k-2} p_i}p_1^{\alpha_1+x_{11}-1} p_2^{\alpha_2+x_{21}-1} \dots (1-\sum_{i=1}^{k-1}p_i)^{\alpha_k+x_{k1}-1} \md p_1 \md p_2 \dots \md p_{k-1}} \\
&=\frac{\mathrm{\Gamma}\left(\sum_{i=1}^{k}(\alpha_i+x_{i1})\right)}{\prod_{i=1}^{k}\mathrm{\Gamma}(\alpha_i+x_{i1})}p_1^{\alpha_1+x_{11}-1} p_2^{\alpha_2+x_{21}-1} \dots p_k^{\alpha_k+x_{k1}-1} \\
&=\mathrm{Dir}(\alpha_1+x_{11},\alpha_2+x_{21},\dots,\alpha_k+x_{k1})
\end{split}
\label{eq:7}
\end{align}
$$

根据$\eqref{eq:7}$易得，$f(p_1,p_2,\dots,p_k \vert X_1, X_2, \dots, X_m) = \mathrm{Dir}(\alpha_1+\sum_{j=1}^{m}x_{1j},\alpha_2+\sum_{j=1}^{m}x_{2j},\dots,\alpha_k+\sum_{j=1}^{m}x_{kj})$，期望为$\hat{p_i}=\frac{\sum_{j=1}^{m}x_{ij}+\alpha_i}{\sum_{i=1}^{k}\sum_{j=1}^{m}x_{ij} + \sum_{i=1}^{k}\alpha_i}$。特别，当$\alpha_1=\alpha_2=\dots=\alpha_k=1$时，$\hat{p_i}=\frac{\sum_{j=1}^{m}x_{ij}+1}{\sum_{i=1}^{k}\sum_{j=1}^{m}x_{ij}+k}$。

### 1.5 泊松分布 ###

泊松分布（Poisson distribution）的概率质量函数为：

$$
\begin{align}
\begin{split}
f(k;\lambda) = \frac{\lambda^k \mathrm{e}^{-\lambda}}{k!}
\end{split}
\label{eq:8}
\end{align}
$$

对于随机变量$X_j \in \\{X_1, X_2, \dots, X_m\\}$易得，MLE为$\hat{\lambda}=\frac{\sum_{i=1}^m k_i}{m}$。

该分布的共轭先验为Gamma分布$\mathrm{Gamma}(\alpha, \beta)$，即对于随机变量$X_i$：

$$
\begin{align}
\begin{split}
f(\lambda|X_i) &= \frac{\frac{1}{k_i!}\lambda^{k_i} \mathrm{e}^{-\lambda}\frac{1}{\mathrm{\Gamma}(\alpha)}\beta^\alpha \lambda^{\alpha-1}\mathrm{e}^{-\beta\lambda}}{f(X_i)} \\
&= \frac{\lambda^{k_i+\alpha-1}\mathrm{e}^{-(\beta+1)\lambda}}{\int_0^\infty \lambda^{k_i+\alpha-1}\mathrm{e}^{-(\beta+1)\lambda} \md \lambda} \\
&= \frac{(\beta+1)^{k_i+\alpha} \lambda^{k_i+\alpha-1}\mathrm{e}^{-(\beta+1)\lambda}}{\mathrm{\Gamma}(k_i+\alpha)} \\
&= \mathrm{Gamma}(k_i+\alpha, \beta+1)
\end{split}
\label{eq:9}
\end{align}
$$

根据$\eqref{eq:9}$易得，$f(\lambda \vert X_1, X_2, \dots, X_m) = \mathrm{Gamma}(\sum_{i=1}^{m}k_i+\alpha, \beta+m)$，期望$\hat{p}=\frac{\sum_{i=1}^{m}k_i+\alpha}{\beta+m}$。

## 2. 连续分布 ##

### 2.1 指数分布 ###

指数分布（exponential distribution）的概率密度函数为：

$$
\begin{align}
\begin{split}
f(x;\lambda) = \lambda \mathrm{e}^{-\lambda x} \quad \mathrm{for} \quad x\ge0
\end{split}
\label{eq:10}
\end{align}
$$

对于随机变量$X_j \in \\{X_1, X_2, \dots, X_m\\}$易得，MLE为$\hat{\lambda}=\frac{m}{\sum_{i=1}^m k_i}$。

该分布的共轭先验为Gamma分布$\mathrm{Gamma}(\alpha, \beta)$。因此，类似于$\eqref{eq:9}$，$f(\lambda \vert X_1, X_2, \dots, X_m) = \mathrm{Gamma}(\alpha+m, \beta+\sum_{i=1}^{m}k_i)$，期望$\hat{p}=\frac{\alpha+m}{\beta+\sum_{i=1}^{m}k_i}$

### 2.2 已知方差的正态分布 ###

正态分布的概率密度函数为：

$$
\begin{align}
\begin{split}
f(x; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi \sigma^2}} \mathrm{e}^{-\frac{(x-\mu)^2}{2\sigma^2}}
\end{split}
\label{eq:11}
\end{align}
$$

对于随机变量$X_j \in \\{X_1, X_2, \dots, X_m\\}$易得，MLE为$\hat{\sigma}^2=\frac{\sum_{i=1}^m(x_i-\mu)^2}{m}$。



## 3. 一些积分证明 ##

### 3.1 多元beta函数 ###

$\int_0^1 \int_0^{1-p_1} \cdots \int_0^{1-\sum_{i=1}^{k-2} p_i}p_1^{\alpha_1-1} p_2^{\alpha_2-1} \dots (1-\sum_{i=1}^{k-1}p_i)^{\alpha_k-1} \md p_1 \md p_2 \dots \md p_{k-1} = \frac{\mathrm{\Gamma}\left(\sum_{i=1}^{k}(\alpha_i)\right)}{\prod_{i=1}^{k}\mathrm{\Gamma}(\alpha_i)} \quad \mathrm{for} \quad \sum_{i=1}^k p_i=1$

令$p_{k-1} = (1-\sum_{i=1}^{k-2}p_i)u$，考察积分：

$$
\begin{align*}
\begin{split}
\int_0^{1-\sum_{i=1}^{k-2} p_i} p_{k-1}^{\alpha_{k-1}-1} (1-\sum_{i=1}^{k-1}p_i)^{\alpha_k-1} \md p_{k-1} &= \left(1-\sum_{i=1}^{k-2}p_i \right)^{\alpha_{k-1}+\alpha_k-2} \int_0^1 u^{\alpha_{k-1}-1}(1-u)^{\alpha_k-1} \md u \\
&= \left(1-\sum_{i=1}^{k-2}p_i \right)^{\alpha_{k-1}+\alpha_k-2} \frac{\mathrm{\Gamma}(\alpha_{k-1})\mathrm{\Gamma}(\alpha_{k})}{\mathrm{\Gamma}(\alpha_{k-1} + \alpha_{k})}
\end{split}
\end{align*}
$$

迭代相乘后即得。

### 3.2 Gamma积分 ###

$\int_0^\infty \frac{x^{\alpha-1}}{\mathrm{e}^{\lambda x}} \md x=\frac{\mathrm{\Gamma}(\alpha)}{\lambda^\alpha}$

令$\lambda x=u$：

$$
\begin{align*}
\begin{split}
\int_0^\infty \frac{x^{\alpha-1}}{\mathrm{e}^{\lambda x}} \md x &= \lambda^{-\alpha} \int_0^\infty \frac{u^{\alpha-1}}{\mathrm{e}^u} \md u \\
&= \frac{\mathrm{\Gamma}(\alpha)}{\lambda^\alpha}
\end{split}
\end{align*}
$$




### 参考资料 ###

* [Chapter 9 The exponential family: Conjugate priors](https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/other-readings/chapter9.pdf)



### 更新记录 ###

2017年7月15日
