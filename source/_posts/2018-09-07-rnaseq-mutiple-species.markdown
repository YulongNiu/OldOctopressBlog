---
layout: post
title: "两物种混合RNA-Seq的基因表达量模型"
date: 2018-09-07 23:29:44 +0800
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
\newcommand{\tildelu}{\widetilde{l_u}}
$$

$$
\newcommand{\tildesv}{\widetilde{l_v}}
$$

$$
\newcommand{\P}{\mathrm{P}}
$$

## 1. 符号表示 ##

物种$A$有$U$个长度为$l_u$的转录序列，构成的转录本集合为$R=\\{r_1, r_2, \dots, r_U\\}$，每个转录序列的相对丰度为$\rho_u$，$\sum\limits_{u=1}^U \rho_u= 1$。物种$B$有$V$个长度为$l_v$的转录序列$s_v$，构成转录本集合$S=\\{s_1, s_2, \dots, s_V\\}$，每个转录序列的相对丰度为$\theta_v$, $\sum\limits_{v=1}^V \theta_v= 1$。两物种的所有转录序列集合构成$T=\\{t_1, t_2, \dots, t_{U+V}\\}$。

<!--more-->

对于某个转录片段$f$，来自于转录序列$r_u$（有效长度为$\tildelu$）的概率为：

$$
\begin{align}
\begin{split}
\P(f \in r_u) &= \frac{\rho_u \tildelu}{\sum\limits_{r \in R} \rho_r \widetilde{l_r}} \\
&= \alpha_u
\end{split}
\label{eq:1}
\end{align}
$$

转录片段$f$，来自转录序列$t_u$的某个位置的概率为：

$$
\begin{align}
\begin{split}
\P(\mathrm{pos}, f \in r_u) = \frac{\alpha_u}{\tildelu}
\end{split}
\label{eq:2}
\end{align}
$$


对于某个转录片段$f$，来自于转录序列$s_v$（有效长度为$\tildesv$）的概率为：

$$
\begin{align}
\begin{split}
\P(f \in s_v) &= \frac{\theta_v \tildesv}{\sum\limits_{s \in S} \theta_s \widetilde{l_s}} \\
&= \beta_v
\end{split}
\label{eq:3}
\end{align}
$$

转录片段$f$，来自转录序列$t_v$的某个位置的概率为：

$$
\begin{align}
\begin{split}
\P(\mathrm{pos}, f \in s_v) = \frac{\beta_v}{\tildesv}
\end{split}
\label{eq:4}
\end{align}
$$


$\eqref{eq:1}$和从$\eqref{eq:2}$也提供了$\rho_u$和$\alpha_u$，$\theta_v$和$\beta_v$转换。单个转录组中，物种$A$和物种$B$的转录片段混合在一起，构成集合$F=\\{f_1, f_2, \dots, f_N\\}$，总转录片段数目为$N=\|F\|$。

## 2. 混合物种模型 ##

所有转录片段的长度为定值，事实上某一个转录片段只能来自于单一物种的单一转录序列。实际上，可观测为：某一个转录片段对应至多个可能物种的多个可能转录序列。不可观测数据为：该转录片段真实来源于哪个物种$Y = \\{y_A, y_B\\}$的哪个转录序列$T=\\{t_1, t_2, \dots, t_{U+V}\\}$。混合物种模型要估计的参数为$\eta=\\{\alpha_1, \dots, \alpha_u, \beta_1, \dots, \beta_v\\}$。

对于转录片段$f$，实际观测数据为$0$和$1$构成的长度为$U+V$的指示向量$M$。元素$m_i=1$表示$f$对应至转录序列$r_i$，反之$0$表示$f$没有比对至$r_i$。同时，还可以构造另一个$0$和$1$构成长度为$2$的指示向量$N=(n_1, n_2)$表示$f$是否对应至某一物种。例如：$(1, 1)$表示$f$对应至物种$A$或$B$，$(0, 1)$表示$f$只对应到物种$B$。

$$
\begin{align}
\begin{split}
\P(\mathrm{pos}|\eta)  &= \P(\mathrm{pos}, f \in A|\eta) + \P(\mathrm{pos}, f \in B|\eta)\\
&= \P(\mathrm{pos} | f \in A, \eta) \P(f \in A) + \P(\mathrm{pos} | f \in B, \eta) \P(f \in B) \\
&= \sum_{i=1}^{U} \P(\mathrm{pos}, f \in r_i|\eta) \P(f \in A) + \sum_{j=U+1}^{U+V} \P(\mathrm{pos}, f \in r_j|\eta) \P(f \in B) \\
&= \sum_{i=1}^{U} n_1 m_i \frac{\alpha_i}{\widetilde{l_i}} \P(f \in A) + \sum_{j=U+1}^{U+V} n_2 m_j \frac{\beta_j}{\widetilde{l_j}} \P(f \in B) \\
\end{split}
\label{eq:5}
\end{align}
$$

其中，$\P(f \in A) + \P(f \in B) = 1$。由于$n_1 m_i = m_i$，即如果观察到转录片段$f$对应至物种$A$，则一定可以观察到$f$对应至物种$A$中的某个转录片段$\\{r_1, \dots, r_U\\}$。同理，$n_2 m_j = m_j$。

$\eqref{eq:5}$可以写为：

$$
\begin{align}
\begin{split}
\P(\mathrm{pos}, f \in r|\eta) &= \sum_{i=1}^{U}m_i \frac{\alpha_i}{\widetilde{l_i}} \P(f \in A) + \sum_{j=U+1}^{U+V}m_j \frac{\beta_j}{\widetilde{l_j}} \P(f \in B) \\
\end{split}
\label{eq:6}
\end{align}
$$

使用EM算法，对于第$n$次迭代：

$$
\begin{align}
\begin{split}
\lambda_i^{(n+1)} &= \frac{\alpha_j^{(n)} \frac{m_i}{\widetilde{l_j}} \P(f \in A)}{\sum\limits_{i=1}^{U}\alpha_i^{(n)} \frac{m_i}{\widetilde{l_i}} \P(f \in A) + \sum\limits_{j=U+1}^{U+V}\beta_j^{(n)} \frac{m_j}{\widetilde{l_j}} \P(f \in B)}\\
\lambda_j^{(n+1)} &= \frac{\beta_j^{(n)} \frac{m_i}{\widetilde{l_j}} \P(f \in A)}{\sum\limits_{i=1}^{U}\alpha_i^{(n)} \frac{m_i}{\widetilde{l_i}} \P(f \in A) + \sum\limits_{j=U+1}^{U+V}\beta_j^{(n)} \frac{m_j}{\widetilde{l_j}} \P(f \in B)}\\
\alpha_i^{(n+1)} &= \frac{\sum\limits_{f \in F} \lambda_i^{(n+1)}}{\sum\limits_{f \in F} \sum\limits_{i=1}^{U} \lambda_i^{(n+1)}} \\
\beta_j^{(n+1)} &= \frac{\sum\limits_{f \in F} \lambda_j^{(n+1)}}{\sum\limits_{f \in F} \sum\limits_{j=U+1}^{U+V} \lambda_j^{(n+1)}} \\
\P(f \in A)^{(n+1)} &= \frac{\sum\limits_{f \in F} \sum\limits_{i=1}^{U} \lambda_i^{(n+1)}}{\sum\limits_{f \in F} \sum\limits_{i=1}^{U} \lambda_i^{(n+1)} + \sum\limits_{f \in F} \sum\limits_{j=U+1}^{U+V} \lambda_j^{(n+1)}} = \frac{\sum\limits_{f \in F} \sum\limits_{i=1}^{U} \lambda_i^{(n+1)}}{N}\\
\P(f \in B)^{(n+1)} &= \frac{\sum\limits_{f \in F} \sum\limits_{j=U+1}^{U+V} \lambda_j^{(n+1)}}{\sum\limits_{f \in F} \sum\limits_{i=1}^{U} \lambda_i^{(n+1)} + \sum\limits_{f \in F} \sum\limits_{j=U+1}^{U+V} \lambda_j^{(n+1)}} = \frac{\sum\limits_{f \in F} \sum\limits_{j=U+1}^{U+V} \lambda_j^{(n+1)}}{N} 
\end{split}
\label{eq:7}
\end{align}
$$


### 更新记录 ###

2018年9月13日
