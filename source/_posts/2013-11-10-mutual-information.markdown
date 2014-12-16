---
layout: post
title: "A brief intruction of mutual information and demonstration with R"
date: 2013-11-10 01:21:16 -0500
comments: true
styles: [data-table]
categories: bioinfor
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

$\newcommand{\entropfrac}[2]{\frac{#1}{#2} \log \left( \frac{#1}{#2} \right)}$


## Mututal Information (MI) ##

### Introduction ###

[Mutual Information (MI)][MI wiki] distance is used to measure the distance between two genes vectors, for example $x_1 = \{1, 0, 1, 1, 1, 1, 0\}$ and $y_1 = \{0, 1, 1, 1, 1, 1, 0\}$. It is easily to transfer the two vectors into a binary table:

------------------

|---------------+---------------+--------------+---------|
|**X/Y**        |**1(Presence)**|**0(Absence)**|**Sum**  |
|:-------------:|:-------------:|:------------:|:-------:|
|**1(Presence)**|a              |b             |a+b      |
|---------------|---------------|--------------|---------|
|**0(Absence)** |c              |d             |c+d      |
|---------------|---------------|--------------|---------|
|**Sum**        |a+c            |b+d           |n=a+b+c+d|
|---------------|---------------|--------------|---------|

----------------------

<!--more-->

Typically, here we give the example of two discrete variables, the mutual information between $x_1$ and $y_1$ is

$$
\begin{align}
\begin{split}
I(X;Y) &= H(X) + H(Y) - H(X,Y)\\
&= -\sum_{x \in \{0, 1\}}p(x)\log(p(x)) - \sum_{y \in \{0, 1\}}p(y)\log(p(y))\\
& \quad -\left( -\sum_{x \in \{0, 1\}}\sum_{y \in \{0, 1\}}p(x,y)\log(p(x,y)) \right)\\
\end{split}
\label{eq:1}
\end{align}
$$

The $\eqref{eq:1}$ is equal to

$$
\begin{equation}
I(X;Y) = \sum_{x \in \{1, 0\}}\sum_{y \in \{1, 0\}} p(x,y) \frac{\log(p(x,y))}{p(x)p(y)}
\label{eq:2}
\end{equation}
$$

$p(x)$ is the probability that a symbol (here is 0 or 1) appears in the gene vector X regardless that what the symbol is in gene vector Y. $p(y)$ has a similar definition of $p(x). $$p(x, y)$ is probability of a symbol combination appears in gene vector X and Y. In this example, there are four kinds of symbol combination $(1, 1)$, $(1, 0)$, $(0, 1)$ and $(1, 1)$.

If we use the binary table to illustrate this equation, the $\eqref{eq:1}$ is:

$$
\begin{align}
\begin{split}
I(X; Y) &= -\left( \entropfrac{a+c}{n} + \entropfrac{b+d}{n} \right)\\
& \quad -\left( \entropfrac{a+b}{n} + \entropfrac{c+d}{n} \right)\\
& \quad - \left( - \left(
\entropfrac{a}{n} + \entropfrac{b}{n} + \entropfrac{c}{n} + \entropfrac{d}{n}
\right) \right)\\
\end{split}
\label{eq:3}
\end{align}
$$

The $\eqref{eq:3}$ is mathmatically equal to:

$$
\begin{align}
\begin{split}
I(X; Y) &= \frac{a}{n}\log\frac{na}{(a+b)(a+c)} + \frac{b}{n}\log\frac{nc}{(a+d)(b+d)}\\
& \quad \frac{c}{n}\log\frac{nc}{(a+b)(a+c)} + \frac{d}{n}\log\frac{nd}{(d+c)(d+b)}
\end{split}
\label{eq:4}
\end{align}
$$

### Example ###

We can use R to directly calculate the MI between two gene vectors mentioned above.

1. Use basic R function

{% codeblock lang:r %}
x1 <- c(1, 0, 1, 1, 1, 1, 0)
y1 <- c(0, 1, 1, 1, 1, 1, 0)
table(x1, y1)
   y1
x1  0 1
  0 1 1
  1 1 4
# calculate MI
4/7 * log(28/25) + 1/7 * log(7/10) + 1/7 * log(7/10) + 1/7 * log(7/4)
[1] 0.04279723
{% endcodeblock %}

2. Use R package <span style="color: blue">bioDist</span>

{% codeblock lang:r %}
library(bioDist)
mutualInfo(rbind(x1, y1))
           x1
y1 0.04279723
{% endcodeblock %}

### Reference ###

* [An example of MI](http://nlp.stanford.edu/IR-book/html/htmledition/mutual-information-1.html)

* [Slices of MI](http://www1.ece.uic.edu/~devroye/courses/ECE534/lectures/ch2.pdf)

* [Entropy and Mutual Information](http://people.cs.umass.edu/~elm/Teaching/Docs/mutInf.pdf)

* [Useful Websites 1](http://www.scholarpedia.org/article/Mutual_information)

* [Slices about MI application in bioinformatics](http://montana.informatics.indiana.edu/LabWebPage/Presentations/Vikas_Nov02_2011.pdf)

* Seung-Seok Choi, Sung-Hyuk Cha, Charles C. Tappert: [A Survey of Binary Similarity and Distance Measures](http://www.iiisci.org/journal/sci/Abstract.asp?var=&id=GS315JG).

* Huynen M, Snel B, Lathe W 3rd, Bork P: Predicting protein function by genomic context: quantitative evaluation and qualitative inferences.Genome Res. 2000;10(8):1204-10.

* Korber BT, Farber RM, Wolpert DH, Lapedes AS: Covariation of mutations in the V3 loop of human immunodeficiency virus type 1 envelope protein: an information theoretic analysis. Proc Natl Acad Sci U S A. 1993;90(15):7176-80.

* Kensche PR, van Noort V, Dutilh BE, Huynen MA: Practical and theoretical advances in predicting the function of a protein by its phylogenetic distribution. J R Soc Interface. 2008;5(19):151-70.

[MI wiki]: http://en.wikipedia.org/wiki/Mutual_information

### Update record ###

12/16/2014
