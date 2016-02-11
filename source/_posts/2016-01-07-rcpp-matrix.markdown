---
layout: post
title: "Rcpp操作矩阵集锦"
date: 2016-01-07 20:50:41 +0800
comments: true
published: false
categories: R
---

使用<span style="color: blue">Rcpp</span>或者<span style="color: blue">RcppArmadillo</span>操作矩阵。


## 1. <span style="color: blue">Rcpp</span> ##

## 2. <span style="color: blue">RcppArmadillo</span> ##

基本类型是`mat`。

Cpp矩阵转时，可以避免拷贝矩阵，以提升效率，比如：`mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false) `{:.language-cpp}。

### <a id="Ref">参考网址</a> ###

* [Armadillo矩阵](http://arma.sourceforge.net/docs.html#adv_constructors_mat)

* [StackOverflow中矩阵提取](https://stackoverflow.com/questions/13038256/subset-of-a-rcpp-matrix-that-matches-a-logical-statement) 




