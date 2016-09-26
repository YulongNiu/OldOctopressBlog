---
layout: post
title: "Rcpp操作矩阵和向量集锦"
date: 2016-01-07 20:50:41 +0800
comments: true
published: true
categories: R
---

收集和记录<span style="color: blue">Rcpp</span>或者<span style="color: blue">RcppArmadillo</span>操作矩阵和向量。


## 1. <span style="color: blue">Rcpp</span> ##

* 可以使用逻辑下标（`LogicalVector`）对向量和列表[取值](http://gallery.rcpp.org/articles/subsetting/)。

## 2. <span style="color: blue">RcppArmadillo</span> ##

基本类型是`mat`、`vec`（`colvec`）和`rowvec`。

* Rcpp矩阵转换为RcppArmadillo矩阵，可以避免拷贝矩阵，以提升效率，`mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false) `{:.language-cpp}。同样道理，可以转化向量。例如：

{% codeblock lang:cpp transfer matrix and vector %}
arma::mat TransferMatArma(Rcpp::NumericMatrix x, Rcpp::NumericVector y) {
    mat tx(x.begin, x.nrow(), x.ncol(), false);
    mat ty(y.begin(), y.size(), false);
    return tx;
}

Rcpp::NumericVector TransferMatRcpp(arma::mat x, arma::vec y) {
    NumericVector ty(y.begin(), y.end());
    return ty;
    
// 不要使用as<IntegerVector>(wrap(y))，会有内存泄露。
}
{% endcodeblock %}




### <a id="Ref">参考网址</a> ###

* [Armadillo矩阵](http://arma.sourceforge.net/docs.html#adv_constructors_mat)

* [StackOverflow中矩阵提取](https://stackoverflow.com/questions/13038256/subset-of-a-rcpp-matrix-that-matches-a-logical-statement)

* [RcppArmadillo和R对照表](https://github.com/petewerner/misc/wiki/RcppArmadillo-cheatsheet)

* [Rcpp Quick Reference Guide](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-quickref.pdf) 



### 更新记录 ###

2016年9月26日




