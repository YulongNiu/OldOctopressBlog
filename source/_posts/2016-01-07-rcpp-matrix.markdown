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

<!--more-->

## 2. <span style="color: blue">RcppArmadillo</span> ##

基本类型是`mat`、`vec`（`colvec`）和`rowvec`。

* 属性
  
    * 对于矩阵，行数：`m.n_rows;`；列数：`m.n_cols;`；维度：`m.size();`或`size(m);`。对于向量，元素数：`v.n_elem;`。

* 特殊向量或矩阵

   * 全是0`ones<mat>(3, 4);`/`vec(10, fill::ones);`/；全是1`zeros<vec>(10);`/`mat(3, 4, fill::zeros);`；全是某个数`mat a(4, 5); a.fill(123.4);`。
   
   * 连续向量，规定长度`linspace<vec>(0, 5, 6);`；连续向量，规定间距`regspace<vec>(0, 2, 9);`。

* 取值

   * 对于向量，连续取值：`v.subvec(1stIdx, lastIdx);`；非连续，可以考虑使用`find()`函数，比如：`v.elem(find(v > 0));`。
   
   * 对于矩阵，连续取值：`m.col(Idx);`/`m.row(Idx);`/`m.cols(Idx);`/`m.rows(Idx);`/`m.submat(1stRowIdx, lastRowIdx, 1stColIdx, lastColIdx);`；非连续，`m.submat(vecRowIdx, vecColIdx);`。

* Rcpp矩阵转换为RcppArmadillo矩阵，可以避免拷贝矩阵，以提升效率，`mat(ptr_aux_mem, n_rows, n_cols, copy_aux_mem = true, strict = false)`{:.language-cpp}。同样道理，可以转化向量。例如：

{% codeblock lang:cpp transfer matrix and vector %}
arma::mat TransferMatArma(Rcpp::NumericMatrix x, Rcpp::NumericVector y) {
    mat tx(x.begin(), x.nrow(), x.ncol(), false);
    vec ty(y.begin(), y.size(), false);
    return tx;
}

Rcpp::NumericVector TransferMatRcpp(arma::mat x, arma::vec y) {
    NumericMatrix tx(x.n_rows, x.n_cols, x.begin());
    NumericVector ty(y.begin(), y.end());
    return ty;
    
// 不要使用as<IntegerVector>(wrap(y))，会有内存泄露。
}
{% endcodeblock %}

* 使用`.each_col()`/`.each_row()`/`.for_each()`替代`apply()`

{% codeblock lang:cpp replace apply() %}
arma::mat TestMat(arma::mat M, double a) {

  M.for_each([a](mat::elem_type& val) {
      val = val > 0 ? val : a;
    });

  M.each_row([](rowvec& r) {
      r /= r.max();
    });

  return M;
}
{% endcodeblock %}

* 使用`sum(M, 0);`和`sum(M, 1);`分别替代`colSums(M)`和`rowSums(M)`。








### <a id="Ref">参考网址</a> ###

* [Armadillo矩阵](http://arma.sourceforge.net/docs.html#adv_constructors_mat)

* [StackOverflow中矩阵提取](https://stackoverflow.com/questions/13038256/subset-of-a-rcpp-matrix-that-matches-a-logical-statement)

* [RcppArmadillo和R对照表](https://github.com/petewerner/misc/wiki/RcppArmadillo-cheatsheet)

* [Rcpp Quick Reference Guide](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-quickref.pdf) 

* [Rcpp note](http://statr.me/rcpp-note/) 



### 更新记录 ###

2017年1月15日




