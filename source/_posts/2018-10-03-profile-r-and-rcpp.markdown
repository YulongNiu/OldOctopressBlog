---
layout: post
title: "R和Rcpp的性能检测"
date: 2018-10-03 19:38:37 +0800
comments: true
categories: r
---

## 1. R性能检测 ##

直接使用<span style="color: blue">profvis</span>包即可，例如[示例](https://rstudio.github.io/profvis/)。

## 2. Rcpp性能检测 ##

### 2.1 安装依赖软件 ###

{% codeblock lang:bash Pre-requested tools %}
$ sudo dnf install gperftools-devel google-perftools graphviz ghostscript kcachegrind
{% endcodeblock %} 

<!--more-->

### 2.2 编译 ###

在包（包名称为`Mypkg`）目录`src`建立如下文件：

{% codeblock lang:c++ Profiler Rcpp code %}
#include <Rcpp.h>
#include <gperftools/profiler.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP start_profiler(SEXP str) {
  ProfilerStart(as<const char*>(str));
  return R_NilValue;
}

// [[Rcpp::export]]
SEXP stop_profiler() {
  ProfilerStop();
  return R_NilValue;
}
{% endcodeblock %} 

包目录`src`的Makevars文件中添加`-lprofile`选项，例如`PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -lprofiler`。之后，安装包，并重新载入。

### 2.3 调试 ###

使用方法为：

{% codeblock lang:r Profiling %}
Mypkg:::start_profiler("/tmp/profile.out")
run_cpp_codes()
Mypkg:::stop_profiler()
{% endcodeblock %} 

查看profile结果：

{% codeblock lang:bash Profiling results %}
## text
pprof --text src/Mypkg.so /tmp/profile.out

## pdf
pprof --pdf src/Mypkg.so /tmp/profile.out > profile.pdf

## kcachegrind
pprof --callgrind src/RNASeqEM.so R/profile.out > profile.res
{% endcodeblock %} 



### <a id="Ref">参考资料</a> ###

*  [Introduction to High-Performance Computing with R by Dr. Dirk Eddelbuettel](https://arxiv.org/abs/1104.3889)

* [Profiling Rcpp packages](https://minimallysufficient.github.io/r%20programming%20cpp/2018/02/16/profiling-rcpp-packages.html)

* [用gperftools对C/C++程序进行profile](https://airekans.github.io/cpp/2014/07/04/gperftools-profile)

### 更新记录 ###

2018年10月4日
