---
layout: post
title: "使用RcppParallel并行计算"
date: 2018-07-19 22:19:24 +0800
comments: true
categories: R
---

在之前的[博文](http://yulongniu.bionutshell.org/blog/2014/06/25/parallel-package/)中，我详细讨论了使用多种R包实现并行计算。其中，提到一个非常重要的问题：

    当循环数很大时（1万以上），`foreach`会变得非常慢。
    
这个问题在Florian Privé的[A guide to parallelism in R](https://privefl.github.io/blog/a-guide-to-parallelism-in-r/)中也提到，解释是`foreach`每次只合并100个循环结果。

<!--more-->

这篇博文里，我尝试使用<span style="color: blue">RcppParallel</span>包调用`C++`的并行方法。结论是：**在循环数很大时，<span style="color: blue">RcppParallel</span>包提供的并行方法优于`foreach`**。

一个简单的测试场景：对一个数值向量的每个元素做平方根运算，结果按原始顺序返回。在[Gist1](https://gist.github.com/YulongNiu/add0d9f066299613b64b8458fd5d741a)和[Gist2](https://gist.github.com/YulongNiu/9331ea0d3ef46f0571c5f2dc061c3f8a)中，分别实现了：

* `SqrtR`：用循环非并行操作每个元素。这种方法在`R`语言编程中不推荐，而应该尽量“向量化”操作。

* `SqrtRPara`：`foreach`并行版本。

* `SqrtCpp`：`C++`非并行版本。

* `SqrtCppPara`：<span style="color: blue">RcppParallel</span>包的`C++`并行版本。

* `sqrt`：R内置的向量化方法，`C`非并行版本。

首先，比较5种实现效率，并行计算调用8个线程（Intel i7-4790 CPU@3.6GHz）。测试结果显示`SqrtRPara`（使用`foreach`）和非向量化的R版本`SqrtR`效率较低。

{% codeblock lang:bash 5 versions %}
tmp1 <- runif(10e3)

all.equal(SqrtCpp(tmp1),
          sqrt(tmp1),
          SqrtR(tmp1),
          SqrtRPara(tmp1),
          SqrtCppPara(tmp1))

## TRUE

microbenchmark(
  SqrtCpp(tmp1),
  sqrt(tmp1),
  SqrtRPara(tmp1),
  SqrtR(tmp1),
  SqrtCppPara(tmp1))
  
## Unit: microseconds
##               expr         min           lq         mean       median
##      SqrtCpp(tmp1)      55.936      70.8940 8.512877e+01      76.3405
##         sqrt(tmp1)      35.801      46.1545 4.984745e+01      49.2785
##    SqrtRPara(tmp1) 1467517.454 1525986.5995 1.574319e+06 1571794.1885
##        SqrtR(tmp1)    3009.722    3092.6530 4.208856e+03    3165.9255
##  SqrtCppPara(tmp1)      23.432      57.3595 9.719465e+01      93.6050
##            uq         max neval
##       82.8530     837.569   100
##       53.2275      96.710   100
##  1607484.0645 1768153.919   100
##     3461.4475   55492.891   100
##      107.2245    1144.543   100
{% endcodeblock %} 

然后，增加循环数，比较效率较高的前三种方法。测试结果显示调用<span style="color: blue">RcppParallel</span>包的`C++`并行版本`SqrtCppPara`胜出。

{% codeblock lang:bash top 3 versions %}
tmp1 <- runif(10e6)

all.equal(SqrtCpp(tmp1),
          sqrt(tmp1),
          SqrtCppPara(tmp1))

## TRUE

microbenchmark(
  SqrtCpp(tmp1),
  sqrt(tmp1),
  SqrtCppPara(tmp1))
  
## Unit: milliseconds
##               expr      min       lq     mean   median       uq       max neval
##      SqrtCpp(tmp1) 76.68263 78.55146 82.51442 79.48709 87.03026 100.56873   100
##         sqrt(tmp1) 52.19705 53.67441 58.16940 54.60512 66.67642  70.94672   100
##  SqrtCppPara(tmp1) 37.10116 38.34199 42.23896 39.17889 42.98785  61.94529   100
{% endcodeblock %} 


### <a id="Ref">参考网址</a> ###

* [Summing a Vector in Parallel with RcppParallel](http://gallery.rcpp.org/articles/parallel-vector-sum/)
