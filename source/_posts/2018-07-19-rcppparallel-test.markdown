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

## 1. 测试 ##

我尝试使用<span style="color: blue">RcppParallel</span>包调用`C++`的并行方法。结论是：**在循环数很大时，<span style="color: blue">RcppParallel</span>包提供的并行方法优于`foreach`**。

一个简单的测试场景：对一个数值向量的每个元素做平方根运算，结果按原始顺序返回。在[Gist1](https://gist.github.com/YulongNiu/add0d9f066299613b64b8458fd5d741a)和[Gist2](https://gist.github.com/YulongNiu/9331ea0d3ef46f0571c5f2dc061c3f8a)中，分别实现了：

* `SqrtR`：用循环非并行操作每个元素。这种方法在`R`语言编程中不推荐，而应该尽量“向量化”操作。

* `SqrtRforeach`：`foreach`并行版本。

* `SqrtRParSapply`: `parSapply`并行版本。

* `SqrtCpp`：`C++`非并行版本。

* `SqrtCppPara`：<span style="color: blue">RcppParallel</span>包的`C++`并行版本。

* `sqrt`：R内置的向量化方法，`C`非并行版本。

首先，比较5种实现效率，并行计算调用8个线程（Intel i7-4790 CPU@3.6GHz）。测试结果显示`SqrtRPara`（使用`foreach`）和非向量化的R版本`SqrtR`效率较低。

{% codeblock lang:bash 5 versions %}
tmp1 <- runif(10e3)

all.equal(SqrtCpp(tmp1),
          sqrt(tmp1),
          SqrtR(tmp1),
          SqrtRforeach(tmp1),
          SqrtRParSapply(tmp1),
          SqrtCppPara(tmp1))

## TRUE

microbenchmark(
    SqrtCpp(tmp1),
    sqrt(tmp1),
    SqrtR(tmp1),
    SqrtRforeach(tmp1),
    SqrtRParSapply(tmp1),
    SqrtCppPara(tmp1))

## Unit: microseconds
##                 expr         min          lq         mean       median
##        SqrtCpp(tmp1)      56.295      72.648 9.338755e+01      82.0335
##           sqrt(tmp1)      36.216      46.074 4.865115e+01      48.3090
##          SqrtR(tmp1)    3030.682    3116.718 4.229971e+03    3947.9380
##   SqrtRforeach(tmp1) 1488851.181 1532937.096 1.561865e+06 1547849.9610
## SqrtRParSapply(tmp1)  954757.348  963478.755 9.701841e+05  969925.9090
##    SqrtCppPara(tmp1)      23.837      79.314 1.069003e+02     104.5975
##           uq         max neval
##      89.0800    1183.279   100
##      52.2995      66.875   100
##    4560.0930   10391.379   100
## 1584297.5760 1750382.995   100
##  974233.5690 1012400.281   100
##     111.9160    1331.442   100

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

## 2. 使用`vector`代替`List` ##

在使用<span style="color: blue">RcppParallel</span>并行计算时，不能在并行循环中调用`Rcpp::List`对象。一个解决办法是：使用`std::vector`替代`Rcpp:List`。例如，`List`中都是数值向量，那么可以建立`std::vector<Rcpp::NumericVector>`对象替代。
[Gist3](https://gist.github.com/YulongNiu/0a11282216162b6e350c9575b68e91cc)中提供了一个例子。这种方法的局限在于`List`中每一个元素的类型需要相同。

## 3. 同步 ##

如果多个线程同时操作某一个共享内存对象，需要在<span style="color: blue">RcppParallel</span>包中使用“锁”。如[Gist4](https://gist.github.com/YulongNiu/5af268df461c8890c73c9640ae9ac754)所示，多个线程都需要操作`estcount`对象。测试代码如下：

{% codeblock lang:r Synchronization %}
library('Rcpp')

sourceCpp('TestSynchron.cpp')

n <- 10
g <- 1000
ecin <- sample(0:9, g*n, replace = TRUE) %>%
  split(1:g)

ecin %>% unlist %>% table

TestShare(ecin, 10)
{% endcodeblock %} 

如果去掉[Gist4](https://gist.github.com/YulongNiu/5af268df461c8890c73c9640ae9ac754)代码中的第`27`和`31`行，可以发现测试结果不正确。

### <a id="Ref">参考网址</a> ###

* [Summing a Vector in Parallel with RcppParallel](http://gallery.rcpp.org/articles/parallel-vector-sum/)

* [Intel TBB Simple Mutex - Example](https://scc.ustc.edu.cn/zlsc/tc4600/intel/2017.0.098/advisor/help/GUID-D98B389E-61B9-414A-9450-D28EF9F61A95.htm)

### 更新记录 ###

2018年10月5日
