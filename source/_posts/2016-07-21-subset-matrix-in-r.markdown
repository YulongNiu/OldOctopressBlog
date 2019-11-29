---
layout: post
title: "为什么不推荐在R语言中随意按照下标操作矩阵"
date: 2016-07-21 19:35:29 +0800
comments: true
categories: R
---

这篇博文的目的是展示R语言中下标操作矩阵的潜在问题。R语言提供了多种方法提取一个矩阵的单个或者部分元素，不同方法对应的效率在Hadley Wickham的[Advance R](http://adv-r.had.co.nz/Performance.html#language-performance)中已有讨论。这些方法中，使用最广泛的是通过下标（行或者列）取值，即操作符`[`{:.language-R}。然而，这种方法存在潜在问题，即内存中会拷贝原始对象。

举例：首先建立一个矩阵，之后取这个矩阵除了第一行之外的部分，接下来操作这个部分矩阵。

{% codeblock lang:R manipulate %}
## step1: build matrix
n <- 8000
tmp1 <- matrix(rnorm(n * n), nrow = n, ncol = n)
gc()

## step2: manipulate a subset of matrix
sink('/dev/null')
apply(tmp1[2:n, ], 1, function(x) x[1])
sink()

## step3: garbage collection
gc()
{% endcodeblock %}

<!--more-->

内存使用情况如下：

<img src="/images/R_apply_memory.png" title="image" alt="UCSC下载rRNA注释">

* 标记1内存上升，因为建立了`tmp1`{:.language-R}的矩阵；

* 标记2内存再次上升，主要因为使用下标取矩阵操作，`tmp1[2:n, ]`{:.language-R}；

* 标记3内存下降，因为手动执行垃圾回收。

可以明显看到内存中多余的垃圾对象。如果使用`for`{:.language-R}循环形式，就可以有效避免内存对象拷贝。虽然，R在内存空间不足时，会自动执行`gc()`{:.language-R}。但是，执行程序时，不能全指望自动垃圾回收，毕竟有时回收得并不及时，而新的对象又相继生成。这种情况下，内存空间不足就成为很严重的问题。

R版本3.3.1。



### 更新记录 ###

2016年7月21日








