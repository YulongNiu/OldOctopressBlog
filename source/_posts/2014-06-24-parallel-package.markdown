---
layout: post
title: "R使用parallel包并行计算"
date: 2014-06-24 22:10:20 -0400
comments: true
categories: R
---

最新版本的R已经内置<span style="color: blue">parallel</span>包，<span style="color: blue">parallel</span>包是从<span style="color: blue">[snow](http://cran.r-project.org/web/packages/snow/index.html)</span>包和<span style="color: blue">[multicore](http://cran.r-project.org/web/packages/multicore/index.html)</span>包继承而来，包含了很多非常好用的函数。<span style="color: blue">parallel</span>包可以通过PVM（<span style="color: blue">rpvm</span>包）、MPI（<span style="color: blue">[Rmpi](http://cran.r-project.org/web/packages/Rmpi/index.html)</span>包）、NetWorkSpaces（<span style="color: blue">[nws](http://cran.r-project.org/web/packages/nws/index.html)</span>包）和raw sockets（如果以上3种都不能使用）平台进行分布计算，支持cluster和多核个人/服务器计算机。在Linux系统上，通常使用[openMPI](http://www.open-mpi.org/)。


## 1. 安装<span style="color: blue">Rmpi</span>包 ##

因为使用openMPI，所以<span style="color: blue">parallel</span>包需要<span style="color: blue">Rmpi</span>包来设定节点，所以首先需要在计算机上安装openMPI。

<!--more-->

### 1.1 Linux系统下安装openMPI环境 ###

~~~ bash
# 安装openmpi环境
# yum install openmpi openmpi-devel

# 配置环境（安装时执行，可能之后运行也要执行）
# ldconfig /usr/lib64/openmpi/lib/
~~~

在`~/.bashrc`{:.language-bash}下写入

~~~ bash
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:}/usr/lib64/openmpi/lib/"
~~~

载入`~/.bashrc`{:.language-bash}

~~~ bash
$ source ~/.bashrc
~~~

### 1.2 安装Rmpi包 ###

在启动的R窗口中输入：

~~~ r
install.packages("Rmpi",
                 configure.args =
                 c("--with-Rmpi-include=/usr/include/openmpi-x86_64/",
                   "--with-Rmpi-libpath=/usr/lib64/openmpi/lib/",
                   "--with-Rmpi-type=OPENMPI"))
~~~

## 2. 使用<span style="color: blue">parallel</span>包 ##

### 2.1 设定节点数 ###

首先，需要设定cluster的节点（nodes）数目

~~~ r
cl <- makeCluster(2, type = "MPI")
~~~

这里对“节点数”设定做一些解释，如果使用cluster，可以直接设定cluster数据即可；如果是在小型服务器或者个人电脑上使用，最大节点数可以设定为“线程数（processor）-1”。比如一个双核四线程计算机，节点数目最大可以设定为3。这是因为<span style="color: blue">snow</span>包（<span style="color: blue">parallel</span>包的主要依赖包）在设计时，总是要保留一个**“主线程”**来处理和整合数据。

在linux系统下，线程数可以通过 `$ nproc`{:.language-bash} 查看。

### 2.2 内置函数 ###

使用<span style="color: blue">parallel</span>包中的内置并行运算函数
比如使用`parApply()`{:.language-r}、`parCapply()`{:.language-r}、`parRapply()`{:.language-r}、`parLapply()`{:.language-r}和`parSapply()`{:.language-r}（如果返回矩阵，使用
`cbind()`{:.language-r}）等函数。其中文档中指出，`parApply()`{:.language-r}函数对于**二维矩阵**的每一个单元进行操作，因此要慢一些。如果可能，使用`parCapply()`{:.language-r}和`parRapply()`{:.language-r}对列和行进行操作，以加快运行速度。

### 2.3 回收节点 ###

~~~ r
stopCluster(cl)
~~~

### 3. 并行计算的包依赖问题 ###

在并行计算过程中，不可避免地会用到其他包辅助。这里涉及到<span style="color: blue">snow</span>包的一个设计原理：并行运算多个R进程，只有一个主进程载入完整的依赖包环境。这就意味着其他并行的R进程中也要载入依赖的包环境。

有两个思路，第一个思路是修改`Rprofile.site`{:.language-bash}文件，让任意R进程在启动时都载入依赖的包。但是，不推荐这种做法，因为这样会增加R载入的速度；并且如果不同的代码用了不同的依赖包，就要不停地修改`Rprofile.site`{:.language-bash}文件。

第二个思路是在新开的R进程中“动态”加载需要的包。所谓**“动态”**，没有什么高深的意思，就是“需要的时候加载即可”。根据需要，可以选择以下两种方法。

* 第一种方法是在直接在启动的R进程中加载包。

这种方法非常直观，推荐。

~~~ r
# 以下代码摘抄自Parallel R，其中packages
# 是一个要选择加载的package列表，
# 比如c('bigmemory', 'foreach', 'doMC')
worker.init <- function(packages) {
  for (p in packages) {
    library(p, character.only=TRUE)
  }
  NULL
}
~~~

* 第二种方法是在调用函数中加入。

这种方法不推荐，因为我们将看到这种方法是“投机”了<span style="color: blue">parallel</span>包的并行`apply`家族函数。原理是：<span style="color: blue">parallel</span>包中最主要的就是`apply`家族函数，比如`parApply(cl = NULL, X, MARGIN, FUN, ...)`{:.language-r}函数，是<span style="color: blue">base</span>包中`apply()`{:.language-r}的并行版本。其中会用到一个`FUN`函数，我们可以在这个函数中加载包，比如写入`require('bigmemory')`{:.language-r}等。这样，并行的R进程就会载入需要的包。举例如下：

~~~ r
Getft <- function(i, arg1, arg2){
  require(packages)
  ...
}

adft <- parSapply(cl, 1:10, Getft, argInput1, argInput2)
~~~


## 4. 与<span style="color: blue">bigmemory</span>包结合 ##

<span style="color: blue">parallel</span>包可以很好地与<span style="color: blue">bigmemory</span>包结合，进而进一步提升R操作大数据的能力。

但是，有一个问题是`parApply()`{:.language-r}、`parCapply()`{:.language-r}、`parRapply()`{:.language-r}函数是不能直接调用<span style="color: blue">bigmemory</span>包的`big.memory`这种S4对象。当然也可以使用`mat[, ]`之类语句引用big.matrix对象。但是这会把矩阵全部载入内存，也就失去了`big.matrix`对象的意义，只有在内存允许的情况下这样操作。

解决办法：

> 1. 将`big.matrix`对象的操作放在一个函数中，函数传入的是`big.matrix`的`description file`（描述文件），而不是`big.matrix`对象本身。
> 
> 2. 把这个函数作为`parLapply()`{:.language-r}和`parSapply()`{:.language-r}的`FUN`，达到分布计算，而又不直接引用`big.matrix`对象的目的。
> 
> 3. 这个思路的前提是：创建的`big.matrix`对象必须是“**内存共享**”的，否则不能将其分布到不同的节点计算。

举一个例子，完整版本见[补充材料：Final version](#final_version)，这个例子中首先创建一个`Getft()`{:.language-r}函数，接受`adAllRowColDesc`和`adMatDesc`两个变量是`big.matrix`对象的描述文件。在这个函数中，`attach.big.matrix()`{:.language-r}通过描述文件引用`big.matrix`对象，并完成相关操作。

~~~ r
Getft <- function(i, adAllRowColDesc, adMatDesc){
  adAllRowColData <- attach.big.matrix(adAllRowColDesc)
  adMatData <- attach.big.matrix(adMatDesc)
  rowIndex <- adAllRowColData[i, 1]
  colIndex <- adAllRowColData[i, 2]
  linkData <- c(rowNames[rowIndex], rowNames[colIndex], adMatData[rowIndex, colIndex])
  return(linkData)
}
~~~
之后，使用`parSapply()`{:.language-r}函数调用`Getft()`{:.language-r}函数，使用`1:nrow(adAllRowCol)`作为“**计数器**”。

~~~r
adft <- parSapply(cl, 1:nrow(adAllRowCol), Getft, adAllRowColDescFile, adMatDescFile)
~~~

如果需要处理的`big.matrix`对象不大，也可以直接使用`parSapply()`{:.language-r}函数，详细参考[补充材料：Bigmatrix direct version](#bigmatrix_direct)。

## 5. 与<span style="color: blue">foreach</span>包比较##

另外一个支持并行计算的包是<span style="color: blue">foreach</span>包，它天生与`big.matrix`对象匹配。所以，我也提供了使用`foreach`{:.language-r}版本，详细参考[补充材料：Foreach version](#foreach)。

通过测试可以发现，在数据量较少时（1000行左右），`foreach`{:.language-r}[版本](#foreach)和`parSapply()`{:.language-r}[版本](#final_version)速度基本持平。但是，数据量增大时（百万行），`foreach`{:.language-r}[版本](#foreach)速度明显减慢。原因是在使用`foreach`{:.language-r}并行计算时，计算开始时候需要建立索引。这个过程在循环数变大时，会变得非常缓慢。

因此，我们可以看到，如果使用<span style="color: blue">foreach</span>包，会减少代码量，而且程序逻辑也非常清晰，但是遇到超大循环数，速度明显减慢。同时，如果使用<span style="color: blue">parallel</span>包，那么需要一些“技巧”才能与`big.matrix`对象有效融合。所以，我们的结论是原生态的R（包括提供的一些包）不适合做并行大数据计算。

### <a id="Ref">参考网址</a> ###

* [R task: High-Performance and Parallel Computing with R](http://cran.r-project.org/web/views/HighPerformanceComputing.html)

* [State of the Art in Parallel Computing with R](http://www.jstatsoft.org/v31/i01/paper)，这篇文章详细介绍了很多R并行计算的平台

* Rmpi安装：[1](http://www.stats.uwo.ca/faculty/yu/Rmpi/), [2](https://www.sharcnet.ca/help/index.php/Using_R_and_MPI), [3](http://www.cybaea.net/Blogs/R-tips-Installing-Rmpi-on-Fedora-Linux.html)

* [snow包介绍](http://www.sfu.ca/~sblay/R/snow.html)

* [The R Package bigmemory: Supporting Efficient Computation and Concurrent Programming with Large Data Sets.](http://www.stat.yale.edu/~mjk56/temp/bigmemory-vignette.pdf)

* 书籍[Parallel R](http://shop.oreilly.com/product/0636920021421.do)

### <a id="appendix">补充材料</a> ###

* <a id="final_version">Final version</a>

{% codeblock combine parSapply and big.matrix lang:r %}
adj2ftBig <- function(adMat, adAllRowCol, n = 2){
  
  # INPUT: 'adMat' should be a bigmatrix. 'adAllRowCol' is the row and column combination, also is a bigmatrix
  
  require(bigmemory)
  require(parallel)
  cl <- makeCluster(n, type = "MPI")

  adMatDescFile <- describe(adMat)
  adAllRowColDescFile <- describe(adAllRowCol)
  
  rowNames <- rownames(adMat)
  colNames <- colnames(adMat)

  ignore <- clusterEvalQ(cl, {library(bigmemory); NULL})

  Getft <- function(i, adAllRowColDesc, adMatDesc){
    adAllRowColData <- attach.big.matrix(adAllRowColDesc)
    adMatData <- attach.big.matrix(adMatDesc)
    rowIndex <- adAllRowColData[i, 1]
    colIndex <- adAllRowColData[i, 2]
    linkData <- c(rowNames[rowIndex], rowNames[colIndex], adMatData[rowIndex, colIndex])
    return(linkData)
  }

  adft <- parSapply(cl, 1:nrow(adAllRowCol), Getft, adAllRowColDescFile, adMatDescFile)

  stopCluster(cl)

  return(adft)
}
{% endcodeblock %}


* <a id="bigmatrix_direct">Bigmatrix direct version</a>

{% codeblock directly use the big.matrix lang:r %}
adj2ftBig3 <- function(adMat, adAllRowCol, n = 2){
  
  # INPUT: 'adMat' is a matrix. 'adAllRowCol' is the row and column combination, also a matrix.
  
  library(parallel)
  cl <- makeCluster(n, type = "MPI")
  
  rowNames <- rownames(adMat)
  colNames <- colnames(adMat)
    
  adft <- parRapply(cl = cl, adAllRowCol, function(x) {
    linkData <- c(rowNames[x[1]], colNames[x[2]], adMat[x[1], x[2]])
    return(linkData)
  })
  
  stopCluster(cl)

  return(adft)
}
{% endcodeblock %}


* <a id="foreach">Foreach version</a>

{% codeblock apply foreach to big.matrix lang:r %}
adj2ftBig2 <- function(adMat, adAllRowCol, n = 4){
  
  # INPUT: 'adMat' should be a bigmatrix. 'adAllRowCol' is the row and column combination, also a bigmatrix.
  
  library(bigmemory)
  library(foreach)
  library(doMC)
  registerDoMC(n)

  rowNames <- rownames(adMat)
  colNames <- colnames(adMat)

  adft <- foreach (i = 1:nrow(adAllRowCol), .combine = rbind, .inorder=TRUE) %dopar% {
    print(paste('It is running ', i, ' in total of ', nrow(adAllRowCol), '.', sep = ''))
    linkData <- c(rowNames[adAllRowCol[i, 1]], colNames[adAllRowCol[i, 2]], adMat[adAllRowCol[i, 1], adAllRowCol[i, 2]])
    return(linkData)
  }
    
  return(adft)
}
{% endcodeblock %}

### 更新记录 ###

2014年7月22日
