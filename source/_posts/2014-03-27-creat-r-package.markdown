---
layout: post
title: "创建R包的最简流程"
date: 2014-03-27 14:53:05 -0400
comments: true
categories: r
---

一个最简单创建R包的流程。

## 1. 载入工具包 ##

安装并载入<span style="color: blue">devtools</span>包和<span style="color: blue">roxygen2</span>包。<span style="color: blue">devtools</span>包提供了一些包的检查、安装和打包的基本工具。<span style="color: blue">roxygen2</span>包则使得书写R帮助文档变得轻松简单。如果习惯使用Emacs，可以结合[ESS](http://ess.r-project.org/)建立R包，可以将R代码和帮助文档有效组合在一起，便于管理。当然，也可以使用[Rstudio](http://www.rstudio.com/)。

{% codeblock lang:r %}
library('devtools')
library('roxygen2')
{% endcodeblock %}

<!--more-->

## 2. <span style="color: blue">Rcpp</span>和相关的包 ##

如果使用了<span style="color: blue">Rcpp</span>或者相关的包，比如<span style="color: blue">RcppArmadillo</span>，需要格外设置。所有cpp代码都写在src文件夹下。

**首先**，执行：

{% codeblock lang:r %}
use_rcpp()
{% endcodeblock %}

**其次**，在`DESCRIPTION`{:.language-r}中添加依赖或者需要链接的包名称，比如：

{% raw %}
```
Encoding: UTF-8
Imports: 
    Rcpp
LinkingTo: 
    Rcpp,
    RcppArmadillo
```
{% endraw %}

**之后**，在包的`R/`{:.language-bash}目录下，添加一个文件`RcppChk.R`{:.language-bash}（文件名称自定），并写入：

{% codeblock lang:r %}
#' @useDynLib my-pkg-name, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
{% endcodeblock %}

注意，修改`my-pkg-name`{:.language-bash}为自己的包名称。

**然后**，在包的`src/`{:.language-bash}目录下，添加一个文件`registerDynamicSymbol.c`{:.language-bash}（文件名称自定），并写入：

{% codeblock lang:c %}
// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_markovchain(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
} 
{% endcodeblock %}

**最后**，可能需要在`~/.R/Makevars`文件下添加：

{% raw %}
```
# Settings from /etc/R/Makeconf with "non-portable flag(s):"
# ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’ replaced by -Wall -pedantic
# and without -fdebug-prefix-map=... 
CFLAGS = -g -O2 -Wall -pedantic -fstack-protector-strong -D_FORTIFY_SOURCE=2 $(LTO)
CXXFLAGS = -g -O2 -Wall -pedantic -fstack-protector-strong -D_FORTIFY_SOURCE=2 $(LTO)
CXX98FLAGS = -g -O2 -Wall -pedantic -fstack-protector-strong -D_FORTIFY_SOURCE=2
CXX11FLAGS = -g -O2 -Wall -pedantic -fstack-protector-strong -D_FORTIFY_SOURCE=2
CXX14FLAGS = -g -O2 -Wall -pedantic -fstack-protector-strong -D_FORTIFY_SOURCE=2
```
{% endraw %}

## 3. 创建DESCRIPTION文件模板 ##

{% codeblock lang:r %}
load_all()
{% endcodeblock %}
其中`import`{:.language-bash}栏目，在源代码中使用了哪些包，需要逐步在import项目中添加和修改。


## 4. 更新文档和检查包 ##

{% codeblock lang:r %}
document()
check()
{% endcodeblock %}

## 5. 安装包 ##

{% codeblock lang:r %}
install()
{% endcodeblock %}

## 6. 生成`.tar.gz`{:.language-bash}压缩文件 ##

{% codeblock lang:r %}
build()
{% endcodeblock %}


### <a id="Ref">参考网址</a> ###

* [创建R包视频](https://www.youtube.com/watch?v=9PyQlbAEujY)

* [R package by Hadley Wickham](http://r-pkgs.had.co.nz/src.html) 

* [Advanced R by Hadley Wickham](http://adv-r.had.co.nz/Rcpp.html) 

* [CRAN submission - R CMD Check warning - compilation flags used](https://stackoverflow.com/questions/50658198/cran-submission-r-cmd-check-warning-compilation-flags-used/52100124#52100124)

* [Warning about UTF-8 with roxygen2](https://stackoverflow.com/questions/51694929/warning-about-utf-8-with-roxygen2)

### 更新记录 ###

2019年04月19日
