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
library(devtools)
library(roxygen2)
{% endcodeblock %}

<!--more-->

## 2. 创建DESCRIPTION文件模板 ##

{% codeblock lang:r %}
load_all()
{% endcodeblock %}
其中`import`{:.language-bash}栏目，在源代码中使用了哪些包，需要逐步在import项目中添加和修改。


## 3. 检查包 ##

{% codeblock lang:r %}
check()
{% endcodeblock %}

## 4. 安装包 ##

{% codeblock lang:r %}
install()
{% endcodeblock %}

## 5. 生成`.tar.gz`{:.language-bash}压缩文件 ##

{% codeblock lang:r %}
build()
{% endcodeblock %}


### <a id="Ref">参考网址</a> ###

* [创建R包视频](https://www.youtube.com/watch?v=9PyQlbAEujY)


### 更新记录 ###

2014年9月5日
