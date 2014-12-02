---
layout: post
title: "Linux安装R语言包"
date: 2010-09-04 10:21:38 -0400
comments: true
categories: R
---

## 1. R包介绍 ##

R的包（package）通常有两种:

* 二进制代码包（Binary package）：这种包属于即得即用型（ready-to-use），但是依赖与平台，比如Windows和Linux平台下不同。

* 源代码包（Source package）: 此类包可以跨平台使用，但用之前需要处理或者编译（compiled）。同时，源代码包可以查看到程序源代码，便于查找、修改和引用。

## 2. R包安装 ##

### 2.1 源代码安装 ###

{% codeblock lang:bash %}
# R CMD INSTALL /.../myPackage.tar.gz
{% endcodeblock %}

使用此方法，需要解决包依赖问题，即安装此包所依赖的包，安装过程有提示。

<!--more-->

### 2.2 内置`install.packages()`{:.language-bash}函数安装###

使用`install.packages()`安装，比较简便，联网即可装，装了就可用。使用这种方法安装包时，R会自动安装依赖的包。如果出现安装报错，可能的原因是缺少依赖的系统文件。同时，需要注意的是，一些R包只能在特定的平台上使用。比如[Rsubread](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)不能在Windows操作系统下使用。

{% codeblock lang:r %}
install.packages('myPackage')
{% endcodeblock %}

同时，可以使用`install.packages()`{:.language-bash}安装本地下载的包，尤其适用于在服务器上安装包

{% codeblock lang:r %}
install.packages(
  c('XML_0.99-5.tar.gz', '../../Interfaces/Perl/RSPerl_0.8-0.tar.gz'),
  repos = NULL,
  configure.args = c(XML = '--with-xml-config=xml-config', RSPerl = '--with-modules= "IO Fcntl"')) 
{% endcodeblock %}

## 3. R包版本查询和更新 ##

R和R包版本查询

{% codeblock lang:r %}
# 在启动的R中执行
R.version

# R包版本
packageVersion('myPackage')

# 查询当前R的详细信息，包括R版本、R包版本、命名空间等
sessionInfo()
{% endcodeblock %}


[CRAN](http://cran.r-project.org/)包更新

{% codeblock lang:r %}
# 可以定期执行以下
update.packages()  
{% endcodeblock %}

[Bioconductor](http://www.bioconductor.org/)的安装和更新方法

{% codeblock lang:r %}
source('http://bioconductor.org/biocLite.R')
biocLite('myPackage')
{% endcodeblock %}

## 4. 卸载R包##

{% codeblock lang:r %}
remove.packages('myPackage')
{% endcodeblock %}

## 5. R包相关函数 ##

{% codeblock lang:r %}
# 查看包的安装目录
.libPaths()

# 查看已经安装的包目录
library()

# 查看已安装包信息
installed.packages()

# 载入myPackage包
library(myPackage)
require(myPackage)

# 查看当前载入的包
search()

# 查看启动R时自动载入的包。
getOption('defaultPackages')
{% endcodeblock %}

## 6. 帮助信息查询 ##

### 6.1 R和R包帮助信息 ###

{% codeblock lang:r %}
# 查询R HOME安装地址
Sys.getenv('R_HOME')

# 查询用户HOME地址
Sys.getenv('HOME')

# 查看某个“函数”或者“方法”的详细内容
?myFunction
?myMethod

# 关键词查询
??myKeyword

# 查看已经安装包的详细HTML文档
help.start()

# 搜索R网站上的“helpinfor”相关信息
RSiteSearch('helpinputinfor')

# 查看“myPackage”的帮助
help(package = 'myPackage')

# 有的包，特别是bioconductor的包有vignette，用函数查看
vignette('myPackage')

# 这个函数也可以查看vignette，更好用一些
openVignette('myPackage')

# 展示一些包中demostration
demo('package')
{% endcodeblock %}

### 6.2 查询对象信息 ###

{% codeblock lang:r %}
# 查看"myPackage"中的所有对象
ls('package:myPackage')

# 查看函数的参数
args(myFunction)

# 自动运行该函数帮助文档中的例子
example(myFunction)

# 查看某个对象的模式（mode）
mode(myObject)

# 查看某个对象的属性（attribute）
attributes(myObject)

# 快速查看某个对象的信息
# 尤其适用于对象有很多行/列
str(myObject)

# 查看某个对象的类
class(myObject)

# 查询某个中某个类的帮助信息，举例如下
class?graph::graph

# 查看某个S3泛型函数中所有的方法或者一个类中所有的方法（S3：S version 3）
methods('myMethods')

# 查看S4类的方法
showMethods(class = 'myClass')

# 查看某个类或者包的具体内容
getClass('class/package')

# 查看某个类的slot
getSlots('class')

# 查看某个对象的slot
slotNames(myObject)

# 访问对象的slot值使用@，可以连续用
Myobject@slotNames
{% endcodeblock %}

## 7. 查看函数源代码 ##

### 7.1 普通函数源代码 ###

直接输入函数名称，不加后面的括号。比如:

{% codeblock lang:r %}
> fivenum
{% endcodeblock %}

### 7.2 查询S3/S4函数源代码 ###

因为如果使用了R的S3和S4方法，可能一个函数名称对应的多个不同的方法（对应不同的对象）。这时，采用的策略是：

1. 查询函数对应的方法

2. 查询具体方法对应的源代码

比如`t.test`{:.language-r}的源代码

{% codeblock lang:r %}
> t.test
function (x, ...) 
UseMethod("t.test")
<bytecode: 0x244a248>
<environment: namespace:stats>

# 查询S3对象的方法
> methods('t.test')
[1] t.test.default* t.test.formula*

   Non-visible functions are asterisked
   

# 查看源代码，三种方法
> getAnywhere('t.test.default')
> getS3method('t.test','default')
> stats:::t.test.default
{% endcodeblock %}

另外一个例子，查询<span style="color: blue">affy</span>包中`pm()`{:.language-r}的源代码：

{% codeblock lang:r %}
> library(affy)
> pm
standardGeneric for "pm" defined from package "affy"

function (object, ...) 
standardGeneric("pm")
<environment: 0x4496f10>
Methods may be defined for arguments: object
Use  showMethods("pm")  for currently available ones.

> findMethods('pm')
{% endcodeblock %}

### <a id="Ref">参考网址</a> ###

* [R包安装](http://yanping.me/cn/blog/2012/01/09/customizing-the-startup-environment/)

* 《R语言编程技术》(The Art of R Programming)


### 更新记录 ###

2014年11月25日
