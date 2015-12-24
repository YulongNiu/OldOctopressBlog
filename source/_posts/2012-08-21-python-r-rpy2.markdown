---
layout: post
title: "Python中使用rpy2模块调用R"
date: 2012-08-21 19:33:47 +0800
comments: true
categories: PPR
---

需要在python中调用R，实在是一种无奈的选择。如果能在一门语言中独立完成课题，是一个比较理想的做法。但是，这种想法也不太现实，毕竟每一种语言都有自己的长处。如果能取长补短，综合使用各种语言，也能起到不错的效果。

现在遇到的问题是，如何在python中调用R？这其中包括了如何调用R的对象（函数和包），R和python的对象如何互相转换，以及如何调用R的脚本（外界参数的输入）。python提供了一个模块[rpy2](http://rpy.sourceforge.net/) ，可以较好地完成这项工作。rpy2只是提供了一个Python调用R的接口，因此也直接继承了所有R的缺点。一个有意思的项目是[renjin](http://www.renjin.org/)，一个基于JVM的R语言解释器。

本文着重记录一些使用过程中的注意事项和小技巧，如有不一致则以[官方文档](http://rpy.sourceforge.net/rpy2/doc-2.5/html/)为准。

## 1. 安装 ##

rpy2作为Python的一个模块，其[安装](http://rpy.sourceforge.net/rpy2/doc-2.5/html/overview.html#installation)非常方便。需要做的准备工作是在本地安装好R。

{% codeblock lang:bash Install rpy2 from pip %}
# pip install rpy2
{% endcodeblock %}

<!--more-->

## 2. python调用R对象 ##

### 2.1 使用`rpy2.robjects`{:.language-python} ###

在rpy2中调用R对象，实际上是开启了一个R的交互进程。主要思路是将R的代码写入一个字符串内，之后执行即可。

对于R代码，将一段R代码写成一行，尽管看起来很丑陋而且不推荐，一样可以执行。但是，反过来，对于Python代码则没有这么简单。因为，Python代码是靠缩进来划分代码的区域，假若一段代码中有两个循环嵌套，如果将代码写成一行，执行起来就要麻烦的多（很可能要依赖空格的多少进行解释）。

调用方法：

{% codeblock lang:python Call R objects %}
from rpy2.robjects import r
{% endcodeblock %}


有三种方式可以选择：

1. 使用`r.obj`{:.language-python}，比如 `r.c(1, 3)`{:.language-python}。

   > 这种方法虽然方便，但是对于名称中有“点号”的函数会出问题，比如 `data.frame`{:.language-r}或者 `read.csv`{:.language-r}等。

2. 使用`r['obj']`{:.language-python}，比如 `r['c'](1, 3)`{:.language-python}。

   > * 这种方法几乎可以调用任何R的函数，而且写法与原始调用很相似，无非是`r['func'](value1, para2 = value2)`{:.language-python}。
   >
   > * 如果一个R函数中的变量名是有“点号”的，不能直接赋值，需要构建一个字典形式，比如 `r['func'](value1, para2 = value, **{para.3: value3})`{:.language-python}。
   >
   > * 如果一个
   
3. 使用`r('obj')`{:.language-python}，比如 `r('c(1, 3)')`{:.language-python}。

   > 这种方法从某种程度上讲是万能的，因为总是可以将任意长度的R代码写成一个Python字符串，之后通过`r('Rcode')`{:.language-python}调用执行。


{% codeblock lang:python Example of Calling R objects from rpy2.objects %}
# import r
>>> from rpy2.robjects import r

# creat an R function
>>> r('''f <- function(r){pi * r}''')
>>> r.f(3)
[9.424778]

# internal function in R
>>> r['ls']()

# two ways of getting 'paste' function in R
# first: direct use R function
>>> print(r['paste'](l, collapse = '-'))
# second: eavl R codes
>>> coder = 'paste(%s, collapse = "-")' % (l.r_repr())
>>> print(r(coder))

# call Python function
>>> l = r['letters']
>>> len(l)
26
>>> dir(l)
{% endcodeblock %}


### 2.2 创建R对象和提取对象的数值 ###

创建向量，可以使用`rpy2.robjects.vectors`{:.language-python}中提供的一系列函数，将一个Python的元组、列表或者字符串转换为R的向量。其中包括 `StrVector()`{:.language-python}、`IntVector()`{:.language-python}、`FloatVector()`{:.language-python}、`FactorVector()`{:.language-python}和`BoolVector()`，分别提供了转换字符、整数、浮点、因子和布尔向量。

创建列表，可以使用`ListVector()`{:.language-python}将Python的字典转换为R的列表。

创建矩阵和数据框（data.frame）建议直接使用R函数`matrix()`{:.language-r}和`data.frame()`{:.language-r}。尽管如此，rpy2提供了 `DataFram()`{:.language-python}将Python的字典转换成R的数据框（列顺序可能与输入不一致，Python字典特性造成），<span style="color: blue">注意</span>会把字符串自动转换成因子。

以上这些构建的R对象，有一系列的属性和方法，比如`names`{:.language-python}，可以访问和赋值。


{% codeblock lang:Python Build R objects %}
>>> from rpy2.robjects import r
>>> from rpy2.robjects.vectors import StrVector, IntVector, ListVector, DataFrame

# build R vector
>>> testVec = IntVector([1, 2, 3])
>>> testVec = StrVector(('a', 'nice', 'day'))
>>> testVec = StrVector('abc')
>>> testVec.names = StrVector(('name1', 'name2', 'name3'))
>>> dir(testVec)

# build List with names
>>> testList = ListVector({'a': 1, 'nice': 2, 'day': 3})
>>> testList.names

# build matrix
>>> m = r.matrix(IntVector([1, 3, 8, 6]), nrow = 2)
>>> m = r.matrix(range(10), nrow = 5)

# build data.frame
# use R function, string vector are automatically transformed as factors
>>> dataf = r['data.frame'](S = StrVector(['x', 'y', 'z']), F = StrVector('acb'))
>>> dataf
<DataFrame - Python:0x7f072235f7e8 / R:0x19fb250>
[FactorVector, FactorVector]
  S: <class 'rpy2.robjects.vectors.FactorVector'>
  <FactorVector - Python:0x7f072235f518 / R:0x1978fa8>
[       1,        2,        3]
  F: <class 'rpy2.robjects.vectors.FactorVector'>
  <FactorVector - Python:0x7f072235f440 / R:0x18e9c20>
[       1,        3,        2]
# use string vector "as it is"
>>> datav = DataFrame({'S': r.I(StrVector(['x', 'y', 'z'])), 'F': StrVector('acb')})
>>> datav
<DataFrame - Python:0x7f072235f638 / R:0x1155ff8>
[StrVector, FactorVector]
  S: <class 'rpy2.robjects.vectors.StrVector'>
  <StrVector - Python:0x7f07222d95f0 / R:0x23b0168>
[str, str, str]
  F: <class 'rpy2.robjects.vectors.FactorVector'>
  <FactorVector - Python:0x7f07222d9488 / R:0x10e7e70>
[       1,        3,        2]
{% endcodeblock %}


在rpy2中构建的对象，可以使用名字、索引或者布尔值进行提取或者修改，<span style="color: blue">注意</span>Python从0开始索引，而R从1开始。R对象的提供了方法`rx()`{:.language-python}相当于R的 `[`{:.language-r}操作，而`rx2()`{:.language-python}相当于`[[`{:.language-r}操作。


{% codeblock lang:python Extract values %}
>>> from rpy2.robjects import r
>>> from rpy2.robjects.vectors import StrVector, IntVector, BoolVector, ListVector, DataFrame

# vector
>>> testVec = IntVector([1, 2, 3])
>>> testVec.names[0]
>>> testVec.names = StrVector(('name1', 'name2', 'name3'))
# assign values
>>> testVec[0] = 20
>>> testVec[1: 3] = IntVector([100, 101])
>>> testVec[testVec.names.index('name2')] = 10
# using the rx() method
>>> testVec.rx('name2')
>>> testVec.rx(3)
>>> testVec.rx(-1)
>>> testVec.rx(IntVector([-1, -3]))
>>> testVec.rx(BoolVector([False, True, False]))
# empty vector returns
>>> testVec.rx(0)
# error because R cannot determine the last element 
>>> testVec[1: ] = IntVector([100, 101])
# error because rx() is a method
>>> testVec.rx(IntVector([1, 3])) = IntVector([100, 101])

# list and matrix
>>> tmp = r("list(a = matrix(1:10, nrow = 2), b = 'Hello')")
>>> print tmp
$a
     [,1] [,2] [,3] [,4] [,5]
[1,]    1    3    5    7    9
[2,]    2    4    6    8   10

$b
[1] "Hello"
>>> tmp.names
>>> tmp.names[1]
>>> tmp.rx('a')
<ListVector - Python:0x8afd86c / R:0x8cf71c0>
[Matrix]
  a: <class 'rpy2.robjects.vectors.Matrix'>
  <Matrix - Python:0x8b013cc / R:0x97de388>
[       1,        2,        3, ...,        8,        9,       10]
>>> tmp.rx2('a')
# same as the former one
>>> tmp.rx(1)
>>> tmp.rx2(1)
<Matrix - Python:0x8b01b4c / R:0x97de388>
[       1,        2,        3, ...,        8,        9,       10]
>>> tmp[tmp.names.index('b')] = 9

# operate matrix
>>> tmpMat = tmp.rx2('a')
# first element of 'a'
>>> tmpMat.rx(1, 1) 
# first row of 'a'
>>> tmpMat.rx(1, True)
>>> tmpMat.rx
>>> tmpMat.colnames = r['letters'][1:6]
>>> tmpMat.rx(True, IntVector([1, 3]))
>>> tmpMat.rx(True, BoolVector((True, False, True, False, False)))
>>> tmpMat.rx(True, 'b')
{% endcodeblock %}


### 2.3 特殊对象 ###

* R中的`NA`{:.language-r}通过`rpy2.robjects.NA_Integer`{:.language-python}或`rpy2.robjects.NA_Character`{:.language-python}等引用。

* 全局变量使用`rpy2.robjects.globalenv`{:.language-python}。

   > 特别是遇到Python找不到一个R对象时（此时R对象可能通过`r('Rcode')`{:.language-python}调用），留意将R对象转变成全局变量。

{% codeblock lang:python Special objects %}
>>> from rpy2.robjects import r
>>> from rpy2.robjects.vectors import StrVector
>>> from rpy2.robjects import globalenv, NA_Character

# test string NA
>>> testNA = StrVector(['a', 'b', 'c'])
>>> testNA[1] = NA_Character

# test global objects
>>> testGlobal = StrVector(['a', 'b', 'c'])
>>> globalenv['testGlobal'] = testGlobal
>>> connectVec = r("paste(testGlobal, collapse = '-')")
{% endcodeblock %}


### 2.4 载入和使用包 ###

使用`rpy2.robjects.packages.importr`{:.language-python}调用R包：

{% codeblock lang:python import R packags %}
>>> from rpy2.robjects.packages import importr

>>> base = importr('base')
>>> stats = importr('stats')
>>> affy = importr('affy')
>>> stats.rnorm(10)
{% endcodeblock %}

如果想引用一个包中的隐变量，也很简单，只要载入包，然后所有r命令化成成字符串，之后引用即可（这种方法是万能的），比如：

{% codeblock lang:python Hidden objects in R packages  %}
>>> from rpy2.robjects.packages import importr

>>> importr('hwriter')
>>> a = r('hwriter:::hwrite.table(matrix(1:10, 2))')
>>> print(a)
[1] "<table border=\"1\">\n<tr>\n<td>1</td><td>3</td><td>5</td><td>7</td><td>9</td></tr>\n<tr>\n<td>2</td><td>4</td><td>6</td><td>8</td><td>10</td></tr>\n</table>\n"
{% endcodeblock %}


## 3. R对象转换成Python对象 ##

推荐使用`tuple()`{:.language-python}或者`list()`{:.language-python}函数，将R对象转换成Python的元组或者列表对象。


{% codeblock lang:python Transform R objects to Python ones %}
>>> a = r('c(1, 2, 3)')
>>> str(a)
'[1] 1 2 3\n'
>>> tuple(a)
(1.0, 2.0, 3.0)
>>> list(a)
[1.0, 2.0, 3.0]
>>> b = r('matrix(1:6, 2, 3)')
>>> print b
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
>>> tuple(b)
(1, 2, 3, 4, 5, 6)
>>> list(b)
[1, 2, 3, 4, 5, 6]
{% endcodeblock %}


## 5. 注意事项 ##

* 分清楚需要调用的R还是Python的对象；

* 善于使用构建好的rpy2对象的属性和方法，比如`rpy2Obj.len`{:.language-python}；

* 如果函数有警告（warnings），在ipython等IDE上能够执行，但是如果是脚本或者与网页服务器交互，则会产生错误。解决办法：

   > 1. 鲁莽的解决很简单，强行忽略R的警告，`options(warn = -1)`{:.language-r}或者R代码放入函数中`suppressWarnings()`{:.language-r}。
   >
   >  2. 第二种办法，如果是自己代码中使用了R的 `warning()`{:.language-r}函数，则将warning信息换成字符串，之后单独输出。

* rpy2不是万能药，它直接继承了R的所有好和不好。



### 参考资料 ###

* [stackoverflow: Python and Rpy2: Calling plot function with options that have “.” in them](http://stackoverflow.com/questions/2125218/python-and-rpy2-calling-plot-function-with-options-that-have-in-them) 


### 更新记录 ###

2015年7月26日



