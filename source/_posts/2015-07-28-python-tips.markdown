---
layout: post
title: "Python使用小技巧"
date: 2015-07-28 14:37:38 +0800
comments: true
categories: PPR
---

收集了一些使用Python过程中的小技巧或者常见错误。

## 1. While--if--break ##

适用情况：**当需要执行一个循环，但是不能事先判断何时终止。** 可以在循环内部加入判断，符合要求时跳出。比如，使用程序在某个目录下新建一个文件夹，如果文件名已经存在，需要重新生成文件名；否则，创建文件夹。

{% codeblock lang:python While--if--break %}
import random, os

while True:
    letter = [chr(i) for i in range(97, 123)]
    folderName = [random.choice(letter) for i in range(5)]
    fn = ''.join(folderName)
    if os.path.exists(fn) is not True:
        # get an unique name 
        os.mkdir(fn)
        break
{% endcodeblock %}


## 2. 链式比较操作 ##

Python支持一个链式的比较操作。

{% codeblock lang:python comparison chain %}
testNum = 10

# 链式比较操作
# 区别if not和if a is not None
if (not 0 < testNum < 5):
    print('It passed!')
{% endcodeblock %}




<!--more-->

### 参考网址 ###

* [Python Wiki](https://wiki.python.org/moin/FrontPage)

* [Wei Shen's Python note](http://blog.shenwei.me/python-note/#more-3951)

* [十个我希望早点知道的python方法](http://www.codefrom.com/c/214)



### 更新记录 ###

2015年7月27日




