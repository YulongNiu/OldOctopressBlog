---
layout: post
title: "C语言学习记录"
date: 2014-12-07 18:23:16 -0500
comments: true
published: true
categories: c
---

## 1. 基本数据类型 ##

### 1.1 数字形式 ###

数字的位数是多少？

* `int`：整数

* `short`：短整数

* `long`：长整数

* `float`：浮点数

* `double`：双浮点数

* `char`：字符

<!--more-->

### 1.2 其他形式 ###

* `arrays`

* `structures`

* `unions`

* `pointers`

* `functions`


## 2. 操作类型和顺序 ##

* 赋值操作`=`：从右向左

* 或操作`||`、且操作`&&`：从左向右，遇到第一个`true`即停止




## 3. 语法注意事项 ##

* `#define LOWER 0`定义常量的语句之后，没有分号`;`

* `EOF`是文档结束的标志（End of File），在`<stdio.h>`{:.language-c}定义为一个整数`-1`，代码如下：

{% codeblock definition of EOF in C lang:c %}
#ifndef EOF
# define EOF (-1)
#endif
{% endcodeblock %}

* 用引号`''`标记的字符串，对应的是ASCII的数字值

* 赋值例如`c = getchar()`{:.language-c}也有值，其值等于等号左侧赋值操作后的值

* `for`循环声明需要有主体，如无，在`for`语句后添加`;`作“无效声明（null statement）”

* 矩阵计数从0开始


## 4. 标准库 ##

### 4.1 `#include <stdio.h>`{:.language-c} ###

{% codeblock lang:c %}
printf("%3.0f %6.1f\n", fahr, celsius);
{% endcodeblock %}

* `%d`：十进制整数

* `%ld`：长十进制整数

* `%6d`：十进制整数，至少6位宽

* `%f`：浮点数

* `%6f`：浮点数，至少6位宽

* `%6.0f`：浮点树，至少6位宽，无小数点而且无小数位

* `%.2f`：浮点数，小数点后两位

* `%6.2f`：浮点数，至少6位宽，小数点后两位

* `%o`：八进制

* `%x`：六进制

* `%c`：字符

* `%s`：字符串

* `%%`：`%`自身


### 参考资料 ###

* BW Kernighan and DM Ritchie: [The C Programming Language (2nd Edition)](http://www.amazon.com/The-Programming-Language-2nd-Edition/dp/0131103628), 1988.

### 更新记录 ###

2014年12月7日
