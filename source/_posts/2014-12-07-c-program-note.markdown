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

<img src="/images/coperator.png" title="image" alt="C operators">

图片取自[参考资料2](#Ref)。


* 赋值操作`=`：从右向左。

* 或操作`||`、且操作`&&`：从左向右，第一个判定成功即终止。



## 3. 语法注意事项 ##

* `#define LOWER 0`定义常量的语句之后，没有分号`;`。

* `EOF`是文档结束的标志（End of File），在`<stdio.h>`{:.language-c}定义为一个整数`-1`，代码如下：

{% codeblock definition of EOF in C lang:c %}
#ifndef EOF
# define EOF (-1)
#endif
{% endcodeblock %}

* 用引号`''`标记的字符串，对应的是ASCII的数字值。

* `for`循环声明需要有主体，如无，在`for`语句后添加`;`作“无效声明（null statement）”。

* 矩阵计数从0开始。

* `i++`{:.language-c}与`++i`{:.language-c}

    * `++i`{:.language-c}马上自增，`i++`{:.language-c}自增则在两个相邻顺序点之间进行。表达式`++i`{:.language-c}值为`i+1`，表达式`i++`{:.language-c}值为`i`。

    * 对于表达式`f(i++)`{:.language-c}，传入的参数值为`i`，但是在函数内部开始执行前，`i`{:.language-c}完成自增。这是因为在**函数的所有参数赋值和函数第一条语句执行之前**有一个顺序点。

    * 在`for`{:.language-c}语句中，使用`for(i = 0; i < 10; i++)`{:.language-c}与`for(i = 0; i < 10; ++i)`{:.language-c}效果一样。

    * 对于现代编译器，`i++`{:.language-c}和`++i`{:.language-c}的执行效率没有区别。所以写代码时，按照自认为最清楚的方式写。


* 副作用

   * 所有赋值运算、`i++`{:.language-c}和`++i`{:.language-c}都有副作用，即改变原始变量的值。

   * 赋值运算的值是赋值操作后，左侧的值；强制转换为左侧值类型；赋值运算左侧必须为“左值（lvalue）”。




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


### <a id="Ref">参考资料</a> ###

* BW Kernighan and DM Ritchie: [The C Programming Language (2nd Edition)](http://www.amazon.com/The-Programming-Language-2nd-Edition/dp/0131103628), 1988.

* KN King: [C Programming: A Modern Approach, 2nd Edition](http://www.amazon.com/Programming-Modern-Approach-2nd-Edition/dp/0393979504), 2008.

* [C各种标准](http://port70.net/~nsz/c/)


### 更新记录 ###

2014年12月7日
