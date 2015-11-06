---
layout: post
title: "C语言学习记录"
date: 2014-12-07 18:23:16 -0500
comments: true
published: true
categories: c
---

## 1. 基本数据类型 ##

### 1.1 算术类型 ###

强制编辑器处理常量为长整数，`1135L`；强制处理为无符号，`1135U`；混合使用，`1135UL`，`U`和`L`顺序和大小写不重要。强制编辑器处理常量为单精度浮点数，`11.3f`；强制为双精度，`11.3L`或者`11.3l`。**避免无符号和有符号整数混用。**

强制类型转化表达式的一般形式为：`(int)floatNumber`。C语言把`(type)`视为一元运算符，所以其等级高于二元运算符。

类型定义一般形式为：`typedef int Newint;`，注意结尾的`；`。

`sizeof`运算符一般形式为：`sizeof(type)`，比如`sizeof(long int)`计算`int`类型占用多少个字节。也可以应用与常量、变量或者表达式。

* `int`：有符号整数（最左边一位表示符号，`0`为正数，`1`为负数），完整形式为`signed int`（singed可以省略）；`unsinged int`：无符号整数。

* `short`：有符号短整数，完整形式为`signed short int`；`unsigned short`：无符号短整数。

* `long`：有符号长整数，完整形式为`signed long int`；`unsigned long`：无符号长整数。

* `float`：单精度浮点数

* `double`：双精度浮点数

* `long double`：扩展精度浮点数

* `char`：字符，C语言把字符当做小整数处理，字符常量使用**单引号**。

<!--more-->

### 1.2 其他形式 ###

* `arrays`

* `structures`

* `unions`

* `pointers`

* `functions`


## 2. 操作符类型和优先级 ##

<img src="/images/coperator.png" title="image" alt="C operators">

图片取自[参考资料2](#Ref)。

* 操作符有两个性质：**结合方向**和**优先级**。**结合方向**决定操作符的执行对象，比如多个同等操作符；而**优先级**决定操作符的结合方式，通俗来讲即谁和谁结合在一起。但是，C没有规定表达式运算的先后顺序。比如对于二元操作符`+`，`a = i + i++;`{:.language-c}，由编译器决定是`i`还是`i++`先执行。

* C89中操作符`/`和`%`，用于负整数操作，结果由编译器决定。

* 运算符`&&`和`||`（从左向右结合），两侧的两个表达式有运算顺序，先左后右。**有可能右侧表达式没有计算，因此不要在右侧放入有副作用的表达式**。



## 3. 表达式 ##

* 条件表达式`i > 0 ? i : f`{:.language-c}，如果`i`和`f`是整数型和浮点型，即使条件判定为真，表达式的值为浮点型。

* `switch`语句最后一个分支，添加`break;`{:.language-c}语句。防止之后修改程序，需要再添加判断条件时，遗漏`break;`{:.language-c}语句。

* 逗号表达式中，`表达式1, 表达式2`{:.language-c}，`表达式1`先计算之后丢弃其值，之后计算`表达式2`。因此，**`表达式1`必须有副作用**。


## 4. 语法注意事项 ##

* `#define LOWER 0`定义常量的语句之后，没有分号`;`。

* `break;`{:.language-c}只能跳出一层循环。

* `continue;`{:.language-c}。一个有意思的应用场景，用`continue;`语句**条件性**忽略一些语句。

* `EOF`是文档结束的标志（End of File），在`<stdio.h>`{:.language-c}定义为一个整数`-1`，代码如下：

{% codeblock definition of EOF in C lang:c %}
#ifndef EOF
# define EOF (-1)
#endif
{% endcodeblock %}

* 用引号`''`标记的字符串，对应的是ASCII的数字值。

* `for`循环声明需要有主体；如无，在`for`语句后添加`;`作“无效声明（null statement）”。**为了避免误解，`;`单独占一行**，或者使用`{}`代替单独一行的`;`。

* 矩阵计数从0开始。

* 顺序点（sequence point）

    * `&&`{:.language-c}、`||`{:.language-c}和comma operators，左边和右边表达式之间。

    * 三元条件操作符`?:`{:.language-c}，在条件判断表达式与第二（第三）表达式之间。

    * 完整表达式结束，包括赋值、`return`{:.language-c}语句、`if`{:.language-c}/`switch`{:.language-c}/`while`{:.language-c}/`do-while`{:.language-c}条件表达式判断结束和`for`{:.language-c}三个表示式。

    * 函数的所有参数赋值和函数第一条语句执行之前（见后举例）。

    * 变量初始化语句结束，比如`int a = 1;`{:.language-c}。如果多个变量初始化（`,`分割），则在每一个`,`结束处，比如 `int a = 1, b = 2;`{:.language-c}。

* `i++`{:.language-c}与`++i`{:.language-c}

    * 大多数情况用于整数操作。
    
    * `++i`{:.language-c}马上自增，`i++`{:.language-c}自增则在两个相邻顺序点之间进行。表达式`++i`{:.language-c}值为`i+1`，表达式`i++`{:.language-c}值为`i`。

    * 表达式`i++ == 0;`{:.language-c}中，`i`使用原始数值，在表达式结束后自增；相反，表达式`++i == 0;`{:.language-c}中，`i`马上自增后与0比较。

    * 对于表达式`f(i++)`{:.language-c}，传入的参数值为`i`，但是在函数内部开始执行前，`i`{:.language-c}完成自增。这是因为在**函数的所有参数赋值和函数第一条语句执行之前**有一个**顺序点**。

    * 在`for`{:.language-c}语句中，使用`for(i = 0; i < 10; i++)`{:.language-c}与`for(i = 0; i < 10; ++i)`{:.language-c}效果一样。

    * 对于现代编译器，`i++`{:.language-c}和`++i`{:.language-c}的执行效率没有区别。所以写代码时，按照自认为最清楚的方式写。


* 副作用（side effect）

   * 所有赋值运算、`i++`{:.language-c}和`++i`{:.language-c}都有副作用，即改变原始变量的值。

   * 赋值运算的值是赋值操作后右侧的值，并且将其强制转换为左侧值类型；赋值运算左侧必须为“左值（lvalue）”。

   * 在同一个表达式，即访问某个变量，同时又修改这个变量，会造成**“未定义行为（undefined behavior）”**。有副作用的操作，会带来隐晦为定义行为。未定义行为会随着不同的编译器，而产生不同的结果。其危险性不仅在于阻碍跨平台使用，而且也会有程序运行失败或者得到意想不到结果。**建议：不在一个表达式中即访问又修改同一个变量**。一些典型的未定义行为的例子：

{% codeblock lang:c Undefined behavior in C %}
# can not decide whether "++", "=", or "+" is the first
a = i + i++;
i = i++;
a[i] = i++;

# "," is not comma operator
printf("%d %d\n", ++i, i);
{% endcodeblock %}



## 5. 标准库 ##

### 5.1 `#include <stdio.h>`{:.language-c} ###

{% codeblock lang:c %}
printf("%3.0f %6.1f\n", fahr, celsius);
{% endcodeblock %}

`%m.pX`：格式串模板，`m`表示要显示的最小字段宽度，`p`表示精度，`X`表示类型。

对于 `printf()`{:.language-c}函数：

* `%d`：十进制整数；`%1d`一个十进制整数，配合`scanf()`实现单个整数输入操作。

* `%6d`：十进制整数，至少6位宽，右对齐；如果数字位数超过6位，则全部显示。

* `%-6d`：同上，区别左对齐。

* `%-6.3d`：同上，如果数字位数少于3位，左侧加0。比如，`printf(".2d%", 5)`{:.language-c}的输出为`05`。

* `%f`：浮点数，默认小数位为6位。

* `%6f`：浮点数，至少6位宽。

* `%6.0f`：浮点树，至少6位宽，无小数点而且无小数位。

* `%.2f`：浮点数，小数点后两位。

* `%6.2f`：浮点数，至少6位宽，小数点后两位。

* `%e`：科学计数。

* `%g`：科学计数或者浮点数。

* `%o`：无符号八进制。`%x`：无符号六进制。

* `%c`：字符。

* `%s`：字符串。

* `%%`：`%`自身。

* 对于`%d`、`%o`、`%x`，`%ld`等表示长整数，`%hd`等表示短整数。

对于`scanf()`{:.language-c}函数：

* 格式串中出现空白字符（空格、水平或者垂直制表符、换页符、换行符），数量无关紧要，因为可以匹配任意数量（包括0个）空白字符。

* **不要输把空白字符（比如`"%d "`）放入格式串的结尾**。原因：`scanf()`{:.language-c}会挂起，直到出现不能匹配的空白字符。

* `%f`、`%e`、`%g`：可以互换使用；double型，使用`%lf`等；long double型，`%Lf`等。

## 6. C99注意事项 ##

* 在`for`语句中声明的计数变量，只能在for循环中使用。比如 `for(int i = 0; i < 10; i++){...}`{:.language-c}，变量`i`只能在循环内起作用。















### <a id="Ref">参考资料</a> ###

* BW Kernighan and DM Ritchie: [The C Programming Language (2nd Edition)](http://www.amazon.com/The-Programming-Language-2nd-Edition/dp/0131103628), 1988.

* KN King: [C Programming: A Modern Approach, 2nd Edition](http://www.amazon.com/Programming-Modern-Approach-2nd-Edition/dp/0393979504), 2008.

* [C各种标准](http://port70.net/~nsz/c/)


### 更新记录 ###

2014年12月7日
