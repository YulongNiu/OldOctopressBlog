---
layout: post
title: "C语言指针记录"
date: 2016-04-17 22:42:25 +0800
comments: true
published: true
styles: [data-table]
categories: c 
---

## 1. C语言指针基础 ##

C语言的指针设计是一致和优雅的。C语言中“指针（pointer）”就是**地址**（所以不能用普通整数储存地址），“指针变量（pointer variable）”是存储地址的变量。一个指针变量，**只能**指向一个特定类型的变量，比如整数、浮点数、字符或者指针。

{% codeblock lang:c Initiate a pointer %}
int tmp1 = 1, tmp2;

/* "=" does not mean "assignment", it just means "initiating" */
/* p is the address of tmp1, *p is equal to the value of tmp1*/
int *p = &tmp1;

int *q;
q = &tmp2;

/* p points to tmp1, q points to tmp2, now the value of tmp2 is 1*/
*p = *q;

/* p and q now both points to tmp1*/
q = p;
{% endcodeblock %}

<!--more-->

## 2. 指针运算有且只有三种 ##

* 指针加一个整数，该表达式值为同类型指针；

* 指针减一个整数，该表达式值为同类型指针；

* 指针与指针相减，该表达式值为整数。

--------------
指针与`++`和`--`结合的表达式

| 表达式   | 意义                                      |
|----------+-------------------------------------------|
| `*p++`   | 表达式值为p指针指向内容，之后指针自增     |
| `(*p)++` | 表达式值为p指针指向内容，之后指向内容自增 |
| `++*p`   | 表达式值为p指针指向内容自增，指针不变     |
| `*++p`   | 表达式值为p指针自增后指向内容             |

-------------

## 3. 指针与数组 ##

* C语言只有一维数组，其中元素可以是数（整数或浮点数）、字符和指针（字符串、其他类型数组或者其他类型指针）。

* 数组地址为第一个元素地址。可以使用数组名作为指向数组第一个元素的指针，但数组名<span style="color: red">不能</span>被修改，例如不能被重新赋值。因此，假如`a`数组，`a+i`等价与`&a[i]`，`*(a+i)`等价与`a[i]`。

* 对于二维数组`a`，`a`表示指向第一行的指针（即指针指向一维数组），`a[0]`表示指向第一行第一个元素的指针。理解`a[0]`：`a[0]`等价与`*(a + 0)`，表示指针`a`指向的内容，即第一行数组；同时，一维数组名表示指向第一个元素的指针。

* “字符串字面量（string literal）”被作为字符数组储存，类型为`char *`，因此对于字符串变量`char s[] = "abc";`和`char *s = "abc";`都合法。但是，`int a[] = {1, 2, 3};`合法，<s>int *a = {1, 2, 3};</s>非法。

------------------

| 数组类型                               | 初始化声明<sup>1</sup>             | 函数形参声明的指针形式<sup>2</sup>   | 第一个元素指针声明                         |
|----------------------------------------+-------------------------+--------------------------+--------------------------------------------|
| 元素为整数的数组                       | `int a[LEN]`            | `int *`                  | `int *p = &a[0]`或`int *p = a`             |
| 元素为整数数组的数组（“二维数组”）     | `int a[ROWNUM][COLNUM]` | `int (*)[COLNUM]`        | `int *p = &a[0]`或`int (*p)[COLNUM] = a`　 |
| 元素为字符的数组（“字符串”）           | `char a[LEN]`           | `char *`                 | `char *p = &a[0]`或`char *p = a`           |
| 元素为字符串指针的数组（“字符串数组”） | `char *a[LEN]`          | `char **`或`char *[LEN]` | `char **p = &a[0]`或`char **p = a`         | 
  
<sup>1</sup>：初始化声明表示在声明同时初始化的形式，比如`int a[3] = {1, 2, 3}`、`char a[] = 'hello'`或者`char a[] = {"hello", "world!"}`。

<sup>2</sup>：在函数中声明形参时，对应的指针类型。形参可以是完整类型或者元素类型，比如，形参`char *a[LEN]`是完整类型，形参`char **a`是元素类型；再比如，形参`int a[ROWNUM][COLNUM]`是完整类型，形参`int (*a)[COLNUM]`是元素类型；再比如，形参`char a[]`是完整类型，形参`char *a`是元素类型。编译器把数组型的形参视为指针。

------------------

## 4. 指针与函数 ##

* C语言传入函数的都是值（数组被当做指针传入），而且形参对应对象的一个副本。

    * 形参为指针，可以改变指向的内容。
    
    * 形参为数组，传入指针（指向第一个元素地址）副本。因此，即便是数组名，也可以修改，即可以把数组名当做一个指针用。如下代码合法。
    
{% codeblock lang:c %}
void TestFun(int const *a) {
  ...
  ++a;
  ...
}
{% endcodeblock %}

## 5. 注意事项 ##

* 留意未初始化指针，修改未初始化指针所指向内容是<span style="color: red">危险</span>的。字符指针<span style="color: red">必须</span>初始化，比如指向已有字符变量、字符串字面量或动态分配字符串。

* 已有数组名<span style="color: red">不能被</span>重新赋值，<span style="color: red">不能</span>指向其他地址。

* <span style="color: red">不能</span>返回指向局部自动变量的指针，因为局部变量和对应指针在返回时销毁。

### 补充材料 ###

* <a id="pointer_array">Pointers and arrays in C</a>

{% codeblock lang:c Using pointers to operate arrays in c %}
#include <stdio.h>

#define N 5

void PrintVal(int *a);
void PrintArray(int *a, int length);
void PrintString(char *a, int length);
void Print2Array(int colnum, int rownum, int (*a)[colnum]);
void Print2Array2(int colnum, int rownum, int **a);
void PrintStringArray(char *a[], int length);
void PrintStringArray2(char **a, int length);

int main(void)
{
  int testVal = 2;
  PrintVal(&testVal);
  printf("\n");

  int testArray[N] = {2, 3, 5};
  PrintArray(testArray, N);
  printf("\n");

  char testString[N] = "hell";
  PrintString(testString, N);
  printf("\n");

  int test2Array[N][N] = {{1, 2, 3, 4, 5}};
  Print2Array(N, N, test2Array);
  int *test2Array2[N] = {test2Array[0], test2Array[1], test2Array[2], test2Array[3], test2Array[4]};
  Print2Array2(N, N, test2Array2);

  char *testStringArray[N] = {"Hello,", "it", "is", "me", "!"};
  PrintStringArray(testStringArray, N);
  PrintStringArray2(testStringArray, N);

  return 0;
}

void PrintVal(int *a) {
  printf("%3d", *a);
}

void PrintArray(int *a, int length) {
  int *p;
  for (p = a; p < a + length; ++p) {
    printf("%3d", *p);
  }
}

void PrintString(char *a, int length) {
  char *p;
  for (p = a; p < a + length; ++p) {
    printf("%c", *p);
  }
}


void Print2Array(int colnum, int rownum, int (*a)[colnum]) {
  int (*p)[colnum];
  for (p = a; p < a + colnum; ++p) {
    for (int *q = *p; q < *p + rownum; ++q) {
      printf("%3d", *q);
    }
    printf("\n");
  }
}

void Print2Array2(int colnum, int rownum, int **a) {
  int **p;
  for (p = a; p < a + colnum; ++p) {
    for (int *q = *p; q < *p + rownum; ++q) {
      printf("%3d", *q);
    }
    printf("\n");
  }
}

void PrintStringArray(char *a[], int length) {
  char **p;
  for (p = a; p < a + length; ++p) {
    printf("%s\n", *p);
  }
}

void PrintStringArray2(char **a, int length) {
  char **p;
  for (p = a; p < a + length; ++p) {
    printf("%s\n", *p);
  }
}
{% endcodeblock %}

输出结果为：

{% raw %}
```
  2
  2  3  5  0  0
hell
  1  2  3  4  5
  0  0  0  0  0
  0  0  0  0  0
  0  0  0  0  0
  0  0  0  0  0
  1  2  3  4  5
  0  0  0  0  0
  0  0  0  0  0
  0  0  0  0  0
  0  0  0  0  0
Hello,
it
is
me
!
Hello,
it
is
me
!
```
{% endraw %}

### 参考资料 ###

* KN King: [C Programming: A Modern Approach, 2nd Edition](http://www.amazon.com/Programming-Modern-Approach-2nd-Edition/dp/0393979504), 2008.

### 更新记录 ###

201６年４月17日
