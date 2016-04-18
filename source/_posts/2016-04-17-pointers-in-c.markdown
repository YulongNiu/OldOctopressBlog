---
layout: post
title: "C语言指针记录"
date: 2016-04-17 22:42:25 +0800
comments: true
published: false
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

1. 指针加一个整数，该表达式值为同类型指针；

2. 指针减一个整数，该表达式值为同类型指针；

3. 指针与指针相减，该表达式值为整数。

## 3. 指针与数组 ##

1. C语言只有一维数组，其中元素可以是数（整数或浮点数）、字符和指针（字符串、其他类型数组或者其他类型指针）。

2. 数组类型为元素类型。

3. 数组地址为第一个元素地址。可以使用数组名作为指向数组第一个元素的指针，但数组名<span style="color: red">不能</span>被修改，例如不能被重新赋值。因此，假如`a`数组，`a+i`等价与`&a[i]`，`*(a+i)`等价与`a[i]`。

------------------

| 数组类型                               | 初始化声明<sup>1</sup>             | 函数形参声明的指针形式<sup>2</sup>   | 第一个元素指针声明                         |
|----------------------------------------+-------------------------+--------------------------+--------------------------------------------|
| 元素为整数的数组                       | `int a[LEN]`            | `int *`                  | `int *p = &a[0]`或`int *p = a`             |
| 元素为整数数组的数组（“二维数组”）     | `int a[ROWNUM][COLNUM]` | `int (*)[COLNUM]`        | `int *p = &a[0]`或`int (*p)[COLNUM] = a`　 |
| 元素为字符的数组（“字符串”）           | `char a[LEN]`           | `char *`                 | `char *p = &a[0]`或`char *p = a`           |
| 元素为字符串指针的数组（“字符串数组”） | `char *a[LEN]`          | `char **`或`char *[LEN]` | `char **p = &a[0]`或`char **p = a`         | 
  
<sup>1</sup>：初始化声明表示在声明同时初始化的形式，比如`int a[3] = {1, 2, 3}`、`char a[] = 'hello'`或者`char a[] = {"hello", "world!"}`。

<sup>1</sup>：在声明形式参数时，对应的指针类型。

------------------






## 4. 指针与函数 ##







## 5. 注意事项 ##

1. 留意未初始化指针，修改未初始化指针所指向内容是<span style="color: red">危险</span>的。

2. 

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




### 更新记录 ###

201６年４月17日
