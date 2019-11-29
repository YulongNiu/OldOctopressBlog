---
layout: post
title: "探索Java基本类型的二进制表示"
date: 2017-05-16 20:17:46 +0800
comments: true
styles: [data-table]
categories: Java
---

## 1. 基本数据类型 ##

Java有[八种](https://docs.oracle.com/javase/specs/jls/se8/html/jls-4.html#jls-4.2)基本数据类型（primitive data type），分别是`boolean`、`char`、`byte`、`short`、`int`、`long`、`float`和`double`，所占用的比特数如下表所示。

<!--more-->

| Data type | Bit                   | Range                       |
|-----------|-----------------------|-----------------------------|
| `boolean` | not precisely defined | `True` or `false`           |
| `char`    | unsigned 16           | $2^{16}$                    |
| `byte`    | signed 8              | $-2^{7}$ ~ $2^7-1$          |
| `short`   | signed 16             | $-2^{15}$ ~ $2^{15}-1$      |
| `int`     | signed 32             | $-2^{31}$ ~ $2^{31}-1$      |
| `long`    | signed 64             | $-2^{63}$ ~ $2^{63}-1$      |
| `float`   | signed 32             | $\pm (1-2^{-24}) * 2^{128}$ 精度$2^{-126}$   |
| `double`  | signed 64             | $\pm (1-2^{-53}) * 2^{1024}$ 精度$2^{-1022}$ |


## 2. 整数二进制表示和范围 ##

以`byte`类型为例，每个整数占用了8个比特，最左边比特位表示正负（0为正，1为负）。例如，$1$表示为$0000 0001_2$。因此，能表示的最大正整数为$0111 1111_2$，即

$$
2^0 + 2^1 + 2^2 + \cdots + 2^6 = 2^7 - 1
$$

由于要满足$-1 + 1 = 0$，所以$-1$表示为$1111 1111_2$，$-2$表示为$1111 1110_2$。类似于正整数，可以得到能表示的最大负整数为$1000 0001_2$，即$-(2^7 - 1)$。

但是，这里出现一个问题：出现了两个0，一个是$+0$（$0000 0000_2$），一个是$-0$（$1000 0000_2$）。因此，规定$1000 0000_2$为最大负整数$-2^7$。

## 3. 浮点数二进制表示和范围 ##

Java使用[IEEE 754](https://en.wikipedia.org/wiki/IEEE_floating_point)标准表示浮点数，其二进制表示分为三个部分：1. 最左边比特位表示正负（0为正，1为负）；2. 指数（`float`有8个比特位，`double`有11比特位）；3. 尾数（`float`有23个比特位，`double`有52比特位）。同时，由于指数需要区分正负，所以`float`指数位转为十进制后需要减去$2^7-1$，而`double`需要减去$2^{10}-1$。

以`float`为例，最大的二进制表示为$0 11111110 11111111111111111111111_2$：

* 最右比特位是0；

* 指数为$1111 1110_2$，即$127$；

* 尾数全为1，即：

$$
2^0 + 2^{-1} + 2^{-2} + \cdots + 2^{-23} = 2-2^{-23}
$$

该数为$(1-2^{-24}) * 2^{128}$。

同理，最小正数为$0 00000001 00000000000000000000000_2$，即$2^{-126}$。

$0 00000000 00000000000000000000000_2$是$+0$，而$1 00000000 00000000000000000000000_2$是$-0$。

$0 11111111 00000000000000000000000_2$是$\infty$，而$1 11111111 00000000000000000000000_2$是$-\infty$。



## 4. 关于2的n次方的有趣事实 ##

| n次方 |             十进制 |
|-------|--------------------|
|     1 |                  2 |
|     2 |                  4 |
|     3 |                  8 |
|     4 |                 16 |
|     5 |                 32 |
|     6 |                 64 |
|     7 |                128 |
|     8 |                256 |
|     9 |                512 |
|    10 |         1024（千） |
|    20 |    1048576（百万） |
|    30 | 1073741824（十亿） |


### 参考资料 ###

1. [The Java Language Specification, Java SE 8 Edition](https://docs.oracle.com/javase/specs/index.html) 

2. [Introduction to Programming in Java](http://introcs.cs.princeton.edu/java/home/)

3. [Wiki Two's complement](https://en.wikipedia.org/wiki/Two's_complement) 

4. [Java Primitive Data Types. Size, Range and Default Value of Basic Data Types](http://cs-fundamentals.com/java-programming/java-primitive-data-types.php) 



### 更新记录 ###

2017年5月16日
