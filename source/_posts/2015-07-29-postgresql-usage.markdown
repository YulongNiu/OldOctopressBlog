---
layout: post
title: "postgresql使用注释"
date: 2015-07-29 18:24:10 +0800
comments: true
categories: Linux
---

## 1. 查询 ##

{% codeblock lang:psql Query %}
-- 选取特定的一列或者多列
SELECT column1, column2
FROM dataset;

-- 选取全部列
SELECT *
FROM dataset;

-- 去除重复，选取一列或多列中唯一元素
-- 如输入多咧，则去除多列组合后的重复
SELECT DISTINCT column1, column2
FROM dataset;

-- 制定输入列数
-- a是输出行数，b是h输出起始行（第一行计为0）
SELECT column
FROM dataset
LIMIT a OFFSET b;
{% endcodeblock %}

<!--more-->

* 使用`;`{:.language-bash}结束一条SQL语句；

* 返回未排序数据；

## 2. 排序 ##

{% codeblock lang:psql Ordering %}
-- 按照一列或者多列升序排序
-- 可以按照“不查询”的列排序
-- 先按照column2排，再按照column5排
SELECT column1, column2, column3
FROM dataset
ORDER BY column2, column5;

-- 降序
-- DESC只作用于最靠近的唯一一列，即column2，不作用于column5
SELECT column1, column2, column3
FROM dataset
ORDER BY column2 DESC, column5;

--使用查询列相对编号
SELECT column1, column2, column3
FROM dataset
ORDER BY 2, column5;

{% endcodeblock %}

* 使用相对列查询时，不查询的列，比如上述例子的column5，不能用数字代替。

* 使用相对列查询有风险，不采用。

* 升序关键字为`ASC`{:.language-psql}，通常升序为默认。

* `ORDER BY`{:.language-psql}必须在制定列和数据集后出现。


## 3. 筛选 ##








### 参考资料 ###

* 《SQL必知必会（SQL in 10 Minutes, Sams Teach Yourself (4th Edition)）》[豆瓣链接](https://book.douban.com/subject/24250054/)

* [Postgresql Wiki](https://wiki.postgresql.org/wiki/9.1%E7%AC%AC%E5%9B%9B%E7%AB%A0) 



### 更新记录 ###

2015年7月29日
