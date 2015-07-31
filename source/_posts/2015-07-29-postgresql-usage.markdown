---
layout: post
title: "postgresql使用指南"
date: 2015-07-29 18:24:10 +0800
comments: true
styles: [data-table]
categories: Linux
---

## 1. 查询 ##

{% codeblock lang:psql Query %}
-- 选取特定的一列或者多列
SELECT column1, column2
FROM dataset

-- 选取全部列
SELECT *
FROM dataset

-- 去除重复，选取一列或多列中唯一元素
-- 如输入多列，则去除多列组合后的重复
SELECT DISTINCT column1, column2
FROM dataset

-- 指定输出列数
-- a是输出行数，b是输出起始行（第一行计为0）
SELECT column
FROM dataset
LIMIT a OFFSET b
{% endcodeblock %}

<!--more-->

* 使用`;`{:.language-bash}结束一条SQL语句，本文中省略；

* 返回未排序数据；

## 2. 排序 ##

{% codeblock lang:psql Ordering %}
-- 按照一列或者多列升序排序
-- 可以按照“不查询”的列排序
-- 先按照column2排，再按照column5排
SELECT column1, column2, column3
FROM dataset
ORDER BY column2, column5

-- 降序
-- DESC只作用于最靠近的唯一一列，即column2，不作用于column5
SELECT column1, column2, column3
FROM dataset
ORDER BY column2 DESC, column5

--使用查询列相对编号
SELECT column1, column2, column3
FROM dataset
ORDER BY 2, column5
{% endcodeblock %}

* 使用相对列查询时，不查询的列，比如上述例子的column5，不能用数字代替。

* 使用相对列查询有风险，不采用。

* 升序关键字为`ASC`{:.language-psql}，通常升序为默认。

* `ORDER BY`{:.language-psql}必须在制定列和数据集后出现。


## 3. 筛选 ##

{% codeblock lang:psql Filter %}
-- ORDER BY语句需要在WHERE之后
-- AND表示“且”， OR表示“或”，第一个条件满足即终止
-- AND和OR可有任意多个
SELECT column1, column2
FROM dataset 
WHERE (column3 = a OR column3 = b) AND column1 = c
ORDER BY column3

-- 多个OR语句使用IN代替，比如上述例子为：
SELECT column1, column2
FROM dataset 
WHERE column3 IN (a, b) AND column1 = c
ORDER BY column3

-- NOT与IN连用
SELECT column1, column2
FROM dataset 
WHERE column3 NOT IN (a, b)

-- BETWEEN a AND b语句，a必须小于等于b
-- 如a等于b，则相当于筛选与a（或者b）相等数值
SELECT column1, column2
FROM dataset
WHERE column1 BETWEEN a AND b

--筛选NULL值使用IS NULL
--筛选非NULL值使用IS NOT NULL
SELECT column1, column2
FROM dataset
WHERE column1 IS NULL
{% endcodeblock %}

* [Postgresql支持的比较符](http://www.postgresql.org/docs/9.4/static/functions-comparison.html#FUNCTIONS-COMPARISON-TABLE)：`<`{:.language-psql}、`<=`{:.language-psql}、`>`{:.language-psql}、`>=`{:.language-psql}、`=`{:.language-psql}和`!=`{:.language-psql}（“不等于”也可以表示为`<>`{:.language-psql}）。

* 筛选字符串条件，需要对筛选串加引号，比如例子中`a`{:.language-psql}为`"testStr"`{:.language-psql}。

* 合理使用括号，强制规定`AND`{:.language-psql}和`OR`{:.language-psql}先后顺序。

* `IN`{:.language-psql}语句执行效率高，并且可以嵌套多层`SELECT`{:.language-psql}语句（每个`SELECT`{:.language-psql}只返回一列数据）。

* 尽量在数据库查询过程，而非自己后续手写，完成数据筛选，因为：1. SQL数据库操作通常比自己手写效率高；2. 便于后续扩展。


## 4. 模式匹配 ##

{% codeblock lang:psql Pattern Matching %}
-- LIKE和NOT LIKE支持对字符串的模式匹配
-- “_”匹配单一一个字符
-- “%”匹配0个或多个字符
SELECT column1, column2
FROM dataset
WHERE column1 LIKE 'F_y%'

-- SIMILAR TO和NOT SIMILAR TO支持正则匹配
-- 正则匹配中，仍然使用“_”和“%”
SELECT column1, column2
FROM dataset
WHERE column1 SIMILAR TO '[^JM]%'
{% endcodeblock %}

* Postgresql支持`ILIKE`{:.language-psql}和`NOT ILIKE`{:.language-psql}忽略大小写敏感搜索，这不是标准SQL语法。

* 注意数据库自动补充的空格，比如 `'F_y'`{:.language-psql}只能匹配“F开头-间隔一个字符-y结尾”的字符串，如果字符串后跟有空格，则不能匹配。

* [Postgresql支持的通配符](http://www.postgresql.org/docs/9.4/static/functions-matching.html)

* 模式匹配效率不高，尽量后置，不要过度使用。


## 5. 函数 ##

{% codeblock lang:psql Function %}
-- 字符串连接
SELECT RTRIM(column1) || ' (' || RTRIM(column1) || ')' AS newName
FROM dataset
ORDER BY column1

--数学计算
SELECT column1, column2 * column3 AS newName
FROM dataset
{% endcodeblock %}


* Postgresql支持的函数：


| 函数名                     | 意义               |
|----------------------------+--------------------|
| `RTRIM()`{:.language-bash} | 删除字符串左侧空格 |
| `LTRIM()`{:.language-bash} | 删除字符串右侧空格 |
| `TRIM()`{:.language-bash}  | 删除双侧空格       |


* 使用`AS`{:.language-psql}即使命名新列。









## 使用建议 ##

1. 数据集的名字为一个单词，比如`priceCustom`{:.language-psql}而不是一个字符串 `'price custom'`{:.language-psql}。同样，命名别名（使用`AS`{:.language-psql}）也如此。







### 参考资料 ###

* 《SQL必知必会（SQL in 10 Minutes, Sams Teach Yourself (4th Edition)）》[豆瓣链接](https://book.douban.com/subject/24250054/)

* [Postgresql Wiki](https://wiki.postgresql.org/wiki/9.1%E7%AC%AC%E5%9B%9B%E7%AB%A0) 



### 更新记录 ###

2015年7月29日
