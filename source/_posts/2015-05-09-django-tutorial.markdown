---
layout: post
title: "Django使用介绍"
date: 2015-05-09 13:42:15 +0800
comments: true
published: true
categories: PPR
---

本文是学习[Django Tutorial](https://docs.djangoproject.com/en/1.8/intro/tutorial01/)的记录，目的为了帮助快速浏览和查找Django使用细节。


## 1. 安装Django ##

本文使用[Python 3.3.6](https://www.python.org/)和[PostgreSQL 9.3.6](http://www.postgresql.org/)学习Django。可以参考[“使用Pyevn控制多个版本Python”](http://yulongniu.bionutshell.org/blog/2015/05/09/python-different-version/)和 [“Fedora安装与使用PostgreSQL”](http://yulongniu.bionutshell.org/blog/2015/05/08/install-use-postgresql/)，安装对应版本Python和PostgreSQL。

{% codeblock lang:bash %}
# 安装Django
$ pip install django

# 安装PostgreSQL支持psycopg2 
$ pip install psycopg2 
{% endcodeblock %}

<!--more-->



### <a id="Ref">参考网址</a> ###

* [Django Tutorial](https://docs.djangoproject.com/en/1.8/intro/tutorial01/)

### 更新记录 ###

2015年5月9日
