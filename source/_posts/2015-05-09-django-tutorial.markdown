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

本文使用[Python 3.3.6](https://www.python.org/)和[PostgreSQL 9.3.6](http://www.postgresql.org/)学习Django。可以参考[“使用Pyenv控制多个版本Python”](http://yulongniu.bionutshell.org/blog/2015/05/09/python-different-version/)和 [“Fedora安装与使用PostgreSQL”](http://yulongniu.bionutshell.org/blog/2015/05/08/install-use-postgresql/)，安装对应版本Python和PostgreSQL。

{% codeblock lang:bash %}
# 安装Django
$ pip install django

# 安装PostgreSQL支持psycopg2 
$ pip install psycopg2 
{% endcodeblock %}

<!--more-->

{% codeblock lang:bash %}
# 检查Django版本
$ python -c "import django; print(django.get_version())"
{% endcodeblock %}


## 2. 创建项目 ##

{% codeblock lang:bash %}
# 创建名为mysite的项目
$ django-admin startproject mysite
{% endcodeblock %}

之后，配制数据库。修改`mysite/settings.py`{:.language-bash}对应位置。

{% codeblock lang:python mysite/settings.py%}
# Database
# https://docs.djangoproject.com/en/1.8/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'testdb',
        'USER': 'userName',
        'PASSWORD': 'passWord',
        'HOST': '/var/run/postgresql',
        'PORT': '5432',
    }
}
{% endcodeblock %}

{% codeblock lang:bash %}
# 配制好数据库后，链接数据库
$ python manage.py migrate
{% endcodeblock %}

开启Django测试网络服务器

{% codeblock lang:bash %}
$ python manage.py runserver
{% endcodeblock %}

## 3. 创建app ##

1. 创建新的app，比如`polls`{:.language-bash}

{% codeblock lang:bash %}
$ python manage.py startapp polls
{% endcodeblock %}

2. 修改app

app文件位置`polls/models.py`{:.language-bash}，之后在项目配制文件`mysite/settings.py`{:.language-bash}中添加app，最后添加app。

{% codeblock lang:bash %}
# 添加app
$ python manage.py makemigrations polls

# 如果有数据库操作，可以打印具体的数据库操作脚本
$ python manage.py sqlmigrate polls 0001

# 也可以检查
$ python manage.py check
{% endcodeblock %}

3. 链接app与数据库

{% codeblock lang:bash %}
$ python manage.py migrate
{% endcodeblock %}













### <a id="Ref">参考网址</a> ###

* [Django Tutorial](https://docs.djangoproject.com/en/1.8/intro/tutorial01/)

### 更新记录 ###

2015年5月9日
