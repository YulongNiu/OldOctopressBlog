---
layout: post
title: "Fedora安装与使用PostgreSQL"
date: 2015-05-08 18:30:54 +0800
comments: true
categories: Linux
---

## 1. 安装和开启postgresql ##

Fedora直接使用`yum`{:.language-bash}安装相关组件。

{% codeblock lang:bash %}
# 安装PostgreSQL
# yum install postgresql postgresql-server postgresql-contrib postgresql-devel pgadmin3
{% endcodeblock %}

<!--more-->


初始化（只需要执行一次）。如果出现类似`Data directory is not empty!`{:.language-bash}报错，可以尝试执行删除操作`rm -f -r /usr/local/pgsql/data`

{% codeblock lang:bash %}
# 初始化
# postgresql-setup initdb

# 开启服务
# service postgresql start

# 查询工作状态
# service postgresql status

# 关闭服务
# service postgresql stop

# 开机启动
# chkconfig postgresql on
{% endcodeblock %}



## 2. 创建用户和数据库 ##

为了方便使用，可以在PostgreSQL中创建一个与当前系统用户名相同的用户，比如目前系统登陆用户为Kitty。

{% codeblock lang:bash %}
# 进入home路径
$ cd \home

# 进入默认postgres用户，需要输入root密码。
# 之后系统命令提示符会变成类似“bash-4.3$”。
$ su postgres

# 创建用户
# -P：创建用户后立即创建密码
# -s：用户为superuser
# -e：打印消息
bash-4.3$ createuser -P -s -e Kitty

# 创建数据库
# -O：制定用户
bash-4.3$ createdb -O Kitty testdb

# 退出
bash-4.3$ exit

# 再创建新的数据库时，不需要进入postgres。
# 直接在当前登录用户下创建即可。
# 由于PostgreSQL用户与系统登录用户相同，不需要指定PostgreSQL用户
$ creatdb testdb2

# 删除数据库
$ dropdb testdb2
{% endcodeblock %}

创建完用户和对应数据库后，可以登录数据库控制台。登录后，系统命令提示符会变成类似“testdb=#”。

{% codeblock lang:bash %}

# 登录testdb数据库控制台
$ psql testdb

# 完整登录命令
# -U：用户名
# -d：数据库
# -h：host，默认为local socket
# -p：端口
$ psql -U Yulong -d testdb -h /var/run/postgresql -p 5432

# 查看PostgreSQL配制文件路径
testdb=# SHOW config_file;

# 查看所有用户
testdb=# \du

# 查看所有数据库
testdb=# \l

# 退出
testdb=# \q
{% endcodeblock %}









### <a id="Ref">参考网址</a> ###

* [PostgreSQL 9.4 Manuals](http://www.postgresql.org/docs/9.4/interactive/index.html)

* [PostgreSQL新手入门](http://www.ruanyifeng.com/blog/2013/12/getting_started_with_postgresql.html)

* [初始化错误](http://www.heatware.net/linux-unix/how-install-postgresql-8-4-centos-5/)



### 更新记录 ###

2015年5月8日
