---
layout: post
title: "Fedora安装与使用postgresql"
date: 2015-05-08 18:30:54 +0800
comments: true
categories: Linux
---

## 安装和开启postgresql ##

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









### <a id="Ref">参考网址</a> ###

* [初始化错误](http://www.heatware.net/linux-unix/how-install-postgresql-8-4-centos-5/)



### 更新记录 ###

2015年5月8日
