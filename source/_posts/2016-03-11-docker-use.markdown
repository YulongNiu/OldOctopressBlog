---
layout: post
title: "Docker使用记录"
date: 2016-03-11 16:53:51 +0800
comments: true
categories: linux
---

关于Docker的安装和使用，有详细的[文档](https://docs.docker.com/)可供参考。本文收集一些有趣和重要的记录。

## 1. 普通用户权限执行Docker ##

创建`docker`用户组并添加普通用户。

{% codeblock lang:bash %}
$ sudo groupadd docker
$ sudo usermod -aG docker myUserName
{% endcodeblock %}

注销后，再次登录即可。

<!--more-->

## 2. Docker命令集锦 ##

{% codeblock lang:bash %}
# 测试Docker
$ docker run hello-world

# Docker镜象列表
$ docker image ls

# Docker容器列表
$ docker container ls
$ docker container ls --all

# 列出容器
$ docker ps
$ docker ps -al
{% endcodeblock %}

## 3. 运行镜象 ##

{% codeblock lang:bash %}
$ docker run -it --rm myDockerImage myCommand
{% endcodeblock %}

## 4. 挂载卷 ##

挂载文件目录至容器，可以挂载多个。

{% codeblock lang:bash %}
$ docker run -it --rm -v /localpath/data:/data -v /localpath/file:/file myDockerImage myCommand
{% endcodeblock %}

### <a id="Ref">参考网址</a> ###

* [Docker官方文档](https://docs.docker.com/)

* [Docker — 从入门到实践](https://www.gitbook.com/book/yeasy/docker_practice/details)

### 更新记录 ###

2018年4月4日
