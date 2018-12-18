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
$ docker image ls --digests

# Docker容器列表
$ docker container ls
$ docker container ls --all

# 列出容器
$ docker ps
$ docker ps -al

# 终止所有容器
$ docker stop $(docker ps -aq)

# 删除镜像
$ docker image rm myImg@sha256:xxx

## 删除虚悬镜像
$ docker ps -a | grep "Exited" | awk '{print $1 }'|xargs docker stop
$ docker ps -a | grep "Exited" | awk '{print $1 }'|xargs docker rm
$ docker images|grep none|awk '{print $3 }'|xargs docker rmi
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

## 5. 保存和本地载入镜像 ##

{% codeblock lang:bash %}
# 查看镜像列表
$ docker images

# 保存镜像
$ docker save -o dockerImage.tar /example/dockerImage

# 载入镜像
$ docker load --input dockerImage.tar
{% endcodeblock %}

## 6. 修改镜像和容器储存位置 ##

查看镜像储存位置，例如`/var/lib/docker`

{% codeblock lang:bash %}
$ docker info | grep "Docker Root Dir"
{% endcodeblock %}

移动镜像和容器存储位置

{% codeblock lang:bash %}
$ mv /var/lib/docker /localpath/docker
$ ln -s /localpath/docker /var/lib/docker
{% endcodeblock %}

### <a id="Ref">参考网址</a> ###

* [Docker官方文档](https://docs.docker.com/)

* [Docker — 从入门到实践](https://www.gitbook.com/book/yeasy/docker_practice/details)

* [删除Docker镜像中为none的镜像](https://www.centos.bz/2017/08/docker-delete-none-images/)

### 更新记录 ###

2018年12月16日
