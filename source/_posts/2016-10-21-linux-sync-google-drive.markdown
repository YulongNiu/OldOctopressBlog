---
layout: post
title: "命令行使用Google Drive"
date: 2016-10-21 18:24:31 +0800
comments: true
categories: linux
---

本文简单介绍[gdrive](https://github.com/prasmussen/gdrive)，它是一个跨多个平台软件，提供命令行操作Google Drive。同类软件还有[drive](https://github.com/odeke-em/drive)和[overGrive](https://www.thefanclub.co.za/overgrive)（Linux下Google Drive界面客户端）。

基本思路：对每一个上传至Google Drive的文件或文件夹都分配一个Id，所以云端操作需要指定Id。

突出优点：自动同步，比较云端和本地文件后，决定上传、删除或者替换；同步加入版本控制，可以下载和删除不同版本文件。

<!--more-->

## 1. 初始化 ##

[下载](https://github.com/prasmussen/gdrive)系统对应版本，执行`gdrive about`，根据提示设置。

## 2. 查找 ##

[查找规则](https://developers.google.com/drive/v3/web/search-parameters)

{% codeblock lang:bash %}
# 查找所有文件夹，不包括垃圾箱，所有者为自己
$ gdrive list --query "mimeType = 'application/vnd.google-apps.folder' and trashed = false and 'me' in owners"

# 加入上层目录Id
gdrive list --query "mimeType = 'application/vnd.google-apps.folder' and trashed = false and 'me' in owners and 'YUlPWWdLcy1mX2c' in parents"

# 查询信息
gdrive info YUlPWWdLcy1mX2c
{% endcodeblock %}


## 3. 建立和删除 ##

{% codeblock lang:bash %}
# 建立文件夹
$ gdrive mkdir newFolder

# 建立下一层文件夹，-p指定上层目录Id
$ gdrive mkdir -p M1h4M1dGYUhpSFE newFolder

# 删除文件（文件Id，非名称）
$ gdrive delete 0BzTeuubJesi

# 删除文件夹（文件Id，非名称）
$ gdrive delete -r 0BzTeuubJesi
{% endcodeblock %}


## 5. 同步 ##
{% codeblock lang:bash %}
# gdrive sync会标记同步文件，因此不要在同步文件夹中使用gdrive upload或者网页上传文件。未标记文件会被忽略。
# 同步列表
$ gdrive sync list

# 同步列表内容
$ gdrive sync content VUxydm5iMnM5LWs

# 上传
$ gdrive sync upload myLocaldir 0BzTeuubJesi

# 下载
$ gdrive sync download 0BzTeuubJesi myLocaldir

# 查询所有版本
gdrive revision list YUlPWWdLcy1mX2c

# 下载某一版本，最后两个Id分别为文件Id和版本Id
gdrive revision download YUlPWWdLcy1mX2c Y3JBWEJ5a0gwZndlR3hzWlZubFlUMWFnaHVnPQ

# 删除某一版本
gdrive revision delete YUlPWWdLcy1mX2c Y3JBWEJ5a0gwZndlR3hzWlZubFlUMWFnaHVnPQ
{% endcodeblock %}

### 更新记录 ###

2016年10月21日




