---
layout: post
title: "两个git的rebase命令应用"
date: 2016-11-11 16:24:39 +0800
comments: true
categories: linux
---

介绍两个`git rebase`的应用场景，一个是合并commits记录，另一个是贡献代码。

## 1. 合并commit记录 ##

假定有多条commits，按照离当前时间从近至远依次为：`c1`、`b2`、`b1`和`a1`。希望合并`c1`、`b2`和`b1`，即只保留`c1`和`a1`。操作流程如下：

<!--more-->

* 打开交互式rebase

{% codeblock lang:bash interactive rebase %}
$ git rebase -i sha1id-of-a1
{% endcodeblock %}

* 标记合并commits

弹出的文本编辑器初始可能为：

{% raw %}
```
pick sha1id-of-c1 c1
pick sha1id-of-b2 b2
pick sha1id-of-b1 b1
...
```
{% endraw %}

修改为：

{% raw %}
```
pick sha1id-of-c1 c1
squash sha1id-of-b2 b2
squash sha1id-of-b1 b1
...
```
{% endraw %}

* 记录合并commits

在弹出的文本编辑器中标记和注释commits

* 提交远程

由于本地和远程记录不一致，需要强制合并。

{% codeblock lang:bash force push %}
$ git push -f origin mybranch
{% endcodeblock %}

## 2. 贡献代码 ##

* Fork项目

Fork在GitHub对应的项目（famous/project，master分支）至自己账户（my/project），克隆至本地并添加远程地址：

{% codeblock lang:bash clone project %}
$ git clone https://github.com/my/project.git
$ git remote add upstream https://github.com/famous/project.git
{% endcodeblock %}

* 建立分支并提交修改

{% codeblock lang:bash commit %}
$ git checkout -b devbranch
$ git commit -a -m 'these commits'
$ git push origin devbranch
{% endcodeblock %}

* 获取和合并最新远程修改

{% codeblock lang:bash commit %}
$ git checkout master
$ git pull upstream master
$ git checkout devbranch

$ git rebase master
## compare and merge the latest commits
$ git rebase --continue

## force push
$ git push -f origin devbranch
{% endcodeblock %}

* 发起pull request。


### 参考资料 ###

1. [合并分支](http://itspg.logdown.com/posts/1731-git-squash-master-commits) 

2. [贡献代码](https://segmentfault.com/a/1190000000736629)


### 更新记录 ###

2016年11月11日


