---
layout: post
title: "使用Pyenv控制多个版本Python"
date: 2015-05-09 00:14:10 +0800
comments: true
categories: PPR 
---


同时在一台电脑上使用多个版本的Python，可以通过[pyenv](https://github.com/yyuu/pyenv)控制和管理。

## 1. 安装方法 ##

{% codeblock lang:bash %}
$ git clone git://github.com/yyuu/pyenv.git ~/.pyenv

# 写入路径信息
$ echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
$ echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
$ echo 'eval "$(pyenv init -)"' >> ~/.bashrc

# 重新载入
$ source ~/.bashrc
{% endcodeblock %}

<!--more-->

## 2. 安装多个版本Python和对应包 ##

{% codeblock lang:bash %}
# 查看可提供的Python版本列表
$ pyenv install --list

# 安装Python所依赖包
# yum install readline readline-devel readline-static openssl openssl-devel openssl-static sqlite-devel bzip2-devel bzip2-libs
# 安装其他版本Python
$ pyenv install 3.4.3
$ pyenv rehash

# 查看已安装Python版本
$ pyenv versions

# 全局切换Python版本
$ pyenv global 3.4.3

# 安装对应版本Python包。
# 每次安装包之后，都要执行rehash。
$ pip install ipython
$ pyenv rehash
{% endcodeblock %}


## 3. 更新pyenv ##

{% codeblock lang:bash %}
$ cd ~/.pyenv
$ git pull
{% endcodeblock %}


## 4. 删除特定版本Python ##

{% codeblock lang:bash %}
# 查找特定版本Python文件夹位置，之后直接删除即可。
$ pyenv prefix 3.4.3
$ rm -rf ~/.pyenv/versions/3.4.3
{% endcodeblock %}




### <a id="Ref">参考网址</a> ###

* [Python多版本共存之pyenv](http://seisman.info/python-pyenv.html)

### 更新记录 ###

2015年5月9日
