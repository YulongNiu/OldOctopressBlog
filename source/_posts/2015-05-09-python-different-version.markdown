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
# 如果命令行下载安装Python太慢，可以将下载的安装包放入~/.pyenv/cache/文件夹中，之后安装
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

如果在使用`pip`{:.language-bash}安装包时，下载速度过慢，可以考虑使用国内源。比如：

* http://pypi.mirrors.ustc.edu.cn/

* http://pypi.douban.com/

使用方法为：

{% codeblock lang:bash Using different pip mirror%}
$ pip install --upgrade numpy -i http://pypi.mirrors.ustc.edu.cn/simple
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

* [使用国内镜像源来加速python pypi包的安装](http://www.leadnt.com/2013/08/%E4%BD%BF%E7%94%A8%E5%9B%BD%E5%86%85%E9%95%9C%E5%83%8F%E6%BA%90%E6%9D%A5%E5%8A%A0%E9%80%9Fpython-pypi%E5%8C%85%E7%9A%84%E5%AE%89%E8%A3%85/) 


### 更新记录 ###

2015年5月23日
