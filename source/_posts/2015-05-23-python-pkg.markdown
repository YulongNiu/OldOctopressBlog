---
layout: post
title: "Python打包和上传PyPI"
date: 2015-05-23 01:03:26 +0800
comments: true
categories: PPR
---

一个简单的Python包创建和上传PyPI流程。

## 1. 打包 ##

{% codeblock lang:bash %}
## 安装相关工具
pip install --upgrade setuptools wheel twine
{% endcodeblock %}

<!--more-->

当Python包代码完成后，在根目录下创建`setup.py`文件，具体参考[sampleproject](https://github.com/pypa/sampleproject/blob/master/setup.py)。

{% codeblock lang:bash %}
## 进入Python包目录
## 测试
python setup.py test

## 打包
python setup.py sdist bdist_wheel
{% endcodeblock %}

## 2. 上传PyPI ##

[PyPI](https://pypi.org/)注册帐号。

{% codeblock lang:bash %}
## 进入Python包目录
## 上传
twine upload dist/*
{% endcodeblock %}

### <a id="Ref">参考网址</a> ###

* [Packaging Python Projects](https://packaging.python.org/tutorials/packaging-projects/)

* [twine package](https://pypi.org/project/twine/)

### 更新记录 ###

2018年6月23日



