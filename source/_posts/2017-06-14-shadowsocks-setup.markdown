---
layout: post
title: "Shadowsocks设置推荐"
date: 2017-06-14 23:49:13 +0800
comments: true
categories: Linux
---

## 1. Shadowsocks客户端 ##

Shadowsocks在主流平台上都有[客户端](https://shadowsocks.org/en/download/clients.html)，下载直接使用即可。对于Linux平台，可能还需要在Chrome浏览器中安装[SwitchyOmega](https://chrome.google.com/webstore/detail/proxy-switchyomega/padekgcemlokbadohgkifijomclgjgif?hl=en)插件。

<!--more-->

## 2. 谷歌设置 ##

进入[网址](https://encrypted.google.com/preferences?hl=zh-CN#languages)，设置自己喜欢的语言，之后重启即可。


## 3. 自己配置Shadowsocks ##

## 3.1 软件准备 ##

{% codeblock lang:bash %}
# apt-get install python-setuptools m2crypto supervisor
# apt-get install python-pip
# pip install --upgrade pip
# pip install --upgrade shadowsocks
{% endcodeblock %}

## 3.2 配置文件 ##

文件位置`/etc/shadowsocks.json`，设置模板：

{% raw %}
```
{
    "server":"***.***.***.***",
    "local_address": "127.0.0.1",
    "port_password":{
     "8381":"******",
     "8382":"******"
    },
    "local_port":1080,
    "timeout":600,
    "method":"aes-256-cfb"
}
```
{% endraw %}

## 3.3 启动 ##

{% codeblock lang:bash %}
# ssserver -c /etc/shadowsocks.json -d start
{% endcodeblock %}



### 参考资料 ###

1. [更改谷歌语言偏好](http://nga.178.com/read.php?tid=8798506)

### 更新记录 ###

2017年6月14日
