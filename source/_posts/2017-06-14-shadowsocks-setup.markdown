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

### 3.1 软件准备 ###

参考[shadowsocks网站](https://github.com/shadowsocks/shadowsocks-libev)安装。

{% codeblock lang:bash %}
$ sudo dnf copr enable librehat/shadowsocks
$ sudo dnf update
$ sudo dnf install shadowsocks-libev
{% endcodeblock %}

### 3.2 配置文件 ###

文件位置`/etc/shadowsocks-libenv/config.json`，设置模板：

{% raw %}
```
{
    "server":"0.0.0.0",
    "port_password":{
     "8381":"******",
     "8382":"******"
    },
    "timeout":600,
    "method":"aes-256-cfb"
}
```
{% endraw %}

### 3.3 启动 ###

{% codeblock lang:bash %}
$ sudo ss-manager -c /etc/shadowsocks-libev/config.json --manager-address 127.0.0.1:8000 -u config.json
{% endcodeblock %}

## 4. 使用TCP BBR加速 ##

{% codeblock lang:bash %}
$ sudo echo "net.core.default_qdisc=fq" >> /etc/sysctl.conf
$ sudo echo "net.ipv4.tcp_congestion_control=bbr" >> /etc/sysctl.conf

$ sudo sysctl -p

$ sudo sysctl net.ipv4.tcp_available_congestion_control
$ sudo sysctl net.ipv4.tcp_congestion_control
{% endcodeblock %}

## 5. 全局配置 ##

使用proxychains全局调用Shadowsocks。

首先，安装和配置`proxychains`：

{% codeblock lang:bash %}
$ sudo dnf install -y proxychains-ng

$ sudo echo 'socks5    127.0.0.1    1080' >> /etc/proxychains.conf
{% endcodeblock %}

之后，打开Shadowsocks后，在需要使用的命令行前加入`proxychains4`，例如：

{% codeblock lang:bash %}
$ proxychains4 git push origin master
{% endcodeblock %}

## 6. privoxy全局代理 ##

首先，安装和配置`privoxy`：

{% codeblock lang:bash %}
$ sudo dnf install -y privoxy

# /etc/privoxy/config修改
# listen-address 127.0.0.1:8118 
# forward-socks5t / 127.0.0.1:1080 
{% endcodeblock %}

配置环境变量并启动：

{% codeblock lang:bash %}
$ export http_proxy="127.0.0.1:8118"
$ export https_proxy="127.0.0.1:8118"
$ export ftp_proxy="127.0.0.1:8118"

$ sudo systemctl restart privoxy
{% endcodeblock %}

### 参考资料 ###

1. [更改谷歌语言偏好](http://nga.178.com/read.php?tid=8798506)

2. [通过TCP BBR为ShadowSocks加速](https://dirtysalt.github.io/blogs/boost-shadowsocks-with-tcp-bbr.html) 

3. [shadowsocks-libev多用户](https://github.com/shadowsocks/shadowsocks-libev/issues/1668)

4. [linux下的ss+privoxy代理配置](https://blog.csdn.net/ypbsyy/article/details/81146866)


### 更新记录 ###

2017年12月19日
