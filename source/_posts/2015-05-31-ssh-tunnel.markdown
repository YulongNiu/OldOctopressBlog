---
layout: post
title: "SSH建立端口转发"
date: 2015-05-31 21:14:02 +0800
comments: true
categories: Linux
---

SSH建立端口转发分为两步：

## 1. 连接远程主机 ##

{% codeblock lang:bash %}
$ ssh -NT -D 8080 usrname@host
{% endcodeblock %}

<!--more-->

其中，`-N`表示只连接远程主机，不打开远程shell；`-T`表示不分配TTY；`-D`表示某端口数据都通过SSH传向远程主机；`8080`可以设置为其他端口。


## 2. 添加Chrome浏览器支持 ##

添加[Proxy SwitchySharp插件](https://chrome.google.com/webstore/detail/proxy-switchysharp/dpplabbmogkhghncfbfdeeokoefdjegm?hl=en)，之后在SOCKS host栏目中填入地址`127.0.0.1`，端口`8080`并启用即可。


### 参考网址 ###

* [SSH隧道翻墙的原理和实现](http://www.pchou.info/linux/2015/11/01/ssh-tunnel.html)

* [SSH原理与运用（二）：远程操作与端口转发](http://www.ruanyifeng.com/blog/2011/12/ssh_port_forwarding.html) 

### 更新记录 ###

2017年5月31日
