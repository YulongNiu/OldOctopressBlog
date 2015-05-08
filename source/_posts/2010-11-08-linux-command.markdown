---
layout: post
title: "Linux常用命令集锦"
date: 2010-11-08 18:00:40 -0400
comments: true
categories: Linux
---

## 1. 压缩与解压缩 ##

`.tar.gz`{:.language-bash}格式文件

{% codeblock lang:bash %}
# 解压 
$ tar -zxvf /filePath/filename.tar.gz

# 压缩 
$ tar -zcvf /filePath/filename.tar.gz /filePath/filename

# 解压到特定文件夹  
$ tar -zxvf /filePath/filename.tar.gz -C /filePath/filename 

# 压缩到特定文件夹  
$ tar -zcvf /filePath/filename.tar.gz -C /filePath/filename 

# 压缩所有txt类型文件 
$ gzip *.txt
{% endcodeblock %}

`.tar.bz2`{:.language-bash}格式文件

{% codeblock lang:bash %}
# 解压 
$ tar -jxvf /filePath/filename.taz.bz2

# 压缩
$ tar -jcvf /filePath/filename.tar.bz2 /filePath/filename
{% endcodeblock %}

<!--more-->

`.rar`{:.language-bash}格式文件

需要安装[rar工具](http://www.rarsoft.com/download.html), 下载对应的linux版本，解压，`make`{:.language-bash}即可。

{% codeblock lang:bash %}
# 解压
$ rar e /filePath/filename.rar /filePath/filename

# 压缩
$ rar a -m5 /filePath/filename.rar /filePath/filename
{% endcodeblock %}

`.zip`{:.language-bash}格式文件

需要安装zip和unzip工具

{% codeblock lang:bash %}
# yum install zip
# yum install unzip
{% endcodeblock %}

{% codeblock lang:bash %}
# 解压
$ unzip filename.zip

# 只打印最简短解压信息，并不解压
$ unzip -tq filename.zip

# 解压其中一个文件
$ unzip filename.zip onefile

# 解压到特定目录
$ unzip filename.zip -d /filePath/filename

# 压缩当前目录
$ zip filename *

# 压缩当前目录包括所有子目录
$ zip -r filename *
{% endcodeblock %}

## 2. 文件夹操作 ##

{% codeblock lang:bash %}
# 创建 
$ mkdir

# 删除 
$ rm  

# 删除整个文件夹 
$ rm -rf

# 复制 
$ cp   

# 复制文件夹  
$ cp -r
{% endcodeblock %}

## 3. 目录 ##

{% codeblock lang:bash %}
# 显示 
$ ls

# 列表显示文件和相关权限 
$ ls -l

# 列表显示文件并按照文件名逆序排列
$ ls -rl

# 显示隐藏文件 
$ ls -a      

# 可以配合使用
$ ls -al

# 查看文件夹大小 
$ ls -lhs

# 查看目录树 
$ tree
{% endcodeblock %}


## 4. 运行 *.sh文件 ##

{% codeblock lang:bash %}
# 将sh文件赋予可执行权限
$ chomd 777 filename.sh
$ sh /filePath/filename.sh
{% endcodeblock %}

## 5. 更改文件名 ##

{% codeblock lang:bash %}
$ mv oldfile newfile
{% endcodeblock %}

## 6. md5报文摘要算法 ##

md5（Message-Digest Algorithm 5）报文摘要，可以用来验证网络文件传输的完整性。

{% codeblock lang:bash %}
$ md5sum file
{% endcodeblock %}

## 7. java类型文件 ##

Java Control Panel位置`/usr/java/jdk1.7.0_45/bin/ControlPanel`{:.language-bash}

{% codeblock lang:bash %}
# 执行.jar格式文件
$ java -jar filename.jar

# 执行.jnlp格式文件
$ javaws filename.jnlp
{% endcodeblock %}

## 8. 查看文件 ##

{% codeblock lang:bash %}
# 将文件直接打印在屏幕上
$ cat filename    
{% endcodeblock %}

## 9. 查看自己系统32位还是64位 ##

{% codeblock lang:bash %}
$ uname -a
{% endcodeblock %}

## 10. 查看当前路径 ##

{% codeblock lang:bash %}
$ pwd
{% endcodeblock %}

## 11. 查看某个程序/库的安装路径 ##

{% codeblock lang:bash %}
# 比如查看R的安装位置
$ which R

# 查看某个文件的位置
$ whereis filename
{% endcodeblock %}

## 12. 查看一个命令的文档帮助 ##

{% codeblock lang:bash %}
# 比如man R
$ man commandname
{% endcodeblock %}

## 13. 批量处理文件 ##

{% codeblock lang:bash %}
# 删除满足条件的数据
$ find ./ -iname '*' | xargs rm -rf 

# 移动大数据量的文件.
$ find ./ -name "*.gif" | xargs -i mv {} /filePath/ 
{% endcodeblock %}

## 14. 修改PATH ##

假设我们程序的绝对路径是`/opt/arbtest/arb`{:.language-bash}
两种方法: 

* 直接命令行运行

{% codeblock lang:bash %}
# export PATH=$PATH:/opt/arbtest
{% endcodeblock %}

这种方法当前有效，重启之后就失效了。

{% codeblock lang:bash %}
# 查看修改该好的PATH
$ export
{% endcodeblock %}
* 修改`/etc/profile`{:.language-bash}（系统设置，任何用户都可使用）或者`~/.bashrc`{:.language-bash}（当前用户）文件。
向这两个文件中添加`export PATH=$PATH:/opt/arbtest`{:.language-bash}

{% codeblock lang:bash %}
# 载入修改好的文件 

# 载入root权限profile
# source /etc/profile

# 载入当前用户.bashrc文件
$ source ~/.bashrc

# 查看修改好的路径
$ echo $PATH
{% endcodeblock %}

## 15. 定向输入输出 ##

禁止屏幕输出，即将屏幕输出导入Linux的无底洞 `/dev/null`{:.language-bash}，比如 

{% codeblock lang:bash %}
$ cat myFile > /dev/null
{% endcodeblock %}

此时，导入的是标准屏幕输出（标号为1）stdout。如果要导入标准错误输出（标号为2）stderr，执行

{% codeblock lang:bash %}
$ cat myFile 2 > /dev/null
{% endcodeblock %}

如果将两种输出全部导入`/dev/null`{:.language-bash}，执行

{% codeblock lang:bash %}
# 最后的"2>$1"表示2的操作等同于1
$ cat myFile > /dev/null 2>$1
{% endcodeblock %}

## 16. 更改文件权限和所有者 ##

使用命令`chmod`{:.language-bash}更改文件权限。

* `u`{:.language-bash}：所有者（user）

* `g`{:.language-bash}：群组（group）

* `o`{:.language-bash}：其他人（others）

* `a`{:.language-bash}：所有人（all）

* `r`{:.language-bash}：表示可读（read），对应数值4

* `w`{:.language-bash}：表示可写（write），对应数值2

* `x`{:.language-bash}：表示可执行（excute），对应数值1

* `-`{:.language-bash}：表示什么操作都不行，对应数值0

{% codeblock lang:bash %}
# 比如rw-rw-r--对应664
$ chomd 664 filename
{% endcodeblock %}

使用命令`chown`{:.language-bash}更改文件所有者

## 17.更改文件时间戳 ##

{% codeblock lang:bash %}
# 更新myfile的存取和修改时间；如果myfile不存在，则创建该文件
$ touch filename
{% endcodeblock %}

## 18. 挂载NTFS分区 ##

CentOS需要手动挂在NTFS分区，下载并安装[NTFS-3G](http://www.tuxera.com/community/ntfs-3g-download/)，安装方法：

{% codeblock lang:bash %}
# ./configure 
# make
# make install 
{% endcodeblock %}

查找NTFS分区路径，此处假定为`/dev/sta1`{:.language-bash}：

{% codeblock lang:bash %}
# 查看硬盘分区
# fdisk -l

# 挂载NTFS分区
# mkdir /mnt/NTFStest
# mount -t ntfs-3g /dev/sta1 /mnt/NTFStest

# 卸载分区
# umount /mnt/NTFStest
{% endcodeblock %} 

如果需要自动挂载ntfs分区，首先要查看分区信息，比如uuid

{% codeblock lang:bash %}
# blkid
{% endcodeblock %}

之后修改`/etc/fstab`{:.language-bash}，添加需要挂载的分区

{% codeblock lang:bash %}
UUID=12D345251F34 /media/D ntfs defaults 0 0
{% endcodeblock %}

## 19. yum服务相关 ##

以下命令都可以配合`grep`{:.language-bash}使用

{% codeblock lang:bash %}
# 终止yum安装
# rm -f /var/run/yum.pid

# 查询包
# yum search pkg

# 重新安装包
# yum reinstall pkg

# 升级包
# yum update pkg

# 卸载包
# yum remove pkg

# 查询已安装包信息
# yum info pkg

# 查看仓库包列表
# yum list pkg*

# 查看已安装的包
# yum list installed
{% endcodeblock %}

## 20. 修改配置文件 ##

可以使用多种文本编辑器，最常用的是[Emacs](http://www.gnu.org/software/emacs/)和[Vim](http://www.vim.org/)，使用方法直`emacs`{:.language-bash}或者`vim`{:.language-bash}和文件名即可。

{% codeblock lang:bash %}
# 修改Apache服务器配制文件
$ vim httpd.conf
$ emacs httpd.conf
{% endcodeblock %}

## 21. 用户管理 ##

{% codeblock lang:bash %}
# 查看用户 
$ w 
$ who
# 超过500为后建用户
$ cat /etc/passwd

# 新建用户 
$ useradd usrname

# 新建用户设置密码 
$ passwd usrname newpasswd

# 删除用户 
$ userdel -r usrname

# 查看用户登陆
$ last
$ last usrname

# 查看当前任务 
$ top
{% endcodeblock %}

## 22. deb和rpm包互转 ##

* 第一种方法是使用[alien](http://joeyh.name/code/alien/)

{% codeblock lang:bash %}
# deb转换为rpm
$ alien -r filename.deb

# rpm转换成deb
$ alien -d filename.rpm
{% endcodeblock %}

* 第二种方法是直接使用`apt`{:.language-bash}，非常方便，配置方法如下

{% codeblock lang:bash %}
# yum install apt
# apt-get update
# apt-get pkg
{% endcodeblock %}

## 23. 查看和终止进程 ##

{% codeblock lang:bash %}
# 查看进程树
$ pstree -p

# 看全部进程
$ ps -A

# 强制终止进程
$ kill -9 7473

# 释放内存
$ free -m
{% endcodeblock %}

## 24. 查看网络有监听的端口 ##

{% codeblock lang:bash %}
$ netstat -lntp
{% endcodeblock %}

## 25. 断开SSH终端，程序后台执行 ##

使用`nohup`{:.language-bash}命令

{% codeblock lang:bash %}
$ nohup /filepath/testScript.py 
{% endcodeblock %}

## 26. 查看某一个库文件的位置 ##

{% codeblock lang:bash %}
$ locate libGLU.so
{% endcodeblock %}

## 27. rpm包 ##

{% codeblock lang:bash %}
# 安装rpm包 
# rpm -ivh pkg.rpm

# 更新rpm包
# rpm -Uvh pkg.rpm

# 查看已经安装的包 
# rpm -qa | grep pkg.rpm

# 卸载rpm包 
# rpm -e pkg
{% endcodeblock %}

## 28. Centos 6 升级 python 2.7 ##

参考网址 [1](https://github.com/0xdata/h2o/wiki/Installing-python-2.7-on-centos-6.3.-Follow-this-sequence-exactly-for-centos-machine-only), [2](http://toomuchdata.com/2012/06/25/how-to-install-python-2-7-3-on-centos-6-2/)


## 29. CPU信息 ##

{% codeblock lang:bash %}
# 查询CPU信息
$ lscpu

# 物理CPU个数
$ grep 'physical id' /proc/cpuinfo | sort -u | wc -l

# 核心个数
$ grep 'core id' /proc/cpuinfo | sort -u | wc -l

# 线程个数
$ grep 'processor' /proc/cpuinfo | sort -u | wc -l
$ nproc
{% endcodeblock %} 

查看CPU温度

{% codeblock lang:bash %}
# 安装lm_sensors工具
# yum install lm_sensors
$ sensors
{% endcodeblock %}

## 30. 查询文件和文件夹大小 ##

{% codeblock lang:bash %}
$ du -h myfile
$ du -h filepath
{% endcodeblock %}

## 31. 查看内存情况 ##

{% codeblock lang:bash %}
# 查看内存使用情况
$ free

# 物理内存
# dmidecode -t memory | grep Size
{% endcodeblock %}

## 32. 查看文件头部和尾部 ##

{% codeblock lang:bash %}
# 头部
$ head -5 file

# 尾部
$ tail -7 file
{% endcodeblock %}

## 33. 下载命令##

{% codeblock lang:bash %}
$ wget -c -t 0 -w 30 httplink
{% endcodeblock %}

* `-c`{:.language-bash}：表示接着下载没下载的文件

* `-t`{:.language-bash}：表示尝试连接次数

* `0`：表示不停尝试

* `-w`{:.language-bash}：表示每两次尝试的时间间隔

## 34. 开机启动 ##

{% codeblock lang:bash %}
# 开机启动httpd
# chkconfig httpd on
# 关闭httpd
# chkconfig httpd off
# 开启启动列表
# chkconfig --list
{% endcodeblock %}

## 35. 系统运行日志 ##

{% codeblock lang:bash %}
# cat /var/log/messages
{% endcodeblock %}






### <a id="Ref">参考网址</a> ###

* [大量文件操作](http://blog.csdn.net/dqswuyundong/article/details/5970004)

* [修改PATH](http://www.cnblogs.com/amboyna/archive/2008/03/08/1096024.html)

*  挂载NTFS分区：[1](http://www.ha97.com/2832.html), [2](http://blog.csdn.net/zouyongjin/article/details/6439232)

* rpm下的apt：[1](http://wiki.debian.org.hk/w/Install_APT_in_Fedora_Linux), [2](http://yinbiao.blog.51cto.com/2765456/511542)

* [Linux使用rar](http://www.liukai.cn/in-linux-setup-rar-for-linux/)

* [后台执行程序](http://www.ibm.com/developerworks/cn/linux/l-cn-nohup/)

* [wget使用](http://www.dayuer.com/freebsd-tooltips/wget_help)

* [zip/unzip](http://www.cyberciti.biz/tips/how-can-i-zipping-and-unzipping-files-under-linux.html)

### 更新记录 ###

2015年5月8日
