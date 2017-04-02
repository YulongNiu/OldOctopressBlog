---
layout: post
title: "Linux常用软件集锦"
date: 2010-07-07 22:09:31 -0500
comments: true
categories: linux
---

本文收集了一些好用或者好玩的Linux软件，安装方法使用[Fedora](https://getfedora.org/)系统示例。


## 1. 添加常用源 ##

* 添加[rpmfusion](http://rpmfusion.org/Configuration)的free和nonfree源。

* 添加[FZUG](http://repo.fdzh.org/)源。一个中文Fedora的源和社区，在Fedora下安装搜狗输入法、有道词典、百度云客户端等。

## 2. 辅助安装工具  ##

* fed

[fedy](http://folkswithhats.org/)是一个辅助在Fedora系统下，更容易安装Skypy、Chrome、Steam等的工具。

<!--more-->

## 3. 影音 ##

* VLC

{% codeblock lang:bash %}
# yum install vlc
{% endcodeblock %}

* Mplayer

{% codeblock lang:bash %}
# yum install mplayer-gui xine smplayer
{% endcodeblock %}

* 录音工具Audacity

{% codeblock lang:bash %}
# yum install audacity
{% endcodeblock %}


* 屏幕录制RecordMyDesktop

Linux下有很多录制屏幕的软件，推荐RecordMyDesktop。可以自己选定屏幕区域和大小，也可以把界面缩小到屏幕的下方，变成一个小的按钮，以方便操作。录制的文件为ogv格式，推荐使用VLC打开。

{% codeblock lang:bash %}
# 安装方法
# yum install recordmydesktop

# 使用方法
$ qt-recordMyDesktop
{% endcodeblock %}

## 4. 办公 ##

* 字典GoldenDict

Linux下曾经风靡一时的字典StarDict，现在有了更加先进和方便的接班人[GoldenDict](http://goldendict.org/) 。GoldenDict的主要特点有：

1. 字典库丰富;

2. 直接查询维基百科和其他网络字典;

3. 支持字典分类、发音（字典库包括发音）、光标取词等。


**安装方法：**

参考GoldenDict的[GitHub安装介绍](https://github.com/goldendict/goldendict)，字典发音需要安装Mplayer，Fedora安装可能依赖的库如下：

{% codeblock lang:bash %}
# dnf install bzip2-devel gcc-c++ git hunspell-devel \
libvorbis-devel libXtst-devel phonon-devel \
qt-devel qtwebkit ffmpeg ffmpeg-devel \
lzo lzo-devel eb eb-devel libao libao-devel \
qtwebkit-devel libtiff-devel lzma-devel

# 安装字典
$ git clone git://github.com/goldendict/goldendict.git
$ cd goldendict
$ qmake-qt4 "CONFIG+=zim_support" && make
# make install
{% endcodeblock %}

**添加离线字典**

* 安装dictd-server

{% codeblock lang:bash %}
# yum install dictd-server
{% endcodeblock %}

* dsl文件处理

GoldenDict只能识别dsl格式的字典文件，所以先进行处理。一般得到dsl文件内容如下：

{% raw %}
```
mydict.dsl
mydict.bmp
mydict.ann
mydict.images.rar
mydict.sounds.rar
```
{% endraw %}

> 压缩dsl文件：

{% codeblock lang:bash %}
# dictzip mydict.dsl 
{% endcodeblock %}

> 解压images和sounds的所有文件到同一个文件夹，之后将其[全部压缩到一个文件夹中](http://forum.ubuntu.org.cn/viewtopic.php?f=48&t=316122&start=0)：

{% codeblock lang:bash %}
# find . -name "*" -print | zip -9 ../mydict.dsl.dz.files.zip -@
{% endcodeblock %}

> 最后，得到如下文件，使用GoldenDict载入即可，注意文件的命名都统一为`mydict`：

{% raw %}
```
mydict.dsl.dz
mydict.dsl.dz.files.zip
mydict.ann
mydict.bmp
```
{% endraw %}


以下是两张在GoldenDict中查询维基百科和大英百科的效果图。

<img src="/images/linux_software_goldendict1.png" width="500" height="700" title="image" alt="images">


<img src="/images/linux_software_goldendict2.png" width="500" height="700" title="image" alt="images">


* chm阅读器

{% codeblock lang:bash %}
# yum install kchmviewer
{% endcodeblock %}

* 文档注释工具Xournal

{% codeblock lang:bash %}
# yum install xournal
{% endcodeblock %}


* ePub文件阅读器

{% codeblock lang:bash %}
# yum install fbreader
{% endcodeblock %}

* djvu文件阅读器

{% codeblock lang:bash %}
# yum install djvulibre
{% endcodeblock %}


* Fedora 安装搜狗输入法

参考[fedora 21 安装sogoupinyin 安装搜狗拼音](http://www.8dlive.com/post/199.html)

* unrar解压工具

{% codeblock lang:bash %}
# yum install unrar
{% endcodeblock %}

* 截屏工具KSnapshot

{% codeblock lang:bash %}
# yum install ksnapshot
{% endcodeblock %}

* 图片转换工具[ImageMagick](http://www.imagemagick.org/)

{% codeblock lang:bash %}
# yum install ImageMagick

# 设置转换图片质量，1质量最低，100质量最高
$ convert -quality 100 input.pdf output.jpg

# 设置像素，比如转换的是500px
$ convert -density 500 input.pdf output.jpg
{% endcodeblock %}

* TeX文本编辑器

[TeXstudio](http://texstudio.sourceforge.net/)支持自动补全、代码高亮、错误提示、文档预览、图片表格公式生成、LaTeX/PDFLaTeX/XeLaTeX。甚至还有一个“放大镜”，放大观察生成文档字体和公式细节。TeXstudio仍然保持持续更新的态势，以下是安装方法和一个阅览图。

{% codeblock lang:bash %}
# yum install texstudio
{% endcodeblock %}

<img src="/images/linux_software_texstudio.jpg" width="500" height="700" title="image" alt="images">

* PDF文件编辑工具

[cpdf](http://community.coherentpdf.com/)用于对PDF文件执行分割、合并和压缩等操作。下载预编译版本使用。

{% codeblock lang:bash Example of using cpdf%}
# 选择第一页输出
$ cpdf pdfFile.pdf 1 -o outPdfFile.pdf
{% endcodeblock %}


## 5. 模拟Windows程序 ##

* Wine

Wine可以尽可能模拟Window软件

{% codeblock lang:bash %}
# yum install wine
{% endcodeblock %}

同时可能需要Winetricks辅助

{% codeblock lang:bash %}
# wget http://www.kegel.com/wine/winetricks
# chmod +x winetricks
# mv winetricks /usr/local/bin
# winetricks mfc42 
{% endcodeblock %}

* PlayOnLinux

{% codeblock lang:bash %}
# wget http://rpm.playonlinux.com/playonlinux.repo
# mv playonlinux.repo /etc/yum.repos.d/
# yum update
# yum install playonlinux
{% endcodeblock %}


## 6. 驱动程序 ##

* Fedora安装NVIDIA显卡驱动

安装过程参考[Fedora 20安装Nvidia 显卡驱动](http://blog.csdn.net/icemango1209/article/details/18766353)。需要注意关闭窗口界面，从而进入完全字符界面操作。

{% codeblock lang:bash %}
# 进入字符界面
# init 3

# 进入图形界面
# init 5
{% endcodeblock %}


## 7. Samsung SSD ##

三星固态硬盘支持Linux版本的[Magician DC](http://www.samsung.com/global/business/semiconductor/minisite/SSD/global/html/support/server_downloads.html)。

{% codeblock lang:bash Samsung SSD Magician %}
# 列出SSD
# magician -L

# 列出指定SSD详细信息
# magician -d 0 -I

# 列出指定SSD的S.M.A.R.T信息
# magician -d 0 -S
# 列出制定SSD的温度
# magician -d 0 -S -t

# trim
# magician -d 0 -T
{% endcodeblock %}

## 8. SD卡录入镜像 ##

推荐使用[Etcher](https://etcher.io/)录入镜像至SD卡。


### <a id="Ref">参考网址</a> ###

* [pdf转换jpg/jpeg文件](http://xmodulo.com/convert-pdf-files-to-jpg-format-on-linux.html)

### 更新记录 ###

2017年4月2日

















