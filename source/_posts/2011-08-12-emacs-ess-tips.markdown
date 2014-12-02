---
layout: post
title: "Emacs和ESS的使用技巧"
date: 2011-08-12 17:20:57 -0400
comments: true
categories: lisp
---

进入GNU Emacs界面之后，输入`M-X R` 即可进入R界面。在这个过程中，会询问是否在当前运行目录下运行R，可以选择在不同目录下运行。


* `M-x R`：

> 1. 启动R。R运行的buffer因为是在Emacs编辑器下运行，所以称为inferior（Emacs文档中称之为iESS buffer）。
> 
> 2. `C-u M-x R RET --no-save RET`：启动R并且不保存。
> 
> 3. `M-x ess-transcript-clean-region`：清理R界面。
> 
> 4. `C-c C-z`：强制停止运行的R进程。

<!--more-->

* `C-c C-x`：

> 1. 代替`ls()`{:.language-r}函数。
> 
> 2. `C-c C-s`：代替`search()`{:.language-r}函数。
> 
> 3. `C-c C-d`：修改已经建立的对象，非常实用。


* `C-c C-n`：

> 1. 把当前行送到R。
> 
> 2. `C-c C-c`：把当前段送到R。
> 
> 3. `C-c C-b`：把当前整个文件送到R。

* `C-c tab`：自动补全R代码。

* `C-x o`：滚动屏幕。

* `C--`：

> 1. 自动给出R语言中特有的 `<-`。
> 
> 2. `C---`：连按两下--，则可以显示下划线。


* `C-c C-o`：

> 1. 在函数头按下会激活ESS对R代码的注释功能。
> 
> 2. `M-x customize-group RET ess RET`：配置默认模板，如果需要插入空行，回车没用，猛戳空格。
> 
> 3. `C-c C-e C-c`：将代码注释为roxygen的`##'`开头格式，特别是注释example的时候很好用。
> 
> 4. `C-c C-e p`：光标跳转到注释段落开头。
> 
> 5. `C-c C-e n`：光标跳转到注释所在函数段落结尾。
> 
> 6. `M-q`：整理roxygen注释，将多行注释压缩整理。

{% img middle /images/Emacs_ESS_snap.jpg 700 700 'Emacs ESS #1' 'a snap of Emacs ESS' %}


### <a id="Ref">参考网址</a> ###

* [像忍者一样写R包](http://cos.name/2011/05/write-r-packages-like-a-ninja/)

* [google code R stype](http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html)

* [记载ESS的博客](http://joysofprogramming.com/install-emacs-ess-el-fedora-rhel/)

* [ESS幻灯片](http://www.damtp.cam.ac.uk/user/sje30/ess11/ess-slides.pdf)

* [ESS文档](http://ess.r-project.org/Manual/ess.html)

### 更新记录 ###

2014年9月10日
