---
layout: post
title: "Emacs和ESS的使用技巧"
date: 2011-08-12 17:20:57 -0400
comments: true
categories: lisp
---

## 1. 安装ESS ##

有两种方法可以安装，一种是直接使用系统自带的包安装系统，比如`yum`：

{% codeblock lang:bash %}
# yum install emacs-ess
{% endcodeblock %}

但是，有时可能不是ESS最新版本。所以，推荐第二种方法，使用Emacs自带的包系统，方便、更新及时，设置方法[参考](http://yulongniu.bionutshell.org/blog/2012/06/24/emacs-extend-skills/)，安装`ess`包。


## 2. 使用ESS ##

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


* `C-c C-o C-o`：

> 1. 在函数头按下会激活ESS对R代码的注释功能。
> 
> 2. `M-x customize-group RET ess RET`：配置默认模板，如果需要插入空行，回车没用，猛戳空格。
> 
> 3. `C-c C-o C-c`：将代码注释为roxygen的`##'`开头格式，特别是注释example的时候很好用。
> 
> 4. `C-c C-o p`：光标跳转到注释段落开头。
> 
> 5. `C-c C-o n`：光标跳转到注释所在函数段落结尾。
> 
> 6. `M-q`：整理roxygen注释，将多行注释压缩整理。

* `/ssh:user@host: RET`:

> 远程连接服务器R。

{% img middle /images/Emacs_ESS_snap.jpg 900 900 'Emacs ESS #1' 'a snap of Emacs ESS' %}


## 3. 高级设置 ##

结合Emacs的自动补全、补全括号、颜色显示等功能，设置更加强大的ESS编辑环境。将[附加环境设置](#c-mode-config)添加到Emacs设置文件，比如`~/.emacs`{:.language-bash}。所需要的包使用`M-x list-package`{:.language-emacs-lisp}安装。

* 自动补全：`auto-complete`包

<img src="/images/emacs_ess_tipes_autocomplete.png" width="900" height="900" title="image" alt="images">

* 代码折叠：`hs-minor-mode`（系统自带）

* 括号补全：`smartparens`

* 括号颜色：`rainbow-delimiters`

* 颜色显示：`rainbow-mode`

<img src="/images/emacs_ess_tipes_parent.png" width="500" height="500" title="image" alt="images">


## 附加内容 ##

<a id="c-mode-config">ESS的Emacs编程环境设置</a>

{% codeblock Emacs configure lang:emacs-lisp %}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; auto-completion
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(require 'auto-complete)
(require 'auto-complete-config)
(ac-config-default)

;;;;;;;;;;;;;;;;;;;;;;;
; ESS set
;;;;;;;;;;;;;;;;;;;;;;;
(require 'ess-site)
(setq ess-use-auto-complete t)

;;;;;;;;;;;;;;;;;;;;;;;;
;open hs-mode
;;;;;;;;;;;;;;;;;;;;;;;
(add-hook 'ess-mode-hook 'hs-minor-mode)

;;;;;;;;;;;;;;;;;;;;;
;smartparens
;;;;;;;;;;;;;;;;;;;;;
(require 'smartparens-config)
(show-smartparens-global-mode t)
(smartparens-global-mode t)

;;;;;;;;;;;;;;;;;;;;;
;rainbow mode 
;;;;;;;;;;;;;;;;;;;;;
(require 'rainbow-mode)
(dolist (hook '(ess-mode-hook inferior-ess-mode-hook))
(add-hook hook 'rainbow-turn-on))

;;;;;;;;;;;;;;
;rainbow-delimiters
;;;;;;;;;;;;;;
(require 'rainbow-delimiters)
(add-hook 'prog-mode-hook 'rainbow-delimiters-mode)
(add-hook 'ess-mode-hook 'rainbow-delimiters-mode)
{% endcodeblock %}


### <a id="Ref">参考资料</a> ###

* [像忍者一样写R包](http://cos.name/2011/05/write-r-packages-like-a-ninja/)

* [google code R stype](http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html)

* [记载ESS的博客](http://joysofprogramming.com/install-emacs-ess-el-fedora-rhel/)

* [ESS幻灯片](http://www.damtp.cam.ac.uk/user/sje30/ess11/ess-slides.pdf)

* [ESS文档](http://ess.r-project.org/Manual/ess.html)

* [YGC auto-complete in ESS](http://ygc.name/2014/12/07/auto-complete-in-ess/) 


### 更新记录 ###

2017年1月6日

