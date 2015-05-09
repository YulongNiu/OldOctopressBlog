---
layout: post
title: "Emacs高级使用技巧"
date: 2012-06-24 19:01:36 -0500
comments: true
categories: lisp
---

## 1. Emacs配置文件位置 ##

在Fedora系统下，Emacs的配置文件位置是`~/.emacs`{:.language-bash}。在Emacs配置文件中添加内容后，使用`M-x eval-buffer`{:.language-emacs-lisp}，使当前配置生效。编译代码以加快加载速度，`M-x byte-compile-file`{:.language-eamcs-lisp}。

使用Eamcs解释器`M-x ielm`{:.language-emacs-lisp}。


## 2. Emacs自带的包管理系统 ##

在Emacs24之后，自带包管理系统，使用`M-x eval-buffer`{:.language-bash}进入。可以添加[MELPA源](http://melpa.org/)、[marmalade源](https://marmalade-repo.org/)、[GNU源](http://elpa.gnu.org/packages/)。

{% codeblock lang:emacs-lisp Add repositories of Emacs packages%}
(require 'package)
(add-to-list 'package-archives
	     '("marmalade" . "http://marmalade-repo.org/packages/"))
(add-to-list 'package-archives
	     '("melpa" . "http://melpa.milkbox.net/packages/")
(add-to-list 'package-archives
	     '("gnu" . "http://elpa.gnu.org/packages/")))
(package-initialize)
{% endcodeblock %}

<!--more-->

这样设置之后，就不需要类似`(add-to-list 'load-path "~/.emacs.d/elpa/popup-20140207.1702/")`{:.language-emacs-lisp}的语句了，因为Emacs会自动识别安装的包。但是，如果需要对某个包进行进一步设置，需要加上`(require 'popup)`{:.language-bash}之类的语句。

## 3. 显示行号 ##

使用`M-x linum-mode`{:.language-eamcs-lisp}添加行号。如果需要永久显示，在Emacs配置文件中添加一下内容。

{% codeblock lang:emacs-lisp Show line number%}
;;;;;;;;;;;;;;;;;;;;;;;;;
;open linum mode
;;;;;;;;;;;;;;;;;;;;;;;;;
(setq linum-format "%4d \u2502")
(add-hook 'prog-mode-hook 'linum-mode)
(add-hook 'ess-mode-hook 'linum-mode)
{% endcodeblock %}

## 4. 进入Shell ##

三种方法：

* `M-x shell`{:.language-eamcs-lisp}

* `M-x ansi-term`{:.language-eamcs-lisp}

* `M-x eshell`{:.language-eamcs-lisp}

## 5. root权限 ##

`C-x C-f`{:.language-eamcs-lisp} 之后输入root密码`/su:root@usrname password`{:.language-eamcs-lisp}

## 6. 移动整体代码块 ##

选中代码块后：
向左移动2个字符：`C-u -2 C-x TAB`{:.language-eamcs-lisp}
向右移动4个字符：`C-u 4 C-x TAB`{:.language-eamcs-lisp}

## 7. 添加彩虹猫 ##

添加`nyan-mode`包，之后在Emacs配置文档中写入：

{% codeblock lang:emacs-lisp Add nyan in modeline%}
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;nyan-mode
;;;;;;;;;;;;;;;;;;;;;;;;;;
(nyan-mode t)
{% endcodeblock %}


## 8. 自定义YASnippet模板 ##

使用`M-x yas-new-snippet`{:.language-emacs-lisp}，打开一个模板。比如，添加Octopress的语言高亮模板

{% raw %}
```
# -*- mode: snippet; require-final-newline: nil -*-
# name: Highlight Language
# key: hl
# binding: direct-keybinding
# contributor: Yulong Niu <yulong.niu@aol.com>
# --

{% codeblock lang:${1:bash} %}
$0
{% endcodeblock %}
```
{% endraw %}

其中 `${1:bash}`{:.languag-emacs-lisp}表示光标跳转位置，编号从1开始。`$0`{:.language-emacs-lisp}表示光标最后停留位置。如果不需要插入空行，在模板中把多余空行去掉。

自定义的模板建议存放与 `~/.emacs.d/snippets/`{:.language-bash}的对应mode文件夹下。比如在markdown模式下使用的模板，存放与 `~/.emacs.d/snippets/markdown-mode/`{:.language-bash}。

### 参考资料 ###

* [YASnippet Tutorial](http://capitaomorte.github.io/yasnippet/)

* [自定义 Yasnippet 模板](http://www.cnblogs.com/ibgo/p/3900317.html)

* [YASnippet添加模板](https://github.com/mad4alcohol/mad4a-blog/blob/master/_posts/2012-08-02-emacs-summary-cont.md)



### 更新记录 ###

2015年5月9日
