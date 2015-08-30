---
layout: post
title: "Emacs配置HTML/JS/CSS编辑环境"
date: 2012-07-30 22:31:13 +0800
comments: true
categories: lisp
---

以下所有Emacs mode都使用[Emacs包安装系统](http://yulongniu.bionutshell.org/blog/2012/06/25/emacs-extend-skills/)。

## 1. web-mode ##

[web-mode](http://web-mode.org/) 提供了很好的wen配件（HTML、JavaScript、CSS、PHP等）的代码缩进、折叠和高亮等出色的功能。

有用技巧：

* `C-c C-n`{:.language-emacs-lisp}：放在HTML标签上，在标签间跳转。

* `C-c C-f`{:.language-emacs-lisp}：放在HTML标签上，在标签折叠。
<!--more-->

## 2. emmet-mode ##

[emmet-mode](https://github.com/smihica/emmet-mode)是[emmet](http://emmet.io/)的Emacs扩展，取代了陈旧的[ZenCoding](http://www.emacswiki.org/emacs/ZenCoding)。只需要输入制定的缩写，这个工具可以自动补全HTML标签。非常方便。

有用技巧：

* `M-x emmet-mode`{:.language-emacs-lisp}：打开emmet-mode。

* `C-j`{:.language-emacs-lisp}：自动补全。

补全缩写规律：

* `a`{:.language-emacs-lisp}：a+href

* `#q`{:.language-emacs-lisp}：div+id(q)

* `.x`{:.language-emacs-lisp}：div+class(x)

* `#q.x`{:.language-emacs-lisp}：div+id(q)+class(x)

## 附加内容 ##

Emacs编程环境设置

{% codeblock lang:emacs-lisp Emacs web config %}
;;;;;;;;;;;;;;
;emmet-mode
;;;;;;;;;;;;;
(require 'emmet-mode)
(add-hook 'sgml-mode-hook 'emmet-mode) ;; Auto-start on any markup modes
(add-hook 'html-mode-hook 'emmet-mode)
(add-hook 'web-mode-hook 'emmet-mode)
(add-hook 'css-mode-hook  'emmet-mode)


;;;;;;;;;;;;;;
;web-mode
;;;;;;;;;;;;;;;
(require 'web-mode)
(add-to-list 'auto-mode-alist '("\\.phtml\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.tpl\\.php\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.[agj]sp\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.as[cp]x\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.erb\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.mustache\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.djhtml\\'" . web-mode))
(add-to-list 'auto-mode-alist '("\\.html?\\'" . web-mode))
(defun my-web-mode-hook ()
  "Hooks for Web mode."
  (setq web-mode-markup-indent-offset 2)
)
(add-hook 'web-mode-hook  'my-web-mode-hook)
{% endcodeblock %}

### 参考资料 ###

* [web-mode说明文档](http://web-mode.org/)




### 更新记录 ###

2015年8月30日

