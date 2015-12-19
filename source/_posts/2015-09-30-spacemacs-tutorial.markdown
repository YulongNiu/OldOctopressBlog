---
layout: post
title: "Spacemacs使用记录"
date: 2015-09-30 16:31:09 +0800
comments: true
categories: lisp 
---

[Spacemacs](https://github.com/syl20bnr/spacemacs)结合了Vim和Emacs，而且定制了各种好用的设置，极大地减轻了Emacs的配置负担。推荐使用`hybrid`模式，这样浏览代码时可以使用Vim快捷键，进入Insert mode后使用Emacs快捷键。`hybrid`模式还有一个好处，编辑时方便汉字输入。Vim模式的先导键是`SPC`，在Emacs模式下是`M-m`。

<!--more-->

## 1. Vim快捷键记录 ##

在Spacemacs中`SPC-h-T`进入Emacs Evil快速入门。

### 1.1 移动 ###

* `h`：前

* `j`：上

* `k`：下

* `l`：后

* `gg`：

    * 文档开始位置。

    * `G`：文档结束位置。

    * `:[number]`：跳转到制定行。

### 1.2 插入和替换 ###

* `i`：

    * 在**光标前**的位置，进入insert mode，进行编辑。

    * `ESC`或者`C-[`退出insert mode。

* `r`：替换**光标所在位置**的单个字符。

* `[number] d object`：

    * `cw`或`ce`：从光标处删除整个单词，并进入insert mode。
    
    * `c$`：从光标处删除整行，并进入insert mode。

### 1.3 搜索和替换 ###

* `\`：

    * `\searchWord`：输入搜索内容，`n`向下搜索，`N`向上搜索。

    * `%`：在各种括号跳转。

    * `:s/old/new`：替换第一个匹配；`:s/old/new/g`：替换当前行匹配；`:#,#s/old/new/g`：替换行（`#`为行号）之间匹配；`:%s/old/new/g `：替换全文匹配。

    * `:%s/old/new/gc`：替换全文匹配，每一个匹配会提示是否匹配（输入`y`表示执行替换，`n`表示跳过匹配）。

### 1.4 删除、剪切和粘贴 ###

* `x`：删除**光标所在位置**的单个字符。

* `[number] d object`：

    * `dw`：从光标处删除整个单词，包括单词后的空格。

    * `de`：从光标处删除整个单词，不包括空格。

    * `d$`：从光标处删除整行。

    * `dd`：删除光标所在的整行。

* `p`：在**光标后**的位置粘贴剪切（`d`类和`x`操作）的内容。

### 1.5 撤销 ###

* `u`：撤销

* `Ctr-R`：反撤销










