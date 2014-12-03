---
layout: post
bootstrap_theme_url: http://bootswatch.com/flatly/bootstrap.min.css
title: "Octopress安装和使用"
date: 2014-07-22 13:45:46 -0400
comments: true
categories: PPR
---

安装和使用Octopress的一些注意事项，详细的内容可以[参考网址](#Ref)。

## 1. 安装 ##

请参考[官网](http://octopress.org/)，其他的博客介绍的安装已经失效或者错误。

常用命令:

{% codeblock lang:bash %}
# 预览，可自动更新。使用Ctrl+c终止。
$ rake preview
{% endcodeblock %}
    
## 2. Ruby版本调整 ##

因为Octopress需要使用Ruby旧版本，推荐使用[RVM](https://rvm.io/)安装Ruby 1.9.3版本。在安转过程中可能会出现`"gpg: Can't check signature: public key not found"`的错误提示，需要执行一下类似命名，添加公用匙。

{% codeblock lang:bash %}
$ gpg2 --keyserver hkp://keys.gnupg.net --recv-keys D39DC0E3
{% endcodeblock %}

使用以下操作设定ruby版本：

{% codeblock lang:bash %}
$ source ~/.rvm/scripts/rvm
$ rvm use 1.9.3
{% endcodeblock %}
<!--more-->    

## 3. 预览错误 ##

如果在预览博文时出现`TCPServer Error: Address already in use - bind(2)`{:.language-bash}的错误，表示端口（Octopress默认4000）被占，解决办法:

{% codeblock lang:bash %}
# 获取占据4000端口程序的PID
$ lsof -wni tcp:4000
$ kill -9 PID
{% endcodeblock %}

## 4. 更新博文 ##

博文放置在`source/_posts`{:.language-bash}目录下。

{% codeblock lang:bash %}
$ rake new_post["new post"]
$ rake generate
$ git add .
$ git commit -m "my comment" 
$ git push origin source
$ rake deploy
{% endcodeblock %}
如果需要在首页显示部分博文，在需要隔断的地方（博文markdown文件）加入：

{% codeblock lang:html%}
<!--more-->
{% endcodeblock %}

## 5. 使用主题 ##

当前博客使用的是[octostrap](http://kaworu.github.io/octopress/)主题。

* [添加Category侧边栏](http://kaworu.github.io/octostrap3/blog/2013/10/03/category-list-aside/)

* [每个页面更换主题](http://kaworu.github.io/octostrap3/blog/2013/10/02/pick-a-theme-for-only-one-page/)

* [选择Bootstrap主题](http://kaworu.github.io/octostrap3/setup/pick-a-theme/)

## 6. 修改标签图标记 ##

可以使用[在线转换工具](http://converticon.com/)，将png格式图片转换为ico格式（比如32*32）。之后将ico图片直接放置在`souce`{:.language-bash}文件夹下。最后，修改对应`source/_includes/head.html`{:.language-bash}文件即可。

## 7. 404公益 ##

在`source`{:.language-bash}文件夹下添加*404.markdown*文件，之后写入

{% codeblock lang:javascript %}
---
layout: page
title: "404 Error"
date: 2014-07-22
comments: false
sharing: false
footer: false
---

<script type="text/javascript" src="http://www.qq.com/404/search_children.js?edition=small" charset="utf-8"></script>
{% endcodeblock %}

## 8. 私密博文 ##

Octopress提供了隐藏博文的功能，即使文章已经推送到了github，也可以不在博客主页显示。具体方法是在每篇markdown文件头部添加：

{% codeblock lang:ruby %}
published: false
{% endcodeblock %}

如果需要公开发表，将其删除或者改为：

{% codeblock lang:ruby %}
published: true
{% endcodeblock %}

## 9. 代码高亮设置 ##

安装`coderay`和`kramdown`

{% codeblock lang:bash %}
$ gem install coderay
$ gem install kramdown
{% endcodeblock %}

在`_config.yml`{:.language-bash}文件中写入：

{% codeblock lang:bash %}
markdown: kramdown
kramdown:
  use_coderay: true
  coderay:
    coderay_line_numbers: table
    coderay_css: class
{% endcodeblock %}

## 10. Kramdown语法小技巧##

* 代码高亮

代码段高亮参考[Octopress codeblock](http://octopress.org/docs/plugins/codeblock/)

行内引用代码使用
`` `source/_includes/custom/head.html`{:.language-bash} ``


* 四个空格或者个一个Tab可以生成一个文本块

## 11. 链接在新的选择卡中打开##

在`source/_includes/custom/head.html`{:.language-bash}中添加：

{% codeblock lang:html%}
<!-- link open with new tab  -->
function addBlankTargetForLinks () {
  $('a[href^="http"]').each(function(){
  $(this).attr('target', '_blank');
  });
}

$(document).bind('DOMNodeInserted', function(event) {
  addBlankTargetForLinks();
});
{% endcodeblock %}

## 12. 修改分页数 ##

修改`_config.yml`{:.language-bash}文件

{% codeblock lang:bash %}
# 每页最多展示的博文数目
paginate: 5

# 分页后博文地址
paginate_path: "posts/:num"
{% endcodeblock %}

## 13. 添加图片 ##

将需要添加的图片移动到目录`source/images/`{:.language-bash}，之后在正文中添加：

{% codeblock lang:bash %}
{% img [class names] /path/to/image [width] [height] [title text [alt text]] %}

# 例子
{% img left /images/testimg.png 350 350 'image' 'images' %}
{% img right http://placekitten.com/300/500 150 250 Place Kitten #3 %}
{% endcodeblock %}


## 14. 修改favicon ##

需要修改的favicon文件，比如`favicon.ico`移动到`source/`{:.language-bash}目录下。之后修改文件`source/_includes/head.html`{:.language-bash}，找到`favicon.png`将其改为`favicon.ico`。


## 15. 添加新页面 ##

首先，添加新的页面：

{% codeblock lang:bash %}
rake new_page[ANewPage]
{% endcodeblock %}

这会生成一个新的文件`source/anewpage/index.markdown`{:.language-bash}。之后，修改`source/_includes/custom/navigation.html`{:.language-bash}文件，根据自己主题，添加如下类似内容

{% codeblock %}
<li {% if page.navbar == 'ANewPage' %}class="active"{% endif %}>
  <a href="{{ root_url }}/anewpage">ANewPage</a>
</li>
{% endcodeblock %}


### <a id="Ref">参考网址</a> ###

* Octopress安装和域名设置：[1](http://tchen.me/posts/2012-12-16-first-blog.html), [2](http://beyondvincent.com/blog/2013/08/03/108-creating-a-github-blog-using-octopress/)

* Octopress其他配制：[1](http://812lcl.com/blog/2013/10/26/octopressce-bian-lan-ji-ping-lun-xi-tong-ding-zhi/), [2](http://cn.soulmachine.me/blog/20130402/)

* 添加多说：[1](http://havee.me/internet/2013-02/add-duoshuo-commemt-system-into-octopress.html), [2](http://kaiimeng.cn/my-first-octopress-blog/), [3](http://cn.soulmachine.me/blog/20130402/)

* 添加Mathjax支持：[1](http://yanping.me/cn/blog/2012/03/10/octopress-with-latex/), [2](http://www.idryman.org/blog/2012/03/10/writing-math-equations-on-octopress/), [Mathjax彩色公式](http://adereth.github.io/blog/2013/11/29/colorful-equations/)

* [Kramdown语法](http://kramdown.gettalong.org/syntax.html)

* [Kramdown演示](http://kramdown.gettalong.org/quickref.html)

* [pygments错误详细提示](http://i.rexdf.org/blog/2014/09/26/octopressbo-ke-geng-xin-ri-zhi/)

* [Octopress highlight language list](http://pygments.org/docs/lexers/)

### 更新记录 ###

2014年7月25日


