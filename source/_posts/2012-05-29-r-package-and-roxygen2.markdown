---
layout: post
bootstrap_theme_url: http://bootswatch.com/readable/bootstrap.css
title: "R包制作和roxygen2使用说明"
date: 2012-05-29 19:46:59 -0400
comments: true
categories: r
---

查看创建R包的各种细节，权威的文献是[Writing R Extensions](http://cran.r-project.org/doc/manuals/r-release/R-exts.html)。

## 1. 创建R包目录##

像盖房子一样，创建R包需要先搭建一个骨架，这个骨架往往是固定的，即一些文件夹（如`R`，`man`等）和文件（如`DESCRIPTION`，`NAMESPACE`）是必须的，而另外一些则可选择性添加。一个典型的R包目录,比如<span style="color: blue">Biobase</span>包，如下图：

<!--more-->

{% img middle /images/r_catalog.png 350 300 'Biobase catalog #1' 'a snap of Biobase catalog' %}


### 1.1 `DESCRIPTION`文件 ###

一个纯文本ASCII文件，其中Package, Version, License, Description, Title, Author和Maintainer是必备条目，一个例子<span style="color: blue">knitr</span>包的`DESCRIPTION`文件：

* Package：由字母、数字和“点”构成的，至少含有两个字符，第一个必须是字母，结束不能用“点”（这也是包名称的命名规范）。

* Version：版本号，经典的命名方式比如0.1-0。

* License：GPL。

* Description：一段话描述包主要的功能。

* Title：包的一句话解释。

* Author和Maintainer：包的作者和维护者，姓+名 （比如Karl Pearson ）。也可以使使用Author@R描述作者，使用R函数`person()`{:.language-r}，其中变量`role`的选项`aut`表示author，`cre`表示creator（维护者），`ctb`表示contributor。一个例子：

{% codeblock lang:r %}
c(person("Hadley", "Wickham", email = "hadley@rice.edu", role ="ctb"),  
person("Yihui", "Xie", email = "xie@yihui.name", role = c("aut", "cre")))
{% endcodeblock %}

* Date（可选）：当前版本包日期，格式yyyy-mm-dd。

* Depends（可选）：依赖的R环境版本和包名称，比如 `R (>= 2.14.1)`。

<span style="color: blue">注意：</span>

> 1. 加上版本号，没有版本号等于没用，因为版本之间可能差异很大；
>
> 2. 对于<span style="color: blue">base</span>等这样的包就不用写了，因为是必然依赖的，同时也是R启动自动载入的。对于依赖的包名称，在R 2.14.0之后完全没有必要写，因为从这个版本后，所有包都有`NAMESPACE`，直接使用`Imports`就可以了。

* Imports（可选）：如果只是使用某些包中类、方法或者（一般）函数，而不用完全载入包，可以在此栏列出包的名称，最好加上版本号（在`R CMD check`{:.language-bash}会检查版本）。在代码中，引用其他包的namespace可以使用`::`或者`:::`操作符。与之对应的，需要在`NAMESPACE`文件中指明引用（见后）。

* Suggests（可选）：如果只是在帮助文档的`examples`，`tests`或者vignettes中用到了一些包，那么没有必要“依赖”或者“引用”，只用“建议”安装即可。版本号同样也要加上，在`R CMD check`{:.language-bash}时会用到。当然，我们要考虑到如果读者也想重现一下`examples/tests/vignettes`{:.language-bash}的例子，最好使用`if(require(pkgname)))`{:.languag-r}的条件句控制：`TRUE`{:.languag-r}执行，`FALSE`{:.languag-r}返回错误。

* Enhances（可选）：这个条目中的包可以被这个正在创建的包“增强”功能，比如我们创建的这个包为（对列表中的包）提供了更好的方法、更优秀的处理手段等。

* URL（可选）：一系列相关的网址，使用“逗号”或者空格隔开。链接可能是作者的主页、包的补充材料等。

* BugReports（可选）：一个网址，用于提交bug，代替了向作者发邮件。一个好的想法是使用[Github](https://github.com/)，在项目的issures版块提交bug。

* Collate（可选）：看不懂。

* KeepSource（可选）：逻辑值，设定为`TRUE`{:.languag-r}，表示包是源码。

* ByteCompile（可选）：逻辑值，目前默认设置为否；如果性能确实可以得到提升，可以选择`TRUE`{:.languag-r}，但这可能增加安装时间和体积。
 
* BuildVignettes（可选）：如果设定位`FALSE`{:.languag-r}，`R CMD build`{:.languag-bash}不会重新编辑vignette，同时`R CMD check`{:.languag-bash}也不会检查。

### 1.2 包目录下的各种文件夹说明 ###

* R文件夹：这个文件夹下放置用R语言写成的源代码，文件最好使用后缀`*.R`。

* man文件夹：man文件夹下放置`*.Rd`结尾的文档，这写文档详细描述的函数的用途。如果使用roxygen2，就不需要建立这个文件夹。

* src文件夹：放置用其他语言写成的代码，比如C/C++/Fortran等。

* data文件夹：放置一些data文件，比如R代码、表格（`*txt`/`*.tab`/`*.csv`）文件或者是用`save()`{:.language-r}函数保存的镜像（`*.RData`/`*.rda`）。这些文件必须是“自给自足”的，即不用调用包也能使用这些data。如果这些data比较大，可以将其压缩为`*.gzip`/`*.bzip2`/`*.xz`格式，或者在保存镜像时，使用`save(, compress = TRUE)`{:.language-r}。这些方法可以有效地减少包的体积和安装时间，以方便用户使用。

* demo文件夹：放置一些包的示例，通过`demo()`函数访问。（需要有一个00Index文件？不能自动生成？）

* exec文件夹：放置脚本文件比如shell、Perl、Python、Tcl等。

* inst文件夹：包在编译安装中，会将这个文件下的文件放入包的根目录下，所以可以放包的更新文件NEWS（或者CHANGELOG）文件、用于引用的CITATION文件等。还有很重要的子文件夹是doc，vignette文件就放在此处，格式为*.Rnw，这个文件会进一步被编程成LaTeX文件，最后生成常见的PDF文档。现在Swear过时了，用<span style="color: blue">knitr</span>吧。

## 2. `NAMESPACE`文件 ##

### 2.1 文件说明 ###

`NAMESPACE`文件非常重要，因为它指明哪些函数是“输出（export）”给用户使用，哪些又是从其他包中“输入（import）”的。如果作者没有指明这个文件，那么在创建包时会被自动创建，此时所有的对象都被输出，在`Imports`和`Depends`下的包全部被输入。

**“输出”格式**：`export(f, g)`{:.language-r}，其中`f`{:.language-r}和`g`{:.language-r}是变量名，一般不用引号括起。也可以使用正则表达式批量输出，比如 `exportPattern("^[^\\.]")`{:.language-r}。

**“输入”格式**：`import(for, bar)`{:.language-r}，其中`for`{:.language-r}和`bar`{:.language-r}是包名称，这样做会导致`for`{:.language-r}和`bar`{:.language-r}中所有的输出对象都被引用。如果我们只想引用其中的某几个对象， 可以使用`importFrom(foo, f, g)`{:.language-r}，其中`f`{:.language-r}和`g`{:.language-r}是`foo`{:.language-r}包的变量名。<span style="color: blue">Base</span>包的命名空间被默认引用。当然，也可以输出其他包的变量，前提是这些变量需要先被输入。

### 2.2 S3类 ###

`S3method(Mymethod, Myclass)`{:.language-r}

### 2.3 S4类 ###

**“输出”格式**：如果创建的S4类需要被外部调用，使用`exportClasses(Myclass1, Myclass2)`{:.language-r}，同样可以使用`exportClassPattern()`{:.language-r}批量输出。对于S4，输出“方法（methods）”也会输出“泛函（generic）”，输出泛函同样也会输出方法。如果一个方法，在要创建的包中没有找到对应的方法，有两种可能的情况：

> 1. 泛函已经从别的包中输入，比如`importMethods()`{:.language-r}；
>
> 2. 从别的包中输入函数，但被转换为默认方法，比如`importFrom("graphics", plot)`{:.language-r}，`plot()`{:.language-r}函数又被分配其他的方法，这种情况下输入自己包的`plot()`{:.language-r}函数，使用`exportMethods(plot)`{:.language-r}。

**“输入”格式**：输入某个包的全部类或者特定的类，`ImportClasses(pkg)`{:.language-r}或者`ImportClassesFrom(pkg, class1)`{:.language-r}；输入某个包全部方法或者特定方法，`ImportMethods(pkg)`{:.language-r}或者`ImportMethodsFrom(pkg, method1)`{:.language-r}。

<span style="color: blue">stats4</span>给出的例子非常好，R 2.15.0之后，对隐式泛函要求严格，必须先输入显示函数：

{% codeblock lang:r %}
  export(mle) 
  importFrom("graphics", plot) 
  importFrom("stats", optim, qchisq) 
  ## For these, we define methods or (AIC, BIC, nobs) an implicit generic: 
  importFrom("stats", AIC, BIC, coef, confint, logLik, nobs, profile, update, vcov)   exportClasses(mle, profile.mle, summary.mle) 
  ## All methods for imported generics: 
  exportMethods(coef, confint, logLik, plot, profile, summary, show, update, vcov)
  ## implicit generics which do not have any methods here 
  export(AIC, BIC, nobs)
{% endcodeblock %}


## 3. 书写R函数帮助文档 ##

书写函数文档的利器当属Emacs+ESS+<span style="color: blue">roxygen2</span>包。Emacs+ESS搭建了一个非常棒的统计学软件（R、SAS、S-Plus等）的代码编辑和运行环境，从ESS命名Emacs Speaks Statistics上就可以感受到它外漏的霸气。而[roxygen2](http://roxygen.org/)是R的一个包，仿照Doxygen构建的工具，实现了文档和源码用一个文件共同管理的想法。

这个想法相当棒：

> 1. 在写代码时，方便“就地”编写文档（函数功能、参数解释、实现功能、注意事项、依赖关系和例子等）；
>
> 2. 在源码修改时，及时地更新文档。

对于编辑R包来说，对应的文档是`*.Rd`（R documentation）格式的文件，以下是一个框架示意图：

{% img middle /images/r_rd.png 600 600 'Rd file #1' 'a snap of Rd file framework' %}


### 3.1 <span style="color: blue">roxygen2</span>使用注释 ###

* `@title`：文档标题，通常是函数名缩写（比如首字母缩写）的解释。出现在该函数文档的标题位置和检索（index）目录。

* `@aliases`：函数别名。这个函数的别名，如果使用，在HTML帮助文档中将出现这个别名，链接的文档还是原来函数文档；同样，`?“别名”`{:.language-r}也能调出帮助文件。具体参考<span style="color: blue">stats</span>包中，“Normal”是`dnorm()`{:.language-r}，`pnorm()`{:.language-r}，`qnorm()`{:.language-r}，`rnorm()`{:.language-r}三个函数的别名。

* `@description`：描述，对函数进行稍微详细的描述，可以讲解函数的功能，不能有空行。出现在文档开始位置。

* `@details`：函数的具体描述，啰嗦各种细节，比如在运行过程中的注意事项、输入结果解释、运用的公式等。

* `@param`：各个变量（argument）的解释。

* `@inheritParams`：用于需要重复解释的变量，防止过多“剪切”“复制”。比如一些变量已经在foo1解释了，而foo2需要用foo1的一些变量，就可以使用这个字段。这里的foo1可以是当前包的函数，也可以是其他包函数，比如pkg::foo1。

* `@return`：函数返回的结果。如果想具体解释返回的结果，格式`\item{name a}{description a}`{:.language-html}。

* `@author`：作者，加上作者邮件地址，格式：`\email{MrR@@stat.cn}`{:.language-html}。

* `@examples`：给出函数使用的例子。例子非常重要，很多用户对例子的关注要远远高于其他的部分，啰啰嗦嗦说了一大堆，不如来两个漂亮的例子。我认为好的例子应该包括（但不限于）：

> * 重要参数的演示；
>
> * 参数之间的协调控制；
>
> * 特殊参数（用语言描述困难）的演示；
>
> * 可能的错误使用；
> 
> * 来几张漂亮图。

这些例子可以直接使用R代码编排、“#”注释、多行也问题。这些例子可以用函数`example()`{:.language-r}运行。如果不想运行某些例子（比如演示错误使用等），格式：`\dontrun{Rcode}`{:.language-html}。如果例子中依赖其他包，使用`if(require(pkg)){ }`{:.language-r}格式。

也可以在包的根目录下建立文件，通过`@example path/relative/to/packge/root`{:.language-r}引用。如果函数不需要输出，则不用添加`@example`，否则在R包检查中会报错。

* `@references`： 参考文献，比如算法出处、函数编写思路、例子的参考文献等。如果想多行出现，只需再用一次`@references`即可。如果要网址，格式为`\url{http://www.google.com}`{:.language-html}。

* `@seealso`：其他可供参考或者有联系的对象，并说明可供参考的原因。如果想在HTML中有超链接，格式为`\code{\link{function}}`{:.language-html}。

* `@family`：为`@seealso`提供了另外一种解决方法，将需要相互引用函数编排成一个“家族（family）”，之后同一个家族的函数自动出现在see-also中。同一个函数可以属于多个不同家族。

* `@note`：添加一些注释信息。

* `@section`：用户自定义一个注释区（包括名称和内容）。

* `@rdname`：解决了多个文档合并的问题。可以将多个函数标记为一个共同的`@rdname（override name）`{:.language-html}，这样在最后的HTML文档中，就可以将多个功能相似或者一个系列的函数总结在一个文档中。如果使用了`@rdname`，那么`@description`和`@details`需要在第一个函数中全部写完，`@title`也只在第一个函数中出现。

* `@keywords internal`：隐藏函数，不在help文档中显示，在其他包中可以正常使用。

* `@aliases`：命名一个函数的别名，将指向同一个文档；`@family`：如果几个函数共用一个family name，它们将在seealso中交叉链接；`@seealso`：其他函数的超链接；`@rdname`：将多个函数合并在一个文档中。这些标记用于在HTML帮助文档中，可能出现这些超级链接或者文档合并。

* 对data文件夹下的R data file进行注释。R data file通常以.RData或者.rda作为文件名扩展。[参考数据注释](http://r-pkgs.had.co.nz/data.html#documenting-data)。


## 4. 写vignittes文件##

推荐使用[<span style="color: blue">knitr</span>包](https://github.com/yihui/knitr)。

## 5. 检查和压缩包 ##

在提交代码之前，首先需要检查自己的代码，运行命名

{% codeblock lang:bash %}
R CMD check Mypkg

{% endcodeblock %}

之后，因为R源码包（source package）的发布都是采用压缩文件（`*.tar.gz`）形式，所以需要将源代码打包，方法是运行命令

{% codeblock lang:bash %}
R CMD build Mypkg

{% endcodeblock %}

当然，有时我们需要发布二进制包（binary packages）,操作方法是

{% codeblock lang:bash %}
R CMD INSTALL --build Mypkg
{% endcodeblock %}

具体执行过程是：首先，在默认目录树下安装包；之后，将成功安装的包转化成（平台依赖的）二进制代码，再将二进制代码打包输出。所以，默认安装目录要有“可写”权限，如果没有，执行`R CMD INSTALL -l location --build Mypkg`{:.language-bash}，指定一个“可写”权限目录（location）安装。

CRAN很体贴地考虑了Window以外平台的作者，通过[网站](http://win-builder.r-project.org/)上传源码包，之后编译成window平台二进制包。

最后是提交代码，依次执行以下命名或者操作:

{% codeblock lang:bash %}
R CMD build Mypkg
# 检查例子（example）、检查（test）和文档输出，如果有报错或者警告，需要仔细阅读提示，原则上消除所有报错和警告
R CMD check --as-cran
{% endcodeblock %}

使用 ftp://CRAN.R-project.org/incoming/ （用户名： anonymous；密码：自己邮箱地址）上传自己的包，并向CRAN@R-project.org发邮件通知自己上传了包（名称、版本等）。


## 注意事项 ##

### 1. 命名 ###

给自己的包起一个漂亮的名字，注意大小写，因为一些操作系统可能对大小写不敏感。因此，为了平台兼容性，包的命名尽量不要用大小写区分；不仅如此，在包源代码文件夹下，文件命名也要注意这个问题。起了名字之后，最好到[R包列表](http://cran.r-project.org/web/packages/available_packages_by_name.html)去看看，检查下自己包的名字是否与其他人的冲突。

### 2. 名称符号 ###

> * 文件名中不允许出现如下符号`"`、 `*`、 `:`、`/`、`<</font>`、 `>`、 `?`、 `\` 和 `|`；
>
> * 文件命名使用的字符包括 `A-Za-z0-9._!#$%&+,;=@^(){}'[]`，但是不推荐使用 `(){}'[]$`；
>
> * `*.Rd`文件可以包括URL，但必须是ASCII，不能出现`%`。

### 3. GNU版本号 ###

GNU版本命名规则：Major_Version_Number.Minor_Version_Number[.Revision_Number[.Build_Number]]。

### 4. 关于“命名空间（namespace）”的讨论 ###

在制作R包时候，有两种方法可以引用其他包的变量：

> **第一种方法**：只需要在`NAMESPACE`中声明引用`import()`或者`importFrom()`，同时在`DESCRIPTION`中`Imports`添加需引用的包即可。这样，在编写代码时不用加上`require(pkg)`{:.language-r}之类载入包的语句。在安装过程中，通过引用可以很方便地获取其他包的变量（比如函数），因而避免了载入整个包；对应的，这很大程度上保证了代码的简洁，一个非常不错的设计。
> 
> 这里需要再多啰嗦两句，`require(pkg)`{:.language-r}和`NAMESPACE`声明`import()`有不同，前一个函数载入包，并将其添加（attach）到搜索目录下；而后一个则只载入包，并有添加搜索目录。这有两点好处：1. 用户只将注意力集中在载入的包，如果import()的包中有与载入包命名冲突的函数，只会调用载入的包；2. 设计简洁。
> 
> **第二种方法**：也可以通过操作符`::`和`:::`获得命名空间的“显变量（exported variables）和“隐变量（internal variables）”。需要注意的是，不推荐使用`:::`，这可能暗示这一个自己设计不足。因为“隐变量”之所以要隐藏起来，作者是有着充分的考虑；而且，这些“隐变量”可能随着作者对包的更新而改名。可以使用函数`loadedNamespaces()`{:.language-r}查看已经载入的命名空间，`unloadNamespace(name)`卸载命名空间。

对于这两种方法，不推荐第二种方法，原因很简单：要多写包名+两个冒号。第二中方法适合在终端操作时，遇到函数名称冲突（见后）使用。

如果变量命名有冲突，R将采取“就近”策略应对，考虑两种特殊情况:

> 1. 假定自己包的名称是Mypkg，如果在`NAMESPACE`中声明输入`export(f)`，但是又从另外一包foo输入同名函数f，比如`importFrom(foo, f)`。这时，“就近”策略发挥作用，即在载入包Mypkg时，默认f函数是自己编写的，即声明输出的函数f，而不是foo中的f函数。总体而言，命名空间的查找顺序是：**包的命名空间 --> import声明 --> Base包的命名空间 --> 最后是搜索路径**。
> 
> 2. 假定pkg1和pkg2两个包，都有命名空间，而且输出同名函数f。这是，R会怎么处理呢？同样是“就近”原则，哪个包后载入，先用哪个包的函数f。不过，如果函数命名有冲突，后载入的包会给出一个警告“marsked from package: ***”。


### <a id="Ref">参考网址</a> ###

* [knitr包的DESCRIPTION文件](https://github.com/yihui/knitr/blob/master/DESCRIPTION)

* 命名空间讨论：[1](https://github.com/hadley/devtools/wiki/Namespaces)、[2](http://www.r-bloggers.com/namespaces-and-name-conflicts/)、[3](http://r-pkgs.had.co.nz/namespace.html)

* 其他帮助文档：[1](http://cos.name/2011/05/write-r-packages-like-a-ninja/)、[2](http://yihui.name/knitr/)、[3](http://cran.r-project.org/doc/manuals/R-exts.html#Dynamic-pages)、[4](http://sjp.co.nz/posts/emacs-ess-knitr/
)、[5](https://github.com/yihui/knitr/issues/252)、[6](http://cos.name/cn/topic/106644)

### 更新记录 ###

2014年9月5日
