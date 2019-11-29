---
layout: post
title: "R ggplot2 notes"
date: 2013-10-02 02:47:11 -0400
comments: true
categories: R
---

# R Package <span style="color: blue">ggplot2</span> Notes #

## 1. Basic grammar ##

### 1.1 Plot types ###

The R package <span style="color: blue">ggplot2</span> is a famous plot tool for high quality scientific figures. The <span style="color: blue">ggplot2</span> style figures are widely seen in papers published in high quality journals such as *PNAS*, *Nature* and *Cell*.

The input data should be in **data frame** form, and it is easily to use the function `as.data.frame()`{:.language-r}. "<span style="color:red">**+**</span>" is used to connect different plot statement. A typical <span style="color: blue">ggplot2</span> plot statement is like:

<!--more-->

{% codeblock lang:r %}
require('ggplot2')
ggplot(data=mpg, mapping=aes(x=cty, y=hwy, colour=factor(year))) +
geom_point() + stat_smooth()
{% endcodeblock %}
`ggplot()`{:.language-r}: **data** is a <span style="color:blue">data.frame</span> class object. **mapping** is an `aes()`{:.language-r} function to specify the X-axis and Y-axis. When a `aes()`{:.language-r} is used, a figure legend will be added. If we do not want the legends appear, use `show_guide = FALSE`{:.language-r} in *geom_XXX* or *stat_XXX*

`geom_point()`{:.language-r}: is used to plot points with the attributes **x**, **y**, alpha, colour, fill, shape, size.

`geom_line()`{:.language-r}: is used to plot points with the attributes **x**, **y**, alpha, colour, fill, shape, size.

`geom_bar()`{:.language-r}: bar plot. `stat = 'identity'`{:.language-r} for draw the identical, `hjust`{:.language-r} and `vjust`{:.language-r} is to adjust the x and y axis distance. `coord_flip()`{:.language-r} to reverse X and Y axis. `position = 'dodge'`{:.language-r} to set position of two bars, can be set as 'dodge', 'stack', 'fill' and 'identity'.

`geom_box()`{:.language-r}: boxplot.

`geom_tile()`{:.language-r}: fill blocks.



{% codeblock lang:r %}
# example
require('ggplot2')
p <- ggplot(mtcars, aes(factor(cyl), mpg))
# basic plot
p + geom_boxplot()
# add colors to boxes
p + geom_boxplot(aes(fill = factor(cyl)))
# change default colors
p + geom_boxplot(aes(fill = factor(cyl))) + scale_fill_manual(values = c('red', 'green', 'blue'))
{% endcodeblock %}
`geom_rect(mapping = NULL, data = NULL, stat = "identity", position = "identity", ...)`{:.language-r}: plot rectangles. 
In `aes()`{:.language-r}, `xmin`{:.language-r}, `xmax`{:.language-r}, `ymin`{:.language-r}, and `ymax`{:.language-r} are necessary. 
`inherit.aes = FALSE`{:.language-r} may be used if new `data`{:.language-r} is applied.

{% codeblock lang:r %}
# example
ggplot(mtcars) +
  geom_density(aes(x=disp, group=cyl, fill=cyl), alpha=0.6, adjust=0.75) + 
  geom_rect(data=mtcars[1,], aes(xmin=100, xmax=200, ymin=0,ymax=Inf), fill="red", alpha=0.2)
{% endcodeblock %}

### 1.2 Statistics ###

`geom_smooth()`{:.language-r}: is used for the add smooth line with the **method** lm, glm, gam, loess and rlm. `se = TRUE`{:.language-r} is to display the confident region. The following aesthetics **x**, **y**, alpha, colour, fill, linetype, size, weight could be added.

`stat_boxplot()`{:.language-r}: plot error lines in boxplot.

### 1.3 Add elements ###

`xlab()`{:.language-r}: change X axis label, set `xlab('')` to remove the X axis label; `ylab()`{:.language-r}: change Y axis label; `ggtitle()`{:.language-r}: add figure title; `scale_y_continuous(limits=c(0, 20))`{:.language-r} and `scale_x_continuous(limits=c(0, 20))`{:.language-r} to adjust range of X and Y axis.

`geom_abline(intercept = 37, slope = -5)`{:.language-r}: to add line.

`geom_hline`{:.language-r} and `geom_vline`{:.language-r}: to add horizontal and vertical lines.

{% codeblock lang:r %}
# example 
p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
geom_vline(xintercept = 1:5, colour="green", linetype = "longdash")
{% endcodeblock %}

{% codeblock lang:r %}
# ggplot2 line type
d <- data.frame(lt=c('blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash', '1F', 'F1', '4C88C488', '12345678'))
ggplot() +
  scale_x_continuous(name='', limits=c(0,1), breaks=NULL) +
  scale_y_discrete(name='linetype') +
  scale_linetype_identity() +
  geom_segment(data=d, mapping=aes(x=0, xend=1, y=lt, yend=lt, linetype=lt))
{% endcodeblock %}

`geom_text()`{:.language-r}: to add text. Set `parse = TRUE` to use expression and greek letters.

`scale_fill_discrete(..., values)`{:.language-r}: change labels. `name`{:.language-r} to reset label names, `labels`{:.language-r} to reset labels.

`scale_shape_manual(..., values)`{:.language-r}: change the shape of points.

`scale_linetype_manual(..., values)`{:.language-r}: change the types of lines. line referring [R plot](http://www.cookbook-r.com/Graphs/Shapes_and_line_types/).`name`, `value`, `labels` are used to change value.

`scale_color_manual`{:.language-r} is used for change the colors. Please refer to [Useful color palette](http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/), [Introduction of ggplot2 colors](http://blog.ggplot2.org/post/24607351280/choosing-colour-palettes-part-ii-educated-choices), and [R Color Chart](http://research.stowers-institute.org/efg/R/Color/Chart/) 
. The default ggplot2 colors are generated from the "*scales*" package, for example the default "hue pallet" could be view as `show_col(hue_pal(h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1)(9))`{:.language-r}

`scale_fill_manual(..., alues)`{:.language-r} to change filled colors.

`theme`{:.language-r} is used for exact control. `legend.position='none'` to remove the side legend.

## 2. Save plot ##

The function `ggsave()`{:.language-r} is used to save the screen plot to file. `print()`{:.language-r} is also applied like:

{% codeblock lang:r %}
pdf('testfile.pdf')
q <- ggplot()
print(q)
dev.off()
{% endcodeblock %}

## 3. Other issues ###

### 3.1 Plot mutiple ggplot2 ###

Use <span style="color: blue">gridExtra</span> package to plot multiple ggplot2 figures in the one figure. 

{% codeblock lang:r %}
# example
require('gridExtra')
# save ggplot object into a list like "plotList"
do.call(grid.arrange, plotList)
{% endcodeblock %}

### <a id="Ref">参考网址</a> ###

* [Better labels](http://directlabels.r-forge.r-project.org/examples.html)
* [ggplot2 doc](http://docs.ggplot2.org/current/)
* [ggplot2 cheatsheet](http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/#working-with-colors)

### 更新记录 ###

2014年8月28日
