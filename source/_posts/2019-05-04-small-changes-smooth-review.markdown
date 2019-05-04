---
layout: post
title: "小的改进能提升文章审稿速度"
date: 2019-05-04 16:03:35 +0200
comments: true
categories: share
---

我最近看到Nature的一篇[报道](https://www.nature.com/articles/d41586-019-01431-z)，讲述一些小的改进能提升文章的升高速度。内容包括了草稿排版、写作和图表等方面，简要摘录关键信息。尽管遵循这些建议不会增加或减少文章被接受的概率，但能更加清楚地向读者和审稿人展现自己文章的核心内容。

1. 草稿排版

    * 宽行间距和大号字体。单列排版时，推荐使用1.5倍行距，每行12-15个单词。
    
    * 连续行编号，方便审稿人指出问题所在位置。
    
    * 图和图例尽量靠近相关文字。

<!--more-->

2. 写作

    * 避免主观词汇。审稿人通常反感主观词汇，例如：unprecedented、 paradigm shift、 amazing、 dramatic、 interesting、 remarkable等。
    
    * 减少使用缩写。当文中出现5次或以上时，再考虑使用缩写。但一些熟知的专业名词，例如DNA、RNA等要使用缩写。
    
    * 避免使用无方向性词汇，influence。应该明确指出怎样影响（增加或减少），最好能有具体的数值表述。例如，改变×××倍提高了×××产量等。
    
    * 避免在所有情况下都使用significant。应为signficant容易和统计检验混淆，应该使用具体数值描述影响。
    
    * 文章起一个陈述性的题目。避免使用暗示行、假设性的题目，而应该具体描述文章的发现。注意不要夸大结果或者模糊描述。“通过避免使用主观词汇、仔细描述文章真实而非潜在暗示的发现，应该能找到一个可读性高、信息量大、有趣的题目，从而不会夸大文章的发现”。
    
    * 摘要写法，Nature提供了一个[模板](https://www.nature.com/documents/nature-summary-paragraph.pdf)。
    
    {% img middle /images/Nature_abstract_template.jpg 900 900 'Emacs ESS #1' 'a snap of Emacs ESS' %}


3. 数据和图标

    * 定义不确定量（统计学）。在图例中描述error bar、盒箱图等不确定量的定义。
    
    * 使用[统计检验](https://www.nature.com/collections/qghhqm)。

    * 展示和提取潜在数据。例如使用散点图（小数据）、盒箱图/小提琴图（大数据）等、提交数据至数据库等。
    
    * 合理上色。避免使用彩虹图、使用[ColorBrewer](http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3)等更加友好的配色方案。比如一种diverge配色方案（热图），中间数据使用白色、两边渐进黄蓝强对比色。
    
    * 简化图。去除3D、阴影或不必要的颜色等。
    
    * 图例中添加小标题。


### 参考资料 ###

* [How small changes to a paper can help to smooth the review process](https://www.nature.com/articles/d41586-019-01431-z)

* [Statistics for Biologists](https://www.nature.com/collections/qghhqm)


### 更新记录 ###

2019年05月04日
