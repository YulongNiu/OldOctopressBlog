---
layout: post
bootstrap_theme_url: http://bootswatch.com/united/bootstrap.css
title: "My first post"
date: 2014-07-04 22:04:38 -0400
comments: true
categories: saysomething
---

This is my first post. I want to say something. Then say something.

This is a piece of html code `<html>`{:.language-html}

{% codeblock lang:ruby %}
def what?
  42
end
{% endcodeblock %}

<!-- more -->

1. Mathjax supports.

    $$
    \begin{align*}
      & \phi(x,y) = \phi \left(\sum_{i=1}^n x_ie_i, \sum_{j=1}^n y_je_j \right)
      = \sum_{i=1}^n \sum_{j=1}^n x_i y_j \phi(e_i, e_j) = \\
      & (x_1, \ldots, x_n) \left( \begin{array}{ccc}
	  \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
	  \vdots & \ddots & \vdots \\
	  \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
	\end{array} \right)
      \left( \begin{array}{c}
	  y_1 \\
	  \vdots \\
	  y_n
	\end{array} \right)
    \end{align*}
    $$

2. Colorful mathjax euqation.

    $$
    \definecolor{energy}{RGB}{114,0,172}
    \definecolor{freq}{RGB}{45,177,93}
    \definecolor{spin}{RGB}{251,0,29}
    \definecolor{signal}{RGB}{18,110,213}
    \definecolor{circle}{RGB}{217,86,16}
    \definecolor{average}{RGB}{203,23,206}
    \color{energy} X_{\color{freq} k} \color{black} =
    \color{average} \frac{1}{N} \sum_{n=0}^{N-1}
    \color{signal}x_n \color{spin}
    e^{\mathrm{i} \color{circle} 2\pi \color{freq}k
    \color{average} \frac{n}{N}}
    $$

    To find <font color="#7200AC">the energy</font>
    <font color="2DB15D">at a particular frequency</font>,
    <font color="#FB001D">spin</font> <font color="#126ED5">your
    signal</font> <font color="#D04400">around a circle</font>
    <font color="2DB15D">at that frequency</font>, and
    <font color="#CB17CE">average a bunch of points along that
    path</font>.

* Show the block.

> This is the block


