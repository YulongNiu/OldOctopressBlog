---
layout: post
title: "模拟退火算法的思考"
date: 2017-07-23 21:11:42 +0800
comments: true
categories: bioinfor
---

模拟退火（simulated annealing，SA）是一种基于蒙特卡洛思想的近似最优化算法，突出的优势是避免陷入局部最优。该算法受启发与冶金学上的退火过程，例如金属、玻璃等达到融点后，温度缓慢降低，以使得内部形成完美晶体，从而获得高质量的材料。达到完美晶体时，全局能量最低，对应了最优化的极值点。
<!--more-->

SA之所以能跳出局部最优解，是因为在随机搜索过程中允许依照概率接受非局部最优数值。这个概率即为：

$$
p = \begin{cases}
1 & \text{if } \Delta E > 0,\\
e^{\frac{\Delta E}{T}} & \text{if } \Delta E < 0
\end{cases}
$$

其中，$\Delta E = E(x_{new}) - E(x_{old})$表示前后两个随机搜索值之差，$T$表示当前温度。

### 参考资料 ###

* [Search and Optimization by Metaheuristics-Techniques and Algorithms Inspired by Nature](https://link.springer.com/book/10.1007/978-3-319-41192-7/page/1)

* 《如何求解问题-现代启发式方法》

* 中文科普[1](https://www.cnblogs.com/ranjiewen/p/6084052.html)、[2](https://www.cnblogs.com/heaad/archive/2010/12/20/1911614.html)

### 更新记录 ###

2017年7月23日

