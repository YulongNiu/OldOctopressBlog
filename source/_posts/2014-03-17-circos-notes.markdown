---
layout: post
title: "Circos Notes"
date: 2014-03-17 15:37:45 -0400
comments: true
categories: Bioinfor
---

Mapping data onto a Circos figure requires that you identify what patterns in your data are (a) likely to be important and (b) likely to be present, and create a figure that exposes such patterns. Remember, if the pattern exists, you can't afford to miss it. If it doesn't exist, you can't afford to be fooled into thinking that it's there, or left wondering whether it's occluded by other data.

## 1. Run circos ##

{% codeblock lang:bash %}
bin/circos
-png
-svg
-conf etc/circos.conf
-outputdir /path/to/your/output/directory
-outputfile yourimage.png
{% endcodeblock %}
<!--more-->

* `-debug`{:.language-bash} is used for debugging

* `-cdump`{:.language-bash} is used for check the parsing file.

{% codeblock lang:bash %}
# store all debug output in circos.debug.txt and display only karyotype
> circos ... -group_debug _all | tee circos.debug.txt | egrep "debuggroup karyotype"
# extract other debug reports
> egrep "debuggroup rule" circos.debug.txt
> egrep "debuggroup (rule|scale)" circos.debug.txt
{% endcodeblock %}

## 2. Syntax and configure files ##

Configuration syntax is like:

{% codeblock lang:html %}
<ideogram>
 thickness = 30p
 fill      = yes
 ...
</ideogram>
{% endcodeblock %}
Usefule syntax: hierarchical structure `<ideogram>`{:.language-html}; all data points in all plot tracks `<links>`{:.language-html} and `<plots>`{:.language-html}; local data points in a given track `<link>`{:.language-html} and `<plot>`{:.language-html}; rules `<rules>`{:.language-html} and `<rule>`{:.language-html}; highlight`<hightlight>`{:.language-html}.

Extermal imports

{% codeblock lang:html %}
<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# colors, fonts and fill patterns
<<include etc/colors_fonts_patterns.conf>>
# system and debug parameters
<<include etc/housekeeping.conf>>
{% endcodeblock %}

Accessing configuration values`parameter2 = conf(parameter1)`{:.language-html} or `parameter2 = conf(block1,block2,parameter1)`{:.language-html} for:

{% codeblock lang:html %}
<block1>
<block2>
parameter1 = ...
</block2>
</block1>
{% endcodeblock %}

Eval

{% codeblock lang:html %}
track_color = blue
track_width = 100
track_start = 0.5

<plots>
<plot>
# color=blue
color = conf(track_color)
# r0 = 0.5r
r0    = eval(conf(track_start) . "r")
# r1 = 0.5r+100p
r1    = eval(conf(track_start) . "r" + conf(track_width) . "p")
</plot>
</plots>
{% endcodeblock %}

Color

Use the color name given, like `red`{:.language-html}, `lred`{:.language-html}, `vlred`{:.language-html}, `vvlred`{:.language-html}, `dred`{:.language-html}, `vdred`{:.language-html}, `vvdred`{:.language-html}. If a pure color is need, use `pred`{:.language-html}. For the transparent color, `red_a1`{:.language-html}, `red_a2`{:.language-html}, and to `red_a5`{:.language-html}. `red_a1`{:.language-html} has a 17% tansparency and `red_a5`{:.language-html} have a 83% transparency.

## 3. Ideogram ##

One of the most common used ideogram is the karyotypes plot. The chromosome is marked as:

{% codeblock lang:perl %} 
chr - id label start end color
{% endcodeblock %}
for example

{% codeblock lang:perl %}
chr - hs1 1 0 249250621 chr1
{% endcodeblock %}

The cytogenetic bands data is like:

{% codeblock lang:perl %}
band hs1 p36.33 p36.33 0 2300000 gneg
{% endcodeblock %}
* Choose a subset of chromosomes

{% codeblock lang:perl %}
chromosomes_display_default = no

chromosomes = hs1;hs2;hs3;h4
chromosomes = /hs[1-4]$/
chromosomes = /hs[1-4]$/;hs10;hs11
chromosomes = /hs[1-4]$/;-hs3
chromosomes = hs1:(-100,120-);hs2;hs3;h4
{% endcodeblock %}

* Order

{% codeblock lang:perl %}
chromosomes_order = hs2,hs3,hs1,hs5,hs4,hs8,hs7,hs6
{% endcodeblock %}

* Scale

{% codeblock lang:perl %}
chromosomes_scale = hs1=0.25,hs2=2.0
# 25% and 50%
chromosomes_scale = hs1=0.25r,hs2=0.50r
# each is 25%
chromosomes_scale = /hs[12]/=0.5rn
# all ideograms distributed evenly within entire figure
chromosomes_scale   = /./=1rn
{% endcodeblock %}

* Rotation

{% codeblock lang:html %}
<image>
angle_orientation* = counterclockwise
<<include etc/image>>
</image>
{% endcodeblock %}

{% codeblock lang:perl %}
chromosomes_reverse = /hs[234]/
{% endcodeblock %}

* Chromosome color

{% codeblock lang:perl %}
chromosomes_color = hs1=red,hs2=orange,hs3=green,hs4=blue
{% endcodeblock %}
Also, the color could be redefined as

{% codeblock lang:perl %}
chr1* = red
chr2* = orange
chr3* = green
chr4* = blue
{% endcodeblock %}

* Position

{% codeblock lang:perl %}
chromosomes_radius = hs4:0.9r
{% endcodeblock %}

* Show bands

{% codeblock lang:perl %}
band_transparency = 0
{% endcodeblock %}

* Add chromosome

For example, add the human mitochondria chromosome.

Add the following to the file `data/karyotype/karyotype.human.txt`{:.language-bash}

{% codeblock lang:perl %}
chr - hsMT MT 0 16569 chrMT
{% endcodeblock %}

Add "hsMT" color in the file `etc/colors.ucsc.conf`{:.language-bash}

{% codeblock lang:perl %}
chrMT = 121,204,61
{% endcodeblock %}


## 4. Highlight ##

### 4.1 Data type ###

Data type of highlight is

{% codeblock lang:perl %}
chr start end
{% endcodeblock %}
for exmaple,

{% codeblock lang:perl %}
hs1 1298972 1300443
{% endcodeblock %}

We can also add the highlight information into the data file, for example:

{% codeblock lang:perl %}
hs1 100433463 100487964 fill_color=red,r0=0.6r,r1=0.6r+50p
hs1 232817594 240828534 fill_color=chr9,z=61,r0=0.4r-78.7058p,r1=0.4r+78.7058p
{% endcodeblock %}

### 4.2 Conf format ###

{% codeblock lang:html %}

<highlights>

<highlight>
...
</highlight>

<highlight>
...
</highlight>

<highlight>
...
</highlight>

</highlights>
{% endcodeblock %}

### 4.3 Parameters###

`r0`: inner radius of highlight

`r1`: outer radius of highlight

`offset`:  an offset applied to both r0 and r1 (useful for overriding default r0,r1 values defined at lower precedence)

`fill_color`: color of the highlight slice

`stroke_color`: color of the highlight border, drawn if stroke_thickness is set

`stroke_thickness`: border thickness, if any, of the highlight slice

`z`: z-depth of the highlight, controlling the order in which highlights are drawn

`ideogram`: toggles the position of the highlight to be within the ideogram extent

### 4.4 Plot hightlight ###

{% codeblock lang:html %}
<plots>
<plot>
type = highlight
file = data/3/chr.highlights.txt
r0   = 0.3r
r1   = 0.35r
z    = 10
</plot>
</plots>
{% endcodeblock %}

## 5. Links ##

### 5.1 Data type ###

{% codeblock lang:perl %}
hs1 100 200 hs2 250 300
{% endcodeblock %}
or

{% codeblock lang:perl %}
segdup00010 hs1 100 200
segdup00010 hs2 250 300
{% endcodeblock %}

with link options

{% codeblock lang:perl %}
hs1 100 200 hs2 250 300 color=blue
hs1 400 550 hs3 500 750 color=red,thickness=5p
hs1 600 800 hs4 150 350 color=black
{% endcodeblock %}

### 5.2 Conf format ###

{% codeblock lang:html %}
<links>
<link>
file          = data/5/segdup.txt
radius        = 0.8r
bezier_radius = 0r
color         = black_a4
thickness     = 2
</link>
</links>
{% endcodeblock %}

### 5.3 Parameters###

`radius`: this is the radial position of the termination of the link; for relative values, if radius < 1 then it is defined in terms of the inner ideogram radius, otherwise it is defined in terms of the outer ideogram radius

`bezier_radius`: the radial position of the third control point (in addition to the two positions defined by the link coordinates) used to draw the Bezier curve; if this parameter is not defined then straight lines will be used

`color`: color of the link line

`thickness`: thickness of the link line (note that this is not stroke_thickness, since the line isn't technically stroked)

`record_limit`: if this is defined, the number of records read from the file is capped; coordinate records are sampled from the start of the file; useful for debugging

In each `<link>`{:.language-html}, `<rules>`{:.language-html} and `<rule>`{:.language-html} can be applied to special ruls. Each `<rule>`{:.language-html} contains three parts: *a condition*, *formatting statements* and *an optional 'flow' statement*.

Marker the priority of rules

{% codeblock lang:html %}
<rule>
# 1st
importance = 10
</rule>
{% endcodeblock %}

Add a tag

{% codeblock lang:html %}
<rule>
flow = goto special_rule if true
...
</rule>

<rule>
tag = special_rule
...
</rule>
{% endcodeblock %}

## 6. Histograms ##

Histogram, line plot, scatter plot, and heat map share the same data format

{% codeblock lang:perl %}
chr start end value [options]
{% endcodeblock %}
for exmaple

{% codeblock lang:perl %}
hs3 196000000 197999999 71.0000 fill_color=blue
{% endcodeblock %}
{% codeblock lang:perl %}
# in data file
hs3 196000000 197999999 71.0000 id=abc

# in rule block
<rule>
condition  = var(id) eq "abc"
fill_color = blue
</rule>
{% endcodeblock %}
Another data format has multiple values

{% codeblock lang:perl %}
chr start end value,value,value,... [options]
{% endcodeblock %}
for example

{% codeblock lang:perl %}
hs3 196000000 197999999 0.0000,7.0000,64.0000,0.0000
{% endcodeblock %}

* Plot histogram

{% codeblock lang:html %}
<plots>
<plot>
type = histogram
file = data/5/segdup.hs1234.hist.txt
r1   = 0.88r
r0   = 0.81r

fill_color = vdgrey
extend_bin = no
</plot>
{% endcodeblock %}

For multiple values, `fill_color`{:.language-bash} is used to set the different colors.

* Orientation

{% codeblock lang:html %}
orientation = in
{% endcodeblock %}

* Rules

{% codeblock lang:html %}
<rules>
<<include exclude.hs1.rule>>
</rules>
{% endcodeblock %}

* Backgrounds and axes

{% codeblock lang:html %}
<plot>
...

<axes>
 <axis>
 ...
 </axis>
 <axis>
 ...
 </axis>
 ...
</axes>

<backgrounds>
 <background>
 ...
 </background>
 <background>
 ...
 </background>
 ...
</backgrounds>

</plot>
{% endcodeblock %}

* text

text data format is

{% codeblock lang:perl %}
chr start end value
{% endcodeblock %}

for exmaple

{% codeblock lang:perl %}
hs1 100425066 100487997 DBT
{% endcodeblock %}

{% codeblock lang:html %}
<plots>
<plot>
type             = text
color            = black
file             = data/phylo/labelGene.txt
</plot>
</plots>
{% endcodeblock %}

apply rules

{% codeblock lang:html %}
<rules>
<rule>
importance = 90
condition  = var(value) eq "GeneName"
color = blue
</rule>
</rules>
{% endcodeblock %}



## Reference ##

* [Ciros grocery of published pictures](http://circos.ca/images/scientific_literature/)

