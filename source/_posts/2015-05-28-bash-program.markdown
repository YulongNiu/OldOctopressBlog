---
layout: post
title: "Bash简易编程"
date: 2015-05-28 21:24:32 +0800
comments: true
published: true
categories: Linux
---

## 1. 循环 ##

`for`{:.language-bash}循环体

{% codeblock lang:bash for in Bash%}
for i in *.zip
do
    echo "$i"
done
{% endcodeblock %} 

<!--more-->

## 2. 创建数组 ##

{% codeblock lang:bash create array %}
declare -a testArray={"element1" "element2"}
echo ${testArray[0]}
for i in "${testArray[@]}"
do
    echo "$i"
done
{% endcodeblock %} 

## 3. 字符串分割 ##

{% codeblock lang:bash Split strings%}
# write in file "testsplit.sh"
IFS=',' read -ra splitArray <<< "This,is,a,test"
for i in "${splitArray[@]}"
do
    echo "$i"
done

$ bash testsplit.sh
This
is
a
test
{% endcodeblock %}

## 4. 屏幕输出存入变量 ##

某个bash命令，比如`ls -l`{:.language-bash}存入变量，之后引用变量。

{% codeblock lang:bash Save Command Output into Variables%}
# 注意等号前后不能加空格
listOutput=`ls -l`
echo "$listOutput"
{% endcodeblock %}

## 5. 文件末尾添加内容 ##

{% codeblock lang:bash append %}
touch testfile
printf "hello\n" > testfile
printf "world\n" >> testfile
{% endcodeblock %} 

## 6. 传入sudo密码 ##

{% codeblock lang:bash append %}
echo "myPassword" | sudo -S myCommond
{% endcodeblock %} 


### 参考网址 ###

* [How do I split a string on a delimiter in Bash?](http://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash)

* [How to set a BASH variable equal to the output from a command?](http://stackoverflow.com/questions/4651437/how-to-set-a-bash-variable-equal-to-the-output-from-a-command) 

### 更新记录 ###

2018年4月7日
