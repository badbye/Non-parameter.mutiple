---
title       : 非参数统计
subtitle    : 多样本检验
author      : yaleidu
job         : 开心就好
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : prettify  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : mathjax            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
github:
  user: badbye
repo: Non-parameter.mutiple
---
  
## 目录
  
3.1 Kruskal-Wallis秩和检验

3.2 正态记分检验

3.3 Jonckheere-Terpstra检验

3.4 完全区组设计: Friedman秩和检验

3.5 Kendall协同系数检验

3.6 完全区组设计: 二元响应的Cochran检验

3.7 完全区组设计: Page检验

3.8 不完全区组设计: Durbin检验

--- .class #id 

## 3.1 Kruskal-Wallis秩和检验
假设数据来自相似的连续分布.

```r
kw.test = function(data, group) {
    ### lengths of each group
    n.vector = tapply(group, group, length)
    R = rank(data)
    N = length(data)
    ### squares of each group's sum group
    sum.vector = tapply(R, group, function(x) sum(x)^2)
    H = 12/N/(N + 1) * sum(sum.vector/n.vector) - 3 * (N + 1)
    ### p value to caculate
}
```


注: "p value to caculate" 是书中没有给出具体是什么分布,故无从计算p值,只计算出统计量.下同.

内置函数: kruskal.test(data, group)

---
## 3.2 正态记分检验
  

```r
norm.mutiple = function(data, group) {
    w = cbind(data, group)
    w = w[order(data), ]
    n.vec = tapply(data, group, length)
    norm.score = qnorm((1:sum(n.vec))/(sum(n.vec) + 1))
    w = cbind(w, 1:sum(n.vec), norm.score)
    T = (sum(n.vec) - 1) * sum(tapply(norm.score, w[, 2], function(x) sum(x)^2/length(x)))/sum(norm.score^2)
    p = pchisq(T, length(unique(group)) - 1)
    list(H1 = "H1: Not all the accumulate functions are identical", p.value = min(p, 
        1 - p))
}
```

---
  
## 3.3 Jonckheere-Terpstra检验
  
假定k个独立样本有同样形状的连续分布.比Kruskal-walliss检验有更强的势,即更容易拒绝原假设.    
可以选择有序的备择假设,判断数据是否有上升或者下降的趋势,如:
  
  $$H_1: \theta_1 \le \theta_2  \le \ldots \le \theta_k $$
  

```r
jt.test = function(data, group) {
    k = length(unique(group))
    U = matrix(0, k, k)
    data.cate = split(data, group)
    for (i in 1:k) {
        for (j in (i + 1):k) {
            U[i, j] = sum(sapply(data.cate[[i]], "-", data.cate[[j]]) > 0)
        }
    }
    J = sum(U)
    ## p.value to caculate
}
```


注: SAGx包中有JT.test函数.

---
  
## 3.4 完全区组设计: Friedman秩和检验
  
每个处理在每个区组恰好有一个观测值.


```r
Friedman = function(data) {
    # 行表示处理,列代表区组
    if ((is.matrix(data) | is.data.frame(data)) == FALSE) 
        stop("input should be a matrix or a dataframe")
    k = nrow(data)
    b = ncol(data)
    r = apply(m, 2, rank)
    Q = 12/b/k/(k + 1) * sum(rowSums(r)^2) - 3 * b * (k + 1)
    W = Q/b/(k - 1)
    ## p.value to caculate
}
```


注: 内置函数friedman.test()

---
  
## 3.5 Kendall协同系数检验
  

```r
kendall.test = function(data) {
    # 行数据表示按某一标准的评判结果 列数据表示同一个体在不同标准下的评判结果
    if ((is.matrix(data) | is.data.frame(data)) == FALSE) 
        stop("input should be a matrix or a dataframe")
    r = colSums(data)
    m = nrow(data)
    n = ncol(data)
    S = sum((r - m * (n + 1)/2)^2)
    W = 12 * S/m^2/(n^3 - n)
    chi.p = pchisq(m * (n - 1) * W, n - 1, low = F)
    F.p = pf((m - 1) * W/(1 - W), n - 1 - 2/m, (m - 1) * (n - 1 - 2/m), low = F)
    list(H1 = "Samples are related to some extent", Chisq.pValue = chi.p, F.pValue = F.p)
}
```

---
  
## 3.6 完全区组设计: 二元响应的Cochran检验
  
二元变量的观测值.


```r
cochran.test = function(data) {
    # 行数据表示按某一标准对所有个体的评判结果
    # 列数据表示同一个体在不同标准下的评判结果
    if ((is.matrix(data) | is.data.frame(data)) == FALSE) 
        stop("input should be a matrix or a dataframe")
    
    n = colSums(data)
    l = rowSums(data)
    N = sum(n)
    k = ncol(data)
    Q = k * (k - 1) * sum((n - mean(n))^2)/(k * N - sum(l^2))
    list(H1 = "Not all the individuals are regarded equally", p.value = pchisq(Q, 
        k - 1, low = F))
}
```


---
## 3.7 完全区组设计: Page检验
  
正态近似,自动处理打结情况.

```r
page.test = function(data) {
    # 行数据表示同一区组的不同处理结果 列数据表示不同区组的同一处理结果
    if ((is.matrix(data) | is.data.frame(data)) == FALSE) 
        stop("input should be a matrix or a dataframe")
    r = apply(data, 1, rank)
    R = rowSums(r)
    L = sum(R * 1:length(R))
    k = ncol(data)
    b = nrow(data)
    u = b * k * (k + 1)^2/4
    tie = as.vector(r)%%1 == 0
    s = ifelse(all(tie), b * (k^3 - k)^2, k * (k^2 - 1) * (b * k * (k^2 - 1) - 
        sum(as.vector(r)[!tie]^3 - as.vector(r)[!tie])))
    Z = (L - u)/sqrt(s/144/(k - 1))
    list(H1 = "Not all the infulences of diffierent processings are the same", 
        p.value = pnorm(Z, low = F))
}
```


---
## 3.8 不完全区组设计: Durbin检验
  
这种数据的存在太反人类,该睡觉的时候不该考虑太复杂的问题.
