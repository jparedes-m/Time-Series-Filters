# Hodrick Prescott Filter (and its modifications)

The Hodrick Prescott Filter is a well known process to decompose an univariate time series. Suppose you have a time serie with no seasonal component,
that can be decomposed additively or multiplicatively as:

#### Additively
```math
y_{t}=\tau_{t}+c_{t}
```


#### Multiplicatively
```math
y_{t} =  \tau_{t} \cdot c_{t}
```

$\tau_{t}$ is the series trend while $c_{t}$ is the series cycle. 

In the context of macroeconomics, say $y_{t}$ is the GDP (output), therefore the trend will be the *potential gdp* while the cycle will be the *output gap*

## The optimization problem 
So Hodrick - Prescott filters looks for a vector $\boldsymbol{\tau}$ that solves the following minimization problem:

```math
\min_{\boldsymbol{\tau}} \sum^{T}_{t=1}(y_{t}-\tau_{t})^{2}+\lambda\sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}
```
To solve this problem:
- We will find a $T \times 1$ vector $\tau$ that will have the trend values of the serie
- We will have $T$ First order conditions, where its coefficients will form the ($T \times T$) $\mathbf{F}$ matrix

```math
\boldsymbol{\tau}_{t} = (\mathbf{I}+\lambda \mathbf{F})^{-1} \mathbf{y}
```
where:
- $\mathbf{I}$ is the $T \times T$ identity matrix
- $\mathbf{y}$ is the original series
- $\mathbf{F}$ is a $T\times T$ pentadiagonal matrix

```math
\mathbf{F}=\left[ \begin{matrix}
1 & -2 & 1 &  &  &  &  \\
-2 & 5 & -4 & 1 &  &  &    \\
1 & -4 & 6 & -4 & 1 &  &     \\ 
  &  \ddots & \ddots & \ddots & \ddots & \ddots &   \\
  &  & 1 & -4 & 6 & -4 & 1  \\
  &  &  & 1 & -4 & 5 & -2 \\
  &  &  &  & 1 & -2 & 1 \\
\end{matrix}
\right]
```

## Results of the code

At the end of the script there is the following code that accounts as an example of a 7 observation univariate serie. 

```
> gdp7 <- filter(gdp, row_number()<=7)
> hp_filter(data= gdp7, lambda = 1600, show_F = TRUE, additional_info = T)
```
That will display the following result: 

```
 Hodrick - Prescott Filter 
 -------------------------
 λ = 1600 
 Series Type: y = trend + cycle 

     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]    1   -2    1    0    0    0    0
[2,]   -2    5   -4    1    0    0    0
[3,]    1   -4    6   -4    1    0    0
[4,]    0    1   -4    6   -4    1    0
[5,]    0    0    1   -4    6   -4    1
[6,]    0    0    0    1   -4    5   -2
[7,]    0    0    0    0    1   -2    1
        date variable    trend        cycle
1 2000-01-01 9.467712 9.476978 -0.009266082
2 2000-04-01 9.485754 9.480357  0.005397061
3 2000-07-01 9.486751 9.483730  0.003021151
4 2000-10-01 9.492677 9.487095  0.005582231
5 2001-01-01 9.489429 9.490451 -0.001021451
6 2001-04-01 9.495624 9.493801  0.001822790
7 2001-07-01 9.491613 9.497148 -0.005535700
```

As we can see we got some information on the variable, the penality $\lambda$, the $\mathbf{F}$ matrix. So assigning it to an object, it returns the dataframe, for the U.S output here are the results:

   ![image](https://user-images.githubusercontent.com/103344273/171779709-1404c923-0afd-44be-9bed-d5723963faac.png)

# Marcet and Ravn approach

Marcet and Ravn in 2003 in their [paper][https://bse.eu/research/working-papers/hp-filter-cross-country-comparisons] point out that mantaining $\lambda$ constant
