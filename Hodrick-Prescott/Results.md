# Results on the Hodrick Prescott Filter

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

$`\tau_{t}`$ is the series trend while $`c_{t}`$ is the series cycle. 

In the context of macroeconomics, say $`y_{t}`$ is the GDP (output), therefore the trend will be the *potential gdp* while the cycle will be the *output gap*

## The optimization problem 
So Hodrick - Prescott filters looks for a vector $`\boldsymbol{\tau}`$ that solves the following minimization problem:

```math
\min_{\boldsymbol{\tau}} \sum^{T}_{t=1}(y_{t}-\tau_{t})^{2}+\lambda\sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}
```
To solve this problem:
- We will find a $`T \times 1`$ vector $`\tau`$ that will have the trend values of the serie
- We will have $`T`$ First order conditions, where its coefficients will form the ($`T \times T`$) $`\mathbf{F}`$ matrix

```math
\boldsymbol{\tau}_{t} = (\mathbf{I}+\lambda \mathbf{F})^{-1} \mathbf{y}
```
where:
- $`\mathbf{I}`$ is the $`T \times T`$ identity matrix
- $`\mathbf{y}`$ is the original series
- $`\mathbf{F}`$ is a $`T\times T`$ pentadiagonal matrix

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
gdp7 <- filter(gdp, row_number()<=7)
hp_filter(data= gdp7, lambda = 1600, var_name = `Gross Domestic Product`, show_F = TRUE)
```
That will display the following result: 

```
 Hoddrick - Prescott Filter 
 -------------------------
 Variable: Gross Domestic Product 
 λ = 1600 
 Series Type: Additive |  y=t+c 

     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]    1   -2    1    0    0    0    0
[2,]   -2    5   -4    1    0    0    0
[3,]    1   -4    6   -4    1    0    0
[4,]    0    1   -4    6   -4    1    0
[5,]    0    0    1   -4    6   -4    1
[6,]    0    0    0    1   -4    5   -2
[7,]    0    0    0    0    1   -2    1
        date variable    trend         cycle
1 1947-01-01 7.617981 7.606784  0.0111964202
2 1947-04-01 7.615310 7.616166 -0.0008557825
3 1947-07-01 7.613243 7.625555 -0.0123115276
4 1947-10-01 7.628765 7.634957 -0.0061918890
5 1948-01-01 7.643695 7.644371 -0.0006761621
6 1948-04-01 7.660067 7.653792  0.0062744931
7 1948-07-01 7.665780 7.663215  0.0025644479
```

As we can see we got some information on the variable, the penality $`\lambda`$, the $`\mathbf{F}`$ matrix (option `show_F=TRUE`) and if we assign it to an object we 
would have the dataframe containing the variable, trend and cycle columns. 

Using my habilities in `ggplot2` and the function developed I created this graph using the data for the U.S. Real GDP for the time-span of Jan - 1947 to Jan 2022. This information is in the Data Folder, furthermore, the data was obtained at the FRED St. Louis 


![image](https://user-images.githubusercontent.com/103344273/168700256-119908a5-3266-4103-8953-78338428bcf6.png)

