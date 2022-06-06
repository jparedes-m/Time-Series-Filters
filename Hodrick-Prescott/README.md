# Hodrick Prescott Filter (and its modifications)

The Hodrick Prescott Filter is a well known process to decompose an univariate time series. Suppose you have a time serie with no seasonal component,
that can be decomposed additively or multiplicatively as:

```math
y_{t}=\tau_{t}+c_{t}
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

## Code Results

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

Marcet and Ravn in 2003 in their [paper](https://bse.eu/research/working-papers/hp-filter-cross-country-comparisons) point out that having a rule of thumb that $\lambda = 1600$ for all countries quarterly gdp series is wrong. They sustain this argument by taking into account the Spanish 70-80s crisis, where having $\lambda = 1600$ for the output gap estimation didn't match the historical view of the crisis. Furthermore, the argument from $\lambda=1600$ comes from Hodrick and Prescott seminal paper where they compute that the 1600 value would be a good fit for the US post-war economy. 

Therefore, Marcet and Ravn come up with two rules:

## Rule one

As stated by Marcet and Ravn (2003): "Rule 1 may be used instead if the researcher believes that deviations from linear trend are larger in some countries considered." Therefore, this rule is good for cross-country comparisons when you try to compare business cycles between a developed country and emerging markets. Formally this rule solves the following minimization problem:

```math
\min_{\tau_{t}}\sum_{t=1}^{T}(y_{t}-\tau_{t})^{2} \quad\text{s.t.}\quad  \frac{\displaystyle \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}}{\displaystyle\sum_{t=1}^{T}(y_{t}-\tau_{t})^{2}} \leq V
```

However, this minimization problem can be arranged as: 

```math
 \min_{\tau_{t}}\sum_{t=1}^{T} (y_{t}-\tau_{t})\quad \text{s.t.} \quad \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2} \leq V\cdot \sum_{t=1}^{T}(y_{t}-\tau_{t})^{2}
```

Therefore the computation for **Rule 1** is to:

1. Assume the $\lambda = 1600$ is right for the U.S. GDP.
2. Find the V value with the following formula:

```math
V = \frac{\displaystyle \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}}{\displaystyle\sum_{t=1}^{T}(y_{t}-\tau_{t})^{2}}
```
3. Now with another dataset (i.e. Ecuadors GDP) find $\lambda$ such that $F(\lambda) = V$. We will refer to this $\lambda$ as $\lambda^{1}$

```math
F(\lambda) = \frac{\displaystyle \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}}{\displaystyle\sum_{t=1}^{T}(y_{t}-\tau_{t})^{2}} = V
```
4. Compute the Hodrick Prescott Filter with the second dataset and $\lambda^{Rule\ 1}$.

### Code results

At the file of `HP Marcet and Ravn.R`, you should have in the environment the `hp_filter` function.

You get the dataframe of the US GDP and another for Ecuador GDP. Then you execute the following code: 

```
HP_ec <- hp_filter_MR(data1 = gdp, data2 = gdp_ec, lambda = 1600, rule = "rule 1", start = 1200, end=3800)
```
In the console it will display the following result:

```
 Root finding: F(λ) - V = 0
 ----------------- 
 V: 0.000127887077934372 
 λ for data2 under rule 1 is: 3276.746 
 (F-V) = -9.56157895914966e-16 
```

With base R plot function I made this graphs for Ecuador GDP using the Rule 1 provided by Marcet and Ravn. 

![image](https://user-images.githubusercontent.com/103344273/172071831-a396eee5-f1ea-496c-8397-5d565853358c.png)

## Rule two

The second rule stated by Marcet and Ravn is to modify the restriction of the main minimization problem so that the constraint restricts the variability of the acceleration in the trend component directly.  This second rule might be used if the researcher believes that the deviation of actual trend from a linear trend is simialr across countries (i.e. UK and US).

Formally, the minimization problem becomes:

```math
\min_{\tau_{t}}\sum_{t=1}^{T}(y_{t}-\tau_{t})^{2} \quad\text{s.t.}\quad  \frac{\displaystyle \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}}{T-2} \leq W
```
By solving this optimization problem, we will find another $\lambda$, we will refer to this parameter as $\lambda^{rule\ 2}$

I solve this problem in the same way as I solved for rule 1, to assume $\lambda = 1600$ is the value for the US, then to solve for W and for series 2, find such $\lambda$ that satisfies $F(\lambda) = W$.

W (with original series)
```math
W = \frac{\displaystyle \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}}{T-2}
```

$F(\lambda)$ (with the second series)
```math
F(\lambda) = \frac{\displaystyle \sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}}{T-2}
```
### Code Results

Here is the input for the UK output gap under rule 2. 

```
HP_uk <- hp_filter_MR(data1 = gdp, data2 = gdp_uk, lambda = 1600, rule = "rule 2", end=3800)
```

This is the console output:
```
 Root finding: F(λ) - W = 0
 ----------------- 
 W: 3.04786220922701e-08 
 λ for data2 under rule 2 is: 4591.0397 
 (F-V) = -3.51968659363763e-19 
```

Here are the results for the UK under rule 2

![image](https://user-images.githubusercontent.com/103344273/172076199-b58362af-1ed6-407e-87e3-ac0fa412e5f6.png)

# References:
Marcet, A., & Morten, R. (2015, septiembre 14). The HP-Filter in Cross-Country Comparisons | Barcelona School of Economics Working Papers. Barcelona School of Economics. https://bse.eu/research/working-papers/hp-filter-cross-country-comparisons

