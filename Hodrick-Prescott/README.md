# Hodrick Prescott Filter (and its modifications)

The Hodrick Prescott Filter is a well known process to decompose an univariate time series. Suppose you have a time serie with no seasonal component,
that can be decomposed additively or multiplicatively as:

```math
y_{t}=\tau_{t}+c_{t}
```

$\tau_{t}$ is the series trend while $c_{t}$ is the series cycle. 

In the context of macroeconomics, say $y_{t}$ is the GDP (output), therefore the trend will be the *potential gdp* while the cycle will be the *output gap*

## The optimization problem 
So Hodrick - Prescott filters look for a vector $\boldsymbol{\tau}$ that solves the following minimization problem:

```math
\min_{\boldsymbol{\tau}} \sum^{T}_{t=1}(y_{t}-\tau_{t})^{2}+\lambda\sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}
```
To solve this problem:
- We will find a $T \times 1$ vector $\tau$ that will have the trend values of the serie
- We will have $T$ First order conditions, where its coefficients will form the $(T \times T)$   $\mathbf{F}$ matrix

```math
\boldsymbol{\tau}_{t} = (\mathbf{I}+\lambda \mathbf{F})^{-1} \mathbf{y}\\

\boldsymbol{\tau}_{t} = H\cdot y, \text{ where } H = (\mathbf{I}+\lambda \mathbf{F})^{-1}
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
The function `hp_filter` needs as inputs:
- a dataframe containing a column of dates and a column of data
- the lambda parameter

As output this function returns a list of three elements: 
- The first element is a dataframe that contains the original series and the estimated trend and cycle for the provided time span.
- The $H$ matrix
- $\lambda$ parameter

Using `plot` we can obtain the following graph:
![image](https://user-images.githubusercontent.com/103344273/194721925-b9c4e25b-8ada-4bc6-80b5-49644085654a.png)


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
This function needs as inputs:
- Data1: A dataframe of date and series that you know what the optimal lambda is
- Data2: A dataframe which you want to know the optimal lambda
- $\lambda$: this parameter is known for Data1
- rule: Which rule you want to accomplish? 

The output of this function is:
- A dataframe that contains the original series and the estimated trend and cycle for the provided time span.
- The $H$ matrix
- $\lambda$ parameter

At the file of `HP Marcet and Ravn.R`, you should have in the environment the `hp_filter` function.

You get the dataframe of the US GDP and another for Ecuador GDP. Then you execute the following code: 

```
> HP_ec <- hp_filter_MR(data1 = gdp_us, data2 = gdp_ec, lambda = 1600, rule = "rule 1")
```
In the console it will display the following result:

```
 Root finding: F(λ) - V = 0
 ----------------- 
 V: 0.000158265640468348 
 λ data2 under rule 1 is: 1966.4815 
 (F-V) = 1.66641873911022e-16
```

With base R plot function I made this graphs for Ecuador GDP using the Rule 1 provided by Marcet and Ravn. 

![image](https://user-images.githubusercontent.com/103344273/194731694-9d245fe5-4de8-48eb-a69a-41b5c213482c.png)

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
> HP_uk <- hp_filter_MR(data1 = gdp_us, data2 = gdp_uk, lambda = 1600, rule = "rule 2")
```

This is the console output:
```
 Root finding: F(λ) - W = 0
 ----------------- 
 W: 4.00653265974394e-08 
 λ for data2 under rule 2 is: 3647.6276 
 (F-V) = -1.41034956184781e-18 
```

Here are the results for the UK under rule 2

![image](https://user-images.githubusercontent.com/103344273/194731698-1bf5455c-7888-4da6-bca3-6909e923df20.png)

# One sided Hodrick Prescott Filter

Since it first appeareance was in Stock and Watson (1999), the one sided Hodrick Prescott filter tries to detrend the series into a trend and cyclical component. However it is Wolf, et. al (2020) who computes it without the Kalman Filter. 

As Wolf et. al (2020) state: 

> "This procedure is equivalent to applying HP-2s recursively on an expanding sample and
keeping, from each recursion step, only the estimate of the trend component for the latest
period. "

So we can easily achieve that with a for loop. However, Hodrick Prescott [2s] filter needs at least 5 observations to compute the filter, the first 4 observations come from the whole sample, and from the fifth observation we start applying the two sided filter recursively. 

As inputs the created function needs: 
- a dataframe with two vectors one of dates and one of data
- $\lambda$ parameter.

This function just has one output:
- a dataframe with 4 columns one for date, original data, trend, and cycle.

### Code results 

Here's the way to use the function. 
```
us_hp <- hp_filter1(data = gdp_us , lambda = 1600)
```

Using the plot function for the US gdp:

![image](https://user-images.githubusercontent.com/103344273/194733166-bcb2a759-e65f-4a91-9834-f12a1b42c230.png)



# References:
Marcet, A., & Morten, R. (2015, septiembre 14). The HP-Filter in Cross-Country Comparisons | Barcelona School of Economics Working Papers. Barcelona School of Economics. https://bse.eu/research/working-papers/hp-filter-cross-country-comparisons

Stock, J. H., & Watson, M. W. (1999). Forecasting inflation [Working Paper]. National Bureau of Economic Research. https://doi.org/10.3386/w7023

Wolf, E., Mokinski, F., & Schüler, Y. S. (2020). On adjusting the one-sided hodrick-prescott filter [SSRN Scholarly Paper]. https://doi.org/10.2139/ssrn.3536248
