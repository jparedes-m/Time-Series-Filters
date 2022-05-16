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

$\tau_{t}$ is the series trend while $c_{t}$ is the series cycle. In the context of macroeconomics, say $y_{t}$ is the GDP (output), therefore the trend will be the
*potential gdp* while the cycle will be the *output gap*

## The optimization problem 
So Hodrick - Prescott filters looks for a vector $\boldsymbol{\tau}$ that solves the following minimization problem:

```math
\min_{\boldsymbol{\tau}} \sum^{T}_{t=1}(y_{t}-\tau_{t})^{2}+\lambda\sum_{t=2}^{T-1}(\tau_{t+1}-2\tau_{t}+\tau_{t-1})^{2}
```
To solve this problem:
- Note that the arguments to minimize this expression are $`\tau_{1}, \cdots, \tau_{T}`$. Therefore you have $`T`$ **First Order Conditions**
- The coefficients of the F.O.C will give the T coefficients to build the following system of equations:

```\math
\boldsymbo{\tau}_{t} = (\mathbf{I}+\lambda \mathbf{F})^{-1} \mathbf{y}
```
where:
- $`\mathbf{I}`$ is the $T \times T$ identity matrix
- $`\mathbf{y}`$ is the original series
- $`\mathbf{F}`$ is a $T\times T$ pentadiagonal matrix

```\math
\mathbf{F}=\left[ \begin{matrix}
1 & -2 & 1 &  &  &  &  \\
-2 & 5 & -4 & 1 &  &  &    \\
1 & -4 & 6 & -4 & 1 &  &     \\ 
  &  \ddots & \ddots & \ddots & \ddots & \ddots &   \\
  &  & 1 & -4 & 6 & -4 & 1  \\
  &  &  & 1 & -4 & 5 & -2 \\
  &  &  &  & 1 & -2 & 1 \\
\end{matrix}
\right].
```
