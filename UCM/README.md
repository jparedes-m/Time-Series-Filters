# Unobserved Components Model

The unobserved components model is a semi-structural model since you give defined structure to the trend and the cycle as time series. For now I will treat two variants of the unobserved components model: constant, and stochastic drift. 

## Constant drift

Consider the following system of equations:

```math
y_{t} = \tau_{t} + c_{t}\\
\tau_{t} = \delta + \tau_{t-1} + \varepsilon_{t}\\
c_{t} = \rho_{1} c_{t-1} + \rho_{2} c_{t-2} + \mu_{t}\\
\varepsilon_{t} \sim \text{ iid } \mathcal{N} (0, \sigma^{2}_{\varepsilon})\\
\mu_{t} \sim \text{ iid }\mathcal{N}(0, \sigma^{2}_{\mu})
```

So I use the State Space representation for estimating with the Kalman Filter. I use the package `FKF`, for the estimation. 

I change just a little the notation from the documentation of the `FKF` package

### State Equation: $\alpha_{t} = d_{t} + T_{t} \alpha_{t-1} + H_{t} \eta_{t}$

```math
\begin{pmatrix}
    \tau_{t} \\  c_{t} \\ c_{t-1}
    \end{pmatrix} = \begin{pmatrix}
    \delta \\0 \\0
    \end{pmatrix} + \begin{pmatrix}
    1 & 0 & 0 \\
    0 & \rho_{1} & \rho_{2}\\
    0 & 1 & 0
    \end{pmatrix}\begin{pmatrix}
    \tau_{t-1} \\ c_{t-1} \\ c_{t-2}+
    \end{pmatrix}+\begin{pmatrix}
    \sigma^{2}_{\epsilon} & 0 & 0 \\
    0 & \sigma^{2}_{\mu} & 0 \\
    0 & 0 & 0
    \end{pmatrix}\begin{pmatrix}
    \epsilon_{t} \\  \mu_{t} \\ 0
    \end{pmatrix}
```

Note that $H_{t}$ matrix should be filled with ones rather than the variances, I do this because the `FKF` package in R assumes that variances are 1, which is not entirely true.

### Measurement Equation: $y_{t} = c_{t} + Z_{t} \alpha_{t} + G_{t} \epsilon_{t}$

```math
 y_{t} = \begin{bmatrix} 0\end{bmatrix} + \begin{bmatrix} 1 & 1 & 0\end{bmatrix}\begin{pmatrix} \tau_{t} \\ c_{t} \\ c_{t-1} \end{pmatrix} + \begin{bmatrix} 0 \end{bmatrix} \epsilon_{t}
 ```
### Initial Guesses
Consider that the initial guesses $\alpha_{0}$, and $P_{0}$ can be computed from the two-sided Hodrick Prescott Filter.

```math
 \alpha_{0} =  \begin{pmatrix}
    \tau_{2} \\ c_{2} \\ c_{1}
    \end{pmatrix}_{HP} \quad 
    P_{0} = \begin{bmatrix} 
    \gamma_{0}(\tau_{t}) & 0 & 0 \\
    0 & \gamma_{0}(c_{t}) & \gamma_{1}(c_{t}, c_{t-1})\\
    0 & \gamma_{1}(c_{t-1}, c_{t}) & \gamma_{0}(c_{t-1})
    \end{bmatrix}_{HP}
```

## Stochastic Drift

Consider the following system:

```math
y_{t} = \tau_{t} + c_{t}\\
\tau_{t} = \delta_{t-1} + \tau_{t-1} + \varepsilon_{t}\\
\delta_{t} = \delta_{t-1} + \omega_{t}\\
c_{t} = \rho_{1} c_{t-1} \rho_{2} c_{t-2} + \mu_{t}\\
\varepsilon_{t} \sim \text{ iid } \mathcal{N}(0, \sigma^{2}_{\varepsilon}\\
\omega_{t} \sim \text{ iid } \mathcal{N}(0, \sigma^{2}_{\omega}\\
\mu_{t} \sim \text{ iid } \mathcal{N}(0, \sigma^{2}_{\mu}
```
Its state space representation can be written as:

### State Equation: $\alpha_{t} = d_{t} + T_{t} \alpha_{t-1} + H_{t}\eta_{t}$

```math
\begin{pmatrix}
    \tau_{t}\\\delta_{t} \\ c_{t} \\ c_{t-1}
    \end{pmatrix} =  \begin{pmatrix}
    0 \\ 0 \\0 \\0
    \end{pmatrix} + \begin{pmatrix}
    1 & 1 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & \rho_{1} & \rho_{2} \\
    0 & 0 & 1 & 0
    \end{pmatrix}\begin{pmatrix}
    \tau_{t-1}\\\delta_{t-1}\\c_{t-1}\\c_{t-2}
    \end{pmatrix}+\begin{pmatrix}
    \sigma^{2}_{\epsilon} & 0 & 0 &0 \\
    0 & \sigma^{2}_{\omega} & 0 &0 \\
    0 & 0 & \sigma^{2}_{\mu} & 0 \\
    0 & 0 & 0 & 0
    \end{pmatrix}\begin{pmatrix}
    \epsilon_{t} \\ \omega_{t}\\ \mu_{t} \\ 0
    \end{pmatrix}
```
Note that $H_{t}$ matrix should be filled with ones rather than the variances, I do this because the `FKF` package in R assumes that variances are 1, which is not entirely true.

### Measurement Equation: $y_{t} = c_{t} + Z_{t}\alpha_{t} + G_{t} \epsilon_{t}$

```math
 y_{t} = \begin{bmatrix} 0 \end{bmatrix} + \begin{bmatrix}1& 0 & 1 &0\end{bmatrix}\begin{pmatrix}
    \tau_{t}\\\delta_{t}\\c_{t}\\c_{t-1}
    \end{pmatrix}+\begin{bmatrix}0\end{bmatrix}\epsilon_{t}
```
