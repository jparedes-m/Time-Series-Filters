# Unobserved Components Model

The unobserved components model is a semi-structural model since you give defined structure to the trend and the cycle as time series. For now I will treat two variants of the unobserved components model: constant, and stochastic drift. 

## Constant drift

Consider the following system of equations:

```math
\displaylines{y_{t} = \tau_{t} + c_{t}\\

\tau_{t} = \delta + \tau_{t-1} + \varepsilon_{t}\\

c_{t} = \rho_{1} c_{t-1} + \rho_{2} c_{t-2} + \mu_{t}\\

\varepsilon_{t} \sim \text{ iid } \mathcal{N} (0, \sigma^{2}_{\varepsilon})\\

\mu_{t} \sim \text{ iid }\mathcal{N}(0, \sigma^{2}_{\mu})}
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

### Code Results 

This function needs just need one compulsory input:
- the dataframe with two columns: date vector and data vector
- optional: `p_estimates = TRUE` if you want to print the estimates of the Kalman Filter. 


The output of the function is a list of three elements:
- a table containing the estimates of the Kalman Filter and the initial guesses (computed with Hodrick Prescott and ML)
- a dataframe containing the data, trend, drift and cycle (data_filter)
- a dataframe containing the smoothed (smoother)
```
us_ucm <- ucm_const(data = gdp_us, p_estimates = TRUE)
```

In the console it will print the following:

```
                 g_hp   kalman_est
phi1     7.803432e-01 9.737668e-01
phi2    -2.438401e-02 1.821226e-02
delta    6.666738e-03 6.660969e-03
sigma_e  4.381146e-06 1.533193e-05
sigma_w  1.056176e-04 1.134649e-04
```

Now we can plot the trend and cycle

![image](https://github.com/jparedes-m/Time-Series-Filters/assets/103344273/0b3c302c-7ad9-4e76-a913-6b399bf5e471)





## Stochastic Drift

Consider the following system:

```math
\displaylines{y_{t} = \tau_{t} + c_{t}\\
\tau_{t} = \delta_{t-1} + \tau_{t-1} + \varepsilon_{t}\\
\delta_{t} = \delta_{t-1} + \omega_{t}\\
c_{t} = \rho_{1} c_{t-1} \rho_{2} c_{t-2} + \mu_{t}\\
\varepsilon_{t} \sim \text{ iid } \mathcal{N}(0, \sigma^{2}_{\varepsilon})\\
\omega_{t} \sim \text{ iid } \mathcal{N}(0, \sigma^{2}_{\omega})\\
\mu_{t} \sim \text{ iid } \mathcal{N}(0, \sigma^{2}_{\mu})}
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

### Initial Guesses:
Consider that the initial guesses $\alpha_0$ and $P_{0}$ can be computed from the two sided Hodrick Prescott Filter.

```math
 \alpha_{0} = \begin{pmatrix}
    \tau_{2}\\\delta_{2} \\ c_{2} \\ c_{1}
    \end{pmatrix}_{HP} \quad \quad P_{0} = \begin{bmatrix}
    \gamma_{0}(\tau_{t}) & \gamma_{1}(\tau_{t}, \delta_{t}) & 0 &0 \\
    \gamma_{1}(\delta_{t}, \tau_{t}) & \gamma_{0}(\delta_{t}) &0 &0\\
    0 & 0 & \gamma_{0}(c_{t}) & \gamma_{1}(c_{t}, c_{t-1})\\
    0 & 0 & \gamma_{1}(c_{t-1}, c_{t}) & \gamma_{0}(c_{t-1})
    \end{bmatrix}_{HP}
```

### Code Results 

This function needs just need one compulsory input:
- the dataframe with two columns: date vector and data vector
- optional: `p_estimates = TRUE` if you want to print the estimates of the Kalman Filter. 


The output of the function is a list of three elements:
- a table containing the estimates of the Kalman Filter and the initial guesses (computed with Hodrick Prescott and ML)
- a dataframe containing the data, trend, drift and cycle (data_filter)
- a dataframe containing the smoothed (smoother)
```
us_ucm <- ucm_stoch(data = gdp_us, p_estimates = TRUE)
```

In the console it will print the following:

```
                 g_hp    kalman_est
phi1     7.803432e-01  9.223821e-01
phi2    -2.438401e-02 -3.000073e-03
sigma_e  4.615888e-07  5.077939e-07
sigma_w  3.505747e-08  8.773087e-08
sigma_u  1.056176e-04  1.056471e-04
```

Now we can plot the trend and cycle

![image](https://github.com/jparedes-m/Time-Series-Filters/assets/103344273/73f9b944-148e-4713-863b-6b0e6aca6ef3)


