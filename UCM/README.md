# Unobserved Components Model

The unobserved components model is a semi-structural model since you give defined structure to the trend and the cycle as time series. For now I will treat two variants of the unobserved components model: constant, and stochastic drift. 

## Constant drift

Consider the following system of equations:

```math
y_{t} = \tau_{t} + c_{t}\\
\tau_{t} = \delta + \tau_{t-1} + \varepsilon_{t}\\
c_{t} = \rho_{1} c_{t-1} + \rho_{2} c_{t-2} + \mu_{t}\\
\varepsilon_{t} \sim \mathcal{N} (0, \sigma^{2}_{\varepsilon})\\
\mu_{t} \sim \mathcal{N}(0, \sigma^{2}_{\mu})
```
