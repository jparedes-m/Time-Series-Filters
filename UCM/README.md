# Unobserved Components Model

The unobserved components model is a semi-structural model since you give defined structure to the trend and the cycle as time series. For now I will treat two variants of the unobserved components model: constant, and stochastic drift. 

## Constant drift

Consider the following system of equations:

```math
y_{t} = \tau_{t} + c_{t}\\
\tau_{t} = \delta + \tau_{t-1} + \varepsilon_{t}\\
c_{t} = \rho_{1} c_{t-1} + \rho_{2} c_{t-2} + \mu_{t}\\
\begin{pmatrix}
\varepsilon_{t} \\ \mu_{t} 
\end{pmatrix} \sim \text{ i.i.d  } \mathcal{N} \left(0, \begin{pmatrix} \sigma^{2}_{\varepsilon} & \\ & \sigma^{2}_{\mu} \end{pmatrix}
```
