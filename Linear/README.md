# Linear Filtering

Supposing you have a series that can be separated additively into trend and cycle:

```math
y_t = \tau_t + c_t
```

This may be the simplest filter there is in the repository. It uses an Ordinary Least Square estimator to get the trend and cycle.

In simple terms, a linear filter (with drift) computes the following regression. 

```math
y_t = \alpha + \beta \cdot t + \varepsilon
```

$\alpha$ is the drift, seen as the intersection between the trend and the vertical axis. $\beta$ is the slope of time ($t$) where time is defined as a sequence from $1$ to $N$, where $N$ is the number of periods. Finally, $\varepsilon$ is the error term from the regression.

The linear filter function has three arguments: `data`, `p_estimates`, and `drift`. 

The drift option 
