# Linear Filtering

Supposing you have a series that can be separated additively into trend and cycle:

```math
y_t = \tau_t + c_t
```

This may be the simplest filter there is in the repository. It uses an Ordinary Least Square estimator for getting the trend and cycle.

In simple terms, linear filter (with drift) computes the following regression. 

```math
y_t = \alpha + \beta \cdot t + \varepsilon
```

Where, $\alpha$ is the drift, seen as the intersection between the trend and the vertical axis. $\beta$ is the slope of time ($t$)

The linear filter function has three arguments: `data`, `p_estimates`, and `drift`. 

The drift option 
