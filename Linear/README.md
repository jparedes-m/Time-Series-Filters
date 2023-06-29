# Linear Filtering

Supposing you have a series that can be separated additively into trend and cycle:

```math
y_t = \tau_t + c_t
```

This may be the simplest filter there is in the repository. It uses an Ordinary Least Square estimator for getting the trend and cycle.

For this filter, I use an OLS estimator. I programmed it to get a drift option (if desired)

So basically you run this regression: 
