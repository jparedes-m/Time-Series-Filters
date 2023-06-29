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

The linear filter function has three arguments: `data`, `p_estimates`, and `drift`.  `data` argument must be a dataframe with two columns (one of class `"Date"`, and the other that is numeric.). `p_estimates` is a boolean (FALSE as default) that if set to TRUE prints the summary of the regression in the console. Finally, `drift` is a boolean (TRUE as default) that if set to FALSE does not account for the constant term. 

This function returns a dataframe with four columns: `date`, `serie`, `trend` which are the predicted values of the regression, and `cycle` which can be seen as the residuals of the regression. For better filtering, in the input assure that the series is seasonally adjusted. 


These are the results of the US GDP filtering from 1970 to 2023:

![image](https://github.com/jparedes-m/Time-Series-Filters/assets/103344273/3fe23274-0063-4b8c-9cb6-605468515043)

Recall that these results are from the LinearFilterExample.R  file in this repo. 
