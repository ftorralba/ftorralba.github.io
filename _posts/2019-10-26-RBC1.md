

```stata
clear all
```

# Real business cycle model

Today I am going to solve the basic Real Business Cycle model with Stata's new(ish) packages dsge and dsgenl. Why Stata? Why not! Yes, most people use Matlab/Octave, perhaps with Dynare, but I wanted to check out the functionality and overall user experience with Stata, while writing a jupyter notebook with the Stata kernel.

The model is a version of [King and Rebelo (1999)](#references). I will follow Stata's [DSGE Reference Manual,  Release 16](https://www.stata.com/bookstore/dynamic-stochastic-general-equilibrium-reference-manual/) (Intro 3b and 3e).

## Description of the model

In this model the representative household maximizes the welfare function
$$E_t \sum_{t=0}^\infty \beta^t \left( \frac{C^{1-\sigma}_t}{1-\sigma} - \frac{H^{1+\phi}_t}{1+\phi} \right)$$
subject to the resource constraint
$$P_t(C_t + I_t) \leq W_t H_t + R_t K_t + \Pi_t               \text{    (1)}$$ 

where $P$ is the general price level, $I$ is investment, $W$ is the nominal wage, $K$ is the capital stock, $R$ is the nominal return on capital, and $\Phi$ are profits. Since firms (below) operate in a perfectly competitive market, $\Pi_t=0$.

Capital accumulates via the usual transition equation
$$K_{t+1} = (1-\delta)K_t + I_t \text{    (2)}$$ 

The household maximizes its welfare function with respect to $C_t$, $H_t$, and $K_{t+1}$, and the solution is given by the Euler equation and the labor supply equation:

$$C^{-\sigma}_t = \beta E_t \left\{C^{-\sigma}_{t+1} \left[(1-\delta) + \frac{R_{t+1}}{P_{t+1}} \right] \right\} \text{    (3)}$$

$$C^{\sigma}_t H^{\phi}_t = \frac{W_t}{P_t} \text{    (4)}$$ 

Firms maximize profits using a Cobb-Douglas production function with labor and capital as inputs, and labor-augmenting (Harrod-neutral) technological progress, which follows a first-order autoregressive process:

$$Y_t = K^{\alpha}_t (Z_tH_t)^{1-\alpha} \text{    (5)}$$

$$Z_{t+1} = \rho_z Z_t + \epsilon^Z_{t+1} \text{    (6)}$$ 

$\epsilon^Z$ is a technology shock of the familiar type: $\epsilon^Z_t \sim N(0,\sigma_Z)$.

Here the example in the Stata Reference Manual departs a little from [King and Rebelo (1999)](#references), who use a random total factor productivity variable, and a deterministic labor-augmenting productivity variable. We can explore the implications of these alternative technology specifications later. 

The firm chooses the level of capital and labor, taking wages and the cost of capital as given, which produces the optimization conditions

$$\frac{R_t}{P_t} = \alpha \frac{Y_t}{K_t} \text{    (7)}$$ 

$$\frac{W_t}{P_t} = (1-\alpha) \frac{Y_t}{H_t} \text{    (8)}$$

The eight equations above completely describe the model. 


## Log-linearization of the model
Stata allows the option to solve the model either in linearized or non-linearized form. I'll try the former first. 

Before I proceed, two things. First, remember that the equations above describe the relationships among the model's variables *at a solution*, hence the levels of output, consumption, investment, etc. in those equations correspond to their steady-state values. Second, we set output as the numeraire, and hence $P=1$. With this normalization we should interpret $W$ and $R$ as the prices of labor and capital relative to output (i.e. the "real" price of factors).

We use [Uhlig's (1999)](#references) tricks to re-write the above equations in log-linearized form as:

$$y_t = \theta_i i_t + \theta_c c_t \text{    (1')}$$

$$k_{t+1} = (1-\delta)k_t + \delta i_t \text{    (2')}$$ 

$$E_t c_{t+1} - c_t = \frac{1-\beta + \beta\delta}{\sigma} E_tr_{t+1} \text{    (3')}$$

$$\sigma c_t + \phi h_t = w_t \text{    (4')}$$

$$y_t = \alpha k_t + (1-\alpha)(z_t + h_t) \text{    (5')}$$

$$z_{t+1} = \rho_z z_t + \epsilon^z_{t+1} \text{    (6')}$$

$$r_t = y_t - k_t \text{    (7')}$$

$$w_t = y_t - h_t \text{    (8')}$$

After the log-linearization, the (now lower-case) variables correspond to deviations from the steady state, $x_t = X_t - X^{SS}$: $c_t$, $r_t$, $w_t$, $k_t$, $i_t$, $y_t$, $h_t$ and $z_t$.

Notice the appearance of two new parameters in the first equation: $\theta_i$ and $\theta_c$. They are the output shares of investment and consumption in the steady state: $\theta_x = X^{SS}/Y^{SS}$.



## Solving the model

The first thing we can do is solve the model by assuming values for the preference and technology parameters. This means we don't estimate any parameters, and thus we make no use of data: "solving" simply means finding the set of state space equations that describe the solution to the model. The appropriate values for those parameters vary (and part of the fun is experimenting with different values and seeing how things change). I will start with values close to the ones chosen by the Stata DSGE Reference Manual:

| parameter | value   |
|-----------|---------|
| $\beta$   |   0.96  |
| $\sigma$  |   1     |
| $\phi$    |   1     |
| $\theta_i$|   0.3   |
| $\theta_c$|   0.7   |
| $\delta$  |   0.025 |
| $\alpha$  |   0.3   |
| $\rho_z$  |   0.8   |


STATA tells us to load a little dataset for this exercise:


```stata
use https://www.stata-press.com/data/r16/usmacro2
```

    (Federal Reserve Economic Data - St. Louis Fed, 2017-01-15)
    


```stata
desc
```

    
    Contains data from https://www.stata-press.com/data/r16/usmacro2.dta
      obs:           244                          Federal Reserve Economic Data - St. Louis Fed, 2017-01-15
     vars:            11                          1 May 2019 17:52
                                                  (_dta has notes)
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                  storage   display    value
    variable name   type    format     label      variable label
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    daten           int     %td                   Numeric (daily) date
    year            int     %9.0g                 Year
    quarter         byte    %9.0g                 Quarter
    dateq           int     %tq                   Date (quarters)
    y               double  %10.0g                Growth rate of real GDP (GDPC96)
    p               double  %10.0g                Growth rate of prices (GDPDEF)
    r               double  %10.0g                Federal funds rate (FEDFUNDS)
    c               double  %10.0g                Growth rate of consumption (PCECC96)
    n               double  %10.0g                Growth rate of hours worked (HOANBS)
    i               double  %10.0g                Corporate bond interest rate (AAA)
    e               double  %10.0g                Percentage change in US exchange rate (TWEXBMTH)
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sorted by: dateq
    


```stata
sum year quarter
```

    
        Variable |        Obs        Mean    Std. Dev.       Min        Max
    -------------+---------------------------------------------------------
            year |        244        1985    17.64301       1955       2015
         quarter |        244         2.5    1.120332          1          4
    

As you can see, the data set contains a handful of U.S. macro variables, at quartertely frequency, from 1955 to 2015.

Next, we write the model in a way that Stata can understand it, using the dsge command (for a linearized model). we write one equation for each variable in the model: $i$, $k_{t+1}$, $c$ (or $E_tc_{t+1}$), $h$, $y$, $w$, $r$, and the exogenous technology variable $z$. (In the Stata command below I have relabelled investment $i$ as $x$ to avoid any confusion with the variable i in the U.S. dataset (see above).) For each variable we tell Stata whether the variable is "unobserved", "state" or nothing (in which case Stata understand that it's an observed variable). A state variable is one whose period-to-period change is defined by a given linear equation. State variables can be deterministic, in which case we add "noshock" to the equation, or have a random component. In this model we specify output, $y$, as the observed variable, the demand components, labor quantity, and the factor prices as unobserved variables, the stock of capital as deterministic state variable, and the technology variable as the state variable with a random element:

dsge ({thetai}*x = y - {thetac}*c, unobserved) ///
     (F.k = {delta}*x+ (1-{delta})*k, state noshock) ///
     (c = F.c - ((1-{beta}+{beta}*{delta})/{sigma})*F.r, unobserved) ///
     ({phi}*h = w - {sigma}*c, unobserved) ///
     (y = (1-{alpha})*(z+h) + {alpha}*k) ///
     (w = y - h, unobserved) ///
     (r = y - k, unobserved) ///   
     (F.z = {rhoz}*z, state), ///
     from(sigma=1 beta=0.96 phi=1 alpha=0.3 delta=0.025 thetai=0.2 thetac=0.6 rhoz=0.8) ///
     solve 

Because we are solving the model, given certain parameters, we have the entire list of parameter values, and we add "solve" at the end.


```stata
dsge ({thetai}*x = y - {thetac}*c, unobserved) ///
     (F.k = {delta}*x+ (1-{delta})*k, state noshock) ///
     (c = F.c - ((1-{beta}+{beta}*{delta})/{sigma})*F.r, unobserved) ///
     ({phi}*h = w - {sigma}*c, unobserved) ///
     (y = (1-{alpha})*(z+h) + {alpha}*k) ///
     (w = y - h, unobserved) ///
     (r = y - k, unobserved) ///   
     (F.z = {rhoz}*z, state), ///
     from(sigma=1 beta=0.96 phi=1 alpha=0.3 delta=0.025 thetai=0.2 thetac=0.6 rhoz=0.8) ///
     solve 
```

    
    
    DSGE model
    
    Sample: 1955q1 - 2015q4                         Number of obs     =        244
    Log likelihood = -1970.8929
    ------------------------------------------------------------------------------
                 |                 OIM
               y |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
    /structural  |
          thetai |         .2          .        .       .            .           .
          thetac |         .6          .        .       .            .           .
           delta |       .025          .        .       .            .           .
            beta |        .96          .        .       .            .           .
           sigma |          1          .        .       .            .           .
             phi |          1          .        .       .            .           .
           alpha |         .3          .        .       .            .           .
            rhoz |         .8          .        .       .            .           .
    -------------+----------------------------------------------------------------
          sd(e.z)|          1          .                             .           .
    ------------------------------------------------------------------------------
    Note: Model solved at specified parameters.
    

The output is not terribly interesting: it's just giving us the parameter values that we had supplied in the first place. Notice that by default the technology shock variable has standard deviation = 1.

### Policy and transition matrices

Given the chosen parameter, we might be interested in knowing how a control variable responds to a unit increase in a state variable, in which case we type:


```stata
estat policy
```

    
    Policy matrix
    
    ------------------------------------------------------------------------------
                 |            Delta-method
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
    x            |
               k |  -1.972861          .        .       .            .           .
               z |   4.313127          .        .       .            .           .
    -------------+----------------------------------------------------------------
    c            |
               k |   .7519892          .        .       .            .           .
               z |   .1882345          .        .       .            .           .
    -------------+----------------------------------------------------------------
    h            |
               k |   -.347684          .        .       .            .           .
               z |   .3936658          .        .       .            .           .
    -------------+----------------------------------------------------------------
    y            |
               k |   .0566212          .        .       .            .           .
               z |   .9755661          .        .       .            .           .
    -------------+----------------------------------------------------------------
    w            |
               k |   .4043052          .        .       .            .           .
               z |   .5819003          .        .       .            .           .
    -------------+----------------------------------------------------------------
    r            |
               k |  -.9433788          .        .       .            .           .
               z |   .9755661          .        .       .            .           .
    ------------------------------------------------------------------------------
    Note: Standard errors reported as missing for constrained policy matrix values.
    

Investment is most sensitive to a unit technology shock, whereas consumption is the least sensitive. The return to capital declines when the capital stock increases, and so does investment (as one would expect in order to bring the system back to steady state). We can also print the transition matrix, which shows the parameters that rule the dynamics of $k_t$ and $z_t$:


```stata
estat transition
```

    
    Transition matrix of state variables
    
    ------------------------------------------------------------------------------
                 |            Delta-method
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
    F.k          |
               k |   .9256785          .        .       .            .           .
               z |   .1078282          .        .       .            .           .
    -------------+----------------------------------------------------------------
    F.z          |
               k |          0  (omitted)
               z |         .8          .        .       .            .           .
    ------------------------------------------------------------------------------
    Note: Standard errors reported as missing for constrained transition matrix values.
    

### Impulse-response functions and comparison of dynamics with different parameter values

I follow the Stata manual and compare the impulse-responses of this model with those of a model where the technology variable is less persistent.


```stata
irf set rbcirf
irf create persistent
irf graph irf, irf(persistent) impulse(z) response(y c x h w z) byopts(yrescale) noci

```

    
    (file rbcirf.irf now active)
    
    (file rbcirf.irf updated)
    


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;222.14&quot; width=&quot;1053.98&quot; height=&quot;1022.08&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;223.87&quot; width=&quot;1050.52&quot; height=&quot;1018.63&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1106.61&quot; x2=&quot;1342.32&quot; y2=&quot;1106.61&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;683.37&quot; x2=&quot;1342.32&quot; y2=&quot;683.37&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;260.13&quot; x2=&quot;1342.32&quot; y2=&quot;260.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;222.14&quot; width=&quot;1053.98&quot; height=&quot;1022.08&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;223.87&quot; width=&quot;1050.52&quot; height=&quot;1018.63&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1111.06&quot; x2=&quot;2583.41&quot; y2=&quot;1111.06&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;898.33&quot; x2=&quot;2583.41&quot; y2=&quot;898.33&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;685.60&quot; x2=&quot;2583.41&quot; y2=&quot;685.60&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;472.86&quot; x2=&quot;2583.41&quot; y2=&quot;472.86&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;260.13&quot; x2=&quot;2583.41&quot; y2=&quot;260.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;222.14&quot; width=&quot;1053.98&quot; height=&quot;1022.08&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;223.87&quot; width=&quot;1050.52&quot; height=&quot;1018.63&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1206.11&quot; x2=&quot;3824.49&quot; y2=&quot;1206.11&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;969.61&quot; x2=&quot;3824.49&quot; y2=&quot;969.61&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;733.12&quot; x2=&quot;3824.49&quot; y2=&quot;733.12&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;496.62&quot; x2=&quot;3824.49&quot; y2=&quot;496.62&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;260.13&quot; x2=&quot;3824.49&quot; y2=&quot;260.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;1400.03&quot; width=&quot;1053.98&quot; height=&quot;1022.08&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;1401.76&quot; width=&quot;1050.52&quot; height=&quot;1018.63&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2384.12&quot; x2=&quot;1342.32&quot; y2=&quot;2384.12&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2164.70&quot; x2=&quot;1342.32&quot; y2=&quot;2164.70&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1945.41&quot; x2=&quot;1342.32&quot; y2=&quot;1945.41&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1726.12&quot; x2=&quot;1342.32&quot; y2=&quot;1726.12&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1506.70&quot; x2=&quot;1342.32&quot; y2=&quot;1506.70&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;1400.03&quot; width=&quot;1053.98&quot; height=&quot;1022.08&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;1401.76&quot; width=&quot;1050.52&quot; height=&quot;1018.63&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2362.96&quot; x2=&quot;2583.41&quot; y2=&quot;2362.96&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2131.66&quot; x2=&quot;2583.41&quot; y2=&quot;2131.66&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1900.49&quot; x2=&quot;2583.41&quot; y2=&quot;1900.49&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1669.32&quot; x2=&quot;2583.41&quot; y2=&quot;1669.32&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1438.02&quot; x2=&quot;2583.41&quot; y2=&quot;1438.02&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;1400.03&quot; width=&quot;1053.98&quot; height=&quot;1022.08&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;1401.76&quot; width=&quot;1050.52&quot; height=&quot;1018.63&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2347.49&quot; x2=&quot;3824.49&quot; y2=&quot;2347.49&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2120.15&quot; x2=&quot;3824.49&quot; y2=&quot;2120.15&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1892.69&quot; x2=&quot;3824.49&quot; y2=&quot;1892.69&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1665.36&quot; x2=&quot;3824.49&quot; y2=&quot;1665.36&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1438.02&quot; x2=&quot;3824.49&quot; y2=&quot;1438.02&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;path d=&quot; M326.33 1206.23 L448.59 838.56 L570.74 595.38 L693.00 448.11 L815.26 373.98 L937.53 355.17 L1059.80 377.57 L1182.06 430.04 L1304.20 504.30&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;path d=&quot; M1567.42 273.62 L1689.68 520.88 L1811.82 712.70 L1934.09 860.71 L2056.35 974.07 L2178.62 1060.08 L2300.88 1124.43 L2423.15 1171.95 L2545.29 1206.23&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;path d=&quot; M2808.51 302.95 L2930.77 475.09 L3053.04 620.50 L3175.18 743.88 L3297.44 849.20 L3419.71 939.54 L3541.97 1017.38 L3664.24 1084.83 L3786.50 1143.73&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;path d=&quot; M326.33 1438.14 L448.59 1674.02 L570.74 1859.15 L693.00 2004.19 L815.26 2117.18 L937.53 2204.80 L1059.80 2272.37 L1182.06 2324.10 L1304.20 2363.33&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;path d=&quot; M1567.42 1466.36 L1689.68 1684.91 L1811.82 1860.27 L1934.09 2000.98 L2056.35 2113.96 L2178.62 2204.92 L2300.88 2277.94 L2423.15 2336.72 L2545.29 2384.12&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;path d=&quot; M2808.51 1438.14 L2930.77 1665.48 L3053.04 1847.40 L3175.18 1992.81 L3297.44 2109.26 L3419.71 2202.32 L3541.97 2276.82 L3664.24 2336.47 L3786.50 2384.12&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1244.22&quot; x2=&quot;288.34&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1106.61&quot; x2=&quot;264.33&quot; y2=&quot;1106.61&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1127.65&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;683.37&quot; x2=&quot;264.33&quot; y2=&quot;683.37&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;704.41&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.25&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;260.13&quot; x2=&quot;264.33&quot; y2=&quot;260.13&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;281.17&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1244.22&quot; x2=&quot;1529.43&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1111.06&quot; x2=&quot;1505.42&quot; y2=&quot;1111.06&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1132.10&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;898.33&quot; x2=&quot;1505.42&quot; y2=&quot;898.33&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;919.37&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.1&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;685.60&quot; x2=&quot;1505.42&quot; y2=&quot;685.60&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;706.63&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;472.86&quot; x2=&quot;1505.42&quot; y2=&quot;472.86&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;493.90&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;260.13&quot; x2=&quot;1505.42&quot; y2=&quot;260.13&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;281.17&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1244.22&quot; x2=&quot;2770.51&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1206.11&quot; x2=&quot;2746.51&quot; y2=&quot;1206.11&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1227.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;969.61&quot; x2=&quot;2746.51&quot; y2=&quot;969.61&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;990.65&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;733.12&quot; x2=&quot;2746.51&quot; y2=&quot;733.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;754.16&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;496.62&quot; x2=&quot;2746.51&quot; y2=&quot;496.62&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;517.66&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;260.13&quot; x2=&quot;2746.51&quot; y2=&quot;260.13&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;281.17&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2422.11&quot; x2=&quot;288.34&quot; y2=&quot;1400.03&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2384.12&quot; x2=&quot;264.33&quot; y2=&quot;2384.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;2405.16&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2164.70&quot; x2=&quot;264.33&quot; y2=&quot;2164.70&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;2185.74&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1945.41&quot; x2=&quot;264.33&quot; y2=&quot;1945.41&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1966.45&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1726.12&quot; x2=&quot;264.33&quot; y2=&quot;1726.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1747.16&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1506.70&quot; x2=&quot;264.33&quot; y2=&quot;1506.70&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1527.74&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2422.11&quot; x2=&quot;1529.43&quot; y2=&quot;1400.03&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2362.96&quot; x2=&quot;1505.42&quot; y2=&quot;2362.96&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;2383.99&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2131.66&quot; x2=&quot;1505.42&quot; y2=&quot;2131.66&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;2152.70&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1900.49&quot; x2=&quot;1505.42&quot; y2=&quot;1900.49&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1921.53&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1669.32&quot; x2=&quot;1505.42&quot; y2=&quot;1669.32&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1690.35&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.8&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1438.02&quot; x2=&quot;1505.42&quot; y2=&quot;1438.02&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1459.06&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2422.11&quot; x2=&quot;2770.51&quot; y2=&quot;1400.03&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2347.49&quot; x2=&quot;2746.51&quot; y2=&quot;2347.49&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;2368.53&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2120.15&quot; x2=&quot;2746.51&quot; y2=&quot;2120.15&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;2141.19&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1892.69&quot; x2=&quot;2746.51&quot; y2=&quot;1892.69&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1913.73&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1665.36&quot; x2=&quot;2746.51&quot; y2=&quot;1665.36&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1686.39&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.8&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1438.02&quot; x2=&quot;2746.51&quot; y2=&quot;1438.02&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1459.06&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2422.11&quot; x2=&quot;1342.32&quot; y2=&quot;2422.11&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;326.33&quot; y1=&quot;2422.11&quot; x2=&quot;326.33&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;326.33&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;570.86&quot; y1=&quot;2422.11&quot; x2=&quot;570.86&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;570.86&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;815.39&quot; y1=&quot;2422.11&quot; x2=&quot;815.39&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;815.39&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;1059.80&quot; y1=&quot;2422.11&quot; x2=&quot;1059.80&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1059.80&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;1304.33&quot; y1=&quot;2422.11&quot; x2=&quot;1304.33&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1304.33&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2422.11&quot; x2=&quot;2583.41&quot; y2=&quot;2422.11&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1567.42&quot; y1=&quot;2422.11&quot; x2=&quot;1567.42&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1567.42&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1811.95&quot; y1=&quot;2422.11&quot; x2=&quot;1811.95&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1811.95&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2056.48&quot; y1=&quot;2422.11&quot; x2=&quot;2056.48&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2056.48&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;2300.88&quot; y1=&quot;2422.11&quot; x2=&quot;2300.88&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2300.88&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;2545.41&quot; y1=&quot;2422.11&quot; x2=&quot;2545.41&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2545.41&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2422.11&quot; x2=&quot;3824.49&quot; y2=&quot;2422.11&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2808.51&quot; y1=&quot;2422.11&quot; x2=&quot;2808.51&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2808.51&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3053.04&quot; y1=&quot;2422.11&quot; x2=&quot;3053.04&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3053.04&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;3297.57&quot; y1=&quot;2422.11&quot; x2=&quot;3297.57&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3297.57&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;3541.97&quot; y1=&quot;2422.11&quot; x2=&quot;3541.97&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3541.97&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;3786.50&quot; y1=&quot;2422.11&quot; x2=&quot;3786.50&quot; y2=&quot;2446.12&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3786.50&quot; y=&quot;2500.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;135.39&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;137.11&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;815.39&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;persistent, z, c&lt;/text&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;135.39&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;137.11&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2056.48&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;persistent, z, h&lt;/text&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;135.39&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;137.11&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3297.57&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;persistent, z, w&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;1313.28&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;1315.00&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;815.39&quot; y=&quot;1369.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;persistent, z, x&lt;/text&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;1313.28&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;1315.00&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2056.48&quot; y=&quot;1369.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;persistent, z, y&lt;/text&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;1313.28&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;1315.00&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3297.57&quot; y=&quot;1369.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;persistent, z, z&lt;/text&gt;
	&lt;text x=&quot;1980.00&quot; y=&quot;2634.61&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;text x=&quot;118.06&quot; y=&quot;2737.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:79.94px;fill:#000000&quot;&gt;Graphs by irfname, impulse, and response&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



    
    
    
    

Now what happens if the technology persistence, $\rho_z$, is 0.6 instead of 0.8?


```stata
quietly dsge ({thetai}*x = y - {thetac}*c, unobserved) ///
     (F.k = {delta}*x+ (1-{delta})*k, state noshock) ///
     (c = F.c - ((1-{beta}+{beta}*{delta})/{sigma})*F.r, unobserved) ///
     ({phi}*h = w - {sigma}*c, unobserved) ///
     (y = (1-{alpha})*(z+h) + {alpha}*k) ///
     (w = y - h, unobserved) ///
     (r = y - k, unobserved) ///   
     (F.z = {rhoz}*z, state), ///
     from(sigma=1 beta=0.96 phi=1 alpha=0.3 delta=0.025 thetai=0.2 thetac=0.6 rhoz=0.6) ///
     solve 
```


```stata
irf create transitory, replace
irf graph irf, irf(transitory) impulse(z) response(y c x h w z) byopts(yrescale)
```

    
    irfname transitory not found in rbcirf.irf
    (file rbcirf.irf updated)
    


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;222.14&quot; width=&quot;1053.98&quot; height=&quot;894.25&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;223.87&quot; width=&quot;1050.52&quot; height=&quot;890.79&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1031.74&quot; x2=&quot;1342.32&quot; y2=&quot;1031.74&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;673.10&quot; x2=&quot;1342.32&quot; y2=&quot;673.10&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;314.33&quot; x2=&quot;1342.32&quot; y2=&quot;314.33&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;222.14&quot; width=&quot;1053.98&quot; height=&quot;894.25&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;223.87&quot; width=&quot;1050.52&quot; height=&quot;890.79&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;983.35&quot; x2=&quot;2583.41&quot; y2=&quot;983.35&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;649.09&quot; x2=&quot;2583.41&quot; y2=&quot;649.09&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;314.83&quot; x2=&quot;2583.41&quot; y2=&quot;314.83&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;222.14&quot; width=&quot;1053.98&quot; height=&quot;894.25&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;223.87&quot; width=&quot;1050.52&quot; height=&quot;890.79&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1078.39&quot; x2=&quot;3824.49&quot; y2=&quot;1078.39&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;805.64&quot; x2=&quot;3824.49&quot; y2=&quot;805.64&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;532.88&quot; x2=&quot;3824.49&quot; y2=&quot;532.88&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;260.13&quot; x2=&quot;3824.49&quot; y2=&quot;260.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;1272.31&quot; width=&quot;1053.98&quot; height=&quot;894.25&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;1274.04&quot; width=&quot;1050.52&quot; height=&quot;890.79&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2091.19&quot; x2=&quot;1342.32&quot; y2=&quot;2091.19&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1830.94&quot; x2=&quot;1342.32&quot; y2=&quot;1830.94&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1570.56&quot; x2=&quot;1342.32&quot; y2=&quot;1570.56&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;1272.31&quot; width=&quot;1053.98&quot; height=&quot;894.25&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;1274.04&quot; width=&quot;1050.52&quot; height=&quot;890.79&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2128.57&quot; x2=&quot;2583.41&quot; y2=&quot;2128.57&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1720.67&quot; x2=&quot;2583.41&quot; y2=&quot;1720.67&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1312.66&quot; x2=&quot;2583.41&quot; y2=&quot;1312.66&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;1272.31&quot; width=&quot;1053.98&quot; height=&quot;894.25&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;1274.04&quot; width=&quot;1050.52&quot; height=&quot;890.79&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2128.57&quot; x2=&quot;3824.49&quot; y2=&quot;2128.57&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1719.44&quot; x2=&quot;3824.49&quot; y2=&quot;1719.44&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1310.31&quot; x2=&quot;3824.49&quot; y2=&quot;1310.31&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_PnHqpJYA&quot; d=&quot;M326.33 1078.39 L326.33 1078.39 L448.59 512.34 L570.86 287.98 L693.12 260.13 L815.39 342.18 L937.53 483.01 L1059.80 652.06 L1182.06 832.00 L1304.33 1012.43 L1304.33 1012.43 L1182.06 832.00 L1059.80 652.06 L937.53 483.01 L815.39 342.18 L693.12 260.13 L570.86 287.98 L448.59 512.34 L326.33 1078.39 L326.33 1078.39 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_akvSgome&quot;&gt;
			&lt;use xlink:href=&quot;#path_PnHqpJYA&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_PnHqpJYA&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_PnHqpJYA&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_akvSgome)&quot;/&gt;
	&lt;path d=&quot; M326.33 1078.51 L448.59 512.34 L570.74 287.98 L693.00 260.13 L815.26 342.30 L937.53 483.01 L1059.80 652.18 L1182.06 832.00 L1304.20 1012.55&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_oQVdrPlv&quot; d=&quot;M1567.42 260.13 L1567.42 260.13 L1689.68 616.29 L1811.95 824.94 L1934.21 945.60 L2056.48 1013.67 L2178.62 1050.67 L2300.88 1069.11 L2423.15 1076.91 L2545.41 1078.39 L2545.41 1078.39 L2423.15 1076.91 L2300.88 1069.11 L2178.62 1050.67 L2056.48 1013.67 L1934.21 945.60 L1811.95 824.94 L1689.68 616.29 L1567.42 260.13 L1567.42 260.13 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_zTtmdMAq&quot;&gt;
			&lt;use xlink:href=&quot;#path_oQVdrPlv&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_oQVdrPlv&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_oQVdrPlv&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_zTtmdMAq)&quot;/&gt;
	&lt;path d=&quot; M1567.42 260.13 L1689.68 616.29 L1811.82 825.07 L1934.09 945.60 L2056.35 1013.79 L2178.62 1050.67 L2300.88 1069.23 L2423.15 1076.91 L2545.29 1078.51&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_mfbZizIr&quot; d=&quot;M2808.51 300.85 L2808.51 300.85 L2930.77 548.35 L3053.04 701.68 L3175.30 797.96 L3297.57 859.84 L3419.83 900.68 L3541.97 928.65 L3664.24 948.57 L3786.50 963.55 L3786.50 963.55 L3664.24 948.57 L3541.97 928.65 L3419.83 900.68 L3297.57 859.84 L3175.30 797.96 L3053.04 701.68 L2930.77 548.35 L2808.51 300.85 L2808.51 300.85 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_oMstIhhd&quot;&gt;
			&lt;use xlink:href=&quot;#path_mfbZizIr&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_mfbZizIr&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_mfbZizIr&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_oMstIhhd)&quot;/&gt;
	&lt;path d=&quot; M2808.51 300.85 L2930.77 548.48 L3053.04 701.68 L3175.18 798.09 L3297.44 859.84 L3419.71 900.68 L3541.97 928.65 L3664.24 948.70 L3786.50 963.67&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_JtqNzGbd&quot; d=&quot;M326.33 1492.22 L326.33 1492.22 L448.59 1761.39 L570.86 1920.66 L693.12 2014.22 L815.39 2068.42 L937.53 2099.24 L1059.80 2116.19 L1182.06 2124.73 L1304.33 2128.57 L1304.33 2128.57 L1182.06 2124.73 L1059.80 2116.19 L937.53 2099.24 L815.39 2068.42 L693.12 2014.22 L570.86 1920.66 L448.59 1761.39 L326.33 1492.22 L326.33 1492.22 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_bWUmcDON&quot;&gt;
			&lt;use xlink:href=&quot;#path_JtqNzGbd&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_JtqNzGbd&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_JtqNzGbd&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_bWUmcDON)&quot;/&gt;
	&lt;path d=&quot; M326.33 1492.22 L448.59 1761.39 L570.74 1920.78 L693.00 2014.22 L815.26 2068.55 L937.53 2099.36 L1059.80 2116.19 L1182.06 2124.85 L1304.20 2128.69&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_LVUXexSO&quot; d=&quot;M1567.42 1310.31 L1567.42 1310.31 L1689.68 1632.31 L1811.95 1825.86 L1934.21 1942.44 L2056.48 2012.61 L2178.62 2055.06 L2300.88 2080.92 L2423.15 2096.64 L2545.41 2106.29 L2545.41 2106.29 L2423.15 2096.64 L2300.88 2080.92 L2178.62 2055.06 L2056.48 2012.61 L1934.21 1942.44 L1811.95 1825.86 L1689.68 1632.31 L1567.42 1310.31 L1567.42 1310.31 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_CrXIAsTj&quot;&gt;
			&lt;use xlink:href=&quot;#path_LVUXexSO&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_LVUXexSO&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_LVUXexSO&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_CrXIAsTj)&quot;/&gt;
	&lt;path d=&quot; M1567.42 1310.31 L1689.68 1632.31 L1811.82 1825.99 L1934.09 1942.44 L2056.35 2012.73 L2178.62 2055.18 L2300.88 2080.92 L2423.15 2096.64 L2545.29 2106.41&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_UxUotZvX&quot; d=&quot;M2808.51 1310.31 L2808.51 1310.31 L2930.77 1637.63 L3053.04 1834.03 L3175.30 1951.85 L3297.57 2022.51 L3419.83 2064.96 L3541.97 2090.45 L3664.24 2105.67 L3786.50 2114.83 L3786.50 2114.83 L3664.24 2105.67 L3541.97 2090.45 L3419.83 2064.96 L3297.57 2022.51 L3175.30 1951.85 L3053.04 1834.03 L2930.77 1637.63 L2808.51 1310.31 L2808.51 1310.31 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_MBBbowPm&quot;&gt;
			&lt;use xlink:href=&quot;#path_UxUotZvX&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_UxUotZvX&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_UxUotZvX&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_MBBbowPm)&quot;/&gt;
	&lt;path d=&quot; M2808.51 1310.31 L2930.77 1637.63 L3053.04 1834.03 L3175.18 1951.85 L3297.44 2022.63 L3419.71 2064.96 L3541.97 2090.45 L3664.24 2105.67 L3786.50 2114.83&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1116.38&quot; x2=&quot;288.34&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1031.74&quot; x2=&quot;264.33&quot; y2=&quot;1031.74&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1052.77&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.14&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;673.10&quot; x2=&quot;264.33&quot; y2=&quot;673.10&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;694.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.16&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;314.33&quot; x2=&quot;264.33&quot; y2=&quot;314.33&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;335.37&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.18&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1116.38&quot; x2=&quot;1529.43&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;983.35&quot; x2=&quot;1505.42&quot; y2=&quot;983.35&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1004.39&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;649.09&quot; x2=&quot;1505.42&quot; y2=&quot;649.09&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;670.13&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;314.83&quot; x2=&quot;1505.42&quot; y2=&quot;314.83&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;335.87&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1116.38&quot; x2=&quot;2770.51&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1078.39&quot; x2=&quot;2746.51&quot; y2=&quot;1078.39&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1099.43&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;805.64&quot; x2=&quot;2746.51&quot; y2=&quot;805.64&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;826.68&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;532.88&quot; x2=&quot;2746.51&quot; y2=&quot;532.88&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;553.92&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;260.13&quot; x2=&quot;2746.51&quot; y2=&quot;260.13&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;281.17&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2166.56&quot; x2=&quot;288.34&quot; y2=&quot;1272.31&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2091.19&quot; x2=&quot;264.33&quot; y2=&quot;2091.19&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;2112.23&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1830.94&quot; x2=&quot;264.33&quot; y2=&quot;1830.94&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1851.98&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1570.56&quot; x2=&quot;264.33&quot; y2=&quot;1570.56&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1591.60&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1310.31&quot; x2=&quot;264.33&quot; y2=&quot;1310.31&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1331.34&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2166.56&quot; x2=&quot;1529.43&quot; y2=&quot;1272.31&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2128.57&quot; x2=&quot;1505.42&quot; y2=&quot;2128.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;2149.60&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1720.67&quot; x2=&quot;1505.42&quot; y2=&quot;1720.67&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1741.71&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;1312.66&quot; x2=&quot;1505.42&quot; y2=&quot;1312.66&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1493.41&quot; y=&quot;1333.70&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2166.56&quot; x2=&quot;2770.51&quot; y2=&quot;1272.31&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2128.57&quot; x2=&quot;2746.51&quot; y2=&quot;2128.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;2149.60&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1719.44&quot; x2=&quot;2746.51&quot; y2=&quot;1719.44&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1740.47&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;1310.31&quot; x2=&quot;2746.51&quot; y2=&quot;1310.31&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2734.50&quot; y=&quot;1331.34&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2166.56&quot; x2=&quot;1342.32&quot; y2=&quot;2166.56&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;326.33&quot; y1=&quot;2166.56&quot; x2=&quot;326.33&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;326.33&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;570.86&quot; y1=&quot;2166.56&quot; x2=&quot;570.86&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;570.86&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;815.39&quot; y1=&quot;2166.56&quot; x2=&quot;815.39&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;815.39&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;1059.80&quot; y1=&quot;2166.56&quot; x2=&quot;1059.80&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1059.80&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;1304.33&quot; y1=&quot;2166.56&quot; x2=&quot;1304.33&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1304.33&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;line x1=&quot;1529.43&quot; y1=&quot;2166.56&quot; x2=&quot;2583.41&quot; y2=&quot;2166.56&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1567.42&quot; y1=&quot;2166.56&quot; x2=&quot;1567.42&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1567.42&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1811.95&quot; y1=&quot;2166.56&quot; x2=&quot;1811.95&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1811.95&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2056.48&quot; y1=&quot;2166.56&quot; x2=&quot;2056.48&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2056.48&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;2300.88&quot; y1=&quot;2166.56&quot; x2=&quot;2300.88&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2300.88&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;2545.41&quot; y1=&quot;2166.56&quot; x2=&quot;2545.41&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2545.41&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;line x1=&quot;2770.51&quot; y1=&quot;2166.56&quot; x2=&quot;3824.49&quot; y2=&quot;2166.56&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2808.51&quot; y1=&quot;2166.56&quot; x2=&quot;2808.51&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2808.51&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3053.04&quot; y1=&quot;2166.56&quot; x2=&quot;3053.04&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3053.04&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;3297.57&quot; y1=&quot;2166.56&quot; x2=&quot;3297.57&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3297.57&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;3541.97&quot; y1=&quot;2166.56&quot; x2=&quot;3541.97&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3541.97&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;3786.50&quot; y1=&quot;2166.56&quot; x2=&quot;3786.50&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3786.50&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;135.39&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;137.11&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;815.39&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;transitory, z, c&lt;/text&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;135.39&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;137.11&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2056.48&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;transitory, z, h&lt;/text&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;135.39&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;137.11&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3297.57&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;transitory, z, w&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;1185.56&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;1187.29&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;815.39&quot; y=&quot;1242.07&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;transitory, z, x&lt;/text&gt;
	&lt;rect x=&quot;1529.43&quot; y=&quot;1185.56&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;1531.15&quot; y=&quot;1187.29&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2056.48&quot; y=&quot;1242.07&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;transitory, z, y&lt;/text&gt;
	&lt;rect x=&quot;2770.51&quot; y=&quot;1185.56&quot; width=&quot;1053.98&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;2772.24&quot; y=&quot;1187.29&quot; width=&quot;1050.52&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3297.57&quot; y=&quot;1242.07&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;transitory, z, z&lt;/text&gt;
	&lt;rect x=&quot;567.64&quot; y=&quot;2443.64&quot; width=&quot;2824.59&quot; height=&quot;186.50&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;570.52&quot; y=&quot;2446.52&quot; width=&quot;2818.83&quot; height=&quot;180.74&quot; style=&quot;fill:none;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;610.83&quot; y=&quot;2486.96&quot; width=&quot;374.34&quot; height=&quot;99.99&quot; style=&quot;fill:#C0C0C0&quot;/&gt;
	&lt;rect x=&quot;615.15&quot; y=&quot;2491.28&quot; width=&quot;365.70&quot; height=&quot;91.35&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;1526.46&quot; y1=&quot;2536.95&quot; x2=&quot;1900.80&quot; y2=&quot;2536.95&quot; style=&quot;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;text x=&quot;1045.19&quot; y=&quot;2571.98&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;95% CI&lt;/text&gt;
	&lt;text x=&quot;1960.82&quot; y=&quot;2571.98&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;impulse-response function (irf)&lt;/text&gt;
	&lt;text x=&quot;1980.00&quot; y=&quot;2379.06&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;text x=&quot;118.06&quot; y=&quot;2737.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:79.94px;fill:#000000&quot;&gt;Graphs by irfname, impulse, and response&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



    
    
    
    

It's hard to compare the impulse-responses across separate sets of charts, so let's overlay the impulse-response functions for, say, consumption:


```stata
irf ograph (persistent z c irf) (transitory z c irf)
```


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;290.81&quot; y=&quot;100.86&quot; width=&quot;3568.21&quot; height=&quot;1992.81&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;293.69&quot; y=&quot;103.74&quot; width=&quot;3562.45&quot; height=&quot;1987.05&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1885.76&quot; x2=&quot;3859.02&quot; y2=&quot;1885.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1311.91&quot; x2=&quot;3859.02&quot; y2=&quot;1311.91&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;738.07&quot; x2=&quot;3859.02&quot; y2=&quot;738.07&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;164.22&quot; x2=&quot;3859.02&quot; y2=&quot;164.22&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;path d=&quot; M354.05 1446.93 L784.33 948.45 L1214.48 618.77 L1644.64 419.03 L2074.92 318.54 L2505.07 292.93 L2935.23 323.37 L3365.51 394.65 L3795.66 495.14&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;path d=&quot; M354.05 2030.43 L784.33 1668.08 L1214.48 1524.52 L1644.64 1506.70 L2074.92 1559.30 L2505.07 1649.39 L2935.23 1757.55 L3365.51 1872.64 L3795.66 1988.23&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2093.67&quot; x2=&quot;290.81&quot; y2=&quot;100.86&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1885.76&quot; x2=&quot;250.84&quot; y2=&quot;1885.76&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1885.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1885.76)&quot; text-anchor=&quot;middle&quot;&gt;.15&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1311.91&quot; x2=&quot;250.84&quot; y2=&quot;1311.91&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1311.91&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1311.91)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;738.07&quot; x2=&quot;250.84&quot; y2=&quot;738.07&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;738.07&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,738.07)&quot; text-anchor=&quot;middle&quot;&gt;.25&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;164.22&quot; x2=&quot;250.84&quot; y2=&quot;164.22&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;164.22&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,164.22)&quot; text-anchor=&quot;middle&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2093.67&quot; x2=&quot;3859.02&quot; y2=&quot;2093.67&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;354.17&quot; y1=&quot;2093.67&quot; x2=&quot;354.17&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;354.17&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1214.48&quot; y1=&quot;2093.67&quot; x2=&quot;1214.48&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1214.48&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2074.92&quot; y1=&quot;2093.67&quot; x2=&quot;2074.92&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2074.92&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;2935.35&quot; y1=&quot;2093.67&quot; x2=&quot;2935.35&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2935.35&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;3795.66&quot; y1=&quot;2093.67&quot; x2=&quot;3795.66&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3795.66&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;2074.92&quot; y=&quot;2333.64&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;1317.32&quot; y=&quot;2418.27&quot; width=&quot;1515.32&quot; height=&quot;326.34&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1320.20&quot; y=&quot;2421.15&quot; width=&quot;1509.56&quot; height=&quot;320.58&quot; style=&quot;fill:none;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1360.51&quot; y1=&quot;2511.46&quot; x2=&quot;1734.85&quot; y2=&quot;2511.46&quot; style=&quot;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;1360.51&quot; y1=&quot;2651.43&quot; x2=&quot;1734.85&quot; y2=&quot;2651.43&quot; style=&quot;stroke:#90353B;stroke-width:8.64&quot;/&gt;
	&lt;text x=&quot;1794.87&quot; y=&quot;2546.49&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; c&lt;/text&gt;
	&lt;text x=&quot;1794.87&quot; y=&quot;2686.46&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; c&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



As you can see, the response of consumption in the case of a persistent shocks is a lot larger and longer-lasting. Let's see all the comparisons of impulse-response functions at once:


```stata
foreach v in c h x k y r w z {
    irf ograph (persistent z `v' irf) (transitory z `v' irf), name(`v') nodraw
}

gr combine c h x k y r w z, rows(2) cols(4)

```

    
    graph c already exists
    graph h already exists
    graph x already exists
    graph k already exists
    graph y already exists
    graph r already exists
    graph w already exists
    graph z already exists
    


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;63.36&quot; y=&quot;63.36&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;65.28&quot; y=&quot;65.28&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;217.31&quot; y=&quot;130.56&quot; width=&quot;737.18&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;219.22&quot; y=&quot;132.48&quot; width=&quot;733.34&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;807.00&quot; x2=&quot;954.48&quot; y2=&quot;807.00&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;595.63&quot; x2=&quot;954.48&quot; y2=&quot;595.63&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;384.13&quot; x2=&quot;954.48&quot; y2=&quot;384.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;172.76&quot; x2=&quot;954.48&quot; y2=&quot;172.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M259.50 645.38 L341.06 461.73 L422.61 340.32 L504.28 266.69 L585.83 229.69 L667.38 220.28 L748.93 231.42 L830.61 257.78 L912.16 294.78&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M259.50 860.34 L341.06 726.81 L422.61 673.96 L504.28 667.40 L585.83 686.71 L667.38 719.88 L748.93 759.85 L830.61 802.17 L912.16 844.74&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;902.41&quot; x2=&quot;217.31&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;807.00&quot; x2=&quot;190.70&quot; y2=&quot;807.00&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.26&quot; y=&quot;807.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.26,807.00)&quot; text-anchor=&quot;middle&quot;&gt;.15&lt;/text&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;595.63&quot; x2=&quot;190.70&quot; y2=&quot;595.63&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.26&quot; y=&quot;595.63&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.26,595.63)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;384.13&quot; x2=&quot;190.70&quot; y2=&quot;384.13&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.26&quot; y=&quot;384.13&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.26,384.13)&quot; text-anchor=&quot;middle&quot;&gt;.25&lt;/text&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;172.76&quot; x2=&quot;190.70&quot; y2=&quot;172.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.26&quot; y=&quot;172.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.26,172.76)&quot; text-anchor=&quot;middle&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;217.31&quot; y1=&quot;902.41&quot; x2=&quot;954.48&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;259.50&quot; y1=&quot;902.41&quot; x2=&quot;259.50&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;259.50&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;422.73&quot; y1=&quot;902.41&quot; x2=&quot;422.73&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;422.73&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;585.83&quot; y1=&quot;902.41&quot; x2=&quot;585.83&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;585.83&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;749.06&quot; y1=&quot;902.41&quot; x2=&quot;749.06&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;749.06&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;912.16&quot; y1=&quot;902.41&quot; x2=&quot;912.16&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;912.16&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;585.83&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;240.32&quot; y=&quot;1118.86&quot; width=&quot;691.14&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;242.24&quot; y=&quot;1120.78&quot; width=&quot;687.30&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;269.16&quot; y1=&quot;1180.98&quot; x2=&quot;431.64&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;269.16&quot; y1=&quot;1287.66&quot; x2=&quot;431.64&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;470.75&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; c&lt;/text&gt;
	&lt;text x=&quot;470.75&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; c&lt;/text&gt;
	&lt;rect x=&quot;1021.68&quot; y=&quot;63.36&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1023.60&quot; y=&quot;65.28&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;1175.50&quot; y=&quot;130.56&quot; width=&quot;737.18&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1177.42&quot; y=&quot;132.48&quot; width=&quot;733.34&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;860.21&quot; x2=&quot;1912.68&quot; y2=&quot;860.21&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;731.14&quot; x2=&quot;1912.68&quot; y2=&quot;731.14&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;602.19&quot; x2=&quot;1912.68&quot; y2=&quot;602.19&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;473.11&quot; x2=&quot;1912.68&quot; y2=&quot;473.11&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;344.16&quot; x2=&quot;1912.68&quot; y2=&quot;344.16&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;215.08&quot; x2=&quot;1912.68&quot; y2=&quot;215.08&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1217.58 223.38 L1299.25 373.24 L1380.80 489.69 L1462.48 579.42 L1544.03 648.10 L1625.58 700.32 L1707.26 739.31 L1788.81 768.14 L1870.48 788.93&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1217.58 172.88 L1299.25 447.87 L1380.80 608.99 L1462.48 702.18 L1544.03 754.77 L1625.58 783.24 L1707.26 797.47 L1788.81 803.53 L1870.48 804.65&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;902.41&quot; x2=&quot;1175.50&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;860.21&quot; x2=&quot;1148.77&quot; y2=&quot;860.21&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1115.33&quot; y=&quot;860.21&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1115.33,860.21)&quot; text-anchor=&quot;middle&quot;&gt;-.1&lt;/text&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;731.14&quot; x2=&quot;1148.77&quot; y2=&quot;731.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1115.33&quot; y=&quot;731.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1115.33,731.14)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;602.19&quot; x2=&quot;1148.77&quot; y2=&quot;602.19&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1115.33&quot; y=&quot;602.19&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1115.33,602.19)&quot; text-anchor=&quot;middle&quot;&gt;.1&lt;/text&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;473.11&quot; x2=&quot;1148.77&quot; y2=&quot;473.11&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1115.33&quot; y=&quot;473.11&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1115.33,473.11)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;344.16&quot; x2=&quot;1148.77&quot; y2=&quot;344.16&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1115.33&quot; y=&quot;344.16&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1115.33,344.16)&quot; text-anchor=&quot;middle&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;215.08&quot; x2=&quot;1148.77&quot; y2=&quot;215.08&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1115.33&quot; y=&quot;215.08&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1115.33,215.08)&quot; text-anchor=&quot;middle&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;1175.50&quot; y1=&quot;902.41&quot; x2=&quot;1912.68&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1217.70&quot; y1=&quot;902.41&quot; x2=&quot;1217.70&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1217.70&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1380.93&quot; y1=&quot;902.41&quot; x2=&quot;1380.93&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1380.93&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;1544.15&quot; y1=&quot;902.41&quot; x2=&quot;1544.15&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1544.15&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;1707.26&quot; y1=&quot;902.41&quot; x2=&quot;1707.26&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1707.26&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;1870.48&quot; y1=&quot;902.41&quot; x2=&quot;1870.48&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1870.48&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;1544.15&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;1198.52&quot; y=&quot;1118.86&quot; width=&quot;691.14&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1200.44&quot; y=&quot;1120.78&quot; width=&quot;687.30&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1227.35&quot; y1=&quot;1180.98&quot; x2=&quot;1389.59&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1227.35&quot; y1=&quot;1287.66&quot; x2=&quot;1389.59&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1428.57&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; h&lt;/text&gt;
	&lt;text x=&quot;1428.57&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; h&lt;/text&gt;
	&lt;rect x=&quot;1980.00&quot; y=&quot;63.36&quot; width=&quot;958.20&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1981.92&quot; y=&quot;65.28&quot; width=&quot;954.36&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;2134.19&quot; y=&quot;130.56&quot; width=&quot;736.81&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2136.11&quot; y=&quot;132.48&quot; width=&quot;732.97&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;822.96&quot; x2=&quot;2871.00&quot; y2=&quot;822.96&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;692.90&quot; x2=&quot;2871.00&quot; y2=&quot;692.90&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;562.83&quot; x2=&quot;2871.00&quot; y2=&quot;562.83&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;432.89&quot; x2=&quot;2871.00&quot; y2=&quot;432.89&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;302.83&quot; x2=&quot;2871.00&quot; y2=&quot;302.83&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;172.76&quot; x2=&quot;2871.00&quot; y2=&quot;172.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2176.39 262.23 L2257.94 401.95 L2339.49 511.85 L2421.05 597.73 L2502.60 664.68 L2584.02 716.66 L2665.57 756.75 L2747.13 787.45 L2828.68 810.59&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2176.39 224.49 L2257.94 493.41 L2339.49 652.55 L2421.05 745.99 L2502.60 800.19 L2584.02 831.01 L2665.57 847.84 L2747.13 856.50 L2828.68 860.34&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;902.41&quot; x2=&quot;2134.19&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;822.96&quot; x2=&quot;2107.46&quot; y2=&quot;822.96&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2074.15&quot; y=&quot;822.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2074.15,822.96)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;692.90&quot; x2=&quot;2107.46&quot; y2=&quot;692.90&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2074.15&quot; y=&quot;692.90&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2074.15,692.90)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;562.83&quot; x2=&quot;2107.46&quot; y2=&quot;562.83&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2074.15&quot; y=&quot;562.83&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2074.15,562.83)&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;432.89&quot; x2=&quot;2107.46&quot; y2=&quot;432.89&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2074.15&quot; y=&quot;432.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2074.15,432.89)&quot; text-anchor=&quot;middle&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;302.83&quot; x2=&quot;2107.46&quot; y2=&quot;302.83&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2074.15&quot; y=&quot;302.83&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2074.15,302.83)&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;172.76&quot; x2=&quot;2107.46&quot; y2=&quot;172.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2074.15&quot; y=&quot;172.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2074.15,172.76)&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;2134.19&quot; y1=&quot;902.41&quot; x2=&quot;2871.00&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2176.39&quot; y1=&quot;902.41&quot; x2=&quot;2176.39&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2176.39&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2339.49&quot; y1=&quot;902.41&quot; x2=&quot;2339.49&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2339.49&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2502.60&quot; y1=&quot;902.41&quot; x2=&quot;2502.60&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2502.60&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;2665.70&quot; y1=&quot;902.41&quot; x2=&quot;2665.70&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2665.70&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;2828.80&quot; y1=&quot;902.41&quot; x2=&quot;2828.80&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2828.80&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;2502.60&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;2157.21&quot; y=&quot;1118.86&quot; width=&quot;690.77&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2159.13&quot; y=&quot;1120.78&quot; width=&quot;686.93&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2186.04&quot; y1=&quot;1180.98&quot; x2=&quot;2349.15&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2186.04&quot; y1=&quot;1287.66&quot; x2=&quot;2349.15&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2388.38&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; x&lt;/text&gt;
	&lt;text x=&quot;2388.38&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; x&lt;/text&gt;
	&lt;rect x=&quot;2938.20&quot; y=&quot;63.36&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2940.12&quot; y=&quot;65.28&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;3092.02&quot; y=&quot;130.56&quot; width=&quot;737.30&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3093.94&quot; y=&quot;132.48&quot; width=&quot;733.46&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;860.21&quot; x2=&quot;3829.32&quot; y2=&quot;860.21&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;645.25&quot; x2=&quot;3829.32&quot; y2=&quot;645.25&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;430.29&quot; x2=&quot;3829.32&quot; y2=&quot;430.29&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;215.21&quot; x2=&quot;3829.32&quot; y2=&quot;215.21&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3134.22 860.34 L3215.77 628.42 L3297.44 460.24 L3378.99 341.56 L3460.67 261.49 L3542.22 211.00 L3623.77 183.28 L3705.45 172.88 L3787.00 175.36&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3134.22 860.34 L3215.77 612.95 L3297.44 482.89 L3378.99 421.88 L3460.67 401.09 L3542.22 403.07 L3623.77 417.79 L3705.45 439.20 L3787.00 463.58&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;902.41&quot; x2=&quot;3092.02&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;860.21&quot; x2=&quot;3065.41&quot; y2=&quot;860.21&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3031.97&quot; y=&quot;860.21&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3031.97,860.21)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;645.25&quot; x2=&quot;3065.41&quot; y2=&quot;645.25&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3031.97&quot; y=&quot;645.25&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3031.97,645.25)&quot; text-anchor=&quot;middle&quot;&gt;.1&lt;/text&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;430.29&quot; x2=&quot;3065.41&quot; y2=&quot;430.29&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3031.97&quot; y=&quot;430.29&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3031.97,430.29)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;215.21&quot; x2=&quot;3065.41&quot; y2=&quot;215.21&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3031.97&quot; y=&quot;215.21&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3031.97,215.21)&quot; text-anchor=&quot;middle&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;3092.02&quot; y1=&quot;902.41&quot; x2=&quot;3829.32&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3134.22&quot; y1=&quot;902.41&quot; x2=&quot;3134.22&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3134.22&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3297.44&quot; y1=&quot;902.41&quot; x2=&quot;3297.44&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3297.44&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;3460.67&quot; y1=&quot;902.41&quot; x2=&quot;3460.67&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3460.67&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;3623.89&quot; y1=&quot;902.41&quot; x2=&quot;3623.89&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3623.89&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;3787.12&quot; y1=&quot;902.41&quot; x2=&quot;3787.12&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3787.12&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;3460.67&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;3115.03&quot; y=&quot;1118.86&quot; width=&quot;691.27&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3116.95&quot; y=&quot;1120.78&quot; width=&quot;687.43&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3143.87&quot; y1=&quot;1180.98&quot; x2=&quot;3306.11&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3143.87&quot; y1=&quot;1287.66&quot; x2=&quot;3306.11&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3345.21&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; k&lt;/text&gt;
	&lt;text x=&quot;3345.21&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; k&lt;/text&gt;
	&lt;rect x=&quot;63.36&quot; y=&quot;1440.00&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;65.28&quot; y=&quot;1441.92&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;217.43&quot; y=&quot;1507.20&quot; width=&quot;737.05&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;219.35&quot; y=&quot;1509.12&quot; width=&quot;733.22&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;2236.85&quot; x2=&quot;954.48&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;2099.73&quot; x2=&quot;954.48&quot; y2=&quot;2099.73&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1962.74&quot; x2=&quot;954.48&quot; y2=&quot;1962.74&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1825.62&quot; x2=&quot;954.48&quot; y2=&quot;1825.62&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1688.50&quot; x2=&quot;954.48&quot; y2=&quot;1688.50&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1551.50&quot; x2=&quot;954.48&quot; y2=&quot;1551.50&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M259.63 1568.21 L341.18 1697.78 L422.73 1801.73 L504.28 1885.27 L585.83 1952.22 L667.51 2006.05 L749.06 2049.36 L830.61 2084.26 L912.16 2112.35&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M259.63 1549.52 L341.18 1820.05 L422.73 1982.66 L504.28 2080.55 L585.83 2139.58 L667.51 2175.22 L749.06 2196.88 L830.61 2210.12 L912.16 2218.16&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;2279.05&quot; x2=&quot;217.43&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;2236.85&quot; x2=&quot;190.82&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.39&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.39,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;2099.73&quot; x2=&quot;190.82&quot; y2=&quot;2099.73&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.39&quot; y=&quot;2099.73&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.39,2099.73)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1962.74&quot; x2=&quot;190.82&quot; y2=&quot;1962.74&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.39&quot; y=&quot;1962.74&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.39,1962.74)&quot; text-anchor=&quot;middle&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1825.62&quot; x2=&quot;190.82&quot; y2=&quot;1825.62&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.39&quot; y=&quot;1825.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.39,1825.62)&quot; text-anchor=&quot;middle&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1688.50&quot; x2=&quot;190.82&quot; y2=&quot;1688.50&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.39&quot; y=&quot;1688.50&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.39,1688.50)&quot; text-anchor=&quot;middle&quot;&gt;.8&lt;/text&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;1551.50&quot; x2=&quot;190.82&quot; y2=&quot;1551.50&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;157.39&quot; y=&quot;1551.50&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 157.39,1551.50)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;217.43&quot; y1=&quot;2279.05&quot; x2=&quot;954.48&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;259.63&quot; y1=&quot;2279.05&quot; x2=&quot;259.63&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;259.63&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;422.85&quot; y1=&quot;2279.05&quot; x2=&quot;422.85&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;422.85&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;585.96&quot; y1=&quot;2279.05&quot; x2=&quot;585.96&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;585.96&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;749.06&quot; y1=&quot;2279.05&quot; x2=&quot;749.06&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;749.06&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;912.16&quot; y1=&quot;2279.05&quot; x2=&quot;912.16&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;912.16&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;585.96&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;240.45&quot; y=&quot;2495.50&quot; width=&quot;691.02&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;242.37&quot; y=&quot;2497.42&quot; width=&quot;687.18&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;269.28&quot; y1=&quot;2557.62&quot; x2=&quot;432.13&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;269.28&quot; y1=&quot;2664.30&quot; x2=&quot;432.13&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;471.24&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; y&lt;/text&gt;
	&lt;text x=&quot;471.24&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; y&lt;/text&gt;
	&lt;rect x=&quot;1021.68&quot; y=&quot;1440.00&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1023.60&quot; y=&quot;1441.92&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;1176.49&quot; y=&quot;1507.20&quot; width=&quot;736.19&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1178.41&quot; y=&quot;1509.12&quot; width=&quot;732.35&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;2008.15&quot; x2=&quot;1912.68&quot; y2=&quot;2008.15&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;1779.46&quot; x2=&quot;1912.68&quot; y2=&quot;1779.46&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;1550.76&quot; x2=&quot;1912.68&quot; y2=&quot;1550.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1218.57 1562.02 L1300.12 1697.78 L1381.55 1802.97 L1463.10 1883.91 L1544.52 1945.66 L1625.95 1992.31 L1707.50 2027.09 L1788.93 2052.58 L1870.48 2070.90&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1218.57 1549.52 L1300.12 1782.67 L1381.55 1918.80 L1463.10 1997.14 L1544.52 2040.95 L1625.95 2064.34 L1707.50 2075.60 L1788.93 2079.93 L1870.48 2080.18&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;2279.05&quot; x2=&quot;1176.49&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;2236.85&quot; x2=&quot;1149.76&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1116.32&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1116.32,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;-.5&lt;/text&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;2008.15&quot; x2=&quot;1149.76&quot; y2=&quot;2008.15&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1116.32&quot; y=&quot;2008.15&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1116.32,2008.15)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;1779.46&quot; x2=&quot;1149.76&quot; y2=&quot;1779.46&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1116.32&quot; y=&quot;1779.46&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1116.32,1779.46)&quot; text-anchor=&quot;middle&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;1550.76&quot; x2=&quot;1149.76&quot; y2=&quot;1550.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1116.32&quot; y=&quot;1550.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1116.32,1550.76)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;1176.49&quot; y1=&quot;2279.05&quot; x2=&quot;1912.68&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1218.69&quot; y1=&quot;2279.05&quot; x2=&quot;1218.69&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1218.69&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1381.67&quot; y1=&quot;2279.05&quot; x2=&quot;1381.67&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1381.67&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;1544.65&quot; y1=&quot;2279.05&quot; x2=&quot;1544.65&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1544.65&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;1707.50&quot; y1=&quot;2279.05&quot; x2=&quot;1707.50&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1707.50&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;1870.48&quot; y1=&quot;2279.05&quot; x2=&quot;1870.48&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1870.48&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;1544.65&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;1199.51&quot; y=&quot;2495.50&quot; width=&quot;690.15&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1201.43&quot; y=&quot;2497.42&quot; width=&quot;686.31&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1228.34&quot; y1=&quot;2557.62&quot; x2=&quot;1392.56&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1228.34&quot; y1=&quot;2664.30&quot; x2=&quot;1392.56&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1432.04&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; r&lt;/text&gt;
	&lt;text x=&quot;1432.04&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; r&lt;/text&gt;
	&lt;rect x=&quot;1980.00&quot; y=&quot;1440.00&quot; width=&quot;958.20&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1981.92&quot; y=&quot;1441.92&quot; width=&quot;954.36&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;2132.95&quot; y=&quot;1507.20&quot; width=&quot;738.05&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2134.87&quot; y=&quot;1509.12&quot; width=&quot;734.21&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;2215.81&quot; x2=&quot;2871.00&quot; y2=&quot;2215.81&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;2082.53&quot; x2=&quot;2871.00&quot; y2=&quot;2082.53&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1949.25&quot; x2=&quot;2871.00&quot; y2=&quot;1949.25&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1815.96&quot; x2=&quot;2871.00&quot; y2=&quot;1815.96&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1682.68&quot; x2=&quot;2871.00&quot; y2=&quot;1682.68&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1549.40&quot; x2=&quot;2871.00&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2175.03 1573.65 L2256.83 1670.68 L2338.50 1752.48 L2420.18 1822.03 L2501.85 1881.43 L2583.65 1932.29 L2665.33 1976.10 L2747.00 2014.22 L2828.68 2047.38&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2175.03 1589.25 L2256.83 1831.19 L2338.50 1980.93 L2420.18 2075.10 L2501.85 2135.50 L2583.65 2175.47 L2665.33 2202.82 L2747.00 2222.25 L2828.68 2236.97&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;2279.05&quot; x2=&quot;2132.95&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;2215.81&quot; x2=&quot;2106.22&quot; y2=&quot;2215.81&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2072.79&quot; y=&quot;2215.81&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2072.79,2215.81)&quot; text-anchor=&quot;middle&quot;&gt;.1&lt;/text&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;2082.53&quot; x2=&quot;2106.22&quot; y2=&quot;2082.53&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2072.79&quot; y=&quot;2082.53&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2072.79,2082.53)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1949.25&quot; x2=&quot;2106.22&quot; y2=&quot;1949.25&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2072.79&quot; y=&quot;1949.25&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2072.79,1949.25)&quot; text-anchor=&quot;middle&quot;&gt;.3&lt;/text&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1815.96&quot; x2=&quot;2106.22&quot; y2=&quot;1815.96&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2072.79&quot; y=&quot;1815.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2072.79,1815.96)&quot; text-anchor=&quot;middle&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1682.68&quot; x2=&quot;2106.22&quot; y2=&quot;1682.68&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2072.79&quot; y=&quot;1682.68&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2072.79,1682.68)&quot; text-anchor=&quot;middle&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;1549.40&quot; x2=&quot;2106.22&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2072.79&quot; y=&quot;1549.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2072.79,1549.40)&quot; text-anchor=&quot;middle&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;2132.95&quot; y1=&quot;2279.05&quot; x2=&quot;2871.00&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2175.15&quot; y1=&quot;2279.05&quot; x2=&quot;2175.15&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2175.15&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2338.50&quot; y1=&quot;2279.05&quot; x2=&quot;2338.50&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2338.50&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2501.98&quot; y1=&quot;2279.05&quot; x2=&quot;2501.98&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2501.98&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;2665.33&quot; y1=&quot;2279.05&quot; x2=&quot;2665.33&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2665.33&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;2828.80&quot; y1=&quot;2279.05&quot; x2=&quot;2828.80&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2828.80&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;2501.98&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;2155.97&quot; y=&quot;2495.50&quot; width=&quot;692.01&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2157.89&quot; y=&quot;2497.42&quot; width=&quot;688.17&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2184.81&quot; y1=&quot;2557.62&quot; x2=&quot;2345.31&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2184.81&quot; y1=&quot;2664.30&quot; x2=&quot;2345.31&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2383.92&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; w&lt;/text&gt;
	&lt;text x=&quot;2383.92&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; w&lt;/text&gt;
	&lt;rect x=&quot;2938.20&quot; y=&quot;1440.00&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2940.12&quot; y=&quot;1441.92&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;3092.26&quot; y=&quot;1507.20&quot; width=&quot;737.06&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3094.18&quot; y=&quot;1509.12&quot; width=&quot;733.22&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;2236.85&quot; x2=&quot;3829.32&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;2099.36&quot; x2=&quot;3829.32&quot; y2=&quot;2099.36&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1961.87&quot; x2=&quot;3829.32&quot; y2=&quot;1961.87&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1824.38&quot; x2=&quot;3829.32&quot; y2=&quot;1824.38&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1686.89&quot; x2=&quot;3829.32&quot; y2=&quot;1686.89&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1549.40&quot; x2=&quot;3829.32&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3134.46 1549.52 L3216.01 1687.01 L3297.57 1797.03 L3379.24 1885.02 L3460.79 1955.31 L3542.34 2011.62 L3623.89 2056.67 L3705.45 2092.80 L3787.00 2121.64&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3134.46 1549.52 L3216.01 1824.50 L3297.57 1989.47 L3379.24 2088.47 L3460.79 2147.87 L3542.34 2183.51 L3623.89 2204.80 L3705.45 2217.67 L3787.00 2225.34&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;2279.05&quot; x2=&quot;3092.26&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;2236.85&quot; x2=&quot;3065.66&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3032.22&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3032.22,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;2099.36&quot; x2=&quot;3065.66&quot; y2=&quot;2099.36&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3032.22&quot; y=&quot;2099.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3032.22,2099.36)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1961.87&quot; x2=&quot;3065.66&quot; y2=&quot;1961.87&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3032.22&quot; y=&quot;1961.87&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3032.22,1961.87)&quot; text-anchor=&quot;middle&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1824.38&quot; x2=&quot;3065.66&quot; y2=&quot;1824.38&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3032.22&quot; y=&quot;1824.38&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3032.22,1824.38)&quot; text-anchor=&quot;middle&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1686.89&quot; x2=&quot;3065.66&quot; y2=&quot;1686.89&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3032.22&quot; y=&quot;1686.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3032.22,1686.89)&quot; text-anchor=&quot;middle&quot;&gt;.8&lt;/text&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;1549.40&quot; x2=&quot;3065.66&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3032.22&quot; y=&quot;1549.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3032.22,1549.40)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;3092.26&quot; y1=&quot;2279.05&quot; x2=&quot;3829.32&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3134.59&quot; y1=&quot;2279.05&quot; x2=&quot;3134.59&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3134.59&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3297.69&quot; y1=&quot;2279.05&quot; x2=&quot;3297.69&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3297.69&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;3460.79&quot; y1=&quot;2279.05&quot; x2=&quot;3460.79&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3460.79&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;3623.89&quot; y1=&quot;2279.05&quot; x2=&quot;3623.89&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3623.89&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;6&lt;/text&gt;
	&lt;line x1=&quot;3787.12&quot; y1=&quot;2279.05&quot; x2=&quot;3787.12&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3787.12&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;8&lt;/text&gt;
	&lt;text x=&quot;3460.79&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;3115.41&quot; y=&quot;2495.50&quot; width=&quot;690.90&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3117.33&quot; y=&quot;2497.42&quot; width=&quot;687.06&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3144.12&quot; y1=&quot;2557.62&quot; x2=&quot;3306.97&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3144.12&quot; y1=&quot;2664.30&quot; x2=&quot;3306.97&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3346.20&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;persistent: irf of z -&amp;gt; z&lt;/text&gt;
	&lt;text x=&quot;3346.20&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;transitory: irf of z -&amp;gt; z&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



    
    
    
    

## Estimating the model

Instead of imposing all the parameter values, we can assume some and estimate others. For example, let's assume, as before $\beta$=0.96, $\sigma$=1, $\alpha$=0.3, $\delta$=0.025, $\phi$=1, $\theta_i$=0.3, and $\theta_c$=0.7. We will estimate $\rho_z$ and $\sigma_z$.

We set up the constraints:


```stata
constraint 1 _b[beta] = 0.96
constraint 2 _b[sigma] = 1
constraint 3 _b[alpha] = 0.3
constraint 4 _b[delta] = 0.025
constraint 5 _b[phi] = 1
constraint 6 _b[thetai] = 0.3
constraint 7 _b[thetac] = 0.7
```

And we rewrite the dsge command, but now to *estimate* the model:


```stata
dsge ({thetai}*x = y - {thetac}*c, unobserved) ///
     (F.k = {delta}*x+ (1-{delta})*k, state noshock) ///
     (c = F.c - ((1-{beta}+{beta}*{delta})/{sigma})*F.r, unobserved) ///
     ({phi}*h = w - {sigma}*c, unobserved) ///
     (y = (1-{alpha})*(z+h) + {alpha}*k) ///
     (w = y - h, unobserved) ///
     (r = y - k, unobserved) ///   
     (F.z = {rhoz}*z, state), ///
     constraint(1/7) 
```

    
    (setting technique to bfgs)
    Iteration 0:   log likelihood = -1489.9913  
    Iteration 1:   log likelihood = -806.32338  (backed up)
    Iteration 2:   log likelihood = -802.14483  
    Iteration 3:   log likelihood = -770.30362  
    Iteration 4:   log likelihood =  -769.3696  
    (switching technique to nr)
    Iteration 5:   log likelihood =  -769.3696  (not concave)
    Iteration 6:   log likelihood = -641.83043  
    Iteration 7:   log likelihood = -639.20494  
    Iteration 8:   log likelihood = -639.17484  
    Iteration 9:   log likelihood = -639.17481  
    
    DSGE model
    
    Sample: 1955q1 - 2015q4                         Number of obs     =        244
    Log likelihood = -639.17481
     ( 1)  [/structural]beta = .96
     ( 2)  [/structural]sigma = 1
     ( 3)  [/structural]alpha = .3
     ( 4)  [/structural]delta = .025
     ( 5)  [/structural]phi = 1
     ( 6)  [/structural]thetai = .3
     ( 7)  [/structural]thetac = .7
    ------------------------------------------------------------------------------
                 |                 OIM
               y |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
    /structural  |
          thetai |         .3  (constrained)
          thetac |         .7  (constrained)
           delta |       .025  (constrained)
            beta |        .96  (constrained)
           sigma |          1  (constrained)
             phi |          1  (constrained)
           alpha |         .3  (constrained)
            rhoz |   .3370701   .0609055     5.53   0.000     .2176975    .4564427
    -------------+----------------------------------------------------------------
          sd(e.z)|   3.204836    .145079                      2.920486    3.489185
    ------------------------------------------------------------------------------
    

The estimated persistence of technology, 0.34, is a *a lot* lower than any of the two values we assumed before: 0.6 and 0.8. 

### Prediction

We can also produce the one-step-ahead prediction of output:


```stata
predict dep*
```

    (option xb assumed; fitted values)
    


```stata
tsline y dep1, legend(col(1))
```


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;290.81&quot; y=&quot;100.86&quot; width=&quot;3568.21&quot; height=&quot;1992.81&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;293.69&quot; y=&quot;103.74&quot; width=&quot;3562.45&quot; height=&quot;1987.05&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1993.67&quot; x2=&quot;3859.02&quot; y2=&quot;1993.67&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1631.45&quot; x2=&quot;3859.02&quot; y2=&quot;1631.45&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1269.34&quot; x2=&quot;3859.02&quot; y2=&quot;1269.34&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;907.12&quot; x2=&quot;3859.02&quot; y2=&quot;907.12&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;545.01&quot; x2=&quot;3859.02&quot; y2=&quot;545.01&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;182.78&quot; x2=&quot;3859.02&quot; y2=&quot;182.78&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;path d=&quot; M354.05 453.81 L368.28 801.55 L382.39 882.37 L396.62 1095.59 L410.73 1381.34 L424.83 1031.49 L439.06 1293.10 L453.17 799.20 L467.40 1083.22 L481.51 1333.70 L495.74 988.55 L509.85 1567.96 L524.08 2030.43 L538.19 1078.89 L552.42 607.14 L566.53 599.46 L580.63 735.35 L594.87 573.23 L608.97 1325.28 L623.21 1155.74 L637.31 631.02 L651.54 1379.48 L665.65 1197.32 L679.88 1622.29 L693.99 1073.07 L708.10 735.72 L722.33 790.29 L736.44 688.07 L750.67 752.18 L764.77 958.84 L779.01 994.73 L793.11 1157.10 L807.35 948.70 L821.45 894.12 L835.56 709.23 L849.79 1063.29 L863.90 649.34 L878.13 926.92 L892.24 878.53 L906.47 1166.88 L920.58 563.45 L934.81 875.93 L948.91 686.83 L963.02 594.02 L977.25 562.59 L991.36 1150.91 L1005.59 1063.66 L1019.70 1022.08 L1033.93 1004.63 L1048.04 1243.97 L1062.27 1020.35 L1076.38 1037.06 L1090.61 686.96 L1104.72 784.23 L1118.82 1061.44 L1133.06 1140.64 L1147.16 822.34 L1161.39 1176.65 L1175.50 1088.54 L1189.73 1396.07 L1203.84 1320.83 L1218.07 1218.11 L1232.18 1013.79 L1246.29 1568.95 L1260.52 503.80 L1274.63 1104.50 L1288.86 1043.74 L1302.96 1184.82 L1317.19 754.16 L1331.30 606.02 L1345.53 1003.52 L1359.64 792.15 L1373.75 564.69 L1387.98 942.76 L1402.09 1427.25 L1416.32 1000.43 L1430.43 1511.03 L1444.66 1193.36 L1458.77 1551.25 L1473.00 1385.55 L1487.10 1621.55 L1501.34 1046.96 L1515.44 794.50 L1529.55 881.38 L1543.78 622.36 L1557.89 1051.29 L1572.12 1122.45 L1586.23 1052.77 L1600.46 934.22 L1614.57 706.26 L1628.80 761.33 L1642.90 1266.37 L1657.01 1168.73 L1671.24 164.22 L1685.35 987.43 L1699.58 883.36 L1713.69 1211.92 L1727.92 1234.32 L1742.03 1062.06 L1756.26 1194.47 L1770.37 1175.79 L1784.47 1862.99 L1798.71 1313.28 L1812.81 737.33 L1827.05 675.94 L1841.15 1481.71 L1855.38 938.80 L1869.49 1609.42 L1883.72 1758.17 L1897.83 1111.93 L1911.94 1373.92 L1926.17 1241.13 L1940.28 892.39 L1954.51 615.68 L1968.62 707.62 L1982.85 677.92 L1996.95 699.33 L2011.18 764.92 L2025.29 985.70 L2039.52 1039.16 L2053.63 982.73 L2067.74 1005.38 L2081.97 822.10 L2096.08 1052.77 L2110.31 1002.04 L2124.42 1136.68 L2138.65 979.14 L2152.76 1119.72 L2166.99 1067.50 L2181.09 945.97 L2195.20 1007.73 L2209.43 795.12 L2223.54 1106.85 L2237.77 889.05 L2251.88 1102.65 L2266.11 888.06 L2280.22 978.52 L2294.45 1042.38 L2308.56 1053.89 L2322.66 1208.09 L2336.89 953.89 L2351.00 1157.59 L2365.23 1262.17 L2379.34 1517.47 L2393.57 1405.60 L2407.68 1045.35 L2421.91 1130.62 L2436.02 1143.11 L2450.13 929.14 L2464.36 951.67 L2478.47 989.04 L2492.70 980.50 L2506.80 1215.39 L2521.03 1097.57 L2535.14 1128.39 L2549.37 885.09 L2563.48 986.69 L2577.71 876.05 L2591.82 1099.06 L2605.93 942.26 L2620.16 1170.46 L2634.27 1168.36 L2648.50 1022.45 L2662.61 1064.41 L2676.84 1079.88 L2690.94 767.89 L2705.18 1002.53 L2719.28 964.91 L2733.39 1049.43 L2747.62 835.59 L2761.73 902.78 L2775.96 1045.47 L2790.07 984.09 L2804.30 989.78 L2818.41 892.76 L2832.64 797.59 L2846.74 1038.91 L2860.85 1031.36 L2875.08 906.99 L2889.19 770.74 L2903.42 1185.44 L2917.53 727.18 L2931.76 1234.57 L2945.87 1105.12 L2960.10 1351.89 L2974.21 1116.14 L2988.32 1361.17 L3002.55 1189.03 L3016.65 1003.77 L3030.89 1109.95 L3044.99 1128.64 L3059.22 1250.90 L3073.33 1119.72 L3087.56 1001.79 L3101.67 788.19 L3115.90 932.61 L3130.01 1103.27 L3144.12 1057.72 L3158.35 1007.23 L3172.45 1019.98 L3186.69 962.19 L3200.79 1118.49 L3215.03 1027.03 L3229.13 1104.26 L3243.36 923.33 L3257.47 1182.96 L3271.58 1243.48 L3285.81 1043.49 L3299.92 1251.52 L3314.15 1048.57 L3328.26 1075.17 L3342.49 1166.13 L3356.59 1467.84 L3370.83 1125.91 L3384.93 1408.57 L3399.04 1888.11 L3413.27 1673.65 L3427.38 1308.57 L3441.61 1174.80 L3455.72 990.28 L3469.95 1144.23 L3484.06 990.77 L3498.29 1074.43 L3512.40 1087.43 L3526.50 1381.46 L3540.74 1059.33 L3554.84 1208.33 L3569.07 944.98 L3583.18 1077.77 L3597.41 1134.45 L3611.52 1234.69 L3625.75 1262.91 L3639.86 1067.38 L3654.09 1213.90 L3668.20 1046.59 L3682.30 988.05 L3696.54 1355.48 L3710.64 987.93 L3724.88 918.62 L3738.98 1103.76 L3753.21 1122.57 L3767.32 1082.60 L3781.55 1126.90 L3795.66 1206.35&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;path d=&quot; M354.05 1050.42 L368.28 848.09 L382.39 965.28 L396.62 992.51 L410.73 1064.65 L424.83 1161.55 L439.06 1043.37 L453.17 1132.10 L467.40 965.03 L481.51 1061.07 L495.74 1145.96 L509.85 1029.38 L524.08 1225.53 L538.19 1382.70 L552.42 1061.31 L566.53 901.42 L580.63 898.33 L594.87 943.87 L608.97 888.68 L623.21 1142.99 L637.31 1085.82 L651.54 908.23 L665.65 1161.31 L679.88 1100.05 L693.99 1244.10 L708.10 1058.59 L722.33 944.37 L736.44 962.43 L750.67 927.66 L764.77 948.94 L779.01 1018.74 L793.11 1030.87 L807.35 1085.82 L821.45 1015.40 L835.56 996.84 L849.79 934.09 L863.90 1053.76 L878.13 913.67 L892.24 1007.36 L906.47 990.90 L920.58 1088.54 L934.81 884.35 L948.91 989.78 L963.02 925.68 L977.25 894.00 L991.36 883.11 L1005.59 1082.10 L1019.70 1052.77 L1033.93 1038.79 L1048.04 1033.10 L1062.27 1114.16 L1076.38 1038.79 L1090.61 1044.48 L1104.72 926.05 L1118.82 958.72 L1133.06 1052.53 L1147.16 1079.50 L1161.39 971.84 L1175.50 1091.76 L1189.73 1062.18 L1203.84 1166.38 L1218.07 1141.26 L1232.18 1106.85 L1246.29 1037.80 L1260.52 1225.78 L1274.63 865.53 L1288.86 1068.49 L1302.96 1047.95 L1317.19 1095.84 L1331.30 950.06 L1345.53 899.57 L1359.64 1033.84 L1373.75 962.31 L1387.98 885.09 L1402.09 1012.68 L1416.32 1176.78 L1430.43 1032.60 L1444.66 1205.61 L1458.77 1098.44 L1473.00 1219.84 L1487.10 1164.28 L1501.34 1244.47 L1515.44 1050.42 L1529.55 964.79 L1543.78 993.87 L1557.89 906.00 L1572.12 1050.79 L1586.23 1074.93 L1600.46 1051.29 L1614.57 1011.19 L1628.80 933.85 L1642.90 952.16 L1657.01 1122.94 L1671.24 1090.15 L1685.35 750.07 L1699.58 1028.02 L1713.69 992.88 L1727.92 1104.01 L1742.03 1111.80 L1756.26 1053.76 L1770.37 1098.69 L1784.47 1092.50 L1798.71 1325.40 L1812.81 1140.02 L1827.05 945.11 L1841.15 923.95 L1855.38 1196.58 L1869.49 1013.05 L1883.72 1240.01 L1897.83 1290.88 L1911.94 1072.70 L1926.17 1161.31 L1940.28 1116.63 L1954.51 998.57 L1968.62 904.52 L1982.85 935.21 L1996.95 924.81 L2011.18 931.62 L2025.29 953.52 L2039.52 1028.02 L2053.63 1046.09 L2067.74 1027.03 L2081.97 1034.58 L2096.08 972.58 L2110.31 1050.55 L2124.42 1033.34 L2138.65 1079.01 L2152.76 1025.80 L2166.99 1073.32 L2181.09 1055.74 L2195.20 1014.66 L2209.43 1035.57 L2223.54 963.42 L2237.77 1068.86 L2251.88 995.23 L2266.11 1067.38 L2280.22 994.86 L2294.45 1025.42 L2308.56 1046.96 L2322.66 1050.92 L2336.89 1103.27 L2351.00 1017.26 L2365.23 1086.19 L2379.34 1121.83 L2393.57 1208.46 L2407.68 1171.08 L2421.91 1049.43 L2436.02 1078.14 L2450.13 1082.48 L2464.36 1010.08 L2478.47 1017.50 L2492.70 1030.00 L2506.80 1027.03 L2521.03 1106.48 L2535.14 1066.76 L2549.37 1077.15 L2563.48 994.86 L2577.71 1029.01 L2591.82 991.52 L2605.93 1066.88 L2620.16 1013.79 L2634.27 1091.01 L2648.50 1090.40 L2662.61 1041.14 L2676.84 1055.25 L2690.94 1060.57 L2705.18 954.88 L2719.28 1034.09 L2733.39 1021.22 L2747.62 1049.80 L2761.73 977.41 L2775.96 1000.06 L2790.07 1048.20 L2804.30 1027.53 L2818.41 1029.38 L2832.64 996.47 L2846.74 964.17 L2860.85 1045.72 L2875.08 1043.25 L2889.19 1001.05 L2903.42 954.88 L2917.53 1095.10 L2931.76 940.16 L2945.87 1111.68 L2960.10 1068.12 L2974.21 1151.78 L2988.32 1072.33 L3002.55 1155.37 L3016.65 1097.45 L3030.89 1034.83 L3044.99 1070.72 L3059.22 1077.03 L3073.33 1118.61 L3087.56 1074.31 L3101.67 1034.33 L3115.90 961.94 L3130.01 1010.57 L3144.12 1068.24 L3158.35 1052.90 L3172.45 1035.82 L3186.69 1040.15 L3200.79 1020.47 L3215.03 1073.32 L3229.13 1042.38 L3243.36 1068.61 L3257.47 1007.36 L3271.58 1095.22 L3285.81 1115.76 L3299.92 1048.32 L3314.15 1118.73 L3328.26 1050.18 L3342.49 1059.09 L3356.59 1089.90 L3370.83 1192.24 L3384.93 1076.78 L3399.04 1172.57 L3413.27 1335.18 L3427.38 1263.28 L3441.61 1140.14 L3455.72 1094.85 L3469.95 1032.35 L3484.06 1084.21 L3498.29 1032.23 L3512.40 1060.32 L3526.50 1064.65 L3540.74 1164.15 L3554.84 1055.25 L3569.07 1105.62 L3583.18 1016.39 L3597.41 1061.19 L3611.52 1080.25 L3625.75 1114.28 L3639.86 1123.81 L3654.09 1057.72 L3668.20 1107.23 L3682.30 1050.67 L3696.54 1030.75 L3710.64 1154.99 L3724.88 1030.62 L3738.98 1007.11 L3753.21 1069.48 L3767.32 1075.92 L3781.55 1062.30 L3795.66 1077.28&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2093.67&quot; x2=&quot;290.81&quot; y2=&quot;100.86&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1993.67&quot; x2=&quot;250.84&quot; y2=&quot;1993.67&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1993.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1993.67)&quot; text-anchor=&quot;middle&quot;&gt;-10&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1631.45&quot; x2=&quot;250.84&quot; y2=&quot;1631.45&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1631.45&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1631.45)&quot; text-anchor=&quot;middle&quot;&gt;-5&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1269.34&quot; x2=&quot;250.84&quot; y2=&quot;1269.34&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1269.34&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1269.34)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;907.12&quot; x2=&quot;250.84&quot; y2=&quot;907.12&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;907.12&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,907.12)&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;545.01&quot; x2=&quot;250.84&quot; y2=&quot;545.01&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;545.01&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,545.01)&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;182.78&quot; x2=&quot;250.84&quot; y2=&quot;182.78&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;182.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,182.78)&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2093.67&quot; x2=&quot;3859.02&quot; y2=&quot;2093.67&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;382.51&quot; y1=&quot;2093.67&quot; x2=&quot;382.51&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;382.51&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;1955q3&lt;/text&gt;
	&lt;line x1=&quot;1232.30&quot; y1=&quot;2093.67&quot; x2=&quot;1232.30&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1232.30&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;1970q3&lt;/text&gt;
	&lt;line x1=&quot;2081.97&quot; y1=&quot;2093.67&quot; x2=&quot;2081.97&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2081.97&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;1985q3&lt;/text&gt;
	&lt;line x1=&quot;2931.76&quot; y1=&quot;2093.67&quot; x2=&quot;2931.76&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2931.76&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2000q3&lt;/text&gt;
	&lt;line x1=&quot;3781.55&quot; y1=&quot;2093.67&quot; x2=&quot;3781.55&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3781.55&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;2015q3&lt;/text&gt;
	&lt;text x=&quot;2074.92&quot; y=&quot;2333.64&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;Date (quarters)&lt;/text&gt;
	&lt;rect x=&quot;996.43&quot; y=&quot;2418.27&quot; width=&quot;2156.96&quot; height=&quot;326.34&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;999.31&quot; y=&quot;2421.15&quot; width=&quot;2151.20&quot; height=&quot;320.58&quot; style=&quot;fill:none;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1039.62&quot; y1=&quot;2511.46&quot; x2=&quot;1414.09&quot; y2=&quot;2511.46&quot; style=&quot;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;1039.62&quot; y1=&quot;2651.43&quot; x2=&quot;1414.09&quot; y2=&quot;2651.43&quot; style=&quot;stroke:#90353B;stroke-width:8.64&quot;/&gt;
	&lt;text x=&quot;1473.99&quot; y=&quot;2546.49&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;Growth rate of real GDP (GDPC96)&lt;/text&gt;
	&lt;text x=&quot;1473.99&quot; y=&quot;2686.46&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;xb  prediction, y, onestep&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



Predicted output matches some of the swings in actual output, but its volatility is much lower. 

## Estimation of the non-linear model

A newer Stata command, dsgenl, allows us to estimate a DSGE directly, without linearizing it first. I'll estimate the same basic RBC as before:

$$C_t + I_t = Y_t              \text{    (1)}$$ 

$$K_{t+1} = (1-\delta)K_t + I_t \text{    (2)}$$ 

$$C^{-\sigma}_t = \beta E_t \left\{C^{-\sigma}_{t+1} \left[(1-\delta) + {R_{t+1}} \right] \right\} \text{    (3)}$$

$$C^{\sigma}_t H^{\phi}_t = {W_t} \text{    (4)}$$ 

$$Y_t = K^{\alpha}_t (Z_tH_t)^{1-\alpha} \text{    (5)}$$

$$ln(Z_{t+1}) = \rho_z ln(Z_t) + \epsilon^Z_{t+1} \text{    (6)}$$ 

$\epsilon^Z$ is a technology shock of the familiar type: $\epsilon^Z_t \sim N(0,\sigma_Z)$.

$${R_t} = \alpha \frac{Y_t}{K_t} \text{    (7)}$$ 

$${W_t} = (1-\alpha) \frac{Y_t}{H_t} \text{    (8)}$$

where I have already imposed the normalization $P_t$=1.


```stata
use https://www.stata-press.com/data/r16/usmacro2, replace
```

    (Federal Reserve Economic Data - St. Louis Fed, 2017-01-15)
    

We define the constraints, as before:


```stata
constraint 1 _b[beta] = 0.96
constraint 2 _b[sigma] = 1
constraint 3 _b[alpha] = 0.3
constraint 4 _b[delta] = 0.025
constraint 5 _b[phi] = 1
```

And we write the model in Stata syntax for dsgenl, which is different from the syntax for dsge:

dsgenl (y = c + x) 
     (F.k = x + (1-{delta})*k) 
     (1/(c^sigma) = {beta}*(1/((F.c)^sigma))*(1-{delta}+r)) 
     (h^phi = w/(c^sigma)) 
     (y = (k^{alpha})*((z*h)^(1-{alpha}))) 
     (ln(F.z)={rhoz}*ln(z)) 
     (r = {alpha}*y/k)    
     (w = (1-{alpha})*y/h), 
     , observed(y) unobserved(c x r w h) exostate(z) endostate(k) 
     constraint(1/5)

As you can see, now instead of indicating the type of variable at the end of each equation, we make a list of each type of variable after we write the equations. As before, the observed variable is output, the endogenous control variables are the demand components, the quantity of labor, and factor prices, the endogenous state variable is capital, which is deterministically determined by its own past value and investment, and the exogenous state is the technology variable.

Let's run the dsgenl command:


```stata
dsgenl (y = c + x) ///
     (F.k = x + (1-{delta})*k) ///
     (1/(c^{sigma}) = {beta}*(1/((F.c)^{sigma}))*(1-{delta}+r)) ///
     (h^{phi} = w/(c^{sigma})) ///
     (y = (k^{alpha})*((z*h)^(1-{alpha}))) ///
     (ln(F.z)={rhoz}*ln(z)) ///
     (r = {alpha}*y/k) ///   
     (w = (1-{alpha})*y/h) ///
     , observed(y) unobserved(c x r w h) exostate(z) endostate(k) ///
     constraint(1/5)
```

    
    Solving at initial parameter vector ...
    Checking identification ...
    (setting technique to bfgs)
    Iteration 0:   log likelihood = -1485.3519  
    Iteration 1:   log likelihood = -1422.3411  
    Iteration 2:   log likelihood = -1397.9517  (backed up)
    Iteration 3:   log likelihood = -1397.9517  (backed up)
    Iteration 4:   log likelihood = -1204.2349  
    Iteration 5:   log likelihood = -1204.2349  (backed up)
    BFGS stepping has contracted, resetting BFGS Hessian
    Iteration 6:   log likelihood = -1204.2349  (backed up)
    Iteration 7:   log likelihood = -1191.1306  (backed up)
    Iteration 8:   log likelihood = -1092.5378  
    Iteration 9:   log likelihood = -1092.5378  (backed up)
    (switching technique to nr)
    Iteration 10:  log likelihood = -1092.5378  (not concave)
    Iteration 11:  log likelihood = -772.52054  (not concave)
    Iteration 12:  log likelihood = -697.50474  
    Iteration 13:  log likelihood = -651.45979  
    Iteration 14:  log likelihood = -639.49151  
    Iteration 15:  log likelihood = -639.35322  
    Iteration 16:  log likelihood = -639.35279  
    Iteration 17:  log likelihood = -639.35279  
    
    First-order DSGE model
    
    Sample: 1955q1 - 2015q4                         Number of obs     =        244
    Log likelihood = -639.35279
     ( 1)  [/structural]beta = .96
     ( 2)  [/structural]sigma = 1
     ( 3)  [/structural]alpha = .3
     ( 4)  [/structural]delta = .025
     ( 5)  [/structural]phi = 1
    ------------------------------------------------------------------------------
                 |                 OIM
               y |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
    /structural  |
           delta |       .025  (constrained)
           sigma |          1  (constrained)
            beta |        .96  (constrained)
             phi |          1  (constrained)
           alpha |         .3  (constrained)
            rhoz |   .2981728   .0619475     4.81   0.000      .176758    .4195876
    -------------+----------------------------------------------------------------
          sd(e.z)|   3.206035   .1453994                      2.921058    3.491013
    ------------------------------------------------------------------------------
    

We can find the steady state values with the command


```stata
estat steady
```

    
    Location of model steady-state
    
    ------------------------------------------------------------------------------
                 |            Delta-method
                 |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
               k |   7.614223          .        .       .            .           .
               z |          1          .        .       .            .           .
               c |   1.501694          .        .       .            .           .
               x |   .1903556          .        .       .            .           .
               r |   .0666667          .        .       .            .           .
               w |   1.333663          .        .       .            .           .
               h |   .8881061          .        .       .            .           .
               y |   1.692049          .        .       .            .           .
    ------------------------------------------------------------------------------
    Note: Standard errors reported as missing for constrained steady-state values.
    

### Impulse-response functions


```stata
irf set rbcirf, replace
irf create est, step(20) replace
irf graph irf, impulse(z) response(y c x h w r z) byopts(yrescale)
```

    
    (file rbcirf.irf created)
    (file rbcirf.irf now active)
    
    irfname est not found in rbcirf.irf
    (file rbcirf.irf updated)
    


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;222.14&quot; width=&quot;1041.11&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;223.87&quot; width=&quot;1037.65&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;692.28&quot; x2=&quot;1329.45&quot; y2=&quot;692.28&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;548.23&quot; x2=&quot;1329.45&quot; y2=&quot;548.23&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;404.18&quot; x2=&quot;1329.45&quot; y2=&quot;404.18&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;260.13&quot; x2=&quot;1329.45&quot; y2=&quot;260.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;1551.58&quot; y=&quot;222.14&quot; width=&quot;1041.23&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1553.31&quot; y=&quot;223.87&quot; width=&quot;1037.78&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;657.26&quot; x2=&quot;2592.81&quot; y2=&quot;657.26&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;539.07&quot; x2=&quot;2592.81&quot; y2=&quot;539.07&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;420.76&quot; x2=&quot;2592.81&quot; y2=&quot;420.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;302.58&quot; x2=&quot;2592.81&quot; y2=&quot;302.58&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;2783.39&quot; y=&quot;222.14&quot; width=&quot;1041.11&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2785.11&quot; y=&quot;223.87&quot; width=&quot;1037.65&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;548.23&quot; x2=&quot;3824.49&quot; y2=&quot;548.23&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;404.18&quot; x2=&quot;3824.49&quot; y2=&quot;404.18&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;260.13&quot; x2=&quot;3824.49&quot; y2=&quot;260.13&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;886.20&quot; width=&quot;1041.11&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;887.93&quot; width=&quot;1037.65&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1356.47&quot; x2=&quot;1329.45&quot; y2=&quot;1356.47&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1248.43&quot; x2=&quot;1329.45&quot; y2=&quot;1248.43&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1140.39&quot; x2=&quot;1329.45&quot; y2=&quot;1140.39&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1032.23&quot; x2=&quot;1329.45&quot; y2=&quot;1032.23&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;924.19&quot; x2=&quot;1329.45&quot; y2=&quot;924.19&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;1551.58&quot; y=&quot;886.20&quot; width=&quot;1041.23&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1553.31&quot; y=&quot;887.93&quot; width=&quot;1037.78&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1327.26&quot; x2=&quot;2592.81&quot; y2=&quot;1327.26&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1193.36&quot; x2=&quot;2592.81&quot; y2=&quot;1193.36&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1059.58&quot; x2=&quot;2592.81&quot; y2=&quot;1059.58&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;925.80&quot; x2=&quot;2592.81&quot; y2=&quot;925.80&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;2783.39&quot; y=&quot;886.20&quot; width=&quot;1041.11&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2785.11&quot; y=&quot;887.93&quot; width=&quot;1037.65&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1356.47&quot; x2=&quot;3824.49&quot; y2=&quot;1356.47&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1248.43&quot; x2=&quot;3824.49&quot; y2=&quot;1248.43&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1140.39&quot; x2=&quot;3824.49&quot; y2=&quot;1140.39&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1032.23&quot; x2=&quot;3824.49&quot; y2=&quot;1032.23&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;924.19&quot; x2=&quot;3824.49&quot; y2=&quot;924.19&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;1658.30&quot; width=&quot;1041.11&quot; height=&quot;508.26&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;1660.03&quot; width=&quot;1037.65&quot; height=&quot;504.80&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2126.46&quot; x2=&quot;1329.45&quot; y2=&quot;2126.46&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2018.92&quot; x2=&quot;1329.45&quot; y2=&quot;2018.92&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1911.38&quot; x2=&quot;1329.45&quot; y2=&quot;1911.38&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1803.84&quot; x2=&quot;1329.45&quot; y2=&quot;1803.84&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1696.29&quot; x2=&quot;1329.45&quot; y2=&quot;1696.29&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_tMbiJpGo&quot; d=&quot;M326.33 549.84 L326.33 549.84 L374.59 405.29 L422.85 401.58 L471.12 426.08 L519.38 452.57 L567.64 476.58 L615.90 497.86 L664.17 516.80 L712.43 533.87 L760.69 549.34 L808.95 563.20 L857.22 575.70 L905.48 587.09 L953.74 597.24 L1002.00 606.52 L1050.14 614.81 L1098.40 622.36 L1146.67 629.16 L1194.93 635.23 L1243.19 640.80 L1291.45 645.75 L1291.45 621.00 L1243.19 613.45 L1194.93 604.91 L1146.67 595.50 L1098.40 585.11 L1050.14 573.60 L1002.00 560.85 L953.74 546.74 L905.48 531.03 L857.22 513.70 L808.95 494.52 L760.69 473.36 L712.43 449.85 L664.17 423.98 L615.90 395.52 L567.64 364.95 L519.38 333.52 L471.12 305.55 L422.85 293.42 L374.59 327.82 L326.33 489.69 L326.33 549.84 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_YnRRMdoA&quot;&gt;
			&lt;use xlink:href=&quot;#path_tMbiJpGo&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_tMbiJpGo&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_tMbiJpGo&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_YnRRMdoA)&quot;/&gt;
	&lt;path d=&quot; M326.33 519.89 L374.59 366.68 L422.85 347.62 L471.12 365.94 L519.38 393.04 L567.64 420.76 L615.78 446.75 L664.04 470.39 L712.30 491.92 L760.57 511.35 L808.83 528.92 L857.09 544.76 L905.36 559.12 L953.62 571.99 L1001.88 583.75 L1050.14 594.27 L1098.40 603.80 L1146.67 612.33 L1194.93 620.13 L1243.19 627.18 L1291.45 633.50&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_mxoQeqHi&quot; d=&quot;M1589.57 325.35 L1589.57 325.35 L1637.83 623.97 L1686.09 685.35 L1734.36 691.41 L1782.62 691.66 L1830.88 692.28 L1879.14 690.55 L1927.41 687.82 L1975.67 684.98 L2023.93 682.38 L2072.19 680.03 L2120.46 677.80 L2168.72 675.82 L2216.98 673.96 L2265.24 672.35 L2313.51 670.87 L2361.77 669.51 L2410.03 668.39 L2458.17 667.28 L2506.43 666.29 L2554.70 665.42 L2554.70 662.58 L2506.43 663.20 L2458.17 663.82 L2410.03 664.56 L2361.77 665.30 L2313.51 666.17 L2265.24 667.16 L2216.98 668.15 L2168.72 669.38 L2120.46 670.62 L2072.19 672.11 L2023.93 673.72 L1975.67 675.45 L1927.41 677.55 L1879.14 679.90 L1830.88 682.75 L1782.62 685.60 L1734.36 677.80 L1686.09 640.18 L1637.83 536.23 L1589.57 260.13 L1589.57 325.35 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_fPFUbRNU&quot;&gt;
			&lt;use xlink:href=&quot;#path_mxoQeqHi&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_mxoQeqHi&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_mxoQeqHi&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_fPFUbRNU)&quot;/&gt;
	&lt;path d=&quot; M1589.57 292.80 L1637.83 580.16 L1686.09 662.83 L1734.36 684.61 L1782.62 688.69 L1830.88 687.58 L1879.14 685.23 L1927.41 682.75 L1975.67 680.28 L2023.81 678.05 L2072.07 676.07 L2120.33 674.21 L2168.59 672.60 L2216.86 671.12 L2265.12 669.76 L2313.38 668.64 L2361.64 667.53 L2409.91 666.54 L2458.17 665.55 L2506.43 664.81 L2554.70 664.06&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_NJzAxSsC&quot; d=&quot;M2821.38 330.18 L2821.38 330.18 L2869.64 548.11 L2917.90 589.93 L2966.16 594.51 L3014.43 597.61 L3062.69 596.12 L3110.95 592.66 L3159.21 588.70 L3207.48 584.98 L3255.74 581.40 L3304.00 578.18 L3352.14 575.33 L3400.40 572.73 L3448.66 570.38 L3496.93 568.15 L3545.19 566.30 L3593.45 564.57 L3641.72 562.96 L3689.98 561.47 L3738.24 560.23 L3786.50 559.12 L3786.50 555.28 L3738.24 556.03 L3689.98 556.89 L3641.72 557.88 L3593.45 558.87 L3545.19 559.99 L3496.93 561.22 L3448.66 562.71 L3400.40 564.19 L3352.14 565.93 L3304.00 567.91 L3255.74 570.01 L3207.48 572.36 L3159.21 574.96 L3110.95 577.93 L3062.69 581.52 L3014.43 585.48 L2966.16 586.84 L2917.90 563.82 L2869.64 488.83 L2821.38 287.60 L2821.38 330.18 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_LFCjluSM&quot;&gt;
			&lt;use xlink:href=&quot;#path_NJzAxSsC&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_NJzAxSsC&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_NJzAxSsC&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_LFCjluSM)&quot;/&gt;
	&lt;path d=&quot; M2821.38 308.89 L2869.64 518.53 L2917.78 576.94 L2966.04 590.80 L3014.30 591.54 L3062.57 588.82 L3110.83 585.36 L3159.09 581.89 L3207.35 578.67 L3255.61 575.70 L3303.88 573.10 L3352.14 570.75 L3400.40 568.53 L3448.66 566.55 L3496.93 564.81 L3545.19 563.20 L3593.45 561.72 L3641.72 560.48 L3689.98 559.24 L3738.24 558.25 L3786.50 557.26&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_CKvdfdFi&quot; d=&quot;M326.33 1005.87 L326.33 1005.87 L374.59 1237.79 L422.85 1293.10 L471.12 1306.10 L519.38 1311.17 L567.64 1315.13 L615.90 1318.84 L664.17 1322.43 L712.43 1325.65 L760.69 1328.62 L808.95 1331.34 L857.22 1333.70 L905.48 1335.92 L953.74 1337.90 L1002.00 1339.76 L1050.14 1341.37 L1098.40 1342.85 L1146.67 1344.09 L1194.93 1345.33 L1243.19 1346.44 L1291.45 1347.43 L1291.45 1342.61 L1243.19 1341.12 L1194.93 1339.39 L1146.67 1337.53 L1098.40 1335.55 L1050.14 1333.32 L1002.00 1330.85 L953.74 1328.13 L905.48 1325.03 L857.22 1321.69 L808.95 1317.85 L760.69 1313.77 L712.43 1309.07 L664.17 1303.87 L615.90 1297.81 L567.64 1290.26 L519.38 1279.49 L471.12 1260.93 L422.85 1222.94 L374.59 1138.66 L326.33 937.19 L326.33 1005.87 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_XaWuRVfp&quot;&gt;
			&lt;use xlink:href=&quot;#path_CKvdfdFi&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_CKvdfdFi&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_CKvdfdFi&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_XaWuRVfp)&quot;/&gt;
	&lt;path d=&quot; M326.33 971.59 L374.59 1188.28 L422.85 1258.08 L471.12 1283.58 L519.38 1295.46 L567.64 1302.76 L615.78 1308.45 L664.04 1313.15 L712.30 1317.48 L760.57 1321.20 L808.83 1324.66 L857.09 1327.76 L905.36 1330.60 L953.62 1333.08 L1001.88 1335.30 L1050.14 1337.41 L1098.40 1339.26 L1146.67 1340.87 L1194.93 1342.48 L1243.19 1343.84 L1291.45 1345.08&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_ykyhkxTV&quot; d=&quot;M1589.57 990.28 L1589.57 990.28 L1637.83 1288.15 L1686.09 1349.91 L1734.36 1356.22 L1782.62 1355.48 L1830.88 1356.47 L1879.14 1355.10 L1927.41 1353.00 L1975.67 1350.65 L2023.93 1348.42 L2072.19 1346.44 L2120.46 1344.59 L2168.72 1342.85 L2216.98 1341.37 L2265.24 1340.01 L2313.51 1338.77 L2361.77 1337.66 L2410.03 1336.67 L2458.17 1335.68 L2506.43 1334.93 L2554.70 1334.19 L2554.70 1331.72 L2506.43 1332.21 L2458.17 1332.71 L2410.03 1333.32 L2361.77 1334.07 L2313.51 1334.69 L2265.24 1335.55 L2216.98 1336.42 L2168.72 1337.41 L2120.46 1338.52 L2072.19 1339.76 L2023.93 1341.12 L1975.67 1342.61 L1927.41 1344.34 L1879.14 1346.44 L1830.88 1348.92 L1782.62 1350.90 L1734.36 1341.12 L1686.09 1302.76 L1637.83 1198.93 L1589.57 924.19 L1589.57 990.28 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_bfhOVfvh&quot;&gt;
			&lt;use xlink:href=&quot;#path_ykyhkxTV&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_ykyhkxTV&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_ykyhkxTV&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_bfhOVfvh)&quot;/&gt;
	&lt;path d=&quot; M1589.57 957.36 L1637.83 1243.60 L1686.09 1326.39 L1734.36 1348.67 L1782.62 1353.25 L1830.88 1352.75 L1879.14 1350.90 L1927.41 1348.79 L1975.67 1346.69 L2023.81 1344.83 L2072.07 1343.10 L2120.33 1341.62 L2168.59 1340.25 L2216.86 1339.02 L2265.12 1337.78 L2313.38 1336.79 L2361.64 1335.92 L2409.91 1335.06 L2458.17 1334.31 L2506.43 1333.57 L2554.70 1332.95&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_nmqhGcYF&quot; d=&quot;M2821.38 1029.26 L2821.38 1029.26 L2869.64 1281.84 L2917.90 1337.41 L2966.16 1346.32 L3014.43 1347.31 L3062.69 1347.56 L3110.95 1348.05 L3159.21 1348.67 L3207.48 1349.41 L3255.74 1350.03 L3304.00 1350.65 L3352.14 1351.27 L3400.40 1351.76 L3448.66 1352.13 L3496.93 1352.63 L3545.19 1353.00 L3593.45 1353.37 L3641.72 1353.62 L3689.98 1353.87 L3738.24 1354.11 L3786.50 1354.36 L3786.50 1353.25 L3738.24 1352.88 L3689.98 1352.51 L3641.72 1352.13 L3593.45 1351.64 L3545.19 1351.14 L3496.93 1350.53 L3448.66 1349.91 L3400.40 1349.16 L3352.14 1348.42 L3304.00 1347.56 L3255.74 1346.57 L3207.48 1345.45 L3159.21 1344.09 L3110.95 1342.36 L3062.69 1339.39 L3014.43 1333.20 L2966.16 1318.60 L2917.90 1282.09 L2869.64 1192.37 L2821.38 965.53 L2821.38 1029.26 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_TSVZrmZN&quot;&gt;
			&lt;use xlink:href=&quot;#path_nmqhGcYF&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_nmqhGcYF&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_nmqhGcYF&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_TSVZrmZN)&quot;/&gt;
	&lt;path d=&quot; M2821.38 997.46 L2869.64 1237.17 L2917.78 1309.81 L2966.04 1332.58 L3014.30 1340.25 L3062.57 1343.47 L3110.83 1345.20 L3159.09 1346.44 L3207.35 1347.56 L3255.61 1348.42 L3303.88 1349.16 L3352.14 1349.91 L3400.40 1350.53 L3448.66 1351.14 L3496.93 1351.64 L3545.19 1352.13 L3593.45 1352.51 L3641.72 1352.88 L3689.98 1353.25 L3738.24 1353.62 L3786.50 1353.87&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;defs&gt;
		&lt;path id=&quot;path_ePWFqRsn&quot; d=&quot;M326.33 1812.38 L326.33 1812.38 L374.59 2067.06 L422.85 2121.14 L471.12 2128.57 L519.38 2128.20 L567.64 2127.33 L615.90 2126.83 L664.17 2126.59 L712.43 2126.46 L760.69 2126.46 L808.95 2126.46 L857.22 2126.46 L905.48 2126.46 L953.74 2126.46 L1002.00 2126.46 L1050.14 2126.46 L1098.40 2126.46 L1146.67 2126.46 L1194.93 2126.46 L1243.19 2126.46 L1291.45 2126.46 L1291.45 2126.46 L1243.19 2126.46 L1194.93 2126.46 L1146.67 2126.46 L1098.40 2126.46 L1050.14 2126.46 L1002.00 2126.46 L953.74 2126.46 L905.48 2126.46 L857.22 2126.46 L808.95 2126.46 L760.69 2126.46 L712.43 2126.34 L664.17 2126.22 L615.90 2125.60 L567.64 2123.99 L519.38 2119.28 L471.12 2106.04 L422.85 2070.53 L374.59 1980.31 L326.33 1751.12 L326.33 1812.38 Z&quot; stroke-linecap=&quot;round&quot;/&gt;
	&lt;/defs&gt;
	&lt;defs&gt;
		&lt;clipPath id=&quot;clipPath_aiGizvPp&quot;&gt;
			&lt;use xlink:href=&quot;#path_ePWFqRsn&quot;/&gt;
		&lt;/clipPath&gt;
	&lt;/defs&gt;
	&lt;use xlink:href=&quot;#path_ePWFqRsn&quot; style=&quot;stroke:none;fill:#C0C0C0&quot;/&gt;
	&lt;use xlink:href=&quot;#path_ePWFqRsn&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:10.37&quot; clip-path=&quot;url(#clipPath_aiGizvPp)&quot;/&gt;
	&lt;path d=&quot; M326.33 1781.81 L374.59 2023.75 L422.85 2095.90 L471.12 2117.43 L519.38 2123.86 L567.64 2125.72 L615.78 2126.34 L664.04 2126.46 L712.30 2126.46 L760.57 2126.46 L808.83 2126.59 L857.09 2126.59 L905.36 2126.59 L953.62 2126.59 L1001.88 2126.59 L1050.14 2126.59 L1098.40 2126.59 L1146.67 2126.59 L1194.93 2126.59 L1243.19 2126.59 L1291.45 2126.59&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.18&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;730.40&quot; x2=&quot;288.34&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;692.28&quot; x2=&quot;264.33&quot; y2=&quot;692.28&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;713.32&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;548.23&quot; x2=&quot;264.33&quot; y2=&quot;548.23&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;569.27&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;404.18&quot; x2=&quot;264.33&quot; y2=&quot;404.18&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;425.22&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;260.13&quot; x2=&quot;264.33&quot; y2=&quot;260.13&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;281.17&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;730.40&quot; x2=&quot;1551.58&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;657.26&quot; x2=&quot;1527.57&quot; y2=&quot;657.26&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;678.29&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;539.07&quot; x2=&quot;1527.57&quot; y2=&quot;539.07&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;560.11&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;420.76&quot; x2=&quot;1527.57&quot; y2=&quot;420.76&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;441.80&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;302.58&quot; x2=&quot;1527.57&quot; y2=&quot;302.58&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;323.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1.5&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;730.40&quot; x2=&quot;2783.39&quot; y2=&quot;222.14&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;692.28&quot; x2=&quot;2759.38&quot; y2=&quot;692.28&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;713.32&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;-2&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;548.23&quot; x2=&quot;2759.38&quot; y2=&quot;548.23&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;569.27&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;404.18&quot; x2=&quot;2759.38&quot; y2=&quot;404.18&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;425.22&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;260.13&quot; x2=&quot;2759.38&quot; y2=&quot;260.13&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;281.17&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1394.46&quot; x2=&quot;288.34&quot; y2=&quot;886.20&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1356.47&quot; x2=&quot;264.33&quot; y2=&quot;1356.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1377.50&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1248.43&quot; x2=&quot;264.33&quot; y2=&quot;1248.43&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1269.47&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1140.39&quot; x2=&quot;264.33&quot; y2=&quot;1140.39&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1161.43&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1032.23&quot; x2=&quot;264.33&quot; y2=&quot;1032.23&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1053.27&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1.5&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;924.19&quot; x2=&quot;264.33&quot; y2=&quot;924.19&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;945.23&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1394.46&quot; x2=&quot;1551.58&quot; y2=&quot;886.20&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1327.26&quot; x2=&quot;1527.57&quot; y2=&quot;1327.26&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;1348.30&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1193.36&quot; x2=&quot;1527.57&quot; y2=&quot;1193.36&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;1214.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1059.58&quot; x2=&quot;1527.57&quot; y2=&quot;1059.58&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;1080.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;20&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;925.80&quot; x2=&quot;1527.57&quot; y2=&quot;925.80&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1515.57&quot; y=&quot;946.84&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;30&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1394.46&quot; x2=&quot;2783.39&quot; y2=&quot;886.20&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1356.47&quot; x2=&quot;2759.38&quot; y2=&quot;1356.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;1377.50&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1248.43&quot; x2=&quot;2759.38&quot; y2=&quot;1248.43&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;1269.47&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1140.39&quot; x2=&quot;2759.38&quot; y2=&quot;1140.39&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;1161.43&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1032.23&quot; x2=&quot;2759.38&quot; y2=&quot;1032.23&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;1053.27&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;924.19&quot; x2=&quot;2759.38&quot; y2=&quot;924.19&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2747.37&quot; y=&quot;945.23&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2166.56&quot; x2=&quot;288.34&quot; y2=&quot;1658.30&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2126.46&quot; x2=&quot;264.33&quot; y2=&quot;2126.46&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;2147.50&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2018.92&quot; x2=&quot;264.33&quot; y2=&quot;2018.92&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;2039.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1911.38&quot; x2=&quot;264.33&quot; y2=&quot;1911.38&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1932.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1803.84&quot; x2=&quot;264.33&quot; y2=&quot;1803.84&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1824.87&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;1696.29&quot; x2=&quot;264.33&quot; y2=&quot;1696.29&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;252.33&quot; y=&quot;1717.33&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;end&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;1551.58&quot; y1=&quot;1394.46&quot; x2=&quot;2592.81&quot; y2=&quot;1394.46&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;1589.57&quot; y1=&quot;1394.46&quot; x2=&quot;1589.57&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1589.57&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1830.88&quot; y1=&quot;1394.46&quot; x2=&quot;1830.88&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1830.88&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;2072.19&quot; y1=&quot;1394.46&quot; x2=&quot;2072.19&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2072.19&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;2313.51&quot; y1=&quot;1394.46&quot; x2=&quot;2313.51&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2313.51&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;2554.70&quot; y1=&quot;1394.46&quot; x2=&quot;2554.70&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2554.70&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;line x1=&quot;2783.39&quot; y1=&quot;1394.46&quot; x2=&quot;3824.49&quot; y2=&quot;1394.46&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;2821.38&quot; y1=&quot;1394.46&quot; x2=&quot;2821.38&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2821.38&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3062.69&quot; y1=&quot;1394.46&quot; x2=&quot;3062.69&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3062.69&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;3304.00&quot; y1=&quot;1394.46&quot; x2=&quot;3304.00&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3304.00&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;3545.19&quot; y1=&quot;1394.46&quot; x2=&quot;3545.19&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3545.19&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;3786.50&quot; y1=&quot;1394.46&quot; x2=&quot;3786.50&quot; y2=&quot;1418.47&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3786.50&quot; y=&quot;1472.42&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;line x1=&quot;288.34&quot; y1=&quot;2166.56&quot; x2=&quot;1329.45&quot; y2=&quot;2166.56&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;line x1=&quot;326.33&quot; y1=&quot;2166.56&quot; x2=&quot;326.33&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;326.33&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;567.64&quot; y1=&quot;2166.56&quot; x2=&quot;567.64&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;567.64&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;808.95&quot; y1=&quot;2166.56&quot; x2=&quot;808.95&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;808.95&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;1050.14&quot; y1=&quot;2166.56&quot; x2=&quot;1050.14&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1050.14&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;1291.45&quot; y1=&quot;2166.56&quot; x2=&quot;1291.45&quot; y2=&quot;2190.57&quot; style=&quot;stroke:#000000;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;1291.45&quot; y=&quot;2244.52&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:60.02px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;135.39&quot; width=&quot;1041.11&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;137.11&quot; width=&quot;1037.65&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;808.95&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, c&lt;/text&gt;
	&lt;rect x=&quot;1551.58&quot; y=&quot;135.39&quot; width=&quot;1041.23&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;1553.31&quot; y=&quot;137.11&quot; width=&quot;1037.78&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2072.19&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, h&lt;/text&gt;
	&lt;rect x=&quot;2783.39&quot; y=&quot;135.39&quot; width=&quot;1041.11&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;2785.11&quot; y=&quot;137.11&quot; width=&quot;1037.65&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3304.00&quot; y=&quot;191.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, r&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;799.45&quot; width=&quot;1041.11&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;801.18&quot; width=&quot;1037.65&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;808.95&quot; y=&quot;855.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, w&lt;/text&gt;
	&lt;rect x=&quot;1551.58&quot; y=&quot;799.45&quot; width=&quot;1041.23&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;1553.31&quot; y=&quot;801.18&quot; width=&quot;1037.78&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;2072.19&quot; y=&quot;855.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, x&lt;/text&gt;
	&lt;rect x=&quot;2783.39&quot; y=&quot;799.45&quot; width=&quot;1041.11&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;2785.11&quot; y=&quot;801.18&quot; width=&quot;1037.65&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;3304.00&quot; y=&quot;855.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, y&lt;/text&gt;
	&lt;rect x=&quot;288.34&quot; y=&quot;1571.55&quot; width=&quot;1041.11&quot; height=&quot;86.75&quot; style=&quot;fill:#D9E6EB&quot;/&gt;
	&lt;rect x=&quot;290.07&quot; y=&quot;1573.28&quot; width=&quot;1037.65&quot; height=&quot;83.30&quot; style=&quot;fill:none;stroke:#D9E6EB;stroke-width:3.46&quot;/&gt;
	&lt;text x=&quot;808.95&quot; y=&quot;1628.06&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:65.96px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;est, z, z&lt;/text&gt;
	&lt;rect x=&quot;567.64&quot; y=&quot;2443.64&quot; width=&quot;2824.59&quot; height=&quot;186.50&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;570.52&quot; y=&quot;2446.52&quot; width=&quot;2818.83&quot; height=&quot;180.74&quot; style=&quot;fill:none;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;610.83&quot; y=&quot;2486.96&quot; width=&quot;374.34&quot; height=&quot;99.99&quot; style=&quot;fill:#C0C0C0&quot;/&gt;
	&lt;rect x=&quot;615.15&quot; y=&quot;2491.28&quot; width=&quot;365.70&quot; height=&quot;91.35&quot; style=&quot;fill:none;stroke:#C0C0C0;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;1526.46&quot; y1=&quot;2536.95&quot; x2=&quot;1900.80&quot; y2=&quot;2536.95&quot; style=&quot;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;text x=&quot;1045.19&quot; y=&quot;2571.98&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;95% CI&lt;/text&gt;
	&lt;text x=&quot;1960.82&quot; y=&quot;2571.98&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;impulse-response function (irf)&lt;/text&gt;
	&lt;text x=&quot;1980.00&quot; y=&quot;2379.06&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;text x=&quot;118.06&quot; y=&quot;2737.89&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:79.94px;fill:#000000&quot;&gt;Graphs by irfname, impulse, and response&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



    
    
    
    

### Sensitivity to alternative parameter values

In order to check the sensitivity of the DSGE model to alternative parameter values, we store the estimated/assume values, we change the parameter we are interested in, and then *solve* (not re-estimate) the model with the new value:


```stata
quietly dsgenl
matrix b = e(b)
```


```stata
matrix b[1,6] = 0.6
dsgenl (y = c + x) ///
     (F.k = x + (1-{delta})*k) ///
     (1/(c^{sigma}) = {beta}*(1/((F.c)^{sigma}))*(1-{delta}+r)) ///
     (h^{phi} = w/(c^{sigma})) ///
     (y = (k^{alpha})*((z*h)^(1-{alpha}))) ///
     (ln(F.z)={rhoz}*ln(z)) ///
     (r = {alpha}*y/k) ///   
     (w = (1-{alpha})*y/h) ///
     , observed(y) unobserved(c x r w h) exostate(z) endostate(k) ///
     from(b) solve noidencheck
```

    
    
    Solving at initial parameter vector ...
    
    First-order DSGE model
    
    Sample: 1955q1 - 2015q4                         Number of obs     =        244
    Log likelihood = -651.84841
    ------------------------------------------------------------------------------
                 |                 OIM
               y |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
    /structural  |
           delta |       .025          .        .       .            .           .
           sigma |          1          .        .       .            .           .
            beta |        .96          .        .       .            .           .
             phi |          1          .        .       .            .           .
           alpha |         .3          .        .       .            .           .
            rhoz |         .6          .        .       .            .           .
    -------------+----------------------------------------------------------------
          sd(e.z)|   3.206035          .                             .           .
    ------------------------------------------------------------------------------
    Note: Skipped identification check.
    Note: Model solved at specified parameters; maximization options ignored.
    

Here's the response of output to a technology shock, under the old and the new values of persistence. Naturally, the response of output is more persistent if the shock is more persistent.


```stata
irf create alt, step(20) replace
irf ograph (est z y irf) (alt z y irf)
```

    
    irfname alt not found in rbcirf.irf
    (file rbcirf.irf updated)
    


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;290.81&quot; y=&quot;100.86&quot; width=&quot;3568.21&quot; height=&quot;1992.81&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;293.69&quot; y=&quot;103.74&quot; width=&quot;3562.45&quot; height=&quot;1987.05&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2030.31&quot; x2=&quot;3859.02&quot; y2=&quot;2030.31&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1563.75&quot; x2=&quot;3859.02&quot; y2=&quot;1563.75&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1097.20&quot; x2=&quot;3859.02&quot; y2=&quot;1097.20&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;630.65&quot; x2=&quot;3859.02&quot; y2=&quot;630.65&quot; style=&quot;stroke:#EAF2F3;stroke-width:8.64&quot;/&gt;
	&lt;path d=&quot; M354.05 480.17 L526.18 1514.99 L698.20 1828.71 L870.33 1926.85 L1042.35 1960.39 L1214.48 1974.12 L1386.49 1981.67 L1558.63 1986.99 L1730.77 1991.45 L1902.78 1995.28 L2074.92 1998.75 L2246.93 2001.72 L2419.07 2004.57 L2591.08 2007.04 L2763.21 2009.27 L2935.23 2011.37 L3107.36 2013.23 L3279.38 2014.84 L3451.51 2016.45 L3623.52 2017.81 L3795.66 2018.92&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;path d=&quot; M354.05 513.45 L526.18 1070.35 L698.20 1409.31 L870.33 1616.97 L1042.35 1745.55 L1214.48 1826.36 L1386.49 1877.96 L1558.63 1911.87 L1730.77 1934.89 L1902.78 1950.98 L2074.92 1962.86 L2246.93 1971.89 L2419.07 1979.07 L2591.08 1984.89 L2763.21 1989.84 L2935.23 1994.05 L3107.36 1997.76 L3279.38 2001.10 L3451.51 2003.95 L3623.52 2006.67 L3795.66 2008.90&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2093.67&quot; x2=&quot;290.81&quot; y2=&quot;100.86&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2030.31&quot; x2=&quot;250.84&quot; y2=&quot;2030.31&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;2030.31&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,2030.31)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1563.75&quot; x2=&quot;250.84&quot; y2=&quot;1563.75&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1563.75&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1563.75)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;1097.20&quot; x2=&quot;250.84&quot; y2=&quot;1097.20&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;1097.20&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,1097.20)&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;630.65&quot; x2=&quot;250.84&quot; y2=&quot;630.65&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;630.65&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,630.65)&quot; text-anchor=&quot;middle&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;164.22&quot; x2=&quot;250.84&quot; y2=&quot;164.22&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;200.73&quot; y=&quot;164.22&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; transform=&quot;rotate(-90 200.73,164.22)&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;290.81&quot; y1=&quot;2093.67&quot; x2=&quot;3859.02&quot; y2=&quot;2093.67&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;354.17&quot; y1=&quot;2093.67&quot; x2=&quot;354.17&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;354.17&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1214.48&quot; y1=&quot;2093.67&quot; x2=&quot;1214.48&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1214.48&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;2074.92&quot; y1=&quot;2093.67&quot; x2=&quot;2074.92&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2074.92&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;2935.35&quot; y1=&quot;2093.67&quot; x2=&quot;2935.35&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2935.35&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;3795.66&quot; y1=&quot;2093.67&quot; x2=&quot;3795.66&quot; y2=&quot;2133.64&quot; style=&quot;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3795.66&quot; y=&quot;2223.62&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;2074.92&quot; y=&quot;2333.64&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;1473.24&quot; y=&quot;2418.27&quot; width=&quot;1203.34&quot; height=&quot;326.34&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1476.12&quot; y=&quot;2421.15&quot; width=&quot;1197.59&quot; height=&quot;320.58&quot; style=&quot;fill:none;stroke:#000000;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1516.56&quot; y1=&quot;2511.46&quot; x2=&quot;1890.90&quot; y2=&quot;2511.46&quot; style=&quot;stroke:#1A476F;stroke-width:8.64&quot;/&gt;
	&lt;line x1=&quot;1516.56&quot; y1=&quot;2651.43&quot; x2=&quot;1890.90&quot; y2=&quot;2651.43&quot; style=&quot;stroke:#90353B;stroke-width:8.64&quot;/&gt;
	&lt;text x=&quot;1950.92&quot; y=&quot;2546.49&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; y&lt;/text&gt;
	&lt;text x=&quot;1950.92&quot; y=&quot;2686.46&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:99.99px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; y&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



    
    
    
    

And now we can do the same for all the other variables:


```stata
graph drop c h x k y r w z
foreach v in c h x k y r w z {
    quietly irf ograph (est z `v' irf) (alt z `v' irf), name(`v') nodraw
}
```


```stata
graph combine c h x k y r w z, rows(2) cols(4)
```


                <iframe frameborder="0" scrolling="no" height="436" width="600"                srcdoc="<html><body>&lt;?xml version=&quot;1.0&quot; encoding=&quot;UTF-8&quot; standalone=&quot;no&quot;?&gt;
&lt;!-- This is a Stata 16.0 generated SVG file (http://www.stata.com) --&gt;

&lt;svg version=&quot;1.1&quot; width=&quot;600px&quot; height=&quot;436px&quot; viewBox=&quot;0 0 3960 2880&quot; xmlns=&quot;http://www.w3.org/2000/svg&quot; xmlns:xlink=&quot;http://www.w3.org/1999/xlink&quot;&gt;
	&lt;desc&gt;Stata Graph - Graph&lt;/desc&gt;
	&lt;rect x=&quot;0&quot; y=&quot;0&quot; width=&quot;3960&quot; height=&quot;2880&quot; style=&quot;fill:#EAF2F3;stroke:none&quot;/&gt;
	&lt;rect x=&quot;0.00&quot; y=&quot;0.00&quot; width=&quot;3959.88&quot; height=&quot;2880.00&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2.88&quot; y=&quot;2.88&quot; width=&quot;3954.12&quot; height=&quot;2874.24&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;rect x=&quot;63.36&quot; y=&quot;63.36&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;65.28&quot; y=&quot;65.28&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;235.25&quot; y=&quot;130.56&quot; width=&quot;719.24&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;237.17&quot; y=&quot;132.48&quot; width=&quot;715.40&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;860.21&quot; x2=&quot;954.48&quot; y2=&quot;860.21&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;688.32&quot; x2=&quot;954.48&quot; y2=&quot;688.32&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;516.55&quot; x2=&quot;954.48&quot; y2=&quot;516.55&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;344.65&quot; x2=&quot;954.48&quot; y2=&quot;344.65&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;172.76&quot; x2=&quot;954.48&quot; y2=&quot;172.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M277.45 654.41 L309.13 471.75 L340.93 448.98 L372.61 470.88 L404.42 503.31 L436.10 536.23 L467.77 567.29 L499.58 595.50 L531.26 621.12 L563.06 644.39 L594.74 665.30 L626.55 684.24 L658.23 701.31 L690.03 716.66 L721.71 730.64 L753.51 743.27 L785.19 754.53 L817.00 764.80 L848.68 774.08 L880.36 782.50 L912.16 790.04&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M277.45 539.32 L309.13 360.50 L340.93 283.03 L372.61 263.47 L404.42 275.97 L436.10 305.42 L467.77 342.92 L499.58 383.39 L531.26 423.73 L563.06 462.59 L594.74 499.10 L626.55 532.88 L658.23 563.82 L690.03 592.16 L721.71 617.90 L753.51 641.17 L785.19 662.33 L817.00 681.51 L848.68 698.84 L880.36 714.43 L912.16 728.66&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;902.41&quot; x2=&quot;235.25&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;860.21&quot; x2=&quot;208.52&quot; y2=&quot;860.21&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.21&quot; y=&quot;860.21&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.21,860.21)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;688.32&quot; x2=&quot;208.52&quot; y2=&quot;688.32&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.21&quot; y=&quot;688.32&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.21,688.32)&quot; text-anchor=&quot;middle&quot;&gt;.2&lt;/text&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;516.55&quot; x2=&quot;208.52&quot; y2=&quot;516.55&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.21&quot; y=&quot;516.55&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.21,516.55)&quot; text-anchor=&quot;middle&quot;&gt;.4&lt;/text&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;344.65&quot; x2=&quot;208.52&quot; y2=&quot;344.65&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.21&quot; y=&quot;344.65&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.21,344.65)&quot; text-anchor=&quot;middle&quot;&gt;.6&lt;/text&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;172.76&quot; x2=&quot;208.52&quot; y2=&quot;172.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.21&quot; y=&quot;172.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.21,172.76)&quot; text-anchor=&quot;middle&quot;&gt;.8&lt;/text&gt;
	&lt;line x1=&quot;235.25&quot; y1=&quot;902.41&quot; x2=&quot;954.48&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;277.45&quot; y1=&quot;902.41&quot; x2=&quot;277.45&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;277.45&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;436.22&quot; y1=&quot;902.41&quot; x2=&quot;436.22&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;436.22&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;594.87&quot; y1=&quot;902.41&quot; x2=&quot;594.87&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;594.87&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;753.51&quot; y1=&quot;902.41&quot; x2=&quot;753.51&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;753.51&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;912.16&quot; y1=&quot;902.41&quot; x2=&quot;912.16&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;912.16&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;594.87&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;258.27&quot; y=&quot;1118.86&quot; width=&quot;673.20&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;260.19&quot; y=&quot;1120.78&quot; width=&quot;669.36&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;287.10&quot; y1=&quot;1180.98&quot; x2=&quot;487.45&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;287.10&quot; y1=&quot;1287.66&quot; x2=&quot;487.45&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;535.71&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; c&lt;/text&gt;
	&lt;text x=&quot;535.71&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; c&lt;/text&gt;
	&lt;rect x=&quot;1021.68&quot; y=&quot;63.36&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1023.60&quot; y=&quot;65.28&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;1193.32&quot; y=&quot;130.56&quot; width=&quot;719.36&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1195.24&quot; y=&quot;132.48&quot; width=&quot;715.52&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;798.71&quot; x2=&quot;1912.68&quot; y2=&quot;798.71&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;595.63&quot; x2=&quot;1912.68&quot; y2=&quot;595.63&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;392.67&quot; x2=&quot;1912.68&quot; y2=&quot;392.67&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;189.71&quot; x2=&quot;1912.68&quot; y2=&quot;189.71&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1235.52 172.88 L1267.20 666.29 L1299.00 808.11 L1330.68 845.73 L1362.49 852.66 L1394.29 850.81 L1425.97 846.72 L1457.78 842.39 L1489.45 838.31 L1521.26 834.47 L1552.94 831.01 L1584.74 827.91 L1616.42 825.07 L1648.23 822.47 L1679.91 820.24 L1711.71 818.14 L1743.51 816.28 L1775.19 814.55 L1807.00 812.94 L1838.68 811.58 L1870.48 810.34&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1235.52 214.59 L1267.20 499.10 L1299.00 664.81 L1330.68 759.85 L1362.49 812.82 L1394.29 841.03 L1425.97 854.64 L1457.78 859.84 L1489.45 860.34 L1521.26 858.11 L1552.94 854.64 L1584.74 850.68 L1616.42 846.48 L1648.23 842.27 L1679.91 838.31 L1711.71 834.72 L1743.51 831.25 L1775.19 828.16 L1807.00 825.44 L1838.68 822.84 L1870.48 820.49&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;902.41&quot; x2=&quot;1193.32&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;798.71&quot; x2=&quot;1166.71&quot; y2=&quot;798.71&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1133.28&quot; y=&quot;798.71&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1133.28,798.71)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;595.63&quot; x2=&quot;1166.71&quot; y2=&quot;595.63&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1133.28&quot; y=&quot;595.63&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1133.28,595.63)&quot; text-anchor=&quot;middle&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;392.67&quot; x2=&quot;1166.71&quot; y2=&quot;392.67&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1133.28&quot; y=&quot;392.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1133.28,392.67)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;189.71&quot; x2=&quot;1166.71&quot; y2=&quot;189.71&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1133.28&quot; y=&quot;189.71&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1133.28,189.71)&quot; text-anchor=&quot;middle&quot;&gt;1.5&lt;/text&gt;
	&lt;line x1=&quot;1193.32&quot; y1=&quot;902.41&quot; x2=&quot;1912.68&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1235.52&quot; y1=&quot;902.41&quot; x2=&quot;1235.52&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1235.52&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1394.29&quot; y1=&quot;902.41&quot; x2=&quot;1394.29&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1394.29&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;1553.06&quot; y1=&quot;902.41&quot; x2=&quot;1553.06&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1553.06&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;1711.71&quot; y1=&quot;902.41&quot; x2=&quot;1711.71&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1711.71&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;1870.48&quot; y1=&quot;902.41&quot; x2=&quot;1870.48&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1870.48&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;1553.06&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;1216.34&quot; y=&quot;1118.86&quot; width=&quot;673.32&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1218.26&quot; y=&quot;1120.78&quot; width=&quot;669.48&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1245.17&quot; y1=&quot;1180.98&quot; x2=&quot;1445.15&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1245.17&quot; y1=&quot;1287.66&quot; x2=&quot;1445.15&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1493.17&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; h&lt;/text&gt;
	&lt;text x=&quot;1493.17&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; h&lt;/text&gt;
	&lt;rect x=&quot;1980.00&quot; y=&quot;63.36&quot; width=&quot;958.20&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1981.92&quot; y=&quot;65.28&quot; width=&quot;954.36&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;2152.26&quot; y=&quot;130.56&quot; width=&quot;718.74&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2154.18&quot; y=&quot;132.48&quot; width=&quot;714.90&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;813.56&quot; x2=&quot;2871.00&quot; y2=&quot;813.56&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;599.96&quot; x2=&quot;2871.00&quot; y2=&quot;599.96&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;386.36&quot; x2=&quot;2871.00&quot; y2=&quot;386.36&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;172.76&quot; x2=&quot;2871.00&quot; y2=&quot;172.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2194.46 223.13 L2226.14 680.03 L2257.82 812.20 L2289.50 847.84 L2321.30 855.14 L2352.98 854.27 L2384.66 851.18 L2416.47 847.84 L2448.15 844.62 L2479.83 841.65 L2511.51 838.93 L2543.31 836.45 L2574.99 834.22 L2606.67 832.24 L2638.47 830.39 L2670.15 828.78 L2701.83 827.29 L2733.64 825.93 L2765.32 824.82 L2797.00 823.71 L2828.68 822.72&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2194.46 259.26 L2226.14 520.88 L2257.82 673.96 L2289.50 762.32 L2321.30 812.20 L2352.98 839.30 L2384.66 853.04 L2416.47 858.85 L2448.15 860.34 L2479.83 859.22 L2511.51 856.87 L2543.31 853.90 L2574.99 850.81 L2606.67 847.59 L2638.47 844.62 L2670.15 841.77 L2701.83 839.17 L2733.64 836.70 L2765.32 834.47 L2797.00 832.49 L2828.68 830.64&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;902.41&quot; x2=&quot;2152.26&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;813.56&quot; x2=&quot;2125.53&quot; y2=&quot;813.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2092.22&quot; y=&quot;813.56&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2092.22,813.56)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;599.96&quot; x2=&quot;2125.53&quot; y2=&quot;599.96&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2092.22&quot; y=&quot;599.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2092.22,599.96)&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;386.36&quot; x2=&quot;2125.53&quot; y2=&quot;386.36&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2092.22&quot; y=&quot;386.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2092.22,386.36)&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;172.76&quot; x2=&quot;2125.53&quot; y2=&quot;172.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2092.22&quot; y=&quot;172.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2092.22,172.76)&quot; text-anchor=&quot;middle&quot;&gt;30&lt;/text&gt;
	&lt;line x1=&quot;2152.26&quot; y1=&quot;902.41&quot; x2=&quot;2871.00&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2194.46&quot; y1=&quot;902.41&quot; x2=&quot;2194.46&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2194.46&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2353.11&quot; y1=&quot;902.41&quot; x2=&quot;2353.11&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2353.11&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;2511.63&quot; y1=&quot;902.41&quot; x2=&quot;2511.63&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2511.63&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;2670.15&quot; y1=&quot;902.41&quot; x2=&quot;2670.15&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2670.15&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;2828.80&quot; y1=&quot;902.41&quot; x2=&quot;2828.80&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2828.80&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;2511.63&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;2175.28&quot; y=&quot;1118.86&quot; width=&quot;672.70&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2177.20&quot; y=&quot;1120.78&quot; width=&quot;668.87&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2204.11&quot; y1=&quot;1180.98&quot; x2=&quot;2405.33&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2204.11&quot; y1=&quot;1287.66&quot; x2=&quot;2405.33&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2453.84&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; x&lt;/text&gt;
	&lt;text x=&quot;2453.84&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; x&lt;/text&gt;
	&lt;rect x=&quot;2938.20&quot; y=&quot;63.36&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2940.12&quot; y=&quot;65.28&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;3109.96&quot; y=&quot;130.56&quot; width=&quot;719.36&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3111.88&quot; y=&quot;132.48&quot; width=&quot;715.52&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;860.21&quot; x2=&quot;3829.32&quot; y2=&quot;860.21&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;631.02&quot; x2=&quot;3829.32&quot; y2=&quot;631.02&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;401.95&quot; x2=&quot;3829.32&quot; y2=&quot;401.95&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3152.04 860.34 L3183.84 543.53 L3215.52 479.79 L3247.32 488.58 L3279.13 516.18 L3310.81 547.12 L3342.61 576.69 L3374.29 604.04 L3406.09 628.79 L3437.78 651.19 L3469.58 671.49 L3501.26 689.80 L3533.06 706.39 L3564.74 721.36 L3596.55 734.73 L3628.23 746.98 L3660.03 757.99 L3691.83 767.89 L3723.51 776.80 L3755.32 784.97 L3787.00 792.27&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3152.04 860.34 L3183.84 562.96 L3215.52 413.34 L3247.32 349.60 L3279.13 334.88 L3310.81 347.38 L3342.61 373.98 L3374.29 407.27 L3406.09 442.92 L3437.78 478.43 L3469.58 512.46 L3501.26 544.39 L3533.06 573.97 L3564.74 601.07 L3596.55 625.82 L3628.23 648.35 L3660.03 668.77 L3691.83 687.33 L3723.51 704.04 L3755.32 719.13 L3787.00 732.87&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;902.41&quot; x2=&quot;3109.96&quot; y2=&quot;130.56&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;860.21&quot; x2=&quot;3083.23&quot; y2=&quot;860.21&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3049.79&quot; y=&quot;860.21&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3049.79,860.21)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;631.02&quot; x2=&quot;3083.23&quot; y2=&quot;631.02&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3049.79&quot; y=&quot;631.02&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3049.79,631.02)&quot; text-anchor=&quot;middle&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;401.95&quot; x2=&quot;3083.23&quot; y2=&quot;401.95&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3049.79&quot; y=&quot;401.95&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3049.79,401.95)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;172.76&quot; x2=&quot;3083.23&quot; y2=&quot;172.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3049.79&quot; y=&quot;172.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3049.79,172.76)&quot; text-anchor=&quot;middle&quot;&gt;1.5&lt;/text&gt;
	&lt;line x1=&quot;3109.96&quot; y1=&quot;902.41&quot; x2=&quot;3829.32&quot; y2=&quot;902.41&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3152.16&quot; y1=&quot;902.41&quot; x2=&quot;3152.16&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3152.16&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3310.93&quot; y1=&quot;902.41&quot; x2=&quot;3310.93&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3310.93&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;3469.58&quot; y1=&quot;902.41&quot; x2=&quot;3469.58&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3469.58&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;3628.35&quot; y1=&quot;902.41&quot; x2=&quot;3628.35&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3628.35&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;3787.12&quot; y1=&quot;902.41&quot; x2=&quot;3787.12&quot; y2=&quot;929.14&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3787.12&quot; y=&quot;989.14&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;3469.58&quot; y=&quot;1062.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;3132.98&quot; y=&quot;1118.86&quot; width=&quot;673.32&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3134.90&quot; y=&quot;1120.78&quot; width=&quot;669.48&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3161.69&quot; y1=&quot;1180.98&quot; x2=&quot;3361.79&quot; y2=&quot;1180.98&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3161.69&quot; y1=&quot;1287.66&quot; x2=&quot;3361.79&quot; y2=&quot;1287.66&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3409.81&quot; y=&quot;1204.36&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; k&lt;/text&gt;
	&lt;text x=&quot;3409.81&quot; y=&quot;1311.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; k&lt;/text&gt;
	&lt;rect x=&quot;63.36&quot; y=&quot;1440.00&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;65.28&quot; y=&quot;1441.92&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;235.50&quot; y=&quot;1507.20&quot; width=&quot;718.99&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;237.42&quot; y=&quot;1509.12&quot; width=&quot;715.15&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;2236.85&quot; x2=&quot;954.48&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;2064.96&quot; x2=&quot;954.48&quot; y2=&quot;2064.96&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;1893.19&quot; x2=&quot;954.48&quot; y2=&quot;1893.19&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;1721.29&quot; x2=&quot;954.48&quot; y2=&quot;1721.29&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M277.57 1665.85 L309.38 2047.01 L341.06 2162.60 L372.86 2198.86 L404.54 2211.11 L436.22 2216.18 L468.02 2219.03 L499.70 2221.01 L531.38 2222.62 L563.19 2223.98 L594.87 2225.22 L626.67 2226.33 L658.35 2227.45 L690.03 2228.31 L721.83 2229.18 L753.51 2229.92 L785.19 2230.54 L817.00 2231.16 L848.68 2231.78 L880.48 2232.27 L912.16 2232.77&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M277.57 1678.10 L309.38 1883.29 L341.06 2008.15 L372.86 2084.63 L404.54 2132.03 L436.22 2161.73 L468.02 2180.79 L499.70 2193.29 L531.38 2201.71 L563.19 2207.65 L594.87 2211.98 L626.67 2215.32 L658.35 2218.04 L690.03 2220.14 L721.83 2222.00 L753.51 2223.61 L785.19 2224.97 L817.00 2226.08 L848.68 2227.20 L880.48 2228.19 L912.16 2229.05&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;2279.05&quot; x2=&quot;235.50&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;2236.85&quot; x2=&quot;208.77&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.33&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.33,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;2064.96&quot; x2=&quot;208.77&quot; y2=&quot;2064.96&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.33&quot; y=&quot;2064.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.33,2064.96)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;1893.19&quot; x2=&quot;208.77&quot; y2=&quot;1893.19&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.33&quot; y=&quot;1893.19&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.33,1893.19)&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;1721.29&quot; x2=&quot;208.77&quot; y2=&quot;1721.29&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.33&quot; y=&quot;1721.29&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.33,1721.29)&quot; text-anchor=&quot;middle&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;1549.40&quot; x2=&quot;208.77&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;175.33&quot; y=&quot;1549.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 175.33,1549.40)&quot; text-anchor=&quot;middle&quot;&gt;4&lt;/text&gt;
	&lt;line x1=&quot;235.50&quot; y1=&quot;2279.05&quot; x2=&quot;954.48&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;277.69&quot; y1=&quot;2279.05&quot; x2=&quot;277.69&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;277.69&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;436.34&quot; y1=&quot;2279.05&quot; x2=&quot;436.34&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;436.34&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;594.99&quot; y1=&quot;2279.05&quot; x2=&quot;594.99&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;594.99&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;753.64&quot; y1=&quot;2279.05&quot; x2=&quot;753.64&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;753.64&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;912.16&quot; y1=&quot;2279.05&quot; x2=&quot;912.16&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;912.16&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;594.99&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;258.51&quot; y=&quot;2495.50&quot; width=&quot;672.95&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;260.43&quot; y=&quot;2497.42&quot; width=&quot;669.11&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;287.22&quot; y1=&quot;2557.62&quot; x2=&quot;488.19&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;287.22&quot; y1=&quot;2664.30&quot; x2=&quot;488.19&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;536.46&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; y&lt;/text&gt;
	&lt;text x=&quot;536.46&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; y&lt;/text&gt;
	&lt;rect x=&quot;1021.68&quot; y=&quot;1440.00&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1023.60&quot; y=&quot;1441.92&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;1194.81&quot; y=&quot;1507.20&quot; width=&quot;717.87&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1196.73&quot; y=&quot;1509.12&quot; width=&quot;714.03&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;2236.85&quot; x2=&quot;1912.68&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;2077.83&quot; x2=&quot;1912.68&quot; y2=&quot;2077.83&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;1918.80&quot; x2=&quot;1912.68&quot; y2=&quot;1918.80&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;1759.78&quot; x2=&quot;1912.68&quot; y2=&quot;1759.78&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;1600.76&quot; x2=&quot;1912.68&quot; y2=&quot;1600.76&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1237.01 1549.52 L1268.56 2012.11 L1300.24 2141.19 L1331.92 2171.63 L1363.60 2173.49 L1395.28 2167.43 L1426.96 2159.75 L1458.64 2152.08 L1490.32 2144.90 L1522.00 2138.47 L1553.68 2132.65 L1585.36 2127.33 L1617.04 2122.50 L1648.72 2118.17 L1680.40 2114.21 L1712.08 2110.75 L1743.76 2107.53 L1775.44 2104.68 L1807.12 2102.08 L1838.80 2099.73 L1870.48 2097.63&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M1237.01 1560.91 L1268.56 1853.83 L1300.24 2021.27 L1331.92 2114.21 L1363.60 2163.09 L1395.28 2186.36 L1426.96 2194.65 L1458.64 2194.65 L1490.32 2190.20 L1522.00 2183.39 L1553.68 2175.59 L1585.36 2167.55 L1617.04 2159.75 L1648.72 2152.33 L1680.40 2145.40 L1712.08 2139.09 L1743.76 2133.27 L1775.44 2127.95 L1807.12 2123.12 L1838.80 2118.79 L1870.48 2114.83&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;2279.05&quot; x2=&quot;1194.81&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;2236.85&quot; x2=&quot;1168.08&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1134.64&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1134.64,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;-1&lt;/text&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;2077.83&quot; x2=&quot;1168.08&quot; y2=&quot;2077.83&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1134.64&quot; y=&quot;2077.83&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1134.64,2077.83)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;1918.80&quot; x2=&quot;1168.08&quot; y2=&quot;1918.80&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1134.64&quot; y=&quot;1918.80&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1134.64,1918.80)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;1759.78&quot; x2=&quot;1168.08&quot; y2=&quot;1759.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1134.64&quot; y=&quot;1759.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1134.64,1759.78)&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;1600.76&quot; x2=&quot;1168.08&quot; y2=&quot;1600.76&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1134.64&quot; y=&quot;1600.76&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 1134.64,1600.76)&quot; text-anchor=&quot;middle&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;1194.81&quot; y1=&quot;2279.05&quot; x2=&quot;1912.68&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1237.01&quot; y1=&quot;2279.05&quot; x2=&quot;1237.01&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1237.01&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;1395.40&quot; y1=&quot;2279.05&quot; x2=&quot;1395.40&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1395.40&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;1553.81&quot; y1=&quot;2279.05&quot; x2=&quot;1553.81&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1553.81&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;1712.08&quot; y1=&quot;2279.05&quot; x2=&quot;1712.08&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1712.08&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;1870.48&quot; y1=&quot;2279.05&quot; x2=&quot;1870.48&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;1870.48&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;1553.81&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;1217.82&quot; y=&quot;2495.50&quot; width=&quot;671.84&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;1219.74&quot; y=&quot;2497.42&quot; width=&quot;668.00&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;1246.66&quot; y1=&quot;2557.62&quot; x2=&quot;1449.73&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;1246.66&quot; y1=&quot;2664.30&quot; x2=&quot;1449.73&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;1498.49&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; r&lt;/text&gt;
	&lt;text x=&quot;1498.49&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; r&lt;/text&gt;
	&lt;rect x=&quot;1980.00&quot; y=&quot;1440.00&quot; width=&quot;958.20&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;1981.92&quot; y=&quot;1441.92&quot; width=&quot;954.36&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;2150.40&quot; y=&quot;1507.20&quot; width=&quot;720.60&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2152.32&quot; y=&quot;1509.12&quot; width=&quot;716.76&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;2236.85&quot; x2=&quot;2871.00&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;2064.96&quot; x2=&quot;2871.00&quot; y2=&quot;2064.96&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;1893.19&quot; x2=&quot;2871.00&quot; y2=&quot;1893.19&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;1721.29&quot; x2=&quot;2871.00&quot; y2=&quot;1721.29&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;1549.40&quot; x2=&quot;2871.00&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2192.60 1624.64 L2224.41 1969.30 L2256.21 2080.30 L2288.01 2120.89 L2319.82 2139.70 L2351.62 2151.46 L2383.43 2160.37 L2415.23 2168.04 L2447.03 2174.73 L2478.84 2180.79 L2510.64 2186.24 L2542.44 2191.19 L2574.25 2195.64 L2606.05 2199.60 L2637.86 2203.19 L2669.66 2206.53 L2701.46 2209.50 L2733.27 2212.10 L2765.07 2214.58 L2796.87 2216.68 L2828.68 2218.66&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M2192.60 1614.00 L2224.41 1783.29 L2256.21 1892.69 L2288.01 1965.34 L2319.82 2015.21 L2351.62 2050.85 L2383.43 2077.33 L2415.23 2097.88 L2447.03 2114.46 L2478.84 2128.20 L2510.64 2139.83 L2542.44 2149.85 L2574.25 2158.76 L2606.05 2166.56 L2637.86 2173.49 L2669.66 2179.80 L2701.46 2185.37 L2733.27 2190.32 L2765.07 2194.90 L2796.87 2198.98 L2828.68 2202.70&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;2279.05&quot; x2=&quot;2150.40&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;2236.85&quot; x2=&quot;2123.80&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2090.36&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2090.36,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;2064.96&quot; x2=&quot;2123.80&quot; y2=&quot;2064.96&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2090.36&quot; y=&quot;2064.96&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2090.36,2064.96)&quot; text-anchor=&quot;middle&quot;&gt;.5&lt;/text&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;1893.19&quot; x2=&quot;2123.80&quot; y2=&quot;1893.19&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2090.36&quot; y=&quot;1893.19&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2090.36,1893.19)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;1721.29&quot; x2=&quot;2123.80&quot; y2=&quot;1721.29&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2090.36&quot; y=&quot;1721.29&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2090.36,1721.29)&quot; text-anchor=&quot;middle&quot;&gt;1.5&lt;/text&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;1549.40&quot; x2=&quot;2123.80&quot; y2=&quot;1549.40&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2090.36&quot; y=&quot;1549.40&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 2090.36,1549.40)&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;2150.40&quot; y1=&quot;2279.05&quot; x2=&quot;2871.00&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2192.60&quot; y1=&quot;2279.05&quot; x2=&quot;2192.60&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2192.60&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;2351.74&quot; y1=&quot;2279.05&quot; x2=&quot;2351.74&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2351.74&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;2510.76&quot; y1=&quot;2279.05&quot; x2=&quot;2510.76&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2510.76&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;2669.78&quot; y1=&quot;2279.05&quot; x2=&quot;2669.78&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2669.78&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;2828.80&quot; y1=&quot;2279.05&quot; x2=&quot;2828.80&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;2828.80&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;2510.76&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;2173.42&quot; y=&quot;2495.50&quot; width=&quot;674.56&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;2175.34&quot; y=&quot;2497.42&quot; width=&quot;670.72&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;2202.26&quot; y1=&quot;2557.62&quot; x2=&quot;2399.64&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;2202.26&quot; y1=&quot;2664.30&quot; x2=&quot;2399.64&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;2447.03&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; w&lt;/text&gt;
	&lt;text x=&quot;2447.03&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; w&lt;/text&gt;
	&lt;rect x=&quot;2938.20&quot; y=&quot;1440.00&quot; width=&quot;958.32&quot; height=&quot;1376.64&quot; style=&quot;fill:#EAF2F3&quot;/&gt;
	&lt;rect x=&quot;2940.12&quot; y=&quot;1441.92&quot; width=&quot;954.48&quot; height=&quot;1372.80&quot; style=&quot;fill:none;stroke:#EAF2F3;stroke-width:3.84&quot;/&gt;
	&lt;rect x=&quot;3110.33&quot; y=&quot;1507.20&quot; width=&quot;718.99&quot; height=&quot;771.85&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3112.25&quot; y=&quot;1509.12&quot; width=&quot;715.15&quot; height=&quot;768.01&quot; style=&quot;fill:none;stroke:#FFFFFF;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;2236.85&quot; x2=&quot;3829.32&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;2022.39&quot; x2=&quot;3829.32&quot; y2=&quot;2022.39&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;1808.04&quot; x2=&quot;3829.32&quot; y2=&quot;1808.04&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;1593.58&quot; x2=&quot;3829.32&quot; y2=&quot;1593.58&quot; style=&quot;stroke:#EAF2F3;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3152.53 1549.65 L3184.21 2032.04 L3215.89 2175.84 L3247.70 2218.66 L3279.38 2231.53 L3311.18 2235.24 L3342.86 2236.48 L3374.54 2236.73 L3406.34 2236.85 L3438.02 2236.85 L3469.70 2236.97 L3501.51 2236.97 L3533.19 2236.97 L3564.99 2236.97 L3596.67 2236.97 L3628.35 2236.97 L3660.15 2236.97 L3691.83 2236.97 L3723.51 2236.97 L3755.32 2236.97 L3787.00 2236.97&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;path d=&quot; M3152.53 1549.52 L3184.21 1824.50 L3215.89 1989.47 L3247.70 2088.47 L3279.38 2147.87 L3311.18 2183.51 L3342.86 2204.80 L3374.54 2217.67 L3406.34 2225.34 L3438.02 2230.04 L3469.70 2232.77 L3501.51 2234.38 L3533.19 2235.37 L3564.99 2235.98 L3596.67 2236.36 L3628.35 2236.60 L3660.15 2236.73 L3691.83 2236.85 L3723.51 2236.85 L3755.32 2236.85 L3787.00 2236.85&quot; stroke-linejoin=&quot;round&quot; style=&quot;fill:none;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;2279.05&quot; x2=&quot;3110.33&quot; y2=&quot;1507.20&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;2236.85&quot; x2=&quot;3083.60&quot; y2=&quot;2236.85&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3050.29&quot; y=&quot;2236.85&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3050.29,2236.85)&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;2022.39&quot; x2=&quot;3083.60&quot; y2=&quot;2022.39&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3050.29&quot; y=&quot;2022.39&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3050.29,2022.39)&quot; text-anchor=&quot;middle&quot;&gt;1&lt;/text&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;1808.04&quot; x2=&quot;3083.60&quot; y2=&quot;1808.04&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3050.29&quot; y=&quot;1808.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3050.29,1808.04)&quot; text-anchor=&quot;middle&quot;&gt;2&lt;/text&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;1593.58&quot; x2=&quot;3083.60&quot; y2=&quot;1593.58&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3050.29&quot; y=&quot;1593.58&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; transform=&quot;rotate(-90 3050.29,1593.58)&quot; text-anchor=&quot;middle&quot;&gt;3&lt;/text&gt;
	&lt;line x1=&quot;3110.33&quot; y1=&quot;2279.05&quot; x2=&quot;3829.32&quot; y2=&quot;2279.05&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3152.53&quot; y1=&quot;2279.05&quot; x2=&quot;3152.53&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3152.53&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;0&lt;/text&gt;
	&lt;line x1=&quot;3311.18&quot; y1=&quot;2279.05&quot; x2=&quot;3311.18&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3311.18&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;5&lt;/text&gt;
	&lt;line x1=&quot;3469.83&quot; y1=&quot;2279.05&quot; x2=&quot;3469.83&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3469.83&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;10&lt;/text&gt;
	&lt;line x1=&quot;3628.47&quot; y1=&quot;2279.05&quot; x2=&quot;3628.47&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3628.47&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;15&lt;/text&gt;
	&lt;line x1=&quot;3787.12&quot; y1=&quot;2279.05&quot; x2=&quot;3787.12&quot; y2=&quot;2305.78&quot; style=&quot;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;text x=&quot;3787.12&quot; y=&quot;2365.78&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;20&lt;/text&gt;
	&lt;text x=&quot;3469.83&quot; y=&quot;2439.04&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot; text-anchor=&quot;middle&quot;&gt;step&lt;/text&gt;
	&lt;rect x=&quot;3133.35&quot; y=&quot;2495.50&quot; width=&quot;672.95&quot; height=&quot;230.92&quot; style=&quot;fill:#FFFFFF&quot;/&gt;
	&lt;rect x=&quot;3135.27&quot; y=&quot;2497.42&quot; width=&quot;669.11&quot; height=&quot;227.08&quot; style=&quot;fill:none;stroke:#000000;stroke-width:3.84&quot;/&gt;
	&lt;line x1=&quot;3162.18&quot; y1=&quot;2557.62&quot; x2=&quot;3363.03&quot; y2=&quot;2557.62&quot; style=&quot;stroke:#1A476F;stroke-width:5.76&quot;/&gt;
	&lt;line x1=&quot;3162.18&quot; y1=&quot;2664.30&quot; x2=&quot;3363.03&quot; y2=&quot;2664.30&quot; style=&quot;stroke:#90353B;stroke-width:5.76&quot;/&gt;
	&lt;text x=&quot;3411.29&quot; y=&quot;2581.00&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;est: irf of z -&amp;gt; z&lt;/text&gt;
	&lt;text x=&quot;3411.29&quot; y=&quot;2687.67&quot; style=&quot;font-family:&#x27;Helvetica&#x27;;font-size:66.70px;fill:#000000&quot;&gt;alt: irf of z -&amp;gt; z&lt;/text&gt;
&lt;/svg&gt;
</body></html>"></iframe>



## References

<a id='references'></a>
King, R. G., and S. T. Rebelo (1999), Resuscitating real business cycles, in: Taylor, J.B., and M. Woodford (eds.), *Handbook of Macroeconomics*, Volume 1A. New York: Elsevier.

Uhlig, H. (1999), A toolkit for analysing nonlinear dynamic stochastic models easily, in: Marion, R., and A. Scott(eds.), *Computational methods for the study of dynamic economies*. New York: Oxford University Press.
