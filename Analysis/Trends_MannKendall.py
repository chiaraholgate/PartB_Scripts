# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:47:34 2019



20/3/19
PROBABLY BETTER TO USE AN R-PACKAGE LIKE THIS ONE:
https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf
SINCE SEASONAL AUTOCORRELATION IS TAKEN CARE OF.

THERE ARE SOME ISSUES WITH THE M-K METHOD BELOW THAT HAVEN'T YET
BEEN FIXED.



This script calculates trends using the Mann-Kendall test. The first option calculates the Mann-Kendall 
statistic and whether there is a trend or not, and if it is significant.

The second part just calculates Kendall's tau, wihch is a non-parameteric measure of association, 
and it's signficance. 

If the Mann-Kendall statistic inidcates there is a trend, then also calculate Theil Sen's slope. I haven't 
done this yet since there haven't been any significant trends.

* Scipy version described here: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.kendalltau.html
* Another option here: https://up-rs-esp.github.io/mkt/

* I'm not sure how scipy calulates the p-value. In any case, a p-value of 0.8 means the probability that your
tau value is wrong is 80%, so you can only be 20% confident in rejecting H0.
* You can calculate a z-value, which is just the area under a bell-curve for a given standard deviation away from
the mean. If you want to be 95% confident in rejecting Ho, you need to have a z-value that is greater than
1.96. Other critical z-values can be found here: https://www.statisticshowto.datasciencecentral.com/probability-and-statistics/find-critical-values/
* I'm not sure how to calculate z-values when you don't want to assume a normal distribution. I think that's 
the t-test? Although below website says that for N>20, the t-distribtuion looks almost identical to normal.
https://www.statisticshowto.datasciencecentral.com/probability-and-statistics/t-distribution/

@author: z3131380
"""
import pandas as pd
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

region = 'MDB'
seasons = ['DJF','MAM','JJA','SON']
variables = ['RR','region_land_mm','P_total','outregion_land_Aus_pct','outregion_land_Aus_mm',\
        'outregion_land_outside_Aus_pct','outregion_land_outside_Aus_mm','outregion_ocean_pct',\
        'outregion_ocean_mm','outregion_ocean_km']

#==============================================================================
# Mann_kendall test from https://michaelpaulschramm.com/2015/08/01/simple-time-series-trend-analysis/
#==============================================================================

def mk_test(x, alpha = 0.05):  
    """   
    Input:
        x:   a vector of data
        alpha: significance level (0.05 default)

    Output:
        trend: tells the trend (increasing, decreasing or no trend)
        h: True (if trend is present) or False (if trend is absence)
        p: p value of the significance test
        z: normalized test statistics 

    Examples
    --------
      >>> x = np.random.rand(100)
      >>> trend,h,p,z = mk_test(x,0.05) 
    """
    n = len(x)

    # calculate S 
    s = 0
    for k in range(n-1):
        for j in range(k+1,n):
            s += np.sign(x[j] - x[k])
            
    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g: # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18
    else: # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(unique_x[i] == x)
        var_s = (n*(n-1)*(2*n+5) + np.sum(tp*(tp-1)*(2*tp+5)))/18

    if s>0:
        z = (s - 1)/np.sqrt(var_s)
    elif s == 0:
            z = 0
    elif s<0:
        z = (s + 1)/np.sqrt(var_s)

    # calculate the p_value
    p = 2*(1-norm.cdf(abs(z))) # two tail test
    h = abs(z) > norm.ppf(1-alpha/2) 

    if (z<0) and h:
        trend = 'decreasing'
    elif (z>0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'

    return trend, h, p, z



###

df = pd.read_csv('/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/Seasonal/'+region+'_seasonal_rainfall_recycling_1979-2013.csv')


for v in variables:
    df_var = pd.DataFrame(df[v])
    # Since sm can't handle missing values, delete first row of df
    df_var.drop(df_var.index[0],inplace=True) # only do this for DJF

    # Convert index to a datetime
    df_var.reset_index(inplace=True)
    df_var['Season_start'] =pd.date_range(start='1979-3-1',end='2013-9-1',freq='3MS')
    df_var = df_var.set_index('Season_start')
    
    res = sm.tsa.seasonal_decompose(df_var[v])
    #res.plot()
    
    # Can access each component
    #residual = res.residual
#        seasonal = res.seasonal
    trend = res.trend
#        plt.plot(trend)
    
    trend_nonans = trend[2:-2]
    test_trend,h,p,z = mk_test(trend_nonans,alpha=0.05) 
    print v,' = ',test_trend, h, z, p 


#==============================================================================
# Kendall's Tau
#==============================================================================

#df = pd.read_csv('/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/Seasonal/'+region+'_seasonal_rainfall_recycling_1979-2013.csv')
#groups = df.groupby('Season')
#Years = np.arange(1979,2014)
#
#for s in seasons: 
#    variable = groups.get_group(s)['RR']
#    tau, p_value = stats.kendalltau(Years,variable)
#    # Calc. z from https://www.statisticshowto.datasciencecentral.com/kendalls-tau/
#    z = (3*tau*np.sqrt(len(variable)*(len(variable)-1)))/np.sqrt(2*(2*len(variable)+5))
#    print region,s, 'tau=',np.round(tau,2), 'prob. that tau is wrong is ',np.round(p_value*100 ,2)
#    if z<1.96:
#        print 'z=',np.round(z,2), 'so you cannot reject Ho'
#    else:
#        print 'z=',np.round(z,2), 'so you can reject Ho'
#    
#plt.plot(Years,variable)