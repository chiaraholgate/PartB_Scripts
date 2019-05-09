# Package outline here https://cran.r-project.org/web/packages/trend/vignettes/trend.pdf
# This script calculates the Mann-Kendall statistic using the correlated seasonal version of M-K.

workdir <- "/home/z3131380/hdrive/PhD/PartB/Scripts/Analysis/"

#20/3/19 package installation not yet working
install.packages("trend",dependencies=TRUE)

data_dir <- "/srv/ccrc/data19/z3131380/PartB/Output/Australia/100parcels/TS10min/exp01/Processed/Rainfall_recycling/Seasonal/"
mydata = read.csv(paste(data_dir,"MDB_seasonal_rainfall_recycling_1979-2013.csv",sep=""))  # read csv file

require(trend)
data(mydata)
P_total <- mydata[,"P_total"]

# Regular M-K test
mk.test(P_total)

# Seasonal M-K test
smk.test(P_total)

# Correlated seasonal M-K test
csmk.test(P_total)