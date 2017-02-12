nhp=new.env()
source('../R/nhpoisp/nhpoisp.R', local=nhp)

events <- c(32000, 38500, 45000, 55000, 62500, 70000, 1000000, 2000000, 4500000) # Cadell
#events <- c(99000)
#events <- c(30000, 50000, 90000, 1000000)
#events <- c(5500, 11500, 32500, 69000, 101000) #Galeen (Roer Valley)
#events <- c(440, 1100, 1800, 2400, 99000) # Rurrand (Roer Valley)
#events <- c(365, 526, 530, 601, 692, 733, 764, 794, 824, 840, 1208, 1237, 1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumatra (Mentawai Segment) from Philibosian et al 2017, Figures 19,20)
#events <- c(1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumatra (Mentawai Segment) from Philibosian et al 2017, Figures 20)
#events <-c(1314, 1350, 1388, 1597, 1613, 1703, 1797, 1833, 2007) #site BLS/PSR#
#events <- events - 1313
# Below is a statistical check on whether the time spacing between 'events'
# is exponential. It does not rely on fitting a model, so is a good 'independent
# check'

# Re-order to make time run the right way
events <- rev(events)
events <- (events - events[1])*(-1)
print(events)

## Suppose the underlying distribution of 'diff(events)' is exponential. Consider the sampling
## distribution of (sample-variance/sample-mean^2). This will not depend on the
## true mean, as illustrated below. 
f<-function(rate=1){
    # Random value of var(data)/mean(data)^2 for exponentially distributed samples
    dat = rexp(length(diff(events)), rate=rate)
    var(dat)/mean(dat)^2
}
# Claim: f() distribution does not depend on the rate
# If true, the qq-plots below should appear as 1:1 lines
pdf('poisson_test.pdf')
random_stat = replicate(10000, f(rate=1))
random_statB = replicate(10000, f(rate=100))
random_statC = replicate(10000, f(rate=0.1))
par(mfrow=c(2,1))
qqplot(random_stat, random_statB)
abline(0,1,col='red')
qqplot(random_stat, random_statC)
abline(0,1,col='red')
##
## Accepting the above theory....
##
## Let's see if the data is 'unusual' given the null hypothesis, based
## on this measure. Note that strongly clustered data should have 'very large' values
## of this statistic, as compared with the null hypothesis samples.
p_value = 1 - mean(random_stat < var(diff(events))/mean(diff(events))^2)
print(p_value)
# ~ 0.23 -- i.e. not unusual.
dev.off()