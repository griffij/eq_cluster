nhp=new.env()
source('../R/nhpoisp/nhpoisp.R', local=nhp)

#events <- c(30000, 50000, 90000, 1000000)
#events <- c(30000, 50000, 90000, 250000, 500000, 530000, 570000, 590000, 850000)
events <- c(20000, 32000, 45000, 54000, 62000, 70000, 1000000, 2000000, 4000000) # Cadell
#events <- c(365, 526, 530, 601, 692, 733, 764, 794, 824, 840, 1208, 1237, 1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumat\
#ra (Mentawai Segment) from Philibosian et al 2017, Figures 19,20)
#events <- events - 364

#thetas1 = seq(0, 1e-05, len=100)
thetas1 = seq(0, 1e-5, len=100)
thetas2 = seq(-1e-03, 1e-03, len=101)
rate_equation = 'theta[1] + theta[2]*exp(-(t-tlast)/14486)'
#thetas2 = seq(-12, -3, len=101)
#thetas2 = seq(-8, -1, len=101)
#rate_equation = 'theta[1] + 10**(theta[2])*exp(-(t-tlast)/30)'

negloglik_grid = matrix(NA, nrow=length(thetas1), ncol=length(thetas2))

for(j in 1:length(thetas2)){
    for(i in 1:length(thetas1)){
        negloglik_grid[i,j] = nhp$negloglik_from_theta(
            theta=c(thetas1[i], thetas2[j]),
            observed_data=events,
            x0 = 0,
            rate_equation=rate_equation,
            integration_dt = 100)
    }
}

pdf('grid_search.pdf')
image(thetas1, thetas2, negloglik_grid, col=rainbow(255)[1:200],
    main=paste0('Negative log likelihood for ', rate_equation, 
        '\n with minimum (point) and 95% bivariate CI') )
# Supposing the minimum is the true min, then an approx bivariate 95% CI is the
# region with likelihoods within qchisq(0.95, 2) (where the '2' comes from the fact
# we have 2 variables). For each variable individually, you can set that to '1'.
min_ind = which(negloglik_grid == min(negloglik_grid), arr.ind=TRUE)
ci_bivariate_95 = qchisq(0.95, 2)/2 # 
contour(thetas1, thetas2, negloglik_grid, col='black', add=TRUE, 
    level=(min(negloglik_grid)+ci_bivariate_95))
points(thetas1[min_ind[1]], thetas2[min_ind[2]], col='black', pch=19, cex=2)
dev.off()