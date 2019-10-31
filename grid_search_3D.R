nhp=new.env()
source('../stormwavecluster/R/nhpoisp/nhpoisp.R', local=nhp)
library('plot3D')
library('rgl')
#par(mar = c(2, 2, 2, 2))
#par(mfrow = c(1, 1))

#events <- c(30000, 50000, 90000, 1000000)
#events <- c(30000, 50000, 90000, 250000, 500000, 530000, 570000, 590000, 850000)
events <- c(20000, 32000, 45000, 54000, 62000, 70000, 1000000, 2000000, 4000000) # Cadell
#events <- c(15100, 18000, 21000, 24000, 55000) # D'n 4 events
#events <- c(15100, 17000,19000, 21000, 22500, 24000, 55000) # D'n 6 events
#events <- c(780, 1300, 10400, 120000) #        Aka
#events <- c(365, 526, 530, 601, 692, 733, 764, 794, 824, 840, 1208, 1237, 1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumat\
#ra (Mentawai Segment) from Philibosian et al 2017, Figures 19,20)
#events <- events - 364


thetas1 = seq(-1e-07, 5e-06, len=50)
thetas2 = seq(-1e-05, 5e-04, len=50)
thetas3 = seq(1,100000, len=50)
rate_equation = 'theta[1] + theta[2]*exp(-(t-tlast)/theta[3])'

negloglik_grid = array(NA, c(length(thetas1), length(thetas2), length(thetas3)))

for(k in 1:length(thetas3)){
    for(j in 1:length(thetas2)){
    	for(i in 1:length(thetas1)){
	    negloglik_grid[i,j,k] = nhp$negloglik_from_theta(
                theta=c(thetas1[i], thetas2[j], thetas3[k]),
            	observed_data=events,
            	x0 = 0,
            	rate_equation=rate_equation,
            	integration_dt = 100)
	}
    }
}

M = mesh(thetas1, thetas2, thetas3)

pdf('grid_search_3D.pdf')
#par(mfrow = c(1, 2))
#M = mesh(thetas1, thetas2, thetas3)
#slice3D(thetas1, thetas2, thetas3, colvar = negloglik_grid, xs = NULL, ys = NULL, zs = 10000)
min_ind = which(negloglik_grid == min(negloglik_grid), arr.ind=TRUE)
ci_bivariate_95 = qchisq(0.95, 3)/3 
isosurf3D(thetas1, thetas2, thetas3, colvar = negloglik_grid, level = (min(negloglik_grid)+ci_bivariate_95), col = "red", ticktype = 'detailed')
grid3d(c('x', 'y', 'z'), col='black', n=10)
#surf3D(thetas1, thetas2, thetas3, negloglik_grid, col=rainbow(255)[1:200],
#    main=paste0('Negative log likelihood for ', rate_equation, 
#        '\n with minimum (point) and 95% bivariate CI') )
# Supposing the minimum is the true min, then an approx bivariate 95% CI is the
# region with likelihoods within qchisq(0.95, 2) (where the '2' comes from the fact
# we have 2 variables). For each variable individually, you can set that to '1'.
##min_ind = which(negloglik_grid == min(negloglik_grid), arr.ind=TRUE)
##ci_bivariate_95 = qchisq(0.95, 2)/2 # 
##contour(thetas1, thetas2, negloglik_grid, col='black', add=TRUE, 
##    level=(min(negloglik_grid)+ci_bivariate_95))
##points(thetas1[min_ind[1]], thetas2[min_ind[2]], col='black', pch=19, cex=2)
dev.off()