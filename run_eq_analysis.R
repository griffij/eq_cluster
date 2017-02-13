nhp=new.env()
source('../stormwavecluster/R/nhpoisp/nhpoisp.R', local=nhp)
nhp

# Read the data
#datafile <- '../Data/eq_data.csv'
#event_data <- read.table(file = datafile, sep =',')
#event_times <- subset(event_data, select = -c(V1))
#print(event_data)
#print(event_times)

# Calculate intervent times
#interevent_times <- c()
#ncols = ncol(event_times)
#for (i in seq(2, ncols)) {
#    interevent_times <-c(interevent_times,c(event_times[i] - event_times[(i-1)]))
#    print('here')
#    print(interevent_times)
#    }
## Put all interevent times into a single vector
#interevent_times <- c(interevent_times, recursive=TRUE)
#interevent_times <- unname(interevent_times)
#interevent_times <- sort(interevent_times, decreasing=FALSE)
#print(interevent_times)

# Define rate equation
rate_equation = 'theta[1] + theta[2]*exp(-(t-tlast)/theta[3])'
#rate_equation = 'theta[1] + theta[2]*exp(-(t-tlast)/14486)'
#print('event_times[1,]')
#print(event_times[1,])
#events = unname(event_times[1,])
#print(events)
#events <- c(30000, 50000, 90000, 1000000)
#events <- c(30000, 50000, 90000, 250000, 500000, 530000, 570000, 590000, 850000)
#events <- c(32000, 38500, 45000, 55000, 62500, 70000, 1000000, 2000000, 4500000) # Cadell
#events <- c(365, 526, 530, 601, 692, 733, 764, 794, 824, 840, 1208, 1237, 1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumatra (Mentawai Segment) from Philibosian et al 2017, Figures 19,20)
#events <- c(1314, 1350, 1388, 1569, 1597, 1613, 1631, 1658, 1703, 1797, 1833, 2007, 2010) #Sumatra (Mentawai Segment) from Philibosian et al 2017, Figures 20)
#events <- events - 1313
# Re-order to make time run the right way
events <- rev(events)
events <- (events - events[1])*(-1)
print(events)
# Call nhpoisp optimisation
model_fit = nhp$fit_nhpoisp(
    events,
    rate_equation=rate_equation,
    minimum_rate=0.,
    initial_theta=c(1/2e6, 1/50000, 30000.),
    ##
    ## The arguments below control details of the optimization
    ## It can be difficult to fit these models with complex rate functions,
    ## so adjustment may be required.
    ## Consult the code documentation and see help on R's "optim"
    ## function for details.
    ##
    integration_dt = 1.0e2,
    number_of_passes=1,
    enforce_nonnegative_theta=TRUE,
    optim_method=c('Nelder-Mead'),
    optim_control=list(maxit = 500, parscale=c(1e-5, 1e-4, 1e4)),
    verbose=TRUE)

print(model_fit$rate_equation)
print(model_fit$par)
present_rate <- model_fit$par[1] + model_fit$par[2]*exp(-(32000)/model_fit$par[3])
#present_rate <- model_fit$par[1] + model_fit$par[2]*exp(-(20000)/14486)
print(present_rate)
nhp$get_fit_standard_errors(model_fit)