nhp=new.env()
source('../stormwavecluster/R/nhpoisp/nhpoisp.R', local=nhp)
nhp

# Define rate equation
rate_equation = 'theta[1] + theta[2]*exp(-(t-tlast)/theta[3])'
#rate_equation = 'theta[1] + theta[2]*exp(-(t-tlast)/14486)'

#events <- c(4000000.0,3500000.0,3496429.0,3492858.0,3489287.0,3485716.0,3482145.0,3478574.0,3475003.0,3471432.0,3467861.0,3464290.0,3460719.0,3457148.0,3453577.0,3450006.0,3446435.0,3442864.0,3439293.0,3435722.0,3432151.0,3428580.0,3425009.0,3421438.0,3417867.0,3414296.0,3410725.0,3407154.0,3403583.0,3400012.0,1680000.0,1673750.0,1667500.0,1661250.0,1655000.0,1648750.0,1642500.0,1636250.0,1630000.0,1623750.0,1617500.0,1611250.0,1605000.0,1598750.0,1592500.0,1586250.0,1070000.0,1066800.0,1063600.0,1060400.0,1057200.0,1054000.0,1050800.0,1047600.0,1044400.0,1041200.0,1038000.0,1034800.0,1031600.0,1028400.0,1025200.0,1022000.0,1018800.0,1015600.0,1012400.0,1009200.0,1006000.0,1002800.0,999600.0,996400.0,993200.0) # Lake George
events <- c(993200.0,996400.0,999600.0,1002800.0,1006000.0,1009200.0,1012400.0,1015600.0,1018800.0,1022000.0,1025200.0,1028400.0,1031600.0,1034800.0,1038000.0,1041200.0,1044400.0,1047600.0,1050800.0,1054000.0,1057200.0,1060400.0,1063600.0,1066800.0,1070000.0,1586250.0,1592500.0,1598750.0,1605000.0,1611250.0,1617500.0,1623750.0,1630000.0,1636250.0,1642500.0,1648750.0,1655000.0,1661250.0,1667500.0,1673750.0,1680000.0,3400012.0,3403583.0,3407154.0,3410725.0,3414296.0,3417867.0,3421438.0,3425009.0,3428580.0,3432151.0,3435722.0,3439293.0,3442864.0,3446435.0,3450006.0,3453577.0,3457148.0,3460719.0,3464290.0,3467861.0,3471432.0,3475003.0,3478574.0,3482145.0,3485716.0,3489287.0,3492858.0,3496429.0,3500000.0,4000000.0) # Lake george
#events <- c(30000, 50000, 90000, 1000000)
#events <- c(30000, 50000, 90000, 250000, 500000, 530000, 570000, 590000, 850000)
#events <- c(32000, 38500, 45000, 55000, 62500, 70000, 1000000, 2000000, 4500000) # Cadell
#events <- c(15100, 18000, 21000, 24000, 55000) # D'n 4 events
#events <- c(15100, 17000,19000, 21000, 22500, 24000, 55000) # D'n 6 events
#events <- c(780, 1300, 10400, 120000) #	Aka
#events <- c(9000, 18000, 20000, 90000, 150000) #Wharekuri estimated
#events <- c(10000,18000,26000,34000,50000,58000,1000000) #Waitangi est
events <- c(440, 1100, 1800, 2400, 99000) # Rurrand (Roer Valley)
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
    initial_theta=c(1/1e6, 1/10000, 10000.),
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
    optim_control=list(maxit = 500, parscale=c(1e-6, 1e-4, 1e3)),
    verbose=TRUE)

print(model_fit$rate_equation)
print(model_fit$par)
present_rate <- model_fit$par[1] + model_fit$par[2]*exp(-(12100)/model_fit$par[3])
#present_rate <- model_fit$par[1] + model_fit$par[2]*exp(-(20000)/14486)
print(present_rate)
# Calculate standard Poisson rate (i.e. mean event rate)
empirical_pois_rate = length(events)/events[length(events)]
print('Empirical Poisson rate')
print(empirical_pois_rate)
nhp$get_fit_standard_errors(model_fit)