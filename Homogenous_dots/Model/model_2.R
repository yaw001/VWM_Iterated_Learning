models = list()

#concentration parameter
alpha = 0.081
# data dimension 
K = 2
scale = 1

#Positive scalar for IW
#encoding noise (measurement)
dot.sd = 30
dot.cov = diag(2)*dot.sd^2


#prior mean
mu_0 = c(0,0)

#cluster
models[['cluster']] = list('basename' = "cluster",
                             #prior mean and cov for the unknown mean: mu_0, sigma_0
                             #priors for the unknown covariance matrix:
                             # (1) scalar degress of freddom: nu 
                             # (2) symmetric psd scale matrix S (uncertainty on the variance)
                             # (3) positive scalar: lambda
                             priors = list('mu_0' = c(0,0),
                                           'lambda' = 0.5,
                                           'nu' = K + scale,
                                           'S' = diag(K)*scale,
                                           'crpalpha' = alpha))

# 
# scale = 10
# S = diag(K)*0.5^2*scale
# nu = K + scale
# 
# scale = 0.001
# S = diag(K)*scale
# nu = K + scale

#line
models[['line']] = list('basename' = 'line',
                        param = list('angle_noise' = 0.59,
                                       'pos_noise' = diag(2)*55^2,
                                       'orth_noise' = 2.5,
                                       'len_noise' = 2,
                                       'kappa' = 1))



