# Conditional posteriors
ConditionalPosterior <- function(dots, priors) {
  posterior = priors
  n <-nrow(dots)
  y.bar = unname(colMeans(dots[,c('x','y')]))
  posterior[['mu']] <- (priors[['lambda']]*priors[['mu_0']] + n*(y.bar))/(priors[['lambda']] + n)
  posterior[['lambda']] <- priors[['lambda']] + n
  posterior[['nu']] <- priors[['nu']] + n
  dev.mu = as.vector(y.bar - priors[['mu_0']])
  posterior[['S']] <- priors[['S']] + (priors[['lambda']]*n / (priors[['lambda']] + n)) * (dev.mu %*% t(dev.mu))
  if(n > 2) {
    S.n <- unname(cov(cbind(dots$x,dots$y)))*(n-1)
    posterior[['S']] <- posterior[['S']] + S.n
  }
  return(posterior)
}

#get the parameters of the Predictive posterior for clusters
#PPD is t distribution for MVN clusters
getPPD.param <- function(posterior.mu,posterior.S, posterior.lambda, posterior.nu, K){
  PPD.param = list()
  PPD.param[['mu']] = posterior.mu
  PPD.param[['cov']] = (posterior.lambda + 1)/(posterior.lambda * (posterior.nu - 2 + 1)) * posterior.S
  PPD.param[['df']] = posterior.nu - 2 + 1
  return(PPD.param)
}
#ll for existed cluster
logl_gaussian_cluster <- function(cluster,dot,priors) {
  n_k = nrow(dots)
  posterior.param <-  ConditionalPosterior(cluster, priors)
  PPD.param <-  getPPD.param(posterior.param[['mu']],
                             posterior.param[['S']],
                             posterior.param[['lambda']], 
                             posterior.param[['nu']], 
                             K)
  #Assuming groups are MVN clusters, PPD is MVT
  #the probability of assigning to occupied table is proportional to the 
  #number of people are already in the table
  return(log(n_k) + dmvt(as.matrix(dot), 
                                PPD.param[['mu']],
                                PPD.param[['cov']], 
                                PPD.param[['df']], 
                                log = T))
}
#ll for new cluster
logl_new_cluster <- function(dot, priors) {
  return(log(priors[['crpalpha']]) + dnorminvwishart(mu = as.matrix(dot),
                                                                  mu0 =priors[['mu_0']],
                                                                  lambda = priors[['lambda']],
                                                                  Sigma = priors[['sigma_0']] + priors[['dot.cov']],
                                                                  priors[['S']],
                                                                  priors[['nu']], log = T))
}

#get center of mass
get_cm <- function(dots){
  return(unname(colMeans(dots[,c('x','y')])))
}

#ll for line
logl_line <- function(dots, dot, line_param) {
  n_k = nrow(dots)
  cm <-  get_cm(dots)
  eigens <- eigen(cov(dots))
  proj.x <- as.matrix(dots)%*%eigens$vector[,1]
  proj.y <- as.matrix(dots)%*%eigens$vector[,2]
  return(log(n_k) + 
           log(1/(2*eigens$values[1])) + 
           dnorm(proj.x,cm[1],sqrt(eigens$values[1]),log=T) + 
           dnorm(proj.y,cm[2],sqrt(eigens$values[2]),log=T))
}
