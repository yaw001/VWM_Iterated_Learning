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
getPPD.param <- function(posterior.mu, posterior.S, posterior.lambda, posterior.nu){
  PPD.param = list()
  PPD.param[['mu']] = posterior.mu
  PPD.param[['cov']] = (posterior.lambda + 1)/(posterior.lambda * (posterior.nu - 2 + 1)) * posterior.S
  PPD.param[['df']] = posterior.nu - 2 + 1
  return(PPD.param)
}

#ll for existed cluster
logl_gaussian_cluster <- function(dots,newdot,priors) {
  n_k = nrow(dots)
  posterior.param <-  ConditionalPosterior(dots, priors)
  PPD.param <-  getPPD.param(posterior.param[['mu']],
                             posterior.param[['S']],
                             posterior.param[['lambda']], 
                             posterior.param[['nu']])
  #Assuming groups are MVN clusters, PPD is MVT
  #the probability of assigning to occupied table is proportional to the 
  #number of people are already in the table
  return(log(n_k) + dmvt(as.matrix(newdot), 
                                PPD.param[['mu']],
                                PPD.param[['cov']], 
                                PPD.param[['df']], 
                                log = T))
}
# #ll for new cluster
# logl_new_cluster <- function(dot, priors) {
#   PPD.param <-  getPPD.param(priors[['mu_0']],
#                              priors[['S']],
#                              priors[['lambda']],
#                              priors[['nu']])
#   return(log(priors[['crpalpha']]) + dmvt(as.matrix(dot),
#        PPD.param[['mu']],
#        PPD.param[['cov']],
#        PPD.param[['df']],
#        log = T))
#   # return(log(priors[['crpalpha']]) + dnorminvwishart(mu = as.matrix(dot),
#   #                                                    mu0 =priors[['mu_0']],
#   #                                                    lambda = priors[['lambda']],
#   #                                                    Sigma = priors[['sigma_0']] + priors[['dot.cov']],
#   #                                                    priors[['S']],
#   #                                                    priors[['nu']], log = T))
# }

logl_new_cluster <- function(dot, priors) {
  return(log(priors[['crpalpha']]) + dmvnorm(as.matrix(dot),
                                             priors[['mu_0']],
                                             priors[['S']],
                                             log = T))
}

cacheFilename = function(seed, chain, iter, idstring){
  return(sprintf('%s.result.seed-%d.chain-%d.iter-%d.RData', idstring, seed, chain, iter))
}

get.CRP.results<-function(chaindata, priors, max.iter){
  result = list()
  for(i in 1:20){  
    dots=filter(chaindata,Iter==i) %>% select(x=x,y=y)
    alpha = 0.5
    priors = list('mu_0' = c(0,0),
                  'lambda' = 0.1,
                  'nu' = 5,
                  'S' = diag(diag(var(dots))),
                  'crpalpha' = alpha)
    c.init = sample(1,15,replace = T)
    result.CRP<-CRP.gibbs(dots, c.init, priors, max.iter)
    a = assignment(result.CRP[['assignment']])
    result[['priors']] = priors
    result[['assignment']] = a
    result[['mean']] = mean_cov_estimate(dots,a,priors)[["mean"]]
    result[['cov']] = mean_cov_estimate(dots,a,priors)[["cov"]]
    filename<-cacheFilename(chaindata$Seed[1], 
                            chaindata$Chain[1], 
                            i, 
                            idstring='Inferred_clusters')
    save(result,file=filename)
  }
}



#ll for line
logl_line_cluster <- function(dots, n_k, newdot, line_param) {
  dot_mean <-  unname(colMeans(dots))
  dot_cov <- unname(cov(dots))
  dot_pc_vectors = eigen(unname(cov(dots)))$vectors
  pc1_proj = as.matrix(dots) %*% dot_pc_vectors[,1]
  newdot_pc1 = as.matrix(newdot) %*% dot_pc_vectors[,1]
  # newdot_pc1 = t(as.matrix(newdot)) %*% dot_pc_vectors[,1]
  # dot_pc1_radians <- acos(dot_pc_vectors[,1][1])
  #angle <- rvonmises(1,dot_pc1_radians, line_param[['kappa']])
  len <- diff(range(pc1_proj))
  deviation <- as.matrix(newdot) %*% dot_pc_vectors[,2]
  return(log(n_k) +
           log(dnorm(newdot_pc1,mean(pc1_proj),sd(pc1_proj))*0.5 +
              ifelse(newdot_pc1 < max(pc1_proj) & newdot_pc1 > min(pc1_proj),
                  1/len,
                  0)*0.5)+
           log(1/len) +
           dnorm(deviation, 0, line_param[['orth_noise']], log=T))
}

#
#ll for line
logl_line_cluster <- function(dots, n_k, newdot, line_param) {
  dot_mean <-  unname(colMeans(dots))
  dot_cov <- unname(cov(dots))
  dot_pc_vectors = eigen(unname(cov(dots)))$vectors
  pc1_proj = as.matrix(dots) %*% dot_pc_vectors[,1]
  newdot_pc1 = as.matrix(newdot) %*% dot_pc_vectors[,1]
  # newdot_pc1 = t(as.matrix(newdot)) %*% dot_pc_vectors[,1]
  # dot_pc1_radians <- acos(dot_pc_vectors[,1][1])
  #angle <- rvonmises(1,dot_pc1_radians, line_param[['kappa']])
  len <- diff(range(pc1_proj))
  deviation <- as.matrix(newdot) %*% dot_pc_vectors[,2]
  return(log(n_k) +
           log(dnorm(newdot_pc1,mean(pc1_proj),sd(pc1_proj))*0.5 +
                 ifelse(newdot_pc1 < max(pc1_proj) & newdot_pc1 > min(pc1_proj),
                        1/len,
                        0)*0.5)+
           log(1/len) +
           dnorm(deviation, 0, line_param[['orth_noise']], log=T))
}

logl_line_gaussian_mixture <- function(dots, newdot, priors){
  # estimate posterior p(line | dots)
  # estimate 
  # p(newdot | line, dots) p(line | dots) + p(newdot | gaussian, dots) p(gaussian | dots)
  # x p(newdot | line, dots) 
  # x p(newdot | gaussian, dots)
  # p(line | dots) + p(gaussian | dots) = 1
  # p(line | dots) = p(dots | line) p(line) / 
  #             (p(dots | line) p(line) + (p(dots | gaussian) (1-p(line)))
  # p(dots | line) = sum_{center, length, angle} p(dots | center, length, angle) * p(center) p(length) p(angle)
  # p(dots | gaussian)  dmvt(priors)
  # p(dots | line) = p(dots | ML center, ML length, ML angle) dnorm(ML center) 1/length 1/2pi
}

# eigendecomposition
# dt(pc1) * dt(pc2)



