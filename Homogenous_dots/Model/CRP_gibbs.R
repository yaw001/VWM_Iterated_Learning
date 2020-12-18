#A tutorial on DPMM
library('MASS')
library('tidyverse')
#CRP_gibbs

crp_gibbs <- function(data, alpha, mu0, lambda, nu, psi, c_init, maxIter) {
  require(mvtnorm)
  data_dim <- ncol(data)
  N <- nrow(data)
  #initialize the CRP gibbs sampler
  z <- c_init
  #initial data counts at each cluster
  n_k <- as.vector(table(z))
  #initial number of clusters
  Nclust <- length(n_k)
  #CRP gibbs sampler
  #cluster membership storage
  res <- matrix(NA,nrow=N,ncol=maxIter)
  pb <- txtProgressBar(min=0,max=maxIter,style=3)
  for(iter in 1:maxIter){
    for(n in 1:N){
      #One data point at a time; when nth customer enters the Chinese restaurant, we need to first un-assign his/her intial cluster membership,
      #then calculate the updated probability of table assignment
      #The nth person's table assignment
      c_i <- z[n]
      #remove the nth person from the table
      n_k[c_i] <- n_k[c_i] - 1 
      #If the table becomes empty when the nth person is removed, then that table is removed
      if(n_k[c_i] == 0){
        #last cluster to replace this empty cluster
        n_k[c_i] <- n_k[Nclust]
        #who are in the last cluster?
        loc_z <- (z == Nclust)
        #move them up to fill just emptied cluster
        z[loc_z] <- c_i
        #Take out the last cluster
        n_k <- n_k[-Nclust]
        #Decrease the total number of clusters by 1
        Nclust <- Nclust -1
      }
      #ensures z[n] doesnt get counted as a cluster
      z[n] <- -1
      #Log probabilities for the clusters
      logp <- rep(NA, Nclust + 1)
      #loop over already occupied tables 1:J and calculate probability
      for(c_i in 1:Nclust) {
        #posterior lambda
        lambda_p = lambda + n_k[c_i]
        #posterior nu
        nu_p = nu + n_k[c_i]
        #find all of the points in this cluster
        loc_z <- which(z == c_i)
        #sum all the points in this cluster
        if(length(loc_z)>1){
          sum_data <- colSums(data[z==c_i, ])
        } else{
          sum_data <- data[z == c_i,]
        }
        #posterior cluster mu
        mean_p = (lambda*mu0 + sum_data)/lambda_p
        #posterior covariance sigma_p
        #sigma_p = psi + variance/covairance of cluster + deviation from prior mean
        #variance/covariance
        if(length(loc_z)>1){
          # the sum of squared deviations from the cluster mean,
          # calculated by first taking the deviations (subtracting the mean from y1 and y2 with the sweep() function) 
          # so that they can be squared and summed using the crossprod() function.
          mean_data = colMeans(data[z==c_i, ])
          y_ss <- sweep(data[z==c_i, ], MARGIN = 2, STATS = mean_data, FUN = "-")
          S = crossprod(y_ss)
          # S = cov(data[z==c_i, ])*(length(loc_z))
        }else{
          S = 0
        }
        #deviation from prior mean
        dev.mu0 = as.vector(mean_data - mu0)
        #posterior covariance amtrix sigma_p
        #lambda*n/lambda_p is a term captures the precision of the prior (lambda less than 1 (small), prior on cluster mean is diffused (play less role)).
        #when lambda greater than 1, when data mean differs more from the prior on cluster mean, the larger uncertainty is incorporated in the posterior sigma.
        #when it is too large, significantly increase the sigma_p => single cluster
        sigma_p = psi + S + (lambda*n/lambda_p) * (dev.mu0%*%t(dev.mu0))
        #Predicitve distribution for already occupied tables to predict the next customer sitting there
        #Predictive distribution is multivariate t distribution with df = nu_p - p + 1, p is dimension of the data
        #predictive distribution mean = mean_p
        #predictive distribution covariance = (lambda_p + 1)/(lambda_p*(nu_p - 2 + 1)) *sigma_p
        #likelihood of the assignment is proportion to the number of customers occupied in the table
        #(lambda_p + 1)/(lambda_p*(nu_p - 2 + 1)) = 1/(nu_p - 2 + 1) + 1/(lambda_p*(nu_p - 2 + 1)),lambda_p>1, large lambda => nu_p primarily scales the sigma_p,
        #small lambda => lambda_p close to n => larger n (more samples) implies more informative/accurate cluster estimate (smaller sigma_p), supporting more and little clusters
        logp[c_i] <- log(n_k[c_i]) + dmvt(data[n,], mean_p, sigma = (lambda_p + 1)/(lambda_p*(nu_p - 2 + 1)) * sigma_p, df = nu_p - 2 + 1, log = T)
      }
      logp[Nclust + 1] <- log(alpha) + dmvt(data[n,], delta = mu0, sigma = psi, df = nu - 2 + 1, log = T)
      #transform unnormalized log probabilities into probabilities
      max_logp <- max(logp)
      logp <- logp - max_logp
      loc_probs <- exp(logp)
      loc_probs <- loc_probs/sum(loc_probs)
      #draw a sample of which cluster this new customer should belong to
      newz <- sample(1:(Nclust+1),1,replace = T, prob = loc_probs)
      #spawn a new cluster if necessary
      if(newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1
      }
      z[n] <- newz
      n_k[newz] <- n_k[newz] + 1
    }
    setTxtProgressBar(pb,iter)
    #cluster membership of N observations
    res[,iter] <- z
  }
  close(pb)
  #return results, n by maxIter matrix
  invisible(res)
}

#Estimate Cluster mean and covariance
cluster_param <- function(dots,assignment,alpha,lambda,mu0,nu,psi){
  mean = list()
  cov = list()
  len = max(assignment)
  for(i in 1:len) {
    n = sum(assignment==i)
    dots_of_i = dots[assignment == i,c('x','y')]
    y.bar = unname(colMeans(dots_of_i))
    posterior.lambda = lambda + n
    posterior.nu <- nu + n
    dev.mu = as.vector(y.bar - mu0)
    posterior.S <- psi + (lambda*n / (lambda + n)) * (dev.mu %*% t(dev.mu))
    S.n <- unname(cov(dots_of_i))*(n)
    posterior.S <- posterior.S + S.n
    mean[[i]] = (lambda*mu0+ n*(y.bar))/(lambda + n)
    #Mean of inverse-wishart psi/ (nu - p - 1), psi (p by p)
    cov[[i]] = posterior.S / (posterior.nu - 2 - 1)
  }
  return(list(mean = mean,
              cov = cov))
}




