# input
#   dots[n x d ] data frame of dot x,y coordinates (and other junk?)
#   c.init [n] vector initial assignment.
#   priors list of prior parameters
#   line_param [priors about lines]
#   max.iter (integer) number of gibbs samples?

CRP.gibbs <- function(dots,c.init, priors, line_param, max.iter) {
  num.dots <- nrow(dots)
  # DELETE:  K = 2
  # DELETE: groups = list()
  #initial membership assignment
  n_k <- as.vector(table(c.init))
  #initial number of clusters
  Nclust <- length(n_k)
  #cluster membership storage [samplse x dots]
  z <- array(NA,dim=c(max.iter, num.dots))
  # z_category <- array(NA,dim=c(max.iter, num.dots))
  z[1,] = c.init
  pb <-  txtProgressBar(min = 0, max = max.iter, style = 3)
  #maxIters prevents endless loops
  for(iter in 1:max.iter) {
    if(iter > 1){
      #The assignment of next iteration initially based on the previous one
      z[iter,] = z[iter-1,]
      n_k <- as.vector(table(z[iter-1,]))
    }
    #one data point at a time
    for(n in 1:num.dots) {
      #when nth data point (customer) enters the CR, we need to first un-assign its initial
      #cluster membership, then use the PPD for occupied table and Prior PD for new table, 
      #update the assignment
      #c_i: the membership assignment of nth person (Data point)
      c_i <- z[iter,n]
      #remove the nth person from the table
      n_k[c_i] <- n_k[c_i] - 1
      #if that cluster becomes empty, then remove the cluster
      if(n_k[c_i] == 0){
        #move last cluster to replace the empty table
        n_k[c_i] <- n_k[Nclust]
        #move the people from the last table to fill the emptied table
        z[iter,(z[iter,] == Nclust)] <- c_i
        #take out the last cluster (now empty)
        n_k <- n_k[-Nclust]
        #drop the total number of clusters by 1
        Nclust <- Nclust - 1
      }
      # unassign the removed dot
      z[iter,n] <- -1
      #update table assignment
      #log probabilities for the clusters, add previously unoccupied table
      logp <- rep(NA, Nclust + 1)
      #loop over already occupied tables 1:J for ll
      for(j in 1:Nclust) {
        if(n_k[j] < 4) {  # TODO: move line threshold number to priors priors[['line_min_count']]
          logp[j] <- logl_existed_cluster(dots[z[iter,]==j,], n_k[j], dots[n,],priors)  # remove n_k since that is redundant with nrow of first argument
        } else {
          # TODO: rename logl_existr_cluster to logl_gaussian_cluster and logl_line to logl_line_cluster
          # TODO: more renaming
          logl_cluster = logl_existed_cluster(dots[z[iter,]==j,],n_k[j],dots[n,],priors)
          logl_line = logl_line(dots[z[iter,]==j,], n_k[j], dots[n,], line_param)
          # TODO: this seems wrong....
          w = logl_cluster/(logl_line+logl_cluster)
          logp[j] <- logl_cluster*w + logl_line*(1-w)
        }
      }
      #ll of unocuupied ("new") table, using prior predictive distribution
      logp[Nclust + 1] <- logl_new_cluster(dots[n,], priors)
      #normalization for the log probabilities
      max_logp <- max(logp)
      #for numerical stability purpose
      logp <- logp - max_logp
      loc_probs <- exp(logp)
      #apply softmax
      loc_probs <- loc_probs / sum(loc_probs)
      #draw a sample of which cluster this new person should belong to
      newz <- sample(1:(Nclust+1), 1, replace = T, prob = loc_probs)
      #spawn a new cluster if neccessary
      if(newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1
      }
      z[iter,n] <- newz
      n_k[newz] <- n_k[newz] + 1
    }
    #update txt progress bar after each iter
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  return(z)
}

# tabulates the counts and assign the dot to most likely cluster membership assignment 
# 
# input: assignment_sample_matrix[s(samples) x n(items/dots)] indicating the cluster-assignment index (integer) for each dot
# output: assign_z [n(items)] indicating the most frequent cluster assignment in the set of samples?
assignment <- function(assignment_sample_matrix) {
  n = ncol(assignment_sample_matrix)
  assign_z = rep(NA, n)
  for(i in 1:n){
    assign_z[i] = unname(which.max(table(assignment_sample_matrix[,i])))
  }
  return(assign_z)
}

#cluster mean and cov estimate
# inputs: 
#   dots[n (dots), d (dimensions)] data frame indicating the x,y coordinates of each item
#   z: assignment vector (integer)
#   priors: list of prior parameters
mean_cov_estimate <- function(dots, z, priors){
  mean = list()
  cov = list()
  len = max(z)
  for(i in 1:len) {
    n = sum(z==i)
    dots_of_i = dots[z == i,c('x','y')]
    y.bar = unname(colMeans(dots_of_i))
    posterior.lambda = priors[['lambda']] + n
    posterior.nu <- priors[['nu']] + n
    dev.mu = as.vector(y.bar - priors[['mu_0']])
    posterior.S <- priors[['S']] + (priors[['lambda']]*n / (priors[['lambda']] + n)) * (dev.mu %*% t(dev.mu))
    S.n <- unname(cov(dots_of_i))*(n-1)
    posterior.S <- posterior.S + S.n
    mean[[i]] = (priors[['lambda']]*priors[['mu_0']] + n*(y.bar))/(priors[['lambda']] + n)
    #Mean of inverse-wishart psi/ (nu - p - 1), psi (p by p)
    cov[[i]] = posterior.S / (posterior.nu - 2 - 1)
  }
  return(list(mean = mean,
              cov = cov))
}

