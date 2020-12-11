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

setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Data')
load('all.data.Rdata')
all.data
all.data$x = all.data$x/300
all.data$y = all.data$y/300

datc = all.data %>% filter(Seed == 2, Chain == 2, Iter==20) %>% select(x=x,y=y) %>% as.matrix()
datc = all.data %>% filter(Seed == 1, Chain == 3, Iter==17) %>% select(x=x,y=y) %>% as.matrix()

plot(datc)

#Prior specification
alpha = 0.08
mu0 <- matrix(rep(0,2),ncol=2,byrow =TRUE)
lambda=0.02
# lambda=0.02
# lambda = 0.05
nu=10
# psi=diag(diag(var(datc)))
psi=diag(2)*mean(diag(var(datc)))/2
# psi=var(datc)/2
# psi=diag(2)*0.0625

c_init <- rep(1,nrow(datc))
results <- crp_gibbs(datc,alpha=alpha,mu0=mu0, lambda=lambda, nu=nu, psi=psi, c_init=c_init, maxIter=2000)
c_model <- apply(results[,-(1:200)], 1, FUN = function(x) { 
  tab <- table(x)
  ans <- names(tab[which.max(tab)]) 
  return(ans) })
# c_true <- rep(1:4, each = 5)
# table(c_true, c_model)
datc_model = data.frame(x=datc[,1],y=datc[,2],assignment = c_model)
datc_model %>% ggplot(aes(x,y,color=assignment))+geom_point()


# lambda=c(0.01,0.05,0.1,0.5,1)
# nu=c(5,10,15)
# psi=diag(diag(var(datc))/2)
# assignment = c()
# for(i in 1:length(lambda)){
#   for(j in 1:length(nu)){
#     results = crp_gibbs(datc,alpha=alpha,mu0=mu0, lambda=lambda[i], nu=nu[j], psi=psi, c_init=c_init, maxIter=2000)
#     c_model <- apply(results, 1, FUN = function(x) { 
#       tab <- table(x)
#       ans <- names(tab[which.max(tab)]) 
#       return(ans) })
#     cat("lambda = ", lambda[i], "nu =", nu[j], "\n")
#     print(table(c_model))
#     assignment = c(assignment,c_model)
#   }
# }

datc = all.data %>% filter(Seed==1,Chain==3)
assignment = c()
for(iter in 1:20) {
  dat_iter=datc %>% filter(Iter==iter) %>% select(x,y) %>% as.matrix()
  psi=diag(2)*mean(diag(var(dat_iter)))/2
  c_init <- rep(1,nrow(dat_iter))
  results <- crp_gibbs(dat_iter,alpha=alpha,mu0=mu0, lambda=lambda, nu=nu, psi=psi, c_init=c_init, maxIter=2000)
  c_model_iter = apply(results[,-(1:200)], 1, FUN = function(x) { 
    tab <- table(x)
    ans <- names(tab[which.max(tab)]) 
    return(ans) })
  assignment = c(assignment,c_model_iter)
}

datc_model = cbind(datc,assignment=assignment)
datc_model$Iter = as.factor(datc_model$Iter)
datc_model$assignment = as.integer(datc_model$assignment)
# datc_model$assignment = as.factor(datc_model$assignment

setwd("/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Results/Test")
save(datc_model, file = "test_result_seed_1_chain_3.Rdata")
library("plotly")
# cols <- c("1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange","5" = "purple")
p=datc_model %>% ggplot(aes(x,y,color=assignment,frame=Iter))+geom_point()+
  scale_color_gradient(low="blue", high="red")
#   scale_color_manual(values=cols)
ggplotly(p)

datc_model %>%
  plot_ly(
    x = ~x, 
    y = ~y, 
    color = ~assignment, 
    frame = ~as.factor(Iter), 
    type = 'scatter',
    mode = 'markers'
  )


#seed 1 chain 3
assignment_1 = assignment
assignment_2 = assignment
assignment_3 = assignment
assignment_4 = assignment
# var(as.matrix(dist(datc))[lower.tri(as.matrix(dist(datc)))])
# 
# var_pair_dist <- function(){
#   # x=runif(15,-1,1)
#   # y=runif(15,-1,1)
#   x=runif(15,min(all.data$x),max(all.data$x))
#   y=runif(15,min(all.data$y),max(all.data$y))
#   dat_random = data.frame(x=x,y=y)
#   plot(dat_random)
#   return(var(as.matrix(dist(dat_random))[lower.tri(as.matrix(dist(dat_random)))]))
# }
# a=var_pair_dist()
# a
# ran_pair_var=replicate(10000,var_pair_dist())
# mean(ran_pair_var)
# sd(ran_pair_var)
# hist(ran_pair_var)
# datc = all.data %>% filter(Seed == 1, Chain == 3, Iter==1) %>% select(x=x,y=y) %>% as.matrix()
# var(as.matrix(dist(datc))[lower.tri(as.matrix(dist(datc)))])

plot(datc)

mean.det_all <- all.data%>%group_by(Seed,Chain,Iter)%>%
  summarise(det.cov=det(cov(cbind(x,y))))%>%
  group_by(Iter)%>%
  summarise(N=n(),mean.det=mean(det.cov),se.mean=sd(det.cov)/sqrt(N))

mean.det_all%>%ggplot(aes(x=Iter,y=mean.det))+
  geom_point()+
  geom_line(size=1.1)+
  geom_errorbar(aes(ymin=mean.det-se.mean,ymax=mean.det+se.mean),width=0.5,size=1.1)+
  labs(x = 'Iteration',
       y = 'Covariance Determinant (Pattern spread)')+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 16, face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = 'none',
        axis.title = element_blank())
