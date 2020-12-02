source('helper2.R')
source('model.R')
source('dpmm_with_gibbs.R')

library('mvtnorm')
library('LaplacesDemon')
library('MCMCpack')
library('tidyverse')
library('MASS')
library('abind')
library('circular')
library("plotly")

#Testing
n <- 5
# Gaussian
m1 <- c(1,1.5)
S1 <- matrix(c(0.3,0.05,0.05,0.3),ncol = 2)
clus1 <- mvrnorm(n=n,mu=m1,Sigma = S1)

# Gaussian
m2 <- c(1,-3)
S2 <- matrix(c(0.2,-0.08,-0.08,0.2),ncol = 2)
clus2 <- mvrnorm(n=n,mu=m2,Sigma = S2)

#Gaussian
m3 <- c(-3,2)
S3 <- matrix(c(0.1,0,0,0.1),ncol = 2)
clus3 <- mvrnorm(n=n,mu=m3,Sigma = S3)

#Line
clus4x <- c(-1.5,-1.8,-1.9,-2,-2.5)
clus4y <- -2+rnorm(1,0,1)
clus4 <- matrix(c(clus4x,clus4y),ncol=2)

m5 <- c(0,0)
S5 <- matrix(c(1,0.05,0.05,1),ncol = 2)
clus5 <- mvrnorm(n=15,mu=m1,Sigma = S1)
plot(clus5)
test_data = data.frame(x=clus5[,1],y=clus5[,2])

test_data = rbind(clus1,clus2,clus3)
plot(test_data)
test_data = data.frame(x=test_data[,1],y=test_data[,2])

alpha = 1
priors = list('mu_0' = c(0,0),
              # 'sigma_0' = diag(2)*3^2,
              # 'line_min_count' = 4,
              #measurement error
              'dot.cov' = diag(2),
              'lambda' = 0.5,
              'nu' = 2,
              'S' = diag(2),
              'crpalpha' = alpha)

c.init = sample(1,15,replace = T)
results = CRP.gibbs(test_data,c.init, priors, max.iter=1000)
a=assignment(results[['assignment']])
test_data_a = cbind(test_data,a)
test_data_a %>%
  ggplot(aes(x = x, y = y, color=LETTERS[a]))+
  geom_point()+
  theme_minimal()

Gaussian_param <- mean_cov_estimate(test_data,a,priors)


