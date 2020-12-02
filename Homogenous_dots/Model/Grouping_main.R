rm(list=ls())

setwd('/Users/young/Desktop/UCSD/Research/Color-Iterated-Learning/Paper submission/Model')

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


all.data <- read.csv('dots_perc.csv')
names(all.data)
all.data <- all.data %>% select(Seed,Chain,Iter,dots,x,y)
save(all.data,file = 'all.data.Rdata')

load('all.data.Rdata')
all.data



setwd('/Users/young/Desktop/UCSD/Research/Color-Iterated-Learning/Paper submission/Inferred_groups')
for(seeds in 1){
  for(chains in 3){
    chaindata <- filter(data, Seed == seeds, Chain == chains)
    get.CRP.results(chaindata,
                    priors,
                    max.iter=1000)
  }
}  

alpha = 0.5
priors = list('mu_0' = c(0,0),
              'lambda' = 0.1,
              'nu' = 5,
              'S' = diag(2)*20000,
              'crpalpha' = alpha)

assignment.iter = c()
for(seeds in 1){
  for(chains in 3){
    for(iter in 1:20){
      group.filename<-cacheFilename(seeds, 
                                    chains, 
                                    iter, 
                                    idstring='Inferred_clusters')
      load(group.filename)
      assignment.iter = c(assignment.iter,result[['assignment']])
    }
  }
}
Seed_1_chain_3 = all.data %>% filter(Seed == 1, Chain == 3) %>% cbind(assignment = as.factor(assignment.iter))

Seed_1_chain_3_iter_20 = all.data %>% filter(Seed == 1, Chain == 3, Iter==20) %>% select(x=x,y=y)
results = CRP.gibbs(Seed_1_chain_3_iter_20,c.init, priors, max.iter=1000)
a=assignment(results[['assignment']])

a

p <- Seed_1_chain_3 %>% group_by(Iter) %>% 
  ggplot(aes(x=x, y=y,color=assignment.iter, frame=Iter))+
  # facet_grid(Chain~Seed)+
  scale_color_gradient(low="blue", high="red")+
  geom_point(size=1)+
  theme_bw()+
  # coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none')

ggplotly(p)
