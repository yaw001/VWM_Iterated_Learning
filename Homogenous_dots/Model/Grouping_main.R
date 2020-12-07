library('mvtnorm')
library('LaplacesDemon')
library('MCMCpack')
library('tidyverse')
library('MASS')
library('abind')
library('circular')
library("plotly")

setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Data')
load('all.data.Rdata')
all.data

setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Model')

source('helper2.R')
source('model_2.R')
source('dpmm_with_gibbs.R')

setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Results/Inferred_cluster')


for(seeds in 1){
  for(chains in 1:10){
    chaindata <- filter(all.data, Seed == seeds, Chain == chains)
    get.CRP.results(chaindata,
                    max.iter=1000)
  }
}  




# assignment.iter = c()

tb.results.group = tibble()
for(seeds in 1){
  for(chains in 1:3){
    for(iter in 1:20){
      group.filename<-cacheFilename(seeds, 
                                    chains, 
                                    iter, 
                                    idstring='Inferred_clusters')
      load(group.filename)
      inferred_group = tibble(prior=list(result$prior),assignment = list(result$assignment), group_mean = list(result$mean),group_cov = list(result$cov))
      tb.results.group = bind_rows(tb.results.group,inferred_group)
      # assignment.iter = c(assignment.iter,result[['assignment']])
    }
  }
}

tb.result_seed_1_3<-all.data%>% filter(Seed==1,Chain==3) %>% 
  group_by(Seed,Chain,Iter)%>% 
  do(tibble(xy = list(data.frame(x=.$x, y=.$y)))) %>% bind_cols(tb.results.group) %>% 
  rowwise() %>% 
  mutate(group_numbers = length(unique(unlist(assignment))))

tb.result_seed_1_to_3 <- tb.result_seed_1_to_3 %>% bind_cols(tb.results.group)
tb.result_seed_1_to_3 %>% filter(Iter == 20)

assignment_1_3 = tb.result_seed_1_3 %>% pull(assignment) %>% unlist()
Seed_1_chain_3 = all.data %>% filter(Seed == 1, Chain == 3) %>% cbind(assignment = as.factor(unlist(assignment_1_3)))

p <- Seed_1_chain_3 %>% group_by(Iter) %>% 
  ggplot(aes(x=x, y=y, color=assignment),frame=Iter)+
  # facet_grid(Iter~)+
  # scale_color_gradient(low="blue", high="red")+
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

Seed_1_chain_3_iter_20 = all.data %>% filter(Seed == 1, Chain == 3, Iter==20) %>% select(x=x,y=y)
p <- all.data %>% filter(Seed == 1, Chain == 3) %>% group_by(Iter) %>% 
  ggplot(aes(x=x, y=y, frame=Iter))+
  # facet_grid(Chain~Seed)+
  # scale_color_gradient(low="blue", high="red")+
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
