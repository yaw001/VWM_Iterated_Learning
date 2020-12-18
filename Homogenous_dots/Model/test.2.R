library("plotly")

setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Data')
load('all.data.Rdata')
all.data
all.data$x = all.data$x/300
all.data$y = all.data$y/300

setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Model')
source('CRP_gibbs.R')
source('cluster_assess.R')
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

datc = all.data %>% filter(Seed==1,Chain==3)
assignment = c()
samples = list()
psi_list = list()
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
  samples[[iter]] = results
  psi_list[[iter]] = psi
}

datc_model = cbind(datc,assignment=assignment)
datc_model$Iter = as.factor(datc_model$Iter)
datc_model$assignment = as.integer(datc_model$assignment)
# datc_model$assignment = as.factor(datc_model$assignment

setwd("/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Homogenous_dots/Results/Test")
save(datc_model, file = "test_result_seed_1_chain_3.Rdata")
save(samples, file = "samples_seed_1_chain_3.Rdata")
save(psi_list, file = "psi_list_seed_1_chain_3.Rdata")
library("plotly")
load("test_result_seed_1_chain_3.Rdata")
load("samples_seed_1_chain_3.Rdata")
# cols <- c("1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange","5" = "purple")
p=datc_model %>% ggplot(aes(x,y,color=assignment,frame=Iter))+geom_point()+
  scale_color_gradient(low="blue", high="red")
#   scale_color_manual(values=cols)
ggplotly(p)

reassign <- function(assignment) {
  groups <- sort(unique(assignment))
  num = length(groups)
  if(sum(groups == 1:num) != 0) {
    for(i in 1:num){
      for(j in 1:15){
        if(assignment[j]==groups[i]){
          assignment[j] = i 
        }
      }
    }
  }
  return(assignment)
}

datc_model$assignment = unlist(datc_model %>% group_by(Iter) %>% 
  group_map(~ {
    .$assignment %>% reassign}))

#cluster parameters
cluster_mean_cov = list()
for(iter in 1:20) {
  dots = filter(datc_model,Iter==iter)
  cluster_mean_cov[[iter]] = cluster_param(dots,dots$assignment,alpha,lambda,mu0,nu,psi_list[[iter]])
}


#silouette coefficient
data = datc_model %>% filter(Iter==13)
sil_coef(data)
sil_coef_vec = c()
for(iter in 1:20) {
  data = datc_model %>% filter(Iter==iter)
  sil_coef_vec = c(sil_coef_vec,sil_coef(data))
}

df = data.frame(iter = 1:20, sil_coef = sil_coef_vec)
df %>% ggplot(aes(x=iter,y=sil_coef)) + geom_point() + scale_y_continuous(limits = c(0,1))


setwd('/Users/young/Desktop/UCSD/Research/VWM_Iterated_Learning/Colored_dots/Data')
load('all.data_color.Rdata')
all.data<-filter(all.data,iter<=25)
all.data = all.data %>% mutate(assignment = ifelse(color=="green", 1, ifelse(color == "blue", 2, 3)))

s_coef = all.data %>% group_by(seed,chain,iter) %>% group_map(~sil_coef(.))
s_coef = s_coef %>% unlist() %>% rep(each=15)
mean_s_coef = all.data %>% select(seed,chain,iter) %>% cbind(s_coef=s_coef) %>% group_by(seed,chain,iter) %>% summarise(s_coef = mean(s_coef)) %>% 
  group_by(iter) %>% 
  summarise(mean_s = mean(s_coef)) 
mean_s_coef %>% ggplot(aes(x=iter,y=mean_s))+geom_point()+scale_y_continuous(limits = c(-1,1))

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
