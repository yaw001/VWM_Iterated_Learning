#Data Analysis and Graphs
library(tidyverse)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggfortify)
library(cluster)
library(reshape)
setwd('/Users/young/Desktop/UCSD/Research/Color-Iterated-Learning/Data/')
load('all.data.Rdata')
all.data<-filter(all.data,iter<=25)
#analysis of the data 
tb.result<-all.data%>%filter(iter<=25) %>% 
  group_by(seed,chain,iter,color)%>% 
  do(tibble(xy = list(data.frame(x=.$x, y=.$y)),
            covariance=list(cov(cbind(.$x,.$y)))))%>%
  rowwise()%>%
  mutate(eigenval = list(eigen(covariance)),
         eigenvec = list(eigen(covariance)$vector),
         det = det(covariance),
         ratio = eigenval$values[1]/eigenval$values[2],
         proj.x = list(as.matrix(xy)%*%eigenval$vectors[,1]),
         proj.y = list(as.matrix(xy)%*%eigenval$vectors[,2]))

# tb.result%>%group_by(iter,color)%>%summarise(mean.ratio = mean(ratio))%>%
#   ggplot(aes(x=iter,y=mean.ratio,color=color))+geom_point()
# partition<-all.data%>%
#   group_by(seed,chain,iter)%>%
#   do(tibble(xy = list(data.frame(x=.$x, y=.$y)),
#             cl=list(kmeans(cbind(.$x,.$y),3))))%>%
#   rowwise()%>%
#   mutate('cl.vec'=list(cl$cluster))

#Determinant of Covariance matrix Contingent on Colors
mean.det<-all.data%>%group_by(seed,chain,iter,color)%>%
  summarise(det.cov=det(cov(cbind(x,y))))%>%
  group_by(iter,color)%>%
  summarise(N=n(),mean.det=mean(det.cov),se.mean=sd(det.cov)/sqrt(N))

# iter_count<-all.data%>%group_by(seed,chain)%>%summarise(max.iter=max(iter))
# min_iter<-min(iter_count$max.iter)

mean.det%>%ggplot(aes(x=iter,mean,y=mean.det,color=color))+
  geom_point()+
  geom_line(size=1.1)+
  scale_color_manual(values = c('red'='red',
                                'green'='green2',
                                'blue'='blue')) +
  geom_errorbar(aes(ymin=mean.det-se.mean,ymax=mean.det+se.mean),width=0.5,size=1.1)+
  labs(x = 'iteration',
       y = 'Mean Determinant of Covariance Matrix')+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = 'none')
  # ggtitle('Mean Determinant of Covariance matrix Contigent on Colors')

mean.det_all<-all.data%>%group_by(seed,chain,iter)%>%
  summarise(det.cov=det(cov(cbind(x,y))))%>%
  group_by(iter)%>%
  summarise(mean.det=mean(det.cov),se.mean=sd(det.cov)/sqrt(100))

mean.det_all %>%  lm(formula=mean.det~iter+I(iter^2)) %>% summary()
mean.det_all%>%ggplot(aes(x=iter,y=mean.det))+
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


#Color Randomization
color.rand<-function(col){
  index<-sample(length(col),replace=F)
  col.rand<-col[index]
  return(col.rand)
}
randcol<-as.vector(replicate(26*100,color.rand(rep(c('red','blue','green'),5))))
dat_randcol<-all.data %>% filter(iter>=0,iter<=25)
dat_randcol['color']=randcol
  
mean.det.rc<-dat_randcol%>%group_by(seed,chain,iter,color)%>%
  summarise(det.cov=det(cov(cbind(x,y))))%>%
  group_by(iter,color)%>%
  summarise(mean.det.rc=mean(det.cov),se.mean.rc=sd(det.cov)/sqrt(100))

df.mean.det<-data.frame(iter=c(mean.det$iter,rep(0:25,each=3)),
                        color=c(mean.det$color,mean.det.rc$color),
                        mean.det=c(mean.det$mean.det,mean.det.rc$mean.det.rc),
                        se.mean=c(mean.det$se.mean,mean.det.rc$se.mean.rc),
                        category=c(rep('Observed',78),rep('Color randomized',78)))
non_rand=all.data%>%group_by(seed,chain,iter,color)%>%
  summarise(det.cov=det(cov(cbind(x,y)))) %>% cbind(category=as.factor(rep(1,7800)))
rand=dat_randcol%>%group_by(seed,chain,iter,color)%>%
  summarise(det.cov=det(cov(cbind(x,y)))) %>% cbind(category=as.factor(rep(2,7800)))
df.det=rbind(non_rand,rand)
df.det$category=df.det$category %>% as.factor()
mod1=df.det %>% filter(color=="red") %>% lm(formula=det.cov~iter*category)
mod2=df.det %>% filter(color=="red") %>% lm(formula=det.cov~iter) 
anova(mod1,mod2)
df.det %>% filter(color=="red") %>% ggplot(aes(x=iter,y=det.cov,color=category))+geom_point() + geom_smooth(method = "lm",se=F)
df.mean.det$category = ifelse(df.mean.det$category=="Observed",0,1) 
df.mean.det %>% filter(color=="blue",category=="Observed") %>% lm(formula=mean.det~iter) %>% summary() 
df.mean.det %>% filter(color=="blue",category=="Color randomized") %>% lm(formula=mean.det~iter) %>% summary()
mod1=df.mean.det %>%  lm(formula=mean.det~iter+I(iter^2))
mod2=df.mean.det %>% lm(formula=mean.det~iter+I(iter^2) + iter:category)
summary(mod1)
summary(mod2)
anova(mod1,mod2)

blue_rand = filter(df.mean.det,category=="Color randomized",color == "blue")
red_rand = filter(df.mean.det,category=="Color randomized",color == "red")
green_rand = filter(df.mean.det,category=="Color randomized",color == "green")
blue_ob = filter(df.mean.det,category=="Observed",color == "blue")
red_ob = filter(df.mean.det,category=="Observed",color == "red")
green_ob = filter(df.mean.det,category=="Observed",color == "green")
ks.test(blue_rand$mean.det,blue_ob$mean.det,alternative = "less")
ks.test(red_rand$mean.det,red_ob$mean.det,alternative = "less")
ks.test(green_rand$mean.det,green_ob$mean.det,alternative = "less")

df.mean.det%>%
  ggplot(aes(x=iter,y=mean.det,color=color))+
  geom_point(aes(colour=color))+
  geom_line(aes(linetype=category,color=color),size=1)+
  scale_color_manual(values = c('red'='red',
                                'green'='green2',
                                'blue'='blue'))+
                                # 'black'='black'))
  scale_linetype_manual(values = c('Observed'="solid",
                                   'Color randomized'="dotted"), name="category")+
  geom_errorbar(aes(ymin=mean.det-se.mean,ymax=mean.det+se.mean),width=0.4,size=0.9)+
  labs(x = 'Iteration',
       y = 'Covariance determinant (group spread)')+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 16, face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.title = element_blank())

# mean.det<-all.data%>%group_by(seed,chain,iter)%>%
#   summarise(det.cov=det(cov(cbind(x,y))))%>%
#   group_by(iter)%>%
#   summarise(mean.det=mean(det.cov),se.mean=sd(det.cov)/sqrt(100))
# mean.det%>%ggplot(aes(x=iter,y=mean.det))+
#   geom_point()+
#   geom_line(size=1.1)+
#   geom_errorbar(aes(ymin=mean.det-se.mean,ymax=mean.det+se.mean),width=0.5,size=1.1)+
#   labs(x= 'iteration',
#        y = 'Mean Determinant of Covariance Matrix')+
#   theme_bw()+
#   theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5))


#ratio. 1st pc/2nd pc.
ratio<-all.data%>%group_by(seed,chain,iter,color)%>%
  summarise(ratio=eigen(cov(cbind(x,y)))$values[1]/eigen(cov(cbind(x,y)))$values[2]) %>% 
  group_by(iter) %>% 
  summarise(mean_ratio = mean(ratio), se.mean = sd(ratio)/sqrt(100))

# ratio %>% ggplot(aes(x=iter,y=mean_ratio)) +
#   geom_point()+
#   geom_line(size=1.1)+
#   geom_errorbar(aes(ymin=mean_ratio-se.mean,ymax=mean_ratio+se.mean),width=0.5,size=1.1)+
#   labs(x= 'iteration',
#        y = 'Mean Determinant of Covariance Matrix')+
#   theme_bw()+
#   theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5))



#0.7 quantile = 10.741309
ratio[which(ratio$ratio<5&ratio$ratio>1),]
#seed 0 chain 0 iter 1  (1.08)
#seed 0 chain 6 iter 3  (3.05)
#seed 1 chain 1 iter 4  (10)
#seed 6 chain 8 iter 13 (20.1)
#seed 2 chain 2 iter 22 (74.4)
#seed 1 chain 3 iter 22 (39634)
# ex1<-filter(all.data,seed==9&chain==0&iter==0)
# ex1%>%ggplot(aes(x=x,y=y,color=color))+
#   geom_point(size=2)+
#   stat_ellipse(type='norm',
#                level=0.95,
#                size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = 'none',
#         panel.border = element_blank())
# 
# ex2<-filter(all.data,seed==9&chain==0&iter==10)
# ex2%>%ggplot(aes(x=x,y=y,color=color))+
#   geom_point(size=2)+
#   stat_ellipse(type='norm',
#                level=0.95,
#                size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = 'none',
#         panel.border = element_blank())
# 
# ex3<-filter(all.data,seed==9&chain==0&iter==15)
# ex3%>%ggplot(aes(x=x,y=y,color=color))+
#   geom_point(size=2)+
#   stat_ellipse(type='norm',
#                level=0.95,
#                size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = 'none',
#         panel.border = element_blank())
# 
# ex4<-filter(all.data,seed==9&chain==0&iter==20)
# ex4%>%ggplot(aes(x=x,y=y,color=color))+
#   geom_point(size=2)+
#   stat_ellipse(type='norm',
#                level=0.95,
#                size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = 'none',
#         panel.border = element_blank())
# 
# ex5<-filter(all.data,seed==9&chain==0&iter==25)
# ex5%>%ggplot(aes(x=x,y=y,color=color))+
#   geom_point(size=2)+
#   stat_ellipse(type='norm',
#                level=0.95,
#                size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = 'none',
#         panel.border = element_blank())
# 

count.ratio.2<-function(ratio,threshold1,threshold2){
  count=length(which(ratio>threshold1&ratio<threshold2))
  return(prop=count/(3*10))
}

ratio<-all.data%>%group_by(seed,chain,iter,color)%>%
  summarise(ratio=eigen(cov(cbind(x,y)))$values[1]/eigen(cov(cbind(x,y)))$values[2])
quantile(ratio$ratio,probs = seq(0,1,0.25))

# 2.6/5/12

# fourth quartile (linear)
threshold1 = 14
threshold2 = 100000000

prop.result.1<-ratio%>%
  group_by(seed,iter)%>%summarise(prop=count.ratio.2(ratio,threshold1,threshold2))%>%
  group_by(iter)%>%summarise(mean.prop=mean(prop),se.prop=sd(prop)/sqrt(10))

df.ratio.1<-data.frame(iter=0:25,mean.prop=prop.result.1$mean.prop,prop.se=prop.result.1$se.prop)
df.ratio.1 <- df.ratio.1 %>% mutate(type = "aalinear(ratio>14)")

# df.ratio.1%>%ggplot(aes(x = iter,y = mean.prop))+geom_point()+
#   # coord_cartesian(ylim = c(0,0.3,0.05))+
#   scale_y_continuous(limits = c(0,0.4), breaks = seq(0,1,0.05))+
#   geom_line(size=1.2)+
#   # geom_hline(aes(yintercept = mean(proportion)),color='red')+
#   # geom_hline(aes(yintercept = mean(proportion)+sd(proportion)),color='blue',linetype='dashed')+
#   # geom_hline(aes(yintercept = mean(proportion)-sd(proportion)),color='blue',linetype='dashed')+
#   theme_bw()+
#   theme(panel.grid = element_blank())+
#   geom_errorbar(aes(ymax=mean.prop+prop.se,ymin=mean.prop-prop.se),size=0.9,width=0.5)+
#   labs(x='iteration',
#        y='proportion')+
#   theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5))

#Isotropic (first quartile)
threshold1=1
threshold2=2.6
prop.result.2<-ratio%>%
  group_by(seed,iter)%>%summarise(prop=count.ratio.2(ratio,threshold1,threshold2))%>%
  group_by(iter)%>%summarise(mean.prop=mean(prop),se.prop=sd(prop)/sqrt(10))
df.ratio.2<-data.frame(iter=0:25,mean.prop=prop.result.2$mean.prop,prop.se=prop.result.2$se.prop)
df.ratio.2<-df.ratio.2 %>% mutate(type = "ccisotropic(1<ratio<2.6)")

df.ratio.2%>%ggplot(aes(x = iter,y = mean.prop))+geom_point()+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  geom_line(size=1.2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=mean.prop+prop.se,ymin=mean.prop-prop.se),size=0.9,width=0.5)+
  labs(x='iteration',
       y='proportion')+
  theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
        panel.grid = element_blank())

# anisotropic (second quartile + third)
threshold1=2.6
threshold2=14
prop.result.3<-ratio%>%
  group_by(seed,iter)%>%summarise(prop=count.ratio.2(ratio,threshold1,threshold2))%>%
  group_by(iter)%>%summarise(mean.prop=mean(prop),se.prop=sd(prop)/sqrt(10))

df.ratio.3<-data.frame(iter=0:25,mean.prop=prop.result.3$mean.prop,prop.se=prop.result.3$se.prop)
df.ratio.3<-df.ratio.3 %>% mutate(type = "bbAnisotropic(2.6<ratio<14")

df.ratio.3%>%ggplot(aes(x = iter,y = mean.prop))+geom_point()+
  scale_y_continuous(limits = c(0.1,0.4), breaks = seq(0.1,0.4,0.05))+
  geom_line(size=1.2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=mean.prop+prop.se,ymin=mean.prop-prop.se),size=0.9,width=0.5)+
  labs(x='iteration',
       y='proportion')+
  theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
        panel.grid = element_blank())

# # anisotropic (third quartile)
# threshold1=5
# threshold2=12
# prop.result.4<-ratio%>%
#   filter(iter>=0)%>%
#   group_by(seed,iter)%>%summarise(prop=count.ratio.2(ratio,threshold1,threshold2))%>%
#   group_by(iter)%>%summarise(mean.prop=mean(prop),se.prop=sd(prop)/sqrt(10))
# 
# df.ratio.4<-data.frame(iter=0:28,mean.prop=prop.result.4$mean.prop,prop.se=prop.result.2$se.prop)
# df.ratio.4<-df.ratio.4 %>% mutate(type = "ratio(5-12)")

ratio.df <- bind_rows(df.ratio.1,df.ratio.2,df.ratio.3)
ggplot(ratio.df, aes(iter,mean.prop)) + geom_area(aes(fill = type))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
  labs(x = "Iteration", y = "Proportion of cluster shapes") + 
  # scale_fill_discrete(labels = c("Linear(ratio>12.4)","Anisotropic(2.6<ratio<12.4)","Isotropic(1<ratio<2.6)"))+
  scale_fill_manual(labels = c("Linear(ratio>14)","Anisotropic(2.6<ratio<14)","Isotropic(1<ratio<2.6)"),
                    values = c("#CC6666", "#9999CC", "#66CC99"))+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
        panel.grid = element_blank(),
        legend.title = element_blank())

prop.result.1 %>% lm(formula = mean.prop~iter) %>% summary()
prop.result.2 %>% lm(formula = mean.prop~iter) %>% summary()
prop.result.3 %>% lm(formula = mean.prop~iter) %>% summary()
# #extreme cases
# threshold1 = 75
# threshold2 = 10000000
# tb.result%>%filter(iter==0&chain==0&threshold1<ratio&ratio<threshold2)%>%
#   group_by(chain,iter,color)
# tb.result%>%filter(seed==7&chain==0&color=='green')%>%
#   ggplot(aes(x=iter,y=ratio))+geom_point()+geom_line()
# #Test the pattern duality hypothesis
# binwidth=0.1
# ratio.breaks<-quantile(ratio$ratio,probs = seq(0,1,binwidth))
# ratio %>% 
#   mutate(ratio.quantile = cut(ratio, 
#                               breaks = ratio.breaks,
#                               labels = seq(0.1,1,0.1))) %>%
#   ggplot(aes(x=iter, fill=ratio.quantile))+
#   geom_bar(position='fill')
# 
# quantile(ratio$ratio,probs = seq(0,1,binwidth))
# 
# ratio_quantiles<-ratio%>%mutate(Ratio.quantile = cut(ratio,breaks = ratio.breaks,
#                                                      labels = seq(binwidth,1,binwidth)))
# 
# count_ratio.3<-function(ratio.quantile){
#   count=length(ratio.quantile)
#   return(prop=count/(3*10))
# }
# 
# prop.result.3<-ratio_quantiles%>%filter(iter>0)%>%
#   group_by(seed,iter,Ratio.quantile)%>%summarise(prop=count_ratio.3(Ratio.quantile))%>%
#   group_by(iter,Ratio.quantile)%>%summarise(mean.prop=mean(prop),se.prop=sd(prop)/sqrt(10))
# 
# df.ratio.3<-as.data.frame(prop.result.3)
# df.ratio.3<-na.omit(df.ratio.3)
# df.ratio.3%>%ggplot(aes(x = iter,y = mean.prop))+geom_point()+
#   geom_line(size=1.2)+
#   theme_bw()+
#   theme(panel.grid = element_blank())+
#   geom_errorbar(aes(ymax=mean.prop+se.prop,ymin=mean.prop-se.prop),size=0.9,width=0.5)+
#   labs(title=sprintf('Proportion of ratios in 10 quantiles'),
#        x='iteration',
#        y='proportion')+
#   theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5))+
#   facet_wrap(~Ratio.quantile,nrow=5)


#Orientation and length



# orientation difference in angle
angle.diff <- function(vec.1, vec.2){
  return(acos(vec.1%*%vec.2)*180/pi)
}

vec.extract <- function(list){
  len = length(list)
  vec = c()
  for(i in 1:len){
    vec = c(vec, list[[i]][,1])
  }
  return(matrix(vec,nrow = 2, ncol = len))
}

angles <- function(mat) {
  angle.list = c()
  if (ncol(mat) == 1) {
    return(0)
  }else{
    for(j in 1:(ncol(mat)-1)){
      for(i in (j+1):ncol(mat)){
        angle.list = c(angle.list, angle.diff(mat[,j],mat[,i]))
      }
    }
    return(mean(angle.list))
  }
}

len.diff = function(len.1,len.2) {
  return(abs(len.1 - len.2))
}

lens <- function(vec) {
  lens.list = c()
  for(j in 1:(length(vec))-1){
    for(i in (j+1):(length(vec))){
      lens.list = c(lens.list, len.diff(vec[j],vec[i]))
    }
  }
  return(mean(lens.list))
}




tb.result = tb.result %>% group_by(seed,chain,iter) %>% mutate(mean.len = lens(len))

linear.result = linear.result %>% rowwise() %>% mutate(len=sqrt((max(proj.x) - min(proj.x))^2 + 
                                                                  (proj.y[which.max(proj.x)] - proj.y[which.min(proj.x)])^2))
linear.result = linear.result %>% group_by(seed,chain,iter) %>% mutate(mean.len = lens(len))


summary.diff = linear.result %>% group_by(seed,chain,iter) %>% summarise(angle.meandiff = unique(mean.angle), len.meandiff = unique(mean.len))

counts.cal <- function(diff,threshold1, threshold2){
  counts = sum(diff>threshold1 & diff<threshold2)
  return(counts)
}

#orientation analysis

quantile(summary.diff$angle.meandiff,probs = seq(0,1,0.25))
#0-5/5-12/12-55/55-170 (12 cutoff 0.75 quantile)
#0-15/15-45/45-90/90-170
summary_counts = summary.diff %>% group_by(iter) %>% count() 

threshold1 = 0
threshold2 = 15
prop1 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(angle.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))
prop1 = prop1 %>% mutate(angle_diff = "4")

threshold1 = 15
threshold2 = 45
prop2 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(angle.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))

prop2 = prop2 %>% mutate(angle_diff = "3")

threshold1 = 45
threshold2 = 90
prop3 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(angle.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))

prop3 = prop3 %>% mutate(angle_diff = "2")

threshold1 = 90
threshold2 = 17000000
prop4 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(angle.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))

prop4 = prop4 %>% mutate(angle_diff = "1")

prop.dist = bind_rows(prop1,prop2,prop3,prop4)




ggplot(prop.dist, aes(iter,prop)) + geom_area(aes(fill = angle_diff))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
  labs(x = "iteration", y = "proportion of angles") + 
  scale_fill_discrete(labels = c(">90°","45°- 90°","15°- 45°","<15°"))+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
        panel.grid = element_blank(),
        legend.title = element_blank())


hist(summary.diff$angle.meandiff)
hist(summary.diff$len.meandiff)

isotropic.result= tb.result %>% filter(iter>=0,ratio < 2.59)

linear.result= tb.result %>% filter(iter>=0,ratio > 14.14)

linear.result <- linear.result %>% group_by(seed,chain,iter) %>% mutate(mean.angle = angles(vec.extract(eigenvec)))

isotropic.result <- isotropic.result %>% group_by(seed,chain,iter) %>% mutate(mean.angle = angles(vec.extract(eigenvec)))

linear.result[which(linear.result$mean.angle!=0),]$iter %>% unique() %>% length()

# angle.sim=table(linear.result$mean.angle)[table(linear.result$mean.angle) == 2] %>% names() %>% as.numeric() 
# linear.result[linear.result$mean.angle%in%angle.sim,] %>% print(n=100)

linear.result = linear.result[which(linear.result$mean.angle!=0),]

isotropic.result = isotropic.result[which(isotropic.result$mean.angle!=0),]

#vector sum
tb.result <- tb.result %>% rowwise() %>% mutate(PC1.vec = list(eigenvec[,1]))
new.pc1 = tb.result$PC1.vec
length(new.pc1)

min_angle <- function(vec.list) {
  len <-  length(vec.list)
  vec.list_all <- c(vec.list,Map("*",vec.list,-1))
  comb = combn(1:(2*len),len)
  if(len == 2) {
    comb = comb[,-(which(comb[2,] - comb[1,] == len))]
  }else{
    comb = comb[,-(which(comb[2,] - comb[1,] == len | comb[3,] - comb[1,] == len | comb[3,] - comb[2,] == len))]
  }
  list_sum = list()
  for(i in 1:ncol(comb)){
    list_sum = c(list_sum,list(Reduce("+",vec.list_all[comb[,i]])))
  }
  vec.sum = lapply(list_sum,norm,"2")
  return(max(unlist(vec.sum)))
}

#aggregate
summary_vec.sum = tb.result %>% 
  group_by(seed,chain,iter) %>% mutate(vec.sum = min_angle(PC1.vec)) %>% 
  group_by(seed,chain,iter) %>% summarise(iter_vec.sum = unique(vec.sum)) %>% group_by(iter) %>% 
  summarise(n=n(),mean_vec.sum = mean(iter_vec.sum),se_vec.sum = sd(iter_vec.sum)/sqrt(n))

summary_vec.sum %>% lm(formula=mean_vec.sum~iter) %>% summary()
summary_vec.sum %>% ggplot(aes(x = iter,y = mean_vec.sum))+geom_point()+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  geom_point(size = 2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=mean_vec.sum+se_vec.sum,ymin=mean_vec.sum-se_vec.sum),size=0.9,width=0.5)+
  # geom_hline(aes(yintercept=random_mean),color='blue',linetype='dashed',size=1.5)+
  # geom_hline(aes(yintercept=3),color='blue',linetype='dashed',size=1.5)+
  labs(x='Iteration',
       y='Orientation similarity')+
  theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
        panel.grid = element_blank())

smoothing <- function(vec.sum,vec.sum.se,smooth_idx) {
  smoothed_se = rep(NA,length(vec.sum))
  for(i in 1:length(vec.sum)){
    if(i == 1){
      vec.sum[i] = mean(vec.sum[i:(i+smooth_idx)])
      smoothed_se[i] = mean(vec.sum.se[i:(i+smooth_idx)])
    }else if(i == length(vec.sum)){
      vec.sum[i] = mean(vec.sum[(i-smooth_idx):i])
      smoothed_se[i] = mean(vec.sum.se[i:(i-smooth_idx)])
    }else{
      vec.sum[i] = mean(vec.sum[(i-smooth_idx/2):(i+smooth_idx/2)])
      smoothed_se[i] = mean(vec.sum.se[(i-smooth_idx/2):(i+smooth_idx/2)])
    }
  }
  return(list(vec.sum=vec.sum,
              smoothed_se = smoothed_se))
}

#linear
summary_vec.sum_linear = linear.result %>% 
  group_by(seed,chain,iter) %>% mutate(vec.sum = min_angle(PC1.vec)) %>% 
  group_by(seed,chain,iter) %>% summarise(iter_vec.sum = unique(vec.sum)) %>% group_by(iter) %>% 
  summarise(n=n(),mean_vec.sum = mean(iter_vec.sum),se_vec.sum = sd(iter_vec.sum)/sqrt(n)) %>% ungroup() %>% 
  mutate(smoothed_mean = smoothing(mean_vec.sum,se_vec.sum,2)[['vec.sum']],smoothed_se = smoothing(mean_vec.sum,se_vec.sum,2)[['smoothed_se']])

iter_0 = data.frame(iter = 0, n = 0, mean_vec.sum = 2.18, se_vec.sum = 0.456, smoothed_mean = mean(c(2.18,1.86)), smoothed_se = mean(0.456,0.066))

summary_vec.sum_linear = rbind(iter_0,summary_vec.sum_linear)
summary_vec.sum_linear %>% ggplot(aes(x = iter,y = smoothed_mean))+geom_point()+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  geom_point(size = 2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=smoothed_mean+smoothed_se,ymin=smoothed_mean-smoothed_se),size=0.9,width=0.5)+
  # geom_hline(aes(yintercept=random_mean),color='blue',linetype='dashed',size=1.5)+
  # geom_hline(aes(yintercept=3),color='blue',linetype='dashed',size=1.5)+
  labs(x='iteration',
       y='Magnitude of directional vector sum')+
  theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
        panel.grid = element_blank())

summary_vec.sum_linear %>% lm(formula=mean_vec.sum~iter) %>% summary()

#isotropic
summary_vec.sum_isotropic = isotropic.result %>% 
  group_by(seed,chain,iter) %>% mutate(vec.sum = min_angle(PC1.vec)) %>% 
  group_by(seed,chain,iter) %>% summarise(iter_vec.sum = unique(vec.sum)) %>% group_by(iter) %>% 
  summarise(n=n(),mean_vec.sum = mean(iter_vec.sum),se_vec.sum = sd(iter_vec.sum)/sqrt(n)) %>% ungroup() %>% 
  mutate(smoothed_mean = smoothing(mean_vec.sum,se_vec.sum,2)[['vec.sum']],smoothed_se = smoothing(mean_vec.sum,se_vec.sum,2)[['smoothed_se']])

summary_vec.sum_isotropic %>% ggplot(aes(x = iter,y = smoothed_mean))+geom_point()+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  geom_point(size = 2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=smoothed_mean+smoothed_se,ymin=smoothed_mean-smoothed_se),size=0.9,width=0.5)+
  # geom_hline(aes(yintercept=random_mean),color='blue',linetype='dashed',size=1.5)+
  # geom_hline(aes(yintercept=3),color='blue',linetype='dashed',size=1.5)+
  labs(x='iteration',
       y='Magnitude of directional vector sum')+
  theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
        panel.grid = element_blank())

summary_vec.sum_isotropic %>% lm(formula=mean_vec.sum~iter) %>% summary()

#isotropic and linear
isotropic_linear = rbind(summary_vec.sum_isotropic,summary_vec.sum_linear) 
summary_isotropic_linear = isotropic_linear %>% mutate(Type = c(rep("Isotropic",26),rep("Linear",26)))

summary_isotropic_linear %>% ggplot(aes(x = iter,y = smoothed_mean,color=Type))+geom_point()+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  scale_color_manual(values = c('Isotropic'='orange',
                                'Linear'='purple'))+
  geom_point(size = 2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=smoothed_mean+smoothed_se,ymin=smoothed_mean-smoothed_se),size=0.9,width=0.5)+
  # geom_hline(aes(yintercept=random_mean),color='blue',linetype='dashed',size=1.5)+
  # geom_hline(aes(yintercept=3),color='blue',linetype='dashed',size=1.5)+
  labs(x='iteration',
       y='Magnitude of directional vector sum')+
  theme(text = element_text(family='Times New Roman', size= 16, face='bold'),
        panel.grid = element_blank(),
        legend.text=element_text(size=16),
        axis.title = element_blank())

acos((2-0.4)/2/1)/(pi/180)
# random
x = runif(30000,-1,1)
sign=sample(c(1,-1),size=30000,replace = T)
y = sqrt(1-x^2)*sign
z = cbind(x,y)
z=split(z,seq(nrow(z)))
random_list = tibble(index = rep(1:10000,each =3), lengths=z)
random_summary_vec.sum = random_list %>% 
  group_by(index) %>% mutate(vec.sum = min_angle(lengths))
random_mean=mean(unique(random_summary_vec.sum$vec.sum))

# tb.result %>% filter(iter>0,iter<=20) %>% rowwise() %>% mutate(vec.sum_mag = norm(vec.sum,"2")) %>% 
#   group_by(seed,chain,iter) %>% 
#   summarise(vec.sum.mag = unique(vec.sum_mag)) %>% ggplot(aes(x=iter, y=))

summary_vec.sum %>% ggplot(aes(x = iter,y = mean_vec.sum))+geom_point()+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  geom_point(size = 2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=mean_vec.sum+se_vec.sum,ymin=mean_vec.sum-se_vec.sum),size=0.9,width=0.5)+
  geom_hline(aes(yintercept=random_mean),color='red',linetype='dashed',size=1.5)+
  geom_hline(aes(yintercept=3),color='blue',linetype='dashed',size=1.5)+
  labs(x='Iteration',
       y='Mean magnitude of vector sum')+
  theme(text = element_text(family='Times New Roman', size= 16, face='bold'),
        panel.grid = element_blank(),
        axis.title = element_blank())

summary(lm(mean_vec.sum~iter,data = summary_vec.sum))

cor.test(summary_vec.sum$mean_vec.sum,summary_vec.sum$iter)

summary_vec.sum_linear %>% ggplot(aes(x = iter,y = mean_vec.sum))+geom_point()+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  geom_point(size = 2)+
  geom_smooth(method="lm",se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=mean_vec.sum+se_vec.sum,ymin=mean_vec.sum-se_vec.sum),size=0.9,width=0.5)+
  # geom_hline(aes(yintercept=random_mean),color='blue',linetype='dashed',size=1.5)+
  # geom_hline(aes(yintercept=3),color='blue',linetype='dashed',size=1.5)+
  labs(x='iteration',
       y='Mean magnitude of vector sum')+
  theme(text = element_text(family='Times New Roman', size= 15, face='bold'),
        panel.grid = element_blank())

#length

#variance of length
tb.result<-all.data%>%filter(iter>=0,iter<=25) %>% 
  group_by(seed,chain,iter,color)%>% 
  do(tibble(xy = list(data.frame(x=.$x, y=.$y)),
            covariance=list(cov(cbind(.$x,.$y)))))%>%
  rowwise()%>%
  mutate(eigenval = list(eigen(covariance)),
         eigenvec = list(eigen(covariance)$vector),
         det = det(covariance),
         ratio = eigenval$values[1]/eigenval$values[2],
         proj.x = list(as.matrix(xy)%*%eigenval$vectors[,1]),
         proj.y = list(as.matrix(xy)%*%eigenval$vectors[,2]))

tb.result = tb.result%>% rowwise() %>% mutate(len=sqrt((max(proj.x) - min(proj.x))^2 + 
                                                         (proj.y[which.max(proj.x)] - proj.y[which.min(proj.x)])^2))
identity = data.frame(seed_1=tb.result$seed, chain_1=tb.result$chain, iter_1 = tb.result$iter)

tb.result_rand = tb.result %>% rowwise() %>% 
  mutate(seed_chain = paste(seed,chain,iter)) %>% 
  group_by(iter) %>% mutate(newchain = sample(seed_chain,size = n())) %>% arrange(newchain)

tb.result_rand=tb.result %>% group_by(iter) %>% sample_n(size=100)

tb.result %>% rowwise() %>% 
  mutate(seed_chain = paste(seed,chain,iter)) %>% 
  group_by(seed,iter,color) %>% mutate(newchain = sample(seed_chain,size = n())) %>% ungroup() %>% arrange(newchain)

tb.result_rand$seed_chain = NULL
tb.result_rand$newchain = NULL

tb.result_rand = tb.result_rand %>% bind_cols(identity)
tb.result_original = tb.result %>% bind_cols(identity)

combined_result = bind_rows(tb.result_original,tb.result_rand)

Type = c(rep("Observed",7800),rep("Randomized",7800))
combined_result = bind_cols(combined_result,as.data.frame(Type))


summary_len.var = combined_result %>% group_by(Type,seed_1,chain_1,iter_1) %>% summarise(var.len = var(len)) %>% group_by(Type,iter_1) %>% 
  summarise(n=n(),mean_var = mean(var.len), se = sd(var.len)/sqrt(n))


# result_len = list()
# for(i in 1:100){
#   tb.result_rand = tb.result %>% rowwise() %>%
#     mutate(seed_chain = paste(seed,chain,iter)) %>%
#     group_by(iter,color) %>% mutate(newchain = sample(seed_chain,size = n())) %>% ungroup() %>% arrange(newchain)
# 
#   # tb.result %>% rowwise() %>%
#   #   mutate(seed_chain = paste(seed,chain,iter)) %>%
#   #   group_by(seed,iter,color) %>% mutate(newchain = sample(seed_chain,size = n())) %>% ungroup() %>% arrange(newchain)
# 
#   tb.result_rand$seed_chain = NULL
#   tb.result_rand$newchain = NULL
# 
#   tb.result_rand = tb.result_rand %>% bind_cols(identity)
#   rand_result = tb.result_rand %>% group_by(seed_1,chain_1,iter_1) %>% summarise(var.len = var(len)) %>% group_by(iter_1) %>%
#     summarise(n=n(),mean_var = mean(var.len), se = sd(var.len))
#   result_len[[i]] = rand_result
# }

identity_rand = data.frame(seed_1=rep(rep(rep(0:9,each=10),3),26), chain_1=rep(rep(0:9,30),26), iter_1 = rep(0:25,each=300))

result_len = list()
for(i in 1:100){
  tb.result_rand=tb.result %>% group_by(iter,color) %>% sample_n(size=100)
  tb.result_rand = tb.result_rand %>% bind_cols(identity_rand)
  rand_result = tb.result_rand %>% group_by(seed_1,chain_1,iter_1) %>% summarise(var.len = var(len)) %>% group_by(iter_1) %>%
      summarise(n=n(),mean_var = mean(var.len), se = sd(var.len))
  result_len[[i]] = rand_result
}

summary_rand_len = bind_rows(result_len) %>% group_by(iter_1) %>% summarise(n=n(),mean_var = mean(mean_var), se = mean(se)/sqrt(n))
summary_original_len = tb.result_original %>% group_by(seed_1,chain_1,iter_1) %>% summarise(var.len = var(len)) %>% group_by(iter_1) %>%
  summarise(n=n(),mean_var = mean(var.len), se = sd(var.len)/sqrt(n))
combined_result = bind_rows(summary_rand_len,summary_original_len)

Type = c(rep("Randomized",26),rep("Observed",26))

combined_result = bind_cols(combined_result,as.data.frame(Type))


identity_rand = data.frame(seed_1=rep(rep(0:9,each=10),3), chain_1=rep(0:9,30), iter_1 = rep(0,each=300))

sim_var=list()
for(j in 1:100){
  mean_var=c()
  for(i in 0:25) {
    rand_iter_result = tb.result %>% filter(iter==i) %>% group_by(color) %>% sample_n(size=100)
    rand_iter_result = rand_iter_result %>% bind_cols(identity_rand)
    rand_var_iter = rand_iter_result %>% group_by(seed_1, chain_1, iter_1) %>% summarise(var.len=var(len))
    mean_var_iter = mean(rand_var_iter$var.len)
    mean_var = c(mean_var,mean_var_iter)
  }
  sim_var[[j]] = data.frame(iter=0:25, mean_var)
}

summary_sim_var=bind_rows(sim_var) %>% group_by(iter) %>% summarise(n=n(),mean_var = mean(mean_var), se = mean(mean_var)/sqrt(n))
summary_data = tb.result %>% group_by(seed,chain,iter) %>% summarise(var.len = var(len)) %>% group_by(iter) %>%
  summarise(n=n(),mean_var = mean(var.len), se = sd(var.len)/sqrt(n))

combined_result = bind_rows(summary_sim_var,summary_data)

Type = c(rep("Randomized",26),rep("Observed",26))

combined_result = bind_cols(combined_result,as.data.frame(Type))


combined_result%>%ggplot(aes(x = iter,y = mean_var, color= Type))+
  # scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.05))+
  # geom_line(size=1.2)+
  geom_point(size=2)+
  scale_color_manual(values = c('Observed'='orange',
                                'Randomized'='purple'))+
  geom_smooth(method = 'lm',se=F)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_errorbar(aes(ymax=mean_var+se,ymin=mean_var-se),size=0.9,width=0.5)+
  labs(x='Iteration',
       y='Variance of lengths between groups')+
  theme(text = element_text(family='Times New Roman', size= 16, face='bold'),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.text=element_text(size=16))

a=c(1,2,3,1,2,3,1,2,3,1.5,1.7,1.9,1.5,1.7,1.9,1.5,1.7,1.9)
var(c(1,2,3))
var(c(1.5,1.7,1.9))

a=c(1,2,3,1,2,3,1,2,3,4,6,8,4,6,8,4,6,8)
var(c(4,6,8))

mean(c(var(a[1:3]),var(a[4:6]),var(a[7:9]),var(a[10:12]),var(a[13:15]),var(a[16:18])))
b=replicate(100,sample(a))
mean(apply(b,2,function(x) mean(c(var(x[1:3]),var(x[4:6]),var(x[7:9]),var(x[10:12]),var(x[13:15]),var(x[16:18])))))

var_seed = c()
for(i in 0:9){
  var_each = tb.result %>% filter(seed==i,chain==0,iter==0) %>% pull(len) %>% var()
  var_seed=c(var_seed,var_each)
}
max(var_seed)
min(var_seed)

mod1=combined_result %>% lm(formula = mean_var~ iter+Type)
summary(mod1)
mod2=combined_result %>% lm(formula = mean_var~iter*Type)
summary(mod2)
anova(mod1,mod2)


quantile(summary.diff$len.meandiff,probs = seq(0,1,0.25))
#0.05/0.1/0.2 (12 or 17)
#0.1/0.18/0.28 (2)
summary_counts = summary.diff %>% group_by(iter) %>% count() 

threshold1 = 0
threshold2 = 0.05
prop_1 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(len.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))
prop_1 = prop_1 %>% mutate(len_diff = "4")

threshold1 = 0.05
threshold2 = 0.12
prop_2 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(len.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))

prop_2 = prop_2 %>% mutate(len_diff = "3")

threshold1 = 0.12
threshold2 = 0.22
prop_3 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(len.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))

prop_3 = prop_3 %>% mutate(len_diff = "2")

threshold1 = 0.22
threshold2 = 17000000
prop_4 = summary.diff %>% group_by(iter) %>% summarise(true.counts = counts.cal(len.meandiff,threshold1,threshold2)) %>% 
  ungroup() %>% mutate(prop = true.counts/summary_counts$n, smooth_prop = smoothing(prop,4))

prop_4 = prop_4 %>% mutate(len_diff = "1")

prop_dist = bind_rows(prop_1,prop_2,prop_3,prop_4)

ggplot(prop_dist, aes(iter,smooth_prop)) + geom_area(aes(fill = len_diff))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
  labs(x = "iteration", y = "proportion of length difference") + 
  scale_fill_discrete(labels = c(">0.22","0.12-0.22","0.05-0.12","<0.05"))+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
        panel.grid = element_blank(),
        legend.title = element_blank())


#Regularity analysis
ratio<-all.data%>%group_by(seed,chain,iter,color)%>%
  summarise(ratio=eigen(cov(cbind(x,y)))$values[1]/eigen(cov(cbind(x,y)))$values[2])
quantile(ratio$ratio,probs = seq(0,1,0.05))

dist.mat <- function(x){
  x <- sort(x)
  x <- (x-min(x))/(max(x)-min(x))
  m <- matrix(0, ncol=length(x), nrow=length(x))
  return(matrix(x[row(m)]-x[col(m)], ncol=length(x)))}

std = dist.mat(1:5)
diag(std[-nrow(std),-1])
upperstd = std[upper.tri(std)]
dist.result = tb.result %>% filter(ratio > 14,iter>0) %>% rowwise() %>% 
  mutate(dist.vec = list(dist.mat(proj.x)[upper.tri(dist.mat(proj.x))]), rss = sum((dist.vec-upperstd)^2)) %>% 
  group_by(iter) %>% 
  summarise(mean_rss = mean(rss), n=n(),se_rss = sd(rss)/sqrt(n))
dist.result %>% lm(formula=mean_rss~iter) %>% summary()
dist.result%>%ggplot(aes(x=iter,y=mean_rss))+
  geom_point()+
  # geom_line(size=1.1)+
  geom_smooth(method = "lm",se=F)+
  geom_errorbar(aes(ymin=mean_rss-se_rss,ymax=mean_rss+se_rss),width=0.5,size=1.1)+
  labs(x= 'Iteration',
       y = 'Deviation from equally spaced dot distribution')+
  theme_bw()+
  theme(text = element_text(family='Times New Roman', size= 12, face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5))

# standard deviations
paired.dist_sd<-tb.result%>%filter(ratio>165)%>%
  .$proj.x%>%lapply(dist.mat)%>%
  lapply(., function(x) diag(x[-nrow(x),-1]))%>%
  lapply(sd)
paired.dist_sd
mean.dat<-mean(unlist(paired.dist_sd))
#Null Hypothesis
uniform_sd<-replicate(100000,list(runif(3)))%>%lapply(.,function(x) c(0,x,1))%>%
  lapply(dist.mat)%>%lapply(., function(x) diag(x[-nrow(x),-1]))%>%
  lapply(sd)
mean.unif<-mean(unlist(uniform_sd))
df.uniform_sd<-unlist(uniform_sd)%>%as.data.frame()
colnames(df.uniform_sd)<-'sd'
df.uniform_sd%>%ggplot(aes(sd))+geom_histogram(binwidth = 0.005, color='black')+
  # geom_vline(aes(xintercept=mean.unif),color='blue',linetype='dashed',size=2)
  geom_vline(aes(xintercept=mean.dat),color='red',linetype='dashed',size=1)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank())+
  # annotate("text",x=0.09,y=420,label='0.0569')+
  # annotate("text",x=0.24,y = 550,label='0.208')
  labs(title='Null distribution of standard deviations of paired distance',
       x = "standard deviations of paired distance",
       y = "counts")+
  theme(text = element_text(family='Times New Roman', size= 10, face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
p.value<-sum(df.uniform_sd$sd<mean.dat)/100000
p.value  


#regularity in isotropic structures
eucl.dist<-function(dots.vec,center){
  center.vec=matrix(rep(center,each=5),nrow=5)
  sqrt(rowSums((dots.vec-center.vec)^2))
}
get.center<-function(xy){
  return(c(center.x=mean(xy$x),
         center.y=mean(xy$y)))
}
paired.dist<-function(xy){
  dist.matrix=as.matrix(dist(xy))
  c(diag(dist.matrix[-5,-1]),dist.matrix[1,5])
}

#1, 1.78, 2.32, 2.96, 3.80, 5.05,7.1,10.74,20.3,74.53,55669.95

tb.result.2<-tb.result%>%filter(ratio>1&ratio<1000000&iter>0&iter<20)%>%
  rowwise()%>%
  mutate(center=list(get.center(xy)),
         cent.dist=list(eucl.dist(as.matrix(xy),melt(center)$value)),
         sd.cent.dist=sd(cent.dist),
         pair.dist=list(paired.dist(xy)),
         sd.pair.dist=sd(pair.dist))

# normalized<-function(dist.vec){
#   dist.vec<-sort(dist.vec)
#   (dist.vec-min(dist.vec))/(max(dist.vec)-min(dist.vec))
# }
# paired.dist_sd<-tb.result.2%>%
#   .$pair.dist%>%lapply(., function(x) dist.mat(c(0,cumsum(sort(x)))))%>%
#   lapply(., function(x) diag(x[-nrow(x),-1]))%>%
#   lapply(sd)
# mean.dat.2<-mean(unlist(paired.dist_sd))
# uniform_sd<-replicate(10000,list(runif(4)))%>%lapply(.,function(x) c(0,x,1))%>%
#   lapply(dist.mat)%>%lapply(., function(x) diag(x[-nrow(x),-1]))%>%
#   lapply(sd)
# df.uniform_sd<-unlist(uniform_sd)%>%as.data.frame()
# colnames(df.uniform_sd)<-'sd'
# df.uniform_sd%>%ggplot(aes(sd))+geom_histogram(binwidth = 0.005, color='black')+
#   geom_vline(aes(xintercept=mean.dat.2),color='red',linetype='dashed',size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.ticks = element_blank())+
#   labs(title='Null distribution of standard deviations of paired distance',
#        x = "standard deviations of paired distance",
#        y = "counts")+
#   theme(text = element_text(family='Times New Roman', size= 10, face='bold'),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())
# p.value<-sum(df.uniform_sd$sd<0.06)/10000
# p.value  

center.dist_sd<-tb.result.2%>%filter(ratio<1.78)%>%
  .$cent.dist%>%lapply(., function(x) dist.mat(cumsum(sort(x))))%>%
  lapply(., function(x) diag(x[-nrow(x),-1]))%>%
  lapply(sd)
mean.dat.center<-mean(unlist(center.dist_sd))
uniform_sd<-replicate(100000,list(runif(5)))%>%
  lapply(dist.mat)%>%lapply(., function(x) diag(x[-nrow(x),-1]))%>%
  lapply(sd)
df.uniform_sd<-unlist(uniform_sd)%>%as.data.frame()
colnames(df.uniform_sd)<-'sd'
df.uniform_sd%>%ggplot(aes(sd))+geom_histogram(binwidth = 0.005, color='black')+
  geom_vline(aes(xintercept=mean.dat.center),color='red',linetype='dashed',size=2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank())+
  labs(title='Null distribution of standard deviations of distance to center',
       x = "standard deviations of distance to center",
       y = "counts")+
  theme(text = element_text(family='Times New Roman', size= 10, face='bold'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())
p.value<-sum(df.uniform_sd$sd<mean.dat.center)/100000
p.value  
ratio  %>%  mutate(group = ifelse(iter < 5, 1, ifelse(iter > 20, 2, 3))) %>% 
  filter(group %in% c(1, 2)) %>% 
  ggplot(aes(x=(ratio)))+geom_histogram()+geom_vline(xintercept = (1.78), color='red')+geom_vline(xintercept = (74.53), color='red')+scale_x_log10()+facet_grid(group~.)

#Polar system analysis

# angle.cal<-function(center,xy){
#   new.xy<-xy - matrix(rep(center,each=5),nrow=5)
#   tan.angle<-new.xy[,2]/new.xy[,1]
#   angle<-atan(tan.angle)
#   angle[angle<0]=2*pi+angle[angle<0]
#   sort.angle<-sort(angle)
#   angles<-(c(sort.angle[2]-sort.angle[1],
#              sort.angle[3]-sort.angle[2],
#              sort.angle[4]-sort.angle[3],
#              2*pi+sort.angle[1]-sort.angle[4]))
#   return(angles)
# }
# 
# tb.result.angle<-tb.result.2%>%filter(ratio>75)
#   rowwise()%>%
#   mutate(angles=list(angle.cal(center,xy)))
#   
# angle_sd<-tb.result.angle%>%
#   .$angles%>%lapply(., function(x) dist.mat(cumsum(sort(x))))%>%
#   lapply(., function(x) diag(x[-nrow(x),-1]))%>%
#   lapply(sd)
# angle_sd<-tb.result.angle%>%
#   .$angles%>%lapply(sd)
# mean.sd.angle<-mean(unlist(angle_sd))
# uniform_sd<-replicate(10000,list(runif(5)))%>%
#   lapply(dist.mat)%>%lapply(., function(x) diag(x[-nrow(x),-1]))%>%
#   lapply(sd)
# uniform_sd<-replicate(10000,list(runif(5,0,2*pi)))%>%lapply(.,sort)%>%lapply(.,diff)%>%
#   lapply(sd)
# df.uniform_sd<-unlist(uniform_sd)%>%as.data.frame()
# colnames(df.uniform_sd)<-'sd'
# df.uniform_sd%>%ggplot(aes(sd))+geom_histogram(binwidth = 0.005, color='black')+
#   geom_vline(aes(xintercept=mean.sd.angle),color='red',linetype='dashed',size=2)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.ticks = element_blank())+
#   labs(title='Null distribution of standard deviations of paired distance',
#        x = "standard deviations of paired distance",
#        y = "counts")+
#   theme(text = element_text(family='Times New Roman', size= 10, face='bold'),
#         panel.grid = element_blank(),
#         plot.title = element_text(hjust=0.5),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())
# p.value<-sum(df.uniform_sd$sd<mean.sd.angle)/10000
# p.value  
