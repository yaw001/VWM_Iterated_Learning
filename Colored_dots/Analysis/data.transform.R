library(rjson)
library(dplyr)
library(ggplot2)
library(gganimate)
library(plotly)
getwd()
setwd('/Users/young/Desktop/UCSD/Research/Color-Iterated-Learning/VSS/')
all.data = data.frame()
for(file.name in list.files(pattern = '*.json')){
  json_file<-fromJSON(file = file.name)
  for(i in 1:length(json_file[['stimuli']])){
    all.data <- rbind(all.data,
                      data.frame(seed = json_file[['seed']],
               chain = json_file[['chain']],
               iter = json_file[['iter']],
               x=json_file[["stimuli"]][[i]][['center']][['x']],
               y=json_file[["stimuli"]][[i]][['center']][['y']],
               color = json_file[["stimuli"]][[i]][['color']]))
  }
}

all.data <- all.data %>% mutate(x=as.numeric(x),
                    y = as.numeric(y))
all.data$color<-as.character(all.data$color)
all.data$color[which(all.data$color=="#bf404d")]<-"red"
all.data$color[which(all.data$color=="#4cbf40")]<-"green"
all.data$color[which(all.data$color=="#404cbf")]<-"blue"

save(all.data, file='all.data.Rdata')
load("all.data.Rdata")
p1 <- all.data %>% filter(iter == 15, seed == 9, chain == 0) %>%
  ggplot(aes(x=x, y=y, color=color))+
  # facet_grid(chain~seed)+
  geom_point(size=4)+
  theme_bw()+
  scale_color_manual(values=c("red"="red",
                              "green"="#00aa00",
                              "blue"="blue"))+
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none')
p1 

p <- all.data %>% filter(iter <= 20) %>%
  ggplot(aes(x=x, y=y, color=color, frame=iter))+
  facet_grid(chain~seed)+
  geom_point(size=2.5)+
  theme_bw()+
  scale_color_manual(values=c("red"="red",
                              "green"="#00aa00",
                              "blue"="blue"))+
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none')

animation::ani.options(ani.res=2000)
gganimate(p, ani.width=600, ani.height=600, interval=0.7)
gganimate(p, "output.html")
gganimate(p, "output.gif")


p <- all.data %>% filter(iter <= 20, seed == 1, chain == 4) %>%
  plot_ly(x = ~x,
          y = ~y,
          color = ~color,
          frame = ~iter,
          type = "scatter",
          mode = "markers")
p


