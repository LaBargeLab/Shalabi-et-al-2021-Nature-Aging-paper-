#### This file provides the code that was used to plot the waterfall plot in figure 1 in the paper 
### after computing the ratio of the AND mask percent area to the KRT19 mask percent area (% of LEps expressing KRT14)
### This info is provided in supplementarry table 13

### The rest of the panles in Figure 1 were produced in GraphPad Prism 8.3.0.

#### import this file

###Load needed packages 
library(ggplot2)


col <- ifelse(data$Risk.State== "AR","#009296", "#BC5A42")
myplot <- ggplot(data) + 
  geom_bar(aes(reorder(Strain.Name, -AND.K19),AND.K19), fill=col, stat="identity", width=0.6)+
  ###change the / to . for the colunm name (AND/K19)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.00))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  geom_errorbar( aes(x=Strain.Name, ymin=data$AND.K19-data$SEM, ymax=data$AND.K19+data$SEM), width=0.2, colour="black", position=position_dodge(.9) , size=0.3) 
####### to remove grids and background 
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
