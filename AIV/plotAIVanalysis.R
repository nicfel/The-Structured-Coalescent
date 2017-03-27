######################################################
######################################################
# plots the mcc tree from the sisco AIV analysis
# also plots the root states tree height etc
######################################################
######################################################
# needed to plot the mcc tree
library(ape)
library(ggplot2)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the mcc tree
tr <- read.nexus(file="./structcoal/sisco_mcc.trees")

# ladderize the tree
tr <- ladderize(tr, right = F)

# get the tip labels
tip_labs = tr$tip.label;

# define the colors for the different states
alaska <- rgb(red=0, green=0.4470, blue=0.7410)
northwest <- rgb(red=0.8500, green=0.3250, blue=0.0980)
northeast <- rgb(red=0.9290, green=0.6940, blue=0.1250)
southeast <- rgb(red=0.4940, green=0.1840, blue=0.5560)
eastcoast <- rgb(red=0.4660, green=0.6740, blue=0.1880)
northmideast <- rgb(red=0.3010, green=0.7450, blue=0.9330)
center <- rgb(red=0.6350, green=0.0780, blue=0.1840)
black <- rgb(red=0, green=0, blue=0)

# initialize the node colors
tip_col = matrix(, nrow = length(tip_labs))
for (i in 1: length(tip_col)){
  
  l1 <- regexpr('alaska', tip_labs[i])
  l2 <- regexpr('northwest', tip_labs[i])
  l3 <- regexpr('northeast', tip_labs[i])
  l4 <- regexpr('southeast', tip_labs[i])
  l5 <- regexpr('eastcoast', tip_labs[i])
  l6 <- regexpr('northmideast', tip_labs[i])
  l7 <- regexpr('center', tip_labs[i])
  
  # assine colors to the tips
  if (l1[1]>0){
    tip_col[i] = alaska
  }
  if (l2[1]>0){
    tip_col[i] = northwest
  }
  if (l3[1]>0){
    tip_col[i] = northeast
  }
  if (l4[1]>0){
    tip_col[i] = southeast
  }
  if (l5[1]>0){
    tip_col[i] = eastcoast
  }
  if (l6[1]>0){
    tip_col[i] = northmideast
  }
  if (l7[1]>0){
    tip_col[i] = center
  }
  
}

# get the time of the root. The nodeHeigths function is in the phytools lubrary
library(phytools)
root_heigth <- max(nodeHeights(tr))

# plot the sisco mcc tree (doesn't really matter which mcc tree is used)
setEPS()
postscript("../text/figures/AIV/sisco_mcc.eps",width=5, height=5)
plot(tr,show.tip.label=F,root.edge = T)
tiplabels(pch = 19, col = tip_col, adj = 0.5, cex = 0.3)
# add a scale with the most recent sample being sampled in 2012.47123287671
axisPhylo(side = 1,root.time = 2012.47123287671-root_heigth, backward = F)
dev.off()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the root state probabilities from a txt file generated in
# matlab that used all treea after a burnin to get the mean root
# state probabilities
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot the inferred root state probabilities
t <- read.table(file="rootStates.txt", header=TRUE, sep="\t")


# complicated way of reading in the data one by one to plot them in ggplot
plot.t1_masco <- data.frame(place="alaska",prob=t$alaska[1],method = "MASCO",x=1)
plot.t1_volz <- data.frame(place="alaska",prob=t$alaska[2],method = "SISCO",x=1)
plot.t2_masco <- data.frame(place="northwest",prob=t$northwest[1],method = "MASCO",x=2)
plot.t2_volz <- data.frame(place="northwest",prob=t$northwest[2],method = "SISCO",x=2)
plot.t3_masco <- data.frame(place="northeast",prob=t$northeast[1],method = "MASCO",x=3)
plot.t3_volz <- data.frame(place="northeast",prob=t$northeast[2],method = "SISCO",x=3)
plot.t4_masco <- data.frame(place="southeast",prob=t$southeast[1],method = "MASCO",x=4)
plot.t4_volz <- data.frame(place="southeast",prob=t$southeast[2],method = "SISCO",x=4)
plot.t5_masco <- data.frame(place="eastcoast",prob=t$eastcoast[1],method = "MASCO",x=5)
plot.t5_volz <- data.frame(place="eastcoast",prob=t$eastcoast[2],method = "SISCO",x=5)
plot.t6_masco <- data.frame(place="northmideast",prob=t$northmideast[1],method = "MASCO",x=6)
plot.t6_volz <- data.frame(place="northmideast",prob=t$northmideast[2],method = "SISCO",x=6)
plot.t7_masco <- data.frame(place="center",prob=t$center[1],method = "MASCO",x=7)
plot.t7_volz <- data.frame(place="center",prob=t$center[2],method = "SISCO",x=7)

plot.t_masco <- rbind(plot.t1_masco,plot.t2_masco,
                plot.t3_masco,plot.t4_masco,
                plot.t5_masco,plot.t6_masco,
                plot.t7_masco)
plot.t_sisco <- rbind(plot.t1_volz,plot.t2_volz,
                plot.t3_volz,plot.t4_volz,
                plot.t5_volz,plot.t6_volz,
                plot.t7_volz)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pie_masco <- ggplot(plot.t_masco, aes(x="", y=prob*100, fill=place)) + geom_bar(stat="identity",width=1) 
pie_masco <- pie_masco + coord_polar(theta = "y") +
  scale_fill_manual("",values = c("alaska" = alaska,
                                  "northwest" = northwest,
                                  "northeast" = northeast,
                                  "southeast" = southeast,
                                  "eastcoast" = eastcoast,
                                  "northmideast" = northmideast,
                                  "center" = center)) +  blank_theme
plot(pie_masco)

pie_sisco <- ggplot(plot.t_sisco, aes(x="", y=prob*100, fill=place)) + geom_bar(stat="identity",width=1) 
pie_sisco <- pie_sisco + coord_polar(theta = "y") +
  scale_fill_manual("",values = c("alaska" = alaska,
                                  "northwest" = northwest,
                                  "northeast" = northeast,
                                  "southeast" = southeast,
                                  "eastcoast" = eastcoast,
                                  "northmideast" = northmideast,
                                  "center" = center)) +  blank_theme
plot(pie_sisco)


# print the masco and sisco pies
ggsave(plot=pie_masco,"../text/figures/AIV/RootStatesMasco.eps",width=5, height=5)

ggsave(plot=pie_sisco,"../text/figures/AIV/RootStatesSisco.eps",width=5, height=5)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the legend for the different sampling states
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# data frame to plot points and colors, the points are not needed,
# but the legend is
d <- data.frame(x=c(1,2,3,4,5,6,7),y=c(1,2,3,4,5,6,7),lab=c("Alaska",
                                            "North West",
                                            "North East",
                                            "South West",
                                            "East Coast",
                                            "North Mid East",
                                            "Center"))

p_leg <- ggplot()+
  geom_point(data=d, aes(x=x,y=y,color=lab),size=5) +
  theme(
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(face = "bold", size = 20),
        legend.position="bottom") + xlab("A B C MASCO SISCO") +
  scale_colour_manual("",values = c("Alaska" = alaska,
                                    "North West" = northwest,
                                    "North East" = northeast,
                                    "South West" = southeast,
                                    "East Coast" = eastcoast,
                                    "North Mid East" = northmideast,
                                    "Center" = center),
                      breaks=c("Alaska",
                                   "North West",
                                   "South West",
                                   "Center",
                                   "East Coast",
                                   "North Mid East",
                                   "North East"),
                      guide=guide_legend(nrow=3, size=5)) 
plot(p_leg)
ggsave(plot=p_leg,"../text/figures/AIV/legend.eps",width=5, height=2)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the distribution of coalescent rates
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Read in all the log files of the 3 different runs
log_masco1 <- read.table(file="./structcoal/individual/AIV_1masco.log", header=TRUE, sep="\t")
log_masco2 <- read.table(file="./structcoal/individual/AIV_2masco.log", header=TRUE, sep="\t")
log_masco3 <- read.table(file="./structcoal/individual/AIV_3masco.log", header=TRUE, sep="\t")
log_sisco1 <- read.table(file="./structcoal/individual/AIV_1sisco.log", header=TRUE, sep="\t")
log_sisco2 <- read.table(file="./structcoal/individual/AIV_2sisco.log", header=TRUE, sep="\t")
log_sisco3 <- read.table(file="./structcoal/individual/AIV_3sisco.log", header=TRUE, sep="\t")

# Combine the legends after a burn in of 10%
log_masco <- rbind(log_masco1[40:length(log_masco3$Sample),],
                   log_masco2[40:length(log_masco3$Sample),],
                   log_masco3[40:length(log_masco3$Sample),])
log_sisco <- rbind(log_sisco1[40:length(log_masco3$Sample),],
                   log_sisco2[40:length(log_masco3$Sample),],
                   log_sisco3[40:length(log_masco3$Sample),])

log_masco$method = "MASCO"
log_sisco$method = "SISCO"

log <- rbind(log_masco, log_sisco)
  
  
# plot the tree height distributions and save as eps
p_violin <- ggplot(data=log) +
  geom_violin(aes(1,coalRates1,color="Alaska")) + 
  geom_violin(aes(2,coalRates2,color="North West")) + 
  geom_violin(aes(3,coalRates3,color="North East")) + 
  geom_violin(aes(4,coalRates4,color="South West")) + 
  geom_violin(aes(5,coalRates5,color="East Coast")) + 
  geom_violin(aes(6,coalRates6,color="North Mid East")) + 
  geom_violin(aes(7,coalRates7,color="Center"))
# Set the seed such that the plotted prior distribution is always exactly the same
set.seed(1)
exprnd <- data.frame(vals=rexp(100000,1))
p_violin <- p_violin +
  geom_violin(data=exprnd, aes(8,vals,color="prior")) +
  facet_grid(method ~.) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 20, colour = black),
        axis.title=element_text(size=20)) + ylab("coalescent rates") +
  scale_colour_manual("",values = c("Alaska" = alaska,
                                    "North West" = northwest,
                                    "North East" = northeast,
                                    "South West" = southeast,
                                    "East Coast" = eastcoast,
                                    "North Mid East" = northmideast,
                                    "Center" = center,
                                    "prior" = black),
                      breaks=c("Alaska",
                               "North West",
                               "South West",
                               "Center",
                               "East Coast",
                               "North Mid East",
                               "North East",
                               "prior")) 
p_violin <- p_violin + annotate("text", x=8, y=0.7, label="prior")
plot(p_violin)

ggsave(plot=p_violin,"../text/figures/AIV/coalescentRates.eps",width=8, height=4)
