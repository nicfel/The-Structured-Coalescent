######################################################
######################################################
# Reads in the simulated and inferred tree heights
# and analyses their distribution. The sampled trees
# are analysed by looking at the logegd tree heights
# while for the master trees, the maximal node height
# is used as the tree height
######################################################
######################################################
library(ggplot2)
# needed to get the root heights
library(phytools)
# needed to read the trees
library(ape)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Set how the base of the filename is made up
filename <- c("treeheight_3_2_1_coal_1.000_2.000_4.000_mig_0.010_0.020_0.001_0.003_0.010_0.010",
               "treeheight_3_2_1_coal_1.000_2.000_4.000_mig_0.100_0.200_0.010_0.030_0.100_0.100",
               "treeheight_3_2_1_coal_1.000_2.000_4.000_mig_1.000_2.000_0.100_0.300_1.000_1.000")

for (f in 1:length(filename)){

  # read in the tree and log files
  master <- read.nexus(sprintf("./out/%s_master.tree", filename[f]))
  esco <- read.table(sprintf("./out/%s_esco.log", filename[f]), header=TRUE, sep="\t")
  lisco <- read.table(sprintf("./out/%s_lisco.log", filename[f]), header=TRUE, sep="\t")
  masco <- read.table(sprintf("./out/%s_masco.log", filename[f]), header=TRUE, sep="\t")
  sisco <- read.table(sprintf("./out/%s_sisco.log", filename[f]), header=TRUE, sep="\t")
  basta <- read.table(sprintf("./out/%s_basta.log", filename[f]), header=TRUE, sep="\t")
  
  if(f==1){overallmigration="slow"}
  if(f==2){overallmigration="medium"}
  if(f==3){overallmigration="fast"}
  
  # Analyse all the simulated and sampled tree files. Skip the first 2001 logs.
  for (i in 1:length(master)){
    tree_name = sprintf("TREE_%d",i-1)
    if (i==1 && f==1){
      height <- data.frame(master=max(nodeHeights(master[[tree_name]])),
                           esco=esco$tree_height[length(esco$Sample)-length(master)-1+i],
                           lisco=lisco$tree_height[length(lisco$Sample)-length(master)-1+i],
                           masco=masco$tree_height[length(masco$Sample)-length(master)-1+i],
                           sisco=sisco$tree_height[length(sisco$Sample)-length(master)-1+i],
                           basta=basta$TreeHeightLog2[length(basta$Sample)-length(master)-1+i],
                           overallmig=overallmigration)
    }else{
      new.height <- data.frame(master=max(nodeHeights(master[[tree_name]])),
                           esco=esco$tree_height[length(esco$Sample)-length(master)-1+i],
                           lisco=lisco$tree_height[length(lisco$Sample)-length(master)-1+i],
                           masco=masco$tree_height[length(masco$Sample)-length(master)-1+i],
                           sisco=sisco$tree_height[length(sisco$Sample)-length(master)-1+i],
                           basta=basta$TreeHeightLog2[length(basta$Sample)-length(master)-1+i],
                           overallmig=overallmigration)
      height <- rbind(height, new.height)
    }
  }
}

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)


# plot the tree height distributions and save as eps
p1 <- ggplot() + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="direct simulation",x=master),size=2, stat="density") +  
  geom_line(data=subset(height,overallmig=="fast"), aes(color="ESCO",x=esco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="MASCO",x=masco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="SISCO",x=sisco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="BASTA",x=basta),size=2,linetype="dashed", stat="density") +
  xlab("tree height") + ylab("density") + 
  scale_colour_manual(guide = FALSE,values = c("direct simulation" = col1, "ESCO" = col3, "MASCO" = col4, "SISCO" = col2, "BASTA" = col0),
                      breaks=c("direct simulation", "ESCO", "MASCO", "SISCO", "BASTA")) +
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlim(c(0,10)) + labs(title = "fast")
p2 <- ggplot() + 
  geom_line(data=subset(height,overallmig=="medium"), aes(color="direct simulation",x=master),size=2, stat="density") +  
  geom_line(data=subset(height,overallmig=="medium"), aes(color="ESCO",x=esco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="medium"), aes(color="MASCO",x=masco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="medium"), aes(color="SISCO",x=sisco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="medium"), aes(color="BASTA",x=basta),size=2,linetype="dashed", stat="density") +
  xlab("tree height") + ylab("density") + 
  scale_colour_manual(guide = FALSE,values = c("direct simulation" = col1, "ESCO" = col3, "MASCO" = col4, "SISCO" = col2, "BASTA" = col0),
                      breaks=c("direct simulation", "ESCO", "MASCO", "SISCO", "BASTA")) +
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlim(c(0,40))+ labs(title = "medium")
p3 <- ggplot() + 
  geom_line(data=subset(height,overallmig=="slow"), aes(color="direct simulation",x=master),size=2, stat="density") +  
  geom_line(data=subset(height,overallmig=="slow"), aes(color="ESCO",x=esco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="slow"), aes(color="MASCO",x=masco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="slow"), aes(color="SISCO",x=sisco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="slow"), aes(color="BASTA",x=basta),size=2,linetype="dashed", stat="density") +
  xlab("tree height") + ylab("density") + 
  scale_colour_manual(guide = FALSE,values = c("direct simulation" = col1, "ESCO" = col3, "MASCO" = col4, "SISCO" = col2, "BASTA" = col0),
                      breaks=c("direct simulation", "ESCO", "MASCO", "SISCO", "BASTA")) +
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlim(c(0,250))+ labs(title = "slow")
library(gridExtra)
library(grid)
grid.arrange(p1, p2, p3, ncol=3)
p <- arrangeGrob(p1, p2, p3, ncol=3) #generates g


ggsave(plot=p,"../../text/figures/TreeHeightTest/treeheight_figure.eps",width=8, height=3)
  
# build the legend and axis labels
pleg <- ggplot() + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="direct simulation",x=master),size=2, stat="density") +  
  geom_line(data=subset(height,overallmig=="fast"), aes(color="ESCO",x=esco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="MASCO",x=masco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="SISCO",x=sisco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=subset(height,overallmig=="fast"), aes(color="BASTA",x=basta),size=2,linetype="dashed", stat="density") +
  xlab("tree height") + ylab("density") + 
  scale_colour_manual("",values = c("direct simulation" = col1, "ESCO" = col3, "MASCO" = col4, "SISCO" = col2, "BASTA" = col0),
                      breaks=c("direct simulation", "ESCO", "MASCO", "SISCO", "BASTA")) + 
  theme(legend.position="top")

  
plot(pleg)
ggsave(plot=pleg,"../../text/figures/TreeHeightTest/treeheight_legend.eps",width=8, height=3)
