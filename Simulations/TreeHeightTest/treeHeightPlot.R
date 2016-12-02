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
filename <- "treeheight_3_2_1_coal_1.000_2.000_4.000_mig_0.010_0.020_0.001_0.003_0.010_0.010"

# read in the tree and log files
master <- read.nexus(sprintf("./out/%s_master.tree", filename))
esco <- read.table(sprintf("./out/%s_esco.log", filename), header=TRUE, sep="\t")
lisco <- read.table(sprintf("./out/%s_lisco.log", filename), header=TRUE, sep="\t")
sisco <- read.table(sprintf("./out/%s_sisco.log", filename), header=TRUE, sep="\t")
basta <- read.table(sprintf("./out/%s_basta.log", filename), header=TRUE, sep="\t")

# Analyse all the simulated and sampled tree files. Skip the first 2001 logs.
for (i in 1:length(master)){
  tree_name = sprintf("TREE_%d",i-1)
  if (i==1){
    height <- data.frame(master=max(nodeHeights(master[[tree_name]])),
                         esco=esco$tree_height[length(esco$Sample)-length(master)-1+i],
                         lisco=lisco$tree_height[length(lisco$Sample)-length(master)-1+i],
                         sisco=sisco$tree_height[length(sisco$Sample)-length(master)-1+i],
                         basta=basta$TreeHeightLog2[length(basta$Sample)-length(master)-1+i])
  }else{
    new.height <- data.frame(master=max(nodeHeights(master[[tree_name]])),
                         esco=esco$tree_height[length(esco$Sample)-length(master)-1+i],
                         lisco=lisco$tree_height[length(lisco$Sample)-length(master)-1+i],
                         sisco=sisco$tree_height[length(sisco$Sample)-length(master)-1+i],
                         basta=basta$TreeHeightLog2[length(basta$Sample)-length(master)-1+i])
    height <- rbind(height, new.height)
  }
}

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)


# plot the tree height distributions and save as eps
p <- ggplot() + 
  geom_line(data=height, aes(color="direct simulation",x=master),size=2, stat="density") +  
  geom_line(data=height, aes(color="ESCO",x=esco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=height, aes(color="LISCO",x=lisco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=height, aes(color="SISCO",x=sisco),size=2,linetype="dashed", stat="density") + 
  geom_line(data=height, aes(color="BASTA",x=basta),size=2,linetype="dashed", stat="density") +
  xlab("tree height") + ylab("density") + xlim(0,200) +
  scale_colour_manual("",values = c("direct simulation" = col1, "ESCO" = col3, "LISCO" = col4, "SISCO" = col2, "BASTA" = col0),
                      breaks=c("direct simulation", "ESCO", "LISCO", "SISCO", "BASTA")) 
plot(p)
ggsave(plot=p,"../../text/figures/TreeHeightTest/treeheight.eps",width=7, height=5)
  

  
    