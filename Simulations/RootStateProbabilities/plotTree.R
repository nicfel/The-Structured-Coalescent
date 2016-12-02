######################################################
######################################################
# plots the tree for the figure in the root state
# probability part. Also plots the trees used for
# the figure in the methods part
######################################################
######################################################
library(ape)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# tree for which the root state probabilities were inferred
tree <- "(((inv3loc_0:2.5,inv4loc_1:8):4,(inv5loc_1:10,inv6loc_1:6):5.8):27,(inv2loc_0:10,inv1loc_0:20):22):0.0;"
tr <- read.tree(text=tree)

setEPS()
postscript("../../text/figures/RootStateProbabilities/tree.eps",width=5, height=5)
plot(tr,show.tip.label=F,root.edge = T,edge.width=8)
tiplabels(pch = 19, col = c("blue", "chocolate4", "chocolate4", "chocolate4", "blue", "blue"), adj = .5, cex = 2)
dev.off()

# Plot trees for the methods figure ---------------------------------

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)

# Plot two times the same tree in diffferent colors
tree <- "((lin3_0:10,lin4_1:10):10,lin6_1:15);"
tr <- read.tree(text=tree)
setEPS()
postscript("../../text/figures/MethodsTree/methodsTree1.eps",width=10, height=5)
plot(tr,show.tip.label=F,root.edge = F,edge.width=8, edge.color = col1)
dev.off()

setEPS()
postscript("../../text/figures/MethodsTree/methodsTree2.eps",width=10, height=5)
plot(tr,show.tip.label=F,root.edge = F,edge.width=8, edge.color = col0)
dev.off()
