######################################################
######################################################
# plots the proability of the root being a state 
# conditioned on the migration rate from state 1 to 
# state 2.
######################################################
######################################################
library(ggplot2)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# read in the root state probabilities from a text file
# created in matlab
stateProbs <- read.table("rootStateProbabilities.txt", header=TRUE, sep="\t")

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)

p <- ggplot(data=stateProbs)+
  geom_line(aes(color="MultiTypeTree",x=migrationrate,y=MTT),size=2) +
  geom_line(aes(color="ESCO",x=migrationrate,y=esco),size=2,linetype="longdash") +
  geom_line(aes(color="MASCO",x=migrationrate,y=Masco),size=2,linetype="dashed") +
  geom_line(aes(color="SISCO",x=migrationrate,y=Sisco),size=2) + 
  xlab("migration rate     to     ") + ylab("P(     )") +
  scale_x_log10(breaks=c(0.001,0.01,0.1)) +
  scale_colour_manual("",values = c("MultiTypeTree"=col1, "ESCO"=col3,"MASCO"=col4,"SISCO"=col2),
                      breaks=c("MultiTypeTree", "ESCO", "MASCO", "SISCO")) 
plot(p)
ggsave(plot=p,"../../text/figures/RootStateProbabilities/rootStateProbabilties.eps",width=4.5, height=2.5)


