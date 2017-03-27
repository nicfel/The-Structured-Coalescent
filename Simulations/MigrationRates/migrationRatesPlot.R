######################################################
######################################################
# This script analyses the output of the structured 
# coalescent for a 8 tip tree sampled homochronously 
# with all samples from the same deme
######################################################
######################################################
library(ggplot2)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get the names of all log files in the esco directory and from there
# derive the names of the log files in the other directories
log_files <- list.files(path="./logs", pattern="*esco.log", full.names = TRUE)

# Read In Data ---------------------------------

migrationRates <- data.frame()
for (i in 1:length(log_files)){
  # Make the filenames for all possible migration rates
  filename <- paste(log_files[i], sep="")
  
  # Get the true migration rate value under which simulation was performed
  splitted <- strsplit(filename, split="_")[[1]][2]
  migration_rate <- as.numeric(splitted)
  
  # Read in the structcoal *.logs
  t_esco <- read.table(filename, header=TRUE, sep="\t")
  t_masco <- read.table(gsub("esco","masco",filename), header=TRUE, sep="\t")
  t_sisco <- read.table(gsub("esco","sisco",filename), header=TRUE, sep="\t")
  
  # get the index of the maximum posterior estimate
  m_esco <- which.max(t_esco$posterior)
  m_masco <- which.max(t_masco$posterior)
  m_sisco <- which.max(t_sisco$posterior)
  
  # check if the maximum posterior(== likelihood since no priors) estimate
  # is either the first or the last value in the list, if this is the case
  # the explored parameter space wasn't large enough for this tree
  if(m_esco>1 && m_masco>1 && m_sisco>1 && 
     m_esco<length(t_esco$posterior) && 
     m_masco<length(t_masco$posterior) &&
     m_sisco<length(t_sisco$posterior)){
    # initialize the data frames
    if (i==1){
      rate <- data.frame(rate=migration_rate, 
                         esco=t_esco$migRates1[m_esco],
                         masco=t_masco$migRates1[m_masco],
                         sisco=t_sisco$migRates1[m_sisco])
    }else{
      new.rate <- data.frame(rate=migration_rate, 
                             esco=t_esco$migRates1[m_esco],
                             masco=t_masco$migRates1[m_masco],
                             sisco=t_sisco$migRates1[m_sisco])
      rate <- rbind(rate, new.rate)
    }
    # print the index of the map estimate if its the last or first
    # for all three methods to see for which the explored space
    # was too small
  }else{
    print(m_esco)
    print(length(t_esco$posterior))
    print(m_masco)
    print(length(t_masco$posterior))
    print(m_sisco)
    print(filename)
  }
}

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


p <- ggplot(data=rate)+
  geom_point(aes(color="ESCO",x=rate,y=esco,shape="ESCO"),size=2) +
  geom_point(aes(color="MASCO",x=rate,y=masco,shape="MASCO"),size=1) +
  geom_point(aes(color="SISCO",x=rate,y=sisco,shape="SISCO"),size=1) +
  xlab("true migration rate") + ylab("estimated migration rate")  +
  scale_colour_manual("",values = c("ESCO" = col3, "MASCO" = col4, "SISCO" = col2)) +
  scale_shape_manual("",values = c("ESCO" = 3, "MASCO" = 16, "SISCO" = 16)) +
  scale_x_log10() + scale_y_log10()  + guides(colour = guide_legend(override.aes = list(size=2))) + theme(legend.position="top")

d <- data.frame(x=c(0.00001,1), y=c(0.00001,1))
p <- p + geom_line(data=d,aes(x=x, y=y), color="red", linetype="dashed")
plot(p)
ggsave(plot=p,"../../text/figures/migrationRates/migrationRateEstimates.eps",width=4, height=4)


