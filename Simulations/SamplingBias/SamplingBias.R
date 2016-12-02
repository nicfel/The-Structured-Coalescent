######################################################
######################################################
# plot the sampling bias output after checking that 
# every combined has an ESS value of over 200
######################################################
######################################################
library(ggplot2)
# needed to calculate ESS values
library(coda)

# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# get the names of all first run esco log files
esco <- list.files(path="./logs", pattern="*1esco.log", full.names = TRUE)

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


# Read In Data ---------------------------------
mig_rate <- data.frame()
for (i in seq(1,length(esco),1)){
  print(i)
  # Make the filenames for all possible migration rates
  filename1 <- paste(esco[i], sep="")
  filename2 <- gsub("1esco","2esco",filename1)
  filename3 <- gsub("1esco","3esco",filename1)
  
  # Read in the asco *.logs and combine all 3 runs that had different initial values
  t_esco1 <- read.table(filename1, header=TRUE, sep="\t")
  t_esco2 <- read.table(filename2, header=TRUE, sep="\t")
  t_esco3 <- read.table(filename3, header=TRUE, sep="\t")
  
  # Read in the lisco *.logs
  t_lisco1 <- read.table(gsub("esco","lisco",filename1), header=TRUE, sep="\t")
  t_lisco2 <- read.table(gsub("esco","lisco",filename2), header=TRUE, sep="\t")
  t_lisco3 <- read.table(gsub("esco","lisco",filename3), header=TRUE, sep="\t")
  
  # Read in the sisco *.logs
  t_sisco1 <- read.table(gsub("esco","sisco",filename1), header=TRUE, sep="\t")
  t_sisco2 <- read.table(gsub("esco","sisco",filename2), header=TRUE, sep="\t")
  t_sisco3 <- read.table(gsub("esco","sisco",filename3), header=TRUE, sep="\t")
  
  # Get the true migration rate value under which simulation was performed
  splitted <- strsplit(filename1, split="_")
  bias <- splitted[[1]][6]
  migration <- splitted[[1]][4]
  samplingBias <- as.numeric(sprintf("%s",bias))
  
  # combine the data after a burn in of 20%
  t_esco1 <- t_esco1[-seq(1,ceiling(length(t_esco1$migRates1)/10)), ]
  t_esco2 <- t_esco2[-seq(1,ceiling(length(t_esco2$migRates1)/10)), ]
  t_esco3 <- t_esco3[-seq(1,ceiling(length(t_esco3$migRates1)/10)), ]
  t_esco = rbind(t_esco1,t_esco2,t_esco3)
  
  t_lisco1 <- t_lisco1[-seq(1,ceiling(length(t_lisco1$migRates1)/10)), ]
  t_lisco2 <- t_lisco2[-seq(1,ceiling(length(t_lisco2$migRates1)/10)), ]
  t_lisco3 <- t_lisco3[-seq(1,ceiling(length(t_lisco3$migRates1)/10)), ]
  t_lisco = rbind(t_lisco1,t_lisco2,t_lisco3) 
  
  t_sisco1 <- t_sisco1[-seq(1,ceiling(length(t_sisco1$migRates1)/10)), ]
  t_sisco2 <- t_sisco2[-seq(1,ceiling(length(t_sisco2$migRates1)/10)), ]
  t_sisco3 <- t_sisco3[-seq(1,ceiling(length(t_sisco3$migRates1)/10)), ]
  t_sisco = rbind(t_sisco1,t_sisco2,t_sisco3) 

  # calculate ess values
  esco_ess <- effectiveSize(t_esco)
  lisco_ess <- effectiveSize(t_lisco)
  sisco_ess <- effectiveSize(t_sisco)
  
  # Check if any ESS value calculated by the coda package is below 200
  # The first entry are the sample number, hence their ESS value is
  # not needed
  if (min(esco_ess[2:6])<200){
    print("esco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(esco_ess[2:6]),filename1))
  }
  if (min(lisco_ess[2:6])<200){
    print("lisco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(lisco_ess[2:6]),filename1))
  }
  if (min(sisco_ess[2:6])<200){
    print("sisco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(sisco_ess[2:6]),filename1))
  }
  

  # if i == 1, initialize the dataframes
  if (i==1){
    mig_rate <- data.frame(bias=as.numeric(bias), migrationrate=as.numeric(migration), 
                           esco_mig1=mean(t_esco$migRates1), esco_mig2=mean(t_esco$migRates2),
                           esco_coal1=mean(t_esco$coalRates1), esco_coal2=mean(t_esco$coalRates2),
                           lisco_mig1=mean(t_lisco$migRates1), lisco_mig2=mean(t_lisco$migRates2),
                           lisco_coal1=mean(t_lisco$coalRates1), lisco_coal2=mean(t_lisco$coalRates2),
                           sisco_mig1=mean(t_sisco$migRates1), sisco_mig2=mean(t_sisco$migRates2),
                           sisco_coal1=mean(t_sisco$coalRates1), sisco_coal2=mean(t_sisco$coalRates2))
  }else{
    new.mig_rate <- data.frame(bias=as.numeric(bias), migrationrate=as.numeric(migration), 
                           esco_mig1=mean(t_esco$migRates1), esco_mig2=mean(t_esco$migRates2),
                           esco_coal1=mean(t_esco$coalRates1), esco_coal2=mean(t_esco$coalRates2),
                           lisco_mig1=mean(t_lisco$migRates1), lisco_mig2=mean(t_lisco$migRates2),
                           lisco_coal1=mean(t_lisco$coalRates1), lisco_coal2=mean(t_lisco$coalRates2),
                           sisco_mig1=mean(t_sisco$migRates1), sisco_mig2=mean(t_sisco$migRates2),
                           sisco_coal1=mean(t_sisco$coalRates1), sisco_coal2=mean(t_sisco$coalRates2))
    mig_rate <- rbind(mig_rate, new.mig_rate)
  }

}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the migration rates
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# lisco data
lisco1 = data.frame(bias = mig_rate$bias*0, rate=mig_rate$lisco_mig1, 
                    facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                    facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
lisco2 = data.frame(bias = mig_rate$bias*0+1, rate=mig_rate$lisco_mig2, 
                    facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                    facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
# esco data
esco1 = data.frame(bias = mig_rate$bias*0+2, rate=mig_rate$esco_mig1, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
esco2 = data.frame(bias = mig_rate$bias*0+3, rate=mig_rate$esco_mig2, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
# sisco data
sisco1 = data.frame(bias = mig_rate$bias*0+4, rate=mig_rate$sisco_mig1, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
sisco2 = data.frame(bias = mig_rate$bias*0+5, rate=mig_rate$sisco_mig2, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
esco = rbind(esco1, esco2)
lisco = rbind(lisco1, lisco2)
sisco = rbind(sisco1, sisco2)

esco$facet_posy <- factor(esco$facet_posy,
                          levels = c("migration rate = 1", "migration rate = 0.1", "migration rate = 0.01"))
lisco$facet_posy <- factor(lisco$facet_posy,
                           levels = c("migration rate = 1", "migration rate = 0.1", "migration rate = 0.01"))
sisco$facet_posy <- factor(sisco$facet_posy,
                          levels = c("migration rate = 1", "migration rate = 0.1", "migration rate = 0.01"))


p_mig <- ggplot()+
  geom_violin(data=lisco, aes(factor(bias),rate,color="LISCO"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=esco, aes(factor(bias),rate,color="ESCO"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=sisco, aes(factor(bias),rate,color="SISCO"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  facet_grid(facet_posy ~ facet_posx,scales = "free_y") + scale_y_log10() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab("") +
  scale_colour_manual("",values = c("ESCO" = col3, "LISCO" = col4, "SISCO" = col2)) 
plot(p_mig)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the coalescent rates
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# lisco data
lisco1 = data.frame(bias = mig_rate$bias*0, rate=mig_rate$lisco_coal1, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
lisco2 = data.frame(bias = mig_rate$bias*0+1, rate=mig_rate$lisco_coal2, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
# esco data
esco1 = data.frame(bias = mig_rate$bias*0+2, rate=mig_rate$esco_coal1, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
esco2 = data.frame(bias = mig_rate$bias*0+3, rate=mig_rate$esco_coal2, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
# sisco data
sisco1 = data.frame(bias = mig_rate$bias*0+4, rate=mig_rate$sisco_coal1, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
sisco2 = data.frame(bias = mig_rate$bias*0+5, rate=mig_rate$sisco_coal2, 
                   facet_posx=sprintf("sample numbers: %d & %d",mig_rate$bias,100-mig_rate$bias), 
                   facet_posy=sprintf("migration rate = %g",mig_rate$migrationrate))
esco = rbind(esco1, esco2)
lisco = rbind(lisco1, lisco2)
sisco = rbind(sisco1, sisco2)

esco$facet_posy <- factor(esco$facet_posy,
                       levels = c("migration rate = 1", "migration rate = 0.1", "migration rate = 0.01"))
lisco$facet_posy <- factor(lisco$facet_posy,
                          levels = c("migration rate = 1", "migration rate = 0.1", "migration rate = 0.01"))
sisco$facet_posy <- factor(sisco$facet_posy,
                          levels = c("migration rate = 1", "migration rate = 0.1", "migration rate = 0.01"))

p_coal <- ggplot() +
  geom_violin(data=lisco, aes(factor(bias),rate,color="LISCO"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=esco, aes(factor(bias),rate,color="ESCO"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  geom_violin(data=sisco, aes(factor(bias),rate,color="SISCO"),draw_quantiles = c(0.25, 0.5, 0.75),adjust = 2) +
  facet_grid(facet_posy ~ facet_posx,scales = "free_y") + scale_y_log10() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab("") +
  scale_colour_manual("",values = c("ESCO" = col3, "LISCO" = col4, "SISCO" = col2)) 
plot(p_coal)

ggsave(plot=p_mig,"../../text/figures/SamplingBias/SamplingBiasMigrationRates.eps",width=10, height=4.5)
ggsave(plot=p_coal,"../../text/figures/SamplingBias/SamplingBiasCoalescentRates.eps",width=10, height=4.5)
