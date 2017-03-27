######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
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

# get the names of all SISCO first (of three) run log files
log <- list.files(path="./logs", pattern="*1sisco.log", full.names = TRUE)

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


# Read In Data ---------------------------------
for (i in seq(1,length(log),1)){
  print(i)
  # Make the filenames for all the three runs
  filename1 <- paste(log[i], sep="")
  filename2 <- gsub("1sisco","2sisco",filename1)
  filename3 <- gsub("1sisco","3sisco",filename1)
  
  # Read in the SISCO *.logs
  t_sisco1 <- read.table(filename1, header=TRUE, sep="\t")
  t_sisco2 <- read.table(filename2, header=TRUE, sep="\t")
  t_sisco3 <- read.table(filename2, header=TRUE, sep="\t")
  
  # Read in the MASCO *.logs
  t_masco1 <- read.table(gsub("sisco","masco",filename1), header=TRUE, sep="\t")
  t_masco2 <- read.table(gsub("sisco","masco",filename2), header=TRUE, sep="\t")
  t_masco3 <- read.table(gsub("sisco","masco",filename3), header=TRUE, sep="\t")
  
  # combine all tree runs after a burn in of 10%
  t_masco1 <- t_masco1[-seq(1,ceiling(length(t_masco1$migRates1)/10)), ]
  t_masco2 <- t_masco2[-seq(1,ceiling(length(t_masco2$migRates1)/10)), ]
  t_masco3 <- t_masco3[-seq(1,ceiling(length(t_masco3$migRates1)/10)), ]
  t_masco = rbind(t_masco1,t_masco2,t_masco3) 
  
  t_sisco1 <- t_sisco1[-seq(1,ceiling(length(t_sisco1$migRates1)/10)), ]
  t_sisco2 <- t_sisco2[-seq(1,ceiling(length(t_sisco2$migRates1)/10)), ]
  t_sisco3 <- t_sisco3[-seq(1,ceiling(length(t_sisco3$migRates1)/10)), ]
  t_sisco = rbind(t_sisco1,t_sisco2,t_sisco3) 
  
  # calculate ess values
  masco_ess <- effectiveSize(t_masco)
  sisco_ess <- effectiveSize(t_sisco)
  
  # Check if any ESS value calculated by the coda package is below 200
  # The first entry are the sample number, hence their ESS value is
  # not needed
  if (min(masco_ess[2:6])<200){
    print("masco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(masco_ess[2:6]),filename1))
  }
  if (min(sisco_ess[2:6])<200){
    print("sisco ESS value to low")
    print(sprintf("ESS value is %f for file %s",min(sisco_ess[2:6]),filename1))
  }
  
  # Get the true asymmetry under which simulation was performed
  # this is encoded in the file name
  splitted <- strsplit(filename1, split="_")
  asymmetry <- as.numeric(splitted[[1]][3])
  
  # get the ratios
  c_ratio_masco = t_masco$coalRates1/t_masco$coalRates2
  c_ratio_sisco = t_sisco$coalRates1/t_sisco$coalRates2
  
  m_ratio_masco = t_masco$migRates1/t_masco$migRates2
  m_ratio_sisco = t_sisco$migRates1/t_sisco$migRates2
  
  
  # if i == 1, initialize the dataframes
  if (i==1){
    
    # check if the coalescent or the migration rates were asymmetric
    if (splitted[[1]][2]=="coal"){
      ratios_tmp1 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(c_ratio_masco,0.5),
                                upper_masco=quantile(c_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(c_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(c_ratio_sisco,0.5),
                                upper_sisco=quantile(c_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(c_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric coalescent rates symmetric migration rates",
                                facet_posy="coalescent rate ratios")
      ratios_tmp2 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(m_ratio_masco,0.5),
                                upper_masco=quantile(m_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(m_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(m_ratio_sisco,0.5),
                                upper_sisco=quantile(m_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(m_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric coalescent rates symmetric migration rates",
                                facet_posy="migration rate ratios")
    }else{
      ratios_tmp1 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(c_ratio_masco,0.5),
                                upper_masco=quantile(c_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(c_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(c_ratio_sisco,0.5),
                                upper_sisco=quantile(c_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(c_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric migration rates symmetric coalescent rates",
                                facet_posy="coalescent rate ratios")
      ratios_tmp2 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(m_ratio_masco,0.5),
                                upper_masco=quantile(m_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(m_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(m_ratio_sisco,0.5),
                                upper_sisco=quantile(m_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(m_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric migration rates symmetric coalescent rates",
                                facet_posy="migration rate ratios")
    }
    ratios <- rbind(ratios_tmp1, ratios_tmp2)
    # if it is not the first round extend the dataframes rather than create a new one
  }else{
    
    if (splitted[[1]][2]=="coal"){
      ratios_tmp1 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(c_ratio_masco,0.5),
                                upper_masco=quantile(c_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(c_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(c_ratio_sisco,0.5),
                                upper_sisco=quantile(c_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(c_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric coalescent rates symmetric migration rates",
                                facet_posy="coalescent rate ratios")
      ratios_tmp2 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(m_ratio_masco,0.5),
                                upper_masco=quantile(m_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(m_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(m_ratio_sisco,0.5),
                                upper_sisco=quantile(m_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(m_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric coalescent rates symmetric migration rates",
                                facet_posy="migration rate ratios")
    }else{
      ratios_tmp1 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(c_ratio_masco,0.5),
                                upper_masco=quantile(c_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(c_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(c_ratio_sisco,0.5),
                                upper_sisco=quantile(c_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(c_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric migration rates symmetric coalescent rates",
                                facet_posy="coalescent rate ratios")
      ratios_tmp2 <- data.frame(asymmetry=as.numeric(splitted[[1]][3]),
                                mean_masco=quantile(m_ratio_masco,0.5),
                                upper_masco=quantile(m_ratio_masco,0.975)[[1]],
                                lower_masco=quantile(m_ratio_masco,0.025)[[1]],
                                mean_sisco=quantile(m_ratio_sisco,0.5),
                                upper_sisco=quantile(m_ratio_sisco,0.975)[[1]],
                                lower_sisco=quantile(m_ratio_sisco,0.025)[[1]],
                                facet_posx="asymmetric migration rates symmetric coalescent rates",
                                facet_posy="migration rate ratios")
    }
    new.ratios <- rbind(ratios_tmp1, ratios_tmp2)
    # add new mean values
    ratios <- rbind(ratios, new.ratios)
  }
  
  
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


p <- ggplot()+
  geom_point(data=ratios, aes(x=asymmetry,y=mean_masco,color="MASCO"), size=0.5) +
  geom_point(data=ratios, aes(x=asymmetry,y=mean_sisco,color="SISCO"), size=0.5) +
  facet_grid(facet_posy ~ facet_posx) + scale_y_log10() + scale_x_log10() +
  ylab("rate ratio") + xlab("asymmetry") +
  scale_colour_manual("",values = c("MASCO" = col4, "SISCO" = col2)) + guides(colour = guide_legend(override.aes = list(size=2)))

# make the data frame for the red dotted line
d <- data.frame(x=c(0.01, 1, 0.01, 1, 0.01, 1, 0.01, 1),
                y=c(0.01, 1, 1, 1, 1, 1, 0.01, 1),
                facet_posx=c("asymmetric coalescent rates symmetric migration rates",
                             "asymmetric coalescent rates symmetric migration rates",
                             "asymmetric migration rates symmetric coalescent rates",
                             "asymmetric migration rates symmetric coalescent rates",
                             "asymmetric coalescent rates symmetric migration rates",
                             "asymmetric coalescent rates symmetric migration rates",
                             "asymmetric migration rates symmetric coalescent rates",
                             "asymmetric migration rates symmetric coalescent rates"),
                facet_posy=c("coalescent rate ratios","coalescent rate ratios",
                             "coalescent rate ratios","coalescent rate ratios",
                             "migration rate ratios","migration rate ratios",
                             "migration rate ratios","migration rate ratios"))

# plot the red dotten line
p <- p + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed") +
  facet_grid(facet_posy ~ facet_posx)
plot(p)
ggsave(plot=p,"../../text/figures/Asymmetry/Asymmetry.eps",width=10, height=6)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# plot the rate ratios with confidence intervals as in indivual
# figures for both MASCO and SISCO
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


p_masco <- ggplot()+
  geom_pointrange(data=ratios, aes(x=asymmetry, y=mean_masco, ymin=lower_masco, ymax=upper_masco), color=col4, size=0.01) +
  geom_point(data=ratios, aes(x=asymmetry,y=mean_masco), color="black", size=0.01) +
  facet_grid(facet_posy ~ facet_posx) + scale_y_log10() + scale_x_log10() +
  ylab("rate ratio") + xlab("asymmetry")

p_sisco <- ggplot()+
  geom_pointrange(data=ratios, aes(x=asymmetry, y=mean_sisco, ymin=lower_sisco, ymax=upper_sisco), color=col2, size=0.01) +
  geom_point(data=ratios, aes(x=asymmetry,y=mean_sisco), color="black", size=0.01) +
  facet_grid(facet_posy ~ facet_posx) + scale_y_log10() + scale_x_log10() +
  ylab("rate ratio") + xlab("asymmetry")

# make the data frame for the red dotted line
d <- data.frame(x=c(0.01, 1, 0.01, 1, 0.01, 1, 0.01, 1),
                y=c(0.01, 1, 1, 1, 1, 1, 0.01, 1),
                facet_posx=c("asymmetric coalescent rates symmetric migration rates",
                             "asymmetric coalescent rates symmetric migration rates",
                             "asymmetric migration rates symmetric coalescent rates",
                             "asymmetric migration rates symmetric coalescent rates",
                             "asymmetric coalescent rates symmetric migration rates",
                             "asymmetric coalescent rates symmetric migration rates",
                             "asymmetric migration rates symmetric coalescent rates",
                             "asymmetric migration rates symmetric coalescent rates"),
                facet_posy=c("coalescent rate ratios","coalescent rate ratios",
                             "coalescent rate ratios","coalescent rate ratios",
                             "migration rate ratios","migration rate ratios",
                             "migration rate ratios","migration rate ratios"))

# plot the red dotten line
p_masco <- p_masco + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed") +
  facet_grid(facet_posy ~ facet_posx)
p_sisco <- p_sisco + geom_line(data=d,aes(x=x,y=y), color="red",linetype="dashed") +
  facet_grid(facet_posy ~ facet_posx)
plot(p_masco)
plot(p_sisco)
ggsave(plot=p_masco,"../../text/figures/Asymmetry/AsymmetryConfidenceLisco.eps",width=10, height=6)
ggsave(plot=p_sisco,"../../text/figures/Asymmetry/AsymmetryConfidenceSisco.eps",width=10, height=6)



