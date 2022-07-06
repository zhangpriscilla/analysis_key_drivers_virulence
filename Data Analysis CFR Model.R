rm(list = ls())

library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(cowplot)
library(MuMIn)
library(corrplot)
library(lme4)
library(scatterpie)
library(reshape2)
library(zoib)
library(ggpubr)
library(tidyr)
library(ggnewscale)
library(ggplot2)
library(cowplot)
library(randomForest)
library(data.table)

library(rpart)				    
library(rpart.plot)	
library(RColorBrewer)
library(caret)
library(ROCR)
library(randomForest)
library(missForest)
library(pdp)
library(pbapply)
library(forcats)

source("C:/Users/prisc/Desktop/Boots Lab/Code/Code From Papers/Mollentze_functions.R")

pathogen_data = read.csv("C:/Users/prisc/Desktop/Boots Lab/Data/final data sets/final_transm_pathogendata.csv")

pathogen_data$Family <- as.factor(pathogen_data$Family)
pathogen_data$Reservoir_order <- as.factor(pathogen_data$Reservoir_order)
pathogen_data$spill_type <- as.factor(pathogen_data$spill_type)
pathogen_data$transmissibility <- as.factor(pathogen_data$transmissibility)
pathogen_data$Cell <- as.factor(pathogen_data$Cell)
pathogen_data$CFR_avg <- as.numeric(pathogen_data$CFR_avg)


gam_cfr <- gam(CFR_avg ~
                  s(Reservoir_order, bs="re")+
                  s(spill_type, bs="re")+
                  s(transmissibility, bs="re") +
                  s(Cell, bs="re")+
                  s(Family, bs = 're'),
                  data = pathogen_data,
                  select= TRUE,
                  family = betar(link = "logit"))

summary(gam_cfr)

gam.check(gam_cfr)
effects_dat <- plot.gam(gam_cfr, residuals=TRUE, select = 5)[[1]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)

#Transmission route and Reservoir Order are significant predictors of the model
#VISUALIZATIONS

#CFRavg ~ cell, colored by spill type
ggplot(data=pathogen_data, aes(x=Cell, y=CFR_avg)) +
  geom_boxplot() +
  geom_point(aes(color=spill_type)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Average human case fataltiy rate") + xlab("Life cycle types")


ggplot(data=pathogen_data, aes(x=transmission_route, y=CFR_avg)) +
  geom_boxplot() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(data=subset(pathogen_data, transmission_route == "via_environment" & Reservoir_order == "ENVIRONMENT"), size=2.3, aes(label=Disease))



# #CFRavg ~ cell, colored by transmissibility
# ggplot(data=pathogen_data, aes(x=Cell, y=CFR_avg, label=Disease)) +
#   geom_boxplot() +
#   geom_point(aes(color=transmissibility)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + 
#   geom_text(data=subset(pathogen_data, CFR_avg > 0.25 & Cell == "Extracellular"), size=2.3, aes(label=Disease)) + 
#   geom_text(data=subset(pathogen_data, CFR_avg > 0.3 & Cell == "Facultative intracellular"), size=2.3, aes(label=Disease))

###############################################################################
#Obligate Intracellular pathogens
intracellular = pathogen_data[pathogen_data$Cell == "Obligate intracellular",]

intracellular_cfr <- gam(CFR_avg ~ 
                         s(Family, bs = 're') +
                         s(Reservoir_order, bs="re") +
                         s(spill_type, bs="re") +
                         #s(transmission_route, bs="re") +
                         s(transmissibility, bs="re") +
                         s(Type, bs="re"),
                         data = intracellular,
                         select=TRUE,
                         family = betar(link = "logit"))

summary(intracellular_cfr)

effects_dat <- plot.gam(intracellular_cfr, residuals=TRUE, select = 5)[[4]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)

###############################################################################
#Bacteria ONLY - nothing significant in the model

bacteria = pathogen_data[pathogen_data$Type == "Bacteria",]

bacteria_cfr <- gam(CFR_avg ~ 
                           s(Family, bs = 're') +
                           s(Reservoir_order, bs="re") +
                           s(spill_type, bs="re") +
                           #s(transmission_route, bs="re") +
                           s(transmissibility, bs="re"),
                           data = bacteria,
                           select=TRUE,
                           family = betar(link = "logit"))

summary(bacteria_cfr)

################################################################################
#Virus ONLY

virus = pathogen_data[pathogen_data$Type == "Virus",]


#######################
# #GET Phylogenetic distance for reservoir order
#
# #Read in Stringent Data
stringent_data = read.csv("C:/Users/prisc/Desktop/Boots Lab/Data/stringent_data.csv")

#Grab "h0rder" and "phylo_dist" columns
phylo = stringent_data[c("hOrder","phylo_dist")]

#Unique h0rder
phylo = phylo %>% distinct(hOrder,.keep_all = TRUE)

#Add in other Reservoirs and their distance
phylo[nrow(phylo) + 1,] = c("Environment", Inf)
phylo[nrow(phylo) + 1,] = c("Unknown", Inf)
phylo[nrow(phylo) + 1,] = c("Fish", 864)
phylo[nrow(phylo) + 1,] = c("Lagomorpha", 180)
phylo[nrow(phylo) + 1,] = c("CINGULATA", 210)
phylo[nrow(phylo) + 1,] = c("Vector", 1594)

#Remove EULIPOTYPHLA
phylo = phylo[!(phylo$hOrder == "EULIPOTYPHLA"),]

#Upper Case all values in h0rder
phylo$hOrder = toupper(phylo$hOrder)

#Ensure that phylo_dist is numeric
phylo$phylo_dist <- as.numeric(phylo$phylo_dist)

#Sort phylo distance from lowest to highest
reservoir_distances = phylo[order(phylo$phylo_dist),]$hOrder
reservoir_distances


#CFRavg ~ Reservoir by phylogenetic distance to PRIMATES
pathogen_data %>%
  mutate(Reservoir_order = fct_relevel(Reservoir_order, reservoir_distances)) %>%
  ggplot(aes(x=Reservoir_order, y=CFR_avg)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + geom_point(aes())