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

#Load in Data
pathogen_data = read.csv("C:/Users/prisc/Desktop/Boots Lab/Data/final data sets/Final_pathogen.csv")
trans = read.csv("C:/Users/prisc/Desktop/Boots Lab/Data/final data sets/updated_pathogens.csv") %>% select("Species", "transmissibility")
colnames(pathogen_data)[2] = "Species"
pathogen_data$transmissibility = NULL
pathogen_data = merge(pathogen_data, trans, by="Species")

##SUBSET DATA BY UNIQUE DISEASE, CFR, RESERVOIR
pathogen_data <- pathogen_data %>% distinct(Disease, CFR_avg, Reservoir_order, spill_type, .keep_all = TRUE)

#Rename Column Names
colnames(pathogen_data)[26] = "transmissibility"
colnames(pathogen_data)[15] = "transmission_route"

#Remove unnecessary columns
pathogen_data$Reservoir_spill_notes = NULL
pathogen_data$spill_reference = NULL
pathogen_data$CFR_low..decimals. = NULL
pathogen_data$CFR_high = NULL
pathogen_data$CFR_low = NULL
pathogen_data$CFR_notes = NULL
pathogen_data$CFR_reference = NULL
pathogen_data$transmissibility_notes = NULL
pathogen_data$transmissibility_reference = NULL
pathogen_data$transmission_route_notes = NULL
pathogen_data$transmission_route_reference = NULL
pathogen_data$cell_notes = NULL
pathogen_data$cell_reference = NULL
pathogen_data$Reservoir_spill_reference = NULL
pathogen_data$Synonym = NULL

#Change categorical columns to factor columns
pathogen_data$Family = as.factor(pathogen_data$Family)
pathogen_data$spill_type = as.factor(pathogen_data$spill_type)
pathogen_data$transmission_route = as.factor(pathogen_data$transmission_route)
pathogen_data$Cell = as.factor(pathogen_data$Cell)
pathogen_data$Type = as.factor(pathogen_data$Type)
pathogen_data$Motility = as.factor(pathogen_data$Motility)
pathogen_data$Spore = as.factor(pathogen_data$Spore)
pathogen_data$Type = as.factor(pathogen_data$Type)
pathogen_data$CFR_avg = as.numeric(pathogen_data$CFR_avg)
pathogen_data$transmissibility <- as.numeric(pathogen_data$transmissibility)

#RECLASSIFYING TRANSMISSION ROUTES

#Change "spores_in_environment" to "via_environment" transmission route
pathogen_data$transmission_route[pathogen_data$transmission_route == "spores_in_environment"] <- "via_environment"

#Change "fecal_oral" to "via_environment"
pathogen_data$transmission_route[pathogen_data$transmission_route == "fecal_oral"] <- "via_environment"

#RECLASSIFYING Reservoir Orders

#Change Artiodactyla to Cetartiodactyla for Reservoir Order
pathogen_data$Reservoir_order[pathogen_data$Reservoir_order == "Artiodactyla"] <- "Cetartiodactyla"

#Change Reservoir_Order to factor type
pathogen_data$Reservoir_order = toupper(pathogen_data$Reservoir_order)
pathogen_data$Reservoir_order = as.factor(pathogen_data$Reservoir_order)

################################################################################
#REMOVING UNKNOWNS
#Remove unknown and NA cell types
pathogen_data = pathogen_data[!(pathogen_data$Cell == "unknown"),]
pathogen_data = pathogen_data[!is.na(pathogen_data$Cell),]


#Remove unknown and NA reservoir
pathogen_data = pathogen_data[!(pathogen_data$Reservoir_order == "UNKNOWN"),]
pathogen_data = pathogen_data[!is.na(pathogen_data$Reservoir_order),]


#Remove unknown and NA transmission route
pathogen_data = pathogen_data[!(pathogen_data$transmission_route == "unknown"),]
pathogen_data = pathogen_data[!is.na(pathogen_data$transmission_route),]

#Remove rows with NA values for transmissibility column
pathogen_data = pathogen_data[!is.na(pathogen_data$transmissibility),]


#DATA CLEANING IS OVER
###############################################################################

#GAMS for all pathogens using automatic double penalty smoothing
gam_trans_global <- gam(transmissibility ~
                          s(Reservoir_order, bs="re") +
                          s(Family, bs = 're') +
                          s(spill_type, bs="re") +
                          s(CFR_avg, bs="tp") +
                          s(Cell, bs = 're' ),
                          data = pathogen_data,
                          select=TRUE,
                          method = 'REML',
                        family = ocat(theta = c(1,2,3,4)))

summary(gam_trans_global)


effects_dat <- plot.gam(gam_trans_global, residuals=TRUE, select = 5)[[1]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)


################################################################################
################################################################################

#### model selection ####
Terms <- list(
  Family = c(NA, "s(Family, bs = 're' )"),
  Reservoir_order = c(NA, "s(Reservoir_order, bs = 're' )"),
  spill_type = c(NA, "s(spill_type, bs='re')"),
  CFR_avg = c(NA,  "s(CFR_avg, bs = 'tp' )"), 
  Cell = c(NA, "s(Cell, bs='re')"))

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)

## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('transmissibility ~', Formula))

## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        data = pathogen_data,
                        select=FALSE,
                        family = ocat(theta = c(1,2,3,4)))))

removeRows <- lapply(CompetingFits$ModelFit, function(x) 'try-error' %in% class(x)) %>% 
  unlist()
FailedFormulas <- CompetingFits[removeRows, ]
stopifnot(nrow(FailedFormulas) == 0)

CompetingFits <- CompetingFits[!removeRows, ]

## Add AIC:
RankedModels <- CompetingFits %>% 
  mutate(AIC = AIC(ModelFit),
         DevianceExplained = summary(ModelFit)$dev.expl) %>% 
  ungroup() %>% 
  arrange(AIC) %>% 
  mutate(DeltaAIC = AIC - AIC[1])

RankedModels$Formula[which.min(RankedModels$AIC)]

gam_tran <- gam(as.formula(RankedModels$Formula[1]),
                data = pathogen_data,
                select = FALSE,
                family = ocat(theta = c(1,2,3,4)))
summary(gam_tran)
#gam.check(gam_tran)
AIC(gam_tran)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## ---- Plot models rankings: -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Find models to display:
#  - Top 15
library(tidyr)

TERM_NAMES <- tribble(~ Term,                   ~PlotLabel,             ~TermCategory,
                      'Family',          'Family',  NA,
                      "Reservoir_order",    'Reservoir_order',  NA,
                      "spill_type",       'spill_type', NA,
                      "CFR_avg", 'CFR_avg', NA,
                      "Cell",           'Cell', NA)




displayModels <- c(1:6, 8:11) 

# Get deviance explained, etc for these models:
rankedSummary <- summarise_ranked_models_tran(RankedModels[displayModels, ], cores = 8)

# Add ranks:
ranks <- RankedModels %>% 
  select(Formula) %>% 
  mutate(Rank = row_number())

rankedSummary <- rankedSummary %>% 
  select(-Rank) %>% 
  left_join(ranks, by = 'Formula')


# Expand summary to include missing terms in each model:
# - full join ensures rows are created for missing terms
termDf <- lapply(displayModels, function(x) data.frame(Rank = x, TERM_NAMES)) %>% 
  bind_rows()

rankedSummary <- rankedSummary %>% 
  full_join(termDf, by = c('termClean' = 'Term',
                           'Rank' = 'Rank'))

# Calculate model labels:
AICs <- rankedSummary %>% 
  distinct(Rank, DeltaAIC, DevianceExplained) %>% 
  na.omit() %>% 
  mutate(DeltaAIC = sprintf("%.2f", DeltaAIC),
         DeltaAIC = str_pad(DeltaAIC, width = 5),
         DevianceExplained = sprintf("%.1f", DevianceExplained*100),
         Label = paste0(' ', DeltaAIC, '       ', DevianceExplained, '%'))

AICvals <- AICs$Label
names(AICvals) <- AICs$Rank

# Plot:
rankedSummary$sig[is.na(rankedSummary$sig)==TRUE] <- "not_app"

tran_models <- ggplot(rankedSummary, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
  geom_tile(aes(colour = sig, width=0.7, height=0.7), size = 1) +
  scale_color_manual(values = c(N = 'grey60', Y = 'black', not_app = 'grey60'), guide = F) +
  scale_y_continuous(breaks = displayModels,
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = displayModels,
                                         labels = AICvals),
                     trans = 'reverse', expand = expansion(add = c(0.1, 0.1))) +
  scale_x_discrete(expand = c(0.055, 0.055)) +
  scale_fill_viridis_c(option = 'viridis', direction = -1, na.value = 'white',
                       breaks = c(1, seq(5, 50, by = 5))) +
  labs(x = NULL, y = NULL, fill = 'Relative deviance\nexplained (%)') + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.ticks.y = element_blank(),
        legend.position = 'top', 
        legend.title = element_text(vjust = 0.8, hjust = 1),
        legend.key.height = unit(0.8, 'lines'),
        legend.key.width = unit(1.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0))


print(tran_models)

#########################################################################
#########################################################################
#########################################################################
#########################################################################

##############################################################################

#VISUALIZATIONS FOR MODEL SELECTION OF ALL PATHOGENS
## Reservoir Order 
plotData <- get_partial_effects(gam_tran,
                                var = 'Reservoir_order')

reservoir_order <- ggplot(plotData$effects, aes(x = Reservoir_order, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(reservoir_order)


## Spill Type Visualization 
plotData <- get_partial_effects(gam_tran,
                                var = 'spill_type')

spill_type <- ggplot(plotData$effects, aes(x = spill_type, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(spill_type)

