rm(list = ls())

#Load in libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(tidyr)

#Import Data
source("C:/Users/prisc/Desktop/Boots Lab/Code/Code From Papers/Mollentze_functions.R")
pathogen_data = read.csv("C:/Users/prisc/Desktop/Boots Lab/Data/final data sets/final_transm_pathogendata.csv")

#Clean up data
pathogen_data$X = NULL
pathogen_data$Family <- as.factor(pathogen_data$Family)
pathogen_data$Reservoir_order <- as.factor(pathogen_data$Reservoir_order)
pathogen_data$spill_type <- as.factor(pathogen_data$spill_type)
pathogen_data$transmissibility <- as.factor(pathogen_data$transmissibility)
pathogen_data$Cell <- as.factor(pathogen_data$Cell)
pathogen_data$transmission_route <- as.factor(pathogen_data$transmission_route)
pathogen_data$Type <- as.factor(pathogen_data$Type)
pathogen_data$CFR_avg <- as.numeric(pathogen_data$CFR_avg)

#Generate GAM using automated term selection by double penalty smoothing
gam_cfr <- gam(CFR_avg ~
                 s(Reservoir_order, bs="re")+
                 s(spill_type, bs="re")+
                 s(transmissibility, bs="re") +
                 s(Cell, bs="re")+
                 s(Family, bs = 're'),
               data = pathogen_data,
               method = 'REML',
               select= TRUE,
               family = betar(link = "logit"))

#summary of the GAM
summary(gam_cfr)

#Determine the effect size of cell life cycle on human CFR
gam.check(gam_cfr)
effects_dat <- plot.gam(gam_cfr, residuals=TRUE, select = 5)[[4]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)


#################################################################################
#MODEL SELECTION using AIC to select the top model

Terms <- list(
  Family = c(NA, "s(Family, bs = 're' )"),
  Reservoir_order = c(NA, "s(Reservoir_order, bs='re')"),
  spill_type = c(NA, "s(spill_type, bs = 're' )"),
  transmissibility = c(NA, "s(transmissibility, bs = 're' )"),#
  Cell = c(NA, "s(Cell, bs='re')"))


## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)


## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('CFR_avg ~', Formula))


## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        family = betar(link = "logit"),
                        data = as.data.frame(pathogen_data), 
                        select = FALSE,    # Not using internal model selection - mixing selection strategies makes presentation akward
                        method = 'REML')))

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

gam_mort <- gam(as.formula(RankedModels$Formula[1]),
                data = pathogen_data,
                family = betar(link = "logit"),
                method = "REML")
summary(gam_mort)
gam.check(gam_mort)
AIC(gam_mort)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[3]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)


#selected model with lowest AIC
gam_mort <- gam(as.formula(RankedModels$Formula[1]),
                data = pathogen_data,
                select=FALSE,
                family = betar(link = "logit"))
summary(gam_mort)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Find models to display:
#  - Top 10  models predicting pathogen characteristics predictive of human CFR


TERM_NAMES <- tribble(~ Term,                   ~PlotLabel,             ~TermCategory,
                      "Family", "Family", "",
                      "Reservoir_order", "Reservoir_order", "",
                      "spill_type", "spill_type", "",
                      "transmissibility", "transmissibility", "",
                      "Cell", "Cell", "")

#Plot models 1-10
displayModels <- c(1:10) 


# Get deviance explained, etc for these models:
rankedSummary <- summarise_ranked_models_mort(RankedModels[displayModels, ], cores = 8)


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

mort_models <- ggplot(rankedSummary, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
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

print(mort_models)


#################################################################################
#VISUALIZATIONS
#CELL EFFECT SIZES 
plotData <- get_partial_effects(gam_mort,
                                var = 'Cell')

Cell <- ggplot(plotData$effects, aes(x = Cell, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = "none") +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = "none") +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on CFR in humans') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(Cell)


#VISUALIZATIONS
##RESERVOIR ORDER EFFECT SIZES 

plotData <- get_partial_effects(gam_mort,
                                var = 'Reservoir_order')

Reservoir_order <- ggplot(plotData$effects, aes(x = Reservoir_order, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#7a378b'), guide = "none") +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#7a378b'), guide = "none") +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on CFR in humans') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(Reservoir_order)


#VISUALIZATIONS
##Spill Type EFFECT SIZES 

plotData <- get_partial_effects(gam_mort,
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
  scale_colour_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = "none") +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#85c874'), guide = "none") +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on CFR in humans') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(spill_type)


