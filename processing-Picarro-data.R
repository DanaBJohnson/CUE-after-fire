##title: "Processing CUE gas data"
##author: "Dana Johnson"
##date: "06/27/2023"

##objective: Collect all the Picarro %C and d13C data from the CUE incubations, 
#   and compile into one place. Make adjustments for blanks and standards. 
#   Calculate g CO2 respired and d13C of C respired for each sample. 
#   !! Want ug 13CO2-C per gram soil for CUE calculations !!

library(ggplot2)
library(dplyr)
library(tidyr)

# setwd('../Box/WhitmanLab/Projects/WoodBuffalo/ModellingFires_DOE_TES/code')

# Import raw data -----------
df.24 <- read.csv("../data/CUE/Picarro/Compiled_data_24HR.csv")
df.48 <- read.csv("../data/CUE/Picarro/Compiled_data_48HR.csv")
df.72 <- read.csv("../data/CUE/Picarro/Compiled_data_72HR.csv")
df.96 <- read.csv("../data/CUE/Picarro/Compiled_data_96HR.csv")
df.std <- read.csv("../data/CUE/Picarro/Compiled_data_STANDARDS.csv")



### Processing Picarro data: 21-March-2024 -----------

# First look at measured vs target values: 
#   Plot target CO2 concentration vs. measured CO2 concentration
p1 <- ggplot(df.std, aes(x=target_CO2_conc_ppm, y = X12CO2_mean + X13CO2_mean, color = sample_run)) + 
  geom_point(size=2) + 
  theme_bw() +
  # geom_abline(slope=1, intercept=0, color = 'grey70')+
  geom_smooth(method = 'glm', se = FALSE, color = 'grey10')

#   Plot target d13C vs. measured d13C
p2 <- ggplot(df.std, aes(x=target_d13C, y = Delta_30s_iCO2_mean, color = sample_run)) + 
  geom_point(size=2, alpha=0.7)+
  theme_bw() +
  geom_smooth(method='glm', se=FALSE, color = 'grey10')

cowplot::plot_grid(p1, p2, ncol=2, 
                   labels = 'AUTO',
                   vjust=-5)

# Calculate difference btw measured std and expected std CO2 concentrations and isotopic compositions -----------
df.std <- subset(df.std, sample_run != 'troubleshooting leak') %>%
  mutate(delta_CO2_ppm = target_CO2_conc_ppm - (X12CO2_mean + X13CO2_mean),
         delta_iCO2 = target_d13C - Delta_30s_iCO2_mean)

# Calculate the target 12CO2 ppm and 13CO2 ppm concentrations:
#   delta 13C = [ (13C/12C)_sample / (13C/12C)_std - 1 ] * 1000
#   (13C/12C)_std = 0.01237 = PBD standard, 0.011100 = VPDB standard
#   target_CO2_conc_ppm = target_12CO2_conc_ppm + target_12CO2_conc_ppm
#   (13C/12C)_sample = target_13CO2_conc_ppm/target_12CO2_conc_ppm

df.std <- df.std %>%
  mutate(y = 0.011100*(target_d13C/1000+1), 
         target_13CO2_conc_ppm = (target_CO2_conc_ppm*y)/(1+y),
         target_12CO2_conc_ppm = target_CO2_conc_ppm-target_13CO2_conc_ppm)

# Plot comparison of barget vs. measured CO2
ggplot(subset(df.std, H2O == 0), 
       aes(x=target_13CO2_conc_ppm, y = X13CO2_mean,
           color = target_d13C)) +
  geom_point(size=2)+
  # facet_wrap(~Rundate, scales = 'free') +
  theme_bw() +
  labs(title = "comparison of target and measured 13CO2 concentrations")+
  geom_smooth(method = 'glm', se = FALSE, color = 'grey')+
  geom_abline(slope=1, intercept = 0, color = 'black', linetype='dotted') 


### Correct for delta 13C and C concentration separately -----------

# Create linear models for d13C, 12C-CO2 conc., and 13C-CO2 conc. 
lm_delta13 <- lm(Delta_30s_iCO2_mean ~ target_d13C, data=df.std) 
lm_conc12 <- lm(X12CO2_mean ~ target_12CO2_conc_ppm , data=df.std) 
lm_conc13 <- lm(X13CO2_mean ~ target_13CO2_conc_ppm , data=df.std) 

# Create linear model for total CO2 conc. -----------
lm_conc <- lm(X12CO2_mean+X13CO2_mean ~ target_CO2_conc_ppm , data=df.std) 

# Extract coefficients from models
conc_intercept <- summary(lm_conc)$coefficients[1]
conc_slope <- summary(lm_conc)$coefficients[2]

conc12_intercept <- summary(lm_conc12)$coefficients[1]
conc12_slope <- summary(lm_conc12)$coefficients[2]

conc13_intercept <- summary(lm_conc13)$coefficients[1]
conc13_slope <- summary(lm_conc13)$coefficients[2]

delta13_intercept  <- summary(lm_delta13)$coefficients[1]
delta13_slope <- summary(lm_delta13)$coefficients[2]



### Combine dataframes, correct for dilutions, correct for blanks -----------

# Remove the methane data for now and combine all the dataframes. 
df.24 <- df.24[, c(1:16, 33:42)]
colnames(df.24)[21:22] <- c('Delta_30s_iCO2_mean', 'Delta_30s_iCO2_stdev')
colnames(df.96)[21:22] <- c('Delta_30s_iCO2_mean', 'Delta_30s_iCO2_stdev')

df <- rbind(df.24, df.48, df.72, df.96) %>%
  subset(sample.id != "")

## Correct concentrations based on linear models
#    y = mx+b --> x = (y-b)/m
df <- df %>%
  mutate(corrected.12CO2.concentration = (X12CO2_mean - conc12_intercept)/conc12_slope,
         corrected.13CO2.concentration = (X13CO2_mean - conc13_intercept)/conc13_slope,
         corrected.total.CO2.concentration = (X12CO2_mean + X13CO2_mean - conc_intercept)/conc_slope,
         # corrected.13CO2.concentration.V2 = (X13CO2_mean - conc_intercept)/conc_slope,
         corrected.delta13C = (Delta_30s_iCO2_mean - delta13_intercept)/delta13_slope)


# Check zero Air readings to see if there are any discrepancies in data here: -----------
#   identify and average zeroair
df.zero <- subset(df, horizon == 'zeroair') %>%
  dplyr::mutate(zero.CO2.total.ppm = (corrected.12CO2.concentration+corrected.13CO2.concentration))

# Isolating zeroair samples from Dec 13-21
for (i in 1:nrow(df.zero)) {
  if (df.zero$Rundate[i] >= 20221213 & df.zero$Rundate[i] <= 20221221) {
    df.zero$timepoint[i] = '24HR'
  }
}

# Plot results:
ggplot(df.zero, aes(x=Rundate+Runtime, y = zero.CO2.total.ppm, color =timepoint)) + 
  geom_point(size=5) + 
  theme_bw()

# Does not look like there are any instances of big leaks leading to high CO2 values in zeroair samples
lm(zero.CO2.total.ppm~(Rundate+Runtime), data=df.zero) %>% summary()



# Correct for dilutions  -----------
#   All samples were diluted 1:1 with zeroAir:sample at initial collection
#   Samples were further diluted for analysis on Picarro if they were organic-rich
for (i in 1:nrow(df)) {
  if (df$dilution.for.analysis.on.picarro..mL.zero..mL.sample.[i]=='') {
    df$initial.12CO2_C_mean_ppm[i] = df$corrected.12CO2.concentration[i] *2
    df$initial.13CO2_C_mean_ppm[i] = df$corrected.13CO2.concentration[i] *2
  } else if (df$dilution.for.analysis.on.picarro..mL.zero..mL.sample.[i]=='1:1') {
    df$initial.12CO2_C_mean_ppm[i] = df$corrected.12CO2.concentration[i] *4
    df$initial.13CO2_C_mean_ppm[i] = df$corrected.13CO2.concentration[i] *4
  } else if (df$dilution.for.analysis.on.picarro..mL.zero..mL.sample.[i]=='5:1') {
    df$initial.12CO2_C_mean_ppm[i] = df$corrected.12CO2.concentration[i] *12
    df$initial.13CO2_C_mean_ppm[i] = df$corrected.13CO2.concentration[i] *12
  }
}

### Correct for Run "sets"  -----------
# Order dataframe by date and time of Picarro run
df <- df[order(df$Rundate, df$Runtime),]

# The inital run "set" equals 1
set = 1
df$Run.set=1

for (i in 2:nrow(df)) {
  if (df$Rundate[i] == df$Rundate[i-1] & df$Runtime[i]-df$Runtime[i-1] < 20000) {
    df$Run.set[i] = set
  } else if (df$Rundate[i] == df$Rundate[i-1] & df$Runtime[i]-df$Runtime[i-1] > 20000){
    set = set+1
    df$Run.set[i] = set
  } else if (df$Rundate[i] != df$Rundate[i-1] & df$Runtime[i]-df$Runtime[i-1] < -234000) {
    df$Run.set[i] = set
  } else if (df$Rundate[i] != df$Rundate[i-1] & df$Runtime[i]-df$Runtime[i-1] > -234000) {
    set = set + 1
    df$Run.set[i] = set
  }
}


# Check that the groupings make sense - yep looks good. 
ggplot(df, aes(x=as.character(Rundate), y = Runtime, color = as.character(Run.set))) + 
  geom_point()+
  scale_colour_manual(values = c('orangered','blue','green','yellow4','purple',
                                 'gold3','black','brown','magenta','grey30','orange',
                                 'blue4','red4','forestgreen'))


# Identify and average blanks for each "run set" -----------
df.blank <- subset(df, horizon == 'blank') %>%
# Group blank values by Run set and calculate average blank
  dplyr::group_by(Run.set) %>%
  mutate(avg.blank.CO2.total.ppm=mean(initial.12CO2_C_mean_ppm+initial.13CO2_C_mean_ppm),
         avg.blank.d13C = mean(Delta_30s_iCO2_mean),
         avg.blank.13CO2.ppm = mean(initial.13CO2_C_mean_ppm)) %>%
  subset(select = c(Run.set, avg.blank.CO2.total.ppm,avg.blank.d13C, avg.blank.13CO2.ppm)) %>%
  unique()


# sanity check
ggplot(df.blank, aes(x=Run.set, y = avg.blank.13CO2.ppm, color = as.character(Run.set))) + 
  geom_point(size=5) + 
  theme_bw() 


# Subtract average blank from measured CO2 values by grouping by date of run 
df.merge <- df  %>%
  subset(carbon.addition %in% c("12CG","12CP","13CG","13CP") & horizon %in% c('O','M')) %>%
  merge(df.blank, by = c('Run.set'), all = TRUE) %>%
  mutate(combined_CO2_ppm = initial.12CO2_C_mean_ppm + initial.13CO2_C_mean_ppm,
         combined_CO2_ppm_blank_corrected = combined_CO2_ppm-avg.blank.CO2.total.ppm,
         blank_corrected_13_CO2_ppm = initial.13CO2_C_mean_ppm-avg.blank.13CO2.ppm)

# Run.set 3 = only blanks
df.merge <- df.merge %>%
  subset(Run.set != 3)


### Calculate gCO2-C starting with CO2 ppm -----------
#   n = total moles in the gas mixture: n = PV/RT
P = 1 # atm
R = 0.08206 # L*atm/mol*K
Temp = 273+25 # K
V = 0.128 # L
#   moles CO2 = PPM/1,000,000 * n
#   g C-CO2 = moles CO2 * 12.01 g C/mol CO2

df.merge <- df.merge %>%
  mutate(n = P*V/(R*Temp),
         moles_13CO2 = blank_corrected_13_CO2_ppm/1000000*n,
         g_13CO2_C = moles_13CO2*12.01,
         mg_13CO2_C = g_13CO2_C *1000,
         mg_CO2_C_total = combined_CO2_ppm_blank_corrected/1000000*n*12.01*1000)
         # g_CO2_C = combined_CO2_ppm_blank_corrected/1000000 * 22.41 * 1/0.128 * 12.01,
         # mg_CO2_C = g_CO2_C*1000) %>%
#        X mol CO2 = X ppm CO2/million mol gas 
#        g CO2-C = X mol CO2 * 22.41 L/mol * 1/0.128 L (jar headspace) * 12.01 g C


### Import moisture data to make correction of CO2-C respired per gram dry soil -----------
df.moisture <- read.csv('../data/CUE/CUE-dry-soil-mass.csv')

# clean up dataframe
df.moisture <- df.moisture %>%
  mutate(carbon.addition = amendment) %>%
  subset(select = c(site, core, horizon, carbon.addition, dry.soil.mass.g)) 

df.moisture$site <- as.character(df.moisture$site)

# merge moisture data with primary dataframe  
df.merge <- merge(df.merge, df.moisture, 
                  by = c('site', 'core', 'horizon', 'carbon.addition'))

# Calculate mg CO2 respired per grams of dry soil -----------
df.merge <- df.merge %>%
  mutate(mg_13CO2_C_per_g_dry_soil = mg_13CO2_C/dry.soil.mass.g) %>%
  mutate(mg_CO2_C_per_g_dry_soil = mg_CO2_C_total/dry.soil.mass.g)

# Plot final results: -----------
p3 = ggplot(subset(df.merge), 
            aes(x=timepoint, y = mg_13CO2_C_per_g_dry_soil)) + 
    geom_boxplot()+
  geom_jitter(aes(color = site), size=2, width=0.1,alpha=0.7) + 
  theme_bw() + 
  facet_grid(horizon~carbon.addition, scales='free') + 
  labs(title = 'With dilution correction') + 
  # ylim(-0.001, 0.08) + 
  geom_hline(yintercept=0, color = 'grey30')

p3

# Save data frame:  -----------
# write.csv(df.merge, '../data/CUE/Picarro/calculated-mg-C-respired.csv')
