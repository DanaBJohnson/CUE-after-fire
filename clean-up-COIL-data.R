# --------------------
# Title: Clean up C isotope data from COIL
# Author: Dana Johnson
# Date: 2024-February-28
#
# OBJECTIVE: 



library(tidyr)
library(ggplot2)


# setwd('../Box/WhitmanLab/Projects/WoodBuffalo/ModellingFires_DOE_TES/data/CUE')

# Soil horizon labels and colors:
Horizon_labels = c('O horizon','Mineral soil')
names(Horizon_labels) = c('O','M')

Horizon_colors = c('green','black')
names(Horizon_colors) = c('O','M')

# Site type labels and colors:
SiteType_labels = c('Organic sites','Sandy sites')
names(SiteType_labels) = c('ORG','SANDY')

SiteType_colors = c('green','purple')
names(SiteType_colors) = c('ORG','SANDY')




df.Pine <- read.csv('COIL/raw-data-Jan2024/Stable_isotope_data_from_combined_runs_Jan2024.csv')
df.Glucose <- read.csv('COIL/raw-data-Oct2023/Stable_isotope_data_from_combined_runs_Oct2023.csv')

head(df.Pine)
head(df.Glucose)

df <- rbind(df.Pine, df.Glucose)

head(df)

# clean up columns names
colnames(df) <-  c('sample.ID', 
                   'weight.mg',
                   'N2.Amp',
                   'percent.N',
                   'atom.percent.15N.14N', 
                   'delta.15N.vs.at.air',
                   'CO2.Amp',
                   'percent.C',
                   'delta.13C.vs.VPDB',
                   'atom.percent.13C.12C')

# separate the sample ID to create more meaningful columns
df <- df %>%
  separate(sample.ID, c("core.hor.id",'amendment','fumigation'), sep = "_", remove = FALSE) 

# create a dataframe of blanks
df.blank <- df[grep('^[B]', df$sample.ID),]


# create a dataframe for samples
df.samples <- df[grep('^[0-9]', df$sample.ID),]



#### Problem: Some of these samples are reruns/duplicates and that messes up the 
#   naming scheme.

# To fix this, Create new column for this information:
for (i in 1:nrow(df.samples)) {
  if (df.samples$amendment[i] != 'rerun') {
    df.samples$rerun[i] = 'no'
  } else if (df.samples$amendment[i] == 'rerun') {
    df.samples$rerun[i] = 'yes'
  } 
  
  if (df.samples$amendment[i] != 'duplicate') {
    df.samples$duplicate[i] = 'no'
  } else if (df.samples$amendment[i] == 'duplicate') {
    df.samples$duplicate[i] = 'yes'
  }
}


# Separate reruns and re-split sample ID
df.reruns <- subset(df.samples,rerun == 'yes') %>%
  separate(sample.ID, c("core.hor.id",'rerun','amendment','fumigation'), sep = "_", remove = FALSE) 

# Separate duplicates and re-split sample ID
df.duplicates <- subset(df.samples, duplicate == 'yes') %>%
  separate(sample.ID, c("core.hor.id",'duplicate','amendment','fumigation'), sep = "_", remove = FALSE) 

# Merge everything back together
df.samples.full <- df.samples %>%
  subset(rerun == 'no' & duplicate == 'no') %>%
  rbind(df.reruns) %>%
  rbind(df.duplicates) %>%
  separate(core.hor.id, c('site','core','horizon'), sep = c(2,4), remove = FALSE)


# Assign soils to site types
for (i in 1:nrow(df.samples.full)) {
  if (df.samples.full$site[i] %in% c('01','02','03','05','08','09')) {
    df.samples.full$site.type[i] = 'SANDY'
  } else {
    df.samples.full$site.type[i] = 'ORG'
  }
}

# Assign substrate to each sample
for (i in 1:nrow(df.samples.full)) {
  if (df.samples.full$amendment[i] %in% c('12CG','13CG')) {
    df.samples.full$substrate[i] = 'glucose'
  } else if(df.samples.full$amendment[i] %in% c('12CP','13CP')) {
    df.samples.full$substrate[i] = 'pine'
  } else {
    df.samples.full$substrate[i] = 'blank'
  }
}

# write.csv(df.samples.full, "merged_CFE_C_concentration_and_atom-percent_data.csv")