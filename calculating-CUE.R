# --------------------
# Title: Calculating CUE
# Author: Dana Johnson
# Date: 2024-February-28
#
# OBJECTIVE: Using data from incubations and chloroform fumigation extraction,
#   calculate CUE

# setwd('../Box/WhitmanLab/Projects/WoodBuffalo/ModellingFires_DOE_TES/code')

### List of figures: -----------------------------------------------------------

# PLOT: Total MBC separated by amendment
# PLOT: Soil- and substrate-derived MBC, combined pine and glucose
# PLOT: Cumulative soil-derived CO2-C
# PLOT: Cumul. soil-derived CO2-C, combined pine and glucose
# PLOT: Soil-derived MBC, combined pine and glucose
# PLOT: CUE
# PLOT: qCO2


### Set up -----------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)
library(minpack.lm)


# Set colors and labels for plots
Soil_colors = c('green', 'brown')
names(Soil_colors) = c('ORG','SANDY') 

Soil_labels = c('Histosol','Gleysol')
names(Soil_labels) = c('ORG','SANDY')

Burn_colors = c('black', 'darksalmon','darkred')
names(Burn_colors) = c('0','30','120') 

Burn_labels = c('Unburned', '30 second exposure','120 second exposure') 
names(Burn_labels) = c('0','30','120')

Horizon_labels = c('O horizon','Mineral soil')
names(Horizon_labels) = c('O','M')

# Site type labels and colors:
SiteType_labels = c('Histosol','Gleysol')
names(SiteType_labels) = c('ORG','SANDY')

SiteType_colors = c('darkgreen','purple')
names(SiteType_colors) = c('ORG','SANDY')



### Import and clean up CFE and incubation data ----------------------------------
df.CFE = read.csv("../data/CUE/COIL/merged_CFE_C_concentration_and_atom-percent_data.csv")
# Notes: df.CFE contains 461 samples, which includes some duplicates

# import metadata for burn treatments
df.meta <- read.csv('../data/metadata.csv')

# clean up dataframe:
df.CFE = df.CFE[,2:21]
colnames(df.CFE)[1] = 'sample.id'

# fix fumigation label typos:
for (i in 1:nrow(df.CFE)) {
  if (df.CFE$fumigation[i] == 'CHCl') {
    df.CFE$fumigation[i] = 'CHCl3'
  } else {
    df.CFE$fumigation[i] = df.CFE$fumigation[i]
  }
}

# df.CFE contains a lot of duplicates and re-runs. Average across duplicates
df.CFE <- df.CFE %>%
  dplyr::group_by(core.hor.id, site, core, horizon, amendment, fumigation, site.type, substrate) %>%
  summarize(average.atom.percent.13C.12C = mean(atom.percent.13C.12C),
                   average.percent.C.DOC = mean(percent.C)) # C (%) comes from COIL data and is the C concentration of salt extracts 
# Notes: df.CFE now contains 432 samples = 54 core.hor x 4 amendments x 2 fumigation treatments



### Save required constants: ---------------------------------------------------

k_EC = 1 # An extraction coefficient of 0.45 suggests that 55% of the isotopic 
#   label is stuck in an unextractable form. This could be in microbial spores
#   or bacterial cell wall, which way be too large to pass through 0.45 um filter?
#   Setting k_EC = 1 because we are not using an extraction coefficient.

# atom % pine
at_perc_unlabelled_pine = 1.06997       ## IRMS data, January 7, 2025
at_perc_labelled_pine = 1.57612   ## IRMS data, January 7,2025

# delta 13C for pine
d13C_unlab_pine = -26.78 ## IRMS data, January 7,2025
d13C_lab_pine = 440.98 ## IRMS data, January 7,2025

# atom % glucose
at_perc_unlabelled_glucose = 1.08694    ## IRMS data, January 7,2025
at_perc_labelled_glucose = 5.16665      ## calculated from IRMS data and "CUE_labeling_calculations.xlsx" 

# delta 13C for glucose
d13C_unlab_glucose= -11.17 ## IRMS data, January 7, 2025
d13C_lab_glucose = 3902.49 ## calculated from IRMS data and "CUE_labeling_calculations.xlsx" 

# at% for Pinus and Picea soils (average from all pinus or picea, respectively,
#   samples from 70 day incubation using days 3-21)
at_perc_sandy_soil =1.07125
at_perc_org_soil = 1.07248
  
# delta 13C for Pinus and Picea samples (average from all pinus or picea, respectively,
#   samples from 70 day incubation using days 3-21)
d13C_sandy_soil = -25.60482
d13C_org_soil = -24.47303





### Import soil mass data and combine with CFE data. ---------------------------

# Import data: this is a dataframe with the soil mass used in the CUE incubations
df.convert <- read.csv('../data/CUE/CFE-soil-mass.csv')
head(df.convert,2)
### Notes: df.convert contains 216 samples = 54 core.hor x 4 amendments 


# Need to create a sample.id that matches the CFE data (Example: 0912O_13CP_CHCl)
df.convert <- df.convert %>%
  separate(sample.id.short, c('site','core','hor'), sep = "-", remove = TRUE) %>%
  unite(sample.id.short, c(site,core,hor), sep = "")

# reformat amendment and fumigation 
for (i in 1:nrow(df.convert)) {
  if (df.convert$amendment[i] == '12C glucose') {
    df.convert$amendment[i] = '12CG'
  } else if (df.convert$amendment[i] == '13C glucose') {
    df.convert$amendment[i] = '13CG'
  } else if (df.convert$amendment[i] == '12C pine') {
    df.convert$amendment[i] = '12CP'
  } else if (df.convert$amendment[i] == '13C pine') {
    df.convert$amendment[i] = '13CP'
  }
}

# Fix column names to match df.CFE
df.convert <- subset(df.convert, select = c(sample.id.short, amendment, core.id,
                                            Mass.with.chloroform..no.X.,mass.without.chloroform..g...X.)) %>%
  unite(sample.id, c(sample.id.short, amendment))

colnames(df.convert)[3:4] <- c('wet.soil.mass.g.fumigated', 'wet.soil.mass.g.not.fumigated')


# Rework CFE dataframe so each row is a sample and contains both fumigated and non-fumigated data
df.CFE.fumigated = subset(df.CFE, fumigation=='CHCl3')
### Notes: df.CFE.fumigated contains 216 samples = 54 core.hor x 4 amendments x 1 fumigation trtmt


# create sample id that includes amendment
df.CFE.fumigated <- df.CFE.fumigated %>%
  unite(sample.id, c(core.hor.id, amendment), sep = "_", remove = FALSE) %>%
  subset(select = c(sample.id, site.type, amendment, substrate, horizon, average.percent.C.DOC, 
                    average.atom.percent.13C.12C)) %>%
  unique()

colnames(df.CFE.fumigated)[6:7] <- c('percent.C.DOC.fumigated','atom.percent.13C.12C.fumigated')


# repeat for un-fumigated dataset
df.CFE.not.fumigated = subset(df.CFE, fumigation == 'nofum')
### Notes: df.CFE.not.fumigated contains 216 samples = 54 core.hor x 4 amendments x 1 fumigation trtmt

# create sample ID that includes amendment
df.CFE.not.fumigated <- df.CFE.not.fumigated %>%
  unite(sample.id, c(core.hor.id, amendment), sep = "_", remove = FALSE) %>%
  subset(select = c(sample.id, site.type, amendment, substrate, horizon, average.percent.C.DOC, 
                    average.atom.percent.13C.12C)) %>%
  unique()

colnames(df.CFE.not.fumigated)[6:7] <- c('percent.C.DOC.notfumigated','atom.percent.13C.12C.notfumigated')



# combine fumigated and not fumigated dataframes by the sample ID (which includes core.hor.amendment)
df.CFE.full <- merge(df.CFE.fumigated, df.CFE.not.fumigated, 
                     by = c('sample.id', 'site.type', 'amendment', 'horizon', 'substrate'))
# Notes: df.CFE.full contains 216 samples = 54 core.hor x 4 amendments
#           fumigation trtmt is now include as columns instead of rows


# Combine soil mass data (df.fum.combined) with CFE data
df.CFE.merge <- merge(df.convert, df.CFE.full, by = 'sample.id')
#    Notes: df.CFE. merge contains 216 samples = 54 core.hor x 4 amendments



### Import moisture data and normalize respiration per gram dry soil -----------

df.CUE.moisture <- read.csv('../data/CUE/incubation-set-up.csv')
# Notes: df.CUE.moisture contains 216 samples = 54 core.hor x 4 amendments

# Notes:
#   wet soil mass (g) = wet.soil.mass.initial..g. - mass.of.empty.jar.and.lid..g.
#   dry soil mass (g) = wet soil mass (g) / (1 + estimated moisture (%) at end of incubation)


# calculate initial soil moisture:
### Initial soil moisture = water / dry soil 
### soil moisture for each incubation = (initial water + mL amendment solution added) / dry soil 
df.CUE.moisture <- df.CUE.moisture %>%
  dplyr::mutate(initial.soil.moisture = (wet.soil.mass.initial..g. - mass.of.empty.jar.and.lid..g. - dry.soil.mass.initial..g.)/dry.soil.mass.initial..g.) %>%
  dplyr::mutate(CG.incub.soil.moisture = (wet.soil.mass.initial..g.+ mL.glucose.solution.to.add- mass.of.empty.jar.and.lid..g. - dry.soil.mass.initial..g.)/dry.soil.mass.initial..g.,
                CP.incub.soil.moisture = (wet.soil.mass.initial..g.+ mL.pine.solution.to.add- mass.of.empty.jar.and.lid..g. - dry.soil.mass.initial..g.)/dry.soil.mass.initial..g.)

# Create dataframe subset with just glucose samples
df.CG.moisture <- df.CUE.moisture %>%
  subset(amendment %in% c('12C glucose', '13C glucose')) %>%
  subset(select = c(core.id, sample.id.short, horizon, amendment, dry.soil.mass.initial..g., CG.incub.soil.moisture, initial.soil.moisture)) 

colnames(df.CG.moisture)[5:6] <- c('dry.soil.mass.g','incub.soil.moisture')

# Create dataframe subset with just pine samples
df.CP.moisture <- df.CUE.moisture %>%
  subset(amendment %in% c('12C pine', '13C pine')) %>%
  subset(select = c(core.id, sample.id.short, horizon, amendment, dry.soil.mass.initial..g., CP.incub.soil.moisture, initial.soil.moisture)) 

colnames(df.CP.moisture)[5:6] <- c('dry.soil.mass.g','incub.soil.moisture')

# Recombine glucose and pine dataframes with new column headings
df.moisture <- rbind(df.CG.moisture, df.CP.moisture)

head(df.moisture)

# Reformat sample IDs
df.moisture <- df.moisture %>%
  separate(sample.id.short, c('site', 'core','hor'), sep = '-', remove = TRUE) %>%
  unite(sample.id, c(site, core, hor), sep = '', remove = FALSE)
  
# Reformat amendment column and create substrate column for later merging
for (i in 1:nrow(df.moisture)) {
  if (df.moisture$amendment[i] == '12C glucose') {
    df.moisture$amendment[i] = '12CG'
    df.moisture$substrate[i] = 'glucose'
  } else if (df.moisture$amendment[i] == '13C glucose') {
    df.moisture$amendment[i] = '13CG'
    df.moisture$substrate[i] = 'glucose'
  } else if (df.moisture$amendment[i] == '12C pine') {
    df.moisture$amendment[i] = '12CP'    
    df.moisture$substrate[i] = 'pine'
  } else if (df.moisture$amendment[i] == '13C pine') {
    df.moisture$amendment[i] = '13CP'
    df.moisture$substrate[i] = 'pine'
  }
}

# write.csv(df.moisture, '../data/CUE/CUE-dry-soil-mass.csv')


### Convert percent C (g C/g salt) data to ug C / g dry soil -------------------

# Now I have the moisture data for the start of the incubation (incub.soil.moisture)
colnames(df.CFE.merge)[1] <- 'sample.id.amend'
df.CFE.merge <- merge(df.CFE.merge, df.moisture, by = c('core.id', 'horizon', 'amendment', 'substrate'))
# Notes: df.CFE.merge contains 216 samples = 54 core.hor x 4 amendments


# Calculate ug C per g soil from the CFE:
#   Using 0.174 g of salt per sample (0.174 g K2SO4 per X g dry soil)
df.CFE.merge <- df.CFE.merge %>%
  dplyr::mutate(g.dry.soil.fumigated = wet.soil.mass.g.fumigated/(1+incub.soil.moisture),
                g.dry.soil.not.fumigated = wet.soil.mass.g.not.fumigated/(1+incub.soil.moisture),
                g.C.DOC.per.g.dry.soil.fumigated = percent.C.DOC.fumigated/100 * 0.174/g.dry.soil.fumigated,
                ### Objective - convert the percent DOC-C to amount of DOC-C (g) per g dry soil
                  # % DOC = (g DOC / g salt) x 100
                  # (g DOC / g salt) x (g salt / g dry soil) = g DOC / g dry soil
                  # Used 0.174 g of salt for each sample (0.174 g K2SO4 per X g dry soil)
                g.C.DOC.per.g.dry.soil.not.fumigated = percent.C.DOC.notfumigated/100 * 0.174/g.dry.soil.not.fumigated,
                ug.C.DOC.per.g.dry.soil.fumigated = g.C.DOC.per.g.dry.soil.fumigated*10^6,
                ug.C.DOC.per.g.dry.soil.not.fumigated = g.C.DOC.per.g.dry.soil.not.fumigated*10^6)

  # dplyr::mutate(ug_C_per_g_soil.fumigated = percent.C.DOC.fumigated/100 * 0.174/wet.soil.mass.g.fumigated*10^6) %>%
    
  


### EQ 1 & 2. Adjust DOC-C concentrations for fum. and nonfum samples ----------

#     C_f_corrected = C_f_sample - C_f_blank      # **all blanks were below the detection limit.
#     C_nf_corrected = C_nf_sample - C_nf_blank   # **all blanks were below the detection limit.

### EQ 3. Calculate MBC --------------------------------------------------------

#     Total MBC = (C_f_corrected - C_nf_corrected)/(k_EC) 
df.calculations  <- df.CFE.merge %>%
  dplyr::mutate(ug.MBC.per.g.dry.soil = (ug.C.DOC.per.g.dry.soil.fumigated - ug.C.DOC.per.g.dry.soil.not.fumigated)/k_EC)
  

df.calculations <- df.calculations %>%
  group_by(core.id, horizon) %>%
  dplyr::mutate(mean.ug.MBC.per.g.dry.soil = mean(ug.MBC.per.g.dry.soil, na.rm=TRUE))
  # Calculates mean MBC for a horizon by averaging across 4 amendments
  # Result is 54 unique MBC values



### Calculate ug glucose, pine, and chloroform added per g soil ----------------

# import C data
df.C.added <- read.csv('../data/CUE/incubation-set-up.csv')
df.CN <- read.csv('../data/soil-properties/total-C-and-N.csv')

# reformat amendment and fumigation
for (i in 1:nrow(df.C.added)) {
  if (df.C.added$amendment[i] == '12C glucose') {
    df.C.added$amendment[i] = '12CG'
  } else if (df.C.added$amendment[i] == '13C glucose') {
    df.C.added$amendment[i] = '13CG'
  } else if (df.C.added$amendment[i] == '12C pine') {
    df.C.added$amendment[i] = '12CP'
  } else if (df.C.added$amendment[i] == '13C pine') {
    df.C.added$amendment[i] = '13CP'
  }
}


df.CN.filter <- df.CN %>%
  separate(sample.id, c('year','project','site','core', 'horizon'), sep = '-', remove = TRUE) %>%
  unite('core.id', c(year, project,site,core), sep='-', remove = TRUE) %>%
# I can't merge the CN data with the CUE data because the CN data = 2 days post-fire
#   and the CUE data = 3 weeks post-fire
  merge(df.meta, by = c('core.id','horizon')) %>%
  subset(select = c(site, horizon, burn.trtmt.duration.seconds, C.percent))


# adjust site IDs...
for (i in 1:nrow(df.CN.filter)){
  if (df.CN.filter$site[i] == 1) {
    df.CN.filter$site[i] = '01'
  } else if (df.CN.filter$site[i] == '2')  {
    df.CN.filter$site[i] = '02'
  } else if (df.CN.filter$site[i] == '3') {
    df.CN.filter$site[i] = '03'
  } else if (df.CN.filter$site[i] == '4') {
    df.CN.filter$site[i] = '04'
  } else if (df.CN.filter$site[i] == '5') {
    df.CN.filter$site[i] = '05'
  } else if (df.CN.filter$site[i] == '6') {
    df.CN.filter$site[i] = '06'
  } else if (df.CN.filter$site[i] == '7') {
    df.CN.filter$site[i] = '07'
  } else if (df.CN.filter$site[i] == '8') {
    df.CN.filter$site[i] = '08'
  } else if (df.CN.filter$site[i] == '9') {
    df.CN.filter$site[i] = '09'
  } else {
    df.CN.filter$site[i] = df.CN.filter$site[i]
  }
}

# Merge C data with meta data and MBC data
df.C.added <- subset(df.C.added,select = c(core.id, horizon, amendment, dry.soil.mass.initial..g.,
                                pine.amendment.conc...mg.mL.,glucose.amendment.conc...mg.mL.,
                                mL.pine.solution.to.add, mL.glucose.solution.to.add)) %>%
  merge(df.calculations, by = c('core.id', 'horizon', 'amendment')) %>%
  merge(subset(df.meta, select = c(core.id, burn.trtmt.duration.seconds))) %>%
  merge(df.CN.filter, by = c('site','horizon','burn.trtmt.duration.seconds')) %>%
  # Calculate ug C added as glucose, pine, or chloroform
  # Glucose = approx. 40% C
  # Pine roots = 42% C
  dplyr::mutate(mg.glucose.added = glucose.amendment.conc...mg.mL.*mL.glucose.solution.to.add,
         mg.pine.added = pine.amendment.conc...mg.mL.*mL.pine.solution.to.add,
         mg.glucose.per.g.dry.soil = mg.glucose.added/dry.soil.mass.initial..g.,
         mg.pine.per.g.dry.soil = mg.pine.added/dry.soil.mass.initial..g.,
         mg.glucose.C.per.g.dry.soil = mg.glucose.per.g.dry.soil*0.4,
         mg.pine.C.per.g.dry.soil = mg.pine.per.g.dry.soil*0.42) %>%
  dplyr::mutate(g.C.per.g.dry.soil = C.percent/100,
                glucose.C.added.as.percent.of.total.C = (mg.glucose.C.per.g.dry.soil/10^3)/g.C.per.g.dry.soil*100,
                pine.C.added.as.percent.of.total.C = (mg.pine.C.per.g.dry.soil/10^3)/g.C.per.g.dry.soil*100) %>%
  # Chloroform-C added = 0.5 mL added * density (1.4832 g/mL)
  dplyr::mutate(mg.chloroform.added.per.g.dry.soil = 0.5*1.4832/g.dry.soil.fumigated*1000,
         # chloroform = CHCl3, molar mass = 119.38 g/mol, C:CHCl3 = 1:1
         mg.chloroform.C.added.per.g.dry.soil = mg.chloroform.added.per.g.dry.soil/119.39*12.01,
         # ug.chloroform.C.added.per.g.dry.soil = g.chloroform.added.per.g.dry.soil/119.39*12.01*10^6,
         chloroform.C.added.as.percent.of.total.C = (mg.chloroform.C.added.per.g.dry.soil/10^3)/(g.C.per.g.dry.soil)*100)



### EQ 4 & 5: Calculate atom% of 13C in total MBC from soils -------------------

# 4. Calculate at% MBC_soil+labelled amendmend
# 5. Calculate at% MBC_soil+unlabelled amendment 
#     at_perc_MBC_sample = ((at_perc_C_f_sample * C_f_corrected) - (at_perc_C_nf_sample * C_nf_corrected))/(C_f_corrected - C_nf_corrected)


df.calculations <- df.calculations %>%
  # Calculate the at% of MBC from soil + amendment (either labelled or unlabelled)
  dplyr::mutate(at_perc_MBC_sample = ((atom.percent.13C.12C.fumigated * ug.C.DOC.per.g.dry.soil.fumigated) - 
                                 (atom.percent.13C.12C.notfumigated * ug.C.DOC.per.g.dry.soil.not.fumigated))/
                              (ug.C.DOC.per.g.dry.soil.fumigated-ug.C.DOC.per.g.dry.soil.not.fumigated))


## Reorganize dataframe to separate at% MBC of trtmt and control -------------

# separate dataframe by amended vs. unamended samples and pine vs glucose
df.calc.pine.sample <- subset(df.calculations, amendment == '13CP') %>%
  separate(sample.id.amend, c('sample.id', 'amendment'), sep = "_", remove = TRUE) %>%
  subset(select = c(sample.id,core.id, site.type, amendment, 
                    mean.ug.MBC.per.g.dry.soil,ug.MBC.per.g.dry.soil,
                    at_perc_MBC_sample, 
                    horizon,dry.soil.mass.g))


df.calc.pine.control <- subset(df.calculations, amendment == '12CP') %>%
  separate(sample.id.amend, c('sample.id', 'amendment'), sep = "_", remove = TRUE) %>%
  subset(select = c(sample.id, at_perc_MBC_sample))

colnames(df.calc.pine.control)[2] <- c('at_perc_MBC_control')


df.calc.glucose.sample <- subset(df.calculations, amendment == '13CG') %>%
  separate(sample.id.amend, c('sample.id', 'amendment'), sep = "_", remove = TRUE) %>%
  subset(select = c(sample.id, core.id, site.type, amendment, 
                    mean.ug.MBC.per.g.dry.soil, ug.MBC.per.g.dry.soil,
                    at_perc_MBC_sample, 
                    horizon,dry.soil.mass.g))

df.calc.glucose.control <- subset(df.calculations, amendment == '12CG') %>%
  separate(sample.id.amend, c('sample.id', 'amendment'), sep = "_", remove = TRUE) %>%
  subset(select = c(sample.id, at_perc_MBC_sample))

colnames(df.calc.glucose.control)[2] <- c('at_perc_MBC_control')


# recombine into pine and glucose dataframes
df.pine <- merge(df.calc.pine.sample, df.calc.pine.control, by = 'sample.id')
df.glucose <- merge(df.calc.glucose.sample, df.calc.glucose.control, by = 'sample.id')


df.pine$horizon <- factor(df.pine$horizon, levels = c('O','M'))
df.glucose$horizon <- factor(df.glucose$horizon, levels = c('O','M'))


### EQ 6. Fraction of total MBC derived from the labelled amendment ------------

#     Option 1: f_MBC,LA = (at_perc_MBC_soil_LA - at_perc_soil) / (at_perc_LA - at_perc_soil)
#                 Where LA = "Labelled amendment"

#     Option 2: f_MBC,LA = (at_perc_MBC_sample - at_perc_MBC_control)/(at_perc_sol - at_perc_unlabelled_amendment)

df.pine <- df.pine %>%
  # Option 1
  # dplyr::mutate(fraction_of_LA_derived_MBC = (at_perc_MBC_sample - at_perc_sandy_soil) / (at_perc_labelled_pine - at_perc_sandy_soil)) %>%
  # Option 2
  dplyr::mutate(fraction_of_LA_derived_MBC = (at_perc_MBC_sample - at_perc_MBC_control)/(at_perc_labelled_pine - at_perc_unlabelled_pine))
  

df.glucose <- df.glucose %>%
  # Option 1
  # dplyr::mutate(fraction_of_LA_derived_MBC = (at_perc_MBC_sample - at_perc_org_soil) / (at_perc_labelled_glucose - at_perc_org_soil)) %>%
  # Option 2
  dplyr::mutate(fraction_of_LA_derived_MBC = (at_perc_MBC_sample - at_perc_MBC_control)/(at_perc_labelled_glucose - at_perc_unlabelled_glucose))


### EQ 7. Calculate substrate-derived MBC, (ug C g^-1 soil) --------------------

#     Substrate-derived MBC = total MBC * fraction_of_LA_derived_MBC
df.pine <- df.pine %>%
  dplyr::mutate(ug_substrate_derived_MBC_per_g_dry_soil = ug.MBC.per.g.dry.soil * fraction_of_LA_derived_MBC) %>%
  # dplyr::mutate(alternative_ug_substrate_derived_MBC_per_g_dry_soil = alternative.mean.ug.MBC.per.g.dry.soil * fraction_of_LA_derived_MBC) %>%
  dplyr::mutate(ug_soil_derived_MBC_per_g_dry_soil  = ug.MBC.per.g.dry.soil * (1-fraction_of_LA_derived_MBC)) %>%
  separate(core.id, c('year','project','site','core'), sep = "-", remove = FALSE) 


df.glucose <- df.glucose %>%
  # use mean ug MBC per g dry soil for glucose samples to deal with negative MBC values
  dplyr::mutate(ug_substrate_derived_MBC_per_g_dry_soil = ug.MBC.per.g.dry.soil * fraction_of_LA_derived_MBC) %>%
  # dplyr::mutate(alternative_ug_substrate_derived_MBC_per_g_dry_soil = alternative.mean.ug.MBC.per.g.dry.soil * fraction_of_LA_derived_MBC) %>%
  dplyr::mutate(ug_soil_derived_MBC_per_g_dry_soil  = ug.MBC.per.g.dry.soil * (1-fraction_of_LA_derived_MBC)) %>%
  separate(core.id, c('year','project','site','core'), sep = "-", remove = FALSE)



### RESPIRATION: Import and clean up resp. data --------------------------------

# Import data
df.CO2 <- read.csv('../data/CUE/Picarro/calculated-mg-C-respired.csv')  
head(df.CO2,3)
colnames(df.CO2)


# reformat column names:
df.CO2.glucose.labelled <- df.CO2 %>%
  subset(carbon.addition == '13CG') %>%
  dplyr::mutate(delta_13C_of_CO2_from_soils_with_labelled_amendment = Delta_30s_iCO2_mean) %>%
  subset(select = c(site,core,horizon, carbon.addition,timepoint,
                    delta_13C_of_CO2_from_soils_with_labelled_amendment))

df.CO2.glucose.unlabelled <- df.CO2 %>%
  subset(carbon.addition == '12CG') %>%
  dplyr::mutate(delta_13C_of_CO2_from_soils_with_unlabelled_amendment = Delta_30s_iCO2_mean)%>%
  subset(select = c(site,core,horizon,timepoint,
                    delta_13C_of_CO2_from_soils_with_unlabelled_amendment))

df.CO2.glucose <- merge(df.CO2.glucose.labelled, df.CO2.glucose.unlabelled,
                        by=c('site','core','horizon','timepoint')) 

## Repeat for Pine:
df.CO2.pine.labelled <- df.CO2 %>%
  subset(carbon.addition == '13CP') %>%
  dplyr::mutate(delta_13C_of_CO2_from_soils_with_labelled_amendment = Delta_30s_iCO2_mean) %>%
  subset(select = c(site,core,horizon, carbon.addition,timepoint,
                    delta_13C_of_CO2_from_soils_with_labelled_amendment))

df.CO2.pine.unlabelled <- df.CO2 %>%
  subset(carbon.addition == '12CP') %>%
  dplyr::mutate(delta_13C_of_CO2_from_soils_with_unlabelled_amendment = Delta_30s_iCO2_mean)%>%
  subset(select = c(site,core,horizon,timepoint,
                    delta_13C_of_CO2_from_soils_with_unlabelled_amendment))

df.CO2.pine <- merge(df.CO2.pine.labelled, df.CO2.pine.unlabelled,
                     by=c('site','core','horizon','timepoint'))



### EQ 8. Fraction of total CO2 derived from labelled substrate ----------------

df.CO2.glucose.sandy <- subset(df.CO2.glucose, site %in% c(1,2,3,5,8,9)) %>%
  dplyr::mutate(fraction_of_CO2_from_labelled_amendment = (delta_13C_of_CO2_from_soils_with_labelled_amendment- d13C_sandy_soil)/(d13C_lab_glucose-d13C_sandy_soil)) %>%
  dplyr::mutate(fraction_of_soil_derived_CO2 = 1-fraction_of_CO2_from_labelled_amendment)

df.CO2.glucose.org <- subset(df.CO2.glucose, site %in% c(4,6,7,10,11,12)) %>%
  dplyr::mutate(fraction_of_CO2_from_labelled_amendment = (delta_13C_of_CO2_from_soils_with_labelled_amendment- d13C_org_soil)/(d13C_lab_glucose-d13C_org_soil)) %>%
  dplyr::mutate(fraction_of_soil_derived_CO2 = 1-fraction_of_CO2_from_labelled_amendment)

df.CO2.glucose.merge <- 
  rbind(subset(df.CO2.glucose.sandy, select = c(site,core,horizon,timepoint,
                                                carbon.addition,fraction_of_CO2_from_labelled_amendment,
                                                fraction_of_soil_derived_CO2)),
        subset(df.CO2.glucose.org, select = c(site,core,horizon,timepoint,
                                              carbon.addition,fraction_of_CO2_from_labelled_amendment,
                                              fraction_of_soil_derived_CO2))) 


# Repeat for pine: 
df.CO2.pine.sandy <- subset(df.CO2.pine, site %in% c(1,2,3,5,8,9)) %>%
  dplyr::mutate(fraction_of_CO2_from_labelled_amendment = (delta_13C_of_CO2_from_soils_with_labelled_amendment- d13C_sandy_soil)/(d13C_lab_pine-d13C_sandy_soil)) %>%
  dplyr::mutate(fraction_of_soil_derived_CO2 = 1-fraction_of_CO2_from_labelled_amendment)

df.CO2.pine.org <- subset(df.CO2.pine, site %in% c(4,6,7,10,11,12)) %>%
  dplyr::mutate(fraction_of_CO2_from_labelled_amendment = (delta_13C_of_CO2_from_soils_with_labelled_amendment- d13C_org_soil)/(d13C_lab_pine-d13C_org_soil)) %>%
  dplyr::mutate(fraction_of_soil_derived_CO2 = 1-fraction_of_CO2_from_labelled_amendment)

df.CO2.pine.merge <- 
  rbind(subset(df.CO2.pine.sandy, select = c(site,core,horizon,timepoint,
                                                carbon.addition,fraction_of_CO2_from_labelled_amendment,
                                             fraction_of_soil_derived_CO2)),
        subset(df.CO2.pine.org, select = c(site,core,horizon,timepoint,
                                              carbon.addition,fraction_of_CO2_from_labelled_amendment,
                                           fraction_of_soil_derived_CO2))) 


### EQ 9. Calculate cumulative substrate-derived CO2 ---------------------------

df.CO2.substrates <- rbind(df.CO2.glucose.merge, df.CO2.pine.merge)

## Multiply substrate-derived CO2 fraction x total CO2 concentration --> substrate-derived CO2
df.CO2.sub.derived <- merge(df.CO2, df.CO2.substrates, by = c('site', 'core', 'horizon', 'carbon.addition', 'timepoint')) %>% # removes 12CG and 12CP amendments
  dplyr::mutate(mg.substrate.derived.CO2 = fraction_of_CO2_from_labelled_amendment*mg_CO2_C_total,
                mg.soil.derived.CO2 = mg_CO2_C_total-mg.substrate.derived.CO2) %>%
  ## Calculate cumulative substrate-derived CO2 across all four timepoints.
  dplyr::group_by(site,core, horizon, carbon.addition) %>%
  dplyr::mutate(cumulative.mg.substrate.derived.CO2 = sum(mg.substrate.derived.CO2),
                cumulative.mg.substrate.derived.CO2.per.g.dry.soil = cumulative.mg.substrate.derived.CO2/dry.soil.mass.g,
                cumulative.ug.substrate.derived.CO2.per.g.dry.soil=cumulative.mg.substrate.derived.CO2.per.g.dry.soil*1000,
                cumulative.mg.soil.derived.CO2 = sum(mg.soil.derived.CO2),
                cumulative.mg.soil.derived.CO2.per.g.dry.soil = cumulative.mg.soil.derived.CO2/dry.soil.mass.g,
                cumulative.ug.soil.derived.CO2.per.g.dry.soil=cumulative.mg.soil.derived.CO2.per.g.dry.soil*1000)


### Clean up dataframe and add treatment data ----------------------------------

df.CO2.clean <- df.CO2.sub.derived %>%
  dplyr::mutate(amendment=carbon.addition) %>%
  separate(sample.id, c('sam', 'amend', 'time'), sep = '_', remove = TRUE) %>%
  dplyr::mutate(sample.id = sam) %>%
  subset(select = c(site,core,horizon,amendment,sample.id,
                    cumulative.ug.substrate.derived.CO2.per.g.dry.soil,
                    cumulative.ug.soil.derived.CO2.per.g.dry.soil,
                    cumulative.mg.soil.derived.CO2.per.g.dry.soil)) %>%
  unique()
  

### PLOT: Figure 2B: Cumul. soil-derived CO2-C, combined pine and glucose: ----------------

# Stats for CO2:
df.CO2.plot <- merge(df.CO2.clean, 
                     subset(df.C.added, select = c(site.type, sample.id,substrate, amendment,
                                                   burn.trtmt.duration.seconds)), 
                     by=c('sample.id', 'amendment')) %>%
  unique()

df.CO2.plot$burn.trtmt.duration.seconds <- as.factor(df.CO2.plot$burn.trtmt.duration.seconds)

### Creating dataframe of significant letters to add to plots. Method 2:
df.CO2.combo.stats = data.frame(burn.trtmt.duration.seconds = 'duration', 
                          site.type = 'soil', 
                          horizon = 'k',
                          max.CO2 = 0,
                          cld = 'cld')

# Now create list of significant letters for ggplot:
for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CO2.plot, site.type == j & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(cumulative.ug.soil.derived.CO2.per.g.dry.soil~as.character(site) + burn.trtmt.duration.seconds, data = x)
      # Interaction between site.type and burn treatment 
      #    (site.type*burn.trtmt.duration.seconds) is not significant so it is 
      #    not including in the model, and we report results without the 
      #    interaction. 
      tukey <- TukeyHSD(anova)
      cld <- multcompLetters4(anova, tukey)
      x.cld <- group_by(subset(df.CO2.plot, 
                               site.type == j & 
                                 horizon == k), burn.trtmt.duration.seconds) %>%
        subset(select = c(burn.trtmt.duration.seconds, site.type, horizon)) %>%
        unique()%>%
        dplyr::arrange((burn.trtmt.duration.seconds)) %>%
        dplyr::mutate(site.type = j)
      cld <- as.data.frame.list(cld$burn.trtmt.duration.seconds)
      cld <- cld[order(as.numeric(row.names(cld))),]
      x.cld$cld <- cld$Letters
      x.cld$max.CO2 = max(x$cumulative.ug.soil.derived.CO2.per.g.dry.soil, na.rm = TRUE)
      df.CO2.combo.stats <- rbind(df.CO2.combo.stats, x.cld)
    }
  }
}

df.CO2.combo.stats <- subset(df.CO2.combo.stats, burn.trtmt.duration.seconds != 'duration')




# Now create list of p-values :
df.CO2.combo.pvalues = data.frame(ID = 'ID',comparison = 'comparison',p.value=0)

for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CO2.plot, site.type == j & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(cumulative.ug.soil.derived.CO2.per.g.dry.soil~as.character(site) + burn.trtmt.duration.seconds, data = x)
      tukey <- TukeyHSD(anova)
      comparison = rownames(tukey$burn.trtmt.duration.seconds)
      p.value = tukey$burn.trtmt.duration.seconds[10:12]
      temp <- data.frame(ID = paste(j,k), comparison = comparison, p.value)
      df.CO2.combo.pvalues <- rbind(df.CO2.combo.pvalues, temp)
    }
  }
}

df.CO2.combo.pvalues <- subset(df.CO2.combo.pvalues, ID != 'ID')
df.CO2.combo.stats$horizon <- factor(df.CO2.combo.stats$horizon, level=c('O','M'))
df.CO2.plot$horizon <- factor(df.CO2.plot$horizon, level=c('O','M'))

df.CO2.combo.stats.sub <- subset(df.CO2.combo.stats, horizon != 'M')

# Plot:
# Create each panel separately to allow for viewing mineral soil panel.
pCO2_combined_ORG_O = ggplot(subset(df.CO2.plot, site.type=='ORG' & horizon == 'O'), 
                             aes(x=(burn.trtmt.duration.seconds), 
                                 y=cumulative.ug.soil.derived.CO2.per.g.dry.soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.CO2.combo.stats.sub, site.type=='ORG' & horizon == 'O'), 
            aes(y=max.CO2+0.1, label = cld), vjust = -1, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       y = expression(atop(`Soil-derived`~`CO`[2]-C,(mu*g~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') + 
  ylim(0,max(subset(df.CO2.plot, site.type=='ORG' & horizon == 'O')$cumulative.ug.soil.derived.CO2.per.g.dry.soil, na.rm = TRUE)*1.25)


pCO2_combined_SANDY_O = ggplot(subset(df.CO2.plot, site.type=='SANDY' & horizon == 'O'), 
                               aes(x=(burn.trtmt.duration.seconds), 
                                   y=cumulative.ug.soil.derived.CO2.per.g.dry.soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.CO2.combo.stats.sub, site.type=='SANDY' & horizon == 'O'),
            aes(y=max.CO2+0.1, label = cld), vjust = -1, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       y = expression(atop(`Soil-derived`~`CO`[2]-C,(mu*g~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') + 
  ylim(0,max(subset(df.CO2.plot, site.type=='SANDY' & horizon == 'O')$cumulative.ug.soil.derived.CO2.per.g.dry.soil, na.rm = TRUE)*1.25)


pCO2_combined_SANDY_M = ggplot(subset(df.CO2.plot, site.type=='SANDY' & horizon == 'M'), 
                               aes(x=(burn.trtmt.duration.seconds), 
                                   y=cumulative.ug.soil.derived.CO2.per.g.dry.soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.CO2.combo.stats.sub, site.type=='SANDY' & horizon == 'M'),
            aes(y=max.CO2+0.1, label = cld), vjust = -1, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       y = expression(atop(`Soil-derived`~`CO`[2]-C,(mu*g~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') +
  ylim(0,max(subset(df.CO2.plot, site.type=='SANDY' & horizon == 'M')$cumulative.ug.soil.derived.CO2.per.g.dry.soil, na.rm = TRUE)*1.05)


cowplot::plot_grid(pCO2_combined_ORG_O, pCO2_combined_SANDY_O, pCO2_combined_SANDY_M, ncol = 3)

# ggsave('CUE_figure_X_soil_derived_CO2', path = "../manuscripts/CUE/figures/",
#        device='png', plot = last_plot(),
#        width = 6.5, height=2.5, units="in")

paste(sprintf("Figure X. Caption details - Soil-derived CO2-C for all samples "),
      sprintf("n=108 total; n=12 for each boxplot (6 cores for each burn trtmt x 2 amendments;"),
      sprintf("stats = ANOVA with TukeyHSD test, letters indicated significant differences at p<0.05 as cutoff"))


# stats for CO2
df.stats.CO2 = data.frame(ID = 'ID',MEAN = 0, SD = 0)

for (j in c('ORG', 'SANDY')) {
  for (k in c('O', 'M')) {
    for (i in c(0,30,120)) {
      x <- subset(df.CO2.plot, site.type == j & horizon == k & burn.trtmt.duration.seconds == i)
      if (dim(x)[1] == 0){
        next
      }
      mean = mean(x$cumulative.ug.soil.derived.CO2.per.g.dry.soil, na.rm=TRUE)
      sd = sd(x$cumulative.ug.soil.derived.CO2.per.g.dry.soil, na.rm=TRUE)
      ID = paste(j,k,i)
      temp <- data.frame(ID = ID, MEAN=mean, SD = sd)
      df.stats.CO2 = rbind(df.stats.CO2, temp)
    }
  }
}

df.stats.CO2

paste(sprintf("Organic, unburned: mean="),sprintf(subset(df.stats.CO2, ID == 'ORG O 0')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'ORG O 0')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.CO2, ID == 'ORG O 30')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'ORG O 30')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.CO2, ID == 'ORG O 120')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'ORG O 120')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.CO2, ID == 'SANDY O 0')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'SANDY O 0')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.CO2, ID == 'SANDY O 30')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'SANDY O 30')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.CO2, ID == 'SANDY O 120')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'SANDY O 120')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.CO2, ID == 'SANDY M 0')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'SANDY M 0')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.CO2, ID == 'SANDY M 30')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'SANDY M 30')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.CO2, ID == 'SANDY M 120')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.CO2, ID == 'SANDY M 120')$SD, fmt='%#.1f'),"ug C per g dry soil")


# Now create list of p-values :
df.substrate.CO2.pvalues = data.frame(ID = 'ID',comparison = 'comparison',
                                      horizon = 'horizon', 
                                      p.value=0)


for (j in c('ORG', 'SANDY')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CO2.plot, site.type == j & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(cumulative.ug.soil.derived.CO2.per.g.dry.soil~as.character(site) + burn.trtmt.duration.seconds, data = x)
      tukey <- TukeyHSD(anova)
      comparison = rownames(tukey$burn.trtmt.duration.seconds)
      p.value = tukey$burn.trtmt.duration.seconds[10:12]
      temp <- data.frame(ID = paste(j,k), horizon = k, 
                         comparison = comparison, p.value)
      df.substrate.CO2.pvalues <- rbind(df.substrate.CO2.pvalues, temp)
  }
}

df.substrate.CO2.pvalues <- subset(df.substrate.CO2.pvalues, ID != 'ID')
df.substrate.CO2.pvalues

# print stats:
paste(sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'ORG O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'ORG O' & comparison == '120-0')$p.value, 
              fmt='%#.5f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'ORG O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'SANDY O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'SANDY O' & comparison == '120-0')$p.value, 
              fmt='%#.3e'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'SANDY O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'SANDY M' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'SANDY M' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.CO2.pvalues, ID == 'SANDY M' & comparison == '120-30')$p.value, 
              fmt='%#.3f'))





### PLOT: Figure 4. Cumulative substrate-derived CO2-C -----

# Stats for CO2:
### Creating dataframe of significant letters to add to plots. 
df.CO2.substrate.stats = data.frame(burn.trtmt.duration.seconds = 'duration', 
                          site.type = 'soil', 
                          horizon = 'k',
                          substrate = 'incubation.period.days',
                          max.CO2 = 0,
                          cld = 'cld')

# Now create list of significant letters for ggplot:
for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CO2.plot, site.type == j & substrate == h & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(cumulative.ug.substrate.derived.CO2.per.g.dry.soil~as.character(site) + burn.trtmt.duration.seconds, data = x)
      # Interaction between site.type and burn treatment 
      #    (site.type*burn.trtmt.duration.seconds) is not significant so it is 
      #    not including in the model, and we report results without the 
      #    interaction. 
      tukey <- TukeyHSD(anova)
      cld <- multcompLetters4(anova, tukey)
      x.cld <- group_by(subset(df.CO2.plot, 
                               site.type == j & 
                                 horizon == k &
                                 substrate == h), burn.trtmt.duration.seconds) %>%
        subset(select = c(burn.trtmt.duration.seconds, site.type, substrate, horizon)) %>%
        unique()%>%
        dplyr::arrange((burn.trtmt.duration.seconds)) %>%
        dplyr::mutate(substrate = h, site.type = j)
      cld <- as.data.frame.list(cld$burn.trtmt.duration.seconds)
      cld <- cld[order(as.numeric(row.names(cld))),]
      x.cld$cld <- cld$Letters
      x.cld$max.CO2 = max(x$cumulative.ug.substrate.derived.CO2.per.g.dry.soil, na.rm = TRUE)
      df.CO2.substrate.stats <- rbind(df.CO2.substrate.stats, x.cld)
    }
  }
}

df.CO2.substrate.stats <- subset(df.CO2.substrate.stats, burn.trtmt.duration.seconds != 'site.type')


# Now create list of p-values :
df.CO2.substrate.pvalues = data.frame(ID = 'ID',comparison = 'comparison',p.value=0)

for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CO2.plot, site.type == j & substrate == h & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(cumulative.ug.substrate.derived.CO2.per.g.dry.soil~as.character(site) + burn.trtmt.duration.seconds, data = x)
      tukey <- TukeyHSD(anova)
      comparison = rownames(tukey$burn.trtmt.duration.seconds)
      p.value = tukey$burn.trtmt.duration.seconds[10:12]
      temp <- data.frame(ID = paste(j,h,k), comparison = comparison, p.value)
      df.CO2.substrate.pvalues <- rbind(df.CO2.substrate.pvalues, temp)
    }
  }
}

df.CO2.substrate.pvalues <- subset(df.CO2.substrate.pvalues, ID != 'ID')
df.CO2.substrate.stats$horizon <- factor(df.CO2.substrate.stats$horizon, level=c('O','M'))
df.CO2.plot$horizon <- factor(df.CO2.plot$horizon, level=c('O','M'))

df.CO2.substrate.stats.sub <- subset(df.CO2.substrate.stats, horizon != 'M')


# Plot total CO2-C respired:
pCO2_substrate_pine = ggplot(subset(df.CO2.plot, substrate == 'pine'),
                   # aes(x=site.type, y = CUE, color = as.character(burn.trtmt.duration.seconds), shape=horizon )) +
                   aes(x=burn.trtmt.duration.seconds,y=cumulative.ug.substrate.derived.CO2.per.g.dry.soil, 
                       fill=burn.trtmt.duration.seconds))+ 
  geom_boxplot(alpha=0.7)+
  # geom_jitter(aes(shape=horizon,color = burn.trtmt.duration.seconds), width = 0.2, alpha=0.7,size=2) + 
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type=SiteType_labels)) +
  # scale_colour_manual(values = c('black','darksalmon','darkred'),
  #                     breaks= c('0', '30','120'))+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  # geom_hline(yintercept=1) + 
  # geom_hline(yintercept=0) + 
  # geom_text(data=subset(df.CO2.substrate.stats.sub, substrate == 'pine'),
            # aes(y=max.CO2*1.1, label = cld), vjust = 0, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       y = expression(atop(`Pine-derived`~`CO`[2]-C,(mu*g~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12),
        axis.text=element_text(size=12), 
        strip.text = element_text(size=12))+
  theme(legend.position ='')+
  ylim(0,max(subset(df.CO2.plot, substrate == 'pine')$cumulative.ug.substrate.derived.CO2.per.g.dry.soil, na.rm = TRUE)*1.2)


pCO2_substrate_glucose = ggplot(subset(df.CO2.plot, substrate == 'glucose'),
                      # aes(x=site.type, y = CUE, color = as.character(burn.trtmt.duration.seconds), shape=horizon )) +
                      aes(x=burn.trtmt.duration.seconds,y=cumulative.ug.substrate.derived.CO2.per.g.dry.soil, 
                          fill=burn.trtmt.duration.seconds))+ 
  geom_boxplot(alpha=0.7)+
  # geom_jitter(aes(shape=horizon,color = burn.trtmt.duration.seconds), width = 0.2, alpha=0.7,size=2) + 
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  # scale_colour_manual(values = c('black','darksalmon','darkred'),
  #                     breaks= c('0', '30','120'))+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  # geom_hline(yintercept=1) + 
  # geom_hline(yintercept=0) + 
  geom_text(data=subset(df.CO2.substrate.stats.sub, substrate == 'glucose'  & horizon == 'O' & site.type == 'SANDY'),
            aes(y=max.CO2*1.1, label = cld), vjust = 0, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       # y = expression(atop(Respired~soil~C~derived~`CO`[2],(mu*g~per~g~dry~soil))),
       y = expression(atop(`Glucose-derived`~`CO`[2]-C,(mu*g~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12),
        axis.text=element_text(size=12), 
        strip.text = element_text(size=12))+
  theme(legend.position ='')+
  ylim(0, max(subset(df.CO2.plot, substrate == 'glucose')$cumulative.ug.substrate.derived.CO2.per.g.dry.soil, na.rm = TRUE)*1.2)


cowplot::plot_grid(pCO2_substrate_pine, pCO2_substrate_glucose, ncol=1)



# stats for CO2
df.stats.substrate.CO2 = data.frame(ID = 'ID',MEAN = 0, SD = 0)

for (j in c('ORG', 'SANDY')) {
  for (k in c('O', 'M')) {
    for (i in c(0,30,120)) {
      for (l in c('glucose','pine')) {
      x <- subset(df.CO2.plot, site.type == j & horizon == k & burn.trtmt.duration.seconds == i & substrate == l)
      if (dim(x)[1] == 0){
        next
      }
      mean = mean(x$cumulative.ug.substrate.derived.CO2.per.g.dry.soil, na.rm=TRUE)
      sd = sd(x$cumulative.ug.substrate.derived.CO2.per.g.dry.soil, na.rm=TRUE)
      ID = paste(j,k,i, l)
      temp <- data.frame(ID = ID, MEAN=mean, SD = sd)
      df.stats.substrate.CO2 = rbind(df.stats.substrate.CO2, temp)
      }
    }
  }
}

df.stats.substrate.CO2

paste('GLUCOSE: ', 
  sprintf("Organic, unburned: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 0 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 0 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 30 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 30 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 120 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 120 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 0 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 0 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 30 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 30 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 120 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 120 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 0 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 0 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 30 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 30 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 120 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 120 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  'PINE: ', 
  sprintf("Organic, unburned: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 0 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 0 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
  sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 30 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 30 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
  sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 120 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'ORG O 120 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 0 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 0 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 30 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 30 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 120 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY O 120 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 0 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 0 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 30 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 30 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
  sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 120 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
  sprintf(subset(df.stats.substrate.CO2, ID == 'SANDY M 120 pine')$SD, fmt='%#.1f'),"ug C per g dry soil")



# print stats: GLUCOSE
paste(sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'ORG glucose O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'ORG glucose O' & comparison == '120-0')$p.value, 
              fmt='%#.5f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'ORG glucose O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY glucose O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY glucose O' & comparison == '120-0')$p.value, 
              fmt='%#.3e'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY glucose O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY glucose M' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY glucose M' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY glucose M' & comparison == '120-30')$p.value, 
              fmt='%#.3f'))

# print stats: Pine
paste(sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'ORG pine O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'ORG pine O' & comparison == '120-0')$p.value, 
              fmt='%#.5f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'ORG pine O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY pine O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY pine O' & comparison == '120-0')$p.value, 
              fmt='%#.3e'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY pine O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY pine M' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY pine M' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CO2.substrate.pvalues, ID == 'SANDY pine M' & comparison == '120-30')$p.value, 
              fmt='%#.3f'))










### PLOT: Figure 2A: Soil-derived MBC, combined pine and glucose --------------------------

df.sub.MBC <- rbind(df.pine, df.glucose) %>%
  merge(subset(df.CO2.plot, select = c(sample.id, horizon, amendment, 
                                       cumulative.ug.substrate.derived.CO2.per.g.dry.soil,
                                       cumulative.ug.soil.derived.CO2.per.g.dry.soil,
                                       site.type, burn.trtmt.duration.seconds, substrate)), 
        by = c('sample.id', 'horizon', 'amendment', 'site.type'))


# Did substrate have an effect on SOC-derived MBC?
anova <- aov(ug_soil_derived_MBC_per_g_dry_soil~substrate, data = df.sub.MBC)
summary(anova)
# Answer: nope. 

# Interaction between substrate and burn treatment on SOC-derived MBC?
anova <- aov(ug_soil_derived_MBC_per_g_dry_soil~substrate*burn.trtmt.duration.seconds, data = df.sub.MBC)
summary(anova)
# Nope

### Creating dataframe of significant letters to add to plots. Method 2:
df.soil.MBC.combined.stats = data.frame(burn.trtmt.duration.seconds = 'duration', 
                                   site.type = 'soil', 
                                   horizon = 'k',
                                   max.soil.MBC = 0,
                                   cld = 'cld')

# Now create list of significant letters for ggplot:
for (j in c('ORG', 'SANDY')) {
  for (k in c('O', 'M')) {
    x <- subset(df.sub.MBC, site.type == j & horizon == k)
    if (dim(x)[1] == 0){
      next
    }
    anova <- aov(ug_soil_derived_MBC_per_g_dry_soil~as.character(site)+ burn.trtmt.duration.seconds, data = x)
    # Interaction between site.type and burn treatment 
    #    (site.type*burn.trtmt.duration.seconds) is not significant so it is 
    #    not including in the model, and we report results without the 
    #    interaction. 
    tukey <- TukeyHSD(anova)
    cld <- multcompLetters4(anova, tukey)
    x.cld <- group_by(subset(df.sub.MBC, 
                             site.type == j & 
                               horizon == k), burn.trtmt.duration.seconds) %>%
      subset(select = c(burn.trtmt.duration.seconds, site.type, horizon)) %>%
      unique()%>%
      dplyr::arrange((burn.trtmt.duration.seconds)) %>%
      dplyr::mutate(site.type = j)
    cld <- as.data.frame.list(cld$burn.trtmt.duration.seconds)
    cld <- cld[order(as.numeric(row.names(cld))),]
    x.cld$cld <- cld$Letters
    x.cld$max.soil.MBC = max(x$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE)
    df.soil.MBC.combined.stats <- rbind(df.soil.MBC.combined.stats, x.cld)
  }
}

df.soil.MBC.combined.stats <- subset(df.soil.MBC.combined.stats, site.type != 'soil')



# Now create list of p-values :
df.soil.MBC.combined.pvalues = data.frame(ID = 'ID',comparison = 'comparison',p.value=0)

for (j in c('ORG', 'SANDY')) {
  for (k in c('O', 'M')) {
    x <- subset(df.sub.MBC, site.type == j & horizon == k)
    if (dim(x)[1] == 0){
      next
    }
    anova <- aov(ug_soil_derived_MBC_per_g_dry_soil~as.character(site)+burn.trtmt.duration.seconds, data = x)
    tukey <- TukeyHSD(anova)
    comparison = rownames(tukey$burn.trtmt.duration.seconds)
    p.value = tukey$burn.trtmt.duration.seconds[10:12]
    temp <- data.frame(ID = paste(j,k), comparison = comparison, p.value)
    df.soil.MBC.combined.pvalues <- rbind(df.soil.MBC.combined.pvalues, temp)
  }
}

df.soil.MBC.combined.pvalues <- subset(df.soil.MBC.combined.pvalues, ID != 'ID')


df.soil.MBC.combined.stats.sub <- subset(df.soil.MBC.combined.stats, horizon != 'M')
df.soil.MBC.combined.stats.sub$horizon <- factor(df.soil.MBC.combined.stats.sub$horizon, level=c('O','M'))
df.sub.MBC$horizon <- factor(df.sub.MBC$horizon, level = c('O','M'))


pMBC_soil_combined = ggplot(subset(df.sub.MBC), 
                       aes(x=(burn.trtmt.duration.seconds), 
                           y=ug_soil_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.soil.MBC.combined.stats.sub), 
            aes(y=max.soil.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') + 
  ylim(min(subset(df.sub.MBC)$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE), 
       max(subset(df.sub.MBC)$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE)*1.45)

pMBC_soil_combined 


# Create each panel separately to allow for viewing mineral soil panel.
pMBC_soil_combined_ORG_O = ggplot(subset(df.sub.MBC, site.type=='ORG' & horizon == 'O'), 
                             aes(x=(burn.trtmt.duration.seconds), 
                                 y=ug_soil_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.soil.MBC.combined.stats.sub, site.type=='ORG' & horizon == 'O'), 
            aes(y=max.soil.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(`Soil-derived`~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') + 
  ylim(min(subset(df.sub.MBC, site.type=='ORG' & horizon == 'O')$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE), 
       max(subset(df.sub.MBC, site.type=='ORG' & horizon == 'O')$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE)*1.3)


pMBC_soil_combined_SANDY_O = ggplot(subset(df.sub.MBC, site.type=='SANDY' & horizon == 'O'), 
                               aes(x=(burn.trtmt.duration.seconds), 
                                   y=ug_soil_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.soil.MBC.combined.stats.sub, site.type=='SANDY' & horizon == 'O'), 
            aes(y=max.soil.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(`Soil-derived`~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') + 
  ylim(min(subset(df.sub.MBC, site.type=='ORG' & horizon == 'O')$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE), 
       max(subset(df.sub.MBC, site.type=='ORG' & horizon == 'O')$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE)*1.15)


pMBC_soil_combined_SANDY_M = ggplot(subset(df.sub.MBC, site.type=='SANDY' & horizon == 'M'), 
                               aes(x=(burn.trtmt.duration.seconds), 
                                   y=ug_soil_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.soil.MBC.combined.stats.sub, site.type=='SANDY' & horizon == 'M'), 
            aes(y=max.soil.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(`Soil-derived`~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '')  +
  # ylim(min(subset(df.sub.MBC, site.type=='SANDY' & horizon == 'M')$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE),
  # max(subset(df.sub.MBC, site.type=='SANDY' & horizon == 'M')$ug_soil_derived_MBC_per_g_dry_soil, na.rm = TRUE)*1.45)
  ylim(-127,120)

cowplot::plot_grid(pMBC_soil_combined_ORG_O, pMBC_soil_combined_SANDY_O, pMBC_soil_combined_SANDY_M, ncol = 3)

# ggsave('CUE_figure_1_Total_MBC', path = "../manuscripts/CUE/figures/",
#        device='png', plot = last_plot(),
#        width = 6.5, height=2.5, units="in")

paste(sprintf("Figure X. Caption details - Soil-derived microbial biomass C for all samples "),
      sprintf("n=108 total; n=12 for each boxplot (6 cores for each burn trtmt x 2 amendments;"),
      sprintf("stats = ANOVA with TukeyHSD test, letters indicated significant differences at p<0.05 as cutoff"))



# stats for MBC
df.stats.soil.MBC = data.frame(ID = 'ID',MEAN = 0, SD = 0)

for (j in c('ORG', 'SANDY')) {
  for (k in c('O', 'M')) {
    for (i in c(0,30,120)) {
      x <- subset(df.sub.MBC, site.type == j & horizon == k & burn.trtmt.duration.seconds == i)
      if (dim(x)[1] == 0){
        next
      }
      mean = mean(x$ug_soil_derived_MBC_per_g_dry_soil, na.rm=TRUE)
      sd = sd(x$ug_soil_derived_MBC_per_g_dry_soil, na.rm=TRUE)
      ID = paste(j,k,i)
      temp <- data.frame(ID = ID, MEAN=mean, SD = sd)
      df.stats.soil.MBC = rbind(df.stats.soil.MBC, temp)
    }
  }
}

df.stats.soil.MBC

paste(sprintf("Organic, unburned: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'ORG O 0')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'ORG O 0')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'ORG O 30')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'ORG O 30')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'ORG O 120')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'ORG O 120')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'SANDY O 0')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'SANDY O 0')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'SANDY O 30')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'SANDY O 30')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'SANDY O 120')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'SANDY O 120')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'SANDY M 0')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'SANDY M 0')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'SANDY M 30')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'SANDY M 30')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.soil.MBC, ID == 'SANDY M 120')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.soil.MBC, ID == 'SANDY M 120')$SD, fmt='%#.1f'),"ug C per g dry soil")


head(df.soil.MBC.combined.pvalues,2)
paste(sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'ORG O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'ORG O' & comparison == '120-0')$p.value, 
              fmt='%#.5f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'ORG O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'SANDY O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'SANDY O' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'SANDY O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'SANDY M' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'SANDY M' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.soil.MBC.combined.pvalues, ID == 'SANDY M' & comparison == '120-30')$p.value, 
              fmt='%#.3f'))






### Stats for substrate-derived MBC -----

anova <- aov(ug_substrate_derived_MBC_per_g_dry_soil~burn.trtmt.duration.seconds*substrate, 
             data=subset(df.sub.MBC))

summary(anova)
# Significant effect of substrate on soil-derived MBC & significant interaction 
#   between substrate and burn treatment

# Effect of substrate on substrate-derived MBC:
anova <- aov(ug_substrate_derived_MBC_per_g_dry_soil~substrate*burn.trtmt.duration.seconds, data = df.sub.MBC)
summary(anova)
# p < 0.001

### Creating dataframe of significant letters to add to plots. Method 2:
df.substrate.MBC.stats = data.frame(burn.trtmt.duration.seconds = 'duration', 
                                        site.type = 'substrate',
                                    substrate = 'substrate',
                                        horizon = 'k',
                                        max.substrate.MBC = 0,
                                        cld = 'cld')

# Now create list of significant letters for ggplot:
for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.sub.MBC, site.type == j & substrate == h & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(ug_substrate_derived_MBC_per_g_dry_soil~as.character(site) + burn.trtmt.duration.seconds, data = x)
      # Interaction between site.type and burn treatment 
      #    (site.type*burn.trtmt.duration.seconds) is not significant so it is 
      #    not including in the model, and we report results without the 
      #    interaction. 
      tukey <- TukeyHSD(anova)
      cld <- multcompLetters4(anova, tukey)
      x.cld <- group_by(subset(df.sub.MBC, 
                               site.type == j & substrate == h & horizon == k), 
                        burn.trtmt.duration.seconds) %>%
        subset(select = c(burn.trtmt.duration.seconds, site.type, horizon, substrate)) %>%
        unique()%>%
        dplyr::arrange((burn.trtmt.duration.seconds)) %>%
        dplyr::mutate(site.type = j)
      cld <- as.data.frame.list(cld$burn.trtmt.duration.seconds)
      cld <- cld[order(as.numeric(row.names(cld))),]
      x.cld$cld <- cld$Letters
      x.cld$max.substrate.MBC = max(x$ug_substrate_derived_MBC_per_g_dry_soil, na.rm = TRUE)
      df.substrate.MBC.stats <- rbind(df.substrate.MBC.stats, x.cld)
    }
  }
}

df.substrate.MBC.stats <- subset(df.substrate.MBC.stats, site.type != 'substrate')



# Now create list of p-values :
df.substrate.MBC.pvalues = data.frame(ID = 'ID',comparison = 'comparison',
                                    substrate = 'substrate', horizon = 'horizon', 
                                    p.value=0)


for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.sub.MBC, site.type == j & substrate == h & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(ug_substrate_derived_MBC_per_g_dry_soil~as.character(site) +burn.trtmt.duration.seconds, data = x)
      tukey <- TukeyHSD(anova)
      comparison = rownames(tukey$burn.trtmt.duration.seconds)
      p.value = tukey$burn.trtmt.duration.seconds[10:12]
      temp <- data.frame(ID = paste(j,k,h), substrate = h, horizon = k, 
                         comparison = comparison, p.value)
      df.substrate.MBC.pvalues <- rbind(df.substrate.MBC.pvalues, temp)
    }
  }
}

df.substrate.MBC.pvalues <- subset(df.substrate.MBC.pvalues, ID != 'ID')
df.substrate.MBC.pvalues

# stats for MBC
df.stats.substrate.MBC = data.frame(ID = 'ID',substrate = 'substrate', horizon = 'horizon',
                                    site.type='site',
                                    MEAN = 0, SD = 0)

for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      for (i in c(0,30,120)) {
        x <- subset(df.sub.MBC, site.type == j & horizon == k & 
                      burn.trtmt.duration.seconds == i & substrate == h)
        if (dim(x)[1] == 0){
          next
        }
        mean = mean(x$ug_substrate_derived_MBC_per_g_dry_soil, na.rm=TRUE)
        sd = sd(x$ug_substrate_derived_MBC_per_g_dry_soil, na.rm=TRUE)
        ID = paste(j,k,i,h)
        temp <- data.frame(substrate = h, horizon = k,site.type=j,
                           ID = ID, MEAN=mean, SD = sd)
        df.stats.substrate.MBC = rbind(df.stats.substrate.MBC, temp)
      }
    }
  }
}

df.stats.substrate.MBC


# Glucose:
paste(sprintf("Glucose: "),
      sprintf("Organic, unburned: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 0 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 0 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 30 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 30 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 120 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 120 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 0 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 0 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 30 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 30 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 120 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 120 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 0 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 0 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 30 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 30 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 120 glucose')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 120 glucose')$SD, fmt='%#.1f'),"ug C per g dry soil")

# Pine:
paste(sprintf("Pine: "),
      sprintf("Organic, unburned: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 0 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 0 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 30 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 30 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 120 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'ORG O 120 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 0 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 0 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 30 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 30 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 120 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY O 120 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 0 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 0 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 30 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 30 pine')$SD, fmt='%#.1f'),"ug C per g dry soil;",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 120 pine')$MEAN, fmt='%#.1f'), "\u00B1", 
      sprintf(subset(df.stats.substrate.MBC, ID == 'SANDY M 120 pine')$SD, fmt='%#.1f'),"ug C per g dry soil")

head(df.substrate.MBC.stats,2)

# Glucose: 
paste(sprintf("Glucose: "),
      sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'ORG O glucose' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'ORG O glucose' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'ORG O glucose' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY O glucose' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY O glucose' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY O glucose' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY M glucose' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY M glucose' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY M glucose' & comparison == '120-30')$p.value, 
              fmt='%#.3f'))

# Pine  : 
paste(sprintf("Pine: "),
      sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'ORG O pine' & comparison == '30-0')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'ORG O pine' & comparison == '120-0')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'ORG O pine' & comparison == '120-30')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY O pine' & comparison == '30-0')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY O pine' & comparison == '120-0')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY O pine' & comparison == '120-30')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY M pine' & comparison == '30-0')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY M pine' & comparison == '120-0')$p.value, 
              fmt='%#.2f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.substrate.MBC.pvalues, ID == 'SANDY M pine' & comparison == '120-30')$p.value, 
              fmt='%#.2f'))


### PLOT: Figure 3: Substrate-derived MBC ------------------------------------------------
df.substrate.MBC.stats$horizon <- factor(df.substrate.MBC.stats$horizon, levels = c('O','M'))


## glucose plots
p_substrate_MBC_glucose_1 = ggplot(subset(df.sub.MBC, substrate == 'glucose' & site.type=='ORG' & horizon == 'O'), 
                            aes(x=(burn.trtmt.duration.seconds), 
                                y=ug_substrate_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.substrate.MBC.stats, substrate=='glucose' & site.type=='ORG' & horizon == 'O'), 
            aes(y=max.substrate.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Glucose-derived~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') + 
  ylim(0,max(subset(df.sub.MBC, substrate=='glucose')$ug_substrate_derived_MBC_per_g_dry_soil, na.rm = TRUE)*1.25)


p_substrate_MBC_glucose_2 = ggplot(subset(df.sub.MBC, substrate == 'glucose' & site.type=='SANDY' & horizon == 'O'), 
                                   aes(x=(burn.trtmt.duration.seconds), 
                                       y=ug_substrate_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  geom_text(data=subset(df.substrate.MBC.stats, substrate=='glucose' & site.type=='SANDY' & horizon == 'O'), 
            aes(y=max.substrate.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Glucose-derived~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') +
  ylim(0,max(subset(df.sub.MBC, substrate=='glucose')$ug_substrate_derived_MBC_per_g_dry_soil, na.rm = TRUE)*1.25)


p_substrate_MBC_glucose_3 = ggplot(subset(df.sub.MBC, substrate == 'glucose' & site.type=='SANDY' & horizon == 'M'), 
                                   aes(x=(burn.trtmt.duration.seconds), 
                                       y=ug_substrate_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  # geom_text(data=subset(df.substrate.MBC.stats, substrate=='glucose'& site.type=='SANDY' & horizon == 'M'), 
            # aes(y=max.substrate.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Glucose-derived~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '')


cowplot::plot_grid(p_substrate_MBC_glucose_1, p_substrate_MBC_glucose_2, p_substrate_MBC_glucose_3, ncol=3)



## Pine plots
p_substrate_MBC_pine_1 = ggplot(subset(df.sub.MBC, substrate == 'pine' & site.type=='ORG' & horizon == 'O'), 
                                aes(x=(burn.trtmt.duration.seconds), 
                                    y=ug_substrate_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  # geom_text(data=subset(df.substrate.MBC.stats, substrate=='pine'), 
  # aes(y=max.substrate.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Pine-derived~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') +
  ylim(-40, 120)

p_substrate_MBC_pine_2 = ggplot(subset(df.sub.MBC, substrate == 'pine' & site.type=='SANDY' & horizon == 'O'), 
                                aes(x=(burn.trtmt.duration.seconds), 
                                    y=ug_substrate_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  # geom_text(data=subset(df.substrate.MBC.stats, substrate=='pine'), 
  # aes(y=max.substrate.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Pine-derived~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') +
  ylim(-15, 45)

p_substrate_MBC_pine_3 = ggplot(subset(df.sub.MBC, substrate == 'pine' & site.type=='SANDY' & horizon == 'M'), 
                                 aes(x=(burn.trtmt.duration.seconds), 
                                     y=ug_substrate_derived_MBC_per_g_dry_soil)) +
  geom_boxplot(aes(fill=burn.trtmt.duration.seconds), alpha=0.6) + 
  # geom_jitter(aes(), color = 'grey20', width = 0.1, alpha=0.7, size=3)+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  # geom_text(data=subset(df.substrate.MBC.stats, substrate=='pine'), 
            # aes(y=max.substrate.MBC+0.1, label = cld), vjust = -1, size=4) +
  geom_hline(yintercept = 0) + 
  labs(x='Burn treatment exposure (seconds)',
       y=expression(atop(Pine-derived~microbial~biomass, (mu*g~C~per~g~dry~soil))),
       fill = 'Burn treatment exposure\n(seconds)') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) + 
  theme(legend.position = '') +
  ylim(-2, 6)

cowplot::plot_grid(p_substrate_MBC_pine_1, p_substrate_MBC_pine_2, p_substrate_MBC_pine_3, ncol=3)






### EQ 10 & 11: Calculate CUE and metabolic quotient, qCO2 ---------------------

#     CUE = substrate-derived MBC / (substrate-derived MBC + substrate-derived CO2)
#     qCO2 = substrate-derived CO2 / substrate-derived MBC

head(df.CO2.clean,2)


df.glucose.merge <- subset(df.CO2.clean, amendment =='13CG') %>%
  subset(select = c(sample.id, amendment,
                    cumulative.ug.substrate.derived.CO2.per.g.dry.soil,
                    cumulative.ug.soil.derived.CO2.per.g.dry.soil)) %>%
  merge(df.glucose, by = c('sample.id','amendment')) %>%
  dplyr::mutate(CUE = ug_substrate_derived_MBC_per_g_dry_soil/(ug_substrate_derived_MBC_per_g_dry_soil+cumulative.ug.substrate.derived.CO2.per.g.dry.soil)) %>%
  dplyr::mutate(qCO2 = cumulative.ug.substrate.derived.CO2.per.g.dry.soil/ug_substrate_derived_MBC_per_g_dry_soil,
                qCO2_soil_derived = cumulative.ug.soil.derived.CO2.per.g.dry.soil/ug_soil_derived_MBC_per_g_dry_soil) %>%
  # dplyr::mutate(alternative_CUE = alternative_ug_substrate_derived_MBC_per_g_dry_soil/(alternative_ug_substrate_derived_MBC_per_g_dry_soil+cumulative.ug.substrate.derived.CO2.per.g.dry.soil))
  dplyr::mutate(substrate = 'glucose')



df.pine.merge <- subset(df.CO2.clean, amendment =='13CP') %>%
  subset(select = c(sample.id, amendment,
                    cumulative.ug.substrate.derived.CO2.per.g.dry.soil,
                    cumulative.ug.soil.derived.CO2.per.g.dry.soil)) %>%
  merge(df.pine, by = c('sample.id','amendment')) %>%
  dplyr::mutate(CUE = ug_substrate_derived_MBC_per_g_dry_soil/(ug_substrate_derived_MBC_per_g_dry_soil+cumulative.ug.substrate.derived.CO2.per.g.dry.soil)) %>%
  dplyr::mutate(qCO2 = cumulative.ug.substrate.derived.CO2.per.g.dry.soil/ug_substrate_derived_MBC_per_g_dry_soil,
                qCO2_soil_derived = cumulative.ug.soil.derived.CO2.per.g.dry.soil/ug_soil_derived_MBC_per_g_dry_soil) %>%
  # dplyr::mutate(alternative_CUE = alternative_ug_substrate_derived_MBC_per_g_dry_soil/(alternative_ug_substrate_derived_MBC_per_g_dry_soil+cumulative.ug.substrate.derived.CO2.per.g.dry.soil))
  dplyr::mutate(substrate = 'pine')

df.CUE <- rbind(df.pine.merge, df.glucose.merge)

### Add burn data --------------------------------------------------------------

df.CUE <- merge(df.CUE, subset(df.meta, select = c(core.id, burn.trtmt.duration.seconds, vegetation)))

df.CUE$burn.trtmt.duration.seconds <- factor(df.CUE$burn.trtmt.duration.seconds, levels = c(0,30,120))
df.CUE$horizon <- factor(df.CUE$horizon, levels = c('O','M'))

dim(df.CUE)
df.CUE <- unique(df.CUE)


df.CUE$horizon <- factor(df.CUE$horizon, levels = c('O','M'))


### PLOT: Figure 5. CUE ------------------------------------------------------------------

# Remove samples where substrate-derived MBC is less than 0
df.CUE.prune <- df.CUE %>%
  subset(ug_substrate_derived_MBC_per_g_dry_soil >0) # removes 10 samples


# # Creating dataframe of significant letters to add to plots. Method 2:
df.CUE.stats = data.frame(burn.trtmt.duration.seconds = 'duration', 
                       site.type = 'soil', 
                       horizon = 'k',
                       substrate = 'incubation.period.days',
                       max.CUE = 0,
                       cld = 'cld')

# Now create list of significant letters for ggplot:
for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CUE.prune, site.type == j & substrate == h & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(CUE~as.character(site)+burn.trtmt.duration.seconds, data = x)
      # Interaction between site.type and burn treatment 
      #    (site.type*burn.trtmt.duration.seconds) is not significant so it is 
      #    not including in the model, and we report results without the 
      #    interaction. 
      tukey <- TukeyHSD(anova)
      cld <- multcompLetters4(anova, tukey)
      x.cld <- group_by(subset(df.CUE.prune, 
                               site.type == j & 
                                 horizon == k &
                              substrate == h), burn.trtmt.duration.seconds) %>%
        subset(select = c(burn.trtmt.duration.seconds, site.type, substrate, horizon)) %>%
        unique()%>%
        dplyr::arrange((burn.trtmt.duration.seconds)) %>%
        dplyr::mutate(substrate = h, site.type = j)
      
      cld <- as.data.frame.list(cld$burn.trtmt.duration.seconds) 
      cld <- cld[order(as.numeric(row.names(cld))),]
      x.cld$cld <- cld$Letters
      x.cld$max.CUE = max(x$CUE, na.rm = TRUE)
      df.CUE.stats <- rbind(df.CUE.stats, x.cld)
    }
  }
}
  
df.CUE.stats <- subset(df.CUE.stats, burn.trtmt.duration.seconds != 'site.type')
  
  

  
# Now create list of p-values :
df.CUE.pvalues = data.frame(ID = 'ID',comparison = 'comparison',p.value=0)
  
for (j in c('ORG', 'SANDY')) {
  for (h in c('glucose','pine')) {
    for (k in c('O', 'M')) {
      x <- subset(df.CUE.prune, site.type == j & substrate == h & horizon == k)
      if (dim(x)[1] == 0){
        next
      }
      anova <- aov(CUE~as.character(site)+burn.trtmt.duration.seconds, data = x)
      tukey <- TukeyHSD(anova)
      comparison = rownames(tukey$burn.trtmt.duration.seconds)
      p.value = tukey$burn.trtmt.duration.seconds[10:12]
      temp <- data.frame(ID = paste(j,h,k), comparison = comparison, p.value)
      df.CUE.pvalues <- rbind(df.CUE.pvalues, temp)
    }
  }
}

df.CUE.pvalues <- subset(df.CUE.pvalues, ID != 'ID')
df.CUE.stats$horizon <- factor(df.CUE.stats$horizon, level=c('O','M'))

df.CUE.stats.sub <- subset(df.CUE.stats, substrate == 'glucose' & horizon != 'M') %>%
  rbind(subset(df.CUE.stats, substrate == 'pine' & site.type == 'SANDY' & horizon != 'M'))


pCUE_pine = ggplot(subset(df.CUE.prune, substrate == 'pine'),
       # aes(x=site.type, y = CUE, color = as.character(burn.trtmt.duration.seconds), shape=horizon )) +
       aes(x=burn.trtmt.duration.seconds,y=CUE, 
           fill=burn.trtmt.duration.seconds))+ 
  geom_boxplot(alpha=0.6)+
  # geom_jitter(aes(color=burn.trtmt.duration.seconds),size=3, width=0.2,alpha=0.7)+
  # geom_jitter(aes(shape=horizon,color = burn.trtmt.duration.seconds), width = 0.2, alpha=0.7,size=2) + 
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type=SiteType_labels)) +
  scale_colour_manual(values = c('black','darksalmon','darkred'),
                      breaks= c('0', '30','120'))+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                      breaks= c('0', '30','120'))+
  # geom_hline(yintercept=1) + 
  # geom_hline(yintercept=0) + 
  geom_text(data=subset(df.CUE.stats.sub, substrate == 'pine' ),
            aes(y=max.CUE+0.1, label = cld), vjust = 0, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       y = 'CUE',
       fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Pine amendment') + 
  theme(axis.title = element_text(size=12),
        axis.text=element_text(size=12), 
        strip.text = element_text(size=12))+
  theme(legend.position='')+
  ylim(0,1)
  # ylim(min(subset(df.CUE.prune, substrate == 'pine')$CUE, na.rm = TRUE), 
       # max(subset(df.CUE.prune, substrate == 'pine')$CUE, na.rm = TRUE)+0.25)


pCUE_glucose = ggplot(subset(df.CUE.prune, substrate == 'glucose'),
                      # aes(x=site.type, y = CUE, color = as.character(burn.trtmt.duration.seconds), shape=horizon )) +
                      aes(x=burn.trtmt.duration.seconds,y=CUE, 
                          fill=burn.trtmt.duration.seconds))+ 
  geom_boxplot(alpha=0.6)+
  # geom_jitter(aes(color=burn.trtmt.duration.seconds),size=3, width=0.2,alpha=0.7)+
  # geom_jitter(alpha=0.7, color='grey')+
  # geom_jitter(aes(shape=horizon,color = burn.trtmt.duration.seconds), width = 0.2, alpha=0.7,size=2) + 
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free', 
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  scale_colour_manual(values = c('black','darksalmon','darkred'),
                      breaks= c('0', '30','120'))+
  scale_fill_manual(values = c('black','darksalmon','darkred'),
                    breaks= c('0', '30','120'))+
  # geom_hline(yintercept=1) + 
  # geom_hline(yintercept=0) + 
  geom_text(data=subset(df.CUE.stats.sub, substrate == 'glucose' ),
            aes(y=max.CUE+0.1, label = cld), vjust = 0, size=4) +
  labs(x='Burn treatment exposure (seconds)',
       y = 'CUE',
       fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12),
        axis.text=element_text(size=12), 
        strip.text = element_text(size=12))+
  theme(legend.position='')+
  ylim(0,1)

cowplot::plot_grid(pCUE_pine, pCUE_glucose, ncol=1)


### Figure caption details:
paste(sprintf("Figure [CUE]. Caption details - CUE for all samples calculated after removing (8) samples with MBC < 0"),
      sprintf("and loss of 02-02-O for lack of data (sandy sites, O horizon, 120 s burn treatment)"),
      sprintf("n=98 total; n=6 for each boxplot (6 cores for each burn trtmt, with exceptions);"),
      sprintf("stats = ANOVA with TukeyHSD test, letters indicated significant differences at p<0.05 as cutoff"))


### CUE data written results:
df.stats.CUE.plot = data.frame(ID = 'ID',substrate='x',MEAN = 0, SD = 0)

for (j in c('ORG', 'SANDY')) {
  for (k in c('O', 'M')) {
    for (i in c(0,30,120)) {
      for (m in c('glucose','pine')) {
        x <- subset(df.CUE.prune, site.type == j & horizon == k & 
                      burn.trtmt.duration.seconds == i & substrate==m)
        if (dim(x)[1] == 0){
          next
        }
        mean = mean(x$CUE, na.rm=TRUE)
        sd = sd(x$CUE, na.rm=TRUE)
        ID = paste(j,k,i)
        temp <- data.frame(ID = ID, substrate=m, MEAN=mean, SD = sd)
        df.stats.CUE.plot = rbind(df.stats.CUE.plot, temp)
      }
    }
  }
}

df.stats.CUE.plot

# GLUCOSE:
df.stats.CUE.glucose <- subset(df.stats.CUE.plot, substrate=='glucose')

paste(sprintf("GLUCOSE: Organic, unburned: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'ORG O 0')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'ORG O 0')$SD, fmt='%#.2f'),";", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'ORG O 30')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'ORG O 30')$SD, fmt='%#.2f'),";", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'ORG O 120')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'ORG O 120')$SD, fmt='%#.2f'),";",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY O 0')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY O 0')$SD, fmt='%#.2f'),";",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY O 30')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY O 30')$SD, fmt='%#.2f'),";",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY O 120')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY O 120')$SD, fmt='%#.2f'),";",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY M 0')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY M 0')$SD, fmt='%#.2f'),";",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY M 30')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY M 30')$SD, fmt='%#.2f'),";",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY M 120')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.glucose, ID == 'SANDY M 120')$SD, fmt='%#.2f'))

# PINE:
df.stats.CUE.pine <- subset(df.stats.CUE.plot, substrate=='pine')

paste(sprintf("PINE: Organic, unburned: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'ORG O 0')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'ORG O 0')$SD, fmt='%#.2f'),";", 
      sprintf("Organic, 30 s: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'ORG O 30')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'ORG O 30')$SD, fmt='%#.2f'),";", 
      sprintf("Organic, 120s: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'ORG O 120')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'ORG O 120')$SD, fmt='%#.2f'),";",
      sprintf("Sandy O hor, unburned: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'SANDY O 0')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'SANDY O 0')$SD, fmt='%#.2f'),";",
      sprintf("Sandy O hor, 30 s: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'SANDY O 30')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'SANDY O 30')$SD, fmt='%#.2f'),";",
      sprintf("Sandy O hor, 120s: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'SANDY O 120')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'SANDY O 120')$SD, fmt='%#.2f'),";",
      sprintf("Sandy M hor, unburned: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'SANDY M 0')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'SANDY M 0')$SD, fmt='%#.2f'),";",
      sprintf("Sandy M hor, 30 s: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'SANDY M 30')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'SANDY M 30')$SD, fmt='%#.2f'),";",
      sprintf("Sandy M hor, 120s: mean="),sprintf(subset(df.stats.CUE.pine, ID == 'SANDY M 120')$MEAN, fmt='%#.2f'), "\u00B1", 
      sprintf(subset(df.stats.CUE.pine, ID == 'SANDY M 120')$SD, fmt='%#.2f'))


# df.CUE.pvalues

# Glucose: 
paste(sprintf("Glucose: "),
      sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'ORG glucose O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'ORG glucose O' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'ORG glucose O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY glucose O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY glucose O' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY glucose O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY glucose M' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY glucose M' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY glucose M' & comparison == '120-30')$p.value, 
              fmt='%#.3f'))

# Pine  : 
paste(sprintf("Pine: "),
      sprintf("Histosol, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'ORG pine O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'ORG pine O' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Histosol, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'ORG pine O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY pine O' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY pine O' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol O, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY pine O' & comparison == '120-30')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 30 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY pine M' & comparison == '30-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, unburned vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY pine M' & comparison == '120-0')$p.value, 
              fmt='%#.3f'), ';',
      sprintf("Gleysol M, 30 vs. 120 s: p-value="),
      sprintf(subset(df.CUE.pvalues, ID == 'SANDY pine M' & comparison == '120-30')$p.value, 
              fmt='%#.3e'))

### PLOT: Figure 6A. Max Temp vs. CUE_glucose, O horizon ----------------------------------

# Import raw data
df.Temp.filtered <- read.csv('../data/burn-simulations/Output-files/Compiled-decimated-temperature-all-cores.csv') %>%
  # Filter time series remove instances of temperature increases or decreases of
  #   >100 units
  group_by(core.id, thermocouple_position) %>%
  mutate(NeighborLag = c(NA,diff(temperature.C,lag=1)),
         NeighborLead = c(diff(temperature.C,lead=1),NA)) %>%
  subset(!(abs(NeighborLag) >100) & !(abs(NeighborLead) >100)) %>%
  # Take one data point per minute (to make it faster to graph)
  filter((run.time.sec/60) %% 1 == 0) %>%
  subset(run.time.sec/3600 < 7)

# import vegetation and burn treatment information
df.duration <- read.csv('../data/burn-simulations/core-prep.csv')

df.duration <- df.duration %>%
  subset(select = c(core.id, vegetation, duration.s)) %>%
  mutate(full.id = core.id) %>%
  separate(core.id, c('year','project','site','core'), sep = '-', remove = TRUE) %>%
  unite('core.id', c(site,core), sep='-', remove = TRUE)

# merge temperature and veg.+ burn treatment data frames:
df.Temp <- merge(df.Temp.filtered, df.duration, by = 'core.id') 

# Re order the thermocouples:
df.Temp$thermocouple_position <- factor(df.Temp$thermocouple_position, levels = c('U','M','L'))


### Organic sites ###
# Max O horizon temp:
df.Temp.max <- df.Temp %>%
  group_by(core.id, thermocouple_position, duration.s) %>%
  mutate(max.temp.C = max(temperature.C)) %>%
  mutate(core.id = full.id) %>%
  subset(select = c(core.id,
                    thermocouple_position,
                    max.temp.C)) %>% 
  unique() 


df.Temp.unburned <- subset(df.CUE, 
                           burn.trtmt.duration.seconds==0 &
                           ug_substrate_derived_MBC_per_g_dry_soil >0) %>%
  subset(select = c(core.id)) %>%
  dplyr::mutate(max.temp.C = 25, 
                thermocouple_position = 'U') 

df.Temp.full <- rbind(df.Temp.max, df.Temp.unburned)

head(df.Temp.full,2)


# Combine T and CUE dataframes:
df.T.CUE <- subset(df.Temp.full, thermocouple_position == 'U') %>%
  # mutate(burn.trtmt.duration.seconds = duration.s) %>%
  subset(select =c(core.id, max.temp.C)) %>%
  merge(subset(df.CUE.prune,ug_substrate_derived_MBC_per_g_dry_soil >0), 
        by = c('core.id'))

# Clean up dataframe for plotting
df.T.CUE$burn.trtmt.duration.seconds <- factor(df.T.CUE$burn.trtmt.duration.seconds,
                                            levels = c('0','30','120'))



## TRYING WITH NLS.LM FUNCTION
getPred = function(params, xx) params$a*exp(-params$b * xx) + params$c
# Predicted CUE = a * exp(-b*max.temp.C)

residFun <- function(params, CUE, max.temp.C) CUE - getPred(params,max.temp.C) # function(params, observed, xx) observed - getPred(p,xx)
# Observation = CUE, getPred = predicted CUE; nls.lm and ResidFun is trying to minimilize the sum of squares between the two vectors. 



## Fitting non-linear regression: exponential decay, GLUCOSE:
df.T.glucose = subset(df.T.CUE, horizon == 'O' & substrate=='glucose')

anova <- aov(CUE~max.temp.C*burn.trtmt.duration.seconds*site.type, data = subset(df.T.glucose))
summary(anova) # No significant difference in CUE across site types

# starting params, glucose: a=0.6, b=0.8, c=0.3
fit = nls.lm(par=list(a=0.6, b=0.05, c=0.3), 
             fn = residFun, 
             max.temp.C = df.T.glucose$max.temp.C,
             CUE = df.T.glucose$CUE,
             lower = c(0,0,0), upper = c(Inf, Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
glucose_a = summary(fit)$coefficients[[1]]
glucose_b = summary(fit)$coefficients[[2]]
glucose_c = summary(fit)$coefficients[[3]] # parameter c sets the minimum for CUE. 
  # i.e. above some maximum T, CUE can't decrease any further
  # Is this valid or would we expect CUE to decrease to 0 at extremely high T?
  # Well, CUE can't equal 0. If all microbes are dead, then respiration = 0 and MBC = 0
  # and CUE doesn't "exist". 


for (i in 1:nrow(df.T.glucose)) {
  # df.T.glucose$NonLinearFit[i] = glucose_a *exp(-glucose_b * df.T.glucose$max.temp.C[i]) + glucose_c
  # df.T.glucose$NonLinearFit[i] = glucose_a - glucose_b * log(df.T.glucose$max.temp.C[i]) 
  df.T.glucose$CUE.fit <- getPred(as.list(coef(fit)),df.T.glucose$max.temp.C)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.T.glucose)
summary(df.lm)
# summary(df.lm)$r.squared
T.glucose.R2adj = summary(df.lm)$adj.r.squared
T.glucose.p.value = summary(df.lm)$coefficients[8]

# df.fit <- data.frame('max.temp.C' = seq(from = 25, to = 533, by = 50))
# 
# for (i in 1:nrow(df.T.glucose)) { 
#   df.fit$CUE.fit <- getPred(as.list(coef(fit)),df.fit$max.temp.C)
# }
# 
# ggplot(df.fit, aes(x=max.temp.C, y= CUE.fit)) + 
#   geom_point()

## PLOT:
p_CUEgl_max_temp_O = ggplot(subset(df.T.glucose, 
                                   substrate == 'glucose' & horizon == 'O' ), 
                     aes(x=(max.temp.C), 
                         y=CUE,
                         color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  geom_line(aes(x=max.temp.C,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = expression(paste("Maximum temperature (",degree*C,") ")), 
       y = 'CUE',
       color='Burn treatment\nexposure',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 
p_CUEgl_max_temp_O

paste(sprintf('Glucose: y = a * exp(-b * max T) + c; where a = '), 
      sprintf(glucose_a, fmt='%#.2f'), ", b = ", 
      sprintf(glucose_b, fmt='%#.2f'), ", and c =  ",
      sprintf(glucose_c, fmt='%#.2f'), ", R2adj. = ", 
      sprintf(T.glucose.R2adj, fmt='%#.2f'), ' and p = ',
      sprintf(T.glucose.p.value, fmt='%#.1e'))




### PLOT: Figure 6A. Max Temp vs. CUE_glucose, Mineral horizon ----------------------------

# Create sample plot for mineral soil using lower thermocouple:
df.Temp.unburned <- subset(df.CUE, 
                           horizon == 'M' & 
                           burn.trtmt.duration.seconds==0 &
                             ug_substrate_derived_MBC_per_g_dry_soil >0) %>%
  subset(select = c(core.id)) %>%
  dplyr::mutate(max.temp.C = 25, 
                thermocouple_position = 'L') 

df.Temp.full <- rbind(df.Temp.max, df.Temp.unburned)

head(df.Temp.full,2)

# Combine T and CUE dataframes:
df.T.CUE <- subset(df.Temp.full, thermocouple_position == 'L') %>%
  # mutate(burn.trtmt.duration.seconds = duration.s) %>%
  subset(select =c(core.id, max.temp.C)) %>%
  merge(subset(df.CUE.prune,ug_substrate_derived_MBC_per_g_dry_soil >0), 
        by = c('core.id'))

# Clean up dataframe for plotting
df.T.CUE$burn.trtmt.duration.seconds <- factor(df.T.CUE$burn.trtmt.duration.seconds,
                                               levels = c('0','30','120'))

## Fitting non-linear regression: exponential decay, GLUCOSE:
df.T.glucose = subset(df.T.CUE, horizon == 'M' & substrate=='glucose')

# starting params, glucose: a=0.6, b=0.8, c=0.3
fit = nls.lm(par=list(a=0.6, b=0.05, c=0.3), 
             fn = residFun, 
             max.temp.C = df.T.glucose$max.temp.C,
             CUE = df.T.glucose$CUE,
             lower = c(0,0,0), upper = c(Inf, Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
glucose_a = summary(fit)$coefficients[[1]]
glucose_b = summary(fit)$coefficients[[2]]
glucose_c = summary(fit)$coefficients[[3]] 

for (i in 1:nrow(df.T.glucose)) {
  df.T.glucose$CUE.fit <- getPred(as.list(coef(fit)),df.T.glucose$max.temp.C)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.T.glucose)
summary(df.lm)
# summary(df.lm)$r.squared
T.glucose.R2adj = summary(df.lm)$adj.r.squared
T.glucose.p.value = summary(df.lm)$coefficients[8]

## PLOT:
p_CUEgl_max_temp_M = ggplot(subset(df.T.glucose, substrate == 'glucose' & horizon == 'M' ), 
       aes(x=(max.temp.C), 
           y=CUE,
           color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  geom_line(aes(x=max.temp.C,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = expression(paste("Maximum temperature (",degree*C,") ")), 
       y = 'CUE',
       color='Burn treatment\nexposure',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 

p_CUEgl_max_temp_M

paste(sprintf('Glucose: y = a * exp(-b * max T) + c; where a = '), 
      sprintf(glucose_a, fmt='%#.2f'), ", b = ", 
      sprintf(glucose_b, fmt='%#.2f'), ", and c =  ",
      sprintf(glucose_c, fmt='%#.2f'), ", R2adj. = ", 
      sprintf(T.glucose.R2adj, fmt='%#.2f'), ' and p = ',
      sprintf(T.glucose.p.value, fmt='%#.2f'))



### PLOT: Figure S1. pH vs. CUE -----------------------------------------------------------
# head(df.CUE,2)
 
df.pH <- read.csv('../data/soil-properties/pH.csv')
# head(df.pH,2)
colnames(df.pH)[1] <- 'core.id'

df.pH.merge <- df.pH %>%
  subset(select = c(core.id, horizon, pH)) %>%
  merge(df.CUE, by = c('core.id', 'horizon')) %>%
  # Remove samples where substrate-derived MBC is less than 0
  subset(ug_substrate_derived_MBC_per_g_dry_soil >0) # removes 10 samples

# head(df.pH.merge,2)

# Effect of pH on CUE:
# Model with interaction btw pH and site; only plotting O horizons here
anova <- aov(CUE~pH*substrate*site.type, data = df.pH.merge)
summary(anova)


## TRYING WITH NLS.LM FUNCTION
getPred = function(params, xx) params$a*exp(-params$b * xx)
# Predicted CUE = a * exp(-b*max.temp.C)

residFun <- function(params, CUE, pH) CUE - getPred(params,pH) # function(params, observed, xx) observed - getPred(p,xx)
# Observation = CUE, getPred = predicted CUE; nls.lm and ResidFun is trying to minimilize the sum of squares between the two vectors. 

### Glucose, HISTOSOL:
## Fitting non-linear regression: exponential decay, GLUCOSE:
df.pH.glucose = subset(df.pH.merge, horizon == 'O' & substrate=='glucose' & site.type == 'ORG')

anova <- aov(CUE~pH*burn.trtmt.duration.seconds, data = subset(df.pH.merge, horizon == 'O' & substrate=='glucose'))
summary(anova) # No significant difference in CUE across site types

# starting params, glucose: a= , b= , c= 
fit = nls.lm(par=list(a=1, b=1), 
             fn = residFun, 
             pH = df.pH.glucose$pH,
             CUE = df.pH.glucose$CUE,
             lower = c(0,0), upper = c(Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
glucose_a = summary(fit)$coefficients[[1]]
glucose_b = summary(fit)$coefficients[[2]]


for (i in 1:nrow(df.pH.glucose)) {
  # df.pH.glucose$NonLinearFit[i] = glucose_a *exp(-glucose_b * df.pH.glucose$max.temp.C[i]) + glucose_c
  # df.pH.glucose$NonLinearFit[i] = glucose_a - glucose_b * log(df.pH.glucose$max.temp.C[i]) 
  df.pH.glucose$CUE.fit <- getPred(as.list(coef(fit)),df.pH.glucose$pH)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.pH.glucose)
summary(df.lm)
# summary(df.lm)$r.squared
pH.glucose.R2adj = summary(df.lm)$adj.r.squared
pH.glucose.p.value = summary(df.lm)$coefficients[8]


## PLOT:
p_pHvCUE_Glucose_Hist = ggplot(subset(df.pH.glucose, substrate == 'glucose' & horizon == 'O' & site.type == 'ORG'), 
       aes(x=(pH), 
           y=CUE,
           color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  geom_line(aes(x=pH,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = 'Soil pH', 
       y = 'CUE',
       color='Burn treatment\nexposure (seconds)',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 


paste(sprintf('Glucose: y = a * exp(-b * max T) ; where a = '), 
      sprintf(glucose_a, fmt='%#.2f'), ", b = ", 
      sprintf(glucose_b, fmt='%#.2f'), ", R2adj. = ", 
      sprintf(pH.glucose.R2adj, fmt='%#.2f'), ' and p = ',
      sprintf(pH.glucose.p.value, fmt='%#.1e'))


### GLUCOSE, GLEYSOL, O horizon:
## Fitting non-linear regression: exponential decay, GLUCOSE:
df.pH.glucose.sandy = subset(df.pH.merge, horizon == 'O' & substrate=='glucose' & site.type == 'SANDY')

# starting params, glucose: a= , b= , c= 
fit = nls.lm(par=list(a=1, b=1), 
             fn = residFun, 
             pH = df.pH.glucose.sandy$pH,
             CUE = df.pH.glucose.sandy$CUE,
             lower = c(0,0), upper = c(Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
glucose_a = summary(fit)$coefficients[[1]]
glucose_b = summary(fit)$coefficients[[2]]


for (i in 1:nrow(df.pH.glucose.sandy)) {
  # df.pH.glucose$NonLinearFit[i] = glucose_a *exp(-glucose_b * df.pH.glucose$max.temp.C[i]) + glucose_c
  # df.pH.glucose$NonLinearFit[i] = glucose_a - glucose_b * log(df.pH.glucose$max.temp.C[i]) 
  df.pH.glucose.sandy$CUE.fit <- getPred(as.list(coef(fit)),df.pH.glucose.sandy$pH)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.pH.glucose)
summary(df.lm)
# summary(df.lm)$r.squared
pH.glucose.R2adj = summary(df.lm)$adj.r.squared
pH.glucose.p.value = summary(df.lm)$coefficients[8]


## PLOT:
p_pHvCUE_Glucose_Gley = ggplot(subset(df.pH.glucose.sandy, substrate == 'glucose' & horizon == 'O' & site.type == 'SANDY'), 
       aes(x=(pH), 
           y=CUE,
           color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_line(aes(x=pH,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~site.type+horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = 'Soil pH', 
       y = 'CUE',
       color='Burn treatment\nexposure (seconds)',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 


paste(sprintf('Glucose: y = a * exp(-b * max T) ; where a = '), 
      sprintf(glucose_a, fmt='%#.2f'), ", b = ", 
      sprintf(glucose_b, fmt='%#.2f'), ", R2adj. = ", 
      sprintf(pH.glucose.R2adj, fmt='%#.2f'), ' and p = ',
      sprintf(pH.glucose.p.value, fmt='%#.1e'))


cowplot::plot_grid(p_pHvCUE_Glucose_Hist,p_pHvCUE_Glucose_Gley, ncol=2)







## Fitting non-linear regression: exponential decay, PINE:
df.pH.pine = subset(df.pH.merge, horizon == 'O' & substrate=='pine' & site.type == 'SANDY')

anova <- aov(CUE~pH*burn.trtmt.duration.seconds, data = subset(df.pH.pine))
summary(anova) # No significant difference in CUE across site types

fit = nls.lm(par=list(a=0.6, b=0.05), 
             fn = residFun, 
             pH = df.pH.pine$pH,
             CUE = df.pH.pine$CUE,
             lower = c(0,0), upper = c(Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
pine_a = summary(fit)$coefficients[[1]]
pine_b = summary(fit)$coefficients[[2]]


for (i in 1:nrow(df.pH.pine)) {
  df.pH.pine$CUE.fit <- getPred(as.list(coef(fit)),df.pH.pine$pH)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.pH.pine)
summary(df.lm)
# summary(df.lm)$r.squared
pH.pine.R2adj = summary(df.lm)$adj.r.squared
pH.pine.p.value = summary(df.lm)$coefficients[8]


## PLOT:
ggplot(subset(df.pH.pine, substrate == 'pine' & horizon == 'O' & site.type == 'SANDY'), 
       aes(x=(pH), 
           y=CUE,
           color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  # geom_line(aes(x=pH,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~site.type +horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = 'Soil pH', 
       y = 'CUE',
       color='Burn treatment\nexposure (seconds)',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Pine amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 





### PLOT: Figure 6B. gene copy number vs. CUE_glucose -------------------------------------

# Import data for weighted mean 16S gene copy number: 
df.WtMean <-  read.csv('../data/sequence/rrnDB_estimates/Feb2025/Weighted-mean-copy-number-FireMod22.csv')
# head(df.WtMean,2)

# Merge with CUE dataset using just O horizon copy number
df.WM.merge <- df.WtMean %>%
  subset(respiration.incubation.length.days==14) %>%
  subset(select = c(core.id, horizon, WtMeanCopyNum)) %>%
  merge(df.CUE, by = c('core.id', 'horizon')) %>%
  # Remove samples where substrate-derived MBC is less than 0
  subset(ug_substrate_derived_MBC_per_g_dry_soil >0) # removes 10 samples


# Effect of WtMeanCopyNum on CUE:
# Model with interaction btw WtMeanCopyNum and site; only plotting O horizons here
anova <- aov(CUE~WtMeanCopyNum*substrate*site.type, data = df.WM.merge)
summary(anova)



## TRYING WITH NLS.LM FUNCTION
getPred = function(params, xx) params$a*exp(-params$b * xx)
# Predicted CUE = a * exp(-b*max.temp.C)

residFun <- function(params, CUE, WtMeanCopyNum) CUE - getPred(params,WtMeanCopyNum) 
# function(params, observed, xx) observed - getPred(p,xx)
# Observation = CUE, getPred = predicted CUE; nls.lm and ResidFun is trying to minimilize the sum of squares between the two vectors. 



## Fitting non-linear regression: exponential decay, GLUCOSE:
df.WM.glucose = subset(df.WM.merge, horizon == 'O' & substrate=='glucose')

anova <- aov(CUE~WtMeanCopyNum*site.type*burn.trtmt.duration.seconds, data = subset(df.WM.glucose))
summary(anova) # No significant difference in CUE across site types

# starting params, glucose: a= , b= , c= 
fit = nls.lm(par=list(a=1, b=1), 
             fn = residFun, 
             WtMeanCopyNum = df.WM.glucose$WtMeanCopyNum,
             CUE = df.WM.glucose$CUE,
             lower = c(0,0), upper = c(Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
glucose_a = summary(fit)$coefficients[[1]]
glucose_b = summary(fit)$coefficients[[2]]


for (i in 1:nrow(df.WM.glucose)) {
  df.WM.glucose$CUE.fit <- getPred(as.list(coef(fit)),df.WM.glucose$WtMeanCopyNum)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.WM.glucose)
summary(df.lm)
# summary(df.lm)$r.squared
WtMeanCopyNum.glucose.R2adj = summary(df.lm)$adj.r.squared
WtMeanCopyNum.glucose.p.value = summary(df.lm)$coefficients[8]

## PLOT:
p_CUEgl_copyNum_O = ggplot(subset(df.WM.glucose, 
                                  substrate == 'glucose' & horizon == 'O'), 
       aes(x=(WtMeanCopyNum), 
           y=CUE,
           color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  geom_line(aes(x=WtMeanCopyNum,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = 'Predicted weighted mean 16S rRNA gene copy number', 
       y = 'CUE',
       color='Burn treatment\nexposure',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 
p_CUEgl_copyNum_O

paste(sprintf('Glucose: y = a * exp(-b * WtMeanCopyNum) ; where a = '), 
      sprintf(glucose_a, fmt='%#.2f'), ", b = ", 
      sprintf(glucose_b, fmt='%#.2f'), ", R2adj. = ", 
      sprintf(WtMeanCopyNum.glucose.R2adj, fmt='%#.2f'), ' and p = ',
      sprintf(WtMeanCopyNum.glucose.p.value, fmt='%#.1e'))



### PLOT: Figure 6. gene copy number, M horizon vs. CUE_pine -----------------------------
df.WM.glucose = subset(df.WM.merge, horizon == 'M' & substrate=='glucose')

## TRYING WITH NLS.LM FUNCTION
getPred = function(params, xx) params$a*exp(-params$b * xx)
# Predicted CUE = a * exp(-b*max.temp.C)

residFun <- function(params, CUE, WtMeanCopyNum) CUE - getPred(params,WtMeanCopyNum) 
# function(params, observed, xx) observed - getPred(p,xx)
# Observation = CUE, getPred = predicted CUE; nls.lm and ResidFun is trying to minimilize the sum of squares between the two vectors. 

# starting params, glucose: a= , b= , c= 
fit = nls.lm(par=list(a=1, b=1), 
             fn = residFun, 
             WtMeanCopyNum = df.WM.glucose$WtMeanCopyNum,
             CUE = df.WM.glucose$CUE,
             lower = c(0,0), upper = c(Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
glucose_a = summary(fit)$coefficients[[1]]
glucose_b = summary(fit)$coefficients[[2]]


for (i in 1:nrow(df.WM.glucose)) {
  df.WM.glucose$CUE.fit <- getPred(as.list(coef(fit)),df.WM.glucose$WtMeanCopyNum)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.WM.glucose)
summary(df.lm)
# summary(df.lm)$r.squared
WtMeanCopyNum.glucose.R2adj = summary(df.lm)$adj.r.squared
WtMeanCopyNum.glucose.p.value = summary(df.lm)$coefficients[8]

## PLOT:
p_CUEgl_copyNum_M = ggplot(subset(df.WM.glucose, 
                                  substrate == 'glucose' & horizon == 'M'), 
                           aes(x=(WtMeanCopyNum), 
                               y=CUE,
                               color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  # geom_line(aes(x=WtMeanCopyNum,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = 'Predicted weighted mean 16S rRNA gene copy number', 
       y = 'CUE',
       color='Burn treatment\nexposure',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Glucose amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 
p_CUEgl_copyNum_M

paste(sprintf('Glucose, M horizon: y = a * exp(-b * WtMeanCopyNum) ; where a = '), 
      sprintf(glucose_a, fmt='%#.2f'), ", b = ", 
      sprintf(glucose_b, fmt='%#.2f'), ", R2adj. = ", 
      sprintf(WtMeanCopyNum.glucose.R2adj, fmt='%#.2f'), ' and p = ',
      sprintf(WtMeanCopyNum.glucose.p.value, fmt='%#.1e'))


### PLOT: Figure S2. gene copy number vs. CUE_pine ----------------------------------------

## Fitting non-linear regression: exponential decay, PINE:
df.WM.pine = subset(df.WM.merge, horizon == 'O' & substrate=='pine')

anova <- aov(CUE~WtMeanCopyNum*site.type*burn.trtmt.duration.seconds, data = subset(df.WM.pine))
summary(anova) # No significant difference in CUE across site types

fit = nls.lm(par=list(a=1, b=1), 
             fn = residFun, 
             WtMeanCopyNum = df.WM.pine$WtMeanCopyNum,
             CUE = df.WM.pine$CUE,
             lower = c(0,0), upper = c(Inf, Inf),
             control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))

summary(fit)
summary(fit)$coefficients
pine_a = summary(fit)$coefficients[[1]]
pine_b = summary(fit)$coefficients[[2]]


for (i in 1:nrow(df.WM.pine)) {
  df.WM.pine$CUE.fit <- getPred(as.list(coef(fit)),df.WM.pine$WtMeanCopyNum)
}

df.lm <- lm(CUE ~ CUE.fit, data = df.WM.pine)
summary(df.lm)
# summary(df.lm)$r.squared
WtMeanCopyNum.pine.R2adj = summary(df.lm)$adj.r.squared
WtMeanCopyNum.pine.p.value = summary(df.lm)$coefficients[8]

## PLOT:
p_CUEpi_copyNum_O = ggplot(subset(df.WM.pine, 
                                  substrate == 'pine' & horizon == 'O'), 
       aes(x=(WtMeanCopyNum), 
           y=CUE,
           color=burn.trtmt.duration.seconds)) +
  geom_point(size=3, alpha=0.7)+
  # geom_smooth(formula = y~x, method='lm', se=FALSE, color = 'grey40') +
  # geom_line(aes(x=WtMeanCopyNum,y=CUE.fit), color = 'grey20') +
  scale_color_manual(values = Burn_colors,
                     labels=Burn_labels)+
  theme_bw() + 
  facet_grid(~horizon, scales='free',
             labeller=labeller(horizon=Horizon_labels,
                               site.type = SiteType_labels)) +
  labs(x = 'Predicted weighted mean 16S rRNA gene copy number', 
       y = 'CUE',
       color='Burn treatment\nexposure',
       shape='Soil horizon',
       # fill = 'Burn treatment exposure\n(seconds)', 
       title = 'Pine amendment') + 
  theme(axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)) 

p_CUEpi_copyNum_O

