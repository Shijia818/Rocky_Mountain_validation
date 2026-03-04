
library(geometry)
library(sp)
library(rgl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)
library(randomForest)
library(circular)
library(bpnreg)
library(lubridate)
library(lme4)
library(lmerTest)
library(purrr)
library(performance)
library(relaimpo)
library(emmeans)
library(tibble)
library(Hmisc)
library(fitdistrplus)
library(quantreg)
library(spaMM)


setwd("E:/Rocky_Mountain")

data.phenology.field <- get(load("./Data/Phenology_data.RData"))
data.phenology.field$species <- gsub(" ","_",data.phenology.field$species)
data.phenology.field <- data.phenology.field %>% mutate(DOY = as.numeric(format(as.Date(date),"%j"))) %>% filter(floralcount > 0)
df.field <- data.phenology.field

data.phenology.herbarium <- read.csv("./Data/data.herbarium.refine.csv")
data.phenology.herbarium <- data.phenology.herbarium %>% mutate(date = make_date(Year, Month, Day)) %>% mutate(DOY = as.numeric(format(as.Date(date),"%j"))) %>% filter(flower > 0)

df.herbarium <- data.phenology.herbarium
df.field <- df.field[which(df.field$species %in% unique(df.herbarium$species)), ]

results.field <- df.field %>% group_by(species) %>% summarise(fit = list(fitdist(doy, "weibull"))) %>% rowwise() %>% mutate(
  shape = fit$estimate["shape"],
  scale = fit$estimate["scale"],
  Q10 = qweibull(0.10, shape, scale),
  Q50 = qweibull(0.50, shape, scale),
  Q90 = qweibull(0.90, shape, scale),
) %>% dplyr::select(-fit)


results.herbarium <- df.herbarium %>% group_by(species) %>% summarise(fit = list(fitdist(julian, "weibull"))) %>% rowwise() %>% mutate(
  shape = fit$estimate["shape"],
  scale = fit$estimate["scale"],
  Q10 = qweibull(0.10, shape, scale),
  Q50 = qweibull(0.50, shape, scale),
  Q90 = qweibull(0.90, shape, scale),
) %>% dplyr::select(-fit)


################################################
################################################

df_T_PRISM <- get(load("./Analyses/Temperature/df_T_PRISM.RData"))
df_T_Field <- get(load("./Analyses/Temperature/df_T_Field.RData"))
df_T_Herbarium <- get(load("./Analyses/Temperature/df_T_Herbarium.RData"))

result.flr.field <- df_T_Field %>% group_by(species, year, plot) %>% summarise(max_flower = max(floralcount, na.rm = TRUE)) %>% ungroup()
result.flr.field$species_year <- paste(result.flr.field$species, result.flr.field$year,sep= "_")
df_T_Field$species_year <- paste(df_T_Field$species,df_T_Field$year,sep= "_")

result.flr.field.match <- result.flr.field[match(df_T_Field$species_year, result.flr.field$species_year), ]
df_T_Field$maxfloral <- result.flr.field.match$max_flower

result.flr.PRISM <- df_T_PRISM %>% group_by(species, year, plot) %>% summarise(max_flower = max(floralcount, na.rm = TRUE)) %>% ungroup()
result.flr.PRISM$species_year <- paste(result.flr.PRISM$species, result.flr.PRISM$year,sep= "_")
df_T_PRISM$species_year <- paste(df_T_PRISM$species,df_T_PRISM$year,sep= "_")

result.flr.PRISM.match <- result.flr.PRISM[match(df_T_PRISM$species_year, result.flr.PRISM$species_year), ]
df_T_PRISM$maxfloral <- result.flr.PRISM.match$max_flower

df_T_Herbarium$species <- gsub("_", " ", df_T_Herbarium$species)
df_T_PRISM <- df_T_PRISM[which(df_T_PRISM$species %in% df_T_Herbarium$species),]
df_T_Field <- df_T_Field[which(df_T_Field$species %in% df_T_Herbarium$species),]

df_T_Herbarium <- df_T_Herbarium[, c("species", "julian", "flower", "T_spring","Year")]
df_T_PRISM <- df_T_PRISM[, c("plot", "species", "doy", "maxfloral", "T_spring","year")]
df_T_Field <- df_T_Field[, c("plot", "species", "doy", "maxfloral", "T_spring", "year")]
colnames(df_T_Herbarium) <- c("species", "doy", "maxfloral", "T_spring", "year")
df_T_Herbarium <- cbind(plot = "NA", df_T_Herbarium)
df_T_Herbarium$Type <- "Herbarium"
df_T_Field$Type <- "Field"
df_T_PRISM$Type <- "PRISM"

df_T_total <- rbind(df_T_Herbarium, df_T_Field, df_T_PRISM)

res.mean.T <- lmerTest::lmer(doy ~ scale(T_spring)* Type + scale(maxfloral) + (1 | species) + (1 | plot) + (0 + scale(T_spring):Type | species), data = df_T_total)


coefs <- coef(res.mean.T)$species
slopes <- coefs %>% as.data.frame() %>% mutate(species = rownames(.), slopes_PRISM = coefs$`scale(T_spring)` + coefs$`scale(T_spring):TypePRISM`, slopes_Field = coefs$`scale(T_spring)` + coefs$`scale(T_spring):TypeField`, slopes_Herbarium = coefs$`scale(T_spring)` + coefs$`scale(T_spring):TypeHerbarium`)%>% dplyr::select(species, slopes_PRISM, slopes_Field, slopes_Herbarium) %>% pivot_longer(col = starts_with("slopes"), names_to = "Type",values_to= "Slope")
slopes$Type <- factor(slopes$Type, levels = c("slopes_Field","slopes_PRISM","slopes_Herbarium"))

slopes_wide <- slopes %>% pivot_wider(names_from = Type, values_from = Slope)
t.test(slopes_wide$slopes_Field,slopes_wide$slopes_Herbarium, paired = T)
t.test(slopes_wide$slopes_PRISM,slopes_wide$slopes_Herbarium, paired = T)
t.test(slopes_wide$slopes_PRISM,slopes_wide$slopes_Field, paired = T)

###########################################
#############################################################

species_list <- unique(df_T_PRISM$species)
results.field.rq <- as.data.frame(matrix(NA, 45, 7))
colnames(results.field.rq) <- c("species","Q10", "Q10_SE","Q50", "Q50_SE", "Q90", "Q90_SE")
results.field.rq$species <- species_list

for(i in 1:length(species_list)){tryCatch({
  
  sp <- species_list[i]
  sp_data <- df_T_PRISM %>% filter(species == sp)
  
  Q10 = rq(scale(doy) ~ T_spring + maxfloral, tau = 0.10, data = sp_data)
  Q50 = rq(scale(doy) ~ T_spring + maxfloral, tau = 0.50, data = sp_data)
  Q90 = rq(scale(doy) ~ T_spring + maxfloral, tau = 0.90, data = sp_data)
  
  results.field.rq[i,2] <- summary(Q10)$coe[2,][1]
  results.field.rq[i,3] <- summary(Q10)$coe[2,][2]
  results.field.rq[i,4] <- summary(Q50)$coe[2,][1]
  results.field.rq[i,5] <- summary(Q50)$coe[2,][2]
  results.field.rq[i,6] <- summary(Q90)$coe[2,][1]
  results.field.rq[i,7] <- summary(Q90)$coe[2,][2]
}, error = function(e){
  message("species", species_list[i], " error，skip")
  results.field.rq[i,2:7] <- NA  
})
  
}

species_list_H <- unique(df_T_Herbarium$species)
results.herbarium.rq <- as.data.frame(matrix(NA, 45, 7))
colnames(results.herbarium.rq) <- c("species","Q10", "Q10_SE", "Q50", "Q50_SE", "Q90", "Q90_SE")
results.herbarium.rq$species <- species_list_H

for(i in 1:length(species_list_H)){
  
  sp <- species_list_H[i]
  sp_data <- df_T_Herbarium %>% filter(species == sp)
  
  Q10 = rq(scale(doy) ~ T_spring + maxfloral, tau = 0.10, data = sp_data)
  Q50 = rq(scale(doy) ~ T_spring + maxfloral, tau = 0.50, data = sp_data)
  Q90 = rq(scale(doy) ~ T_spring + maxfloral, tau = 0.90, data = sp_data)
  
  results.herbarium.rq[i,2] <- summary(Q10)$coe[2,][1]
  results.herbarium.rq[i,3] <- (summary(Q10)$coe[2,][3] - summary(Q10)$coe[2,][2])/(2*1.96)
  results.herbarium.rq[i,4] <- summary(Q50)$coe[2,][1]
  results.herbarium.rq[i,5] <- (summary(Q50)$coe[2,][3] - summary(Q50)$coe[2,][2])/(2*1.96)
  results.herbarium.rq[i,6] <- summary(Q90)$coe[2,][1]
  results.herbarium.rq[i,7] <- (summary(Q90)$coe[2,][3] - summary(Q90)$coe[2,][2])/(2*1.96)
}

results.field.rq.match <- results.field.rq[match(results.herbarium.rq$species, results.field.rq$species), ]
results.field.rq.match <- results.field.rq.match[which(results.field.rq.match$species !="Viola adunca"), ]
results.herbarium.rq <- results.herbarium.rq[which(results.herbarium.rq$species !="Viola adunca"), ]


######################
######################  snow density ########################

df_May_PRISM <- get(load("./Snow_density/df_May_PRISM.RData"))
df_May_Field <- get(load("./Snow_density/df_May_Field.RData"))
df_May_Herbarium <- get(load("./Snow_density/df_May_Herbarium.RData"))

df_May_Herbarium <- df_May_Herbarium %>% filter(floralcount != 0)

df_May_PRISM <- df_May_PRISM[which(df_May_PRISM$year != "1973" & df_May_PRISM$year != "1974"), ]
df_May_Herbarium <- df_May_Herbarium[which(df_May_Herbarium$year != "1973" & df_May_Herbarium$year != "1974"), ]

result.flr.field <- df_May_Field %>% group_by(species, year, plot) %>% summarise(max_flower = max(floralcount, na.rm = TRUE)) %>% ungroup()
result.flr.field$species_year <- paste(result.flr.field$species, result.flr.field$year,sep= "_")
df_May_Field$species_year <- paste(df_May_Field$species,df_May_Field$year,sep= "_")

result.flr.field.match <- result.flr.field[match(df_May_Field$species_year, result.flr.field$species_year), ]
df_May_Field$maxfloral <- result.flr.field.match$max_flower

result.flr.PRISM <- df_May_PRISM %>% group_by(species, year, plot) %>% summarise(max_flower = max(floralcount, na.rm = TRUE)) %>% ungroup()
result.flr.PRISM$species_year <- paste(result.flr.PRISM$species, result.flr.PRISM$year,sep= "_")
df_May_PRISM$species_year <- paste(df_May_PRISM$species,df_May_PRISM$year,sep= "_")

result.flr.PRISM.match <- result.flr.PRISM[match(df_May_PRISM$species_year, result.flr.PRISM$species_year), ]
df_May_PRISM$maxfloral <- result.flr.PRISM.match$max_flower

df_May_Herbarium$species_year <- paste(df_May_Herbarium$species, df_May_Herbarium$year, sep = "_")
df_May_Herbarium$maxfloral <- df_May_Herbarium$floralcount

df_May_PRISM <- df_May_PRISM[which(df_May_PRISM$species %in% df_May_Herbarium$species),]
df_May_Field <- df_May_Field[which(df_May_Field$species %in% df_May_Herbarium$species),]

df_May_PRISM$Type <- "PRISM"
df_May_Field$Type <- "Field"
df_May_Herbarium$Type <- "Herbarium"


df_snow_total <- rbind(df_May_PRISM, df_May_Field, df_May_Herbarium)

res.mean.snow <- lmerTest::lmer(doy ~  scale(snow.density.bb)*Type + scale(maxfloral) + (1 | species) + (1 | plot) + (0 + scale(snow.density.bb):Type | species), data = df_snow_total)

coefs <- coef(res.mean.snow)$species
slopes <- coefs %>% as.data.frame() %>% mutate(species = rownames(.), slopes_PRISM = coefs$`scale(snow.density.bb)` + coefs$`scale(snow.density.bb):TypePRISM`, slopes_Field = coefs$`scale(snow.density.bb)` + coefs$`scale(snow.density.bb):TypeField`, slopes_Herbarium = coefs$`scale(snow.density.bb)` + coefs$`scale(snow.density.bb):TypeHerbarium`)%>% dplyr::select(species, slopes_PRISM, slopes_Field, slopes_Herbarium) %>% pivot_longer(col = starts_with("slopes"), names_to = "Type",values_to= "Slope")
slopes$Type <- factor(slopes$Type, levels = c("slopes_Field","slopes_PRISM","slopes_Herbarium"))
slopes_wide <- slopes %>% pivot_wider(names_from = Type, values_from = Slope)

t.test(slopes_wide$slopes_Herbarium,slopes_wide$slopes_PRISM, paired = T)
t.test(slopes_wide$slopes_Herbarium,slopes_wide$slopes_Field, paired = T)
t.test(slopes_wide$slopes_PRISM,slopes_wide$slopes_Field, paired = T)

###############################
#######################################################

species_list <- unique(df_May_PRISM$species)
results.field.rq <- as.data.frame(matrix(NA, 45, 7))
colnames(results.field.rq) <- c("species","Q10", "Q10_SE", "Q50", "Q50_SE", "Q90", "Q90_SE")
results.field.rq$species <- species_list

for(i in 1:length(species_list)){tryCatch({
  
  sp <- species_list[i]
  sp_data <- df_May_PRISM %>% filter(species == sp)
  
  Q10 = rq(scale(doy) ~ snow.density.bb + maxfloral, tau = 0.10, data = sp_data)
  Q50 = rq(scale(doy) ~ snow.density.bb + maxfloral, tau = 0.50, data = sp_data)
  Q90 = rq(scale(doy) ~ snow.density.bb + maxfloral, tau = 0.90, data = sp_data)
  
  results.field.rq[i,2] <- summary(Q10)$coe[2,][1]
  results.field.rq[i,3] <- summary(Q10)$coe[2,][2]
  results.field.rq[i,4] <- summary(Q50)$coe[2,][1]
  results.field.rq[i,5] <- summary(Q50)$coe[2,][2]
  results.field.rq[i,6] <- summary(Q90)$coe[2,][1]
  results.field.rq[i,7] <- summary(Q90)$coe[2,][2]
}, error = function(e){
  message("species ", species_list[i], " error，skip")
  results.field.rq[i,2:7] <- NA  
})
  
}

species_list_H <- unique(df_May_Herbarium$species)
results.herbarium.rq <- as.data.frame(matrix(NA, 45, 7))
colnames(results.herbarium.rq) <- c("species","Q10", "Q10_SE","Q50", "Q50_SE", "Q90", "Q90_SE")
results.herbarium.rq$species <- species_list_H

for(i in 1:length(species_list_H)){
  
  sp <- species_list_H[i]
  sp_data <- df_May_Herbarium %>% filter(species == sp)
  
  Q10 = rq(scale(doy) ~ snow.density.bb + maxfloral, tau = 0.10, data = sp_data)
  Q50 = rq(scale(doy) ~ snow.density.bb + maxfloral, tau = 0.50, data = sp_data)
  Q90 = rq(scale(doy) ~ snow.density.bb + maxfloral, tau = 0.90, data = sp_data)
  
  results.herbarium.rq[i,2] <- summary(Q10)$coe[2,][1]
  results.herbarium.rq[i,3] <- (summary(Q10)$coe[2,][3] - summary(Q10)$coe[2,][2])/(2*1.96)
  results.herbarium.rq[i,4] <- summary(Q50)$coe[2,][1]
  results.herbarium.rq[i,5] <- (summary(Q50)$coe[2,][3] - summary(Q50)$coe[2,][2])/(2*1.96)
  results.herbarium.rq[i,6] <- summary(Q90)$coe[2,][1]
  results.herbarium.rq[i,7] <- (summary(Q90)$coe[2,][3] - summary(Q90)$coe[2,][2])/(2*1.96)
}

results.field.rq.match <- results.field.rq[match(results.herbarium.rq$species, results.field.rq$species), ]
results.field.rq.match <- results.field.rq.match[which(results.field.rq.match$species !="Viola adunca"), ]
results.herbarium.rq <- results.herbarium.rq[which(results.herbarium.rq$species !="Viola adunca"), ]

###############################
############################# sensitivity analyses ##################


df_T_PRISM <- get(load("./Analyses/Temperature/df_T_PRISM.RData"))
df_T_Field <- get(load("./Analyses/Temperature/df_T_Field.RData"))
df_T_Herbarium <- get(load("./Analyses/Temperature/df_T_Herbarium.RData"))

df_T_mean_Field_match <- df_T_Field %>% group_by(species, year, plot) %>% summarise(max_flower = max(floralcount, na.rm = TRUE), mean_deg = mean(circular(doy, units = "degrees", template = "clock24", modulo = "2pi")) %>% as.numeric(), mean_T = mean(T_spring, na.rm = TRUE), .groups = "drop") %>% mutate(mean_doy = round(mean_deg * 365 / 360)) %>% ungroup()


df_T_mean_PRISM_match <- df_T_PRISM %>% group_by(species, year, plot) %>% summarise(max_flower = max(floralcount, na.rm = TRUE), mean_deg = mean(circular(doy, units = "degrees", template = "clock24", modulo = "2pi")) %>% as.numeric(), mean_T = mean(T_spring, na.rm = TRUE), .groups = "drop") %>% mutate(mean_doy = round(mean_deg * 365 / 360)) %>% ungroup()



df_T_Herbarium$species <- gsub("_", " ", df_T_Herbarium$species)
df_T_mean_PRISM_match <- df_T_mean_PRISM_match[which(df_T_mean_PRISM_match$species %in% df_T_Herbarium$species),]
df_T_mean_Field_match <- df_T_mean_Field_match[which(df_T_mean_Field_match$species %in% df_T_Herbarium$species),]

df_T_Herbarium <- df_T_Herbarium[, c("species", "julian", "flower", "T_spring","Year")]
df_T_mean_PRISM_match <- df_T_mean_PRISM_match[, c("plot", "species", "mean_doy", "max_flower", "mean_T","year")]
df_T_mean_Field_match <- df_T_mean_Field_match[, c("plot", "species", "mean_doy", "max_flower", "mean_T","year")]
colnames(df_T_Herbarium) <- c("species", "mean_doy", "max_flower", "mean_T", "year")
df_T_Herbarium <- cbind(plot = "NA", df_T_Herbarium)
df_T_Herbarium$Type <- "Herbarium"
df_T_mean_Field_match$Type <- "Field"
df_T_mean_PRISM_match$Type <- "PRISM"

df_T_total <- rbind(df_T_mean_PRISM_match, df_T_mean_Field_match, df_T_Herbarium)

res.T.model <- lmerTest::lmer(mean_doy ~ scale(mean_T) * Type + (1 | species) + (1 | plot) + (0 + scale(mean_T):Type | species), data = df_T_total)

model.T.mean <- fitme(mean_doy ~ scale(mean_T) * Type + (1 | species) + (0 + scale(mean_T):Type | species) + AR1(1 | year), data = df_T_total, family = gaussian(), method = "REML")












