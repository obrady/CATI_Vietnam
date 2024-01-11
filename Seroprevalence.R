

library(mapview)
library(tidyverse)
library(sf)
library(raster)
library(rgdal)
library(mgcv)
library(broom)
library(SpatialDengue)
library(Information)
library(ClustOfVar)
library(reshape2)
library(caret)
library(pROC)
library(ggpubr)
library(binom)

path = "C:/Users/User1/Desktop/Dengue_Vietnam/"

nha.trang.survey.ras <-raster(paste0(path, "Data/Maps/nha.trang.survey.ras.tif"))
nha.trang.study.ras <-raster(paste0(path, "Data/Maps/nha.trang.study.ras.tif"))
nha.trang.survey.shp <- readOGR(paste0(path, "Data/Maps/nha.trang.survey.shp.shp"))
nha.trang.study.shp <- readOGR(paste0(path, "Data/Maps/nha.trang.study.shp.shp" ))

#import sero-survey
surv1 <- read.csv(paste0(path, 'Data/se.csv'), fileEncoding="latin1")
surv<-surv1 %>%
  mutate(agegroup2 = ifelse(commune == 'DB', as.numeric(agegroup) + 5, as.numeric(agegroup))) %>%
  mutate(agegroup2 = as.factor(agegroup2)) %>%
  mutate(commune2 = ifelse(commune == "DA", 1, 2)) %>%
  mutate(educ = ifelse(highschool == 1, 2, ifelse(postsecond == 1, 3, ifelse(noschool ==1, 1, NA)))) %>%
  mutate(educ = as.factor(educ)) %>%
  #any education
  mutate(educ2 = ifelse(noschool== 0, 1, ifelse(highschool==1 | postsecond ==1, 0, NA))) %>%
  mutate(occup = ifelse(OfficeProfessional ==1, 1, ifelse(Worker==1,2 , ifelse(Other==1, 3, NA)))) %>%
  mutate (occup = as.factor(occup)) %>%
  mutate(occup2 = ifelse(OfficeProfessional==1, 1, ifelse(Worker==1 | Other==1, 0, NA))) %>%
  mutate(head = ifelse(head_father==1, 1,
                       ifelse(head_mother==1, 2, ifelse(head_gparent ==1, 3,
                                                        ifelse(head_other2==1, 4, NA))))) %>%
  mutate(head2 = ifelse(head==1 | head==2, 1, ifelse(head>2, 0, NA))) %>%
  mutate(head=as.factor(head)) %>%
  mutate(hhsize = ifelse(hlb05 < 0 , NA, hlb05)) %>%
  mutate(new_IgGPos = ifelse(IgG_panbio_units >=2.2, 1, 0), new_IgMPos = ifelse(IgM_panbio_units>= 8.1, 1, 0)) %>%
  #mutate(sero = ifelse(new_IgGPos == 1 | new_IgMPos == 1, 1,0 )) %>%
  mutate(sero = ifelse(new_IgGPos == 1, 1,0 )) %>%
  mutate(smoke.in = ifelse(hlb09_1d ==1 |hlb09_2d == 1 | hlb09_3d == 1 | hlb09_4d == 1, 1, 0))


gis_full2 <- read.csv(paste0(path, "/Data/gis_full2.csv"), fileEncoding = "latin1" )
gis_inerested <- gis_full2 %>%
  dplyr::select(hlperno, final_X, final_Y, lat, long ) %>%
  mutate(long1=lat, lat1=long) %>%
  dplyr::select(-lat, -long) %>%
  mutate(long=long1, lat=lat1)

surv <- left_join(surv, gis_inerested, by=c("hlperno"="hlperno"))

#survey covariates
covariates <- stack(paste0(path, "Data/Maps/survey_covariates.tif"))

slope <- raster(paste0(path, "Data/Maps/vnm_srtm_slope_100m.tif"))
slope_2 <- slope
slope_resampled <- resample(slope, nha.trang.survey.ras, method = "bilinear" ) #bilinear because longitudinal data

slope.survey<- crop(slope_resampled, extent(nha.trang.survey.shp))
slope.survey <- mask(slope.survey, nha.trang.survey.shp)
#plot(slope.survey)


covariates <- stack(covariates, slope.survey)

#study covariates

names(covariates) <-  c("pop", "built", "elevation", "lights", "waterway", "nature", "slope")



#convert pop raster to a dataframe with patchID
nt_unipix <- make.unipix(nha.trang.survey.ras)
#head(covariates_unipix)

#match patchID from pop dataframe to survey data
surv$patchID <- apply(cbind(surv$long, surv$lat), 1, pix.id.find, nt_unipix)
#head(survey)

covariates_df<- as.data.frame(covariates)

covariates_df_2 <- covariates_df %>%
  filter(!is.na(pop))


covariates_unipix = data.frame(pixID = (1:length(as.vector(covariates$pop)))[!is.na(as.vector(covariates$pop))], 
                               pop = as.vector(covariates$pop)[!is.na(as.vector(covariates$pop))])
covariates_unipix = covariates_unipix[covariates_unipix$pop > 0, ]
covariates_unipix = data.frame(patchID = 1:nrow(covariates_unipix), covariates_unipix, coordinates(covariates)[covariates_unipix$pixID, 
                                                                                                               ])

covariates_unipix_2 = cbind(covariates_unipix, built = covariates_df_2$built, 
                            elevation = covariates_df_2$elevation, lights = covariates_df_2$lights,
                            waterway= covariates_df_2$waterway, nature = covariates_df_2$nature,
                            slope = covariates_df_2$slope)

surv<- left_join(surv, covariates_unipix_2, by =c("patchID" = "patchID")) 

surv <- surv %>%
  mutate(built = ifelse(built<0.5, 0, ifelse(built >=-.5, 1, built)))



surv_sf <- st_as_sf(surv, coords = c("long", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#now need to match the cases to survey area


surv_int <- st_intersection(surv_sf, st_as_sf(nha.trang.survey.shp))


survey <- surv_int
survey_sf <- st_as_sf(survey, coords = c("long1", "lat1"),crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" ) %>%
  mutate(sero1 = ifelse(sero ==1, "Positive", "Negative"))

#survey$patchID <- apply(cbind(survey$long, survey$lat), 1, pix.id.find, covariates_unipix)


#sero_data_cov <- left_join(survey, covariates_unipix_2, by =c("patchID" = "patchID")) %>%
  #mutate(built = ifelse(built<0.5, 0, ifelse(built >=-.5, 1, built)))

as.data.frame(survey_sf) %>%
  dplyr::select(pop, elevation, lights,
          waterway, nature, slope, sero, long1, lat1) %>%
  cor(method = "spearman")





# survey_sf_int <- st_intersection(survey_sf, st_as_sf(nha.trang.survey.shp))

mapview(nha.trang.survey.ras)
mapview(nha.trang.survey.ras, na.color = NA, layer.name= "nha.trang.survey.ras") + 
  mapview(survey_sf, zcol = "sero1", cex = 3 , alpha = 0.5, layer.name = "Case")

survey_sf$weight_2 = survey_sf$weight/2



nt_unipix2 <- nt_unipix %>%
  rename(long1 = x, lat1=y)


nt_unipix2 %>% 
  bind_cols(predict.gam(object = gam_coor,
                        newdata = ., 
                        type = "response",
                        se.fit = F) %>%
              tibble(value = .)) %>%
  na.omit %>%
  ggplot(data = ., aes(x=long1, y=lat1)) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  theme_bw() +
  ylab("Latitude")  + xlab("Longitude") +
  scale_fill_viridis_c(limits = c(0.8,1),
                       name = "Previous exposure to dengue") +
  theme(legend.position = "bottom") +
  ggtitle("Predicted previous exposure to dengue")


nt_unipix2 %>% 
  bind_cols(predict.gam(object = gam_coor1,
                        newdata = ., 
                        type = "response",
                        se.fit = F) %>%
              tibble(value = .)) %>%
  na.omit %>%
  ggplot(data = ., aes(x=long1, y=lat1)) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  theme_bw() +
  ylab("Latitude")  + xlab("Longitude") +
  scale_fill_viridis_c(
                       name = "Previous exposure to dengue") +
  theme(legend.position = "bottom") +
  ggtitle("Predicted previous exposure to dengue")



nt_unipix2 %>% 
  bind_cols(predict.gam(object = gam_coor1,
                        newdata = ., 
                        type = "response",
                        se.fit = F) %>%
              tibble(value = .)) %>%
  summary(.)






nt_unipix2 %>% 
  bind_cols(predict.gam(object = gam_coor,
                        newdata = ., 
                        type = "response",
                        se.fit = F) %>%
              tibble(value = .)) %>%
  summary(.)





nt_unipix_full <- make.unipix(nha.trang.study.ras) %>%
rename(long1 = x, lat1=y)










formula.all <- as.formula(factor(sero) ~ s(pop) + s(elevation) + s(nature) + s(waterway) +  s(lights) + s(slope))
formula.all_sp <- as.formula(factor(sero) ~ s(pop) +  s(elevation) + s(nature) + s(waterway) +  s(lights) + s(slope) + s(long1, lat1, bs = "gp"))


formula.all_pop <- as.formula(factor(sero) ~ s(elevation) + s(nature) + s(waterway) +  s(lights) + s(slope))
formula.all_pop_sp <- as.formula(factor(sero) ~ s(elevation) + s(nature) + s(waterway) +  s(lights) + s(slope) + s(long1, lat1, bs = "gp"))


formula.sp<-as.formula(factor(sero) ~ s(long1, lat1, bs = "gp"))
formula.all.lm <- as.formula(factor(sero) ~ pop + elevation + nature + waterway +  lights + slope)
formula.all.lm.sp <- as.formula(factor(sero) ~ pop + elevation + nature + waterway +  lights + slope + s(long1, lat1, bs = "gp"))

gam_all_train <- gam(formula = formula.all, data = seroTrain, weights =weight_2, family = binomial, method = "REML", select = TRUE)
a_all_train <- gam_all_train$aic
  
  
gam_all_sp_train <- gam(formula = formula.all_sp, data = seroTrain, weights = weight_2, family = binomial, method= "REML", select = TRUE)
a_all_sp_train <- gam_all_sp_train$aic


gam_all_pop_train <- gam(formula = formula.all_pop, data = seroTrain, weights = weight_2, family = binomial, method = "REML", select = TRUE)
a_all_pop_train <- gam_all_pop_train$aic


gam_all_pop_sp_train <- gam(formula = formula.all_pop_sp, data = seroTrain, weights = weight_2,family = binomial, method= "REML", select = TRUE)
a_all_pop_sp_train <- gam_all_pop_sp_train$aic


gam_sp_train <- gam(formula = formula.sp, data = seroTrain, weights = weight2, family = binomial, method = "REML", select = TRUE)
a_sp_train <- gam_sp_train$aic



gam_lm_train <- gam(formula = formula.all.lm, data = seroTrain, weights = weight_2, family = binomial, method= "REML", select = TRUE)
a_lm_train <- gam_lm_train$aic


gam_lm_sp_train <- gam(formula = formula.all.lm.sp, data = seroTrain, weights = weight_2, family = binomial, method = "REML", select = TRUE)
a_lm_sp_train <- gam_lm_sp_train$aic




gam_null <- gam(factor(sero) ~ 1, weights = weight_2, data = seroTrain,  family = binomial, method = "REML")

a_null <- gam_null$aic





aic_values <- data.frame(a_all_train, a_all_sp_train, a_all_pop_train ,a_all_pop_sp_train, a_sp_train, a_lm_train, a_lm_sp_train, a_null)

min_aic <- aic_values[which.min(aic_values)]


#gam_null had lowest aic so select this model
  
  
  
  
  #STUDY DATA
  study.covariates <- stack(paste0(path, "Data/Maps/study_covariates.tif"))
  
  names(study.covariates) <- c("pop", "built", "elevation", "lights", "waterway", "nature")
  
  slope <- raster(paste0(path, "Data/Maps/vnm_srtm_slope_100m.tif"))
  slope_2 <- slope
  slope_resampled <- resample(slope, nha.trang.study.ras, method = "bilinear" ) #bilinear because longitudinal data
  
  slope.study<- crop(slope_resampled, extent(nha.trang.study.shp))
  slope.study <- mask(slope.study, nha.trang.study.shp)
  #plot(slope.study)
  
  
  study.covariates <- stack(study.covariates, slope.study)
  names(study.covariates) <- c("pop", "built", "elevation", "lights", "waterway", "nature", "slope")
  
  
  #convert study.covariates to data, to  use for prediction
  
  study.covariates_df <- as.data.frame(study.covariates)
  
  study.covariates_df_2 <-  study.covariates_df %>%
    filter(!is.na(pop))
  
  
  study.covariates_unipix = data.frame(pixID = (1:length(as.vector(study.covariates$pop)))[!is.na(as.vector(study.covariates$pop))], 
                                       pop = as.vector(study.covariates$pop)[!is.na(as.vector(study.covariates$pop))])
  study.covariates_unipix = study.covariates_unipix[study.covariates_unipix$pop > 0, ]
  study.covariates_unipix = data.frame(patchID = 1:nrow(study.covariates_unipix), study.covariates_unipix, coordinates(study.covariates)[study.covariates_unipix$pixID, 
                                                                                                                                         ])
  
  
  
  study.covariates_unipix_2 = cbind(study.covariates_unipix, built = study.covariates_df_2$built, 
                                    elevation = study.covariates_df_2$elevation,lights = study.covariates_df_2$lights,
                                    waterway= study.covariates_df_2$waterway, nature = study.covariates_df_2$nature, slope = study.covariates_df_2$slope) %>%
    mutate(built = ifelse(built<0.5, 0, ifelse(built >=-.5, 1, built)))
  
  
  #, 
  
  
  study.covariates_unipix_2 <- study.covariates_unipix_2 %>%
    rename(long1 = x, lat1=y)
  
  
  
# study.covariates_unipix_2$lat <-study.covariates_unipix_2$x
# study.covariates_unipix_2$long <- study.covariates_unipix_2$y




predict_study_basic <- predict.gam(gam_none, newdata = study.covariates_unipix_2, type = "response")


new_study_data <- cbind( study.covariates_unipix_2, predict_study_basic)



ggplot(new_study_data) +
  geom_tile(aes(x=x, y=y, fill = predict_study_basic))
quantile(as.numeric(as.vector(predict_study_all)), c(0.025, 0.5, 0.95))







make.unipix.modified<- function (popras){
unipix = data.frame(pixID = (1:length(as.vector(popras))))
unipix = data.frame(patchID = 1:nrow(unipix), unipix, coordinates(popras)[unipix$pixID,])
}



nt.study.unipix <- make.unipix(nha.trang.study.ras)



nt.study.pred <- make.unipix.modified(nha.trang.study.ras)

nt.study.pred <- left_join(nt.study.pred, new_study_data, by =("pixID" = "pixID"))
       
nha.trang.study.ras_sero <- nha.trang.study.ras

values(nha.trang.study.ras_sero)<- nt.study.pred$predict_study_basic
plot(nha.trang.study.ras_sero, xlab= "Longitude", ylab ="Latitude")


#writeRaster(nha.trang.study.ras_sero, file = paste0(path, "Data/Output/10092020_sero_newcutoff.tiff"), format="GTiff", overwrite = TRUE)

