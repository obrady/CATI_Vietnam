
#recoded to fit whole data fit.

rm(list = ls())

require(ggplot2)
require(SpatialDengue)
require(geosphere)
require(automap)
require(raster)
library(tidyverse)
require(tmap)
require(sf)
require(RColorBrewer)





# load in final predictions

path = "C:/Users/User1/Desktop/Dengue_Vietnam/"


vndat <- read.csv(paste0(path, "Data/vndata_without_patch_2.csv"), header = TRUE) %>%
  dplyr::select(-X) 
# %>%
#   mutate(Week = Week - 1)


nha.trang.study.ras<-raster(paste0(files, "nha.trang.study.ras.tif"))
nha.trang.study.ras <- round(nha.trang.study.ras, 0)
unipix <- make.unipix(nha.trang.study.ras)
pixdistmat <- distm(cbind(unipix$x, unipix$y))

stim <-raster("10092020_sero_newcutoff.tif")
stim <- as.vector(stim)[unipix$pixID]

load("seasonal_vector.RData")
seasonal_vector <- incidence_est_2

seasonal_start <- 259







col_df <- matrix(NA, ncol = 2, nrow = 3)
colnames(col_df) <- c("median_dev_pred3","mean_dev_pred1_pred3" )
rownames(col_df) <- c("exponential", "gravity", "radiation")

pdf(paste0(files, "Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/fit_figures/8space/Cluster_overall_fits_16102020.pdf"),
    height = 15,
    width = 15)
par(mfrow = c(3,1))


#/Volumes/Work/H/Documents/Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/ABC_rounds/8space/Mov_model_1_fulldat_Round1_092020_trim.RData"

#change k from 3 to 2 as only 2 models area done
# for(i in 1:3){ # cluster loop
  for(k in 1:3){ # movement model loop
    
    Movement_types = c("Exponential","Gravity", "Radiation")
    # round 1
    load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round1_16102020.RData"))
    # round 2
    load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round2_16102020.RData"))
    # round 3
    load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round3_16102020.RData"))
    
    # recreate total devs
    Preds1$totaldevs = rowSums(Preds1[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
    Preds2$totaldevs = rowSums(Preds2[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
    Preds3$totaldevs = rowSums(Preds3[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
    Preds1$Preds_n = 1
    Preds2$Preds_n = 2
    Preds3$Preds_n = 3
    
    Preds <- rbind(Preds1, Preds2, Preds3) %>%
      mutate(Preds_n = as.factor(Preds_n))
    
    ##############
    ## Figure A: overall deviation improvement by SMC round for each cluster movement model combo
    ##############
    
    # limits
    x_range = c(rowSums(Preds1[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE),
                rowSums(Preds2[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE),
                rowSums(Preds3[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE))
                
    nbreaks = 20
    
    
    # ggplot(data = Preds, aes(totaldevs, fill = Preds_n, colour = Preds_n)) +
    #   geom_histogram( alpha=0.25, position = 'identity', bins = 20) +
    #   scale_fill_manual(values=c("red", "green","blue")) +
    #   ylim(0, 2500) +
    #   theme_classic() +
    #   xlab ("Deviation from data") + 
    #   ylab ("") +
    #   geom_vline(aes(xintercept = median(Preds1$totaldevs), col =" red", alpha = 0.25)) +
    #   geom_vline(aes(xintercept = median(Preds2$totaldevs), col ="green", alpha = 0.25)) +
    #   geom_vline(aes(xintercept = median(Preds3$totaldevs), col ="blue", alpha = 0.25)) +
    #   theme(legend.position = "none")
    #   
    #   
    #  geom_rangeframe(color = 'black') +
    #   theme_tufte()
    
    #red
    hist(Preds1$totaldevs, 
         col = rgb(1, 0, 0, 0.25),
         border = rgb(1, 0, 0, 0.25),
         xlab = "Deviation from data", 
         ylab = "",
         #main = paste0( "Movement model: ", Movement_types[[k]] ),
         main="",
         ylim = c(0, 2000), 
         #xlim = c(min(x_range), max(x_range)), 
         xlim = c(0,2.5),
         breaks = nbreaks)
    abline(v = median(Preds1$totaldevs), col = rgb(1, 0, 0, 0.25), lwd =2)
    text(2, 1000, 
         paste0("d = ", round(median(Preds1$totaldevs), 2)),
         col = rgb(1, 0, 0, 0.75))
    
    #green
    hist(Preds2$totaldevs,
         col = rgb(0, 1, 0, 0.25),
         border = rgb(0, 1, 0, 0.25),
         ylab = "",
         add = T,
         breaks = nbreaks)
    abline(v = median(Preds2$totaldevs), col = rgb(0, 1, 0, 0.25), lwd =2)
    text(2, 900,
         paste0("d = ", round(median(Preds2$totaldevs), 2)),
         col = rgb(0, 1, 0, 0.75))
    
    #blue
    hist(Preds3$totaldevs,
         col = rgb(0, 0, 1, 0.25),
         border = rgb(0, 0, 1, 0.25),
         ylab = "",
         add = T,
         breaks = nbreaks)
    abline(v = median(Preds3$totaldevs), col = rgb(0, 0, 1, 0.25), lwd =2)
    text(2, 800,
         paste0("d = ", round(median(Preds3$totaldevs), 2)),
         col = rgb(0, 0, 1, 0.75))

    # summary metrics
    #% reduction between round 1 and 3
    s_devs = round(100 *(mean(Preds1[, "spacedevs_R"], na.rm = TRUE) - mean(Preds3[, "spacedevs_R"], na.rm = TRUE)) / mean(Preds1[, "spacedevs_R"], na.rm = TRUE), 1)
    t_devs = round(100 *(mean(Preds1[, "timedevs_R"], na.rm = TRUE) - mean(Preds3[, "timedevs_R"], na.rm = TRUE)) / mean(Preds1[, "timedevs_R"], na.rm = TRUE), 1)
    o_devs = round(100 *(mean(Preds1$totaldevs) - mean(Preds3$totaldevs)) / mean(Preds1$totaldevs), 1)

    

    
    text(2, 500, paste0("Space = ", s_devs, "%"), xpd = T)
    text(2, 400, paste0("Time = ", t_devs, "%"), xpd = T)
    text(2, 300, paste0("Overall = ", o_devs, "%"), xpd = T)
    
    # collect
    col_df[k, 1] = median(Preds3$totaldevs)
    col_df[k, 2] = 100 *(mean(Preds1$totaldevs) - mean(Preds3$totaldevs)) / mean(Preds1$totaldevs)
  }
#}

mtext("A: Exponential movement", outer=TRUE,  cex=2, line=-3, adj = 0)
mtext("B: Gravity movement", outer=TRUE,  cex=2, line=-41,  adj = 0 )
mtext("C: Radiation movement", outer=TRUE,  cex=2, line=-79,  adj = 0)


dev.off()

write.csv(col_df, file = paste0(files, "Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/fit_figures/8space/Model_selection_table_16102020.csv"))

check <-     load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round1_16102020_8space_trim.RData"))



##############
## Figure B: parameter values
##############

pdf(paste0(files, "Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/fit_figures/8space/Parameter_overall_fits_16102020.pdf"),
    height = 15,
    width = 22)

par(mfrow = c(3,5))


#for(i in 1:3){ # cluster loop
    for(k in 1:3){ # movement model loop
 
    #     pdf(paste0(files, "Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/fit_figures/8space/Movmod_", k, "_04062020.pdf"),
    #      height = 3,
    #      width = 15)
    # par(mfrow = c(1,5))

    # round 1
    load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round1_16102020.RData"))
    # round 2
    load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round2_16102020.RData"))
    # round 3
    load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round3_16102020.RData"))


    nbreaks = 20
    
    Movement_types = c("Exponential", "Gravity", "Radiation")
    
    print(Movement_types[[k]])

    # hometime
    hist(Preds1$hometime,
         col = rgb(1, 0, 0, 0.25),
         border = rgb(1, 0, 0, 0.25),
         xlab = "Time spent at home",
         main = "",
         ylab = "",
         ylim = c(0, 800),
         breaks = nbreaks)
    hist(Preds2$hometime,
         col = rgb(0, 1, 0, 0.25),
         border = rgb(0, 1, 0, 0.25),
         add = T,
         breaks = nbreaks)
    hist(Preds3$hometime,
         col = rgb(0, 0, 1, 0.25),
         border = rgb(0, 0, 1, 0.25),
         add = T,
         breaks = nbreaks)
    #lines(seq(0, 1, length.out = 4000), rep(200, 4000), lty = 2)
    #4000/nbreaks (20) = 200
    
    

    # mospdeath
    hist(Preds1$mospdeath,
         col = rgb(1, 0, 0, 0.25),
         border = rgb(1, 0, 0, 0.25),
         xlab = "Daily mosquito mortality probability",
         main = "",
         ylab = "",
         ylim = c(0, 800),
         breaks = nbreaks)
    hist(Preds2$mospdeath,
         col = rgb(0, 1, 0, 0.25),
         border = rgb(0, 1, 0, 0.25),
         add = T,
         breaks = nbreaks)
    hist(Preds3$mospdeath,
         col = rgb(0, 0, 1, 0.25),
         border = rgb(0, 0, 1, 0.25),
         add = T,
         breaks = nbreaks)
    mospdeath_func <- function(x) rbeta(x, shape1 = 11.93328, shape2 = 107.3995)
    set.seed(1234)
    #lines(seq(0, 0.2, length.out = 1000),
    #      20*dbeta(seq(0, 0.2, length.out = 1000), shape1 = 11.93328, shape2 = 107.3995))

    # betaEnv mean
    hist(Preds1$betaEnv_mean,
         col = rgb(1, 0, 0, 0.25),
         border = rgb(1, 0, 0, 0.25),
         xlab = "Transmission coefficient (mean)",
         main = "",
         ylab = "",
         ylim = c(0, 800),
         breaks = nbreaks)
    hist(Preds2$betaEnv_mean,
         col = rgb(0, 1, 0, 0.25),
         border = rgb(0, 1, 0, 0.25),
         add = T,
         breaks = nbreaks)
    hist(Preds3$betaEnv_mean,
         col = rgb(0, 0, 1, 0.25),
         border = rgb(0, 0, 1, 0.25),
         add = T,
         breaks = nbreaks)
    #lines(seq(0.1, 10, length.out = 1000), rep(100, 1000), lty = 2)


    # betaEnv cor
    hist(Preds1$betaEnv_cor,
         col = rgb(1, 0, 0, 0.25),
         border = rgb(1, 0, 0, 0.25),
         xlab = "Transmission coefficient (correlation)",
         main = "",
         ylab = "",
         ylim = c(0, 800),
         breaks = nbreaks)
    hist(Preds2$betaEnv_cor,
         col = rgb(0, 1, 0, 0.25),
         border = rgb(0, 1, 0, 0.25),
         add = T,
         breaks = nbreaks)
    hist(Preds3$betaEnv_cor,
         col = rgb(0, 0, 1, 0.25),
         border = rgb(0, 0, 1, 0.25),
         add = T,
         breaks = nbreaks)
    #lines(seq(-1, 1, length.out = 50), rep(50, 50), lty = 2)


    # pdetec
    # pdetec_func <- function(x) rbeta(x, shape1 = 6.756883, shape2 = 100.2028) / 8
    # set.seed(1234)
    # 
    # hist(pdetec_func(1000),
    #     col = rgb(0, 0, 0, 0.25),
    #     border = rgb(0, 0, 0, 0.25),
    #     xlab = "Time spent at home",
    #    main = "",
    #    ylab = "",
    #     ylim = c(0, 50),
    #     xlim = c(0, 1))

    Preds1$pdetec[Preds1$pdetec > 1] = 1
    Preds2$pdetec[Preds2$pdetec > 1] = 1
    Preds3$pdetec[Preds3$pdetec > 1] = 1

    hist((Preds1$pdetec*8),
         col = rgb(1, 0, 0, 0.25),
         border = rgb(1, 0, 0, 0.25), #red
         main = "",
         ylab = "",
         xlab = "Probability of reporting",
         ylim =c(0, 800),
         xlim = c(0.07, 0.09),
         breaks = nbreaks)
    hist(Preds2$pdetec*8,
         col = rgb(0, 1, 0, 0.25), #green
         border = rgb(0, 1, 0, 0.25),
         add = T,
         breaks = nbreaks)
    hist(Preds3$pdetec*8,
         col = rgb(0, 0, 1, 0.25), #blue
         border = rgb(0, 0, 1, 0.25),
         add = T,
         breaks = nbreaks)
    
    pdetec_func <- function(x) rlnorm (x, meanlog = -2.28263583, sdlog = 0.01930149 ) / 8
    set.seed(1234)
    #lines(seq(0, 0.2, length.out = 500),
    #      5*dlnorm(seq(0, 0.2, length.out = 500), meanlog = -2.28263583, sdlog = 0.01930149))
    

    
 # }
    #dev.off()
}

mtext("A: Exponential movement", outer=TRUE,  cex=2, line=-3, adj = 0)
mtext("B: Gravity movement", outer=TRUE,  cex=2, line=-41,  adj = 0 )
mtext("C: Radiation movement", outer=TRUE,  cex=2, line=-79,  adj = 0)



dev.off()




#combine three pdfs with title






#summary of Preds3

k=3

load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round3_102020.RData"))

round(100*mean(Preds3$hometime),1)
round(100*quantile(Preds3$hometime, c(0.025, 0.95 )), 1)

round(100*mean(Preds3$mospdeath), 1)
round(100*quantile(Preds3$mospdeath, c(0.025, 0.95 )), 1)

round(mean(Preds3$betaEnv_mean), 2)
round(quantile(Preds3$betaEnv_mean, c(0.025, 0.95 )), 2)

round(mean(Preds3$betaEnv_cor), 2)
round(quantile(Preds3$betaEnv_cor, c(0.025, 0.95 )), 2)

8*round(100*mean(Preds3$pdetec), 1)
8*round(100*quantile(Preds3$pdetec, c(0.025, 0.95 )), 1)







##############
## Figure C: spatiotemporal plots of best fit for each cluster
##############

#best clsuter fits are:

mov_mod_type_opts = c(1)

source("DENSpatial_fitting_functions.R")



###########



unipix <- make.unipix(nha.trang.study.ras)
unipix<-round(unipix,0)
pixdistmat <- distm(cbind(unipix$x, unipix$y))


# load parameter values:
load(paste0("ABC_rounds/8space/Mov_model_", k, "_fulldat_Round3_102020_trim.RData"))
Preds3b$totaldevs = rowSums(Preds3b[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
Preds3b<- Preds3b[order(Preds3b$totaldevs),]
Preds3b <- Preds3b[1:1000,]
set.seed(42)
rows <- sample(nrow(Preds3b))
Preds3b <- Preds3b[rows, ]


#summary(Preds3b)
#

# load back in fitting data
#
#mov_mod_type = mov_mod_type_opts[i]
mov_mod_type = k

fitdat <-  vndat[vndat$Week>=100 & vndat$Week<=126,]

#pastdat <- vndat[vndat$Week<=79,]


fitdat <- patch_aggregate_cluster(fitdat)
colnames(fitdat)[4:5] = c("Longitude", "Latitude")

# pastdat <- patch_aggregate_cluster(pastdat)
# colnames(pastdat)[4:5] = c("Longitude", "Latitude")

weekdates <- c(min(fitdat$Week), max(fitdat$Week))

# implementing model with final parameter sets
nruns = 100
set.seed(123)


hold <- multi.model.run(Preds3b) # ran on cluster
save(hold, file = "analysis_save/temp_mod_1_simsave_8space_100runs.RData")

hold_2 <- multi.model.run(Preds3b) # ran on cluster
save(hold_2, file = "analysis_save/temp_mod_2_simsave_8space_100runs.RData")

hold_3 <- multi.model.run(Preds3b) #ran on cluster
save(hold_3, file = "analysis_save/temp_mod_3_simsave_8space_100runs.RData")

#####the above multi.model.run was run on cluster
####now look at results

fitdat_plot2 <- aggregate(fitdat$Number_of_cases, by = list(fitdat$patchID), FUN = sum)
fitdat_plot2$Longitude = fitdat$Longitude[match(fitdat_plot2$Group.1, fitdat$patchID)]
fitdat_plot2$Latitude = fitdat$Latitude[match(fitdat_plot2$Group.1, fitdat$patchID)]



fitdat_start <- fitdat[fitdat$Week>=81 & fitdat$Week<=85,]
fitdat_start2 <- aggregate(fitdat_start$Number_of_cases, by = list(fitdat_start$patchID), FUN = sum)
fitdat_start2$Longitude = fitdat_start$Longitude[match(fitdat_start2$Group.1, fitdat_start$patchID)]
fitdat_start2$Latitude = fitdat_start$Latitude[match(fitdat_start2$Group.1, fitdat_start$patchID)]



#name of this data is hold_1


load(paste0(path, "Output/hold_mov_1_simsave_8space_200runs.RData"))




ts_mat_1 <- hold_1[[1]]
space_mat_1 <- hold_1[[2]]
space_mat_peak_1 <- hold_1[[3]]


#plot.state(apply(space_mat, 1, median), sgpop, unipix)
p1 <- Sing.plot.FES(apply(space_mat_1, 1, median), nha.trang.study.ras, unipix, plotPoints = fitdat_plot2)


# plotPoints_start = st_as_sf(fitdat_start2, coords = c("Longitude", "Latitude"))
# plotPoints_start = st_set_crs(plotPoints_start, 4326)
# colnames(plotPoints_start)[2] = "Initial cases"
  

# p1 <- p1 + tm_shape(plotPoints_start) +
#   tm_symbols("Initial cases", col = "blueviolet", size = 0.5)

p1
tmap_save(p1, filename = "FES_mov1_seasonal_day226.png")

# r2:
cor_check <- data.frame(unipix$patchID,
                        apply(space_mat, 1, median),
                        fitdat_plot2$x[match(unipix$patchID, fitdat_plot2$Group.1)])
cor_check[is.na(cor_check[, 3]), 3] = 0
cor(cor_check[, 2], cor_check[, 3])^2 # 0.225787



patch_peak = rep(NA, nrow(unipix))
for(r in 1:length(unique(fitdat$patchID))){
  fitdat_f <- fitdat[fitdat$patchID == unique(fitdat$patchID)[r], ]
  # then aggregate and find the peak
  if(any(duplicated(fitdat_f$Week))){
    fitdat_f = aggregate(fitdat_f$Number_of_cases, by = as.list(fitdat_f$Week), FUN = sum)
  }
  
  patch_peak[unique(fitdat$patchID)[r]] = fitdat_f$Week[which.max(fitdat_f$Number_of_cases)]
}

# assembing values and building a unified scale
space_mat_peak[space_mat_peak == 1] = NA
vals <- apply(space_mat_peak, 1, median, na.rm = T) / 7 + min(fitdat$Week)


save_plots <- Sing.plot.PEAK(Observed = patch_peak,
                             Predicted = vals,
                             pras = nha.trang.study.ras,
                             unipix = unipix,
                             byscale = 20,
                             min_bracket = 0)

tmap_save(save_plots[[1]], filename = "analysis_save/8space/PEAK_OBS_Mod_1_1000runs_seasonal_day156.png")
tmap_save(save_plots[[2]], filename = "analysis_save/8space/PEAK_PRED_Mod_1_1000runs_seasonal_day156.png")




# weekly_fit <- fitdat %>%
#   group_by(Week) %>%
#   summarize(Numer_of_cases = sum(Number_of_cases))


#select from day 8 onwards as first week is the week used for initialization
daily_hold_1 <- ts_mat_1 %>%
  as.data.frame() %>%
  slice(8:196)

daily_hold_1_b <- ts_mat_1 %>%
  as.data.frame() %>%
  slice(8:196)

daily_hold_1$median =  apply(daily_hold_1_b, 1, median)
daily_hold_1$lower =  apply(daily_hold_1_b, 1, quantile, 0.025)
daily_hold_1$upper =  apply(daily_hold_1_b, 1, quantile, 0.975)



days_date <- read.csv(paste0(path, 'Data/days_date.csv'), header =T) %>%
  mutate(Date = as.Date(Date))



cases_daily <- vndat[vndat$Week>=100 & vndat$Week<=126,] %>%
  group_by(Date) %>%
  mutate(Number_of_cases=sum(Number_of_cases)) %>%
  ungroup() %>%
  distinct(Date, .keep_all=T) %>%
  mutate(Date = as.Date(Date)) %>%
  select(Date, Number_of_cases)


days_date_1 <- days_date %>%
  left_join(., cases_daily, by =c ("Date")) %>%
  mutate(Number_of_cases = if_else(is.na(Number_of_cases), 0, Number_of_cases)) %>%
  mutate(Estimate =daily_hold_1$median,
         Upper = daily_hold_1$upper,
         Lower = daily_hold_1$lower) %>%
  mutate(Mov = "Exponential")

#c("exponential", "gravity", "radiation")

#MODEL 2
load(paste0(path, "Output/hold_mov_2_simsave_8space_200runs.RData"))



ts_mat_2<- hold_2[[1]]
space_mat_2 <- hold_2[[2]]
space_mat_peak_2 <- hold_2[[3]]


#plot.state(apply(space_mat, 1, median), sgpop, unipix)
p2 <- Sing.plot.FES(apply(space_mat_2, 1, median), nha.trang.study.ras, unipix, plotPoints = fitdat_plot2)




# plotPoints_start = st_as_sf(fitdat_start2, coords = c("Longitude", "Latitude"))
# plotPoints_start = st_set_crs(plotPoints_start, 4326)
# colnames(plotPoints_start)[2] = "Initial cases"


# p1 <- p1 + tm_shape(plotPoints_start) +
#   tm_symbols("Initial cases", col = "blueviolet", size = 0.5)

p2
tmap_save(p2, filename = "FES_mov2_seasonal_day226.png")



#select from day 8 onwards as first week is the week used for initialization
daily_hold_2 <- ts_mat_2 %>%
  as.data.frame() %>%
  slice(8:196)

daily_hold_2_b <- ts_mat_2 %>%
  as.data.frame() %>%
  slice(8:196)

daily_hold_2$median =  apply(daily_hold_2_b, 1, median)
daily_hold_2$lower =  apply(daily_hold_2_b, 1, quantile, 0.025)
daily_hold_2$upper =  apply(daily_hold_2_b, 1, quantile, 0.975)



days_date_2 <- days_date %>%
  left_join(., cases_daily, by =c ("Date")) %>%
  mutate(Number_of_cases = if_else(is.na(Number_of_cases), 0, Number_of_cases)) %>%
  mutate(Estimate =daily_hold_2$median,
         Upper = daily_hold_2$upper,
         Lower = daily_hold_2$lower) %>%
 mutate(Mov = "Gravity")


 



###MODEL 3



load(paste0(path, "Output/hold_mov_3_simsave_8space_200runs.RData"))



ts_mat_3 <- hold_3[[1]]
space_mat_3 <- hold_3[[2]]
space_mat_peak_3 <- hold_3[[3]]


ts_mat_4 <- hold_all[[1]]
space_mat_4 <- hold_all[[2]]
space_mat_peak_4 <- hold_all[[3]]

#plot.state(apply(space_mat, 1, median), sgpop, unipix)
p3 <- Sing.plot.FES(apply(space_mat_3, 1, median), nha.trang.study.ras, unipix, plotPoints = fitdat_plot2)


p3
tmap_save(p3, filename = "FES_mov3_seasonal_day266.png")

# r2:
cor_check <- data.frame(unipix$patchID,
                        apply(space_mat, 1, median),
                        fitdat_plot2$x[match(unipix$patchID, fitdat_plot2$Group.1)])
cor_check[is.na(cor_check[, 3]), 3] = 0
cor(cor_check[, 2], cor_check[, 3])^2 # 0.225787






patch_peak = rep(NA, nrow(unipix))
for(r in 1:length(unique(fitdat$patchID))){
  fitdat_f <- fitdat[fitdat$patchID == unique(fitdat$patchID)[r], ]
  # then aggregate and find the peak
  if(any(duplicated(fitdat_f$Week))){
    fitdat_f = aggregate(fitdat_f$Number_of_cases, by = as.list(fitdat_f$Week), FUN = sum)
  }

  patch_peak[unique(fitdat$patchID)[r]] = fitdat_f$Week[which.max(fitdat_f$Number_of_cases)]
}

# assembing values and building a unified scale
space_mat_peak[space_mat_peak == 1] = NA
vals <- apply(space_mat_peak, 1, median, na.rm = T) / 7 + min(fitdat$Week)


save_plots <- Sing.plot.PEAK(Observed = patch_peak,
                             Predicted = vals,
                             pras = nha.trang.study.ras,
                             unipix = unipix,
                             byscale = 20,
                             min_bracket = 0)

tmap_save(save_plots[[1]], filename = "analysis_save/8space/PEAK_OBS_Mod_3_1000runs_seasonal_day156.png")
tmap_save(save_plots[[2]], filename = "analysis_save/8space/PEAK_PRED_Mod_3_1000runs_seasonal_day156.png")




#select from day 8 onwards as first week is the week used for initialization
daily_hold_3 <- ts_mat_3 %>%
  as.data.frame() %>%
  slice(8:196)

daily_hold_3_b <- ts_mat_3 %>%
  as.data.frame() %>%
  slice(8:196)

daily_hold_3$median =  apply(daily_hold_3_b, 1, median)
daily_hold_3$lower =  apply(daily_hold_3_b, 1, quantile, 0.025)
daily_hold_3$upper =  apply(daily_hold_3_b, 1, quantile, 0.975)



days_date_3 <- days_date %>%
  left_join(., cases_daily, by =c ("Date")) %>%
  mutate(Number_of_cases = if_else(is.na(Number_of_cases), 0, Number_of_cases)) %>%
  mutate(Estimate =daily_hold_3$median,
         Upper = daily_hold_3$upper,
         Lower = daily_hold_3$lower) %>%
  mutate(Mov = "Radiation")



days_date_tot = bind_rows(days_date_1, days_date_2, days_date_3)

ggplot(days_date_tot, aes(x = Date)) +
  geom_point(aes(y=Number_of_cases)) +
   geom_ribbon(aes(min=Lower, max=Upper, group = Mov, fill=Mov), alpha= 0.5) +
  geom_line(aes(y=Estimate, group = Mov, colour=Mov))
  
  


three_mov <- ggplot(days_date_tot, aes(x = Date)) +
  geom_ribbon(aes(min=Lower, max=Upper, fill=Mov), alpha= 0.3) +
  geom_line(aes(y=Estimate,  col = Mov), size = 1) +
  geom_point(aes(y=Number_of_cases), alpha =0.5) +
  scale_x_date(date_labels="%b %y",date_breaks  ="1 month") +
  facet_wrap(.~Mov, ncol=1) +
  scale_fill_discrete(name = "Movement type") +
  scale_colour_discrete(name = "Movement type") +
  labs(y="Number of cases") +
  theme_minimal()




ggsave(three_mov, file=paste0(path, "Plots/three_mov.jpg"))



#################################


mean(Preds3$hometime)
quantile(Preds3$hometime, c(0.025, 0.975))

mean(Preds3$mospdeath)
quantile(Preds3$mospdeath, c(0.025, 0.975))

mean(Preds3$betaEnv_mean)
quantile(Preds3$betaEnv_mean, c(0.025, 0.975))

mean(Preds3$betaEnv_cor)
quantile(Preds3$betaEnv_cor, c(0.025, 0.975))

mean(Preds3$pdetec) *8
quantile(Preds3$pdetec, c(0.025, 0.975)) *8


##########################################
# assigning parameter distributions back to the model
##########################################


# movement model 1
load(paste0("ABC_rounds/8space/Mov_model_", 1, "_fulldat_Round3_16102020_trim.RData"))
Preds3b$totaldevs = rowSums(Preds3b[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
Preds3b<- Preds3b[order(Preds3b$totaldevs),]
Preds3b <- Preds3b[1:1000,]

finalWeights = Preds3b

# movement model 2
load(paste0("ABC_rounds/8space/Mov_model_", 2, "_fulldat_Round3_16102020_trim.RData"))
Preds3b$totaldevs = rowSums(Preds3b[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
Preds3b<- Preds3b[order(Preds3b$totaldevs),]
Preds3b <- Preds3b[1:1000,]

finalWeights = rbind(finalWeights, Preds3b)

# movement model 3
load(paste0("ABC_rounds/8space/Mov_model_", 3, "_fulldat_Round3_16102020_trim.RData"))
Preds3b$totaldevs = rowSums(Preds3b[, c("timedevs_R", "spacedevs_R")], na.rm = TRUE)
Preds3b<- Preds3b[order(Preds3b$totaldevs),]
Preds3b <- Preds3b[1:1000,]

finalWeights = rbind(finalWeights, Preds3b)

# rescale finalweight total deviations sos they represent probabilities
finalWeights$totaldevs = 1 - finalWeights$totaldevs / max(finalWeights$totaldevs)

save(finalWeights, file = "finalWeights/16102020/finalWeights.rda")



# assignign to the DenSpatial package
setwd("FINAL_fits/")
library(devtools)
library(roxygen2)
use_data(finalWeights, overwrite = T)

# then reinstall
setwd("..")
install("SpatialDengue")

# and check
rm(list = ls())

### restart R ####

require(SpatialDengue)
data("finalWeights")
dim(finalWeights)

#summary(Preds3$pdetec)
#quantile(Preds3$pdetec, c(0.025, 0.975))

