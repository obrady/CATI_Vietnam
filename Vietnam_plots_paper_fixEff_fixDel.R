
library(tidyverse)
library(ggplot2)
library(raster)
library(geosphere)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(scales)
library(patchwork)
 

#######

path = "C:/Users/User1/Desktop/Dengue_Vietnam/"

vndat <- read.csv(paste0(path,"Data/vndata_without_patch.csv"), header = TRUE) %>%
  dplyr::select(-X) %>%
  mutate(Week = Week - 1)

fitdat <-  vndat[vndat$Week>79,]
pastdat <- vndat[vndat$Week<=79,]



base_259 <- get(load(paste0(path, "Output/base_day259.RData")))

steprun = 365

newD <- data.frame( median=base_259[[3]]$newD,
                    low = base_259[[1]]$newD,
                    high = base_259[[5]]$newD, 
                    days = c(1:365))

days_outbreak <- ((max(fitdat$Week) - min(fitdat$Week)) + 1)  * 7

newD_out <- newD %>%
  filter(days <= days_outbreak )



newD_hist <- newD %>%
  filter(days <= days_outbreak )%>%
group_by(x = ceiling(row_number()/7)) %>%
  summarise(median_cases = sum(median),
            low_cases = sum(low),
            high_cases = sum(high)) 
 





base_plot <- ggplot(newD, aes(x = days, y = median)) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, colour = NA, fill = "red") +
  theme(legend.position = "none") +
  theme_bw() +
  scale_x_continuous(limits = c(1, 365), breaks = c(1, 50, 100, 150, 200, 250, 300, 350)) +
  expand_limits(x = 1, y = 0) + 
  labs(x= "Days", y ="Median cases") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  geom_vline(xintercept=1, linetype = "dashed", alpha = 0.7) +
  geom_text(aes(x=1, label="Eary outbreak\n", y=10),angle=90, size = 3,check_overlap = TRUE) +
  geom_vline(xintercept=52, linetype = "dashed", alpha = 0.7) +
  geom_text(aes(x=52, label="Mid outbreak\n", y=10), angle=90, size = 3, check_overlap = TRUE)+
  geom_vline(xintercept=92, linetype = "dashed", alpha = 0.7) +
  geom_text(aes(x=92, label="Late outbreak\n", y=10), angle=90, size = 3, check_overlap = TRUE)

base_plot

#this goes to appendix
ggsave(base_plot, file = paste0(path, "Plots/base_plot.png"))





base_152 <- get(load(paste0(path, "Output/base_day152.RData")))

steprun = 365

newD_152 <- data.frame( median =base_152[[3]]$newD,
                    low = base_152[[1]]$newD,
                    high = base_152[[5]]$newD, 
                    days = c(1:365))



ggplot(newD_152, aes(x = days, y = median)) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, colour = NA, fill = "red") +
  theme(legend.position = "none") +
  theme_bw()






base_335 <- get(load(paste0(path, "Output/base_day335.RData")))

steprun = 365

newD_335 <- data.frame( median =base_335[[3]]$newD,
                        low = base_335[[1]]$newD,
                        high = base_335[[5]]$newD, 
                        days = c(1:365))



ggplot(newD_335, aes(x = days, y = median)) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, colour = NA, fill = "red") +
  theme(legend.position = "none") +
  theme_bw()




################### FROM HERE, PLOTS IN BODY OF PAPER #################################
### MDA drug with fix effecive coverage at 0.81 (with 7d delay)

MDA_drugs_str0_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start0_259_MDA_fixEff.RData")))
cases_MDA_drugs_str0_fixEff <- MDA_drugs_str0_fixEff[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 0, Radius = "MDA")

cases_MDA_drugs_str0<- rbind(cases_MDA_drugs_str0_fixEff %>% filter(Delays == 7), (cases_MDA_drugs_str0_fixEff %>% filter(Delays == 7))[rep(1, 3), ])
cases_MDA_drugs_str0$Length <- c(30, 60, 90, 180)




MDA_drugs_str52_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start52_259_MDA_fixEff.RData")))
cases_MDA_drugs_str52_fixEff <- MDA_drugs_str52_fixEff[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 52, Radius = "MDA")

cases_MDA_drugs_str52<- rbind(cases_MDA_drugs_str52_fixEff %>% filter(Delays == 7), (cases_MDA_drugs_str52_fixEff %>% filter(Delays == 7))[rep(1, 3), ])
cases_MDA_drugs_str52$Length <- c(30, 60, 90, 180)

MDA_drugs_str92_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start92_259_MDA_fixEff.RData"))) 
cases_MDA_drugs_str92_fixEff <- MDA_drugs_str92_fixEff[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 92, Radius = "MDA")

cases_MDA_drugs_str92<- rbind(cases_MDA_drugs_str92_fixEff %>% filter(Delays == 7), (cases_MDA_drugs_str92_fixEff %>% filter(Delays == 7))[rep(1, 3), ])
cases_MDA_drugs_str92$Length <- c(30, 60, 90, 180)



MDA_VC_str0_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_VC_start0_259_MDA_fixEff.RData"))) 
cases_MDA_VC_str0_fixEff <- MDA_VC_str0_fixEff[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 0, Radius = "MDA")

cases_MDA_VC_str0<- rbind(cases_MDA_VC_str0_fixEff %>% filter(Delays == 7), (cases_MDA_VC_str0_fixEff %>% filter(Delays == 7))[rep(1, 3), ])
cases_MDA_VC_str0$Length <- c(30, 60, 90, 180)




MDA_VC_str52_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_VC_start52_259_MDA_fixEff.RData"))) 
cases_MDA_VC_str52_fixEff <- MDA_VC_str52_fixEff[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 52, Radius = "MDA")

cases_MDA_VC_str52<- rbind(cases_MDA_VC_str52_fixEff %>% filter(Delays == 7), (cases_MDA_VC_str52_fixEff %>% filter(Delays == 7))[rep(1, 3), ])
cases_MDA_VC_str52$Length <- c(30, 60, 90, 180)

MDA_VC_str92_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_VC_start92_259_MDA_fixEff.RData"))) 
cases_MDA_VC_str92_fixEff <- MDA_VC_str92_fixEff[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 92, Radius = "MDA")

cases_MDA_VC_str92<- rbind(cases_MDA_VC_str92_fixEff %>% filter(Delays == 7), (cases_MDA_VC_str92_fixEff %>% filter(Delays == 7))[rep(1, 3), ])
cases_MDA_VC_str92$Length <- c(30, 60, 90, 180)







## MDA fixed del


MDA_drugs_str0_fixDel <- get(load(paste0(path,"Output/seasonal_CATI_drugs_start0_259_MDA_fixDel.RData"))) 
cases_MDA_drugs_str0_fixDel <- MDA_drugs_str0_fixDel[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 0, Radius = "MDA")



MDA_drugs_str52_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start52_259_MDA_fixDel.RData")))
cases_MDA_drugs_str52_fixDel <- MDA_drugs_str52_fixDel[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 52, Radius = "MDA")


MDA_drugs_str92_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start92_259_MDA_fixDel.RData"))) 
cases_MDA_drugs_str92_fixDel <- MDA_drugs_str92_fixDel[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 92, Radius = "MDA")




MDA_VC_str0_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_VC_start0_259_MDA_fixDel.RData"))) 
cases_MDA_VC_str0_fixDel <- MDA_VC_str0_fixDel[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 0, Radius = "MDA")



MDA_VC_str52_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_VC_start52_259_MDA_fixDel.RData"))) 
cases_MDA_VC_str52_fixDel <- MDA_VC_str52_fixDel[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 52, Radius = "MDA")


MDA_VC_str92_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_VC_start92_259_MDA_fixDel.RData"))) 
cases_MDA_VC_str92_fixDel <- MDA_VC_str92_fixDel[[2]] %>%
  filter(Radius == 1000) %>%
  mutate(start = 92, Radius = "MDA")
























base <- get(load(paste0(path, "Output/seasonal_CATI_drugs_Length_start0_259.RData")))
base <- base[[3]]

#LENGTH



CATI_drugs_str0 <- get(load(paste0(path, "Output/seasonal_CATI_drugs_Length_start0_259.RData")))
cases_drugs_str0 <- CATI_drugs_str0[[2]] %>%
  mutate(start = 0) %>%
  mutate(Radius = as.character(Radius))

CATI_drugs_str52 <- get(load(paste0(path, "Output/seasonal_CATI_drugs_Length_star52_259.RData")))
cases_drugs_str52 <- CATI_drugs_str52[[2]] %>%
  mutate(start = 52) %>%
  mutate(Radius = as.character(Radius))

CATI_drugs_str92 <- get(load(paste0(path, "Output/seasonal_CATI_drugs_Length_start92_259.RData")))
cases_drugs_str92 <- CATI_drugs_str92[[2]] %>%
  mutate(start = 92) %>%
  mutate(Radius = as.character(Radius))


CATI_VC_str0 <- get(load(paste0(path, "Output/seasonal_CATI_VC_Length_start0_259.RData")))
cases_VC_str0 <- CATI_VC_str0[[2]] %>%
  mutate(start = 0) %>%
  mutate(Radius = as.character(Radius))

CATI_VC_str52<- get(load(paste0(path, "Output/seasonal_CATI_VC_Length_star52_259.RData")))
cases_VC_str52 <- CATI_VC_str52[[2]] %>%
  mutate(start = 52) %>%
  mutate(Radius = as.character(Radius))

CATI_VC_str92 <- get(load(paste0(path, "Output/seasonal_CATI_VC_Length_start92_259.RData")))
cases_VC_str92 <- CATI_VC_str92[[2]] %>%
  mutate(start = 92) %>%
  mutate(Radius = as.character(Radius))






#Fixed effect

CATI_drugs_str0_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start0_259_fixEff.RData")))
cases_drugs_str0_fixEff <- CATI_drugs_str0_fixEff[[2]] %>%
  mutate(start = 0) %>%
  mutate(Radius = as.character(Radius))

CATI_drugs_str52_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start52_259_fixEff.RData")))
cases_drugs_str52_fixEff <- CATI_drugs_str52_fixEff[[2]] %>%
  mutate(start = 52) %>%
  mutate(Radius = as.character(Radius))

CATI_drugs_str92_fixEff <- get(load(paste0(path,"Output/seasonal_CATI_drugs_start92_259_fixEff.RData")))
cases_drugs_str92_fixEff <- CATI_drugs_str92_fixEff[[2]] %>%
  mutate(start = 92) %>%
  mutate(Radius = as.character(Radius))



CATI_VC_str0_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_VC_start0_259_fixEff.RData")))
cases_VC_str0_fixEff <- CATI_VC_str0_fixEff[[2]] %>%
  mutate(start = 0) %>%
  mutate(Radius = as.character(Radius))

CATI_VC_str52_fixEff <- get(load(paste0(path, "Output/seasonal_CATI_VC_start52_259_fixEff.RData")))
cases_VC_str52_fixEff <- CATI_VC_str52_fixEff[[2]] %>%
  mutate(start = 52) %>%
  mutate(Radius = as.character(Radius))

CATI_VC_str92_fixEff <- get(load(paste0(path,"Output/seasonal_CATI_VC_start92_259_fixEff.RData")))
cases_VC_str92_fixEff <- CATI_VC_str92_fixEff[[2]] %>%
  mutate(start = 92) %>%
  mutate(Radius = as.character(Radius))



##Fix del

CATI_drugs_str0_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start0_259_fixDel.RData")))
cases_drugs_str0_fixDel <- CATI_drugs_str0_fixDel[[2]] %>%
  mutate(start = 0) %>%
  mutate(Radius = as.character(Radius))

CATI_drugs_str52_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_drugs_start52_259_fixDel.RData")))
cases_drugs_str52_fixDel <- CATI_drugs_str52_fixDel[[2]] %>%
  mutate(start = 52) %>%
  mutate(Radius = as.character(Radius))

CATI_drugs_str92_fixDel <- get(load(paste0(path,"Output/seasonal_CATI_drugs_start92_259_fixDel.RData")))
cases_drugs_str92_fixDel <- CATI_drugs_str92_fixDel[[2]] %>%
  mutate(start = 92) %>%
  mutate(Radius = as.character(Radius))



CATI_VC_str0_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_VC_start0_259_fixDel.RData")))
cases_VC_str0_fixDel <- CATI_VC_str0_fixDel[[2]] %>%
  mutate(start = 0) %>%
  mutate(Radius = as.character(Radius))

CATI_VC_str52_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_VC_start52_259_fixDel.RData")))
cases_VC_str52_fixDel <- CATI_VC_str52_fixDel[[2]] %>%
  mutate(start = 52) %>%
  mutate(Radius = as.character(Radius))

CATI_VC_str92_fixDel <- get(load(paste0(path, "Output/seasonal_CATI_VC_start92_259_fixDel.RData")))
cases_VC_str92_fixDel <- CATI_VC_str92_fixDel[[2]] %>%
  mutate(start = 92) %>%
  mutate(Radius = as.character(Radius))






###########

cases_MDA_drugs_str0 <- cases_MDA_drugs_str0 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)  


cases_MDA_drugs_str52 <- cases_MDA_drugs_str52 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)  


cases_MDA_drugs_str92 <- cases_MDA_drugs_str92 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)  



cases_MDA_VC_str0 <- cases_MDA_VC_str0 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)  


cases_MDA_VC_str52 <- cases_MDA_VC_str52 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)  


cases_MDA_VC_str92 <- cases_MDA_VC_str92 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)  











  cases_drugs_str0 <- cases_drugs_str0 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_drugs_str52 <- cases_drugs_str52 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_drugs_str92 <- cases_drugs_str92 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)






cases_VC_str0 <- cases_VC_str0 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_VC_str52 <- cases_VC_str52 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_VC_str92 <- cases_VC_str92 %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_length <- bind_rows(bind_rows(cases_MDA_drugs_str0,
                                    cases_MDA_drugs_str52,
                                    cases_MDA_drugs_str92,
                                    cases_drugs_str0,
                                cases_drugs_str52,
                                cases_drugs_str92) %>% 
                            mutate(type = "Drug control"),
                          bind_rows(cases_MDA_VC_str0,
                                    cases_MDA_VC_str52,
                                    cases_MDA_VC_str92,
                                    cases_VC_str0,
                                    cases_VC_str52,
                                    cases_VC_str92) %>%
                            mutate(type = "Vector control"))


# print(MDA_drugs_fixEff_266[[3]])
#[1] 621.000 323.900 777.025


cases_length$Cases_median= 100*((CATI_drugs_str0[[3]][1] - cases_length$Cases_count_median) / CATI_drugs_str0[[3]][1])
cases_length$Cases_low= 100*((CATI_drugs_str0[[3]][2] - cases_length$Cases_count_low) / CATI_drugs_str0[[3]][2])
cases_length$Cases_high= 100*((CATI_drugs_str0[[3]][3] - cases_length$Cases_count_high) / CATI_drugs_str0[[3]][3])




head(cases_length)


cases_length$start <-factor(cases_length$start,
                                  levels = c(0, 52, 92),
                                  labels = c("Early outbreak", "Mid outbreak", "Late outbreak")) 




cases_length$Radius <-as.character(cases_length$Radius)




#Fix effect - 0.81




cases_MDA_drugs_str0_fixEff <- cases_MDA_drugs_str0_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_MDA_drugs_str52_fixEff <- cases_MDA_drugs_str52_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_MDA_drugs_str92_fixEff <- cases_MDA_drugs_str92_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_MDA_VC_str0_fixEff <- cases_MDA_VC_str0_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_MDA_VC_str52_fixEff <- cases_MDA_VC_str52_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_MDA_VC_str92_fixEff <- cases_MDA_VC_str92_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)




cases_drugs_str0_fixEff <- cases_drugs_str0_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_drugs_str52_fixEff <- cases_drugs_str52_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_drugs_str92_fixEff <- cases_drugs_str92_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)






cases_VC_str0_fixEff <- cases_VC_str0_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_VC_str52_fixEff <- cases_VC_str52_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_VC_str92_fixEff <- cases_VC_str92_fixEff %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_fixEff <- bind_rows(bind_rows(cases_MDA_drugs_str0_fixEff,
                                    cases_MDA_drugs_str52_fixEff,
                                    cases_MDA_drugs_str92_fixEff,
                                    cases_drugs_str0_fixEff,
                                    cases_drugs_str52_fixEff,
                                    cases_drugs_str92_fixEff) %>% 
                            mutate(type = "Drug control"),
                          bind_rows(cases_MDA_VC_str0_fixEff,
                                    cases_MDA_VC_str52_fixEff,
                                    cases_MDA_VC_str92_fixEff,
                                    cases_VC_str0_fixEff,
                                    cases_VC_str52_fixEff,
                                    cases_VC_str92_fixEff) %>%
                            mutate(type = "Vector control"))


# print(MDA_drugs_fixEff_266[[3]])
#[1] 621.000 323.900 777.025


cases_fixEff$Cases_median= 100*((CATI_drugs_str0[[3]][1] - cases_fixEff$Cases_count_median) / CATI_drugs_str0[[3]][1])
cases_fixEff$Cases_low= 100*((CATI_drugs_str0[[3]][2] - cases_fixEff$Cases_count_low) / CATI_drugs_str0[[3]][2])
cases_fixEff$Cases_high= 100*((CATI_drugs_str0[[3]][3] - cases_fixEff$Cases_count_high) / CATI_drugs_str0[[3]][3])

head(cases_length)


cases_fixEff$start <-factor(cases_fixEff$start,
                            levels = c(0, 52, 92),
                            labels = c("Early outbreak", "Mid outbreak", "Late outbreak")) 




cases_fixEff$Radius <-as.character(cases_fixEff$Radius)




###### fixed del at 7 days




cases_MDA_drugs_str0_fixDel <- cases_MDA_drugs_str0_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_MDA_drugs_str52_fixDel <- cases_MDA_drugs_str52_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_MDA_drugs_str92_fixDel <- cases_MDA_drugs_str92_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)






cases_MDA_VC_str0_fixDel <- cases_MDA_VC_str0_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_MDA_VC_str52_fixDel <- cases_MDA_VC_str52_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_MDA_VC_str92_fixDel <- cases_MDA_VC_str92_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)




cases_drugs_str0_fixDel <- cases_drugs_str0_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_drugs_str52_fixDel <- cases_drugs_str52_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_drugs_str92_fixDel <- cases_drugs_str92_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)






cases_VC_str0_fixDel <- cases_VC_str0_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_VC_str52_fixDel <- cases_VC_str52_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)


cases_VC_str92_fixDel <- cases_VC_str92_fixDel %>%
  rename(Cases_count_median = Cases_median,
         Cases_count_low = Cases_low,
         Cases_count_high = Cases_high)



cases_fixDel <- bind_rows(bind_rows(cases_MDA_drugs_str0_fixDel,
                                    cases_MDA_drugs_str52_fixDel,
                                    cases_MDA_drugs_str92_fixDel,
                                    cases_drugs_str0_fixDel,
                                    cases_drugs_str52_fixDel,
                                    cases_drugs_str92_fixDel) %>% 
                            mutate(type = "Drug control"),
                          bind_rows(cases_MDA_VC_str0_fixDel,
                                    cases_MDA_VC_str52_fixDel,
                                    cases_MDA_VC_str92_fixDel,
                                    cases_VC_str0_fixDel,
                                    cases_VC_str52_fixDel,
                                    cases_VC_str92_fixDel) %>%
                            mutate(type = "Vector control"))


#print(MDA_drugs_fixDel_266[[3]])
#[1] 621.000 323.900 777.025


cases_fixDel$Cases_median= 100*((CATI_drugs_str0[[3]][1] - cases_fixDel$Cases_count_median) / CATI_drugs_str0[[3]][1])
cases_fixDel$Cases_low= 100*((CATI_drugs_str0[[3]][2] - cases_fixDel$Cases_count_low) / CATI_drugs_str0[[3]][2])
cases_fixDel$Cases_high= 100*((CATI_drugs_str0[[3]][3] - cases_fixDel$Cases_count_high) / CATI_drugs_str0[[3]][3])

head(cases_fixDel)


cases_fixDel$start <-factor(cases_fixDel$start,
                            levels = c(0, 52, 92),
                            labels = c("Early outbreak", "Mid outbreak", "Late outbreak")) 












cases_fixDel









####

update_df_2 <- function(pb){
  
  
  vals = pb[, 15:17]
  pb$Cases_median = apply(vals, 1, median)
  pb$Cases_low = apply(vals, 1, min)
  pb$Cases_high = apply(vals, 1, max)
  
  pb <- pb %>%
    mutate(Cases_low = ifelse(Cases_low <0, 0, Cases_low))
  
  pb$Radius <- factor(pb$Radius,
                      levels = c("50","100","200", "500", "1000", "MDA"),
                      labels = c("50m", "100m", "200m", "500m", "1000m", "MA"))
  

  
  return(pb)
}




update_df_3 <- function(pb){
  
  
  vals = pb[, 23:25]
  pb$Cases_median = apply(vals, 1, median)
  pb$Cases_low = apply(vals, 1, min)
  pb$Cases_high = apply(vals, 1, max)
  
  pb <- pb %>%
    mutate(Cases_low = ifelse(Cases_low <0, 0, Cases_low))
  
  pb$Radius <- factor(pb$Radius,
                      levels = c("50","100","200", "500", "1000", "MDA"),
                      labels = c("50m", "100m", "200m", "500m", "1000m", "MDA"))
  
  
  
  return(pb)
}


cases_length_2 <- cases_length %>%
  mutate(Cases_median = ifelse(Cases_median <0, 0, Cases_median)) %>%
  mutate(Cases_low=ifelse(Cases_low<0, 0, Cases_low)) %>%
  mutate(Cases_high = ifelse(Cases_high<0, 0, Cases_high))


cases_length_2 <- update_df_3(cases_length_2)


#colour codes
#lightest to darkest
#94E8B4
#72BDA3
#5E8C61
#4E6151
#3B322C

myColours <- c("#94E8B4", "#72BDA3", "#5E8C61", "#4E6151", "#3B322C", "orange")

names(myColours) <- levels(cases_length_2$Radius)

colScale <- scale_colour_manual(name = "Radius",values = myColours)
colScale_2 <- scale_fill_manual(name = "Radius",values = myColours)

colpal = brewer.pal(length(radvec), "Accent")
colpalpoly = addalpha(brewer.pal(length(radvec), "Accent"), 0.3)


# length_intervention<- ggplot(cases_length_2, aes(x = Length, y = Cases_median, group = Radius, color = Radius, fill = Radius)) +
#   geom_point() +
#   scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
#   geom_errorbar(aes(ymin = Cases_low, ymax = Cases_high, group = Radius, color = Radius), alpha = 0.2) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.position = "top") +
#   labs(x = "Length of intervention (days)",
#        y = "% of cases averted") +
#   facet_grid(type ~ start) +
#   # colScale +
#   # colScale_2 +
# guides(color=guide_legend(nrow=1))
# 
# 


length_intervention_2<- ggplot(cases_length_2,
                               aes(x = Length, y = Cases_median, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin = Cases_low, ymax = Cases_high, group = Radius, color = Radius), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "% of cases averted") +
  facet_grid(type ~ start)+
  # colScale +
  # colScale_2 +
 guides(color=guide_legend(nrow=1))



length_intervention_2


ggsave(length_intervention_2, file=paste0(path, "Plots/2023_length_delay7_noVC_MA.png"))



#### Fixed Effect
# 0.81


cases_fixEff_2 <- cases_fixEff%>%
  mutate(Cases_median = ifelse(Cases_median <0, 0, Cases_median)) %>%
  mutate(Cases_low=ifelse(Cases_low<0, 0, Cases_low)) %>%
  mutate(Cases_high = ifelse(Cases_high<0, 0, Cases_high))


cases_fixEff_2 <- update_df_3(cases_fixEff_2)


fixEff_intervention<- ggplot(cases_fixEff_2, aes(x = as.numeric(Delays), y = Cases_median, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 2, 4, 7, 14)) +
  geom_errorbar(aes(ymin = Cases_low, ymax = Cases_high, group = Radius, color = Radius), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Delay in intervention (days)",
       y = "% of cases averted") +
  facet_grid(type ~ start)+
  # colScale +
  # colScale_2 +
 guides(color=guide_legend(nrow=1))


fixEff_intervention

ggsave(fixEff_intervention, file=paste0(path, "Plots/2023_fixEff_0.81.png"))



fixEff_intervention_2<- ggplot(cases_fixEff_2, aes(x = as.numeric(Delays), y = Cases_median, group = Radius, color = Radius, fill = Radius)) +
                                 geom_point() +
                                 scale_x_continuous(breaks = c(0, 2, 4, 7, 14)) +
                                 geom_errorbar(aes(ymin = Cases_low, ymax = Cases_high, group = Radius, color = Radius), alpha = 0.2) +
                                 theme_bw() +
                                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                                 theme(legend.position = "top") +
                                 labs(x = "Delay in intervention (days)",
                                      y = "% of cases averted") +
                                 facet_grid(type ~ start)+
                                  # colScale +
                                  # colScale_2 +
                                  guides(color=guide_legend(nrow=1))


fixEff_intervention_2

ggsave(fixEff_intervention_2, file=paste0(path, "Plots/2023_fixEff_0.81.png"))







####fix delay




cases_fixDel_2 <- cases_fixDel%>%
  mutate(Cases_median = ifelse(Cases_median <0, 0, Cases_median)) %>%
  mutate(Cases_low=ifelse(Cases_low<0, 0, Cases_low)) %>%
  mutate(Cases_high = ifelse(Cases_high<0, 0, Cases_high))


cases_fixDel_2 <- update_df_3(cases_fixDel_2)

fixDel_intervention_2<- ggplot(cases_fixDel_2, aes(x = as.numeric(CrudeEFFcov), y = Cases_median, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8)) +
  geom_errorbar(aes(ymin = Cases_low, ymax = Cases_high, group = Radius, color = Radius), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Effective coverage",
       y = "% of cases averted") +
  facet_grid(type ~ start)+
  # colScale +
  # colScale_2 +
  guides(color=guide_legend(nrow=1))


fixDel_intervention_2

ggsave(fixDel_intervention_2, file=paste0(path, "Plots/2023_fixDel_7days_noVC_MA.png"))




###### Patches treated


cases_length_3 <- cases_length_2 %>%
  mutate(Patches_median = ifelse(type == "Drug control", Patches_median_drug, Patches_median_vec)) %>%
  mutate(Patches_low = ifelse(type == "Drug control", Patches_low_drug, Patches_low_vec)) %>%
  mutate(Patches_high = ifelse(type == "Drug control", Patches_high_drug, Patches_high_vec))


patches_treated_plot <- ggplot(cases_length_3, aes(x = as.numeric(Length), y = Patches_median, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin = Patches_low, ymax = Patches_high, group = Radius, color = Radius), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "# patches treated") +
  facet_grid(type  ~ start)+
  colScale +
  colScale_2 +
  guides(color=guide_legend(nrow=1))

patches_treated_plot

ggsave(patches_treated_plot, file=paste0(path, "Plots/2023_patches_treated_length_delay7.png"))



early_out_patches <- ggplot(cases_length_3 %>% filter(type=="Drug control" & start == "Early outbreak"), aes(x = as.numeric(Length), y = Patches_median, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin = Patches_low, ymax = Patches_high, group = Radius, color = Radius), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "# patches treated") +
  facet_grid(type  ~ start)+
  # colScale +
  # colScale_2 +
  guides(color=guide_legend(nrow=1))

ggsave(early_out_patches, file=paste0(path, "Plots/2023_early_out_patches_delay7.png"))






# 
# 
# 
#   pop_treated_plot <- ggplot(cases_drugs_length_2, aes(x = as.numeric(Length), y = Poptot_median, group = Radius, color = Radius, fill = Radius)) +
#     geom_line() +
#     scale_x_continuous(breaks = c(0, 30, 60, 90, 180, 360)) +
#     geom_ribbon(aes(ymin = Poptot_low, ymax = Poptot_high), alpha = 0.2, colour = NA) +
#     theme_bw() +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#     theme(legend.position = "top") +
#     labs(x = "Length of intervention (days)",
#          y = "# people treated") +
#     
#     
#     facet_grid(~ start) +
#     guides(color=guide_legend(nrow=1))
# 
#   pop_treated_plot
# 
#   ggsave(pop_treated_plot, file="pop_treated_length_delay7.png")
#   
# 
# 
# comb_fixEff<-ggarrange(drugs_fixEff, VC_fixEff,
#                        drugs_fixdel, VC_fixdel,
#                        labels = c("A", "B", "C", "D"),
#                        ncol = 2, nrow = 2,
#                        common.legend = TRUE, legend = "bottom"
#                        )
# 
# 
# 
# 
# ggexport(comb_fixEff, filename = "/Volumes/Work/H/Documents/Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/analysis_save/Simulations/sim_comparison.pdf",
#          width = 8, height = 6)
# 
# ggexport(comb_fixEff, filename = "/Volumes/Work/H/Documents/Yalda_dengue_data/Oliver_model/vietnam run/To run on cluster/analysis_save/Simulations/sim_comparison.png",
#          width = 1202, height = 930, res = 175)
# 




#calculate ratio of cases averted to patches treated


CATI_drugs_str0[[3]][1] 


cases_length_4 <- cases_length_2 %>%
  mutate(Poptot_median = ifelse(type == "Drug control", Poptot_median_drug, Poptot_median_vec)) %>%
  mutate(Poptot_low = ifelse(type == "Drug control", Poptot_low_drug, Poptot_low_vec)) %>%
  mutate(Poptot_high = ifelse(type == "Drug control", Poptot_high_drug, Poptot_high_vec))


cases_length_4$Ratio_median <-(CATI_drugs_str0[[3]][1]- cases_length_4$Cases_count_median) /cases_length_4$Poptot_median   
cases_length_4$Ratio_low <- (CATI_drugs_str0[[3]][2]- cases_length_4$Cases_count_low) /cases_length_4$Poptot_low
cases_length_4$Ratio_high <- (CATI_drugs_str0[[3]][3]- cases_length_4$Cases_count_high) /cases_length_4$Poptot_high




vals = cases_length_4[, 29:31]
cases_length_4$Ratio_median = apply(vals, 1, median)
cases_length_4$Ratio_low = apply(vals, 1, min)
cases_length_4$Ratio_high = apply(vals, 1, max)



cases_length_4$Ratio_high <- ifelse(cases_length_4$Ratio_high<0, 0, cases_length_4$Ratio_high)
cases_length_4$Ratio_low <- ifelse(cases_length_4$Ratio_low<0, 0, cases_length_4$Ratio_low)
cases_length_4$Ratio_median <- ifelse(cases_length_4$Ratio_median<0, 0, cases_length_4$Ratio_median)




cases_length_4$Ratio_median2 <-cases_length_4$Poptot_median / (CATI_drugs_str0[[3]][1]- cases_length_4$Cases_count_median) 
cases_length_4$Ratio_low2 <- cases_length_4$Poptot_low / (CATI_drugs_str0[[3]][2]- cases_length_4$Cases_count_low) 
cases_length_4$Ratio_high2 <- cases_length_4$Poptot_high / (CATI_drugs_str0[[3]][3]- cases_length_4$Cases_count_high) 



vals = cases_length_4[, 32:34]
cases_length_4$Ratio_median2 = apply(vals, 1, median)
cases_length_4$Ratio_low2 = apply(vals, 1, min)
cases_length_4$Ratio_high2 = apply(vals, 1, max)



cases_length_4$Ratio_high2 <- ifelse(cases_length_4$Ratio_high2<0, 0, cases_length_4$Ratio_high2)
cases_length_4$Ratio_low2 <- ifelse(cases_length_4$Ratio_low2<0, 0, cases_length_4$Ratio_low2)
cases_length_4$Ratio_median2 <- ifelse(cases_length_4$Ratio_median2<0, 0, cases_length_4$Ratio_median2)









# 
# ggplot(cases_length_4 %>% 
#          filter(type == "Drug control", 
#        start == "Early outbreak")
#        , aes(x = as.numeric(Length), y = Ratio_median, group = Radius, color = Radius, fill = Radius)) +
#   geom_line() +
#   scale_x_continuous(breaks = c(0, 30, 60, 90, 180, 360)) +
#   geom_ribbon(aes(ymin = Ratio_low, ymax = Ratio_high), alpha = 0.2, colour = NA) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.position = "top") +
#   labs(x = "Length of intervention (days)",
#        y = "cases averted/ dose") +
#   facet_grid(type  ~ start)+
#   colScale +
#   colScale_2 +
#   guides(color=guide_legend(nrow=1))
# 






ppl_treated_averted <- ggplot(cases_length_4 %>%
         filter(type == "Drug control", 
                start == "Early outbreak") , aes(x = as.numeric(Length), y = Ratio_median2, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin = Ratio_low2, ymax = Ratio_high2), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "Treatment epidoses / cases averted") +
  # facet_grid(type  ~ start)+
  # colScale +
  # colScale_2 +
  guides(color=guide_legend(nrow=1))


ppl_treated_averted

ggsave(ppl_treated_averted, file=paste0(path_new, "Plots/2023_ratio_ppl_treated_avertedtreated_length_delay7_drug_early_v2.png"))


ppl_treated <- ggplot(cases_length_4 %>%
         filter(type == "Drug control", 
                start == "Early outbreak") , aes(x = as.numeric(Length), y =Poptot_median , group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin =Poptot_low, ymax = Poptot_high), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "Treatment episodes") +
  # colScale +
  # colScale_2 +
 scale_y_continuous(labels = comma) +  guides(color=guide_legend(nrow=1))





ppl_treated

ggsave(ppl_treated, file=paste0(path_new, "Plots/2023_ppl_treated_length_delay7_drug_v2.png"))


ppl_treated_averted_2 <- ggplot(cases_length_4 %>%
                                filter(type == "Drug control", 
                                       start == "Early outbreak") , aes(x = as.numeric(Length), y = Ratio_median2, group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin = Ratio_low2, ymax = Ratio_high2), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "Treatment epidoses / cases averted") +
  # colScale +
  # colScale_2 +
  guides(color=guide_legend(nrow=1))

  # ggtitle("B") +
  

  # facet_grid(type  ~ start)+

#guides(color=guide_legend(nrow=1))

ppl_treated_averted_2

ppl_treated_2 <- ggplot(cases_length_4 %>%
                        filter(type == "Drug control", 
                               start == "Early outbreak") , aes(x = as.numeric(Length), y =Poptot_median , group = Radius, color = Radius, fill = Radius)) +
  geom_point() +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 180)) +
  geom_errorbar(aes(ymin =Poptot_low, ymax = Poptot_high), alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "top") +
  labs(x = "Length of intervention (days)",
       y = "Treatment episodes") +
  scale_y_continuous(labels = comma)+
  # colScale +
  # colScale_2 +
  # ggtitle("A") +
  guides(color=guide_legend(nrow=1))
    

patch <- ppl_treated_2 + ppl_treated_averted_2 + plot_annotation(tag_levels = "A") +
plot_layout(guide = "collect")  & theme(legend.position = "top")

patch

ggsave(patch, file = paste0(path, "Plots/2023_pop_comp.png"))




