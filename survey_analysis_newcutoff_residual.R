
# lb <- "C:/R_packages" 
# .libPaths(lb)


library(tidyverse)
library(lubridate)
library(rjags)
library(R2jags)
library(runjags)
library(lattice)
library(coda)
library(survey)
library(splines)
library(sf)
library(rgdal)
library(raster)
library(SpatialDengue)
library(lrtest)


path = "C:/Users/User1/Desktop/Dengue_Vietnam/"


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
  mutate(head = ifelse(hlb02==1, 1,
                       ifelse(hlb02==2, 2, ifelse(hlb02==3 | hlb02 ==4, 3,
                                                        ifelse(hlb02==5 | hlb02 ==6, 4, ifelse(hlb02==7, 5,NA)))))) %>%
  mutate(head2 = ifelse(head==1 | head==2, 1, ifelse(head>2, 0, NA))) %>%
  mutate(head=as.factor(head)) %>%
  mutate(hhsize = ifelse(hlb05 < 0 , NA, hlb05)) %>%
  mutate(new_IgGPos = ifelse(IgG_panbio_units >=2.2, 1, 0), new_IgMPos = ifelse(IgM_panbio_units>= 8.1, 1, 0)) %>%
  #mutate(sero = ifelse(new_IgGPos == 1 | new_IgMPos == 1, 1,0 )) %>%
  mutate(sero = ifelse(new_IgGPos == 1, 1,0 )) %>%
  mutate(smoke.in = ifelse(hlb09_1d ==1 |hlb09_2d == 1 | hlb09_3d == 1 | hlb09_4d == 1, 1, 0)) %>%
  mutate(JE = ifelse(hla39==1 | hla39==2 | hla39==3, 1, ifelse(hla39==0, 0, ifelse(hla39==99, NA, NA)))) %>%
  mutate(chronic = ifelse(hla20 ==1, 1, ifelse(hla20 ==2, 0, NA))) %>%
  mutate(smoke = ifelse(hlb09a == 0, 1, ifelse(hlb09a==1, 0, NA))) %>%
  mutate(hosp = ifelse(hla21 == 1, 1, ifelse(hla21 ==2, 0, NA)))


gis_full2 <- read.csv(paste0(path, "/Data/gis_full2.csv"), fileEncoding = "latin1" )
gis_inerested <- gis_full2 %>%
  dplyr::select(hlperno, final_X, final_Y, lat, long ) %>%
  mutate(long1=lat, lat1=long) %>%
  dplyr::select(-lat, -long) %>%
  mutate(long=long1, lat=lat1)

surv <- left_join(surv, gis_inerested, by=c("hlperno"="hlperno"))



surv$risk = 1 - (exp(-1.254e-01 *surv$ageyears_f))  
surv$residual = surv$new_IgGPos - surv$risk

surv_sf <- st_as_sf(surv, coords = c("long", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#now need to match the cases to survey area

nha.trang.survey.ras <-raster(paste0(path, "Data/Maps/nha.trang.survey.ras.tif"))
nha.trang.survey.shp <- readOGR(paste0(path, "Data/Maps/nha.trang.survey.shp.shp"))


surv_int <- st_intersection(surv_sf, st_as_sf(nha.trang.survey.shp))

surv_int_df <- as.data.frame(surv_int)

write.csv(surv_int_df, file = paste0(path, "Data/surv_int.csv"))

des <- svydesign(id=~1,strata=~agegroup2, data=surv_int, fpc=~pop.x)

#for the two communes, 1 for each.
des.1 <- svydesign(id=~1,strata=~agegroup2, data=filter(surv_int, zoneid==1), fpc=~pop.x)
des.2 <- svydesign(id=~1,strata=~agegroup2, data=filter(surv_int, zoneid==2), fpc=~pop.x)


#chose JKn as the help file says use JKn for stratified samples
#jdes<- svydesign(id=~1,strata=~agegroup2, data=surv_int, fpc=~pop.x, type = "JKn")


#jdes <- as.svrepdesign(des, type="JKn")
bdes2 <- as.svrepdesign(des, type = "bootstrap", replicates = 10000)  

#for the two communes, 1 for each
bdes2.1 <- as.svrepdesign(des.1, type = "bootstrap", replicates = 10000)  
bdes2.2 <- as.svrepdesign(des.2, type = "bootstrap", replicates = 10000)  
#zone1 = phuoc
#zone 2 = hai


#zone 1 IgM and IgG population prevalence

100*svymean(~(new_IgMPos==1 & new_IgGPos==0), bdes2.1, na.rm = T)
100*confint(svymean(~(new_IgMPos==1 & new_IgGPos==0), bdes2.1, na.rm = T))

prev <- svymean(~(new_IgMPos==1 & new_IgGPos==1), bdes2.1, na.rm = T)
100*prev
100*confint(prev)

prev <- svymean(~(new_IgMPos==0 & new_IgGPos==1), bdes2.1, na.rm = T)
100*prev
100*confint(prev)

prev <- svymean(~(new_IgMPos==1 | new_IgGPos==1), bdes2.1, na.rm = T)
100*prev
100*confint(prev)



#zone 2 IgM and IgG population prevalence
prev <- svymean(~(new_IgMPos==1 & new_IgGPos==0), bdes2.2, na.rm = T)
100*prev
100*confint(svymean(~(new_IgMPos==1 & new_IgGPos==0), bdes2.2, na.rm = T))

prev <- svymean(~(new_IgMPos==1 & new_IgGPos==1), bdes2.2, na.rm = T)
100*prev
100*confint(prev)

prev <- svymean(~(new_IgMPos==0 & new_IgGPos==1), bdes2.2, na.rm = T)
100*prev
100*confint(prev)
#new_IgMPos == 0 & new_IgGPos == 1TRUE  76.968 0.0266
#new_IgMPos == 0 & new_IgGPos == 1TRUE  71.75743 82.17863


prev <- svymean(~(new_IgMPos==1 | new_IgGPos==1), bdes2.2, na.rm = T)
100*prev
100*confint(prev)





#putting in the bootstrap option directly in the svydesign argument works.
#bdes <- svydesign(id=~1,strata=~agegroup2, data=surv_int, fpc=~pop.x, type = "bootstrap", replicates=1000)

#bdes_1 <- svydesign(id=~1,strata=~agegroup2, data=filter(surv_int, commune == "DA"), fpc=~pop.x, type = "bootstrap")
#bdes_2 <- svydesign(id=~1,strata=~agegroup2, data=filter(surv_int, commune == "DB"), fpc=~pop.x, type = "bootstrap", replicates=10000)


sero<- svymean(~sero, bdes2,na.rm=TRUE)
100*sero
100*confint(sero)


#https://www.datasciencecentral.com/profiles/blogs/resampling-methods-comparison


hh_distribution <- data.frame(svytable(~hlb01, design = bdes))
write.csv(hh_distribution, "Output/hh_distribution.csv")

###

svyciprop(~sero, bdes2, method="logit")

  
  svyciprop(~I(sero==1), bdes2, method="logit") #method li produces NA


svyciprop(~sero, bdes2.1, method="logit")
svyciprop(~sero, bdes2.2, method="logit")

svyciprop(~I(agegroup==1), bdes2, method="li")

str( svyciprop(~I(agegroup %in% "E"), bdes2, method="lo") )

sapply( levels(factor(bdes2$variables$agegroup)),
        function(x){ 
          form <- as.formula( substitute( ~I(agegroup %in% x), list(x=x)))
          z <- svyciprop(form, bdes, method="lo", df=degf(bdes))
          c( 100*z, 100*c(attr(z,"ci")) )}  )


svyciprop(~I(agegroup %in% 1), bdes2, method="logit")

#non-weight adjusted


survey_mean(~I(agegroup %in% 1), bdes2)



100*svymean(~(factor(agegroup)), bdes2, na.rm = T)
100*confint(svymean(~(factor(agegroup)), bdes2, na.rm = T))


















#risk factors by positivity

female_rr <- svyglm(residual ~ female, design = bdes2, family = gaussian(link = "identity"))

female_rr
confint(female_rr)

binom.confint(sum(surv_int$female==0), dim(surv_int)[1],method="exact")

model_female_OR <-(svyglm(sero ~ female, design = des, family = quasibinomial()))
summary(model_female_OR)
exp(coef(model_female_OR))
exp(confint(model_female_OR))



m.female <- svycoxph(Surv(time=rep(1, 508), sero)~female, design = bdes2) 

100*svymean(~(factor(female)), bdes2, na.rm = T)
100*confint(svymean(~(factor(female)), bdes2, na.rm = T))

# model_female <-(svyglm(sero ~ female, design = bdes, family = quasibinomial()))
# summary(model_female)
# exp(coef(model_female))
# exp(confint(model_female))


# chronic condition

# table(surv_int$sero, surv_int$hla20)
# svyby (~hla20==1, by = ~sero, design = bdes, FUN = svymean, vartype = 'ci')
# 

binom.confint(sum(surv_int$chronic==1), dim(surv_int)[1],method="exact")


model_chronic <-(svyglm(sero ~ hla20, design = bdes2, family = quasibinomial()))

# summary(model_chronic)
# exp(coef(model_chronic))[-1]
# exp(confint(model_chronic))[-1, ]

m.chronic <- svyglm(residual ~ chronic, design = bdes2, family = gaussian(link = "identity"))
summary(m.chronic)
confint(m.chronic)



#ever hospitalized

# model_hosp <-(svyglm(sero ~ hla21, design = bdes, family = quasibinomial()))
# summary(model_hosp)
# exp(coef(model_hosp))[-1]
# exp(confint(model_hosp))[-1, ]

binom.confint(sum(surv_int_df$hosp==1), dim(surv_int_df)[1],method="exact")

m.hosp <- svyglm(residual ~ hosp, design = bdes2, family = gaussian(link = "identity"))
summary(m.hosp)
confint(m.hosp)

100*svymean(~(factor(hosp)), bdes2, na.rm = T)
100*confint(svymean(~(factor(hosp)), bdes2, na.rm = T))

#edu of head of household
# educaiton of head of housheold
# svyby (~noschool==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')
# svyby (~highschool==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')
# svyby (~postsecond==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')


# model_educ <-(svyglm(sero ~ educ, design = bdes, family = quasibinomial()))
# summary(model_educ)
# exp(coef(model_educ))
# exp(confint(model_educ))

binom.confint(sum(surv_int_df$educ==1), dim(surv_int_df)[1],method="exact")
binom.confint(sum(surv_int_df$educ==2), dim(surv_int_df)[1],method="exact")
binom.confint(sum(surv_int_df$educ==3), dim(surv_int_df)[1],method="exact")

m.educ <- svyglm(residual ~ as.factor(educ), design = bdes2, family = gaussian(link = "identity"))

summary(m.educ)
(confint(m.educ))
#no school is reference


# 
# model_educ2 <-(svyglm(sero ~ educ2, design = bdes, family = quasibinomial()))
# summary(model_educ2)
# exp(coef(model_educ2))
# exp(confint(model_educ2))

100*svymean(~(factor(educ)), bdes2, na.rm = T)

100*confint(svymean(~(factor(educ)), bdes2, na.rm = T))


#profession
#not much variation sooo not informative
# 
# svyby (~OfficeProfessional==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci', na.rm=TRUE)
# svyby (~Worker==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci', na.rm = TRUE)
# svyby (~Other==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci', na.rm=TRUE)

# model_occup <-(svyglm(sero ~ occup, design = bdes, family = quasibinomial()))
# summary(model_occup)
# exp(coef(model_occup))
# exp(confint(model_occup))
# 
# model_occup2 <-(svyglm(sero ~ occup2, design = bdes, family = quasibinomial()))
# summary(model_occup2)
# exp(coef(model_occup2))
# exp(confint(model_occup2))

binom.confint(sum(surv_int_df$occup==1, na.rm=TRUE), length(surv_int_df$occup[!is.na(surv_int_df$occup)]),method="exact")
binom.confint(sum(surv_int_df$occup==2, na.rm=TRUE), length(surv_int_df$occup[!is.na(surv_int_df$occup)]),method="exact")
binom.confint(sum(surv_int_df$occup==3, na.rm=TRUE), length(surv_int_df$occup[!is.na(surv_int_df$occup)]),method="exact")




m.occ <- svyglm(residual ~ as.factor(occup), design = bdes2, family = gaussian(link = "identity"))
summary(m.occ)
(confint(m.occ))

100*confint(svymean(~(occup), bdes2, na.rm = T))



# doses of JE vaccine
#japanese encephalitis vaccine


svyby (~hla39>=1, by = ~sero, design = bdes, FUN = svymean, vartype = 'ci', na.rm = TRUE)
 
# model_je <-(svyglm(sero ~ hla39, design = bdes, family = quasibinomial()))
# summary(model_je)
# exp(coef(model_je))[-1]
# exp(confint(model_je))[-1, ]

binom.confint(sum(surv_int_df$JE==1, na.rm=TRUE), length(surv_int_df$JE[!is.na(surv_int_df$JE)]),method="exact")


m.je <- svyglm(residual ~ JE, design = bdes2, family = gaussian(link = "identity"))
summary(m.je)
confint(m.je)


100*confint(svymean(~factor((JE)), bdes2))



# head of housheold
binom.confint(sum(surv_int_df$head==1, na.rm=TRUE), length(surv_int_df$head[!is.na(surv_int_df$head)]),method="exact")
binom.confint(sum(surv_int_df$head==2, na.rm=TRUE), length(surv_int_df$head[!is.na(surv_int_df$head)]),method="exact")
binom.confint(sum(surv_int_df$head==3, na.rm=TRUE), length(surv_int_df$head[!is.na(surv_int_df$head)]),method="exact")
binom.confint(sum(surv_int_df$head==4, na.rm=TRUE), length(surv_int_df$head[!is.na(surv_int_df$head)]),method="exact")
binom.confint(sum(surv_int_df$head==5, na.rm=TRUE), length(surv_int_df$head[!is.na(surv_int_df$head)]),method="exact")


svyby (~head_father==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')
svyby (~head_mother==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')
svyby (~head_gparent==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')
svyby (~head_other2==1, by = ~sero, design = des, FUN = svymean, vartype = 'ci')


# model_head <-(svyglm(sero ~ head, design = bdes, family = quasibinomial()))
# summary(model_head)
# exp(coef(model_head))[-1]
# exp(confint(model_head))[-1, ]


m.head <- svyglm(residual ~ as.factor(head), design = bdes2, family = gaussian(link = "identity"))
summary(m.head)
(confint(m.head))

100*confint(svymean(~factor((head)), bdes2))

svymean(~factor(head), design = bdes2)


#head
#1-father
#2-mother
#3-gparent
#4-other

#did not use head2 because combining grandparent and other removes the protective effect of gparent.




# mean log income

mean(surv_int_df$logincome)
quantile(surv_int_df$logincome, c(0.025, 0.975))

svyby (~logincome, by = ~sero, design = bdes, FUN = svymean, vartype = 'ci')

# model_income<-(svyglm(sero ~logincome, design = bdes, family = quasibinomial()))
# summary(model_income)
# exp(coef(model_income))[-1]
# exp(confint(model_income))[-1, ]

m.income <- svyglm(residual ~ logincome, design = bdes2, family = gaussian(link = "identity"))
summary(m.income)
(confint(m.income))

svymean(~((logincome)), bdes2, na.rm = T)

# hhsize
mean(log(surv_int_df$hhsize, na.rm=TRUE))
quantile(log(surv_int_df$hhsize), na.rm=TRUE, c(0.025, 0.975))


svyby (~hhsize, by = ~sero, design = des, FUN = svymean, vartype = 'ci', na.rm=TRUE)

# model_size<-(svyglm(sero ~ log(hhsize), design = des, family = quasibinomial()))
# summary(model_size)
# exp(coef(model_size))[-1]
# exp(confint(model_size))[-1, ]


m.hh <- svyglm(residual ~ log(hhsize), design = bdes2, family = gaussian(link = "identity"))
summary(m.hh)
(confint(m.hh))

svymean(~(log(hhsize)), bdes2, na.rm = T)
confint(svymean(~(log(hhsize)), bdes2, na.rm = T))


#kitchen inside or out -> only 1 outside so not usefeul



m.kit <- svycoxph(Surv(time=rep(1, 508), sero)~ log(hlb07), design = bdes2) 

#anyone smoke in house
#1 = no, so reverse tom ake 1 = yes
binom.confint(sum(surv_int_df$smoke==1, na.rm=TRUE), length(surv_int_df$smoke[!is.na(surv_int_df$smoke)]),method="exact")


m.smoke <- svyglm(residual ~ fct_rev(as.factor(hlb09a)), design = bdes2, family = gaussian(link = "identity"))
summary(m.smoke)
(confint(m.smoke))


str(surv$hlb09a)


m.smoke.in <- svycoxph(Surv(time=rep(1, 508), sero)~ smoke.in, design = bdes2) 
m.smoke.in
exp(confint(m.smoke.in))  

100*svymean(~ fct_rev(as.factor(hlb09a)), bdes2, na.rm = T)

100*confint(svymean(~ fct_rev(as.factor(hlb09a)), bdes2, na.rm = T))


# female, JE vaccine, head of household had an effect, now add spatial terms

multi <- svycoxph(Surv(time=rep(1, 508), sero)~ female + head + JE, design = bdes2) 
fem.head <- svycoxph(Surv(time=rep(1, 508), sero)~ female + head , design = bdes2) 

lrtest(multi, m.female)
anova.svycoxph

lrtest(fem.head, m.female)
fem.head
exp(confint(fem.head))  

summary(fem.head)

(fem.head, m.female)
  
  #now, assess spatial terms
# "pop"       "built"     "elevation" "lights"    "waterway"  "nature"   


m.pop <- svycoxph(Surv(time=rep(1, 508), sero)~ pop, design = bdes2) 
m.pop
exp(confint(m.pop))  


m.built <- svycoxph(Surv(time=rep(1, 508), sero)~ built, design = bdes2) 
m.built
exp(confint(m.built))  

m.elev<- svycoxph(Surv(time=rep(1, 508), sero)~ elevation, design = bdes2) 
m.elev
exp(confint(m.elev))  

m.lights<- svycoxph(Surv(time=rep(1, 508), sero)~ lights, design = bdes2) 
m.lights
exp(confint(m.lights))

m.wat<- svycoxph(Surv(time=rep(1, 508), sero)~ waterway, design = bdes2) 
m.wat
exp(confint(m.wat))  

m.nat<- svycoxph(Surv(time=rep(1, 508), sero)~ nature, design = bdes2) 
m.nat
exp(confint(m.nat))  

m.slope<- svycoxph(Surv(time=rep(1, 508), sero)~ slope, design = bdes2) 
m.slope
exp(confint(m.slope))  

m.spat<- svycoxph(Surv(time=rep(1, 508), sero)~ long1 + lat1, design = bdes2) 
m.spat
exp(confint(m.spat))  


m.m1 <- svycoxph(Surv(time=rep(1, 508), sero)~ female + long1 + lat1, design = bdes2) 
m.m1
exp(confint(m.m1))

anova.svycoxph(m.m1, m.female)
lrtest(m.m1, m.female)





;