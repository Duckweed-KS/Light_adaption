#summarise temp alt var
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
cov <- read.csv("Temp_weath_climate_alt_data.csv")
names(cov)

library(dplyr)

#check classes of all cols
sapply(cov, class)

cov$Sum2020_A.temp <- as.numeric(cov$Sum2020_A.temp)
cov$Sum2020_weath <- as.numeric(cov$Sum2020_weath)
cov$Aut2020_A.temp <- as.numeric(cov$Aut2020_A.temp)
cov$Aut2020_weath <- as.numeric(cov$Aut2020_weath)
cov$Sum2021_weath <- as.numeric(cov$Sum2021_weath)
cov$Aut2021_weath <- as.numeric(cov$Aut2021_weath)
cov$Spr2022_weath <- as.numeric(cov$Spr2022_weath)
cov$Sum2022_weath <- as.numeric(cov$Sum2022_weath)
cov$Annual.mean.temp <- as.numeric(cov$Annual.mean.temp)
cov$Diurnal.range <- as.numeric(cov$Diurnal.range)
cov$Isothermality <- as.numeric(cov$Isothermality)
cov$Temp.seasonality <- as.numeric(cov$Temp.seasonality)
cov$Max.temp <- as.numeric(cov$Max.temp)
cov$Min.temp <- as.numeric(cov$Min.temp)
cov$Temp.annual.range <- as.numeric(cov$Temp.annual.range)
cov$Mean.temp.wettest <- as.numeric(cov$Mean.temp.wettest)
cov$Mean.temp.driest <- as.numeric(cov$Mean.temp.driest)
cov$Mean.temp.warmest <- as.numeric(cov$Mean.temp.warmest)
cov$Mean.temp.coldest <- as.numeric(cov$Mean.temp.coldest)
cov$Annual.precip <- as.numeric(cov$Annual.precip)
cov$Precip.wettest <- as.numeric(cov$Precip.wettest)
cov$Precip.seasonality <- as.numeric(cov$Precip.seasonality)
cov$Precip.seaonsality <- as.numeric(cov$Precip.seasonality)
cov$Precip.wettest.quarter <- as.numeric(cov$Precip.wettest.quarter)
cov$Precip.seasonality.quarter <- as.numeric(cov$Precip.seasonality.quarter)
cov$Precip.warmest.quarter <- as.numeric(cov$Precip.warmest.quarter)
cov$Precip.coldest.quarter <- as.numeric(cov$Precip.coldest.quarter)
cov$Altitude <- as.numeric(cov$Altitude)


colMeans(cov[, 3:21], na.rm = TRUE)
colMeans(cov[, 24:44], na.rm = TRUE) 

rowMeans(cov[, 3:21], na.rm = TRUE) 
rowMeans(cov[, 24:44], na.rm = TRUE) 

rowSums(cov[, 3:21], na.rm = TRUE)
rowSums(cov[, 24:44], na.rm = TRUE)

#add light level forenv
#add grouping for light levels FOR 24 observations
cov$EnvLight <- c("dLL", "dLL", "dHL", "dLL", "dLL", "dLL", "dHL",
                   "dHL", "dHL", "dLL", "dHL", "dLL", "dLL", "dHL",
                   "dHL", "dLL", "dLL", "dLL", "dLL", "dHL", "dHL", "dHL",
                   "dLL", "dLL")

cov$EnvLight <- as.factor(cov$EnvLight)

names(cov)

#boxplots for hl ll sites
boxplot(Sum2020_A.temp~EnvLight,data=cov)
boxplot(Aut2020_A.temp~EnvLight,data=cov)
boxplot(Spr2021_A.temp~EnvLight,data=cov)
boxplot(Sum2021_A.temp~EnvLight,data=cov)
boxplot(Aut2021_A.temp~EnvLight,data=cov)
boxplot(Spr2022_A.temp~EnvLight,data=cov)
boxplot(Sum2022_A.temp~EnvLight,data=cov)

boxplot(Aut2020_W.temp~EnvLight,data=cov)
boxplot(Spr2021_W.temp~EnvLight,data=cov)
boxplot(Sum2021_W.temp~EnvLight,data=cov)
boxplot(Aut2021_W.temp~EnvLight,data=cov)
boxplot(Spr2022_W.temp~EnvLight,data=cov)
boxplot(Sum2022_W.temp~EnvLight,data=cov)

boxplot(Annual.mean.temp~EnvLight,data=cov)
boxplot(Diurnal.range~EnvLight,data=cov)
boxplot(Isothermality~EnvLight,data=cov)
boxplot(Temp.seasonality~EnvLight,data=cov)
boxplot(Max.temp~EnvLight,data=cov)
boxplot(Min.temp~EnvLight,data=cov)
boxplot(Temp.annual.range~EnvLight,data=cov)
boxplot(Mean.temp.wettest~EnvLight,data=cov)
boxplot(Mean.temp.driest~EnvLight,data=cov)
boxplot(Mean.temp.warmest~EnvLight,data=cov)
boxplot(Mean.temp.coldest~EnvLight,data=cov)
boxplot(Annual.precip~EnvLight,data=cov)
boxplot(Precip.wettest~EnvLight,data=cov)
boxplot(Precip.seasonality~EnvLight,data=cov)
boxplot(Precip.seaonsality~EnvLight,data=cov)
boxplot(Precip.wettest.quarter~EnvLight,data=cov)
boxplot(Precip.seasonality.quarter~EnvLight,data=cov)
boxplot(Precip.warmest.quarter~EnvLight,data=cov)
boxplot(Precip.coldest.quarter~EnvLight,data=cov)
boxplot(Altitude~EnvLight,data=cov)

#t tests
t.test(Sum2020_A.temp~EnvLight,data=cov)
t.test(Aut2020_A.temp~EnvLight,data=cov)
t.test(Spr2021_A.temp~EnvLight,data=cov)
t.test(Sum2021_A.temp~EnvLight,data=cov)
t.test(Aut2021_A.temp~EnvLight,data=cov)
t.test(Spr2022_A.temp~EnvLight,data=cov)
t.test(Sum2022_A.temp~EnvLight,data=cov)
#all n.s

t.test(Aut2020_W.temp~EnvLight,data=cov)
t.test(Spr2021_W.temp~EnvLight,data=cov)
t.test(Sum2021_W.temp~EnvLight,data=cov)
t.test(Aut2021_W.temp~EnvLight,data=cov)
t.test(Spr2022_W.temp~EnvLight,data=cov)
t.test(Sum2022_W.temp~EnvLight,data=cov)
#all n.s

t.test(Annual.mean.temp~EnvLight,data=cov)
t.test(Diurnal.range~EnvLight,data=cov)
t.test(Isothermality~EnvLight,data=cov)
t.test(Temp.seasonality~EnvLight,data=cov)
t.test(Max.temp~EnvLight,data=cov)
t.test(Min.temp~EnvLight,data=cov)
t.test(Temp.annual.range~EnvLight,data=cov)
t.test(Mean.temp.wettest~EnvLight,data=cov)
t.test(Mean.temp.driest~EnvLight,data=cov)
t.test(Mean.temp.warmest~EnvLight,data=cov)
t.test(Mean.temp.coldest~EnvLight,data=cov)
t.test(Annual.precip~EnvLight,data=cov)
t.test(Precip.wettest~EnvLight,data=cov)
t.test(Precip.seasonality~EnvLight,data=cov)
t.test(Precip.seaonsality~EnvLight,data=cov)
t.test(Precip.wettest.quarter~EnvLight,data=cov)
t.test(Precip.seasonality.quarter~EnvLight,data=cov)
t.test(Precip.warmest.quarter~EnvLight,data=cov)
t.test(Precip.coldest.quarter~EnvLight,data=cov)
t.test(Altitude~EnvLight,data=cov)
#all n.s

Summary <- cov %>%
  group_by(EnvLight) %>%
  summarise(Sum2020_A.temp_mean = mean(Sum2020_A.temp, na.rm = TRUE), Sum2020_A.temp_stdev = sd(Sum2020_A.temp, na.rm = TRUE), Sum2020_A.temp_maximum = max(Sum2020_A.temp, na.rm = TRUE))
Summary

Summary1 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Aut2020_A.temp_mean = mean(Aut2020_A.temp, na.rm = TRUE), Aut2020_A.temp_stdev = sd(Aut2020_A.temp, na.rm = TRUE), Aut2020_A.temp_maximum = max(Aut2020_A.temp, na.rm = TRUE))
Summary1

Summary2 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Spr2021_A.temp_mean = mean(Spr2021_A.temp, na.rm = TRUE), Spr2021_A.temp_stdev = sd(Spr2021_A.temp, na.rm = TRUE), Spr2021_A.temp_maximum = max(Spr2021_A.temp, na.rm = TRUE))
Summary2

Summary3 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Sum2021_A.temp_mean = mean(Sum2021_A.temp, na.rm = TRUE), Sum2021_A.temp_stdev = sd(Sum2021_A.temp, na.rm = TRUE), Sum2021_A.temp_maximum = max(Sum2021_A.temp, na.rm = TRUE))
Summary3

Summary4 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Aut2021_A.temp_mean = mean(Aut2021_A.temp, na.rm = TRUE), Aut2021_A.temp_stdev = sd(Aut2021_A.temp, na.rm = TRUE), Aut2021_A.temp_maximum = max(Aut2021_A.temp, na.rm = TRUE))
Summary4

Summary5 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Spr2022_A.temp_mean = mean(Spr2022_A.temp, na.rm = TRUE), Spr2022_A.temp_stdev = sd(Spr2022_A.temp, na.rm = TRUE), Spr2022_A.temp_maximum = max(Spr2022_A.temp, na.rm = TRUE))
Summary5

Summary6 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Sum2022_A.temp_mean = mean(Sum2022_A.temp, na.rm = TRUE), Sum2022_A.temp_stdev = sd(Sum2022_A.temp, na.rm = TRUE), Sum2022_A.temp_maximum = max(Sum2022_A.temp, na.rm = TRUE))
Summary6

Summary7 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Aut2020_W.temp_mean = mean(Aut2020_W.temp, na.rm = TRUE), Aut2020_W.temp_stdev = sd(Aut2020_W.temp, na.rm = TRUE), Aut2020_W.temp_maximum = max(Aut2020_W.temp, na.rm = TRUE))
Summary7

Summary8 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Spr2021_W.temp_mean = mean(Spr2021_W.temp, na.rm = TRUE), Spr2021_W.temp_stdev = sd(Spr2021_W.temp, na.rm = TRUE), Spr2021_W.temp_maximum = max(Spr2021_W.temp, na.rm = TRUE))
Summary8

Summary9 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Sum2021_W.temp_mean = mean(Sum2021_W.temp, na.rm = TRUE), Sum2021_W.temp_stdev = sd(Sum2021_W.temp, na.rm = TRUE), Sum2021_W.temp_maximum = max(Sum2021_W.temp, na.rm = TRUE))
Summary9

Summary10 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Aut2021_W.temp_mean = mean(Aut2021_W.temp, na.rm = TRUE), Aut2021_W.temp_stdev = sd(Aut2021_W.temp, na.rm = TRUE), Aut2021_W.temp_maximum = max(Aut2021_W.temp, na.rm = TRUE))
Summary10

Summary11 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Spr2022_W.temp_mean = mean(Spr2022_W.temp, na.rm = TRUE), Spr2022_W.temp_stdev = sd(Spr2022_W.temp, na.rm = TRUE), Spr2022_W.temp_maximum = max(Spr2022_W.temp, na.rm = TRUE))
Summary11

Summary12 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Sum2022_W.temp_mean = mean(Sum2022_W.temp, na.rm = TRUE), Sum2022_W.temp_stdev = sd(Sum2022_W.temp, na.rm = TRUE), Sum2022_W.temp_maximum = max(Sum2022_W.temp, na.rm = TRUE))
Summary12

Summary13 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Annual.mean.temp_mean = mean(Annual.mean.temp, na.rm = TRUE), Annual.mean.temp_stdev = sd(Annual.mean.temp, na.rm = TRUE), Annual.mean.temp_maximum = max(Annual.mean.temp, na.rm = TRUE))
Summary13

Summary14 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Diurnal.range_mean = mean(Diurnal.range, na.rm = TRUE),Diurnal.range_stdev = sd(Diurnal.range, na.rm = TRUE), Diurnal.range_maximum = max(Diurnal.range, na.rm = TRUE))
Summary14

Summary15 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Isothermality_mean = mean(Isothermality, na.rm = TRUE),Isothermality_stdev = sd(Isothermality, na.rm = TRUE), Isothermality_maximum = max(Isothermality, na.rm = TRUE))
Summary15

Summary16 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Temp.seasonality_mean = mean(Temp.seasonality, na.rm = TRUE),Temp.seasonality_stdev = sd(Temp.seasonality, na.rm = TRUE), Temp.seasonality_maximum = max(Temp.seasonality, na.rm = TRUE))
Summary16

Summary17 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Max.temp_mean = mean(Max.temp, na.rm = TRUE),Max.temp_stdev = sd(Max.temp, na.rm = TRUE), Max.temp_maximum = max(Max.temp, na.rm = TRUE))
Summary17

Summary18 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Min.temp_mean = mean(Min.temp, na.rm = TRUE),Min.temp_stdev = sd(Min.temp, na.rm = TRUE), Min.temp_maximum = max(Min.temp, na.rm = TRUE))
Summary18

Summary19 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Temp.annual.range_mean = mean(Temp.annual.range, na.rm = TRUE),Temp.annual.range_stdev = sd(Temp.annual.range, na.rm = TRUE), Temp.annual.range_maximum = max(Temp.annual.range, na.rm = TRUE))
Summary19

Summary20 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Mean.temp.wettest_mean = mean(Mean.temp.wettest, na.rm = TRUE),Mean.temp.wettest_stdev = sd(Mean.temp.wettest, na.rm = TRUE), Mean.temp.wettest_maximum = max(Mean.temp.wettest, na.rm = TRUE))
Summary20

Summary21 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Mean.temp.driest_mean = mean(Mean.temp.driest, na.rm = TRUE),Mean.temp.driest_stdev = sd(Mean.temp.driest, na.rm = TRUE), Mean.temp.driest_maximum = max(Mean.temp.driest, na.rm = TRUE))
Summary21

Summary22 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Mean.temp.warmest_mean = mean(Mean.temp.warmest, na.rm = TRUE),Mean.temp.warmest_stdev = sd(Mean.temp.warmest, na.rm = TRUE), Mean.temp.warmest_maximum = max(Mean.temp.warmest, na.rm = TRUE))
Summary22

Summary23 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Mean.temp.coldest_mean = mean(Mean.temp.coldest, na.rm = TRUE),Mean.temp.coldest_stdev = sd(Mean.temp.coldest, na.rm = TRUE), Mean.temp.coldest_maximum = max(Mean.temp.coldest, na.rm = TRUE))
Summary23

Summary24 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Annual.precip_mean = mean(Annual.precip, na.rm = TRUE),Annual.precip_stdev = sd(Annual.precip, na.rm = TRUE), Annual.precip_maximum = max(Annual.precip, na.rm = TRUE))
Summary24

Summary25 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.wettest_mean = mean(Precip.wettest, na.rm = TRUE),Precip.wettest_stdev = sd(Precip.wettest, na.rm = TRUE), Precip.wettest_maximum = max(Precip.wettest, na.rm = TRUE))
Summary25

Summary26 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.driest_mean = mean(Precip.driest, na.rm = TRUE),Precip.driest_stdev = sd(Precip.driest, na.rm = TRUE), Precip.driest_maximum = max(Precip.driest, na.rm = TRUE))
Summary26

Summary27 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.seasonality_mean = mean(Precip.seasonality, na.rm = TRUE),Precip.seasonality_stdev = sd(Precip.seasonality, na.rm = TRUE), Precip.seasonality_maximum = max(Precip.seasonality, na.rm = TRUE))
Summary27

Summary28 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.wettest.quarter_mean = mean(Precip.wettest.quarter, na.rm = TRUE),Precip.wettest.quarter_stdev = sd(Precip.wettest.quarter, na.rm = TRUE), Precip.wettest.quarter_maximum = max(Precip.wettest.quarter, na.rm = TRUE))
Summary28

Summary29 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.driest.quarter_mean = mean(Precip.driest.quarter, na.rm = TRUE),Precip.driest.quarter_stdev = sd(Precip.driest.quarter, na.rm = TRUE), Precip.driest_maximum = max(Precip.driest.quarter, na.rm = TRUE))
Summary29

Summary30 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.warmest.quarter_mean = mean(Precip.warmest.quarter, na.rm = TRUE),Precip.warmest.quarter_stdev = sd(Precip.warmest.quarter, na.rm = TRUE), Precip.warmest.quarter_maximum = max(Precip.warmest.quarter, na.rm = TRUE))
Summary30

Summary31 <- cov %>%
  group_by(EnvLight) %>%
  summarise(Precip.coldest.quarter_mean = mean(Precip.coldest.quarter, na.rm = TRUE),Precip.coldest.quarter_stdev = sd(Precip.coldest.quarter, na.rm = TRUE), Precip.coldest_maximum = max(Precip.coldest.quarter, na.rm = TRUE))
Summary31

Summ_TempClim <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5,
                      Summary6,Summary7,Summary8,Summary9,Summary10,Summary11,
                      Summary12,Summary13,Summary14,Summary15,Summary16,Summary17,
                      Summary18,Summary19,Summary20,Summary21,Summary22,Summary23,
                      Summary24,Summary25,Summary26,Summary27,Summary28,Summary29,
                      Summary30,Summary31)

names(cov)
