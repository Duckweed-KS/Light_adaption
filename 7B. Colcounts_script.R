#script to explore hl ll dw growth rate data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Biomass - start weights.csv")
grod3 <- read.csv("Col_countsRGR_nospREM.csv") #just to use
#without 4 species included

#dont need to normalise as doing col gain in this script
#need to take start bio away from end if working with fw

tail(gr)
head(gr)
names(gr)
#names
#[1] "Accession"  "T0_cols"    "T7_cols"    "T14_cols"   "T21_cols"   "FW..mg."   
#[7] "FDW..mg."   "Start_date" "Rep"        "Treatment"  "Species"    "X"         
#[13] "X.1"        "X.2"  

#check classes of all cols
sapply(gr, class)
#make numeric
gr$T0_cols <- as.numeric(as.character(gr$T0_cols))
gr$T7_cols <- as.numeric(as.character(gr$T7_cols))
gr$T14_cols <- as.numeric(as.character(gr$T14_cols))
gr$T21_cols <- as.numeric(as.character(gr$T21_cols))
class(gr$T0_cols)

#need to remove last 3 cols as empty
#gr <- gr[-c(12:14)]

library(dplyr)
library(stringr)
library(ggplot2)

gr %>% select(Accession, Treatment) #just displays them
unique(gr$Accession, gr$Treatment) 
length(unique(gr$Accession, gr$Treatment))
length(unique(gr$Accession,order = ascending))
#24 accessions total
length(unique(gr$Species,order = ascending))
#4 species total

as.factor(gr$Accession)
as.factor(gr$Species)
as.factor(gr$Treatment)
as.factor(gr$Rep)
as.factor(gr$LightIntensity)

#COL COUNTS
boxplot(T0_cols~Species,data=gr)
boxplot(T0_cols~Accession,data=gr)
boxplot(T0_cols~Rep,data=gr)
#always 1 col, no normalisation needed

boxplot(T21_cols~Species*Treatment,data=gr)
boxplot(T21_cols~Accession*Treatment,data=gr) 
boxplot(T21_cols~Accession,data=gr) 
boxplot(T21_cols~Rep,data=gr) 
boxplot(T21_cols~Species,data=gr) #spoly makes less colonies at t21 than lemnas

boxplot(T14_cols~Species*Treatment,data=gr) #lminor nore cols at t14 ll than hl
boxplot(T14_cols~Accession*Treatment,data=gr)
boxplot(T14_cols~Accession,data=gr)
boxplot(T14_cols~Species,data=gr)

#RUN THIS IF PUTTING 4 ACCESSIONS INTO CODE
#create new variables for differences during log period
gr %>% mutate()
grmod <- gr %>% group_by(Accession) %>% mutate(Col_RGR = (log(T21_cols-T14_cols)/7))
grmod <- grmod %>% group_by(Accession) %>% mutate(Col_RGRnonlog = ((T21_cols-T14_cols)/7))

#pull out all that have na, normally because t14 higher than t21
#RGR_Na <- grmod[is.na(grmod$Col_RGR),]
#create new df for cols missing in RGR

#grmod$Col_RGR

#write.csv(grmod, "RGR.csv")

#grmod2 <- na.omit(grmod)
grod3 <- grmod[complete.cases(grmod$Col_RGR), ]   
grod3 <- grmod[complete.cases(grmod$Col_RGRnonlog), ]  
#remove all NA cases for RGR column in data frame = 20 rem

range(grod3$Col_RGR)
#between 0 -0.86
range(grod3$Col_RGRnonlog)
#between 0.14 - 60.42 no of cols per day

write.csv(grod3, "Col_countsRGR.csv")


Cols_aov <- aov(Col_RGR ~ Treatment*Accession, data = grod3)
Cols_aov <- aov(ColRGR_nonlog ~ Treatment*Accession, data = grod3)
summary(Cols_aov)

#raw data for all with RGR_col and RGR_nonlog col

#subsetting
#remove rep 1 as sig difference between other LLs
REPLL1 <- grod3 %>% filter(!Rep %in% "LL1")
grod3 <- REPLL1

#work out difference in Col_RGR between LL and HL treats per accession
RGR_new <- grod3 %>% arrange(desc(Col_RGRnonlog))

#expect L. minuta to have more cols per day as less no of fronds in isolates

#SUMMARISING
sum <- grod3 %>% group_by(Accession, Treatment) %>% summarise_all(mean)
#now have summarised RGR rate per access and treatment and other cols
sum$Col_RGR
sum$Col_RGRnonlog

plot(Col_RGR~Accession*Treatment,data=grod3)
plot(Col_RGR~Treatment,data=grod3)
plot(Col_RGR~LightIntensity,data=grod3)
plot(Col_RGR~Accession,data=grod3)
plot(Col_RGRnonlog~Accession*Treatment,data=grod3)
plot(Col_RGRnonlog~Treatment,data=grod3)
plot(Col_RGRnonlog~Accession,data=grod3)
#range of light levels for two treatments
plot(Col_RGRnonlog~Rep,data=grod3)

#detailed summary
Summary <- grod3 %>%
  group_by(Treatment,Accession) %>%
  summarise(RGRnonlog_mean = mean(Col_RGRnonlog), RGRnonlog_stdev = sd(Col_RGRnonlog), RGRnonlog_n= n(), RGRnonlog_maximum = max(Col_RGRnonlog))
Summary

#more detailed summary
Summary1 <- grod3 %>%
  group_by(Treatment,Accession) %>%
  summarise(RGR_mean = mean(Col_RGR), RGR_stdev = sd(Col_RGR), RGR_n= n(), RGR_maximum = max(Col_RGR))
Summary1

Summ <- cbind(Summary, Summary1)

length(Summ)

#to split into hl and ll manuually and read in
cols <- read.csv("Colcounts_RGR_Summary_HLLLsep_nospREM.csv")

#jump down to hl_llmeans


try1 <- split(Summ,cumsum(1:nrow(Summ)%in%25)) #made 2 lists 0 and 1
str(try1)
#call the two lists
try1$`0`
try1$`1`
try1$`0`$Treatment
try1$`1`$Treatment
backtogether <- cbind(try1$`0`, try1$`1`) #need to be the same no args
#puts all hl cols together first then all ll cols after
#need to rename colnames
names(backtogether)

write.csv(backtogether, "Colcounts_RGR_Summary_new.csv")

#find top 10 performing accession RGR values
#6A top of HL, LY01B top of LL
#12 bottom of HL and LL - slowest growing


#look at differences between col nos in treatments
Cols_aov <- aov(Col_RGRnonlog ~ Treatment, data = grod3) 
Cols_aov <- aov(Col_RGRnonlog ~ Accession, data = grod3) 
Cols_aov <- aov(Col_RGRnonlog ~ Treatment*Accession, data = grod3)
Cols_aov <- aov(Col_RGRnonlog ~ Treatment*Species, data = grod3) 
Cols_aov <- aov(Col_RGRnonlog ~ Species, data = grod3) 
Cols_aov <- aov(Col_RGRnonlog ~ Rep, data = grod3)
Cols_aov <- aov(Col_RGRnonlog ~ Treatment*Rep, data = grod3)
summary(Cols_aov) #treatment significant also when do treat*access and treat*species
#rep significant
#accession or species not sig

Cols_aov <- aov(Col_RGR ~ Treatment*Accession, data = grod3)
summary(Cols_aov) 
grod3$Col_RGRnonlog
grod3$Col_RGR

t.test(Col_RGRnonlog ~ Treatment, data = grod3) # not sig
#col counts not affected by HL or LL but lots of variation between reps
t.test(Col_RGR ~ Treatment, data = grod3)

#Tukey works when not means data summary
tuk_out <- TukeyHSD(Cols_aov, "Rep", conf.level=.95)
tuk_out <- TukeyHSD(Cols_aov, "Treatment", conf.level=.95)
tuk_out <- TukeyHSD(Cols_aov, "Accession", conf.level=.95)
tuk_out <- TukeyHSD(Cols_aov, "Species", conf.level=.95)
str(tuk_out)
tuk_out
plot(tuk_out , las=1 , col="brown") #not useful to visualise
#treat significant with or without accession or species included
#rep significant LL1-HL1 LL1-HL2 LL1-HL3 LL1-HL4 LL1-LL2 LL1-LL3 LL1-LL4 LL1-LL5
#accession not sig with treat or not
#species not sig with treat or not
#rep still sig when include treatment - HL all similar, LL1 different from others

#boxplot faceted by treatment
RGRcol_boxplot <- ggplot(grod3,aes(x = reorder(Accession,Col_RGRnonlog, FUN = median),y=Col_RGRnonlog)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 3) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) RGR_col"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Col_RGRnonlog))+
  theme(axis.title.x = element_text(color="#000000", 
                                    size=12, angle=0),
        axis.title.y = element_text(color="#000000",
                                    size=14),
        axis.text.x = element_text(color="#000000", 
                                   size=10, angle = 45, hjust = 1),
        axis.text.y = element_text(color="#000000",
                                   size=14, angle=0),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"))+
  theme(legend.position = "bottom")
RGRcol_boxplot

#quite messy on graph when use log
#change to parameter you want
RGRcol2_boxplot <- ggplot(grod3,aes(x = reorder(Accession,Col_RGR, FUN = median),y=Col_RGR)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 3) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) RGR_col2"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Col_RGR))+
  theme(axis.title.x = element_text(color="#000000", 
                                    size=12, angle=0),
        axis.title.y = element_text(color="#000000",
                                    size=14),
        axis.text.x = element_text(color="#000000", 
                                   size=10, angle = 45, hjust = 1),
        axis.text.y = element_text(color="#000000",
                                   size=14, angle=0),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"))+
  theme(legend.position = "bottom")
RGRcol2_boxplot

#why ll1 so different? producing more cols - lower light? or fresh from site?


#remove LL1 from analysis
cat <- cbind(HL1, HL2)
rep_edit <- (grod3 - Replist$LL1)
#cant do this as unequal numbers

#remove LL1 from analysis manually and read back in
rep_edit <- read.csv("Biomass_norepLL1.csv")

#need to recalculate RGR_cols
rep_edit %>% mutate()
reped3 <- rep_edit %>% group_by(Accession) %>% mutate(Col_RGR = (log(T21_cols-T14_cols)/7))
reped4 <- reped3 %>% group_by(Accession) %>% mutate(Col_RGRnonlog = ((T21_cols-T14_cols)/7))

#remove NAs
reped4_na <- reped4[complete.cases(reped4$Col_RGR), ] 

#repeat boxplots for reps
plot(Col_RGRnonlog~Rep,data=reped4_na)
plot(Col_RGRnonlog~Treatment,data=reped4_na)
plot(Col_RGRnonlog~Accession,data=reped4_na)

#look at differences between col nos now rep removed
Cols_aov <- aov(Col_RGRnonlog ~ Treatment, data = reped4_na) 
Cols_aov <- aov(Col_RGRnonlog ~ Accession, data = reped4_na) 
Cols_aov <- aov(Col_RGRnonlog ~ Treatment*Accession, data = reped4_na)
Cols_aov <- aov(Col_RGRnonlog ~ Treatment*Species, data = reped4_na) 
Cols_aov <- aov(Col_RGRnonlog ~ Species, data = reped4_na) 
Cols_aov <- aov(Col_RGRnonlog ~ Rep, data = reped4_na)
Cols_aov <- aov(Col_RGRnonlog ~ Treatment*Rep, data = reped4_na)
summary(Cols_aov)
#nothing sig

HL <- reped4_na %>% filter(Treatment == "HL") %>% select(Col_RGRnonlog)
HL %>% arrange(Accession) #order all HL RGR values by accession
reped4_na %>% arrange(desc(Col_RGRnonlog)) %>% select(Treatment, Accession) %>% top_n(10)
reped4_na %>% arrange(desc(Col_RGR)) %>% select(Treatment, Accession) %>% top_n(10)

#reflects gr data better in terms of order but treat no longer significant

#how to group RGR_diff into HL or LL but keep non averaged

#look at difference between RGR col means in treatments
cat <- cbind(avg_HL, avg_LL, max_HL, max_LL)
cat <- cat[,-3]
cat <- cat[,-4]
cat <- cat[,-5]
#col 2 = HL mean(Col_RGRnonlog), col 3= LLmean(Col_RGRnonlog).1
names(cat)

avg_HL <- Summary %>% filter(Treatment == "HL") %>% group_by(Accession) %>% summarise(mean(mean))
avg_LL <- Summary %>% filter(Treatment == "LL") %>% group_by(Accession) %>% summarise(mean(mean))
max_HL <- Summary %>% filter(Treatment == "HL") %>% group_by(Accession) %>% summarise(max(maximum))
max_LL <- Summary %>% filter(Treatment == "LL") %>% group_by(Accession) %>% summarise(max(maximum))

colnames(cat)[2]  <- "mean.HL"
colnames(cat)[3]  <- "mean.LL"
colnames(cat)[4]  <- "maximum.HL"
colnames(cat)[5]  <- "maximum.LL"

hl_llmeans <- cat
hl_llmeans <- backtogether

#make summary long version so mean and max per treatment
#called spread from tidyr package
library(tidyr)
#non of these working
#sum_spread <-spread(Summary, Treatment, maximum)
#str(sum_spread)
#sum_sep <- separate(Summary,
#         col = "mean",
#         into = c("Mean HL", "Mean LL"),
#         sep = "Treatment")

raw <- read.csv("Col_countsRGR_nospREM.csv")

RGR_new <- raw

#remove NAs
RGR_new <- RGR_new[complete.cases(RGR_new$Col_RGR), ] 

Summary <- RGR_new %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Col_RGR), stdev = sd(Col_RGR), n= n(), maximum = max(Col_RGR))
Summary

#calc se for ggplot std error
my_sum <- RGR_new %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Col_RGR),
    sd=sd(Col_RGR)
  ) %>%
  mutate( se=sd/sqrt(n))

#USE THIS FOR SUPPLEMENTERY
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25) +
  theme_bw() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90))+
  theme(axis.text.y=element_text(angle = 90))+
  ylab("RGRlog by col gain (cols per day)") +
  xlab("Ecotype")
#facet_wrap(~Treatment) harder to read accession var


#read in new file with grouped HL LL means and max
hl_llmeans <- read.csv("Colcounts_RGR_Summary_HLLLsep.csv")
#this doesnt include the differences
hl_llmeans <- read.csv("Colcounts_RGR_Summary_new_splitHL_LL.csv") #rem sp


#find cols above
hl_llmeans <- cols

#USE 4 OLD ACCESSIONS
hl_llmeans <- Summ #NEEDS TO BE SPLIT TO HL/LL THEN USE

names(hl_llmeans)

#NEWEST VERSION 16/12/2022
library(dplyr)
#take ll - hl
hl_llmeans <- hl_llmeans %>% mutate(ColRGR_diff = (ColRGR_HL_mean - ColRGR_LL_mean)) %>% 
  mutate(Colmax_diff = (ColRGR_HL_maximum - ColRGR_LL_maximum))
hl_llmeans <- hl_llmeans %>% mutate(ColRGRnonlog_diff = (RGRnonlog_HL_mean - RGRnonlog_LL_mean)) %>% 
  mutate(Colmaxnonlog_diff = (RGRnonlog_HL_maximum - RGRnonlog_LL_maximum))
#RGR prop diff
hl_llmeans <- hl_llmeans %>% mutate(ColRGR_propdiff = (ColRGR_diff/ColRGR_LL_mean)*100)
hl_llmeans %>% arrange(desc(ColRGR_propdiff)) %>% select(Accession) %>% top_n(24)
hl_llmeans <- hl_llmeans %>% mutate(ColRGRnonlog_propdiff = (ColRGRnonlog_diff/ColRGR_LL_mean)*100)
hl_llmeans %>% arrange(desc(ColRGRnonlog_propdiff)) %>% select(Accession) %>% top_n(24)

#24 obs
hl_llmeans$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                         "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                         "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                         "HL", "HL", "LL", "LL")

library(dplyr)
library(ggplot2)
#not working to call accessions
#hl_llmeans[which(hl_llmeans$Accession) | hl_llmeans$RGR_diff < 0, ]
hl_llmeans %>% arrange(desc(ColRGR_HL_mean)) %>% select(Accession, ColRGR_HL_mean) %>% top_n(24)
hl_llmeans %>% arrange(desc(ColRGR_LL_mean)) %>% select(Accession, ColRGR_LL_mean) %>% top_n(24)
hl_llmeans %>% arrange(desc(RGRnonlog_HL_mean)) %>% select(Accession, RGRnonlog_HL_mean) %>% top_n(24)
hl_llmeans %>% arrange(desc(RGRnonlog_LL_mean)) %>% select(Accession, RGRnonlog_LL_mean) %>% top_n(24)

hl_llmeans %>% arrange(desc(ColRGR_diff)) %>% select(Accession, ColRGR_diff) %>% top_n(24)
hl_llmeans %>% arrange(desc(ColRGRnonlog_diff)) %>% select(Accession, ColRGRnonlog_diff) %>% top_n(24)
hl_llmeans %>% arrange(desc(ColRGR_propdiff)) %>% select(Accession, ColRGR_propdiff) %>% top_n(24)

names(hl_llmeans)

#effect of treatment growth RGR
#edit this graph to color by light level
par(mfrow = c(2,5)) 
ggplot(hl_llmeans, aes(y=ColRGR_diff, x=reorder(Accession, ColRGR_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90))+
  theme(axis.text.y=element_text(angle = 90))+
  ylab("Difference growth rate (HL - LL col gain per day)") +
  xlab("Ecotype")

write.csv(hl_llmeans, "Col_RGR__oppcalc_propchange_nospREM.csv")

hl_llmeans %>% arrange(desc(ColRGR_HL_mean)) %>% select(Accession, ColRGR_HL_mean) %>% top_n(24)
hl_llmeans %>% arrange(desc(ColRGR_LL_mean)) %>% select(Accession, ColRGR_LL_mean) %>% top_n(24)

#OLD
#take hl away from ll
hl_llmeans <- hl_llmeans %>% mutate(RGR_diff = (RGR_mean1 - RGR_mean)) %>% 
  mutate(max_diff = (RGR_maximum1 - RGR_maximum))
#take ll away from hl
#hl_llmeans <- hl_llmeans %>% mutate(RGR_diff = (mean.HL - mean.LL)) %>% 
#  mutate(max_diff = (maximum.HL - maximum.LL))
#take ll away from hl
#as says online to do % proportion increase
hl_llmeans <- hl_llmeans %>% mutate(RGR_propdiff = (RGR_diff/RGR_mean1)*100)
hl_llmeans %>% arrange(desc(RGR_propdiff)) %>% select(Accession) %>% top_n(28)
#RGR prop diff increase in HL compared to decrease at bottom

write.csv(hl_llmeans, "Col_RGR_propchange.csv")

#plot difference as scatterplot

hl_llmeans <- read.csv("Col_RGR_propchange.csv")
hl_llmeans <- read.csv("Col_RGR_propchange_nospREM.csv")
plot(max_diff ~ RGR_diff, hl_llmeans)
plot(RGR_diff ~ max_diff, hl_llmeans)
#does the maximum growing one also have biggest difference in growth rate?
#RGR_diff	maximum.HL	0.663433863	+	1.27E-06
#max_diff	RGR_diff	0.795616888	+	4.85E-09


range(hl_llmeans$ RGR_maximum)
#non log [1]  0.7142857 24.4285714
# [1] 0.2299197 0.7345234
range(hl_llmeans$ RGR_maximum1)
#[1]  1.00000 17.28571
#[1] 0.2779872 0.6851129

#group data into + or - RGR diff to see which direction performed better
hl_llmeans$RGR_diff > 0

which(hl_llmeans$RGR_diff > 0) #19 accessions positive
#[1]  2  4  5  9 11 16

#grew best at high light NON LOG
#KS03, KS06A, KS06B, KS14, KS16, KS22

which(hl_llmeans$RGR_diff < 0) #9 accessions negative
#[1]  1  3  5  6  8 12 17 26 28
#APP2 KS03 KS06A KS06B KS12 KS17 MOOR1 SEL1
#grew best at low light

range(hl_llmeans$RGR_diff)
#[1] -6.285714  6.535714
#[1] -0.2462871  0.2407263 

hl_llmeans$Accession
hl_llmeans %>% arrange(desc(RGR_diff)) %>% select(Accession) %>% top_n(24)
#shows same thing - which had highest values for difference

#effect of treatment growth col gain
ggplot(hl_llmeans, aes(y=RGR_diff, x=reorder(Accession, RGR_diff))) +
  geom_col(position='dodge', color='black') +
  ylab("RGR diff col gain") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)

hl_llmeans$Accession

#same graph but color by light
#28 OBS
hl_llmeans$Light <- c("HL","LL", "LL", "HL", "LL", "LL", "LL",
                      "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                      "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                      "HL", "HL", "LL", "LL", "HL", "LL", "HL")

#24 obs
hl_llmeans$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                         "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                         "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                         "HL", "HL", "LL", "LL")

par(mfrow = c(2,5)) 
ggplot(hl_llmeans, aes(y=RGR_diff, x=reorder(Accession, RGR_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("RGR difference colony gain") +
  xlab("Ecotype")

#automate correlations
#cant include 4 and 5
lms <- expand.grid(10:13, 10:13)
lms_names <- expand.grid(names(hl_llmeans)[10:13], names(hl_llmeans)[10:13])
out <- vector(mode = "list", length = nrow(lms))
for(i in 1:nrow(lms)){
  lms_col_2 <- lms$Var2[i]
  lms_col_1 <- lms$Var1[i]
  plot_name <- paste0(stringr::str_pad(i, width = 3, pad = "0"), " ", lms_names$Var2[i], " vs. ", lms_names$Var1[i], ".png")
  png(plot_name, width = 500, height = 500, type = "cairo")
  plot(hl_llmeans[, lms_col_1], hl_llmeans[, lms_col_2], main = paste0(lms_names$Var2[i], " vs. ", lms_names$Var1[i]), type = "n")
  text(hl_llmeans[, lms_col_1], hl_llmeans[, lms_col_2], row.names(hl_llmeans), cex = 0.8)
  abline(lm(hl_llmeans[, lms_col_2] ~ hl_llmeans[, lms_col_1]))
  #plots per column and output as png with correspinding names
  dev.off()
  out[[i]] <- data.frame(r.squared = summary(lm(hl_llmeans[, lms_col_2] ~ hl_llmeans[, lms_col_1]))$r.squared,
                         pos_neg = ifelse(summary(lm(hl_llmeans[, lms_col_2] ~ hl_llmeans[, lms_col_1]))$coef[2, 1] > 0, "+", "-"),
                         p.value = summary(lm(hl_llmeans[, lms_col_2] ~ hl_llmeans[, lms_col_1]))$coefficients[2, 4])
}
#select r squared and p val from table and make df output
#use estimate col as gradient to make + or - col

#write to files
library(dplyr)
(all_output <- bind_cols(lms_names, do.call(rbind, out)))
all_output <- all_output[all_output$r.squared != 1, ]
all_output <- all_output[all_output$p.value < 0.05, ]
all_output$r.squared <- ifelse(all_output$pos_neg == "+", all_output$r.squared, -all_output$r.squared)

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Colcounts")
write.csv(all_output, "hl_ll_correlation.csv")


#work out std error using base r by treatment accession combo
aggregate(Col_RGRnonlog ~ Accession + Treatment, data = grod3, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

aggregate(RGR ~ Species + Treatment, data = RGR_new, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- grod3 %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Col_RGRnonlog),
    sd=sd(Col_RGRnonlog)
  ) %>%
  mutate( se=sd/sqrt(n))

ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)


aov <- aov(mean ~ Treatment * Accession, data = my_sum)
aov <- aov(mean ~ Treatment, data = my_sum)
aov <- aov(mean ~ Accession, data = my_sum)
summary(aov)
#neither treatment or accession significant for  col counts
#too much variation in col counts to determine effects
#no point doing tukey

#BIOMASS
#fw variable biomass after 6 weeks
boxplot(FW..mg.~Species*Treatment,data=gr) 
boxplot(FW..mg.~Accession*Treatment,data=gr) 
boxplot(FW..mg.~Accession,data=gr) 
boxplot(FW..mg.~Species,data=gr) #biomass fairly equal all sp, few l.minors make more
boxplot(FW..mg.~Rep,data=gr) #hl1 and ll1 higher than others
boxplot(FW..mg.~Treatment*Rep,data=gr) #hl1 and ll2 higher than other reps
##ll high avg biomass than hl

Biomass_aov <- aov(FW..mg. ~ Treatment, data = gr) 
Biomass_aov <- aov(FW..mg. ~ Accession, data = gr) 
Biomass_aov <- aov(FW..mg. ~ Treatment*Accession, data = gr) 
Biomass_aov <- aov(FW..mg. ~ Treatment*Species, data = gr) 
Biomass_aov <- aov(FW..mg. ~ Species, data = gr) 
Biomass_aov <- aov(FW..mg. ~ Rep, data = gr)
summary(Biomass_aov)#treatment sig, rep very significant
#species and accession not sig. species sig when add treatment term

#Tukey works when not means data summary
tuk_out <- TukeyHSD(Biomass_aov, "Rep", conf.level=.95)
tuk_out <- TukeyHSD(Biomass_aov, "Treatment", conf.level=.95)
tuk_out <- TukeyHSD(Biomass_aov, "Accession", conf.level=.95)
tuk_out <- TukeyHSD(Biomass_aov, "Species", conf.level=.95)
str(tuk_out)
tuk_out
plot(tuk_out , las=1 , col="brown") #not useful to visualise
#diffs biomass between treatments with sig p values
#diffs biomass between reps with sig p values

#diffs between species but not sig p values
#diffs between accessions but not sig p values

#cant see any accession or species differences due to variation in repeats/treats

#subsetting
#remove rep 1 and 2 pre 2021
REPLL1 <- gr %>% filter(!Rep %in% "LL1")
REPHL1 <- REPLL1 %>% filter(!Rep %in% "HL1")
# REPHL1 has LL1 and HL1 reps missing

boxplot(FW..mg.~Rep,data=REPHL1) 

#remove t5, t6 and t9 cols NOT RELEVANT AS NOT INCLUDED
#gr1 <- subset(gr, select=-c(T5.area, T6.area, T9.area))

Summary <- reped4_na %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Col_RGRnonlog), stdev = sd(Col_RGRnonlog), n= n(), maximum = max(Col_RGRnonlog))
Summary

#shows same data but messier
Summary <- grod3 %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Col_RGR), stdev = sd(Col_RGR), n= n(), maximum = max(Col_RGR))
Summary

ggplot(Summary, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), position = position_dodge(0.9), width = 0.25)

library(forcats)

as.factor(Summary$Accesssion)
as.factor(Summary$Treatment)

class(Summary$mean)

#try to seperate data by HL and LL to make graphs
p1 <- Summary %>% filter(!Treatment %in% "LL") #HL
p2 <- Summary %>% filter(Treatment %in% "LL") #LL

p1 <- p1 %>% arrange(desc(mean)) 

# Reorder following the value of another column:
Summary %>% mutate(mean = fct_reorder(mean, Treatment))
Summ <- Summary %>% filter(Treatment) %>% group_by(Accession) %>% arrange(desc(mean))
Summary %>% select(mean, Accession, Treatment) %>% filter(Treatment) %>% arrange(desc(mean)) 


#bar graph for HL
ggplot(p1, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), position = position_dodge(0.9), width = 0.25)

#bar graph for LL
ggplot(p2, aes(y=mean, x=Accession)) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), position = position_dodge(0.9), width = 0.25)

#hl ll treatment split graphs
library(stringr)

RGR_lag$Accession <- as.factor(RGR_lag$Accession)
RGR_lag$Treatment <- as.factor(RGR_lag$Treatment)

#remove if water
#param$Accession == "Water"

##labelling you levels  You can put in your light here :) 

RGR_lag <- grod3
levels(RGR_lag$Treatment)
levels(RGR_lag$Treatment) <- c("HL", "LL")  

# numbers of reps not working
library(dplyr)
RGR_lag %>% group_by(Treatment, Accession) %>% dplyr::summarize(n = n()) %>% print(n = Inf)
RGR_lag$Treatment

names(RGR_lag)

#need to tidy up variance using graphs + growth rate data
#correlate to gr data

#t tests and anova summary stats
Sum2 <- t(Summary)
Sum2 <- c-(colnames(Summary$Accession))
aov <- aov(mean ~ Treatment * Accession, data = Summary) #same result
aov <- aov(mean ~ Accession * Treatment, data = Summary) #same result
aov <- aov(mean ~ Treatment, data = Summary) #sig when just treatment used
t.test(mean ~ Treatment, data = Summary) # t test difference 
# cant use t test for accession as more than 2 
summary(aov) #accession and access*treat interaction
tuk_out <- TukeyHSD(aov, conf.level=.95)
(tuk_out)
#interaction most significant, need to look at each accession, treatment pair individually


#summarise by mean and sd
Summary_lag <- RGR_lag %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(RGR7), stdev = sd(RGR), n= n(), maximum = max(RGR7))
Summary_lag

#save to csv
write.csv(Summary, file="Summary_RGR.csv")
write.csv(Summary_lag, file="Summary_RGR.csv")

##labelling you levels  You can put in your light here :) 

levels(RGR_lag$Treatment)
levels(RGR_lag$Treatment) <- c("HL", "LL")  

# numbers of reps not working
library(dplyr)
RGR_lag %>% group_by(Treatment, Accession) %>% dplyr::summarize(n = n()) %>% print(n = Inf)
RGR_lag$Treatment

names(RGR_lag)


#change to parameter you want
RGR7_boxplot <- ggplot(RGR_lag,aes(x = reorder(Accession,RGR7, FUN = median),y=RGR7)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 3) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) RGR_lag"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(RGR7))+
  theme(axis.title.x = element_text(color="#000000", 
                                    size=12, angle=0),
        axis.title.y = element_text(color="#000000",
                                    size=14),
        axis.text.x = element_text(color="#000000", 
                                   size=10, angle = 45, hjust = 1),
        axis.text.y = element_text(color="#000000",
                                   size=14, angle=0),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"))+
  theme(legend.position = "bottom")
RGR7_boxplot

boxplot(RGR ~ Rep, data = RGR_lag)
boxplot(RGR ~ Rep * Accession, data = RGR_lag)

levels(RGR_lag$Rep)
levels(RGR_lag$Rep) <- c("LL1", "LL2", "LL3", "LL4", "LL5", "HL1", "HL2", "HL3", "HL4",
                         "SLL1", "SLL2", "SLL3", "SHL1", "SHL2", "SHL3")  

#remove last six reps
#RGR_lag$Rep == "SLL1"#, "SLL2", "SLL3", "SHL1", "SHL2", "SHL3"

#see how RGR changes between reps
RGRreps_boxplot <- ggplot(RGR_lag,aes(x = reorder(Accession,RGR, FUN = median),y=RGR)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Rep, scales = "free_x", ncol = 3) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) RGR_log"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(RGR))+
  theme(axis.title.x = element_text(color="#000000", 
                                    size=12, angle=0),
        axis.title.y = element_text(color="#000000",
                                    size=14),
        axis.text.x = element_text(color="#000000", 
                                   size=10, angle = 45, hjust = 1),
        axis.text.y = element_text(color="#000000",
                                   size=14, angle=0),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"))+
  theme(legend.position = "bottom")
RGRreps_boxplot

#spearmans rank
access.ranked <- data.frame(cbind(rank(RGR_lag$RGR, ties.method = 'average'),
                                  rank(RGR_lag$Rep, ties.method = 'average')))
colnames(access.ranked) <- c('RGR', 'Rep')

corr <- cor.test(x=RGR_lag$RGR, y=RGR_lag$Rep, method = 'spearman')
corr

#chi sq
#library("gplots")
# 1. convert the data as a table
#RGRchi <- as.table(as.matrix(RGR_lag))
#balloonplot(t(RGRchi), main ="housetasks", xlab ="", ylab="",
#            label = FALSE, show.margins = FALSE)
#chisq <- chisq.test(RGR)
#chisq

chisq.test(RGR_lag$Treatment, RGR_lag$RGR, correct=FALSE)
table(RGR_lag$Treatment, RGR_lag$RGR)
chisq.test(RGR_lag$Treatment, RGR_lag$RGR)
table(RGR_lag$Rep, RGR_lag$RGR)
chisq.test(RGR_lag$Rep, RGR_lag$RGR)
table(RGR_lag$Accession, RGR_lag$RGR)
chisq.test(RGR_lag$Accession, RGR_lag$RGR)

#wilcox test
#needs 2 levels to test difference
wilcox.test()
wilcox.test(RGR_lag$RGR~RGR_lag$Treatment)

library(reshape2)
melted <- melt(RGR_lag)
trans <- t(RGR_lag)
class(trans)
colnames(trans)
RGR_lag$Accession
trans <- colnames(RGR_lag$Accession)

#remove 221-240 when arranged by RGR
#RGR_newd <- RGR_new[is.na(RGR_new),] #removes NA values
#remove NA values
#RGR_sumn <- RGR_sum[!is.na(RGR_sum$RGR),] #removes NA values
#remove rows 41-56 all NAs

#cant summarise until removed NAs as doesnt calc properly
#removed rows manually on excel and re-read in file above _noNArow
#summarise all cols by accession and treatment
library(dplyr)
#done on data without LL?
RGR_sum <- RGR_new %>% group_by(Treatment, Accession) %>% summarise_all(mean, na.rm = TRUE)
RGR_sum3 <- reped4_na %>% group_by(Treatment, Accession) %>% summarise_all(mean, na.rm = TRUE)
#same data

#summarised plot
barplot(RGR_sum$RGR, 
        names.arg = RGR_sum$Accession,
        horiz = T, las = 1,
        xlim = c(0, 1.5))
box()

#summarised plot
ggplot(RGR_sum, aes(x=Accession, y=RGR, fill=Treatment, group=Treatment)) +
  geom_col(position='dodge', color='black')

RGR_sumsp <- RGR_new %>% group_by(Treatment, Species) %>% summarise_all(mean) %>% arrange(desc(RGR)) #add desc as argument to do descending

ggplot(RGR_sumsp, aes(x=Species, y=RGR, fill=Treatment, group=Treatment)) +
  geom_col(position='dodge', color='black')

#get RGR for HL and LL and find the difference between them
# not working
HL_RGR <- RGR_sum %>% filter(Treatment == "HL") %>% select(Accession, RGR)
LL_RGR <- RGR_sum %>% filter(Treatment == "LL") %>% select(Accession, RGR)
comb <- rbind(HL_RGR, LL_RGR)
comb %>% group_by(Accession, Treatment) %>% mutate(RGR_diff = (RGR))
diff <- (LL_RGR$RGR - HL_RGR$RGR)


#split by HL LL to produce bar plots

#split by HL LL to produce line plot overtime
#produce individual plots per rep, summarise as average
#get rid of non numeric columns
#transpose so accessions are col names
#change row names to just numbers of days

#ADAPTED FOR COL NUMBERS
accessions <- (RGR_sum[,2])
#RUN THIS IN ORDER AND HASH OUT RELEVANT VERSIONS
accessions <- (HL[,2])
#accessions <- (LL[,2])

names(grod3)
names(reped4_na)

#percencov
HL <- RGR_sum %>% filter(Treatment == "HL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)
#LL <- RGR_sum %>% filter(Treatment == "LL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)


#to see hl and ll seperately need to split by treatment here
#HL <- grod3 %>% filter(Treatment == "HL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)
LL <- grod3 %>% filter(Treatment == "LL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)

#without rep LL1
#HL <- reped4_na %>% filter(Treatment == "HL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)
LL <- reped4_na %>% filter(Treatment == "LL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)
#HL$Treatment <- NULL
LL$Treatment <- NULL

#sub <- (RGR_sum[,17:22]) #just t0 to t42 columns to show all data
#to see hl and ll seperately need to split by treatment 
#sa_comb <- sub
sa_comb <- HL
#sa_comb <- LL
names(sa_comb)

sa_comb$Accession <- NULL
sa_comb$Treatment <- NULL

#row.names(sa_comb) <- sa_comb$Accession
#sa_comb$Accession <- NULL
#take col names from first col and give new names
sa_comb <- t(sa_comb) #become a matrix during transposition
sa_comb
#turn back into dataframe again
sadf <- data.frame(sa_comb)
class(sadf)

row.names(sadf)

#row names defined col
#null col 

time_vec <- numeric(length = nrow(sadf)) #make vector same name as row names

library(stringr)

#replace characters in a string
row.names(sadf) <- str_replace_all(row.names(sadf), "T0_", "0")
row.names(sadf) <- str_replace_all(row.names(sadf), "_cols", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "T", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "cols", "")
#looks for pattern and give replacement

row.names(sadf)

time_vec <- row.names(sadf)
time_vec <- as.numeric(time_vec) #coerce to numeric

class(sadf) 
accessions <- t(accessions)
colnames(sadf) <- (accessions)

sadf$time <- time_vec

str(sadf)

class(sadf)



library(tidyr)

sadf$KS02
sadf$time

#par(mfrow = c(1,1)) #no of rows and cols in plot display
#plot(KS02 ~ time, sadf, type = "n", xlab = "Time (days)", ylab = "RGR col gain T14-T21/7 (gain per day)")
#lines(KS03 ~ time, sadf) # one at a time
#lines(SEL1 ~ time, sadf) 
#lines(KS04 ~ time, sadf) # one at a time
#lines(LY03 ~ time, sadf) 

#do something to every i
plot(KS02 ~ time, sadf, type = "n", 
     main = "Growth rate in HL",
     xlab = "Time (days)", 
     ylab = expression(paste("Col gain per day (gain per day)")), 
     axes = F,
     xlim = c(0, 22),
     ylim = c(0, 100))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:24){
  lines(sadf$time, sadf[, i], col = line_colours[i], lty = sort(line_type)[i], lwd = 2)
  
}
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
legend("bottomright", pt.cex = 150, cex = 0.3, col = line_colours, lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()

line_colours <- rep(1:24, length.out = ncol(sadf)-1) #8 colors
line_colours
line_type <- rep(1:2, length.out =  ncol(sadf)-1)
line_type

#redo plot for hl with 1 color
plot(KS02 ~ time, sadf, type = "n", 
     main = "Growth rate in HL",
     xlab = "Time (days)", 
     ylab = "Colony gain per day", 
     axes = F,
     xlim = c(0, 21),
     ylim = c(0, 120))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:24){
  lines(sadf$time, sadf[, i], col = "#000000", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 2], col = "#FF0000", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 11], col = "#006400", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 16], col = "#006400", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 7], col = "#81007F", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 5], col = "#0000FF", lty = 1, lwd = 2)
}
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
#legend("bottomright", pt.cex = 150, cex = 0.3, col = "#00000033", lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()


#redo plot for ll with 1 color
plot(KS02 ~ time, sadf, type = "n", 
     main = "Growth rate in LL",
     xlab = "Time (days)", 
     ylab = "Colony gain per day", 
     axes = F,
     xlim = c(0, 21),
     ylim = c(0, 120))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:24){
  lines(sadf$time, sadf[, i], col = "#000000", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 2], col = "#FF0000", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 11], col = "#006400", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 16], col = "#006400", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 7], col = "#81007F", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 5], col = "#0000FF", lty = 1, lwd = 2)
}
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
#legend("bottomright", pt.cex = 150, cex = 0.3, col = "#00000033", lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()
