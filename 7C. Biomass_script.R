#BIOMASS SCRIPT
#script to explore hl ll dw growth rate data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Biomass.csv")
gr <- read.csv("Biomass - start weights.csv")

tail(gr)
head(gr)
names(gr)
#names
#[1] "Accession"  "T0_cols"    "T7_cols"    "T14_cols"   "T21_cols"   "FW..mg."   
#[7] "FDW..mg."   "Start_date" "Rep"        "Treatment"  "Species"    "X"         
#[13] "X.1"        "X.2"  

#check classes of all cols
sapply(gr, class)

gr<- subset(gr, select=-c(T0_cols, T7_cols, T14_cols, T21_cols))

library(dplyr)
library(stringr)
library(ggplot2)

gr %>% select(Accession, Treatment) #just displays them
unique(gr$Accession, gr$Treatment) 
length(unique(gr$Accession, gr$Treatment))
length(unique(gr$Accession,order = ascending))
#25 accessions total
length(unique(gr$Species,order = ascending))
#5 species total

as.factor(gr$Accession)
as.factor(gr$Species)
as.factor(gr$Treatment)
as.factor(gr$Rep)
as.factor(gr$LightIntensity)


#BIOMASS
#fw variable biomass after 6 weeks
boxplot(FW..mg.~Treatment,data=gr) #sig
boxplot(FW..mg.~Species*Treatment,data=gr) #only sig when data not normalised to start
#biomass, l. minor higher mean ll and more outliers, l. turion higher mean hl all can say
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
Biomass_aov <- aov(FW..mg. ~ Treatment*Rep, data = gr)
summary(Biomass_aov)#treatment sig, rep very significant
#species and accession not sig. species sig when add treatment term

t.test(FW..mg. ~ Treatment, data = gr) # very sig
#data:  FW..mg. by Treatment
#t = -7.4042, df = 138.99, p-value = 1.152e-11
#more growth in LL
#when remove ll1 and hl2 species sig - rep still sig

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

#edit - when remove ll1 hl2, species treat - L.turion / L.minor sig diff,
#L.minuta-L.minor sig diff
#species still sig when normalise start values of starting weight
#L. minuta to L. minor, L. turion to L. minor
#treat still significant
#treat accession visual diff on graph but only treat sig

#subsetting
#reps 1 and 2 sig different from other like repeats in treatment
#remove rep 1 and 2 pre 2021
REPLL1 <- gr %>% filter(!Rep %in% "LL1")
REPHL1 <- REPLL1 %>% filter(!Rep %in% "HL1")
# REPHL1 has LL1 and HL1 reps missing

boxplot(FW..mg.~Rep,data=REPHL1) 

#replace to remove reps LL1 and LL2
gr <- REPHL1

#to see if species affect is due to diff start weights, need to include
#as a variable and do mutate to take away from end weight

gr <- read.csv("Biomass - start weights.csv")

names(gr)

#measured start weights once, replicated for all rows
#could improve by getting exact 3 fr and measuring more times

#subsetting
gr<- subset(gr, select=-c(T0_cols, T7_cols, T14_cols, T21_cols))
#reps 1 and 2 sig different from other like repeats in treatment
#remove rep 1 and 2 pre 2021
REPLL1 <- gr %>% filter(!Rep %in% "LL1")
REPHL1 <- REPLL1 %>% filter(!Rep %in% "HL1")
# REPHL1 has LL1 and HL1 reps missing

gr <- REPHL1

#create new variable for normalised start biomass value
gr %>% mutate()
grmod <- gr %>% group_by(Accession) %>% mutate(FW_mod = (FW..mg.-Start.weights..mg.))
#grmod <- grmod %>% group_by(Accession) %>% mutate(FDW = (FDW..mg./Start.weights..mg.))

#grmod2 <- na.omit(grmod)
grod3 <- grmod[complete.cases(grmod$FW_mod), ]   
#grod4 <- grmod[complete.cases(grmod$FDW..mg.), ]  
#remove all NA cases for RGR column in data frame = 20 rem

#to replace old column and re run boxplots and aovs
#gr$FW..mg. <- grmod$FW_mod

range(grod3$FW_mod)
#between -15.9 5213.7
range(grod3$FDW..mg.)
#between 0 785.9

write.csv(grod3, "Biomass_minusstart.csv")

#work out difference in RGR_mod between LL and HL treats per accession
RGR_new <- grod3 %>% arrange(desc(FW_mod))

#expect S.poly to have more biomass as weigh more to start with
#even when normalised 

#SUMMARISING
sum <- grod3 %>% group_by(Accession, Treatment) %>% summarise_all(mean)
#now have summarised RGR rate per access and treatment and other cols
sum$FW_mod

#range of light levels for two treatments
plot(FW_mod~LightIntensity,data=grod3) #good spread

#all RGR values for all independent observations in HL
#only 2 reps HL as removed first
HL <- grod3 %>% filter(Treatment == "HL") %>% select(FW_mod)
HL %>% arrange(Accession) #order all HL RGR values by accession
grod3 %>% arrange(desc(FW_mod)) %>% select(Treatment, Accession) %>% top_n(10)

LL <- grod3 %>% filter(Treatment == "LL") %>% select(FW_mod)
LL %>% arrange(Accession) #order all HL RGR values by accession
grod3 %>% arrange(desc(FW_mod)) %>% select(Treatment, Accession) %>% top_n(10)

#Summarised AVERAGED RGR values for all in HL and LL
avg_HL <- grod3 %>% filter(Treatment == "HL") %>% group_by(Accession) %>% summarise(mean(FW_mod))
avg_LL <- grod3 %>% filter(Treatment == "LL") %>% group_by(Accession) %>% summarise(mean(FW_mod))

#more detailed summary
Summary <- grod3 %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(FW_mod), stdev = sd(FW_mod), n= n(), maximum = max(FW_mod))
Summary

Summary1 <- grod3 %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(FDW..mg.), stdev = sd(FDW..mg.), n= n(), maximum = max(FDW..mg.))
Summary1

Summary <- cbind(Summary, Summary1)

write.csv(Summary, "Biomass_RGR_Summary.csv")

#find top 10 performing accession RGR values
#6A top of HL, LY01B top of LL
#12 bottom of HL and LL - slowest growing

grod3 %>% 
  ggplot(aes(FW_mod, Accession, col = Treatment)) + 
  geom_point(alpha = 0.8) +
  facet_wrap(~Treatment)

#boxplot faceted by treatment
RGRbiom_boxplot <- ggplot(grod3,aes(x = reorder(Accession, FW_mod, FUN = median),y=FW_mod)) +
  
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
  labs(title =expression("(a) FW"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(FW_mod))+
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
RGRbiom_boxplot 

#boxplot faceted by treatment
RGRFDW_boxplot <- ggplot(grod3,aes(x = reorder(Accession, FDW..mg., FUN = median),y=FDW..mg.)) +
  
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
  labs(title =expression("(a) FDW"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(FDW..mg.))+
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
RGRFDW_boxplot 


#by species
RGRbiomsp_boxplot <- ggplot(grod3,aes(x = reorder(Species,FW_mod, FUN = median),y=FW_mod)) +
  
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
  labs(title =expression("(a) FW"),subtitle = ("Sort by Species"),x=expression(),
       y=expression(FW_mod))+
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
RGRbiomsp_boxplot

library(tidyr)

#not run so far

#read in new file with grouped HL LL means and max
hl_llmeans <- read.csv("Biomass_RGR_Summary_HLLLsep.csv")

names(hl_llmeans)

#changed calc
#hll-ll
hl_llmeans <- hl_llmeans %>% mutate(FW_mean_diff = (FW.mean.HL - FW.mean.LL)) %>% 
  mutate(FW_max_diff = (FW.maximum.HL - FW.maximum.LL))
hl_llmeans <- hl_llmeans %>% mutate(FDW_mean_diff = (FDW.mean.HL - FDW.mean.LL)) %>% 
  mutate(FDW_max_diff = (FDW.maximum.HL - FDW.maximum.LL))

#hl_llmeans <- hl_llmeans %>% mutate(FW_mean_diff = (FW.mean.LL - FW.mean.HL)) %>% 
#  mutate(FW_max_diff = (FW.maximum.LL - FW.maximum.HL))
#hl_llmeans <- hl_llmeans %>% mutate(FDW_mean_diff = (FDW.mean.LL - FDW.mean.HL)) %>% 
#  mutate(FDW_max_diff = (FDW.maximum.LL - FDW.maximum.HL))


#plot difference as scatterplot
plot(FW_mean_diff ~ FDW_mean_diff, hl_llmeans)
plot(FW_max_diff ~ FDW_max_diff, hl_llmeans)
#any correlation between FDW and FW?

plot(FW_mean_diff ~ FW.mean.LL, hl_llmeans) #looks positive
plot(FW_mean_diff ~ FW.mean.HL, hl_llmeans) #looks negative
plot(FW_mean_diff ~ FDW_mean_diff, hl_llmeans)

plot(FDW_mean_diff ~ FDW.mean.LL, hl_llmeans) #looks positive
plot(FDW_mean_diff ~ FDW.mean.HL, hl_llmeans) #looks negative


range(hl_llmeans$FW_mean_diff)
#[1]  -860.200 2119.083
range(hl_llmeans$FW_max_diff)
#[1]  -770.5 4627.3
range(hl_llmeans$FDW_mean_diff)
#[1]  -55.4675 135.7667
range(hl_llmeans$FDW_max_diff)
#[1]  -52.1 687.8

#group data into + or - RGR diff to see which direction performed better
hl_llmeans$FW_mean_diff > 0

which(hl_llmeans$FW_mean_diff > 0) #25 accessions positive
which(hl_llmeans$FW_mean_diff < 0) #3 accessions negative
#more biomass at highlight
# 1 7 27 
#APP2 KS09 NUFF1

hl_llmeans$FDW_mean_diff > 0

which(hl_llmeans$FDW_mean_diff > 0) #21 accessions positive
which(hl_llmeans$FDW_mean_diff < 0) #7 accessions negative
# 1 4 11 17 24 27 28
#APP2 KS04 KS15 KS22 LY02 NUFF1

hl_llmeans$Accession

#as says online to do % proportion increase
hl_llmeans <- hl_llmeans %>% mutate(FW_propdiff = (FW_mean_diff/FW.mean.LL)*100)
hl_llmeans %>% arrange(desc(FW_propdiff)) %>% select(Accession) %>% top_n(28)
#TOP = DID BETTER IN HL
hl_llmeans <- hl_llmeans %>% mutate(FDW_propdiff = (FDW_mean_diff/FDW.mean.LL)*100)
hl_llmeans %>% arrange(desc(FDW_propdiff)) %>% select(Accession) %>% top_n(28)
#TOP = DID BETTER IN HL

write.csv(hl_llmeans, "Biomass_RGR_propchange.csv")
write.csv(hl_llmeans, "Biomass_RGR_propchange_oppcalc.csv")#NEW
names(hl_llmeans)

ggplot(hl_llmeans, aes(y=FW_mean_diff, x=reorder(Accession, FW_mean_diff))) +
  geom_col(position='dodge', color='black') +
  ylab("Fresh weight difference (mg)") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

#BIOMASS FW
ggplot(hl_llmeans, aes(y=FW_mean_diff, x=reorder(Accession, FW_mean_diff))) +
  geom_col(position='dodge', color='black') +
  ylab("Fresh weight difference (mg)") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=FDW_mean_diff, x=reorder(Accession, FDW_mean_diff))) +
  geom_col(position='dodge', color='black') +
  ylab("Freeze dry weight difference (mg)") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

Biomass_aov <- aov(FW_mod ~ Treatment*Accession, data = grod3) 
Biomass_aov <- aov(FDW..mg. ~ Treatment*Accession, data = grod3) 
summary(Biomass_aov)

#work with relative water content
rwc <- read.csv("Biomass_minusstart+RWC_remsp.csv")
names(rwc)
library(dplyr)
#cut any values >100 and below 0 so not included in averages
which(rwc$RWC > 0) #keep all above 0, just 37 fails
which(rwc$RWC < 100) #keep all below, fails 64, 73, 132

rwc1 <- rwc %>% filter(RWC > 0)
rwc2 <- rwc1 %>% filter(RWC < 100)
rwc <- rwc2

#check
which(rwc$RWC < 0) #0
which(rwc$RWC > 100) #0

(rwc$Treatment)
length(rwc$Treatment)#just matches number of obs

as.factor(rwc$Treatment) #2 LEVELS

levels(rwc$Treatment) #2 LEVELS
levels(rwc$Treatment) <- c("LL", "HL")
which(rwc$Treatment == "") #0

#measure affect of species, treatment, accession
rwc_aov <- aov(RWC ~ Treatment*Accession, data = rwc) #all sig, esp treat
rwc_aov <- aov(RWC ~ Treatment*Species, data = rwc) #treat, sp sig
summary(rwc_aov)

tuk_out <- TukeyHSD(rwc_aov, "Treatment", conf.level=.95) #7% diff ll - hl
tuk_out <- TukeyHSD(rwc_aov, "Accession", conf.level=.95) #no p adj <0.05
tuk_out <- TukeyHSD(rwc_aov, "Species", conf.level=.95) #L. minuta and L. turion
#both less than L. minor? -13.8 Lt Lm <0.00001. Lmu Lmi -5.4 0.02.
str(tuk_out)
tuk_out

boxplot(RWC~Species*Treatment,data=rwc)
boxplot(RWC~Accession*Treatment,data=rwc)
boxplot(RWC~Treatment,data=rwc)

#higher water content in LL
#L TURION lower in ll
#L minor higher in hl

#summarise data into accessions and treatment
Summary <- rwc %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(RWC), stdev = sd(RWC), n= n(), maximum = max(RWC))
Summary

#check summary length
length(Summary) #variable no
#obs= 48

try1 <- split(Summary,cumsum(1:nrow(Summary)%in%25)) #made 2 lists 0 and 1
str(try1)
#call the two lists
try1$`0`
try1$`1`
backtogether <- cbind(try1$`0`, try1$`1`) 

#change names of columns
colnames(backtogether)
colnames(backtogether) <- c("Treatment", "Accession", "mean.HL_RWC", 
                            "stdev", "n", "max.HL_RWC", "Treatment1",
                            "Accession1", "mean.LL_RWC", 
                            "stdev1", "n1", "max.LL_RWC")

#add column for light env
backtogether$EnvLight <- c("dLL", "dLL", "dHL", "dLL", "dLL", "dLL",
                  "dHL", "dHL", "dHL", "dLL", "dHL", "dLL", "dLL",
                  "dHL", "dHL", "dLL", "dLL", "dLL", "dLL", "dHL", "dHL",
                  "dHL", "dLL", "dLL")

#see if rwc varies between sites came from
boxplot(mean.HL_RWC~EnvLight,data=backtogether) #hl
boxplot(mean.LL_RWC~EnvLight,data=backtogether) #ll

rwc_aov <- aov(mean.HL_RWC ~ EnvLight, data = backtogether) #not quite sig
rwc_aov <- aov(mean.LL_RWC ~ EnvLight, data = backtogether) #ns 0.745
summary(rwc_aov)

t.test(mean.HL_RWC ~ EnvLight, data = backtogether) # not sig 0.09 but dll higher 79, 87
t.test(mean.LL_RWC ~ EnvLight, data = backtogether) #not sig, similar 91

#20 14 17 9 least water content in HL
#6b 16 25 LL