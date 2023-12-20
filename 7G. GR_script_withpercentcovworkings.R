#script to explore hl ll dw growth rate data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("GR_Oct2022.csv")
gr <- read.csv("RGR_norm.csv") #normalised

tail(gr)
head(gr)
names(gr)
#names
# [1] "Accession"      "T0.area"        "T5.area"        "T6.area"       
#[5] "T7.area"        "T9.area"        "T14.area"       "T21.area"      
#[9] "T28.area"       "T42.area"       "Start_date"     "Rep"           
#[13] "Treatment"      "Species"        "LightIntensity" "SpecIntensity" 

#check classes of all cols
sapply(gr, class)
#make numeric
gr$T0.area <- as.numeric(as.character(gr$T0.area))
class(gr$T0.area)

library(dplyr)
library(stringr)
library(ggplot2)

gr %>% select(Accession, Treatment) #just displays them
unique(gr$Accession, gr$Treatment) 
length(unique(gr$Accession, gr$Treatment))
length(unique(gr$Accession,order = ascending))
#28 accessions total
length(unique(gr$Species,order = ascending))
#5 species total

#see how variable RGRlog is for summer 2022 accessions
#response to treatment compared to others
bp <- ggplot(gr, aes(x=RGRlog, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

as.factor(gr$Accession)
as.factor(gr$Species)
as.factor(gr$Treatment)
as.factor(gr$LightIntensity)

#gr before normalisation, 
boxplot(T0.area~Species,data=gr) #varies per species
boxplot(T0.area~Accession,data=gr) #varies per accession
#need to normalise all values to t0? how make it zero
boxplot(T0.area~Rep,data=gr) #varies per accession
#hl most consistant, but smaller areas than ll

boxplot(T42.area~Species*Treatment,data=gr) #varies not all within sp finished by t 42
boxplot(T42.area~Accession*Treatment,data=gr) #varies per accession
boxplot(T42.area~Accession,data=gr) #ks03 always finish, ks16 ly01b never
boxplot(T42.area~Species,data=gr)

gr %>% mutate()
grmod <- gr %>% group_by(Accession) %>% mutate(RGR = (log(T21.area-T14.area)/7))
#gr %>% group_by(Accession, Treatment) %>% mutate(PRI = (log(T21.area-T14.area)/7))

#pull out all that have na, normally because t14 higher than t21
RGR_Na <- grmod[is.na(grmod$RGR),]
#create new df for cols missing in RGR

grmod$RGR

#write.csv(grmod, "RGR.csv")

#grmod2 <- na.omit(grmod)
grod3 <- grmod[complete.cases(grmod$RGR), ]    
#remove all NA cases for RGR column in data frame = 20 rem

range(grod3$RGR)
#how get a minus RGR?

#subsetting
#remove rep 1 and 2 pre 2021
REPLL1 <- gr %>% filter(!Rep %in% "LL1")
REPLL2 <- REPLL1 %>% filter(!Rep %in% "LL2")
gr <- REPLL2
#remove t5, t6 and t9 cols
gr1 <- subset(gr, select=-c(T5.area, T6.area, T9.area))

#SUMMARISING
sum <- grod3 %>% group_by(Accession, Treatment) %>% summarise_all(mean)
#now have summarised RGR rate per access and treatment and other cols
sum$RGR

plot(RGR~Accession*Treatment,data=grod3)
plot(RGR~Treatment,data=grod3)
plot(RGR~Accession,data=grod3)
plot(RGR~LightIntensity,data=grod3)
boxplot(LightIntensity~Treatment,data=grod3)
boxplot(LightIntensity~Rep,data=grod3)
#range of light levels for two treatments
plot(RGR~Rep,data=grod3)


HL <- grod3 %>% filter(Treatment == "HL") %>% select(RGR)
HL %>% arrange(Accession) #order all HL RGR values by accession
grod3 %>% arrange(desc(RGR)) %>% select(Treatment, Accession) %>% top_n(10)
#find top 10 performing accession RGR values

avg_HL <- grod3 %>% filter(Treatment == "HL") %>% group_by(Accession) %>% summarise(mean(RGR))
avg_LL <- grod3 %>% filter(Treatment == "LL") %>% group_by(Accession) %>% summarise(mean(RGR))

grod3 %>% 
  ggplot(aes(RGR, Accession, col = Treatment)) + 
  geom_point(alpha = 0.8) +
  facet_wrap(~Treatment)

#affect on RGR
RGR_aov <- aov(RGR ~ Treatment, data = grod3) 
#not sig treatment if remove ll1 and ll2
RGR_aov <- aov(RGR ~ Accession, data = grod3) 
#accession sig with and without ll1 nd ll2
RGR_aov <- aov(RGR ~ Rep, data = grod3) 
#rep sig when all included and also when ll1 and ll2 not included
RGR_aov <- aov(RGR ~ LightIntensity, data = grod3) 
#light intensity not sig on RGR without ll1 ll2 = no treatment affect
RGR_aov <- aov(RGR ~ LightIntensity, data = grod3)
RGR_aov <- aov(RGR ~ Accession * Treatment, data = grod3)
RGR_aov <- aov(RGR ~ Accession * Treatment * Rep, data = grod3) #nothing sig, all 3 terms
str(RGR_aov)
summary(RGR_aov) #accession and treatment sig

T0_aov <- aov(T0.area ~ Treatment, data = gr) #not sig due to treat no ll1 ll2
T0_aov <- aov(T0.area ~ Species, data = gr) #sig area due to species no ll1 ll2
T0_aov <- aov(T0.area ~ Accession, data = gr) #sig area due to accession no ll1 ll2
T0_aov <- aov(T0.area ~ Rep, data = gr) #sig area due to rep with no ll1 ll2
T0_aov <- aov(T0.area ~ Treatment * Accession, data = gr) #no affect treat no ll1 ll2
T0_aov <- aov(T0.area ~ Treatment * Species, data = gr) #no affect treat no ll1 ll2
str(T0_aov)
summary(T0_aov) 

T0_aov <- aov(T0.area ~ Treatment, data = gr) #sig due to treat with ll1 ll2
T0_aov <- aov(T0.area ~ Species, data = gr) #sig area due to species with or without ll1 ll2
T0_aov <- aov(T0.area ~ Accession, data = gr) #sig area due to accession with or without ll1 ll2
T0_aov <- aov(T0.area ~ Rep, data = gr) #sig area due to rep with or without ll1 ll2
T0_aov <- aov(T0.area ~ Treatment * Accession, data = gr) #affect treat and accession with ll1 ll2
T0_aov <- aov(T0.area ~ Treatment * Species, data = gr) #affect treat and accession with ll1 ll2
T0_aov <- aov(T0.area ~ Treatment * Species * Rep, data = gr) #all 3 have affect + species:rep interaction
str(T0_aov)
summary(T0_aov) 

light_aov <- aov(LightIntensity ~ Treatment, data = gr)
light_aov <- aov(LightIntensity ~ Treatment * Rep, data = gr)
str(light_aov)
summary(light_aov) #sig light between treat and rep

#find average T0 per accession per treatment, remove rep HL4
#REPHL4 <- gr %>% filter(!Rep %in% "HL4")
#t0_sum <- REPHL4 %>% filter(Treatment == "HL") %>% group_by(Accession) %>% summarise(mean(T0.area))

#t42 significant by anything?
T42_aov <- aov(T42.area ~ Treatment, data = gr) #sig due to treat with ll1 ll2
T42_aov <- aov(T42.area ~ Species, data = gr) #sig due to species with ll1 ll2
T42_aov <- aov(T42.area ~ Accession, data = gr) #sig due to accession with ll1 ll2
T42_aov <- aov(T42.area ~ Rep, data = gr) #sig due to rep with ll1 ll2
T42_aov <- aov(T42.area ~ Treatment * Accession, data = gr) #affect treat and accession with ll1 ll2
T42_aov <- aov(T42.area ~ Treatment * Species, data = gr) #affect treat and species with ll1 ll2
T42_aov <- aov(T42.area ~ Treatment * Species * Rep, data = gr) #all 3 have affect
str(T42_aov)
summary(T42_aov) 

#re-read in data with rows removed without T14, T21
gr <- read.csv("GR_Oct2022_noNArow.csv")

#subsetting
#remove rep 1 and 2 pre 2021
REPLL1 <- gr %>% filter(!Rep %in% "LL1")
REPLL2 <- REPLL1 %>% filter(!Rep %in% "LL2")
gr <- REPLL2
#remove t5, t6 and t9 cols
gr1 <- subset(gr, select=-c(T5.area, T6.area, T9.area))

#NORMALISE DATA TO START VALUE 0 becomes RGR_new
adjust <- gr %>% 
  group_by(T0.area) %>% 
  mutate(T0_new = T0.area - T0.area[1]) %>% 
  ungroup() #%>% 
adjust1 <- adjust %>% 
  group_by(T0.area) %>% 
  mutate(T7_new = T7.area - T0.area[1]) %>% 
  ungroup() #%>% 
adjust2 <- adjust1 %>% 
  group_by(T0.area) %>% 
  mutate(T14_new = T14.area - T0.area[1]) %>% 
  ungroup() #%>% 
adjust3 <- adjust2 %>% 
  group_by(T0.area) %>% 
  mutate(T21_new = T21.area - T0.area[1]) %>% 
  ungroup() #%>% 
adjust4 <- adjust3 %>% 
  group_by(T0.area) %>% 
  mutate(T28_new = T28.area - T0.area[1]) %>% 
  ungroup() #%>% 
adjust5 <- adjust4 %>% 
  group_by(T0.area) %>% 
  mutate(T42_new = T42.area - T0.area[1]) %>% 
  ungroup() #%>% 

#looked at normalised data for minus values
#recalculate RGR values
adjust5 %>% mutate()
RGR_new <- adjust5 %>% mutate(RGRlog = (log(T21_new-T14_new)/7))
RGR_new %>% arrange(desc(RGRlog))

#create difference between time 0 and 7
RGR_lag <- RGR_new %>% mutate(RGRlag = (log(T7_new-T0_new)/7))
RGR_lag %>% arrange(desc(RGRlag))

RGR_lag$RGRlag

#remove NAs
RGR_lag <- RGR_lag[complete.cases(RGR_lag$RGRlag), ]  
#RGR_lag RGR_log avgs
lag <- RGR_lag %>% group_by(Accession, Treatment) %>% summarise(mean(RGRlag))

#order by treatment
lag <- lag[order(lag$Treatment),]

#to add lag to density need to summarise per treatment and set as cols
try1 <- split(lag,cumsum(1:nrow(lag)%in%29)) #made 2 lists 0 and 1
str(try1)
#call the two lists
try1$`0`
try1$`1`
backtogether <- cbind(try1$`0`, try1$`1`)
backtogether = subset(backtogether, select = -c(2,4,5))
colnames(backtogether) <- c("Accession", "mean.HL_RGRlag", "mean.LL_RGRlag")
RGRlag <- backtogether

#save and add RGRlag onto density data

#write.csv(RGRlag, "RGRlag.csv")
#write.csv(RGR_lag, "RGR_norm.csv")

#do without 4 accessions, barplot by accession groups
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
raw <- read.csv("RGR_norm_nospREM.csv")
RGR_new <- raw

#for working out percent cov
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
percencov <- read.csv("RGR_norm_nospREM_forpercentagecoverage.csv")
library(dplyr)

names(percencov)
pcov <- percencov %>% mutate(T0_area = (T0_new / Area)*100) %>% 
  mutate(T7_area = (T7_new / Area)*100) %>%
    mutate(T14_area = (T14_new / Area)*100) %>%
    mutate(T21_area = (T21_new / Area)*100) %>%
           mutate(T28_area = (T28_new / Area)*100) %>%
             mutate(T42_area = (T42_new / Area)*100)

#new rgrlog and rgr area for %area data
pcov1 <- pcov %>% mutate(RGRpercovlog = (log(T21_area-T14_area)/7))
pcov1 %>% arrange(desc(RGRpercovlog))

#create difference between time 0 and 7
pcov2 <- pcov1 %>% mutate(RGRpercovlag = (log(T7_area-T0_area)/7))
pcov2 %>% arrange(desc(RGRpercovlag))

#remove NAs
RGR_new <- RGR_new[complete.cases(RGR_new$RGRlog), ] 

Summary <- RGR_new %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(RGRlog), stdev = sd(RGRlog), n= n(), maximum = max(RGRlog))
Summary

#before add letters for sig tukey groups
data_summ <- Summary %>% arrange(desc(mean))

aov <- aov(RGRlog ~ Treatment, data = RGR_new) #treatment affects
aov <- aov(RGRlog ~ Accession, data = RGR_new) #accession affects
aov <- aov(RGRlog ~ Treatment * Accession, data = RGR_new)#accession and treat have affect
summary(aov)
#Tukey works when not means data summary
tuk_out <- TukeyHSD(aov, "Accession", conf.level=.95)
tuk_out <- TukeyHSD(aov, "Treatment", conf.level=.95)
str(tuk_out)
tuk_out
plot(tuk_out , las=1 , col="brown") #not useful to visualise

#prepare name for letter production in multcompView
Tukey = TukeyHSD(aov)
print(Tukey)
print(tuk_out)

#install.packages("multcompView")
library("multcompView")

#add letters for significance of each group

Cld=multcompLetters4(aov, Tukey)
print(Cld)

Cld= as.data.frame.list(Cld$`Treatment:Accession`)
data_summ$Cld=Cld$Letters
print(data_summ)
View(data_summ)
#99 % conf and 95% same result ks04 ll own group, followed by sel1 hl, ks03 ll, ly03 ll, ks02 ll
#ks16 own group hl and ll generally slowest both conditions, of all accessions
#other accessions and treatments in same mid group
#sel1 hl only hl to match RGR of ll

#remove 4 species - ks03, ks04 best in ll, ks16 in hl ll own group slowest

#not working to plot letters on ggplot
#not plotting all the data, apparently NAs
RGR_new <- RGR_new[complete.cases(RGR_new$RGRlog), ] 

#not working to plot numbers of groups
qq=ggplot(RGR_new, aes(x=Accession, y=RGRlog, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(data = data_summ, aes(label=Cld,x=Accession, y=mean +stdev), position=position_dodge2(0.75), vjust = -6) + ylim(0,0.6) +
  labs(x="Ecotypes", y = "RGR_log (area mm2 per day)") + facet_grid(.~Treatment)
qq

#remove NAs
RGR_new <- RGR_new[complete.cases(RGR_new$RGR), ] 

Summary <- RGR_new %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(RGR), stdev = sd(RGR), n= n(), maximum = max(RGR))
Summary

#before add letters for sig tukey groups
data_summ <- Summary %>% arrange(desc(mean))

#work out std error using base r by treatment accession combo
aggregate(RGR ~ Accession + Treatment, data = RGR_new, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

aggregate(RGR ~ Species + Treatment, data = RGR_new, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- RGR_new %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(RGRlog),
    sd=sd(RGRlog)
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
  ylab("RGRlog by area (mm2 per day)") +
  xlab("Ecotype")
  #facet_wrap(~Treatment) harder to read accession var

#aov with summarised
aov <- aov(mean ~ Treatment * Accession, data = my_sum)
#aov with raw
aov <- aov(RGRlog ~ Treatment * Accession, data = raw)
summary(aov)
#accession and treat access interaction sig

#as there are significant differences between treatments, split by this
#LORNA
mylist1 <- split(grod3, grod3$Treatment) # split the data by levels within treat
HL <- mylist1[[1]] #two sets of data, one for each treat
LL<- mylist1[[2]]
#LORNA
#split data by rep
Replist <- split(RGR_lag, RGR_lag$Rep) # split the data by levels within rep
HL1 <- Replist[[1]] #one for each rep
HL2<- Replist[[2]]
HL3<- Replist[[3]]
HL4<- Replist[[4]]
LL1 <- Replist[[5]]
LL2<- Replist[[6]]
LL3<- Replist[[7]]
LL4<- Replist[[8]]

plot(LL$RGR ~ HL$RGR) #cant do as different no of observations

plot(HL1$RGR ~ HL2$RGR)
plot(HL2$RGR ~ HL3$RGR)
plot(HL3$RGR ~ HL4$RGR)

#dont know how to seperate means by LL and HL
#manually pasted columns into new columns belonging to observations

hl_llmeans <- read.csv("Summary_RGR_HL_LLsep.csv")
names(hl_llmeans)
#ll - hl
#hl_llmeans <- hl_llmeans %>% mutate(RGR_diff = (mean.LL - mean.HL)) %>% 
#  mutate(max_diff = (maximum.LL - maximum.HL))

#hl - ll
hl_llmeans <- hl_llmeans %>% mutate(RGR_diff = (mean.HL - mean.LL)) %>% 
  mutate(max_diff = (maximum.HL - maximum.LL))

aov <- aov(RGR_diff ~ Accession, data = hl_llmeans)
summary(aov)
#RGR_diff sig per accession

aov <- aov(RGR_diff ~ Accession, data = hl_llmeans)#accession and treat have affect
summary(aov)

#plot difference as scatterplot
plot(max_diff ~ RGR_diff, hl_llmeans)
#does the maximum growing one also have biggest difference in growth rate?
#RGR_diff	max_diff	+0.251229737 0.006585657
#weak positive
plot(maximum.HL ~ max_diff, hl_llmeans)
#maximum.HL	max_diff	0.650556249	+	2.19E-07
#strong positive and significant - HL having more of an affect in diff than LL
lm(maximum.HL ~ max_diff, hl_llmeans) 
str(lm(maximum.HL ~ max_diff, hl_llmeans))

range(hl_llmeans$maximum.HL)
#[1] 0.5800584 1.0876191
range(hl_llmeans$maximum.LL)
#[1] 0.7968253 1.1786437
#MORE VAR IN MAXIMUM HL


#probing differences in RGR between HL LL

hl_llmeans <- hl_llmeans %>% mutate(propchange1 = (mean.HL-mean.LL)/mean.HL)
hl_llmeans %>% arrange(desc(propchange1)) %>% select(Accession) %>% top_n(28)
hl_llmeans <- hl_llmeans %>% mutate(propchange2 = (mean.HL/mean.LL)/mean.HL)
hl_llmeans %>% arrange(desc(propchange2)) %>% select(Accession) %>% top_n(28)
hl_llmeans <- hl_llmeans %>% mutate(propchange3 = (mean.HL-mean.LL)/maximum.HL)
hl_llmeans %>% arrange(desc(propchange3)) %>% select(Accession) %>% top_n(28)
#as says online to do % proportion increase
hl_llmeans <- hl_llmeans %>% mutate(propchange4 = (RGR_diff/mean.LL)*100)
hl_llmeans %>% arrange(desc(propchange4)) %>% select(Accession) %>% top_n(28)
#can see the biggest increase (HL), vs ones staying same, vs ones biggest declines
#same order as absolute but easier to quantify

#probing differences in RGR between HL LL

hl_llmeans <- hl_llmeans %>% mutate(propchange1 = (mean.HL-mean.LL)/mean.HL)
hl_llmeans %>% arrange(desc(propchange1)) %>% select(Accession) %>% top_n(28)
hl_llmeans <- hl_llmeans %>% mutate(propchange2 = (mean.HL/mean.LL)/mean.HL)
hl_llmeans %>% arrange(desc(propchange2)) %>% select(Accession) %>% top_n(28)
hl_llmeans <- hl_llmeans %>% mutate(propchange3 = (mean.HL-mean.LL)/maximum.HL)
hl_llmeans %>% arrange(desc(propchange3)) %>% select(Accession) %>% top_n(28)
#as says online to do % proportion increase
hl_llmeans <- hl_llmeans %>% mutate(propchange4 = (RGR_diff/mean.LL)*100)
hl_llmeans %>% arrange(desc(propchange4)) %>% select(Accession) %>% top_n(28)
#can see the biggest increase (HL), vs ones staying same, vs ones biggest declines
#same order as absolute but easier to quantify

plot(propchange1 ~ propchange2, data=hl_llmeans)
plot(propchange1 ~ RGR_diff, data=hl_llmeans)
plot(propchange1 ~ propchange3, data=hl_llmeans) #similar
plot(propchange2 ~ propchange3, data=hl_llmeans)

#write.csv(hl_llmeans, "GR_RGR_propchange_oppcalc.csv")
#write.csv(hl_llmeans, "GR_RGR_propchange.csv")

#re read in version with 4 accessions
hl_llmeans <- read.csv("GR_RGR_propchange_oppcalc_spREM.csv") #minuses for prop
hl_llmeans <- read.csv("GR_RGR_propchange_spREM.csv")

#see highest performing in each condition
hl_llmeans %>% arrange(desc(mean.HL)) %>% select(Accession, mean.HL) %>% top_n(24)
hl_llmeans %>% arrange(desc(mean.LL)) %>% select(Accession, mean.LL) %>% top_n(24)


#group data into + or - RGR diff to see which direction performed better
hl_llmeans$RGR_diff > 0

which(hl_llmeans$RGR_diff > 0) #6 opp calc, 18
which(hl_llmeans$RGR_diff < 0) #18 opp calc, 6 otherwise

range(hl_llmeans$RGR_diff)
#[1] -0.3586310  0.2798193

#subset group as 3 groups based on RGR_diff values
nearzero <- hl_llmeans$RGR_diff[hl_llmeans$RGR_diff >= -0.1 & hl_llmeans$RGR_diff <= 0.05] #8 near zero
low <- hl_llmeans$RGR_diff[hl_llmeans$RGR_diff < -0.1] #15 low
high <- hl_llmeans$RGR_diff[hl_llmeans$RGR_diff > 0.05] #5 high
#length is 28 so included all in subset, 3 groups
#makes values of three group intergers but not associated to accessions

#try to include accession list with filtered RGR_diff
nearzero <- hl_llmeans$RGR_diff[hl_llmeans$RGR_diff >= -0.1 & hl_llmeans$RGR_diff <= 0.05] & select(hl_llmeans$Accession)

#2 groups basic
high <- which(hl_llmeans$RGR_diff < 0)
length(which(hl_llmeans$RGR_diff < 0))
#20 have loweer values hl than ll
low <- which(hl_llmeans$RGR_diff > 0)
length(which(hl_llmeans$RGR_diff > 0))
#8 have higher values hl than ll

hl_llmeans$Accession
#6 vs 18 with 24 obs

#add column for light level
#24 obs
hl_llmeans$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                      "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                      "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                      "HL", "HL", "LL", "LL")

#for 28 obs
hl_llmeans$Light <- c("HL","LL", "LL", "HL", "LL", "LL", "LL",
                      "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                      "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                      "HL", "HL", "LL", "LL", "HL", "LL", "HL")

#> 0 is APP2, KS06A, KS06B, KS15, KS16, KS22, LY02, SEL1 =8 PERFORMED BETTER IN HL
#finding those pretty much similar too
#ks15 and ly02 also higher NPQ than others

library(dplyr)
library(ggplot2)
#not working to call accessions
#hl_llmeans[which(hl_llmeans$Accession) | hl_llmeans$RGR_diff < 0, ]
hl_llmeans %>% arrange(desc(RGR_diff)) %>% select(Accession, RGR_diff) %>% top_n(24)
hl_llmeans %>% arrange(desc(RGR_diff)) %>% select(Accession) %>% top_n(28)
#shows same thing - which had highest values for difference

names(hl_llmeans)

#effect of treatment growth RGR
#edit this graph to color by light level
par(mfrow = c(2,5)) 
  ggplot(hl_llmeans, aes(y=RGR_diff, x=reorder(Accession, RGR_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
    theme_classic() +
    theme(axis.text.x=element_text(angle = 90))+
    theme(axis.text.y=element_text(angle = 90))+
    ylab("Difference growth rate (RGRlog HL - LL mm2)") +
    xlab("Ecotype")
  
  #if do ll - hl shows mostly negative affect of hl, only 5 inc gr hl
  #how to add error bars to this? need to use non avgd data
  
  # TO USE RAW DATA = PROBLEM AS NOT SAME NO OF REPS SO NOT
  #ALL INCLUDED IF MINUS ALL LIKE FOR LIKE TIMEPOINTS HL OR LL
  #GROWH TOGETHER
  raw <- read.csv("RGR_norm_allHL_LLrawRGRlog.csv")
  library(tidyr)
  
  lapply(raw, class)
  #raw$RGR_diff <-as.numeric(raw$RGR_diff) #changing values? doesnt work
  raw[,13] <- as.numeric(as.character(raw[,13]))
  #drop_na(raw) #doesnt work
raw_avgs <- raw %>%
    group_by(Accession) %>%
    summarise(avgRGR_diff = mean(RGR_diff, na.rm= TRUE), stdevRGR_diff = sd(RGR_diff, na.rm= TRUE))
  raw_avgs
  
  #calc se for ggplot std error
RAW_AVGS <- raw %>%
    group_by(Accession, Treatment) %>%
    summarise( 
      n=n(),
      mean=mean(RGR_diff, na.rm=TRUE),
      sd=sd(RGR_diff, na.rm=TRUE)
    ) %>%
    mutate( se=sd/sqrt(n))

#24 obs
RAW_AVGS$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                         "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                         "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                         "HL", "HL", "LL", "LL")

#std err graphs 
  ggplot(RAW_AVGS, aes(y=mean, x=reorder(Accession, mean), fill=EnvLight)) +
    geom_col(position='dodge', color='black') +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25) +
    theme_bw() +
    theme_classic() +
    theme(axis.text.x=element_text(angle = 90))+
    theme(axis.text.y=element_text(angle = 90))+
    ylab("Difference growth rate (RGRlog HL - LL mm2)") +
    xlab("Ecotype")
    
 #std dev
ggplot(RAW_AVGS, aes(y=mean, x=reorder(Accession, mean))) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)
    
  
  #graph fr raw
  par(mfrow = c(2,5)) 
  ggplot(raw, aes(y=RGR_diff, x=reorder(Accession, RGR_diff),  fill=Accession)) +
    geom_col(position='dodge') +
    theme_bw() +
    theme_classic() +
    theme(axis.text.x=element_text(angle = 90))+
    theme(axis.text.y=element_text(angle = 90))+
    ylab("Difference growth rate (RGRlog HL - LL mm2)") +
    xlab("Ecotype")
  
  
ggplot(hl_llmeans, aes(x=Light, y=RGR_diff, fill = Light)) +
  geom_boxplot(outlier.shape = NA) +
  ylab("Difference RGR log area gain (mm) between HL/LL")

t.test(RGR_diff ~ Light, data = hl_llmeans) # t test difference 

aov <- aov(RGR_diff ~ Light, data = hl_llmeans) #same result
summary(aov)

ggplot(hl_llmeans, aes(x=Light, y=propchange4, fill = Light)) +
  geom_boxplot(outlier.shape = NA)

#effect of treatment growth col gain
par(mfrow=c(2,2))
par(mar=c(0, 0, 0, 0.))
g <- ggplot(hl_llmeans, aes(y=propchange4, x=reorder(Accession, propchange4), fill=EnvLight))+
  geom_bar(stat="identity", color="black")+
  theme(axis.text.x=element_text(angle = 45), size = 20)+
  theme(axis.text.y=element_text(angle = 90), size = 10)+
  ylab("% Growth change RGRlog area") +
  xlab("Ecotype") +
  theme_minimal()
g+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))

ggplot(hl_llmeans, aes(y=propchange4, x=reorder(Accession, propchange4),  fill=EnvLight))+
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()

#automate correlations
lms <- expand.grid(2:11, 2:11)
lms_names <- expand.grid(names(hl_llmeans)[2:11], names(hl_llmeans)[2:11])
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

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Area_RGR")
#write.csv(all_output, "hl_ll_correlation.csv")

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
#write.csv(Summary, file="Summary_RGR.csv")
#write.csv(Summary_lag, file="Summary_RGR.csv")


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

ggplot(Summary, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), position = position_dodge(0.9), width = 0.25)

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
  
levels(RGR_lag$Treatment)
levels(RGR_lag$Treatment) <- c("HL", "LL")  
  
  # numbers of reps not working
library(dplyr)
RGR_lag %>% group_by(Treatment, Accession) %>% dplyr::summarize(n = n()) %>% print(n = Inf)
RGR_lag$Treatment
  
names(RGR_lag)
  
  
  #change to parameter you want
RGR_boxplot <- ggplot(RGR_lag,aes(x = reorder(Accession,RGR, FUN = median),y=RGR)) +
    
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
  RGR_boxplot
  
  
  RGR_lag$Accession <- as.factor(RGR_lag$Accession)
  RGR_lag$Treatment <- as.factor(RGR_lag$Treatment)
  
  #remove if water
  #param$Accession == "Water"
  
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
RGR_sum <- RGR_new %>% group_by(Treatment, Accession) %>% summarise_all(mean) %>% arrange(desc(RGRlog)) #add desc as argument to do descending

#summarised plot
barplot(RGR_sum$RGRlog, 
        names.arg = RGR_sum$Accession,
        horiz = T, las = 1,
        xlim = c(0, 1.5))
box()

#summarised plot
ggplot(RGR_sum, aes(x=Accession, y=RGRlog, fill=Treatment, group=Treatment)) +
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
#accessions <- (RGR_sum[,2])

#RUN THIS IN ORDER AND HASH OUT RELEVANT VERSIONS
#accessions <- (HL[,1])
accessions <- (LL[,1])

#to see hl and ll seperately need to split by treatment here
#HL <- RGR_sum %>% filter(Treatment == "HL") %>% select(Accession, T0_new, T7_new, T14_new, T21_new, T28_new, T42_new)
LL <- RGR_sum %>% filter(Treatment == "LL") %>% select(Accession, T0_new, T7_new, T14_new, T21_new, T28_new, T42_new)
#HL$Treatment <- NULL
LL$Treatment <- NULL

#RUN TO DO % coverage
percovavg1 <- pcov2 %>% group_by(Treatment, Accession) %>% summarise_all(mean, na.rm = TRUE) #add desc as argument to do descending

HL <- percovavg1 %>% filter(Treatment == "HL") %>% select(Accession, T0_area, T7_area, T14_area, T21_area, T28_area, T42_area)
#LL <- percovavg1 %>% filter(Treatment == "LL") %>% select(Accession, T0_area, T7_area, T14_area, T21_area, T28_area, T42_area)

#HL$Treatment <- NULL
#LL$Treatment <- NULL

#for percent cov
accessions <- (HL[,2])
#accessions <- (LL[,2])

sub <- (percovavg1[,27:32]) #just t0 to t42 area columns to show all data
#to see hl and ll seperately need to split by treatment 
sa_comb <- sub

#sub <- (RGR_sum[,17:22]) #just t0 to t42 columns to show all data
#to see hl and ll seperately need to split by treatment 
#sa_comb <- sub
sa_comb <- HL
#sa_comb <- LL
names(sa_comb)

#for percen cov
sa_comb <- sa_comb[-1]
#sa_comb <- sa_comb[-1]
#sa_comb <- HL[-1]
#sa_comb <- sa_comb[-1]

sa_comb$Accession <- NULL

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

#for working with percen area
#replace characters in a string
row.names(sadf) <- str_replace_all(row.names(sadf), "T0_", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "_area", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "T", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "area", "0")
#looks for pattern and give replacement

#replace characters in a string
row.names(sadf) <- str_replace_all(row.names(sadf), "T0_", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "_new", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "T", "")
row.names(sadf) <- str_replace_all(row.names(sadf), "new", "0")
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

sadf$SEL1
sadf$KS02
sadf$time

par(mfrow = c(1,1)) #no of rows and cols in plot display
plot(KS02 ~ time, sadf, type = "n", xlab = "Time (days)", ylab = "RGR frond area T14-T21/7 (mm2)")
lines(KS03 ~ time, sadf) # one at a time
lines(KS04 ~ time, sadf) 
lines(KS06A ~ time, sadf) # one at a time
lines(KS06B ~ time, sadf) 

#do something to every i
plot(KS02 ~ time, sadf, type = "n", 
     main = "Growth rate in HL",
     xlab = "Time (days)", 
     ylab = expression(paste("Frond area coverage (mm"^2*")")), 
     axes = F,
     xlim = c(0, 50),
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

#graphs to save
tiff('Growthrate_LL_lineplot_cov.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(2,2)) 
#redo plot for ll with 1 main color, colors for ones of interest
plot(KS02 ~ time, sadf, type = "n", 
     main = "Growth rate in LL",
     xlab = "Time (days)", 
     ylab = "Percentage frond area coverage", 
     axes = F,
     xlim = c(0, 50),
     ylim = c(0, 100))
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
dev.off()
#dev.new()


#dev.new()

tiff('Growthrate_HL_lineplot_cov.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(2,2)) 
#redo plot for hl with 1 color
plot(KS02 ~ time, sadf, type = "n", 
     main = "Growth rate in HL",
     xlab = "Time (days)", 
     ylab = "Percentage frond area coverage", 
     axes = F,
     xlim = c(0, 50),
     ylim = c(0, 100))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
#text(15.8,23.8, "RGRlog", srt=45) #place txt and rotate it
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
dev.off()
