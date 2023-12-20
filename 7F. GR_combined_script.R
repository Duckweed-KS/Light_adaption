#combine GR, biomass and Col count data
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

#install.packages("dplyr")
#install.packages("tidyr")
library(dplyr)
library(tidyr)

#NOT BEEN CHANGED TO NEW CALCS FOR DIFFS RGR, BIOMASS,COLS

#growth rate data done manuallyf
gr <- read.csv("Combo_data.csv")
gr <- read.csv("Combo_data_spREM.csv")

#NEW SPECIES L JP AND ENV LIGHT INC
gr <- read.csv("Combo_data_spREM_NEWsp.csv")

#JOIN TOTAL DATASETS
#gr <- read.csv("RGR_norm_4stick.csv")
#cc <- read.csv("Col_countsRGR_4stick.csv")
#bio <- read.csv("Biomass_minusstart_4stick.csv")

#join 2 data sets together (by accession)
#names(gr)
#names(gr)[13] <- "Rep" #name col for observation names
#names(cc)[10] <- "Rep"
#names(bio)[6] <- "Rep"
#gr %>% str()
#cc %>% str()
#gr$Rep <- as.numeric(gr$Rep) #growth rate
#cc$Rep <- as.numeric(cc$Rep) # col counts
#bio$Rep <- as.numeric(bio$Rep) # biomass FW FDW
#View(gr$Sample) #named numeric 1-28
#tail(data.frame(gr), 20) #show last 20 rows
#gr <- gr[1:28, ] #not done anything

#TRIALLING JOINING DATA TOGETHER 
#joined <- gr %>% inner_join(cc, by = c("Rep"))
#joined <- joined %>% inner_join(bio, by = c("Rep"))
#worked but not done exactly as code says
#joined2 <- gr %>% full_join(cc, by = c("Rep"))
#joined2 <- joined2 %>% full_join(bio, by = c("Rep"))

#df = merge(x=gr,y=cc,by="Rep",all=TRUE)
#df

#left join
#df1 = merge(x=gr,y=cc,by="Rep",all.x=TRUE)
#df1

#should be ~400 observations not 5000?

#join <- left_join(gr, cc) %>%
#  distinct(Rep, Accession, Treatment, .keep_all =TRUE)
#not working to add onto
#join <- left_join(join, bio) %>%
#  distinct(Rep, Accession, Treatment, .keep_all =TRUE)

#names(joined)
#unique(joined$Rep) #12 unique
#unique(joined2$Rep) #17 unique
#length(joined$Rep)
#length(join$Rep)
#unique(join$Rep)

names(gr)

#lag area with all
plot(RGRlag ~ RGRlog, gr) #good cor
#plot(FW..mg. ~ RGRlag, gr)
plot(FDW..mg. ~ RGRlag, gr)
plot(FW_mod ~ RGRlag, gr)

#correlation between FW and FDW
corr <- cor.test(gr$RGRlag, gr$RGRlog,
                 method = "pearson"
)
corr$p.value
corr$estimate

corr <- cor.test(gr$RGRlag, gr$FW_mod,
                 method = "pearson"
)
corr$p.value
corr$estimate

#log area with all
plot(FW_mod ~ RGRlog, gr)
plot(FDW..mg. ~ RGRlog, gr)
#plot(FW..mg. ~ RGRlog, gr)

corr <- cor.test(gr$RGRlog, gr$FW_mod,
                 method = "pearson"
)
corr$p.value
corr$estimate

corr <- cor.test(gr$RGRlag, gr$FDW..mg.,
                 method = "pearson"
)
corr$p.value
corr$estimate #not good

plot(FW_mod ~ FDW..mg., gr)

corr <- cor.test(gr$FW_mod, gr$FDW..mg.,
                 method = "pearson"
)
corr$p.value
corr$estimate #not good

#col lag with all
#doesnt work well as mostly low var
plot(Col_lag_RGR ~ RGRlog, gr)
plot(Col_lag_RGR ~ RGRlag, gr)
plot(FW_mod ~ Col_lag_RGR, gr)
plot(FDW..mg. ~ Col_lag_RGR, gr)
plot(Col_lag_RGR ~ Col_RGR, gr)
plot(Col_lag_RGR ~ Col_RGRnonlog, gr)

#col log with all
plot(Col_RGR ~ RGRlog, gr)
plot(Col_RGR ~ RGRlag, gr)
plot(FW_mod ~ Col_RGR, gr)
plot(FDW..mg. ~ Col_RGR, gr)
plot(Col_RGR ~ Col_RGRnonlog, gr) #arc relationship

#col diff calc with all
plot(Col_RGRnonlog ~ RGRlog, gr)
plot(Col_RGRnonlog ~ RGRlag, gr)
plot(FW_mod ~ Col_RGRnonlog, gr)
plot(FDW..mg. ~ Col_RGRnonlog, gr)

#with 4 accessions included
#does treatment affect the data or just accession?
aov <- aov(RGRlog ~ Treatment*Accession, data = gr)  #treat and access
aov <- aov(Col_RGR ~ Treatment*Accession, data = gr) #treat and access
aov <- aov(FW_mod ~ Treatment*Accession, data = gr)  #not sig
aov <- aov(FDW..mg. ~ Treatment*Accession, data = gr) #not sig
summary(aov)

#without 4
#does treatment affect the data or just accession?
aov <- aov(RGRlog ~ Treatment*Accession, data = gr)  #treat and access
aov <- aov(Col_RGR ~ Treatment*Accession, data = gr) #treat and access (more treat)
aov <- aov(FW_mod ~ Treatment*Accession, data = gr)  #not sig
aov <- aov(FDW..mg. ~ Treatment*Accession, data = gr) #not sig
summary(aov)

#does sp affect data?
aov <- aov(RGRlog ~ Treatment*Species, data = gr)  #treat and species
aov <- aov(Col_RGR ~ Treatment*Species, data = gr) #treat and sp and interaction
aov <- aov(FW_mod ~ Treatment*Species, data = gr)  #just sp sig
aov <- aov(FDW..mg. ~ Treatment*Species, data = gr) #not sig
summary(aov)

#does habitat affect data?
aov <- aov(RGRlog ~ EnvLight*Species, data = gr)  #treat and species
aov <- aov(Col_RGR ~ EnvLight*Species, data = gr) #treat and sp and interaction
aov <- aov(FW_mod ~ EnvLight*Species, data = gr)  #just sp sig
aov <- aov(FDW..mg. ~ EnvLight*Species, data = gr) #not sig
summary(aov)

#does habitat affect data?
aov <- aov(RGRlog ~ EnvLight*Accession, data = gr)  #treat and species
aov <- aov(Col_RGR ~ EnvLight*Accession, data = gr) #treat and sp and interaction
aov <- aov(FW_mod ~ EnvLight*Accession, data = gr)  #just sp sig
aov <- aov(FDW..mg. ~ EnvLight*Accession, data = gr) #not sig
summary(aov)

#does sp affect data?
aov <- aov(RGRlog ~ Treatment*Species, data = gr)  #treat and species
aov <- aov(Col_RGR ~ Treatment*Species, data = gr) #treat and sp and interaction
aov <- aov(FW_mod ~ Treatment*Species, data = gr)  #just sp sig
aov <- aov(FDW..mg. ~ Treatment*Species, data = gr) #not sig
summary(aov)

boxplot(RGRlog ~ Species, data = gr)
boxplot(RGRlog ~ EnvLight, data = gr)
boxplot(RGRlog ~ Treatment*Species, data = gr)
boxplot(RGRlog ~ EnvLight*Species, data = gr)
boxplot(Col_RGR ~ Species, data = gr)
boxplot(Col_RGR ~ Treatment*Species, data = gr)
boxplot(FW_mod ~ Species, data = gr) 

#tukey tests
tuk_out <- TukeyHSD(aov, "Species", conf.level=.95)
tuk_out

#using 24 obs not avgd
t.test(RGRlog ~ Treatment, data = gr)  #sig ll higher
t.test(Col_RGR ~ Treatment, data = gr) #sig ll higher
t.test(FW_mod ~ Treatment, data = gr) #ll higher ns
t.test(FDW..mg. ~ Treatment, data = gr) #ll higher ns

#consider as whole rather than treatment parts
#nicer plot
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)
ggplot(gr, aes(RGRlag, RGRlog, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("RGR lag (mm2 per day)") +
  ylab("RGR log (mm2 per day)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

ggplot(gr, aes(Col_RGR, RGRlog, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Col gain per day") +
  ylab("RGR log (mm2 per day)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

#include in paper #remove labels
myggp <- ggplot(gr, aes(Col_RGR, RGRlog, label = Accession)) +    # ggplot2 plot with labels
  geom_point(pch=3) +
  geom_smooth(method='lm') +
  xlab("RGR log col gain per day") +
  ylab("RGR area log (mm2 per day)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 2)
myggp + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))

ggplot(gr, aes(Col_RGR, RGRlog, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Col gain per day") +
  ylab("RGR lag (mm2 per day)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

ggplot(gr, aes(FW_mod, RGRlog, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Fresh weight (mg)") +
  ylab("RGR log (mm2 per day)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

ggplot(gr, aes(FDW..mg., RGRlog, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("RGR log (mm2 per day)") +
  ylab("Freeze dry weight (mg)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

#doesnt work well with col lag
ggplot(gr, aes(Col_lag_RGR, RGRlag, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Col_lag_RGR (col gain per day)") +
  ylab("RGR lag (mm2 per day)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

library(dplyr)
#compare ranking data for col gain and rgr area all raw
Col_log_RGR_ord <- gr %>% arrange(desc(Col_RGR)) %>% select(Accession, Col_RGR)
RGRlog_ord <- gr %>% arrange(desc(RGRlog)) %>% select(Accession, RGRlog)
rankings <- cbind(Col_log_RGR_ord, RGRlog_ord)
#top 4 accessions ly03, sel1, nuff1, moor1

#RGR_sum <- rowMeans(gr[, 29], na.rm = TRUE)
#to average data based on accession
gr$Accession <- as.factor(gr$Accession)
by_access <- gr %>% group_by(Accession)
RGR_sum <- by_access %>% summarise_each(funs(mean(., na.rm = TRUE)))

#do again when summarise data
#compare ranking data for col gain and rgr area
Col_log_RGR_ord <- RGR_sum %>% arrange(desc(Col_RGR)) %>% select(Accession, Col_RGR)
RGRlog_ord <- RGR_sum %>% arrange(desc(RGRlog)) %>% select(Accession, RGRlog)
rankings <- cbind(Col_log_RGR_ord, RGRlog_ord)

#plot avg gr data color by light level
RGR_sum$Light <- c("LL", "LL", "HL", "LL", "LL", "LL",
                   "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                   "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                   "HL", "HL", "LL", "LL")

#28 obs
RGR_sum$Light <- c("HL","LL", "LL", "HL", "LL", "LL", "LL",
                      "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                      "HL", "HL", "LL", "LL", "LL", "LL", "HL",
                      "HL", "HL", "LL", "LL", "HL", "LL", "HL")

#looking at performance without splitting into treatments
par(mfrow = c(1,2)) 
ggplot(RGR_sum, aes(y=Col_RGR, x=reorder(Accession, Col_RGR),  fill=Light)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("mean log RGR (colony gain)") +
  xlab("Ecotype")
  

par(mfrow = c(2,5)) 
ggplot(RGR_sum, aes(y=RGRlog, x=reorder(Accession, RGRlog),  fill=Light)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("RGR area gain (mm2)") +
  xlab("Ecotype")

#L. japonicas fastest, L. minutas and S. poly slower. L min middle

par(mfrow = c(2,5)) 
ggplot(gr, aes(y=RGRlog, x=reorder(Species, RGRlog))) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("RGR area gain (mm2)") +
  xlab("Ecotype")

#cant see species difference without splitting into LL and HL responses

#use gr to make boxplots for RGRlog and Col RGR split by treatment and species
#read in GR_Combo_data_spREM_NEWsp_splitbytreat.csv