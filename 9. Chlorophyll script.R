#script to explore hl ll dw growth rate data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Chlorophyll_Oct2022.csv")
gr <- read.csv("Chlorophyll_Oct2022_nogenera.csv")
gr <- read.csv("Chlorophyll_Oct2022_nogenera_nospecies.csv")
#new ljp included
gr <- read.csv("Chlorophyll_Oct2022_nogenera_nospecies_NEWsp.csv")

tail(gr)
head(gr)
names(gr)
#names
#[1] "Accession"      "Within_Rep"     "Treatment"     
#[4] "Chl.a..mg.g."   "Chl.b..mg.g."   "Carot..mg.g."  
#[7] "Start_date"     "Rep"            "Species"       
#[10] "LightIntensity" "SpecIntensity"  "X"             
#[13] "X.1"

#check classes of all cols
sapply(gr, class)

#need to remove last 3 cols as empty
#gr <- gr[-c(12:13)]
#gr <- gr[-c(451:462),]

library(dplyr)
library(stringr)
library(ggplot2)

gr %>% select(Accession, Treatment) #just displays them
unique(gr$Accession, gr$Treatment) 
unique(gr$Accession) #32 #now 28 without genera
unique(gr$Rep, gr$Treatment) #10 reps
unique(gr$Rep) #no LL2 or LL6 or HL 4 in this set
length(unique(gr$Accession, gr$Treatment)) #54 inc genera
#49 when no genera
length(unique(gr$Accession,order = ascending))
#32 includes genera #28 no genera
length(unique(gr$Species,order = ascending))
#8 as includes genera #5 no genera

as.factor(gr$Accession)
as.factor(gr$Species)
as.factor(gr$Treatment)
as.factor(gr$Rep)
as.factor(gr$LightIntensity)

#levels(gr$Species) <- c("S. polyrhiza", "L. minor", 
#                     "L. minuta", "L. gibba", "L. turionifera",
#                      "S. int", "L. punc", "W. arr")
gr$Accession <- factor(gr$Accession,levels = c("KS02", "KS03", "KS04", "KS06A","KS06B",
                        "KS09", "KS12", "KS13", "KS14", "KS15", "KS16",
                        "KS17", "KS18", "KS20", "KS21", "KS22", "KS25",
                        "KS27", "KS28", "KS29", "LY01A", "LY01B", "LY02",
                        "LY03", "MOOR1", "APP2", "NUFF1", "SEL1",
                        "S. int", "L. punc", "L. min", "W. arr"))

#no genera
#levels(gr$Species) <- c("S. polyrhiza", "L. minor", 
#                       "L. minuta", "L. gibba", "L. turionifera")
#levels(gr$Accession) <- c("KS02", "KS03", "KS04", "KS06A","KS06B",
#                          "KS09", "KS12", "KS13", "KS14", "KS15", "KS16",
#                          "KS17", "KS18", "KS20", "KS21", "KS22", "KS25",
#                          "KS27", "KS28", "KS29", "LY01A", "LY01B", "LY02",
#                          "LY03", "MOOR1", "APP2", "NUFF1", "SEL1")

#with genera
#filter by treatment = LL 
gr %>% filter(Treatment == "LL")
ggplot(gr, aes(x = Accession, y = Chl.a..mg.g.,)) + geom_bar(stat = "identity")
ggplot(gr, aes(y=Chl.a..mg.g., x=Accession)) +
  geom_col(position='dodge', color='black')

names(gr)
#affects of other variables
aov <- aov(Chl.a..mg.g.~Rep,data=gr) #sig for all chl a, b and carot
aov <- aov(Carot..mg.g.~Within_Rep,data=gr) #not sig apart from for carot
summary(aov)

tukey <- TukeyHSD(aov, "Rep", conf.level=.95)
str(tukey)
tukey
#LL4-LL1 -1.6427338 -2.51150980 -0.7739577526 0.0000008
#LL5-LL1 -4.2904438 -5.97173559 -2.6091520275 0.0000000
#LL4-LL3 -0.8680980 -1.73687399  0.0006780582 0.0503282
#LL5-LL3 -3.5158080 -5.19709977 -1.8345162167 0.0000000
#LL5-LL4 -2.6477100 -4.31883408 -0.9765859846 0.0000739
#disclude LL1 or LL3

boxplot(Chl.a..mg.g.~Accession*Rep,data=gr)

#how different reps are
boxplot(Chl.a..mg.g.~Within_Rep,data=gr)
boxplot(Chl.a..mg.g.~Rep,data=gr)
boxplot(Chl.a..mg.g.~Rep*Within_Rep,data=gr)
boxplot(Chl.b..mg.g.~Rep,data=gr)
boxplot(Carot..mg.g.~Rep,data=gr) #big difference in LL4

#look at HL LL seperately as looks like LL declined over time, between reps
#for chl a and chl b

#CHL A plots
boxplot(Chl.a..mg.g.~Species,data=gr)
boxplot(Chl.a..mg.g.~Accession,data=gr)
boxplot(Chl.a..mg.g.~Rep,data=gr)

boxplot(Chl.a..mg.g.~Species*Treatment,data=gr) #more chl a in LL
boxplot(Chl.a..mg.g.~Accession*Treatment,data=gr) 
boxplot(Chl.a..mg.g.~Accession,data=gr) 
boxplot(Chl.a..mg.g.~Rep,data=gr) 
boxplot(Chl.a..mg.g.~Species,data=gr)

#CHL B Plots
boxplot(Chl.b..mg.g.~Species,data=gr)
boxplot(Chl.b..mg.g.~Accession,data=gr)
boxplot(Chl.b..mg.g.~Rep,data=gr)

boxplot(Chl.b..mg.g.~Species*Treatment,data=gr) #more chl b in LL
boxplot(Chl.b..mg.g.~Accession*Treatment,data=gr) 
boxplot(Chl.b..mg.g.~Accession,data=gr) 
boxplot(Chl.b..mg.g.~Rep,data=gr) 
boxplot(Chl.b..mg.g.~Species,data=gr)

#Car Plots
boxplot(Carot..mg.g.~Species,data=gr)
boxplot(Carot..mg.g.~Accession,data=gr)
boxplot(Carot..mg.g.~Rep,data=gr)

boxplot(Carot..mg.g.~Species*Treatment,data=gr)
boxplot(Carot..mg.g.~Accession*Treatment,data=gr) 
boxplot(Carot..mg.g.~Accession,data=gr) 
boxplot(Carot..mg.g.~Rep,data=gr) 
boxplot(Carot..mg.g.~Species,data=gr)


#adding new variables
#gr %>% mutate()
#grmod <- gr %>% group_by(Accession) %>% mutate(Col_RGR = (log(T21_cols-T14_cols)/7))
#gr %>% group_by(Accession, Treatment) %>% mutate(PRI = (log(T21.area-T14.area)/7))
#grmod <- grmod %>% group_by(Accession) %>% mutate(Col_RGRnonlog = ((T21_cols-T14_cols)/7))

#remove na cases
#grmod2 <- na.omit(grmod)
#gr <- gr[complete.cases(gr$Carot..mg.g), ]   
#grod3 <- grmod[complete.cases(grmod$Col_RGRnonlog), ]  
#remove all NA cases for RGR column in data frame = 20 rem

range(gr$Chl.a..mg.g.)
#between 0.60490 10.99875
range(gr$Chl.b..mg.g.)
#between -0.02722  5.76273
#0.02722 5.76273
range(gr$Carot..mg.g.)
#between -1.257858  3.040272
#0.000281187 3.040272146

#minus values to remove
#which(gr$Chl.b..mg.g. < 0) #just 103 is negative - KS04 HL
#which(gr$Carot..mg.g. < 0) #128, 211, 280 negative
#LY01A KS09, KS21 LL
#removed and resaved

#add new cols for a+b total and a/b ratio, car per chl

gr <- gr %>% mutate(Chla_b= (Chl.a..mg.g. + Chl.b..mg.g.)) %>% 
  mutate(Chlab_rat = (Chl.a..mg.g. / Chl.b..mg.g.)) %>%
  mutate(Car_rat = (Carot..mg.g./Chla_b))

#look at differences between col nos in treatments
ChlA_aov <- aov(Chl.a..mg.g. ~ Treatment, data = gr) 
ChlA_aov <- aov(Chl.a..mg.g. ~ Accession, data = gr) 
ChlA_aov <- aov(Chl.a..mg.g. ~ Treatment*Accession, data = gr)
ChlA_aov <- aov(Chl.a..mg.g. ~ Treatment*Species, data = gr) 
ChlA_aov <- aov(Chl.a..mg.g. ~ Species, data = gr) 
ChlA_aov <- aov(Chl.a..mg.g. ~ Rep, data = gr)
ChlA_aov <- aov(Chl.a..mg.g. ~ Treatment*Rep, data = gr)
ChlA_aov <- aov(Chl.a..mg.g. ~ Treatment*Species*Rep, data = gr) 
summary(ChlA_aov)

lm <- lm(Chl.a..mg.g. ~ Treatment*Species, data = gr) 
lm <- lm(Chl.b..mg.g. ~ Treatment*Species, data = gr) 
lm <- lm(Carot..mg.g. ~ Treatment*Species, data = gr) 
lm <- lm(Chl_ab ~ Treatment*Species, data = gr) 
lm <- lm(Chlab_rat ~ Treatment*Species, data = gr) #sp *
lm <- lm(Car_rat ~ Treatment*Species, data = gr) #~***, ***,*
summary(lm)

ChlB_aov <- aov(Chl.b..mg.g. ~ Treatment*Accession, data = gr) 

ChlAB_rat_aov <- aov(Chlab_rat ~ Treatment*Accession, data = gr) #ns
ChlAB_rat_aov <- aov(Chlab_rat ~ Treatment*Species, data = gr) #now * sig species
summary(ChlAB_rat_aov)

TotChl_aov <- aov(Chla_b ~ Treatment*Accession, data = gr)
TotChl_aov <- aov(Chla_b ~ Treatment*Species, data = gr)
TotChl_aov <- aov(Chla_b ~ EnvLight*Species, data = gr) 
TotChl_aov <- aov(Chla_b ~ Accession*Treatment, data = gr)
TotChl_aov <- aov(Chla_b ~ Species*Treatment, data = gr)
TotChl_aov <- aov(Chla_b ~ Species*EnvLight, data = gr) 
summary(TotChl_aov)

Car_aov <- aov(Carot..mg.g. ~ Treatment*Accession, data = gr) 
Car_aov <- aov(Carot..mg.g. ~ Treatment*Species, data = gr) 
Car_aov <- aov(Carot..mg.g. ~ EnvLight*Species, data = gr) 
summary(Car_aov)

Car_rat_aov <- aov(Car_rat ~ Treatment*Accession, data = gr)
Car_rat_aov <- aov(Car_rat ~ Treatment*Species, data = gr)
Car_rat_aov <- aov(Car_rat ~ EnvLight*Species, data = gr) #ns

summary(TotChl_aov)
summary(ChlAB_rat_aov)
summary(Car_rat_aov)

summary(ChlB_aov)
summary(Car_aov)

#SIGNIFICANCE
#TREATMENT ***, ACCESSION **, TREAT * ACCESS ALL ***,
#SPECIES ***, REP ***, TREAT * REP *BOTH ***

#Tukey works when not means data summary
tuk_out <- TukeyHSD(ChlA_aov, "Rep", conf.level=.95)
tuk_out <- TukeyHSD(ChlA_aov, "Treatment", conf.level=.95)
tuk_out <- TukeyHSD(ChlA_aov, "Accession", conf.level=.95)
tuk_out <- TukeyHSD(ChlA_aov, "Species", conf.level=.95)
str(tuk_out)
tuk_out
plot(tuk_out , las=1 , col="brown") #not useful to visualise
#treat significant with or without accession or species included
#rep significant LL1-HL1 LL1-HL2 LL1-HL3 LL1-HL4 LL1-LL2 LL1-LL3 LL1-LL4 LL1-LL5
#accession not sig with treat or not
#species not sig with treat or not
#rep still sig when include treatment - HL all similar, LL1 different from others


aov <- aov(Chl.a..mg.g. ~ Treatment*Species, data = gr) 
aov <- aov(Chl.b..mg.g. ~ Treatment*Species, data = gr) 
aov <- aov(Carot..mg.g. ~ Treatment*Species, data = gr) 
aov <- aov(Chl_ab ~ Treatment*Species, data = gr) 
aov <- aov(Chlab_rat ~ Treatment*Species, data = gr) #sp *
aov <- aov(Car_rat ~ Treatment*Species, data = gr) #~***, ***,*
summary(aov)
tukey <- TukeyHSD(aov, conf.level=.95)
tukey
#Chl A L minor, L. minu and L jap all diff between treats. LL - L mino diff S poly, L. jap diff to L. minu and S poly
#Chl B L mino, L. minu and L jap all diff between. LL - L jap and L mino both diff S poly
#ChlAB L mino, L jap diff between
#ChlAB rat LL - Diff between L jap and L minu
#Carot L jap and L mino diff between
#Car rat L jap L mino L minu L turio diff between, LL - s poly diff l turio, l jap diff s poly, l mino dif s poly, l mino diff l minu

aov <- aov(Chl.a..mg.g. ~ Treatment*EnvLight, data = gr) #treat env light int
aov <- aov(Chl.b..mg.g. ~ Treatment*EnvLight, data = gr) # treat env light int
aov <- aov(Carot..mg.g. ~ Treatment*EnvLight, data = gr) #just treat
aov <- aov(Chl_ab ~ Treatment*EnvLight, data = gr) #treat
aov <- aov(Chlab_rat ~ Treatment*EnvLight, data = gr) #nothing
aov <- aov(Car_rat ~ Treatment*EnvLight, data = gr) #treat
summary(aov)
tukey <- TukeyHSD(aov, conf.level=.95)
tukey

try <- read.csv("Fluo+Chl_stitched+LL_HL_sep.csv")
try$EnvLight <- factor(try$EnvLight, levels = c("dLL", "dHL"))

#FOR FLUO PARAM
#incorp species USING THESE FOR PLOTS
aov <- aov(QY_max_LL ~ Species*EnvLight, data = try) #ns
aov <- aov(QY_max_HL ~ Species*EnvLight, data = try) #species, env light int
aov <- aov(Fv.Fm_L3_LL ~ Species*EnvLight, data = try) #ns
aov <- aov(Fv.Fm_L3_HL ~ Species*EnvLight, data = try) #species, env light int
aov <- aov(Fv.Fm_L5_LL ~ Species*EnvLight, data = try) # ns
aov <- aov(Fv.Fm_L5_HL ~ Species*EnvLight, data = try) # ns
aov <- aov(Fv.Fm_L11_LL ~ Species*EnvLight, data = try) # ns
aov <- aov(Fv.Fm_L11_HL ~ Species*EnvLight, data = try) # ns species***
aov <- aov(NPQ_L3_LL ~ Species*EnvLight, data = try) # species ***
aov <- aov(NPQ_L3_HL ~ Species*EnvLight, data = try) # species ***
aov <- aov(NPQ_L5_LL ~ Species*EnvLight, data = try) # species **
aov <- aov(NPQ_L5_HL ~ Species*EnvLight, data = try) # species ***, int *
aov <- aov(NPQ_L11_LL ~ Species*EnvLight, data = try) # species *
aov <- aov(NPQ_L11_HL ~ Species*EnvLight, data = try) # species ***, species envlight **
summary(aov)

#FOR PIGMENT EXT
#incorp species USING THESE FOR PLOTS
aov <- aov(Chl.a..mg.g._LL ~ Species*EnvLight, data = try) #ns 0.12
aov <- aov(Chl.a..mg.g._HL ~ Species*EnvLight, data = try) #** 0.002
aov <- aov(Chl.b..mg.g._LL ~ Species*EnvLight, data = try) #* 0.04
aov <- aov(Chl.b..mg.g._HL ~ Species*EnvLight, data = try) #** 0.005
aov <- aov(TotChl_LL ~ Species*EnvLight, data = try) # ns 0.08
aov <- aov(TotChl_HL ~ Species*EnvLight, data = try) # ** 0.002
aov <- aov(ChlAB_rat_LL ~ Species*EnvLight, data = try) # ns 0.2
aov <- aov(ChlAB_rat_HL ~ Species*EnvLight, data = try) # ns 0.3
aov <- aov(Carot..mg.g._LL ~ Species*EnvLight, data = try) # ns 0.9
aov <- aov(Carot..mg.g._HL ~ Species*EnvLight, data = try) #** 0.006
aov <- aov(Car_rat_LL ~ Species*EnvLight, data = try) # ns 0.3
aov <- aov(Car_rat_HL ~ Species*EnvLight, data = try) # ns 0.07
summary(aov)

#single factor
aov <- aov(Chl.a..mg.g._LL ~ EnvLight, data = try) #* 0.05
aov <- aov(Chl.a..mg.g._HL ~ EnvLight, data = try) #** 0.006
aov <- aov(Chl.b..mg.g._LL ~ EnvLight, data = try) #* 0.01
aov <- aov(Chl.b..mg.g._HL ~ EnvLight, data = try) #* 0.01
aov <- aov(TotChl_LL ~ EnvLight, data = try) # * 0.03
aov <- aov(TotChl_HL ~ EnvLight, data = try) # ** 0.006
aov <- aov(ChlAB_rat_LL ~ EnvLight, data = try) # ns 0.2
aov <- aov(ChlAB_rat_HL ~ EnvLight, data = try) # ns 0.3
aov <- aov(Carot..mg.g._LL ~ EnvLight, data = try) # ns 0.9
aov <- aov(Carot..mg.g._HL ~ EnvLight, data = try) #** 0.009
aov <- aov(Car_rat_LL ~ EnvLight, data = try) # ns 0.1
aov <- aov(Car_rat_HL ~ EnvLight, data = try) # ns 0.2
summary(aov)
tuk_out <- TukeyHSD(aov, conf.level=.95)
tuk_out

par(mfrow = c(1, 1))
#tiff('ChlCarot_boxplots_EnvLight_sigraw.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('ChlCarot_boxplots_EnvLight_sigraw.pdf', width=14, height=12)
#par(mfrow = c(3, 3),  mar=c(2,4.5,2,2))
#plot boxplots by species and treatment
par(mfrow = c(3, 4))
boxplot(Chl.a..mg.g._LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 13), xlab = "",
        ylab = "Chl a (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 12.90, "LL", cex=2.75)
text(1, 11.5, "n.s", cex=2.75)
text(2.1, 12.90, "p = 0.12", cex=2)
title("A.", cex.main=2, adj=0)
stripchart(Chl.a..mg.g._LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl.a..mg.g._HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 13), xlab = "",
        ylab = "Chl a (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 12.90, "HL", cex=2.75)
text(1, 10.2, "**", cex=2.75)
text(2.1, 11.90, "p = 0.002", cex=2)
#title("A.", cex.main=2, adj=0)
stripchart(Chl.a..mg.g._HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl.b..mg.g._LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 4), xlab = "",
        ylab = "Chl b (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 3.90, "LL", cex=2.75)
text(2, 3.4, "*", cex=2.75)
text(2.1, 3.90, "p = 0.04", cex=2)
title("B.", cex.main=2, adj=0)
stripchart(Chl.b..mg.g._LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl.b..mg.g._HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 4), xlab = "",
        ylab = "Chl b (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 3.90, "HL", cex=2.75)
text(1, 3.4, "**", cex=2.75)
text(2.1, 3.90, "p = 0.005", cex=2)
#title("A.", cex.main=2, adj=0)
stripchart(Chl.b..mg.g._HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(TotChl_LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 16), xlab = "",
        ylab = "Total Chl a+b (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 15.90, "LL", cex=2.75)
text(1, 14.7, "n.s", cex=2.75)
text(2.1, 15.90, "p = 0.08", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(TotChl_LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(TotChl_HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 16), xlab = "",
        ylab = "Total Chl a+b (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 15.90, "HL", cex=2.75)
text(1, 11.90, "**", cex=2.75)
text(2.1, 15.90, "p = 0.002", cex=2)
#title("A.", cex.main=2, adj=0)
stripchart(TotChl_HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(ChlAB_rat_LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 10), xlab = "",
        ylab = "Chl a:b",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 9.90, "LL", cex=2.75)
text(1, 8.8, "n.s", cex=2.75)
text(2.1, 9.90, "p = 0.2", cex=2)
title("D.", cex.main=2, adj=0)
stripchart(ChlAB_rat_LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(ChlAB_rat_HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 10), xlab = "",
        ylab = "Chl a:b",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 9.90, "HL", cex=2.75)
text(1, 8.6, "n.s", cex=2.75)
text(2.1, 9.90, "p = 0.3", cex=2)
#title("A.", cex.main=2, adj=0)
stripchart(ChlAB_rat_HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Carot..mg.g._LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 4), xlab = "",
        ylab = "Total carotenoids (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 3.90, "LL", cex=2.75)
text(1, 3.1, "n.s", cex=2.75)
text(2.1, 3.90, "p = 0.9", cex=2)
title("E.", cex.main=2, adj=0)
stripchart(Carot..mg.g._LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Carot..mg.g._HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 4), xlab = "",
        ylab = "Total carotenoids (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 3.90, "HL", cex=2.75)
text(1, 3.2, "**", cex=2.75)
text(2.1, 3.90, "p = 0.006", cex=2)
#title("A.", cex.main=2, adj=0)
stripchart(Carot..mg.g._HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Car_rat_LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 0.6), xlab = "",
        ylab = "Car:Chl",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 0.58, "LL", cex=2.75)
text(1, 0.48, "n.s", cex=2.75)
text(2.1, 0.58, "p = 0.3", cex=2)
title("F.", cex.main=2, adj=0)
stripchart(Car_rat_LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Car_rat_HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 0.6), xlab = "",
        ylab = "Car:Chl",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 0.58, "HL", cex=2.75)
text(1, 0.45, "n.s", cex=2.75)
text(2.1, 0.58, "p = 0.07", cex=2)
#title("A.", cex.main=2, adj=0)
stripchart(Car_rat_HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

dev.off()

#PLOT RAW BOXPLOT DATA FOR SPECIES EACH PIGMENT PARAMTER
#change order species
gr <- read.csv("Chlorophyll_Oct2022_nogenera_nospecies_NEWsp.csv")
gr$Species <- factor(gr$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))
gr$Treatment <- factor(gr$Treatment, levels = c("LL", "HL"))

par(mfrow = c(1, 1))
#tiff('ChlCarot_boxplots_Species_sigraw.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('ChlCarot_boxplots_Species_sigraw.pdf', width=14, height=12)
#par(mfrow = c(3, 3),  mar=c(2,4.5,2,2))
#plot boxplots by species and treatment
par(mfrow = c(2, 3))
boxplot(Chl.a..mg.g. ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 13), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 12.90, "a", cex=2)
text(2, 12.90, "a", cex=2)
text(1, 12.99, "***", cex=2)
text(2, 12.99, "***", cex=2)
text(4, 12.99, "***", cex=2)
text(3, 12.90, "ab", cex=2)
text(4, 12.90, "b", cex=2)
text(5, 12.90, "b", cex=2)
text(7, 12.90, "a", cex=2)
text(8, 12.90, "a", cex=2)
text(9, 12.90, "a", cex=2)
text(10, 12.90, "a", cex=2)
text(11, 12.90, "a", cex=2)
#Chl A L minor, L. minu and L jap all diff between treats. LL - L mino diff S poly, L. jap diff to L. minu and S poly
title("A.", cex.main=2, adj=0)
#FvFm_L3  HL - L jap L mino, L jap L turion, L minu L jap, between L jap with itself
stripchart(Chl.a..mg.g. ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl.b..mg.g. ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 6), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 5.90, "a", cex=2)
text(2, 5.90, "a", cex=2)
text(1, 5.98, "***", cex=2)
text(2, 5.98, "***", cex=2)
text(4, 5.98, "***", cex=2)
text(3, 5.90, "ab", cex=2)
text(4, 5.90, "ab", cex=2)
text(5, 5.90, "b", cex=2)
text(7, 5.90, "a", cex=2)
text(8, 5.90, "a", cex=2)
text(9, 5.90, "a", cex=2)
text(10, 5.90, "a", cex=2)
text(11, 5.90, "a", cex=2)
#Chl B L mino, L. minu and L jap all diff between. LL - L jap and L mino both diff S poly
stripchart(Chl.b..mg.g. ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chla_b ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 15), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 14.90, "a", cex=2)
text(2, 14.90, "a", cex=2)
text(1, 14.90, "***", cex=2)
text(2, 14.90, "***", cex=2)
text(3, 14.90, "a", cex=2)
text(4, 14.90, "a", cex=2)
text(5, 14.90, "a", cex=2)
text(7, 14.90, "a", cex=2)
text(8, 14.90, "a", cex=2)
text(9, 14.90, "a", cex=2)
text(10, 14.90, "a", cex=2)
text(11, 14.90, "a", cex=2)
#ChlAB L mino, L jap diff between
stripchart(Chla_b ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chlab_rat ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 5), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 4.90, "ab", cex=2)
text(2, 4.90, "a", cex=2)
#text(2, 4.90, "***", cex=2)
text(3, 4.90, "ab", cex=2)
text(4, 4.90, "b", cex=2)
text(5, 4.90, "ab", cex=2)
text(7, 4.90, "a", cex=2)
text(8, 4.90, "a", cex=2)
text(9, 4.90, "a", cex=2)
text(10, 4.90, "a", cex=2)
text(11, 4.90, "a", cex=2)
title("A.", cex.main=2, adj=0)
#ChlAB rat LL - Diff between L jap and L minu
stripchart(Chlab_rat ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Carot..mg.g. ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 4), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 3.90, "a", cex=2)
text(2, 3.90, "a", cex=2)
text(1, 3.98, "***", cex=2)
text(2, 3.98, "***", cex=2)
text(3, 3.90, "a", cex=2)
text(4, 3.90, "a", cex=2)
text(5, 3.90, "a", cex=2)
text(7, 3.90, "a", cex=2)
text(8, 3.90, "a", cex=2)
text(9, 3.90, "a", cex=2)
text(10, 3.90, "a", cex=2)
text(11, 3.90, "a", cex=2)
title("A.", cex.main=2, adj=0)
#Carot L jap and L mino diff between
stripchart(Carot..mg.g. ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Car_rat ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 0.5), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 0.47, "b", cex=2)
text(2, 0.47, "b", cex=2)
text(1, 0.47, "***", cex=2)
text(2, 0.47, "***", cex=2)
text(3, 0.47, "***", cex=2)
text(4, 0.47, "***", cex=2)
text(3, 0.47, "b", cex=2)
text(4, 0.47, "a", cex=2)
text(5, 0.47, "a", cex=2)
text(7, 0.47, "a", cex=2)
text(8, 0.47, "a", cex=2)
text(9, 0.47, "a", cex=2)
text(10, 0.47, "a", cex=2)
text(11, 0.47, "a", cex=2)
title("A.", cex.main=2, adj=0)
#Car rat L jap L mino L minu L turio diff between, LL - s poly diff l turio, l jap diff s poly, l mino dif s poly, l mino diff l minu
stripchart(Car_rat ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

dev.off()
#no species differences at hl, differences due to envlight? previous adaption

#subsetting
#remove genera
#gr_sp <- gr[-c(64:69), ] #remove L. min, L. punc
#gr <- gr_sp
#gr_spn <- gr[-c(81:86), ] #not working to remove

#save without genera

#remove last 3 reps
#gr$Rep == "SLL1", "SLL2", "SLL3"
#gr_rep <- gr %>% filter(!Rep %in% "SLL1", "SLL2", "SLL3")


#work out difference in Chl A between LL and HL treats per accession
#Chla_new <- gr %>% arrange(desc(Chl.a..mg.g.))

names(gr)
  
#carotenoid ratio worked out by total carotenoid / total chloro

library(ggplot2)

bp <- ggplot(gr, aes(x=Chl.a..mg.g., y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=Chl.b..mg.g., y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=Carot..mg.g., y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=Car_rat, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=Chlab_rat, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)


bp <- ggplot(gr, aes(x=Chla_b, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

#SUMMARISING BY ACCESSION AND TREATMENT
sum <- gr %>% group_by(Accession, Treatment) %>% summarise_all(mean)
#now have summarised RGR rate per access and treatment and other cols
sum$Chl.a..mg.g.
sum$Chl.a..mg.g.
sum$Carot..mg.g.
sum$Chla_b
sum$Chlab_rat
sum$Car_rat

plot(Chl.a..mg.g.~LightIntensity,data=gr)
plot(Chla_b~LightIntensity,data=sum)
plot(Chla_b~Treatment,data=sum) #total chl higher in LL clear graph
plot(Chlab_rat~Treatment,data=sum)
plot(Car_rat~Treatment,data=sum) #car ratio much higher in HL clear graph

Chla_b_aov <- aov(Chla_b~Treatment,data=gr) #sig
Chla_b_aov <- aov(Chla_b~Accession,data=gr) #accession 1 star sig
Chla_b_aov <- aov(Chla_b~Treatment*Accession,data=gr) #all 3 sig
Chla_b_aov <- aov(Chla_b~Treatment*Species,data=gr) #sig seperately
Chlab_rat_aov <- aov(Chlab_rat~Treatment*Accession,data=gr) #not sig
Car_rat_aov <- aov(Car_rat~Treatment*Accession,data=gr) #both treat + access sig
summary(Chla_b_aov)
summary(Chlab_rat_aov) #ratio doesnt change in HL or LL
summary(Car_rat_aov) 

#all Chl A values for all independent observations in HL
# not working
#HL <- gr %>% filter(Treatment == "HL") %>% select(Chl.a..mg.g.)
#HL %>% arrange(Accession) #order all HL by accession
#gr %>% arrange(desc(Chl.a..mg.g.)) %>% select(Treatment, Accession) %>% top_n(10)

#LL <- gr %>% filter(Treatment == "LL") %>% select(Chl.a..mg.g.)
#LL %>% arrange(Accession) #order all HL RGR values by accession
#gr %>% arrange(desc(Chl.a..mg.g.)) %>% select(Treatment, Accession) %>% top_n(10)

#Summarised AVERAGED Chl A for all in HL and LL
#avg_HL <- gr %>% filter(Treatment == "HL") %>% group_by(Accession) %>% summarise(mean(Chl.a..mg.g.))
#avg_LL <- gr %>% filter(Treatment == "LL") %>% group_by(Accession) %>% summarise(mean(Chl.a..mg.g.))

#more detailed summary
Summary <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Chl.a..mg.g.), stdev = sd(Chl.a..mg.g.), n= n(), maximum = max(Chl.a..mg.g.))
Summary

Summary1 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Chl.b..mg.g.), stdev = sd(Chl.b..mg.g.), n= n(), maximum = max(Chl.b..mg.g.))
Summary1

Summary2 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Carot..mg.g.), stdev = sd(Carot..mg.g.), n= n(), maximum = max(Carot..mg.g.))
Summary2

Summary3 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Chla_b), stdev = sd(Chla_b), n= n(), maximum = max(Chla_b))
Summary3

Summary4 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Chlab_rat), stdev = sd(Chlab_rat), n= n(), maximum = max(Chlab_rat))
Summary4

Summary5 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(mean = mean(Car_rat), stdev = sd(Car_rat), n= n(), maximum = max(Car_rat))
Summary5

#combine 5 summaries manually

Summ <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5)

#write.csv(Summ, "Chloro_Summary.csv")
#write.csv(Summ, "Chloro_Summary_incSp.csv") #given colnames manually

range(Summ$maximum4)
range(Summ$maximum3)

#find top 10 performing accession RGR values
#6A top of HL, LY01B top of LL
#12 bottom of HL and LL - slowest growing

gr %>% 
  ggplot(aes(Chl.a..mg.g., Accession, col = Treatment)) + 
  geom_point(alpha = 0.8) +
  facet_wrap(~Treatment)


t.test(Chl.a..mg.g. ~ Treatment, data = gr) # sig
t.test(Chl.b..mg.g. ~ Treatment, data = gr) # sig
t.test(Carot..mg.g. ~ Treatment, data = gr) # sig
t.test(Chla_b ~ Treatment, data = gr) # sig
t.test(Chlab_rat ~ Treatment, data = gr) # not sig
t.test(Car_rat ~ Treatment, data = gr) # sig

#boxplot faceted by treatment
ChlA_boxplot <- ggplot(gr,aes(x = reorder(Accession,Chl.a..mg.g., FUN = median),y=Chl.a..mg.g.)) +
  
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
  labs(title =expression("(a) Chl A (mg/g)"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Chl.a..mg.g.g))+
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
ChlA_boxplot

#chl b
ChlB_boxplot <- ggplot(gr,aes(x = reorder(Accession,Chl.b..mg.g., FUN = median),y=Chl.b..mg.g.)) +
  
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
  labs(title =expression("(a) Chl B (mg/g)"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Chl.b..mg.g.))+
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
ChlB_boxplot

#chl b
Car_boxplot <- ggplot(gr,aes(x = reorder(Accession,Carot..mg.g., FUN = median),y=Carot..mg.g.)) +
  
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
  labs(title =expression("(a) Carot (mg/g)"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Carot..mg.g.))+
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
Car_boxplot

#plot by rep
#boxplot faceted by treatment
ChlArep_boxplot <- ggplot(gr,aes(x = reorder(Accession,Chl.a..mg.g., FUN = median),y=Chl.a..mg.g.)) +
  
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
  labs(title =expression("(a) Chl A (mg/g)"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Chl.a..mg.g.g))+
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
ChlArep_boxplot

#plot scatter plot to show relationship with rep
#doesnt work as datapoints belong to that rep
library(ggplot2)
ggplot(gr, aes(Chl.a..mg.g., Rep, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)


#for inc Sp manually sorted into HL LL and created summary
#to add lag to density need to summarise per treatment and set as cols

#order by treatment
Summ <- read.csv("Chloro_Summary_incSp.csv") 
Summ <- Summ[order(Summ$Treatment),]
try1 <- split(Summ,cumsum(1:nrow(Summ)%in%25)) #made 2 lists 0 and 1
str(try1)
#call the two lists
try1$`0`
try1$`1`
backtogether <- cbind(try1$`0`, try1$`1`) #not working as no HL for some accessions
#backtogether = subset(backtogether, select = -c(2,4,5))
#colnames(backtogether) <- c("Accession", "mean.HL_RGRlag", "mean.LL_RGRlag")
#RGRlag <- backtogether

#read in new file with grouped HL LL means and max
hl_llmeans <- read.csv("Chloro_Summary_SepHLLL.csv")

names(hl_llmeans)
library(dplyr)
#CHANGED THIS TO OTHER WAY ROUND SO match other data
hl_llmeans %>% mutate()
#take hl away from ll calculated paramters ChlAB tot ChlAB rat Car rat
hl_llmeans <- hl_llmeans %>% mutate(ChlABtot_diff = (Chl_ab_HL_mean - Chl_ab_LL_mean)) %>% 
  mutate(ChlABrat_diff = (Chlab_rat_HL_mean - Chlab_rat_LL_mean)) %>%
  mutate(Car_rat_diff = (Car_rat_HL_mean - Car_rat_LL_mean))
#take hl away from ll ChlA
hl_llmeans <- hl_llmeans %>% mutate(ChlB_diff = (Chl_b_HL_mean - Chl_b_LL_mean)) %>% 
  mutate(ChlB_max_diff = (Chl_b_HL_maximum - Chl_b_LL_maximum))
#take hl away from ll ChlB
hl_llmeans <- hl_llmeans %>% mutate(ChlA_diff = (Chl_a_HL_mean - Chl_a_LL_mean)) %>% 
  mutate(ChlA_max_diff = (Chl_a_HL_maximum - Chl_a_LL_maximum))
#take hl away from ll Carot
hl_llmeans <- hl_llmeans %>% mutate(Carot_diff = (Carot_HL_mean - Carot_LL_mean)) %>% 
  mutate(Carot_max_diff = (Carot_HL_maximum - Carot_LL_maximum))

hl_llmeans <- hl_llmeans %>% mutate(prop_ChlAB = (ChlABtot_diff/Chl_ab_HL_mean)*100)
hl_llmeans <- hl_llmeans %>% mutate(prop_ChlA = (ChlA_diff/Chl_a_HL_mean)*100)
hl_llmeans <- hl_llmeans %>% mutate(prop_ChlB = (ChlB_diff/Chl_b_HL_mean)*100)
hl_llmeans <- hl_llmeans %>% mutate(prop_Chlabrat = (ChlABrat_diff/Chlab_rat_HL_mean)*100)
hl_llmeans <- hl_llmeans %>% mutate(prop_Carrat = (Car_rat_diff/Car_rat_HL_mean)*100)
hl_llmeans <- hl_llmeans %>% mutate(prop_Car = (Carot_diff/Carot_HL_mean)*100)

write.csv(hl_llmeans, "Chloro_Summary_sepHLLL_withdiffs.csv")

#chlab_rat
names(hl_llmeans)

#plot difference as scatterplot
plot(ChlA_max_diff ~ ChlA_diff, hl_llmeans)
plot(ChlA_diff ~ ChlA_max_diff, hl_llmeans)
plot(ChlB_max_diff ~ ChlB_diff, hl_llmeans)
plot(ChlB_diff ~ ChlB_max_diff, hl_llmeans)
plot(Carot_max_diff ~ Carot_diff, hl_llmeans)
plot(Carot_diff ~ Carot_max_diff, hl_llmeans)

#range(hl_llmeans$maximum.HL)
#[1]  0.7142857 24.4285714
#range(hl_llmeans$maximum.LL)
#[1]  1.00000 17.28571
#MORE VAR IN MAXIMUM HL

#group data into + or - RGR diff to see which direction performed better
hl_llmeans$ChlA_diff > 0
hl_llmeans$ChlB_diff > 0
hl_llmeans$Carot_diff > 0

hl_llmeans$ChlABtot_diff > 0 #23 negative = decreased AB content
hl_llmeans$ChlABrat_diff > 0 #some pos some neg

which(hl_llmeans$ChlABtot_diff > 0) #23 accessions positive
which(hl_llmeans$ChlABtot_diff < 0) #1 accession negative
#[1] 8 KS13 gained more chlAB total in HL than LL
#other acessions more chlorosis

which(hl_llmeans$ChlABrat_diff > 0) #13 accessions positive
which(hl_llmeans$ChlABrat_diff < 0) #11 accessions negative
#[1]  2  4  5  7 13 14 15 17 18 21 23
#KS03 KS06A KS06B KS12 KS18 KS20 KS21 KS25 KS27 LY01A LY02

which(hl_llmeans$Car_rat_diff < 0) #all 24 accession negative
#all gain carotenoids in HL 

#BAR PLOTS
library(ggplot2)

#work out std error using base r by treatment accession combo
aggregate(Car_rat ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Car_rat),
    sd=sd(Car_rat)
  ) %>%
  mutate( se=sd/sqrt(n))

#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25) +
  theme_bw()
#car higher in all accessions in HL, KS12 similar

#CHLA B ratios
#22 28 LY03 15 decreased, KS03 increased
#work out std error using base r by treatment accession combo
aggregate(Chlab_rat ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Chlab_rat),
    sd=sd(Chlab_rat)
  ) %>%
  mutate( se=sd/sqrt(n))

#plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#ks15, ks22, ks28, ly03 higher ab ratio in LL. KS03 and 06B higher in HL
#different responses between accessions, quantify differences

#CHL A B TOT
#work out std error using base r by treatment accession combo
aggregate(Chla_b ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#how much chl a and b lost due to chlorosis, biggest decline = worst affected

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Chla_b),
    sd=sd(Chla_b)
  ) %>%
  mutate( se=sd/sqrt(n))

#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#ks 13 only one higher total in HL
#ks04, 09, 15, LY02 highest chlab in LL

range(hl_llmeans$ChlABtot_diff)
# -0.5412347  6.8266250
range(hl_llmeans$ChlABrat_diff)
# -0.3508654  1.0041849

hl_llmeans$Accession

which(hl_llmeans$ChlA_diff > 0) #23 accessions positive
#[1]  1  2  3  4  5  6  7  9 10 11 12 13 14 15 16 17 18 19 20 21
#[21] 22 23 24
which(hl_llmeans$ChlB_diff > 0)
#23 positive [1]  1  2  3  4  5  6  7  9 10 11 12 13 14 15 16 17 18 19 20 21
#[21] 22 23 24
which(hl_llmeans$Carot_diff > 0)
#3 positive 3, 9, 11
#KS14, KS04, KS16 

which(hl_llmeans$ChlA_diff < 0) #1 accession negative
#[8] KS13
which(hl_llmeans$ChlB_diff < 0) #1 accession negative
#[8] KS13
which(hl_llmeans$Carot_diff < 0) #21 accession negative
#[1]  1  2  4  5  6  7  8 10 12 13 14 15 16 17 18 19 20 21 22 23
#[21] 24
which(hl_llmeans$ChlABrat_diff < 0) #13
which(hl_llmeans$Car_rat_diff < 0) #0

range(hl_llmeans$ChlA_diff)
#[1] -0.371073  5.309238
range(hl_llmeans$ChlB_diff)
#[1] -0.1701613  1.5232723
range(hl_llmeans$Carot_diff)
#[1] -0.9970552  0.9359410

hl_llmeans$Accession

#add EnvLight factor in for 24 accessions
hl_llmeans$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                     "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                     "HL", "HL", "LL", "LL", "LL", "LL", "HL", "HL",
                     "HL", "LL", "LL")

library(ggplot2)
library(ggpubr)
#effect of treatment Chl A
ggplot(hl_llmeans, aes(y=ChlA_diff, x=reorder(Accession, ChlA_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in Chl A (mg/g)") +
  xlab("Ecotype")

ggplot(hl_llmeans, aes(y=ChlB_diff, x=reorder(Accession, ChlB_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in Chl B (mg/g)") +
  xlab("Ecotype")

#chl ab tot difference
#for paper
ggplot(hl_llmeans, aes(y=ChlABtot_diff, x=reorder(Accession, ChlABtot_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90))+
  theme(axis.text.y=element_text(angle = 90))+
  ylab("Difference in Total Chl (mg/g)") +
  xlab("Ecotype")

ggplot(hl_llmeans, aes(y=ChlABtot_diff, x=reorder(Accession, ChlABtot_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size=14)) +
  ylab("Difference in Total Chl (mg/g)") +
  xlab("Ecotype")

#chl a b ratio diff
ggplot(hl_llmeans, aes(y=ChlABrat_diff, x=reorder(Accession, ChlABrat_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in ratio Chl A/B (mg/g)") +
  xlab("Ecotype")

#diff in carot ratio
ggplot(hl_llmeans, aes(y=Car_rat_diff, x=reorder(Accession, Car_rat_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in ratio Carot/Chl (mg/g)") +
  xlab("Ecotype")

  #effect of treatment Carot
ggplot(hl_llmeans, aes(y=Carot_diff, x=reorder(Accession, Carot_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90))+
  theme(axis.text.y=element_text(angle = 90))+
  ylab("Difference in Carotenoids (mg/g)") +
  xlab("Ecotype")

#automate correlations
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Chloro")
#cant include 4 and 5
lms <- expand.grid(2:31, 2:31)
lms_names <- expand.grid(names(hl_llmeans)[2:31], names(hl_llmeans)[2:31])
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

write.csv(all_output, "hl_ll_correlation.csv")


#work out std error using base r by treatment accession combo
aggregate(Carot.mg.g. ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

aggregate(Carot..mg.g. ~ Species + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Carot..mg.g.),
    sd=sd(Carot..mg.g.)
  ) %>%
  mutate( se=sd/sqrt(n))

ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)


aov <- aov(mean ~ Treatment * Accession, data = my_sum)
aov <- aov(mean ~ Treatment, data = my_sum)
aov <- aov(mean ~ Accession, data = my_sum)
summary(aov)
#just treatment effect for chl a when done alone

#spearmans rank
access.ranked <- data.frame(cbind(rank(RGR_lag$RGR, ties.method = 'average'),
                                  rank(RGR_lag$Rep, ties.method = 'average')))
colnames(access.ranked) <- c('RGR', 'Rep')

corr <- cor.test(x=RGR_lag$RGR, y=RGR_lag$Rep, method = 'spearman')
corr

names(hl_llmeans)
#show Tot chla+b and Carot diff as bar plot
par(mfrow = c(2,5)) 
ggplot(hl_llmeans, aes(y=ChlABtot_diff, x=reorder(Accession, ChlABtot_diff))) +
  geom_col(position='dodge', color='black') +
  theme_bw() +
  theme_classic() +
  ylab("Total Chlorophyll decline in HL mg/g") +
  xlab("Ecotype")

#show Tot chla+b and Carot diff as bar plot
par(mfrow = c(2,5)) 
ggplot(hl_llmeans, aes(y=Car_rat_diff, x=reorder(Accession, Car_rat_diff))) +
  geom_col(position='dodge', color='black') +
  theme_bw() +
  theme_classic() +
  ylab("Total Carotenoids increase in HL (mg/g)") +
  xlab("Ecotype")

hl_llmeans <- hl_llmeans %>% mutate(Car_diff = (Carot_HL_mean/Carot_LL_mean))

par(mfrow = c(2,5)) 
ggplot(hl_llmeans, aes(y=Car_diff, x=reorder(Accession, Car_diff))) +
  geom_col(position='dodge', color='black') +
  theme_bw() +
  theme_classic() +
  ylab("Total Carotenoids increase in HL (mg/g)") +
  xlab("Ecotype")

#addition of light classification groupings for sites
hl_llmeans$Accession

#t test affect env light on differences
aov <- aov(Car_rat_diff ~ EnvLight * Accession, data = hl_llmeans)
aov <- aov(Car_rat_diff ~ EnvLight * Accession, data = hl_llmeans)
t.test(Car_rat_diff ~ EnvLight, data = hl_llmeans) # not sig
t.test(Carot_diff ~ EnvLight, data = hl_llmeans) # not sig
t.test(ChlABrat_diff ~ EnvLight, data = hl_llmeans) # not sig
t.test(ChlABtot_diff ~ EnvLight, data = hl_llmeans) #not sig

t.test(Chl_ab_HL_mean ~ EnvLight, data = hl_llmeans) #sig
t.test(Chl_ab_LL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Chlab_rat_HL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Chlab_rat_LL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Car_rat_HL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Car_rat_LL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Carot_HL_mean ~ EnvLight, data = hl_llmeans) #sig
t.test(Carot_LL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Chl_a_HL_mean ~ EnvLight, data = hl_llmeans) #sig
t.test(Chl_a_LL_mean ~ EnvLight, data = hl_llmeans) #not sig
t.test(Chl_b_HL_mean ~ EnvLight, data = hl_llmeans) #sig
t.test(Chl_b_LL_mean ~ EnvLight, data = hl_llmeans) #not sig

#only hl data sig, no diff data or ratio data

#light not significant by t. tests for parameters when summarised diff
summary(aov)

#look at rankings for variables
names(hl_llmeans)
hl_llmeans %>% arrange(desc(ChlA_diff)) %>% select(Accession, ChlA_diff) %>% top_n(24)
#shows same thing - which had highest values for difference
#only KS13 gained chl A in HL treat, rest declined
hl_llmeans$Accession
hl_llmeans %>% arrange(desc(ChlB_diff)) %>% select(Accession, ChlB_diff) %>% top_n(24)
#shows same thing - which had highest values for difference
#only KS13 gained chl B in HL treat, rest declined
hl_llmeans$Accession
hl_llmeans %>% arrange(desc(Carot_diff)) %>% select(Accession, Carot_diff) %>% top_n(24)
#shows same thing - which had highest values for difference
#some gained carot, 3 decreased in HL
hl_llmeans %>% arrange(desc(ChlABtot_diff)) %>% select(Accession, ChlABtot_diff) %>% top_n(24)
#shows same thing - which had highest values for difference
#KS13 gained tot chl
hl_llmeans %>% arrange(desc(ChlABrat_diff)) %>% select(Accession, ChlABrat_diff) %>% top_n(24)
#shows same thing - which had highest values for difference
#some pos some negative for chl a/b ratio
hl_llmeans %>% arrange(desc(Car_rat_diff)) %>% select(Accession, Car_rat_diff) %>% top_n(24)
#shows same thing - which had highest values for difference
#some gained carot, 3 decreased in H