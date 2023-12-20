#script to explore fluor data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Fluorcam.csv")
gr <- read.csv("Fluorcam+fqfm.csv") #including new variables
gr <- read.csv("Fluorcam+fqfm_REMsum2022_spr2021.csv") #rem obs
gr <- read.csv("Fluorcam+fqfm_REMsum2022_spr2021_REMoutliers.csv") #rem obs
#new l jp species + env light
gr <- read.csv("Fluorcam+fqfm_REMsum2022_spr2021_REMoutliers_NEWsp.csv") #rem obs
tail(gr)
head(gr)
names(gr)

#check classes of all cols
sapply(gr, class)

#need to remove last 3 cols as empty
#gr <- gr[-c(12:13)]
#gr <- gr[-c(451:462),]

#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("plyr")
#install.packages("stats")

library(dplyr)
library(stringr)
library(ggplot2)

gr %>% select(Accession, Treatment) #just displays them
unique(gr$Accession, gr$Treatment) 
unique(gr$Accession) #31
unique(gr$Rep, gr$Treatment) #16
unique(gr$Rep) #16
length(unique(gr$Accession, gr$Treatment)) #69
length(unique(gr$Accession,order = ascending)) #31
length(unique(gr$Species,order = ascending)) #5 species

as.factor(gr$Accession)
as.factor(gr$Species)
as.factor(gr$Treatment)
as.factor(gr$Rep)
as.factor(gr$Within_Rep)
as.factor(gr$Plate)

#from lorna
library(plyr)
library(stats)
mu <- ddply(gr, "Rep", summarise, grp.mean=mean(QY_max))
head(mu)
mu <- ddply(gr, "Rep", summarise, grp.mean=mean(NPQ_L11))
head(mu)
#for more than 1 condition
OneWay_A <- aov(QY_max~Rep, data = gr) #heading +12 during the day
summary(OneWay_A)
#used for 2 conditions
res <- wilcox.test(QY_max ~ Treatment, data = gr,
                   exact = FALSE)
res
#W = 18546, p-value = 7.486e-09

#levels(gr$Species) <- c("S. polyrhiza", "L. minor", 
#                     "L. minuta", "L. gibba", "L. turionifera",
#                      "S. int", "L. punc", "W. arr")
#levels(gr$Accession) <- c("KS02", "KS03", "KS04", "KS06A","KS06B",
#                        "KS09", "KS12", "KS13", "KS14", "KS15", "KS16",
#                        "KS17", "KS18", "KS20", "KS21", "KS22", "KS25",
#                        "KS27", "KS28", "KS29", "LY01A", "LY01B", "LY02",
#                        "LY03", "MOOR1", "APP2", "NUFF1", "SEL1",
#                        "S. int", "L. punc", "L. min", "W. arr")

#no genera
#levels(gr$Species) <- c("S. polyrhiza", "L. minor", 
#                       "L. minuta", "L. gibba", "L. turionifera")
#levels(gr$Accession) <- c("KS02", "KS03", "KS04", "KS06A","KS06B",
#                          "KS09", "KS12", "KS13", "KS14", "KS15", "KS16",
#                          "KS17", "KS18", "KS20", "KS21", "KS22", "KS25",
#                          "KS27", "KS28", "KS29", "LY01A", "LY01B", "LY02",
#                          "LY03", "MOOR1", "APP2", "NUFF1", "SEL1")

names(gr)

gr <- flx
#to run without excess species

#aovs and tukey
aov <- aov(QY_max~Species*Treatment,data=gr) #all sig
aov <- aov(Fv.Fm_L3~Species*Treatment,data=gr) #all sig
aov <- aov(Fv.Fm_L5~Species*Treatment,data=gr) #just sp and int
aov <- aov(Fv.Fm_L11~Species*Treatment,data=gr) #all sig
aov <- aov(NPQ_L3~Species*Treatment,data=gr) #sp treat sig
aov <- aov(NPQ_L5~Species*Treatment,data=gr) #all sig
aov <- aov(NPQ_L11~Species*Treatment,data=gr) #all sig
summary(aov)
tukey <- TukeyHSD(aov, conf.level=.95) # to look at interaction
tukey <- TukeyHSD(aov, "Treatment", conf.level=.95)
tukey 
#QYmax #LL- l minu and l mino, s poly and l mino, l minu l jap, 
#between l jap with itself, 
#HL l minu l mino
#FvFm_L3  HL - L jap L mino, L jap L turion, L minu L jap, between L jap with itself
#FvFm_L5 HL - L jap L mino, L jap L turion
#FvFm_L11 LL - L mino and L minu, L jap with itself, HL - L jap with L turi, S poly, L mino
#NPQ L3 LL - L minu and L jap, s poly and L jap. HL - L jap and L mino, L jap and L turi, L jap and L minu
#NPQ L5 LL - LL - L jap and L minu, L minu with itself, HL - HL - L jap with L mino, L turio and L minu
#NPQ L11 LL - HL - L jap and L turi, L jap and L minu, L jap and S poly

#change order species
gr$Species <- factor(gr$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))
gr$Treatment <- factor(gr$Treatment, levels = c("LL", "HL"))

par(mfrow = c(1, 1))
#tiff('FloChl_boxplots_Species_sig_raw.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('FloChl_boxplots_Species_sigraw.pdf', width=14, height=12)
par(mfrow = c(3, 3),  mar=c(2,4.5,2,2))
#plot boxplots by species and treatment
boxplot(Fv.Fm_L3 ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 1), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 0.94, "a", cex=2)
text(2, 0.94, "a", cex=2)
text(2, 0.98, "***", cex=2)
text(3, 0.94, "a", cex=2)
text(4, 0.94, "a", cex=2)
text(5, 0.94, "a", cex=2)
text(7, 0.94, "a", cex=2)
text(8, 0.94, "b", cex=2)
text(9, 0.94, "a", cex=2)
text(10, 0.94, "ab", cex=2)
text(11, 0.94, "ab", cex=2)
title("A.", cex.main=2, adj=0)
#FvFm_L3  HL - L jap L mino, L jap L turion, L minu L jap, between L jap with itself
stripchart(Fv.Fm_L3 ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Fv.Fm_L5 ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 1), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75) #to edit based on results
text(1, 0.94, "a", cex=2)
text(2, 0.94, "a", cex=2)
#text(2, 0.86, "***", cex=2)
text(3, 0.94, "a", cex=2)
text(4, 0.94, "a", cex=2)
text(5, 0.94, "a", cex=2)
text(7, 0.94, "a", cex=2)
text(8, 0.94, "b", cex=2)
text(9, 0.94, "a", cex=2)
text(10, 0.94, "ab", cex=2)
text(11, 0.94, "ab", cex=2)
#FvFm_L3  #FvFm_L5 HL - L jap L mino, L jap L turion
stripchart(Fv.Fm_L5 ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Fv.Fm_L11 ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 1), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75) #to edit based on results
text(1, 0.94, "a", cex=2)
text(2, 0.94, "ab", cex=2)
text(2, 0.96, "***", cex=2)
text(3, 0.94, "ab", cex=2)
text(4, 0.94, "b", cex=2)
text(5, 0.94, "ab", cex=2)
text(7, 0.94, "a", cex=2)
text(8, 0.94, "b", cex=2)
text(9, 0.94, "a", cex=2)
text(10, 0.94, "ab", cex=2)
text(11, 0.94, "a", cex=2)
##FvFm_L11 LL - L mino and L minu, L jap with itself, HL - L jap with L turi, S poly, L mino
stripchart(Fv.Fm_L11 ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(NPQ_L3 ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 3), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 2.9, "ab", cex=2)
text(2, 2.9, "a", cex=2)
#text(2, 0.88, "***", cex=2)
text(3, 2.9, "ab", cex=2)
text(4, 2.9, "b", cex=2)
text(5, 2.9, "b", cex=2)
text(7, 2.9, "b", cex=2)
text(8, 2.9, "a", cex=2)
text(9, 2.9, "b", cex=2)
text(10, 2.9, "b", cex=2)
text(11, 2.9, "ab", cex=2)
title("B. NPQ", cex.main=2, adj=0)
#NPQ L3 LL - L minu and L jap, s poly and L jap. HL - L jap and L mino, L jap and L turi, L jap and L minu
stripchart(NPQ_L3 ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(NPQ_L5 ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 3), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75) #to edit based on results
text(1, 2.8, "ab", cex=2)
text(2, 2.8, "a", cex=2)
text(4, 2.9, "***", cex=2)
text(3, 2.8, "ab", cex=2)
text(4, 2.8, "b", cex=2)
text(5, 2.8, "ab", cex=2)
text(7, 2.8, "b", cex=2)
text(8, 2.8, "a", cex=2)
text(9, 2.8, "b", cex=2)
text(10, 2.8, "b", cex=2)
text(11, 2.8, "ab", cex=2)
#NPQ L5 LL - LL - L jap and L minu, L minu with itself, HL - HL - L jap with L mino, L turio and L minu
stripchart(NPQ_L5 ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(NPQ_L11 ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 5), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75) #to edit based on results
text(1, 4.84, "a", cex=2)
text(2, 4.84, "a", cex=2)
#text(2, 0.86, "***", cex=2)
text(3, 4.84, "a", cex=2)
text(4, 4.84, "a", cex=2)
text(5, 4.84, "a", cex=2)
text(7, 4.84, "ab", cex=2)
text(8, 4.84, "a", cex=2)
text(9, 4.84, "b", cex=2)
text(10, 4.84, "b", cex=2)
text(11, 4.84, "b", cex=2)
#NPQ L11 LL - HL - L jap and L turi, L jap and L minu, L jap and S poly
stripchart(NPQ_L11 ~ Species*Treatment, gr,  at =c(1,2,3,4,5, 7,8,9,10,11),           # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(QY_max ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0.6, 0.85), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
#text(0.9, 340, "HL", cex=2.75)
text(1, 0.84, "a", cex=2)
text(2, 0.84, "a", cex=2)
text(2, 0.85, "***", cex=2)
text(3, 0.84, "ab", cex=2)
text(4, 0.84, "b", cex=2)
text(5, 0.84, "b", cex=2)
text(7, 0.84, "a", cex=2)
text(8, 0.84, "ab", cex=2)
text(9, 0.84, "ab", cex=2)
text(10, 0.84, "b", cex=2)
text(11, 0.84, "ab", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(QY_max ~ Species*Treatment, gr, at =c(1,2,3,4,5, 7,8,9,10,11),            # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

dev.off()
dev.off()

#SPOT MEASUREMENTS OF INTEREST
#QY max plots
boxplot(QY_max~Accession*Rep,data=gr)
boxplot(QY_max~Accession,data=gr)
boxplot(QY_max~Species,data=gr)
boxplot(QY_max~Treatment,data=gr)
boxplot(QY_max~Accession*Treatment,data=gr)
boxplot(QY_max~Species*Treatment,data=gr)
boxplot(QY_max~Rep*Within_Rep,data=gr)
boxplot(QY_max~Rep,data=gr)
boxplot(QY_max~Within_Rep,data=gr) #not useful

#QY LOW
boxplot(Fv.Fm_L3~Accession*Rep,data=gr)
boxplot(Fv.Fm_L3~Accession,data=gr)
boxplot(Fv.Fm_L3~Species,data=gr)
boxplot(Fv.Fm_L3~Treatment,data=gr)
boxplot(Fv.Fm_L3~Accession*Treatment,data=gr)
boxplot(Fv.Fm_L3~Species*Treatment,data=gr)
boxplot(Fv.Fm_L3~Rep*Within_Rep,data=gr)

#QY MID
boxplot(Fv.Fm_L5~Accession*Rep,data=gr)
boxplot(Fv.Fm_L5~Accession,data=gr)
boxplot(Fv.Fm_L5~Species,data=gr)
boxplot(Fv.Fm_L5~Treatment,data=gr)
boxplot(Fv.Fm_L5~Accession*Treatment,data=gr)
boxplot(Fv.Fm_L5~Species*Treatment,data=gr)
boxplot(Fv.Fm_L5~Rep*Within_Rep,data=gr)

#QY HIGH
boxplot(Fv.Fm_L11~Accession*Rep,data=gr)
boxplot(Fv.Fm_L11~Accession,data=gr)
boxplot(Fv.Fm_L11~Species,data=gr)
boxplot(Fv.Fm_L11~Treatment,data=gr)
boxplot(Fv.Fm_L11~Accession*Treatment,data=gr)
boxplot(Fv.Fm_L11~Species*Treatment,data=gr)
boxplot(Fv.Fm_L11~Rep*Within_Rep,data=gr)

#NPQ LOW plots
boxplot(NPQ_L1~Accession*Rep,data=gr)
boxplot(NPQ_L1~Accession,data=gr)
boxplot(NPQ_L1~Species,data=gr)
boxplot(NPQ_L1~Treatment,data=gr)
boxplot(NPQ_L1~Accession*Treatment,data=gr)
boxplot(NPQ_L1~Species*Treatment,data=gr)
boxplot(NPQ_L1~Rep*Within_Rep,data=gr)

#NPQ LOW plots
boxplot(NPQ_L3~Accession*Rep,data=gr)
boxplot(NPQ_L3~Accession,data=gr)
boxplot(NPQ_L3~Species,data=gr)
boxplot(NPQ_L3~Treatment,data=gr)
boxplot(NPQ_L3~Accession*Treatment,data=gr)
boxplot(NPQ_L3~Species*Treatment,data=gr)
boxplot(NPQ_L3~Rep*Within_Rep,data=gr)

#NPQ MID plots
boxplot(NPQ_L5~Accession*Rep,data=gr)
boxplot(NPQ_L5~Accession,data=gr)
boxplot(NPQ_L5~Species,data=gr)
boxplot(NPQ_L5~Treatment,data=gr)
boxplot(NPQ_L5~Accession*Treatment,data=gr)
boxplot(NPQ_L5~Species*Treatment,data=gr)
boxplot(NPQ_L5~Rep*Within_Rep,data=gr)

#NPQ HIGH plots
boxplot(NPQ_L11~Accession*Rep,data=gr)
boxplot(NPQ_L11~Accession,data=gr)
boxplot(NPQ_L11~Species,data=gr)
boxplot(NPQ_L11~Treatment,data=gr)
boxplot(NPQ_L11~Accession*Treatment,data=gr)
boxplot(NPQ_L11~Species*Treatment,data=gr)
boxplot(NPQ_L11~Rep*Within_Rep,data=gr)

#which data to cut when 31 obs? which reps?
boxplot(NPQ_L11~Rep,data=gr) #hl1 low
boxplot(fqfm_L1~Rep,data=gr) #ll1 low
#only 3 reps of hl so cant cut one out
#hl1 earliest exp, not more prepared for hl,
#ll1 no different to others so doesn't naturally have higher
#npq when coming into lab

#24 obs
boxplot(NPQ_L11~Rep,data=gr) #hl1 low
boxplot(fqfm_L1~Rep,data=gr)
boxplot(NPQ_L3~Rep,data=gr) #ll4 high
boxplot(NPQ_L5~Rep,data=gr)

library(ggplot2)
bp <- ggplot(gr, aes(x=NPQ_L11, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
las=1
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=NPQ_L3, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=NPQ_L5, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)


bp <- ggplot(gr, aes(x=fqfm_L3, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=fqfm_L1, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
bp
bp + facet_grid(Treatment ~ .)

bp <- ggplot(gr, aes(x=QY_max, y=Accession, group=Accession)) + 
  geom_boxplot(aes(fill=Accession))
  theme_classic()
bp
bp + facet_grid(Treatment ~ .)

#see relationships between photosynthetic variables
plot(NPQ_L1 ~ QY_max, gr)
plot(NPQ_L11 ~ QY_max, gr) #CAN HAVE HIGH QYMAX AND HIGH NPQ L11
plot(NPQ_L1 ~ NPQ_L5, gr)
plot(NPQ_L3 ~ NPQ_L5, gr) #INTERESTING GRAPH 2 RESPONSES
plot(NPQ_L3 ~ Fv.Fm_L3, gr)
plot(NPQ_L5 ~ Fv.Fm_L5, gr)

#affects of other variables 32 obs
aov <- aov(QY_max~Accession,data=gr) #sig accession
aov <- aov(QY_max~Treatment,data=gr) # sig treatment
aov <- aov(QY_max~Species,data=gr) # sig species
aov <- aov(QY_max~Accession*Treatment,data=gr) #all sig
aov <- aov(QY_max~Species*Treatment,data=gr) #all sig
aov <- aov(QY_max~EnvLight*Species,data=gr) #all sig
aov <- aov(QY_max~Rep,data=gr) #sig affect of rep
aov <- aov(QY_max~Plate,data=gr) #not sig
aov <- aov(QY_max~Treatment*Rep,data=gr) #sig affect of rep and treat
aov <- aov(QY_max~Within_Rep,data=gr) #not sig
summary(aov)

aov <- aov(fqfm_L1~Accession*Treatment,data=gr) #acession + treat
aov <- aov(NPQ_L1~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L11~Accession*Treatment,data=gr) #all sig
aov <- aov(NPQ_L11~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L3~Accession*Treatment,data=gr) #access and acctreat affect
aov <- aov(fqfm_L5~Accession*Treatment,data=gr) #all sig
aov <- aov(NPQ_L5~Accession*Treatment,data=gr) #all sig
summary(aov)

aov <- aov(fqfm_L3~EnvLight*Species,data=gr) #ENV NOT SIG
aov <- aov(fqfm_L5~EnvLight*Species,data=gr) #ENV NOT SIG
aov <- aov(fqfm_L11~EnvLight*Species,data=gr) #ENV NOT SIG
aov <- aov(NPQ_L5~EnvLight*Species,data=gr) #ENV NOT SIG
aov <- aov(NPQ_L3~EnvLight*Species,data=gr) #ENV NOT SIG
aov <- aov(NPQ_L11~EnvLight*Species,data=gr) #ENV NOT SIG, INT SIG
summary(aov)

aov <- aov(fqfm_L3~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L3~Species*Treatment,data=gr) #seperate treat and sp
aov <- aov(fqfm_L5~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L11~Species*Treatment,data=gr) #

aov <- aov(NPQ_L3~Accession*Treatment,data=gr) #all sig
aov <- aov(NPQ_L3~Species*Treatment,data=gr) #seperate treat and sp
aov <- aov(NPQ_L5~Accession*Treatment,data=gr) #all sig
aov <- aov(NPQ_L5~Species*Treatment,data=gr) #all sig
aov <- aov(NPQ_L11~Species*Treatment,data=gr) #all sig
summary(aov)
tukey <- TukeyHSD(aov, "Species", conf.level=.95) # npq l3
tukey <- TukeyHSD(aov, "Treatment", conf.level=.95) #ll - hl 0.16 npq l3
tukey <- TukeyHSD(aov, "Accession", conf.level=.95)
tukey

#look like npq decline in hl not inc
boxplot(NPQ_L3~Species*Treatment,data=gr) #l minor higher npq l3 ll,
#l turion lowest npq l3 hl
boxplot(NPQ_L5~Species*Treatment,data=gr)

# npq increased in hl l.minor comparable other sp
boxplot(NPQ_L11~Species*Treatment,data=gr)

#maintainance of phi psii
# phi psii increaed in l turion hl
boxplot(fqfm_L11~Species*Treatment,data=gr) # l turion, s poly higher hl

boxplot(fqfm_L3~Species*Treatment,data=gr) # l turion, s poly higher hl
boxplot(fqfm_L5~Species*Treatment,data=gr) # l turion s poly higher hl

boxplot(fqfm_L3~Treatment,data=gr) # l turion, s poly higher hl
boxplot(fqfm_L5~Treatment,data=gr) # l turion s poly higher hl


aov <- aov(fqfm_L3~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L3~Species*Treatment,data=gr) #seperate treat and sp
aov <- aov(fqfm_L5~Accession*Treatment,data=gr) #seperate treat, less accession
aov <- aov(fqfm_L5~Species*Treatment,data=gr) #treat and sp seperate
summary(aov)

aov <- aov(NPQ_L11~Species*Treatment,data=gr) #all sig
aov <- aov(fqfm_L11~Species*Treatment,data=gr) #species sig
summary(aov)

t.test(fqfm_L1 ~ Treatment, data = gr) #higher ll, sig
t.test(fqfm_L3 ~ Treatment, data = gr) #higher ll, sig
t.test(fqfm_L5 ~ Treatment, data = gr) #higher hl, sig
t.test(fqfm_L11 ~ Treatment, data = gr) #higher hl, sig

t.test(QY_max ~ Treatment, data = gr) #higher ll, sig
t.test(NPQ_L11 ~ Treatment, data = gr) #higher hl, sig


aov <- aov(NPQ_L3 ~ Species*Treatment, data = gr) #higher hl, sig
aov <- aov(NPQ_L5 ~ Species*Treatment, data = gr) #higher hl, sig
aov <- aov(fqfm_L3 ~ Accession*Treatment, data = gr) #accession *, treat ***
aov <- aov(fqfm_L5 ~ Accession*Treatment, data = gr) #accession *, treat ***
summary(aov)

aov <- aov(NPQ_L3 ~ Species, data = gr) #sig
aov <- aov(NPQ_L5 ~ Species, data = gr) #sig
aov <- aov(NPQ_L5 ~ Species*Treatment, data = gr) #sig
aov <- aov(fqfm_L3 ~ Species*Treatment, data = gr) #sig
aov <- aov(fqfm_L5 ~ Species*Treatment, data = gr) #sig
aov <- aov(fqfm_L3 ~ Accession*Treatment, data = gr) #accession *, treat ***
aov <- aov(fqfm_L5 ~ Accession*Treatment, data = gr) #accession *, treat ***
tukey <- TukeyHSD(aov, "Species", conf.level=.95)
tukey <- TukeyHSD(aov, "Treatment", conf.level=.95)
tukey <- TukeyHSD(aov, "Accession", conf.level=.95)
tukey


boxplot(fqfm_L11~Species,data=gr)
boxplot(QY_max~Treatment,data=gr)

boxplot(NPQ_L3~Species,data=gr)
boxplot(NPQ_L5~Species,data=gr)
boxplot(NPQ_L5~Species*Treatment,data=gr)

aov <- aov(fqfm_L11 ~ Species, data = gr) #sig
tukey <- TukeyHSD(aov, "Species", conf.level=.95)
tukey

#24 obs
aov <- aov(QY_max~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L1~Accession*Treatment,data=gr) #all sig
aov <- aov(NPQ_L1~Accession*Treatment,data=gr) #treat, access*treat int sig
aov <- aov(fqfm_L11~Accession*Treatment,data=gr) #access sig
aov <- aov(NPQ_L11~Accession*Treatment,data=gr) #non sig
aov <- aov(fqfm_L3~Accession*Treatment,data=gr) #all sig, esp treat
aov <- aov(NPQ_L3~Accession*Treatment,data=gr) #all sig
aov <- aov(fqfm_L5~Accession*Treatment,data=gr) # treat sig, bit of access affect
aov <- aov(NPQ_L5~Accession*Treatment,data=gr) #all sig
summary(aov)

#new t tests for 24 obs
t.test(fqfm_L1 ~ Treatment, data = gr) #higher ll, sig
t.test(fqfm_L3 ~ Treatment, data = gr) #higher ll, sig
t.test(fqfm_L5 ~ Treatment, data = gr) #higher hl, sig
t.test(fqfm_L11 ~ Treatment, data = gr) #not sig
#most accessions had higher photosynthetic efficiency in
#low light when given low light, more had higher photo-
#synthhetic efficiency when given hl - flexible
#no affect at highest light - need to acclimate

t.test(QY_max ~ Treatment, data = gr) #higher ll, sig
t.test(NPQ_L1 ~ Treatment, data = gr) #higher hl, sig
t.test(NPQ_L3 ~ Treatment, data = gr) #higher ll, sig
t.test(NPQ_L5 ~ Treatment, data = gr) #higher ll, sig
t.test(NPQ_L11 ~ Treatment, data = gr) #higher hl, sig
#highest qy max in ll - preference ll (but access affect)
#higher npq in ll treated samples for mid values
#higher npq at lowest light and max light in hl

par(mfrow = c(1, 4))
boxplot(NPQ_L1 ~ Treatment, data = gr, ylim = c(0, 5)) #higher hl
boxplot(NPQ_L3 ~ Treatment, data = gr, ylim = c(0, 5)) #higher ll
boxplot(NPQ_L5 ~ Treatment, data = gr, ylim = c(0, 5)) #higher ll
boxplot(NPQ_L11 ~ Treatment, data = gr, ylim = c(0, 5)) #higher hl

#plot npq and fqfm variables
#install.packages("vioplot")
library("vioplot")

#practice jitter
par(mfrow = c(2, 4), mar = c(2.1, 4.1, 4.1, 2.1))
plot(NPQ_L1 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,5), xaxt="n", xlab = "", ylab = "NPQ_L1", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(NPQ_L1 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(NPQ_L3 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,5), xaxt="n", xlab = "", ylab = "NPQ_L3", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(NPQ_L3 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(NPQ_L5 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,5), xaxt="n", xlab = "", ylab = "NPQ_L5", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(NPQ_L5 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(NPQ_L11 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,5), xaxt="n", xlab = "", ylab = "NPQ_L11", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(NPQ_L11 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(QY_max ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,0.83), xaxt="n", xlab = "", ylab = "Max Quantum yield (Fv'Fm')", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(QY_max ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(fqfm_L3 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,0.83), xaxt="n", xlab = "", ylab = "Fqfm_L3", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(fqfm_L3 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(fqfm_L5 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,0.83), xaxt="n", xlab = "", ylab = "Fqfm_L5", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(fqfm_L5 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)
plot(fqfm_L11 ~ as.factor(Treatment), col=Treatment, gr, ylim = c(0,0.83), xaxt="n", xlab = "", ylab = "Fqfm_L11", whisklty = 1, staplelty = 0, boxwex = 0.6, outline = F)
stripchart(fqfm_L11 ~ as.factor(Treatment), gr, col = "#00000055", method = "jitter", add = T, vertical = T, pch=1)

#basic plots in matrix
par(mfrow = c(2, 4))
boxplot(NPQ_L1 ~ Treatment, data = gr, ylim = c(0, 5)) #higher hl
boxplot(NPQ_L3 ~ Treatment, data = gr, ylim = c(0, 5)) #higher ll
boxplot(NPQ_L5 ~ Treatment, data = gr, ylim = c(0, 5)) #higher ll
boxplot(NPQ_L11 ~ Treatment, data = gr, ylim = c(0, 5)) #higher hl
boxplot(QY_max ~ Treatment, data = gr, ylim = c(0, 0.83)) #higher ll
boxplot(fqfm_L3 ~ Treatment, data = gr, ylim = c(0, 0.83)) #higher ll
boxplot(fqfm_L5 ~ Treatment, data = gr, ylim = c(0, 0.83)) #higher ll
boxplot(fqfm_L11 ~ Treatment, data = gr, ylim = c(0, 0.83)) #higher hl

#find % difference between NPQ see which access has biggest inc
 gr2 <- gr %>% mutate(NPQ_diff_1 = ((NPQ_L2 - NPQ_L1)/NPQ_L2)*100) %>% 
   mutate(NPQ_diff_2 = ((NPQ_L3 - NPQ_L2)/NPQ_L3)*100) %>%
   mutate(NPQ_diff_3 = ((NPQ_L4 - NPQ_L3)/NPQ_L4)*100) %>%
   mutate(NPQ_diff_4 = ((NPQ_L5 - NPQ_L4)/NPQ_L5)*100) %>%
   mutate(NPQ_diff_5 = ((NPQ_L6 - NPQ_L5)/NPQ_L6)*100) %>%
   mutate(NPQ_diff_6 = ((NPQ_L7 - NPQ_L6)/NPQ_L7)*100) %>%
   mutate(NPQ_diff_7 = ((NPQ_L8 - NPQ_L7)/NPQ_L8)*100) %>%
   mutate(NPQ_diff_8 = ((NPQ_L9 - NPQ_L8)/NPQ_L9)*100) %>%
   mutate(NPQ_diff_9 = ((NPQ_L10 - NPQ_L9)/NPQ_L10)*100) %>%
   mutate(NPQ_diff_10 = ((NPQ_L11 - NPQ_L10)/NPQ_L11)*100)

 
 bp <- ggplot(gr2, aes(x=NPQ_diff_2, y=Accession, group=Accession)) + 
   geom_boxplot(aes(fill=Accession))
 bp
 bp + facet_grid(Treatment ~ .)
 
 bp <- ggplot(gr2, aes(x=NPQ_diff_3, y=Accession, group=Accession)) + 
   geom_boxplot(aes(fill=Accession))
 bp
 bp + facet_grid(Treatment ~ .)
 
 bp <- ggplot(gr2, aes(x=NPQ_diff_5, y=Accession, group=Accession)) + 
   geom_boxplot(aes(fill=Accession))
 bp
 bp + facet_grid(Treatment ~ .)
 
 bp <- ggplot(gr2, aes(x=NPQ_diff_10, y=Accession, group=Accession)) + 
   geom_boxplot(aes(fill=Accession))
 bp
 bp + facet_grid(Treatment ~ .)
 
 t.test(NPQ_diff_5 ~ Treatment, data = gr2) #higher hl, sig
 t.test(NPQ_diff_2 ~ Treatment, data = gr2) #higher hl, sig
 t.test(NPQ_diff_3 ~ Treatment, data = gr2) #higher ll, sig
 t.test(NPQ_diff_10 ~ Treatment, data = gr2) #higher hl, sig
 t.test(NPQ_diff_1 ~ Treatment, data = gr2) #higher hl, sig
 
gr %>% arrange(desc(QY_max)) %>% select(Accession, QY_max) %>% top_n(24)

gr2 %>% arrange(desc(NPQ_diff_10)) %>% select(Accession, NPQ_diff_10, Treatment) %>% top_n(24)
gr2 %>% arrange(desc(NPQ_diff_2)) %>% select(Accession, NPQ_diff_2, Treatment) %>% top_n(24)
gr2 %>% arrange(desc(NPQ_diff_5)) %>% select(Accession, NPQ_diff_5, Treatment) %>% top_n(24)

tukey <- TukeyHSD(aov, "Rep", conf.level=.95) #hl2 diff from shl1, shl2
#hl1 diff from shl1 shl2
tukey <- TukeyHSD(aov, "Treatment", conf.level=.95) 
tukey <- TukeyHSD(aov, "Accession", conf.level=.95) 
str(tukey)
tukey

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

range(gr$QY_max)
#between 0.08 0.83
which(gr$QY_max < 0.6)
#47, 58, 82, 176, 186, 190, 251, 312 to remove outliers ly02, ks15,
#ks15, 13, 27 and ly01b, 16, ly01a
gr$Accession
gr %>% arrange(desc(QY_max)) %>% select(Accession, QY_max, Treatment) %>% top_n(36)
gr %>% arrange(QY_max) %>% select(Accession, QY_max) %>% top_n(36)

#remove outliers to mod
#gr_sp <- gr[-c(64:69), ]

range(gr$NPQ_L11)
which(gr$NPQ_L11 > 5)
gr %>% arrange(desc(NPQ_L11)) %>% select(Accession, NPQ_L11) %>% top_n(10)

#SUMMARISING BY ACCESSION AND TREATMENT
gr <- read.csv("Fluorcam+fqfm.csv") #including new variables
sum <- gr %>% group_by(Accession, Treatment) %>% summarise_all(mean)
#now have summarised RGR rate per access and treatment and other cols
sum$QY_max
sum$Fv.Fm_L1
sum$Fv.Fm_L3
sum$Fv.Fm_L11
sum$NPQ_L1
sum$NPQ_L3
sum$NPQ_L11
sum$NPQ_L5
sum$Fv.Fm_L5


sum %>% arrange(desc(NPQ_L11)) %>% select(Accession, Treatment, NPQ_L11) %>% top_n(10)

#more detailed summary
Summary <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(QY_max_mean = mean(QY_max), QY_max_stdev = sd(QY_max), QY_max_n= n(), QY_max_maximum = max(QY_max))
Summary

Summary1 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(Fv.Fm_L3_mean = mean(Fv.Fm_L3), Fv.Fm_L3_stdev = sd(Fv.Fm_L3), Fv.Fm_L3_n= n(), Fv.Fm_L3_maximum = max(Fv.Fm_L3))
Summary1

Summary2 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(Fv.Fm_L11_mean = mean(Fv.Fm_L11), Fv.Fm_L11_stdev = sd(Fv.Fm_L11), Fv.Fm_L11_n= n(), Fv.Fm_L11_maximum = max(Fv.Fm_L11))
Summary2

Summary3 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(NPQ_L1_mean = mean(NPQ_L1), NPQ_L1_stdev = sd(NPQ_L1), NPQ_L1_n= n(), NPQ_L1_maximum = max(NPQ_L1))
Summary3

Summary4 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(NPQ_L3_mean = mean(NPQ_L3), NPQ_L3_stdev = sd(NPQ_L3), NPQ_L3_n= n(), NPQ_L3_maximum = max(NPQ_L3))
Summary4

Summary5 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(NPQ_L11_mean = mean(NPQ_L11), NPQ_L11_stdev = sd(NPQ_L11), NPQ_L11_n= n(), NPQ_L11_maximum = max(NPQ_L11))
Summary5

Summary6 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(NPQ_L5_mean = mean(NPQ_L5), NPQ_L5_stdev = sd(NPQ_L5), NPQ_L5_n= n(), NPQ_L5_maximum = max(NPQ_L5))
Summary6

Summary7 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(Fv.Fm_L5_mean = mean(Fv.Fm_L5), Fv.Fm_L5_stdev = sd(Fv.Fm_L5), Fv.Fm_L5_n= n(), Fv.Fm_L5_maximum = max(Fv.Fm_L5))
Summary7

Summary8 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(fqfm_L1_mean = mean(fqfm_L1), fqfm_L1_stdev = sd(fqfm_L1), fqfm_L1_n= n(), fqfm_L1_maximum = max(fqfm_L1))
Summary8

Summary9 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(fqfm_L3_mean = mean(fqfm_L3), fqfm_L3_stdev = sd(fqfm_L3), fqfm_L3_n= n(), fqfm_L3_maximum = max(fqfm_L3))
Summary9

Summary10 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(fqfm_L5_mean = mean(fqfm_L5), fqfm_L5_stdev = sd(fqfm_L5), fqfm_L5_n= n(), fqfm_L5_maximum = max(fqfm_L5))
Summary10

Summary11 <- gr %>%
  group_by(Treatment,Accession) %>%
  summarise(fqfm_L11_mean = mean(fqfm_L11), fqfm_L11_stdev = sd(fqfm_L11), fqfm_L11_n= n(), fqfm_L11_maximum = max(fqfm_L11))
Summary11

#combine 5 summaries manually

Summ1 <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5, Summary6, Summary7, Summary8, Summary9, Summary10, Summary11)

#write.csv(Summ, "Fluo_Summary_new.csv")
write.csv(Summ1, "Fluo_Summary_REMsp_REMoutliers.csv") #includes fqfm col summaries

#EXERCISE TO NAME COLS, SPLIT COLUMNS BY ROW NUMBER AND PUT BACK TOGETHER AS
#SEPERATE COLUMNS -- WORKS --
try <- read.csv("Fluo_Summary_REMsp_REMoutliers.csv")
#names(try)
#length(try)
#give columns names
#colnames(try) <- c("X", "Treatment", "Accession",
#                        "fqfm_L1_mean", "fqfm_L1_stddev", "fqfm_L1_n", "fqfm_L1_max", 
#                        "Treatment1", "Accession1", "fqfm_L3_mean", "fqfm_L3_stddev", 
#                        "fqfm_L3_n", "fqfm_L3_max", "Treatment2", "Accession2", "fqfm_L5_mean", 
#                        "fqfm_L5_stddev", "fqfm_L5_n", "fqfm_L5_max", "Treatment3", "Accession3", 
#                        "fqfm_L11_mean", "fqfm_L11_stddev", "fqfm_L11_n", "fqfm_L11_max")
                        

try1 <- split(try,cumsum(1:nrow(try)%in%25)) #made 2 lists 0 and 1
str(try1)
#call the two lists
try1$`0`
try1$`1`
#check correct number of hl and ll and split properly
try1$`0`$Treatment
try1$`1`$Treatment
backtogether <- cbind(try1$`0`, try1$`1`) #need to be the same no args
#puts all hl cols together first then all ll cols after
#need to rename colnames
names(backtogether)
#remove extra columns
backtogether2 = subset(backtogether, select = -c(1,8,9,14,15,20,21,26,27,32,33,38,39,44,45,50,51,56,57,62,63,68,69,74,81,82,87,88,93,94,99,100,105,106,111,112,117,118,123,124,129,130,135,136,141,142))
write.csv(backtogether2, "Fluo_Summary1_REMsp_REMout.csv")


#rename col names for either hl or ll

Summ <- read.csv("Fluo_Summary_sepHL_LL.csv")

names(Summ)
range(Summ$mean.LL.QYmax)
#[1] 0.6212500 0.7916667
range(Summ$mean.HL.QYmax)
#[1] 0.5583333 0.7771429

#find top 10 performing accession RGR values
#6A top of HL, LY01B top of LL
#12 bottom of HL and LL - slowest growing

gr %>% 
  ggplot(aes(QY_max, Accession, col = Treatment)) + 
  geom_point(alpha = 0.8) +
  facet_wrap(~Treatment)


t.test(QY_max ~ Treatment, data = gr) # sig
t.test(NPQ_L11 ~ Treatment, data = gr) # sig #highest in HL
t.test(NPQ_L1 ~ Treatment, data = gr) # sig #highest in HL
t.test(NPQ_L3 ~ Treatment, data = gr) # sig #highest in LL
t.test(NPQ_L5 ~ Treatment, data = gr) # sig #highest in LL


#boxplot faceted by treatment
QYmax_boxplot <- ggplot(gr,aes(x = reorder(Accession,QY_max, FUN = median),y=QY_max)) +
  
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
  labs(title =expression("(a) QYmax"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(QY_max))+
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
QYmax_boxplot

#L3 FVFM LL
Fv.Fm_L3_boxplot <- ggplot(gr,aes(x = reorder(Accession,Fv.Fm_L3, FUN = median),y=Fv.Fm_L3)) +
  
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
  labs(title =expression("(a) Fv.Fm_L3"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Fv.Fm_L3))+
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
Fv.Fm_L3_boxplot

#FVFVM L5 HL
Fv.Fm_L5_boxplot <- ggplot(gr,aes(x = reorder(Accession,Fv.Fm_L5, FUN = median),y=Fv.Fm_L5)) +
  
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
  labs(title =expression("(a) Fv.Fm_L5"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(Fv.Fm_L5))+
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
Fv.Fm_L5_boxplot

#NPQ LL
NPQ_L1_boxplot <- ggplot(gr,aes(x = reorder(Accession,NPQ_L1, FUN = median),y=NPQ_L1)) +
  
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
  labs(title =expression("(a) NPQ_L1"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(NPQ_L1))+
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
NPQ_L1_boxplot

#NPQ LL
NPQ_L3_boxplot <- ggplot(gr,aes(x = reorder(Accession,NPQ_L3, FUN = median),y=NPQ_L3)) +
  
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
  labs(title =expression("(a) NPQ_L3"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(NPQ_L3))+
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
NPQ_L3_boxplot

#NPQ HL
NPQ_L5_boxplot <- ggplot(gr,aes(x = reorder(Accession,NPQ_L5, FUN = median),y=NPQ_L5)) +
  
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
  labs(title =expression("(a) NPQ_L5"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(NPQ_L5))+
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
NPQ_L5_boxplot

#NPQ ExtremeHL
NPQ_L11_boxplot <- ggplot(gr,aes(x = reorder(Accession,NPQ_L11, FUN = median),y=NPQ_L11)) +
  
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
  labs(title =expression("(a) NPQ_L11"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(NPQ_L11))+
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
NPQ_L11_boxplot

#plot by rep
#boxplot faceted by treatment
NPQ_L11rep_boxplot <- ggplot(gr,aes(x = reorder(Accession,NPQ_L11, FUN = median),y=NPQ_L11)) +
  
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
  labs(title =expression("(a) NPQ_L11"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(NPQ_L11))+
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
NPQ_L11rep_boxplot

#read in new file with grouped HL LL means and max
hl_llmeans <- read.csv("Fluo_Summary_sepHL_LL.csv")
hl_llmeans <- read.csv("Fluo_Summary1_REMsp_REMout_HL_LLrename.csv")
#new version without species, outliers
#wont work with code until use renames

names(hl_llmeans)

#take hl away from ll calculated paramters NPQ and fvfm main parameters
hl_llmeans <- hl_llmeans %>% mutate(meanQYmax_diff = (QY_max_HL_mean - QY_max_LL_mean)) %>% 
  mutate(maxQYmax_diff = (QY_max_HL_maximum - QY_max_LL_maximum)) %>%
  mutate(meanNPQ_L11_diff = (NPQ_L11_HL_mean - NPQ_L11_LL_mean)) %>%
  mutate(meanNPQ_L5_diff = (NPQ_L5_HL_mean - NPQ_L5_LL_mean)) %>%
  mutate(meanNPQ_L3_diff = (NPQ_L3_HL_mean - NPQ_L3_LL_mean)) %>%
  mutate(meanNPQ_L1_diff = (NPQ_L1_HL_mean - NPQ_L1_LL_mean)) %>%
  mutate(meanfqfm_L11_diff = (fqfm_L11_HL_mean - fqfm_L11_LL_mean)) %>%
  mutate(meanfqfm_L5_diff = (fqfm_L5_HL_mean - fqfm_L5_LL_mean)) %>%
  mutate(meanfqfm_L3_diff = (fqfm_L3_HL_mean - fqfm_L3_LL_mean)) %>%
  mutate(meanfqfm_L1_diff = (fqfm_L1_HL_mean - fqfm_L1_LL_mean))

#plot difference as scatterplot
plot(meanQYmax_diff ~ maxQYmax_diff, hl_llmeans) #not neccessrily a trade off
plot(meanQYmax_diff ~ meanNPQ_L11_diff, hl_llmeans)
plot(meanNPQ_L3_diff ~ meanNPQ_L5_diff, hl_llmeans)
plot(meanNPQ_L3_diff ~ meanNPQ_L11_diff, hl_llmeans)
plot(meanNPQ_L5_diff ~ meanNPQ_L11_diff, hl_llmeans)
plot(meanNPQ_L1_diff ~ meanNPQ_L11_diff, hl_llmeans)
plot(meanNPQ_L1_diff ~ meanQYmax_diff, hl_llmeans) #not neccessrily trade off
plot(meanfqfm_L3_diff ~ meanfqfm_L5_diff, hl_llmeans)
plot(meanfqfm_L3_diff ~ meanfqfm_L11_diff, hl_llmeans)
plot(meanfqfm_L5_diff ~ meanfqfm_L11_diff, hl_llmeans)
plot(meanfqfm_L1_diff ~ meanfqfm_L11_diff, hl_llmeans)
plot(meanQYmax_diff ~ meanNPQ_L11_diff, hl_llmeans)

#want to know if qy max diff or npq max diff corr other factors

range(hl_llmeans$meanQYmax_diff)
#[1] -0.040  0.212
range(hl_llmeans$maxQYmax_diff)
#[1]  -0.06  0.17
range(hl_llmeans$meanNPQ_L11_diff)
#[1] -1.885  2.465
range(hl_llmeans$meanNPQ_L3_diff)
# [1] -0.6078788  1.0850000
range(hl_llmeans$meanNPQ_L5_diff)
# [1] -0.4915152  0.7371591
#some minus some positive in QYmax and NPQ diffs

#group data into + or - RGR diff to see which direction performed better
hl_llmeans$meanQYmax_diff > 0
hl_llmeans$maxQYmax_diff > 0
hl_llmeans$meanNPQ_L11_diff > 0

which(hl_llmeans$meanQYmax_diff > 0) #23 accessions positive
which(hl_llmeans$meanQYmax_diff < 0) #8 accessions neg
which(hl_llmeans$maxQYmax_diff > 0) #23 accession positive
which(hl_llmeans$maxQYmax_diff < 0) #3 accessions neg
which(hl_llmeans$meanNPQ_L11_diff > 0) #4 accession positive
which(hl_llmeans$meanNPQ_L11_diff < 0) #27 accessions neg

hl_llmeans$Accession
hl_llmeans %>% arrange(desc(meanQYmax_diff)) %>% select(Accession) %>% top_n(31)
#shows same thing - which had highest values for difference
hl_llmeans$Accession
hl_llmeans %>% arrange(desc(maxQYmax_diff)) %>% select(Accession) %>% top_n(31)
#shows same thing - which had highest values for difference
hl_llmeans$Accession
hl_llmeans %>% arrange(desc(meanNPQ_L11_diff)) %>% select(Accession) %>% top_n(31)
#shows same thing - which had highest values for difference
hl_llmeans %>% arrange(desc(meanNPQ_L5_diff)) %>% select(Accession) %>% top_n(31)
#shows same thing - which had highest values for difference
hl_llmeans %>% arrange(desc(meanNPQ_L3_diff)) %>% select(Accession) %>% top_n(31)
#shows same thing - which had highest values for difference
hl_llmeans %>% arrange(desc(meanfqfm_L11_diff)) %>% select(Accession) %>% top_n(31)
#shows same thing - which had highest values for difference

#as says online to do % proportion increase
hl_llmeans <- hl_llmeans %>% mutate(prop_NPQ_L3 = (meanNPQ_L3_diff/NPQ_L3_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_NPQ_L3)) %>% select(Accession, prop_NPQ_L3) %>% top_n(31)
hl_llmeans <- hl_llmeans %>% mutate(prop_NPQ_L5 = (meanNPQ_L5_diff/NPQ_L5_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_NPQ_L5)) %>% select(Accession, prop_NPQ_L5) %>% top_n(31)
hl_llmeans <- hl_llmeans %>% mutate(prop_NPQ_L11 = (meanNPQ_L11_diff/NPQ_L11_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_NPQ_L11)) %>% select(Accession, prop_NPQ_L11) %>% top_n(31)

hl_llmeans <- hl_llmeans %>% mutate(prop_fqfm_L3 = (meanfqfm_L3_diff/fqfm_L3_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_fqfm_L3)) %>% select(Accession, prop_fqfm_L3) %>% top_n(31)
hl_llmeans <- hl_llmeans %>% mutate(prop_fqfm_L5 = (meanfqfm_L5_diff/fqfm_L5_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_fqfm_L5)) %>% select(Accession, prop_fqfm_L5) %>% top_n(31)
hl_llmeans <- hl_llmeans %>% mutate(prop_fqfm_L11 = (meanfqfm_L11_diff/fqfm_L11_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_fqfm_L11)) %>% select(Accession, prop_fqfm_L11) %>% top_n(31)

hl_llmeans <- hl_llmeans %>% mutate(prop_QYmax = (meanQYmax_diff/QY_max_HL_mean)*100)
hl_llmeans %>% arrange(desc(prop_QYmax)) %>% select(Accession, prop_QYmax) %>% top_n(31)

write.csv(hl_llmeans, "Fluo_means_diffs_props.csv")

#effect of treatment 
ggplot(hl_llmeans, aes(y=meanQYmax_diff, x=reorder(Accession, meanQYmax_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)
#negative left = higher QYmax in HL. mid = equal. right = positive, higher QYmax in LL
#MOOR1 most highly affected by the HL
ggplot(hl_llmeans, aes(y=maxQYmax_diff, x=reorder(Accession, maxQYmax_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)

ggplot(hl_llmeans, aes(y=meanNPQ_L11_diff, x=reorder(Accession, meanNPQ_L11_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)
#NPQ at highest light: left NPQ highest in HL, right - severerely negatively
#NPQ affected by HL - already high in APP2 and SEL1, gone down in HL,
#already prepared for high light as taken from env?

ggplot(hl_llmeans, aes(y=meanNPQ_L3_diff, x=reorder(Accession, meanNPQ_L3_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)

ggplot(hl_llmeans, aes(y=meanNPQ_L5_diff, x=reorder(Accession, meanNPQ_L5_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)

#difference in photosynthetic efficiency in HL vs LL
ggplot(hl_llmeans, aes(y=meanfqfm_L11_diff, x=reorder(Accession, meanfqfm_L11_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)
#left higher photosynthetic efficiency in HL, right more negatively affected
#by HL treatment really high light here, still any efficiency important?

ggplot(hl_llmeans, aes(y=meanfqfm_L3_diff, x=reorder(Accession, meanfqfm_L3_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)


ggplot(hl_llmeans, aes(y=meanfqfm_L5_diff, x=reorder(Accession, meanfqfm_L5_diff))) +
  geom_col(position='dodge', color='black')
geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position = position_dodge(0.9), width = 0.25)
#left higher in HL, right = more affected by HL, mid = hardly any diff

#BAR PLOTS WITH STD ERROR

#work out std error using base r by treatment accession combo
#QY max
aggregate(QY_max ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(QY_max),
    sd=sd(QY_max)
  ) %>%
  mutate( se=sd/sqrt(n))

#QY max plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#FvFm_L3
#work out std error using base r by treatment accession combo
aggregate(Fv.Fm_L3 ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Fv.Fm_L3),
    sd=sd(Fv.Fm_L3)
  ) %>%
  mutate( se=sd/sqrt(n))

#FVFM L3
#plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#FVFML5
#work out std error using base r by treatment accession combo
aggregate(Fv.Fm_L5 ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#how much chl a and b lost due to chlorosis, biggest decline = worst affected

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(Fv.Fm_L5),
    sd=sd(Fv.Fm_L5)
  ) %>%
  mutate( se=sd/sqrt(n))

#FVFM L5
#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#NPQ L1
#work out std error using base r by treatment accession combo
aggregate(NPQ_L1 ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#how much chl a and b lost due to chlorosis, biggest decline = worst affected

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(NPQ_L1),
    sd=sd(NPQ_L1)
  ) %>%
  mutate( se=sd/sqrt(n))

#NPQ_L1
#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#NPQ L3
#work out std error using base r by treatment accession combo
aggregate(NPQ_L3 ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#how much chl a and b lost due to chlorosis, biggest decline = worst affected

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(NPQ_L3),
    sd=sd(NPQ_L3)
  ) %>%
  mutate( se=sd/sqrt(n))

#NPQ_L3
#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#NPQ L5
#work out std error using base r by treatment accession combo
aggregate(NPQ_L5 ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(NPQ_L5),
    sd=sd(NPQ_L5)
  ) %>%
  mutate( se=sd/sqrt(n))

#NPQ_L5
#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#NPQ L11
#work out std error using base r by treatment accession combo
aggregate(NPQ_L11 ~ Accession + Treatment, data = gr, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))

#how much chl a and b lost due to chlorosis, biggest decline = worst affected

#calc se for ggplot std error
my_sum <- gr %>%
  group_by(Accession, Treatment) %>%
  summarise( 
    n=n(),
    mean=mean(NPQ_L11),
    sd=sd(NPQ_L11)
  ) %>%
  mutate( se=sd/sqrt(n))

#NPQ_L11
#Car ration plots per HL and LL
ggplot(my_sum, aes(y=mean, x=Accession, fill=Treatment)) +
  geom_col(position='dodge', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position = position_dodge(0.9), width = 0.25)

#get fqfm values

flx <- read.csv("Fluorcam.csv")
flx <- read.csv("Fluorcam+fqfm_REMsum2022_spr2021_REMoutliers.csv") #rem obs

library(ggplot2)
library(dplyr)
flx$fqfm_L1 <- (flx$Fq_L1 / flx$Fm_L1) #divide and create new col
flx$fqfm_L2 <- (flx$Fq_L2 / flx$Fm_L2) #divide and create new col
flx$fqfm_L3 <- (flx$Fq_L3 / flx$Fm_L3) #divide and create new col
flx$fqfm_L4 <- (flx$Fq_L4 / flx$Fm_L4) #divide and create new col
flx$fqfm_L5 <- (flx$Fq_L5 / flx$Fm_L5) #divide and create new col
flx$fqfm_L6 <- (flx$Fq_L6 / flx$Fm_L6) #divide and create new col
flx$fqfm_L7 <- (flx$Fq_L7 / flx$Fm_L7) #divide and create new col
flx$fqfm_L8 <- (flx$Fq_L8 / flx$Fm_L8) #divide and create new col
flx$fqfm_L9 <- (flx$Fq_L9 / flx$Fm_L9) #divide and create new col
flx$fqfm_L10 <- (flx$Fq_L10 / flx$Fm_L10) #divide and create new col
flx$fqfm_L11 <- (flx$Fq_L11 / flx$Fm_L11) #divide and create new col

sapply(flx, class) #check data class, make sure your factors are factors and not integers  # used before

flx_ind <- which(names(flx) == "fqfm_L1")

flx$Species <- as.factor(flx$Species)
flx$Treatment <- as.factor(flx$Treatment)
flx$Accession <- as.factor(flx$Accession)

flx %>% ggplot(aes(Accession, fqfm_L1, col = Species)) + geom_boxplot() + facet_wrap(~Treatment) + coord_flip()

flx[flx$Treatment == "HL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, NPQ_L11), NPQ_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "High light (HL)", x = "Accession", y = "NPQ") + theme_minimal() + theme(text = element_text(size = 18))
flx[flx$Treatment == "LL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, NPQ_L11), NPQ_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "Low light (LL)", x = "Accession", y = "NPQ") + theme_minimal() + theme(text = element_text(size = 18))

flx[flx$Treatment == "HL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, fqfm_L11), fqfm_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "High light (HL)", x = "Accession", y = "fqfm") + theme_minimal() + theme(text = element_text(size = 18))
flx[flx$Treatment == "LL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, fqfm_L11), fqfm_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "Low light (LL)", x = "Accession", y = "fqfm") + theme_minimal() + theme(text = element_text(size = 18))

#lines through scatter plots according to species or treat
flx %>% ggplot(aes(fqfm_L11, NPQ_L11, col = Species)) + geom_point() + stat_smooth(method = "lm") +
  theme_bw() +
  theme_classic()

flx %>% ggplot(aes(fqfm_L11, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
    xlab("Photosynthetic efficiency (fqfm) max light") +
    ylab("NPQ max light") +
    theme_bw() +
  theme_classic()

flx %>% ggplot(aes(fqfm_L1, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
  xlab("Photosynthetic efficiency low light") +
  ylab("NPQ max light") +
  theme_bw() +
  theme_classic()

flx %>% ggplot(aes(QY_max, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
  xlab("Photosynthetic efficiency max light") +
  ylab("NPQ max light") +
  theme_bw() +
  theme_classic()

flx %>% ggplot(aes(fqfm_L1, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
theme_bw() +
  theme_classic()

flx %>% ggplot(aes(NPQ_L3, NPQ_L5, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
  xlab("NPQ (100 umol m-2 s-1)") +
  ylab("NPQ (350 umol m-2 s-1)") +
  theme_bw() +
  theme_classic()

flx %>% ggplot(aes(NPQ_L3, fqfm_L3, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
  xlab("NPQ (100 umol m-2 s-1)") +
  ylab("fqfm (100 umol m-2 s-1)") +
  theme_bw() +
  theme_classic()

flx %>% ggplot(aes(NPQ_L5, fqfm_L5, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
  xlab("NPQ (350 umol m-2 s-1)") +
  ylab("fqfm (350 umol m-2 s-1)") +
  theme_bw() +
  theme_classic()

library(ggpubr)
colors <- c("#ff0000", "#262828")
class(flx$Treatment)
levels(flx$Treatment) <- c("HL", "LL")
names(colors) <- levels(flx$Treatment)
#plot rate of change between fqfm L3 and L5
pdf('PSIIQY_LL_HLrateofchange_ind_dep_remslopes.pdf', width=6.5, height=6)
flx %>% ggplot(aes(fqfm_L3, fqfm_L5, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
  scale_color_manual(values=colors)+
  xlab("fqfm (130 umol m-2 s-1)") +
  ylab("fqfm (365 umol m-2 s-1)") +
  theme_bw() +
  theme_classic()#+
#  stat_cor(
#    aes(label = paste(..rr.label..)))+
      #  (label.y = c(0.5,0.55)) +
 # stat_regline_equation(mapping = NULL, label.y = c(0.6,0.65))
  #stat_cor(mapping=NULL)+
  #stat_regline_equation(mapping= NULL)
#  stat_regline_equation(
#    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
#    formula = formula) 
  #stat_regline_equation(mapping = NULL)
dev.off()
#y = -0.043 + 0.83 x 0.52 #HL
#y = 0.092 + 0.57 x r2 = 0.87 #LL
#HL started from higher values but slope less extreme

#check
fit <- lm(formula = fqfm_L3 ~ fqfm_L5 + Treatment, data=flx)
summary(fit)

reshl <- flx[flx$Treatment=="HL",]
resll <- flx[flx$Treatment=="LL",]
fit_hl <- lm(fqfm_L5 ~ fqfm_L3, data=reshl)
fit_ll <- lm(fqfm_L5 ~ fqfm_L3, data=resll)
fit_hl$coefficients
fit_ll$coefficients

plot(fqfm_L5 ~ fqfm_L3, data=flx, col=as.integer(Treatment), xlab="Drinks/week", 
     ylab="Party hours/week", pch=19, las=1)
abline(fit_hl$coefficients, col=1, lwd=2)
abline(fit_ll$coefficients, col=2, lwd=2)

as.integer(flx$Treatment) #1 = HL, 2 = LL
  
#plot rate of change between fqfm L1 and L2
flx %>% ggplot(aes(fqfm_L1, fqfm_L2, col = Treatment)) + geom_point() + stat_smooth(method = "lm") +
xlab("fqfm (0 umol m-2 s-1)") +
ylab("fqfm (20 umol m-2 s-1)") +
  theme_bw() +
  theme_classic()

flx %>% ggplot(aes(NPQ_L1, NPQ_L3, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(NPQ_L1, NPQ_L5, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(QY_max, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm")

#call species
my_spp <- unique(flx[, 144])

#summarised by accession originally but doesnt include select
flx2 <- flx %>% select(Accession, Species, Treatment, NPQ_L11, fqfm_L11) %>% group_by(Accession, Species, Treatment) %>% summarise(mean_npq = mean(NPQ_L11), mean_fqfm = mean(fqfm_L11))
flx2 <- flx %>% select(Accession, Treatment, NPQ_L11, fqfm_L11) %>% group_by(Accession, Treatment) %>% summarise(mean_npq = mean(NPQ_L11), mean_fqfm = mean(fqfm_L11))
flx2$Species <- my_spp 
flx2$Treatment <- my_tr
flx2 <- flx2[flx2$Species %in% c("L. minor", "L. minuta", "L. gibba", "L. turionifera", "S. polyrhiza"), ]
flx2 <- flx2[flx2$Treatment %in% c("LL", "HL"), ]

#plots line plots for ind species and summarises model as intercept and mean_fqfm
flx2 %>% ggplot(aes(mean_fqfm, mean_npq, col = flx$Treatment)) + geom_point() + stat_smooth(method = "lm")
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Treatment == "LL", ]))
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Treatment == "HL", ]))

flx2 %>% ggplot(aes(mean_fqfm, mean_npq, col = Species)) + geom_point() + stat_smooth(method = "lm")
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Species == "L. minuta", ]))
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Species == "L. minor", ]))
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Species == "L. gibba", ]))
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Species == "L. turionifera", ]))
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Species == "S. polyrhiza", ]))

flx2 %>% ggplot(aes(mean_fqfm, mean_npq, col = Treatment)) + geom_point() + stat_smooth(method = "lm") 
summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Treatment == "HL", ])) +
  summary(lm(mean_npq ~ mean_fqfm, flx2[flx2$Treatment == "LL", ])) +
  theme.bw()

write.csv(flx, "Fluorcam+fqfm.csv")

#AUTOMATE TUKEYS AND LINEAR MODELS
#install.packages("metR")

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

#moved species to start
flall <- read.csv("Fluorcam+fqfm.csv")
str(flall)

#manually moved species to start as automation fails if in middle col

library(dplyr)

#run this to save to correct place
#setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Fluorcam\\R")


#ANOVA AUTOMATION FROM MATT K
#starts from col no, #accession and treat in pos 
my_formula <- as.formula(paste(colnames(flall)[10], "~", paste(colnames(flall)[c(2, 3, 5, 7)], collapse = "*")))

anova_output2 <- list()
for(i in 10:156){
  my_col <- colnames(flall[i]) 
  my_formula <- as.formula(paste(my_col, "~", paste(colnames(flall)[c(2, 3, 5, 7)], collapse = "*")))
  sum_aov <- summary(aov(my_formula, flall))
  sum_tab <- sum_aov[[1]] %>% as.data.frame()
  sum_tab$ind <- my_col
  anova_output2[[i]] <- sum_tab
  print(my_formula)
}

(anova_output2)


all_output2 <- do.call(rbind, anova_output2)

all_output2

#write.csv(all_output2, "fluo_HL_LL_all_anovaoutput.csv")

#takes longer when add extra variables, from 2 up to 4
#halves the residual number by adding rep and within rep but also all sig
model=lm(flall$NPQ_L11 ~ flall$Accession * flall$Treatment * flall$Plate * flall$Within_Rep)
#model=lm(flall$NPQ_L11 ~ flall$Accession * flall$Treatment)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'flall$Accession', conf.level=0.95)
summary(TUKEY)
plot(TUKEY , las=1 , col="brown")

flx <- flall

##t tests and plots
#ALL SEEM TO OVERLAP HL LL NO AFFECT? 
par(mfrow = c(1, 2))
plot(NPQ_L1~Treatment, flx, notch = T)
plot(NPQ_L11~Treatment, flx, notch = T)
plot(NPQ_L3~Treatment, flx, notch = T)
plot(NPQ_L5~Treatment, flx, notch = T)
plot(Fv.Fm_L1~Treatment, flx, notch = T)
plot(Fv.Fm_L11~Treatment, flx, notch = T)
plot(Fv.Fm_L3~Treatment, flx, notch = T)
plot(Fv.Fm_L5~Treatment, flx, notch = T)
plot(fqfm_L1~Treatment, flx, notch = T)
plot(fqfm_L11~Treatment, flx, notch = T)
plot(fqfm_L3~Treatment, flx, notch = T)
plot(fqfm_L5~Treatment, flx, notch = T)
plot(QY_max~Treatment, flx, notch = T)

#group by light and see sig
names(hl_llmeans)
hl_llmeans$meanfqfm_L11_diff
hl_llmeans$meanNPQ_L11_diff
hl_llmeans$Accession
hl_llmeans$Light <- c("HL", "LL", "LL", "HL", "LL", "LL", "LL",
              "HL", "HL", "HL", "LL", "HL", "LL", "LL",
              "HL", "HL", "LL", "LL", "LL", "LL", "HL", "HL",
              "HL", "HL", "HL", "HL", "LL", "LL", "HL", "LL", "HL")
t.test(meanNPQ_L11_diff~Light, hl_llmeans)
t.test(meanQYmax_diff~Light, hl_llmeans)
boxplot(meanNPQ_L11_diff~Light,data=hl_llmeans)
boxplot(meanQYmax_diff~Light,data=hl_llmeans)

t.test(Fv.Fm_L1~Treatment, flx)
t.test(Fv.Fm_L11~Treatment, flx)

t.test(NPQ_L1~Treatment, flx)
t.test(NPQ_L11~Treatment, flx)
t.test(NPQ_L3~Treatment, flx)
t.test(NPQ_L5~Treatment, flx)

t.test(fqfm_L1~Treatment, flx)
t.test(fqfm_L11~Treatment,flx)
t.test(fqfm_L3~Treatment, flx)
t.test(fqfm_L5~Treatment,flx)

t.test(QY_max~Treatment, flx)

t.test(QY_max~Treatment, gr)

#ANOVA for Species and treatment
source("SBStat3.R")
mod1<-lm(Fv.Fm_L1~Treatment*Accession*Plate,data=flx)
summa(mod1$residuals)
#skewness low, octile skewness also low?
summaplot(mod1$residuals)
plot(mod1$fitted.values,mod1$residuals,pch=16)

anova(mod1) # treat and access sig, plate not sig, rep sig
#only use 3 var in model as crashes otherwise

mod2<-lm(Fv.Fm_L5~Treatment*Accession*Plate,data=flx)
summa(mod2$residuals)
#skewness low, octile skewness also low?
summaplot(mod2$residuals)
plot(mod2$fitted.values,mod2$residuals,pch=16)

anova(mod2) # accession and access treat interaction

mod3<-lm(NPQ_L1~Treatment*Accession*Plate,data=flx)
summa(mod3$residuals)
#skewness low, octile skewness also low?
summaplot(mod3$residuals)
plot(mod3$fitted.values,mod3$residuals,pch=16)

anova(mod3) # plate sig

mod4<-lm(NPQ_L11~Treatment*Accession*Plate,data=flx)
summa(mod4$residuals)
#skewness low, octile skewness also low?
summaplot(mod4$residuals)
plot(mod4$fitted.values,mod4$residuals,pch=16)

anova(mod4) # plate sig

mod5<-lm(fqfm_L1~Treatment*Accession*Plate,data=flx)
summa(mod5$residuals)
#skewness low, octile skewness also low?
summaplot(mod5$residuals)
plot(mod5$fitted.values,mod5$residuals,pch=16)

anova(mod5) # plate not sig but treatment interaction

mod6<-lm(fqfm_L11~Treatment*Accession*Plate,data=flx)
summa(mod6$residuals)
#skewness low, octile skewness also low?
summaplot(mod6$residuals)
plot(mod6$fitted.values,mod6$residuals,pch=16)

anova(mod6) # sig

mod7<-lm(QY_max~Treatment*Accession*Plate,data=flx)
summa(mod7$residuals)
#skewness low, octile skewness also low?
summaplot(mod7$residuals)
plot(mod7$fitted.values,mod7$residuals,pch=16)

anova(mod7) # treat and accession sig, plate not sig

mod8<-lm(fqfm_L5~Treatment*Accession*Plate,data=flx)
summa(mod8$residuals)
#skewness low, octile skewness also low?
summaplot(mod8$residuals)
plot(mod8$fitted.values,mod8$residuals,pch=16)

anova(mod8) #treatment accession sig, interaction

mod9<-lm(fqfm_L3~Treatment*Accession*Plate,data=flx)
summa(mod9$residuals)
#skewness low, octile skewness also low?
summaplot(mod9$residuals)
plot(mod9$fitted.values,mod9$residuals,pch=16)

anova(mod9) # more accession affect, small treat affect

mod10<-lm(NPQ_L3~Treatment*Accession*Plate,data=flx)
summa(mod10$residuals)
#skewness low, octile skewness also low?
summaplot(mod10$residuals)
plot(mod10$fitted.values,mod10$residuals,pch=16)

anova(mod10) # treat, access, interaction, small plate affect

mod11<-lm(NPQ_L5~Treatment*Accession*Plate,data=flx)
summa(mod11$residuals)
#skewness low, octile skewness also low?
summaplot(mod11$residuals)
plot(mod11$fitted.values,mod11$residuals,pch=16)

anova(mod11) # treat, access and interaction sig

#make QYmax by treatment plot for paper
gr$Treatment <- factor(gr$Treatment, levels=c("LL", "HL"))
#need to configure distances so less gaps between plots
#tiff('QYmax_boxplot_light_new.tiff', units="cm", width=10, height=15, res=300, compression = 'lzw')
pdf('QYmax_boxplot_light_new.pdf', width=10, height=15)
#par(mfrow = c(2, 4),  mar=c(2,4.5,2,2))
plot(QY_max ~ Treatment, gr,
     col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, ylab = "Fv/Fm",
        ylim=c(0.62,0.84), xlab = "",
        cex.lab=1.5, cex.axis=1.5
)
text(1, 0.84, "***", cex=3)
text(2.1, 0.84, "p = <0.0001", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
#title("C. QYmax", cex.main=1, adj=0)
stripchart(QY_max ~ Treatment, gr,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.3,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
dev.off()

#old version
plot(QY_max ~ Treatment, gr,
     col=c("#FF6347", "#999999"),
     las=1, whisklty=1, staplewex=0, ylab="QYmax",
     ylim=c(min(QY_max),max(QY_max)), xlab = "")
text(2, 0.835, "***", cex=2.5)

#boxplots with sig stars
#boxplot graphs per treatment
#sig diff per treatment, add stars
plot(QY_max ~ Treatment, flall,
     col=c("#FF6347", "#999999"),
     las=1, whisklty=1, staplewex=0,
     ylim=c(min(QY_max),max(QY_max) *1.1), xlab = "")
text(2, 0.90, "***", cex=2.5)

#subset and make line graph for NPQ and fqfm
#can do this per rep, per accession or as an average

str(flall)
names(flall)
#use summaries of variables
names(sum)
HL <- sum %>% filter(Treatment == "HL") %>% select(Accession, NPQ_L1, NPQ_L2, NPQ_L3, NPQ_L4, NPQ_L5, NPQ_L6,
                                                     NPQ_L7, NPQ_L8, NPQ_L9, NPQ_L10, NPQ_L11)
LL <- sum %>% filter(Treatment == "LL") %>% select(Accession, NPQ_L1, NPQ_L2, NPQ_L3, NPQ_L4, NPQ_L5, NPQ_L6,
                                                   NPQ_L7, NPQ_L8, NPQ_L9, NPQ_L10, NPQ_L11)
#HL <- sum %>% filter(Treatment == "HL") %>% select(Accession, fqfm_L1, fqfm_L2, fqfm_L3, fqfm_L4, fqfm_L5, fqfm_L6,
#                                                   fqfm_L7, fqfm_L8, fqfm_L9, fqfm_L10, fqfm_L11)
#LL <- sum %>% filter(Treatment == "LL") %>% select(Accession, fqfm_L1, fqfm_L2, fqfm_L3, fqfm_L4, fqfm_L5, fqfm_L6,
#                                                   fqfm_L7, fqfm_L8, fqfm_L9, fqfm_L10, fqfm_L11)
#ALL HL
#HL <- flall %>% filter(Treatment == "HL") %>% select(Accession, fqfm_L1, fqfm_L2, fqfm_L3, fqfm_L4, fqfm_L5, fqfm_L6,
#                                                     fqfm_L7, fqfm_L8, fqfm_L9, fqfm_L10, fqfm_L11)
#SPLIT PER REP
#HL <- flall %>% filter(Treatment == "HL", Rep == "HL1") %>% select(Accession, fqfm_L1, fqfm_L2, fqfm_L3, fqfm_L4, fqfm_L5, fqfm_L6,
#                                                                   fqfm_L7, fqfm_L8, fqfm_L9, fqfm_L10, fqfm_L11)
#LL <- flall %>% filter(Treatment == "LL") %>% select(Accession, fqfm_L1, fqfm_L2, fqfm_L3, fqfm_L4, fqfm_L5, fqfm_L6,
#                                                     fqfm_L7, fqfm_L8, fqfm_L9, fqfm_L10, fqfm_L11)

#sub <- (flall[,145:155]) #just fqfm cols
#to see hl and ll seperately need to split by treatment 
#sa_comb <- sub
sa_comb <- HL #all HL and gives accession column back
#sa_comb <- LL
names(sa_comb)

sa_comb$Accession <- NULL

#row.names(sa_comb) <- sa_comb$Accession
sa_comb$Accession <- NULL #removes accession column again
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
row.names(sadf) <- str_replace_all(row.names(sadf), "NPQ_L", "")
#row.names(sadf) <- str_replace_all(row.names(sadf), "fqfm_L", "")
#looks for pattern and give replacement
#check its worked
row.names(sadf)

time_vec <- row.names(sadf)
time_vec <- as.numeric(time_vec) #coerce to numeric

class(sadf) 

accessions <- (HL[,1]) #define Accessions variable
#accessions <- (LL[,1]) #define Accessions variable
accessions <- t(accessions) #make transverse so correct no of cols
colnames(sadf) <- (accessions) #put names back on matrix as col names

sadf$time <- time_vec #make time column

str(sadf)

class(sadf)

library(tidyr)

sadf$APP2
sadf$time

#par(mfrow = c(1,1)) #no of rows and cols in plot display
#plot(KS02 ~ time, sadf, type = "n", xlab = "Time (days)", ylab = "RGR col gain T14-T21/7 (gain per day)")
#lines(KS03 ~ time, sadf) # one at a time
#lines(SEL1 ~ time, sadf) 
#lines(KS04 ~ time, sadf) # one at a time
#lines(LY03 ~ time, sadf) 

#do something to every i
par(mfrow = c(1, 1))
plot(APP2 ~ time, sadf, type = "n", 
     main = "Photosynthetic efficiency in increasing light (HL treatment)",
     xlab = "Steps", 
     ylab = expression(paste("FqFm")), 
     axes = F,
     xlim = c(0, 12),
     ylim = c(0, 1))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:31){
  lines(sadf$time, sadf[, i], col = line_colours[i], lty = sort(line_type)[i], lwd = 2)
  
}
axis(1, las = 1) #hash out so doesnt add x labels as defined
#axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
legend("bottomright", pt.cex = 150, cex = 0.3, col = line_colours, lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()

line_colours <- rep(1:31, length.out = ncol(sadf)-1) #8 colors
line_colours
line_type <- rep(1:3, length.out =  ncol(sadf)-1)
line_type

#NPQ
#do something to every i
par(mfrow = c(1, 1))
plot(APP2 ~ time, sadf, type = "n", 
     main = "NPQ in increasing light (LL treatment)",
     xlab = "Steps", 
     ylab = expression(paste("NPQ")), 
     axes = F,
     xlim = c(0, 12),
     ylim = c(0, 5))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:31){
  lines(sadf$time, sadf[, i], col = line_colours[i], lty = sort(line_type)[i], lwd = 2)
  
}
axis(1, las = 1) #hash out so doesnt add x labels as defined
#axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
legend("bottomright", pt.cex = 150, cex = 0.3, col = line_colours, lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()

line_colours <- rep(1:31, length.out = ncol(sadf)-1) #8 colors
line_colours
line_type <- rep(1:3, length.out =  ncol(sadf)-1)
line_type

