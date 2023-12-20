#script to explore hl ll dw growth rate data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Xanthophyll_data_working.csv")
gr <- read.csv("Xanthophyll_data_working_addspNEW+envlight.csv")

tail(gr)
head(gr)
names(gr)
#names
#[1] "Accession"      "Treatment"      "chl.a.b.ratio"  "car..chl.ratio"
#[5] "DEPs"           "XC.pool"        "X.XC.pool.Car"  "X.N.Car"       
#[9] "X.V.Car"        "X.A.Car"        "X.L.Car"        "X.Z.Car"       
#[13] "X...car.car"    "N.L"            "XC.N"           "XC.L"          
#[17] "N.Xan"          "XC.xan"         "L.xan"          "V.N"           
#[21] "Z.N" 

#check classes of all cols
sapply(gr, class)

library(dplyr)
library(stringr)
library(ggplot2)

gr %>% select(Accession, Treatment) #just displays them
unique(gr$Accession, gr$Treatment) 
unique(gr$Accession) #32 #now 28 without genera
length(unique(gr$Accession, gr$Treatment)) #22
length(unique(gr$Accession,order = ascending))
#20

as.factor(gr$Accession)
as.factor(gr$Treatment)

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
ggplot(gr, aes(x = Accession, y = chl.a.b.ratio)) + geom_bar(stat = "identity")
ggplot(gr, aes(y=chl.a.b.ratio, x=Accession)) +
  geom_col(position='dodge', color='black')

names(gr)
#affects of other variables
aov <- aov(chl.a.b.ratio~Treatment,data=gr) #sig
aov <- aov(car..chl.ratio~Treatment,data=gr) #not sig
aov <- aov(DEPs~Treatment,data=gr) #sig
aov <- aov(XC.pool~Treatment,data=gr) #not sig
aov <- aov(X.XC.pool.Car~Treatment,data=gr) #sig
summary(aov)

#affects of species and envlight
#affects of species and envlight
aov <- aov(chl.a.b.ratio~Treatment*Species,data=gr) #treat *
aov <- aov(car..chl.ratio~Treatment*Species,data=gr) #ns
aov <- aov(DEPs~Treatment*Species,data=gr) #treat ***
aov <- aov(XC.pool~Treatment*Species,data=gr) #**species, int treat species
aov <- aov(X.XC.pool.Car~Treatment*Species,data=gr) #treat ***
summary(aov)
tuk_out <- TukeyHSD(aov, conf.level=.95)
tuk_out

aov <- aov(chl.a.b.ratio~Species,data=gr) #ns
aov <- aov(car..chl.ratio~Species,data=gr) #ns
aov <- aov(DEPs~Species,data=gr) #ns
aov <- aov(XC.pool~Species,data=gr) #* species
aov <- aov(X.XC.pool.Car~Species,data=gr) #ns
summary(aov)

aov <- aov(chl.a.b.ratio~EnvLight,data=gr) #ns
aov <- aov(car..chl.ratio~EnvLight,data=gr) #ns
aov <- aov(DEPs~EnvLight,data=gr) #ns
aov <- aov(XC.pool~EnvLight,data=gr) #ns
aov <- aov(X.XC.pool.Car~EnvLight,data=gr) #ns
summary(aov)

aov <- aov(chl.a.b.ratio~Species*EnvLight,data=gr) #ns
aov <- aov(car..chl.ratio~Species*EnvLight,data=gr) #ns
aov <- aov(DEPs~Species*EnvLight,data=gr) #ns
aov <- aov(XC.pool~Species*EnvLight,data=gr) #species *
aov <- aov(X.XC.pool.Car~Species*EnvLight,data=gr) #ns
summary(aov)

t.test(chl.a.b.ratio~Treatment,data=gr) #sig
t.test(car..chl.ratio~Treatment,data=gr) #not sig
t.test(DEPs~Treatment,data=gr) #sig
t.test(XC.pool~Treatment,data=gr) #not sig
t.test(X.XC.pool.Car~Treatment,data=gr) #sig

barplot(gr$chl.a.b.ratio)
barplot(gr$car..chl.ratio)
barplot(gr$DEPs)
barplot(gr$XC.pool)
barplot(gr$X.XC.pool.Car)

#with accession
aov <- aov(chl.a.b.ratio~Treatment+Accession,data=gr) #treat sig, access not
aov <- aov(car..chl.ratio~Treatment+Accession,data=gr) #ns
aov <- aov(DEPs~Treatment+Accession,data=gr) #*** treat, ns access
aov <- aov(XC.pool~Treatment+Accession,data=gr) #ns
aov <- aov(X.XC.pool.Car~Treatment+Accession,data=gr) #*** treat, ns access
summary(aov)

aov <- aov(X.N.Car~Treatment,data=gr) #not sig
aov <- aov(X.V.Car~Treatment,data=gr) #sig*
aov <- aov(X.A.Car~Treatment,data=gr) #sig
aov <- aov(X.L.Car~Treatment,data=gr) #sig
aov <- aov(X.Z.Car~Treatment,data=gr) #sig
summary(aov)

t.test(X.N.Car~Treatment,data=gr) #not sig
t.test(X.V.Car~Treatment,data=gr) #sig higher ll
t.test(X.A.Car~Treatment,data=gr) #sig higher hl
t.test(X.L.Car~Treatment,data=gr) #sig higher ll
t.test(X.Z.Car~Treatment,data=gr) #sig higher hl
t.test(X...car.car~Treatment,data=gr) #ns

#with accession
aov <- aov(X.N.Car~Treatment+Accession,data=gr) #not sig
aov <- aov(X.V.Car~Treatment+Accession,data=gr) #treat *
aov <- aov(X.A.Car~Treatment+Accession,data=gr) #treat *
aov <- aov(X.L.Car~Treatment+Accession,data=gr) #** treatment
aov <- aov(X.Z.Car~Treatment+Accession,data=gr) #treat ***
summary(aov)

aov <- aov(X...car.car~Treatment,data=gr) #not sig
aov <- aov(N.L~Treatment,data=gr) #not sig
aov <- aov(XC.N~Treatment,data=gr) #sig
aov <- aov(XC.L~Treatment,data=gr) #sig
aov <- aov(N.Xan~Treatment,data=gr) #not sig
aov <- aov(XC.xan~Treatment,data=gr) #sig
aov <- aov(L.xan~Treatment,data=gr) #sig
aov <- aov(V.N~Treatment,data=gr) #not sig
aov <- aov(Z.N~Treatment,data=gr) #sig
summary(aov)

tukey <- TukeyHSD(aov, "Treatment", conf.level=.95)
str(tukey)
tukey

boxplot(chl.a.b.ratio~Accession*Treatment,data=gr)

par(mfrow = c(2,5)) 
boxplot(chl.a.b.ratio~Treatment,data=gr)
boxplot(car..chl.ratio~Treatment,data=gr)
boxplot(XC.pool~Treatment,data=gr)
boxplot(X.XC.pool.Car~Treatment,data=gr)
boxplot(DEPs~Treatment,data=gr)
boxplot(X.N.Car~Treatment,data=gr) 
boxplot(X.L.Car~Treatment,data=gr) #sig
boxplot(X.V.Car~Treatment,data=gr) 
boxplot(X.A.Car~Treatment,data=gr) 
boxplot(X.Z.Car~Treatment,data=gr) #sig

boxplot(XC.pool~Treatment*Species,data=gr)

t.test(X...car.car~Treatment,data=gr) #not sig
t.test(N.L~Treatment,data=gr) #not sig
t.test(XC.N~Treatment,data=gr) #sig
t.test(XC.L~Treatment,data=gr) #sig
t.test(N.Xan~Treatment,data=gr) #not sig
t.test(XC.xan~Treatment,data=gr) #sig
t.test(L.xan~Treatment,data=gr) #sig
t.test(V.N~Treatment,data=gr) #not sig
t.test(Z.N~Treatment,data=gr) #sig

boxplot(X...car.car~Treatment,data=gr) #not sig
boxplot(N.L~Treatment,data=gr) #not sig
boxplot(XC.N~Treatment,data=gr) #sig
boxplot(XC.L~Treatment,data=gr) #sig
boxplot(N.Xan~Treatment,data=gr) #not sig
boxplot(XC.xan~Treatment,data=gr) #sig
boxplot(L.xan~Treatment,data=gr) #sig
boxplot(V.N~Treatment,data=gr) #not sig
boxplot(Z.N~Treatment,data=gr) #sig

#plot boxplots in style for paper
par(mfrow = c(3, 5),  mar=c(2,4.5,2,2))
boxplot(chl.a.b.ratio ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 3.2), xlab = "",
        ylab = "Chl a:b",
        cex.lab=1.5, cex.axis=1.0
)
text(1, 3.1, "**", cex=1.75)
boxplot(car..chl.ratio ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 3), xlab = "",
        ylab = "Car:chl",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 3, substitute(paste(italic("n.s"))), cex=1.75)
boxplot(XC.pool ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 100), xlab = "",
        ylab = "XC.pool",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 100, substitute(paste(italic("n.s"))), cex=1.75)
boxplot(X.XC.pool.Car ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 50), xlab = "",
        ylab = "XC.pool:Car",
        cex.lab=1.5, cex.axis=1.0
)
text(1, 50, "***", cex=1.75)
boxplot(DEPs ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 100), xlab = "",
        ylab = "DEPs",
        cex.lab=1.5, cex.axis=1.0
)
text(1, 100, "***", cex=1.75)
boxplot(X.N.Car ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 60), xlab = "",
        ylab = "% Neo:Car",
        cex.lab=1.5, cex.axis=1.0
)
text(1, 60, substitute(paste(italic("n.s"))), cex=1.75)
boxplot(X.L.Car ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 60), xlab = "",
        ylab = "% Lut:Car",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 60, "***", cex=1.75)
boxplot(X.V.Car ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 60), xlab = "",
        ylab = "% Vio:Car",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 60, "*", cex=1.75)
boxplot(X.A.Car ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 60), xlab = "",
        ylab = "% Ant:Car",
        cex.lab=1.5, cex.axis=1.0
)
text(1, 60, "*", cex=1.75)
boxplot(X.Z.Car ~ Treatment, gr,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 60), xlab = "",
        ylab = "% Zea:Car",
        cex.lab=1.5, cex.axis=1.0
)
text(1, 60, "***", cex=1.75)

#boxplots with accessions
#change to parameter you want
chl.a.b.ratio_boxplot <- ggplot(gr,aes(x = reorder(Accession,chl.a.b.ratio, FUN = median),y=chl.a.b.ratio)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 2) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) chl.a.b.ratio"),subtitle = ("Sort by Accession"),x=expression(),
       y=expression(chl.a.b.ratio))+
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
chl.a.b.ratio_boxplot

ggsave("HL_LL_chl.a.b.ratio_boxplot.tiff", dpi = 600, width = 20, height = 15, units = "cm")


car..chl.ratio_boxplot <- ggplot(gr,aes(x = reorder(Accession,car..chl.ratio, FUN = median),y=car..chl.ratio)) +
  
  #chane light to treatment
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 2) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) Car:Chl ratio"),subtitle = ("Arranged by Accession"),x=expression(),
       y=expression(car..chl.ratio))+
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
car..chl.ratio_boxplot

ggsave("HL_LL_car..chl.ratio_boxplot.tiff", dpi = 600, width = 20, height = 15, units = "cm")

#DEPs looks neatest greenness parameter
#change to parameter you want
DEPs_boxplot <- ggplot(gr,aes(x = reorder(Accession,DEPs, FUN = median),y=DEPs)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #           linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 2) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) DEPs "),subtitle = ("Arranged by Accession"),x=expression(),
       y=expression(DEPs))+
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
DEPs_boxplot

ggsave("HL_LL_DEPs_boxplot.tiff", dpi = 600, width = 20, height = 15, units = "cm")


#change to parameter you want
XC.pool_boxplot <- ggplot(gr,aes(x = reorder(Accession,XC.pool, FUN = median),y=XC.pool)) +
  
  #geom_hline(aes(yintercept=MaxNPQ_mean),                                                                          #WHAT IS THIS DASHED LINE? 
  #           linetype="dashed", color = "red", size = 1)+ 
  geom_boxplot() +
  facet_wrap(. ~ Treatment, scales = "free_x", ncol = 2) + 
  scale_fill_brewer(palette = "Paired")+			 
  theme_classic() +
  theme(axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  theme(panel.spacing=unit(0,"pt"),panel.border=element_rect(colour="black", fill=NA),
        panel.background = element_rect(size = 1, linetype = "solid"))+
  labs(title =expression("(a) XC.pool "),subtitle = ("Arranged by Accession"),x=expression(),
       y=expression(XC.pool))+
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
XC.pool_boxplot

ggsave("HL_LL_XC.pool_boxplot.tiff", dpi = 600, width = 20, height = 15, units = "cm")

#get as 'averages' so one per treat/access

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


#Tukey works when not means data summary
tuk_out <- TukeyHSD(ChlA_aov, "Rep", conf.level=.95)
tuk_out <- TukeyHSD(ChlA_aov, "Treatment", conf.level=.95)
tuk_out <- TukeyHSD(ChlA_aov, "Accession", conf.level=.95)
tuk_out <- TukeyHSD(ChlA_aov, "Species", conf.level=.95)
str(tuk_out)
tuk_out
plot(tuk_out , las=1 , col="brown") #not useful to visualise

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
sum$chl.a.b.ratio
sum$car..chl.ratio
sum$DEPs
sum$XC.pool
sum$X.XC.pool.Car
sum$X.N.Car

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
  summarise(chl.a.b.ratio = mean(chl.a.b.ratio), car..chl.ratio = mean(car..chl.ratio),
            DEPs = mean(DEPs), XC.pool = mean(XC.pool),
            X.XC.pool.Car = mean(X.XC.pool.Car), X.N.Car = mean(X.N.Car),
            X.V.Car = mean(X.V.Car), X.A.Car = mean(X.A.Car),
            X.L.Car = mean(X.L.Car), X.Z.Car = mean(X.Z.Car),
            X...car.car = mean(X...car.car), N.L = mean(N.L),
            XC.N = mean(XC.N), XC.L = mean(XC.L),
            N.Xan = mean(N.Xan), XC.xan = mean(XC.xan),
            L.xan = mean(L.xan), V.N = mean(V.N),
            Z.N = mean(Z.N))
Summary


#for inc Sp manually sorted into HL LL and created summary
#to add lag to density need to summarise per treatment and set as cols

#order by treatment
try1 <- split(Summary,cumsum(1:nrow(Summary)%in%18)) #made 2 lists 0 and 1
str(try1)
#call the two lists
try1$`0`
try1$`1`
backtogether <- cbind(try1$`0`, try1$`1`)
write.csv(backtogether, "Xanthophyll_Summary.csv") #given colnames manually

hl_llmeans <- read.csv("Xanthophyll_Summary_completepairs_x.csv")

#add envlight
hl_llmeans$EnvLight <- c("dLL", "dLL", "dHL", "dLL", "dLL", "dHL",
                         "dHL", "dLL", "dLL", "dLL", "dLL", "dHL", "dHL",
                         "dLL")

#ks09 ll
hl_llmeans$EnvLight <- c("dLL", "dLL", "dHL", "dLL", "dLL", "dLL",
                         "dHL", "dLL", "dLL", "dLL", "dLL", "dHL", "dHL",
                         "dLL")
gr <- hl_llmeans

t.test(chl.a.b.ratio_HL~EnvLight,data=gr) #not sig
t.test(car..chl.ratio_HL~EnvLight,data=gr) #not sig
t.test(DEPs_HL~EnvLight,data=gr) #not sig
t.test(XC.pool_HL~EnvLight,data=gr) #not sig
t.test(X.XC.pool.Car_HL~EnvLight,data=gr) #not sig
t.test(X.N.Car_HL~EnvLight,data=gr) #not sig
t.test(X.V.Car_HL~EnvLight,data=gr) #not sig
t.test(X.A.Car_HL~EnvLight,data=gr) #not sig
t.test(X.L.Car_HL~EnvLight,data=gr) #not sig
t.test(X.Z.Car_HL~EnvLight,data=gr) #not sig

t.test(X...car.car_HL~EnvLight,data=gr) #not sig
t.test(N.L_HL~EnvLight,data=gr) #not sig
t.test(XC.N_HL~EnvLight,data=gr) #not sig
t.test(XC.L_HL~EnvLight,data=gr) #not sig
t.test(N.Xan_HL~EnvLight,data=gr) #not sig
t.test(XC.xan_HL~EnvLight,data=gr) #not sig
t.test(L.xan_HL~EnvLight,data=gr) #not sig
t.test(V.N_HL~EnvLight,data=gr) #not sig
t.test(Z.N_HL~EnvLight,data=gr) #not sig

boxplot(chl.a.b.ratio_HL~EnvLight,data=gr) #not sig
boxplot(car..chl.ratio_HL~EnvLight,data=gr) #not sig
boxplot(DEPs_HL~EnvLight,data=gr) #not sig
boxplot(XC.pool_HL~EnvLight,data=gr) #not sig
boxplot(X.XC.pool.Car_HL~EnvLight,data=gr) #not sig
boxplot(X.N.Car_HL~EnvLight,data=gr) #not sig
boxplot(X.V.Car_HL~EnvLight,data=gr) #not sig
boxplot(X.A.Car_HL~EnvLight,data=gr) #not sig
boxplot(X.L.Car_HL~EnvLight,data=gr) #not sig
boxplot(X.Z.Car_HL~EnvLight,data=gr) #not sig

boxplot(X...car.car_HL~EnvLight,data=gr) #not sig
boxplot(N.L_HL~EnvLight,data=gr) #not sig
boxplot(XC.N_HL~EnvLight,data=gr) #not sig
boxplot(XC.L_HL~EnvLight,data=gr) #not sig
boxplot(N.Xan_HL~EnvLight,data=gr) #not sig
boxplot(XC.xan_HL~EnvLight,data=gr) #not sig
boxplot(L.xan_HL~EnvLight,data=gr) #not sig
boxplot(V.N_HL~EnvLight,data=gr) #not sig
boxplot(Z.N_HL~EnvLight,data=gr) #not sig

#LL
t.test(chl.a.b.ratio_LL~EnvLight,data=gr) #not sig
t.test(car..chl.ratio_LL~EnvLight,data=gr) #not sig
t.test(DEPs_LL~EnvLight,data=gr) #not sig
t.test(XC.pool_LL~EnvLight,data=gr) #not sig
t.test(X.XC.pool.Car_LL~EnvLight,data=gr) #not sig
t.test(X.N.Car_LL~EnvLight,data=gr) #not sig
t.test(X.V.Car_LL~EnvLight,data=gr) #not sig
t.test(X.A.Car_LL~EnvLight,data=gr) #not sig
t.test(X.L.Car_LL~EnvLight,data=gr) #not sig
t.test(X.Z.Car_LL~EnvLight,data=gr) #not sig

boxplot(chl.a.b.ratio_LL~EnvLight,data=gr) #not sig
boxplot(car..chl.ratio_LL~EnvLight,data=gr) #not sig
boxplot(DEPs_LL~EnvLight,data=gr) #not sig
boxplot(XC.pool_LL~EnvLight,data=gr) #not sig
boxplot(X.XC.pool.Car_LL~EnvLight,data=gr) #not sig
boxplot(X.N.Car_LL~EnvLight,data=gr) #not sig
boxplot(X.V.Car_LL~EnvLight,data=gr) #not sig
boxplot(X.A.Car_LL~EnvLight,data=gr) #not sig
boxplot(X.L.Car_LL~EnvLight,data=gr) #not sig
boxplot(X.Z.Car_LL~EnvLight,data=gr) #not sig

t.test(X...car.car_LL~EnvLight,data=gr) #not sig
t.test(N.L_LL~EnvLight,data=gr) #not sig
t.test(XC.N_LL~EnvLight,data=gr) #not sig
t.test(XC.L_LL~EnvLight,data=gr) #not sig
t.test(N.Xan_LL~EnvLight,data=gr) #not sig
t.test(XC.xan_LL~EnvLight,data=gr) #not sig
t.test(L.xan_LL~EnvLight,data=gr) #not sig
t.test(V.N_LL~EnvLight,data=gr) #not sig
t.test(Z.N_LL~EnvLight,data=gr) #not sig

boxplot(X...car.car_LL~EnvLight,data=gr) #not sig
boxplot(N.L_LL~EnvLight,data=gr) #not sig
boxplot(XC.N_LL~EnvLight,data=gr) #not sig
boxplot(XC.L_LL~EnvLight,data=gr) #not sig
boxplot(N.Xan_LL~EnvLight,data=gr) #not sig
boxplot(XC.xan_LL~EnvLight,data=gr) #not sig
boxplot(L.xan_LL~EnvLight,data=gr) #not sig
boxplot(V.N_LL~EnvLight,data=gr) #not sig
boxplot(Z.N_LL~EnvLight,data=gr) #not sig

#add new cols for differences
library(dplyr)
# hl - ll to get differences
hl_llmeans <- hl_llmeans %>% mutate(chl.a.b.ratio_diff = (chl.a.b.ratio_HL - chl.a.b.ratio_LL))
hl_llmeans <- hl_llmeans %>% mutate(car..chl.ratio_diff = (car..chl.ratio_HL - car..chl.ratio_LL))
hl_llmeans <- hl_llmeans %>% mutate(DEPs_diff = (DEPs_HL - DEPs_LL))
hl_llmeans <- hl_llmeans %>% mutate(XC.pool_diff = (XC.pool_HL - XC.pool_LL))
hl_llmeans <- hl_llmeans %>% mutate(X.XC.pool.Car_diff = (X.XC.pool.Car_HL - X.XC.pool.Car_LL))
hl_llmeans <- hl_llmeans %>% mutate(X.N.Car_diff = (X.N.Car_HL - X.N.Car_LL))
hl_llmeans <- hl_llmeans %>% mutate(X.V.Car_diff = (X.V.Car_HL - X.V.Car_LL))
hl_llmeans <- hl_llmeans %>% mutate(X.A.Car_diff = (X.A.Car_HL - X.A.Car_LL))
hl_llmeans <- hl_llmeans %>% mutate(X.L.Car_diff = (X.L.Car_HL - X.L.Car_LL))
hl_llmeans <- hl_llmeans %>% mutate(X.Z.Car_diff = (X.Z.Car_HL - X.Z.Car_LL))
hl_llmeans <- hl_llmeans %>% mutate(X...car.car_diff = (X...car.car_HL - X...car.car_LL))

#get proportional differences
hl_llmeans <- hl_llmeans %>% mutate(chl.a.b.ratio_propdiff = (chl.a.b.ratio_diff/chl.a.b.ratio_LL)*100)
hl_llmeans %>% arrange(desc(chl.a.b.ratio_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(car..chl.ratio_propdiff = (car..chl.ratio_diff/car..chl.ratio_LL)*100)
hl_llmeans %>% arrange(desc(car..chl.ratio_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(DEPs_propdiff = (DEPs_diff/DEPs_LL)*100)
hl_llmeans %>% arrange(desc(DEPs_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(XC.pool_propdiff = (XC.pool_diff/XC.pool_LL)*100)
hl_llmeans %>% arrange(desc(XC.pool_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X.XC.pool.Car_propdiff = (X.XC.pool.Car_diff/X.XC.pool.Car_LL)*100)
hl_llmeans %>% arrange(desc(X.XC.pool.Car_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X.N.Car_propdiff = (X.N.Car_diff/X.N.Car_LL)*100)
hl_llmeans %>% arrange(desc(X.N.Car_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X.V.Car_propdiff = (X.V.Car_diff/X.V.Car_LL)*100)
hl_llmeans %>% arrange(desc(X.V.Car_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X.A.Car_propdiff = (X.A.Car_diff/X.A.Car_LL)*100)
hl_llmeans %>% arrange(desc(X.A.Car_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X.L.Car_propdiff = (X.L.Car_diff/X.L.Car_LL)*100)
hl_llmeans %>% arrange(desc(X.L.Car_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X.Z.Car_propdiff = (X.Z.Car_diff/X.Z.Car_LL)*100)
hl_llmeans %>% arrange(desc(X.Z.Car_propdiff)) %>% select(Accession) %>% top_n(14)
hl_llmeans <- hl_llmeans %>% mutate(X...car.car_propdiff = (X...car.car_diff/X...car.car_LL)*100)
hl_llmeans %>% arrange(desc(X...car.car_propdiff)) %>% select(Accession) %>% top_n(14)

#plot differences as bar plots
library(ggplot2)

par(mfrow = c(1,2)) 
ggplot(hl_llmeans, aes(y=chl.a.b.ratio_propdiff, x=reorder(Accession, chl.a.b.ratio_propdiff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in Chl a:b") +
  xlab("Ecotype")

ggplot(hl_llmeans, aes(y=car..chl.ratio_propdiff, x=reorder(Accession, car..chl.ratio_propdiff), fill=EnvLight)) +
  geom_col(position='dodge') +
  ylab("Difference in Car:chl") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=DEPs_propdiff, x=reorder(Accession, DEPs_propdiff), fill=EnvLight)) +
  geom_col(position='dodge') +
  ylab("Difference in DEPs") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=XC.pool_propdiff, x=reorder(Accession, XC.pool_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in XC pool") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X.XC.pool.Car_propdiff, x=reorder(Accession, X.XC.pool.Car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in XC pool/Car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X.N.Car_propdiff, x=reorder(Accession, X.N.Car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in N/Car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X.V.Car_propdiff, x=reorder(Accession, X.V.Car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in V/car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X.A.Car_propdiff, x=reorder(Accession, X.A.Car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in A/car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X.L.Car_propdiff, x=reorder(Accession, X.L.Car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in L/Car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X.Z.Car_propdiff, x=reorder(Accession, X.Z.Car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in Z/car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

ggplot(hl_llmeans, aes(y=X...car.car_propdiff, x=reorder(Accession, X...car.car_propdiff), fill=EnvLight)) +
  geom_col(position='dodge', color='black') +
  ylab("Difference in B-carotene/Car") +
  xlab("Ecotypes") +
  theme_bw() +
  theme_classic()

write.csv(hl_llmeans, "Xanthophyll_Summary_completepairs_x+diffs.csv")
