setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
#Lightpercent_dhl <- read.csv("Light_dhl_%_for_bargraphs.csv", header = TRUE)
#names(Lightpercent_dhl)[c(1,2,3)]<- c("Light", "Recording", "Time")
#Lightpercent_dll <- read.csv("Light_dll_%_for_bargraphs.csv", header = TRUE)
#names(Lightpercent_dll)[c(1,2,3)]<- c("Light", "Recording", "Time")

#old %
#Lightpercent_dhl <- read.csv("Light_dhl_%_for_bargraphs_ADDPFD.csv", header = TRUE)
#names(Lightpercent_dhl)[c(1,2,3)]<- c("Light", "Recording", "Time")
#Lightpercent_dll <- read.csv("Light_dll_%_for_bargraphs_ADDPFD.csv", header = TRUE)
#names(Lightpercent_dll)[c(1,2,3)]<- c("Light", "Recording", "Time")

#RAW NEW WITH PFD
Light_dhl <- read.csv("Light_dHL_for_bargraphs_PFDADD.csv", header = TRUE)
names(Light_dhl)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")

Light_dll <- read.csv("Light_dLL_for_bargraphs_PFDADD.csv", header = TRUE)
names(Light_dhl)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")

#new % with sds
Lightpercent_dhl <- read.csv("Light_dhl_%_for_bargraphs_ADDPFD_raw.csv", header = TRUE)
names(Light_dhl)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")
Lightpercent_dll <- read.csv("Light_dll_%_for_bargraphs_ADDPFD_raw.csv", header = TRUE)
names(Light_dll)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")

#RAW
#Light_dhl <- read.csv("Light_dHL_for_bargraphs.csv", header = TRUE)
#names(Light_dhl)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")
# Graph
#Light_dll <- read.csv("Light_dLL_for_bargraphs.csv", header = TRUE)
#names(Light_dll)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")
# Graph

#just frf
Light_rfr_dhl <- read.csv("Light_RFR_dHL_for_bargraphs.csv", header = TRUE)
names(Light_rfr_dhl)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")
Light_rfr_dll <- read.csv("Light_RFR_dll_for_bargraphs.csv", header = TRUE)
names(Light_rfr_dll)[c(1,2,3,4)]<- c("Light", "Recording", "Time", "SD")

class(Lightpercent_dhl$Light)
class(Lightpercent_dhl$Recording)
class(Lightpercent_dhl$Time)
class(Lightpercent_dll$Light)
class(Lightpercent_dll$Recording)
class(Lightpercent_dll$Time)

class(Light_dhl$Light)
class(Light_dhl$Recording)
class(Light_dhl$Time)

library(ggplot2)
library(plyr)
Lightpercent_dhlSumm <- ddply(Lightpercent_dhl, c("Light", "Time"),
                             summarise,
                             N = length(Recording),
                             mean = mean(Recording),
                             sd = sd(Recording),
                             se = sd/(sqrt(N)))

class(Lightpercent_dhl$mean)

#not done per site so dont have sd se N
#not using sum just using raw as no err bars
ggplot(data = Lightpercent_dhlSumm, 
       aes(x = as.factor(Light), y = mean, fill = as.factor(Time)))

#re arrange order of time
Lightpercent_dhl$Time <- factor(Lightpercent_dhl$Time, levels = c("Spring 2021", "Summer 2021",
                                               "Autumn 2021", "Spring 2022",
                                               "Summer 2022"))
par(mfrow = c(1, 2))
p6 <- ggplot(data = Lightpercent_dhl, 
       aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin = Recording - SD, ymax = Recording+SD), width = 0.2,
  position = position_dodge(0.7))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(0, 100))+
  theme(legend.position = "none")+
  ggtitle("dHL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("%PPFD", "%FR", "%R", "%G", "%B", "%UV"),
                   labels = c("PPFD%", "FR%", 
                              "R%", "G%", "B%", "UV%"))+
  labs(x="Light",
       y="Proportion of light (%)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p6

#dll graph
Lightpercent_dll$Time <- factor(Lightpercent_dll$Time, levels = c("Spring 2021", "Summer 2021",
                                                                  "Autumn 2021", "Spring 2022",
                                                                  "Summer 2022"))
p5 <- ggplot(data = Lightpercent_dll, 
       aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin = Recording - SD, ymax = Recording+SD), width = 0.2,
  position = position_dodge(0.7))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(0, 100))+
  theme(legend.position = "none")+
  ggtitle("dLL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("%PPFD", "%FR", "%R", "%G", "%B", "%UV"),
                   labels = c("PPFD%", "FR%", 
                              "R%", "G%", "B%", "UV%"))+
  labs(x="Light",
       y="Proportion of light (%)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
  
  p5

#test grid function
#install.packages("gridExtra")
library(gridExtra)

grid.arrange(p5, p6, nrow = 1)

#DHL RAW
Light_dhl$Time <- factor(Light_dhl$Time, levels = c("Spring 2021", "Summer 2021",
                                                                  "Autumn 2021", "Spring 2022",
                                                                  "Summer 2022"))
#par(mfrow = c(1, 2)) only works with base r
p2 <- ggplot(data = Light_dhl, 
       aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
  position = position_dodge(0.7))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1500))+
  theme(legend.position = "none")+
  ggtitle("dHL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("PFD", "PPFD", "FR", "R", "G", "B", "UV"),
                   labels = c("PFD", "PPFD", "FR", "R", "G", "B", "UV"))+
  labs(x="Light",
       y="Light intensity (??mol m^-2~s^-1)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p2

#LAMBDA DATA ALONE
Light_dhl$Time <- factor(Light_dhl$Time, levels = c("Spring 2021", "Summer 2021",
                                                    "Autumn 2021", "Spring 2022",
                                                    "Summer 2022"))
par(mfrow = c(1, 2))
p4 <- ggplot(data = Light_dhl, 
       aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.4,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
                position = position_dodge(0.4))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0), limits = c(0, 900))+
  theme(legend.position = "none")+
  ggtitle("dHL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("LambdaP", "LambdaD"),
                   labels = c("??p", "??d"))+
  labs(x="Light",
       y="Wavelength (nm)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p4

#raw LL
Light_dll$Time <- factor(Light_dll$Time, levels = c("Spring 2021", "Summer 2021",
                                                    "Autumn 2021", "Spring 2022",
                                                    "Summer 2022"))
par(mfrow = c(1, 2))
p1 <- ggplot(data = Light_dll, 
       aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
                position = position_dodge(0.7))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1500))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  ggtitle("dLL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("PFD", "PPFD", "FR", "R", "G", "B", "UV"),
                   labels = c("PFD", "PPFD", "FR", "R", "G", "B", "UV"))+
  labs(x="Light",
       y=bquote("Light intensity"~(mu~mol~ m^-2~s^-1)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p1

#LAMBDA DATA ALONE
Light_dll$Time <- factor(Light_dll$Time, levels = c("Spring 2021", "Summer 2021",
                                                    "Autumn 2021", "Spring 2022",
                                                    "Summer 2022"))

p3 <- ggplot(data = Light_dll, 
       aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.4,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
                position = position_dodge(0.4))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 900))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  ggtitle("dLL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("LambdaP", "LambdaD"),
                   labels = c("??p", "??d"))+
  labs(x="Light",
       y="Wavelength (nm)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

p3

#rfr playing
Light_rfr_dhl$Time <- factor(Light_rfr_dhl$Time, levels = c("Spring 2021", "Summer 2021",
                                                                        "Autumn 2021", "Spring 2022",
                                                                        "Summer 2022"))

p8 <- ggplot(data = Light_rfr_dhl, 
       aes(x = Light, y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.3,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
                position = position_dodge(0.3))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 2))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")+
  ggtitle("dHL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("R:FR"),
                   labels = c("R:FR"))+
  labs(x="Light",
       y="R:FR")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#rfr playing
Light_rfr_dll$Time <- factor(Light_rfr_dll$Time, levels = c("Spring 2021", "Summer 2021",
                                                            "Autumn 2021", "Spring 2022",
                                                            "Summer 2022"))
par(mfrow = c(1, 2))
p7 <- ggplot(data = Light_rfr_dll, 
       aes(x = Light, y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.3,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
                position = position_dodge(0.3))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 2))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  ggtitle("dLL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("R:FR"),
                   labels = c("R:FR"))+
  labs(x="Light",
       y="R:FR")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#test grid
library(ggpubr) #needed for common legend
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2, bottom = legend)

#more complex layout
lay <- rbind(c(1,1,2,2,3,4),
             c(5,5,6,6,7,8))
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix = lay)

#try share legend
#graph with legend
ggp1_legend <- ggplot(data = Light_dll, 
             aes(x = factor(Light), y = Recording, fill = factor(Time)))+
  geom_bar(stat = "identity", width = 0.7,
           position = position_dodge(), col = "black")+
  geom_errorbar(aes(ymin=Recording-SD, ymax=Recording+SD), width = 0.2,
                position = position_dodge(0.7))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1500))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  ggtitle("dLL")+
  scale_fill_manual(values = c("pink2", "grey",
                               "orange", "cadetblue", "blue"))+
  scale_x_discrete(limits = c("PFD", "PPFD", "FR", "R", "G", "B", "UV"),
                   labels = c("PFD", "PPFD", "FR", "R", "G", "B", "UV"))+
  labs(x="Light",
       y="Light intensity (??mol m^-2~s^-1)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggp1_legend

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(ggp1_legend)

# Draw plots with shared legend
grid.arrange(arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4),
             shared_legend, nrow = 2, heights = c(10, 1))

#works with lay
#best so far
grid.arrange(arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix = lay),
             shared_legend, nrow = 2, heights = c(10, 1))

#save by giving it object name?
#g <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, nrow=2) #generates g

#best way to save, without legend here
tiff('Light_barpolots.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix = lay)
dev.off()

#save with legend YES - BEST SO FAR
#works with lay
tiff('Light_barpolots_leg.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
grid.arrange(arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, layout_matrix = lay),
             shared_legend, nrow = 2, heights = c(10, 1))
dev.off()
#low quality
#dev.copy(pdf,"whatever.pdf")
#dev.off()

#season significant differences
#summary data
gr <- read.csv("Light var for season significance.csv")
#raw data
gr <- read.csv("Light var for season significance raw.csv")
names(gr)
gr$Light <- gr$EnvLight
#PERCENTAGES
aov <- aov(B.~Time+Light,data=gr) #time
summary(aov)
#seasons no sig diff tukey
aov <- aov(R.~Time+Light,data=gr) #timeand light
summary(aov)
#seasons no sig diff tukey
aov <- aov(FR.~Time+Light,data=gr) #time and light
summary(aov)
#fr tukey seasonal summer autumn spring diffs, sping and autumn
#same and so is summer summer autumn autumn spring spring
aov <- aov(G.~Time+Light,data=gr) #time and light
summary(aov)
#spring and summ like aut. summ and spring diff
aov <- aov(PPFD.~Time+Light,data=gr) #time and light
summary(aov)
#seasons summ spring diff, sum aut diff between seasons
aov <- aov(UV.~Time+Light,data=gr) #light sig
summary(aov)

tukey <- TukeyHSD(aov, "Time", conf.level=.95)
tukey <- TukeyHSD(aov, "Light", conf.level=.95)
str(tukey)
tukey

#UV% high in Autumn 2021 when use summary

#RAW time and light sig
#summary light sig
aov <- aov(B~Time+Light,data=gr) #light sig
summary(aov)
#summ and spring diff to aut
aov <- aov(R~Time+Light,data=gr) #light sig
summary(aov)
#summ diff to aut
aov <- aov(FR~Time+Light,data=gr) #light sig
summary(aov)
#summ diff to aut
aov <- aov(R.FR~Time+Light,data=gr)#light sig
summary(aov)
#summ aut spring, summ 2022 diff to all
aov <- aov(G~Time+Light,data=gr) #light sig
summary(aov)
#aut diff to summ both and spr 2022
aov <- aov(PPFD~Time+Light,data=gr) #light sig
summary(aov)
#aut diff to bot summ and spr 2022
aov <- aov(UV~Time+Light,data=gr) #light sig
summary(aov)
#aut diff to both summ, spr 2022

tukey <- TukeyHSD(aov, "Time", conf.level=.95)
tukey <- TukeyHSD(aov, "Light", conf.level=.95)
str(tukey)
tukey

names(gr)
aov <- aov(Spr2021_LambdaD~Time+Light,data=gr) #light sig
summary(aov)
aov <- aov(Spr2021_LambdaP~Time+Light,data=gr) #both sig
summary(aov)
aov <- aov(PFD~Time+Light,data=gr) #both sig
summary(aov)


