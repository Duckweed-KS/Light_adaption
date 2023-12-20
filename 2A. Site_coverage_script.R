#script to explore enviro coverage data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Site_coverage_visit.csv")
#doesnt work with dates
#gr <- read.csv("Site_coverage_date.csv")

tail(gr)
head(gr)
names(gr)

#check classes of all cols
sapply(gr, class)

#gr$Visit <- as.factor(gr$Visit)

library(dplyr)
library(stringr)
library(ggplot2)

boxplot(Avg.cov~Site,data=gr)
boxplot(Avg.cov~Site*Visit,data=gr)
boxplot(Avg.cov~Site*Date,data=gr)

#need to split observations and show over time
#plot as line graph

sa_comb <- gr #all HL and gives accession column back
names(sa_comb)

sa_comb$Site <-NULL
#sa_comb$Date <- NULL
#sa_comb$SD.cov <- NULL
#sa_comb$SE.cov <- NULL

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
row.names(sadf) <- str_replace_all(row.names(sadf), "X", "")
#looks for pattern and give replacement
#check its worked
row.names(sadf)

time_vec <- row.names(sadf)
time_vec <- as.numeric(time_vec) #coerce to numeric

class(sadf) 

accessions <- (gr[,1]) #define Site variable
accessions <- t(accessions) #make transverse so correct no of cols
colnames(sadf) <- (accessions) #put names back on matrix as col names

sadf$time <- time_vec #make time column

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
par(mfrow = c(1, 1))
plot(KS02 ~ time, sadf, type = "n", 
     main = "Duckweed coverage at UK sites",
     xlab = "Visit", 
     ylab = expression(paste("% Site Coverage")), 
     axes = F,
     xlim = c(0, 9),
     ylim = c(0, 100))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:24){
  lines(sadf$time, sadf[, i], col = line_colours[i], lty = sort(line_type)[i], lwd = 2)
  
}
axis(1, las = 1) #hash out so doesnt add x labels as defined
#axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
legend("bottomright", pt.cex = 150, cex = 0.3, col = line_colours, lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()

line_colours <- rep(1:24, length.out = ncol(sadf)-1) #8 colors
line_colours
line_type <- rep(1:3, length.out =  ncol(sadf)-1)
line_type

gr <- read.csv("Site_coverage.csv")
#doesnt work with dates
#gr <- read.csv("Site_coverage_date.csv")

tail(gr)
head(gr)
names(gr)

#check classes of all cols
sapply(gr, class)

#gr$Visit <- as.factor(gr$Visit)

library(dplyr)
library(stringr)
library(ggplot2)

boxplot(Avg.cov~Site,data=gr)
#13, 17, 20, 21, 27, 28, 29 always low
boxplot(Avg.cov~Site*Visit,data=gr)
boxplot(Avg.cov~Date,data=gr)
#mar 21 lowest avgs, jul 22 highest
boxplot(Avg.cov~Visit,data=gr)

aov <- aov(Avg.cov~Site*Visit,data=gr) #sig accession
summary(aov)
tukey <- TukeyHSD(aov, "Site", conf.level=.95) 
tukey <- TukeyHSD(aov, "Visit", conf.level=.95) 
str(tukey)
tukey

gr <- gr[complete.cases(gr$Avg.cov), ]   

Jun20_ord <- gr %>% filter(Date == "Jun-20") %>% arrange(desc(Avg.cov))
Jul20_ord <- gr %>% filter(Date == "Jul-20") %>% arrange(desc(Avg.cov))
Oct20_ord <- gr %>% filter(Date == "Oct-20") %>% arrange(desc(Avg.cov))
Mar21_ord <- gr %>% filter(Date == "Mar-21") %>% arrange(desc(Avg.cov))
Jul21_ord <- gr %>% filter(Date == "Jul-21") %>% arrange(desc(Avg.cov))
Oct21_ord <- gr %>% filter(Date == "Oct-21") %>% arrange(desc(Avg.cov))
Mar22_ord <- gr %>% filter(Date == "Mar-22") %>% arrange(desc(Avg.cov))
Jul22_ord <- gr %>% filter(Date == "Jul-22") %>% arrange(desc(Avg.cov))
gr %>% group_by(Date) %>% summarise(mean_pop = mean(Avg.cov))
gr %>% group_by(Date) %>% summarise(min_pop = min(Avg.cov))
gr %>% group_by(Date) %>% summarise(max_pop = max(Avg.cov))

#cbind not working, dif number observations
rankings <- merge(Jun20_ord, Jul20_ord, Oct20_ord, Mar21_ord, Jul21_ord,
                  Oct21_ord, Mar22_ord, Jul22_ord, by="Site", all = T)
