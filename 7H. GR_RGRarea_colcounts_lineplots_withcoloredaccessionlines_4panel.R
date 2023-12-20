#growth rate rgr area and rgr col together as panel

#area prep
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

#RUN TO DO % coverage
percovavg1 <- pcov2 %>% group_by(Treatment, Accession) %>% summarise_all(mean, na.rm = TRUE) #add desc as argument to do descending

LL <- percovavg1 %>% filter(Treatment == "LL") %>% select(Accession, T0_area, T7_area, T14_area, T21_area, T28_area, T42_area)

#HL$Treatment <- NULL
#LL$Treatment <- NULL

#for percent cov
accessions <- (LL[,2])
sub <- (LL[,3:8]) #just t0 to t42 area columns to show all data
#to see hl and ll seperately need to split by treatment 
sa_comb <- sub

#sub <- (RGR_sum[,17:22]) #just t0 to t42 columns to show all data
#to see hl and ll seperately need to split by treatment 
#sa_comb <- sub
sa_comb <- LL
names(sa_comb)

#for percen cov
sa_comb <- sa_comb[-1] #remove treatment so just numeric
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

#running area graphs
#graphs to save
tiff('GR_Comb_lineplotpanels.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(2, 2),  mar=c(4,4.5,2,2))#run to make a panel of 4
#try to reduce spacing between plots
#redo plot for ll with 1 main color, colors for ones of interest
plot(KS02 ~ time, sadf, type = "n", 
     main = "LL",
     xlab = "", 
     ylab = "Percentage frond area coverage", 
     axes = F,
     xlim = c(0, 50),
     ylim = c(0, 100))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
title("A.", cex.main=2, adj=0)
for(i in 1:24){
  lines(sadf$time, sadf[, i], col = "#000000", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 2], col = "#FF0000", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 11], col = "#006400", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 16], col = "#006400", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 7], col = "#81007F", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 5], col = "#0000FF", lty = 1, lwd = 2)
  lines(sadf$time, sadf[, 6], col = "#0000FF", lty = 1, lwd = 2)
  
}
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf$time[-5], labels = sadf$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
#legend("bottomright", pt.cex = 150, cex = 0.3, col = "#00000033", lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()
#dev.off()

#prep for hl
HL <- percovavg1 %>% filter(Treatment == "HL") %>% select(Accession, T0_area, T7_area, T14_area, T21_area, T28_area, T42_area)
accessions <- (HL[,2])
sub <- (HL[,3:8])
sa_comb1 <- HL
names(sa_comb1)

#for percen cov
sa_comb1 <- sa_comb1[-1] #remove treatment so just numeric
#sa_comb <- sa_comb[-1]
#sa_comb <- HL[-1]
#sa_comb <- sa_comb[-1]

sa_comb1$Accession <- NULL

#row.names(sa_comb) <- sa_comb$Accession
#sa_comb$Accession <- NULL
#take col names from first col and give new names
sa_comb1 <- t(sa_comb1) #become a matrix during transposition
sa_comb1
#turn back into dataframe again
sadf1 <- data.frame(sa_comb1)
class(sadf1)

row.names(sadf1)

#row names defined col
#null col 

time_vec <- numeric(length = nrow(sadf1)) #make vector same name as row names

library(stringr)

#for working with percen area
#replace characters in a string
row.names(sadf1) <- str_replace_all(row.names(sadf1), "T0_", "")
row.names(sadf1) <- str_replace_all(row.names(sadf1), "_area", "")
row.names(sadf1) <- str_replace_all(row.names(sadf1), "T", "")
row.names(sadf1) <- str_replace_all(row.names(sadf1), "area", "0")
#looks for pattern and give replacement

row.names(sadf1)

time_vec <- row.names(sadf1)
time_vec <- as.numeric(time_vec) #coerce to numeric

class(sadf1) 
accessions <- t(accessions)
colnames(sadf1) <- (accessions)

sadf1$time <- time_vec

str(sadf1)

class(sadf1)

library(tidyr)

#tiff('Growthrate_HL_lineplot_cov.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
#redo plot for hl with 1 color
plot(KS02 ~ time, sadf1, type = "n", 
     main = "HL",
     xlab = "", 
     ylab = "", 
     axes = F,
     xlim = c(0, 50),
     ylim = c(0, 100))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
#text(15.8,23.8, "RGRlog", srt=45) #place txt and rotate it
for(i in 1:24){
  lines(sadf1$time, sadf1[, i], col = "#000000", lty = 1, lwd = 2)
  lines(sadf1$time, sadf1[, 2], col = "#FF0000", lty = 1, lwd = 2)
  lines(sadf1$time, sadf1[, 11], col = "#006400", lty = 1, lwd = 2)
  lines(sadf1$time, sadf1[, 16], col = "#006400", lty = 1, lwd = 2)
  lines(sadf1$time, sadf1[, 7], col = "#81007F", lty = 1, lwd = 2)
  lines(sadf1$time, sadf1[, 5], col = "#0000FF", lty = 1, lwd = 2)  
  lines(sadf1$time, sadf1[, 6], col = "#0000FF", lty = 1, lwd = 2)
  
  }
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf1$time[-5], labels = sadf1$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
#legend("bottomright", pt.cex = 150, cex = 0.3, col = "#00000033", lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()
#dev.off()

#col counts prep
raw <- read.csv("Col_countsRGR_nospREM.csv")

RGR_new <- raw
RGR_sum <- RGR_new %>% group_by(Treatment, Accession) %>% summarise_all(mean, na.rm = TRUE)
#percencov
LL2 <- RGR_sum %>% filter(Treatment == "LL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)

sa_comb3 <- LL2
accessions2 <- (LL2[,2])
names(sa_comb3)

sa_comb3$Accession <- NULL
sa_comb3$Treatment <- NULL

#take col names from first col and give new names
sa_comb3 <- t(sa_comb3) #become a matrix during transposition
sa_comb3
#turn back into dataframe again
sadf3 <- data.frame(sa_comb3)
class(sadf3)

row.names(sadf3)

#row names defined col
#null col 

time_vec <- numeric(length = nrow(sadf3)) #make vector same name as row names

library(stringr)

#replace characters in a string
row.names(sadf3) <- str_replace_all(row.names(sadf3), "T0_", "0")
row.names(sadf3) <- str_replace_all(row.names(sadf3), "_cols", "")
row.names(sadf3) <- str_replace_all(row.names(sadf3), "T", "")
row.names(sadf3) <- str_replace_all(row.names(sadf3), "cols", "")
#looks for pattern and give replacement

row.names(sadf3)

time_vec <- row.names(sadf3)
time_vec <- as.numeric(time_vec) #coerce to numeric

class(sadf3) 
accessions2 <- t(accessions2)
colnames(sadf3) <- (accessions2)

sadf3$time <- time_vec

str(sadf3)

class(sadf3)

library(tidyr)


#col counts graphs
#redo plot for ll with 1 color
plot(KS02 ~ time, sadf3, type = "n", 
     main = "LL",
     xlab = "Time (days)", 
     ylab = "Colony gain per day", 
     axes = F,
     xlim = c(0, 21),
     ylim = c(0, 120))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
title("B.", cex.main=2, adj=0)
for(i in 1:24){
  lines(sadf3$time, sadf3[, i], col = "#000000", lty = 1, lwd = 2)
  lines(sadf3$time, sadf3[, 2], col = "#FF0000", lty = 1, lwd = 2)
  lines(sadf3$time, sadf3[, 11], col = "#006400", lty = 1, lwd = 2)
  lines(sadf3$time, sadf3[, 16], col = "#006400", lty = 1, lwd = 2)
  lines(sadf3$time, sadf3[, 7], col = "#81007F", lty = 1, lwd = 2)
  lines(sadf3$time, sadf3[, 5], col = "#0000FF", lty = 1, lwd = 2)
  lines(sadf3$time, sadf3[, 6], col = "#0000FF", lty = 1, lwd = 2)
  
  }
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf3$time[-5], labels = sadf3$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
#legend("bottomright", pt.cex = 150, cex = 0.3, col = "#00000033", lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()

#prep for hl

HL2 <- RGR_sum %>% filter(Treatment == "HL") %>% select(Accession, T0_cols, T7_cols, T14_cols, T21_cols)
sa_comb4 <- HL2
accessions2 <- (HL2[,2])
sa_comb4$Accession <- NULL
sa_comb4$Treatment <- NULL
#row.names(sa_comb) <- sa_comb$Accession
#sa_comb$Accession <- NULL
#take col names from first col and give new names
sa_comb4 <- t(sa_comb4) #become a matrix during transposition
sa_comb4
#turn back into dataframe again
sadf4 <- data.frame(sa_comb4)
class(sadf4)

row.names(sadf4)

#row names defined col
#null col 

time_vec <- numeric(length = nrow(sadf4)) #make vector same name as row names

library(stringr)

#replace characters in a string
row.names(sadf4) <- str_replace_all(row.names(sadf4), "T0_", "0")
row.names(sadf4) <- str_replace_all(row.names(sadf4), "_cols", "")
row.names(sadf4) <- str_replace_all(row.names(sadf4), "T", "")
row.names(sadf4) <- str_replace_all(row.names(sadf4), "cols", "")
#looks for pattern and give replacement

row.names(sadf4)

time_vec <- row.names(sadf4)
time_vec <- as.numeric(time_vec) #coerce to numeric

class(sadf4) 
accessions2 <- t(accessions2)
colnames(sadf4) <- (accessions2)

sadf4$time <- time_vec

str(sadf4)

class(sadf4)

library(tidyr)

#redo plot for hl with 1 color
plot(KS02 ~ time, sadf4, type = "n", 
     main = "HL",
     xlab = "Time (days)", 
     ylab = "", 
     axes = F,
     xlim = c(0, 21),
     ylim = c(0, 120))
rect(xleft = 14, ybottom = -50, xright = 21, ytop = 5000, col ="#00000033", border = "#FFFFFFFF")
for(i in 1:24){
  lines(sadf4$time, sadf4[, i], col = "#000000", lty = 1, lwd = 2)
  lines(sadf4$time, sadf4[, 2], col = "#FF0000", lty = 1, lwd = 2)
  lines(sadf4$time, sadf4[, 11], col = "#006400", lty = 1, lwd = 2)
  lines(sadf4$time, sadf4[, 16], col = "#006400", lty = 1, lwd = 2)
  lines(sadf4$time, sadf4[, 7], col = "#81007F", lty = 1, lwd = 2)
  lines(sadf4$time, sadf4[, 5], col = "#0000FF", lty = 1, lwd = 2)
  lines(sadf4$time, sadf4[, 6], col = "#0000FF", lty = 1, lwd = 2)
  
  }
#axis(1, las = 1) hash out so doesnt add x labels as defined
axis(1, at = sadf4$time[-5], labels = sadf4$time[-5]) #at is where ticks going, labels printed
axis(2, las = 1)
#legend("bottomright", pt.cex = 150, cex = 0.3, col = "#00000033", lty = sort(line_type), lwd = 2, legend = names(sadf[-ncol(sadf)])) #include all except final column
box()
dev.off()
#dev.new()

