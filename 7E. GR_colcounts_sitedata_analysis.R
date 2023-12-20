#combine GR and Col count data
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
#grouping variables by hl or ll

library(dplyr)
library(tidyr)

#read in averaged data old versions
#gr <- read.csv("GR_RGR_propchange.csv")
#cc <- read.csv("Col_RGR_propchange.csv")
#bio <- read.csv("Biomass_RGR_propchange.csv")

#new versions
gr <- read.csv("GR_RGR_propchange_oppcalc_spREM.csv")
#nb doesnt include 4 accessions
cc <- read.csv("Col_RGR__oppcalc_propchange_nospREM.csv")
#nb doesnt include 4 accessions
bio <- read.csv("Biomass_RGR_propchange_oppcalc.csv")
#includes 4 accessions tot 28

#join 2 data sets together (by accession)
names(gr)

joined <- gr %>% inner_join(cc, by = c("Accession"))

#ADD IN WHEN WORKING WITH BIO DATA
#new add
joined <- joined %>% inner_join(bio, by = c("Accession"))

#measured these before added bio
#relationships between col gain and area light responses
#both worked out using LL - HL 
names(joined)

plot(RGR_diff.x ~ RGR_diff.y, joined)
cor(joined$RGR_diff.x, joined$RGR_diff.y,
    method = "pearson"
    )
#0.74
#when diff area increases, diff cols increases

#correlation between FW and FDW
corr <- cor.test(joined$FW_mean_diff, joined$FDW_mean_diff,
                 method = "pearson"
)
corr$estimate
corr$p.value
#0.3 p value 0.15

#correlation area and FW
corr <- cor.test(joined$FW_mean_diff, joined$RGR_diff.x,
                 exact = FALSE,
                 method = "pearson"
)
corr$estimate
corr$p.value
#0.63 p value 0.00009

#correlation cols and FW
corr <- cor.test(joined$FW_mean_diff, joined$RGR_diff.y,
                 exact = FALSE,
                 method = "pearson"
)
corr$estimate
corr$p.value
#0.58 p value 0.002

#correlation area and FDW
corr <- cor.test(joined$FDW_mean_diff, joined$RGR_diff.x,
                 exact = FALSE,
                 method = "pearson"
)
corr$estimate
corr$p.value
#0.39 p value 0.18

#correlation cols and FDW
corr <- cor.test(joined$FDW_mean_diff, joined$RGR_diff.y,
                 exact = FALSE,
                 method = "pearson"
)
corr$estimate
corr$p.value
#-0.06 p value 0.7 NO RELATIONSHIP

#area and cols
corr <- cor.test(joined$RGR_diff.x, joined$RGR_diff.y,
                 exact = FALSE,
                 method = "pearson"
                 )
corr$estimate
corr$p.value
#0.74 p value <0.00001

#area and cols correlate
#area and FW correlate
#cols and fw

#weaker
#fdw and area
#fdw and cols
#fdw and fd

corr <- cor.test(joined$mean.HL, joined$ColRGR_HL_mean,
                 method = "pearson")
corr$estimate
corr$p.value
#weak 047 and 0.01
#higher mean area in HL, higher cols in HL

corr <- cor.test(joined$mean.LL, joined$ColRGR_LL_mean,
                 method = "pearson")
corr$estimate
corr$p.value
#stronger 0.59 0.002
#higher mean area in LL, higher cols in LL

corr <- cor.test(joined$mean.HL, joined$mean.LL,
                 method = "pearson")
corr$estimate
corr$p.value
#strong 0.62, 0.001

corr <- cor.test(joined$ColRGR_HL_mean, joined$ColRGR_LL_mean,
                 method = "pearson")
corr$estimate
corr$p.value
#weak 0.34, 0.09

corr <- cor.test(joined$RGR_diff.x, joined$RGR_diff.y,
                 method = "pearson")
corr$estimate
corr$p.value
#weak 0.74, <0.0001

library(corrplot)
#basic corr plot, first check class, select only numeric, then plot matrix out
sapply(joined, class)
sapply(joined, is.numeric)
num <- select_if(joined, is.numeric)
cornum <- cor(num) #displays r2

par(mfrow = c(1, 1))
corrplot(cor(num),
         method = "number",
         type = "upper" # show only upper side
)

names(joined)
plot(mean.HL ~ mean.HL, joined)
plot(mean.LL ~ mean.LL, joined) #no real corr
plot(RGR_diff.x ~ RGR_diff.y, joined) # looks good
plot(maximum.HL ~ RGR_HL_maximum, joined)
plot(maximum.LL ~ RGR_LL_maximum, joined)

#change names of cols and remove excess
#joined$Accession <- joined$Accession1
joined$X.x <- NULL
joined$X.y <- NULL
joined$X <- NULL
joined$X.1 <- NULL
joined$Treatment <- NULL

names(joined)


library(ggplot2)
library(ggpubr)
#plot scatter plots with labels
#NICER PLOTS FOR SIG AREA WITH COL GAIN AND AREA WITH FW
#difference in area plot difference in colonies

#for paper
myggp <- ggplot(joined, aes(RGR_diff.x, RGR_diff.y, label = Accession)) +    # ggplot2 plot with labels
    geom_point(pch=3) +
    geom_smooth(method='lm') +
    xlab("RGR diff area (mm2)") +
    ylab("RGR diff col gain") +
    theme_bw() +
    theme_classic() +
    stat_cor(method = "pearson", label.x = 0.08, label.y = 0.5)
myggp + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))

#difference in area plot difference in colonies
ggplot(joined, aes(RGR_diff.x, FW_mean_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("RGR diff area (mm2)") +
  ylab("RGR diff fresh weight (mg)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.25, label.y = 2000)

#mean HL area vs mean HL cols
ggplot(joined, aes(mean.HL.x, mean.HL.y, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#max HL area vs mean HL cols
ggplot(joined, aes(maximum.HL.x, maximum.HL.y, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#area vs FW diff
ggplot(joined, aes(RGR_diff.x, FW_mean_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#col vs FW diff
ggplot(joined, aes(RGR_diff.y, FW_mean_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#area vs FDW diff
ggplot(joined, aes(RGR_diff.x, FDW_mean_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#col vs FDW diff
ggplot(joined, aes(RGR_diff.y, FDW_mean_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#not useful
#ggplot(joined, aes(propchange4, RGR_propdiff, label = Accession)) +    # ggplot2 plot with labels
#  geom_point() +
#  geom_text(aes(label = Accession), hjust = - 0.5)

#stronger treatment effect between methods if same accessions top/bottom

#FOR RANKING INDIVIDUALS FOR EACH VARIABLE
names(joined)

#x is area
#y is col counts

#BY HL
mean.HL.x_ord <- joined %>% arrange(desc(mean.HL)) %>% select(Accession, mean.HL) %>% top_n(28)
maximum.HL.x_ord <- joined %>% arrange(desc(maximum.HL)) %>% select(Accession, maximum.HL) %>% top_n(28)
mean.HL.y_ord <- joined %>% arrange(desc(ColRGR_HL_mean)) %>% select(Accession, ColRGR_HL_mean) %>% top_n(28)
maximum.HL.y_ord <- joined %>% arrange(desc(RGR_HL_maximum)) %>% select(Accession, RGR_HL_maximum) %>% top_n(28)
mean.HL.FW_ord <- joined %>% arrange(desc(FW.mean.HL)) %>% select(Accession, FW.mean.HL) %>% top_n(28)
maximum.HL.FW_ord <- joined %>% arrange(desc(FW.maximum.HL)) %>% select(Accession, FW.maximum.HL) %>% top_n(28)
mean.HL.FDW_ord <- joined %>% arrange(desc(FDW.mean.HL)) %>% select(Accession, FDW.mean.HL) %>% top_n(28)
maximum.HL.FDW_ord <- joined %>% arrange(desc(FDW.maximum.HL)) %>% select(Accession, FDW.maximum.HL) %>% top_n(28)

#BY LL
mean.LL.x_ord <- joined %>% arrange(desc(mean.LL)) %>% select(Accession, mean.LL) %>% top_n(28)
maximum.LL.x_ord <- joined %>% arrange(desc(maximum.LL)) %>% select(Accession, maximum.LL) %>% top_n(28)
mean.LL.y_ord <- joined %>% arrange(desc(ColRGR_LL_mean)) %>% select(Accession, ColRGR_LL_mean) %>% top_n(28)
maximum.LL.y_ord <- joined %>% arrange(desc(RGR_LL_maximum)) %>% select(Accession, RGR_LL_maximum) %>% top_n(28)
mean.LL.FW_ord <- joined %>% arrange(desc(FW.mean.LL)) %>% select(Accession, FW.mean.LL) %>% top_n(28)
maximum.LL.FW_ord <- joined %>% arrange(desc(FW.maximum.LL)) %>% select(Accession, FW.maximum.LL) %>% top_n(28)
mean.LL.FDW_ord <- joined %>% arrange(desc(FDW.mean.LL)) %>% select(Accession, FDW.mean.LL) %>% top_n(28)
maximum.LL.FDW_ord <- joined %>% arrange(desc(FDW.maximum.LL)) %>% select(Accession, FDW.maximum.LL) %>% top_n(28)

#by differences
RGR_diff.x_ord <- joined %>% arrange(desc(RGR_diff.x)) %>% select(Accession, RGR_diff.x) %>% top_n(28)
max.diff.x_ord <- joined %>% arrange(desc(max_diff.x)) %>% select(Accession, max_diff.x) %>% top_n(28)
propchange4.x_ord <- joined %>% arrange(desc(propchange4)) %>% select(Accession, propchange4) %>% top_n(28)
RGR_diff.y_ord <- joined %>% arrange(desc(RGR_diff.y)) %>% select(Accession, RGR_diff.y) %>% top_n(28)
max.diff.y_ord <- joined %>% arrange(desc(max_diff.y)) %>% select(Accession, max_diff.y) %>% top_n(28)
RGR_propdiff_ord <- joined %>% arrange(desc(RGR_propdiff)) %>% select(Accession, RGR_propdiff) %>% top_n(28)
FW_mean_diff_ord <- joined %>% arrange(desc(FW_mean_diff)) %>% select(Accession, FW_mean_diff) %>% top_n(28)
FW_max_diff_ord <- joined %>% arrange(desc(FW_max_diff)) %>% select(Accession, FW_max_diff) %>% top_n(28)
FDW_mean_diff_ord <- joined %>% arrange(desc(FDW_mean_diff)) %>% select(Accession, FDW_mean_diff) %>% top_n(28)
FDW_max_diff_ord <- joined %>% arrange(desc(FDW_max_diff)) %>% select(Accession, FDW_max_diff) %>% top_n(28)
FW_prop_diff_ord <- joined %>% arrange(desc(FW_propdiff)) %>% select(Accession, FW_propdiff) %>% top_n(28)
FDW_prop_diff_ord <- joined %>% arrange(desc(FDW_propdiff)) %>% select(Accession, FDW_propdiff) %>% top_n(28)

#add density variabl
#density = divide mass per area
density <- joined %>% mutate(adensity_HL = (FW.mean.HL/mean.HL))
density %>% arrange(desc(adensity_HL)) %>% select(Accession) %>% top_n(28)
density <- density %>% mutate(adensity_LL = (FW.mean.LL/mean.LL))
density %>% arrange(desc(adensity_LL)) %>% select(Accession) %>% top_n(28)

#density = divide mass per col
density <- density %>% mutate(Cdensity_HL = (FW.mean.HL/ColRGR_HL_mean))
density %>% arrange(desc(Cdensity_HL)) %>% select(Accession) %>% top_n(28)
density <- density %>% mutate(Cdensity_LL = (FW.mean.LL/ColRGR_LL_mean))
density %>% arrange(desc(Cdensity_LL)) %>% select(Accession) %>% top_n(28)

#combo
adensity_HL_ord <- density %>% arrange(desc(adensity_HL)) %>% select(Accession, adensity_HL) %>% top_n(28)
adensity_LL_ord <- density %>% arrange(desc(adensity_LL)) %>% select(Accession, adensity_LL) %>% top_n(28)
Cdensity_HL_ord <- density %>% arrange(desc(Cdensity_HL)) %>% select(Accession, Cdensity_HL) %>% top_n(28)
Cdensity_LL_ord <- density %>% arrange(desc(Cdensity_LL)) %>% select(Accession, Cdensity_LL) %>% top_n(28)


#stick all together with rankings next to them - see consistency
class(mean.HL.x_ord) #df
rankings_hl <- cbind(mean.HL.x_ord,maximum.HL.x_ord,mean.HL.y_ord,maximum.HL.y_ord,
                  mean.HL.FW_ord,maximum.HL.FW_ord,mean.HL.FDW_ord,maximum.HL.FDW_ord,
                  adensity_HL_ord,Cdensity_HL_ord)
                 
rankings_ll <- cbind(mean.LL.x_ord,maximum.LL.x_ord,mean.LL.y_ord,maximum.LL.y_ord,
                     mean.LL.FW_ord,maximum.LL.FW_ord,mean.LL.FDW_ord,maximum.LL.FDW_ord,
                     adensity_LL_ord,Cdensity_LL_ord)

#by differences
rankings_diff <- cbind(RGR_diff.x_ord,max.diff.x_ord,propchange4.x_ord,
                       RGR_diff.y_ord,max.diff.y_ord,RGR_propdiff_ord,
                       FW_mean_diff_ord, FW_max_diff_ord,FDW_mean_diff_ord,
                       FDW_max_diff_ord,FW_prop_diff_ord,FDW_prop_diff_ord)

rankings <- cbind(rankings_hl, rankings_ll, rankings_diff)

write.csv(rankings, "GR_HL_LL_diff_rankings.csv")

class(rankings)

#colnames(rankings) <- c("mean.HL.x_ord", "mean.LL.x_ord", "maximum.HL.x_ord", "maximum.LL.x_ord",
#                            "RGR_diff.x_ord", "max.diff.x_ord", "propchange4.x_ord", "mean.HL.y_ord", 
#                            "mean.LL.y_ord", "maximum.HL.y_ord", "maximum.LL.y_ord",
#                            "RGR_diff.y_ord", "max.diff.y_ord", "RGR_propdiff_ord",
#                        "FW_mean_diff_ord", "FW_max_diff_ord", "FDW_mean_diff_ord", "FDW_max_diff_ord",
#                        "adensity_HL_ord", "adensity_LL_ord", "Cdensity_HL_ord", "Cdensity_LL_ord")
#need to add cols for 1-28 between each column

#density = divide mass per area
names(joined)
#density upgraded version of joined - extra variables to correlate

#density (biomass per col) in hl and ll
ggplot(density, aes(Cdensity_HL, Cdensity_LL, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

ggplot(density, aes(adensity_HL, adensity_LL, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

names(density)
density <- density %>% mutate(SLA_HL = (FW.mean.HL/ColRGR_HL_mean))
density %>% arrange(desc(SLA_HL)) %>% select(Accession) %>% top_n(28)
density <- density %>% mutate(SLA_LL = (FW.mean.LL/ColRGR_LL_mean))
density %>% arrange(desc(SLA_LL)) %>% select(Accession) %>% top_n(28)

#areas divided by colonoies
ggplot(density, aes(SLA_HL, SLA_LL, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

library(dplyr)
density <- density%>% mutate(SLA_propdiff = (SLA_HL/SLA_LL)*100)
density <- density%>% mutate(adensity_propdiff = (adensity_HL/adensity_LL)*100)

#correlate density data
#dat <- density
dat <- joined


str(dat)
numeric_cols <- lapply(dat, is.numeric) %>% unlist() %>% unname()
dat2 <- dat[, numeric_cols]
library(corrplot)
?corrplot
cor(dat2)

dat2 <- dat2[complete.cases(dat2), ]

#create matrix

M<-cor(dat2,use='complete.obs')
M[M < 0.5 & M > -0.5] <- 0


col_classes <- unname(sapply(M, class)) #do I need to do this to define col_classes to then be able to do something else with this? where is 'class' coming from?
Msig <- which(col_classes == 0.5:-0.5)
replace(Msig, NA)

#Run correlation, 'complete.obs' allows for na values:
#matr <- read.csv("physio_roots_ionome_sites_num.csv", stringsAsFactors=FALSE)
#str(matr)
#head(matr)
#tail(matr)
#M<-cor(matr,use='complete.obs')

#M<-cor(paramnum,use='complete.obs')

library(dplyr)
library(corrplot)

### Computing p-value of correlations. mat : is a matrix of data.
#This is creating a function to compute p-values, you just need to run this.
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Create a matrix of the p-value of the correlation:
p.mat <- cor.mtest(dat2)


plotcolour<- colorRampPalette(c("darkorange2", "gray88", "darkcyan"))(20) 

corrplot(M, type="upper", order="hclust",method = "number",
         p.mat = p.mat, sig.level = 0.1005, insig = "blank",tl.cex = 0.5, tl.col = 'black',col=plotcolour)

corrplot(M, col=plotcolour, type="upper", order="hclust", # Add coefficient of correlation
         tl.col="black", tl.cex = 0.5, tl.offset = 0.5,
         #Text label color and rotation
         # Combine with significance
         p.mat = p.mat,
         sig.level = c(.001, .01, .05), pch.cex = 0.8, #add signficance stars
         insig = "label_sig", pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

plotcolour<- colorRampPalette(c("darkorange2", "gray88", "darkcyan"))(20)
#This is where I've created a colour scheme of my choosing. 

write.csv(M, "GR_col_bio_mat.csv")

#confirms RGR_x and RGR_y (area and col)
#confirms RGR_x and FW_mean

library(basicTrendline)

#nicer plot for correlation
#need to add labels, pvalue, rsquared, change scales
names(joined)
RGR_diff_area <- (joined[ ,10])
RGR_diff_cols <- (joined[ ,28])
par(pty = "s") #"m" normal #square plots
plot(RGR_diff_area,RGR_diff_cols)
trendline(RGR_diff_area, RGR_diff_cols, model = "line2P", Pvalue.corrected = TRUE, 
          linecolor = "black", lty = 1, lwd = 1, show.equation = TRUE,
          show.Rpvalue = TRUE, Rname = 1, Pname = 1, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, ePos.x = 0.3, ePos.y = -400,
          text.col = "black", eDigit = 3, eSize = 1, CI.fill = TRUE,
          CI.level = 0.95, CI.color = "lightblue", CI.alpha = 1, CI.lty = 1,
          CI.lwd = 1, las = 1, main = "Area gain vs Col gain", cex.main=1,
          xlab = "log RGR difference HL/LL area (mm2)", 
          ylab = "log RGR difference HL/LL col gain", cex.lab=0.9, cex.axis=0.8)

RGR_diff_area <- (joined[ ,10])
FW_mean_diff <- (joined[ ,57])
par(pty = "s") #"m" normal #square plots
plot(RGR_diff_area,FW_mean_diff)
trendline(RGR_diff_area, FW_mean_diff, model = "line2P", Pvalue.corrected = TRUE, 
          linecolor = "red", lty = 1, lwd = 1, show.equation = TRUE,
          show.Rpvalue = TRUE, Rname = 1, Pname = 1, xname = "x", yname = "y",
          yhat = FALSE, summary = TRUE, ePos.x = 0.3, ePos.y = -400,
          text.col = "black", eDigit = 3, eSize = 1, CI.fill = TRUE,
          CI.level = 0.95, CI.color = "grey", CI.alpha = 1, CI.lty = 1,
          CI.lwd = 1, las = 1, main = "Area gain vs Fresh weight", cex.main=1,
          xlab = "log RGR difference HL/LL area (mm2)", 
          ylab = "Biomass difference HL/LL (mg)", cex.lab=0.9, cex.axis=0.8)

#resave so can do sapply func
#write.csv(density, "Growth_combined.csv")
write.csv(density, "Growth_combined_oppcalc.csv")
#manual rename
new <- read.csv("Growth_combined_renamed.csv") #old
new <- read.csv("Growth_combined_oppcalc_renamed.csv") #new
names(new)

lapply(new, class)

#add other files to
cov <- read.csv("Site_coverage_date.csv")
names(cov)
Seasonavg <- colMeans(cov[, 2:8], na.rm = TRUE) #highest cov avg in July 21
colMeans(cov[, 2:4], na.rm = TRUE) #mean 2020
colMeans(cov[, 5:7], na.rm = TRUE) #mean 2021
#colSums(cov[, 2:8], na.rm = TRUE)
Monthavg <- rowMeans(cov[, 2:8], na.rm = TRUE) #avg coverage for each accession
rowSums(cov[, 2:8], na.rm = TRUE)

cov$Accession

cov <- cbind(cov, Monthavg)
#cov <- rbind(cov, Seasonavg)

#write.csv(cov, "Site_coverage_date_summ.csv")

#run from here, rest joining
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
library(dplyr)
library(corrplot)
#grouping variables by hl or ll
#add coverage and light data to growth rate
#new <- read.csv("Growth_combined_renamed.csv") #old
new <- read.csv("Growth_combined_oppcalc_reduced+site+light+altitude.csv") #use this?
new <- read.csv("Growth_combined_oppcalc_renamed.csv") #new RENAME
#new has new calculations for rgr and col gain diffs + proportions

cov <- read.csv("Site_coverage_date_summ.csv")
env <- read.csv("Env-light_accessnames.csv")
lag <- read.csv("RGRlag.csv") #to change

lag$X <- NULL

#still need to avg fluo and chl data to match hs param avged
fluo <- read.csv("Fluo_Summary_sepHL_LL.csv")
#no moor sel etc
hs <- read.csv("HS_Summ_HL_LLseperated.csv")
chl <- read.csv("Chloro_Summary_sepHLLL.csv")

hs$X <- NULL
hs$Treatment <- NULL
hs$Treatment9 <- NULL
hs$Accession9 <- NULL

#env summaries
colMeans(env[, 2:46], na.rm = TRUE) #avg all variables
rowMeans(env[, 2:46], na.rm = TRUE) #avg all accessions for all variables
#avgs not that useful

#28 observations
comb <- new %>% inner_join(cov, by = c("Accession"))
comb <- comb %>% inner_join(lag, by = c("Accession"))
comb <- comb %>% inner_join(fluo, by = c("Accession"))
#fluo has 66a 77a 78a no others do at mo
#reduce observations to 24
comb <- comb %>% inner_join(env, by = c("Accession"))
#to do #reduce observations to 20
#has 24 if dont include env
comb <- comb %>% inner_join(hs, by = c("Accession"))
comb <- comb %>% inner_join(chl, by = c("Accession"))

#summarise coverage data
colMeans(comb[, 63:65], na.rm = TRUE) #mean 2020, 27, 31, 34%
colMeans(comb[, 66:68], na.rm = TRUE) #mean 2021 7, 44, 36%
colMeans(comb[, 69:70], na.rm = TRUE) #mean 2022 20, 43%
#colSums(comb[, 2:8], na.rm = TRUE)
rowMeans(comb[, 63:70], na.rm = TRUE) #avg comberage for each accession
rowSums(comb[, 63:70], na.rm = TRUE)
#highest overall coverage ks03, maintains 67% year and season round, followed
#by 6a. lowest = 18, 15, 14, = KS27, KS21, KS20. 
comb$Accession

#boxplots for accessions changes over seasons
boxplot(Jun.20~Accession,data=comb) #ks06 highest
boxplot(Jul.20~Accession,data=comb) #ks03 highest
boxplot(Oct.20~Accession,data=comb) #ks18 highest
boxplot(Mar.21~Accession,data=comb) 
#KS03 60% march, others mostly 0, 30%, 20%, 5%. cold tolerance
boxplot(Jul.21~Accession,data=comb) #ks03, ks04, ks12 highest
boxplot(Oct.21~Accession,data=comb) #ks25 highest

boxplot(Mar.22~Accession,data=comb) 
#more with higher % this year, KS02 and KS03 highest
boxplot(Jul.22~Accession,data=comb) 
#KS18 highest

#need to group coverage data into seasons to show var as boxplot?
#subset <- comb[c(2,63:70)]
#SEE SITE COVERAGE SCRIPT FOR MORE ANALYSIS OF THIS

#has all RGR groupings + environment light data + coverage data

comb$Accession

#add grouping for light levels FOR 28 observations
comb$EnvLight <- c("HL", "LL", "LL", "HL", "LL", "LL", "LL",
                   "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                   "HL", "HL", "LL", "LL", "LL", "LL", "HL", "HL",
                   "HL", "LL", "LL", "HL", "LL", "HL")

#add grouping for light levels FOR 24 observations
comb$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL", "HL",
                   "HL", "HL", "LL", "HL", "LL", "LL", "HL",
                   "HL", "LL", "LL", "LL", "LL", "HL", "HL", "HL",
                   "LL", "LL")

#add grouping for light levels FOR 24 observations rename
comb$EnvLight <- c("dLL", "dLL", "dHL", "dLL", "dLL", "dLL", "dHL",
                   "dHL", "dHL", "dLL", "dHL", "dLL", "dLL", "dHL",
                   "dHL", "dLL", "dLL", "dLL", "dLL", "dHL", "dHL", "dHL",
                   "dLL", "dLL")

#add grouping for species FOR 24 observations rename
comb$Species <- c("L. minor", "L. japonica", "L. japonica", "L. minuta", "L. minuta", "L. minor", "S. polyrhiza",
                   "L. minor", "L. japonica", "L. japonica", "L. turionifera", "L. japonica", "L. japonica", "L. minuta",
                   "L. japonica", "L. turionifera", "L. minuta", "L. minor", "L. japonica", "L. minor", "L. japonica", "L. minuta",
                   "L. japonica", "L. japonica")



#add grouping for light levels FOR 24 observations INCORRECT core 24
#comb$EnvLight <- c("HL", "LL", "LL", "HL", "LL", "LL", "LL",
#                      "HL", "HL", "HL", "LL", "HL", "LL", "LL",
#                      "HL", "HL", "LL", "LL", "LL", "LL", "HL", "HL",
#                      "LL", "HL")

#for 20 observations
comb$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                   "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                   "HL", "HL", "LL", "LL", "LL", "LL", "HL")
names(comb)
#DIFF SIG DEPENDING ON NO OF OBSERVATIONS INCLUDED
#CURRENT READ OUT FOR 24 - NO ENV INCLUDED
#rgr area
t.test(RGR_diff~EnvLight,comb) #not sig
t.test(RGR_mean.HL~EnvLight,comb) #sig (dll higher) 0.027 *
t.test(RGR_mean.LL~EnvLight,comb) #not sig

#col rgr
t.test(ColRGRnonlog_HL_mean~EnvLight,comb) #sig (ll higher) 0.04 *
t.test(ColRGRnonlog_LL_mean~EnvLight,comb) #not sig
t.test(ColRGR_HL_mean~EnvLight,comb) #sig (ll higher) 0.005 **
t.test(ColRGR_LL_mean~EnvLight,comb) #not sig
t.test(ColRGR_diff~EnvLight,comb) #not sig
t.test(ColRGR_propdiff~EnvLight,comb) #not sig

#fw
t.test(FW.mean.HL~EnvLight,comb) #sig (ll double) * 0.03
t.test(FW.mean.LL~EnvLight,comb) #not sig
t.test(FW_propdiff~EnvLight,comb) #sig

#fdw
t.test(FDW.mean.HL~EnvLight,comb) #sig (ll higher) 0.03 *
t.test(FDW.mean.LL~EnvLight,comb) #not sig
t.test(FDW_propdiff~EnvLight,comb) #sig

#density
t.test(adensity_LL~EnvLight,comb) #not sig
t.test(adensity_HL~EnvLight,comb) #sig (ll higher)
t.test(Cdensity_LL~EnvLight,comb) #not sig
t.test(Cdensity_HL~EnvLight,comb) #not sig

#NPQ
t.test(mean.LL.NPQL1~EnvLight, comb)#not sig
t.test(mean.LL.NPQL3~EnvLight, comb) #not sig
t.test(mean.LL.NPQL5~EnvLight, comb) #not sig
t.test(mean.LL.NPQL11~EnvLight, comb) #not sig #0.36
t.test(mean.HL.NPQL1~EnvLight, comb) #not sig
t.test(mean.HL.NPQL3~EnvLight, comb) #not sig
t.test(mean.HL.NPQL5~EnvLight, comb) #not sig
t.test(mean.HL.NPQL11~EnvLight, comb) #SIG #not sig 24 obs #0.53
t.test(fqfm_L1_mean_HL~EnvLight, comb) #NOT SIG #0.36

#for fig
t.test(mean.LL.NPQL11~EnvLight, comb) #not sig #0.36
t.test(mean.HL.NPQL11~EnvLight, comb) #SIG #not sig 24 obs #0.53
t.test(mean.LL.QYmax~EnvLight, comb) #Not sig 0.40
t.test(mean.HL.QYmax~EnvLight, comb) #NOT SIG #0.29
t.test(fqfm_L1_mean_LL~EnvLight, comb) #Not sig 0.34
t.test(fqfm_L1_mean_HL~EnvLight, comb) #NOT SIG #0.36

#NPQ NOT SIG BY ENVLIGHT SETTING
par(mfrow = c(2, 3))
boxplot(mean.LL.QYmax~EnvLight,data=comb) 
boxplot(mean.HL.QYmax~EnvLight,data=comb) 
boxplot(fqfm_L11_mean_HL~EnvLight,data=comb) 
boxplot(fqfm_L1_mean_LL~EnvLight,data=comb) 
boxplot(mean.LL.NPQL11~EnvLight,data=comb) 
boxplot(mean.HL.NPQL11~EnvLight,data=comb) 

#PHOTOSYNTHESIS BOXPLOTS BY ENV LIGHT
par(mfrow = c(4, 6),  mar=c(2,4.5,2,2))
boxplot(mean.LL.QYmax ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 0.83), xlab = "",
        ylab = "QYmax",
        cex.lab=1.5, cex.axis=1.0)
text(1.0, 0.83, "LL", cex=1)

boxplot(mean.HL.QYmax ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.83), xlab = "",
        ylab = "",
        cex.lab=1.5, cex.axis=1.0
)
text(1.0, 0.83, "HL", cex=1)

boxplot(fqfm_L11_mean_LL ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.6), xlab = "",
        ylab = "Phi PSII",
        cex.lab=1.5, cex.axis=1.0
)
text(4.8, 1.45, "LL", cex=1)

boxplot(fqfm_L11_mean_HL ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.6), xlab = "",
        ylab = "",
        cex.lab=1.5, cex.axis=1.0
)
text(4.8, 1.45, "HL", cex=1)

boxplot(mean.LL.NPQL11 ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 4), xlab = "",
        ylab = "NPQ <1000 PAR",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(mean.HL.NPQL11 ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 4), xlab = "",
        ylab = "",
        cex.lab=1.5, cex.axis=1.0
)


par(mfrow = c(1, 2))
boxplot(mean.HL.NPQL11~EnvLight,data=comb) #sig
boxplot(fqfm_L1_mean_HL~EnvLight,data=comb) #NOT QUITE sig
#boxplot(ColRGR_propdiff~EnvLight,data=comb) not sig

dev.new()

#USING RAW DATA FOR GROWTH PLOTS
#do growth plots with raw data per treatment and look at species
try <- read.csv("GR_Combo_data_spREM_NEWsp_splitbytreat.csv")

#COLOR GROWTH BOXPLOTS BY SPECIES
try$Species <- factor(try$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))
try$EnvLight <- factor(try$EnvLight, levels = c("dLL", "dHL"))

#do statistics species and gr responses, tukey
aov <- aov(RGRlog_LL ~ Species*EnvLight, data = try)  #sp ***, env ns 0.6
aov <- aov(RGRlog_HL ~ Species*EnvLight, data = try)  #sp *, env *, 0.01, env*sp *
aov <- aov(Col_RGR_LL ~ Species*EnvLight, data = try) #sp *, ns env 0.7
aov <- aov(Col_RGR_HL ~ Species*EnvLight, data = try) #sp *, ennv ** 0.01
aov <- aov(FW_mod_LL ~ Species*EnvLight, data = try)  #sp **, env 0.9 ns
aov <- aov(FW_mod_HL ~ Species*EnvLight, data = try) #ns env 0.3
aov <- aov(FDW..mg.LL ~ Species*EnvLight, data = try) #ns env 0.8
aov <- aov(FDW..mg._HL ~ Species*EnvLight, data = try) #ns env 0.5
aov <- aov(FW_mod_LL ~ EnvLight, data = try)  #ns
aov <- aov(FW_mod_HL ~ EnvLight, data = try) #ns
aov <- aov(FDW..mg.LL ~ EnvLight, data = try) #ns
aov <- aov(FDW..mg._HL ~ EnvLight, data = try) #ns

summary(aov)

#tukey tests
tuk_out <- TukeyHSD(aov, "EnvLight", conf.level=.95)
tuk_out
#tukey tests
tuk_out <- TukeyHSD(aov, "Species", conf.level=.95)
tuk_out
#RGRlog LL - l turi and l mino, l turi and l jap, l turi and l minu
#RGRlog HL - l minu and l jap
#Col LL - s poly and l jap
#Col HL - ns
#FW LL - l turi and l mino, l turi and l jap
#others ns
summary(aov)

species <-list("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza")

tiff('Growthrate_boxplots_Species_sig+envlight_raw.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('Growthrate_boxplots_Species_sig+envlight_raw.pdf', width=14, height=12)
par(mfrow = c(4, 4),  mar=c(2,4.5,2,2))
boxplot(RGRlog_LL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 1.5), xlab = "", xaxt="n",
        ylab = "RGRlog area",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 1.45, "LL", cex=2.75)
text(1, 1.3, "a", cex=2)
text(2, 1.3, "a", cex=2)
text(3, 1.3, "b", cex=2)
text(4, 1.3, "a", cex=2)
text(5, 1.3, "ab", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("A.", cex.main=2, adj=0)
stripchart(RGRlog_LL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(RGRlog_HL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1.5), xlab = "", xaxt="n",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 1.45, "HL", cex=2.75)
text(1, 1.25, "ab", cex=2)
text(2, 1.25, "a", cex=2)
text(3, 1.25, "ab", cex=2)
text(4, 1.25, "b", cex=2)
text(5, 1.25, "ab", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(RGRlog_HL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Col_RGR_LL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 1.2), xlab = "", xaxt="n",
        ylab = "Col RGRlog",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 1.15, "LL", cex=2.75)
text(1, 1, "ab", cex=2)
text(2, 1, "a", cex=2)
text(3, 1, "ab", cex=2)
text(4, 1, "ab", cex=2)
text(5, 1, "b", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("B.", cex.main=2, adj=0)
stripchart(Col_RGR_LL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Col_RGR_HL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1.2), xlab = "", xaxt="n",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 1.15, "HL", cex=2.75)
text(1, 1, "a", cex=2)
text(2, 1, "a", cex=2)
text(3, 1, "a", cex=2)
text(4, 1, "a", cex=2)
text(5, 1, "a", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(Col_RGR_HL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW_mod_LL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0, 
        ylim=c(0, 6000), xlab = "",
        ylab = "FW",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 5900, "LL", cex=2.75)
text(1, 5200, "a", cex=2)
text(2, 5200, "a", cex=2)
text(3, 5200, "b", cex=2)
text(4, 5200, "ab", cex=2)
text(5, 5200, "ab", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(FW_mod_LL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW_mod_HL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 6000), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 5900, "HL", cex=2.75)
text(1, 5200, "a", cex=2)
text(2, 5200, "a", cex=2)
text(3, 5200, "a", cex=2)
text(4, 5200, "a", cex=2)
text(5, 5200, "a", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(FW_mod_HL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW..mg.LL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0, 
        ylim=c(0, 350), xlab = "",
        ylab = "FDW",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 340, "LL", cex=2.75)
text(1, 280, "a", cex=2)
text(2, 280, "a", cex=2)
text(3, 280, "a", cex=2)
text(4, 280, "a", cex=2)
text(5, 280, "a", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("D.", cex.main=2, adj=0)
stripchart(FDW..mg.LL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW..mg._HL ~ Species, try,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 350), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 340, "HL", cex=2.75)
text(1, 280, "a", cex=2)
text(2, 280, "a", cex=2)
text(3, 280, "a", cex=2)
text(4, 280, "a", cex=2)
text(5, 280, "a", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(FDW..mg._HL ~ Species, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(RGRlog_LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 1.5), xlab = "",
        ylab = "RGRlog area",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 1.45, "LL", cex=2.75)
text(2.1, 1.45, "p = 0.57", cex=2)
title("A.", cex.main=2, adj=0)
stripchart(RGRlog_LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(RGRlog_HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1.5), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
text(1, 1.2, "*", cex=3.5)
text(0.9, 1.45, "HL", cex=2.75)
text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(RGRlog_HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Col_RGR_LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1), xlab = "",
        ylab = "Col RGRlog",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 0.96, "LL", cex=2.75)
text(2.1, 0.96, "p = 0.7", cex=2)
title("B.", cex.main=2, adj=0)
stripchart(Col_RGR_LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Col_RGR_HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
text(1, 0.8, "**", cex=3.5)
text(0.8, 0.96, "HL", cex=2.75)
text(2.1, 0.96, "p = 0.01", cex=2)
stripchart(Col_RGR_HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW_mod_LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 5500), xlab = "",
        ylab = "FW (mg)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 5430, "LL", cex=2.75)
text(2.1, 5430, "p = 0.9", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(FW_mod_LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW_mod_HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 5500), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1800, "*", cex=3.5)
text(0.8, 5430, "HL", cex=2.75)
text(2.1, 5430, "p = 0.3", cex=2)
stripchart(FW_mod_HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW..mg.LL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 200), xlab = "",
        ylab = "FDW (mg)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 190, "LL", cex=2.75)
text(2.1, 190, "p = 0.75", cex=2)
title("D.", cex.main=2, adj=0)
stripchart(FDW..mg.LL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW..mg._HL ~ EnvLight, try,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 200), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 120, "*", cex=3.5)
text(0.8, 190, "HL", cex=2.75)
text(2.1, 190, "p = 0.5", cex=2)
stripchart(FDW..mg._HL ~ EnvLight, try,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

dev.off()
dev.off()

#WORKING WITH COMB DATA, SPECIES AFFECTS?
#do statistics species and gr responses, tukey
aov <- aov(RGR_mean.LL ~ Species*EnvLight, data = comb)  #sp ***, env ns
aov <- aov(RGR_mean.HL ~ Species*EnvLight, data = comb)  #sp **, env **, env*sp *
aov <- aov(ColRGR_LL_mean ~ Species*EnvLight, data = comb) #sp *, ns env
aov <- aov(ColRGR_HL_mean ~ Species*EnvLight, data = comb) #sp ns, env *
aov <- aov(FW.mean.LL ~ Species*EnvLight, data = comb)  #sp **
aov <- aov(FW.mean.HL ~ Species*EnvLight, data = comb) #ns
aov <- aov(FDW.mean.LL ~ Species*EnvLight, data = comb) #ns
aov <- aov(FDW.mean.HL ~ Species*EnvLight, data = comb) #ns
aov <- aov(FW.mean.LL ~ EnvLight, data = comb)  #sp **
aov <- aov(FW.mean.HL ~ EnvLight, data = comb) #*
aov <- aov(FDW.mean.LL ~ EnvLight, data = comb) #ns
aov <- aov(FDW.mean.HL ~ EnvLight, data = comb) #*

summary(aov)

#tukey tests
tuk_out <- TukeyHSD(aov, "Species", conf.level=.95)
tuk_out
#RGR LL - L jap and L. minu, L. jap and L. turi, L mino L minu,
#L turio and L mino, L turio and L minu
#RGR HL - L minu and L jap
#Col LL - S poly L jap
#Col HL - ns
#FW LL - L turio and L jap, L minu L mino, L turio L mino,
#FW HL - ns
#FDW LL - ns
#FDW HL - ns

#for paper colored with stars
#GROWTH
comb$Species <- factor(comb$Species, levels = c("L. minor", "L. japonica", "L. turionifera", "L. minuta", "S. polyrhiza"))
comb$EnvLight <- factor(comb$EnvLight, levels=c("dLL", "dHL"))
#need to configure distances so less gaps between plots
tiff('Growthrate_boxplots_Envlight+species.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(2, 4),  mar=c(2,4.5,2,2))
boxplot(RGR_mean.LL ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 1.5), xlab = "", xaxt="n",
        ylab = "RGRlog area",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 1.45, "LL", cex=2.75)
text(1, 1.3, "b", cex=2)
text(2, 1.3, "a", cex=2)
text(3, 1.3, "c", cex=2)
text(4, 1.3, "c", cex=2)
text(5, 1.3, "abc", cex=2)
#RGR LL - L jap and L. minu, L. jap and L. turi, L mino L minu,
#L turio and L mino, L turio and L minu
#text(2.1, 1.45, "p = 0.57", cex=2)
title("A.", cex.main=2, adj=0)
stripchart(RGR_mean.LL ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(RGR_mean.HL ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1.5), xlab = "", xaxt="n",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 1.45, "HL", cex=2.75)
text(1, 1.25, "ab", cex=2)
text(2, 1.25, "a", cex=2)
text(3, 1.25, "ab", cex=2)
text(4, 1.25, "b", cex=2)
text(5, 1.25, "ab", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(RGR_mean.HL ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(ColRGR_LL_mean ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 1.2), xlab = "", xaxt="n",
        ylab = "Col RGRlog",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 1.15, "LL", cex=2.75)
text(1, 1, "ab", cex=2)
text(2, 1, "a", cex=2)
text(3, 1, "ab", cex=2)
text(4, 1, "ab", cex=2)
text(5, 1, "b", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("B.", cex.main=2, adj=0)
stripchart(ColRGR_LL_mean ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(ColRGR_HL_mean ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1.2), xlab = "", xaxt="n",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 1.15, "HL", cex=2.75)
text(1, 1, "a", cex=2)
text(2, 1, "a", cex=2)
text(3, 1, "a", cex=2)
text(4, 1, "a", cex=2)
text(5, 1, "a", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(ColRGR_HL_mean ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW.mean.LL ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0, 
        ylim=c(0, 6000), xlab = "",
        ylab = "FW",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 5900, "LL", cex=2.75)
text(1, 5200, "a", cex=2)
text(2, 5200, "a", cex=2)
text(3, 5200, "b", cex=2)
text(4, 5200, "ab", cex=2)
text(5, 5200, "ab", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(FW.mean.LL ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW.mean.HL ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 6000), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 5900, "HL", cex=2.75)
text(1, 5200, "a", cex=2)
text(2, 5200, "a", cex=2)
text(3, 5200, "a", cex=2)
text(4, 5200, "a", cex=2)
text(5, 5200, "a", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(FW.mean.HL ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW.mean.LL ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0, 
        ylim=c(0, 350), xlab = "",
        ylab = "FDW",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 340, "LL", cex=2.75)
text(1, 280, "a", cex=2)
text(2, 280, "a", cex=2)
text(3, 280, "a", cex=2)
text(4, 280, "a", cex=2)
text(5, 280, "a", cex=2)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("D.", cex.main=2, adj=0)
stripchart(FDW.mean.LL ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW.mean.HL ~ Species, comb,
        col=c("#E69F00", "#242424", "#CC79A7", "#009E73", "#D55E00"),
        las=2, whisklty=1, staplewex=0,
        ylim=c(0, 350), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 1.05, "*", cex=3.5)
text(0.9, 340, "HL", cex=2.75)
text(1, 280, "a", cex=2)
text(2, 280, "a", cex=2)
text(3, 280, "a", cex=2)
text(4, 280, "a", cex=2)
text(5, 280, "a", cex=2)
#text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(FDW.mean.HL ~ Species, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(RGR_mean.LL ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
     las=1, whisklty=1, staplewex=0, 
     ylim=c(0, 1.5), xlab = "",
     ylab = "RGRlog area",
     cex.lab=2.0, cex.axis=1.5
)
text(0.8, 1.45, "LL", cex=2.75)
text(2.1, 1.45, "p = 0.57", cex=2)
title("A.", cex.main=2, adj=0)
stripchart(RGR_mean.LL ~ EnvLight, comb,             # Data
            method = "jitter", # Random noise
            pch = 19,          # Pch symbols
            cex = 0.8,
            col = 1,           # Color of the symbol
            vertical = TRUE,   # Vertical mode
            add = TRUE)

boxplot(RGR_mean.HL ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 1.5), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
text(1, 1.05, "*", cex=3.5)
text(0.9, 1.45, "HL", cex=2.75)
text(2.1, 1.45, "p = 0.02", cex=2)
stripchart(RGR_mean.HL ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(ColRGR_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.6), xlab = "",
        ylab = "Col RGRlog",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 0.6, "LL", cex=2.75)
text(2.1, 0.6, "p = 0.19", cex=2)
title("B.", cex.main=2, adj=0)
stripchart(ColRGR_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(ColRGR_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.6), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
text(1, 0.56, "**", cex=3.5)
text(0.8, 0.6, "HL", cex=2.75)
text(2.1, 0.6, "p = 0.0025", cex=2)
stripchart(ColRGR_HL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW.mean.LL ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2500), xlab = "",
        ylab = "FW (mg)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 2500, "LL", cex=2.75)
text(2.1, 2500, "p = 0.75", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(FW.mean.LL ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FW.mean.HL ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2500), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
text(1, 1800, "*", cex=3.5)
text(0.8, 2500, "HL", cex=2.75)
text(2.1, 2500, "p = 0.03", cex=2)
stripchart(FW.mean.HL ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW.mean.LL ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 200), xlab = "",
        ylab = "FDW (mg)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 200, "LL", cex=2.75)
text(2.1, 200, "p = 0.75", cex=2)
title("D.", cex.main=2, adj=0)
stripchart(FDW.mean.LL ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(FDW.mean.HL ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 200), xlab = "",
        ylab = NULL,
        cex.lab=2.0, cex.axis=1.5
)
text(1, 120, "*", cex=3.5)
text(0.8, 200, "HL", cex=2.75)
text(2.1, 200, "p = 0.03", cex=2)
stripchart(FDW.mean.HL ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

dev.off()

#basic plots
boxplot(RGR_mean.LL~EnvLight,comb, ylim = c(0, 1.5)) #not sig
boxplot(RGR_mean.HL~EnvLight,comb, ylim = c(0, 1.5)) #sig
boxplot(ColRGR_LL_mean~EnvLight,comb, ylim = c(0, 0.6)) #not sig
boxplot(ColRGR_HL_mean~EnvLight,comb, ylim = c(0, 0.6)) #sig
boxplot(FW.mean.LL~EnvLight,comb, ylim = c(0, 2500)) #not sig
boxplot(FW.mean.HL~EnvLight,comb, ylim = c(0, 2500)) #sig
boxplot(FDW.mean.LL~EnvLight,comb, ylim = c(0, 200)) #not sig
boxplot(FDW.mean.HL~EnvLight,comb, ylim = c(0, 200)) #sig
#boxplot(adensity_LL~EnvLight,comb) #not sig
#boxplot(adensity_HL~EnvLight,comb) #sig

#PHOTOSYNTHESIS
par(mfrow = c(1, 2))
boxplot(fqfm_L1_mean_LL~EnvLight,data=comb)
boxplot(mean.LL.NPQL11~EnvLight,data=comb)

t.test(mean.LL.NPQL11~EnvLight, comb) #not sig
t.test(fqfm_L1_mean_LL~EnvLight, comb) #not sig

#aov <- aov(mean.LL.NPQL1~EnvLight,data=comb) #not sig
t.test(fqfm_L3_mean_LL~EnvLight, comb) #NOT SIG
t.test(fqfm_L1_mean_LL~EnvLight, comb) #NOT SIG
t.test(fqfm_L11_mean_LL~EnvLight, comb) #NOT SIG
t.test(fqfm_L3_mean_HL~EnvLight, comb) #NOT SIG #sig 24 obs
t.test(fqfm_L1_mean_HL~EnvLight, comb) #NOT SIG
t.test(fqfm_L11_mean_HL~EnvLight, comb) #NOT SIG #sig 24 obs

#fqfm not sig between env light at any light level, treat

comb$Accession1 <- NULL
comb$Treatment.y <- NULL
comb$X <- NULL

names(comb)

#for paper
par(mfrow = c(1, 2)
myggp <- ggplot(comb, aes(PRI_LL_mean, Car_rat_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point(pch=3) +
  geom_smooth(method='lm') +
  xlab("PRI value LL") +
  ylab("Carotenoid ratio LL") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.29)
myggp + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))

#for paper
myggp <- ggplot(comb, aes(PRI_HL_mean, Car_rat_HL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point(pch=3) +
  geom_smooth(method='lm') +
  xlab("PRI value HL") +
  ylab("Carotenoid ratio HL") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.01, label.y = 0.33)
myggp + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))

myggp <- ggplot(comb, aes(PRI_LL_mean, Carot_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point(pch=3) +
  geom_smooth(method='lm') +
  xlab("PRI value LL") +
  ylab("Carotenoid content LL") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.25)
myggp + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))


myggp <- ggplot(comb, aes(PRI_HL_mean, Carot_HL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point(pch=3) +
  geom_smooth(method='lm') +
  xlab("PRI value HL") +
  ylab("Carotenoid content HL") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.25)
myggp + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))


#automate anovas
#moved species to start
flall <- comb
str(flall)

#manually moved species to start as automation fails if in middle col

library(dplyr)

#run this to save to correct place
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Comb")


#ANOVA AUTOMATION FROM MATT K
#starts from col no, #accession and treat in pos 
my_formula <- as.formula(paste(colnames(flall)[190], "~", paste(colnames(flall)[c(191)], collapse = "*")))

anova_output2 <- list()
for(i in 2:190){
  my_col <- colnames(flall[i]) 
  my_formula <- as.formula(paste(my_col, "~", paste(colnames(flall)[c(191)], collapse = "*")))
  sum_aov <- summary(aov(my_formula, flall))
  sum_tab <- sum_aov[[1]] %>% as.data.frame()
  sum_tab$ind <- my_col
  anova_output2[[i]] <- sum_tab
  print(my_formula)
}

(anova_output2)


all_output2 <- do.call(rbind, anova_output2)

all_output2

write.csv(all_output2, "combo_all_anovaoutput.csv")

#takes longer when add extra variables, from 2 up to 4
#halves the residual number by adding rep and within rep but also all sig
model=lm(flall$RGR_mean.HL ~ flall$EnvLight)
#model=lm(flall$NPQ_L11 ~ flall$Accession * flall$Treatment)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'flall$EnvLight', conf.level=0.95)
summary(TUKEY)
plot(TUKEY , las=1 , col="brown")

#from output 24 OBS (not chl/ hs?)
t.test(maximum.LL.FvFm3~EnvLight, comb) #SIG #not sig
t.test(mean.HL.NPQL11~EnvLight, comb) #SIG #not sig
t.test(fqfm_L3_max_LL~EnvLight, comb) #SIG #not sig
t.test(ColRGR_LL_maximum~EnvLight, comb) #SIG #not sig

#from output 28 OBS
t.test(maximum.LL.FvFm3~EnvLight, comb) #SIG
t.test(mean.HL.NPQL11~EnvLight, comb) #NOT SIG
t.test(fqfm_L3_max_LL~EnvLight, comb) #SIG
t.test(ColRGR_LL_maximum~EnvLight, comb) #NOT SIG #sig 24 obs

#from output 20 OBS ALL NOT SIG
t.test(maximum.LL.FvFm3~EnvLight, comb) #NOT SIG
t.test(mean.HL.NPQL11~EnvLight, comb) #NOT SIG
t.test(fqfm_L3_max_LL~EnvLight, comb) #NOT SIG
t.test(ColRGR_LL_maximum~EnvLight, comb) #NOT SIG #sig 24 obs

par(mfrow = c(2, 2))
boxplot(maximum.LL.FvFm3~EnvLight, comb) #fvfm higher in ll 100 umol in accessions from ll
boxplot(mean.HL.NPQL11~EnvLight, comb) #npq higher in ll sites, decline in hl sites
boxplot(fqfm_L3_max_LL~EnvLight, comb) #same as fvfm3 ll
boxplot(ColRGR_LL_maximum~EnvLight, comb) #col gains in ll higher than hl sites
#more that ll sites better equiped ll - faster growth, higher max fqfm
#in hl treated npq higher in ll derived plants

#does it explain anything seen in nature?

plot(maximum.LL.FvFm3~Summer2022_light, comb)
plot(maximum.LL.FvFm3~ColRGR_LL_maximum, comb)
plot(mean.HL.NPQL11~ColRGR_LL_maximum, comb) #most sig

#ranking of npq and photosynthetic efficiency accessions
names(comb)
#hl response
NPQL11_HL_ord <- comb %>% arrange(desc(mean.HL.NPQL11)) %>% select(Accession, mean.HL.NPQL11) %>% top_n(24)
NPQL3_HL_ord <- comb %>% arrange(desc(mean.HL.NPQL3)) %>% select(Accession, mean.HL.NPQL3) %>% top_n(24)
NPQL5_HL_ord <- comb %>% arrange(desc(mean.HL.NPQL5)) %>% select(Accession, mean.HL.NPQL5) %>% top_n(24)
NPQL1_HL_ord <- comb %>% arrange(desc(mean.HL.NPQL1)) %>% select(Accession, mean.HL.NPQL1) %>% top_n(24)
fqfm_L11_HL_ord <- comb %>% arrange(desc(fqfm_L11_mean_HL)) %>% select(Accession, fqfm_L11_mean_HL) %>% top_n(24)
fqfm_L3_HL_ord <- comb %>% arrange(desc(fqfm_L3_mean_HL)) %>% select(Accession, fqfm_L3_mean_HL) %>% top_n(24)
fqfm_L5_HL_ord <- comb %>% arrange(desc(fqfm_L5_mean_HL)) %>% select(Accession, fqfm_L5_mean_HL) %>% top_n(24)
fqfm_L1_HL_ord <- comb %>% arrange(desc(fqfm_L1_mean_HL)) %>% select(Accession, fqfm_L1_mean_HL) %>% top_n(24)
QYmax_HL_ord <- comb %>% arrange(desc(mean.HL.QYmax)) %>% select(Accession, mean.HL.QYmax) %>% top_n(24)

#ll response
NPQL11_LL_ord <- comb %>% arrange(desc(mean.LL.NPQL11)) %>% select(Accession, mean.LL.NPQL11) %>% top_n(24)
NPQL3_LL_ord <- comb %>% arrange(desc(mean.LL.NPQL3)) %>% select(Accession, mean.LL.NPQL3) %>% top_n(24)
NPQL5_LL_ord <- comb %>% arrange(desc(mean.LL.NPQL5)) %>% select(Accession, mean.LL.NPQL5) %>% top_n(24)
NPQL1_LL_ord <- comb %>% arrange(desc(mean.LL.NPQL1)) %>% select(Accession, mean.LL.NPQL1) %>% top_n(24)
fqfm_L11_LL_ord <- comb %>% arrange(desc(fqfm_L11_mean_LL)) %>% select(Accession, fqfm_L11_mean_LL) %>% top_n(24)
fqfm_L3_LL_ord <- comb %>% arrange(desc(fqfm_L3_mean_LL)) %>% select(Accession, fqfm_L3_mean_LL) %>% top_n(24)
fqfm_L5_LL_ord <- comb %>% arrange(desc(fqfm_L5_mean_LL)) %>% select(Accession, fqfm_L5_mean_LL) %>% top_n(24)
fqfm_L1_LL_ord <- comb %>% arrange(desc(fqfm_L1_mean_LL)) %>% select(Accession, fqfm_L1_mean_LL) %>% top_n(24)
QYmax_LL_ord <- comb %>% arrange(desc(mean.LL.QYmax)) %>% select(Accession, mean.LL.QYmax) %>% top_n(24)

#DIFFERENCES (MUTATED BELOW)
QYmax_diff_ord <- comb %>% arrange(desc(meanQYmax_diff)) %>% select(Accession, meanQYmax_diff) %>% top_n(24)
NPQ_L11_diff_ord <- comb %>% arrange(desc(meanNPQ_L11_diff)) %>% select(Accession, meanNPQ_L11_diff) %>% top_n(24)
NPQ_L5_diff_ord <- comb %>% arrange(desc(meanNPQ_L5_diff)) %>% select(Accession, meanNPQ_L5_diff) %>% top_n(24)
NPQ_L3_diff_ord <- comb %>% arrange(desc(meanNPQ_L3_diff)) %>% select(Accession, meanNPQ_L3_diff) %>% top_n(24)
NPQ_L1_diff_ord <- comb %>% arrange(desc(meanNPQ_L1_diff)) %>% select(Accession, meanNPQ_L1_diff) %>% top_n(24)
fqfm_L3_diff_ord <- comb %>% arrange(desc(meanfqfm_L3_diff)) %>% select(Accession, meanfqfm_L3_diff) %>% top_n(24)
fqfm_L5_diff_ord <- comb %>% arrange(desc(meanfqfm_L5_diff)) %>% select(Accession, meanfqfm_L5_diff) %>% top_n(24)
fqfm_L11_diff_ord <- comb %>% arrange(desc(meanfqfm_L11_diff)) %>% select(Accession, meanfqfm_L11_diff) %>% top_n(24)

fluo_rankings_hl <-cbind(NPQL11_HL_ord,NPQL3_HL_ord,NPQL5_HL_ord,NPQL1_HL_ord,
                         fqfm_L11_HL_ord,fqfm_L3_HL_ord,fqfm_L5_HL_ord,fqfm_L1_HL_ord,QYmax_HL_ord)

fluo_rankings_ll <-cbind(NPQL11_LL_ord,NPQL3_LL_ord,NPQL5_LL_ord,NPQL1_LL_ord,
                         fqfm_L11_LL_ord,fqfm_L3_LL_ord,fqfm_L5_LL_ord,fqfm_L1_LL_ord,QYmax_LL_ord)

fluo_rankings_diff <-cbind(QYmax_diff_ord,NPQ_L11_diff_ord,NPQ_L5_diff_ord,NPQ_L3_diff_ord,
                           NPQ_L1_diff_ord,fqfm_L3_diff_ord,fqfm_L5_diff_ord,fqfm_L11_diff_ord)

fluo_rankings <- cbind(fluo_rankings_hl, fluo_rankings_ll, fluo_rankings_diff)

class(fluo_rankings)       

write.csv(fluo_rankings, "Fluo_HL_LL_diff_rankings.csv")

#chlorophyll rankings hl
Chla_HL_ord <- comb %>% arrange(desc(Chl_a_HL_mean)) %>% select(Accession, Chl_a_HL_mean) %>% top_n(24)
Chlb_HL_ord <- comb %>% arrange(desc(Chl_b_HL_mean)) %>% select(Accession, Chl_b_HL_mean) %>% top_n(24)
Chlab_HL_ord <- comb %>% arrange(desc(Carot_HL_mean)) %>% select(Accession, Carot_HL_mean) %>% top_n(24)
Carot_HL_ord <- comb %>% arrange(desc(Chl_ab_HL_mean)) %>% select(Accession, Chl_ab_HL_mean) %>% top_n(24)
Chlabrat_HL_ord <- comb %>% arrange(desc(Chlab_rat_HL_mean)) %>% select(Accession, Chlab_rat_HL_mean) %>% top_n(24)
Car_rat_HL_ord <- comb %>% arrange(desc(Car_rat_HL_mean)) %>% select(Accession, Car_rat_HL_mean) %>% top_n(24)

#chlorophyll rankings ll
Chla_LL_ord <- comb %>% arrange(desc(Chl_a_LL_mean)) %>% select(Accession, Chl_a_LL_mean) %>% top_n(24)
Chlb_LL_ord <- comb %>% arrange(desc(Chl_b_LL_mean)) %>% select(Accession, Chl_b_LL_mean) %>% top_n(24)
Chlab_LL_ord <- comb %>% arrange(desc(Carot_LL_mean)) %>% select(Accession, Carot_LL_mean) %>% top_n(24)
Carot_LL_ord <- comb %>% arrange(desc(Chl_ab_LL_mean)) %>% select(Accession, Chl_ab_LL_mean) %>% top_n(24)
Chlabrat_LL_ord <- comb %>% arrange(desc(Chlab_rat_LL_mean)) %>% select(Accession, Chlab_rat_LL_mean) %>% top_n(24)
Car_rat_LL_ord <- comb %>% arrange(desc(Car_rat_LL_mean)) %>% select(Accession, Car_rat_LL_mean) %>% top_n(24)

chl_rankings_hl <- cbind(Chla_HL_ord,Chlb_HL_ord,Chlab_HL_ord,
                         Carot_HL_ord,Chlabrat_HL_ord,Car_rat_HL_ord)

chl_rankings_ll <- cbind(Chla_LL_ord,Chlb_LL_ord,Chlab_LL_ord,
                         Carot_LL_ord,Chlabrat_LL_ord,Car_rat_LL_ord)

chl_rankings <- cbind(chl_rankings_hl, chl_rankings_ll)

class(chl_rankings)       

write.csv(chl_rankings, "Chl_HL_LL_rankings.csv")


                
#GR CORRELATIONS
corr <- cor.test(comb$RGR_mean.HL, comb$ColRGR_HL_mean,
                 method = "pearson"
)
corr$estimate #weak pos
corr$p.value # sig

corr <- cor.test(comb$RGR_mean.HL, RGR_mean.LL,
                 method = "pearson"
)
corr$estimate #weak pos
corr$p.value # sig

corr <- cor.test(comb$RGR_mean.LL, comb$ColRGR_LL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$RGR_mean.HL, comb$FW.mean.HL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$RGR_mean.LL, comb$FW.mean.LL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$FW.mean.HL, comb$ColRGR_HL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$FW.mean.LL, comb$ColRGR_LL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig


corr <- cor.test(comb$FW.mean.HL, comb$FDW.mean.HL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$FW.mean.LL, comb$FDW.mean.LL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig


corr <- cor.test(comb$FDW.mean.HL, comb$ColRGR_HL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$FDW.mean.LL, comb$ColRGR_LL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$FDW.mean.HL, comb$RGR_mean.HL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig

corr <- cor.test(comb$FDW.mean.LL, comb$RGR_mean.LL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value # sig


#OLD
corr <- cor.test(comb$maximum.LL.FvFm3, comb$Summer2022_light,
                 method = "pearson"
)
corr$estimate #weak neg
corr$p.value #sig

corr <- cor.test(comb$maximum.LL.FvFm3, comb$ColRGR_LL_maximum,
                 method = "pearson"
)
corr$estimate #no relationship
corr$p.value #not sig

#adapted for ll
corr <- cor.test(comb$mean.HL.NPQL11, comb$ColRGR_LL_maximum,
                 method = "pearson"
)
corr$estimate #weak positive
corr$p.value #sig

plot(mean.HL.NPQL11~Jul.22, comb)

corr <- cor.test(comb$mean.HL.NPQL11, comb$Jul.22,
                 method = "pearson"
)
corr$estimate #weak neg
corr$p.value

#none of coverage values linked to gr data
plot(ColRGR_LL_maximum~Jul.22, comb)
plot(ColRGR_HL_maximum~Jul.22, comb)

#go together
plot(mean.HL.NPQL11~Summer2022_light, comb)
plot(ColRGR_LL_maximum~Summer2022_light, comb)

#not quite sig #not sig 24 core access #not sig
t.test(maximum.LL.FvFm5~EnvLight, comb)
t.test(mean.HL.QYmax~EnvLight, comb)
t.test(fqfm_L1_mean_HL~EnvLight, comb)
t.test(fqfm_L5_max_LL~EnvLight, comb)
t.test(maximum.HL.NPQL5~EnvLight, comb)
t.test(RGRnonlog_LL_maximum~EnvLight, comb)

#look at rankings for ColRGR and RGRlog per treatment
#compare ranking data for col gain and rgr area
Col_RGR_HL_ord <- comb %>% arrange(desc(ColRGR_HL_mean)) %>% select(Accession, ColRGR_HL_mean)
Col_RGR_LL_ord <- comb %>% arrange(desc(ColRGR_LL_mean)) %>% select(Accession, ColRGR_LL_mean)
RGR_mean.HL_ord <- comb %>% arrange(desc(RGR_mean.HL)) %>% select(Accession, RGR_mean.HL)
RGR_mean.LL_ord <- comb %>% arrange(desc(RGR_mean.LL)) %>% select(Accession, RGR_mean.LL)
rankings <- cbind(Col_RGR_HL_ord, Col_RGR_LL_ord, RGR_mean.HL_ord, RGR_mean.LL_ord)
#top 4 accessions ly03, sel1, nuff1, moor1

#redo with more gr param and no extra accessions to see rankings/consistancy
names(comb)
Col_RGR_HL_ord <- comb %>% arrange(desc(ColRGR_HL_mean)) %>% select(Accession, ColRGR_HL_mean)
Col_RGR_LL_ord <- comb %>% arrange(desc(ColRGR_LL_mean)) %>% select(Accession, ColRGR_LL_mean)
RGR_mean.HL_ord <- comb %>% arrange(desc(RGR_mean.HL)) %>% select(Accession, RGR_mean.HL)
RGR_mean.LL_ord <- comb %>% arrange(desc(RGR_mean.LL)) %>% select(Accession, RGR_mean.LL)
FW_mean.HL_ord <- comb %>% arrange(desc(FW.mean.HL)) %>% select(Accession, FW.mean.HL)
FW_mean.LL_ord <- comb %>% arrange(desc(FW.mean.LL)) %>% select(Accession, FW.mean.LL)
FDW_mean.HL_ord <- comb %>% arrange(desc(FDW.mean.HL)) %>% select(Accession, FDW.mean.HL)
FDW_mean.LL_ord <- comb %>% arrange(desc(FDW.mean.LL)) %>% select(Accession, FDW.mean.LL)
rankings <- cbind(Col_RGR_HL_ord, Col_RGR_LL_ord, RGR_mean.HL_ord, RGR_mean.LL_ord, FW_mean.HL_ord, FW_mean.LL_ord,
                  FDW_mean.HL_ord, FDW_mean.LL_ord)

#OLD
Col_RGR_HL_ord <- new %>% arrange(desc(ColRGR_HL_mean)) %>% select(Accession, ColRGR_HL_mean)
Col_RGR_LL_ord <- new %>% arrange(desc(ColRGR_LL_mean)) %>% select(Accession, ColRGR_LL_mean)
RGR_mean.HL_ord <- new %>% arrange(desc(RGR_mean.HL)) %>% select(Accession, RGR_mean.HL)
RGR_mean.LL_ord <- new %>% arrange(desc(RGR_mean.LL)) %>% select(Accession, RGR_mean.LL)
rankings <- cbind(Col_RGR_HL_ord, Col_RGR_LL_ord, RGR_mean.HL_ord, RGR_mean.LL_ord)
#top 4 accessions ly03, sel1, nuff1, moor1

names(new)

#see order for diff
RGR_diff_ord <- new %>% arrange(desc(RGR_diff)) %>% select(Accession, RGR_diff)
RGRpropchange4_ord <- new %>% arrange(desc(RGRpropchange4)) %>% select(Accession, RGRpropchange4)
ColRGR_diff_ord <- new %>% arrange(desc(ColRGR_diff)) %>% select(Accession, ColRGR_diff)
rankings <- cbind(RGR_diff_ord, RGRpropchange4_ord, ColRGR_diff_ord)

names(comb)
comb$Accession

#T TESTS
comb$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                     "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                     "HL", "HL", "LL", "LL", "LL", "LL", "HL", "HL",
                     "HL", "LL", "LL")

#all not sig obs 24
t.test(PRI_HL_mean~EnvLight, comb) #not sig 0.46

t.test(Chl_a_LL_mean~EnvLight, comb) #not SIG 0.98
t.test(Chl_a_HL_mean~EnvLight, comb) #SIG * #0.04
t.test(Chl_b_LL_mean~EnvLight, comb) #not SIG 0.71
t.test(Chl_b_HL_mean~EnvLight, comb) #SIG * #0.04

t.test(Chl_ab_LL_mean~EnvLight, comb) #not sig 0.94
t.test(Chl_ab_HL_mean~EnvLight, comb) #SIG 24 obs * 0.04
t.test(Chlab_rat_HL_mean~EnvLight, comb) #NOT SIG 0.21
t.test(Chlab_rat_LL_mean~EnvLight, comb) #NOT SIG 0.18

t.test(Carot_LL_mean~EnvLight, comb) #not SIG =0.87
t.test(Carot_HL_mean~EnvLight, comb) #SIG * #0.03

t.test(Car_rat_HL_mean~EnvLight, comb) #not sig 0.63
t.test(Car_rat_LL_mean~EnvLight, comb) #not sig 0.57

#for paper
#old paper version colored with stars
#CHLOROPHYLL AND CARPTENOIDS
#need to configure distances so less gaps between plots
par(mfrow = c(3, 4),  mar=c(2,4.5,2,2))
boxplot(Chl_a_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 10), xlab = "",
        ylab = "Chl a (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)
text(0.8, 9.5, "LL", cex=1.75)

boxplot(Chl_a_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 10), xlab = "",
        ylab = "Chl a (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 8, "*", cex=2.5)
text(0.8, 9.5, "HL", cex=1.75)
boxplot(Chl_b_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "Chl b (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Chl_b_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "Chl b (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 2, "*", cex=2.5)
boxplot(Carot_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "Total Carotenoids (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Carot_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "Total Carotenoids (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 2.3, "*", cex=2.5)

boxplot(Chl_ab_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 12), xlab = "",
        ylab = "Chl a+b (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Chl_ab_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 12), xlab = "",
        ylab = "Chl a+b (mg/g)",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 10, "*", cex=2.5)

boxplot(Chlab_rat_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(3, 4.5), xlab = "",
        ylab = "Chl a:b",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Chlab_rat_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(3, 4.5), xlab = "",
        ylab = "Chl a:b",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Car_rat_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.5), xlab = "",
        ylab = "Car:Chl",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Car_rat_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.5), xlab = "",
        ylab = "Car:Chl",
        cex.lab=1.5, cex.axis=1.0
)

#format for may 2023 version
comb$EnvLight <- factor(comb$EnvLight, levels=c("dLL", "dHL"))
#need to configure distances so less gaps between plots
tiff('ChlCar_boxplots_Envlight.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
par(mfrow = c(3, 4),  mar=c(2,4.5,2,2))
boxplot(Chl_a_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 10), xlab = "",
        ylab = "Chl a (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 10, "LL", cex=2)
text(2, 10, "p = 0.98", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("A.", cex.main=2, adj=0)
stripchart(Chl_a_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl_a_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 10), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
text(1, 9, "*", cex=3.5)
text(0.8, 10, "HL", cex=2)
text(2, 10, "p = 0.04", cex=1.5)
stripchart(Chl_a_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl_b_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 3), xlab = "",
        ylab = "Chl b (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 3, "LL", cex=2)
text(2, 3, "p = 0.71", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("B.", cex.main=2, adj=0)
stripchart(Chl_b_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl_b_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 3), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
text(1, 2.6, "*", cex=3.5)
text(0.8, 3, "HL", cex=2)
text(2, 3, "p = 0.04", cex=1.5)
stripchart(Chl_b_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Carot_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 2.5), xlab = "",
        ylab = "Total carotenoids (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 2.5, "LL", cex=2)
text(2, 2.5, "p = 0.87", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("C.", cex.main=2, adj=0)
stripchart(Carot_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Carot_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
text(1, 2.3, "*", cex=3.5)
text(0.8, 2.5, "HL", cex=2)
text(2, 2.5, "p = 0.03", cex=1.5)
stripchart(Carot_HL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl_ab_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0, 12), xlab = "",
        ylab = "Total Chla+b (mg/g)",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 12, "LL", cex=2)
text(2, 12, "p = 0.94", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("D.", cex.main=2, adj=0)
stripchart(Chl_ab_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chl_ab_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 12), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
text(1, 10, "*", cex=3.5)
text(0.8, 12, "HL", cex=2)
text(2, 12, "p = 0.04", cex=1.5)
stripchart(Chl_ab_HL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chlab_rat_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(3,4.5), xlab = "",
        ylab = "Chla:b",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 4.5, "LL", cex=2)
text(2, 4.5, "p = 0.21", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("E.", cex.main=2, adj=0)
stripchart(Chlab_rat_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Chlab_rat_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(3, 4.5), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 10, "*", cex=3.5)
text(0.8, 4.5, "HL", cex=2)
text(2, 4.5, "p = 0.18", cex=1.5)
stripchart(Chlab_rat_HL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Car_rat_LL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0, 
        ylim=c(0,0.5), xlab = "",
        ylab = "Car:Chl",
        cex.lab=2.0, cex.axis=1.5
)
text(0.8, 0.5, "LL", cex=2)
text(2, 0.5, "p = 0.63", cex=1.5)
#text(2.1, 1.45, "p = 0.57", cex=2)
title("F.", cex.main=2, adj=0)
stripchart(Car_rat_LL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

boxplot(Car_rat_HL_mean ~ EnvLight, comb,
        col=c("#999999", "#FF6347"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 0.5), xlab = "",
        ylab = "",
        cex.lab=2.0, cex.axis=1.5
)
#text(1, 10, "*", cex=3.5)
text(0.8, 0.5, "HL", cex=2)
text(2, 0.5, "p = 0.56", cex=1.5)
stripchart(Car_rat_HL_mean ~ EnvLight, comb,             # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           cex = 0.8,
           col = 1,           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)
dev.off()

#JUST SHOW 2 SUMMARISE SAVE SPACE
par(mfrow = c(2, 4))
boxplot(Carot_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "Total Carotenoids (mg/g) (LL)",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Carot_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 2.5), xlab = "",
        ylab = "Total Carotenoids (mg/g) (HL)",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 2.3, "*", cex=2.5)
boxplot(Chl_ab_LL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 12), xlab = "",
        ylab = "Total ChlA+B (mg/g) (LL)",
        cex.lab=1.5, cex.axis=1.0
)

boxplot(Chl_ab_HL_mean ~ EnvLight, comb,
        col=c("#FF6347", "#999999"),
        las=1, whisklty=1, staplewex=0,
        ylim=c(0, 12), xlab = "",
        ylab = "Total ChlA+B (mg/g) (HL)",
        cex.lab=1.5, cex.axis=1.0
)
text(2, 10, "*", cex=2.5)

#basic plots
#chl vars between treatments and env light
par(mfrow = c(1, 4))
#NEED TO MAKE SCALES SAME
boxplot(Chl_a_LL_mean~EnvLight, comb, ylim = c(0, 10))
boxplot(Chl_a_HL_mean~EnvLight, comb, ylim = c(0, 10))
boxplot(Chl_b_LL_mean~EnvLight, comb, ylim = c(0, 2.5))
boxplot(Chl_b_HL_mean~EnvLight, comb, ylim = c(0, 2.5))
boxplot(Carot_LL_mean~EnvLight, comb, ylim = c(0, 2.5))
boxplot(Carot_HL_mean~EnvLight, comb, ylim = c(0, 2.5))
boxplot(Chl_ab_LL_mean~EnvLight, comb, ylim = c(0, 12))
boxplot(Chl_ab_HL_mean~EnvLight, comb, ylim = c(0, 12))

boxplot(Chlab_rat_HL_mean~EnvLight, comb) #sig with 24 obs
#ll accessions making more chl a, b, tot and carot in hl than hl

#test
plot(Car_rat_HL_mean~PRI_HL_mean, data=comb)

corr <- cor.test(comb$Car_rat_HL_mean, comb$PRI_HL_mean,
                 method = "pearson"
)
corr$p.value #sig
corr$estimate #weak neg

corr <- cor.test(comb$Carot_HL_mean, comb$PRI_HL_mean,
                 method = "pearson"
)
corr$p.value #sig
corr$estimate #weak neg
#not as good as using ratio obs 24, better for obs 20

#nicer plot
library(ggplot2)
library(ggpubr)
library(basicTrendline)

#PRI CORRELATES WELL WITH CAR RAT
par(mfrow = c(2, 1))
ggplot(comb, aes(PRI_HL_mean, Car_rat_HL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in HL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.01, label.y = 0.4)

ggplot(comb, aes(PRI_LL_mean, Car_rat_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.4)

#try by carot mean
ggplot(comb, aes(PRI_HL_mean, Carot_HL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in HL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.01, label.y = 0.4)

ggplot(comb, aes(PRI_LL_mean, Carot_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in HL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.01, label.y = 0.4)


#try by carot mean
ggplot(comb, aes(PRI_LL_mean, Chl_ab_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.4)

ggplot(comb, aes(PRI_HL_mean, Chl_ab_HL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.01, label.y = 0.4)

ggplot(comb, aes(PRI_HL_mean, mean.HL.QYmax, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = -0.01, label.y = 0.4)

#try by carot mean
ggplot(comb, aes(PRI_LL_mean, Carot_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("PRI value") +
  ylab("Ratio carotenoids in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0, label.y = 0.4)

#good highest
ggplot(comb, aes(NDVI_LL_mean, FW.mean.LL, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value LL") +
  ylab("Fresh weight biomass in LL (mg)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.72, label.y = 100)


#good pos
ggplot(comb, aes(NDVI_HL_mean, fqfm_L5_mean_HL, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value HL") +
  ylab("Photosynthetic efficiency HL (365 umol m-2 s-1)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.73, label.y = 0.2)

#good neg
ggplot(comb, aes(NDVI_LL_mean, RGR_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value") +
  ylab("RGR diff") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.7, label.y = 0.2)

#shit
ggplot(comb, aes(NDVI_HL_mean, FW.mean.HL, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value") +
  ylab("Total chlorophyll in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.72, label.y = 250)

#neg
ggplot(comb, aes(NDVI_LL_mean, RGR_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value") +
  ylab("Total chlorophyll in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.72, label.y = 1)

#shit
ggplot(comb, aes(NDVI_HL_mean, RGR_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value") +
  ylab("Total chlorophyll in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.72, label.y = 1)

#low correlation
ggplot(comb, aes(NDVI_LL_mean, Chl_ab_LL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value") +
  ylab("Total chlorophyll in LL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 1, label.y = 12)

#not useful
ggplot(comb, aes(NDVI_HL_mean, Chlab_rat_HL_mean, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("NDVI value") +
  ylab("Total ratio chlorophyll A/B in HL (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 1.5, label.y = 4)


plot(Car_rat_LL_mean~PRI_LL_mean, data=comb)

corr <- cor.test(comb$PRI_LL_mean, comb$Car_rat_LL_mean,
                 method = "pearson"
)
corr$p.value #sig
corr$estimate #weak neg

plot(NDVI_LL_mean~Chl_a_LL_mean, data=comb)
plot(NDVI_LL_mean~Chl_b_LL_mean, data=comb)

corr <- cor.test(comb$Chl_a_LL_mean, comb$NDVI_LL_mean,
                 method = "pearson"
)
corr$p.value
corr$estimate #not strong

write.csv(comb, "combined_24obs_oppcalcs.csv")

#fluo analysis per treatment
names(comb)
comb$EnvLight
#take hl away from ll calculated paramters NPQ and fvfm main parameters
comb <- comb %>% mutate(meanQYmax_diff = (mean.HL.QYmax - mean.LL.QYmax)) %>% 
  mutate(maxQYmax_diff = (maximum.HL.QYmax - maximum.LL.QYmax)) %>%
  mutate(meanNPQ_L11_diff = (mean.HL.NPQL11 - mean.LL.NPQL11)) %>%
  mutate(meanNPQ_L5_diff = (mean.HL.NPQL5 - mean.LL.NPQL5)) %>%
  mutate(meanNPQ_L3_diff = (mean.HL.NPQL3 - mean.LL.NPQL3)) %>%
  mutate(meanNPQ_L1_diff = (mean.HL.NPQL1 - mean.LL.NPQL1)) %>%
  mutate(meanfqfm_L11_diff = (fqfm_L11_mean_HL - fqfm_L11_mean_LL)) %>%
  mutate(meanfqfm_L5_diff = (fqfm_L5_mean_HL- fqfm_L5_mean_LL)) %>%
  mutate(meanfqfm_L3_diff = (fqfm_L3_mean_HL - fqfm_L3_mean_LL)) %>%
  mutate(meanfqfm_L1_diff = (fqfm_L1_mean_HL - fqfm_L1_mean_LL))

t.test(meanQYmax_diff~EnvLight, comb) #NOT SIG
t.test(meanNPQ_L11_diff~EnvLight, comb) #NOT SIG
t.test(meanNPQ_L5_diff~EnvLight, comb) #NOT SIG
t.test(meanNPQ_L3_diff~EnvLight, comb) #NOT SIG
t.test(meanNPQ_L1_diff~EnvLight, comb) #NOT SIG
t.test(meanfqfm_L11_diff~EnvLight, comb) #NOT SIG
t.test(meanfqfm_L5_diff~EnvLight, comb) #NOT SIG
t.test(meanfqfm_L3_diff~EnvLight, comb) #NOT SIG
t.test(meanfqfm_L1_diff~EnvLight, comb) #NOT SIG

#non of differences sig by env light

#effect of treatment QY max
ggplot(comb, aes(y=meanQYmax_diff, x=reorder(Accession, meanQYmax_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in QY max (HL - LL)") +
  xlab("Ecotype")

#effect of treatment NPQ max
ggplot(comb, aes(y=meanNPQ_L11_diff, x=reorder(Accession, meanNPQ_L11_diff),  fill=EnvLight)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme_classic() +
  ylab("Difference in NPQ max (HL - LL)") +
  xlab("Ecotype")

#possibly not useful as not just treatment effect poss keep seperate


#manually moved species to start as automation fails if in middle col

library(dplyr)

#run this to save to correct place
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Comb")

#save 24 obs
write.csv(comb, "combined_24obs.csv")
#save 20 obs and move end col to 2nd pos then run
write.csv(flall2, "combined_20obs.csv")
flall2 <- read.csv("combined_20obs_renamed.csv")

#ANOVA AUTOMATION FROM MATT K
#starts from col no, #accession and treat in pos 
#flall2 <- comb
flall2$Accession

#for 24 obs
flall2$EnvLight <- c("LL", "LL", "HL", "LL", "LL", "LL",
                   "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                   "HL", "HL", "LL", "LL", "LL", "LL", "HL", "HL",
                   "HL", "LL", "LL")

#not working properly just law column done against gr data
my_formula <- as.formula(paste(colnames(flall2)[312], "~", paste(colnames(flall2)[c(313)], collapse = "*")))

anova_output2 <- list()
for(i in 2:312){
  my_col <- colnames(flall2[i]) 
  my_formula <- as.formula(paste(my_col, "~", paste(colnames(flall2)[c(313)], collapse = "*")))
  sum_aov <- summary(aov(my_formula, flall2))
  sum_tab <- sum_aov[[1]] %>% as.data.frame()
  sum_tab$ind <- my_col
  anova_output2[[i]] <- sum_tab
  print(my_formula)
}

(anova_output2)


all_output2 <- do.call(rbind, anova_output2)

all_output2

#write.csv(all_output2, "combo_all+fluo+hs_noenv_anovaoutput.csv")
write.csv(all_output2, "combo_all+fluo+hs_env_anovaoutput.csv")

model=lm(flall2$RGR_mean.HL ~ flall2$EnvLight)
#model=lm(flall2$NPQ_L11 ~ flall2$Accession * flall2$Treatment)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'flall2$EnvLight', conf.level=0.95)
summary(TUKEY)
plot(TUKEY , las=1 , col="brown")

