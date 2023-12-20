#combine GR and Col count data
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

library(dplyr)
library(tidyr)

#read in averaged data
gr <- read.csv("GR_RGR_propchange.csv")
cc <- read.csv("Col_RGR_propchange.csv")
bio <- read.csv("Biomass_RGR_propchange.csv")

#join 2 data sets together (by accession)
names(gr)
names(gr)[2] <- "Sample" #name col for observation names
names(cc)[2] <- "Sample"
names(bio)[2] <- "Sample"
gr %>% str()
cc %>% str()
gr$Sample <- as.numeric(gr$Sample) #growth rate
cc$Sample <- as.numeric(cc$Sample) # col counts
bio$Sample <- as.numeric(bio$Sample) # biomass FW FDW
#View(gr$Sample) #named numeric 1-28
tail(data.frame(gr), 20) #show last 20 rows
gr <- gr[1:28, ] #not done anything

joined <- gr %>% inner_join(cc, by = c("Sample"))
#worked but not done exactly as code says

#ADD IN WHEN WORKING WITH BIO DATA
#new add
joined <- joined %>% inner_join(bio, by = c("Sample"))

#measured these before added bio
#relationships between col gain and area light responses
#both worked out using LL - HL 
names(joined)

plot(RGR_diff.x ~ RGR_diff.y, joined)
cor(joined$RGR_diff.x, joined$RGR_diff.y,
    method = "spearman"
    )
#0.626 
#when diff area increases, diff cols increases

sapply(class, joined)
names(joined)

#correlation between FW and FDW
corr <- cor.test(joined$FW_mean_diff, joined$FDW_mean_diff,
                 method = "spearman"
)
corr$p.value
#0.414 p value 0.02

#correlation area and FW
corr <- cor.test(joined$FW_mean_diff, joined$RGR_diff.x,
                 exact = FALSE,
                 method = "spearman"
)
corr$p.value
#0.646 p value 0.00019

#correlation cols and FW
corr <- cor.test(joined$FW_mean_diff, joined$RGR_diff.y,
                 exact = FALSE,
                 method = "spearman"
)
corr$p.value
#0.530 p value 0.003

#correlation area and FDW
corr <- cor.test(joined$FDW_mean_diff, joined$RGR_diff.x,
                 exact = FALSE,
                 method = "spearman"
)
corr$p.value
#0.432 p value 0.02

#correlation cols and FDW
corr <- cor.test(joined$FDW_mean_diff, joined$RGR_diff.y,
                 exact = FALSE,
                 method = "spearman"
)
corr$p.value
#0.063 p value 0.747

corr <- cor.test(joined$RGR_diff.x, joined$RGR_diff.y,
                 exact = FALSE,
                 method = "spearman"
                 )
corr$p.value
#0.626 p value 0.0003

#area and cols correlate
#area and FW correlate
#weaker
#cols and fw
#fdw and cols
#fdw and area
#fdw and fd

corr <- cor.test(joined$mean.HL.x, joined$mean.LL.y,
                 method = "spearman")
corr$p.value
#0.512 and 0.005
#higher mean area in HL, higher area in HL

library(corrplot)
#basic corr plot, first check class, select only numeric, then plot matrix out
sapply(joined, class)
sapply(joined, is.numeric)
num <- select_if(joined, is.numeric)
cornum <- cor(num) #displays r2

corrplot(cor(num),
         method = "number",
         type = "upper" # show only upper side
)

plot(mean.HL.x ~ mean.HL.y, joined)
plot(mean.LL.x ~ mean.LL.y, joined) #no real corr
plot(RGR_diff.x ~ RGR_diff.y, joined)
plot(maximum.HL.x ~ maximum.HL.y, joined)
plot(maximum.LL.x ~ maximum.LL.y, joined)

joined$Accession <- joined$Accession.x
#joined$Accession.y

library(ggplot2)
library(ggpubr)
#plot scatter plots with labels
#NICER PLOTS FOR SIG AREA WITH COL GAIN AND AREA WITH FW
#difference in area plot difference in colonies
ggplot(joined, aes(RGR_diff.x, RGR_diff.y, label = Accession)) +    # ggplot2 plot with labels
    geom_point() +
    geom_text(aes(label = Accession), hjust = - 0.5) +
    geom_smooth(method='lm') +
    xlab("RGR diff area (mm2)") +
    ylab("RGR diff col gain") +
    theme_bw() +
    theme_classic() +
    stat_cor(method = "spearman", label.x = 0.08, label.y = 6)

#difference in area plot difference in colonies
ggplot(joined, aes(RGR_diff.x, FW_mean_diff, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("RGR diff area (mm2)") +
  ylab("RGR diff fresh weight (mg)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "spearman", label.x = -0.25, label.y = 2000)

#mean HL area vs mean HL cols
ggplot(joined, aes(mean.HL.x, mean.HL.y, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#max HL area vs mean HL cols
ggplot(joined, aes(maximum.HL.x, maximum.HL.y, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5)

#area vs FW diff
ggplot(joined, aes(RGR_diff.x, FW_mean_diff, label = Accession.y)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession.y), hjust = - 0.5)

#col vs FW diff
ggplot(joined, aes(RGR_diff.y, FW_mean_diff, label = Accession.y)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession.y), hjust = - 0.5)

#area vs FDW diff
ggplot(joined, aes(RGR_diff.x, FDW_mean_diff, label = Accession.y)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession.y), hjust = - 0.5)

#col vs FDW diff
ggplot(joined, aes(RGR_diff.y, FDW_mean_diff, label = Accession.y)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession.y), hjust = - 0.5)

#not useful
#ggplot(joined, aes(propchange4, RGR_propdiff, label = Accession)) +    # ggplot2 plot with labels
#  geom_point() +
#  geom_text(aes(label = Accession), hjust = - 0.5)

#stronger treatment effect between methods if same accessions top/bottom

#FOR RANKING INDIVIDUALS FOR EACH VARIABLE
names(joined)
#x is area
mean.HL.x_ord <- joined %>% arrange(desc(mean.HL.x)) %>% select(Accession) %>% top_n(28)
mean.LL.x_ord <- joined %>% arrange(desc(mean.LL.x)) %>% select(Accession) %>% top_n(28)
maximum.HL.x_ord <- joined %>% arrange(desc(maximum.HL.x)) %>% select(Accession) %>% top_n(28)
maximum.LL.x_ord <- joined %>% arrange(desc(maximum.LL.x)) %>% select(Accession) %>% top_n(28)
RGR_diff.x_ord <- joined %>% arrange(desc(RGR_diff.x)) %>% select(Accession) %>% top_n(28)
max.diff.x_ord <- joined %>% arrange(desc(max_diff.x)) %>% select(Accession) %>% top_n(28)
propchange4.x_ord <- joined %>% arrange(desc(propchange4)) %>% select(Accession) %>% top_n(28)

#y is col counts
mean.HL.y_ord <- joined %>% arrange(desc(mean.HL.y)) %>% select(Accession) %>% top_n(28)
mean.LL.y_ord <- joined %>% arrange(desc(mean.LL.y)) %>% select(Accession) %>% top_n(28)
maximum.HL.y_ord <- joined %>% arrange(desc(maximum.HL.y)) %>% select(Accession) %>% top_n(28)
maximum.LL.y_ord <- joined %>% arrange(desc(maximum.LL.y)) %>% select(Accession) %>% top_n(28)
RGR_diff.y_ord <- joined %>% arrange(desc(RGR_diff.y)) %>% select(Accession) %>% top_n(28)
max.diff.y_ord <- joined %>% arrange(desc(max_diff.y)) %>% select(Accession) %>% top_n(28)
RGR_propdiff_ord <- joined %>% arrange(desc(RGR_propdiff)) %>% select(Accession) %>% top_n(28)

#bio data
FW_mean_diff_ord <- joined %>% arrange(desc(FW_mean_diff)) %>% select(Accession.y) %>% top_n(28)
FW_max_diff_ord <- joined %>% arrange(desc(FW_max_diff)) %>% select(Accession.y) %>% top_n(28)
FDW_mean_diff_ord <- joined %>% arrange(desc(FDW_mean_diff)) %>% select(Accession.y) %>% top_n(28)
FDW_max_diff_ord <- joined %>% arrange(desc(FDW_max_diff)) %>% select(Accession.y) %>% top_n(28)

#combo
adensity_HL_ord <- density %>% arrange(desc(adensity_HL)) %>% select(Accession.y) %>% top_n(28)
adensity_LL_ord <- density %>% arrange(desc(adensity_LL)) %>% select(Accession.y) %>% top_n(28)
Cdensity_HL_ord <- density %>% arrange(desc(Cdensity_HL)) %>% select(Accession.y) %>% top_n(28)
Cdensity_LL_ord <- density %>% arrange(desc(Cdensity_LL)) %>% select(Accession.y) %>% top_n(28)

#stick all together with rankings next to them - see consistency
class(mean.HL.x_ord) #df
rankings <- cbind(mean.HL.x_ord, mean.LL.x_ord, maximum.HL.x_ord, maximum.LL.x_ord,
      RGR_diff.x_ord, max.diff.x_ord, propchange4.x_ord, mean.HL.y_ord, 
      mean.LL.y_ord, maximum.HL.y_ord, maximum.LL.y_ord,
      RGR_diff.y_ord, max.diff.y_ord, RGR_propdiff_ord, FW_mean_diff_ord, 
      FW_max_diff_ord, FDW_mean_diff_ord, FDW_max_diff_ord, adensity_HL_ord,
      adensity_LL_ord, Cdensity_HL_ord, Cdensity_LL_ord)

class(rankings)

colnames(rankings) <- c("mean.HL.x_ord", "mean.LL.x_ord", "maximum.HL.x_ord", "maximum.LL.x_ord",
                            "RGR_diff.x_ord", "max.diff.x_ord", "propchange4.x_ord", "mean.HL.y_ord", 
                            "mean.LL.y_ord", "maximum.HL.y_ord", "maximum.LL.y_ord",
                            "RGR_diff.y_ord", "max.diff.y_ord", "RGR_propdiff_ord",
                        "FW_mean_diff_ord", "FW_max_diff_ord", "FDW_mean_diff_ord", "FDW_max_diff_ord",
                        "adensity_HL_ord", "adensity_LL_ord", "Cdensity_HL_ord", "Cdensity_LL_ord")
#need to add cols for 1-28 between each column

#density = divide mass per area
density <- joined %>% mutate(adensity_HL = (FW.mean.HL/mean.HL.x))
density %>% arrange(desc(adensity_HL)) %>% select(Accession.y) %>% top_n(28)
density <- density %>% mutate(adensity_LL = (FW.mean.LL/mean.LL.x))
density %>% arrange(desc(adensity_LL)) %>% select(Accession.y) %>% top_n(28)

#density = divide mass per col
density <- density %>% mutate(Cdensity_HL = (FW.mean.HL/mean.HL.y))
density %>% arrange(desc(Cdensity_HL)) %>% select(Accession.y) %>% top_n(28)
density <- density %>% mutate(Cdensity_LL = (FW.mean.LL/mean.LL.y))
density %>% arrange(desc(Cdensity_LL)) %>% select(Accession.y) %>% top_n(28)

#density upgraded version of joined - extra variables to correlate

#density (biomass per col) in hl and ll
ggplot(density, aes(Cdensity_HL, Cdensity_LL, label = Accession.y)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession.y), hjust = - 0.5)

ggplot(density, aes(adensity_HL, adensity_LL, label = Accession.y)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession.y), hjust = - 0.5)

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
RGR_diff_area <- (joined[ ,11])
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

RGR_diff_area <- (joined[ ,11])
FW_mean_diff <- (joined[ ,52])
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

