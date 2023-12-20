setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Comb")

#diff number of observations test as corr plots + pca
#20 obs including environment + hl or ll classification
forpca <- read.csv("combined_24obs_REM.csv") #with accessions remove std devs, n
dat <- read.csv("combined_20obs_renamed.csv") # just env grouping
dat <- read.csv("combined_24obs_renamed_REMCOV.csv") # just env grouping, no env data
dat <- read.csv("combined_24obs_renamed_REMCOV_subset.csv") #after analysis just including sig
dat <- read.csv("combined_24obs_oppcalcs_reduced.csv") #after analysis just including sig
#newest
dat <- read.csv("combined_24obs_opp_calc_subset.csv") #after analysis of new, same var as above but fw, fdw extra
#FOR CORR
dat <- read.csv("combined_24obs_opp_calc_subset_gr+photosyn.csv") #gr + photosyn
dat <- read.csv("combined_24obs_opp_calc_subset_gr+photosyn_rename.csv") #gr + photosyn QYmax name change

names(dat) 
#View(dat)

#NAs in env variables in 20 obs

#corr plot

str(dat)
sapply(dat, class)

#gr <- dat

gr$Mar2021_LambdaD <- as.numeric(gr$Mar2021_LambdaD)
gr$Mar2021_LambdaP <- as.numeric(gr$Mar2021_LambdaP)
gr$Sum2021_LambdaD <- as.numeric(gr$Sum2021_LambdaD)
gr$Sum2021_LambdaP <- as.numeric(gr$Sum2021_LambdaP)
gr$Aut2021_LambdaP <- as.numeric(gr$Aut2021_LambdaP)
gr$Aut2021_LambdaD <- as.numeric(gr$Aut2021_LambdaD)
gr$Spr2022_LambdaP <- as.numeric(gr$Spr2022_LambdaP)
gr$Spr2022_LambdaD <- as.numeric(gr$Spr2022_LambdaD)
gr$Summer2022_LambdaP <- as.numeric(gr$Summer2022_LambdaP)
gr$Summer2022_LambdaD <- as.numeric(gr$Summer2022_LambdaD)
gr$Sum2021_PFD.B <- as.numeric(gr$Sum2021_PFD.B)
gr$Sum2021_PFD.B

library(dplyr)

numeric_cols <- lapply(dat, is.numeric) %>% unlist() %>% unname()
dat2 <- dat[, numeric_cols]
library(corrplot)
#?corrplot
cor(dat2)

#dat2 <- dat2[complete.cases(dat2), ] #severely reduced no obs

#create matrix

M<-cor(dat2,use='complete.obs')
#M<-cor(dat2)
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

#slightly better plot
plotcolour<- colorRampPalette(c("darkred", "gray88", "darkcyan"))(20) 

corrplot(M, type="upper", order="hclust",method = "circle",
         p.mat = p.mat, sig.level = 0.05, insig = "blank",tl.cex = 0.5, tl.col = 'black',col=plotcolour)

#to run nicer plot
plotcolour<- colorRampPalette(c("darkred", "gray88", "darkcyan"))(20)
#This is where I've created a colour scheme of my choosing. 
tiff('Physiology_rename_corrplot.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
corrplot(M, col=plotcolour, type="upper", order="hclust", # Add coefficient of correlation
         tl.col="black", tl.cex = 0.8, tl.offset = 0.5,
         #Text label color and rotation
         # Combine with significance
         p.mat = p.mat,
         sig.level = c(.001, .01, .05), pch.cex = 0.7, #add signficance stars
         insig = "label_sig", pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )
dev.off()
#tl.cex changes txt size and pch.cex in plot size

#SAVE MATRIX AS EXCEL FILE
#write.csv(M, "Comb_20obs_mat.csv")
#write.csv(M, "Comb_24obs_mat_REMCOV.csv")
write.csv(M, "Comb_24obs_mat_REMCOV_subset.csv")
write.csv(M, "Comb_24obs_mat_opp_calc.csv") #new
write.csv(M, "Comb_24obs_mat_opp_calc+GR+PHOTOSYN.csv") #new

#GROWTH RATE SIGNIFICANT RELATIONSHIPS
comb <-dat
#correlations as signified by var pca plot
#BEST CORRELATIN FOR PRI HL?
corr <- cor.test(comb$PRI_HL_mean, comb$mean.HL.QYmax,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$PRI_HL_mean, comb$Chl_ab_HL_mean,
                 method = "pearson"
)
corr$estimate #weak positive
corr$p.value #sig

corr <- cor.test(comb$PRI_HL_mean, comb$Car_rat_HL_mean,
                 method = "pearson"
)
corr$estimate #weak neg
corr$p.value #sig

corr <- cor.test(comb$PRI_HL_mean, comb$RGR_mean.HL,
                 method = "pearson"
)
corr$estimate #not as pos
corr$p.value #sig

corr <- cor.test(comb$PRI_LL_mean, comb$FW.mean.LL,
                 method = "pearson"
)
corr$estimate #weak positive
corr$p.value #sig

#good relationship with pigments

corr <- cor.test(comb$PRI_LL_mean, comb$Chl_ab_LL_mean,
                 method = "pearson"
)
corr$estimate #weak positive
corr$p.value #sig

corr <- cor.test(comb$PRI_LL_mean, comb$Car_rat_LL_mean,
                 method = "pearson"
)
corr$estimate #weak neg
corr$p.value #sig


#NDVI CORRELATED WELL WITH LL GR

corr <- cor.test(comb$NDVI_LL_mean, comb$FW.mean.LL,
                 method = "pearson"
)
corr$estimate #godd positive
corr$p.value #sig



#apparent high correlation from matrix
corr <- cor.test(comb$RGR_mean.LL, comb$Car_rat_LL_mean,
                 method = "pearson"
)
corr$estimate #weak neg
corr$p.value #sig

corr <- cor.test(comb$RGR_diff, comb$NDVI_LL_mean,
                 method = "pearson"
)
corr$estimate #neg
corr$p.value #sig

corr <- cor.test(comb$RGR_propdiff, comb$NDVI_LL_mean,
                 method = "pearson"
)
corr$estimate #strong neg
corr$p.value #sig

#rgr area due to greenness or NDVI at lowest light

corr <- cor.test(comb$ColRGR_HL_mean, comb$Carot_HL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$ColRGR_HL_mean, comb$Chl_ab_HL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$ColRGR_HL_mean, comb$Carot_HL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

#hl performance col nos due to amount of chl and carot and ratio

corr <- cor.test(comb$PRI_LL_mean, comb$Chl_ab_LL_mean,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$PRI_HL_mean, comb$mean.HL.QYmax,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$NDVI_HL_mean, comb$PRI_LL_mean,
                 method = "pearson"
)
corr$estimate #neg
corr$p.value #sig

corr <- cor.test(comb$NDVI_HL_mean, comb$fqfm_L11_mean_HL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$NDVI_HL_mean, comb$fqfm_L5_mean_HL,
                 method = "pearson"
)
corr$estimate #pos
corr$p.value #sig

corr <- cor.test(comb$NDVI_HL_mean, comb$fqfm_L3_mean_HL,
                 method = "pearson"
)
corr$estimate #not as pos
corr$p.value #sig


#pca for corr sig subset
#pca
names(dat)
ions_all <- dat #5 includes species in col 13
# turn accession into factor and classify data to include
sapply(ions_all, class)

ions_all$EnvLight <- as.factor(ions_all$EnvLight)
options(ggrepel.max.overlaps = Inf)
library(FactoMineR)
library(factoextra)

ions_acc <- ions_all[ ,c(2:43)] 
ions_acc.pca <- PCA(ions_acc, quali.sup=42) #cant exceed max col number
ions_acc <- ions_all[ ,c(2:47)] # selecting columns from csv
ions_acc.pca <- PCA(ions_acc, quali.sup=46) #cant exceed max col number
print(ions_acc.pca)
ions_acc.pca$eig #24 and 18% 1 and 2
#ions_acc.pca$quali.sup
head(ions_acc.pca) #values for each row in pc
print(summary(ions_acc.pca)) #shows pc contirbutions of eigens
write.csv(ions_acc.pca$x, "obs24_pca_data.csv") #doesnt save anything
write.csv(ions_acc.pca$x, "obs24_pca_dataopp_calc.csv") #doesnt save anything


#do this first to get pc columns, then visualise variable
#ions_acc.pca<- prcomp(ions_all[,c(3:152)], center=T,scale=T) #list can include mix mat and df
ions_acc.pca<- prcomp(ions_acc, center=T,scale=T) #list can include mix mat and df
#tutorial
#tutorial
str(ions_acc.pca)
library(vegan) #used for ellipses

mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

#FOR ACCESSION POINTS ON LANDSCAPE
png("24obs_sig_var-PCA_oppcalc.png", type="cairo", width = 650, height = 650) #save as png

plot(mypc[,1], mypc[,2], col = ions_all$EnvLight,
     las=1, xlab="PC1 (24%)", ylab="PC2 (18%)",
     pch=16, cex=1.5, ylim=c(-4,5)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("topright", pch=16, col=1:2, cex=0.6, c("LL", "HL")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
text(x=mypc[,1],y=mypc[,2], labels =ions_all$EnvLight, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom

dev.off() #needs to be off to save it
#dev.new()

#needed to change red to black
reordered_groups <- factor(ions_all$EnvLight, levels = c("LL","HL"))

#basic plot
biplot(ions_acc.pca)

fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=reordered_groups, palette=c("green2", "gold"), addEllipses=TRUE, ellipse.type="confidence")

#FOR VARIABLE COS 2 GRAPH
#plot cos 2 as bar graph #high = good representation on pc
fviz_cos2(ions_acc.pca, choice = "var", axes = 1:2)

fviz_pca_var(ions_acc.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#FOR 43 VAR FOR PAPER WRONG DIRECTION TO ABOVE
#FOR PCA
dat <- read.csv("combined_24obs_opp_calc_subset_FORPAPER.csv")
dat <- read.csv("combined_24obs_opp_calc_subset_FORPAPER_rename.csv") #qymax fvfm

#pca for corr sig subset
#pca
names(dat)

ions_all <- dat #5 includes species in col 13
# turn accession into factor and classify data to include
sapply(ions_all, class)

ions_all$EnvLight <- as.factor(ions_all$EnvLight)
options(ggrepel.max.overlaps = Inf)
#install.packages("FactoMineR")
#install.packages("factoextra")
library(FactoMineR)
library(factoextra)

ions_acc <- ions_all[ ,c(2:43)] # selecting columns from csv
ions_acc.pca <- PCA(ions_acc, quali.sup=42) #cant exceed max col number
print(ions_acc.pca)
ions_acc.pca$eig #24 and 19% 1 and 2
#ions_acc.pca$quali.sup
head(ions_acc.pca) #values for each row in pc
print(summary(ions_acc.pca)) #shows pc contirbutions of eigens

ions_acc.pca<- prcomp(ions_acc, center=T,scale=T) #list can include mix mat and df
#tutorial
str(ions_acc.pca)
#install.packages("vegan")
library(vegan) #used for ellipses

mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

#basic plot
biplot(ions_acc.pca)

fviz_pca_biplot(ions_acc.pca, repel=TRUE, pointsize=6, pointshape=21, col.var="black", arrowsize=0.6, labelsize=5, col.ind=reordered_groups, palette=c("green2", "gold"), addEllipses=TRUE, ellipse.type="confidence")

#FOR VARIABLE COS 2 GRAPH
#plot cos 2 as bar graph #high = good representation on pc
fviz_cos2(ions_acc.pca, choice = "var", axes = 1:2)

fviz_pca_var(ions_acc.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE# Avoid text overlapping
)

#try to adjust easthetics
#rename column names and make txt smaller so all ledgible

library(ggplot2)
pdf('24obs_PCA_variables_renamedqymax_toedit.pdf', width=14, height=12)
plot1 <- fviz_pca_var(ions_acc.pca, col.var = "cos2",
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE)
plot1 + xlab("PC1 (24%)") + ylab("PC2 (19%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     plot.title = element_blank())
dev.off()
#remove title and grids

#save plot 1
ggsave("24obs_biplot.png", dpi = 300, width = 15, height = 15 , units = "cm") 
ggsave("24obs_biplot_renameQYmax_toedit.tiff", dpi = 300, width = 15, height = 15 , units = "cm")

#for accession labels
png("24obs_sig_var_accessions-PCA_opp_calc_FORPAPER_rename.png", type="cairo", width = 650, height = 650) #save as png

plot(mypc[,1], mypc[,2], col = reordered_groups,
     las=1, xlab="PC1 (24%)", ylab="PC2 (18%)",
     pch=16, cex=1.5, ylim=c(-5,5)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend('topright', legend = levels(ions_all$EnvLight), col = 2:1, cex = 1.2, pch = 16)
text(x=mypc[,1],y=mypc[,2], labels =forpca$Accession, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom

dev.off() #needs to be off to save it
#dev.new()

#pca with centroids
par(mfrow = c(1, 1))
library(dplyr)
#iris_pca <- prcomp(sub[, 1:26], center = T, scale. = T)
scoresx <- ions_acc.pca$x
scores2 <- scoresx %>% data.frame() #change to dataframe

scores2$EnvLight <- ions_all$EnvLight #append envlight variable
centroids <- scores2 %>% group_by(EnvLight) %>% summarise_all("mean")

tiff('24obs_PCA_accessions_centroids_renamedqymax.tiff', units="in", width=14, height=12, res=300, compression = 'lzw')
pdf('24obs_PCA_accessions_centroids_renamedqymax_toedit.pdf', width=14, height=12)
#PLOT PCA WITH CENTROIDS
plot(scores2[,1], scores2[,2],
     las=1, xlab="PC1 (24%)", ylab="PC2 (18%)",
     pch=16, cex=3, xlim=c(-10,10), ylim=c(-10,10))
loadings = TRUE
abline(h = 0, col = "lightgrey", lty = 2)
abline(v = 0, col = "lightgrey", lty = 2)
text(x=scores2[,1],y=scores2[,2], labels =forpca$Accession, cex=2, pos=2)

points(centroids[centroids$EnvLight == "LL", 2], centroids[centroids$EnvLight == "LL", 3], pch = 3, lwd = 4, cex = 5)
points(centroids[centroids$EnvLight == "HL", 2], centroids[centroids$EnvLight == "HL", 3], pch = 3, lwd = 4, cex = 5, col = "red")

points(scores2[scores2$EnvLight == "LL", 1], scores2[scores2$EnvLight == "LL", 2], pch = 16, cex = 3)
points(scores2[scores2$EnvLight == "HL", 1], scores2[scores2$EnvLight == "HL", 2], pch = 16, cex = 3, col = "red")

legend("topright", pch=16, col=c("black", "red"), cex=3.0, c("dLL", "dHL"), bty="n")

dev.off()
#remake pca plots with new xanthophyll data
#see comb_join_xanth_data_to_summary_14_script