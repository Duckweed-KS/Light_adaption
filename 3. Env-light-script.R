#script to explore enviro light data

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

gr <- read.csv("Env-light.csv")
#gr <- read.csv("Env-light+coverage.csv")

tail(gr)
head(gr)
names(gr)

#check classes of all cols
sapply(gr, class)

gr$Mar2021_LambdaD <- as.numeric(gr$Mar2021_LambdaD)
gr$Mar2021_LambdaP <- as.numeric(gr$Mar2021_LambdaP)
gr$Sum2021_LambdaD <- as.numeric(gr$Sum2021_LambdaD)
gr$Sum2021_LambdaP <- as.numeric(gr$Sum2021_LambdaP)
gr$Spr2021_LambdaP <- as.numeric(gr$Spr2021_LambdaP)
gr$Spr2021_LambdaD <- as.numeric(gr$Spr2021_LambdaD)
gr$Aut2021_LambdaP <- as.numeric(gr$Spr2022_LambdaP)
gr$Spr2022_LambdaD <- as.numeric(gr$Spr2022_LambdaD)
gr$Summer2022_LambdaP <- as.numeric(gr$Summer2022_LambdaP)
gr$Summer2022_LambdaD <- as.numeric(gr$Summer2022_LambdaD)

names(gr)
length(gr)

gr$Summer2022_LambdaD
gr$Summer2022_LambdaP

#SUMMARY STATS
gr$Duckweed.collection
gr$Light <- c("LL","LL", "LL", "HL", "LL", "LL",
                      "HL", "HL", "HL", "LL", "HL", "LL", "LL",
                      "HL", "HL", "LL", "LL", "LL", "LL", "HL")
as.factor(gr$Light)

library(dplyr)
names(gr)
#not working
gr <- gr[complete.cases(gr), ]  
gr %>% select(-Duckweed.collection, -Name) %>% filter(Light %in% "HL") %>% summarise_all(range)
agg_df <- aggregate(gr$Summer2022_PFD.FR, by=list(gr$Light), FUN=range)
agg_df1 <- aggregate(gr$Summer2022_light, by=list(gr$Light), FUN=range)
agg_df2 <- aggregate(gr$Summer2022_PFD.R, by=list(gr$Light), FUN=range)
agg_df3 <- aggregate(gr$Summer2022_PFD.G, by=list(gr$Light), FUN=range)
agg_df4 <- aggregate(gr$Summer2022_PFD.B, by=list(gr$Light), FUN=range)
agg_df5 <- aggregate(gr$Summer2022_PFD.UV, by=list(gr$Light), FUN=range)
agg_df6 <- aggregate(gr$Summer2022_LambdaD, by=list(gr$Light), FUN=range)
agg_df7 <- aggregate(gr$Summer2022_LambdaP, by=list(gr$Light), FUN=range)
agg_df8 <- aggregate(gr$Summer2022_PPFD, by=list(gr$Light), FUN=range)

dob <- range(gr$Summer2022_light)

#plot as PCA
#plot gr RGR data as PCA
#works for RGR but doesnt really show much - NEED MORE DIMENSIONS?
ions_all <- gr[, c(1:47)] #5 includes species in col 13

sapply(ions_all, class)

ions_all$Duckweed.collection <- as.factor(ions_all$Duckweed.collection)
ions_all$Name <- as.factor(ions_all$Name)

#not working to group by species 
class(ions_all$Duckweed.collection)

library(FactoMineR)
#install.packages("missMDA")
library(missMDA)

ions_acc <- ions_all[ ,c(3:47)] # selecting columns from csv
#just numeric variabeles

#need to remove nas from dataset
#data(gr)
nb = estim_ncpPCA(ions_acc,ncp.max=5)
res.comp = imputePCA(ions_acc,ncp=2)
res.pca = PCA(res.comp$completeObs)

ions_acc.pca <- PCA(ions_acc, quali.sup=45) #cant exceed max col number
print(ions_acc.pca)
head(ions_acc.pca)
print(summary(ions_acc.pca)) #shows pc contirbutions of eigens
#dim 1 54% dim 2 15% 
#write.csv(ions_acc.pca, ions_app_pca_data.csv) # doesnt work

print(res.pca)
head(res.pca)
print(summary(res.pca))
#write.csv(res.pca, ions_app_pca_data.csv) # doesnt work
#dim 1 63.4% dim 2 12.4% explains 75% variance

#do this first to get pc columns, then visualise variable
#ions_acc.pca<-prcomp(ions_all[ ,3:47],center=T,scale=T) #list can include mix mat and df
#error here as NAs in dataset
ions_acc.pca<-prcomp(na.omit(ions_acc), center = TRUE, scale = TRUE)
ions_acc.pca<-prcomp(center = TRUE, scale = TRUE)

ions_all$Duckweed.collection
ions_all$Duckweed.collection <- factor(ions_all$Duckweed.collection, levels=c("KS01", "KS02", "KS03", "KS04", "KS06", "KS09", "KS12", "KS13", "KS14", "KS15",
                                                          "KS16", "KS17", "KS18", "KS20", "KS21", "KS22", "KS25", "KS27",
                                                          "KS28", "KS29"))


#tutorial
str(ions_acc.pca)
ions_acc.pca$x
mypc <- ions_acc.pca$x #define new variable x = pcs just need to plot these on xy graph

png("ionomebyplace.png", type="cairo", width = 900, height = 900) #save as png

plot(mypc[,1], mypc[,2], col = ions_all$Light,
     las=1, xlab="PC1 (63.4%)", ylab="PC2 (12.4%)",
     pch=16, cex=1.5, ylim=c(-5,5)) #pcs as columns, produce xy plot, las is rotation of axis numbers, pch plot shape, ylim expand out so legend room
abline(v=0, lty=2, col="lightgrey") #draw line, lty is segmented
abline(h=0, lty=2, col="lightgrey") #0 lines dashed
legend("topright", pch=16, col=1:2, cex=0.6, c("HL", "LL")) #customise legend seperately cex=txtsize cols 1-4 stnd, concat order of places as want to display
text(x=mypc[,1],y=mypc[,2], labels =ions_all$Duckweed.collection, pos=2) #txt and pos to define labels on plot, 2 = top of plot, 1 bottom

dev.off() #needs to be off to save it
#dev.new()

#PLOT AS DENDROGRAPH
gr <- read.csv("Env-light_noKS01.csv")
gr <- read.csv("Env-light_noKS01+pfd.csv")
ions_all <- gr[, c(1:47)] #5 includes species in col 13
ions_all <- gr[, c(1:52)] #5 includes species in col 13 ??

#remove ks01 manually - still gives 20 levels
ions_all = ions_all[-1,]
ions_all$Duckweed.collection <- factor(ions_all$Duckweed.collection, levels=c("KS02", "KS03", "KS04", "KS06", "KS09", "KS12", "KS13", "KS14", "KS15",
                                                                              "KS16", "KS17", "KS18", "KS20", "KS21", "KS22", "KS25", "KS27",
                                                                              "KS28", "KS29"))
ions_acc <- ions_all[ ,c(3:47)] # selecting num columns from csv
ions_acc <- ions_all[ ,c(3:52)] #??

cluster <- ions_all
clusterno <- ions_acc
str(cluster)

typeof(cluster)
class(cluster)
as.matrix(cluster)

cluster$Duckweed.collection

#install.packages('tidyverse')
#install.packages('cluster')
#install.packages('factoextra')
#install.packages('dendextend')

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

library (dplyr)

clust <- t(cluster)

rownames(cluster) <- cluster$Duckweed.collection
distance <- get_dist(cluster)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

dist.obj <- get_dist(cluster, method = "euclidean", stand = FALSE)
fviz_dist(
        dist.obj,
        order = TRUE,
        show_labels = TRUE,
        lab_size = NULL,
        gradient = list(low = "blue", mid = "white", high = "red")
)

clusterom <- na.omit(clusterno) # removes nas
clustersc <- scale(clusterno) # standardize variables

clustersc

#par(mfrow = c(1, 2))
# Dissimilarity matrix using euclidean
clusterd <- dist(clustersc, method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(clusterd, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, rownames(cluster))
plot(hc1) #basic
summary(hc1) #check labels
#4 now in ll group according to this method

# Dissimilarity matrix using euclidean
clusterd <- dist(clustersc, method = "euclidean")
# Hierarchical clustering using ward Linkage
hc1a <- hclust(clusterd, method = "ward.D" )
# Plot the obtained dendrogram
plot(hc1a, cex = 0.6, hang = -1, rownames(cluster))

#same groupings for euclidean + complete and ward.D linkage

# Dissimilarity matrix using euclidean
clusterd <- dist(clustersc, method = "euclidean")
# Hierarchical clustering using ward2 Linkage
hc1a2 <- hclust(clusterd, method = "ward.D2" )
# Plot the obtained dendrogram
plot(hc1a2, cex = 0.7, hang = -1, rownames(cluster))
rect.hclust(hc1a2, k = 2, which = NULL, x = NULL, h = NULL,
            border = c("blue", "red"), cluster = NULL)


cutree(hc1a2, k = 2, h = NULL) # vector for group memberships

#REPEAT AND CHANGE TO COMPLETE
# Dissimilarity matrix using euclidean
clusterd <- dist(clustersc, method = "euclidean")
# Hierarchical clustering using ward2 Linkage
hc1a2 <- hclust(clusterd, method = "complete" )
# Plot the obtained dendrogram
plot(hc1a2, cex = 0.7, hang = -1, rownames(cluster))
rect.hclust(hc1a2, k = 2, which = NULL, x = NULL, h = NULL,
            border = c("blue", "red"), cluster = NULL)


cutree(hc1a2, k = 2, h = NULL) # vector for group memberships

# other methods to define dissimilarity
pars <- par()
par(lwd=3, mar=c(3.2,4,2,0))
# Dissimilarity matrix using manhattan
clusterd2 <- dist(clustersc, method = "manhattan")
# Hierarchical clustering using Complete Linkage
hc2 <- hclust(clusterd2, method = "complete" )
hc2$labels <- cluster$Duckweed.collection
# Plot the obtained dendrogram
#not a ggplot so not saving as one
plot(hc2, cex = 1.0, hang = -1, rownames(cluster))
rect.hclust(hc2, k = 2, which = NULL, x = NULL, h = NULL,
            border = c("blue", "red"), cluster = NULL)

#try to plot using ggplot
#install.packages("ggdendro")
library(ggplot2)
library(ggdendro)
#best plot for pub
ggdendrogram(hc2, rotate = FALSE, size = 42, 
             labels = TRUE, segments = TRUE) +
        theme(axis.text=element_text(size=16)) +
        annotate("rect", xmin = 0.7, xmax = 11.1, ymin = -5, ymax = 105,
                 color = "blue", size = 1.5, alpha = .3,fill = "white") +
        annotate("rect",xmin = 11.2, xmax = 19.2, ymin = -5, ymax = 105,
                 color = "red", size = 1.5, alpha= .3, fill = "white") +
        #geom_rect(aes(xmin = 11.2, xmax = 19.2, ymin = -5, ymax = 105,
        #              color = "red", alpha= .3, fill = "white")) +
        theme(legend.position = "none")
ggsave("hc2.png", dpi = 300, width = 15, height = 10 , units = "cm") 

#high quality save

cutree(hc2, k = 2, h = NULL) # vector for group memberships

#no borders
par(mfrow = c(1,1)) 
pars <- par()
par(lwd=3, mar=c(3.2,4,2,0))
# Dissimilarity matrix using manhattan
clusterd2 <- dist(clustersc, method = "manhattan")
# Hierarchical clustering using Complete Linkage
hc2 <- hclust(clusterd2, method = "complete" )
# Plot the obtained dendrogram
plot(hc2, cex = 1.0, hang = -1, rownames(cluster))



#same 2 groups using manhattan

# Dissimilarity matrix using maximum
clusterd3 <- dist(clustersc, method = "maximum")
# Hierarchical clustering using Complete Linkage
hc3 <- hclust(clusterd3, method = "complete" )
# Plot the obtained dendrogram
plot(hc3, cex = 0.6, hang = -1, rownames(cluster))

# Dissimilarity matrix using canberra
clusterd4 <- dist(clustersc, method = "canberra")
# Hierarchical clustering using Complete Linkage
hc4 <- hclust(clusterd4, method = "complete" )
# Plot the obtained dendrogram
plot(hc4, cex = 0.6, hang = -1, rownames(cluster))

# Dissimilarity matrix using binary
clusterd5 <- dist(clustersc, method = "binary")
# Hierarchical clustering using Complete Linkage
hc5 <- hclust(clusterd5, method = "complete" )
# Plot the obtained dendrogram
plot(hc5, cex = 0.6, hang = -1, rownames(cluster))

# Dissimilarity matrix using minkowski
clusterd6 <- dist(clustersc, method = "minkowski")
# Hierarchical clustering using Complete Linkage
hc6 <- hclust(clusterd6, method = "complete" )
# Plot the obtained dendrogram
plot(hc6, cex = 0.6, hang = -1, rownames(cluster))
rotate(hc6)

#to try method to editable title
hclust(d, method = "complete", members = NULL)
# S3 method for hclust
plot(hc1, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram ss",
     sub = NULL, xlab = NULL, ylab = "Height aa")

#used one complete cluster
plot(hc1, cex = 0.7, hang = -1, rownames(cluster))
rect.hclust(hc1, k = 2, which = NULL, x = NULL, h = NULL,
            border = c("blue","red"), cluster = NULL)


cutree(hc1, k = 2, h = 11) # vector for group memberships


# load package ape; remember to install it: install.packages('ape')
#install.packages('ape')
library(ape)
# plot basic tree
plot(hc1, type = "phylogram", cex = 0.9, label.offset = 1, rownames(cluster))

#horiz plot instead
#not working with labels
plot(as.phylo(hc1), cex = 0.9, rownames(cluster))

#not worked to attach together using R
#read in file with coverage data attached to light
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")
gr <- read.csv("Env-light+coverage.csv")
gr <- read.csv("Env-light_noKS01.csv")
names(gr)
library(dplyr)

#ll and hl organised sites
gr$Light <- c("LL","LL", "LL", "HL", "LL", "LL",
              "HL", "HL", "HL", "LL", "HL", "LL", "LL",
              "HL", "HL", "LL", "LL", "LL", "LL", "HL")

#no ks01
gr$Light <- c("LL","LL", "HL", "LL", "LL",
              "HL", "HL", "HL", "LL", "HL", "LL", "LL",
              "HL", "HL", "LL", "LL", "LL", "LL", "HL")

as.factor(gr$Light)

#work out r:fr ratio
#2 diff files (with without KS01)
gr <- gr %>% mutate(Sum2022RFR = (Sum2022_PFD.R/Sum2022_PFD.FR))
t.test(Sum2022RFR ~ Light, data=gr)

#no ks01
gr <- gr %>% mutate(Sum2022RFR = (Summer2022_PFD.R/Summer2022_PFD.FR))
t.test(Sum2022RFR ~ Light, data=gr) #s
gr <- gr %>% mutate(Spr2022RFR = (Spr2022_PFD.R/Spr2022_PFD.FR))
t.test(Spr2022RFR ~ Light, data=gr) #s
gr <- gr %>% mutate(Aut2021RFR = (Aut2021_PFD.R/Aut2021_PFD.FR))
t.test(Aut2021RFR ~ Light, data=gr) #s
gr <- gr %>% mutate(Sum2021RFR = (Sum2021_PFD.R/Sum2021_PFD.FR))
t.test(Sum2021RFR ~ Light, data=gr) #s
gr <- gr %>% mutate(Spr2021RFR = (Spr2021_PFD.R/Spr2021_PFD.FR))
t.test(Spr2021RFR ~ Light, data=gr) #s

boxplot(Sum2022RFR ~ Light, gr)
boxplot(Spr2022RFR ~ Light, gr)
boxplot(Aut2021RFR ~ Light, gr)
boxplot(Sum2021RFR ~ Light, gr)
boxplot(Spr2021RFR ~ Light, gr)

library(corrplot)

#rename cols to match light season/month

#methods to correlate variables

corr <- cor.test(x=gr$Sum2022_cov, y=gr$Sum2022_PFD.FR, method = 'pearson')
corr
corr <- cor.test(x=gr$Sum2022_cov, y=gr$Sum2022_PFD.UV, method = 'pearson')
corr
plot(Sum2022_cov ~ Sum2022_PFD.FR, gr) #fair correlation
plot(Sum2022_cov ~ Sum2022_PFD.UV, gr) #fair correlation
plot(Sum2022_cov ~ Sum2022_light, gr)
plot(Sum2022_cov ~ Sum2022_LambdaP, gr) #doesnt work as 0/1 type data
plot(Sum2022_cov ~ Sum2022_LambdaD, gr)

#noticed from corr plot
corr <- cor.test(x=gr$Sum2022_cov, y=gr$Aut2020_cov, method = 'pearson')
corr$p.value
corr$estimate
#HIGH

corr <- cor.test(x=gr$Sum2022_cov, y=gr$Sum2022_PFD.FR, method = 'pearson')
corr$p.value
corr$estimate
#low

corr <- cor.test(x=gr$Sum2022_cov, y=gr$Sum2022_PFD.UV, method = 'pearson')
corr$p.value
corr$estimate
#low

corr <- cor.test(x=gr$Sum2022_cov, y=gr$Sum2022_light, method = 'pearson')
corr$p.value
corr$estimate
#low

#sig of uv in autumn 2021
corr <- cor.test(x=gr$Aut2021_cov, y=gr$Aut2021_PFD.UV, method = 'pearson')
corr$p.value
corr$estimate
#SIG
#only UV significant in 2021
#EVIDENCE:

corr <- cor.test(x=gr$Aut2021_cov, y=gr$Aut2021_light, method = 'pearson')
corr$p.value
corr$estimate
#low #in sig

#aut 2021
corr <- cor.test(x=gr$Aut2021_cov, y=gr$Aut2021_PFD.FR, method = 'pearson')
corr$p.value
corr$estimate
#low #in sig

corr <- cor.test(x=gr$Aut2021_cov, y=gr$Aut2021_PFD.B, method = 'pearson')
corr$p.value
corr$estimate
#low #in sig

corr <- cor.test(x=gr$Aut2021_cov, y=gr$Aut2021_PFD.R, method = 'pearson')
corr$p.value
corr$estimate
#LOW #in sig

corr <- cor.test(x=gr$Aut2021_cov, y=gr$Aut2021_PFD.G, method = 'pearson')
corr$p.value
corr$estimate
#LOW #in sig

#aut 2020 NOT WORKING DIFF NOS X AND Y
gr$Aut2020_cov #20
gr$Aut2020_light #20 with 1 NA
#THESE DONT EXIST AS DIDNT TAKE SPECTROMETER
gr$Aut2020_PFD.UV
gr$Aut2020_PFD.FR
gr$Aut2020_PFD.R
gr$Aut2020_PFD.G
gr$Aut2020_PFD.B

corr <- cor.test(x=gr$Aut2020_cov, y=gr$Aut2020_light, method = 'pearson')
corr$p.value
corr$estimate
#LOW #in sig

plot(Sum2022_cov ~ Aut2020_cov, gr) #fair corr
plot(Sum2020_cov ~ Aut2020_cov, gr)
plot(Sum2022_LambdaP ~ Sum2022_PFD.UV, gr)
plot(Sum2022_LambdaP ~ Mar2021_PFD.FR, gr)
plot(Aut2021_light ~ Aut2020_cov, gr)
plot(Aut2021_PFD.B ~ Aut2020_cov, gr)
plot(Aut2021_PFD.R ~ Aut2020_cov, gr)
plot(Aut2021_PFD.G ~ Aut2020_cov, gr)
plot(Aut2021_PPFD ~ Aut2020_cov, gr)
plot(Sum2021_cov ~ Sum2022_cov, gr)
plot(Aut2021_cov ~ Sum2022_cov, gr)
plot(Aut2021_cov ~ Aut2021_PFD.UV, gr)
plot(Sum2022_LambdaP ~ Spr2021_LambdaP, gr) #doesnt work

gr$Mar2021_LambdaP

#plot correlation matrix
gr$X <- NULL
gr$X.1 <- NULL
gr$X.2 <- NULL

grn <- gr
#corr plot
as.matrix(grn)

sapply(grn, class) #check if numeric
grnnum <- as.data.frame(lapply(grn, as.numeric)) #make numeric
str(grnnum)
matr <- grnnum

M<-cor(matr,use='complete.obs')
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
cor.mtest <- function(matr, ...) {
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
p.mat <- cor.mtest(matr)


plotcolour<- colorRampPalette(c("darkorange2", "gray88", "darkcyan"))(20) 

corrplot(M, type="upper", order="hclust",method = "number",
         p.mat = p.mat, sig.level = 0.1005, insig = "blank",tl.cex = 0.5, tl.col = 'black',col=plotcolour)

corrplot(M, col=plotcolour, type="upper", order="hclust", # Add coefficient of correlation
         tl.col="black", tl.cex = 0.5, tl.offset = 0.5,
         #Text label color and rotation
         # Combine with significance
         p.mat = p.mat,
         sig.level = c(.001, .01, .05), pch.cex = 0.6, #add signficance stars
         insig = "label_sig", pch.col = "black",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

plotcolour<- colorRampPalette(c("darkorange2", "gray88", "darkcyan"))(20)
#This is where I've created a colour scheme of my choosing. 

setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Env")
#automate linear models
#automate correlations
lms <- expand.grid(3:55, 3:55)
lms_names <- expand.grid(names(gr)[3:55], names(gr)[3:55])
out <- vector(mode = "list", length = nrow(lms))
for(i in 1:nrow(lms)){
        lms_col_2 <- lms$Var2[i]
        lms_col_1 <- lms$Var1[i]
        plot_name <- paste0(stringr::str_pad(i, width = 3, pad = "0"), " ", lms_names$Var2[i], " vs. ", lms_names$Var1[i], ".png")
        png(plot_name, width = 500, height = 500, type = "cairo")
        plot(gr[, lms_col_1], gr[, lms_col_2], main = paste0(lms_names$Var2[i], " vs. ", lms_names$Var1[i]), type = "n")
        text(gr[, lms_col_1], gr[, lms_col_2], row.names(gr), cex = 0.8)
        abline(lm(gr[, lms_col_2] ~ gr[, lms_col_1]))
        #plots per column and output as png with correspinding names
        dev.off()
        out[[i]] <- data.frame(r.squared = summary(lm(gr[, lms_col_2] ~ gr[, lms_col_1]))$r.squared,
                               pos_neg = ifelse(summary(lm(gr[, lms_col_2] ~ gr[, lms_col_1]))$coef[2, 1] > 0, "+", "-"),
                               p.value = summary(lm(gr[, lms_col_2] ~ gr[, lms_col_1]))$coefficients[2, 4])
}
#select r squared and p val from table and make df output
#use estimate col as gradient to make + or - col

#write to files
library(dplyr)
(all_output <- bind_cols(lms_names, do.call(rbind, out)))
all_output <- all_output[all_output$r.squared != 1, ]
all_output <- all_output[all_output$p.value < 0.05, ]
all_output$r.squared <- ifelse(all_output$pos_neg == "+", all_output$r.squared, -all_output$r.squared)

write.csv(all_output, "Envlight+cov_correlation.csv")

names(gr)

#nicer plots
library(ggplot2)
ggplot(gr, aes(Aut2021_cov, Aut2021_PFD.UV, label = Duckweed.collection)) +
        geom_point() +# ggplot2 plot with labels
        geom_text(aes(label = Duckweed.collection), hjust = - 0.5) +
        geom_smooth(method='lm') +
        xlab("Duckweed % coverage (Autumn)") +
        ylab("UV (umol m-2 s-1)") +
        theme_bw() +
        theme_classic()

ggplot(gr, aes(Sum2022_cov, Sum2022_PFD.UV, label = Duckweed.collection)) +    # ggplot2 plot with labels
        geom_point() +
        geom_text(aes(label = Duckweed.collection), hjust = - 0.5) +
        geom_smooth(method='lm') +
        xlab("Coverage % (Jul 2022)") +
        ylab("UV (umol m-2 s-1)") +
        theme_bw() +
        theme_classic()

ggplot(gr, aes(Sum2021_cov, Sum2021_PFD.UV, label = Duckweed.collection)) +    # ggplot2 plot with labels
        geom_point() +
        geom_text(aes(label = Duckweed.collection), hjust = - 0.5) +
        geom_smooth(method='lm') +
        xlab("Coverage % (Jul 2021)") +
        ylab("UV (umol m-2 s-1)") +
        theme_bw() +
        theme_classic()



#pointless
ggplot(gr, aes(Sum2022_cov, Aut2020_cov, label = Duckweed.collection)) +    # ggplot2 plot with labels
        geom_point() +
        geom_text(aes(label = Duckweed.collection), hjust = - 0.5) +
        geom_smooth(method='lm') +
        xlab("% Coverage (Sum 2022)") +
        ylab("% Coverage (Aut 2020)") +
        theme_bw() +
        theme_classic()

#to read in with r:fr to make boxplots and t tests for env light
gr <- read.csv("Env-light+coverage+RFRrat+UVprop.csv")

#new data included per site, new rfr calcs, % calcs
gr <- read.csv("Env-light+coverage+moredata.csv")

#boxplots for ll and hl organised sites
gr$Light <- c("LL","LL", "LL", "HL", "LL", "LL",
              "HL", "HL", "HL", "LL", "HL", "LL", "LL",
              "HL", "HL", "LL", "LL", "LL", "LL", "HL")

#no ks01
gr$Light <- c("LL", "LL", "HL", "LL", "LL",
              "HL", "HL", "HL", "LL", "HL", "LL", "LL",
              "HL", "HL", "LL", "LL", "LL", "LL", "HL")

as.factor(gr$Light)

#difference in coverage with light grouping
par(mfrow = c(2, 4))
boxplot(Spr2020_cov ~ Light, gr)
boxplot(Sum2020_cov ~ Light, gr)
boxplot(Aut2020_cov ~ Light, gr)
boxplot(Spr2021_cov ~ Light, gr)
boxplot(Sum2021_cov ~ Light, gr)
boxplot(Aut2021_cov ~ Light, gr)
boxplot(Spr2022_cov ~ Light, gr)
boxplot(Sum2022_cov ~ Light, gr)
#comparable spring, LL improving in 2022. Summer generally higher in HL,
#2022 more comparable
#Aut 2020 coverage higher in LL, Aut 2021 HL higher than LL.

#see if light factors sig diff when group by HL or LL sites

#SUMMER 2022
par(mfrow = c(2, 4))
boxplot(Sum2022_PPFD ~ Light, gr)
boxplot(Sum2022_PFD.FR ~ Light, gr)
boxplot(Sum2022_PFD.R ~ Light, gr)
boxplot(Sum2022_RFR ~ Light, gr)
boxplot(Sum2022_propUV ~ Light, gr)
boxplot(Sum2022_PFD.G ~ Light, gr)
boxplot(Sum2022_PFD.B ~ Light, gr)
boxplot(Sum2022_PFD.UV ~ Light, gr)
boxplot(Sum2022_LambdaD ~ Light, gr)
boxplot(Sum2022_LambdaP ~ Light, gr)

#SPRING 2022
par(mfrow = c(2, 3))
boxplot(Spr2022_PPFD ~ Light, gr)
boxplot(Spr2022_PFD.FR ~ Light, gr)
boxplot(Spr2022_PFD.R ~ Light, gr)
boxplot(Spr2022_RFR ~ Light, gr)
boxplot(Spr2022_propUV ~ Light, gr)
boxplot(Spr2022_PFD.G ~ Light, gr)
boxplot(Spr2022_PFD.B ~ Light, gr)
boxplot(Spr2022_PFD.UV ~ Light, gr)
boxplot(Spr2022_LambdaP ~ Light, gr)
boxplot(Spr2022_LambdaD ~ Light, gr)

names(gr)
#AVERAGES PER LIGHT LEVEL
t.test(Sum2022_PPFD ~ Light, data=gr) #**
t.test(Sum2022_PFD ~ Light, data=gr) #**
t.test(Sum2022_PFD.FR ~ Light, data=gr) #**
t.test(Sum2022_PFD.R ~ Light, data=gr) #**
t.test(Sum2022_RFR ~ Light, data=gr) #***
t.test(Sum2022_PFD.G ~ Light, data=gr) #**
t.test(Sum2022_PFD.B ~ Light, data=gr) #**
t.test(Sum2022_PFD.UV ~ Light, data=gr) #**
t.test(Sum2022_LambdaP ~ Light, data=gr) #**
t.test(Sum2022_LambdaD ~ Light, data=gr) #*

t.test(Sum2022_propFR ~ Light, data=gr) #***
t.test(Sum2022_propR ~ Light, data=gr) #***
t.test(Sum2022_propG ~ Light, data=gr) #***
t.test(Sum2022_propB ~ Light, data=gr) #NS
t.test(Sum2022_propUV ~ Light, data=gr) #NS
t.test(Sum2022_propPPFD ~ Light, data=gr) #**

#AVERAGES PER LIGHT LEVEL
t.test(Sum2021_PPFD ~ Light, data=gr) #**
t.test(Sum2021_PFD ~ Light, data=gr) #**
t.test(Sum2021_PFD.FR ~ Light, data=gr) #**
t.test(Sum2021_PFD.R ~ Light, data=gr) #**
t.test(Sum2021_RFR ~ Light, data=gr) #***
t.test(Sum2021_PFD.G ~ Light, data=gr)#**
t.test(Sum2021_PFD.B ~ Light, data=gr) #**
t.test(Sum2021_PFD.UV ~ Light, data=gr) #**
t.test(Sum2021_LambdaP ~ Light, data=gr) #***
#all sig
t.test(Sum2021_LambdaD ~ Light, data=gr) #n.s

t.test(Sum2021_propFR ~ Light, data=gr) #***
t.test(Sum2021_propR ~ Light, data=gr) #***
t.test(Sum2021_propG ~ Light, data=gr) #***
t.test(Sum2021_propB ~ Light, data=gr) #NS
t.test(Sum2021_propUV ~ Light, data=gr) #**
t.test(Sum2021_propPPFD ~ Light, data=gr) #***

#AVERAGES PER LIGHT LEVEL
t.test(Spr2022_PPFD ~ Light, data=gr) #ns
t.test(Spr2022_PFD ~ Light, data=gr) #ns
t.test(Spr2022_PFD.FR ~ Light, data=gr) #ns
t.test(Spr2022_PFD.R ~ Light, data=gr) #ns
t.test(Spr2022_RFR ~ Light, data=gr) #*
t.test(Spr2022_PFD.G ~ Light, data=gr) #ns
t.test(Spr2022_PFD.B ~ Light, data=gr) #ns
t.test(Spr2022_PFD.UV ~ Light, data=gr) #ns
t.test(Spr2022_LambdaP ~ Light, data=gr) #ns

t.test(Spr2022_propFR ~ Light, data=gr) #ns
t.test(Spr2022_propR ~ Light, data=gr) #***
t.test(Spr2022_propG ~ Light, data=gr) #ns
t.test(Spr2022_propB ~ Light, data=gr) #*
t.test(Spr2022_propUV ~ Light, data=gr) #*
t.test(Spr2022_propPPFD ~ Light, data=gr) #ns

#AVERAGES PER LIGHT LEVEL
t.test(Spr2021_PPFD ~ Light, data=gr) #**
t.test(Spr2021_PFD ~ Light, data=gr) #**
t.test(Spr2021_PFD.FR ~ Light, data=gr) #**
t.test(Spr2021_PFD.R ~ Light, data=gr) #**
t.test(Spr2021_RFR ~ Light, data=gr) #*
t.test(Spr2021_PFD.G ~ Light, data=gr) #**
t.test(Spr2021_PFD.B ~ Light, data=gr) #**
t.test(Spr2021_PFD.UV ~ Light, data=gr) #**
t.test(Spr2021_LambdaP ~ Light, data=gr) #***
t.test(Spr2021_LambdaD ~ Light, data=gr) #*

t.test(Spr2021_propFR ~ Light, data=gr) #ns
t.test(Spr2021_propR ~ Light, data=gr) #***
t.test(Spr2021_propG ~ Light, data=gr) #*
t.test(Spr2021_propB ~ Light, data=gr) #**
t.test(Spr2021_propUV ~ Light, data=gr) #***
t.test(Spr2021_propPPFD ~ Light, data=gr) #ns

#AVERAGES PER LIGHT LEVEL
t.test(Aut2021_PPFD ~ Light, data=gr) #**
t.test(Aut2021_PFD ~ Light, data=gr) #**
t.test(Aut2021_PFD.FR ~ Light, data=gr) #**
t.test(Aut2021_PFD.R ~ Light, data=gr) #**
t.test(Aut2021_RFR ~ Light, data=gr) #*
t.test(Aut2021_PFD.G ~ Light, data=gr) #**
t.test(Aut2021_PFD.B ~ Light, data=gr) #**
t.test(Aut2021_PFD.UV ~ Light, data=gr) #*
t.test(Aut2021_LambdaP ~ Light, data=gr) #ns
t.test(Aut2021_LambdaD ~ Light, data=gr)
#most n.s

t.test(Aut2021_propFR ~ Light, data=gr) #*
t.test(Aut2021_propR ~ Light, data=gr) #ns
t.test(Aut2021_propG ~ Light, data=gr) #**
t.test(Aut2021_propB ~ Light, data=gr) #ns
t.test(Aut2021_propUV ~ Light, data=gr) #ns
t.test(Aut2021_propPPFD ~ Light, data=gr) #**

library(dplyr)
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PPFD_mean = mean(Sum2022_PPFD, na.rm = TRUE), Sum2022_PPFD_stdev = sd(Sum2022_PPFD, na.rm = TRUE), Sum2022_PPFD_maximum = max(Sum2022_PPFD, na.rm = TRUE))
Summary
Summaryx <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PFD_mean = mean(Sum2022_PFD, na.rm = TRUE), Sum2022_PFD_stdev = sd(Sum2022_PFD, na.rm = TRUE), Sum2022_PFD_maximum = max(Sum2022_PFD, na.rm = TRUE))
Summaryx
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PFD.FR_mean = mean(Sum2022_PFD.FR, na.rm = TRUE), Sum2022_PFD.FR_stdev = sd(Sum2022_PFD.FR, na.rm = TRUE), Sum2022_PFD.FR_maximum = max(Sum2022_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PFD.R_mean = mean(Sum2022_PFD.R, na.rm = TRUE), Sum2022_PFD.R_stdev = sd(Sum2022_PFD.R, na.rm = TRUE), Sum2022_PFD.R_maximum = max(Sum2022_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PFD.G_mean = mean(Sum2022_PFD.G, na.rm = TRUE), Sum2022_PFD.G_stdev = sd(Sum2022_PFD.G, na.rm = TRUE), Sum2022_PFD.G_maximum = max(Sum2022_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PFD.B_mean = mean(Sum2022_PFD.B, na.rm = TRUE), Sum2022_PFD.B_stdev = sd(Sum2022_PFD.B, na.rm = TRUE), Sum2022_PFD.B_maximum = max(Sum2022_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_PFD.UV_mean = mean(Sum2022_PFD.UV, na.rm = TRUE), Sum2022_PFD.UV_stdev = sd(Sum2022_PFD.UV, na.rm = TRUE), Sum2022_PFD.UV_maximum = max(Sum2022_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_LambdaP_mean = mean(Sum2022_LambdaP, na.rm = TRUE), Sum2022_LambdaP_stdev = sd(Sum2022_LambdaP, na.rm = TRUE), Sum2022_LambdaP_maximum = max(Sum2022_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2022_LambdaD_mean = mean(Sum2022_LambdaD, na.rm = TRUE), Sum2022_LambdaD_stdev = sd(Sum2022_LambdaD, na.rm = TRUE), Sum2022_LambdaD_maximum = max(Sum2022_LambdaD, na.rm = TRUE))
Summary7

Summ <- cbind(Summary, Summaryx, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)

#spring 2022
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PPFD_mean = mean(Spr2022_PPFD, na.rm = TRUE), Spr2022_PPFD_stdev = sd(Spr2022_PPFD, na.rm = TRUE), Spr2022_PPFD_maximum = max(Spr2022_PPFD, na.rm = TRUE))
Summary
Summaryx <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD_mean = mean(Spr2022_PFD, na.rm = TRUE), Spr2022_PFD_stdev = sd(Spr2022_PFD, na.rm = TRUE), Spr2022_PFD_maximum = max(Spr2022_PFD, na.rm = TRUE))
Summaryx
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.FR_mean = mean(Spr2022_PFD.FR, na.rm = TRUE), Spr2022_PFD.FR_stdev = sd(Spr2022_PFD.FR, na.rm = TRUE), Spr2022_PFD.FR_maximum = max(Spr2022_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.R_mean = mean(Spr2022_PFD.R, na.rm = TRUE), Spr2022_PFD.R_stdev = sd(Spr2022_PFD.R, na.rm = TRUE), Spr2022_PFD.R_maximum = max(Spr2022_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.G_mean = mean(Spr2022_PFD.G, na.rm = TRUE), Spr2022_PFD.G_stdev = sd(Spr2022_PFD.G, na.rm = TRUE), Spr2022_PFD.G_maximum = max(Spr2022_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.B_mean = mean(Spr2022_PFD.B, na.rm = TRUE), Spr2022_PFD.B_stdev = sd(Spr2022_PFD.B, na.rm = TRUE), Spr2022_PFD.B_maximum = max(Spr2022_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.UV_mean = mean(Spr2022_PFD.UV, na.rm = TRUE), Spr2022_PFD.UV_stdev = sd(Spr2022_PFD.UV, na.rm = TRUE), Spr2022_PFD.UV_maximum = max(Spr2022_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_LambdaP_mean = mean(Spr2022_LambdaP, na.rm = TRUE), Spr2022_LambdaP_stdev = sd(Spr2022_LambdaP, na.rm = TRUE), Spr2022_LambdaP_maximum = max(Spr2022_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_LambdaD_mean = mean(Spr2022_LambdaD, na.rm = TRUE), Spr2022_LambdaD_stdev = sd(Spr2022_LambdaD, na.rm = TRUE), Spr2022_LambdaD_maximum = max(Spr2022_LambdaD, na.rm = TRUE))
Summary7

Summ <- cbind(Summary, Summaryx, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)

#aut 2021
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PPFD_mean = mean(Aut2021_PPFD, na.rm = TRUE), Aut2021_PPFD_stdev = sd(Aut2021_PPFD, na.rm = TRUE), Aut2021_PPFD_maximum = max(Aut2021_PPFD, na.rm = TRUE))
Summary
Summaryx <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD_mean = mean(Aut2021_PFD, na.rm = TRUE), Aut2021_PFD_stdev = sd(Aut2021_PFD, na.rm = TRUE), Aut2021_PFD_maximum = max(Aut2021_PFD, na.rm = TRUE))
Summaryx
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.FR_mean = mean(Aut2021_PFD.FR, na.rm = TRUE), Aut2021_PFD.FR_stdev = sd(Aut2021_PFD.FR, na.rm = TRUE), Aut2021_PFD.FR_maximum = max(Aut2021_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.R_mean = mean(Aut2021_PFD.R, na.rm = TRUE), Aut2021_PFD.R_stdev = sd(Aut2021_PFD.R, na.rm = TRUE), Aut2021_PFD.R_maximum = max(Aut2021_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.G_mean = mean(Aut2021_PFD.G, na.rm = TRUE), Aut2021_PFD.G_stdev = sd(Aut2021_PFD.G, na.rm = TRUE), Aut2021_PFD.G_maximum = max(Aut2021_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.B_mean = mean(Aut2021_PFD.B, na.rm = TRUE), Aut2021_PFD.B_stdev = sd(Aut2021_PFD.B, na.rm = TRUE), Aut2021_PFD.B_maximum = max(Aut2021_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.UV_mean = mean(Aut2021_PFD.UV, na.rm = TRUE), Aut2021_PFD.UV_stdev = sd(Aut2021_PFD.UV, na.rm = TRUE), Aut2021_PFD.UV_maximum = max(Aut2021_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_LambdaP_mean = mean(Aut2021_LambdaP, na.rm = TRUE), Aut2021_LambdaP_stdev = sd(Aut2021_LambdaP, na.rm = TRUE), Aut2021_LambdaP_maximum = max(Aut2021_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_LambdaD_mean = mean(Aut2021_LambdaD, na.rm = TRUE), Aut2021_LambdaD_stdev = sd(Aut2021_LambdaD, na.rm = TRUE), Aut2021_LambdaD_maximum = max(Aut2021_LambdaD, na.rm = TRUE))
Summary7

Summ <- cbind(Summary, Summaryx, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)

#spr 2021
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PPFD_mean = mean(Spr2021_PPFD, na.rm = TRUE), Spr2021_PPFD_stdev = sd(Spr2021_PPFD, na.rm = TRUE), Spr2021_PPFD_maximum = max(Spr2021_PPFD, na.rm = TRUE))
Summary
Summaryx <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD_mean = mean(Spr2021_PFD, na.rm = TRUE), Spr2021_PFD_stdev = sd(Spr2021_PFD, na.rm = TRUE), Spr2021_PFD_maximum = max(Spr2021_PFD, na.rm = TRUE))
Summaryx
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.FR_mean = mean(Spr2021_PFD.FR, na.rm = TRUE), Spr2021_PFD.FR_stdev = sd(Spr2021_PFD.FR, na.rm = TRUE), Spr2021_PFD.FR_maximum = max(Spr2021_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.R_mean = mean(Spr2021_PFD.R, na.rm = TRUE), Spr2021_PFD.R_stdev = sd(Spr2021_PFD.R, na.rm = TRUE), Spr2021_PFD.R_maximum = max(Spr2021_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.G_mean = mean(Spr2021_PFD.G, na.rm = TRUE), Spr2021_PFD.G_stdev = sd(Spr2021_PFD.G, na.rm = TRUE), Spr2021_PFD.G_maximum = max(Spr2021_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.B_mean = mean(Spr2021_PFD.B, na.rm = TRUE), Spr2021_PFD.B_stdev = sd(Spr2021_PFD.B, na.rm = TRUE), Spr2021_PFD.B_maximum = max(Spr2021_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.UV_mean = mean(Spr2021_PFD.UV, na.rm = TRUE), Spr2021_PFD.UV_stdev = sd(Spr2021_PFD.UV, na.rm = TRUE), Spr2021_PFD.UV_maximum = max(Spr2021_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_LambdaP_mean = mean(Spr2021_LambdaP, na.rm = TRUE), Spr2021_LambdaP_stdev = sd(Spr2021_LambdaP, na.rm = TRUE), Spr2021_LambdaP_maximum = max(Spr2021_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_LambdaD_mean = mean(Spr2021_LambdaD, na.rm = TRUE), Spr2021_LambdaD_stdev = sd(Spr2021_LambdaD, na.rm = TRUE), Aut2021_LambdaD_maximum = max(Spr2021_LambdaD, na.rm = TRUE))
Summary7

Summ <- cbind(Summary, Summaryx, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)

#sum 2021
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PPFD_mean = mean(Sum2021_PPFD, na.rm = TRUE), Sum2021_PPFD_stdev = sd(Sum2021_PPFD, na.rm = TRUE), Sum2021_PPFD_maximum = max(Sum2021_PPFD, na.rm = TRUE))
Summary
Summaryx <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD_mean = mean(Sum2021_PFD, na.rm = TRUE), Sum2021_PFD_stdev = sd(Sum2021_PFD, na.rm = TRUE), Sum2021_PFD_maximum = max(Sum2021_PFD, na.rm = TRUE))
Summaryx
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.FR_mean = mean(Sum2021_PFD.FR, na.rm = TRUE), Sum2021_PFD.FR_stdev = sd(Sum2021_PFD.FR, na.rm = TRUE), Sum2021_PFD.FR_maximum = max(Sum2021_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.R_mean = mean(Sum2021_PFD.R, na.rm = TRUE), Sum2021_PFD.R_stdev = sd(Sum2021_PFD.R, na.rm = TRUE), Sum2021_PFD.R_maximum = max(Sum2021_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.G_mean = mean(Sum2021_PFD.G, na.rm = TRUE), Sum2021_PFD.G_stdev = sd(Sum2021_PFD.G, na.rm = TRUE), Sum2021_PFD.G_maximum = max(Sum2021_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.B_mean = mean(Sum2021_PFD.B, na.rm = TRUE), Sum2021_PFD.B_stdev = sd(Sum2021_PFD.B, na.rm = TRUE), Sum2021_PFD.B_maximum = max(Sum2021_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.UV_mean = mean(Sum2021_PFD.UV, na.rm = TRUE), Sum2021_PFD.UV_stdev = sd(Sum2021_PFD.UV, na.rm = TRUE), Sum2021_PFD.UV_maximum = max(Sum2021_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_LambdaP_mean = mean(Sum2021_LambdaP, na.rm = TRUE), Sum2021_LambdaP_stdev = sd(Sum2021_LambdaP, na.rm = TRUE), Sum2021_LambdaP_maximum = max(Sum2021_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_LambdaD_mean = mean(Sum2021_LambdaD, na.rm = TRUE), Sum2021_LambdaD_stdev = sd(Sum2021_LambdaD, na.rm = TRUE), Aut2021_LambdaD_maximum = max(Sum2021_LambdaD, na.rm = TRUE))
Summary7

Summ <- cbind(Summary, Summaryx, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)

#AVERAGES PER LIGHT LEVEL
t.test(Spr2022_PPFD ~ Light, data=gr)
t.test(Spr2022_PFD.FR ~ Light, data=gr)
t.test(Spr2022_PFD.R ~ Light, data=gr)
t.test(Spr2022_PFD.G ~ Light, data=gr)
t.test(Spr2022_PFD.B ~ Light, data=gr)
t.test(Spr2022_PFD.UV ~ Light, data=gr)
#all above sig
t.test(Spr2022_LambdaP ~ Light, data=gr) #n.s
t.test(Spr2022_LambdaD ~ Light, data=gr) #doesnt exist

aov <- aov(Sum2022_PPFD ~ Light, data=gr)
summary(aov)
tukey <- TukeyHSD(aov, "Light", conf.level=.95) 
str(tukey)
tukey

#see if differences are year round
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PPFD_mean = mean(Spr2022_PPFD, na.rm = TRUE), Spr2022_PPFD_stdev = sd(Spr2022_PPFD, na.rm = TRUE), Spr2022_PPFD_maximum = max(Spr2022_PPFD, na.rm = TRUE))
Summary
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.FR_mean = mean(Spr2022_PFD.FR, na.rm = TRUE), Spr2022_PFD.FR_stdev = sd(Spr2022_PFD.FR, na.rm = TRUE), Spr2022_PFD.FR_maximum = max(Spr2022_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.R_mean = mean(Spr2022_PFD.R, na.rm = TRUE), Spr2022_PFD.R_stdev = sd(Spr2022_PFD.R, na.rm = TRUE), Spr2022_PFD.R_maximum = max(Spr2022_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.G_mean = mean(Spr2022_PFD.G, na.rm = TRUE), Spr2022_PFD.G_stdev = sd(Spr2022_PFD.G, na.rm = TRUE), Spr2022_PFD.G_maximum = max(Spr2022_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.B_mean = mean(Spr2022_PFD.B, na.rm = TRUE), Spr2022_PFD.B_stdev = sd(Spr2022_PFD.B, na.rm = TRUE), Spr2022_PFD.B_maximum = max(Spr2022_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_PFD.UV_mean = mean(Spr2022_PFD.UV, na.rm = TRUE), Spr2022_PFD.UV_stdev = sd(Spr2022_PFD.UV, na.rm = TRUE), Spr2022_PFD.UV_maximum = max(Spr2022_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_LambdaP_mean = mean(Spr2022_LambdaP, na.rm = TRUE), Spr2022_LambdaP_stdev = sd(Spr2022_LambdaP, na.rm = TRUE), Spr2022_LambdaP_maximum = max(Spr2022_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2022_LambdaD_mean = mean(Spr2022_LambdaD, na.rm = TRUE), Spr2022_LambdaD_stdev = sd(Spr2022_LambdaD, na.rm = TRUE), Spr2022_LambdaD_maximum = max(Spr2022_LambdaD, na.rm = TRUE))
Summary7

Summ_Spr2022 <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)


#see if differences are year round
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PPFD_mean = mean(Aut2021_PPFD, na.rm = TRUE), Aut2021_PPFD_stdev = sd(Aut2021_PPFD, na.rm = TRUE), Aut2021_PPFD_maximum = max(Aut2021_PPFD, na.rm = TRUE))
Summary
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.FR_mean = mean(Aut2021_PFD.FR, na.rm = TRUE), Aut2021_PFD.FR_stdev = sd(Aut2021_PFD.FR, na.rm = TRUE), Aut2021_PFD.FR_maximum = max(Aut2021_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.R_mean = mean(Aut2021_PFD.R, na.rm = TRUE), Aut2021_PFD.R_stdev = sd(Aut2021_PFD.R, na.rm = TRUE), Aut2021_PFD.R_maximum = max(Aut2021_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.G_mean = mean(Aut2021_PFD.G, na.rm = TRUE), Aut2021_PFD.G_stdev = sd(Aut2021_PFD.G, na.rm = TRUE), Aut2021_PFD.G_maximum = max(Aut2021_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.B_mean = mean(Aut2021_PFD.B, na.rm = TRUE), Aut2021_PFD.B_stdev = sd(Aut2021_PFD.B, na.rm = TRUE), Aut2021_PFD.B_maximum = max(Aut2021_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_PFD.UV_mean = mean(Aut2021_PFD.UV, na.rm = TRUE), Aut2021_PFD.UV_stdev = sd(Aut2021_PFD.UV, na.rm = TRUE), Aut2021_PFD.UV_maximum = max(Aut2021_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_LambdaP_mean = mean(Aut2021_LambdaP, na.rm = TRUE), Aut2021_LambdaP_stdev = sd(Aut2021_LambdaP, na.rm = TRUE), Aut2021_LambdaP_maximum = max(Aut2021_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Aut2021_LambdaD_mean = mean(Aut2021_LambdaD, na.rm = TRUE), Aut2021_LambdaD_stdev = sd(Aut2021_LambdaD, na.rm = TRUE), Aut2021_LambdaD_maximum = max(Aut2021_LambdaD, na.rm = TRUE))
Summary7

Summ_Aut2021 <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5,
                      Summary6,Summary7)

#sum 2021
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PPFD_mean = mean(Sum2021_PPFD, na.rm = TRUE), Sum2021_PPFD_stdev = sd(Sum2021_PPFD, na.rm = TRUE), Sum2021_PPFD_maximum = max(Sum2021_PPFD, na.rm = TRUE))
Summary
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.FR_mean = mean(Sum2021_PFD.FR, na.rm = TRUE), Sum2021_PFD.FR_stdev = sd(Sum2021_PFD.FR, na.rm = TRUE), Sum2021_PFD.FR_maximum = max(Sum2021_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.R_mean = mean(Sum2021_PFD.R, na.rm = TRUE), Sum2021_PFD.R_stdev = sd(Sum2021_PFD.R, na.rm = TRUE), Sum2021_PFD.R_maximum = max(Sum2021_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.G_mean = mean(Sum2021_PFD.G, na.rm = TRUE), Sum2021_PFD.G_stdev = sd(Sum2021_PFD.G, na.rm = TRUE), Sum2021_PFD.G_maximum = max(Sum2021_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.B_mean = mean(Sum2021_PFD.B, na.rm = TRUE), Sum2021_PFD.B_stdev = sd(Sum2021_PFD.B, na.rm = TRUE), Sum2021_PFD.B_maximum = max(Sum2021_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_PFD.UV_mean = mean(Sum2021_PFD.UV, na.rm = TRUE), Sum2021_PFD.UV_stdev = sd(Sum2021_PFD.UV, na.rm = TRUE), Sum2021_PFD.UV_maximum = max(Sum2021_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_LambdaP_mean = mean(Sum2021_LambdaP, na.rm = TRUE), Sum2021_LambdaP_stdev = sd(Sum2021_LambdaP, na.rm = TRUE), Sum2021_LambdaP_maximum = max(Sum2021_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Sum2021_LambdaD_mean = mean(Sum2021_LambdaD, na.rm = TRUE), Sum2021_LambdaD_stdev = sd(Sum2021_LambdaD, na.rm = TRUE), Sum2021_LambdaD_maximum = max(Sum2021_LambdaD, na.rm = TRUE))
Summary7

Summ_Sum2021 <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5,
              Summary6,Summary7)

#spr 2021
Summary <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PPFD_mean = mean(Spr2021_PPFD, na.rm = TRUE), Spr2021_PPFD_stdev = sd(Spr2021_PPFD, na.rm = TRUE), Spr2021_PPFD_maximum = max(Spr2021_PPFD, na.rm = TRUE))
Summary
Summary1 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.FR_mean = mean(Spr2021_PFD.FR, na.rm = TRUE), Spr2021_PFD.FR_stdev = sd(Spr2021_PFD.FR, na.rm = TRUE), Spr2021_PFD.FR_maximum = max(Spr2021_PFD.FR, na.rm = TRUE))
Summary1
Summary2 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.R_mean = mean(Spr2021_PFD.R, na.rm = TRUE), Spr2021_PFD.R_stdev = sd(Spr2021_PFD.R, na.rm = TRUE), Spr2021_PFD.R_maximum = max(Spr2021_PFD.R, na.rm = TRUE))
Summary2
Summary3 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.G_mean = mean(Spr2021_PFD.G, na.rm = TRUE), Spr2021_PFD.G_stdev = sd(Spr2021_PFD.G, na.rm = TRUE), Spr2021_PFD.G_maximum = max(Spr2021_PFD.G, na.rm = TRUE))
Summary3
Summary4 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.B_mean = mean(Spr2021_PFD.B, na.rm = TRUE), Spr2021_PFD.B_stdev = sd(Spr2021_PFD.B, na.rm = TRUE), Spr2021_PFD.B_maximum = max(Spr2021_PFD.B, na.rm = TRUE))
Summary4
Summary5 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_PFD.UV_mean = mean(Spr2021_PFD.UV, na.rm = TRUE), Spr2021_PFD.UV_stdev = sd(Spr2021_PFD.UV, na.rm = TRUE), Spr2021_PFD.UV_maximum = max(Spr2021_PFD.UV, na.rm = TRUE))
Summary5
Summary6 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_LambdaP_mean = mean(Spr2021_LambdaP, na.rm = TRUE), Spr2021_LambdaP_stdev = sd(Spr2021_LambdaP, na.rm = TRUE), Spr2021_LambdaP_maximum = max(Spr2021_LambdaP, na.rm = TRUE))
Summary6
Summary7 <- gr %>%
        group_by(Light) %>%
        summarise(Spr2021_LambdaD_mean = mean(Spr2021_LambdaD, na.rm = TRUE), Spr2021_LambdaD_stdev = sd(Spr2021_LambdaD, na.rm = TRUE), Spr2021_LambdaD_maximum = max(Spr2021_LambdaD, na.rm = TRUE))
Summary7

Summ_Spr2021 <- cbind(Summary, Summary1, Summary2, Summary3, Summary4, Summary5,
                      Summary6,Summary7)
