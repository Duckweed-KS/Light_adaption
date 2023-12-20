library(vcfR)
library(ggplot2)
library(ggrepel)

setwd("C:\\Users\\kelli\\OneDrive\\Documents\\Duckweed\\adgenet\\2023")

#original
#vcf <- read.vcfR("Cochlearia_4dg_106.purged_pruned.ann.vcf") #Filtered and LD-pruned VCF file

#ks vcf with old 4fds
#variants = 43,0000
#vcf <- read.vcfR("uklight.vcf.gz")

#ks vcf uklight with old 4fds + filtering prune
#variants = 6,219
vcf <- read.vcfR("4fold_mis05_maf001_ld200_pruned.vcf")

#uk light with diff params 8025 variants
#mis too low according to levi
vcf <- read.vcfR("4fold_ld_pruned_mis0.75_maf0.025_50_50_0.2.vcf")
#no mis included, just do in pca script
#7288 variants
vcf <- read.vcfR("4fold_ld_pruned_maf0.05_50_50_0.1.vcf")

#ks vcf uklight with old 4fds + filtering prune for flav data
#variants = 5,947
vcf <- read.vcfR("4fold_flav_mis05_maf001_200_ld_pruned.vcf")

#diff for flav 5.268 no maf param yet
vcf <- read.vcfR("4fold_flav_ld_pruned_maf0.05_100_50_0.1.vcf")

#Transform VCF to numeric genotypes
df <- extract.gt(vcf)
df[df == "0|0"] <- 0
df[df == "0|1"] <- 1
df[df == "1|0"] <- 1
df[df == "1|1"] <- 2
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/0"] <- 1
df[df == "1/1"] <- 2
df[df == "0/0/0/0"] <- 0
df[df == "0/0/0/1"] <- 1
df[df == "0/0/1/1"] <- 2
df[df == "0/1/1/1"] <- 3
df[df == "1/1/1/1"] <- 4
df[df == "0/0/0/0/0/0"] <- 0
df[df == "0/0/0/0/0/1"] <- 1
df[df == "0/0/0/0/1/1"] <- 2
df[df == "0/0/0/1/1/1"] <- 3
df[df == "0/0/1/1/1/1"] <- 4
df[df == "0/1/1/1/1/1"] <- 5
df[df == "1/1/1/1/1/1"] <- 6
df <- data.frame(apply(df,2,function(x)as.numeric(as.character(x))))

#Remove samples with > 50% missing data
#mis <- apply(df,2,function(x)sum(is.na(x))/length(x))
#df <- df[,mis <= 0.5]

#Remove samples with > 80% missing data
mis <- apply(df,2,function(x)sum(is.na(x))/length(x))
df <- df[,mis <= 0.8]

#Calculate allele frequencies
ploidy <- apply(df,2,max,na.rm=T)
p <- apply(df,1,function(x)sum(x,na.rm=T)/sum(ploidy[!is.na(x)]))

#Removing individuals can change allele frequencies, so we make sure that maf > 0.05
#df <- df[p >= 0.05 & p <= 0.95,]
#p <- p[p >= 0.05 & p <= 0.95]

#df colnames rename
#colnames(df) #gives read out of names
#change names
#colnames(df) <- c('S. intermedia 8410',  'L. minuta 7760', 'L. minor 9441',
#                  'L. minor 7868', 'S. polyrhiza 9504', 'L. turionifera 9434',
#                  'L. minuta 9581', 'L. minuta 9484', 'L. minuta 7612', 'L. minuta 6717',
#                  'KS09',  'KS12',  'KS15', 'KS21', 'KS02', 'KS16',
#                  'KS17', 'KS18', 'KS20', 'KS22', 'KS25', 'KS27',
#                  'KS28', 'KS03', 'KS29', 'LY01A', 'LY01B', 'LY02', 
#                  'LY03', 'KS04', 'KS06A', 'KS06B', 'KS13', 'KS14')

#Estimate a covariance matrix
n <- ncol(df)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
	for(j in 1:i){
		x <- mean(c(ploidy[i],ploidy[j]))
		cov[i,j] <- mean((df[,i]-x*p)*(df[,j]-x*p)/(x*p*(1-p)),na.rm=T)
		cov[j,i] <- cov[i,j]
	}	
}

#Do PCA on the matrix and plot

#missing out cov for accessions with NAs in?

#work ploidy colours in sometime: c("#1e90ff", "#ffa500", "#7cfc00"))

pc <- prcomp(cov,scale=T)
xlab <- paste0("PC1 (",round(summary(pc)$importance[2]*100),"%)")
ylab <- paste0("PC2 (",round(summary(pc)$importance[5]*100),"%)")
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],id=colnames(df),ploidy=ploidy)
colors <- pc <- prcomp(cov,scale=T)
xlab <- paste0("PC1 (",round(summary(pc)$importance[2]*100),"%)")
ylab <- paste0("PC2 (",round(summary(pc)$importance[5]*100),"%)")
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],id=colnames(df),ploidy=ploidy)
colors <- c("#1e90ff", "#ffa500", "#7cfc00")

#not changing names on graph
#cant not do 0.8 as wont take for pca
#for uklight 0.8 missing stays the same doesnt matter which input prune used?
pcs$id #27 acccessions left
pcs$id <- c('L. minuta 7760', 'L. minor 9441', 'L. minor 7868',
                  'KS09',  'KS12',  'KS15', 'KS21', 'KS02', 'KS16',
                  'KS17', 'KS18', 'KS20', 'KS22', 'KS25', 'KS27',
                  'KS28', 'KS03', 'KS29', 'LY01A', 'LY01B', 'LY02', 
                  'LY03', 'KS04', 'KS06A', 'KS06B', 'KS13', 'KS14')

pcs$Species <- c('L. minuta', 'L. minor', 'L. minor',
                 'L. minor',  'S. polyrhiza',  'L. minor', 'L. minor', 'L. minor', 'L. turionifera',
                 'L. minor', 'L. minor', 'L. minuta', 'L. turionifera', 'L. minuta', 'L. minor',
                 'L. minor', 'L. minor', 'L. minor', 'L. minor', 'L. minuta', 'L. minor', 
                 'L. minor', 'L. minor', 'L. minuta', 'L. minuta', 'L. minor', 'L. minor')

#mod for flavour data
#not mattered if included maf in pruning, maf 0.8 here same output no species
pcs$id #31 acccessions left
pcs$id <- c('KS66A', 'KS77A', 'KS78A', 'KSNUFF2', 
            'L. minuta 7760', 'L. minor 9441', 'L. minor 7868',
            'KS09',  'KS12',  'KS15', 'KS21', 'KS02', 'KS16',
            'KS17', 'KS18', 'KS20', 'KS22', 'KS25', 'KS27',
            'KS28', 'KS03', 'KS29', 'LY01A', 'LY01B', 'LY02', 
            'LY03', 'KS04', 'KS06A', 'KS06B', 'KS13', 'KS14')

pcs$Species <- c('L. minor', 'S. polyrhiza', 'S. polyrhiza', 'L. minor', 'L. minuta', 'L. minor', 'L. minor',
                 'L. minor',  'S. polyrhiza',  'L. minor', 'L. minor', 'L. minor', 'L. turionifera',
                 'L. minor', 'L. minor', 'L. minuta', 'L. turionifera', 'L. minuta', 'L. minor',
                 'L. minor', 'L. minor', 'L. minor', 'L. minor', 'L. minuta', 'L. minor', 
                 'L. minor', 'L. minor', 'L. minuta', 'L. minuta', 'L. minor', 'L. minor')


#need to manually change as some accessions missing? 20/34
ggplot(pcs, aes(PC1, PC2, color=as.factor(Species)))+
  geom_point(size=4)+
  labs(x=xlab, y=ylab, color="Species")+
  geom_text_repel(aes(label=id), size=4, force=20, max.overlaps=Inf, color="black")+
  theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		plot.title=element_text(size=14, color="black", hjust = 0.5),
		legend.text=element_text(size=11, color="black"),
		legend.title=element_text(size=12, color="black"),
		legend.key=element_blank(),
		aspect.ratio=0.5)


