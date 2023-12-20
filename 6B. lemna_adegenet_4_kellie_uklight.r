# LOCAL CONDA ENVIRONMENT
# conda activate ade_env

# RUN COMMAND
# Rscript adegenet_PCA_4_kellie.R

options(warn=1)
library(ade4)
library(adegenet)
library(adegraphics)
library(MASS)
library(pegas)
library(StAMPP)
library(vcfR)
library(dplyr)
rm(list = ls())


##################### 
# MODIFIED FUNCTIONS

# a function for conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}


# FAST PCA CODE
## -------------------------------------------------------
### a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
    # if(!is.null(locNames(x))){
        rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
        # rownames(res$loadings) <- paste(locNames(x)) # proper allele names
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
# ---------------------------------------------------------






setwd("C:\\Users\\kelli\\OneDrive\\Documents\\Duckweed\\adgenet\\2023")

#vcf_name <- "filtered.F4.4fds.vcf.gz" # short set for fast testing
#vcf_name <- "uk20_filtered.F4.4fds.vcf" # short set for fast testing
vcf_name <- read.vcfR("uklight.vcf.gz")

# X AXIS PC
PCx = 1
# Y AXIS PC
PCy = 2
# AXIS RANGES
limit = c(-400,400)
# COLORS
bfd_c  <- '#F8A08D'
cwall_c  <- '#2F87B6'
hull_c  <- '#43871D'
york_c  <- '#9628A9'
brist_c <- '#9628A9'
shading <- c(bfd_c,cwall_c,hull_c,york_c, brist_c)

# IMPORT SNP data from VCF
vcf <- read.vcfR(vcf_name)   #read in all data
vcf <- vcf_name

# convert to genlight 	
aa.genlight_x <- vcfR2genlight.tetra(vcf) ## use the modified function vcfR2genlight.tetra at the end of the file

toRemove <- is.na(glMean(aa.genlight_x, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
aa.genlight <- aa.genlight_x[, !toRemove]

aa.genlight$ind.names
length(aa.genlight$ind.names)

labels  = list(
    'ERR3957957'  = c('S. intermedia 8410',  'WW',  'Sintermedia'),
    'L.min'  = c('L. minuta 7760',  'WW',  'Lminuta'),
    'SRR10958765'  = c('L. minor 9441',  'WW',  'Lminor'),
    'SRR10958787'  = c('L. minor 7868',  'WW',  'Lminor'),
    'SRR11472010'  = c('S. polyrhiza 9504',  'WW',  'Spolyrhiza'),
    'SRR8291590'  = c('L. turionifera 9434',  'WW',  'Lturionifera'),
    'SRR8291593'  = c('L. minuta 9581',  'WW',  'Lminuta'),
    'SRR8291594'  = c('L. minuta 9484',  'WW',  'Lminuta'),
    'SRR8291595'  = c('L. minuta 7612',  'WW',  'Lminuta'),
    'SRR8291596'  = c('L. minuta 6717',  'WW',  'Lminuta'),
    'a3'  = c('KS09',  'Bfd',  'Lminor'),
    'a4'  = c('KS12',  'Bfd',  'Spolyrhiza'),
    'a6'  = c('KS15',  'Bfd',  'Lminor'),
    'a7'  = c('KS21',  'Hull',  'Lminor'),
    'onea1'  = c('KS02',  'Bfd',  'Lminor'),
    'onea11' = c('KS16',  'Hull',  'Lturionifera'),
    'onea12' = c('KS17',  'Hull',  'Lminor'),
    'onea13' = c('KS18',  'Hull',  'Lminor'),
    'onea14' = c('KS20',  'Hull',  'Lminuta'),
    'onea16' = c('KS22',  'Hull',  'Lturionifera'),
    'onea17' = c('KS25',  'York',  'Lminuta'),
    'onea18' = c('KS27',  'York',  'Lminor'),
    'onea19' = c('KS28',  'York',  'Lminor'),
    'onea2'  = c('KS03',  'Bfd',  'Lminor'),
    'onea20' = c('KS29',  'York',  'Lminor'),
    'onea21' = c('LY01A', 'Cwall', 'Lminor'),
    'onea22' = c('LY01B', 'Cwall', 'Lminuta'),
    'onea23' = c('LY02',  'Cwall', 'Lminor'),
    'onea24' = c('LY03',  'Cwall', 'Lminor'),
    'onea3'  = c('KS04',  'Bfd',  'Lminor'),
    'onea4'  = c('KS06A', 'Bfd',  'Lminuta'),
    'onea5'  = c('KS06B', 'Bfd',  'Lminuta'),
    'onea8'  = c('KS13',  'Bfd',  'Lminor'),
    'onea9'  = c('KS14',  'Bfd',  'Lminor'))


# extract relevant sample names / populations
samples <- character()
populations <- character()
species <- character()
for (vcf_ind in indNames(aa.genlight)){
    stored = labels[[vcf_ind]]
    samples <- c(samples, stored[1])
    populations <- c(populations, stored[2])
    species <- c(species, stored[3])
}

# update genlight details
indNames(aa.genlight) <- samples
#pop(aa.genlight) <- populations
pop(aa.genlight) <- species

# proportion of explained variance by each axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis
pca.1$eig[4]/sum(pca.1$eig) # proportion of variation explained by 4th axis
pca.1$eig[5]/sum(pca.1$eig) # proportion of variation explained by 5th axis
pca.1$eig[6]/sum(pca.1$eig) # proportion of variation explained by 6th axis
pca.1$eig[7]/sum(pca.1$eig) # proportion of variation explained by 7th axis
#1 and 2 = 43% and 25%, 1 and 3= 43% and 6%


# run PCA

pca.1 <- glPcaFast(aa.genlight, nf=300) # use the modified function glPcaFast at the end of the file

# X AXIS PC
PCx = 1
# Y AXIS PC
PCy = 2
#PCy = 3

#  OUTPUTS
contribution = pca.1$loadings^2
#write.table(contribution[,1:3], sep='\t', file = "contributions")

PCx_var = format(round(pca.1$eig[PCx]/sum(pca.1$eig)*100,2), nsmall=2)
PCy_var = format(round(pca.1$eig[PCy]/sum(pca.1$eig)*100,2), nsmall=2)
PCx_label = paste('PC',PCx,' (',PCx_var,'%)',sep='')
PCy_label = paste('PC',PCy,' (',PCy_var,'%)',sep='')

label_1 = paste('PCA_all_SNPs','.PC',PCx,'.PC',PCy,'.pdf', sep='')
pdf (label_1, width=14, height=7)

#tiff('PCA_species_genotyping.tiff', units="in", width=13, height=7, res=300, compression = 'lzw')
# PLOT; INDIVIDUALS
g0 <- s.label(
  pca.1$scores,
  xax         = PCx,
  yax         = PCy,
  xlim        = limit,
  ylim        = limit, 
  ppoints.col = "red", 
  ppoints.cex = 1,
  plegend.drawKey     = TRUE, 
  plabels     = list(box=list(draw=FALSE), optim=FALSE), 
  paxes.draw  = TRUE, 
  pgrid.draw  = FALSE, 
  plabels.cex = 1, 
  plot        = FALSE, 
  xlab        = PCx_label, 
  ylab        = PCy_label
  )


# PLOT; POPULATIONS  + EIGEN VALUES
g1 <- s.class(
  pca.1$scores, 
  as.factor(as.vector(pop(aa.genlight))), 
  plegend.drawKey     = TRUE, 
  ellipseSize         = 2.5,
  col                 = shading, 
  pellipses.col       = c('transparent'), 
  xax                 = PCx, 
  yax                 = PCy, 
  xlim                = limit, 
  ylim                = limit,  
  ppoints.cex         = 1, 
  paxes.draw          = TRUE, 
  pgrid.draw          = FALSE, 
  plab.cex            = 1.5, 
  plabels     = list(box=list(draw=FALSE), optim=FALSE),
  pellipses.axes.draw = FALSE, 
  xlab                = PCx_label, 
  ylab                = PCy_label, 
  plot                = FALSE
  )

ge <- plotEig(
  pca.1$eig, 
  pca.1$nf, 
  xax         = PCx, 
  yax         = PCy, 
  psub        = list(text = "Eigenvalues"), 
  col.plot    = 'dimgray', 
  pbackground = list(box=TRUE), 
  plot        = FALSE
  )

g1 <- insert(ge, g1, posi=c(0.75,0.0825,0.955,0.25), ratio = 0.2, plot=FALSE)

# SUBPLOT
ADEgS(c(g0, g1), layout = c(1, 2))
dev.off()
#dev.new()
