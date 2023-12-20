
library(raster)

setwd("~/OneDrive - The University of Nottingham/Essex/Arabidopsis natural variation/Merged datasets")

ecotypes <- read.csv("Ecotypes V2.csv", header = T)

setwd("~/OneDrive - The University of Nottingham/Essex/Arabidopsis natural variation/Climate analyses/")

r <- getData("worldclim",var="bio",res=10)
coordinates <- as.data.frame(cbind(ecotypes$Longitude, ecotypes$Latitude))
colnames(coordinates) <- c("Longitude", "Latitude")
climatedata_ecotypes <- extract(r, coordinates[,1:2])

ecotypes_2 <- cbind(ecotypes, climatedata_ecotypes)

write.csv(ecotypes_2, "Ecotypes V3.csv", row.names = F)
