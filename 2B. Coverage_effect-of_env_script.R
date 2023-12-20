setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS\\Comb")
library(dplyr)

new <- read.csv("LatLongCoveragedata.csv")

#remove row 25
new <- new[-25,]

names(new)
#find if EnvLight sig affects coverage at each season
boxplot(Spr2020.cov~EnvLight,data=new)
boxplot(Sum2020.cov~EnvLight,data=new)
boxplot(Aut2020.cov~EnvLight,data=new)
boxplot(Spr2021.cov~EnvLight,data=new)
boxplot(Sum2021.cov~EnvLight,data=new)
boxplot(Aut2021.cov~EnvLight,data=new)
boxplot(Spr2022.cov~EnvLight,data=new)
boxplot(Sum2022.cov~EnvLight,data=new)

t.test(Spr2020.cov~EnvLight,new) #non sig but dll higher
t.test(Sum2020.cov~EnvLight,new) #non sig but dll higher
t.test(Aut2020.cov~EnvLight,new) #non sig but dll higher
t.test(Spr2021.cov~EnvLight,new) #non sig but dll higher
t.test(Sum2021.cov~EnvLight,new) #non sig but dll higher
t.test(Aut2021.cov~EnvLight,new) #non sig but dll higher
t.test(Aut2020.cov~EnvLight,new) #non sig but dll higher
t.test(Spr2022.cov~EnvLight,new) #non sig but dll higher
t.test(Sum2022.cov~EnvLight,new) #non sig but dll higher
