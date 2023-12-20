#combine GR and Col count data
setwd("C:\\Users\\kelli\\OneDrive - The University of Nottingham\\Nottingham\\Year 2\\UK Lemna - Bradford - Hull\\HL LL analysis\\ALL REPS")

library(dplyr)
library(tidyr)

#read in averaged data
gr <- read.csv("Fluo+Chl+HS_stitched.csv")

#TotChl
plot(TotChl ~ fqfm_L1, gr)
plot(TotChl ~ fqfm_L11, gr)
plot(TotChl ~ ChlAB_rat, gr)
plot(TotChl ~ Car_rat, gr) #best so far?
plot(TotChl ~ NPQ_L1, gr)
plot(TotChl ~ NPQ_L11, gr)

plot(Car_rat ~ NPQ_L1, gr)
plot(Car_rat ~ NPQ_L11, gr)
plot(Carot..mg.g. ~ NPQ_L1, gr)
plot(Carot..mg.g. ~ NPQ_L11, gr)
plot(Carot..mg.g. ~ fqfm_L1, gr)
plot(Carot..mg.g. ~ fqfm_L11, gr)
plot(TotChl ~ QY_max, gr)
plot(Chl.a..mg.g. ~ QY_max, gr)

plot(Light_intensity ~ NPQ_L1, gr)
plot(Light_intensity ~ NPQ_L11, gr)
plot(Light_intensity ~ fqfm_L1, gr)
plot(Light_intensity ~ fqfm_L11, gr)

plot(fqfm_L3 ~ fqfm_L5, gr) # GOOD COR
plot(NPQ_L3 ~ NPQ_L5, gr) #GOOD COR

#correlation between FW and FDW
corr <- cor.test(gr$fqfm_L1, gr$TotChl,
                 method = "pearson"
)
corr$p.value
corr$estimate
#low

corr <- cor.test(gr$fqfm_L11, gr$TotChl,
                 method = "pearson"
)
corr$p.value
corr$estimate
#low

corr <- cor.test(gr$NPQ_L1, gr$Car_rat,
                 method = "pearson"
)
corr$p.value
corr$estimate
#low

corr <- cor.test(gr$NPQ_L11, gr$Car_rat,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$ChlAB_rat, gr$TotChl,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Car_rat, gr$TotChl,
                 method = "pearson"
)
corr$p.value
corr$estimate 
# -0.5145569 8.728677e-32

corr <- cor.test(gr$Car_rat, gr$NPQ_L1,
                 method = "pearson"
)
corr$p.value
corr$estimate 
# low

corr <- cor.test(gr$Car_rat, gr$NPQ_L11,
                 method = "pearson"
)
corr$p.value
corr$estimate 
# low

corr <- cor.test(gr$Carot..mg.g., gr$NPQ_L1,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Carot..mg.g., gr$NPQ_L11,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Carot..mg.g., gr$fqfm_L11,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Carot..mg.g., gr$fqfm_L1,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$TotChl, gr$QY_max,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Chl.a..mg.g., gr$QY_max,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Light_intensity, gr$NPQ_L1,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Light_intensity, gr$NPQ_L11,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Light_intensity, gr$fqfm_L1,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

corr <- cor.test(gr$Light_intensity, gr$fqfm_L11,
                 method = "pearson"
)
corr$p.value
corr$estimate 
#low

flx <- gr
flx_ind <- which(names(flx) == "fqfm_L1")

flx$Species <- as.factor(flx$Species)
flx$Treatment <- as.factor(flx$Treatment)
flx$Accession <- as.factor(flx$Accession)

flx$Treatment

flx %>% ggplot(aes(Accession, fqfm_L1, col = Species)) + geom_boxplot() + facet_wrap(~Treatment) + coord_flip()

flx[flx$Treatment == "HL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, NPQ_L11), NPQ_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "High light (HL)", x = "Accession", y = "NPQ") + theme_minimal() + theme(text = element_text(size = 18))
flx[flx$Treatment == "LL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, NPQ_L11), NPQ_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "Low light (LL)", x = "Accession", y = "NPQ") + theme_minimal() + theme(text = element_text(size = 18))

flx[flx$Treatment == "HL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, fqfm_L11), fqfm_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "High light (HL)", x = "Accession", y = "fqfm") + theme_minimal() + theme(text = element_text(size = 18))
flx[flx$Treatment == "LL", ] %>% ggplot(aes(forcats::fct_reorder(Accession, fqfm_L11), fqfm_L11, col = Species)) + geom_boxplot() + coord_flip() + labs(title = "Low light (LL)", x = "Accession", y = "fqfm") + theme_minimal() + theme(text = element_text(size = 18))

#lines through scatter plots according to species or treat
flx %>% ggplot(aes(fqfm_L11, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(NPQ_L3, NPQ_L5, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(NPQ_L1, NPQ_L3, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(NPQ_L1, NPQ_L5, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(QY_max, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(QY_max, fqfm_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(QY_max, TotChl, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(Car_rat, TotChl, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(QY_max, ChlAB_rat, col = Treatment)) + geom_point() + stat_smooth(method = "lm")
flx %>% ggplot(aes(Car_rat, NPQ_L11, col = Treatment)) + geom_point() + stat_smooth(method = "lm")


#nicer plot
library(ggplot2)
library(ggpubr)
ggplot(gr, aes(Car_rat, NPQ_L11, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Carotenoid ratio (mg/g)") +
  ylab("NPQ high light") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

ggplot(gr, aes(TotChl, fqfm_L1, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Tot Chl A+B (mg/g)") +
  ylab("Photosynthetic efficiency QY") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

ggplot(gr, aes(fqfm_L5, ChlAB_rat, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Photosynthetic rate 350 umol") +
  ylab("Chl A/B ratio (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

#doesnt work well with col lag
ggplot(gr, aes(Car_rat, TotChl, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Carotenoid ratio (mg/g)") +
  ylab("ChlA+B (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)

#carotenoid chlorophyllA+B relationship
ggplot(gr, aes(Carot..mg.g., TotChl, label = Accession)) +    # ggplot2 plot with labels
  geom_point() +
  geom_text(aes(label = Accession), hjust = - 0.5) +
  geom_smooth(method='lm') +
  xlab("Carotenoids (mg/g)") +
  ylab("ChlA+B (mg/g)") +
  theme_bw() +
  theme_classic() +
  stat_cor(method = "pearson", label.x = 0.08, label.y = 6)
