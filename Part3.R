rm(list=ls())
gc()

setwd("wd") #Put working directory

# From part 1 and 2
load("Input.RData")
load("LMM_input.RData")

############################################
###      Mixed Effects Model - lmer      ###
############################################                                  
# Packages -----
library(lme4)
library(MASS)
library(car)
library(MuMIn)

# 1. Species richness -----
# Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS (GLMsp5)
SRdata <- cbind(SRdata,Year,SB)

(LMMsp1 <- lmer(Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS+(1|Year)+
                     (1|Subbasin), data = SRdata))
(LMMsp2 <- lmer(Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS+(1|Year),
                     data = SRdata))
(LMMsp3 <- lmer(Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS+(1|Subbasin),
                     data = SRdata))
(LMMsp4 <- lmer(Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS+(Subbasin|Year),
                     data = SRdata))
(LMMsp5 <- lmer(Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS+(1|Year:Subbasin),
                     data = SRdata))

# Overdispersion checks -----
dfun <- function(x) {
  with(x,(sum(residuals(x)^2))/df.residual(x))
} # Ben Bolker, 2017
as.data.frame(matrix(unlist(list(
  dfun(LMMsp1),dfun(LMMsp2),dfun(LMMsp3),dfun(LMMsp4),dfun(LMMsp5)))))

# qAICc (because of overdispersion models)
QAICc(LMMsp1, LMMsp2,LMMsp3,LMMsp4,LMMsp5, 
      chat = deviance(LMMsp1) / df.residual(LMMsp1))

# Evaluations -----
model.sel(LMMsp1,LMMsp2,LMMsp3,LMMsp4,LMMsp5) 
#AICc = 12534 (125 have close AICc i.e. <2)
#qAICc = 25314 (253 have close qAICc i.e. <2)
BIC(LMMsp1,LMMsp2,LMMsp3,LMMsp4,LMMsp5) #BIC
# 25134 (251 have super close BIC i.e. <2)

# Marginal and conditional r2 (Nagakawa, 2013)
# Trigamma-estimate is recommended whenever available 
# (Nakagawa et al., 2014; Johnson, 2014). 
round(r.squaredGLMM(LMMsp2),3) #0.433 0.573
round(r.squaredGLMM(LMMsp5),3) #0.446 0.569
round(r.squaredGLMM(LMMsp3),3) #0.559 0.587
round(r.squaredGLMM(LMMsp1),3) #0.344 0.574
round(r.squaredGLMM(LMMsp4),3) #0.310 0.579

# Stick to model 3 - LMMsp3
LMMSpRch <- LMMsp3

write.table(cbind(
  QAICc(LMMsp1, LMMsp2,LMMsp3,LMMsp4,LMMsp5, chat = deviance(LMMsp1)/df.residual(LMMsp1)),
  BIC(LMMsp1,LMMsp2,LMMsp3,LMMsp4,LMMsp5),
  model.sel(LMMsp1,LMMsp2,LMMsp3,LMMsp4,LMMsp5)),"LMM_SpRch.txt")

# 2. Overall abundance ----
# Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH
Abdata <- cbind(Abdata, Year, SB)

(LMMab1 <- lmer(Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH+(1|Year)+
                   (1|Subbasin), data = Abdata))
(LMMab2 <- lmer(Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH+(1|Year),
                 data = Abdata))
(LMMab3 <- lmer(Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH+(1|Subbasin),
                 data = Abdata))
(LMMab4 <- lmer(Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH+(Subbasin|Year), 
                 data = Abdata))
(LMMab5 <- lmer(Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH+(1|Year:Subbasin),
                 data = Abdata))

model.sel(LMMab1,LMMab2,LMMab3,LMMab4,LMMab5) #AICc 53124
# 3 and 5 are very close (dAICc~2)
BIC(LMMab1,LMMab2,LMMab3,LMMab4,LMMab5) # BIC 53124
# 3 and 5 are very close (dBIC <2)

# Evaluation -----
# Marginal and conditional r2 (Nagakawa, 2013)
round(r.squaredGLMM(LMMab5),2) #0.48 0.61
round(r.squaredGLMM(LMMab3),2) #0.44 0.54

LMMAbun <- LMMab5   # Overall abundance final model

write.table(cbind(model.sel(LMMab1,LMMab2,LMMab3,LMMab4,LMMab5),
                  BIC(LMMab1,LMMab2,LMMab3,LMMab4,LMMab5)), "LMM_Abun.txt")
#####

# 3. Overall biomass ----
# Biomass ~ NO3+TSS+H55+FR+WATR+PO4
Biodata <- cbind(Biodata, Year, SB)

(LMMbio1 <- lmer(Biomass ~ NO3+TSS+H55+FR+WATR+PO4+(1|Year)+
                   (1|Subbasin), data = Biodata))
(LMMbio2 <- lmer(Biomass ~ NO3+TSS+H55+FR+WATR+PO4+(1|Year),
                 data = Biodata))
(LMMbio3 <- lmer(Biomass ~ NO3+TSS+H55+FR+WATR+PO4+(1|Subbasin),
                 data = Biodata))
(LMMbio4 <- lmer(Biomass ~ NO3+TSS+H55+FR+WATR+PO4+(Subbasin|Year), 
                 data = Biodata))
(LMMbio5 <- lmer(Biomass ~ NO3+TSS+H55+FR+WATR+PO4+(1|Year:Subbasin),
                 data = Biodata))

model.sel(LMMbio1,LMMbio2,LMMbio3,LMMbio4,LMMbio5) #AICc 25314
BIC(LMMbio1,LMMbio2,LMMbio3,LMMbio4,LMMbio5) # BIC 25314

# Marginal and conditional r2 (Nagakawa, 2013)
round(r.squaredGLMM(LMMbio1),3) # 0.411 0.437
round(r.squaredGLMM(LMMbio2),3) # 0.411 0.437
round(r.squaredGLMM(LMMbio3),3) # 0.422 0.422
round(r.squaredGLMM(LMMbio4),3) # 0.43 0.483
round(r.squaredGLMM(LMMbio5),3) # 0.421 0.427

LMMBiom <- LMMbio2   # Overall Biomass final model

write.table(cbind(model.sel(LMMbio1,LMMbio2,LMMbio3,LMMbio4,LMMbio5),
                  BIC(LMMbio1,LMMbio2,LMMbio3,LMMbio4,LMMbio5)), "LMM_bio.txt")
#####

### LMM residual plots-----
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(gridExtra)

pdf("LMMSpResidual.pdf", width = 12, height = 8)
grid.arrange(
  plot_model(LMMsp1, type = "diag", title = "",
             terms = c("Richness"))[[4]],
  plot_model(LMMsp2, type = "diag", title = "LMMsp2",
             terms = c("Richness"))[[4]],
  plot_model(LMMsp3, type = "diag", title = "LMMsp3",
             terms = c("Richness"))[[4]],
  plot_model(LMMsp4, type = "diag", title = "LMMsp4",
             terms = c("Richness"))[[4]],
  plot_model(LMMsp5, type = "diag", title = "LMMsp5",
             terms = c("Richness"))[[4]],
  nrow = 3, ncol = 2, top = "Species richness LMM residual plot")
dev.off()

pdf("LMMabResidual.pdf", width = 12, height = 8)
grid.arrange(
  plot_model(LMMab1, type = "diag", terms = c("Abundance"))[[4]],
  plot_model(LMMab2, type = "diag", terms = c("Abundance"))[[4]],
  plot_model(LMMab3, type = "diag", terms = c("Abundance"))[[4]],
  plot_model(LMMab4, type = "diag", terms = c("Abundance"))[[4]],
  plot_model(LMMab5, type = "diag", terms = c("Abundance"))[[4]],
  nrow = 3, ncol = 2, top = "Overall abundance LMM residual plot")
dev.off()

pdf("LMMbioResidual.pdf", width = 12, height = 8)
grid.arrange(
  plot_model(LMMbio1, type = "diag", terms = c("Biomass"))[[4]],
  plot_model(LMMbio2, type = "diag", terms = c("Biomass"))[[4]],
  plot_model(LMMbio3, type = "diag", terms = c("Biomass"))[[4]],
  plot_model(LMMbio4, type = "diag", terms = c("Biomass"))[[4]],
  plot_model(LMMbio5, type = "diag", terms = c("Biomass"))[[4]],
  nrow = 3, ncol = 2, top = "Overall biomass LMM residual plot")
dev.off()

### LMM/Final model plots -----
pdf("LMM plots.pdf", height = 6, width = 12)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(gridExtra)
grid.arrange(
  plot_model(LMMSpRch, show.values = TRUE, value.offset = .3, 
             title = "Species richness"),
  plot_model(LMMAbun, show.values = TRUE, value.offset = .3,
             title = "Overall abundance"),
  plot_model(LMMBiom, show.values = TRUE, value.offset = .3,
             title = "Overall biomass"),
  nrow = 1, top = "Model coefficients", 
  bottom = "* = conf. int >0.95; *** = conf. int >0.99 ")

library(coefplot2)
par(mar=c(2,2,2,2))
par(mfrow=c(1,3))
coefplot2(LMMSpRch)
mtext("LMM - Species richness",side = 3, line = 1, cex = 0.8)
coefplot2(LMMAbun)
mtext("LMM - Overall abundance",side = 3, line = 1, cex = 0.8)
coefplot2(LMMBiom)
mtext("LMM - Overall biomass",side = 3, line = 1, cex = 0.8)
dev.off()

save(SRdata,Abdata, Biodata, file = "LMMdata.RData")
save(LMMSpRch,LMMAbun,LMMBiom, file = "LMM.RData")