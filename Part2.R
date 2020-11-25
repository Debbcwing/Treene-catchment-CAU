rm(list=ls())
gc()

setwd("wd") #Put working directory

# From Part 1
load("Input.RData")
##########################################################################
###                      Fixed Effects Model - glm                     ###
###         Testing fixed effect terms via ML is recommended           ###
###  GLM1 - Global model (Complete fixed var.)                         ###
###  GLM2 - Regional model (Only important parameters)                 ###
###  GLM3 - Local interaction model                                    ###
###         (Important parameters with interactions)                   ###
###  GLM4 - Local model (First 3 important parameters)                 ###
###  GLM5 - Selected significant parameters from global model (GLM1)   ###
###  GLM6 - Selected significant parameters from regional model (GLM2) ###
##########################################################################
# Packages -----
library(MuMIn)
library(MASS)
library(rsq)

# 1. Species richness ----
SRdata <- as.data.frame(Resp.var$Richness)
colnames(SRdata) <- c("Richness")
SRdata <- cbind(SRdata,Av)

# pH+H55+NPR+WPAS+DEPTH+PO4+H29+WIDTH+SO4+H41+H21+H36
(GLMsp1 <- glm(Richness ~ pH+DO+PO4+NH4+NO3+NO2+NPR+Cl+SO4+TSS+H11+H12+H13+
                 H20+H21+H29+H36+H41+H55+WIDTH+DEPTH+VELOCITY+FRSE+FR+RNGE+
                 UIDU+UR+WATR+WETL+WPAS, 
               family = poisson(link = "log"), data = SRdata))
(GLMsp2 <- glm(Richness ~ pH+H55+NPR+WPAS+DEPTH+PO4+H29+WIDTH+SO4+H41+H21+H36,
               family = poisson(link = "log"), data = SRdata))
(GLMsp3 <- glm(Richness ~ pH+H55+NPR+(pH:H55)+(pH:NPR)+(H55:NPR),
               family = poisson(link = "log"), data = SRdata))
(GLMsp4 <- glm(Richness ~ pH+H55+NPR, family = poisson(link = "log"),
               data = SRdata))

# Overdispersion check -----
theta <- function(x){x$deviance/x$df.residual}
as.data.frame(matrix(unlist(list(
  theta(GLMsp1),theta(GLMsp2),theta(GLMsp3),theta(GLMsp4)))))

# Over-dispersion --> negative binomial GLM (Harrison, 2014)
(GLMsp1 <- glm.nb(Richness ~ pH+DO+PO4+NH4+NO3+NO2+NPR+Cl+SO4+TSS+H11+H12+
                    H13+H20+H21+H29+H36+H41+H55+WIDTH+DEPTH+VELOCITY+FRSE+
                    FR+RNGE+UIDU+UR+WATR+WETL+WPAS, data = SRdata))
(GLMsp2 <- glm.nb(Richness ~ pH+H55+NPR+WPAS+DEPTH+PO4+H29+WIDTH+SO4+H41+H21+H36,
                  data = SRdata))
(GLMsp3 <- glm.nb(Richness ~ pH+H55+NPR+(pH:H55)+(pH:NPR)+(H55:NPR),
                  data = SRdata))
(GLMsp4 <- glm.nb(Richness ~ pH+H55+NPR, data = SRdata))

# Overdispersion check again
as.data.frame(matrix(unlist(list(
  theta(GLMsp1),theta(GLMsp2),theta(GLMsp3),theta(GLMsp4))))) # No overdispersion
model.sel(GLMsp1,GLMsp2,GLMsp3,GLMsp4) #2134

# Mod5 = selecting significant var. (<0.05) from GLMsp1 (i.e. global mod) as predictors
summary(GLMsp1) # --> pH,DO,PO4,NPR,SO4,H20,DEPTH,FR,WPAS
(GLMsp5 <- glm.nb(Richness ~ pH+DO+PO4+NPR+SO4+H20+DEPTH+FR+WPAS,
                  data = SRdata))
# Mod6 = selecting significant var. from the GLMsp2 as predictors
summary(GLMsp2) # --> pH,NPR,WPAS,PO4,WIDTH,H41
(GLMsp6 <- glm.nb(Richness ~ pH+NPR+WPAS+PO4+WIDTH+H41, data = SRdata))

# Select the best model
GLMsp <- as.data.frame(model.sel(GLMsp1,GLMsp2,GLMsp3,GLMsp4,GLMsp5,GLMsp6)) #562134
as.data.frame(matrix(unlist(list(
  theta(GLMsp1),theta(GLMsp2),theta(GLMsp3),theta(GLMsp4),
  theta(GLMsp5),theta(GLMsp6)))))
# no overdispersion model

# Conclusion
# Model 5 was ranked the best subsequently, with favour to important
# variable selected by the GLM anova. There was no competitive models.

# 2. Overall abundance ----
Abdata <- as.data.frame(Resp.var$Abundance)
colnames(Abdata) <- c("Abundance")
Abdata <- cbind(Abdata, Av)

# WATR+PO4+UR+FRSE+SO4+WPAS+NPR+WETL+H20+WIDTH+UIDU+NO3+H21
(GLMab1 <- glm(Abundance ~ pH+DO+PO4+NH4+NO3+NO2+NPR+Cl+SO4+TSS+H11+H12+
                 H13+H20+H21+H29+H36+H41+H55+WIDTH+DEPTH+VELOCITY+
                 FRSE+FR+RNGE+UIDU+UR+WATR+WETL+WPAS,data = Abdata, 
               family = gaussian(link = "identity")))
(GLMab2 <- glm(Abundance ~ WATR+PO4+UR+FRSE+SO4+WPAS+NPR+WETL+H20+WIDTH+
                 UIDU+NO3+H21,family = gaussian(link = "identity"),
               data = Abdata))
(GLMab3 <- glm(Abundance ~ WATR+PO4+UR+(WATR:PO4)+(WATR:UR)+(PO4:UR),
               data = Abdata,family = gaussian(link = "identity")))
(GLMab4 <- glm(Abundance ~ WATR+PO4+UR,data = Abdata,
               family = gaussian(link = "identity")))

# Gaussian models do not require overdispersion adjustment because they do 
# not assume a specific mean-variance relationship like other distributions 
# (Harrison et al., 2018)

summary(GLMab1) #DO,PO4,H36,H55,WATR
(GLMab5 <- glm(Abundance ~ DO+PO4+H36+H55+WATR,data = Abdata,
               family = gaussian(link = "identity")))
summary(GLMab2) # WATR,PO4,UR, WIDTH
(GLMab6 <- glm(Abundance ~ WATR+PO4+UR+WIDTH,data = Abdata,
               family = gaussian(link = "identity")))
model.sel(GLMab1,GLMab2,GLMab3,GLMab4,GLMab5,GLMab6) #654231
# Conclusion
# Mod 6 has been rated as the best but mod 5 has comparable performance
# with dAICc < 2 and 1 more predictor. They have similar logLik so more
# evaluation has to be made.
# Mod 7 is finally generated as a combination of Mod 5 and Mod 6
(GLMab7 <- glm(Abundance ~ DO+PO4+H36+H55+WATR+UR+WIDTH, data = Abdata,
              family = gaussian(link="identity")))

GLMab <- as.data.frame(model.sel(GLMab1,GLMab2,GLMab3,GLMab4,
                                 GLMab5,GLMab6, GLMab7)) #7546231 

# 3. Overall biomass ----
Biodata <- as.data.frame(Resp.var$Biomass)
colnames(Biodata) <- c("Biomass")
Biodata <- cbind(Biodata, Av)

# PO4+NH4+NO3+NPR+TSS+H20+FRSE+UIDU+UR+WATR+WPAS
(GLMbio1 <- glm(Biomass ~ pH+DO+PO4+NH4+NO3+NO2+NPR+Cl+SO4+TSS+H11+H12+
                  H13+H20+H21+H29+H36+H41+H55+WIDTH+DEPTH+VELOCITY+FRSE+
                  FR+RNGE+UIDU+UR+WATR+WETL+WPAS,data = Biodata, 
                family = gaussian(link = "identity")))
(GLMbio2 <- glm(Biomass ~ PO4+NH4+NO3+NPR+TSS+H20+FRSE+UIDU+UR+WATR+WPAS,
                family = gaussian(link = "identity"),data = Biodata))
(GLMbio3 <- glm(Biomass ~ PO4+NH4+NO3+(PO4:NH4)+(PO4:NO3)+(NH4:NO3),
                data = Biodata,family = gaussian(link = "identity")))
(GLMbio4 <- glm(Biomass ~ PO4+NH4+NO3,data = Biodata,
                family = gaussian(link = "identity")))

# Gaussian models do not require overdispersion adjustment because they do 
# not assume a specific mean-variance relationship like other distributions (Harrison et al., 2018)

summary(GLMbio1) #NO3+TSS+H55+FR+WATR
(GLMbio5 <- glm(Biomass ~ NO3+TSS+H55+FR+WATR,data = Biodata,
                family = gaussian(link = "identity")))
summary(GLMbio2) # WATR,TSS
(GLMbio6 <- glm(Biomass ~ WATR+TSS+PO4,data = Biodata,
                family = gaussian(link = "identity")))

model.sel(GLMbio1,GLMbio2,GLMbio3,GLMbio4,GLMbio5,GLMbio6) #562143
rse(GLMbio5) # 1.13
rse(GLMbio6) # 1.15

rsq(GLMbio5, adj = T) # 0.3763317
rsq(GLMbio6, adj = T) # 0.3533497

(GLMbio7 <- glm(Biomass ~ NO3+TSS+H55+FR+WATR+PO4,data = Biodata,
                family = gaussian(link = "identity")))
model.sel(GLMbio1,GLMbio2,GLMbio3,GLMbio4,GLMbio5,GLMbio6,GLMbio7) #7562143

# Conclusion
# Mod 7 has been rated as the best.
GLMbio <- as.data.frame(model.sel(GLMbio1,GLMbio2,GLMbio3,GLMbio4,GLMbio5,GLMbio6,GLMbio7)) #562143

# Performance result output -----
write.table(rbind(
  GLMsp[c("df","logLik","AICc","delta","weight")],
  GLMab[c("df","logLik","AICc","delta","weight")],
  GLMbio[c("df","logLik","AICc","delta","weight")]),"GLMresult.txt")

# Other evaluation for competitive models (RSE,deviance, adjusted r2) ----
RSE  <- as.data.frame(matrix(unlist(list(
  rse(GLMsp1),rse(GLMsp2), rse(GLMsp3), rse(GLMsp4), rse(GLMsp5), rse(GLMsp6),
  rse(GLMab1), rse(GLMab2),rse(GLMab3),rse(GLMab4),rse(GLMab5),rse(GLMab6),rse(GLMab7),
  rse(GLMbio1), rse(GLMbio2),rse(GLMbio3),rse(GLMbio4),rse(GLMbio5),rse(GLMbio6),rse(GLMbio7))), 
  nrow = 20, byrow = T))
colnames(RSE) <- c("RSE")
rownames(RSE) <- c("GLMsp1","GLMsp2", "GLMsp3", "GLMsp4", "GLMsp5", 
                   "GLMsp6","GLMab1", "GLMab2","GLMab3","GLMab4",
                   "GLMab5","GLMab6","GLMab7","GLMbio1","GLMbio2",
                   "GLMbio3","GLMbio4","GLMbio5","GLMbio6","GLMbio7")

dev  <- as.data.frame(matrix(unlist(list(
  deviance(GLMsp1),deviance(GLMsp2),deviance(GLMsp3),deviance(GLMsp4),deviance(GLMsp5), deviance(GLMsp6),
  deviance(GLMab1),deviance(GLMab2),deviance(GLMab3),deviance(GLMab4),deviance(GLMab5), deviance(GLMab6),
  deviance(GLMab7),deviance(GLMbio1),deviance(GLMbio2),deviance(GLMbio3),deviance(GLMbio4),
  deviance(GLMbio5),deviance(GLMbio6),deviance(GLMbio7))),
  nrow = 20, byrow = T))
colnames(dev) <- c("deviance")

aR2 <- function(x){(rsq(x,adj = T))}
adjR2  <- as.data.frame(matrix(unlist(list(
  aR2(GLMsp1),aR2(GLMsp2),aR2(GLMsp3),aR2(GLMsp4),aR2(GLMsp5),aR2(GLMsp6),
  aR2(GLMab1),aR2(GLMab2),aR2(GLMab3),aR2(GLMab4),aR2(GLMab5),aR2(GLMab6),
  aR2(GLMab7),aR2(GLMbio1),aR2(GLMbio2),aR2(GLMbio3),aR2(GLMbio4),aR2(GLMbio5),
  aR2(GLMbio6),aR2(GLMbio7))),
  nrow = 20, byrow = T))
colnames(adjR2) <- c("adjR2")

write.table(cbind(RSE,dev,adjR2),"glmEval.txt")

# Residual plots -----
GLMrp <- function(x){
  plot(x$residuals~x$fitted.values)
  abline(0,0,lty = 2, col = "red")}

pdf("GLMrp.pdf",width = 8, height = 12)
par(mar=c(2,2,2,2))
par(mfrow=c(10,2))
GLMrp(GLMsp1)
mtext("GLMsp1", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMsp2)
mtext("GLMsp2", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMsp3)
mtext("GLMsp3", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMsp4)
mtext("GLMsp4", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMsp5)
mtext("GLMsp5", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMsp6)
mtext("GLMsp6", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab1)
mtext("GLMab1", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab2)
mtext("GLMab2", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab3)
mtext("GLMab3", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab4)
mtext("GLMab4", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab5)
mtext("GLMab5", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab6)
mtext("GLMab6", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMab7)
mtext("GLMab7", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio1)
mtext("GLMbio1", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio2)
mtext("GLMbio2", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio3)
mtext("GLMbio3", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio4)
mtext("GLMbio4", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio5)
mtext("GLMbio5", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio6)
mtext("GLMbio6", side = 3,cex = 1, outer = F, line = 1)
GLMrp(GLMbio7)
mtext("GLMbio7", side = 3,cex = 1, outer = F, line = 1)

dev.off()
#####

save(SRdata,Abdata,Biodata,file = "LMM_input.RData")
