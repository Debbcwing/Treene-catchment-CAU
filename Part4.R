rm(list=ls())
gc()

setwd("wd") #Put working directory

# From part 1 - 3,
load("Input.RData")
load("LMMdata.RData")
load("LMM.RData")

#### Selected GLMM models ####
# LMMSpRch = Richness ~ pH + NPR + PO4 + WPAS + DEPTH + H20 + (1 | Year:Subbasin)
# LMMAbun = Abundance ~ DO + PO4 + H36 + H55 + FR + UR + WATR + (1 | Year:Subbasin)
# LMMBiom = Biomass ~ NO3 + TSS + H55 + FR + WATR + PO4 + (1 | Year)
#####

######################################
#### Scenarios study (Projection) ####
######################################
# Scen 1: 15% increment of PO4----
# creating scen1 data
scen1data <- Av
scen1data$PO4 <- Av$PO4*1.15

# Species richness ----
scen1SRdata <- as.data.frame(cbind(Resp.var$Richness,scen1data,SB,Year))
(scen1SR <- as.data.frame(predict(LMMSpRch, newdata = scen1SRdata, 
                                  type = "response", allow.new.level=T)))
scen1SR <- cbind(scen1SR,Resp.var$Richness)
colnames(scen1SR) <- c("Predicted","Observed")
scen1SR$Site <- 1:118
#####
# Overall abundance ----
scen1ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen1data,SB,Year))
(scen1AB <- as.data.frame(predict(LMMAbun, newdata = scen1ABdata, 
                                  type = "response", allow.new.level=T)))
scen1AB <- cbind(scen1AB,Resp.var$Abundance)
colnames(scen1AB) <- c("Predicted","Observed")
scen1AB$Site <- 1:118
#####

# Scen 2: 30% increment of PO4----
# creating scen2 data
scen2data <- Av
scen2data$PO4 <- Av$PO4*1.3

# Species richness ----
scen2SRdata <- as.data.frame(cbind(Resp.var$Richness,scen2data,SB,Year))
(scen2SR <- as.data.frame(predict(LMMSpRch, newdata = scen2SRdata, 
                                  type = "response", allow.new.level=T)))
scen2SR <- cbind(scen2SR,Resp.var$Richness)
colnames(scen2SR) <- c("Predicted","Observed")
scen2SR$Site <- 1:118
#####
# Overall abundance ----
scen2ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen2data,SB,Year))
(scen2AB <- as.data.frame(predict(LMMAbun, newdata = scen2ABdata, 
                                  type = "response", allow.new.level=T)))
scen2AB <- cbind(scen2AB,Resp.var$Abundance)
colnames(scen2AB) <- c("Predicted","Observed")
scen2AB$Site <- 1:118
#####

# Scen 3: 10% increase WPAS; 10% decrease NO3; 10% decrease PO4; 
# 10% decrease AGRL (Not available in variables)-----
scen3data <- Av
scen3data$WPAS <- Av$WPAS*1.1
scen3data$PO4 <- Av$PO4*0.9
scen3data$NO3 <- Av$NO3*0.9

# Species richness ----
scen3SRdata <- as.data.frame(cbind(Resp.var$Richness,scen3data,SB,Year))
(scen3SR <- as.data.frame(predict(LMMSpRch, newdata = scen3SRdata, 
                                  type = "response", allow.new.level=T)))
scen3SR <- cbind(scen3SR,Resp.var$Richness)
colnames(scen3SR) <- c("Predicted","Observed")
scen3SR$Site <- 1:118
#####
# Overall abundance ----
scen3ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen3data,SB,Year))
(scen3AB <- as.data.frame(predict(LMMAbun, newdata = scen3ABdata, 
                                  type = "response", allow.new.level=T)))
scen3AB <- cbind(scen3AB,Resp.var$Abundance)
colnames(scen3AB) <- c("Predicted","Observed")
scen3AB$Site <- 1:118
#####

# Scen 4: 20% increase WPAS; 20% decrease NO3; 20% decrease PO4;
# 10% decrease AGRL (Not available in variables);10% decrease WATR-----
scen4data <- Av
scen4data$WPAS <- Av$WPAS*1.2
scen4data$NO3 <- Av$NO3*0.8
scen4data$WATR <- Av$WATR*0.9

# Species richness ----
scen4SRdata <- as.data.frame(cbind(Resp.var$Richness,scen4data,SB,Year))
(scen4SR <- as.data.frame(predict(LMMSpRch, newdata = scen4SRdata, 
                                  type = "response", allow.new.level=T)))
scen4SR <- cbind(scen4SR,Resp.var$Richness)
colnames(scen4SR) <- c("Predicted","Observed")
scen4SR$Site <- 1:118
#####
# Overall abundance ----
scen4ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen4data,SB,Year))
(scen4AB <- as.data.frame(predict(LMMAbun, newdata = scen4ABdata, 
                                  type = "response", allow.new.level=T)))
scen4AB <- cbind(scen4AB,Resp.var$Abundance)
colnames(scen4AB) <- c("Predicted","Observed")
scen4AB$Site <- 1:118
#####

# Scen 5: 5% reduction precipitation (The hydrological indices reduced correspondingly)----
scen5data <- Av
scen5data$H20 <- Av$H20*0.95    # Weekly skewness
scen5data$H36 <- Av$H36*0.95    # Monthly skewness
scen5data$H55 <- Av$H55*0.95    # Rate of change - difference in flow median

# Species richness ----
scen5SRdata <- as.data.frame(cbind(Resp.var$Richness,scen5data,SB,Year))
(scen5SR <- as.data.frame(predict(LMMSpRch, newdata = scen5SRdata, 
                                  type = "response", allow.new.level=T)))
scen5SR <- cbind(scen5SR,Resp.var$Richness)
colnames(scen5SR) <- c("Predicted","Observed")
scen5SR$Site <- 1:118
#####
# Overall abundance ----
scen5ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen5data,SB,Year))
(scen5AB <- as.data.frame(predict(LMMAbun, newdata = scen5ABdata, 
                                  type = "response", allow.new.level=T)))
scen5AB <- cbind(scen5AB,Resp.var$Abundance)
colnames(scen5AB) <- c("Predicted","Observed")
scen5AB$Site <- 1:118
#####

# Scen 6: 15% reduction precipitation----
scen6data <- Av
scen6data$H20 <- Av$H20*0.85    # Weekly skewness
scen6data$H36 <- Av$H36*0.85    # Monthly skewness
scen6data$H55 <- Av$H55*0.85    # Rate of change - difference in flow median

# Species richness ----
scen6SRdata <- as.data.frame(cbind(Resp.var$Richness,scen6data,SB,Year))
(scen6SR <- as.data.frame(predict(LMMSpRch, newdata = scen6SRdata, 
                                  type = "response", allow.new.level=T)))
scen6SR <- cbind(scen6SR,Resp.var$Richness)
colnames(scen6SR) <- c("Predicted","Observed")
scen6SR$Site <- 1:118
#####
# Overall abundance ----
scen6ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen6data,SB,Year))
(scen6AB <- as.data.frame(predict(LMMAbun, newdata = scen6ABdata, 
                                  type = "response", allow.new.level=T)))
scen6AB <- cbind(scen6AB,Resp.var$Abundance)
colnames(scen6AB) <- c("Predicted","Observed")
scen6AB$Site <- 1:118
#####

# Scen 7: 5% increase precipitation----
scen7data <- Av
scen7data$H20 <- Av$H20*1.05    # Weekly skewness
scen7data$H36 <- Av$H36*1.05    # Monthly skewness
scen7data$H55 <- Av$H55*1.05    # Rate of change - difference in flow median

# Species richness ----
scen7SRdata <- as.data.frame(cbind(Resp.var$Richness,scen7data,SB,Year))
(scen7SR <- as.data.frame(predict(LMMSpRch, newdata = scen7SRdata, 
                                  type = "response", allow.new.level=T)))
scen7SR <- cbind(scen7SR,Resp.var$Richness)
colnames(scen7SR) <- c("Predicted","Observed")
scen7SR$Site <- 1:118
#####
# Overall abundance ----
scen7ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen7data,SB,Year))
(scen7AB <- as.data.frame(predict(LMMAbun, newdata = scen7ABdata, 
                                  type = "response", allow.new.level=T)))
scen7AB <- cbind(scen7AB,Resp.var$Abundance)
colnames(scen7AB) <- c("Predicted","Observed")
scen7AB$Site <- 1:118
#####

# Scen 8: 15% increase precipitation------
scen8data <- Av
scen8data$H20 <- Av$H20*1.15    # Weekly skewness
scen8data$H36 <- Av$H36*1.15    # Monthly skewness
scen8data$H55 <- Av$H55*1.15    # Rate of change - difference in flow median

# Species richness ----
scen8SRdata <- as.data.frame(cbind(Resp.var$Richness,scen8data,SB,Year))
(scen8SR <- as.data.frame(predict(LMMSpRch, newdata = scen8SRdata, 
                                  type = "response", allow.new.level=T)))
scen8SR <- cbind(scen8SR,Resp.var$Richness)
colnames(scen8SR) <- c("Predicted","Observed")
scen8SR$Site <- 1:118
#####
# Overall abundance ----
scen8ABdata <- as.data.frame(cbind(Resp.var$Abundance,scen8data,SB,Year))
(scen8AB <- as.data.frame(predict(LMMAbun, newdata = scen8ABdata, 
                                  type = "response", allow.new.level=T)))
scen8AB <- cbind(scen8AB,Resp.var$Abundance)
colnames(scen8AB) <- c("Predicted","Observed")
scen8AB$Site <- 1:118
#####

#### 1. Statistical test for significant difference with observed data----
scenSR <- as.data.frame(cbind(
  Resp.var$Richness,scen1SR$Predicted, scen2SR$Predicted,scen3SR$Predicted, 
  scen4SR$Predicted, scen5SR$Predicted,scen6SR$Predicted, scen7SR$Predicted, 
  scen8SR$Predicted))
colnames(scenSR) <- c("Observed","Scen1","Scen2","Scen3","Scen4","Scen5",
                          "Scen6","Scen7","Scen8")

scenAB <- as.data.frame(cbind(
  Resp.var$Abundance,scen1AB$Predicted, scen2AB$Predicted,scen3AB$Predicted, 
  scen4AB$Predicted, scen5AB$Predicted,scen6AB$Predicted, scen7AB$Predicted, 
  scen8AB$Predicted))
colnames(scenAB) <- c("Observed","Scen1","Scen2","Scen3","Scen4","Scen5",
                      "Scen6","Scen7","Scen8")
### Wilcoxon test----
wtestSR <- sapply(
  lapply(scenSR[-1], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))
(wtestSR <- as.data.frame(t(wtestSR)))

wtestAB <- sapply(
  lapply(scenAB[-1], function(x) {
  wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))
(testAB <- as.data.frame(t(wtestAB)))

### Kruskal-wallis test ----
ktestSR <- sapply(
  lapply(scenSR[-1], function(x) {
  kruskal.test(x~scenSR$Observed)}), `[`, c("statistic","p.value"))
(ktestSR <- as.data.frame(t(ktestSR)))

ktestAB <- sapply(
  lapply(scenAB[-1], function(x) {
  kruskal.test(x~scenAB$Observed)}), `[`, c("statistic","p.value"))
(ktestAB <- as.data.frame(t(ktestAB)))

#####

#### No significant difference between the observed and simulated data

### 2. Statistical test for significant difference with observed data 
### (within subbasins)-----
### Wilcoxon test----
scenSR <- cbind(scenSR, SB)
wtest_Tr <- sapply(
  lapply(scenSR[ which(scenSR$Subbasin == 'Tr'),2:9], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Sa <- sapply(
  lapply(scenSR[ which(scenSR$Subbasin == 'Sa'),2:9], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Je <- sapply(
  lapply(scenSR[ which(scenSR$Subbasin == 'Je'),2:9], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Bo <- sapply(
  lapply(scenSR[ which(scenSR$Subbasin == 'Bo'),2:9], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Ki <- sapply(
  lapply(scenSR[ which(scenSR$Subbasin == 'Ki'),2:9], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Ju <- sapply(
  lapply(scenSR[ which(scenSR$Subbasin == 'Ju'),2:9], function(x) {
    wilcox.test(x,scenSR$Observed, data = scenSR, conf.level = 0.95)}),
  `[`, c("p.value"))

write.table(cbind(wtest_Tr, wtest_Sa, wtest_Ki, wtest_Ju, wtest_Je, 
                  wtest_Bo), "wtest_SB.txt")

scenAB <- cbind(scenAB, SB)
wtest_Tr <- sapply(
  lapply(scenAB[ which(scenAB$Subbasin == 'Tr'),2:9], function(x) {
    wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Sa <- sapply(
  lapply(scenAB[ which(scenAB$Subbasin == 'Sa'),2:9], function(x) {
    wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Je <- sapply(
  lapply(scenAB[ which(scenAB$Subbasin == 'Je'),2:9], function(x) {
    wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Bo <- sapply(
  lapply(scenAB[ which(scenAB$Subbasin == 'Bo'),2:9], function(x) {
    wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Ki <- sapply(
  lapply(scenAB[ which(scenAB$Subbasin == 'Ki'),2:9], function(x) {
    wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))

wtest_Ju <- sapply(
  lapply(scenAB[ which(scenAB$Subbasin == 'Ju'),2:9], function(x) {
    wilcox.test(x,scenAB$Observed, data = scenAB, conf.level = 0.95)}),
  `[`, c("p.value"))

write.table(cbind(wtest_Tr, wtest_Sa, wtest_Ki, wtest_Ju, wtest_Je, 
                  wtest_Bo), "wtestAB_SB.txt")

#####

### 2. Boxplot------
pdf("Scen_boxplot.pdf", width = 9, height = 5)
par(mar=c(3,3,3,3))
par(mfrow=c(2,1))
boxplot(scenSR[1:9], main = "Predicted species richness in 8 scenarios",
        xlab = "Scenarios", ylab = "Species richness")
boxplot(scenAB[1:9], main = "Predicted abundance in 8 scenarios",
        xlab = "Scenarios", ylab = "Abundance")
dev.off()

write.table(cbind(scenSR,SB),"scenSR.txt")
write.table(cbind(scenAB,SB),"scenAB.txt")
