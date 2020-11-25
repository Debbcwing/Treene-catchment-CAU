rm(list=ls())     # cleaning data
gc()

setwd("wd") #Put working directory

# Installing required packages ----
#install.packages(c("readr","ggplot2","ggpubr","usdm",
#"PerformanceAnalytics","corrplot","RColorBrewer","randomForest",
#"MuMIn","MASS","sjstats","lmtest","rsq","lme4","MASS","car","r2glmm",
#"sjPlot","sjlabelled","sjmisc","gridExtra"))

### Loading packages -----
library(readr)
library(ggplot2)
library(ggpubr)
# Species richness and abundance data-----
options(scipen = 999)
Resp.var <- read.csv(paste0(getwd(),"/InputData.csv"),
                     sep = ",", header = T)
Resp.var$Biomass = log(Resp.var$Biomass)
# Log-transformation to cell density for better model-fitting

#################################
### Resp.var data exploration ###
#################################
# Shapiro test on response variables
Resp.var.Shapiro <- lapply(Resp.var, shapiro.test)
Resp.var.Shapiro <- sapply(Resp.var.Shapiro, `[`, 
                           c("statistic","p.value"))
(as.data.frame(t(Resp.var.Shapiro)))
rm(Resp.var.Shapiro)

# vif
library(usdm)
(vifstep(Resp.var, th = 5))

# Boxplot
pdf("Boxplot_RespVar.pdf", width = 5,height = 5)
par(mar=c(3,3,3,3))
par(mfrow=c(1,3))
boxplot(Resp.var$Richness, main = "Species richness", cex.main = 1.5,
        xlab = NULL)
boxplot(Resp.var$Abundance, main = "Relative abundance", cex.main = 1.5,
        xlab = NULL)
boxplot(Resp.var$Biomass, main = "Overall biomass", cex.main = 1.5,
        xlab = NULL)
dev.off()
 
#######################################
### Predictor var. data exploration ###
#######################################
# Reading data (n=86) ------
Av <- read.csv(paste0(getwd(),"/Av_input.csv"),
               sep = ",", header = T)   
# Av = Ev + Hv (118 obs. of 86 variables)
Av <- Av[c(2:87)]

### VIF_all predictor variables (VIF = 5) ###
{
  # Calculating VIF
  # VIF_step (threshold value = 5) -----      
  vifstep(Av, th = 5)     
  
  ### RESULTS ----      
  #  56 variables from the 86 input variables have collinearity problem: 
  #  AGRL FRSD URLD H51 H42 H38 H45 H40 H39 H49 H53 DIN H06 H35 H01 H30 H14 H23 
  #  H27 H22 H07 H15 H08 H03 H09 H24 H02 H57 H16 H33 H31 H17 H26 H04 H52 H05 H18 
  #  H25 H10 H43 H19 H48 H37 H56 H28 H50 URMD H34 H47 TP H54 H32 FRST H46 H44 WT 
  
  #  After excluding the collinear variables, the linear correlation coef 
  #  ranges between: 
  #  min correlation ( FR ~ SO4 ):  -0.0009779467 
  #  max correlation ( NPR ~ NO3 ):  0.6863922 
  
  # VIFs of the remained variables
  #    Variables      VIF
  #  1         pH 2.053689
  #  2         DO 2.256263
  #  3        PO4 2.406707
  #  4        NH4 2.645954
  #  5        NO3 4.806506
  #  6        NO2 2.739314
  #  7        NPR 4.693686
  #  8         Cl 1.658632
  #  9        SO4 3.204664
  #  10       TSS 1.617709
  #  11       H11 2.436272
  #  12       H12 1.403784
  #  13       H13 1.560900
  #  14       H20 3.282168
  #  15       H21 3.132296
  #  16       H29 2.798156
  #  17       H36 2.716718
  #  18       H41 2.225555
  #  19       H55 3.514445
  #  20     WIDTH 3.766191
  #  21     DEPTH 2.355598
  #  22  VELOCITY 1.946403
  #  23      FRSE 2.050814
  #  24        FR 1.764138
  #  25      RNGE 1.690957
  #  26      UIDU 3.068952
  #  27        UR 2.340011
  #  28      WATR 2.799631
  #  29      WETL 1.630201
  #  30      WPAS 2.244379      
  
  ###########      
}

## Subsetting selected Av (remaining = 30nos.) -----
Av[c("AGRL","FRSD","URLD","H51","H42","H38","H45","H40","H39","H49",
     "H53","DIN","H06","H35","H01","H30","H14","H23","H27","H22","H07",
     "H15","H08","H03","H09","H24","H02","H57","H16","H33","H31","H17",
     "H26","H04","H52","H05","H18","H25","H10","H43","H19","H48","H37",
     "H56","H28","H50","URMD","H34","H47","TP","H54","H32","FRST","H46",
     "H44","WT")] <- NULL

## Plotting Av distribution and correlation -----
library(lattice)
pdf("Av_ClevPlot.pdf", width = 10, height = 8)
dotplot(as.matrix(Av), groups = FALSE,
        main = "Cleveland dotplot for predictor variables", 
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data")
dev.off()

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
# matrix of the p-value of the correlation
p.mat <- cor.mtest(Av)
head(p.mat[, 1:5])
pdf("Corrplot_Av2.pdf",width = 15,height = 13)
library(corrplot)
library(RColorBrewer)
corrplot(cor(Av), method = "color", col=brewer.pal(n=10, name="RdYlBu"),
         type="upper", order="hclust", 
         addCoef.col = T, number.cex = 0.8, 
         tl.col="black", tl.srt=45, tl.cex = 0.8,
         p.mat = p.mat, sig.level = 0.01, insig = "blank", diag = F)
mtext(text = "Correlation plot of Av",side = 3, outer = F, line = 3,
      cex = 1.5)
dev.off()
# Random effects (Year and Subbasin)------
t = matrix(data = NA, nrow = 1, ncol = 118) 
t[c(1:59)] = c("2014")
t[c(60:118)] = c("2015")
Year = as.data.frame(t(t))
colnames(Year) = c("Year")
rm(t)

s = matrix(data = NA, nrow = 1, ncol = 118)
s[c(1:22)] = c("Tr")
s[c(23:34)] = c("Bo")
s[c(35:44)] = c("Je")
s[c(45:51)] = c("Ki")
s[c(52:55)] = c("Sa")
s[c(56:59)] = c("Ju")
s[c(60:118)] = s[c(1:59)]
SB = as.data.frame(t(s))
colnames(SB) = c("Subbasin")
rm(s)

#####
################################
#### Random Forest --> VIMP ####
################################
{
  library(randomForest)
  
  # 1. Species richness (500) -----
  rfsprch <- Resp.var$Richness
  rfsprch = cbind(rfsprch, Av)
  
  set.seed(1234)
  (randomForest(rfsprch~.,data=rfsprch, ntree=8000, importance = T))
  #          Mean of squared residuals: 104.6131
  #                   % Var explained: 51.32
  
  set.seed(1234)
  (randomForest(rfsprch~.,data=rfsprch, ntree=5000, importance = T))
  #            Mean of squared residuals: 104.4229
  #                   % Var explained: 51.41
  
  set.seed(1234)
  (randomForest(rfsprch~.,data=rfsprch, ntree=1000, importance = T))
  #            Mean of squared residuals: 103.1878
  #                   % Var explained: 51.99
  
  set.seed(1234)
  (randomForest(rfsprch~.,data=rfsprch, ntree=500, importance = T))
  #            Mean of squared residuals: 106.3382
  #                   % Var explained: 52.09
  
  #### 500 is best-fitted
  ### VImportance -----
  set.seed(1234)
  VIMPsr <- as.data.frame(
    randomForest::importance(
      randomForest(rfsprch~.,data=rfsprch,ntree=500,importance = T),
      type = 1))
  ###################

  # 2. Overall abundance (1000) -----
  rfAbun <- Resp.var$Abundance
  rfAbun = cbind(rfAbun, Av)
  
  set.seed(1234)
  (randomForest(rfAbun~.,data=rfAbun, ntree=8000, importance = T)) 
  # Mean of squared residuals: 1.865665
  # % Var explained: 64.4
  
  set.seed(1234)
  (randomForest(rfAbun~.,data=rfAbun, ntree=5000, importance = T)) 
  # Mean of squared residuals: 1.861052
  # % Var explained: 64.49
  
  set.seed(1234)
  (randomForest(rfAbun~.,data=rfAbun, ntree=1000, importance = T)) 
  # Mean of squared residuals: 1.858746
  # % Var explained: 64.53
  
  set.seed(1234)
  (randomForest(rfAbun~.,data=rfAbun, ntree=500, importance = T)) 
  #  Mean of squared residuals: 1.877572
  # % Var explained: 64.17
  
  
  #### 1000 is best-fitted
  ## VImportance ------
  set.seed(1234)
  VIMPab <- as.data.frame(
    randomForest::importance(
      randomForest(rfAbun~.,data=rfAbun, ntree=1000,importance = T),
      type = 1))
  ##################
  
  # 3. Overall Biomass (1000) -----
  rfBiom <- Resp.var$Biomass
  rfBiom = cbind(rfBiom, Av)
  
  set.seed(1234)
  (randomForest(rfBiom~.,data=rfBiom, ntree=8000, importance = T)) 
  # Mean of squared residuals: 1.207764
  # % Var explained: 40.49
  
  set.seed(1234)
  (randomForest(rfBiom~.,data=rfBiom, ntree=5000, importance = T)) 
  # Mean of squared residuals: 1.205002
  # % Var explained: 40.63
  
  set.seed(1234)
  (randomForest(rfBiom~.,data=rfBiom, ntree=1000, importance = T)) 
  # Mean of squared residuals: 1.192034
  # % Var explained: 41.27
  
  set.seed(1234)
  (randomForest(rfBiom~.,data=rfBiom, ntree=500, importance = T)) 
  # Mean of squared residuals: 1.208602
  # % Var explained: 40.45
  
  #### 1000 is best-fitted
  ## VImportance ------
  set.seed(1234)
  VIMPbio <- as.data.frame(
    randomForest::importance(
      randomForest(rfBiom~.,data=rfBiom, ntree=1000,importance = T),
      type = 1))
  ##################
  
  # Plotting VIMP -----
  pdf("VIMP_new.pdf", width = 10, height = 5)
  par(mar=c(2,2,2,2))
  par(mfrow=c(1,3))
  set.seed(1234)
  varImpPlot(randomForest(rfsprch~.,data=rfsprch, ntree=500, importance = T), 
             type = 1,main = "VIMP - Species richness") 
  set.seed(1234)
  varImpPlot(randomForest(rfAbun~., data=rfAbun,ntree=1000, importance = T),
             type = 1, main = "VIMP - Overall abundance")
  set.seed(1234)
  varImpPlot(randomForest(rfBiom~., data=rfBiom,ntree=1000, importance = T),
             type = 1, main = "VIMP - Overall biomass")
  # Shouldn't use MDI because it is based towards some predictor variables (Strobl et al.,2007)
  dev.off()
  #####
  
  # Writting VIMP ------
  write.table(round(cbind(VIMPsr,VIMPab,VIMPbio),2), "VIMP_new.txt")
}
# Cleaning----
rm(rfAbun,rfsprch,rfBiom)

save(Resp.var,Av,SB,Year, file = "Input.RData")
