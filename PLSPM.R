
#####################################################################################
#
# Structural equation modeling (PLS-SEM) for Water Use Efficiency (WUE) modeling
# with stratified bootstrapping
#
# Author: Javier Lopatin | javierlopatin@gmail.com
#
#####################################################################################

library(plspm)

setwd('/home/javier/Documents/temp/IWUE_Chiloe')

### load data
forest <- read.csv('Forest_monthly.csv')[,1:25]
peatland <- read.csv('Peatland_monthly.csv')[,1:25]

# clean unuse variables and keep same order
names(  forest[-c(1,3,4,5,12,14,17,18,20,21,22)] )
names( peatland[-c(1,3,4,5,12,14,17,18,20,21,22)] )

forest <- forest[-c(1,3,4,5,12,14,17,18,20,21,22)]
peatland <- peatland[-c(1,3,4,5,12,14,17,18,20,21,22)]

names(forest) == names(peatland)

forest$Month <- factor(forest$Month)
peatland$Month <- factor(peatland$Month)


################################################################################
### Miscelaneos functions
################################################################################

#------------------------------------------------
## stratify sampling. Classes need to be factor
#------------------------------------------------

stratifySampling <- function(data, classes){
  TRAIN <- list()
  VAL <- list()

  for (i in 1:length( levels(classes) )){
    x = grep( levels(classes)[i], classes  )
    x <- data[x, ]
    # random sampling
    if ( length(x[, 1])==1 ){
      TRAIN[[i]] <- x
      VAL[[i]] <- x
    }
    else {
      idx = sample(1:length(x[,1]), length(x[,1]), replace=TRUE)
      TRAIN[[i]] <- x[idx, ]
      VAL[[i]] <- x[-idx, ]
    }
  }

  # unlist
  TRAIN2 <- do.call("rbind", TRAIN)
  VAL2   <- do.call("rbind", VAL)

  # prepare exit
  output <- list(TRAIN2, VAL2)
  names(output) <- c("train", "validation")
  output

}

#------------------------------------------------
## Run PLSPM with block bootstrap validation for Forest and Peatand
#------------------------------------------------

runSEMBlock <- function(obsVar, inner, outer, modes, models){

  # store values
  modelf <- list()
  modelp <- list()
  predf <- list()
  predp <- list()
  r2f <- list()
  r2p <- list()
  nrmsef <- list()
  nrmsep <- list()
  path_rg_f <- list()
  path_ta_f <- list()
  path_vpd_f <- list()
  path_u_f <- list()
  path_p_f <- list()
  path_swc_f <- list()
  path_wtd_f <- list()
  path_gs_f <- list()
  path_rg_p <- list()
  path_ta_p <- list()
  path_vpd_p <- list()
  path_u_p <- list()
  path_p_p <- list()
  path_swc_p <- list()
  path_wtd_p <- list()
  path_gs_p <- list()

  for(i in 1:iter){
    # Run PLSPM forest
    semf = plspm(na.omit(trainf[[i]]), inner, c(outer, obsVar), modes,
                 boot.val = F, br = 500, scheme = "factor", scaled = T)
    semp = plspm(na.omit(trainp[[i]]), inner, c(outer, obsVar), modes,
                 boot.val = F, br = 500, scheme = "factor", scaled = T)

    scoresf = semf$scores
    scoresp = semp$scores

    # store models
    modelf[[i]] <- semf
    modelp[[i]] <- semp
    # r2
    r2f[[i]] <- semf$inner_summary$R2[length(outer)+1]
    r2p[[i]] <- semp$inner_summary$R2[length(outer)+1]
    # predicted
    if ( models==1 || models==3){
      pred1 <- scoresf[,1]*semf$path_coefs[8,1] + scoresf[,2]*semf$path_coefs[8,2] + scoresf[,3]*semf$path_coefs[8,3] +
        scoresf[,4]*semf$path_coefs[8,4] + scoresf[,5]*semf$path_coefs[8,5] + scoresf[,6]*semf$path_coefs[8,6] +
        scoresf[,7]*semf$path_coefs[8,7]
      pred2 <- scoresp[,1]*semp$path_coefs[8,1] + scoresp[,2]*semp$path_coefs[8,2] + scoresp[,3]*semp$path_coefs[8,3] +
        scoresp[,4]*semp$path_coefs[8,4] + scoresp[,5]*semp$path_coefs[8,5] + scoresp[,6]*semp$path_coefs[8,6] +
        scoresp[,7]*semp$path_coefs[8,7]
      predf[[i]] <- pred1
      predp[[i]] <- pred2
    } else if ( models==2 ){
      pred1 <- scoresf[,1]*semf$path_coefs[9,1] + scoresf[,2]*semf$path_coefs[9,2] + scoresf[,3]*semf$path_coefs[9,3] +
        scoresf[,4]*semf$path_coefs[9,4] + scoresf[,5]*semf$path_coefs[9,5] + scoresf[,6]*semf$path_coefs[9,6] +
        scoresf[,7]*semf$path_coefs[9,7] + scoresf[,8]*semf$path_coefs[9,8]
      pred2 <- scoresp[,1]*semp$path_coefs[9,1] + scoresp[,2]*semp$path_coefs[9,2] + scoresp[,3]*semp$path_coefs[9,3] +
        scoresp[,4]*semp$path_coefs[9,4] + scoresp[,5]*semp$path_coefs[9,5] + scoresp[,6]*semp$path_coefs[9,6] +
        scoresp[,7]*semp$path_coefs[9,7] + scoresp[,8]*semp$path_coefs[9,8]
    }

    predf[[i]] <- pred1
    predp[[i]] <- pred2

    # nRMSE
    if ( models==1 || models==3 ){
      rmse = sqrt(mean((scoresf[,8] - pred1)^2))
      nrmsef[[i]] <- (rmse/(max(scoresf[,8])-min(scoresf[,8])))*100
      rmse = sqrt(mean((scoresp[,8] - pred2)^2))
      nrmsep[[i]] <- (rmse/(max(scoresp[,8])-min(scoresp[,8])))*100
    } else {
      rmse = sqrt(mean((scoresf[,9] - pred1)^2))
      nrmsef[[i]] <- (rmse/(max(scoresf[,9])-min(scoresf[,8])))*100
      rmse = sqrt(mean((scoresp[,9] - pred2)^2))
      nrmsep[[i]] <- (rmse/(max(scoresp[,9])-min(scoresp[,8])))*100
    }
    # coefficients
    if ( models==1 ){
      path_rg_f[[i]]  <- semf$path_coefs[8,1]
      path_ta_f[[i]]  <- semf$path_coefs[8,2]
      path_vpd_f[[i]] <- semf$path_coefs[8,3]
      path_u_f[[i]]   <- semf$path_coefs[8,4]
      path_p_f[[i]]   <- semf$path_coefs[8,5]
      path_swc_f[[i]] <- semf$path_coefs[8,6]
      path_wtd_f[[i]] <- semf$path_coefs[8,7]
      path_rg_p[[i]]  <- semp$path_coefs[8,1]
      path_ta_p[[i]]  <- semp$path_coefs[8,2]
      path_vpd_p[[i]] <- semp$path_coefs[8,3]
      path_u_p[[i]]   <- semp$path_coefs[8,4]
      path_p_p[[i]]   <- semp$path_coefs[8,5]
      path_swc_p[[i]] <- semp$path_coefs[8,6]
      path_wtd_p[[i]] <- semp$path_coefs[8,7]
    } else if ( models==2 ){
      path_rg_f[[i]]  <- semf$path_coefs[9,1]
      path_ta_f[[i]]  <- semf$path_coefs[9,2]
      path_vpd_f[[i]] <- semf$path_coefs[9,3]
      path_u_f[[i]]   <- semf$path_coefs[9,4]
      path_p_f[[i]]   <- semf$path_coefs[9,5]
      path_swc_f[[i]] <- semf$path_coefs[9,6]
      path_wtd_f[[i]] <- semf$path_coefs[9,7]
      path_gs_f[[i]]  <- semf$path_coefs[9,8]
      path_rg_p[[i]]  <- semp$path_coefs[9,1]
      path_ta_p[[i]]  <- semp$path_coefs[9,2]
      path_vpd_p[[i]] <- semp$path_coefs[9,3]
      path_u_p[[i]]   <- semp$path_coefs[9,4]
      path_p_p[[i]]   <- semp$path_coefs[9,5]
      path_swc_p[[i]] <- semp$path_coefs[9,6]
      path_wtd_p[[i]] <- semp$path_coefs[9,7]
      path_gs_p[[i]]  <- semp$path_coefs[9,8]
    } else if ( models==3 ){
      path_rg_f[[i]]  <- semf$path_coefs[8,1]
      path_ta_f[[i]]  <- semf$path_coefs[8,2]
      path_u_f[[i]]   <- semf$path_coefs[8,3]
      path_p_f[[i]]   <- semf$path_coefs[8,4]
      path_swc_f[[i]] <- semf$path_coefs[8,5]
      path_wtd_f[[i]] <- semf$path_coefs[8,6]
      path_gs_f[[i]]  <- semp$path_coefs[8,7]
      path_rg_p[[i]]  <- semp$path_coefs[8,1]
      path_ta_p[[i]]  <- semp$path_coefs[8,2]
      path_u_p[[i]]   <- semp$path_coefs[8,3]
      path_p_p[[i]]   <- semp$path_coefs[8,4]
      path_swc_p[[i]] <- semp$path_coefs[8,5]
      path_wtd_p[[i]] <- semp$path_coefs[8,6]
      path_gs_p[[i]]  <- semp$path_coefs[8,7]
    }
  }

  # prepare output
  out_model <- list(modelf, modelp); names(out_model) <- c('Forest', 'Peatland')
  out_r2 <- list(unlist(r2f), unlist(r2p)); names(out_r2) <- c('Forest', 'Peatland')
  out_nrmse <- list(unlist(nrmsef), unlist(nrmsep)); names(out_nrmse) <- c('Forest', 'Peatland')
  if (models==1){
    out_pathf <- list(unlist(path_rg_f), unlist(path_ta_f), unlist(path_vpd_f), unlist(path_u_f),
                      unlist(path_p_f), unlist(path_swc_f), unlist(path_wtd_f))
    names(out_pathf) <- c('Rg', 'Ta', 'VPD', 'U', 'P', 'SWC', 'WTD')
  } else if  (models==2){
    out_pathf <- list(unlist(path_rg_f), unlist(path_ta_f), unlist(path_vpd_f), unlist(path_u_f),
                      unlist(path_p_f), unlist(path_swc_f), unlist(path_wtd_f), unlist(path_gs_f))
    names(out_pathf) <- c('Rg', 'Ta', 'VPD', 'U', 'P', 'SWC', 'WTD', 'Gs')
  } else if  (models==3){
    out_pathf <- list(unlist(path_rg_f), unlist(path_ta_f), unlist(path_u_f),
                      unlist(path_p_f), unlist(path_swc_f), unlist(path_wtd_f), unlist(path_gs_f))
    names(out_pathf) <- c('Rg', 'Ta', 'U', 'P', 'SWC', 'WTD', 'Gs')
  }


  if (models==1){
    out_pathp <- list(unlist(path_rg_p), unlist(path_ta_p), unlist(path_vpd_p), unlist(path_u_p),
                      unlist(path_p_p), unlist(path_swc_p), unlist(path_wtd_p))
    names(out_pathp) <- c('Rg', 'Ta', 'VPD', 'U', 'P', 'SWC', 'WTD')

  } else if (models==2){
    out_pathp <- list(unlist(path_rg_p), unlist(path_ta_p), unlist(path_vpd_p), unlist(path_u_p),
                      unlist(path_p_p), unlist(path_swc_p), unlist(path_wtd_p), unlist(path_gs_p))
    names(out_pathp) <- c('Rg', 'Ta', 'VPD', 'U', 'P', 'SWC', 'WTD', 'Gs')

  } else if  (models==3){
    out_pathp <- list(unlist(path_rg_p), unlist(path_ta_p), unlist(path_u_p),
                      unlist(path_p_p), unlist(path_swc_p), unlist(path_wtd_p), unlist(path_gs_p))
    names(out_pathp) <- c('Rg', 'Ta', 'U', 'P', 'SWC', 'WTD', 'Gs')
  }

  # output all
  out_all <- list(out_model, out_r2, out_nrmse, out_pathf, out_pathp)
  names(out_all) <- c('Models', 'R2', 'nRMSE', 'PathF', 'PathP')

  return(out_all)

}

#------------------------------------------------
## Get significant path coefficients
#------------------------------------------------
signif_paths <- function(sem, outer, models){
  if (models==1 || models==3){
    out <- matrix(nrow = 2, ncol = 7)
  } else{
    out <- matrix(nrow = 2, ncol = 8)
  }
  rownames(out) <- c('Forest', 'Peatland')
  colnames(out) <- unlist(outer)
  for(i in 1:length(sem$PathF)){
    p <- quantile(sem$PathF[[i]], c(0.05, 0.95))
    out[1,i] <- sign(p[1])==sign(p[2])
    p <- quantile(sem$PathP[[i]], c(0.05, 0.95))
    out[2,i] <- sign(p[1])==sign(p[2])
  }

  return(out)
}

#------------------------------------------------
## One-sided bootstrapping significance test
#------------------------------------------------
significanceTest <- function(sem, models){
  ### r²
  r2 <- sem$R2$Forest - sem$R2$Peatland
  ### %RMSE
  nRMSE <- sem$nRMSE$Forest - sem$nRMSE$Peatland
  ### pathcoefficients
  rg <- sem$PathF$Rg - sem$PathP$Rg
  ta <- sem$PathF$Ta - sem$PathP$Ta
  u <- sem$PathF$U - sem$PathP$U
  p <- sem$PathF$P - sem$PathP$P
  swc <- sem$PathF$SWC - sem$PathP$SWC
  wtd <- sem$PathF$WTD - sem$PathP$WTD
  if (models==1){
    vpd <- sem$PathF$VPD - sem$PathP$VPD
  } else if (models==2){
    vpd <- sem$PathF$VPD - sem$PathP$VPD
    Gs <- sem$PathF$Gs - sem$PathP$Gs
  } else if (models==3){
    Gs <- sem$PathF$Gs - sem$PathP$Gs
  }
  # prepare output
  if (models==1){
    output <- list(r2, nRMSE, rg, ta, vpd, u, p, swc, wtd)
  } else if (models==2){
    output <- list(r2, nRMSE, rg, ta, vpd, u, p, swc, wtd, Gs)
  } else if (models==3){
    output <- list(r2, nRMSE, rg, ta, u, p, swc, wtd, Gs)
  }
  # matrix of significances
  if (models==1){
    a = matrix(nrow = 9, ncol = 3)
    colnames(a) <- c("0.1", "0.05", "0.001")
    rownames(a) <- c("r2", "nRMSE", "Rg",  "Ta", "VPD", "U", "P", "SWC", "WTD")
  } else if (models==2){
  a = matrix(nrow = 10, ncol = 3)
  colnames(a) <- c("0.1", "0.05", "0.001")
  rownames(a) <- c("r2", "nRMSE", "Rg",  "Ta", "VPD", "U", "P", "SWC", "WTD", "Gs")
  } else if (models==3){
    a = matrix(nrow = 9, ncol = 3)
    colnames(a) <- c("0.1", "0.05", "0.001")
    rownames(a) <- c("r2", "nRMSE", "Rg",  "Ta", "U", "P", "SWC", "WTD", "Gs")
  }

  for (i in 1:nrow(a)){
    # 0.1
    if ( sign(quantile(output[[i]], probs=c(0.1))) == sign(quantile(output[[i]], probs=c(0.9))) ){
      a[i,1] = "True"
    } else{
      a[i,1] = "False"
    }
    # 0.05
    if ( sign(quantile(output[[i]], probs=c(0.05))) == sign(quantile(output[[i]], probs=c(0.95))) ){
      a[i,2] = "True"
    } else{
      a[i,2] = "False"
    }
    # 0.001
    if ( sign(quantile(output[[i]], probs=c(0.001))) == sign(quantile(output[[i]], probs=c(0.995))) ){
      a[i,3] = "True"
    } else{
      a[i,3] = "False"
    }
  }
  return(a)
}

################################################################################
### END LOAD FUNCTIONS
################################################################################

################################################################################
### SETTING GENERAL PLSPM PARAMETERS
################################################################################

# =======================================================
# Set the inner model
# =======================================================

# Set the rows of the inner model matrix.
# 1 = recive an interaction
# 0 = do not recive an interaction

# predictors
# for GPP, ET
Rg  = c(0, 0, 0, 0, 0, 0, 0, 0)
Ta  = c(0, 0, 0, 0, 0, 0, 0, 0)
VPD = c(0, 0, 0, 0, 0, 0, 0, 0)
U   = c(0, 0, 0, 0, 0, 0, 0, 0)
P   = c(0, 0, 0, 0, 0, 0, 0, 0)
SWC = c(0, 0, 0, 0, 0, 0, 0, 0)
WTD = c(0, 0, 0, 0, 0, 0, 0, 0)
Gs  = c(0, 0, 0, 0, 0, 0, 0, 0)
# for WUE
Rg2  = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
Ta2  = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
VPD2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
U2   = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
P2   = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
SWC2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
WTD2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
Gs2  = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
# observed variables
GPP  = c(1, 1, 1, 1, 1, 1, 1, 0)
ET   = c(1, 1, 1, 1, 1, 1, 1, 0)
WUE  = c(1, 1, 1, 1, 1, 1, 1, 1, 0)
IWUE = c(1, 1, 1, 1, 1, 1, 1, 0)
iWUE = c(1, 1, 1, 1, 1, 1, 1, 0)

# matrix created by row binding. Example if GPP
innerGPP = rbind(Rg, Ta, VPD, U, P, SWC, WTD, GPP); colnames(innerGPP) = rownames(innerGPP)
innerET = rbind(Rg, Ta, VPD, U, P, SWC, WTD, ET); colnames(innerET) = rownames(innerET)
innerWUE = rbind(Rg2, Ta2, VPD2, U2, P2, SWC2, WTD2, Gs2, WUE); colnames(innerWUE) = rownames(innerWUE)
innerIWUE = rbind(Rg, Ta, U, P, SWC, WTD, Gs, IWUE); colnames(innerIWUE) = rownames(innerIWUE)
inneriWUE = rbind(Rg, Ta, VPD, U, P, SWC, WTD, iWUE); colnames(inneriWUE) = rownames(inneriWUE)

# plot the inner matrix
innerplot(inneriWUE)

# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
modes1 = rep("A", 8)
modes2 = rep("A", 9)

# =======================================================
# Set the outer model
# =======================================================

outer1 = list('Rg', 'Ta', 'VPD', 'U', 'P', 'SWC', 'WTD')
outer2 = list('Rg', 'Ta', 'VPD', 'U', 'P', 'SWC', 'WTD', 'Gs')
outer3 = list('Rg', 'Ta', 'U', 'P', 'SWC', 'WTD', 'Gs')

# ==============================================================
# Create train/val data for block validation with bootstrapping
# ==============================================================

##########################################################
### RUN ONLY ONCE TO KEEP CONSISTENCY IF REVISIT IS NEEDED
##########################################################

# number of bootstrapping iteration
iter = 1000
# store samples per season
trainf <- list()
trainp <- list()
validationf <- list()
validationp <- list()
for(i in 1:iter){
  set.seed(i)
  idf <- stratifySampling(forest, forest$Month)
  set.seed(i)
  idp <- stratifySampling(peatland, forest$Month)
  trainf[[i]] <- idf$train
  trainp[[i]] <- idp$train
  validationf[[i]] <- idf$validation
  validationp[[i]] <- idp$validation
}

save(trainf, trainp, validationf, validationp, file='iterationData.RData')

###############################################################################
### GPP; model type 1
###############################################################################

###########################
## LOAD ITERATION DATA
load('iterationData.RData')

# models
semGPP <- runSEMBlock('GPP', innerGPP, outer1, modes1, 1)
boxplot(semGPP$R2$Forest, semGPP$R2$Peatland)

# significance of path coefficients
sig_GPP <- signif_paths(semGPP, outer, 1); sig_GPP

# signimicance differences between Forest and Peatland
sig_betw_GPP <- significanceTest(semGPP, 1); sig_betw_GPP

# median distribution values
lapply(semGPP$R2, median)
lapply(semGPP$nRMSE, median)
lapply(semGPP$PathF, median)
lapply(semGPP$PathP, median)

###############################################################################
### ET ; model type 1
###############################################################################

# models
semET<- runSEMBlock('ET', innerET, outer, modes1, 1)
boxplot(semET$R2$Forest, semET$R2$Peatland)

# significance of path cefficients
sig_ET<- signif_paths(semET, outer1, 1); sig_ET

# signimicance differences between Forest and Peatland
sig_betw_ET <- significanceTest(semET, 1); sig_betw_ET

# median distribution values
lapply(semET$R2, median)
lapply(semET$nRMSE, median)
lapply(semET$PathF, median); sig_ET[1,]
lapply(semET$PathP, median); sig_ET[2,]

###############################################################################
### WUE ; model type 2
###############################################################################

# models
semWUE<- runSEMBlock('WUE', innerWUE, outer2, modes2, 2)
boxplot(semWUE$R2$Forest, semWUE$R2$Peatland)

# significance of path cefficients
sig_WUE <- signif_paths(semWUE, outer2, 2); sig_WUE

# signimicance differences between Forest and Peatland
sig_betw_WUE <- significanceTest(semWUE, 2); sig_betw_WUE

# median distribution values
lapply(semWUE$R2, median)
lapply(semWUE$nRMSE, median)
lapply(semWUE$PathF, median);sig_WUE[1,]
lapply(semWUE$PathP, median);sig_WUE[2,]

###############################################################################
### IWUE ; model type 3
###############################################################################

# models
semIWUE<- runSEMBlock('IWUEn', innerIWUE, outer3, modes1, 3)
boxplot(semIWUE$R2$Forest, semIWUE$R2$Peatland)

# significance of path cefficients
sig_IWUE <- signif_paths(semIWUE, outer3, 3); sig_WUE

# signimicance differences bWUEween Forest and Peatland
sig_betw_IWUE <- significanceTest(semIWUE, 3); sig_betw_IWUE

# median distribution values
lapply(semIWUE$R2, median)
lapply(semIWUE$nRMSE, median)
lapply(semIWUE$PathF, median); sig_IWUE[1,]
lapply(semIWUE$PathP, median); sig_IWUE[2,]

###############################################################################
### iWUE ; model type 1
###############################################################################

# models
semiWUE<- runSEMBlock('iWUE', inneriWUE, outer1, modes1, 1)
boxplot(semiWUE$R2$Forest, semiWUE$R2$Peatland)

# significance of path cefficients
sig_iWUE <- signif_paths(semiWUE, outer1, 1); sig_iWUE

# signimicance differences bWUEween Forest and Peatland
sig_betw_iWUE <- significanceTest(semiWUE, 1); sig_betw_iWUE

# median distribution values
lapply(semiWUE$R2, median)
lapply(semiWUE$nRMSE, median)
lapply(semiWUE$PathF, median); sig_iWUE[1,]
lapply(semiWUE$PathP, median); sig_iWUE[2,]

save.image('PLSPM.RData')
