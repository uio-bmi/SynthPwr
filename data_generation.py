import os, sys
import re
import numpy as np
import pandas as pd
import time
import warnings
import rpy2.robjects as robjects
from matplotlib import pyplot as plt
from rpy2.robjects import pandas2ri
import seaborn as sns
from power_evaluation import plot_power_by_sample, plot_power_simulations_by_sample

# Define static directory names
pandas2ri.activate()
warnings.simplefilter(action='ignore', category=FutureWarning)
starttime = time.time()
dirname = os.sep+"power_experiments"+os.sep
figuredirname = os.sep+"figures"+os.sep
summarycalcdirname = os.sep+"summary_stats"+os.sep

###################################################
# rpy2 R object for extracting a real-world beta_matrix using R's GeoQuery
###################################################

robjects.r('''
            chooseCRANmirror(ind = 1)
            if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos='https://cloud.r-project.org/')
            install.packages("GEOquery", repos='https://cloud.r-project.org/')
            install.packages("R.utils")
            library(BiocManager)
            library(GEOquery)
            library(Biobase)
            get_betamatrix <- function(r, verbose=FALSE) {
            GEO_real_world_data <- "GSE67170"
            gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE)
            beta_matrix <- as.data.frame(exprs(gset[[1]]))
            control_samples <- c('GSM1641098', 'GSM1641099', 'GSM1641100', 'GSM1641101', 'GSM1641102', 'GSM1641103', 'GSM1641104', 'GSM1641105', 'GSM1641106', 'GSM1641107', 'GSM1641108', 'GSM1641109', 'GSM1641110', 'GSM1641111', 'GSM1641112', 'GSM1641113', 'GSM1641114', 'GSM1641115','GSM1641116', 'GSM1641117')
            beta_matrix <- beta_matrix[,control_samples]
            sample_names <- colnames(beta_matrix)
            print(sample_names)
            cpg_names <- rownames(beta_matrix)
            list_df <- list(beta_matrix, sample_names, cpg_names)
            }
            get_betamatrix()
            ''')


def calculate_shape_parameters(beta_matrix):
    """
    Calculate the statistical mean and variance of a beta_matrix
    :param beta_matrix: a matrix of beta_values to use as a reference for generating shape parameters
    :rtype methPara: a dataframe containing a column for means and variances for each CpG in the beta_matrix
    """
    list_of_means = [row for row in range(0, len(beta_matrix.index))]
    list_of_vars = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_means:
        list_of_means[i] = np.nanmean(beta_matrix.iloc[i])
        list_of_vars[i] = np.nanvar(beta_matrix.iloc[i])
    methPara = pd.DataFrame(list(zip(list_of_means, list_of_vars)), index=beta_matrix.index, columns=["mu","var"])
    return methPara

def get_power_calculation(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims):
    """
    Emulation of the R package pwrEWAS (https://rdrr.io/bioc/pwrEWAS/) using rpy2, with the extension of selecting
    in which CpGs the signal is to be introduced based on the methylation state (e.g., <= 0.1, >= 0.9, <= 0.4 and >= 0.6)
    :param targetDmCpGs: Target number of DM CpGs.
    :param methPara: a shape parameter dataframe containing a column for means and variances
    :param detectionLimit: Smallest detectable difference in DNAm (default: 0.01).
    :param J: Number of CpGs tested/simulated (default: 100000).
    :param CpGonArray: a matrix of beta_values to use as a reference for generating shape parameters
    :param targetDelta: targetDelta Target maximum difference in mean DNAm. (Either 'targetDelta' or 'deltaSD' should be specified)
    :param DMmethod: Method of Differential Methylation analysis: "limma" (default), "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc".
    :param minTotSampleSize: Minimum total sample size.
    :param maxTotSampleSize: Maximum total sample size.
    :param SampleSizeSteps: Sample size increments.
    :param FDRcritVal: FDRcritVal (default: 0.05).
    :param NcntPer: Percentage sample group 1 (control group) (NcntPer = 0.5 indicates a balanced design).
    :param core: Number of threads for multi-threading (default: 1).
    :param sims: Number of simulated data sets (default: 50).
    :rtype synthPwr: a list containing: powerArray (power-sample calculations per sim), metrics (marTypeI, classicalPower, FDR, and FDC measures) and FDR values
    """
    robjects.r('''
                install.packages("doSNOW")
                install.packages("doParallel")
                library(doSNOW)
                library(doParallel)
                library(truncnorm)
                library(limma)
                library(CpGassoc)
                library(genefilter)
                
                getTau <- function(targetDmCpGs, targetDelta, methPara, detectionLimit, J, CpGonArray){
                out <- NULL
                tau <- 1
                tauSteps <- 1
                lookForTau <- TRUE
                cnt <- 0
                maxCnt <- 100
    
                while(cnt < maxCnt & lookForTau){
                percentile <- NULL
                for(i in seq_len(100)){
                cpgIdx4Tau <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE) 
                delta <- truncnorm::rtruncnorm(1, mean = 0, sd = tau, 
                a=0.5 - methPara$mu[cpgIdx4Tau] - 
                    sqrt(0.25-methPara$var[cpgIdx4Tau]), 
                b=0.5 - methPara$mu[cpgIdx4Tau] + 
                    sqrt(0.25-methPara$var[cpgIdx4Tau]))
                percentile[i] <- stats::quantile(abs(delta),0.9999, na.rm = TRUE)
                }
                targetDelta <- as.numeric(targetDelta)
                if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau >= 1){
                tau <- tau + 1
                } else if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau < 1){
                tauSteps <- 0.5 * tauSteps
                tau <- tau + tauSteps
                } else if(mean(percentile) > targetDelta + 0.5*detectionLimit){
                tauSteps <- 0.5 * tauSteps
                tau <- tau - tauSteps
                } else {
                lookForTau <- FALSE
                }
                cnt <- cnt + 1
                }
                if(cnt == maxCnt) stop("Max iterations reached")
                truelyDMperc <- mean(abs(delta) > detectionLimit)
                out$tau <- tau
                targetK <- round(1/truelyDMperc * targetDmCpGs)
                out$K <- ifelse(targetK > J, J, targetK)
                return(out)
                }

                getK <- function(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, tau){
                cpgIdx4Tau <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE) 
                delta <- truncnorm::rtruncnorm(1, mean = 0, sd = tau, 
                a=0.5 - methPara$mu[cpgIdx4Tau] - sqrt(0.25 - methPara$var[cpgIdx4Tau]), 
                b=0.5 - methPara$mu[cpgIdx4Tau] + sqrt(0.25 - methPara$var[cpgIdx4Tau]))
                truelyDMperc <- mean(abs(delta) > detectionLimit)
                targetK <- round(1/truelyDMperc * targetDmCpGs)
                K <- ifelse(targetK > J, J, targetK)
                return(K)
                }
                
                getAlphBet <- function(myMean, myVar){
                alpha <- myMean^2 * ((1-myMean)/myVar - 1/myMean) 
                beta <- alpha * (1/myMean - 1)
                return( list(alpha = alpha, beta = beta) )
                }

                getMeanVar <- function(myAlpha, myBeta){
                mean <- myAlpha / (myAlpha + myBeta)
                var <- myAlpha * myBeta / ((myAlpha + myBeta)^2 * (myAlpha + myBeta + 1))
                return( list(mean = mean, var = var) )
                }

                beta2Mvalue <- function(beta){ # beta to m-value
                return( log2(beta / (1-beta)) )
                }

                ttestSlow <- function(g1Beta,g2Beta,rCnt,rTx,paired){
                mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
                ttest <- apply(mvals,1,  function (x) 
                stats::t.test(x[seq_len(rCnt)], x[(rCnt+1):(rCnt+rTx)], paired = FALSE, var.equal = TRUE))
                temp <- NULL
                temp$pval <- unlist(lapply(ttest, function(x) x$p.value))
                temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
                return(temp)
                }

                ttestFast <- function(g1Beta,g2Beta,rCnt,rTx){
                mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
                ttest <- genefilter::rowttests(mvals, 
                fac = factor(c(rep("g1",rCnt),rep("g2",rTx))), tstatOnly = FALSE)
                temp <- NULL
                temp$pval <- ttest$p.value
                temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
                return(temp)
                }

                Wilcox <- function(g1Beta, g2Beta, rCnt,rTx){
                mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
                WRS <- apply(mvals,1, function (x) stats::wilcox.test(x[seq_len(rCnt)] - x[(rCnt+1):(rCnt+rTx)], correct=TRUE))
                temp <- NULL
                temp$pval <- unlist(lapply(WRS, function(x) x$p.value))
                temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
                return(temp)
                }

                limma <- function(g1Beta, g2Beta, rCnt,rTx){
                mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
                design <- stats::model.matrix(~ c(rep("g1",rCnt),rep("g2",rTx)))
                limmaFit <- limma::lmFit(mvals, design)
                temp <- NULL
                temp$pval <- eBayes(limmaFit)$p.value[,2]
                temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
                return(temp)
                }

                CPGassoc <- function(g1Beta, g2Beta, rCnt,rTx){
                assoc <- CpGassoc::cpg.assoc(cbind(g1Beta, g2Beta), c(rep("g1",rCnt),rep("g2",rTx)))
                temp <- NULL
                temp$pval <- assoc$results$P.value
                temp$fdr <- assoc$results$FDR
                return(temp)
                }
                
                combine_tau <- function(listA, listB){
                if(is.null(listA)) return(listB)
                if(!methods::is(listA[["power"]], "array") & !methods::is(listA[["power"]], "matrix")) listA[["power"]] <- matrix(listA[["power"]])
                if(!methods::is(listB[["power"]], "array") & !methods::is(listB[["power"]], "matrix")) listB[["power"]] <- matrix(listB[["power"]])
                if(!methods::is(listA[["metric"]]$marTypeI, "array") & !methods::is(listA[["metric"]]$marTypeI, "matrix")) listA[["metric"]]$marTypeI <- matrix(listA[["metric"]]$marTypeI)
                if(!methods::is(listB[["metric"]]$marTypeI, "array") & !methods::is(listB[["metric"]]$marTypeI, "matrix")) listB[["metric"]]$marTypeI <- matrix(listB[["metric"]]$marTypeI)
                if(!methods::is(listA[["metric"]]$classicalPower, "array") & !methods::is(listA[["metric"]]$classicalPower, "matrix")) listA[["metric"]]$classicalPower <- matrix(listA[["metric"]]$classicalPower)
                if(!methods::is(listB[["metric"]]$classicalPower, "array") & !methods::is(listB[["metric"]]$classicalPower, "matrix")) listB[["metric"]]$classicalPower <- matrix(listB[["metric"]]$classicalPower)
                if(!methods::is(listA[["metric"]]$FDR, "array") & !methods::is(listA[["metric"]]$FDR, "matrix")) listA[["metric"]]$FDR <- matrix(listA[["metric"]]$FDR)
                if(!methods::is(listB[["metric"]]$FDR, "array") & !methods::is(listB[["metric"]]$FDR, "matrix")) listB[["metric"]]$FDR <- matrix(listB[["metric"]]$FDR)
                if(!methods::is(listA[["metric"]]$FDC, "array") & !methods::is(listA[["metric"]]$FDC, "matrix")) listA[["metric"]]$FDC <- matrix(listA[["metric"]]$FDC)
                if(!methods::is(listB[["metric"]]$FDC, "array") & !methods::is(listB[["metric"]]$FDC, "matrix")) listB[["metric"]]$FDC <- matrix(listB[["metric"]]$FDC)
                if(!methods::is(listA[["metric"]]$probTP, "array") & !methods::is(listA[["metric"]]$probTP, "matrix")) listA[["metric"]]$probTP <- matrix(listA[["metric"]]$probTP)
                if(!methods::is(listB[["metric"]]$probTP, "array") & !methods::is(listB[["metric"]]$probTP, "matrix")) listB[["metric"]]$probTP <- matrix(listB[["metric"]]$probTP)
                if(!methods::is(listA[["delta"]], "list")) listA[["delta"]] <- list(listA[["delta"]])
                returnList <- list()
                returnList[["power"]] <- abind::abind(listA[["power"]], listB[["power"]], along = 3)
                returnList[["delta"]] <- listA[["delta"]]
                returnList[["delta"]][[length(listA[["delta"]])+1]] <- listB[["delta"]]
                returnList[["metric"]]$marTypeI <- abind::abind(listA[["metric"]]$marTypeI, listB[["metric"]]$marTypeI, along = 3)
                returnList[["metric"]]$classicalPower <- abind::abind(listA[["metric"]]$classicalPower, listB[["metric"]]$classicalPower, along = 3)
                returnList[["metric"]]$FDR <- abind::abind(listA[["metric"]]$FDR, listB[["metric"]]$FDR, along = 3)
                returnList[["metric"]]$FDC <- abind::abind(listA[["metric"]]$FDC, listB[["metric"]]$FDC, along = 3)
                returnList[["metric"]]$probTP <- abind::abind(listA[["metric"]]$probTP, listB[["metric"]]$probTP, along = 3)
                return(returnList)
                }
    
                combine_totSampleSizes <- function(listA, listB){
                if(is.null(listA)) return(listB)
                returnList <- list()
                returnList[["power"]] <- cbind(listA[["power"]], listB[["power"]])
                returnList[["delta"]] <- cbind(listA[["delta"]], listB[["delta"]])
                returnList[["metric"]]$marTypeI <- cbind(listA[["metric"]]$marTypeI, listB[["metric"]]$marTypeI)
                returnList[["metric"]]$classicalPower <- cbind(listA[["metric"]]$classicalPower, listB[["metric"]]$classicalPower)
                returnList[["metric"]]$FDR <- cbind(listA[["metric"]]$FDR, listB[["metric"]]$FDR)
                returnList[["metric"]]$FDC <- cbind(listA[["metric"]]$FDC, listB[["metric"]]$FDC)
                returnList[["metric"]]$probTP <- cbind(listA[["metric"]]$probTP, listB[["metric"]]$probTP)
                return(returnList)
                }
                
                power_calc <- function(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims, verbose=FALSE) {
                output <- NULL
                totSampleSizes <- seq(minTotSampleSize, maxTotSampleSize, SampleSizeSteps)
                deltaSD <- NULL
                if(is.null(deltaSD)){
                K <- NULL 
                tau <- NULL
                for(d in seq_along(targetDelta)){
                myTau <- getTau(targetDmCpGs, targetDelta[d], methPara, detectionLimit, J, CpGonArray)
                tau[d] <- myTau$tau
                K[d] <- myTau$K 
                }
                print("the following Taus were chosen:")
                print(tau)
                write.table(as.data.frame(tau),file=paste('/Users/malwash/PycharmProjects/Synthetic_Power/summary_stats/taus.csv'),sep=",",row.names=F)
                } else {
                tau <- deltaSD 
                K <- NULL
                for(d in seq_along(tau)){
                K[d] <- getK(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, tau)
                }
                }
                cl <- parallel::makeCluster(core, outfile="")
                doSNOW::registerDoSNOW(cl)
                Ntot <- NULL
                
                multiThreadOut <- foreach(d = seq_along(tau), .combine = combine_tau, .packages = c("truncnorm", "limma", "CpGassoc", "genefilter"), .export = c("getAlphBet", "getMeanVar", "beta2Mvalue", "limma", "ttestSlow", "ttestFast", "Wilcox", "CPGassoc")) %:%
                foreach(Ntot = totSampleSizes, .combine = combine_totSampleSizes) %dopar% { 
                Ncnt <- round(Ntot * NcntPer)
                Ntx <- Ntot - Ncnt
                marPower <- NULL
                deltaSim <- NULL
                marTypeI <- NULL
                FDR <- NULL
                classicalPower <- NULL
                FDC <- NULL
                probTP <- NULL
                for(sim in seq_len(sims)){
                cpgIdx <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE)
                cpgIdxName <- paste(seq_len(J), "_", rownames(methPara)[cpgIdx], sep = "")
                changedCpgsIdx <- sample(x = cpgIdx, size = K[d])
                changedCpgsIdxName <- cpgIdxName[match(changedCpgsIdx, cpgIdx)]
                
                delta <- truncnorm::rtruncnorm(1, mean = 0, sd = as.numeric(tau[d]), 
                a=0.5 - methPara$mu[changedCpgsIdx] - sqrt(0.25-methPara$var[changedCpgsIdx]), 
                b=0.5 - methPara$mu[changedCpgsIdx] + sqrt(0.25-methPara$var[changedCpgsIdx]))
                deltaSim <- c(deltaSim, delta)
                meaningfulDM <- (abs(delta) >= detectionLimit)
                meaningfulDMName <- changedCpgsIdxName[meaningfulDM]
                
                muToBeSimuUNchanged <- methPara$mu[cpgIdx]
                muToBeSimuChanged  <- methPara$mu[cpgIdx]
                
                site_mean_check <- muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] >= 0.4 & muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] <= 0.6
                muToBeSimuChanged[site_mean_check] <- muToBeSimuChanged[site_mean_check] + delta

                params_unchanged <- getAlphBet(myMean = muToBeSimuUNchanged, myVar = methPara$var[cpgIdx])
                alpha_unchanged <- params_unchanged$alpha
                beta_unchanged <- params_unchanged$beta
                params_changed <- getAlphBet(myMean = muToBeSimuChanged, myVar = methPara$var[cpgIdx])
                alpha_changed <- params_changed$alpha
                beta_changed <- params_changed$beta
                
                g1Beta <- NULL
                g2Beta <- NULL
                g1Beta <- matrix(stats::rbeta(J*Ncnt, rep(alpha_unchanged, each = Ncnt), rep(beta_unchanged, each = Ncnt)), ncol = Ncnt, byrow = TRUE) 
                g2Beta <- matrix(stats::rbeta(J*Ntx, rep(alpha_changed, each = Ntx), rep(beta_changed, each = Ntx)), ncol = Ntx, byrow = TRUE) 
                
                g1Beta[g1Beta == 1] <- max(g1Beta[g1Beta != 1])
                g2Beta[g2Beta == 1] <- max(g2Beta[g2Beta != 1])
                g1Beta[g1Beta == 0] <- min(g1Beta[g1Beta != 0])
                g2Beta[g2Beta == 0] <- min(g2Beta[g2Beta != 0])
                rownames(g1Beta) <- rownames(g2Beta) <- paste(seq_len(J),"_",names(alpha_unchanged),sep = "")
                
                DMtest <- switch(DMmethod,
                    "t-test (unequal var)" = ttestSlow(g1Beta,g2Beta,Ncnt,Ntx,paired=FALSE),
                    "t-test (equal var)" = ttestFast(g1Beta,g2Beta,Ncnt,Ntx), 
                    "CPGassoc" = CPGassoc(g1Beta, g2Beta, Ncnt,Ntx), 
                    "Wilcox rank sum" = Wilcox(g1Beta, g2Beta, Ncnt,Ntx),
                    "limma" = limma(g1Beta, g2Beta, Ncnt,Ntx),
                    stop("Test not found")
                )

                notDM  <- cpgIdxName[!(cpgIdxName %in% changedCpgsIdxName)]
                DM_negligible <- changedCpgsIdxName[!(changedCpgsIdxName %in% meaningfulDMName)]
                DM_meaningful <- changedCpgsIdxName[changedCpgsIdxName %in% meaningfulDMName]
                FP  <- intersect(cpgIdxName[DMtest$fdr<FDRcritVal], notDM)
                NP <- intersect(cpgIdxName[DMtest$fdr<FDRcritVal], DM_negligible)
                TP <- intersect(cpgIdxName[DMtest$fdr<FDRcritVal], DM_meaningful)
                detectedCpGs <- cpgIdxName[DMtest$fdr<FDRcritVal]
                TN <- intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], notDM)
                NN <- intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], DM_negligible)
                FN <- intersect(cpgIdxName[!(DMtest$fdr<FDRcritVal)], DM_meaningful)
                
                marPower[sim] <- ifelse(length(DM_meaningful) > 0, length(TP)/length(DM_meaningful), NA)
                marTypeI[sim] <- ifelse(length(notDM) > 0, length(FP)/length(notDM), NA) 
                FDR[sim] <- ifelse(length(detectedCpGs) > 0, length(FP)/length(detectedCpGs), NA)
                FDC[sim] <- ifelse(length(TP) > 0, (length(FP))/length(TP), NA)
                classicalPower[sim] <- (length(NP)+length(TP))/(length(DM_negligible)+length(DM_meaningful))
                probTP[sim] <- ifelse(length(TP)>0, 1, 0)
                
                variance_g1Beta_column <- apply(g1Beta, 2, var)
                variance_g2Beta_column <- apply(g2Beta, 2, var)
                mean_g1Beta_row <- rowMeans(g1Beta)
                mean_g2Beta_row <- rowMeans(g2Beta)
                write.table(variance_g1Beta_column,file=paste('/Users/malwash/PycharmProjects/Synthetic_Power/Output/variances/varg1col','_tau',as.character(tau[d]),'_sim',sim,'_samplesize',Ntx,'.csv'),sep=",",  col.names=FALSE)
                write.table(variance_g2Beta_column,file=paste('/Users/malwash/PycharmProjects/Synthetic_Power/Output/variances/varg2col','_tau',as.character(tau[d]),'_sim',sim,'_samplesize',Ntx,'.csv'),sep=",",  col.names=FALSE)
                write.table(mean_g1Beta_row,file=paste('/Users/malwash/PycharmProjects/Synthetic_Power/Output/means/meang1row','_tau',as.character(tau[d]),'_sim',sim,'_samplesize',Ntx,'.csv'),sep=",",  col.names=FALSE)
                write.table(mean_g2Beta_row,file=paste('/Users/malwash/PycharmProjects/Synthetic_Power/Output/means/meang2row','_tau',as.character(tau[d]),'_sim',sim,'_samplesize',Ntx,'.csv'),sep=",",  col.names=FALSE)
                }

                outSim <- list() 
                outSim[["power"]] <- marPower 
                outSim[["delta"]] <- deltaSim
                outSim[["metric"]]$marTypeI <- marTypeI
                outSim[["metric"]]$FDR <- FDR
                outSim[["metric"]]$classicalPower <- classicalPower
                outSim[["metric"]]$FDC <- FDC
                outSim[["metric"]]$probTP <- probTP
                outSim
                } 
                parallel::stopCluster(cl)
                if(is.null(targetDelta)) targetDelta <- tau
    
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$meanPower <- matrix(mean(multiThreadOut[["power"]], na.rm=TRUE))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$meanPower <- matrix(apply(multiThreadOut[["power"]], 2, mean, na.rm=TRUE), ncol = 1)
                if(length(targetDelta) > 1)                                 output$meanPower <- apply(multiThreadOut[["power"]], c(2,3), mean, na.rm=TRUE)
                rownames(output$meanPower) <- totSampleSizes
                colnames(output$meanPower) <- targetDelta
    
                if(length(targetDelta) == 1) output$powerArray <- array(data = multiThreadOut[["power"]], dim = c(sims, length(totSampleSizes), length(targetDelta)))
                if(length(targetDelta) > 1 & length(totSampleSizes) == 1) output$powerArray <- multiThreadOut[["power"]]
                if(length(targetDelta) > 1) output$powerArray <- multiThreadOut[["power"]]
                dimnames(output$powerArray) <- list(seq_len(sims), totSampleSizes, targetDelta)  
                
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$deltaArray <- list(matrix(multiThreadOut[["delta"]]))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$deltaArray <- list(multiThreadOut[["delta"]])
                if(length(targetDelta) > 1 & length(totSampleSizes) == 1)   output$deltaArray <- lapply(multiThreadOut[["delta"]],as.matrix)
                if(length(targetDelta) > 1 & length(totSampleSizes) > 1)    output$deltaArray <- multiThreadOut[["delta"]]
                names(output$deltaArray) <- targetDelta
                for(d in seq_along(targetDelta)){
                colnames(output$deltaArray[[d]]) <- totSampleSizes
                }
    
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$marTypeI <- matrix(mean(multiThreadOut[["metric"]]$marTypeI, na.rm=TRUE))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$marTypeI <- matrix(apply(multiThreadOut[["metric"]]$marTypeI, 2, mean, na.rm=TRUE), ncol = 1)
                if(length(targetDelta) > 1)                                 output$metric$marTypeI <- apply(multiThreadOut[["metric"]]$marTypeI, c(2,3), mean, na.rm=TRUE)
                rownames(output$metric$marTypeI) <- totSampleSizes
                colnames(output$metric$marTypeI) <- targetDelta
    
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$classicalPower <- matrix(mean(multiThreadOut[["metric"]]$classicalPower, na.rm=TRUE))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$classicalPower <- matrix(apply(multiThreadOut[["metric"]]$classicalPower, 2, mean, na.rm=TRUE), ncol = 1)
                if(length(targetDelta) > 1)                                 output$metric$classicalPower <- apply(multiThreadOut[["metric"]]$classicalPower, c(2,3), mean, na.rm=TRUE)
                rownames(output$metric$classicalPower) <- totSampleSizes
                colnames(output$metric$classicalPower) <- targetDelta
    
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$FDR <- matrix(mean(multiThreadOut[["metric"]]$FDR, na.rm=TRUE))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$FDR <- matrix(apply(multiThreadOut[["metric"]]$FDR, 2, mean, na.rm=TRUE), ncol = 1)
                if(length(targetDelta) > 1)                                 output$metric$FDR <- apply(multiThreadOut[["metric"]]$FDR, c(2,3), mean, na.rm=TRUE)
                rownames(output$metric$FDR) <- totSampleSizes
                colnames(output$metric$FDR) <- targetDelta
    
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$FDC <- matrix(mean(multiThreadOut[["metric"]]$FDC, na.rm=TRUE))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$FDC <- matrix(apply(multiThreadOut[["metric"]]$FDC, 2, mean, na.rm=TRUE), ncol = 1)
                if(length(targetDelta) > 1)                                 output$metric$FDC <- apply(multiThreadOut[["metric"]]$FDC, c(2,3), mean, na.rm=TRUE)
                rownames(output$metric$FDC) <- totSampleSizes
                colnames(output$metric$FDC) <- targetDelta
    
                if(length(targetDelta) == 1 & length(totSampleSizes) == 1)  output$metric$probTP <- matrix(mean(multiThreadOut[["metric"]]$probTP, na.rm=TRUE))
                if(length(targetDelta) == 1 & length(totSampleSizes) > 1)   output$metric$probTP <- matrix(apply(multiThreadOut[["metric"]]$probTP, 2, mean, na.rm=TRUE), ncol = 1)
                if(length(targetDelta) > 1)                                 output$metric$probTP <- apply(multiThreadOut[["metric"]]$probTP, c(2,3), mean, na.rm=TRUE)
                rownames(output$metric$probTP) <- totSampleSizes
                colnames(output$metric$probTP) <- targetDelta
                return(list(output$powerArray, output$metric, multiThreadOut[["metric"]]$FDR))
                }
                ''')
    power_result = robjects.r['power_calc']
    output_powerarray, output_metrics, output_fdr = power_result(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims)
    sample_steps = np.arange(minTotSampleSize, maxTotSampleSize + SampleSizeSteps, SampleSizeSteps)

    for idx, metric in enumerate(output_metrics):
        if idx == 0:
            marTypeI = pd.DataFrame(metric, columns=targetDelta)
            marTypeI.insert(len(marTypeI.columns), column="Sample_Size", value=sample_steps)
            marTypeI.to_csv(os.getcwd() + figuredirname + 'marTypeI.csv')
        elif idx == 1:
            classicalPower = pd.DataFrame(metric, columns=targetDelta)
            classicalPower.insert(len(classicalPower.columns), column="Sample_Size", value=sample_steps)
            classicalPower.to_csv(os.getcwd() + figuredirname + 'classicalPower.csv')
        elif idx == 2:
            fdr = pd.DataFrame(metric, columns=targetDelta)
            fdr.insert(len(fdr.columns), column="Sample_Size", value=sample_steps)
            fdr.to_csv(os.getcwd() + figuredirname + 'fdr.csv')
        elif idx == 3:
            fdc = pd.DataFrame(metric, columns=targetDelta)
            fdc.insert(len(fdc.columns), column="Sample_Size", value=sample_steps)
            fdc.to_csv(os.getcwd() + figuredirname + 'fdc.csv')
        elif idx == 4:
            probTP = pd.DataFrame(metric, columns=targetDelta)
            probTP.insert(len(probTP.columns), column="Sample_Size", value=sample_steps)
            probTP.to_csv(os.getcwd() + figuredirname + 'probTP.csv')

    list_of_sim_results = []
    for idx, sim in enumerate(output_powerarray):
        sim_run = pd.DataFrame(columns=targetDelta)
        for sampleSize in sim:
            sim_run.loc[len(sim_run)] = sampleSize.tolist()
        sim_run.insert(len(sim_run.columns), column="Sample_Size", value=sample_steps)
        list_of_sim_results.append(sim_run)
    plot_power_simulations_by_sample(list_of_sim_results, targetDelta, sample_steps, classicalPower)

    taus = pd.read_csv('/Users/malwash/PycharmProjects/Synthetic_Power/summary_stats/taus.csv', header=0)
    g1col = {row['tau']: [] for index, row in taus.iterrows()}
    g2col = {row['tau']: [] for index, row in taus.iterrows()}

    for sim_iter in range(1,sims+1):
        for step in sample_steps:
            for index, row in taus.iterrows():
                g1colpd = pd.read_csv('/Users/malwash/PycharmProjects/Synthetic_Power/Output/variances/varg1col' + ' _tau ' + str(row['tau']) + ' _sim ' + str(sim_iter) + ' _samplesize ' + str(int(step / 2)) + ' .csv',index_col=0, header=None)
                g2colpd = pd.read_csv('/Users/malwash/PycharmProjects/Synthetic_Power/Output/variances/varg2col' + ' _tau ' + str(row['tau']) + ' _sim ' + str(sim_iter) + ' _samplesize ' + str(int(step / 2)) + ' .csv',index_col=0, header=None)
                g1col[row['tau']].extend(g1colpd.iloc[:,0].values.tolist())
                g2col[row['tau']].extend(g2colpd.iloc[:,0].values.tolist())

    plt.clf()
    plt.figure(figsize=(10, 10))
    plt.title("Distribution of variance across sample columns")
    plt.xlabel('Variance (σ2)')
    plt.ylabel('Frequency')

    for index, row in taus.iterrows():
        sns.kdeplot(g1col[row['tau']], label='Case - τ:'+str(row['tau']))
        sns.kdeplot(g2col[row['tau']], label='Control - τ:'+str(row['tau']))
    plt.legend()
    plt.savefig(os.getcwd() + figuredirname + "dist_samples.png",dpi=300)
    plt.show()

    plt.clf()
    plt.figure(figsize=(10, 10))
    plt.title("Distribution of CpG means")
    plt.xlabel('Average % of Methylation values')
    plt.ylabel('Frequency')

    g1row = {row['tau']: [] for index, row in taus.iterrows()}
    g2row = {row['tau']: [] for index, row in taus.iterrows()}
    for sim_iter in range(1,sims+1):
        for step in sample_steps:
            for index, row in taus.iterrows():
                g1rowpd = pd.read_csv('/Users/malwash/PycharmProjects/Synthetic_Power/Output/means/meang1row' + ' _tau ' + str(row['tau']) + ' _sim ' + str(sim_iter) + ' _samplesize ' + str(int(step / 2)) + ' .csv',index_col=0, header=None)
                g2rowpd = pd.read_csv('/Users/malwash/PycharmProjects/Synthetic_Power/Output/means/meang2row' + ' _tau ' + str(row['tau']) + ' _sim ' + str(sim_iter) + ' _samplesize ' + str(int(step / 2)) + ' .csv',index_col=0, header=None)
                g1row[row['tau']].extend(g1rowpd.iloc[:,0].values.tolist())
                g2row[row['tau']].extend(g2rowpd.iloc[:,0].values.tolist())

    for index, row in taus.iterrows():
        sns.kdeplot(g1row[row['tau']], label='Case - τ:' + str(row['tau']),color='green',fill=True, common_norm=False, palette="crest",alpha=.5, linewidth=0)
        sns.kdeplot(g2row[row['tau']], label='Control - τ:' + str(row['tau']),color='blue',fill=True, common_norm=False, palette="crest",alpha=.5, linewidth=0)
    plt.legend()
    plt.savefig(os.getcwd() + figuredirname + "dist_cpgs.png", dpi=300)
    plt.show()

def synthPwr(minTotSampleSize,maxTotSampleSize,SampleSizeSteps,NcntPer,targetDelta,deltaSD=None,J=100000,targetDmCpGs=100,tissueType="GSE67170",detectionLimit=0.01,DMmethod=["limma", "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc"],FDRcritVal=0.05,core=4,sims=50):
    """
    Main application to retrieve a reference EWAS and use this as a basis to generate a simulation and run synthPwr
    :param minTotSampleSize: Minimum total sample size.
    :param maxTotSampleSize: Maximum total sample size.
    :param SampleSizeSteps: Sample size increments.
    :param NcntPer: Percentage sample group 1 (control group) (NcntPer = 0.5 indicates a balanced design).
    :param targetDelta: targetDelta Target maximum difference in mean DNAm. (Either 'targetDelta' or 'deltaSD' should be specified)
    :param deltaSD: Standard deviation of simulated differences. (Either 'targetDelta' or 'deltaSD' should be specified)
    :param J: Number of CpGs tested/simulated (default: 100000).
    :param targetDmCpGs: Target number of DM CpGs.
    :param tissueType: Select a tissue from a GEO Accession number.
    :param methPara: a shape parameter dataframe containing a column for means and variances
    :param detectionLimit: Smallest detectable difference in DNAm (default: 0.01).
    :param DMmethod: Method of Differential Methylation analysis: "limma" (default), "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc".
    :param FDRcritVal: FDRcritVal (default: 0.05).
    :param core: Number of threads for multi-threading (default: 1).
    :param sims: Number of simulated data sets (default: 50).
    :rtype synthPwr: a list containing: powerArray (power-sample calculations per sim), metrics (marTypeI, classicalPower, FDR, and FDC measures) and FDR values
    """
    beta_matrix_pull = robjects.r['get_betamatrix']
    beta_matrix = pd.DataFrame(beta_matrix_pull()[0]).transpose()
    print("The reference data - Beta matrix:")
    print(beta_matrix.shape)
    print(beta_matrix)
    methPara = calculate_shape_parameters(beta_matrix)
    CpGonArray = len(methPara['mu'])
    get_power_calculation(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims)

###################################################
# Choose parameterisation to run synthPwr
###################################################

if __name__ == '__main__':
    input = {}
    input['Nmin'] = 50
    input['Nmax'] = 950
    input['NCntPer'] = 0.5
    input['Nsteps'] = 150
    input['J'] = 10000
    input['targetDmCpGs'] = 10
    input['targetDeltaString'] = "0.01, 0.1, 0.5"
    input['tauString'] = "0.01, 0.03"
    input['targetDelta'] = list(map(float, re.split(', ', input['targetDeltaString'])))
    input['deltaSDString'] = "0.01, 0.03"
    input['deltaSD'] = list(map(float, re.split(', ', input['deltaSDString'])))
    input['method'] = "limma"
    input['detectionLimit'] = 0.01
    input['FDRcritVal'] = 0.05
    input['cores'] = 4
    input['sim'] = 50
    input['tissueType'] = "Saliva"
    synthPwr(minTotSampleSize = input['Nmin'],maxTotSampleSize = input['Nmax'],SampleSizeSteps = input['Nsteps'],NcntPer = input['NCntPer'],targetDelta = input['targetDelta'],J = input['J'],targetDmCpGs = input['targetDmCpGs'],tissueType = input['tissueType'],detectionLimit = input['detectionLimit'],DMmethod = input['method'],FDRcritVal = input['FDRcritVal'],core = input['cores'],sims = input['sim'])
    print("Time taken: ", time.time() - starttime)