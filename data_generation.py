import csv
import itertools
import math
import os
from zipfile import ZipFile
import numpy as np
import pandas as pd
import random
import time
import warnings
from multiprocessing import Pool
from itertools import product
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
from numpy import dtype
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import DataFrame, StrVector
from rpy2.robjects.packages import importr

#pd.options.display.max_rows
warnings.simplefilter(action='ignore', category=FutureWarning)

starttime = time.time()
dirname = os.sep+"power_experiments"+os.sep
figuredirname = os.sep+"figures"+os.sep
summarycalcdirname = os.sep+"summary_stats"+os.sep

#GSE67170 controls
#control_samples <- c('GSM1641098', 'GSM1641099', 'GSM1641100', 'GSM1641101', 'GSM1641102', 'GSM1641103', 'GSM1641104', 'GSM1641105', 'GSM1641106', 'GSM1641107', 'GSM1641108', 'GSM1641109', 'GSM1641110', 'GSM1641111', 'GSM1641112', 'GSM1641113', 'GSM1641114', 'GSM1641115','GSM1641116', 'GSM1641117')
#control_samples < - c('GSM1641098_9283265129_R01C01', 'GSM1641099_9283265129_R01C02', 'GSM1641100_9283265129_R02C01',
#                      'GSM1641101_9283265129_R02C02', 'GSM1641102_9283265129_R03C01', 'GSM1641103_9283265129_R03C02',
#                      'GSM1641104_9283265129_R04C01', 'GSM1641105_9283265129_R04C02', 'GSM1641106_9283265129_R05C01',
#                      'GSM1641107_9283265129_R06C01', 'GSM1641108_9283265129_R05C02', 'GSM1641109_9283265129_R06C02',
#                      'GSM1641110_9283265149_R01C01', 'GSM1641111_9283265149_R01C02', 'GSM1641112_9283265149_R02C01',
#                      'GSM1641113_9283265149_R02C02', 'GSM1641114_9283265149_R03C01', 'GSM1641115_9283265149_R04C01',
#                      'GSM1641116_9283265149_R05C01', 'GSM1641117_9283265149_R06C01')

#GSE42861 controls
#control_samples <- c('GSM1051533_7800246024_R03C02' , 'GSM1051534_7800246024_R04C02', 'GSM1051535_7800246024_R05C02' , 'GSM1051536_7800246024_R06C02', 'GSM1051537_7800246085_R01C01' , 'GSM1051538_7800246085_R02C01', 'GSM1051539_7800246085_R03C01' , 'GSM1051540_7800246085_R04C01', 'GSM1051541_7800246085_R05C01' , 'GSM1051542_7800246085_R06C01', 'GSM1051543_7800246085_R01C02' , 'GSM1051544_7800246085_R02C02', 'GSM1051545_7800246085_R03C02' , 'GSM1051549_7800246087_R01C01' , 'GSM1051550_7800246087_R02C01', 'GSM1051551_7800246087_R03C01' , 'GSM1051552_7800246087_R04C01', 'GSM1051553_7800246087_R05C01' , 'GSM1051555_7800246087_R01C02' , 'GSM1051556_7800246087_R02C02', 'GSM1051557_7800246087_R03C02' , 'GSM1051558_7800246087_R04C02', 'GSM1051559_7800246087_R05C02' , 'GSM1051560_7800246087_R06C02', 'GSM1051561_7800246123_R01C01' , 'GSM1051562_7800246123_R02C01', 'GSM1051563_7800246123_R03C01' , 'GSM1051564_7800246123_R04C01', 'GSM1051565_7800246123_R05C01' , 'GSM1051566_7800246123_R06C01', 'GSM1051567_7800246123_R01C02' , 'GSM1051568_7800246123_R02C02', 'GSM1051569_7800246123_R03C02' , 'GSM1051570_7800246123_R04C02', 'GSM1051571_7800246123_R05C02' , 'GSM1051572_7800246123_R06C02', 'GSM1051573_7800246132_R01C01' , 'GSM1051574_7800246132_R02C01', 'GSM1051575_7800246132_R04C01' , 'GSM1051576_7800246132_R05C01', 'GSM1051577_7800246132_R06C01' , 'GSM1051578_7800246132_R01C02', 'GSM1051579_7800246132_R02C02' , 'GSM1051580_7800246132_R03C02', 'GSM1051581_7800246132_R04C02' , 'GSM1051582_7800246132_R05C02', 'GSM1051583_7800246132_R06C02' , 'GSM1051584_7800246137_R01C01', 'GSM1051585_7800246137_R02C01' , 'GSM1051586_7800246137_R03C01', 'GSM1051587_7800246137_R04C01' , 'GSM1051588_7800246137_R05C01', 'GSM1051589_7800246137_R06C01' , 'GSM1051590_7800246137_R02C02', 'GSM1051591_7800246137_R03C02' , 'GSM1051592_7800246137_R04C02', 'GSM1051593_7800246137_R05C02' , 'GSM1051594_7800246137_R06C02', 'GSM1051595_7800246140_R01C01' , 'GSM1051596_7800246140_R02C01', 'GSM1051597_7800246140_R03C01' , 'GSM1051599_7800246140_R05C01' , 'GSM1051601_7800246140_R01C02' , 'GSM1051602_7800246140_R02C02', 'GSM1051603_7800246140_R03C02' , 'GSM1051604_7800246140_R04C02', 'GSM1051605_7800246140_R05C02' , 'GSM1051606_7800246140_R06C02', 'GSM1051607_7800246188_R01C01' , 'GSM1051608_7800246188_R02C01', 'GSM1051609_7800246188_R03C01' , 'GSM1051610_7800246188_R04C01', 'GSM1051618_7796814016_R01C01', 'GSM1051619_7796814016_R03C01' , 'GSM1051620_7796814016_R04C01', 'GSM1051621_7796814016_R05C01' , 'GSM1051622_7796814016_R06C01', 'GSM1051623_7796814016_R01C02' , 'GSM1051624_7796814016_R02C02', 'GSM1051625_7796814016_R03C02' , 'GSM1051626_7796814016_R04C02', 'GSM1051627_7796814016_R05C02' , 'GSM1051628_7796814016_R06C02', 'GSM1051629_7800246006_R01C01' , 'GSM1051630_7800246006_R02C01', 'GSM1051631_7800246006_R03C01' , 'GSM1051632_7800246006_R04C01', 'GSM1051633_7800246006_R05C01' , 'GSM1051634_7800246006_R06C01', 'GSM1051635_7800246006_R01C02' , 'GSM1051638_7800246006_R04C02', 'GSM1051639_7800246006_R05C02' , 'GSM1051641_7800246018_R01C01' , 'GSM1051642_7800246018_R02C01', 'GSM1051643_7800246018_R03C01' , 'GSM1051644_7800246018_R04C01', 'GSM1051645_7800246018_R05C01' , 'GSM1051646_7800246018_R06C01', 'GSM1051647_7800246018_R01C02' , 'GSM1051648_7800246018_R02C02', 'GSM1051649_7800246018_R03C02' , 'GSM1051650_7800246018_R04C02', 'GSM1051651_7800246018_R05C02' , 'GSM1051652_7800246018_R06C02', 'GSM1051653_7800246043_R01C01' , 'GSM1051654_7800246043_R02C01', 'GSM1051655_7800246043_R03C01' , 'GSM1051656_7800246043_R05C01', 'GSM1051657_7800246043_R06C01' , 'GSM1051658_7800246043_R01C02', 'GSM1051659_7800246043_R03C02' , 'GSM1051660_7800246043_R04C02', 'GSM1051661_7800246043_R05C02' , 'GSM1051662_7800246043_R06C02', 'GSM1051663_7800246054_R01C01' , 'GSM1051664_7800246054_R02C01', 'GSM1051665_7800246054_R03C01' , 'GSM1051666_7800246054_R04C01', 'GSM1051667_7800246054_R05C01' , 'GSM1051668_7800246054_R06C01', 'GSM1051669_7800246054_R01C02' , 'GSM1051671_7800246054_R03C02' , 'GSM1051672_7800246054_R04C02', 'GSM1051673_7800246054_R05C02' , 'GSM1051674_7800246054_R06C02', 'GSM1051675_7800246055_R01C01' , 'GSM1051687_7800246071_R01C01' , 'GSM1051688_7800246071_R02C01', 'GSM1051689_7800246071_R03C01' , 'GSM1051690_7800246071_R04C01', 'GSM1051691_7800246071_R05C01' , 'GSM1051692_7800246071_R06C01', 'GSM1051693_7800246071_R01C02' , 'GSM1051710_7800246034_R01C01', 'GSM1051712_7800246034_R03C01', 'GSM1051713_7800246034_R04C01' , 'GSM1051714_7800246034_R05C01', 'GSM1051715_7800246034_R06C01' , 'GSM1051716_7800246034_R01C02', 'GSM1051717_7800246034_R02C02' , 'GSM1051726_7800246040_R05C01', 'GSM1051727_7800246040_R06C01' , 'GSM1051730_7800246040_R03C02', 'GSM1051731_7800246040_R04C02' , 'GSM1051732_7800246040_R05C02', 'GSM1051733_7800246040_R06C02' , 'GSM1051742_7800246041_R03C02', 'GSM1051744_7800246041_R05C02', 'GSM1051745_7800246041_R06C02' , 'GSM1051746_7800246044_R01C01', 'GSM1051749_7800246044_R04C01' , 'GSM1051751_7800246044_R06C01' , 'GSM1051752_7800246044_R01C02', 'GSM1051754_7800246044_R03C02', 'GSM1051755_7800246044_R04C02' , 'GSM1051756_7800246044_R05C02', 'GSM1051757_7800246046_R01C01' , 'GSM1051758_7800246046_R02C01', 'GSM1051759_7800246046_R03C01' , 'GSM1051760_7800246046_R04C01', 'GSM1051761_7800246046_R05C01' , 'GSM1051762_7800246046_R06C01', 'GSM1051763_7800246046_R01C02' , 'GSM1051764_7800246046_R02C02', 'GSM1051765_7800246046_R03C02' , 'GSM1051766_7800246046_R04C02', 'GSM1051767_7800246046_R05C02' , 'GSM1051768_7800246046_R06C02', 'GSM1051769_7800246057_R01C01' , 'GSM1051770_7800246057_R02C01', 'GSM1051771_7800246057_R03C01' , 'GSM1051772_7800246057_R04C01', 'GSM1051773_7800246057_R05C01' , 'GSM1051774_7800246057_R06C01', 'GSM1051775_7800246057_R01C02' , 'GSM1051777_7800246057_R03C02' , 'GSM1051781_7800246068_R01C01' , 'GSM1051782_7800246068_R02C01', 'GSM1051784_7800246068_R04C01', 'GSM1051786_7800246068_R06C01', 'GSM1051787_7800246068_R03C02' , 'GSM1051788_7800246068_R04C02', 'GSM1051789_7800246068_R05C02' , 'GSM1051790_7800246068_R06C02', 'GSM1051791_7800246079_R01C01' , 'GSM1051792_7800246079_R02C01', 'GSM1051793_7800246079_R04C01' , 'GSM1051802_7512560084_R05C02', 'GSM1051803_7512560084_R06C02' , 'GSM1051804_7512560103_R01C01', 'GSM1051805_7512560103_R02C01' , 'GSM1051806_7512560103_R03C01', 'GSM1051807_7512560103_R04C01' , 'GSM1051808_7512560103_R05C01', 'GSM1051809_7512560103_R06C01' , 'GSM1051810_7512560103_R01C02', 'GSM1051811_7512560103_R02C02' , 'GSM1051812_7512560103_R03C02', 'GSM1051813_7512560103_R04C02' , 'GSM1051814_7512560103_R05C02', 'GSM1051815_7512560103_R06C02' , 'GSM1051816_7512560104_R01C01', 'GSM1051817_7512560104_R02C01' , 'GSM1051820_7512560104_R05C01', 'GSM1051821_7512560104_R06C01' , 'GSM1051822_7512560104_R01C02', 'GSM1051823_7512560104_R02C02' , 'GSM1051824_7512560104_R03C02', 'GSM1051825_7512560104_R05C02' , 'GSM1051826_7512560104_R06C02', 'GSM1051827_7512560115_R01C01' , 'GSM1051828_7512560115_R02C01', 'GSM1051831_7512560115_R05C01' , 'GSM1051832_7512560115_R06C01', 'GSM1051833_7512560115_R01C02' , 'GSM1051834_7512560115_R02C02', 'GSM1051835_7512560115_R03C02' , 'GSM1051836_7512560115_R04C02', 'GSM1051837_7512560115_R05C02' , 'GSM1051838_7512560115_R06C02', 'GSM1051839_7512560124_R01C01' , 'GSM1051840_7512560124_R02C01', 'GSM1051841_7512560124_R03C01' , 'GSM1051842_7512560124_R04C01', 'GSM1051843_7512560124_R05C01' , 'GSM1051844_7512560124_R06C01', 'GSM1051845_7512560124_R01C02' , 'GSM1051846_7512560124_R02C02', 'GSM1051847_7512560124_R03C02' , 'GSM1051848_7512560124_R04C02', 'GSM1051850_7512560124_R06C02', 'GSM1051851_7512560128_R01C01' , 'GSM1051852_7512560128_R02C01', 'GSM1051853_7512560128_R04C01' , 'GSM1051854_7512560128_R05C01', 'GSM1051855_7512560128_R06C01' , 'GSM1051857_7512560128_R02C02' , 'GSM1051859_7512560128_R04C02' , 'GSM1051860_7512560128_R05C02', 'GSM1051862_7766130158_R01C01', 'GSM1051863_7766130158_R02C01' , 'GSM1051864_7766130158_R03C01', 'GSM1051865_7766130158_R04C01', 'GSM1051873_7766130166_R01C01' , 'GSM1051874_7766130166_R02C01', 'GSM1051876_7766130166_R04C01', 'GSM1051883_7766130166_R05C02' , 'GSM1051884_7766130166_R06C02', 'GSM1051901_5730192053_R01C01' , 'GSM1051902_5730192053_R02C01', 'GSM1051903_5730192053_R03C01' , 'GSM1051904_5730192053_R05C01', 'GSM1051905_5730192053_R06C01' , 'GSM1051906_5730192053_R01C02', 'GSM1051907_5730192053_R02C02' , 'GSM1051923_5730504014_R03C02' , 'GSM1051924_5730504014_R05C02', 'GSM1051925_5730504014_R06C02' , 'GSM1051926_5730504015_R01C01', 'GSM1051927_5730504015_R02C01' , 'GSM1051928_5730504015_R03C01', 'GSM1051929_5730504015_R04C01' , 'GSM1051949_5765205058_R01C01' , 'GSM1051950_5765205058_R02C01', 'GSM1051951_5765205058_R04C01' , 'GSM1051952_5765205058_R05C01', 'GSM1051953_5765205058_R06C01' , 'GSM1051954_5765205058_R01C02', 'GSM1051955_5765205058_R02C02' , 'GSM1051962_5765205059_R05C01', 'GSM1051963_5765205059_R06C01' , 'GSM1051964_5765205059_R01C02', 'GSM1051965_5765205059_R02C02' , 'GSM1051966_5765205059_R04C02', 'GSM1051967_5765205059_R05C02' , 'GSM1051968_5765205059_R06C02', 'GSM1052015_5730053006_R03C02', 'GSM1052016_5730053006_R04C02' , 'GSM1052017_5730053006_R05C02', 'GSM1052018_5730053006_R06C02' , 'GSM1052019_5730053010_R01C01', 'GSM1052020_5730053010_R02C01' , 'GSM1052021_5730053010_R03C01', 'GSM1052022_5730053010_R04C01' , 'GSM1052034_5730053011_R03C02' , 'GSM1052035_5730053011_R04C02', 'GSM1052036_5730053011_R05C02' , 'GSM1052037_5730053011_R06C02', 'GSM1052038_5730053027_R01C01' , 'GSM1052039_5730053027_R02C01', 'GSM1052040_5730053027_R03C01' , 'GSM1052041_5730053027_R04C01', 'GSM1052047_5730053037_R01C01', 'GSM1052048_5730053037_R04C01' , 'GSM1052049_5730053037_R06C01', 'GSM1052050_5730053037_R01C02' , 'GSM1052051_5730053037_R02C02', 'GSM1052099_5730192036_R05C01', 'GSM1052100_5730192036_R06C01' , 'GSM1052101_5730192036_R01C02', 'GSM1052102_5730192036_R02C02' , 'GSM1052103_5730192036_R03C02', 'GSM1052104_5730192036_R04C02' , 'GSM1052122_5730192048_R01C01' , 'GSM1052123_5730192048_R02C01', 'GSM1052124_5730192048_R03C01' , 'GSM1052125_5730192048_R04C01', 'GSM1052126_5730192048_R05C01' , 'GSM1052127_5730192048_R06C01', 'GSM1052128_5730192048_R01C02' , 'GSM1052129_5730192048_R02C02', 'GSM1052130_5730192048_R03C02' , 'GSM1052131_5730192048_R04C02', 'GSM1052132_5730192048_R05C02' , 'GSM1052133_5730192048_R06C02', 'GSM1052134_5730192057_R01C01' , 'GSM1052135_5730192057_R02C01', 'GSM1052136_5730192057_R03C01' , 'GSM1052137_5730192057_R04C01', 'GSM1052188_5730053041_R05C01' , 'GSM1052189_5730053041_R06C01', 'GSM1052190_5730053041_R01C02' , 'GSM1052191_5730053041_R02C02', 'GSM1052192_5730053041_R03C02' , 'GSM1052193_5730053041_R04C02', 'GSM1052194_5730053041_R05C02' , 'GSM1052195_5730053041_R06C02', 'GSM1052206_5730053048_R05C01' , 'GSM1052207_5730053048_R06C01', 'GSM1052208_5730053048_R01C02' , 'GSM1052209_5730053048_R02C02', 'GSM1052211_5730053048_R04C02' , 'GSM1052212_5730053048_R05C02', 'GSM1052213_5730053048_R06C02')

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
            gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE, GSElimits = c(1, 25))
            beta_matrix <- as.data.frame(exprs(gset[[1]]))
            sample_names <- colnames(beta_matrix)
            print(sample_names)
            cpg_names <- rownames(beta_matrix)
            list_df <- list(beta_matrix, sample_names, cpg_names)
            }
            get_betamatrix()
            ''')

#robjects.r('''
#            chooseCRANmirror(ind = 1)
#            if (!require("BiocManager", quietly = TRUE))
#            install.packages("BiocManager", repos='https://cloud.r-project.org/')
#            install.packages("GEOquery", repos='https://cloud.r-project.org/')
#            install.packages("R.utils")
#            library(BiocManager)
#            library(GEOquery)
#            library(Biobase)
#            library(RnBeads)
#            library(here)
#            library(grid)
#            options(timeout = max(300, getOption("timeout")))
#            options(download.file.method.GEOquery = "wget")
#            options("fftempdir" = here::here("Temp"))
#            data.dir = here::here("Data_DNAm")
#            idat.dir = file.path(data.dir, "idats")
#            setwd("/Users/malwash/")
#            pheno.file = file.path(data.dir, "sampleSheet_withBarcode.csv")
#            analysis.dir = here::here("Output/RnBeads")
#            report.dir = file.path(analysis.dir, "diffMeth_report")
#            rnb.initialize.reports(report.dir)
#            get_betamatrix <- function(r, verbose=FALSE) {
#            rnb.options(identifiers.column = "barcode",
#            import.idat.platform = "probes450",
#            filtering.greedycut.pvalue.threshold = 0.05,
#            filtering.sex.chromosomes.removal = TRUE,
#            filtering.snp = "3",
#            filtering.cross.reactive = TRUE,
#            normalization.method = "none",
#            normalization.background.method = "none",
#            exploratory = FALSE,
#            differential = FALSE
#            )
#            data.source = c(idat.dir, pheno.file)
#            result = rnb.run.import(data.source = data.source, data.type = "infinium.idat.dir", dir.reports = report.dir)
#            rnb.set = result$rnb.set
#            unfiltered_rnb.set = rnb.set
#            filtered_results = rnb.run.preprocessing(unfiltered_rnb.set, dir.reports = report.dir)
#            filtered_rnb.set = filtered_results$rnb.set
#            betas_RnBeads = as.data.frame(meth(filtered_rnb.set, row.names = TRUE))
#            print(colnames(betas_RnBeads))
#            pheno_RnBeads = pheno(filtered_rnb.set)
#            control_samples <- c('GSM1641098_9283265129_R01C01', 'GSM1641099_9283265129_R01C02', 'GSM1641100_9283265129_R02C01','GSM1641101_9283265129_R02C02', 'GSM1641102_9283265129_R03C01', 'GSM1641103_9283265129_R03C02','GSM1641104_9283265129_R04C01', 'GSM1641105_9283265129_R04C02', 'GSM1641106_9283265129_R05C01','GSM1641107_9283265129_R06C01', 'GSM1641108_9283265129_R05C02', 'GSM1641109_9283265129_R06C02','GSM1641110_9283265149_R01C01', 'GSM1641111_9283265149_R01C02', 'GSM1641112_9283265149_R02C01','GSM1641113_9283265149_R02C02', 'GSM1641114_9283265149_R03C01', 'GSM1641115_9283265149_R04C01','GSM1641116_9283265149_R05C01', 'GSM1641117_9283265149_R06C01')
#            betas_RnBeads <- betas_RnBeads[,control_samples]
#            list_df <- list(betas_RnBeads)
#            }
#            ''')

beta_matrix_pull = robjects.r['get_betamatrix']
#list_of_samples = beta_matrix_pull()[1]
#list_of_controls = ['GSM1051533','GSM1051534','GSM1051535','GSM1051536','GSM1051537','GSM1051538','GSM1051539','GSM1051540','GSM1051541','GSM1051542','GSM1051543','GSM1051544','GSM1051545','GSM1051549','GSM1051550','GSM1051551','GSM1051552','GSM1051553','GSM1051555','GSM1051556','GSM1051557','GSM1051558','GSM1051559','GSM1051560','GSM1051561','GSM1051562','GSM1051563','GSM1051564','GSM1051565','GSM1051566','GSM1051567','GSM1051568','GSM1051569','GSM1051570','GSM1051571','GSM1051572','GSM1051573','GSM1051574','GSM1051575','GSM1051576','GSM1051577','GSM1051578','GSM1051579','GSM1051580','GSM1051581','GSM1051582','GSM1051583','GSM1051584','GSM1051585','GSM1051586','GSM1051587','GSM1051588','GSM1051589','GSM1051590','GSM1051591','GSM1051592','GSM1051593','GSM1051594','GSM1051595','GSM1051596','GSM1051597','GSM1051599','GSM1051601','GSM1051602','GSM1051603','GSM1051604','GSM1051605','GSM1051606','GSM1051607','GSM1051608','GSM1051609','GSM1051610','GSM1051618','GSM1051619','GSM1051620','GSM1051621','GSM1051622','GSM1051623','GSM1051624','GSM1051625','GSM1051626','GSM1051627','GSM1051628','GSM1051629','GSM1051630','GSM1051631','GSM1051632','GSM1051633','GSM1051634','GSM1051635','GSM1051638','GSM1051639','GSM1051641','GSM1051642','GSM1051643','GSM1051644','GSM1051645','GSM1051646','GSM1051647','GSM1051648','GSM1051649','GSM1051650','GSM1051651','GSM1051652','GSM1051653','GSM1051654','GSM1051655','GSM1051656','GSM1051657','GSM1051658','GSM1051659','GSM1051660','GSM1051661','GSM1051662','GSM1051663','GSM1051664','GSM1051665','GSM1051666','GSM1051667','GSM1051668','GSM1051669','GSM1051671','GSM1051672','GSM1051673','GSM1051674','GSM1051675','GSM1051687','GSM1051688','GSM1051689','GSM1051690','GSM1051691','GSM1051692','GSM1051693','GSM1051710','GSM1051712','GSM1051713','GSM1051714','GSM1051715','GSM1051716','GSM1051717','GSM1051726','GSM1051727','GSM1051730','GSM1051731','GSM1051732','GSM1051733','GSM1051742','GSM1051744','GSM1051745','GSM1051746','GSM1051749','GSM1051751','GSM1051752','GSM1051754','GSM1051755','GSM1051756','GSM1051757','GSM1051758','GSM1051759','GSM1051760','GSM1051761','GSM1051762','GSM1051763','GSM1051764','GSM1051765','GSM1051766','GSM1051767','GSM1051768','GSM1051769','GSM1051770','GSM1051771','GSM1051772','GSM1051773','GSM1051774','GSM1051775','GSM1051777','GSM1051781','GSM1051782','GSM1051784','GSM1051786','GSM1051787','GSM1051788','GSM1051789','GSM1051790','GSM1051791','GSM1051792','GSM1051793','GSM1051802','GSM1051803','GSM1051804','GSM1051805','GSM1051806','GSM1051807','GSM1051808','GSM1051809','GSM1051810','GSM1051811','GSM1051812','GSM1051813','GSM1051814','GSM1051815','GSM1051816','GSM1051817','GSM1051820','GSM1051821','GSM1051822','GSM1051823','GSM1051824','GSM1051825','GSM1051826','GSM1051827','GSM1051828','GSM1051831','GSM1051832','GSM1051833','GSM1051834','GSM1051835','GSM1051836','GSM1051837','GSM1051838','GSM1051839','GSM1051840','GSM1051841','GSM1051842','GSM1051843','GSM1051844','GSM1051845','GSM1051846','GSM1051847','GSM1051848','GSM1051850','GSM1051851','GSM1051852','GSM1051853','GSM1051854','GSM1051855','GSM1051857','GSM1051859','GSM1051860','GSM1051862','GSM1051863','GSM1051864','GSM1051865','GSM1051873','GSM1051874','GSM1051876','GSM1051883','GSM1051884','GSM1051901','GSM1051902','GSM1051903','GSM1051904','GSM1051905','GSM1051906','GSM1051907','GSM1051923','GSM1051924','GSM1051925','GSM1051926','GSM1051927','GSM1051928','GSM1051929','GSM1051949','GSM1051950','GSM1051951','GSM1051952','GSM1051953','GSM1051954','GSM1051955','GSM1051962','GSM1051963','GSM1051964','GSM1051965','GSM1051966','GSM1051967','GSM1051968','GSM1052015','GSM1052016','GSM1052017','GSM1052018','GSM1052019','GSM1052020','GSM1052021','GSM1052022']

#print(beta_matrix_pull()[0])
#print(pd.DataFrame(beta_matrix_pull()[0]))
beta_matrix = pd.DataFrame(beta_matrix_pull()[0]).transpose()
#beta_matrix.columns = list_of_samples

#beta_matrix = beta_matrix[:,list_of_controls]

print("The reference data - Beta matrix:")
print(beta_matrix.shape)
print(beta_matrix)

def calculate_mean(beta_matrix):
    """
    Returns an array of the length of CpGs containing mean measures (non-NaN) across samples.
    :param beta_matrix: an n (sample) by n (CpG) matrix of 0-1 methylation intensities
    :rtype list_of_means: shape mean parameters based on the beta_matrix
    """
    list_of_means = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_means:
        list_of_means[i] = np.nanmean(beta_matrix.iloc[i])
    return list_of_means

def calculate_variance(beta_matrix):
    """
    Returns an array of the length of CpGs containing variance measures (non-NaN) across samples.
    :param beta_matrix: an n (sample) by n (CpG) matrix of 0-1 methylation intensitiesh
    :rtype list_of_means: shape variance parameters based on the beta_matrix
    """
    list_of_vars = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_vars:
        list_of_vars[i] = np.nanvar(beta_matrix.iloc[i])
    return list_of_vars

def calculate_stds(beta_matrix):
    """
    Returns an array of the length of CpGs containing std measures (non-NaN) across samples.
    :param beta_matrix: an n (sample) by n (CpG) matrix of 0-1 methylation intensities
    :rtype list_of_stds: shape std parameters based on the beta_matrix
    """
    list_of_stds = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_stds:
        list_of_stds[i] = np.nanstd(beta_matrix.iloc[i])
    return list_of_stds

# Returns a matrix of ID's for CpGs and their associated mean values to generate desirable number of CpGs in the simulated data.
def user_specified_num_elements(means, user_specified_n_CpGs):
    """
    Returns an array of indices to induce differences chosen at random with replacement.
    :param means: an n (CpG) length array of shape means on which to select locations from
    :param user_specified_n_CpGs: number of CpG sites to induce differences
    :rtype sampling_indices: array of randomly selected CpG sites to induce differences
    """
    indicies_of_means = [num for num in range(0, len(means))]
    sampling_indices = random.choices(indicies_of_means,k=user_specified_n_CpGs)
    return sampling_indices

def induce_group_differnces(num_true_modified, vector_of_ref_means, effect_size):
    """
    Induce differences (i.e., motivate hypermethylation) by changing the mean shape parameters of the affected group at the locations of selected for modification.
    :param num_true_modified: user-specified input of number of CpG sites to change
    :param vector_of_ref_means: array of reference shape means to base simulation of methylation values
    :param effect_size: standardised factor to introduce effect signal within affected group
    :rtype vector_of_affected_means, truly_different_indices: Returns an array of affected means with introduced differences for simulation and the ground truth of CpG indices that were truly modified
    """
    means_of_control_population = [num for num in range(0, len(vector_of_ref_means))]
    truly_different_indices = random.sample(means_of_control_population, k=num_true_modified) #indicies to change
    vector_of_affected_means = np.array(vector_of_ref_means, dtype='f') #Copy and replicate control and seperate for case
    vector_of_affected_means[truly_different_indices] = vector_of_ref_means[truly_different_indices] + float(effect_size) #Induce the difference
    vector_of_affected_means = np.minimum(1, vector_of_affected_means)
    return vector_of_affected_means, truly_different_indices

#get number of samples per group based on the control-case balance (healthy proportion)
def get_group_number(healthy_proportion, total_num_samples):
    """
    Retrieve the number of samples per group based on the specified healthy proportion factor.
    :param healthy_proportion: proportion factor between case-control groups
    :param total_num_samples: sum of samples across groups
    :rtype group_nums: Returns a dataframe containing the number of samples in reference and disease-state (i.e., case-control) groups
    """
    group1_number_samples = healthy_proportion*total_num_samples
    group2_number_samples = total_num_samples-group1_number_samples
    groups = {
        "Reference": [group1_number_samples],
        "Diseased": [group2_number_samples]
    }
    group_nums = pd.DataFrame(groups)
    return group_nums

def sample_distribution(means_vector, stds_vector, number_of_samples_in_group):
    """
    Simulate beta-values for every CpG across all samples and truncate any value above 1 and below 0.
    :param means_vector: shape parameter array of beta-value means
    :param stds_vector: shape parameter array of beta-value stds
    :param stds_vector: user-specified number of samples to generate for a distribution
    :rtype all_cpgs: beta-matrix of CpG by samples for a group
    """
    all_cpgs = pd.DataFrame()
    #cpg_name = pd.DataFrame()
    means_vector_df = pd.DataFrame(means_vector)
    stds_vector_df = pd.DataFrame(stds_vector)
    for index in range(0, len(means_vector_df.index)):
        cpg_i = np.random.normal(loc=means_vector_df.iloc[index], scale=stds_vector_df.iloc[index], size=int(number_of_samples_in_group))
        #cpg_i = cpg_i + 0.2
        cpg_i_corrected = np.maximum(0, np.minimum(1, cpg_i)) # Remove < 0 and > 1
        all_cpgs = all_cpgs.append(pd.DataFrame(cpg_i_corrected.tolist()).T)
        #cpg_name = cpg_name.append(pd.Series(str("cpg") + str(index)), ignore_index=True)
    return all_cpgs

def generate_cpgs_for_groups(g1_means_vector, g2_means_vector, vector_std, g1_number_of_samples, g2_number_of_samples):
    """
    Simulate methylation data for every group based on reference and affected shape parameters and combine case-control group data
    :param g1_means_vector: shape parameter array of reference beta-value means
    :param g2_means_vector: shape parameter array of affected beta-value means
    :param vector_std: shape parameter array of beta-value stds
    :param g1_number_of_samples: number of samples in control
    :param g2_number_of_samples: number of samples in case
    :rtype combined_groups_cpgs: beta-matrix of CpG by samples for all groups
    """
    g1_cpgs = sample_distribution(g1_means_vector, vector_std, g1_number_of_samples)
    g2_cpgs = sample_distribution(g2_means_vector, vector_std, g2_number_of_samples)
    combined_groups_cpgs = pd.concat([g1_cpgs, g2_cpgs], axis=1)
    return(combined_groups_cpgs)

#Generate group-specific sample names
def generate_col_names(g1_number_of_samples, g2_number_of_samples):
    """
    Generate simulated group-specfic sample names
    :param g1_number_of_samples: number of samples in control
    :param g2_number_of_samples: number of samples in case
    :rtype all_sample_names: array of group-specific sample names
    """
    g1_sample_names = ["Gr1_Samp" + str(num) for num in range(0, int(g1_number_of_samples))]
    control_sample_names = pd.DataFrame(g1_sample_names).T
    g2_sample_names = ["Gr2_Samp" + str(num) for num in range(0, int(g2_number_of_samples))]
    case_sample_names = pd.DataFrame(g2_sample_names).T
    all_sample_names = pd.concat([control_sample_names, case_sample_names], axis=1)
    return all_sample_names

#Create all combinations (cartisian product) of the input parameters (5 user-defined parameters)
def get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs):
    """
    Generate all combinations (i.e., cartisian product) based on user-specified input parameters
    :param total_num_samples_vector: user-specified list of total number of samples
    :param effect_size_vector: user-specified list of experiment effect sizes
    :param healthy_proportion: user-specified list of healthy_proportions
    :param num_true_modified: user-specified list of CpG sites to truly modify in simulations
    :param user_specified_n_CpGs: user-specified list of number of CpGs to generate in simulations
    :rtype all_combinations_df: dataframe containing the product of all input parameters
    """
    all_combinations = list(itertools.product(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified,user_specified_n_CpGs))
    param_columns = ['n_samples', 'effect_size', 'healthy_proportion', 'n_true_modified_CpGs','n_CpGs']
    all_combinations_df = pd.DataFrame(all_combinations, columns=param_columns)
    workflows = ["workflow_" + str(workflow_num) for workflow_num in range(0, len(all_combinations_df))]
    all_combinations_df["workflow"] = workflows
    all_combinations_df.to_csv(os.getcwd()+summarycalcdirname+'all_combinations.csv', header=True)
    return all_combinations_df

def simulator(total_num_samples, effect_size, healthy_proportion, num_true_modified, user_specified_n_CpGs, workflow_num, num_simulations, age_group):
    """
    Run the procedural steps for drawing shape parameters from a reference dataset, perform for each specified number of simulations, the inducing of differences and the generation of beta-matrices
    :param total_num_samples: user-specified total number of samples for single experiment
    :param effect_size: user-specified experimental effect size
    :param healthy_proportion: user-specified healthy_proportion factor
    :param num_true_modified: user-specified CpG sites to truly modify in simulations
    :param user_specified_n_CpGs: user-specified number of CpGs to generate in simulations
    :param workflow_num: auto-generated number corresponding to the experimental setup within the product of all experimental setups
    :param num_simulations: user-specified number of simulations to run each experimental setup
    """
    print("Calculating shape parameters for simulation...")
    list_of_simulated_data = []
    list_of_truly_different_indices = []

    means_real_world = pd.DataFrame(calculate_mean(beta_matrix))
    stds_real_world = pd.DataFrame(calculate_stds(beta_matrix))
    shape_parameter_real_world = pd.concat([means_real_world, stds_real_world], axis=1)

    for sim_iter in range(0, num_simulations):
        #list_of_true_sample_ages = []
        indices = user_specified_num_elements(shape_parameter_real_world.iloc[:,0], user_specified_n_CpGs) #indices to sample mean/stds using sampling with replacement
        means_stds_by_indices_sample = shape_parameter_real_world.iloc[indices,:]
        vector_of_ref_means = means_stds_by_indices_sample.iloc[:, 0]
        vector_of_ref_stds = means_stds_by_indices_sample.iloc[:, 1]

        #Inducing the difference between groups
        g1_number_of_samples = get_group_number(healthy_proportion, total_num_samples).iloc[0, 0]
        g2_number_of_samples = get_group_number(healthy_proportion, total_num_samples).iloc[0, 1]
        #if age_group != 0:
        #    for _ in range(0, total_num_samples):
        #        if age_group == 1:
        #            list_of_true_sample_ages.append(random.randint(15, 25))
        #        elif age_group == 2:
        #            list_of_true_sample_ages.append(random.randint(25, 35))
        #        elif age_group == 3:
        #            list_of_true_sample_ages.append(random.randint(35, 45))
        #        elif age_group == 4:
        #            list_of_true_sample_ages.append(random.randint(45, 55))

        vector_of_affected_means, truly_different_indices = induce_group_differnces(num_true_modified,np.array(vector_of_ref_means, dtype='f'), effect_size)

        print("Generating cpgs for groups (", workflow_num, ", num_sim ",sim_iter,")")
        simulated_data = generate_cpgs_for_groups(vector_of_ref_means,vector_of_affected_means, vector_of_ref_stds, g1_number_of_samples, g2_number_of_samples)

        simulated_data_columns = generate_col_names(g1_number_of_samples, g2_number_of_samples).values.tolist()
        simulated_data.columns = simulated_data_columns
        list_of_truly_different_indices.append(truly_different_indices)
        list_of_simulated_data.append(simulated_data.to_numpy())
    file_name_user_parameters = "User_Parameters"+"_"+workflow_num+".csv"
    file_name_truly_modified = "truly_different_sites_indices"+"_"+workflow_num+".npz"
    #file_name_true_sample_ages = "true_sample_ages" + "_" + workflow_num + ".npz"
    file_name_simulated_data = "Simulated_data"+"_"+workflow_num+".npz"
    params = [['Total number of samples', total_num_samples], ['User-specified number of CpGs', user_specified_n_CpGs],
              ['Healthy proportion', healthy_proportion], ['Effect size', effect_size],
              ['Number of true modified CpG sites', num_true_modified]]#, ['Age group', age_group]]
    params_summary = pd.DataFrame(params, columns=['Parameter', 'Value'])

    if os.path.isdir(os.getcwd()+dirname+workflow_num+os.sep) == False:
        os.mkdir(os.getcwd()+dirname+workflow_num+os.sep)
    params_summary.to_csv(os.getcwd()+dirname+workflow_num+os.sep+file_name_user_parameters, header=True)
    np.savez_compressed(os.getcwd()+dirname+workflow_num+os.sep+file_name_truly_modified, list_of_truly_different_indices)
    #np.savez_compressed(os.getcwd() + dirname + workflow_num + os.sep + file_name_true_sample_ages,list_of_true_sample_ages)
    np.savez_compressed(os.getcwd()+dirname+workflow_num+os.sep+file_name_simulated_data, list_of_simulated_data)

def experimental_setups(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs, num_simulations, num_processes, age_group, variable_varied):
    if os.path.isdir(os.getcwd() + dirname) == False:
        os.mkdir(os.getcwd() + dirname)
    if os.path.isdir(os.getcwd() + summarycalcdirname) == False:
        os.mkdir(os.getcwd() + summarycalcdirname)
    if os.path.isdir(os.getcwd() + figuredirname) == False:
        os.mkdir(os.getcwd() + figuredirname)
    combination_df = get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs)
    list_of_workflows = [num for num in range(0, len(combination_df.index))]
    num_workflows = len(list_of_workflows)
    config_params = [['num_simulations', num_simulations], ['num_workflows', num_workflows], ['num_processes', num_processes], ['variable_varied', variable_varied], ['age_group', age_group]]
    config_params_summary = pd.DataFrame(config_params, columns=['Parameter', 'Value'])
    config_params_summary.to_csv(os.getcwd()+summarycalcdirname+'env_inputparams.csv', header=True)
    pool = Pool(processes=num_processes)
    result = pool.map(simulator_worker,list_of_workflows)
    pool.close()

def simulator_worker(workflow):
    combination_df = pd.read_csv(os.getcwd()+summarycalcdirname+"all_combinations.csv", header=0, index_col=0)
    env_params_df = pd.read_csv(os.getcwd() + summarycalcdirname + "env_inputparams.csv", header=0, index_col=0)
    num_simulations = int(env_params_df.iloc[0, 1])
    age_group = int(env_params_df.iloc[4, 1])
    print("Running environmental setup for workflow: ", workflow)
    simulator(combination_df.loc[workflow, 'n_samples'], combination_df.loc[workflow, 'effect_size'],
              combination_df.loc[workflow, 'healthy_proportion'], combination_df.loc[workflow, 'n_true_modified_CpGs'],
              combination_df.loc[workflow, 'n_CpGs'], combination_df.loc[workflow, 'workflow'], num_simulations, age_group)

if __name__ == '__main__':
    total_num_samples_vector = [50, 100, 200, 350, 500, 650, 800, 950]
    effect_size_vector = [0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1]  # [0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032]
    healthy_proportion = [0.5]
    num_true_modified = [15]  # [5,10,15,20,35,50,95,125]#[4,37,371,3713,18568,37137,148550,315670]#,148550,315670]#[5,10,15,20,35,50,95,125]#[4,37,371,3713,18568,37137,148550,315670]
    user_specified_n_CpGs = [1000]
    num_simulations = 50
    num_processes = 8
    age_group = 0  # 15-25, 2#25-30, 3#35-45, 4#45-60
    variable_varied = "effect_size"  # Important to change this for select x variable which is varied # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
    experimental_setups(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs, num_simulations, num_processes, age_group, variable_varied)
    print("Time taken: ",time.time() - starttime)