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

pandas2ri.activate()
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

# Scenario 1 - Mixing controls GSE30870 and GSE124366
#GEO_real_world_data < - "GSE30870"
#gset < - getGEO(GEO_real_world_data, GSEMatrix=TRUE, getGPL=FALSE, GSElimits=c(1, 25))
#beta_matrix < - as.data.frame(exprs(gset[[1]]))
#control_samples < - c('GSM765860', 'GSM765862', 'GSM765863', 'GSM765864', 'GSM765865', 'GSM765866', 'GSM765867',
#                      'GSM765868', 'GSM765869', 'GSM765870', 'GSM765871', 'GSM765872', 'GSM765873', 'GSM765874',
#                      'GSM765875', 'GSM765876', 'GSM765877', 'GSM765878', 'GSM765879', 'GSM765880')
#beta_matrix < - beta_matrix[0:485000, control_samples]
#GEO_real_world_data < - "GSE124366"
#gset < - getGEO(GEO_real_world_data, GSEMatrix=TRUE, getGPL=FALSE, GSElimits=c(1, 25))
#second_beta_matrix < - as.data.frame(exprs(gset[[1]]))
#second_beta_matrix < - second_beta_matrix[0:485000, ]
#beta_matrix < - cbind(beta_matrix, second_beta_matrix)

# Scenario 2 - GSE30870 Nonagenarians and newborns
# control_samples <- c('GSM765860', 'GSM765862', 'GSM765863', 'GSM765864', 'GSM765865', 'GSM765866', 'GSM765867', 'GSM765868', 'GSM765869', 'GSM765870', 'GSM765871', 'GSM765872', 'GSM765873', 'GSM765874', 'GSM765875', 'GSM765876', 'GSM765877', 'GSM765878', 'GSM765879', 'GSM765880')
# control_samples <- c('GSM765861', 'GSM765881', 'GSM765882', 'GSM765883', 'GSM765884', 'GSM765885', 'GSM765886', 'GSM765887', 'GSM765888', 'GSM765889', 'GSM765890', 'GSM765891', 'GSM765892', 'GSM765893', 'GSM765894', 'GSM765895', 'GSM765896', 'GSM765897', 'GSM765898', 'GSM765899')
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
            GEO_real_world_data <- "GSE30870"
            gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE, GSElimits = c(1, 25))
            beta_matrix <- as.data.frame(exprs(gset[[1]]))
            control_samples <- c('GSM765861', 'GSM765881', 'GSM765882', 'GSM765883', 'GSM765884', 'GSM765885', 'GSM765886', 'GSM765887', 'GSM765888', 'GSM765889', 'GSM765890', 'GSM765891', 'GSM765892', 'GSM765893', 'GSM765894', 'GSM765895', 'GSM765896', 'GSM765897', 'GSM765898', 'GSM765899')
            beta_matrix <- beta_matrix[,control_samples]
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
#            list_df <- list(betas_RnBeads)
#            }
#            ''')

def calculate_shape_parameters(beta_matrix):
    list_of_means = [row for row in range(0, len(beta_matrix.index))]
    list_of_vars = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_means:
        list_of_means[i] = np.nanmean(beta_matrix.iloc[i])
        list_of_vars[i] = np.nanvar(beta_matrix.iloc[i])
    methPara = pd.DataFrame(list(zip(list_of_means, list_of_vars)), index=beta_matrix.index, columns=["mu","var"])
    return methPara

#list_of_controls = ['GSM1051533','GSM1051534','GSM1051535','GSM1051536','GSM1051537','GSM1051538','GSM1051539','GSM1051540','GSM1051541','GSM1051542','GSM1051543','GSM1051544','GSM1051545','GSM1051549','GSM1051550','GSM1051551','GSM1051552','GSM1051553','GSM1051555','GSM1051556','GSM1051557','GSM1051558','GSM1051559','GSM1051560','GSM1051561','GSM1051562','GSM1051563','GSM1051564','GSM1051565','GSM1051566','GSM1051567','GSM1051568','GSM1051569','GSM1051570','GSM1051571','GSM1051572','GSM1051573','GSM1051574','GSM1051575','GSM1051576','GSM1051577','GSM1051578','GSM1051579','GSM1051580','GSM1051581','GSM1051582','GSM1051583','GSM1051584','GSM1051585','GSM1051586','GSM1051587','GSM1051588','GSM1051589','GSM1051590','GSM1051591','GSM1051592','GSM1051593','GSM1051594','GSM1051595','GSM1051596','GSM1051597','GSM1051599','GSM1051601','GSM1051602','GSM1051603','GSM1051604','GSM1051605','GSM1051606','GSM1051607','GSM1051608','GSM1051609','GSM1051610','GSM1051618','GSM1051619','GSM1051620','GSM1051621','GSM1051622','GSM1051623','GSM1051624','GSM1051625','GSM1051626','GSM1051627','GSM1051628','GSM1051629','GSM1051630','GSM1051631','GSM1051632','GSM1051633','GSM1051634','GSM1051635','GSM1051638','GSM1051639','GSM1051641','GSM1051642','GSM1051643','GSM1051644','GSM1051645','GSM1051646','GSM1051647','GSM1051648','GSM1051649','GSM1051650','GSM1051651','GSM1051652','GSM1051653','GSM1051654','GSM1051655','GSM1051656','GSM1051657','GSM1051658','GSM1051659','GSM1051660','GSM1051661','GSM1051662','GSM1051663','GSM1051664','GSM1051665','GSM1051666','GSM1051667','GSM1051668','GSM1051669','GSM1051671','GSM1051672','GSM1051673','GSM1051674','GSM1051675','GSM1051687','GSM1051688','GSM1051689','GSM1051690','GSM1051691','GSM1051692','GSM1051693','GSM1051710','GSM1051712','GSM1051713','GSM1051714','GSM1051715','GSM1051716','GSM1051717','GSM1051726','GSM1051727','GSM1051730','GSM1051731','GSM1051732','GSM1051733','GSM1051742','GSM1051744','GSM1051745','GSM1051746','GSM1051749','GSM1051751','GSM1051752','GSM1051754','GSM1051755','GSM1051756','GSM1051757','GSM1051758','GSM1051759','GSM1051760','GSM1051761','GSM1051762','GSM1051763','GSM1051764','GSM1051765','GSM1051766','GSM1051767','GSM1051768','GSM1051769','GSM1051770','GSM1051771','GSM1051772','GSM1051773','GSM1051774','GSM1051775','GSM1051777','GSM1051781','GSM1051782','GSM1051784','GSM1051786','GSM1051787','GSM1051788','GSM1051789','GSM1051790','GSM1051791','GSM1051792','GSM1051793','GSM1051802','GSM1051803','GSM1051804','GSM1051805','GSM1051806','GSM1051807','GSM1051808','GSM1051809','GSM1051810','GSM1051811','GSM1051812','GSM1051813','GSM1051814','GSM1051815','GSM1051816','GSM1051817','GSM1051820','GSM1051821','GSM1051822','GSM1051823','GSM1051824','GSM1051825','GSM1051826','GSM1051827','GSM1051828','GSM1051831','GSM1051832','GSM1051833','GSM1051834','GSM1051835','GSM1051836','GSM1051837','GSM1051838','GSM1051839','GSM1051840','GSM1051841','GSM1051842','GSM1051843','GSM1051844','GSM1051845','GSM1051846','GSM1051847','GSM1051848','GSM1051850','GSM1051851','GSM1051852','GSM1051853','GSM1051854','GSM1051855','GSM1051857','GSM1051859','GSM1051860','GSM1051862','GSM1051863','GSM1051864','GSM1051865','GSM1051873','GSM1051874','GSM1051876','GSM1051883','GSM1051884','GSM1051901','GSM1051902','GSM1051903','GSM1051904','GSM1051905','GSM1051906','GSM1051907','GSM1051923','GSM1051924','GSM1051925','GSM1051926','GSM1051927','GSM1051928','GSM1051929','GSM1051949','GSM1051950','GSM1051951','GSM1051952','GSM1051953','GSM1051954','GSM1051955','GSM1051962','GSM1051963','GSM1051964','GSM1051965','GSM1051966','GSM1051967','GSM1051968','GSM1052015','GSM1052016','GSM1052017','GSM1052018','GSM1052019','GSM1052020','GSM1052021','GSM1052022']
#beta_matrix = beta_matrix[:,list_of_controls]

def get_power_calculation(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims):
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
                percentile[i] <- stats::quantile(abs(delta),0.9999)
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
                } else {
                tau <- deltaSD 
                K <- NULL
                for(d in seq_along(tau)){
                K[d] <- getK(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, tau)
                }
                }
                
                cl <- parallel::makeCluster(core)
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
                muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] <- muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] + delta
                
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
                }
            
                outSim <- list() 
                outSim[["power"]] <- marPower 
                outSim[["delta"]] <- deltaSim 
                outSim[["metric"]]$marTypeI <- marTypeI
                outSim[["metric"]]$FDR <- FDR
                outSim[["metric"]]$classicalPower <- classicalPower
                outSim[["metric"]]$FDC <- FDC
                outSim[["metric"]]$probTP <- probTP
                outSim[["metric"]]$g1Beta <- g1Beta
                outSim[["metric"]]$g2Beta <- g2Beta
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
    sample_steps = np.arange(minTotSampleSize, maxTotSampleSize + SampleSizeSteps, SampleSizeSteps)
    output_powerarray, output_metrics, output_fdr = power_result(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims)

    #fdr_df = pd.DataFrame(index=range(0,7),columns=targetDelta)
    #print(output_g1var)
    #print(pd.DataFrame(output_g1var))
    #print(output_g2var)
    #print(pd.DataFrame(output_g2var))
    fdr_effectsizes = {item: [] for item in targetDelta}
    for delta_row in range(0, len(sample_steps)):
        for sim in output_fdr:
            for idx, effectsize in enumerate(fdr_effectsizes):
                fdr_effectsizes[effectsize].append(sim[delta_row, idx])
    for effectsize in fdr_effectsizes:
        fdr_effectsizes[effectsize].sort(reverse=True)
        sns.kdeplot(fdr_effectsizes[effectsize]).set(title='Frequency Distribution of adjusted p-values')
    plt.xlabel("FDR (q-value)")

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

def synthPwr(minTotSampleSize,maxTotSampleSize,SampleSizeSteps,NcntPer,targetDelta,deltaSD=None,J=100000,targetDmCpGs=100,tissueType="GSE67170",detectionLimit=0.01,DMmethod=["limma", "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc"],FDRcritVal=0.05,core=4,sims=50):
    beta_matrix_pull = robjects.r['get_betamatrix']
    beta_matrix = pd.DataFrame(beta_matrix_pull()[0]).transpose()
    #list_of_samples = beta_matrix_pull()[1]
    #list_of_cpgs = beta_matrix_pull()[2]
    #beta_matrix.columns = list_of_samples
    #beta_matrix.index = list_of_cpgs
    print("The reference data - Beta matrix:")
    print(beta_matrix.shape)
    print(beta_matrix)

    methPara = calculate_shape_parameters(beta_matrix)
    CpGonArray = len(methPara['mu'])
    get_power_calculation(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, targetDelta, DMmethod, minTotSampleSize, maxTotSampleSize, SampleSizeSteps, FDRcritVal, NcntPer, core, sims)

if __name__ == '__main__':
    input = {}
    input['Nmin'] = 50
    input['Nmax'] = 950
    input['NCntPer'] = 0.5
    input['Nsteps'] = 150
    input['J'] = 10000
    input['targetDmCpGs'] = 10
    input['targetDeltaString'] = "0.1, 0.3, 0.5"
    #input['tauString'] = "0.01, 0.03"
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