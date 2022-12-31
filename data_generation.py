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
            GEO_real_world_data <- "GSE77718"
            control_samples <- c('GSM2057626','GSM2057627','GSM2057628','GSM2057629','GSM2057630','GSM2057631','GSM2057632','GSM2057633','GSM2057634','GSM2057635','GSM2057636','GSM2057637','GSM2057638','GSM2057639','GSM2057640','GSM2057641','GSM2057642','GSM2057643','GSM2057662','GSM2057663','GSM2057664','GSM2057665','GSM2057666','GSM2057667','GSM2057668','GSM2057669','GSM2057670','GSM2057671','GSM2057672','GSM2057673','GSM2057674','GSM2057675','GSM2057676','GSM2057677','GSM2057678','GSM2057679','GSM2057680','GSM2057681','GSM2057682','GSM2057683','GSM2057684','GSM2057685','GSM2057686','GSM2057687','GSM2057688','GSM2057689','GSM2057690','GSM2057691','GSM2057692','GSM2057693','GSM2057694','GSM2057695','GSM2057696','GSM2057697','GSM2057698','GSM2057699','GSM2057700','GSM2057701','GSM2057702','GSM2057703','GSM2057704','GSM2057705','GSM2057706','GSM2057707','GSM2057708','GSM2057709','GSM2057758','GSM2057760','GSM2057762','GSM2057764','GSM2057766','GSM2057768','GSM2057770','GSM2057772','GSM2057774','GSM2057776','GSM2057778','GSM2057780','GSM2057782','GSM2057784','GSM2057786','GSM2057788','GSM2057790','GSM2057792','GSM2057794','GSM2057796','GSM2057798','GSM2057800','GSM2057802','GSM2057804','GSM2057806','GSM2057808','GSM2057810','GSM2057812','GSM2057814','GSM2057816')
            gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE)
            beta_matrix <- as.data.frame(exprs(gset[[1]]))
            beta_matrix <- beta_matrix[,control_samples]
            sample_names <- colnames(beta_matrix)
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
#            parallel.setup(8)
#            options(timeout = max(300, getOption("timeout")))
#            options(download.file.method.GEOquery = "wget")
#            options("fftempdir" = here::here("Temp"))
#            data.dir = here::here("Data_DNAm")
#            idat.dir = file.path(data.dir, "idats")
#            print(idat.dir)
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
#            normalization.method = "bmiq",
#            normalization.background.method = "enmix.oob",
#            exploratory = FALSE,
#            differential = FALSE
#            )
#            data.source = c(idat.dir, pheno.file)
#            result = rnb.run.import(data.source = data.source, data.type = "infinium.idat.dir", dir.reports = report.dir)
#            rnb.set = result$rnb.set
#            unfiltered_rnb.set = rnb.set
#            filtered_results = rnb.run.preprocessing(unfiltered_rnb.set, dir.reports = report.dir)
#            filtered_rnb.set = filtered_results$rnb.set
#            betas_RnBeads = as.matrix(as.data.frame(meth(filtered_rnb.set, row.names = TRUE)))
#            print(betas_RnBeads)
#            pheno_RnBeads = pheno(filtered_rnb.set)
#            print(pheno_RnBeads)
#            list_df <- list(betas_RnBeads)
#            }
#            get_betamatrix()
#            ''')
#GEO_real_world_data <- "GSE68777"
#gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE)
#idat_paths <- getGEOSuppFiles("GSE68777", filter_regex = "GSE68777_RAW.tar")
#variables <- names(idat_paths)
#

beta_matrix_pull = robjects.r['get_betamatrix']
list_of_samples = beta_matrix_pull()[1]
beta_matrix = pd.DataFrame(beta_matrix_pull()[0]).transpose()
beta_matrix.columns = list_of_samples

print("The reference data - Beta matrix:")
print(beta_matrix.shape)
print(beta_matrix)

def calculate_mean(beta_matrix):
    """
    Returns an array of the length of CpGs containing mean estimates (non-NaN) for CpGs across samples.
    :param beta_matrix: an n (sample) by n (CpG) matrix of 0-1 methylation intensities
    :rtype list_of_means: shape mean parameters based on the beta_matrix
    """
    list_of_means = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_means:
        list_of_means[i] = np.nanmean(beta_matrix.iloc[i])
    return list_of_means

def calculate_variance(beta_matrix):
    list_of_vars = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_vars:
        list_of_vars[i] = np.nanvar(beta_matrix.iloc[i])
    return list_of_vars

def calculate_stds(beta_matrix):
    list_of_stds = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_stds:
        list_of_stds[i] = np.nanstd(beta_matrix.iloc[i])
    return list_of_stds

# Returns a matrix of ID's for CpGs and their associated mean values to generate desirable number of CpGs in the simulated data.
def user_specified_num_elements(means, user_specified_n_CpGs):
    indicies_of_means = [num for num in range(0, len(means))]
    sampling_indices = random.choices(indicies_of_means,k=user_specified_n_CpGs)
    return(sampling_indices)

# Induce differences (i.e., hypermetalytation) between groups by changing mean beta-values of selected CpGs. The modified CpGs are sampled from mean values of the reference group (control group). The output then recordes indices of the truly modified CpGs and vector of means with introduced corresponding differences.
# This function operates based on indices between a control and case
def induce_group_differnces(num_true_modified, vector_of_ref_means, effect_size, sim_iter):
    means_of_control_population = [num for num in range(0, len(vector_of_ref_means))]
    truly_different_indices = random.sample(means_of_control_population, k=num_true_modified) #indicies to change
    vector_of_affected_means = np.array(vector_of_ref_means, dtype='f') #Copy and replicate control and seperate for case
    vector_of_affected_means[truly_different_indices] = vector_of_ref_means[truly_different_indices] + float(effect_size) #Induce the difference
    vector_of_affected_means = np.minimum(1, vector_of_affected_means)
    return vector_of_affected_means, truly_different_indices

#get number of samples per group based on the control-case balance (healthy proportion)
def get_group_number(healthy_proportion, total_num_samples):
    group1_number_samples = healthy_proportion*total_num_samples
    group2_number_samples = total_num_samples-group1_number_samples
    groups = {
        "Reference": [group1_number_samples],
        "Diseased": [group2_number_samples]
    }
    group_nums = pd.DataFrame(groups)
    return group_nums

#Simulate beta values for every CpG across all samples. Due to the range of the real beta values (0 and 1), truncate/adjust the values within the range
def sample_distribution(means_vector, stds_vector, number_of_samples_in_group):
    all_cpgs = pd.DataFrame()
    #cpg_name = pd.DataFrame()
    means_vector_df = pd.DataFrame(means_vector)
    stds_vector_df = pd.DataFrame(stds_vector)
    for index in range(0, len(means_vector_df.index)):
        cpg_i = np.random.normal(loc=means_vector_df.iloc[index], scale=stds_vector_df.iloc[index], size=int(number_of_samples_in_group))
        cpg_i_corrected = np.maximum(0, np.minimum(1, cpg_i)) # Remove < 0 and > 1
        all_cpgs = all_cpgs.append(pd.DataFrame(cpg_i_corrected.tolist()).T)
        #cpg_name = cpg_name.append(pd.Series(str("cpg") + str(index)), ignore_index=True)
    return all_cpgs

#Simulate methylation data for every group based on the shape parameters specific to the group (different means and common standard deviation). The output is data frame with all samples and all CpGs.
def generate_cpgs_for_groups(g1_means_vector, g2_means_vector, vector_std, g1_number_of_samples, g2_number_of_samples):
    g1_cpgs = sample_distribution(g1_means_vector, vector_std, g1_number_of_samples)
    g2_cpgs = sample_distribution(g2_means_vector, vector_std, g2_number_of_samples)
    combined_groups_cpgs = pd.concat([g1_cpgs, g2_cpgs], axis=1)
    return(combined_groups_cpgs)

#Generate group-specific sample names
def generate_col_names(g1_number_of_samples, g2_number_of_samples):
    g1_sample_names = ["Gr1_Samp" + str(num) for num in range(0, int(g1_number_of_samples))]
    control_sample_names = pd.DataFrame(g1_sample_names).T
    g2_sample_names = ["Gr2_Samp" + str(num) for num in range(0, int(g2_number_of_samples))]
    case_sample_names = pd.DataFrame(g2_sample_names).T
    all_sample_names = pd.concat([control_sample_names, case_sample_names], axis=1)
    return all_sample_names

#Create all combinations (cartisian product) of the input parameters (5 user-defined parameters)
def get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs):
    all_combinations = list(
        itertools.product(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified,
                          user_specified_n_CpGs))

    param_columns = ['n_samples', 'effect_size', 'healthy_proportion', 'n_true_modified_CpGs','n_CpGs']
    all_combinations_df = pd.DataFrame(all_combinations, columns=param_columns)
    workflows = ["workflow_" + str(workflow_num) for workflow_num in range(0, len(all_combinations_df))]
    all_combinations_df["workflow"] = workflows
    all_combinations_df.to_csv(os.getcwd()+summarycalcdirname+'all_combinations.csv', header=True)
    return all_combinations_df

def simulator(total_num_samples, effect_size, healthy_proportion, num_true_modified, user_specified_n_CpGs, workflow_num, num_simulations):
    print("Calculating shape parameters for simulation...")
    list_of_simulated_data = []
    list_of_truly_different_indices = []

    means_real_world = pd.DataFrame(calculate_mean(beta_matrix))
    stds_real_world = pd.DataFrame(calculate_stds(beta_matrix))
    shape_parameter_real_world = pd.concat([means_real_world, stds_real_world], axis=1)

    for sim_iter in range(0, num_simulations):
        indices = user_specified_num_elements(shape_parameter_real_world.iloc[:,0], user_specified_n_CpGs) #indices to sample mean/stds using sampling with replacement
        means_stds_by_indicies_sample = shape_parameter_real_world.iloc[indices,:]
        vector_of_ref_means = means_stds_by_indicies_sample.iloc[:, 0]
        vector_of_ref_stds = means_stds_by_indicies_sample.iloc[:, 1]

        print("Inducing differences between groups...")
        vector_of_affected_means, truly_different_indices = induce_group_differnces(num_true_modified,np.array(vector_of_ref_means, dtype='f'), effect_size, sim_iter)
        g1_number_of_samples = get_group_number(healthy_proportion, total_num_samples).iloc[0, 0]
        g2_number_of_samples = get_group_number(healthy_proportion, total_num_samples).iloc[0, 1]

        print("Generating cpgs for groups (", workflow_num, ", num_sim ",sim_iter,")")
        simulated_data = generate_cpgs_for_groups(vector_of_ref_means,vector_of_affected_means, vector_of_ref_stds, g1_number_of_samples, g2_number_of_samples)

        simulated_data_columns = generate_col_names(g1_number_of_samples, g2_number_of_samples).values.tolist()
        simulated_data.columns = simulated_data_columns
        list_of_truly_different_indices.append(truly_different_indices)
        list_of_simulated_data.append(simulated_data.to_numpy())
    file_name_user_parameters = "User_Parameters"+"_"+workflow_num+".csv"
    file_name_truly_modified = "truly_different_sites_indices"+"_"+workflow_num+".npz"
    file_name_simulated_data = "Simulated_data"+"_"+workflow_num+".npz"
    params = [['Total number of samples', total_num_samples], ['User-specified number of CpGs', user_specified_n_CpGs],
              ['Healthy proportion', healthy_proportion], ['Effect size', effect_size],
              ['Number of true modified CpG sites', num_true_modified]]
    params_summary = pd.DataFrame(params, columns=['Parameter', 'Value'])

    if os.path.isdir(os.getcwd()+dirname+workflow_num+os.sep) == False:
        os.mkdir(os.getcwd()+dirname+workflow_num+os.sep)
    params_summary.to_csv(os.getcwd()+dirname+workflow_num+os.sep+file_name_user_parameters, header=True)
    np.savez_compressed(os.getcwd()+dirname+workflow_num+os.sep+file_name_truly_modified, list_of_truly_different_indices)
    np.savez_compressed(os.getcwd()+dirname+workflow_num+os.sep+file_name_simulated_data, list_of_simulated_data)

def multi_simMethyl():
    if os.path.isdir(os.getcwd() + dirname) == False:
        os.mkdir(os.getcwd() + dirname)
    if os.path.isdir(os.getcwd() + summarycalcdirname) == False:
        os.mkdir(os.getcwd() + summarycalcdirname)
    if os.path.isdir(os.getcwd() + figuredirname) == False:
        os.mkdir(os.getcwd() + figuredirname)

    total_num_samples_vector = [50,100,200,350,500,650,800,950]
    effect_size_vector = [0.01]#,0.02,0.03,0.05,0.07,0.08,0.09,0.1]  # [0.01,0.04,0.07,0.1,0.13,0.16,0.19,0.22]
    healthy_proportion = [0.5]
    num_true_modified = [4,37,371,3713,18568,37137,148550,315670]#,148550,315670]#[5,10,15,20,35,50,95,125]#[4,37,371,3713,18568,37137,148550,315670]
    user_specified_n_CpGs = [371377] # [1000]
    num_simulations = 5
    combination_df = get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs)
    list_of_workflows = [num for num in range(0, len(combination_df.index))]
    num_workflows = len(list_of_workflows)

    #Important to change this for select variable which is varied
    variable_varied = "n_modified_CpGs"#"n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
    config_params = [['num_simulations', num_simulations], ['num_workflows', num_workflows], ['variable_varied', variable_varied]]
    config_params_summary = pd.DataFrame(config_params, columns=['Parameter', 'Value'])
    config_params_summary.to_csv(os.getcwd()+summarycalcdirname+'env_inputparams.csv', header=True)
    pool = Pool(processes=64) # user specified CPUs e.g., processes=8
    result = pool.map(simulator_pool,list_of_workflows)
    pool.close()

def simulator_pool(workflow):
    total_num_samples_vector = [50,100,200,350,500,650,800,950]
    effect_size_vector = [0.01]#,0.02,0.03,0.05,0.07,0.08,0.09,0.1]#[0.01,0.04,0.07,0.1,0.13,0.16,0.19,0.22]
    healthy_proportion = [0.5]
    num_true_modified = [4,37,371,3713,18568,37137,148550,315670]  # [5,10,15,20,35,50,95,125]#[4,37,371,3713,18568,37137,148550,315670]
    user_specified_n_CpGs = [371377]  # [1000]
    num_simulations = 5

    combination_df = get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion,num_true_modified, user_specified_n_CpGs)
    print("Running environmental setup for workflow: ", workflow)
    simulator(combination_df.loc[workflow, 'n_samples'], combination_df.loc[workflow, 'effect_size'],
              combination_df.loc[workflow, 'healthy_proportion'], combination_df.loc[workflow, 'n_true_modified_CpGs'],
              combination_df.loc[workflow, 'n_CpGs'], combination_df.loc[workflow, 'workflow'], num_simulations)

if __name__ == '__main__':
    multi_simMethyl()
    print("Time taken: ",time.time() - starttime)