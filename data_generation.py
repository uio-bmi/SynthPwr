import csv
import itertools
import os
from zipfile import ZipFile
import numpy as np
import pandas as pd
import random
import time
import warnings
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

#from power_evaluation import limma_test

pd.options.display.max_rows
warnings.simplefilter(action='ignore', category=FutureWarning)

starttime = time.time()
dirname = "/experiments/"
#GeoQuery emulated retrieval
robjects.r('''
            library(GEOquery)
            library(Biobase)
            get_betamatrix <- function(r, verbose=FALSE) {
            GEO_real_world_data <- "GSE145714"
            gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE)
            beta_matrix <- as.data.frame(exprs(gset[[1]]))
            sample_names <- colnames(beta_matrix)
            cpg_names <- rownames(beta_matrix)
            list_df <- list(beta_matrix, sample_names, cpg_names)
            }
            get_betamatrix()
            ''')
beta_matrix_pull = robjects.r['get_betamatrix']
list_of_samples = beta_matrix_pull()[1]
#print("sample_names:")
#print([beta_matrix_pull()[1]])
#print("cpg_names:")
#print([beta_matrix_pull()[2]])
#list_of_cpgs = beta_matrix_pull()[2]

beta_matrix = pd.DataFrame(beta_matrix_pull()[0]).transpose()
beta_matrix.columns = list_of_samples
#beta_matrix.index = list_of_cpgs


print("The reference data - Beta matrix:")
#print(beta_matrix.shape)
#print(beta_matrix)


#GeoParse retrieval
#import GEOparse
#gse = GEOparse.get_GEO(geo="GSE145714", destdir="./", include_data=True)
#for gsm_name, gsm in gse.gsms.items():
#    print("Name: ", gsm_name)
#    print("Metadata: ",)
#    for key, value in gsm.metadata.items():
#        print(" - %s : %s" % (key, ", ".join(value)))
#    print("Table data: ",)
#    print(gsm.table)

#for gpl_name, gpl in gse.gpls.items():
#    print("Name: ", gpl_name)
#    print("Metadata:",)
#    for key, value in gpl.metadata.items():
#        print(" - %s : %s" % (key, ", ".join(value)))
#    print("Table data:",)
#    print(gpl.table.head())

# Function for calculating the mean of CpGs across all samples.
def calculate_mean(beta_matrix):
    list_of_means = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_means:
        list_of_means[i] = np.mean(beta_matrix.iloc[i])
    return list_of_means

def calculate_variance(beta_matrix):
    list_of_vars = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_vars:
        list_of_vars[i] = np.var(beta_matrix.iloc[i])
    return list_of_vars

def calculate_stds(beta_matrix):
    list_of_stds = [row for row in range(0, len(beta_matrix.index))]
    for i in list_of_stds:
        list_of_stds[i] = np.std(beta_matrix.iloc[i])
    return list_of_stds

# Returns a matrix of ID's for CpGs and their associated mean values to generate desirable number of CpGs in the simulated data.
def user_specified_num_elements(means, user_specified_n_CpGs):
    indicies_of_means = [num for num in range(0, len(means))]
    sampling_indices = random.choices(indicies_of_means,k=user_specified_n_CpGs) #sampling with replacement
    return(sampling_indices)

# Induce differences (i.e., hypermetalytation) between groups by changing mean beta-values of selected CpGs. The modified CpGs are sampled from mean values of the reference group (control group). The output then recordes indices of the truly modified CpGs and vector of means with introduced corresponding differences.
# This function operates based on indices between a control and case
def induce_group_differnces(num_true_modified, vector_of_ref_means, effect_size):
    means_of_control_population = [num for num in range(0, len(vector_of_ref_means))]
    truly_different_indices = random.sample(means_of_control_population, k=num_true_modified) #indicies to change
    with open(os.getcwd()+dirname+'truly_different_sites_indices.txt', 'w') as f:
        for r in truly_different_indices:
            f.write(str(r)+" ")

    vector_of_affected_means = np.array(vector_of_ref_means, dtype='f') #Copy and replicate control and seperate from case
    vector_of_affected_means[truly_different_indices] = vector_of_ref_means[truly_different_indices] + float(effect_size) #Induce the difference
    vector_of_affected_means = np.minimum(1, vector_of_affected_means)
    return vector_of_affected_means

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
    cpg_name = pd.DataFrame()
    means_vector_df = pd.DataFrame(means_vector)
    stds_vector_df = pd.DataFrame(stds_vector)

    for index, row in enumerate(means_vector.tolist()):
        cpg_i = np.random.normal(loc=means_vector_df.iloc[index], scale=stds_vector_df.iloc[index], size=int(number_of_samples_in_group))
        all_cpgs = all_cpgs.append(pd.DataFrame(cpg_i.tolist()).T)
        cpg_name = cpg_name.append(pd.Series(str("cpg") + str(index)), ignore_index=True)
    for index, row in all_cpgs.iterrows():
        for j, value in all_cpgs.iloc[index].items():
            if float(value) > float(1.0):
                all_cpgs.iloc[index, j] = max(beta_matrix.iloc[index])
            elif float(value) < float(0):
                all_cpgs.iloc[index, j] = min(beta_matrix.iloc[index])
    all_cpgs.index = cpg_name.iloc[:,0].tolist()
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
    #Create a cartesian product of parameter combinations
    all_combinations = list(
        itertools.product(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified,
                          user_specified_n_CpGs))

    param_columns = ['n_samples', 'effect_size', 'healthy_proportion', 'n_true_modified_CpGs',
                     'n_CpGs']
    all_combinations_df = pd.DataFrame(all_combinations, columns=param_columns)
    workflows = ["SimMethyl_run_" + str(workflow_num) for workflow_num in range(0, len(all_combinations_df))]
    all_combinations_df["workflow"] = workflows

    if os.path.isdir(os.getcwd() + dirname) == False:
        os.mkdir(os.getcwd() + dirname)

    all_combinations_df.to_csv(os.getcwd()+dirname+'all_combinations.csv', header=True)
    return all_combinations_df

def simulator(total_num_samples, effect_size, healthy_proportion, num_true_modified, user_specified_n_CpGs, workflow_num):
    print("Calculating shape parameters for simulation...")
    means_real_world = pd.DataFrame(calculate_mean(beta_matrix))
    stds_real_world = pd.DataFrame(calculate_stds(beta_matrix))
    shape_parameter_real_world = pd.concat([means_real_world, stds_real_world], axis=1)

    indices = user_specified_num_elements(shape_parameter_real_world.iloc[:,0], user_specified_n_CpGs) #indices to sample mean/stds
    means_stds_by_indicies_sample = shape_parameter_real_world.iloc[indices,:]
    vector_of_ref_means = means_stds_by_indicies_sample.iloc[:, 0]
    vector_stds = means_stds_by_indicies_sample.iloc[:, 1]
    print("Inducing differences between groups...")
    vector_of_affected_means = induce_group_differnces(num_true_modified,np.array(vector_of_ref_means, dtype='f'), effect_size)
    g1_number_of_samples = get_group_number(healthy_proportion, total_num_samples).iloc[0, 0]
    g2_number_of_samples = get_group_number(healthy_proportion, total_num_samples).iloc[0, 1]
    print("Generating cpgs for groups...")
    simulated_data = generate_cpgs_for_groups(vector_of_ref_means,vector_of_affected_means, vector_stds, g1_number_of_samples, g2_number_of_samples)

    simulated_data_columns = generate_col_names(g1_number_of_samples, g2_number_of_samples).values.tolist()
    simulated_data.columns = simulated_data_columns
    file_name_simulated_data = "Simulated_data.txt"

    print("The synthetic data - Beta matrix:")
    print(simulated_data)
    simulated_data.to_csv(os.getcwd()+dirname+file_name_simulated_data,index=True,header=True, sep=' ')

    file_name_user_parameters = "User_Parameters.csv"
    file_name_truly_modified_indices = "truly_different_sites_indices.txt"
    params = [['Total number of samples', total_num_samples], ['User-spesified number of CpGs', user_specified_n_CpGs], ['Healthy proportion', healthy_proportion], ['Effect size', effect_size], ['Number of true modified CpG sites', num_true_modified]]
    params_summary = pd.DataFrame(params, columns=['Parameter', 'Value'])
    params_summary.to_csv(os.getcwd()+dirname+file_name_user_parameters, header=True)

    with ZipFile(workflow_num+str('.zip'), 'w') as zipObj:
        zipObj.write(os.getcwd()+dirname+file_name_simulated_data)
        zipObj.write(os.getcwd()+dirname+file_name_user_parameters)
        zipObj.write(os.getcwd()+dirname+file_name_truly_modified_indices)

def multi_simMethyl(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs):
    combination_df = get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs)
    for i in range(len(combination_df)):
        print("Running environmental setup using parameters:")
        print(combination_df.iloc[i])
        simulator(combination_df.loc[i, 'n_samples'], combination_df.loc[i, 'effect_size'],
                  combination_df.loc[i, 'healthy_proportion'], combination_df.loc[i, 'n_true_modified_CpGs'],
                  combination_df.loc[i, 'n_CpGs'], combination_df.loc[i, 'workflow'])


healthy_proportion = [0.5]
num_true_modified = [1,4,7,10,13]
user_specified_n_CpGs = [10000]
total_num_samples_vector = [50, 100,200,350,500,650,800,950]
effect_size_vector = [0.05]
multi_simMethyl(total_num_samples_vector, effect_size_vector, healthy_proportion,num_true_modified,user_specified_n_CpGs)

totaltime = round((time.time() - starttime), 2)
print("Total time in seconds: ",str(totaltime))