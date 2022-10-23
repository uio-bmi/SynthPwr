import csv
import itertools
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
user_specified_n_CpGs = [371377]

pd.options.display.max_rows
warnings.simplefilter(action='ignore', category=FutureWarning)

starttime = time.time()
dirname = os.sep+"dgp_experiments"+os.sep
#GeoQuery emulated retrieval
robjects.r('''
            BiocManager::install("GEOquery", update=FALSE)
            library(GEOquery)
            library(Biobase)
            get_betamatrix <- function(r, verbose=FALSE) {
            GEO_real_world_data <- "GSE77718"
            gset <- getGEO(GEO_real_world_data, GSEMatrix =TRUE, getGPL=FALSE)
            beta_matrix <- as.data.frame(exprs(gset[[1]]))
            control_samples <- c('GSM2057626','GSM2057627','GSM2057628','GSM2057629','GSM2057630','GSM2057631','GSM2057632','GSM2057633','GSM2057634','GSM2057635','GSM2057636','GSM2057637','GSM2057638','GSM2057639','GSM2057640','GSM2057641','GSM2057642','GSM2057643','GSM2057662','GSM2057663','GSM2057664','GSM2057665','GSM2057666','GSM2057667','GSM2057668','GSM2057669','GSM2057670','GSM2057671','GSM2057672','GSM2057673','GSM2057674','GSM2057675','GSM2057676','GSM2057677','GSM2057678','GSM2057679','GSM2057680','GSM2057681','GSM2057682','GSM2057683','GSM2057684','GSM2057685','GSM2057686','GSM2057687','GSM2057688','GSM2057689','GSM2057690','GSM2057691','GSM2057692','GSM2057693','GSM2057694','GSM2057695','GSM2057696','GSM2057697','GSM2057698','GSM2057699','GSM2057700','GSM2057701','GSM2057702','GSM2057703','GSM2057704','GSM2057705','GSM2057706','GSM2057707','GSM2057708','GSM2057709','GSM2057758','GSM2057760','GSM2057762','GSM2057764','GSM2057766','GSM2057768','GSM2057770','GSM2057772','GSM2057774','GSM2057776','GSM2057778','GSM2057780','GSM2057782','GSM2057784','GSM2057786','GSM2057788','GSM2057790','GSM2057792','GSM2057794','GSM2057796','GSM2057798','GSM2057800','GSM2057802','GSM2057804','GSM2057806','GSM2057808','GSM2057810','GSM2057812','GSM2057814','GSM2057816')
            beta_matrix <- beta_matrix[,control_samples]
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
print(beta_matrix.shape)
print(beta_matrix)

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

    for index in range(0, len(means_vector_df.index)):
        cpg_i = np.random.normal(loc=means_vector_df.iloc[index], scale=stds_vector_df.iloc[index], size=int(number_of_samples_in_group))
        all_cpgs = all_cpgs.append(pd.DataFrame(cpg_i.tolist()).T)
        cpg_name = cpg_name.append(pd.Series(str("cpg") + str(index)), ignore_index=True)
    for index, row in all_cpgs.iterrows():
        for j, value in all_cpgs.iloc[index].items():
            if float(value) > float(1.0):
                all_cpgs.iloc[index, j] = np.nanmax(beta_matrix.iloc[index])
            elif float(value) < float(0):
                all_cpgs.iloc[index, j] = np.nanmin(beta_matrix.iloc[index])
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

    with ZipFile(os.getcwd()+dirname+workflow_num+str('.zip'), 'w') as zipObj:
        zipObj.write(os.getcwd()+dirname+file_name_simulated_data, file_name_simulated_data)
        zipObj.write(os.getcwd()+dirname+file_name_user_parameters, file_name_user_parameters)
        zipObj.write(os.getcwd()+dirname+file_name_truly_modified_indices, file_name_truly_modified_indices)

def multi_simMethyl(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs):
    combination_df = get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion, num_true_modified, user_specified_n_CpGs)
    list_of_workflows = [num for num in range(0, len(combination_df.index))]
    pool = Pool(processes=40) # user specified CPUs e.g., processes=8
    result = pool.map(simulator_pool,list_of_workflows)
    pool.close()
    #print(result)
    # Serial Version
    #for i in range(len(combination_df)):
    #    print("Running environmental setup using parameters:")
    #    print(combination_df.iloc[i])
    #    simulator(combination_df.loc[i, 'n_samples'], combination_df.loc[i, 'effect_size'],
    #              combination_df.loc[i, 'healthy_proportion'], combination_df.loc[i, 'n_true_modified_CpGs'],
    #              combination_df.loc[i, 'n_CpGs'], combination_df.loc[i, 'workflow'])

def simulator_pool(workflow):
    healthy_proportion = [0.5]
    num_true_modified = [50,100,150,200,350,500,950,1250]#[5,10,15,20,35,50,95,125]
    user_specified_n_CpGs = [371377]
    total_num_samples_vector = [50,100,200,350,500,650,800,950]
    effect_size_vector = [0.01,0.02,0.03,0.05,0.07,0.09,0.11,0.13]
    combination_df = get_all_combinations(total_num_samples_vector, effect_size_vector, healthy_proportion,
                                          num_true_modified, user_specified_n_CpGs)
    print("Running environmental setup using parameters:")
    print(workflow)
    simulator(combination_df.loc[workflow, 'n_samples'], combination_df.loc[workflow, 'effect_size'],
              combination_df.loc[workflow, 'healthy_proportion'], combination_df.loc[workflow, 'n_true_modified_CpGs'],
              combination_df.loc[workflow, 'n_CpGs'], combination_df.loc[workflow, 'workflow'])

if __name__ == '__main__':
    healthy_proportion = [0.5]
    num_true_modified = [50,100,150,200,350,500,950,1250]#[5,10,15,20,35,50,95,125]
    user_specified_n_CpGs = [371377]#[1000]
    total_num_samples_vector = [50,100,200,350,500,650,800,950]
    effect_size_vector = [0.01,0.02,0.03,0.05,0.07,0.09,0.11,0.13]
    multi_simMethyl(total_num_samples_vector, effect_size_vector, healthy_proportion,num_true_modified,user_specified_n_CpGs)
    print("Time taken: ",time.time() - starttime)