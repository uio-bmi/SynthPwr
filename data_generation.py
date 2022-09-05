import csv

import numpy as np
import pandas as pd
from random import random

import GEOparse
gse = GEOparse.get_GEO(geo="GSE145714", destdir="./")
beta_matrix = None
# Potential problem sources
#pwrEWAS uses variance over sd

# Function for calculating the mean of CpGs across all samples.
def calculate_mean(beta_matrix):
    df = pd.DataFrame(data=beta_matrix)
    cpg_means = df.mean(axis=1, skipna=True)
    return cpg_means

def calculate_sd(beta_matrix):
    df = pd.DataFrame(data=beta_matrix)
    cpg_stds = df.std(axis=0, skipna=True)
    return cpg_stds

# Returns a matrix of ID's for CpGs and their associated mean values to generate desirable number of CpGs in the simulated data.
set.seed(42)
def user_specified_num_elements(means, user_specified_n_CpGs):
    CpG_rows = [num for num in range(1, means.shape[0])]
    sampling_indices = random.choices(CpG_rows,user_specified_n_CpGs) #sampling with replacement
    print(sampling_indices)
    return(sampling_indices)


# Induce differences (i.e., hypermetalytation) between groups by changing mean beta-values of selected CpGs. The modified CpGs are sampled from mean values of the reference group (control group). The output then recordes indices of the truly modified CpGs and vector of means with introduced corresponding differences.
# This function operates based on indices between a control and case
def induce_group_differnces(num_true_modified, vector_of_ref_means, effect_size, user_path):
    set.seed(42)
    means_of_control_population = [num for num in range(1, len(vector_of_ref_means))]
    truly_different_indices = np.array(random.sample(means_of_control_population, num_true_modified)) #indicies to change

    #Save the truly changed positions to file
    with open('truly_different_sites_indices.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        [writer.writerow(r) for r in truly_different_indices]

    vector_of_affected_means = np.array(vector_of_ref_means) #Copy and replicate control with potential case
    vector_of_affected_means[truly_different_indices] = vector_of_ref_means[truly_different_indices] + effect_size #Induce the difference
    vector_of_bool_affected_means = [num for num in range(1, len(vector_of_affected_means))]
    for row in vector_of_affected_means:
        if vector_of_affected_means[row] > 1:
            vector_of_affected_means[row] = 1
            vector_of_bool_affected_means[row] = True
        else:
            vector_of_bool_affected_means[row] = False
    return vector_of_bool_affected_means

#get number of samples per group based on the control-case balance (healthy proportion)
def get_group_number(healthy_proportion, total_num_samples):
    group1_number_samples = healthy_proportion*total_num_samples
    group2_number_samples = total_num_samples-group1_number_samples
    groups = {
        "Group": ["Reference", "Diseased"],
        "Num": [group1_number_samples, group2_number_samples]
    }
    group_nums = pd.DataFrame(groups)
    return group_nums

#Simulate beta values for every CpG across all samples. Due to the range of the real beta values (0 and 1), truncate/adjust the values within the range
def sample_distribution(means_vector, sds_vector, number_of_samples_in_group):
    cpg_name = [None] * len(means_vector)
    all_cpgs = [None] * len(means_vector)
    for i in means_vector:
        cpg_i = np.random.normal(number_of_samples_in_group, means_vector[i], sds_vector[i])
        cpg_name[i] = "cpg"+str(i)
        all_cpgs = pd.concat([all_cpgs, cpg_i])

    for j in len(all_cpgs):
        if all_cpgs[j] > 1:
            all_cpgs[j] = max(beta_matrix)
        elif all_cpgs[j] < 0:
            all_cpgs[j] = min(beta_matrix)

    #rownames(all_cpgs) = cpg_name
    return all_cpgs

