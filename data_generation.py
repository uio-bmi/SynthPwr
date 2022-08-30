import csv
from random import random

import GEOparse
gse = GEOparse.get_GEO(geo="GSE145714", destdir="./")

#pwrEWAS uses variance over sd

# Returns a dataframe with shape parameters, column 1 mean, column 2 sd, rows are cpgs
def calculate_shape_parameters():
    calculate_mean() # input beta matrix, output vector of mean values
    calculate_sd() # input beta matrix, output vector of sd values

def calculate_mean():
    return 2

def calculate_sd():
    return 2

# Returns a matrix of ID's for CpGs and their associated mean values to generate desirable number of CpGs in the simulated data.
set.seed(42)
def user_specified_num_elements(means, user_specified_n_CpGs):
    CpG_rows = [num for num in range(1, means.shape[0])]
    sampling_indices = random.choices(CpG_rows,user_specified_n_CpGs)
    print(sampling_indices)
    return(sampling_indices)


# Induce differences between groups by changing mean beta-values of selected CpGs. The modified CpGs are sampled from mean values of the reference group (control group). The output is the recorded indices of the truly modified CpGs and vector of means with introduced corresponding differences.
# This function operates based on indices between a control and case
# ask knut about how to perserve the oriinginal control with keeping updates
#a = np.array([1,2,5,4,6,8,4,5])
#b = np.array([1,3,5])
#a[b]
#array([2, 4, 8])
#a[b] = a[b]+10
#a[b]
#array([12, 14, 18])

#Expected 11 2 15 4 16 8 4 5

#hypermetalytation
def induce_group_differnces(num_true_modified, vector_of_ref_means, effect_size, user_path):
    set.seed(42)
    means_of_control_population = [num for num in range(1, len(vector_of_ref_means))]
    truly_different_indices = random.sample(means_of_control_population, num_true_modified)
    filename = '/truly_different_sites_indices.txt'

    #Save the truly changed positions
    with open('test.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        [writer.writerow(r) for r in truly_different_indices]

    vector_of_affected_means = vector_of_ref_means #Copy and replicate control with potential case
    vector_of_affected_means[truly_different_indices] = vector_of_ref_means[truly_different_indices] + effect_size #Induce the difference
    if vector_of_affected_means > 1:
        vector_of_affected_means = 1
    result_list = [vector_of_affected_means]
    result_list.append(truly_different_indices)
    return result_list