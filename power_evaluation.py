import os
import random

import scipy
from scipy import stats, optimize, interpolate
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import time
import matplotlib.pyplot as plt
#matplotlib.use('TkAgg')
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, FactorVector
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import DataFrame, StrVector
from multiprocessing import Pool
from scipy.stats import ttest_ind, gaussian_kde

starttime = time.time()
pandas2ri.activate()
dirname = os.sep+"power_experiments"+os.sep
figuredirname = os.sep+"figures"+os.sep
summarycalcdirname = os.sep+"summary_stats"+os.sep
base = importr('base')
packnames = ('stats', 'caret', 'ENmix')
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    robjects.r.options(download_file_method='curl')
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    utils.install_packages(StrVector(names_to_install))

stats = importr('stats')
caret = importr('caret')
ENmix = importr('ENmix')

def calculate_delta_beta(group1_means_vector,group2_means_vector):
    """
    Measure the comparative difference between the mean methylation signatures of two groups
    :param group1_means_vector: the methylation values of group1
    :param group1_means_vector: the methylation values of group2
    :rtype list_of_delta: the delta_beta for all CpGs
    """
    list_of_delta = [row for row in range(0, len(group1_means_vector.index))]
    list_of_delta_scores = []
    for i in list_of_delta:
        list_of_delta_scores.append(np.abs((np.nanmean(group1_means_vector.iloc[i,:]) - np.nanmean(group2_means_vector.iloc[i,:]))))
    return list_of_delta

def unpack_workflow(workflow_num):
    """
    Unpack the resources pertaining to a single environmental setup (i.e., workflow)
    :param workflow_num: the environmental identifier amongst all combinatorial environmental setups
    :rtype result: the simulated_data, ground truth CpG locations modified and environmental parameters for a workflow as a list
    """
    list_of_simulated_data_df = np.load(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"Simulated_data"+"_workflow_"+str(workflow_num)+".npz")
    list_of_true_CpG_locations = np.load(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"truly_different_sites_indices"+"_workflow_"+str(workflow_num)+".npz")
    user_params_df = pd.read_csv(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"User_Parameters"+"_workflow_"+str(workflow_num)+".csv", index_col=0)
    result = [list_of_simulated_data_df, list_of_true_CpG_locations, user_params_df]
    return result

def unpack_workflow_params_only(workflow_num):
    """
    Unpack only the environmental parameters for a single environmental setup (i.e., workflow)
    :param workflow_num: the environmental identifier amongst all combinatorial environmental setups
    :rtype result: the environmental parameters for a workflow
    """
    user_params_df = pd.read_csv(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"User_Parameters"+"_workflow_"+str(workflow_num)+".csv", index_col=0)
    return user_params_df

#targets simulated data characteristics to determine num samples
def get_num_sample_per_group(user_params, simulated_df):
    """
    Retreive the split in sample groups based on simulated data and environmental parameters
    :param user_params: user-specified environmental parameters
    :param simulated_df: synthetic dataset of Sample by CpG methylation values
    :rtype result: the number representing how many samples within each group
    """
    healthy_proportion = user_params.iloc[2,1]
    n_Group1 = len(simulated_df.columns)*healthy_proportion
    n_Group2 = len(simulated_df.columns)*(1-healthy_proportion)
    result = [n_Group1,n_Group2]
    return result

def tradtional_t_test(simulated_df, n_group1):
    """
    Perform a student's t-test on a simulated dataset
    :param simulated_df: synthetic dataset of Sample by CpG methylation values
    :param n_group1: number of samples in group1
    :rtype result: a p-value for every CpG in the simulated_data
    """
    simulated_df = pd.DataFrame(logit_transformation(simulated_df))
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]

    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        t_test_statistic, t_test_pvalue = ttest_ind(gr1,gr2, equal_var=True)
        all_p_val[i] = t_test_pvalue
    return all_p_val

def w_test(simulated_df, n_group1):
    """
    Perform a non-parametric Wilcox rank sum test on a simulated dataset
    :param simulated_df: synthetic dataset of Sample by CpG methylation values
    :param n_group1: number of samples in group1
    :rtype result: a p-value for every CpG in the simulated_data
    """
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        w_test_statistic, w_test_pvalue = scipy.stats.ranksums(gr1, gr2)
        all_p_val[i] = w_test_pvalue
    return all_p_val

def limma_test(simulated_df, n_group1, n_group2):
    """
    Fit a linear model using R's limma package on a simulated dataset
    :param simulated_df: synthetic dataset of Sample by CpG length methylation values
    :param n_group1: number of samples in group1
    :param n_group2: number of samples in group2
    :rtype result: a p-value for every CpG in the simulated_data
    """
    simulated_df = pd.DataFrame(logit_transformation(simulated_df))
    n_group1 = int(n_group1)
    n_group2 = int(n_group2)
    Group1 = [1 for num in range(0, n_group1)] + ([0 for num in range(0, n_group2)])
    Group2 = [0 for num in range(0, n_group1)] + ([1 for num in range(0, n_group2)])
    design_matrix = pd.concat([pd.DataFrame(Group1), pd.DataFrame(Group2)], axis=1)
    design_matrix.columns = ['Group1','Group2']
    robjects.r('''
            library(limma)
            limma_test <- function(simulated_df, design_matrix, verbose=FALSE) {
            contrast_matrix = makeContrasts(Diff = Group2 - Group1, levels = design_matrix)
            fit1 <- lmFit(simulated_df, design_matrix)
            fit2 <- contrasts.fit(fit1, contrast_matrix)
            fit3 <- eBayes(fit2)
            limma_output = as.data.frame(fit3)
            limma_Pvalues = as.data.frame(limma_output$p.value)
            }
            ''')
    limma_result = robjects.r['limma_test']
    all_p_val = limma_result(simulated_df, design_matrix)
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_from_r_df = ro.conversion.rpy2py(all_p_val)
    return pd_from_r_df.iloc[:,0].tolist()

def logit_transformation(beta_matrix):
    """
    Perform a logit transformation using R's ENmix package to normalise beta-methylation values into M-methylation values to satisfy statistical assumptions
    :param beta_matrix: a 2D matrix of decimal values between 0 and 1
    :rtype logit_result: the converted matrix as M-methylation values
    """
    logit_result = ENmix.B2M(beta_matrix)
    return logit_result

#Running all the tests to compare the two groups for every CpG site based on either transformed or original values. The output is p-values for every CpG site and for every test
def run_all_test(n_group1, n_group2, sim_data):
    """
    Perform all hypothesis testing on the simulated_data
    :param n_group1: number of samples in group1
    :param n_group2: number of samples in group2
    :param sim_data: synthetic dataset of Sample by CpG length methylation values
    :rtype logit_result: a dict containing the results of all statistical testing
    """
    t_test_output = tradtional_t_test(sim_data, n_group1)
    limma_test_output = limma_test(sim_data, n_group1, n_group2)
    w_test_output = w_test(sim_data, n_group1)

    test_dict = {"T-test": t_test_output,"limma-test": limma_test_output,"W-test": w_test_output}
    all_test_output = pd.DataFrame(test_dict)
    return all_test_output

def getadjustPval(p_val_vector):
    """
    Perform a single p-value test adjustment
    :param p_val_vector: p-values to be adjusted
    :rtype p_val_vector: p-values that have been corrected
    """
    method = "fdr"
    p_val_vector = stats.p_adjust(p_val_vector, method=method)
    return p_val_vector

def multitest_p_adjust(all_test_df):
    """
    Perform multiple testing correction (i.e., adjust p-values) that were derived from multiple statistical tests to control for increases in false positives
    :param all_test_df: a dataframe containing p-values which is of the shape of tests being done and by CpG length
    :rtype store_p_adjust: adjusted p-values corrected for each test
    """
    store_p_adjust = pd.DataFrame(index=all_test_df.index, columns=all_test_df.columns)
    for cpg in range(0, len(all_test_df.index)):
        for p_val_test in range(0, len(all_test_df.columns)):
            store_p_adjust.iloc[cpg,p_val_test] = getadjustPval(all_test_df.iloc[cpg, p_val_test])
    store_p_adjust.columns = ["T-test", "limma", "W-test"]
    return store_p_adjust

def precompute_average_dataframe(list_of_data_df, num_workflows, num_simulations, num_tests):
    average_df = pd.DataFrame(index=range(num_workflows * num_tests), columns=range(2))
    average_df.columns = ['Power', 'Power_se']
    for i in range(0, num_workflows * num_tests):
        power_scores = []
        for sim_iter in range(0, num_simulations):
            power_scores.append(pd.DataFrame(list_of_data_df[sim_iter]).iloc[i, 0])
        row_mean = np.mean(power_scores)
        row_se = scipy.stats.sem(power_scores)
        average_df.iloc[i, 0] = row_mean
        average_df.iloc[i, 1] = row_se
    return average_df

def plot_line(data_df, num_workflows, num_simulations, num_tests, y_parameter_string, varied_parameter_string, p_adjust_method_string):
    average_df = precompute_average_dataframe(data_df, num_workflows, num_simulations, num_tests)
    data_df = data_df[0]
    data_df['Power'] = average_df['Power']
    data_df['Power_se'] = average_df['Power_se']
    selected_parameters_df = data_df.loc[:, ['Power', 'Power_se','Test', varied_parameter_string, y_parameter_string]]
    y_parameter = y_parameter_string
    varied_parameter = varied_parameter_string
    for idx, val in enumerate(selected_parameters_df[varied_parameter].unique()):
        plt.figure()
        selected_parameters = selected_parameters_df.loc[selected_parameters_df[varied_parameter] == val,:]
        varied_param_i = val
        plt.title("Sample Vs Power ("+varied_parameter+" = "+str(val)+")")
        plt.xlabel('Sample')
        plt.ylabel('Power')
        line_plot_fig = sns.lineplot(data=selected_parameters, x=y_parameter, y='Power', hue='Test', errorbar=None)
        plt.errorbar(x=selected_parameters['n_samples'], y=selected_parameters['Power'],yerr=selected_parameters['Power_se'], fmt='none', c='black', capsize=2)
        fig = line_plot_fig.get_figure()
        fig.savefig(os.getcwd()+figuredirname+"LinePlot_" + "vs_Power_" + varied_parameter + "_step_" +str(varied_param_i) + "_" + p_adjust_method_string + ".png", dpi=100)

def plot_every_effect_line(data_df, num_workflows, num_simulations, num_tests, nsample_string, effect_string, p_adjust_method_string):
    average_df = precompute_average_dataframe(data_df, num_workflows, num_simulations, num_tests)
    data_df = data_df[0]
    data_df['Power'] = average_df['Power']
    data_df['Power_se'] = average_df['Power_se']
    selected_parameters_df = data_df.loc[:, ['Power', 'Power_se','Test',effect_string,nsample_string]]
    y_parameter = nsample_string
    for idx, val in enumerate(selected_parameters_df['Test'].unique()):
        plt.figure()
        selected_parameters = selected_parameters_df.loc[selected_parameters_df['Test'] == val,:]
        varied_param_i = val
        plt.title("Sample Vs Power (all effect sizes for "+val+")")
        plt.xlabel('Sample')
        plt.ylabel('Power')
        plt.axhline(0.8, linestyle="--", c="black", alpha=.75)
        line_plot_fig = sns.lineplot(data=selected_parameters, x=y_parameter, y='Power', hue='effect_size', errorbar=None, legend='full')
        plt.errorbar(x=selected_parameters['n_samples'], y=selected_parameters['Power'],yerr=selected_parameters['Power_se'], fmt='none', c='black', capsize=2)

        fig = line_plot_fig.get_figure()
        fig.savefig(os.getcwd()+figuredirname+"LinePlot_" + "vs_Power_" + "alleffectsizes" + "_test_" +str(varied_param_i) + "_" + p_adjust_method_string + ".png", dpi=100)

def heat_map(data_df, num_workflows, num_simulations, num_tests, x_parameter_string, y_parameter_string, p_adjust_method_string, test_string):
    plt.figure()
    average_df = pd.DataFrame(index=range(num_workflows*num_tests), columns=range(1))
    average_df.columns = ['Power']
    for i in range(0, num_workflows*num_tests):
        power_scores = []
        for sim_iter in range(0, num_simulations):
            power_scores.append(pd.DataFrame(data_df[sim_iter]).iloc[i,0])
        row_mean = np.mean(power_scores)
        average_df.iloc[i,0] = row_mean
    data_df = data_df[0]
    data_df['Power'] = average_df['Power']

    selected_parameters_df = data_df.loc[:, ['Power', 'Test', x_parameter_string, y_parameter_string]]
    y_parameter = y_parameter_string
    x_parameter = x_parameter_string
    selected_parameters_df = selected_parameters_df.loc[selected_parameters_df['Test'] == test_string,:]
    selected_parameters_df.to_csv(os.getcwd() + summarycalcdirname + "PowerCalc_avg_" + test_string + '.csv', sep=",", index='ID',header=True)
    df_m = selected_parameters_df.reset_index().pivot_table(index=y_parameter,columns=x_parameter, values='Power')
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.title("HeatMap of "+y_parameter+ " and "+ x_parameter)
    plt.xlabel(x_parameter)
    plt.ylabel('Sample')

    ax = sns.heatmap(df_m, cmap='crest', annot=True)
    fig = ax.get_figure()
    fig.savefig(os.getcwd()+figuredirname+"HeatMap_" + y_parameter + "_vs_" +x_parameter +"_fill_Power_"+"test_"+test_string+"_padjust_"+p_adjust_method_string+".png", dpi=100)

def plot_power_by_sample(data_df, targetDelta):
    plt.figure()
    plt.title("Sample Vs Power (all effect sizes)")
    plt.xlabel('Sample')
    plt.ylabel('Power')
    plt.axhline(0.8, linestyle="--", c="black", alpha=.75)
    for effect_size in targetDelta:
        plt.plot(data_df['Sample_Size'], data_df[effect_size], label = str(effect_size))
    #plt.errorbar(x=selected_parameters['n_samples'], y=selected_parameters['Power'],yerr=selected_parameters['Power_se'], fmt='none', c='black', capsize=2)
    plt.legend()
    plt.show()
    plt.savefig(os.getcwd()+figuredirname+"LinePlot_" + "Sample_vs_Power_" + "alleffectsizes.png", dpi=100)

def plot_power_simulations_by_sample(list_of_data_df, targetDeltas, sample_steps, classicalPower):
    average_df = pd.DataFrame(columns=targetDeltas)
    errors_df = pd.DataFrame(columns=targetDeltas)
    zscore_df = pd.DataFrame(columns=targetDeltas)
    for delta_row in range(0, len(list_of_data_df[0].index)): # for every deltaBeta
        first = []
        second = []
        third = []
        for sim_iter in range(0, len(list_of_data_df)):  # for every sim
            first.append(list_of_data_df[sim_iter].iloc[delta_row, 0])
            second.append(list_of_data_df[sim_iter].iloc[delta_row, 1])
            third.append(list_of_data_df[sim_iter].iloc[delta_row, 2])
        avg_scores = [np.mean(first), np.mean(second), np.mean(third)]
        error_scores = [scipy.stats.sem(first), scipy.stats.sem(second), scipy.stats.sem(third)]
        average_df.loc[len(average_df)] = avg_scores
        errors_df.loc[len(errors_df)] = error_scores
    average_df.insert(len(average_df.columns), column="Sample_Size", value=sample_steps)
    errors_df.insert(len(errors_df.columns), column="Sample_Size", value=sample_steps)
    zscore_df.insert(len(zscore_df.columns), column="Sample_Size", value=sample_steps)
    plt.figure()
    plt.title("Sample Vs Power (all effect sizes) with errorbars")
    plt.xlabel('Sample')
    plt.ylabel('Power')
    plt.axhline(0.8, linestyle="--", c="black", alpha=.75)
    for effect_size in targetDeltas:
        plt.plot(average_df['Sample_Size'], average_df[effect_size], linewidth=2.0,label=str(effect_size)+" Empirical")
        plt.errorbar(x=average_df['Sample_Size'], y=average_df[effect_size],yerr=errors_df[effect_size], fmt='none', c='black', capsize=2)
        plt.plot(classicalPower['Sample_Size'], classicalPower[effect_size], linewidth=1.0, linestyle='dashed', label=str(effect_size)+" Classic")
        zscore_df[effect_size] = scipy.stats.zscore(average_df[effect_size], axis=0, ddof=0, nan_policy='omit')
    plt.legend()
    plt.savefig(os.getcwd() + figuredirname + "LinePlot_" + "Sample_vs_Power_" + "alleffectsizes" + "_errorbars.png",dpi=100)
    plt.show()
    average_df.to_csv(os.getcwd() + figuredirname + 'meanPower.csv')
    errors_df.to_csv(os.getcwd() + figuredirname + 'semError.csv')
    sns.displot(zscore_df)

if __name__ == '__main__':
    endtime = time.time()
    print("Time taken: ", endtime - starttime)
