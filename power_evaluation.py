import csv
import os
import zipfile
from zipfile import ZipFile

import scipy
from scipy import stats, optimize, interpolate
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import time
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, FactorVector
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import DataFrame, StrVector
from multiprocessing import Pool
from scipy.stats import ttest_ind

starttime = time.time()
pandas2ri.activate()
dirname = os.sep+"power_experiments"+os.sep
base = importr('base')
packnames = ('stats', 'caret')
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    robjects.r.options(download_file_method='curl')
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    utils.install_packages(StrVector(names_to_install))

stats = importr('stats')
caret = importr('caret')

def calculate_delta_beta(group1_means_vector,group2_means_vector):
    list_of_delta = [row for row in range(0, len(group1_means_vector.index))]
    for i in list_of_delta:
        list_of_delta[i] = np.abs(np.mean(group1_means_vector.iloc[i,:]) - np.mean(group2_means_vector.iloc[i,:]))
    return list_of_delta

def unzip_workflow(workflow_num):
    zf = zipfile.ZipFile(os.getcwd() + dirname + "SimMethyl_run_" + str(workflow_num) + ".zip")
    zip_directory = "SimMethyl_run_" + str(workflow_num)
    if os.path.isdir(os.getcwd() + dirname + zip_directory +os.sep) == False:
        os.mkdir(os.getcwd() + dirname + zip_directory + os.sep)
    zf.extractall(os.getcwd() + dirname + zip_directory + os.sep)
    simulated_data_df = pd.read_csv(os.getcwd() + dirname + zip_directory + os.sep + 'Simulated_data.csv', header=0, sep=',')
    txt_file = open(os.getcwd() + dirname + zip_directory + os.sep + 'truly_different_sites_indices.txt', "r")
    true_CpG_locations = txt_file.read().split()
    txt_file.close()
    user_params_df = pd.read_csv(os.getcwd() + dirname + zip_directory + os.sep + 'User_Parameters.csv', index_col=0)
    result = [simulated_data_df, true_CpG_locations, user_params_df]
    return result

def unzip_workflow_params_only(workflow_num):
    zf = zipfile.ZipFile(os.getcwd() + dirname + "SimMethyl_run_" + str(workflow_num) + ".zip")
    zip_directory = "SimMethyl_run_" + str(workflow_num)
    if os.path.isdir(os.getcwd() + dirname + zip_directory + os.sep) == False:
        os.mkdir(os.getcwd() + dirname + zip_directory + os.sep)
    zf.extractall(os.getcwd() + dirname + zip_directory + os.sep)
    user_params_df = pd.read_csv(os.getcwd() + dirname + zip_directory + os.sep + 'User_Parameters.csv', index_col=0)
    user_params_df = user_params_df.T
    return user_params_df

def get_all_zip(n_zip):
    all_enviornments_list = []
    for i in range(0, n_zip):
        run = unzip_workflow(i)
        all_enviornments_list.append(run)
    return all_enviornments_list

def extract_data_from_zip(all_runs_zip_num, workflow_num, file_name):
    file_name_to_extract = -1
    if file_name=="Simulated_data":
        file_name_to_extract = 0
    elif file_name=="truly_different_sites_indices":
        file_name_to_extract = 1
    elif file_name=="User_Parameters":
        file_name_to_extract = 2
    else:
        print("something very strange has occured")
    all_zip = get_all_zip(all_runs_zip_num)
    extracted_data = all_zip[workflow_num][file_name_to_extract]
    return extracted_data

#targets simulated data characteristics to determine num samples
def get_num_sample_per_group(user_params, simulated_df):
    healthy_proportion = user_params.iloc[2,1]
    n_Group1 = len(simulated_df.columns)*healthy_proportion
    n_Group2 = len(simulated_df.columns)*(1-healthy_proportion)
    result = [n_Group1,n_Group2]
    return result

def tradtional_t_test(simulated_df, n_group1):
    t_test_label = "Welch Two Sample t-test"
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        t_test_statistic, t_test_pvalue = ttest_ind(gr1.dropna(),gr2.dropna(), equal_var=True)
        all_p_val[i] = t_test_pvalue
    return all_p_val

def ks_test(simulated_df, n_group1):
    ks_test_label = "Two-sample Kolmogorov-Smirnov test"
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        k_test_statistic, k_test_pvalue = scipy.stats.ks_2samp(gr1.dropna(), gr2.dropna())
        #k_test_pvalue = dgof.ks_test(gr1,gr2,alternative = "two.sided", exact = True)
        all_p_val[i] = k_test_pvalue
    return all_p_val

def w_test(simulated_df, n_group1):
    w_test_label = "Wilcoxon rank-sum test"
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        w_test_statistic, w_test_pvalue = scipy.stats.ranksums(gr1.dropna(), gr2.dropna())
        all_p_val[i] = w_test_pvalue
    return all_p_val

def limma_test(simulated_df, n_group1, n_group2):
    n_group1 = int(n_group1)
    n_group2 = int(n_group2)

    Group1 = [1 for num in range(0, n_group1)] + ([0 for num in range(0, n_group2)])
    Group2 = [0 for num in range(0, n_group1)] + ([1 for num in range(0, n_group2)])
    design_matrix = pd.concat([pd.DataFrame(Group1).dropna(), pd.DataFrame(Group2).dropna()], axis=1)
    design_matrix.columns = ['Group1','Group2']

    robjects.r('''
            chooseCRANmirror(ind = 1)
            install.packages('limma', repos='https://cloud.r-project.org/')
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
    robjects.r('''
                BiocManager::install("ENmix", update=False)
                library(ENmix)
                ENmix_transform <- function(beta_matrix, verbose=FALSE) {
                m_matrix = as.data.frame(B2M(beta_matrix))
                }
                ''')
    logit_result = robjects.r['ENmix_transform']
    m_matrix = logit_result(beta_matrix)
    return m_matrix

#Running all the tests to compare the two groups for every CpG site based on either transformed or original values. The output is p-values for every CpG site and for every test
def run_all_test(n_group1, n_group2, user_params, sim_data, beta_or_M_string):
    if beta_or_M_string == "M":
        #print("M values")
        sim_data = pd.DataFrame(logit_transformation(sim_data))
        sim_data = sim_data.T
    elif beta_or_M_string == "beta":
        # print("B values")
        sim_data = sim_data
    else:
        print("invalid input for beta_or_M_string in run_all_test, use either M or beta")

    sim_data.fillna(0)
    t_test_output = tradtional_t_test(sim_data, n_group1)
    ks_test_output = ks_test(sim_data, n_group1)
    limma_test_output = limma_test(sim_data, n_group1, n_group2)
    w_test_output = w_test(sim_data, n_group1)
    # all_test_output = pd.concat([pd.concat([pd.DataFrame(t_test_output), pd.DataFrame(ks_test_output)], axis=1, ignore_index=True), pd.DataFrame(limma_test_output)], axis=1, ignore_index=True)

    test_dict = {"T-test": t_test_output,
                 "KS-test": ks_test_output,
                 "limma-test": limma_test_output,
                 "W-test": w_test_output}
    all_test_output = pd.DataFrame(test_dict)
    print("Multiple tests")
    print(all_test_output)
    return all_test_output

def getadjustPval(p_val_vector):
    method = "fdr"
    p_val_vector = stats.p_adjust(p_val_vector, method=method)
    return p_val_vector

def multitest_p_adjust(all_test_df):
    store_p_adjust = pd.DataFrame(index=all_test_df.index, columns=all_test_df.columns)
    for cpg in range(0, len(all_test_df.index)):
        for p_val_test in range(0, len(all_test_df.columns)):
            store_p_adjust.iloc[cpg,p_val_test] = getadjustPval(all_test_df.iloc[cpg, p_val_test])
    store_p_adjust.columns = ["T-test",
                              "KS-test",
                            "limma", "W-test"]
    return store_p_adjust

# prior to confusion matrix construction, CpGs should be labeled based on two distinction characteristics: whether they are truly modified, and if their adjusted p-value is less then 0.05
# The input is the vector with observed values within the simulated data. The output is a boolean which represents if the observation is signficantly different between the two groups (0) and not significant (1)
def create_predicted_vector(p_val_vector, group1_means_vector,group2_means_vector, p_cut):
    delta_beta_cpg = calculate_delta_beta(group1_means_vector,group2_means_vector)
    boolean_predicted_vector = [False for row in range(0, len(p_val_vector.index))]
    for i in range(0, len(delta_beta_cpg)):
        if delta_beta_cpg[i] >= 0.01 and p_val_vector.iloc[i] < 0.05:
            boolean_predicted_vector[i] = True
        elif delta_beta_cpg[i] >= 0.01 and p_val_vector.iloc[i] > 0.05:
            boolean_predicted_vector[i] = False
    return boolean_predicted_vector

def create_expected_val_vector(p_val_df, truly_different):
    truly_different = [int(numeric_string) for numeric_string in truly_different]
    indices_p_val_vector = [num for num in range(0, len(p_val_df.index))]
    boolean_truly_modifed = np.isin(indices_p_val_vector, truly_different, invert=False)
    return boolean_truly_modifed.tolist()

# Create a confusion matrix using observed and predicted values for every CpG. Create confusion matrix for every test. Could use PyCaret
def create_confusion_matrix(all_test_df, group1_means_vector,group2_means_vector,truly_different, p_cut):
    workflows_confusion_matrix = []
    expected = create_expected_val_vector(all_test_df, truly_different)
    for i in range(0, len(all_test_df.columns)):
        predicted = create_predicted_vector(all_test_df.iloc[:,i], group1_means_vector,group2_means_vector, p_cut)
        result = caret.confusionMatrix(data=FactorVector(robjects.BoolVector(predicted)),reference=FactorVector(robjects.BoolVector(expected)))
        workflows_confusion_matrix.append(result)
    return workflows_confusion_matrix

def calc_empirical_marg_power(workflows_confusion_matrix):
    df_all_test = pd.DataFrame(index=range(4), columns=range(2))
    #df_all_test = pd.DataFrame(index=range(4), columns=range(2))

    df_all_test.columns = ['Power', 'Test']
    df_all_test['Test'] = ['T_test',
                           'KS_test',
                            'Limma_test', 'W_test']

    for i in range(0, len(workflows_confusion_matrix)):
        confusion_matrix = workflows_confusion_matrix[i]
        #print(confusion_matrix[1])
        true_positive = confusion_matrix[1][1][1]
        false_positive = confusion_matrix[1][1][0]
        false_negative = confusion_matrix[1][0][1]
        true_negative = confusion_matrix[1][0][0]
        power_calc_val = int(true_positive)/(int(true_positive)+int(false_negative))
        #print("False Positive proportion: ",int(false_positive)/(int(false_positive)+int(true_positive)))
        df_all_test.iloc[i, 0] = power_calc_val
    #print(df_all_test)
    return df_all_test

def power_calc_multiple_runs(workflow):
    df_all_test = pd.DataFrame()
    print("Power Calculation for Workflow: ",workflow)
    workflow_files = unzip_workflow(workflow)
    sim_data = workflow_files[0]
    truly_different = workflow_files[1]
    user_params = workflow_files[2]

    n_samples_per_gr = get_num_sample_per_group(user_params, sim_data)
    n_group1 = n_samples_per_gr[0]
    n_group2 = n_samples_per_gr[1]

    # Defining default parameters inside temporarly due to Pool signature interaction
    beta_or_M_string = "M"
    p_cut = 0.05
    all_test = run_all_test(n_group1, n_group2, user_params,sim_data, beta_or_M_string)
    p_adjust_all_test = multitest_p_adjust(all_test)
    confusion_matrix = create_confusion_matrix(p_adjust_all_test,sim_data.iloc[:,0:int(n_group1)],sim_data.iloc[:,int(n_group1):], truly_different, p_cut)
    calc_power_value = pd.DataFrame(calc_empirical_marg_power(confusion_matrix))
    calc_power_value['ID'] = np.repeat(workflow, 4).tolist()#i
    df_all_test = df_all_test.append(calc_power_value)
    df_all_test.columns = ['Power', 'Test', 'ID']
    #print("for multiple test runs")
    #print(df_all_test.to_string())
    return df_all_test

def get_all_parameters(all_zips):
    environmental_tests_df = pd.DataFrame(index=range(0), columns=range(6))
    environmental_tests_df.columns = ['n_samples','n_CpGs','healthy_proportion','effect_size','n_modified_CpGs', 'ID']
    for i in range(0, all_zips):
        print("unpacking zip: ",i)
        user_params = unzip_workflow(i)[2]#extract_data_from_zip(all_zips,i,"User_Parameters")
        user_params = user_params.T
        environmental_tests_df.loc[len(environmental_tests_df.index)] = user_params.iloc[1].tolist() + [i]
    return environmental_tests_df

def merge_data(all_zips, p_adjust_method_string, beta_or_M_string, p_cut):
    pool = Pool(processes=64) # user specified CPUs e.g., processes=8
    list_of_workflows = [num for num in range(all_zips)]
    environmental_tests_df = pd.DataFrame(index=range(0), columns=range(6))
    environmental_tests_df.columns = ['n_samples', 'n_CpGs', 'healthy_proportion', 'effect_size', 'n_modified_CpGs','ID']
    result = pool.map(power_calc_multiple_runs, list_of_workflows)
    parameters = pool.map(unzip_workflow_params_only, list_of_workflows)
    pool.close()
    # Serial Version
    # power_calc_val = power_calc_multiple_runs(all_zips,p_adjust_method_string,beta_or_M_string,p_cut)
    #parameters = get_all_parameters(all_zips)
    output_of_map = pd.DataFrame()
    for value in result:
        output_of_map = output_of_map.append(value, ignore_index=True)
    i = 0
    for param_combo in parameters:
        environmental_tests_df.loc[len(environmental_tests_df.index)] = param_combo.iloc[1].tolist() + [i]
        i += 1
    #print(output_of_map)
    output_df = output_of_map.merge(environmental_tests_df, how='inner', on='ID')
    print("The output for analysis: ")
    print(output_df.to_string())
    return output_df

def plot_line(data_df, y_parameter_string, varied_parameter_string, p_adjust_method_string, beta_or_M_string):
    selected_parameters_df = data_df.loc[:,['Power','Test',varied_parameter_string,y_parameter_string]]
    y_parameter = y_parameter_string
    varied_parameter = varied_parameter_string
    for idx, val in enumerate(selected_parameters_df[varied_parameter].unique()):
        plt.figure()
        selected_parameters = selected_parameters_df.loc[selected_parameters_df[varied_parameter] == val,:]
        varied_param_i = val
        #print("One unique value of varied parameter: ", val)
        #print(selected_parameters)
        plt.title("Sample Vs Power ("+varied_parameter+" = "+str(val)+")")
        plt.xlabel('Sample')
        plt.ylabel('Power')
        line_plot_fig = sns.lineplot(data=selected_parameters, x=y_parameter, y='Power', hue='Test', errorbar=None)
        fig = line_plot_fig.get_figure()
        fig.savefig(str("LinePlot_") + str("vs_Power_") + varied_parameter + str("_step_") +
                    str(varied_param_i) + "_" + p_adjust_method_string + "_" + beta_or_M_string + ".png", dpi=100)
        #scatter_fig = sns.scatterplot(data=selected_parameters_df, x=y_parameter, y=selected_parameters_df['Power'], hue=selected_parameters_df['Test'])

def plot_every_effect_line(data_df, nsample_string, effect_string, p_adjust_method_string, beta_or_M_string):
    selected_parameters_df = data_df.loc[:,['Power','Test',effect_string,nsample_string]]
    y_parameter = nsample_string
    for idx, val in enumerate(selected_parameters_df['Test'].unique()):
        plt.figure()
        selected_parameters = selected_parameters_df.loc[selected_parameters_df['Test'] == val,:]
        varied_param_i = val
        print("One unique value of varied parameter: ", val)
        print(selected_parameters.to_string())
        plt.title("Sample Vs Power (all effect sizes for "+val+")")
        plt.xlabel('Sample')
        plt.ylabel('Power')
        #print("extra testing")
        #print(selected_parameters['effect_size'])
        line_plot_fig = sns.lineplot(data=selected_parameters, x=y_parameter, y='Power', hue='effect_size', errorbar=None, legend='full')
        fig = line_plot_fig.get_figure()
        fig.savefig(str("LinePlot_") + str("vs_Power_") + "alleffects" + str("_steptest_") +
                    str(varied_param_i) + "_" + p_adjust_method_string + "_" + beta_or_M_string + ".png", dpi=100)

def heat_map(data_df, x_parameter_string, y_parameter_string, p_adjust_method_string, test_string, beta_or_M_string):
    plt.figure()
    selected_parameters_df = data_df.loc[:, ['Power', 'Test', x_parameter_string, y_parameter_string]]
    y_parameter = y_parameter_string
    x_parameter = x_parameter_string
    selected_parameters_df = selected_parameters_df.loc[selected_parameters_df['Test'] == test_string,:]
    plt.figure()
    df_m = selected_parameters_df.reset_index().pivot_table(index=y_parameter,columns=x_parameter, values='Power')
    #print("heat map 2D grid: ", df_m)
    df_m = df_m.fillna(0)
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.title("HeatMap of "+y_parameter+ " and "+ x_parameter)
    plt.xlabel(x_parameter)
    plt.ylabel('Sample')
    ax = sns.heatmap(df_m, cmap='crest', annot=True)
    fig = ax.get_figure()
    fig.savefig(str("HeatMap_") + y_parameter + "_vs_" +x_parameter +"_fill_Power_"+"test_"+test_string+"_padjust_"+p_adjust_method_string+"_"+beta_or_M_string + ".png", dpi=100)

#Aggregate function that calculates marginal power for all the tests across samples based on user-defined inputs. As an output it generates line plots and heat maps that shows the relationship between power and two varying parameters and three constant parameters
def PowerCalc(num_zips, p_adjust_method_string, beta_or_M_string, y_parameter_string, x_parameter_string, color_string, test_vector_string,p_cut):
    all_test_all_zips = merge_data(num_zips,p_adjust_method_string,beta_or_M_string, p_cut)
    all_test_all_zips.to_csv(os.getcwd()+dirname+"PowerCalc_"+beta_or_M_string+"_"+p_adjust_method_string+'.csv', sep=",", index='ID',header=True)

    plot_line(all_test_all_zips, y_parameter_string, x_parameter_string,p_adjust_method_string,beta_or_M_string)

    for test in test_vector_string:
        heat_map(all_test_all_zips, x_parameter_string, y_parameter_string, p_adjust_method_string,test,beta_or_M_string)
    plot_every_effect_line(all_test_all_zips, "n_samples", "effect_size",p_adjust_method_string,beta_or_M_string)

if __name__ == '__main__':
    num_zips = 16 # number of simulated data environmental workflows
    p_adjust_method_string = "fdr" #"none"/"BH"/"bonferroni"/"fdr"
    beta_or_M_string = "M" #"beta"/"M"
    y_parameter_string = "n_samples" # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
    x_parameter_string = "effect_size" # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
    color_string = "Power"
    test_vector_string = ["T_test","Limma_test",
                          "KS_test",
                          "W_test"]
    p_cut = 0.05 #0.01, 0.1, 0.05

    PowerCalc(num_zips,p_adjust_method_string,beta_or_M_string,y_parameter_string,x_parameter_string,color_string,test_vector_string,p_cut)
    endtime = time.time()
    print("Time taken: ", endtime - starttime)
