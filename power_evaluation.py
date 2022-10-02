import csv
import os
import zipfile
from re import sub
from zipfile import ZipFile

import scipy
from scipy import stats, optimize, interpolate
import pingouin as pg
import pandas as pd
import math
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import rpy2
import rpy2.robjects as robjects
from numpy import log2, float64
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector, FactorVector
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import DataFrame, StrVector
from scipy.stats import ttest_ind

pandas2ri.activate()
dirname = "/experiments/"
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
    #if workflow_num == -1:
    #    workflow_num = 0
    zf = zipfile.ZipFile(os.getcwd()+dirname+"SimMethyl_run_"+str(workflow_num)+".zip")
    zip_directory = "SimMethyl_run_"+str(workflow_num)
    if os.path.isdir(os.getcwd() + dirname+zip_directory+"/") == False:
        os.mkdir(os.getcwd() + dirname+zip_directory+"/")
    zf.extractall(os.getcwd() + dirname + zip_directory+"/")
    simulated_data_df = pd.read_csv(os.getcwd() + dirname + zip_directory+"/"+'Simulated_data.txt', index_col=0, header=0, sep=' ')
    txt_file = open(os.getcwd()+dirname+ zip_directory+"/"+'truly_different_sites_indices.txt', "r")
    true_CpG_locations = txt_file.read().split()
    txt_file.close()
    user_params_df = pd.read_csv(os.getcwd()+dirname+ zip_directory+"/"+'User_Parameters.csv', index_col=0)
    result = [simulated_data_df, true_CpG_locations, user_params_df]
    return result

def get_all_zip(n_zip):
    all_enviornments_list = []
    for i in range(0, n_zip):
        run = unzip_workflow(i)
        all_enviornments_list.append(run)
    return all_enviornments_list

#not only simulated data, based on the permater of text file name, give associated data.
#get targeted experiment results
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
    #print("A",all_zip)
    extracted_data = all_zip[workflow_num][file_name_to_extract]
    #print("Extracting content: ",extracted_data)
    return extracted_data

#targets simulated data characteristics to determine num samples
def get_num_sample_per_group(all_zips, workflow_num, simulated_df):
    user_parameters = extract_data_from_zip(all_zips, workflow_num, "User_Parameters")
    healthy_proportion = user_parameters.iloc[2,1]
    n_Group1 = len(simulated_df.columns)*healthy_proportion
    n_Group2 = len(simulated_df.columns)*(1-healthy_proportion)
    result = [n_Group1,n_Group2]
    return result

def tradtional_t_test(simulated_df, n_group1):
    #user minifi
    t_test_label = "Welch Two Sample t-test"
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        #t_test_pvalue = pg.ttest(gr1, gr2, correction=False)
        t_test_statistic, t_test_pvalue = ttest_ind(gr1.dropna(),gr2.dropna(), equal_var=True)
        all_p_val[i] = t_test_pvalue#t_test_pvalue['p-val'].values[0]
    return all_p_val

def ks_test(simulated_df, n_group1):
    ks_test_label = "Two-sample Kolmogorov-Smirnov test"
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        k_test_statistic, k_test_pvalue = scipy.stats.ks_2samp(gr1.dropna(), gr2.dropna())
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
                BiocManager::install("ENmix", update=FALSE)
                library(ENmix)
                ENmix_transform <- function(beta_matrix, verbose=FALSE) {
                m_matrix = as.data.frame(B2M(beta_matrix))
                }
                ''')
    logit_result = robjects.r['ENmix_transform']
    m_matrix = logit_result(beta_matrix)
    return m_matrix

#Running all the tests to compare the two groups for every CpG site based on either transformed or original values. The output is p-values for every CpG site and for every test
def run_all_test(all_zips, workflow_num, beta_or_M_string):
    simulated_data_frame = extract_data_from_zip(all_zips,workflow_num,"Simulated_data")
    #print(simulated_data_frame)
    if beta_or_M_string == "M":
        #print("M values")
        simulated_data_frame = pd.DataFrame(logit_transformation(simulated_data_frame))
        simulated_data_frame = simulated_data_frame.T
    elif beta_or_M_string == "beta":
        #print("B values")
        simulated_data_frame = simulated_data_frame
    else:
        print("invalid input for beta_or_M_string in run_all_test, use either M or beta")

    # Array probe design

    simulated_data_frame.fillna(0)
    n_samples_per_gr = get_num_sample_per_group(all_zips, workflow_num, simulated_data_frame)
    n_group1 = n_samples_per_gr[0]
    n_group2 = n_samples_per_gr[1]
    t_test_output = tradtional_t_test(simulated_data_frame, n_group1)
    ks_test_output = ks_test(simulated_data_frame, n_group1)
    limma_test_output = limma_test(simulated_data_frame,n_group1,n_group2)
    w_test_output = w_test(simulated_data_frame, n_group1)
    #all_test_output = pd.concat([pd.concat([pd.DataFrame(t_test_output), pd.DataFrame(ks_test_output)], axis=1, ignore_index=True), pd.DataFrame(limma_test_output)], axis=1, ignore_index=True)

    test_dict = {"T-test": t_test_output, "KS-test": ks_test_output,"limma-test": limma_test_output, "W-test": w_test_output}
    all_test_output = pd.DataFrame(test_dict)
    return all_test_output

def getadjustPval(p_val_vector, method):
    p_val_vector = stats.p_adjust(p_val_vector, method=method)
    return p_val_vector

def multitest_p_adjust(all_test_df, p_adjust_method_string):
    store_p_adjust = pd.DataFrame(index=all_test_df.index, columns=all_test_df.columns)

    for cpg in range(0, len(all_test_df.index)):
        for p_val_test in range(0, len(all_test_df.columns)):
            store_p_adjust.iloc[cpg,p_val_test] = getadjustPval(all_test_df.iloc[cpg, p_val_test], p_adjust_method_string)
    store_p_adjust.columns = ["T-test", "KS-test", "limma", "W-test"]
    return store_p_adjust

# prior to confusion matrix construction, CpGs should be labeled based on two distinction characteristics: whether they are truly modified, and if their adjusted p-value is less then 0.05
# The input is the vector with observed values within the simulated data. The output is a boolean which represents if the observation is signficantly different between the two groups (0) and not significant (1)
def create_predicted_vector(p_val_vector, group1_means_vector,group2_means_vector, p_cut):
    bool_signif_p_val = calculate_delta_beta(group1_means_vector,group2_means_vector)
    for i in range(0, len(bool_signif_p_val)):
        print("sfias")
        print(bool_signif_p_val[i])
        print(p_val_vector.iloc[i])

        if bool_signif_p_val[i] >= 0.01 and p_val_vector.iloc[i] < 0.05:
            bool_signif_p_val[i] = True
        elif bool_signif_p_val[i] >= 0.01 and p_val_vector.iloc[i] > 0.05:
            bool_signif_p_val[i] = False
    return bool_signif_p_val

def create_expected_val_vector(p_val_df, all_runs_zip_num,workflow_num):
    indices_truly_different = extract_data_from_zip(all_runs_zip_num, workflow_num,"truly_different_sites_indices")
    #print(indices_truly_different)
    indices_truly_different = [int(numeric_string) for numeric_string in indices_truly_different]
    indices_p_val_vector = [num for num in range(0, len(p_val_df.index))]
    boolean_truly_modifed = np.isin(indices_p_val_vector, indices_truly_different, invert=False)
    return boolean_truly_modifed.tolist()

# Create a confusion matrix using observed and predicted values for every CpG. Create confusion matrix for every test. Could use PyCaret
def create_confusion_matrix(all_test_df, group1_means_vector,group2_means_vector,all_runs_zip_num, workflow_num, p_cut):
    workflows_confusion_matrix = []
    expected = create_expected_val_vector(all_test_df, all_runs_zip_num, workflow_num)
    for i in range(0, len(all_test_df.columns)):
        predicted = create_predicted_vector(all_test_df.iloc[:,i], group1_means_vector,group2_means_vector, p_cut)
        result = caret.confusionMatrix(data=FactorVector(robjects.BoolVector(predicted)),reference=FactorVector(robjects.BoolVector(expected)))
        workflows_confusion_matrix.append(result)
    return workflows_confusion_matrix

def calc_empirical_marg_power(workflows_confusion_matrix):
    df_all_test = pd.DataFrame(index=range(4), columns=range(2))
    #df_all_test = pd.DataFrame()
    #df_all_test.columns = ['Power']
    df_all_test.columns = ['Power', 'Test']
    df_all_test['Test'] = ['T_test', 'KS_test', 'Limma_test', 'W_test']

    for i in range(0, len(workflows_confusion_matrix)):
        confusion_matrix = workflows_confusion_matrix[i]
        print(confusion_matrix[1])
        true_positive = confusion_matrix[1][1][1]
        false_positive = confusion_matrix[1][1][0]
        false_negative = confusion_matrix[1][0][1]
        true_negative = confusion_matrix[1][0][0]
        power_calc_val = int(true_positive)/(int(true_positive)+int(false_negative))
        #df_all_test = df_all_test.append(pd.Series(power_calc_val), ignore_index=True)
        df_all_test.iloc[i, 0] = power_calc_val
    print(df_all_test)
    return df_all_test

def power_calc_multiple_runs(all_zips, p_adjust_method_string, beta_or_M_string, p_cut):
    df_all_test = pd.DataFrame()

    for i in range(0, all_zips):
        sim_data = extract_data_from_zip(all_zips, i, "Simulated_data")
        n_samples_per_gr = get_num_sample_per_group(all_zips, i, sim_data)
        n_group1 = n_samples_per_gr[0]
        n_group2 = n_samples_per_gr[1]

        all_test = run_all_test(all_zips,i,beta_or_M_string)
        p_adjust_all_test = multitest_p_adjust(all_test,p_adjust_method_string)
        confusion_matrix = create_confusion_matrix(p_adjust_all_test,sim_data.iloc[:,0:int(n_group1)],sim_data.iloc[:,int(n_group1):], all_zips,i,p_cut) #delta cutoff
        calc_power_value = pd.DataFrame(calc_empirical_marg_power(confusion_matrix))
        calc_power_value['ID'] = np.repeat(i, 4).tolist()#i
        df_all_test = df_all_test.append(calc_power_value)
    df_all_test.columns = ['Power', 'Test', 'ID']
    #df_all_test.columns = ['Power', 'ID']
    print(df_all_test)
    return df_all_test


# fills a matrix with the parameters as the columns and the
def get_all_parameters(all_zips):
    environmental_tests_df = pd.DataFrame(index=range(0), columns=range(6))
    environmental_tests_df.columns = ['n_samples','n_CpGs','healthy_proportion','effect_size','n_modified_CpGs', 'ID']
    for i in range(0, all_zips):
        user_params = extract_data_from_zip(all_zips,i,"User_Parameters")
        user_params = user_params.T
        environmental_tests_df.loc[len(environmental_tests_df.index)] = user_params.iloc[1].tolist() + [i]
    return environmental_tests_df

def merge_data(all_zips, p_adjust_method_string, beta_or_M_string, p_cut):
    power_calc_val = power_calc_multiple_runs(all_zips,p_adjust_method_string,beta_or_M_string,p_cut)
    parameters = get_all_parameters(all_zips)
    output_df = power_calc_val.merge(parameters, how='inner', on='ID')
    print("The output for analysis: ")
    print(output_df)
    return output_df

def plot_line(data_df, y_parameter_string, varied_parameter_string, p_adjust_method_string, beta_or_M_string):
    #print("selected_df before")
    #print(data_df)
    selected_parameters_df = data_df.loc[:,['Power','Test',varied_parameter_string,y_parameter_string]]
    #selected_parameters_df = data_df.loc[:, ['Power', varied_parameter_string, y_parameter_string]]

    #print("selected_df after")

    y_parameter = y_parameter_string
    varied_parameter = varied_parameter_string
    for idx, val in enumerate(selected_parameters_df[varied_parameter].unique()):
        plt.figure()
        selected_parameters = selected_parameters_df.loc[selected_parameters_df[varied_parameter] == val,:]
        varied_param_i = val
        print("this is one unique")
        print(selected_parameters)
        plt.title("Sample Vs Power ("+varied_parameter+" = "+str(val)+")")
        plt.xlabel('Sample')
        plt.ylabel('Power')
        line_plot_fig = sns.lineplot(data=selected_parameters, x=y_parameter, y='Power', hue='Test')
        fig = line_plot_fig.get_figure()
        fig.savefig(str("LinePlot_") + str("vs_Power_") + varied_parameter + str("_step_") +
                    str(varied_param_i) + "_" + p_adjust_method_string + "_" + beta_or_M_string + ".png", dpi=100)
        #scatter_fig = sns.scatterplot(data=selected_parameters_df, x=y_parameter, y=selected_parameters_df['Power'], hue=selected_parameters_df['Test'])

def heat_map(data_df, x_parameter_string, y_parameter_string, p_adjust_method_string, test_string, beta_or_M_string):
    plt.figure()
    selected_parameters_df = data_df.loc[:, ['Power', 'Test', x_parameter_string, y_parameter_string]]
    #selected_parameters_df = data_df.loc[:, ['Power', x_parameter_string, y_parameter_string]]
    #print("in heat", selected_parameters_df)
    y_parameter = y_parameter_string
    x_parameter = x_parameter_string
    print(test_string)
    #print(selected_parameters_df.loc[selected_parameters_df['Test'] == test_string,:])
    #print(selected_parameters_df['Test'] == test_string)
    selected_parameters_df = selected_parameters_df.loc[selected_parameters_df['Test'] == test_string,:]

    print("just check")
    print(selected_parameters_df)
    #selected_parameters_df = selected_parameters_df.drop(columns=['healthy_proportion', 'n_CpGs', 'ID', 'Test', 'n_modified_CpGs'],axis=1)
    #selected_parameters_df = pd.unique(selected_parameters_df)
    #selected_parameters_df = selected_parameters_df.sort_values(by=['n_samples', 'Power'])
    df_m = selected_parameters_df.reset_index().pivot_table(index=y_parameter,columns=x_parameter, values='Power', aggfunc='mean')
    print("matrix", df_m)
    df_m = df_m.fillna(0)

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.title("HeatMap of "+y_parameter+ " and "+ x_parameter)
    plt.xlabel(x_parameter)
    plt.ylabel('Sample')
    ax = sns.heatmap(df_m, cmap='crest', annot=True)
    fig = ax.get_figure()
    fig.savefig(str("HeatMap_") + y_parameter + "_vs_" +x_parameter+ "_fill_Power_"+"_test_"+test_string+"_padjust_"+p_adjust_method_string+"_"+beta_or_M_string + ".png", dpi=100)

#Aggregate function that calculates marginal power for all the tests across samples based on user-defined inputs. As an output it generates line plots and heat maps that shows the relationship between power and two varying parameters and three constant parameters
def PowerCalc(num_zips, p_adjust_method_string, beta_or_M_string, y_parameter_string, x_parameter_string, color_string, test_vector_string,p_cut):
    all_test_all_zips = merge_data(num_zips,p_adjust_method_string,beta_or_M_string, p_cut)
    all_test_all_zips.to_csv(os.getcwd()+dirname+"PowerCalc_"+beta_or_M_string+"_"+p_adjust_method_string+'.csv', sep=",", index='ID',header=True)

    plot_line(all_test_all_zips, y_parameter_string, x_parameter_string,p_adjust_method_string,beta_or_M_string)

    for test in test_vector_string:
        heat_map(all_test_all_zips, x_parameter_string, y_parameter_string, p_adjust_method_string,test,beta_or_M_string)

num_zips = 64 # number of simulated data environmental workflows
p_adjust_method_string = "fdr" #"none"/"BH"/"bonferroni"/"fdr"
beta_or_M_string = "M" #"beta"/"M"
y_parameter_string = "n_samples" # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
x_parameter_string = "n_modified_CpGs" # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
color_string = "Power"
test_vector_string = ["T_test","Limma_test", "KS_test", "W_test"]
p_cut = 0.05 #0.01, 0.1, 0.05

PowerCalc(num_zips,p_adjust_method_string,beta_or_M_string,y_parameter_string,x_parameter_string,color_string,test_vector_string,p_cut)