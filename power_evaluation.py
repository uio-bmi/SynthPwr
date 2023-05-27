import os
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
from scipy.stats import ttest_ind

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
    list_of_delta = [row for row in range(0, len(group1_means_vector.index))]
    list_of_delta_scores = []
    for i in list_of_delta:
        list_of_delta_scores.append(np.abs((np.nanmean(group1_means_vector.iloc[i,:]) - np.nanmean(group2_means_vector.iloc[i,:]))))
    return list_of_delta

def unpack_workflow(workflow_num):
    list_of_simulated_data_df = np.load(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"Simulated_data"+"_workflow_"+str(workflow_num)+".npz")
    list_of_true_CpG_locations = np.load(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"truly_different_sites_indices"+"_workflow_"+str(workflow_num)+".npz")
    user_params_df = pd.read_csv(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"User_Parameters"+"_workflow_"+str(workflow_num)+".csv", index_col=0)
    result = [list_of_simulated_data_df, list_of_true_CpG_locations, user_params_df]
    return result

def unpack_workflow_params_only(workflow_num):
    user_params_df = pd.read_csv(os.getcwd()+dirname+"workflow_"+str(workflow_num)+os.sep+"User_Parameters"+"_workflow_"+str(workflow_num)+".csv", index_col=0)
    return user_params_df

#targets simulated data characteristics to determine num samples
def get_num_sample_per_group(user_params, simulated_df):
    healthy_proportion = user_params.iloc[2,1]
    n_Group1 = len(simulated_df.columns)*healthy_proportion
    n_Group2 = len(simulated_df.columns)*(1-healthy_proportion)
    result = [n_Group1,n_Group2]
    return result

def tradtional_t_test(simulated_df, n_group1):
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
    n_group1 = int(n_group1)
    all_p_val = [num for num in range(0, len(simulated_df.index))]
    for i in range(0, len(simulated_df.index)):
        gr1 = simulated_df.iloc[i, 0:n_group1]
        gr2 = simulated_df.iloc[i, n_group1:]
        w_test_statistic, w_test_pvalue = scipy.stats.ranksums(gr1, gr2)
        all_p_val[i] = w_test_pvalue
    return all_p_val

def limma_test(simulated_df, n_group1, n_group2):
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
    logit_result = ENmix.B2M(beta_matrix)
    return logit_result

#Running all the tests to compare the two groups for every CpG site based on either transformed or original values. The output is p-values for every CpG site and for every test
def run_all_test(n_group1, n_group2, sim_data):
    t_test_output = tradtional_t_test(sim_data, n_group1)
    limma_test_output = limma_test(sim_data, n_group1, n_group2)
    w_test_output = w_test(sim_data, n_group1)

    test_dict = {"T-test": t_test_output,"limma-test": limma_test_output,"W-test": w_test_output}
    all_test_output = pd.DataFrame(test_dict)
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
    store_p_adjust.columns = ["T-test", "limma", "W-test"]
    return store_p_adjust

# prior to confusion matrix construction, CpGs should be labeled based on two distinction characteristics: whether they are truly modified, and if their adjusted p-value is less then 0.05
# The input is the vector with observed values within the simulated data. The output is a boolean which represents if the observation is signficantly different between the two groups (0) and not significant (1)
def create_predicted_vector(p_val_vector, group1_means_vector,group2_means_vector, p_cut):
    delta_beta_cpg = calculate_delta_beta(group1_means_vector,group2_means_vector)
    boolean_predicted_vector = [False for row in range(0, len(p_val_vector.index))]
    for i in range(0, len(delta_beta_cpg)):
        if delta_beta_cpg[i] >= 0.01 and p_val_vector.iloc[i] < p_cut:
            boolean_predicted_vector[i] = True
        else:
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
    df_all_test = pd.DataFrame(index=range(3), columns=range(2))
    df_all_test.columns = ['Power', 'Test']
    df_all_test['Test'] = ['T_test','Limma_test', 'W_test']

    for i in range(0, len(workflows_confusion_matrix)):
        confusion_matrix = workflows_confusion_matrix[i]
        true_positive = confusion_matrix[1][1][1]
        false_positive = confusion_matrix[1][1][0]
        false_negative = confusion_matrix[1][0][1]
        true_negative = confusion_matrix[1][0][0]
        power_calc_val = int(true_positive)/(int(true_positive)+int(false_negative))
        df_all_test.iloc[i, 0] = power_calc_val
    return df_all_test

def power_calc_multiple_runs(workflow):
    print("Power Calculation for Workflow: ", workflow)
    env_params_df = pd.read_csv(os.getcwd() + summarycalcdirname + "env_inputparams.csv", index_col=0)
    num_simulation = int(env_params_df.iloc[0, 1])
    p_cut = 0.05
    df_workflow_all_test_simiter_list = []
    workflow_files = unpack_workflow(workflow)
    list_of_sim_data = workflow_files[0]
    list_of_truly_different = workflow_files[1]
    user_params = workflow_files[2]

    list_of_sim_data = list_of_sim_data['arr_0']
    list_of_truly_different = list_of_truly_different['arr_0']
    for sim_iter in range(0, num_simulation):
        simulated_data = pd.DataFrame(list_of_sim_data[sim_iter])
        truly_different = list_of_truly_different[sim_iter]
        n_samples_per_gr = get_num_sample_per_group(user_params, simulated_data)
        n_group1 = n_samples_per_gr[0]
        n_group2 = n_samples_per_gr[1]

        all_test = run_all_test(n_group1, n_group2, simulated_data)
        p_adjust_all_test = multitest_p_adjust(all_test)
        confusion_matrix = create_confusion_matrix(p_adjust_all_test,simulated_data.iloc[:,0:int(n_group1)],simulated_data.iloc[:,int(n_group1):], truly_different, p_cut)
        calc_power_value_workflow = pd.DataFrame(calc_empirical_marg_power(confusion_matrix))
        calc_power_value_workflow['ID'] = np.repeat(workflow, 3).tolist()
        calc_power_value_workflow.columns = ['Power', 'Test', 'ID']
        df_workflow_all_test_simiter_list.append(calc_power_value_workflow)
    return df_workflow_all_test_simiter_list

def merge_data(num_workflows, num_simulations, num_processes):
    pool = Pool(processes=num_processes)
    list_of_workflows = [num for num in range(num_workflows)]
    environmental_tests_df = pd.DataFrame(index=range(num_workflows), columns=range(6))
    environmental_tests_df.columns = ['n_samples', 'n_CpGs', 'healthy_proportion', 'effect_size', 'n_modified_CpGs','ID']
    power_calc_result = pool.map(power_calc_multiple_runs, list_of_workflows)
    parameters = pool.map(unpack_workflow_params_only, list_of_workflows)
    pool.close()
    output_of_map_list = list(range(num_simulations))
    for sim_run in range(0, num_simulations):
        output_of_map_list[sim_run] = pd.DataFrame()
    for workflow_run in range(0, len(power_calc_result)):
        workflow = power_calc_result[workflow_run]
        for sim_run in range(0, len(workflow)):
            output_of_map_list[sim_run] = output_of_map_list[sim_run].append(workflow[sim_run], ignore_index=True)
    i = 0
    for param_combo in parameters:
        param_combo = param_combo.T
        environmental_tests_df.iloc[i] = param_combo.iloc[1].tolist() + [i]
        i += 1
    simiter_output = []
    for sim_run in range(0, num_simulations):
        output_df = output_of_map_list[sim_run].merge(environmental_tests_df, how='inner', on='ID')
        simiter_output.append(output_df)
    print("The output for analysis: ")
    print(simiter_output)
    return simiter_output

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

#Aggregate function that calculates marginal power for all the tests across samples based on user-defined inputs. As an output it generates line plots and heat maps that shows the relationship between power and two varying parameters and three constant parameters
def PowerCalc(num_workflows, num_simulations,num_processes,p_adjust_method_string, y_parameter_string, x_parameter_string, test_vector_string):
    all_test_all_zips = merge_data(num_workflows, num_simulations,num_processes)
    for sim_iter in range(0, num_simulations):
        all_test_all_zips[sim_iter].to_csv(os.getcwd()+summarycalcdirname+"PowerCalc_iter"+str(sim_iter)+'.csv', sep=",", index='ID',header=True)

    plot_line(all_test_all_zips, num_workflows, num_simulations, len(test_vector_string),y_parameter_string, x_parameter_string,p_adjust_method_string)
    plot_every_effect_line(all_test_all_zips, num_workflows, num_simulations, len(test_vector_string),"n_samples", "effect_size", p_adjust_method_string)
    for test in test_vector_string:
        heat_map(all_test_all_zips, num_workflows, num_simulations, len(test_vector_string), x_parameter_string, y_parameter_string, p_adjust_method_string,test)

if __name__ == '__main__':
    env_params_df = pd.read_csv(os.getcwd()+summarycalcdirname+"env_inputparams.csv", index_col=0)
    num_simulations = int(env_params_df.iloc[0, 1])
    num_workflows = int(env_params_df.iloc[1, 1])
    num_processes = int(env_params_df.iloc[2, 1])
    p_adjust_method_string = "fdr" #"none"/"BH"/"bonferroni"/"fdr"
    y_parameter_string = "n_samples" # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
    x_parameter_string = env_params_df.iloc[3, 1] # "n_samples"/"n_CpGs"/"healthy_proportion"/"effect_size"/"n_modified_CpGs"
    test_vector_string = ["T_test","Limma_test","W_test"]
    p_cut = 0.05 #0.01, 0.1, 0.05
    #age_group = env_params_df.iloc[3, 1]

    PowerCalc(num_workflows,num_simulations,num_processes,p_adjust_method_string,y_parameter_string,x_parameter_string,test_vector_string)
    endtime = time.time()
    print("Time taken: ", endtime - starttime)
