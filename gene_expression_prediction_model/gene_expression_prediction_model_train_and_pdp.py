#!/usr/bin/env python
# coding: utf-8

# In[29]:


import json
import multiprocessing
import operator
import os
import pickle
from itertools import product

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotx
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import scipy
import seaborn as sns
import shap
import xgboost as xgb
from adjustText import adjust_text
from sklearn.cluster import KMeans
from sklearn.datasets import load_svmlight_file
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import partial_dependence
from sklearn.linear_model import LinearRegression
from sklearn.metrics import root_mean_squared_error, silhouette_score
from sklearn.model_selection import (
    GridSearchCV,
    KFold,
    cross_val_score,
    cross_validate,
    train_test_split,
)
from sklearn.pipeline import Pipeline
import shap
import json
import operator
import scipy
import joblib
import multiprocessing
import argparse

def map_target_names(target_pair_list, target_pair_mapping_df, from_col="targets", to_col="shorthand" ):
    cur_names = target_pair_mapping_df.loc[:,from_col].to_list()
    new_names = target_pair_mapping_df.loc[:,to_col].to_list()
    print()
    result = target_pair_list
    for i in range(len(cur_names)):
        cur_name = cur_names[i]
        new_name = new_names[i]
        result = [target_pair.replace(cur_name, new_name) for target_pair in result]
    return(result)
def column_to_rownames(wgc, var="pos"):
    pos_list_full = wgc[var]
    wgc = wgc.set_index(var)
    return wgc
def remove_model_in_param_grid(best_params):
    results = {}
    for param in best_params.keys():
        new_param = param.replace("model__", "")
        results[new_param] = best_params[param]
    return results
# Instantiate the parser
parser = argparse.ArgumentParser(description='Model train and pdp')
# Required positional argument
parser.add_argument('--gene_select_name', type=str,
                    help='coding_cpg or coding_non_cpg or coding_all')

# Optional positional argument
parser.add_argument('--model_design', type=str, 
                    help='rnaseq_vs_hiplex or df_vs_df')
parser.add_argument('--random_seed', type=int, 
                    help='seed')

args = parser.parse_args()

print("num cpu: " + str(multiprocessing.cpu_count()), flush=True)

# In[2]:

print("start", flush=True)
gene_select_dir = "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island"
gene_select_filename = "gene_categories.json"
gene_select_dir_filename = os.path.join(gene_select_dir, gene_select_filename)
with open(gene_select_dir_filename) as fp:
    gene_select_dict = json.load(fp)
gene_select_dict["coding_all"] = gene_select_dict["coding_cpg"] + gene_select_dict["coding_non_cpg"]
# gene_select_name_list = ["coding_cpg"]
gene_select_name = args.gene_select_name
print(gene_select_name, flush=True)
gene_select_list = gene_select_dict[gene_select_name]

def safe_minmax(df):
    df_scaled = df.copy()
    for col in df.columns:
        col_min = df[col].min()
        col_max = df[col].max()
        if col_max != col_min:
            df_scaled[col] = (df[col] - col_min) / (col_max - col_min)
        else:
            df_scaled[col] = 0.0  # or leave as is
    return df_scaled

def column_to_rownames(wgc, var="pos"):
    pos_list_full = wgc[var]
    wgc = wgc.set_index(var)
    return wgc

with open("/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_hiplex_cutoff.json", 'r') as f:
    rnaseq_cutoffs = json.load(f)

# In[5]:
model_design = args.model_design

random_seed = args.random_seed
print(model_design, flush=True)

if model_design == "rnaseq_vs_hiplex_rm_outlier_log":
    rnaseq_dir = "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq"
    rnaseq_filename = "RNA_seq_TPM_all.csv"
    rnaseq_dir_filename = os.path.join(rnaseq_dir, rnaseq_filename)
    rnaseq = pd.read_csv(rnaseq_dir_filename, sep=",", header=0, index_col=0)
    # rnaseq.loc[:, "sqrt_V"] = rnaseq.loc[:, "sqrt_V1V2"]
    # rnaseq.loc[:, "sqrt_V"] = rnaseq.loc[:, "sqrt_V"]**2
    rnaseq.loc[:, "sqrt_V"] = np.log10(rnaseq.loc[:, "V1V2"]+1)
    
    wgc_root_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc"
    frag_type = "mixed"
    bin_size = "promoter_-1000-1000"
    scen = "V"
    target_qc_type = "all"
    post_process = "libnorm"
    wgc_dir = os.path.join(wgc_root_dir, frag_type, bin_size)
    # /dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/promoter_-1000-1000/V_mixed_promoter_-1000-1000_colQC-all_libnorm.featherX
    wgc_filename = f"{scen}_{frag_type}_{bin_size}_colQC-{target_qc_type}_{post_process}.feather"
    wgc_dir_filename = os.path.join(wgc_dir, wgc_filename)
    wgc_raw = feather.read_feather(wgc_dir_filename)
    with open('/dcs05/hongkai/data/next_cutntag/script/utils/filtered_target_pairs.json', 'r') as file:
        features = json.load(file)
    wgc_vals = wgc_raw[features]
    zero_hiplex_genes = wgc_raw["pos"][wgc_vals.sum(axis=1) == 0].to_list()
    wgc_raw[features] = np.log10(wgc_raw[features]+1)
    wgc_raw[features] = (wgc_raw[features] - wgc_raw[features].min()) / (wgc_raw[features].max() - wgc_raw[features].min())
    
    features = ["pos"]+features
    wgc_raw = wgc_raw.loc[:, features]
    
    wgc_raw_list = wgc_raw.loc[:, "pos"].values.tolist()
    overlap_gene_list = list(set(wgc_raw_list) & set(gene_select_list))
    rnaseq_avail = rnaseq.loc[overlap_gene_list, ["gene_id", "sqrt_V"]]
    # filter outlier
    q99 = rnaseq_cutoffs[model_design]
    rnaseq_avail = rnaseq_avail.loc[rnaseq_avail.loc[:, "sqrt_V"] <= q99,:]
    zero_rnaseq_genes = rnaseq_avail["gene_id"][rnaseq_avail["sqrt_V"]==0].to_list()
    zero_all_genes = list(set(zero_hiplex_genes) & set(zero_rnaseq_genes))
    # rnaseq_avail = rnaseq_avail.loc[rnaseq_avail.loc[:, "sqrt_V"] > 0,:]
    rnaseq_wgc_raw = rnaseq_avail.merge(wgc_raw, left_on="gene_id", right_on="pos", how="inner")
    rnaseq_wgc_raw = rnaseq_wgc_raw.set_index("gene_id")
    # print(rnaseq_wgc_raw.shape)
    rnaseq_wgc_raw = rnaseq_wgc_raw.loc[~rnaseq_wgc_raw.index.isin(zero_all_genes), :]
    # print(rnaseq_wgc_raw.shape)
    rnaseq_wgc = rnaseq_wgc_raw.drop(columns=["pos"], inplace=False)
    rnaseq_wgc_train, rnaseq_wgc_test = train_test_split(rnaseq_wgc, test_size=0.2, random_state=42)
    rnaseq_wgc_train_X = rnaseq_wgc_train.drop(columns=["sqrt_V"], inplace=False)
    rnaseq_wgc_train_y = rnaseq_wgc_train.loc[:, "sqrt_V"].values
    rnaseq_wgc_test_X = rnaseq_wgc_test.drop(columns=["sqrt_V"], inplace=False)
    rnaseq_wgc_test_y = rnaseq_wgc_test.loc[:, "sqrt_V"].values
    rnaseq_wgc_all_X = rnaseq_wgc.drop(columns=["sqrt_V"], inplace=False)
    rnaseq_wgc_all_y = rnaseq_wgc.loc[:, "sqrt_V"].values

all_models = {}

all_models["coding_all"] = {}
model_params_file = f'/dcs05/hongkai/data/next_cutntag/script/explainability/{model_design}/{gene_select_name}.json'
with open(model_params_file) as fp:
    model_params = json.load(fp)
for i in range(len(model_params)):
    model_name = model_params[i]['model_name']
    if model_name == "rf":
        model = RandomForestRegressor(random_state=random_seed)
    else:
        model = LinearRegression()
    params = model_params[i]['best params']

    if len(params) > 0:
        params = remove_model_in_param_grid(params)
        print(params)
        model.set_params(**params)
    all_models["coding_all"][model_name] = model

print(rnaseq_wgc_all_X.shape)
save_dir = f"/dcs05/hongkai/data/next_cutntag/bulk/explainability/{model_design}/{random_seed}"
os.makedirs(save_dir, exist_ok=True)
results = []
for frag_type in ['mixed']:
    print(frag_type)
    # rnaseq_wgc_all_X = data[gene_select_name][model_design][frag_type]["rnaseq_wgc_all_X"]
    # rnaseq_wgc_all_y = data[gene_select_name][model_design][frag_type]["rnaseq_wgc_all_y"]
    models = {"rf": all_models[gene_select_name]["rf"], 
              "lr": all_models[gene_select_name]["lr"]}
    for model_name, model in models.items():
        model.fit(rnaseq_wgc_all_X, rnaseq_wgc_all_y)
        rnaseq_wgc_pred_y = model.predict(rnaseq_wgc_all_X)
        rmse = root_mean_squared_error(rnaseq_wgc_all_y, rnaseq_wgc_pred_y)
        cor = scipy.stats.pearsonr(rnaseq_wgc_all_y, rnaseq_wgc_pred_y)[0]**2
        print(model_name)
        print(f'Root Mean Squared Error: {rmse}')
        print(f'Pearson^2 for all data: {cor}')
        all_models[gene_select_name][model_name] = model
        results.append({
            "model_name": model_name,
            "RMSE": rmse,
            "Pearson^2": cor
        })
results_df = pd.DataFrame(results)
results_df.to_csv(f"{save_dir}/model_performance_results.csv", index=False)


    
top_features_dict ={}
# gene_select_name ="coding_all"
top_features_dict[gene_select_name] = {}

print(gene_select_name)
model = all_models[gene_select_name]["rf"]
# rnaseq_wgc_all_X = data[gene_select_name]['rnaseq_vs_hiplex']["mixed"]["rnaseq_wgc_all_X"]
feature_importance = pd.DataFrame({
    'feature': rnaseq_wgc_all_X.columns,
    'importance': model.feature_importances_
}).sort_values(by='importance', ascending=False)

feature_importance.to_csv(f"{save_dir}/top_features.csv", index=False)

top_features_dict[gene_select_name][model_design] = feature_importance["feature"].to_list()
print(feature_importance)
with open(f"{save_dir}/top_feature_importance.json", "w") as f:
    json.dump(top_features_dict, f, indent=4)
# gene_select_names = ["coding_cpg", "coding_non_cpg"]
pdp_lines_result = {}

# gene_select_name ="coding_all"
pdp_lines_result[gene_select_name] = {}

# we dont need to apply cluster id anymore at this step.
# cluster_result = pd.read_csv(f"/dcs05/hongkai/data/next_cutntag/bulk/explainability/leave_one_out/V/all/all_clusters=5.csv", index_col=0)
# cluster_result = cluster_result.set_index("gene_id")
# rnaseq_wgc_all_X = data[gene_select_name][model_design]["mixed"]["rnaseq_wgc_all_X"]
# rnaseq_wgc_all_X["cluster"] = cluster_result.loc[rnaseq_wgc_all_X.index, "cluster"]

rnaseq_wgc_all_X.to_csv(f"{save_dir}/rnaseq_wgc_all_X.csv")
model = all_models[gene_select_name]["rf"]
pdp_lines_dict = {}
target_pairs = top_features_dict[gene_select_name][model_design]
remove_TV = False
if remove_TV:
    target_pairs = [s for s in target_pairs if not s.startswith(("T_", "V_"))]

print(target_pairs)
for target_pair in target_pairs:
    cluster_result_i = rnaseq_wgc_all_X
    rnaseq_wgc_all_X = cluster_result_i.loc[:,~cluster_result_i.columns.isin(["gene_id", "sqrt_V", "pos", "cluster"])]
    rnaseq_wgc_all_X = rnaseq_wgc_all_X.reset_index(drop=True)
    pdp_lines = partial_dependence(model, rnaseq_wgc_all_X, 
        features = [target_pair], grid_resolution = 85, kind = "both", percentiles = (0,1)
    )
    pdp_lines_result[gene_select_name][target_pair] = pdp_lines
with open(f"{save_dir}/top_features_pdp.pkl", 'wb') as f:
    pickle.dump(pdp_lines_result, f)

# %%
