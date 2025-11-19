#!/usr/bin/env python
# coding: utf-8

# In[29]:


import os
import pyarrow.feather as feather
import pandas as pd
import numpy as np
from sklearn.datasets import load_svmlight_file
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score, cross_validate, GridSearchCV
from sklearn.inspection import partial_dependence
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import xgboost as xgb
import shap
import json
import operator
import scipy
import joblib
import multiprocessing
import argparse

# Instantiate the parser
parser = argparse.ArgumentParser(description='Model paramters search')
# Required positional argument
parser.add_argument('--gene_select_name', type=str,
                    help='coding_cpg or coding_non_cpg, coding_all')

# Optional positional argument
parser.add_argument('--model_design', type=str, 
                    help='rnaseq_vs_hiplex_rm_outlier_log or more')

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


print("start training", flush=True)
rf_model = RandomForestRegressor(random_state = 42)
# xgb_model = xgb.XGBRegressor(objective='reg:squarederror', random_state=42)
lr_model = LinearRegression()

param_grid = [
    {
        'model': [rf_model],
        'model__bootstrap': [True, False],
        'model__max_depth': [10, 20, 30, None],
        'model__max_features': ['log2', 'sqrt', 1.0],
        'model__min_samples_leaf': [1, 2, 4],
        'model__min_samples_split': [2, 5, 10],
        'model__n_estimators': [50, 100, 150]
    },
    {
        'model': [lr_model],
    }
]
result = []
model_names = ["rf", "lr"]
for i in range(len(param_grid)):
    params = param_grid[i]
    model_name = model_names[i]
    #classifier
    model = params['model'][0]

    print("start training model x", flush=True)
    #getting arguments by
    #popping out classifier
    params.pop('model')

    #pipeline
    steps = [('model', model)]

    #cross validation using
    #Grid Search
    grid = GridSearchCV(Pipeline(steps), param_grid=params, cv=5, scoring="neg_root_mean_squared_error", n_jobs=-1)
    grid.fit(rnaseq_wgc_train_X, rnaseq_wgc_train_y)

    #storing result
    result.append\
    (
        {
            'model_name': model_name,
            'best score': grid.best_score_,
            'best params': grid.best_params_,
        }
    )
result = sorted(result, key=operator.itemgetter('best score'),reverse=True)
# grid = result[0]['grid']


# joblib.dump(grid, "/dcs05/hongkai/data/next_cutntag/bulk/explainability/non_linear_model/best_model.pkl") 
print("model saved")
print(result)
# In[39]:

os.makedirs(f'/dcs05/hongkai/data/next_cutntag/script/explainability/{model_design}/', exist_ok=True)
with open(f'/dcs05/hongkai/data/next_cutntag/script/explainability/{model_design}/{gene_select_name}.json', 'w') as f:
    json.dump(result, f)



