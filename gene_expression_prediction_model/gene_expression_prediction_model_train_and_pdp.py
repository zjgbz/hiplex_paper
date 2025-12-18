#!/usr/bin/env python
# coding: utf-8

"""
Model training and partial dependence plot generation for gene expression prediction.
"""

import argparse
import json
import multiprocessing
import os
import pickle

import numpy as np
import pandas as pd
import pyarrow.feather as feather
import scipy.stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import partial_dependence
from sklearn.linear_model import LinearRegression
from sklearn.metrics import root_mean_squared_error
from sklearn.model_selection import train_test_split


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Model train and pdp')
    parser.add_argument('--random_seed', type=int, required=True, help='Random seed')
    return parser.parse_args()


def load_gene_categories(filepath):
    """Load gene selection categories from JSON file."""
    with open(filepath) as fp:
        gene_dict = json.load(fp)
    gene_dict["coding_all"] = gene_dict["coding_cpg"] + gene_dict["coding_non_cpg"]
    return gene_dict


def load_json(filepath):
    """Load JSON file."""
    with open(filepath, 'r') as f:
        return json.load(f)


def column_to_rownames(df, var="pos"):
    """Set a column as the dataframe index."""
    return df.set_index(var)


def remove_model_prefix(params_dict):
    """Remove 'model__' prefix from parameter keys."""
    return {key.replace("model__", ""): value for key, value in params_dict.items()}


def load_rnaseq_data(filepath):
    """Load and preprocess RNA-seq data."""
    rnaseq = pd.read_csv(filepath, sep=",", header=0, index_col=0)
    rnaseq.loc[:, "sqrt_V"] = np.log10(rnaseq.loc[:, "V1V2"] + 1)
    return rnaseq


def load_wgc_data(wgc_dir_filename, features):
    """Load and preprocess WGC data."""
    wgc_raw = feather.read_feather(wgc_dir_filename)
    wgc_vals = wgc_raw[features]
    zero_genes = wgc_raw["pos"][wgc_vals.sum(axis=1) == 0].to_list()
    
    # Log transform and normalize
    wgc_raw[features] = np.log10(wgc_raw[features] + 1)
    wgc_raw[features] = (wgc_raw[features] - wgc_raw[features].min()) / \
                        (wgc_raw[features].max() - wgc_raw[features].min())
    
    wgc_raw = wgc_raw.loc[:, ["pos"] + features]
    return wgc_raw, zero_genes


def prepare_dataset(rnaseq, wgc_raw, overlap_genes, q99_cutoff, zero_all_genes):
    """Prepare merged RNA-seq and WGC dataset."""
    rnaseq_avail = rnaseq.loc[overlap_genes, ["gene_id", "sqrt_V"]]
    rnaseq_avail = rnaseq_avail.loc[rnaseq_avail.loc[:, "sqrt_V"] <= q99_cutoff, :]
    
    rnaseq_wgc = rnaseq_avail.merge(wgc_raw, left_on="gene_id", right_on="pos", how="inner")
    rnaseq_wgc = rnaseq_wgc.set_index("gene_id")
    rnaseq_wgc = rnaseq_wgc.loc[~rnaseq_wgc.index.isin(zero_all_genes), :]
    rnaseq_wgc = rnaseq_wgc.drop(columns=["pos"], inplace=False)
    
    return rnaseq_wgc


def split_data(rnaseq_wgc):
    """Split data into features and target."""
    X = rnaseq_wgc.drop(columns=["sqrt_V"], inplace=False)
    y = rnaseq_wgc.loc[:, "sqrt_V"].values
    return X, y


def load_models(model_params_file, gene_select_name, random_seed):
    """Load and configure models based on saved parameters."""
    with open(model_params_file) as fp:
        model_params = json.load(fp)
    
    models = {}
    for param_set in model_params:
        model_name = param_set['model_name']
        
        if model_name == "rf":
            model = RandomForestRegressor(random_state=random_seed)
        else:
            model = LinearRegression()
        
        params = param_set['best params']
        if params:
            params = remove_model_prefix(params)
            model.set_params(**params)
        
        models[model_name] = model
    
    return models


def train_and_evaluate(models, X, y):
    """Train models and evaluate performance."""
    results = []
    
    for model_name, model in models.items():
        model.fit(X, y)
        y_pred = model.predict(X)
        
        rmse = root_mean_squared_error(y, y_pred)
        r_squared = scipy.stats.pearsonr(y, y_pred)[0] ** 2
        
        print(f"\n{model_name}:")
        print(f"  RMSE: {rmse:.4f}")
        print(f"  R²: {r_squared:.4f}")
        
        results.append({
            "model_name": model_name,
            "RMSE": rmse,
            "Pearson^2": r_squared
        })
    
    return pd.DataFrame(results)


def extract_feature_importance(model, feature_names):
    """Extract and sort feature importances from Random Forest model."""
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': model.feature_importances_
    }).sort_values(by='importance', ascending=False)
    
    return feature_importance


def compute_partial_dependence(model, X, target_pairs):
    """Compute partial dependence plots for top features."""
    pdp_results = {}
    
    for target_pair in target_pairs:
        X_clean = X.loc[:, ~X.columns.isin(["gene_id", "sqrt_V", "pos", "cluster"])]
        X_clean = X_clean.reset_index(drop=True)
        
        pdp_lines = partial_dependence(
            model, X_clean,
            features=[target_pair],
            grid_resolution=85,
            kind="both",
            percentiles=(0, 1)
        )
        pdp_results[target_pair] = pdp_lines
    
    return pdp_results


def main():
    """Main execution function."""
    args = parse_arguments()
    random_seed = args.random_seed
    
    print(f"Number of CPUs: {multiprocessing.cpu_count()}")
    print("Starting analysis...\n")
    
    # Configuration
    model_design = "rnaseq_vs_hiplex_rm_outlier_log"
    gene_select_name = "coding_all"
    
    # Paths
    gene_select_dir = "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island"
    rnaseq_dir = "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq"
    wgc_root_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc"
    save_dir = f"/dcs05/hongkai/data/next_cutntag/bulk/explainability/{model_design}/{random_seed}"
    
    os.makedirs(save_dir, exist_ok=True)
    
    # Load gene categories
    gene_dict = load_gene_categories(
        os.path.join(gene_select_dir, "gene_categories.json")
    )
    gene_select_list = gene_dict[gene_select_name]
    print(f"Gene category: {gene_select_name}")
    
    # Load cutoffs and features
    rnaseq_cutoffs = load_json(
        "/dcs05/hongkai/data/next_cutntag/bulk/explainability/rnaseq_hiplex_cutoff.json"
    )
    features = load_json(
        '/dcs05/hongkai/data/next_cutntag/script/utils/filtered_target_pairs.json'
    )
    
    # Load RNA-seq data
    rnaseq = load_rnaseq_data(os.path.join(rnaseq_dir, "RNA_seq_TPM_all.csv"))
    
    # Load WGC data
    frag_type = "mixed"
    bin_size = "promoter_-1000-1000"
    scen = "V"
    target_qc_type = "all"
    post_process = "libnorm"
    
    wgc_filename = f"{scen}_{frag_type}_{bin_size}_colQC-{target_qc_type}_{post_process}.feather"
    wgc_dir_filename = os.path.join(wgc_root_dir, frag_type, bin_size, wgc_filename)
    wgc_raw, zero_hiplex_genes = load_wgc_data(wgc_dir_filename, features)
    
    # Prepare dataset
    overlap_genes = list(set(wgc_raw.loc[:, "pos"].tolist()) & set(gene_select_list))
    zero_rnaseq_genes = rnaseq.loc[overlap_genes, "gene_id"][
        rnaseq.loc[overlap_genes, "sqrt_V"] == 0
    ].tolist()
    zero_all_genes = list(set(zero_hiplex_genes) & set(zero_rnaseq_genes))
    
    q99_cutoff = rnaseq_cutoffs[model_design]
    rnaseq_wgc = prepare_dataset(rnaseq, wgc_raw, overlap_genes, q99_cutoff, zero_all_genes)
    X_all, y_all = split_data(rnaseq_wgc)
    
    print(f"Dataset shape: {X_all.shape}\n")
    
    # Load models
    model_params_file = f'/dcs05/hongkai/data/next_cutntag/script/explainability/{model_design}/{gene_select_name}.json'
    models = load_models(model_params_file, gene_select_name, random_seed)
    
    # Train and evaluate models
    results_df = train_and_evaluate(models, X_all, y_all)
    results_df.to_csv(f"{save_dir}/model_performance_results.csv", index=False)
    print(f"\nResults saved to {save_dir}/model_performance_results.csv")
    
    # Extract feature importance from Random Forest
    rf_model = models["rf"]
    feature_importance = extract_feature_importance(rf_model, X_all.columns)
    feature_importance.to_csv(f"{save_dir}/top_features.csv", index=False)
    
    top_features_dict = {
        gene_select_name: {
            model_design: feature_importance["feature"].tolist()
        }
    }
    with open(f"{save_dir}/top_feature_importance.json", "w") as f:
        json.dump(top_features_dict, f, indent=4)
    
    print("\nTop 10 features:")
    print(feature_importance.head(10))
    
    # Save data
    X_all.to_csv(f"{save_dir}/rnaseq_wgc_all_X.csv")
    
    # Compute partial dependence plots
    print("\nComputing partial dependence plots...")
    target_pairs = top_features_dict[gene_select_name][model_design]
    pdp_results = {
        gene_select_name: compute_partial_dependence(rf_model, X_all, target_pairs)
    }
    
    with open(f"{save_dir}/top_features_pdp.pkl", 'wb') as f:
        pickle.dump(pdp_results, f)
    
    print(f"\nPDP results saved to {save_dir}/top_features_pdp.pkl")
    print("Analysis complete!")


if __name__ == "__main__":
    main()