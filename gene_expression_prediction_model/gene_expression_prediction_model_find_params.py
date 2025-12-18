#!/usr/bin/env python
# coding: utf-8

"""
Hyperparameter search for gene expression prediction models.
Grid search for Random Forest and Linear Regression models.
"""

import argparse
import json
import multiprocessing
import operator
import os

import numpy as np
import pandas as pd
import pyarrow.feather as feather
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.pipeline import Pipeline


def load_json(filepath):
    """Load JSON file."""
    with open(filepath, 'r') as f:
        return json.load(f)


def load_gene_categories(filepath):
    """Load gene selection categories from JSON file."""
    with open(filepath) as fp:
        gene_dict = json.load(fp)
    gene_dict["coding_all"] = gene_dict["coding_cpg"] + gene_dict["coding_non_cpg"]
    return gene_dict


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


def create_param_grid():
    """Create parameter grid for GridSearchCV."""
    rf_model = RandomForestRegressor(random_state=42)
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
    
    model_names = ["rf", "lr"]
    
    return param_grid, model_names


def run_grid_search(param_grid, model_names, X_train, y_train):
    """Run grid search for each model in parameter grid."""
    results = []
    
    for i, params in enumerate(param_grid):
        model_name = model_names[i]
        model = params['model'][0]
        
        print(f"\nTraining {model_name} model...", flush=True)
        
        # Remove model from params for GridSearchCV
        params_copy = params.copy()
        params_copy.pop('model')
        
        # Create pipeline
        steps = [('model', model)]
        pipeline = Pipeline(steps)
        
        # Grid search with cross-validation
        grid = GridSearchCV(
            pipeline,
            param_grid=params_copy,
            cv=5,
            scoring="neg_root_mean_squared_error",
            n_jobs=-1
        )
        grid.fit(X_train, y_train)
        
        # Store results
        results.append({
            'model_name': model_name,
            'best score': grid.best_score_,
            'best params': grid.best_params_,
        })
        
        print(f"  Best score: {grid.best_score_:.4f}")
        print(f"  Best params: {grid.best_params_}")
    
    # Sort by best score (descending)
    results = sorted(results, key=operator.itemgetter('best score'), reverse=True)
    
    return results


def main():
    """Main execution function."""
    print(f"Number of CPUs: {multiprocessing.cpu_count()}")
    print("Starting parameter search...\n")
    
    # Configuration
    model_design = "rnaseq_vs_hiplex_rm_outlier_log"
    gene_select_name = "coding_all"
    
    # Paths
    gene_select_dir = "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island"
    rnaseq_dir = "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq"
    wgc_root_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc"
    output_dir = f'/dcs05/hongkai/data/next_cutntag/script/explainability/{model_design}'
    
    os.makedirs(output_dir, exist_ok=True)
    
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
    
    # Split data
    rnaseq_wgc_train, rnaseq_wgc_test = train_test_split(
        rnaseq_wgc, test_size=0.2, random_state=42
    )
    
    X_train = rnaseq_wgc_train.drop(columns=["sqrt_V"], inplace=False)
    y_train = rnaseq_wgc_train.loc[:, "sqrt_V"].values
    
    print(f"Training set shape: {X_train.shape}")
    print(f"Test set shape: {rnaseq_wgc_test.shape}\n")
    
    # Create parameter grid
    param_grid, model_names = create_param_grid()
    
    # Run grid search
    results = run_grid_search(param_grid, model_names, X_train, y_train)
    
    # Save results
    output_file = os.path.join(output_dir, f"{gene_select_name}.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=4)
    
    print(f"\n{'='*60}")
    print("Parameter search complete!")
    print(f"Results saved to: {output_file}")
    print(f"{'='*60}")
    print("\nBest models (sorted by score):")
    for result in results:
        print(f"\n{result['model_name'].upper()}:")
        print(f"  Score: {result['best score']:.4f}")
        print(f"  Params: {result['best params']}")


if __name__ == "__main__":
    main()