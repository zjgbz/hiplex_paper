import pyarrow.feather as feather
import matplotlib.pyplot as plt
from itertools import combinations
import pandas as pd
import numpy as np
from random import shuffle
import random
import copy
import sys
import os

def coexpr_generation(df_1, df_2, cooccur_scen, coexp_col=["0|0", "1|0", "0|1", "1|x", "x|1", "1|1"]):
    coexp_df = pd.DataFrame(columns=coexp_col, index=list(df_1.index))
    coexp_df.loc[:, "1|x"] = (df_1 != 0).sum(axis=1)
    coexp_df.loc[:, "x|1"] = (df_2 != 0).sum(axis=1)
    df_1and2 = df_1 + df_2
    df_1x2 = df_1.multiply(df_2)
    coexp_df.loc[:, "0|0"] = (df_1and2 == 0).sum(axis=1)
    coexp_df.loc[:, "1|1"] = (df_1x2 != 0).sum(axis=1)
    coexp_df.loc[:, "1|0"] = coexp_df.loc[:, "1|x"] - coexp_df.loc[:, "1|1"]
    coexp_df.loc[:, "0|1"] = coexp_df.loc[:, "x|1"] - coexp_df.loc[:, "1|1"]
    coexp_array = coexp_df.loc[:, "1|1"].values
    total_num = coexp_df.loc[:, "1|0"].sum() + coexp_df.loc[:, "0|1"].sum() + np.sum(coexp_array)
    coexp_array_binary = copy.deepcopy(coexp_array)
    coexp_array_binary[coexp_array != 0] = 1
    coexp_cell_num = np.sum(coexp_array)
    coexp_region_num = np.sum(coexp_array_binary)
    beta = coexp_region_num / df_1.shape[0]
    if cooccur_scen == "cell":
        coexp_num = coexp_cell_num
    elif cooccur_scen == "region":
        coexp_num = coexp_region_num
    return coexp_df, total_num, coexp_num, beta

argv = sys.argv
tissue_name = argv[1]
frag_type = argv[2]
bin_size = argv[3]
process = argv[4]
region_selection = argv[5]
cell_list_dir_filename = argv[6]
cell_reads_thres = int(argv[7])
region_cell_dir = argv[8]
cooccur_scen = argv[9]
row_cluster_sel = argv[10]
coexp_num_thres = int(argv[11])
observation_all_dir = argv[12]
observation_avail_dir = argv[13]
observation_summary_dir = argv[14]

row_cluster_dir = "/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/data_summary"
row_cluster_filename = "complexheatmap_cluster_row_km-15_row_label.tsv"
row_cluster_dir_filename = os.path.join(row_cluster_dir, row_cluster_filename)
row_cluster = pd.read_csv(row_cluster_dir_filename, sep="\t", header=0, index_col=None)

cell_list_df = pd.read_csv(cell_list_dir_filename, sep="\t", header=None, index_col=None)
cell_list = cell_list_df.loc[:, 0].values.tolist()

if (row_cluster_sel == "all"):
    row_sel_array = row_cluster.loc[:, "pos"].values
else:
    row_sel_array = row_cluster.loc[row_cluster.loc[:, "label"] == row_cluster_sel, "pos"].values

wgc_dir = "/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/data_summary"
col_culster_filename = "V1V2_col_cluster_row_outsource_km-15_column_km-16_K562_available.tsv"
col_cluster_dir_filename = os.path.join(wgc_dir, col_culster_filename)
col_cluster = pd.read_csv(col_cluster_dir_filename, sep='\t', header=0, index_col=None)
target_pair_list = col_cluster.loc[:, "target_pair"].values.tolist()
target_pair_pair_list = list(combinations(target_pair_list, 2))

target_pair_pair_num = len(target_pair_pair_list)
target_pair_pair_df = pd.DataFrame(target_pair_pair_list, columns=["target_pair_1", "target_pair_2"])
target_pair_pair_df.loc[:, "total_num"] = 0
target_pair_pair_df.loc[:, "coexp_num"] = 0

for target_pair_pair_idx in range(target_pair_pair_num):
    target_pair_tuple = target_pair_pair_list[target_pair_pair_idx]
    target_pair_1 = target_pair_tuple[0]
    target_pair_2 = target_pair_tuple[1]

    wgc_1_filename = f"{tissue_name}_{target_pair_1}_{process}_{frag_type}_{bin_size}_{region_selection}.feather"
    wgc_1_dir_filename = os.path.join(region_cell_dir, wgc_1_filename)
    wgc_1 = feather.read_feather(wgc_1_dir_filename)
    wgc_1 = wgc_1.set_index(keys="pos", drop=True, inplace=False)
    wgc_1 = wgc_1.loc[:, cell_list]
    wgc_1_cell_list = list(wgc_1.columns)

    wgc_2_filename = f"{tissue_name}_{target_pair_2}_{process}_{frag_type}_{bin_size}_{region_selection}.feather"
    wgc_2_dir_filename = os.path.join(region_cell_dir, wgc_2_filename)
    wgc_2 = feather.read_feather(wgc_2_dir_filename)
    wgc_2 = wgc_2.set_index(keys="pos", drop=True, inplace=False)
    wgc_2 = wgc_2.loc[:, cell_list]
    wgc_2_cell_list = list(wgc_2.columns)
    wgc_2_cell_list_pre_perm = copy.deepcopy(wgc_2_cell_list)

    wgc_cell_list = list(set(wgc_1_cell_list) & set(wgc_2_cell_list))
    wgc_1_sel = wgc_1.loc[row_sel_array, wgc_cell_list]
    wgc_2_sel = wgc_2.loc[row_sel_array, wgc_cell_list]

    _, total_num, coexp_num, _ = coexpr_generation(wgc_1_sel, wgc_2_sel, cooccur_scen=cooccur_scen)
    target_pair_pair_df.loc[target_pair_pair_idx, "total_num"] = total_num
    target_pair_pair_df.loc[target_pair_pair_idx, "coexp_num"] = coexp_num

target_pair_avail_df = target_pair_pair_df.loc[target_pair_pair_df.loc[:, "coexp_num"] >= coexp_num_thres, ["target_pair_1", "target_pair_2"]]
coexp_num_summary_df = pd.DataFrame(target_pair_pair_df.loc[:, "coexp_num"].value_counts(), columns=["count"])
coexp_num_summary_df = coexp_num_summary_df.reset_index(drop=False, inplace=False, names="coexp_num")
total_num_summary_df = pd.DataFrame(target_pair_pair_df.loc[:, "total_num"].value_counts(), columns=["count"])
total_num_summary_df = total_num_summary_df.reset_index(drop=False, inplace=False, names="total_num")

target_pair_pair_df_filename = f"coexp_total_num_cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{row_cluster_sel}.tsv"
target_pair_pair_df_dir_filename = os.path.join(observation_all_dir, target_pair_pair_df_filename)
target_pair_pair_df.to_csv(target_pair_pair_df_dir_filename, sep="\t", header=True, index=False)

target_pair_avail_df_filename = f"target_pair_avail_cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{row_cluster_sel}_coexp_num_thres-{coexp_num_thres}.tsv"
target_pair_avail_df_dir_filename = os.path.join(observation_avail_dir, target_pair_avail_df_filename)
target_pair_avail_df.to_csv(target_pair_avail_df_dir_filename, sep="\t", header=True, index=False)

coexp_num_summary_df_filename = f"coexp_num_summary_cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{row_cluster_sel}.tsv"
observation_summary_coexp_num_dir = os.path.join(observation_summary_dir, "coexp_num")
os.makedirs(observation_summary_coexp_num_dir, exist_ok=True) 
coexp_num_summary_df_dir_filename = os.path.join(observation_summary_coexp_num_dir, coexp_num_summary_df_filename)
coexp_num_summary_df.to_csv(coexp_num_summary_df_dir_filename, sep="\t", header=True, index=False)

total_num_summary_df_filename = f"total_num_summary_cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{row_cluster_sel}.tsv"
observation_summary_total_num_dir = os.path.join(observation_summary_dir, "total_num")
os.makedirs(observation_summary_total_num_dir, exist_ok=True)
total_num_summary_df_dir_filename = os.path.join(observation_summary_total_num_dir, total_num_summary_df_filename)
total_num_summary_df.to_csv(total_num_summary_df_dir_filename, sep="\t", header=True, index=False)
