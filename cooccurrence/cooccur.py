import pyarrow.feather as feather
import matplotlib.pyplot as plt
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
perm_num = int(argv[9])
cooccur_scen = argv[10]
row_cluster_sel = argv[11]
target_pair_1 = argv[12]
target_pair_2 = argv[13]
test_result_dir = argv[14]
region_dir = argv[15]
fig_dir = argv[16]

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
    
target_pair_coexpress_df, _, coexp_num, beta = coexpr_generation(wgc_1_sel, wgc_2_sel, cooccur_scen=cooccur_scen)
target_pair_coexpress_df_filename = f"cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{target_pair_1}_{target_pair_2}_{row_cluster_sel}_{cooccur_scen}_cooccur_region.tsv"
target_pair_coexpress_df_dir_filename = os.path.join(region_dir, target_pair_coexpress_df_filename)
target_pair_coexpress_df.to_csv(target_pair_coexpress_df_dir_filename, sep="\t", header=True, index=True)
    
coexp_num_perm_list = []
for perm_i in range(perm_num):
    wgc_2_cell_list_perm = copy.deepcopy(wgc_2_cell_list_pre_perm)
    random.seed(perm_i)
    shuffle(wgc_2_cell_list_perm)
    wgc_2_sel_perm = wgc_2_sel.loc[:, wgc_2_cell_list_perm]
    wgc_2_sel_perm.columns = wgc_2_cell_list
    _, _, coexp_num_perm, _ = coexpr_generation(wgc_1_sel, wgc_2_sel_perm, cooccur_scen=cooccur_scen)
    coexp_num_perm_list.append(coexp_num_perm)
    
perm_mean = np.mean(coexp_num_perm_list)
perm_sum = np.sum(coexp_num_perm_list)
coexp_num_perm_array = np.asarray(coexp_num_perm_list)
p = np.sum(coexp_num_perm_array >= coexp_num) / perm_num
alter_event_num = np.sum(coexp_num_perm_array >= coexp_num)
log10_p = -np.log10(p)
coexp_summary_array = np.array([[target_pair_1, target_pair_2, coexp_num, perm_sum, perm_mean, alter_event_num, p, log10_p, beta]])
coexp_summary_columns = ["target_pair_1", "target_pair_2", f"1|1_{cooccur_scen}", f"1|1_{cooccur_scen}_sim_sum", f"1|1_{cooccur_scen}_sim", "alt_num", "pval", "-log10(pval)", "beta"]
coexp_summary_df = pd.DataFrame(coexp_summary_array, columns=coexp_summary_columns)
coexp_summary_filename = f"cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{target_pair_1}_{target_pair_2}_{row_cluster_sel}_{cooccur_scen}.tsv"
coexp_summary_dir_filename = os.path.join(test_result_dir, coexp_summary_filename)
coexp_summary_df.to_csv(coexp_summary_dir_filename, sep="\t", header=True, index=False)

plt.figure()
plt.hist(coexp_num_perm_list, bins=30)
plt.axvline(x=coexp_num, color='r', label=f"1|1, {coexp_num}")
plt.axvline(x=perm_mean, color='b', label=f"null hyp, {perm_mean:.1f}")
plt.legend()
plt.title(f"{target_pair_1} & {target_pair_2}, {row_cluster_sel}, {cooccur_scen}\n"
          f"p-value = {p:.2e}, beta = {beta:.2e}, x {coexp_num / perm_mean:.2f}")

fig_filename = f"cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{target_pair_1}_{target_pair_2}_{row_cluster_sel}_{cooccur_scen}_cooccur_test_result.pdf"
fig_dir_filename = os.path.join(fig_dir, fig_filename)
plt.savefig(fig_dir_filename)

