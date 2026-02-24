import pandas as pd
import numpy as np
import os
import pyarrow.feather as feather
import copy
import sys

argv = sys.argv
tissue_name = argv[1]
target_pair_name = argv[2]
process = argv[3]
sc_root_dir = argv[4]
target_reads_thres = argv[5]
bin_size = argv[6]
target_qc_type = argv[7]

# target_pair_qc_dir = os.path.join(sc_root_dir, "data_summary")
# tissue_name = "brain"
# cell_reads_thres = 100
# bin_size = 800
frag_type = "mixed"
# target_qc_type = f"all-36x-qc-{target_reads_thres}"
# target_qc_type = "all"
# process = "libnorm"
# target_pair_name = "H3K4me3-H3K4me3"
wgc_dir = os.path.join(sc_root_dir, "wgc", tissue_name, frag_type, str(bin_size))
region_cell_dir = os.path.join(sc_root_dir, "wgc", tissue_name, frag_type, "region_cell", bin_size)

cell_list_dir = os.path.join(sc_root_dir, tissue_name, "cell_list")
# cell_list_filename = f"cell_list_36x_qc-{cell_reads_thres}.txt"
cell_list_filename = "cell_list.txt"
cell_list_dir_filename = os.path.join(cell_list_dir, cell_list_filename)
cell_list_df = pd.read_csv(cell_list_dir_filename, sep="\t", header=None, index_col=None)
cell_list = cell_list_df.loc[:, 0].values.tolist()

cell_name_initial = cell_list[0]
wgc_initial_filename = f"{cell_name_initial}_{bin_size}_colQC-{target_qc_type}_{process}.feather"
wgc_initial_dir_filename = os.path.join(wgc_dir, wgc_initial_filename)
wgc_initial = feather.read_feather(wgc_initial_dir_filename)
wgc_initial = wgc_initial.set_index(keys="pos", inplace=False)
cell_gb = pd.DataFrame(np.zeros((wgc_initial.shape[0], len(cell_list))), columns=cell_list, index=list(wgc_initial.index))

for cell_name_i in cell_list:
    wgc_tmp_filename = f"{cell_name_i}_{bin_size}_colQC-{target_qc_type}_{process}.feather"
    wgc_tmp_dir_filename = os.path.join(wgc_dir, wgc_tmp_filename)
    wgc_tmp = feather.read_feather(wgc_tmp_dir_filename)
    wgc_tmp = wgc_tmp.set_index(keys="pos", inplace=False)
    cell_gb.loc[:, cell_name_i] = wgc_tmp.loc[:, target_pair_name]
    print(cell_name_i)

# pos_list_raw = list(cell_gb.index)
# pos_list = [i.replace('_', ':', 1) for i in pos_list_raw]
# cell_gb.index = pos_list

cell_gb_noallzero_row = cell_gb[(cell_gb.T != 0).any()]
print(cell_gb_noallzero_row.shape)
cell_gb_noallzero = cell_gb_noallzero_row.loc[:, (cell_gb_noallzero_row != 0).any(axis=0)]
print(cell_gb_noallzero.shape)
cell_gb_noallzero = cell_gb_noallzero.reset_index(drop=False, inplace=False)
cell_gb_noallzero = cell_gb_noallzero.rename(columns={"index": "pos"}, inplace=False)
cell_gb_noallzero_filename = f"{tissue_name}_{target_pair_name}_{process}_noAllZero_{frag_type}_{bin_size}.feather"
cell_gb_noallzero_dir_filename = os.path.join(region_cell_dir, cell_gb_noallzero_filename)
feather.write_feather(cell_gb_noallzero, cell_gb_noallzero_dir_filename)

cell_gb = cell_gb.reset_index(drop=False, inplace=False)
cell_gb = cell_gb.rename(columns={"index": "pos"}, inplace=False)
cell_gb_filename = f"{tissue_name}_{target_pair_name}_{process}_{frag_type}_{bin_size}.feather"
cell_gb_dir_filename = os.path.join(region_cell_dir, cell_gb_filename)
feather.write_feather(cell_gb, cell_gb_dir_filename)
print(cell_gb.shape)
