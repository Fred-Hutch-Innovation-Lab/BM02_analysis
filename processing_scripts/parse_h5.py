import h5py
import argparse
import pandas as pd
# import numpy as np
import json
# import glob
import os

kits = ['GEMX', 'NextGEM', 'flex']
samples = ['F1A', 'F1B', 'F5A', 'F5B']

json_data_out = {}
for kit in kits:
    i = 1
    for sample in samples:
        if kit == 'GEMX' :
            h5_path = os.path.join('/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/',
                               kit.lower() + '_3p', 'downsampled_runs',
                               kit + '_' + sample, 'outs/per_sample_outs',
                               kit + '_' + sample, 'count/sample_molecule_info.h5')
            f = h5py.File(h5_path, "r")
            json_data = json.loads(f['metrics_json'][()])
            # print(kit, sample)
            # print(json_data['libraries']['0']['feature_read_pairs'])
            # print(json_data['libraries']['0']['usable_read_pairs'])
            json_data_out[(kit, sample)] = json_data['libraries']['0']
        elif kit == 'NextGEM':
            h5_path = os.path.join('/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/',
                               kit.lower() + '_3p', 'downsampled_runs',
                               kit + '_' + sample, 'outs/per_sample_outs',
                               kit + '_' + sample, 'count/sample_molecule_info.h5')
            f = h5py.File(h5_path, "r")
            json_data = json.loads(f['metrics_json'][()])
            # print(kit, sample)
            # print(json_data['libraries']['0']['feature_read_pairs'])
            # print(json_data['libraries']['0']['usable_read_pairs'])
            json_data_out[(kit, sample)] = json_data['libraries']['0']
            # json_data_out = pd.concat([json_data_out, pd.DataFrame([json_data['libraries']['0']], index=[kit + '_' + sample])], ignore_index=True)
        else:
            h5_path = os.path.join('/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/',
                               kit.lower(), 'full_data_runs/per_sample_outs/', 
                               'B2-FL-' + sample + '_BC0' + str(i),
                               'count/sample_molecule_info.h5')
        
            f = h5py.File(h5_path, "r")
            json_data = json.loads(f['metrics_json'][()])
            # print(kit, sample)
            # print(json_data['libraries']['0']['feature_read_pairs'])
            # print(json_data['libraries']['0']['usable_read_pairs'])
            # print(json_data['libraries']['0']['on_target_usable_read_pairs'])
            json_data_out[(kit, sample)] = json_data['libraries']['0']
            # json_data_out = pd.concat([json_data_out, pd.DataFrame([json_data['libraries']['0']], index=[kit + '_' + sample])], ignore_index=True)
            i += 1

df = pd.DataFrame.from_dict(json_data_out, orient='index')
print(df)
df.to_csv('/fh/fast/_IRC/FHIL/grp/BM_paper/analysis/data/3p/10x_read_utilization.csv')

kits = ['GEMX', 'NextGEM']
samples = ['F1', 'F2', 'F4', 'F5']

json_data_out = {}
for kit in kits:
    i = 1
    for sample in samples:
        if kit == 'GEMX' :
            h5_path = os.path.join('/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/',
                               kit.lower() + '_5p', 'downsampled',
                               kit + '_' + sample, 'outs/per_sample_outs',
                               kit + '_' + sample, 'count/sample_molecule_info.h5')
            f = h5py.File(h5_path, "r")
            json_data = json.loads(f['metrics_json'][()])
            # print(kit, sample)
            # print(json_data['libraries']['0']['feature_read_pairs'])
            # print(json_data['libraries']['0']['usable_read_pairs'])
            json_data_out[(kit, sample)] = json_data['libraries']['0']
        elif kit == 'NextGEM':
            h5_path = os.path.join('/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/',
                               kit.lower() + '_5p', 'downsampled',
                               kit + '_' + sample, 'outs/per_sample_outs',
                               kit + '_' + sample, 'count/sample_molecule_info.h5')
            f = h5py.File(h5_path, "r")
            json_data = json.loads(f['metrics_json'][()])
            # print(kit, sample)
            # print(json_data['libraries']['0']['feature_read_pairs'])
            # print(json_data['libraries']['0']['usable_read_pairs'])
            json_data_out[(kit, sample)] = json_data['libraries']['0']
            # json_data_out = pd.concat([json_data_out, pd.DataFrame([json_data['libraries']['0']], index=[kit + '_' + sample])], ignore_index=True)


df = pd.DataFrame.from_dict(json_data_out, orient='index')
print(df)
df.to_csv('/fh/fast/_IRC/FHIL/grp/BM_paper/analysis/data/5p/10x_read_utilization.csv')
