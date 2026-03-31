"""
Prep input metadata.
"""

import os
import pandas as pd
import matplotlib
import plotting_utils as plu
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics import pairwise_distances
matplotlib.use('macOSX')
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/manas_heart'
path_data = os.path.join(path_main, 'data')
path_input = os.path.join(path_main, 'data/input')
path_figures = os.path.join(path_main, 'figures')

# Read all files

# Original LCM metadata
df_coords = pd.read_csv(os.path.join(path_data, 'Heart_final_coorindates_135.csv'))
df_LCM = pd.read_csv(os.path.join(path_data, 'Heart_metadata.csv'))


## 


# Septum:
# Left:'Left IVS'
# Right: 'Right IVS'
# 
# Centre: 'Centre IVS' 'Centre tip IVS' 'Centre base IVS'
# 
# Right Ventricle:
# Right ventricle
# 
# Left Ventricle:
# Left ventricle

# 1.Left septum
# 2.Right septum
# 3.Centre septum
# 
# 4.Right ventricle
# 5.Left ventricle

histo_map = {
    'Left IVS' : 'Left_septum',
    'Right IVS' : 'Right_septum',
    'Centre IVS' : 'Centre_septum',
    'Centre tip IVS' : 'Centre_septum',
    'Centre base IVS' : 'Centre_septum',
    'Right ventricle' : 'Right_Ventricle',
    'Left ventricle' : 'Left_Ventricle'
}

meta = df_LCM[['Sample_ID', 'Histo']].drop_duplicates()
histo_of_interest = list(histo_map.keys())
meta = meta.query('Histo in @histo_of_interest')
meta['regions'] = meta['Histo'].map(histo_map)
meta['chunk'] = meta['regions'].copy()

# Counts
meta['regions'].value_counts()


# Viz
fig , ax = plt.subplots(figsize=(3.5,3.5))
plu.counts_plot(meta, 'regions', ax=ax)
plu.format_ax(ax, xlabel='', ylabel='n samples', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'samples_per_region.pdf'))


##


# Chunk: Left_Ventricle
samples = meta.query('regions == "Left_Ventricle"')['Sample_ID'].values
xyz = df_coords.query('name in @samples').set_index('name').loc[samples]

D = pairwise_distances(xyz)
order = leaves_list(linkage(D))
# plt.imshow(D[order][:, order], cmap='viridis')
# plt.show()

chunks = [
    xyz.index[order][:4],
    xyz.index[order][4:9],
    xyz.index[order][9:14],
    xyz.index[order][14:],
]
for i,chunk in enumerate(chunks):
    meta.loc[meta['Sample_ID'].isin(chunk), 'chunk'] = f'Left_Ventricle.{i+1}'

# Chunk: Left_septum
samples = meta.query('regions == "Left_septum"')['Sample_ID'].values
xyz = df_coords.query('name in @samples').set_index('name').loc[samples]

D = pairwise_distances(xyz)
order = leaves_list(linkage(D))
# plt.imshow(D[order][:, order], cmap='viridis')
# plt.show()

chunks = [
    xyz.index[order][:5],
    xyz.index[order][5:10],
    xyz.index[order][10:15],
    xyz.index[order][15:],
]
for i,chunk in enumerate(chunks):
    meta.loc[meta['Sample_ID'].isin(chunk), 'chunk'] = f'Left_septum.{i+1}'

# Chunk: Right_septum
samples = meta.query('regions == "Right_septum"')['Sample_ID'].values
xyz = df_coords.query('name in @samples').set_index('name').loc[samples]

D = pairwise_distances(xyz)
order = leaves_list(linkage(D))
# plt.imshow(D[order][:, order], cmap='viridis')
# plt.show()

chunks = [
    xyz.index[order][:4].tolist() + xyz.index[order][-2:].tolist(),
    xyz.index[order][4:8],
    xyz.index[order][8:13],
]
for i,chunk in enumerate(chunks):
    meta.loc[meta['Sample_ID'].isin(chunk), 'chunk'] = f'Right_septum.{i+1}'

# Chunk: Right_Ventricle
samples = meta.query('regions == "Right_Ventricle"')['Sample_ID'].values
xyz = df_coords.query('name in @samples').set_index('name').loc[samples]

D = pairwise_distances(xyz)
order = leaves_list(linkage(D))
# plt.imshow(D[order][:, order], cmap='viridis')
# plt.show()

chunks = [
    xyz.index[order][:4],
    xyz.index[order][4:]
]
for i,chunk in enumerate(chunks):
    meta.loc[meta['Sample_ID'].isin(chunk), 'chunk'] = f'Right_Ventricle.{i+1}'


##


# Final check
meta['regions'].value_counts()
meta['chunk'].value_counts()

# Save
meta.to_csv(os.path.join(path_input, 'heart_samples.csv'), index=False)


##



# Placenta
old_meta = pd.read_csv(os.path.join(path_input, 'All_Samples_Fetal_Natalie.csv'))
old_meta = old_meta.query('Histo == "Trophoblasts"')

old_meta['Site'].value_counts()
samples = [
    *old_meta['Sample'].loc[old_meta['Site'].loc[lambda x: x=='Top right'].sample(2).index].to_list(),
    *old_meta['Sample'].loc[old_meta['Site'].loc[lambda x: x=='Bottom right'].sample(2).index].to_list(),
    *old_meta['Sample'].loc[old_meta['Site'].loc[lambda x: x=='Top left'].index].to_list(),
    *old_meta['Sample'].loc[old_meta['Site'].loc[lambda x: x=='Bottom left'].index].to_list()
]
regions = [
    'Top_right',
    'Top_right',
    'Bottom_right',
    'Bottom_right',
    'Top_left',
    'Top_left',
    'Bottom_left'
]
meta = pd.DataFrame({'Sample_ID': samples, 'regions': regions}).assign(chunk='Placenta')
meta['Histo'] = meta['regions']
(
    meta[['Sample_ID', 'Histo', 'regions', 'chunk']]
    .to_csv(os.path.join(path_input, 'placenta_samples.csv'), index=False)
)


##


# Visualize
meta = pd.read_csv(os.path.join(path_input, 'heart_samples.csv'))
meta_p = pd.read_csv(os.path.join(path_input, 'placenta_samples.csv'))
meta = pd.concat([meta, meta_p], ignore_index=True)

fig , ax = plt.subplots(figsize=(3.5,3.5))
plu.counts_plot(meta, 'chunk', ax=ax)
plu.format_ax(ax, xlabel='', ylabel='n samples', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'samples_per_chunk.pdf'))


##
