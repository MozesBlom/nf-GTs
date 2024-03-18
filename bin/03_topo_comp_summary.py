#!/usr/bin/env python3

#########################
## Required modules
#########################
try:
    import os
    import sys
    import argparse
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.patches as mpatches
except ImportError:
    sys.exit("One of the required modules can't be found...")

#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--topo_freqs_tsv_list", help="List with TSV file paths, each summarising gene tree frequencies", action = "store", nargs='+', required=True)
parser.add_argument("-i", "--chromo_info_tsv_fn", help="TSV file with chromo information, incl. following columns; chromo | start | end | mid | start_rel | end_rel | mid_rel", action = "store", required=True)
parser.add_argument("-s", "--sex_chromo_list", help="List with sex chromosomes, needs to treated differently for plots", action = "store", nargs='+', required=False, default=[])
parser.add_argument("-o", "--out_dir", help="Directory to store output files with gene tree frequencies, plots, etc.", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
topo_tsv_list = args.topo_freqs_tsv_list
chromo_length_info_fn = args.chromo_info_tsv_fn
sex_chromo_list = args.sex_chromo_list
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()

#########################
## Create dataframes
#########################

## Differentiate between sex and autosomes
chromo_info_df = pd.read_csv(chromo_length_info_fn, header=0, sep="\t")

chromo_type_dict = {}
for chromo in chromo_info_df["chromo"]:
    if chromo not in sex_chromo_list:
        chromo_type_dict[chromo] = 'autosomes'
    else:
        chromo_type_dict[chromo] = 'sex_chromosomes'


## Read in all the summary dataframes and create one large pandas summary DF
dfs_list = []

for df_fn in topo_tsv_list:
    df = pd.read_csv(df_fn, header=0, sep="\t")
    dfs_list.append(df)

dfs_WIDE = pd.concat(dfs_list)
dfs_WIDE.insert(1, "Chromo_type", dfs_WIDE['Chromosome'].map(chromo_type_dict))


## Create a df for sex and autosomes seperately
dfs_WIDE_auto = dfs_WIDE[dfs_WIDE['Chromo_type']=='autosomes']
dfs_WIDE_sex = dfs_WIDE[dfs_WIDE['Chromo_type']=='sex_chromosomes']

## Group by combination
dfs_WIDE_auto_grouped = dfs_WIDE_auto.groupby(('Combination'))
dfs_WIDE_sex_grouped = dfs_WIDE_sex.groupby(('Combination'))

## Finally, iterate over each dataframe (combination) and calculate the gene tree frequencies
# Autosomes
dfs_sum_auto_list = []

for num, (combo, combo_df) in enumerate(dfs_WIDE_auto_grouped):
    gt_number = combo_df.shape[0]
    d_tmp = {'combination': [combo, combo, combo], 'topology': ['BBAA', 'ABBA', 'BABA'], 'number': [combo_df['BBAA'].sum(), combo_df['ABBA'].sum(), combo_df['BABA'].sum()], 'percentage': [(combo_df['BBAA'].sum()/gt_number)*100, (combo_df['ABBA'].sum()/gt_number)*100, (combo_df['BABA'].sum()/gt_number)*100]}
    dfs_sum_auto_list.append(pd.DataFrame(data=d_tmp))

dfs_sum_auto = pd.concat(dfs_sum_auto_list)

# Sex chromos
dfs_sum_sex_list = []

for num, (combo, combo_df) in enumerate(dfs_WIDE_sex_grouped):
    gt_number = combo_df.shape[0]
    d_tmp = {'combination': [combo, combo, combo], 'topology': ['BBAA', 'ABBA', 'BABA'], 'number': [combo_df['BBAA'].sum(), combo_df['ABBA'].sum(), combo_df['BABA'].sum()], 'percentage': [(combo_df['BBAA'].sum()/gt_number)*100, (combo_df['ABBA'].sum()/gt_number)*100, (combo_df['BABA'].sum()/gt_number)*100]}
    dfs_sum_sex_list.append(pd.DataFrame(data=d_tmp))

dfs_sum_sex = pd.concat(dfs_sum_sex_list)


#########################
## Plot a summary of the observed topologies
#########################

## Absolute numbers:
# Create a figure with two columns
sns.set_style("whitegrid", {'axes.grid' : False})
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=(10,10))
fig.subplots_adjust(hspace=0.2)

# First for autosomes
sns.stripplot(data=dfs_sum_auto,
              x="topology",
              y="number",
              jitter=True,
              alpha=0.5,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              zorder=2,
              ax=axes[0])

sns.pointplot(data=dfs_sum_auto,
              x="topology",
              y="number",
              join=False,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              markers="*",
              scale=1.2,
              ci=None,
              zorder=2,
              ax=axes[0])

# Then for sex chromosomes
sns.stripplot(data=dfs_sum_sex,
              x="topology",
              y="number",
              jitter=True,
              alpha=0.5,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              zorder=2,
              ax=axes[1])

sns.pointplot(data=dfs_sum_sex,
              x="topology",
              y="number",
              join=False,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              markers="*",
              scale=1.2,
              ci=None,
              zorder=2,
              ax=axes[1])

## Format labels and titles
axes[0].set_title("Autosomes", fontweight="bold", size=16)
axes[1].set_title("Sex chromosomes", fontweight="bold", size=16)
axes[0].set_xticklabels(["BBAA", "ABBA", "BABA"], fontsize=12, rotation=45)
axes[1].set_xticklabels(["BBAA", "ABBA", "BABA"], fontsize=12, rotation=45)
axes[0].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)
axes[1].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)

axes[0].set_ylabel("Gene tree number", fontweight="bold", fontsize=12, labelpad=20)
axes[1].set_ylabel("", fontweight="bold", fontsize=12, labelpad=20)

## Remove top and right borders
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

# Save output
plt.savefig(os.path.join(output_path, ('Summary_across_comps_GT_NUMBER.png')), bbox_inches = "tight", dpi=800)


## In percentages:
## Create a figure with two columns
sns.set_style("whitegrid", {'axes.grid' : False})
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(10,10))
fig.subplots_adjust(hspace=0.2)

# First for autosomes
sns.stripplot(data=dfs_sum_auto,
              x="topology",
              y="percentage",
              jitter=True,
              alpha=0.5,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              zorder=2,
              ax=axes[0])

sns.pointplot(data=dfs_sum_auto,
              x="topology",
              y="percentage",
              join=False,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              markers="*",
              scale=1.2,
              ci=None,
              zorder=2,
              ax=axes[0])

# Then for sex chromosomes
sns.stripplot(data=dfs_sum_sex,
              x="topology",
              y="percentage",
              jitter=True,
              alpha=0.5,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              zorder=2,
              ax=axes[1])

sns.pointplot(data=dfs_sum_sex,
              x="topology",
              y="percentage",
              join=False,
              palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
              order=["BBAA", "ABBA", "BABA"],
              markers="*",
              scale=1.2,
              ci=None,
              zorder=2,
              ax=axes[1])

## Format labels and titles
axes[0].set_title("Autosomes", fontweight="bold", size=16)
axes[1].set_title("Sex chromosomes", fontweight="bold", size=16)
axes[0].set_xticklabels(["BBAA", "ABBA", "BABA"], fontsize=12, rotation=45)
axes[1].set_xticklabels(["BBAA", "ABBA", "BABA"], fontsize=12, rotation=45)
axes[0].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)
axes[1].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)

axes[0].set_ylim([0,100])

axes[0].set_ylabel("Gene tree percentage", fontweight="bold", fontsize=12, labelpad=20)
axes[1].set_ylabel("", fontweight="bold", fontsize=12, labelpad=20)

## Remove top and right borders
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

# Save output
plt.savefig(os.path.join(output_path, ('Summary_across_comps_GT_PERCENTAGE.png')), bbox_inches = "tight", dpi=800)