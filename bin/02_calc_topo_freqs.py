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
    import ete3
    from ete3 import Tree
except ImportError:
    sys.exit("One of the required modules can't be found...")

#########################
## Input - parse arguments
#########################

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--combo", help="P1_P2_P3_O, character string, will use order in string to specify P1/P2/P3!!", action = "store", required=True)
parser.add_argument("-p", "--phy_list", help="List with paths to phylogenies to be included", action = "store", nargs='+', required=True)
parser.add_argument("-i", "--chromo_info_fn", help="TSV file with chromo information, incl. following columns; chromo | start | end | mid | start_rel | end_rel | mid_rel", action = "store", required=True)
parser.add_argument("-s", "--sex_chromo_list", help="List with sex chromosomes, needs to treated differently for plots", action = "store", nargs='+', required=False, default=[])
parser.add_argument("-b", "--bin_size", help="To visualise the ratio of topologies across chromosomes, we will calculate topo frequency per bin of size", action = "store", required=False, default=1000000)
parser.add_argument("-o", "--out_dir", help="Directory to store output files with gene tree frequencies, plots, etc.", default='current', action = "store", required=False)

args = parser.parse_args()

## Store variables
p1_p2_p3_o = args.combo
phy_list = args.phy_list
chromo_length_info_fn = args.chromo_info_fn
sex_chromo_list = args.sex_chromo_list
bin_size = int(args.bin_size)
output_path = args.out_dir
if output_path != 'current':
	pass
else:
	output_path = os.getcwd()

#########################
## Get topologies and information about chromosomes
#########################

## P1/P2/P3 and O, as frequently referred to with abba-baba
p1 = p1_p2_p3_o.rsplit('_')[0]
p2 = p1_p2_p3_o.rsplit('_')[1]
p3 = p1_p2_p3_o.rsplit('_')[2]
out = p1_p2_p3_o.rsplit('_')[3]

## All rooted trees possible with the four taxa above
all_trees = dict()
# BBAA
t1 = Tree("(((%s,%s),%s),%s);" % (p1, p2, p3, out))
# ABBA
t2 = Tree("(((%s,%s),%s),%s);" % (p2, p3, p1, out))
# BABA
t3 = Tree("(((%s,%s),%s),%s);" % (p1, p3, p2, out))

# Root all possible trees on the outgroup and add to a list of trees
t1.set_outgroup(t1&"%s" %(out))
all_trees['BBAA'] = t1
t2.set_outgroup(t2&"%s" %(out))
all_trees['ABBA'] = t2
t3.set_outgroup(t3&"%s" %(out))
all_trees['BABA'] = t3

chromo_info_df = pd.read_csv(chromo_length_info_fn, header=0, sep="\t")
chromo_length_dict = dict(zip(chromo_info_df["chromo"],chromo_info_df["end"]))
# Changing the relative position of each chromo/scaffold with one binsize to improve plotting clarity (i.e. no overlapping bins)
chromo_start_rel_dict_tmp = dict(zip(chromo_info_df["chromo"],chromo_info_df["start_rel"]))
chromo_start_rel_dict = {}
count = 0
for chromo in chromo_start_rel_dict_tmp:
	chromo_start_rel_dict[chromo] = chromo_start_rel_dict_tmp[chromo] + count
	count += (bin_size * 1.5)

#########################
## Estimate GT frequency distributions
#########################

## Create a list of phylogenies 
phy_fns = []
for phy in phy_list:
	phy_fns.append(phy)


## Now loop over phylogenies, parse the gene trees and identify topology
# Create a pandas dataframe to store the relevant information for each topology
df_WIDE = pd.DataFrame(columns=['Chromosome','Start', 'End', 'Combination', 'P1', 'P2', 'P3', 'Out', 'BBAA', 'ABBA', 'BABA'])
for tree_fn in phy_fns:
    # Meta information for phylogeny can be found in phylogeny fn, format: chr15_10100000_10110000_filt_indivs_aln.treefile
    basename = os.path.basename(tree_fn)
    chromo = basename.rsplit('_')[0]
    window_start = basename.rsplit('_')[1]
    window_end = basename.rsplit('_')[2]
    # parse tree
    gt_string = open(tree_fn).readline().rstrip()
    gt = Tree(gt_string, format = 1)
    gt.set_outgroup(gt&"%s" %(out))
    match = []
    for topo in all_trees:
        rf = all_trees[topo].robinson_foulds(gt, unrooted_trees = False)
        if (rf[0] == 0):
            match.append(topo)
        else:
            pass
    if (len(match) == 0):
        print(tree_fn)
        sys.exit("One of the genetrees does not match one of the known topologies...")
    elif (len(match) > 1):
        print(tree_fn)
        sys.exit("One of the genetrees matches more than one of the known topologies...")
    else:
        # Score topology hit +1 
        if match[0] == 'BBAA':
            topo_df = pd.DataFrame([[chromo, window_start, window_end, p1_p2_p3_o, p1, p2, p3, out, 1, 0, 0]], columns=['Chromosome','Start', 'End', 'Combination', 'P1', 'P2', 'P3', 'Out', 'BBAA', 'ABBA', 'BABA'])
        elif match[0] == 'ABBA':
            topo_df = pd.DataFrame([[chromo, window_start, window_end, p1_p2_p3_o, p1, p2, p3, out, 0, 1, 0]], columns=['Chromosome','Start', 'End', 'Combination', 'P1', 'P2', 'P3', 'Out', 'BBAA', 'ABBA', 'BABA'])
        elif match[0] == 'BABA':
            topo_df = pd.DataFrame([[chromo, window_start, window_end, p1_p2_p3_o, p1, p2, p3, out, 0, 0, 1]], columns=['Chromosome','Start', 'End', 'Combination', 'P1', 'P2', 'P3', 'Out', 'BBAA', 'ABBA', 'BABA'])
        else:
            print(match[0])
            sys.exit("Unknown topology?!")
        df_WIDE = pd.concat([df_WIDE, topo_df], ignore_index=True)

## Store output as a tab-delim file
tsv_fn = os.path.join(output_path, (p1_p2_p3_o + '_gene_tree_freqs_WIDE.tsv'))
df_WIDE.to_csv(tsv_fn, sep='\t', index=False)

## Finally add some additional information for downstream plotting purposes (relative positions, chromosome type)
chromo_type_dict = {}
for chromo in chromo_info_df["chromo"]:
    if chromo not in sex_chromo_list:
        chromo_type_dict[chromo] = 'autosomes'
    else:
        chromo_type_dict[chromo] = 'sex_chromosomes'

df_WIDE.insert(1, "Chromo_type", df_WIDE['Chromosome'].map(chromo_type_dict))

df_WIDE['Start']= df_WIDE['Start'].astype('int64')
df_WIDE['End']= df_WIDE['End'].astype('int64')

df_WIDE.insert (4, "Start_rel", df_WIDE['Start'] + df_WIDE['Chromosome'].map(chromo_start_rel_dict))
df_WIDE.insert (5, "Mid_rel", (((df_WIDE['End'] - df_WIDE['Start'])/2) + df_WIDE['Start']) + df_WIDE['Chromosome'].map(chromo_start_rel_dict))
df_WIDE.insert (6, "End_rel", df_WIDE['End'] + df_WIDE['Chromosome'].map(chromo_start_rel_dict))

## In addition to the wide dataframe format, it is also helpful to have the same dataframe in a 'LONG' format (where topo names are used rather than tabular score info)
df_tmp = pd.melt(df_WIDE, id_vars=['Chromosome', 'Chromo_type', 'Start', 'End', 'Start_rel', 'Mid_rel', 'End_rel','Combination', 'P1', 'P2','P3','Out'], value_vars=['BBAA', 'ABBA', 'BABA'])
df_LONG = df_tmp[df_tmp['value'] > 0].sort_values(by=['Chromosome', 'Start']).rename(columns={"variable": "Topology"}).drop(['value'], axis=1)

## To differentiate between sex and autosomes create two separate df's
df_WIDE_auto = df_WIDE[df_WIDE['Chromo_type']=='autosomes']
df_WIDE_sex = df_WIDE[df_WIDE['Chromo_type']=='sex_chromosomes']

## Let's create two additional DFs that are grouped by chromosome
df_WIDE_grouped = df_WIDE.groupby(('Chromosome'))
df_LONG_grouped = df_LONG.groupby(('Chromosome'))


#########################
## Plot GT topology number and percentages
## 
## Both in terms of absolute number as well as a percentage
## Sex and autosomes separately
#########################

## Summarise the topology numbers and percentages in a summary dataframe (for plotting purposes)
## Bar chart with the distribution of parental alleles by genotype:
gt_number = df_WIDE_auto.shape[0]
d_tmp = {'topology': ['BBAA', 'ABBA', 'BABA'], 'number': [df_WIDE_auto['BBAA'].sum(), df_WIDE_auto['ABBA'].sum(), df_WIDE_auto['BABA'].sum()], 'percentage': [(df_WIDE_auto['BBAA'].sum()/gt_number)*100, (df_WIDE_auto['ABBA'].sum()/gt_number)*100, (df_WIDE_auto['BABA'].sum()/gt_number)*100]}
df_sum_auto = pd.DataFrame(data=d_tmp)
gt_number = df_WIDE_sex.shape[0]
d_tmp = {'topology': ['BBAA', 'ABBA', 'BABA'], 'number': [df_WIDE_sex['BBAA'].sum(), df_WIDE_sex['ABBA'].sum(), df_WIDE_sex['BABA'].sum()], 'percentage': [(df_WIDE_sex['BBAA'].sum()/gt_number)*100, (df_WIDE_sex['ABBA'].sum()/gt_number)*100, (df_WIDE_sex['BABA'].sum()/gt_number)*100]}
df_sum_sex = pd.DataFrame(data=d_tmp)

## Create a figure with two columns
sns.set_style("whitegrid", {'axes.grid' : False})
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=(10,10))
fig.subplots_adjust(hspace=0.2)

## Specify a custom colour palette for consistency across plots
bar_palette = {'BBAA': '#9b59b6',
               'ABBA': '#3498db',
               'BABA': '#2ecc71'
              }

## Create a list of the xlabs, to be able to use the set_xticklabels call
xlabs_auto = [df_sum_auto.topology[0], df_sum_auto.topology[1], df_sum_auto.topology[2]]
xlabs_sex = [df_sum_sex.topology[0], df_sum_sex.topology[1], df_sum_sex.topology[2]]
col_width = 0.8

## Plot gene tree number for both auto and sex chromosomes
df_sum_auto.plot.bar(x="topology",
				y="number",
				color=df_sum_auto['topology'].replace(bar_palette),
				ax=axes[0],
				width = col_width,
				legend=False)

df_sum_sex.plot.bar(x="topology",
				y="number",
				color=df_sum_sex['topology'].replace(bar_palette),
				ax=axes[1],
				width = col_width,
				legend=False)

## Format labels and titles
axes[0].set_title("Autosomes", fontweight="bold", size=16)
axes[1].set_title("Sex chromosomes", fontweight="bold", size=16)
axes[0].set_xticklabels(xlabs_auto, fontsize=12, rotation=45)
axes[1].set_xticklabels(xlabs_sex, fontsize=12, rotation=45)
axes[0].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)
axes[1].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)

axes[0].set_ylabel("Gene tree number", fontweight="bold", fontsize=12, labelpad=20)

## Remove top and right borders
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

## Create a character string for each topology, to include a legend with individual ID
topo1 = "(((%s,%s),%s),%s)" % (p1, p2, p3, out)
topo2 = "(((%s,%s),%s),%s)" % (p2, p3, p1, out)
topo3 = "(((%s,%s),%s),%s)" % (p1, p3, p2, out)

# Include custom legend
BBAA = mpatches.Patch(color='#9b59b6', label=topo1)
ABBA = mpatches.Patch(color='#3498db', label=topo2)
BABA = mpatches.Patch(color='#2ecc71', label=topo3)
axes[1].legend(handles=[BBAA, ABBA, BABA], bbox_to_anchor=(1.05, 0.8), fontsize=14, frameon=False)

# Save output
plt.savefig(os.path.join(output_path, (p1_p2_p3_o + '_GT_NUMBER.png')), bbox_inches = "tight", dpi=800)


#########################
## Now repeat for percentages rather than absolute counts
#########################
## Create a figure with two columns
sns.set_style("whitegrid", {'axes.grid' : False})
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(10,10))
fig.subplots_adjust(hspace=0.2)

## Plot gene tree percentages for both auto and sex chromosomes
df_sum_auto.plot.bar(x="topology",
				y="percentage",
				color=df_sum_auto['topology'].replace(bar_palette),
				ax=axes[0],
				width = col_width,
				legend=False)

df_sum_sex.plot.bar(x="topology",
				y="percentage",
				color=df_sum_auto['topology'].replace(bar_palette),
				ax=axes[1],
				width = col_width,
				legend=False)

## Format labels and titles
axes[0].set_title("Autosomes", fontweight="bold", size=16)
axes[1].set_title("Sex chromosomes", fontweight="bold", size=16)
axes[0].set_xticklabels(xlabs_auto, fontsize=12, rotation=45)
axes[1].set_xticklabels(xlabs_sex, fontsize=12, rotation=45)
axes[0].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)
axes[1].set_xlabel("Topology", fontweight="bold", fontsize=12, labelpad=20)

axes[0].set_ylim([0, 100])

axes[0].set_ylabel("Gene tree percentage", fontweight="bold", fontsize=12, labelpad=20)

## Remove top and right borders
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

# Include custom legend (handles have been defined above)
axes[1].legend(handles=[BBAA, ABBA, BABA], bbox_to_anchor=(1.05, 0.8), fontsize=14, frameon=False)

# Save output
plt.savefig(os.path.join(output_path, (p1_p2_p3_o + '_GT_PERCENTAGE.png')), bbox_inches = "tight", dpi=800)


#########################
##
## Gene tree distribution plots
##
## Visualise the distribution of gene tree topologies across chromosomes.
## The challenge is that it can be difficult to read when there are so many data points.
##
## An option is to do some form of stacked binning
#########################

# Create a figure with three panels
sns.set_style("whitegrid", {'axes.grid' : False})
fig, axes = plt.subplots(nrows=3, sharex=False, figsize=(19.6,19.6))
fig.subplots_adjust(hspace=0.2)

# Seperate the chromos so they are kind of equal lenght in plot
set_0 = ['chr1A', 'chr1', 'chr2']
set_1 = ['chr3', 'chr4A', 'chr4', 'chr5', 'chr6', 'chr7']

# I will store the name of the chromo, the position along the x-axis and if need be the start/end position so I can include shading
x_labels_0 = []
x_labels_1 = []
x_labels_2 = []
x_labels_pos_0 = []
x_labels_pos_1 = []
x_labels_pos_2 = []
shade_coords_0 = []
shade_coords_1 = []
shade_coords_2 = []

for num, (chromo, chromo_df) in enumerate(df_LONG_grouped):
    if chromo in set_0:
        sns.histplot(data=chromo_df, 
                    x="Mid_rel", hue="Topology", stat="count", binwidth = bin_size, hue_order = ["BABA", "ABBA","BBAA"],
                    multiple="stack", palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
                    ax=axes[0], legend=False, zorder=2)
        x_labels_0.append(chromo)
        x_labels_pos_0.append(int(((chromo_df['End_rel'].iloc[-1] - chromo_df['Start_rel'].iloc[0])/2) + chromo_df['Start_rel'].iloc[0]))
        shade_coords_0.append([chromo_df['Start_rel'].iloc[0], chromo_df['End_rel'].iloc[-1]])
    elif chromo in set_1:
        sns.histplot(data=chromo_df, 
                    x="Mid_rel", hue="Topology", stat="count", binwidth = bin_size, hue_order = ["BABA", "ABBA","BBAA"],
                    multiple="stack", palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
                    ax=axes[1], legend=False, zorder=2)
        x_labels_1.append(chromo)
        x_labels_pos_1.append(int(((chromo_df['End_rel'].iloc[-1] - chromo_df['Start_rel'].iloc[0])/2) + chromo_df['Start_rel'].iloc[0]))
        shade_coords_1.append([chromo_df['Start_rel'].iloc[0], chromo_df['End_rel'].iloc[-1]])
    else:
        sns.histplot(data=chromo_df, 
                    x="Mid_rel", hue="Topology", stat="count", binwidth = bin_size, hue_order = ["BABA", "ABBA","BBAA"],
                    multiple="stack", palette=dict(BBAA="#9b59b6", ABBA="#3498db", BABA="#2ecc71"),
                    ax=axes[2], legend=False, zorder=2)
        x_labels_2.append(chromo)
        x_labels_pos_2.append(int(((chromo_df['End_rel'].iloc[-1] - chromo_df['Start_rel'].iloc[0])/2) + chromo_df['Start_rel'].iloc[0]))
        shade_coords_2.append([chromo_df['Start_rel'].iloc[0], chromo_df['End_rel'].iloc[-1]])

# If needed, include shaded boxes that represent chromosomes 
# for l in shade_coords_0[0::2]:
#     axes[0].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.25, zorder=1)
# for l in shade_coords_1[0::2]:
#     axes[1].axvspan(l[0], l[1], facecolor='#dbd9d9', alpha=0.25, zorder=1)
# for l in shade_coords_2[0::2]:
#     axes[2].axvspan(l[0], l[1],  facecolor='#dbd9d9', alpha=0.25, zorder=1)

# The main plot is now finished but I want to format the image
axes[0].set_xticks(x_labels_pos_0)
axes[1].set_xticks(x_labels_pos_1)
axes[2].set_xticks(x_labels_pos_2)
axes[0].set_xticklabels(x_labels_0, fontweight="bold", fontsize=12, rotation=45)
axes[1].set_xticklabels(x_labels_1, fontweight="bold", fontsize=12, rotation=45)
axes[2].set_xticklabels(x_labels_2, fontweight="bold", fontsize=12, rotation=45)
sns.despine()

# Format labels
axes[0].set_ylabel("")
axes[1].set_ylabel('Topology count \n per ' + str(bin_size) + 'bp', fontweight="bold", fontsize=18, labelpad=20)
axes[2].set_ylabel("")
axes[0].set_xlabel("")
axes[1].set_xlabel("")
axes[2].set_xlabel('Chromosomes', fontweight="bold", fontsize=18, labelpad=20)

# Include custom legend (handles were already created above, since we always use the same colour scheme this shouldn't change)
axes[1].legend(handles=[BBAA, ABBA, BABA], bbox_to_anchor=(1.05, 0.8), fontsize=18, frameon=False)

# Save output
plt.savefig(os.path.join(output_path, (p1_p2_p3_o + '_GT_ratio_DISTRIBUTION.png')), bbox_inches = "tight", dpi=800)
