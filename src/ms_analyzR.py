# - - - - - - - - - - - - Modules - - - - - - - - - - - - #
#                    Importing modules                    #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
import os
import sys
import re
import io
import random
import shutil
import subprocess
import warnings
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.core import indexing
from pandas.core.frame import DataFrame
from pandas.core.indexes.base import Index
import subprocess

# - - - - - - - - - - - - INIT - - - - - - - - - - - - #
#    Initializing the input arguments of ms_analyzR    #
# - - - - - - - - - - - - - - - - - - - - - - - - - -- #

# Function definition:
def ms_analyzR() -> tuple:
    ap = argparse.ArgumentParser()
    # ap.add_argument("-in", "--inputdir", help="path to your input files")
    ap.add_argument("-out", "--outputdir", help="path to your output files")
    ap.add_argument("-rr", "--roary", help="path to ROARY gene_presence_abscence.csv file", type=str)
    ap.add_argument("-cds", "--coding", help="MOTIF_CDS.gff file")
    ap.add_argument("-ncds", "--noncoding", help="MOTIF_nCDS.gff file")
    ap.add_argument("-inter", "--intergenic", help="MOTIF_true_intergenic.gff file")
    ap.add_argument("-ups", "--upstream", help="MOTIF_upstream.gff file")
    ap.add_argument("-bed", "--make_bed", help="Write in [OUT] tabular per-feature file ready for RCircos", action="store_true", default=False)
    ap.add_argument("-chr", "--make_chrom", help="Rearrange your input GFFs for chromosomes", action="store_true", default=False)
    # ap.add_argument("-e", "--show_error", help="if indexer is out-of-bonds error occurs, by returning TRUE this flag you can see which gene-id is wrong and manually erase it from your GFF", action="store_true", default=False)
    return ap.parse_args()

args = ms_analyzR()

# File-object definition:
# input_dir = args.inputdir
output_dir = args.outputdir
roary_csv = args.roary
coding_ms = args.coding
noncoding_ms = args.noncoding
intergenic_ms = args.intergenic
upstream_ms = args.upstream
evo_gff = args.make_bed
split_flag = args.make_chrom
# show_error = args.show_error

# coding_compath = input_dir + coding_ms
# noncoding_compath = input_dir + noncoding_ms
# intergenic_compath = input_dir + intergenic_ms
# upstream_compath = input_dir + upstream_ms

# Output directory creation:
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    basename = os.path.basename(output_dir)
    print("\n")
    print(f"Creating the {basename} folder ..")
    # print("\n")
else:
    sys.exit("[WARNING]: The chosen directory already exists, please change folder's name")

print("Starting the analysis ..")
print("\n")

# - - - - - - - - - - - - PARSING - - - - - - - - - - - - #
#         Start file-parsing and columns-indexing         #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# "<motif>_CDS.gff" file parsing:
with open(coding_ms, 'r') as cds:
    # Init a dataframe
    cds_df = pd.read_table(cds, header=None)
    cds_df.columns=["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    cds_toremove_list = [0]
    cds_df = cds_df[~cds_df['F'].isin(cds_toremove_list)]
    # print(cds_df)
    
    # Init I-column: geneID
    gene_count_serie = cds_df['I'].value_counts()
    gene_count = gene_count_serie.to_frame(name='geneID_occurency')
    gene_number = len(gene_count.index)
    
    cds_occurrency=gene_count['geneID_occurency'].value_counts()


    #Highest Methylation Number (HMN)
    gene_HMN = gene_count['geneID_occurency'].max()
    if gene_count_serie.empty:
        sys.exit("Sorry, the analysis has been stopped: you have not enough methylations")
    cds_seqID_HMN = gene_count_serie.iloc[[0]]
    #Name of seqID with Highest Methylation Number (HMN)
    cds_seqID_HMN = cds_seqID_HMN.keys()[0]

    tot_meth_cds = gene_count.sum(axis=0)
    tot_meth_cds_df = tot_meth_cds.to_frame().T
    # print(cds_df)
    # print(gene_count)

# coding_greed_path = output_dir + '/' + "coding_greed.gff"
# cds_df.to_csv(coding_greed_path, sep="\t")

# "<motif>_nCDS.gff" file parsing:
with open(noncoding_ms, 'r') as ncds:
    # Init a dataframe
    non_cds_df = pd.read_table(ncds, header=None)
    non_cds_df.columns = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    non_cds_toremove_list = [0]
    non_cds_df = non_cds_df[~non_cds_df['F'].isin(non_cds_toremove_list)]

    # Init I-column: geneID
    ncds_count_serie = non_cds_df['I'].value_counts()
    ncds_count = ncds_count_serie.to_frame(name='geneID_occurency')
    ncds_number = len(ncds_count.index)

    ncds_occurrency = ncds_count['geneID_occurency'].value_counts()

    ncds_HMN = ncds_count['geneID_occurency'].max() # Highest Methylation Number (HMN)
    if ncds_count_serie.empty:
        sys.exit("Sorry, the analysis has been stopped: you have not enough methylations")

    ncds_seqID_HMN = ncds_count_serie.iloc[[0]]
    # Name of seqID with Highest Methylation Number (HMN)
    ncds_seqID_HMN = ncds_seqID_HMN.keys()[0]

    tot_meth_ncds = ncds_count.sum(axis=0)
    tot_meth_ncds_df = tot_meth_ncds.to_frame().T
    # print(non_cds_df)
    # print(ncds_count)


# "<motif>_intergenic.gff" file parsing:
with open(intergenic_ms, 'r') as intergenic_file:
    # Init a dataframe
    intergenic_df = pd.read_table(intergenic_file, header=None)
    intergenic_df.columns = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    intergenic_toremove_list = [0]
    intergenic_df = intergenic_df[~intergenic_df['F'].isin(intergenic_toremove_list)]

    # Init I-column: geneID
    intergenic_count_serie = intergenic_df['I'].value_counts()
    intergenic_count = intergenic_count_serie.to_frame(name='geneID_occurency')
    intergenic_number = len(intergenic_count.index)

    intergenic_occurrency = intergenic_count['geneID_occurency'].value_counts()

    # Highest Methylation Number (HMN)
    intergenic_HMN = intergenic_count['geneID_occurency'].max()
    if intergenic_count_serie.empty:
        sys.exit("Sorry, the analysis has been stopped: you have not enough methylations")

    intergenic_seqID_HMN = intergenic_count_serie.iloc[[0]]
    # Name of seqID with Highest Methylation Number (HMN)
    intergenic_seqID_HMN = intergenic_seqID_HMN.keys()[0]

    tot_meth_intergenic = intergenic_count.sum(axis=0)
    tot_meth_intergenic_df = tot_meth_intergenic.to_frame().T
    # print(upstream_df)
    # print(upstream_count


# "<motif>_upstream.gff" file parsing:
with open(upstream_ms, 'r') as upstream_file:
    # Init a dataframe
    upstream_df = pd.read_table(upstream_file, header=None)
    upstream_df.columns = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    upstream_toremove_list = [0]
    upstream_df = upstream_df[~upstream_df['F'].isin(upstream_toremove_list)]

    # Init I-column: geneID
    upstream_count_serie = upstream_df['I'].value_counts()
    upstream_count = upstream_count_serie.to_frame(name='geneID_occurency')
    upstream_number = len(upstream_count.index)

    ups_occurrency = upstream_count['geneID_occurency'].value_counts()

    upstream_HMN = upstream_count['geneID_occurency'].max() # Highest Methylation Number (HMN)
    if upstream_count_serie.empty:
        sys.exit("Sorry, the analysis has been stopped: you have not enough methylations")
    ups_seqID_HMN = upstream_count_serie.iloc[[0]]
    # Name of seqID with Highest Methylation Number (HMN)
    ups_seqID_HMN = ups_seqID_HMN.keys()[0]
    # print(ups_seqID_HMN)
    
    tot_meth_upstream = upstream_count.sum(axis=0)
    tot_meth_upstream_df = tot_meth_upstream.to_frame().T
    # print(upstream_df)
    # print(upstream_count)
##########################################################################################################################################################################
if not roary_csv:
    coding_EVO_df = cds_df[['A', 'D', 'E', 'I']].copy()
    coding_EVO_df = coding_EVO_df.drop_duplicates(subset='I').assign(
        J=coding_EVO_df.groupby('I')['I'].transform('count'))
    coding_EVO_df.columns = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score']
    coding_EVO_df_dropped = coding_EVO_df.reset_index(drop=True)
    coding_EVO_path = output_dir + '/' + "evo_CDS.bed"
    coding_EVO_df_dropped.to_csv(coding_EVO_path, sep="\t", index=False)

    # non_CDS
    non_coding_EVO_df = non_cds_df[['A', 'D', 'E', 'I']].copy()
    non_coding_EVO_df = non_coding_EVO_df.drop_duplicates(subset='I').assign(
        J=non_coding_EVO_df.groupby('I')['I'].transform('count'))
    non_coding_EVO_df.columns = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score']
    non_coding_EVO_df_dropped = non_coding_EVO_df.reset_index(drop=True)
    non_coding_EVO_path = output_dir + '/' + "evo_nCDS.bed"
    non_coding_EVO_df_dropped.to_csv(
        non_coding_EVO_path, sep="\t", index=False)

    # INTERGENIC
    inter_EVO_df = intergenic_df[['A', 'D', 'E', 'I']].copy()
    inter_EVO_df = inter_EVO_df.drop_duplicates(subset='I').assign(
        J=inter_EVO_df.groupby('I')['I'].transform('count'))
    inter_EVO_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    inter_EVO_df_dropped = inter_EVO_df.reset_index(drop=True)
    inter_EVO_path = output_dir + '/' + "evo_intergenic.bed"
    inter_EVO_df_dropped.to_csv(inter_EVO_path, sep="\t", index=False)

    # UPS
    ups_EVO_df = upstream_df[['A', 'D', 'E', 'I']].copy()
    ups_EVO_df = ups_EVO_df.drop_duplicates(subset='I').assign(
        J=ups_EVO_df.groupby('I')['I'].transform('count'))
    ups_EVO_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    ups_EVO_df_dropped = ups_EVO_df.reset_index(drop=True)
    ups_EVO_path = output_dir + '/' + "evo_upstream.bed"
    ups_EVO_df_dropped.to_csv(ups_EVO_path, sep="\t", index=False)

    # CDS categorical-plot
    cds_index = gene_count.index
    cds_index = cds_index.to_numpy()
    cds_values = gene_count.values
    cds_indexlist = cds_index.tolist()
    cds_valuelist = cds_values.tolist()

    cds_valueslist = []
    for sublist in cds_valuelist:
        for item in sublist:
            cds_valueslist.append(item)

    cds_plot_df = pd.DataFrame(list(zip(cds_indexlist, cds_valueslist)), columns=[
                            'gene ID', 'N. of genes'])
    cds_plot_df = cds_plot_df.head(15)
    cds_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
    # plt.yticks("")
    plt.savefig(output_dir + '/' + basename + '_CDS_scatterplot.png')

    # nCDS categorical-plot
    ncds_index = ncds_count.index
    ncds_index = ncds_index.to_numpy()
    ncds_values = ncds_count.values
    ncds_indexlist = ncds_index.tolist()
    ncds_valuelist = ncds_values.tolist()

    ncds_valueslist = []
    for sublist in ncds_valuelist:
        for item in sublist:
            ncds_valueslist.append(item)

    ncds_plot_df = pd.DataFrame(list(zip(ncds_indexlist, ncds_valueslist)), columns=[
        'gene ID', 'N. of genes'])
    ncds_plot_df = ncds_plot_df.head(15)
    ncds_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
    plt.savefig(output_dir + '/' + basename + '_nCDS_scatterplot.png')

    # Intergenic categorical-plot
    intergenic_index = intergenic_count.index
    intergenic_index = intergenic_index.to_numpy()
    intergenic_values = intergenic_count.values
    intergenic_indexlist = intergenic_index.tolist()
    intergenic_valuelist = intergenic_values.tolist()

    intergenic_valueslist = []
    for sublist in intergenic_valuelist:
        for item in sublist:
            intergenic_valueslist.append(item)

    intergenic_plot_df = pd.DataFrame(list(zip(intergenic_indexlist, intergenic_valueslist)), columns=[
        'gene ID', 'N. of genes'])
    intergenic_plot_df = intergenic_plot_df.head(15)
    intergenic_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
    plt.savefig(output_dir + '/' + basename + '_tIG_scatterplot.png')

    # Upstream categorical-plot
    ups_index = upstream_count.index
    ups_index = ups_index.to_numpy()
    ups_values = upstream_count.values
    ups_indexlist = ups_index.tolist()
    ups_valuelist = ups_values.tolist()

    ups_valueslist = []
    for sublist in ups_valuelist:
        for item in sublist:
            ups_valueslist.append(item)

    ups_plot_df = pd.DataFrame(list(zip(ups_indexlist, ups_valueslist)), columns=[
        'gene ID', 'N. of genes'])
    ups_plot_df = ups_plot_df.head(15)
    ups_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
    plt.savefig(output_dir + '/' + basename + '_UPS_scatterplot.png')

    print("\n")
    print("Analysis completed")
    sys.exit()
#########################################################################################################################################################

# ROARY file parsing:
with open(roary_csv, 'r') as roary_file:
    roary_df = pd.read_csv(roary_file, error_bad_lines=False, index_col=False, dtype='unicode')
    tot_genes = roary_df.shape[0]
    # print(cds_seqID_HMN)

    # Hide the 'elementwise comparison failed' [WARNING]
    warnings.simplefilter(action='ignore', category=FutureWarning)
    cds_BLASTX = roary_df.loc[(roary_df == cds_seqID_HMN).any(1), 'Annotation']
    cds_BLASTX = cds_BLASTX.to_frame(name='Gene')
    cds_BLASTX = cds_BLASTX['Gene'].iloc[0]

    cds_core = roary_df.loc[(roary_df == cds_seqID_HMN).any(1), 'No. isolates']
    cds_core = cds_core.to_frame(name='No. isolates')
    cds_core = cds_core['No. isolates'].iloc[0]

    core_count_serie = roary_df['No. isolates'].max()
    if cds_core == core_count_serie:
        cds_core_value = 'Core'
    else:
        cds_core_value = 'Dispensable'

    ncds_BLASTX = roary_df.loc[(roary_df == ncds_seqID_HMN).any(1), 'Annotation']
    ncds_BLASTX = ncds_BLASTX.to_frame(name='Gene')
    ncds_BLASTX = ncds_BLASTX['Gene'].iloc[0]

    ncds_core = roary_df.loc[(roary_df == ncds_seqID_HMN).any(1), 'No. isolates']
    ncds_core = ncds_core.to_frame(name='No. isolates')
    ncds_core = ncds_core['No. isolates'].iloc[0]

    if ncds_core == core_count_serie:
        ncds_core_value = 'Core'
    else:
        ncds_core_value = 'Dispensable'

    # print(intergenic_seqID_HMN)
    inter_BLASTX = roary_df.loc[(roary_df == intergenic_seqID_HMN).any(1), 'Annotation']
    inter_BLASTX = inter_BLASTX.to_frame(name='Gene')
    inter_BLASTX = inter_BLASTX['Gene'].iloc[0]

    inter_core = roary_df.loc[(roary_df == intergenic_seqID_HMN).any(1), 'No. isolates']
    inter_core = inter_core.to_frame(name='No. isolates')
    inter_core = inter_core['No. isolates'].iloc[0]

    if inter_core == core_count_serie:
        inter_core_value = 'Core'
    else:
        inter_core_value = 'Dispensable'

    ups_BLASTX = roary_df.loc[(roary_df == ups_seqID_HMN).any(1), 'Annotation']
    ups_emptylist=[]

    ups_BLASTX = ups_BLASTX.to_frame(name='Gene')
    ups_BLASTX = ups_BLASTX['Gene'].iloc[0]

    ups_core = roary_df.loc[(roary_df == ups_seqID_HMN).any(1), 'No. isolates']
    ups_core = ups_core.to_frame(name='No. isolates')
    ups_core = ups_core['No. isolates'].iloc[0]

    if ups_core == core_count_serie:
        ups_core_value = 'Core'
    else:
        ups_core_value = 'Dispensable'


# - - - - - - - - - - - - RATIO - - - - - - - - - - - - #
#                  Doing some calculus..                #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - #
abs_gene_ratio = gene_number / tot_genes * 100
abs_ncds_ratio = ncds_number / tot_genes * 100
abs_intergenic_ratio = intergenic_number / tot_genes * 100
abs_upstream_ratio = upstream_number / tot_genes * 100

gene_ratio = round(abs_gene_ratio, 2)
ncds_ratio = round(abs_ncds_ratio, 2)
intergenic_ratio = round(abs_intergenic_ratio, 2)
upstream_ratio = round(abs_upstream_ratio, 2)


# - - - - - - - - - - - - LOG - - - - - - - - - - - - #
#       Logging some info to standard output..        #
# - - - - - - - - - - - - - - - - - - - - - - - - - - #
print(f"Total number of FOUND genes: {tot_genes}")
print("Parsing -cds file..")
print("\t", f"Total number of METHYLATED genes: {gene_number} ({gene_ratio} %)")
print("\t", f"Total number of METHYLATIONS found: {tot_meth_cds_df.iloc[0]['geneID_occurency']}")
print("\t", f"Max number of methylations found: {gene_HMN} on gene {cds_seqID_HMN}")
print("\t", f"Roary annotation: {cds_seqID_HMN} = {cds_BLASTX} -> {cds_core_value} genome")
print("\n")
print("Parsing -ncds file..")
print("\t", f"Total number of METHYLATED genes: {ncds_number} ({ncds_ratio} %)")
print("\t", f"Total number of METHYLATIONS found: {tot_meth_ncds_df.iloc[0]['geneID_occurency']}")
print("\t", f"Max number of methylations found: {ncds_HMN} on gene {ncds_seqID_HMN}")
print("\t", f"Roary annotation: {ncds_seqID_HMN} = {ncds_BLASTX} -> {ncds_core_value} genome")
print("\n")
print("Parsing -inter file..")
print("\t", f"Total number of METHYLATED regions: {intergenic_number} ({intergenic_ratio} %)")
print("\t", f"Total number of METHYLATIONS found: {tot_meth_intergenic_df.iloc[0]['geneID_occurency']}")
print("\t", f"Max number of methylations found: {intergenic_HMN} on gene {intergenic_seqID_HMN}")
print("\t", f"Roary annotation: {intergenic_seqID_HMN} = {inter_BLASTX} -> {inter_core_value} genome")
print("\n")
print("Parsing -ups file..")
print("\t", f"Total number of METHYLATED regions: {upstream_number} ({upstream_ratio} %)")
print("\t", f"Total number of METHYLATIONS found: {tot_meth_upstream_df.iloc[0]['geneID_occurency']}")
print("\t", f"Max number of methylations found: {upstream_HMN} up to gene {ups_seqID_HMN}")
print("\t", f"Roary annotation: {ups_seqID_HMN} = {ups_BLASTX} -> {ups_core_value} genome")
print("\n")



# - - - - - - - - - - - - GRAPH - - - - - - - - - - - - #
#                   From data to plots                  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# CDS categorical-plot
cds_index = gene_count.index  # pandas.core.indexes.base.Index
cds_index = cds_index.to_numpy() # numpy.ndarray (len 4368)
cds_values = gene_count.values  # numpy.ndarray (len 4368)

cds_indexlist = cds_index.tolist()
cds_valuelist = cds_values.tolist()

cds_valueslist = []
for sublist in cds_valuelist:
    for item in sublist:
        cds_valueslist.append(item)

cds_plot_df = pd.DataFrame(list(zip(cds_indexlist, cds_valueslist)), columns=['gene ID', 'N. of genes'])
cds_plot_df=cds_plot_df.head(15)
cds_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
# plt.yticks("")
plt.savefig(output_dir + '/' + basename + '_cds_scatterplot.png')

# nCDS categorical-plot
ncds_index = ncds_count.index
ncds_index = ncds_index.to_numpy()
ncds_values = ncds_count.values

ncds_indexlist = ncds_index.tolist()
ncds_valuelist = ncds_values.tolist()

ncds_valueslist = []
for sublist in ncds_valuelist:
    for item in sublist:
        ncds_valueslist.append(item)

ncds_plot_df = pd.DataFrame(list(zip(ncds_indexlist, ncds_valueslist)), columns=['gene ID', 'N. of genes'])
ncds_plot_df = ncds_plot_df.head(15)
ncds_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
# plt.yticks("")
plt.savefig(output_dir + '/' + basename + '_ncds_scatterplot.png')

# Intergenic categorical-plot
intergenic_index = intergenic_count.index
intergenic_index = intergenic_index.to_numpy()
intergenic_values = intergenic_count.values

intergenic_indexlist = intergenic_index.tolist()
intergenic_valuelist = intergenic_values.tolist()

intergenic_valueslist = []
for sublist in intergenic_valuelist:
    for item in sublist:
        intergenic_valueslist.append(item)

intergenic_plot_df = pd.DataFrame(list(zip(intergenic_indexlist, intergenic_valueslist)), columns=['gene ID', 'N. of genes'])
intergenic_plot_df = intergenic_plot_df.head(15)
intergenic_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
# plt.yticks("")
plt.savefig(output_dir + '/' + basename + '_intergenic_scatterplot.png')

# Upstream categorical-plot
ups_index = upstream_count.index
ups_index = ups_index.to_numpy()
ups_values = upstream_count.values

ups_indexlist = ups_index.tolist()
ups_valuelist = ups_values.tolist()

ups_valueslist = []
for sublist in ups_valuelist:
    for item in sublist:
        ups_valueslist.append(item)

ups_plot_df = pd.DataFrame(list(zip(ups_indexlist, ups_valueslist)), columns=['gene ID', 'N. of genes'])
ups_plot_df = ups_plot_df.head(15)
ups_plot_df.plot(x='N. of genes', y='gene ID', kind='scatter')
# plt.yticks("")
plt.savefig(output_dir + '/' + basename + '_upstream_scatterplot.png')


# - - - - - - - - - - - - WRITE - - - - - - - - - - - - #
#                   Write files to the                  #
#                    "outdir" folder                    # 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# LOGfile.log
logfile_path = output_dir + '/' + basename + '.txt'
logfile = open(logfile_path, 'w')
print(f"Total number of FOUND genes: {tot_genes}", file=logfile)
print("Parsing -cds file..", file=logfile)
print("\t", f"Total number of METHYLATED genes: {gene_number} ({gene_ratio} %)", file=logfile)
print("\t", f"Total number of METHYLATIONS found: {tot_meth_cds_df.iloc[0]['geneID_occurency']}", file=logfile)
print("\t", f"Max number of methylations found: {gene_HMN} on gene {cds_seqID_HMN}", file=logfile)
print("\t", f"Roary annotation: {cds_seqID_HMN} = {cds_BLASTX} -> {cds_core_value} genome", file=logfile)
print("\n", file=logfile)
print("Parsing -ncds file..", file=logfile)
print("\t", f"Total number of METHYLATED genes: {ncds_number} ({ncds_ratio} %)", file=logfile)
print("\t", f"Total number of METHYLATIONS found: {tot_meth_ncds_df.iloc[0]['geneID_occurency']}", file=logfile)
print("\t", f"Max number of methylations found: {ncds_HMN} on gene {ncds_seqID_HMN}", file=logfile)
print("\t", f"Roary annotation: {ncds_seqID_HMN} = {ncds_BLASTX} -> {ncds_core_value} genome", file=logfile)
print("\n", file=logfile)
print("Parsing -inter file..", file=logfile)
print("\t", f"Total number of METHYLATED regions: {intergenic_number} ({intergenic_ratio} %)", file=logfile)
print("\t", f"Total number of METHYLATIONS found: {tot_meth_intergenic_df.iloc[0]['geneID_occurency']}", file=logfile)
print("\t", f"Max number of methylations found: {intergenic_HMN} on gene {intergenic_seqID_HMN}", file=logfile)
print("\t", f"Roary annotation: {intergenic_seqID_HMN} = {inter_BLASTX} -> {inter_core_value} genome", file=logfile)
print("\n", file=logfile)
print("Parsing -ups file..", file=logfile)
print("\t", f"Total number of METHYLATED regions: {upstream_number} ({upstream_ratio} %)", file=logfile)
print("\t", f"Total number of METHYLATIONS found: {tot_meth_upstream_df.iloc[0]['geneID_occurency']}", file=logfile)
print("\t", f"Max number of methylations found: {upstream_HMN} up to gene {ups_seqID_HMN}", file=logfile)
print("\t", f"Roary annotation: {ups_seqID_HMN} = {ups_BLASTX} -> {ups_core_value} genome", file=logfile)
print("\n", file=logfile)
logfile.close()


# .num FILES

# cds_classpath = output_dir + '/' + basename + '_cds.num.txt'
# cds_occurrency.to_csv(path_or_buf=cds_classpath, sep="\t")
# ncds_classpath = output_dir + '/' + basename + '_ncds.num.txt'
# ncds_occurrency.to_csv(path_or_buf=ncds_classpath, sep="\t")
# intergenic_classpath = output_dir + '/' + basename + '_intergenic.num.txt'
# intergenic_occurrency.to_csv(path_or_buf=intergenic_classpath, sep="\t")
# ups_classpath = output_dir + '/' + basename + '_ups.num.txt'
# ups_occurrency.to_csv(path_or_buf=ups_classpath, sep="\t")

# - - - - - - - - - - READY-TO-USE-GFF - - - - - - - - - - -  - #
# Write in the output directory a per-feature-GFF that is       #
# more user friendly and ready to be used by the Circos RScript #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if evo_gff == True:

    # CDS
    coding_EVO_df = cds_df[['A','D','E','I']].copy()
    coding_EVO_df = coding_EVO_df.drop_duplicates(subset='I').assign(J=coding_EVO_df.groupby('I')['I'].transform('count'))
    coding_EVO_df.columns=['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    coding_EVO_df_dropped = coding_EVO_df.reset_index(drop=True)

    # print(coding_EVO_df_dropped)
    # print(coding_EVO_df_dropped['ID'])

    # evogenes_list = []
    anno_list = []

    for element in coding_EVO_df_dropped['name']:
        if not roary_df.loc[(roary_df == element).any(1), 'Annotation'].empty:
        # print(element)
        # evogenes_list.append(element)
            annocol = roary_df.loc[(roary_df == element).any(1), 'Annotation']
            annocol = annocol.to_frame(name='Gene')
            annocol = annocol['Gene'].iloc[0] #IndexError: single positional indexer is out-of-bounds -> annocol has an ID that is not present in roary_csv
        # print(annocol) # in fact, if you print annocol, it will show few IDs and then will stop when it finds the non-matched ID .. 
            anno_list.append(annocol)
        else:
            pass



    anno_df = pd.DataFrame(anno_list)
    anno_df.columns=['protein']
    anno_df_dropped = anno_df.reset_index(drop=True)

    coding_EVO_df_dropped = coding_EVO_df_dropped.join(anno_df_dropped)
    # print(coding_EVO_df_dropped)
    coding_EVO_path = output_dir + '/' + "evo_CDS.bed"
    coding_EVO_df_dropped.to_csv(coding_EVO_path, sep="\t", index = False)

    # non_CDS
    non_coding_EVO_df = non_cds_df[['A','D','E','I']].copy()
    non_coding_EVO_df = non_coding_EVO_df.drop_duplicates(subset='I').assign(J=non_coding_EVO_df.groupby('I')['I'].transform('count'))
    non_coding_EVO_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    non_coding_EVO_df_dropped = non_coding_EVO_df.reset_index(drop=True)

    # non_evogenes_list = []
    non_anno_list = []

    for nonelement in non_coding_EVO_df_dropped['name']:
        if not roary_df.loc[(roary_df == nonelement).any(1), 'Annotation'].empty:
            non_annocol = roary_df.loc[(roary_df == nonelement).any(1), 'Annotation']
            non_annocol = non_annocol.to_frame(name='Gene')
            non_annocol = non_annocol['Gene'].iloc[0]
            non_anno_list.append(non_annocol)
        else:
            pass

    non_anno_df = pd.DataFrame(non_anno_list)
    non_anno_df.columns = ['Genes']
    non_anno_df_dropped = non_anno_df.reset_index(drop=True)

    non_coding_EVO_df_dropped = non_coding_EVO_df_dropped.join(non_anno_df_dropped)
    non_coding_EVO_path = output_dir + '/' + "evo_nCDS.bed"
    non_coding_EVO_df_dropped.to_csv(non_coding_EVO_path, sep="\t", index = False)

    # INTERGENIC

    inter_EVO_df = intergenic_df[['A', 'D', 'E', 'I']].copy()
    inter_EVO_df = inter_EVO_df.drop_duplicates(subset='I').assign(J=inter_EVO_df.groupby('I')['I'].transform('count'))
    inter_EVO_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    inter_EVO_df_dropped = inter_EVO_df.reset_index(drop=True)

    inter_anno_list = []

    for interelement in inter_EVO_df_dropped['name']:
        if not roary_df.loc[(roary_df == interelement).any(1), 'Annotation'].empty:

            inter_annocol = roary_df.loc[(roary_df == interelement).any(1), 'Annotation']
            inter_annocol = inter_annocol.to_frame(name='Gene')
            inter_annocol = inter_annocol['Gene'].iloc[0]
            inter_anno_list.append(inter_annocol)
        else:
            pass

    inter_anno_df = pd.DataFrame(inter_anno_list)
    inter_anno_df.columns = ['Genes']
    inter_anno_df_dropped = inter_anno_df.reset_index(drop=True)

    inter_EVO_df_dropped = inter_EVO_df_dropped.join(inter_anno_df_dropped)
    inter_EVO_path = output_dir + '/' + "evo_intergenic.bed"
    inter_EVO_df_dropped.to_csv(inter_EVO_path, sep="\t", index = False)

    # UPS

    ups_EVO_df = upstream_df[['A','D', 'E', 'I']].copy()
    ups_EVO_df = ups_EVO_df.drop_duplicates(subset='I').assign(J=ups_EVO_df.groupby('I')['I'].transform('count'))
    ups_EVO_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score']
    ups_EVO_df_dropped = ups_EVO_df.reset_index(drop=True)

    ups_anno_list = []

    for elementup in ups_EVO_df_dropped['name']:
        if not roary_df.loc[(roary_df == elementup).any(1), 'Annotation'].empty:
            ups_annocol = roary_df.loc[(roary_df == elementup).any(1), 'Annotation']
            ups_annocol = ups_annocol.to_frame(name='Gene')
            ups_annocol = ups_annocol['Gene'].iloc[0]
            ups_anno_list.append(ups_annocol)
        else:
            pass

    ups_anno_df = pd.DataFrame(ups_anno_list)
    ups_anno_df.columns = ['Genes']
    ups_anno_df_dropped = ups_anno_df.reset_index(drop=True)

    ups_EVO_df_dropped = ups_EVO_df_dropped.join(ups_anno_df_dropped)
    ups_EVO_path = output_dir + '/' + "evo_upstream.bed"
    ups_EVO_df_dropped.to_csv(ups_EVO_path, sep="\t", index = False)

    # subprocess.call(f"Rscript --vanilla {input_dir}", shell=True)
    # subprocess.call(["Rscript", "--vanilla", input_dir])

    print("\n")
    print("Analysis completed")



# - - - - - - - - - - CHROMOSPLIT - - - - - - - - - - - #
# Write in the output directory chromosome-level sorted #
# GFFs that will be used as input for ms_chromoStat.py  #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - #

if split_flag == True:
    
    print("Ok, now splitting your input files to CHROMOSOME level")
    contigs_list_output = output_dir + '/' + "contigs.txt"
    contig_list = f"cut -f 1 {coding_ms} | uniq > {contigs_list_output}"
    os.system(contig_list)
    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output, error = process.communicate()

    print("Contigs list saved.")
    print("Converting your GFFs ..")

    with open(contigs_list_output, 'r') as ctgs:
        contigs = ctgs.read().splitlines()
        for line in contigs:
            contig_out = output_dir + '/' + f"{line}"
            with open(contig_out.format(line), "w") as f:
                commandBash_1 = f"grep '{line}' {coding_ms} > {contig_out}"
                os.system(commandBash_1)

# if show_error == True:

#     print(f"This is the error in CDS: {cds_seqID_HMN}")
#     print(f"This is the error in nCDS: {ncds_seqID_HMN}")
#     print(f"This is the error in UP: {intergenic_seqID_HMN}")
#     print(f"This is the error in tIG: {cds_seqID_HMN}")

    
print("\n")
print("Done!")



#                                                           *(((//(.
#               /////.           ,,***//////**,,.            *((//////,
#               //(////,/*/*,,**/*,,,,,***/*,,,,,,,****(*     /((///**//
#                /((***/,,,,**/*,,,,,,*/*,,,,,,***/,,,,,,**(/.(////////(
#                /*,*/*,,****(**,,****/,,,,,,****,,,,,,,**/*,,,*(/////(
#               ,***/*,*//*,,*/**/(**/***,****/*,,,,,**/*,,,*,**//////.
#               **,***/,,,,,**(*,,,,,,*/*//,,,*/**,**//*,,,**/**,*///
#               /,,*/,**,,,*/*,,,,,,,**(*,,,,,,,(*/,,,***,*/*,*%%%
#             ,((*,,*****,*(,,,,,,,,,*/*,,,,,,**/*,,,,,****/*,%%%%%%(
#      ,//((((//(((****/**(**,,,,,,,*/*,,,,,,**/,,,,,,,****///%%%%%########(###
#   .////*//*///////(//*,*/**,,*,,**(**,,,,,,*(*,,,,,,*/**/,*#(##/*/*/###(#((#/((
#   /(//////////(///(   /**/,*//***/,,,,,,,,,***,,,,*,**/,,((//  *,,,,,***/#%#///.
#                          ,/****,,**/*/*****/******/*/,*(/***        .******.
#                                 (/*/,,,/*,***,,,,*((((((.
#                                                .((((((((*
#                                               (((((/(///(
#                                             *((((((/(///,
#                                             ((GREED//**/
#                                            (((/(////(*
#                                          (((////*/(*
#                                        /((///(*(.
#                                     .(((((//,
#                                    .((#.
