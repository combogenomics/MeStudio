## - Libraries - ##
import os
import sys
import argparse
import subprocess
import numpy as np
import pandas as pd
from operator import ge
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.core import indexing
from pandas.core.frame import DataFrame
from pandas.core.indexes.base import Index

## - Functions - ##


def bed_builder(msgff3, annogff3):
    # Init the files
    path_anno = os.path.abspath(annogff3)
    path_ms = os.path.abspath(msgff3)
    # Get the feature (CDS, nCDS, tIG, US)
    motif_name = os.path.basename(msgff3.split("_")[0])
    feature_name = os.path.basename(msgff3.split("_")[1])
    # Define the dataframes
    anno_df = pd.read_table(path_anno, header=None)
    anno_df.columns = ["seqid", "source", "type", "start",
                       "end", "score", "strand", "phase", "attribute"]
    # - - - - - - - - - - - - - - #
    ms_df = pd.read_table(path_ms, header=None)
    ms_df.columns = ["seqid", "source", "type", "start",
                     "end", "score", "strand", "phase", "attribute"]
    ms_df = ms_df[['seqid', 'start', 'end', 'attribute']]
    # Define 2 new dataframes: the gene-id and the gene-annotation
    attribute_df = anno_df["attribute"].str.split("product=", expand=True)
    gene_id = attribute_df[0].str.split(";", expand=True)[0]
    gene_id = gene_id.str.split("ID=", expand=True)[1]
    gene_id.columns = ["gene-id"]  # gene-id dataframe
    gene_anno = attribute_df[1]  # gene-anno dataframe
    # Define a raw dataframe
    raw_df = anno_df[['seqid', 'start', 'end']].copy()
    # Define the BED dataframe as merge of 'raw' and 'gene-id' and 'gene-anno'
    bed_df = pd.concat([raw_df, gene_id, gene_anno], axis=1)
    bed_df.columns = ['seqid', 'start', 'end',
                      'gene-id', 'gene-anno']  # rename 'em all
    # Add a 5th column 'score' which takes the number of methylations x gene (attribute)
    ms_df = ms_df.drop_duplicates(subset='attribute').assign(
        score=ms_df.groupby('attribute')['attribute'].transform('count'))
    ms_df = ms_df.reset_index(drop=True)

    match_list = []  # initi a list for storing the results
    for gene in ms_df["attribute"]:
        if not bed_df.loc[(bed_df == gene).any(1), 'gene-anno'].empty:
            match = str(bed_df.loc[(bed_df == gene).any(
                1), 'gene-anno'].reset_index(drop=True).to_list()[0])
            match_list.append(match)
        else:
            match = "Not found"
            match_list.append(match)
    # Add 'gene-annotation' column to BED file
    ms_df["gene-annotation"] = match_list
    ms_df = ms_df.sort_values(by=['score'], ascending=False)
    # Save the result
    if os.path.basename(path_ms).endswith("_CDS.gff"):
        path_outfile = output_dir + "/" + f"{motif_name}_CDS.bed"
        ms_df.to_csv(path_outfile, index=False, sep="\t")

    if os.path.basename(path_ms).endswith("_nCDS.gff"):
        path_outfile = output_dir + "/" + f"{motif_name}_nCDS.bed"
        ms_df.to_csv(path_outfile, index=False, sep="\t")

    if os.path.basename(path_ms).endswith("_tIG.gff") or os.path.basename(path_ms).endswith("_true_intergenic.gff"):
        path_outfile = output_dir + "/" + f"{motif_name}_tIG.bed"
        ms_df.to_csv(path_outfile, index=False, sep="\t")

    if os.path.basename(path_ms).endswith("_US.gff") or os.path.basename(path_ms).endswith("_upstream.gff"):
        path_outfile = output_dir + "/" + f"{motif_name}_US.bed"
        ms_df.to_csv(path_outfile, index=False, sep="\t")


def roary_reader(rr_csv):
    path_infile = os.path.abspath(rr_csv)
    roary_df = pd.read_csv(path_infile, index_col=False,
                           dtype='unicode', on_bad_lines='skip', sep=",")

    # get the corresponding value for "core" genes
    core_value = roary_df["No. isolates"].max()
    core_df = roary_df[roary_df["No. isolates"].str.contains(
        core_value)]  # get the core genes dataframe
    cloud_df = roary_df[~roary_df["No. isolates"].str.contains(
        core_value)]  # get the dispensable genes dataframe

    # Get some occurrences
    n_tot = roary_df.shape[0]  # total number of genes
    n_core = core_df.shape[0]  # number of core genes
    n_cloud = cloud_df.shape[0]  # number of dispensable genes
    perc_core = (n_core / n_tot)*100  # core genes percentage over totals
    perc_core = round(perc_core, 1)  # and round it
    perc_cloud = (n_cloud / n_tot)*100  # cloud genes percentage over totals
    perc_cloud = round(perc_cloud, 1)  # and round it

    # Print some info to the STDOUT
    print("Parsing the pangenome")
    print(f"{n_tot} found genes")
    print(f"Number of core genes: {n_core} ({perc_core}%)")
    print(f"Number of dispensable genes: {n_cloud} ({perc_cloud}%)")


def roary_miner(rr_csv, id):
    path_infile = os.path.abspath(rr_csv)
    roary_df = pd.read_csv(path_infile, index_col=False,
                           dtype='unicode', on_bad_lines='skip', sep=",")
    # get the corresponding value for "core" genes
    core_value = roary_df["No. isolates"].max()
    # Find the annotation of the gene (without pandas indexes) and return it as string
    match = str(roary_df.loc[(roary_df == id).any(
        1), 'Annotation'].reset_index(drop=True).to_list()[0])
    core_match = str(roary_df.loc[(roary_df == id).any(
        1), 'No. isolates'].reset_index(drop=True).to_list()[0])
    if core_value == core_match:
        res = " (Core)"
        return match + res
    else:
        res = " (Dispensable)"
        return match + res


def stdout_logger(motif, feature, ratio, hmg, hmo, n_motifs, n_memotifs):
    print("\n")
    print(f"Parsing the {motif} on {feature}")
    print(f"Number of {motif}s found: {n_motifs}")
    print(f"Number of methylated {motif}s found: {n_memotifs}")
    print(f"Ratio: {ratio}")
    print(f"Highest Methylated Gene: {hmg} with {hmo} methylated {motif}s")


def ms_gff(gff3):
    path_infile = os.path.abspath(gff3)
    name_motif = os.path.basename(gff3.split("_")[0])
    name_feature = os.path.basename(gff3.split("_")[1])
    with open(path_infile, 'r') as g:
        # Init the dataframe
        gff_df = pd.read_table(g, header=None)
        gff_df.columns = ["seqid", "source", "type", "start",
                          "end", "score", "strand", "phase", "attribute"]
        # Operate to the dataframe:
        # Filter out '0' values from 'score'
        gff_cdf = gff_df[gff_df['score'] != 0]
        # Get some numbers:
        n_motifs = gff_df.shape[0]  # total number of found motifs
        # total number of methylated found motifs
        n_memotifs = gff_cdf.shape[0]
        # Check -s flag:
        if strict == True:
            # Occurrence check: if no methylations have been found, stop the analysis.
            if n_memotifs == 0:
                sys.exit(
                    f"Sorry, the analysis has been stopped: you have not enough methylations for {name_feature} ")

        me_ratio = n_memotifs / n_motifs  # fraction of methylated motifs over the total
        me_ratio = round(me_ratio, 4)  # round it at the 4th decimal
        occurrences_count = gff_cdf['attribute'].value_counts().to_frame(
            name='occurrences')  # dataframe with N* methylations x attribute
        # H-ighest M-ethylation O-ccurence
        hmo = occurrences_count['occurrences'].max()
        hmg = gff_cdf['attribute'].value_counts().iloc[[0]].keys()[
            0]  # H-ighest M-ethylated G-ene

        # Printing some info to the STDOUT
        if os.path.basename(path_infile).endswith("_CDS.gff"):
            stdout_logger(name_motif, "CDS", me_ratio,
                          hmg, hmo, n_motifs, n_memotifs)

        if os.path.basename(path_infile).endswith("_nCDS.gff"):
            stdout_logger(name_motif, "nCDS", me_ratio,
                          hmg, hmo, n_motifs, n_memotifs)

        if os.path.basename(path_infile).endswith("_tIG.gff") or os.path.basename(path_infile).endswith("_intergenic.gff"):
            stdout_logger(name_motif, "tIG", me_ratio,
                          hmg, hmo, n_motifs, n_memotifs)

        if os.path.basename(path_infile).endswith("_US.gff") or os.path.basename(path_infile).endswith("_upstream.gff"):
            stdout_logger(name_motif, "US", me_ratio,
                          hmg, hmo, n_motifs, n_memotifs)


# Arg-parser:
def ms_analyzR() -> tuple:
    ap = argparse.ArgumentParser(
        prog="MeStudio Analyzer",
        description='The post processing tool',
        epilog='Version 2.0 (@greed)'
    )
    requiredNamed = ap.add_argument_group('mandatory arguments')
    requiredNamed.add_argument(
        "-o", "--outputdir", help="path to your output files")
    requiredNamed.add_argument("-CDS", "--coding", help="CDS.gff file")
    requiredNamed.add_argument("-nCDS", "--noncoding", help="nCDS.gff file")
    requiredNamed.add_argument("-tIG", "--intergenic", help="tIG.gff file")
    requiredNamed.add_argument("-US", "--upstream", help="US.gff file")
    requiredNamed.add_argument(
        "-anno", "--annotation", help="path to file produced by genomic annotator [GFF-file]")
    # Optional arguments
    ap.add_argument("-rpl", "--replicon", help="Rearrange your input GFFs for chromosomes",
                    action="store_true", default=False)
    ap.add_argument("-s", "--strict", help="If used, this flag stops the analysis when 0 methylations are found in the GFF3.",
                    action="store_true", default=False)
    ap.add_argument(
        "-r", "--roary", help="path to 'gene_presence_abscence.csv' file produced by Roary", action="store_true")
    ap.add_argument("-q", "--quiet", help="Running ms analyzR quietly (no stdout)",
                    action="store_true", default=False)
    return ap.parse_args()


args = ms_analyzR()

# Mandatory
output_dir = args.outputdir
cds_ms = args.coding
ncds_ms = args.noncoding
tig_ms = args.intergenic
us_ms = args.upstream
annotation_file = args.annotation

# Optionals
replicon = args.replicon
strict = args.strict
roary = args.roary
quiet = args.quiet

# Path definition
path_out = os.path.abspath(output_dir)
path_cds = os.path.abspath(cds_ms)
path_ncds = os.path.abspath(ncds_ms)
path_tig = os.path.abspath(tig_ms)
path_us = os.path.abspath(us_ms)
path_anno = os.path.abspath(annotation_file)
## - - - - - - - - - - - - - - - - - - - - - ##
name_anno = os.path.basename(path_anno)

# Output directory creation:
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    basename = os.path.basename(output_dir)
    print("\n")
    print(f"Creating the {basename} folder ..")
else:
    sys.exit(
        "[WARNING]: The chosen directory already exists, please change folder's name")

if roary == True:
    roary_reader(roary)

if quiet == True:
    bed_builder(cds_ms, annotation_file)
    bed_builder(ncds_ms, annotation_file)
    bed_builder(tig_ms, annotation_file)
    bed_builder(us_ms, annotation_file)
else:
    ms_gff(cds_ms)
    ms_gff(ncds_ms)
    ms_gff(tig_ms)
    ms_gff(us_ms)
    print("\n")
    bed_builder(cds_ms, annotation_file)
    bed_builder(ncds_ms, annotation_file)
    bed_builder(tig_ms, annotation_file)
    bed_builder(us_ms, annotation_file)

if replicon == True:
    print("Replicon level flag: ON")
    print("Splitting your files to replicon level")
    replicon_list_output = output_dir + '/' + "replicons.txt"
    replicon_command = f"cut -f 1 {cds_ms} | uniq > {replicon_list_output}"

    print("\t", "Saving the found replicons to 'replicons.txt ..'")
    os.system(replicon_command)
    print("\t", "\t", "[Done]")

    print("\t", "Now converting GFF3 from feature level to replicon level ..")
    with open(replicon_list_output, 'r') as ctgs:
        contigs = ctgs.read().splitlines()
        for line in contigs:
            contig_out = output_dir + '/' + f"{line}.gff"
            with open(contig_out.format(line), "w") as f:
                saving_command = f"grep '{line}' {cds_ms} > {contig_out}"
                os.system(saving_command)
    print("\t", "\t", "[Done]")
    print("GFF3 files split by replicon are now saved in the output directory.")

print("\n")
print("Post processing has been successfully termined.")
print("Thank you for using MeStudio AnalyzR.")
print("\n")
#
#
#
#
#
#
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
