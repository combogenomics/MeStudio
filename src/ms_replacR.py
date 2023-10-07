#!/usr/bin/env python3
# -*- coding: utf-8 -*

## - Libaries Import - ##
from curses.ascii import isspace
import os
from re import L
import sys
import random
import argparse
from pathlib import Path
import pandas as pd
import gffpandas.gffpandas as gffpd

## - Functions Definition - ##

def idQV(gff3):
    ## Init the GFF3 file
    path_gff = os.path.abspath(gff3)
    with open(path_gff, 'r') as pgff:
        df_gff = pd.read_table(pgff, header=None, sep= "\t")
        # for value in df_gff["iQV"]:
        #     if value == 'nan':
        #         print(value)

## Take 2 GFF3 files and check if the SeqID is in the same order for both
def gff_organizR(file_path, file_a, file_b):
    p = Path(file_path)

    with open (file_a) as f:
        for line in f:
            if not line.startswith('#') or line.startswith('##FASTA'):
                next(f)
                idx = {line.split()[0]: i for i, line in enumerate(f)}
    for fn in p.glob(file_b):
        with open(fn, 'r+') as f:
            dat = f.readlines()
            f.seek(0)
            for line in [dat[0]]+sorted(dat[1:],key=lambda l: idx.get(l.split()[0], 0)):
                f.write(line)

    print("\n","[OK] your files are now with the same order, proceed.")

## Create a tmp string with random numbers
def temporary():
    while True:
        my_string = 'tmp' + str(random.randint(1,1000000000))
        if not os.path.exists(my_string):
            return my_string

## Based on the tmp numeric string of temporary(), create tmp files (FASTA / GFF3)
def tempfile_generator(infile):
    random_id = temporary()
    path_infile = os.path.abspath(infile)

    if path_infile.endswith('.fasta') or path_infile.endswith('.fna'):
        fasta_temporary = open(path_out + '/' + 'genomic_' + random_id + '.fasta', 'w')
        with open(path_fasta, 'r+') as f:
            for line in f:
                line = line.replace('|', '_')
                fasta_temporary.write(line)
        fasta_temporary.close()

    if path_infile.endswith('_anno.gff') or path_infile.endswith('_annotation.gff'):
        anno_temporary = open(path_out + '/' + 'anno_'+random_id + '.gff', 'w')
        with open (path_anno, 'r+') as a:
            for line in a:
                if '\t' in line and not '#' in line:
                    line = line.replace("|", "_")
                    anno_temporary.write(line)
        anno_temporary.close()
    
    if path_infile.endswith('_smart.gff'):
        smart_temporary = open(path_out + '/' + 'smart_'+random_id + '.gff', 'w')
        with open (path_smart, 'r+') as s:
            for line in s:
                if line.startswith('#'):
                    pass
                elif line.startswith('>'):
                    pass
                else:
                    line = line.replace("|", "_")
                    smart_temporary.write(line)
        smart_temporary.close()

## Print to stdout a warning message if pipe symbol is found in a list of SeqIDs 
def warning_message(ilist):
    print("\n")
    for element in ilist:
        for letter in element:
            if letter == "|":
                print("\t", f"[WARNING]: illegal character detected in {element}. Pipe symbol will be modified.")

## Print to stdout name and number of SeqIDs found in the input files (FASTA / GFF3)
def gff_logger(infile):
    path_infile = os.path.abspath(infile)

    if path_infile.endswith('smart.gff'):
        with open(path_infile, 'r') as g:
            genomic_seqids = []
            first = ""
            current = ""
            for line in g:
                if not line.startswith('#') and '\t' in line:
                    current = line.split('\t')[0]
                    if current != first:
                        genomic_seqids.append(current)
                        first = current
                    else:
                        first = current
        contigs_occurrence = len(genomic_seqids)
        print("Parsing the METHYLATION [-Me] file")
        print(f"There are {contigs_occurrence} found headers:")
        for element in genomic_seqids:
            print("\t", element)
        warning_message(genomic_seqids)

    if path_infile.endswith('anno.gff'):
        with open(path_infile, 'r') as g:
            genomic_seqids = []
            first = ""
            current = ""
            for line in g:
                if not line.startswith('#') and '\t' in line:
                    current = line.split('\t')[0]
                    if current != first:
                        genomic_seqids.append(current)
                        first = current
                    else:
                        first = current
        contigs_occurrence = len(genomic_seqids)
        print("Parsing the ANNOTATION [-anno] file")
        print(f"There are {contigs_occurrence} found headers:")
        for element in genomic_seqids:
            print("\t", element)
        warning_message(genomic_seqids)
    
    if path_infile.endswith('.fasta') or path_infile.endswith('.fna'):
        with open(path_fasta, 'r') as f: 
            data_fasta = f.readlines()
            fasta_headers = []
            for line in data_fasta:
                if line.startswith('>'):
                    line = line.split(' ')[0]
                    line = line.replace('>', '')
                    line=line.rstrip()
                    fasta_headers.append(line)

        contigs_number = len(fasta_headers)
        print("Parsing the FASTA [-f] file")
        print(f"There are {contigs_number} found headers:")
        for element in fasta_headers:
            print("\t", element)
        warning_message(fasta_headers)

## - Main Function Definition - ##
def ms_replacR() -> tuple:
    ap = argparse.ArgumentParser()
    ap.add_argument("-out", "--outputdir", help="path to new files directory")
    ap.add_argument("-anno", "--annotation", help="path to file produced by genomic annotator [GFF-file]")
    ap.add_argument("-f", "--fasta", help="path to genome file [FASTA/FNA-file]")
    ap.add_argument("-smart", "--smartlink", help="path to file produced by the sequencer [GFF-file]")
    ap.add_argument("-iqv", "--identificationQV", help="If used, this flag clears the methylations with no identificationQV",action="store_true", default=False)
    return ap.parse_args()

args = ms_replacR()

# Defining file-paths:
output_dir = args.outputdir
genomic_file = args.annotation
fasta_file = args.fasta
methylation_file = args.smartlink
identification_qvalue = args.identificationQV
## - - - - - - - - - - - - - - - - - - - - - ##
path_out = os.path.abspath(output_dir)
path_anno = os.path.abspath(genomic_file)
path_fasta = os.path.abspath(fasta_file)
path_smart = os.path.abspath(methylation_file)
## - - - - - - - - - - - - - - - - - - - - - ##
name_out = os.path.basename(path_out)
name_anno = os.path.basename(path_anno)
name_fasta = os.path.basename(path_fasta)
name_smart = os.path.basename(path_smart)
# name_out = os.path.basename(path_out.split('.')[0])
# name_anno = os.path.basename(path_anno.split('.')[0])
# name_fasta = os.path.basename(path_fasta.split('.')[0])
# name_smart = os.path.basename(path_smart.split('.')[0])

## - Files extensions check - ##

# Checking annotation file
# if not genomic_file.endswith('_anno.gff') or not genomic_file.endswith('_annotation.gff'):
#    sys.exit("[WARNING] The annotation file is not correctly named. Please change it by adding '_anno.gff' or '_annotation.gff' as extension, or make sure you added the correct file.")
# else:
#     print("[OK] The annotation file (-anno, --annotation) is correctly named.")
# # Checking fasta file
# if not fasta_file.endswith('fasta') or not fasta_file.endswith('fna'):
#    sys.exit("[WARNING] The genomic file is not correctly named. Please change it by adding '.fasta' or '.fna' as extension, or make sure you added the correct file.")
# else:
#     print("[OK] The genomic file (-f, --fasta) is correctly named.")
# # Checking smartlink file
# if not methylation_file.endswith('_smart.gff'):
#    sys.exit("[WARNING] The SmartLink file is not correctly named. Please change it by adding '_smart.gff' as extension, or make sure you added the correct file.")
# else:
#     print("[OK] The SmartLink file (-smart, --smartlink) is correctly named.")

# Create -out's folder and log some info
if not os.path.exists(path_out):
    os.makedirs(path_out)
    basename = os.path.basename(path_out)
    print("\n")
    print(f"Creating {basename} folder")
else:
    sys.exit("[WARNING] The chosen directory already exists, please change folder's name")

## - STDOUT printing information - ##

# INFO: FASTA file:
print("\n")
gff_logger(path_fasta)
#INFO: ANNOTATION file
print("\n")
gff_logger(path_anno)
# INFO: METHYLATION file
print("\n")
gff_logger(path_smart)

## - Temporary files creation - ##
print("Writing and saving temporary files")
tempfile_generator(path_fasta)
tempfile_generator(path_anno)
tempfile_generator(path_smart)

## Initialize temporary objects
temporary_genomic = ''
temporary_annotation = ''
temporary_methylation = ''

for element in os.listdir(path_out):
    if element.startswith("genomic"):
        temporary_genomic = element
    if element.startswith("anno"):
        temporary_annotation = element
    if element.startswith("smart"):
        temporary_methylation = element
## Assign a complete path to each temporary file
temporary_genomic = path_out + '/' + temporary_genomic
temporary_annotation = path_out + '/' + temporary_annotation
temporary_methylation = path_out + '/' + temporary_methylation
print("[OK]")
## - Temporary files ordering - ##
print("Indexing and ordering your temporary files")
tm = open(temporary_methylation)
for line in tm:
    if not line.startswith('#') or line.startswith('##FASTA'):
        next(tm)
        idx = {line.split()[0]: i for i, line in enumerate(tm)}

    with open(temporary_annotation, 'r+') as ta:
        dat = ta.readlines()
        ta.seek(0)
        for line in sorted(dat[0:],key=lambda l: idx.get(l.split()[0], 0)):
            ta.write(line)
tm.close() 
print("[OK]")

# if identification_qvalue == True:
#     idQV(temporary_methylation)

## Rename the temporary files as the originals
print("Renaming and saving your files to the output directory")
os.rename(temporary_genomic, path_out + '/' + name_fasta)
os.rename(temporary_annotation, path_out + '/' + name_anno)
os.rename(temporary_methylation, path_out + '/' + name_smart)
print("[OK]")
print("Pre-filtering has been successfully termined.")
print("Thank you for using MeStudio ReplacR, enjoy the rest of the pipeline.")
print("\n")



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
