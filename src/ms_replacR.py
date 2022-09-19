#!/usr/bin/env python3
# -*- coding: utf-8 -*

import os
import sys
import random
import shutil
import argparse
from pathlib import Path

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

def temp_generator():
    while True:
        my_string = 'tmp' + str(random.randint(1,1000000000))
        if not os.path.exists(my_string):
            return my_string

def warning_message(urlist):
    print("\n")
    for element in urlist:
        for letter in element:
            #if letter == in_char:
            if letter == "|":
                print("\t", f"[WARNING]: illegal character detected in {element}")



def ms_replacR() -> tuple:
    ap = argparse.ArgumentParser()
    #ap.add_argument("-in", "--inputdir",help="path to your files directory")
    ap.add_argument("-out", "--outputdir", help="path to new files directory")
    ap.add_argument("-g", "--genomic", help="path to file produced by genomic annotator [GFF-file]")
    ap.add_argument("-f", "--fasta", help="path to genome file [FASTA/FNA-file]")
    ap.add_argument("-Me", "--methylation", help="path to file produced by the sequencer [GFF-file]")
    #ap.add_argument("-i", "--input_word",help="Element to delete [SYMBOL/STRING]")
    #ap.add_argument("-o", "--output_word", help="Element to insert [SYMBOL/STRING]")
    return ap.parse_args()

args = ms_replacR()

# Defining file-objects:
# input_dir = args.inputdir
output_dir = args.outputdir
genomic_file = args.genomic
fasta_file = args.fasta
methylation_file = args.methylation
#in_char = args.input_word
#out_char = args.output_word

# - - - - - - - - - - - - - - DIR - - - - - - - - - - - - - -- #
# Create and populate the -out directory with given files      #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- #

# Create -out's folder and log some info
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    basename = os.path.basename(output_dir)
    print("\n")
    print(f"Creating the {basename} folder ..")
else:
    sys.exit("[WARNING]: The chosen directory already exists, please change folder's name")

# Copy files from -input_dir to -output_dir

# os.system('cp ' + input_dir + genomic_file + ' ' + output_dir)
# os.system('cp ' + input_dir + fasta_file + ' ' + output_dir)
# os.system('cp ' + input_dir + methylation_file + ' ' + output_dir)

# os.system('cp ' + genomic_file + ' ' + output_dir)
# os.system('cp ' + fasta_file + ' ' + output_dir)
# os.system('cp ' + methylation_file + ' ' + output_dir)

# - - - - - - - - - - - - - - LOG - - - - - - - - - - - - - - #
# Printing all replicon-info contained inside the input files #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # 

# INFO: FASTA file
print("\n")
print ("Parsing FASTA [-f] file ...", "\n")

# Parsing FASTA file:
with open(fasta_file, 'r') as f: 
    data_fasta = f.readlines()
    fasta_headers = []
    for line in data_fasta:
        if line.startswith('>'):
            line = line.split(' ')[0]
            line = line.replace('>', '')
            line=line.rstrip()
            fasta_headers.append(line)
f.close()


fasta_len = len(fasta_headers)
print(f"These are the {fasta_len} FASTA headers:", "\n")
for element in fasta_headers:
    print("\t", element)
warning_message(fasta_headers)

#INFO: GENOMIC file
print("\n")
print("Parsing GENOMIC [-g] file ...", "\n")

# Parsing GENOMIC file:
with open(genomic_file, 'r') as g:
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
g.close()

gen_len = len(genomic_seqids)
print(f"These are the {gen_len} GENOMIC headers:", "\n")
for element in genomic_seqids:
    print("\t", element)
warning_message(genomic_seqids)

# INFO: METHYLATION file
print("\n")
print("Parsing METHYLATION [-Me] file ...", "\n")

# # Parsing METHYLATION file:
with open(methylation_file, 'r') as m:
    methylation_seqids = []
    first = ""
    current = ""
    for line in m:
        if not line.startswith('#') and '\t' in line:
            current = line.split('\t')[0]
            if current != first:
                methylation_seqids.append(current)
                first = current
            else:
                first = current
m.close()

meth_len = len(methylation_seqids)
print(f"These are the {meth_len} METHYLATION headers:", "\n")
for element in methylation_seqids:
    print("\t", element)
warning_message(methylation_seqids)

print("\n")
#print(f"Replacing in: '{in_char}' with out: '{out_char}' ")
print(f"Replacing in: '|' with out: '_' ")
print("Doing some abracadabras ...")
print("\n")

# - - - - - - - - - - - - - TMP - - - - - - - - - - - - - - - #
# Create and modify TMP files that you will use as MAIN later #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Initializing randomic numerical index
randomID = temp_generator()

# Create fasta_tmp_<n_index>
fasta_temp_file = open(output_dir + '/' + 'fasta_'+randomID + '.fasta', 'w')
fasta_completepath = output_dir + '/' + fasta_file
fasta_name_path = os.path.basename(fasta_file)
with open(fasta_file, 'r+') as f:
    for line in f:
        line = line.replace("|", "_")
        #line = line.replace(in_char, out_char)
        fasta_temp_file.write(line)
fasta_temp_file.close()


# Create genomic_tmp_<n_index>
genomic_temp_file = open(output_dir + '/' + 'genomic_'+randomID + '.gff', 'w')
genomic_completepath = output_dir + '/' + genomic_file
with open (genomic_file, 'r+') as g:
    for line in g:
        if not line.startswith('#') and not line.startswith('##FASTA') and not line.startswith('>') and not \
                line.startswith('A') and not line.startswith('T') and not line.startswith('C') and not line.startswith('G'):
            line = line.replace("|", "_", 1)
            #line = line.replace(in_char, out_char, 1)
            genomic_temp_file.write(line)
genomic_temp_file.close()


# Create methylation_tmp_<n_index>
methylation_temp_file = open(output_dir + '/' + 'metylation_'+randomID + '.gff', 'w')
methylation_completepath = output_dir + '/' + methylation_file
with open(methylation_file, 'r+') as me:
    for line in me:
        if not line.startswith('#') and not line.startswith('##FASTA') and not line.startswith('>') and not \
                line.startswith('A') and not line.startswith('T') and not line.startswith('C') and not line.startswith('G'):
            #line = line.replace(in_char, out_char, 1)
            line = line.replace("|", "_", 1)
            methylation_temp_file.write(line)
methylation_temp_file.close()


# - - - - - - - - - - - - tmp INIT - - - - - - - - - - - - - - - #
# Initialize the tmp_files in order to read and check them later #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- #

tempfile_genomic = output_dir + '/' + 'genomic_'+randomID + '.gff'
tempfile_methylation = output_dir + '/' + 'metylation_'+randomID + '.gff'
tempfile_fasta = output_dir + '/' + 'fasta_' + randomID + '.fasta'

with open(tempfile_fasta, 'r') as f:
    data_fasta_tmp = f.readlines()
    fasta_headers_tmp = []
    for line in data_fasta_tmp:
        if line.startswith('>'):
            line = line.split(' ')[0]
            line = line.replace('>', '')
            line = line.rstrip()
            fasta_headers_tmp.append(line)
f.close()

with open(tempfile_genomic, 'r') as g:
    genomic_seqids_tmp = []
    first = ""
    current = ""
    for line in g:
        if not line.startswith('#') and '\t' in line:
            current = line.split('\t')[0]
            if current != first:
                genomic_seqids_tmp.append(current)
                first = current
            else:
                first = current
g.close()

with open(tempfile_methylation, 'r') as me:
    methylation_seqids_tmp = []
    first = ""
    current = ""
    for line in me:
        if not line.startswith('#'):
            current = line.split('\t')[0]
            if current != first:
                methylation_seqids_tmp.append(current)
                first = current
            else:
                first = current
me.close()

# - - - - - - - - - - - SECURITY - - - - - - - - - - - #
# Checking if your output files are valid..            #
# - - - - - - - - - - - - - - - - - - - - - - - - - -- #

if os.path.getsize(tempfile_fasta) == 0:
    sys.exit("[ERROR]: something went wront, your FASTA is not valid. Please remove your folder.")
if os.path.getsize(tempfile_genomic) == 0:
    sys.exit("[ERROR]: something went wront, your GENOMIC is not valid. Please remove your folder.")
if os.path.getsize(tempfile_methylation) == 0:
    sys.exit("[ERROR]: something went wront, your METHYLATION is not valid. Please remove your folder.")
# if fasta_headers_tmp[0] != genomic_seqids_tmp[0] or fasta_headers_tmp[0] != methylation_seqids_tmp[0]:
#     sys.exit("[ERROR]: seqids are not coincidents. FILES not valid. Please remove your folder.")

# print(fasta_headers_tmp)
# print(genomic_seqids_tmp)
# print(methylation_seqids_tmp)

id_check = fasta_headers_tmp == genomic_seqids_tmp == methylation_seqids_tmp

if id_check == False:
    sys.exit("[ERROR]: your input file headers are different.")
else:
    print("Your are ready in 3,2,1 ..")
# - - - - - - - - - - - - MODs - - - - - - - - - - - - #
# Use a certain file as footprint in order to set the  #
# GFF's core sequence and avoid tidy-array-disorders   #
# - - - - - - - - - - - - - - - - - - - - - - - - - -- #

p = Path(output_dir)
met_f = open(tempfile_methylation)

for line in met_f:
    if not line.startswith('#') or line.startswith('##FASTA'):
        next(met_f)
        idx = {line.split()[0]: i for i, line in enumerate(met_f)}

with open(tempfile_genomic, 'r+') as tg_f:
    dat = tg_f.readlines()
    tg_f.seek(0)
    for line in sorted(dat[0:],key=lambda l: idx.get(l.split()[0], 0)):
        tg_f.write(line)
met_f.close()

print("\n","[OK] your files are now with the same order, proceed.")


# Rename all the files
old_fasta_path = output_dir + '/' + fasta_file
old_genomic_path = output_dir + '/' + genomic_file
old_methylation_path = output_dir + '/' + methylation_file

os.rename(tempfile_fasta, fasta_file)
os.rename(tempfile_genomic, genomic_file)
os.rename(tempfile_methylation, methylation_file)

os.system('cp ' + genomic_file + ' ' + output_dir)
os.system('cp ' + fasta_file + ' ' + output_dir)
os.system('cp ' + methylation_file + ' ' + output_dir)

# print(methylation_seqids_tmp, genomic_seqids_tmp, fasta_headers_tmp)





# - - - - - - - - - - - - README - - - - - - - - - - - #
# Write a README.txt file with all the informations to #
# proceed the MeStudio processing pipeline.            #
# - - - - - - - - - - - - - - - - - - - - - - - - - -- #

# readmepath = output_dir + '/' + 'README.txt'
# readme_file = open(readmepath,'w')
# text = 'Hello, World!'
# readme_file.write(text)
# readme_file.close()



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
