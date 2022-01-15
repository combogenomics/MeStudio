# MeStudio
MeStudio is a tool which allows us to analyse and combine genome-wide methylation profiles with genomic features.
MeStudio connects several procedural software components that can be run either individually or as part of a pipeline.

*MeStudio version 1.0*

## Table of Contents

- [Quick start](#Quick-start)
- [MeStudio ReplacR](#MeStudio-ReplacR)
- [MeStudio Core](#MeStudio-Core)
- [MeStudio AnalyzR](#MeStudio-AnalyzR)
- [Usage](#Usage)
- [Results](#Results)
- [Reference](#Reference)


## Quick start

Hey sorcerer's apprentice, we know.. You got a lot to work so you may have no time to read all the sorcerer's book.
Please, feel free to directly click on "Usage" and check the *bard-way* to run MeStudio. All's gonna work in a wand flight. 

Inside the [toyset](/toyset/) folder you are going to find three files needed to start the analysis:

- FSMMA_genomic.fna
- FSMMA_genomic.gff
- FSMMA_methylation.gff

```FSMMA_genomic.fna``` is the the genome of FSMMA strain of *Sinorhizobium meliloti*. You can get more informations [here](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/382/?utm_source=None&utm_medium=referral&utm_campaign=KnownItemSensor:taxname).
```FSMMA_genomic.gff``` is the product of the genomic annotation (executed via [Prokka](https://github.com/tseemann/prokka) annotator). The file is reported in GFF3 format.
```FSMMA_methylation.gff``` is the GFF3 file with methylation positions obtained through the sequencer. We performed the analysis using the PacBio RT-SMRT sequencing.

N.B. please note that Prokka and Roary (see *Usage*) can be easily found to [Galaxy](https://usegalaxy.org), the installation of these tools is not strictly needed.


## MeStudio ReplacR
In order to properly run *MeStudio Core*, a pre-processing python-based script named *ms_replacR* has been implemented and is highly suggested to be used. You can find the source code here.
*ms_replacR* expects six arguments:
1. *output directory*, in which results and log files will be written
2. *genomic annotation*, in the GFF3 format
3. *methylation annotation*, a sequencer-produced modified base calls in the GFF3 format
4. *genomic sequence*, in fasta or fna file format
5. *input* string type field
6. *output* string type field


The last flags (.5 and .6) are used to communicate a certain character or string to replace.
As a matter of fact MeStudio requires consistent formatting as far as the sequence identifiers are concerned but sometimes the annotation process can lead to
genomic headers alterations for what concern “seqname” fields identity and special characters (e.g. a pipe symbol replaced by the underscore).

An example is reported here down below:

```FSMMA_genomic.fna```

> 000000F|arrow
> 
> GCCGGTCCAGCGCAAAACCCTCGCTCGGCGTGATCGAGAGTATGCGCTGCGAGCCGAGGT
> CGGGCCAGAAGAGCTTCGAATTCACGAGCCGGAAATGCGGTGCGACGATAACGCGTTCGA

```FSMMA_genomic.gff```


| 000000F_arrow | Prodigal:2.6 | CDS | 95 | 538 | . | + | 0 | ID=JPHAALHC_00001 |
|---------------|--------------|-----|----|-----|---|---|---|-------------------|

As you can see, the header of the fasta file has a pipe as delimeter while the *seqname* column of the GFF3 file as the underscore as delimeter.
This difference in the formatting syntax can create some troubles and here the last flags are crucial to fix the problem.

More over, depending on the annotator, we noticed that sometimes we can find different order of the contigs between fasta and GFF3 files.
*ms_replacR* also fix this kind of anomaly in order to avoid biased results.


###### How do I run ms_replacR?

Before running *ms_replacR*, please make sure to have python3.8 (or above) properly installed on your computer. [Here](https://phoenixnap.com/kb/how-to-install-python-3-ubuntu) you find an "how to install python on Ubuntu" tutorial ;)

After the installation, in order to check the *ms_replacR* usage, you can run:

```
python3.8 ms_replacR.py --help

usage: ms_replacR.py [-h] [-out OUTPUTDIR] [-g GENOMIC] [-f FASTA] [-Me METHYLATION] [-i INPUT_WORD] [-o OUTPUT_WORD]

optional arguments:
  -h, --help            show this help message and exit
  
  -out OUTPUTDIR, --outputdir OUTPUTDIR
                        path to new files directory
  -g GENOMIC, --genomic GENOMIC
                        path to file produced by genomic annotator [GFF-file]
  -f FASTA, --fasta FASTA
                        path to genome file [FASTA/FNA-file]
  -Me METHYLATION, --methylation METHYLATION
                        path to file produced by the sequencer [GFF-file]
  -i INPUT_WORD, --input_word INPUT_WORD
                        Element to delete [SYMBOL/STRING]
  -o OUTPUT_WORD, --output_word OUTPUT_WORD
                        Element to insert [SYMBOL/STRING]
 ```

Now to process the files contained in the [toyset](/toyset/) folder, please check **Usage** at *Wizard Level*.

## MeStudio Core
To add.
## MeStudio AnalyzR
MeStudio also implements a post-processing python-based script named ms_analyzR which has to be run on the four GFF3 MeStudio-deriving files. 
In addition, for comparative genomic analyses a “gene_presence_abscence.csv” file produced by [Roary](https://sanger-pathogens.github.io/Roary/) can be used to
define the methylation level and patterns of core and dispensable genome fractions, as well as annotating the genes-coded proteins.
Here's the mandatory and optional fields required:

```
usage: ms_analyzR.py [-h] [-out OUTPUTDIR] [-rr ROARY] [-cds CODING] [-ncds NONCODING] [-inter INTERGENIC] [-ups UPSTREAM] [-prt] [-split]

optional arguments:
  -h, --help            show this help message and exit
  -out OUTPUTDIR, --outputdir OUTPUTDIR
                        path to your output files
  -rr ROARY, --roary ROARY
                        path to ROARY gene_presence_abscence.csv file (OPTIONAL)
  -cds CODING, --coding CODING
                        MOTIF_CDS.gff file
  -ncds NONCODING, --noncoding NONCODING
                        MOTIF_nCDS.gff file
  -inter INTERGENIC, --intergenic INTERGENIC
                        MOTIF_true_intergenic.gff file
  -ups UPSTREAM, --upstream UPSTREAM
                        MOTIF_upstream.gff file
  -prt, --prt_bed       Write in [OUT] tabular per-feature file ready for RCircos (OPTIONAL)
  -split, --split_features
                        Rearrange your input GFFs for chromosomes (OPTIONAL)
```
The `split` flag saves into the output directory the GFFs at “chromosome level” rather than “feature level”. Each GFF produced will be characterized not for
feature (CDS, nCDS, true intergenic and upstream) but by chromosomes (or contigs), maintaining the MeStudio Core derived contents and layout. 

The `prt` flag produces a BED file for each feature in which is reported: 
- the chrom column, with the name of each chromosome or contig
- start of the feature
- end of the feature
- the name of the ID found in that interval
- the number of methylations found for ID
- the protein product of the ID
 
As well as being significant, the information contained in BED files are directly related to an R script (see [src](/src/)) which plots the distribution of the
methylation density for each feature analysed making use of the r-package [circlize](https://jokergoo.github.io/circlize_book/book/). 

Now to process the files contained in the [toyset](/toyset/) folder, please check **Usage** at *Wizard Level*.

## Usage

###### Bard Level

###### Wizard Level
From [src](/src/) folder you can download all the scripts needed to perform the analysis to wizard level.
Here down below you have a step by step process assuming that you are using the files contained into the [toyset](/toyset/) folder and a Linux-based OS.

`ms_replacR`
```
python3.8 ms_replacR.py -out "/path/to/output_folder" -g "FSMMA_genomic.gff" -f "FSMMA_genomic.fna" -Me "FSMMA_methylation.gff" -i "|" -o "_"
```
`MeStudio Core`

```
mscheck -g "FSMMA_genomic.gff" -f "$FSMMA_genomic.fna" -m "$FSMMA_methylation.gff" -o path/to/output_dir --mo "motifs.txt" --cr "circular.txt"
msmine path/to/output_dir/params.ms
msfasta path/to/output_dir/params.ms
msmatch path/to/output_dir/params.ms
msx path/to/output_dir/params.ms
```

`ms_analyzR`
```
python3.8 ms_analyzR.py -out "/path/to/output_folder" -cds "CDS.gff" -ncds "nCDS.gff" -inter "true_intergenic.gff" -ups "upstream.gff" -rr "gene_presence_absence.csv" -split -evo
```

## Results
## Reference
