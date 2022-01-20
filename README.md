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

Inside the [toyset](/toyset/) folder you are going to find the three files needed to start the analysis:

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

MeStudio core is composed of five software components. These components match the nucleotide motifs to the genomic sequence and maps them to the corresponding category, which are extracted from the annotation file. Categories are defined as follows:
- protein-coding gene with accordant (sense) strand (CDS)
- discordant (antisense) strand (nCDS)
- regions that fall between annotated genes (true intergenic, tIG)
- regions upstream to the reading frame of a gene, with accordant strand (US)
 
The components, which in the paper we refer to as executables, are called internally by the pipeline version of MeStudio, in the following order: mscheck, msmine, msfasta, msmatch, msx. There is one additional executable, msread, which gets never called and serves only for debugging purpose for the wizard user. Following is a detailed explanation of what each executable does.

###### mscheck

mscheck requires four mandatory existing files, i.e., the genomic annotation, the genomic sequence, the methylated base calls, and a list containing the motifs to be scanned. It also expects a non-existing directory in which it will dump all the output. Through two optional arguments it is possible to specify the feature type that is to be extracted from the gff (default: gene) and the limit of nucleotides to include during the upstream category calculation. Note that 'gene' is the default feature type and it will ensure that all CDS fall within the reading frame of the gene, and that the intergenic regions will consist in regions that are between two neighboring and non-overlapping genes. It is highly recommended that if you use MeStudio core outside of a pipeline, you select a non existing name for the output argument in your current working directory. Among others, mscheck will produce a binary file called params.ms in the output directory. This file will serve as input for all the other executables that follow. Note that you cannot pass the directory tree to another user since params.ms stores the absolute paths. Run '''msread params.ms''' to see what's inside!

###### msmine

msmine is intended to "mine" the genomic annotation and compute the four categories reported in the introductory paragraph. In order to do so it calculates the length of each input chromosome, and then uses this information to calculate the range of the upstream and/or intergenic regions that are flanking the first and last feature type (hopefully gene, as we keep default). The rest of the data is taken directly from the annotation file, which however was rendered in a binary version by mscheck to minimize the number of type castings and dynamic memory allocations. msmine populates the output directory with binary files which resemble GFF3 tables and contain nCDS, tIG and upstream information which will be used downstream by the other executables.

###### msfasta

msfasta also requires the params.ms binary parameters file produced by mscheck, but it is only in charge of rendering the genomic sequence into a better computer-parsable format. Once again a binary file is written in the output directory.

###### msmatch

The current implementation of msmatch uses a naive matching algorithm to map motif sequences to the reference genome. Albeit the number of comparisons in the worst case is O(m*(n-m+1)), we assumed that the amount of input meaningful motifs seldom exceeds a few dozen. Therefore, we did not make use of more optimized transforms or algorithms for pattern finding, and opted for the most straightforward one. During the naive matching phase, each replicon or chromosome gets loaded into memory in the form of an array of chars (or null terminated C-string) one at a time as both strands are scanned for the presence of the motif sequences, which can hold ambiguity characters. After you run msmatch you will find as many new directories as there were motifs, having each directory containing a binary file which holds motif matching information (such as start, end, strand and exact pattern that was matched). This information, together with the pre-computed categories will be utilized by the last executable in the pipeline, msx, to produce the final output.

###### msx

The last-but-not-least executable is the one that checks whether the modified bases fall within patterns on the genomic sequence, and bins them by category. It performs several nested loop calculations and outputs a number of GFF3 files which report the methylation status in every category. A detailed explanation of the output produced by msx is reported in the following paragraph.

###### MeStudio core output

msx produces GFF3 files that report the methylation status in every category. Following is a description of each of the nine fields:
1. SequenceID		Replicon/Chromosome identifier
2. Source		Normally this would report the algorithm or the procedure that generated this feature. MeStudio reports the category that holds the information (CDS, nCDS, true_intergenic or upstream)
3. Feature Type		Describes what the feature is (default: gene)
4. Feature Start	Start
5. Feature End		End
6. Score		Typically E-values for sequence similarity and P-values for predictions. MeStudio actually reports the distance from the beginning or end of the feature in its reading frame. When a motif is found but no methylation is present on the strand of interest, a 0 is placed instead of a positive integer
7. Strand		The strand on which the methylation lies on. If no methylation is reported, it displays the strand the feature lies on
8. Phase		The name of the (standard) methylated base found in within that motif. When a motif is found but no methylation is present on the strand of interest, a "." is placed instead of a letter
9. Attributes		MeStudio only reports the locus tag of the region in which the motif is found. For CDS and nCDS category this would be the locus tag of the gene that stretches between Feature Start and Feature End. For upstream category this reports the locus tag of the gene to which the region of interest is upstream to, whereas for a true_intergenic category it reports its rightmost gene.

Unless you're interested in your own parsing of the data, you normally wouldn't worry about the binaries and text files that MeStudio core produces, since its output is processed by ms_analyzR. Just so you know, however, we'd like to clarify what columns 6 and 8 report. Column 6 (Score) reports an unsigned integer that represents the distance of the methylated base from the start of the gene it associates with, in its correct reading frame. For example, if you're looking at a CDS GFF3, any methylated base that falls within a gene stretch on the positive strand, and both the gene and the methylation are on the "Plus" strand, you would see the distance from the end of such gene (Feature End). Instead, if the methylation is found on the reverse strand, and also the gene it was found on has its reading frame on the "Minus" strand, then you would see the distance calculated from the beginning of the gene (Feature Start), since it "makes sense" when in its reverse complement form. The distance calculation can be summed up in the following table:


| Feature | Strand | Methylation | distance |
|---------|--------|-------------|----------|
| CDS | + |	+ | from End |
| CDS | + |	-	|	0 |
| CDS | - |	+	|	0 |
| CDS | - |	-	|	from Start |
| nCDS | + | + | 0 |
| nCDS | + | - | from End |
| nCDS | - | + | from Start |
| nCDS | - | - | 0 |
| upstream | + | + | from End |
| upstream | + | - | 0 |
| upstream | - | + | 0 |
| upstream | - | - | from End |
| true_intergenic | * | * | from Start |

But again, this is an internal detail that's taken care of by MeStudio, under the hood. Lastly, the 8th column reports the modified base found within the reading frame of the analyzed motif. When the eighth column displays a dot, it means that the motif returned a positive match in that category but it did not contain any modified base in the correct reading frame, or any modified base at all.


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
