# MeStudio
MeStudio is a tool which allows us to analyse and combine genome-wide methylation profiles with genomic features.
MeStudio connects several procedural software components that can be run either individually or as part of a pipeline.

When you run MeStudio as a pipeline (see *Installing MeStudio*, `./mestudio -h`) it automatically calls all the scripts and executables it needs to perform genome-wide methylation profiles analysis. Running MeStudio this way is fairly easy:

```
Usage: mestudio -f <str> -g <str> -Me <str> -mo <str> -out <str> [-rr <str>]

Mandatory arguments
-f    <str>           genomic sequence file
-g    <str>           genomic annotation file
-Me   <str>           methylated base calls file
-mo   <str>           newline delimited motifs list
-out  <str>           output directory

Optional arguments
-rr   <str>          "gene_presence_absence.csv" file produced by Roary
```

If you are using *MeStudio*, please cite out work https://doi.org/10.1101/2022.03.23.485463

*MeStudio version 1.0*

## Table of Contents

- [Installing MeStudio](#Installing-MeStudio)
- [MeStudio ReplacR](#MeStudio-ReplacR)
- [MeStudio Core](#MeStudio-Core)
- [MeStudio AnalyzR](#MeStudio-AnalyzR)
- [Results](#Results)
- [Reference](#Reference)

## Installing MeStudio

You have to render `install` executable by typing:
```
chmod +x install
```
Then you can install it with the following command:
```
sudo ./install
```
If you don't have the administrator privileges, please add the directory in which executables are present to the `$PATH`. [Here](https://linuxize.com/post/how-to-add-directory-to-path-in-linux/) you can find how to do it.

Ok, now you're ready to launch `mestudio`. Simply type:
```
mestudio
```
And you'll get:
```
Usage: mestudio -f <str> -g <str> -Me <str> -mo <str> -out <str> [-rr <str>]

Mandatory arguments
-f    <str>           genomic sequence file
-g    <str>           genomic annotation file
-Me   <str>           methylated base calls file
-mo   <str>           newline delimited motifs list
-out  <str>           output directory

Optional arguments
-rr   <str>          "gene_presence_absence.csv" file produced by Roary
```
Please feel free to use the files contained in the [dataset](/dataset/) folder to start a trial analysis or you can directly use your own.
Inside the [dataset](/dataset/) folder you are going to find the three files needed to start the analysis:

```FSMMA_genomic.fna``` is the the genome of FSMMA strain of *Sinorhizobium meliloti*. You can get more information [here](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/382/?utm_source=None&utm_medium=referral&utm_campaign=KnownItemSensor:taxname).

```FSMMA_genomic.gff``` is the product of the genomic annotation (executed via [Prokka](https://github.com/tseemann/prokka) annotator). The file is reported in GFF3 format.

```FSMMA_methylation.gff``` is the GFF3 file with methylation positions obtained through the sequencer. We performed the analysis using the PacBio RT-SMRT sequencing.

N.B. please note that Prokka and Roary (see *Installing MeStudio*) can be easily found to [Galaxy](https://usegalaxy.org), the installation of these tools is not strictly needed.

Once this is done, you are going to find all the results inside the directory you started the analysis from. Here's reported a "tree visualization" of the directories hierarchy, assuming that you downloaded the files from the [dataset](/dataset/). 

```
├── FSMMA_genomic.fasta
├── FSMMA_genomic.gff
├── FSMMA_methylation.gff
├── gene_presence_absence.csv
├── motifs.txt
└── replout
    ├── core
    │   ├── GANTC
    │   │   ├── GANTC_CDS.gff
    │   │   ├── GANTC.ms
    │   │   ├── GANTC_nCDS.gff
    │   │   ├── GANTC_true_intergenic.gff
    │   │   ├── GANTC_upstream.gff
    │   │   └── results
    │   │       ├── output_CDS.bed
    │   │       ├── output_intergenic.bed
    │   │       ├── output_nCDS.bed
    │   │       ├── output_upstream.bed
    │   │       ├── results_cds_scatterplot.png
    │   │       ├── results_intergenic_scatterplot.png
    │   │       ├── results_ncds_scatterplot.png
    │   │       ├── results_upstream_scatterplot.png
    │   │       └── Rplots.pdf
    │   ├── genomic_fasta.ms
    │   ├── genomic.ms
    │   ├── matches.ms
    │   ├── nCDS.ms
    │   ├── params.ms
    │   ├── sequencer.ms
    │   ├── true_intergenic.ms
    │   └── upstream.ms
    ├── ms_circ.R

```
Where:
`FSMMA_genomic.fasta`, `FSMMA_genomic.gff`, `FSMMA_methylation.gff` are the files containing respectevely the genome, genome's annotation and methylation data. 

`gene_presence_absence.csv` is the file produced by Roary with proteins annotation.

`motifs.txt` is the text file containing all the motifs you want to look for in the genome. In the above example we used GANTC motif only.

`replout`, `core`, `results` are the directories created respectively from `ms_replacR`, `ms_core` and `ms_analyzR`. Inside these folders you can find the newly produced tabular files as GFFs and BEDs.

All the files with `.ms` extension are produced by `ms_core`.

Jump to the *Results* chapter to see the results produced by *MeStudio*.


## MeStudio ReplacR
In order to properly run *MeStudio Core*, a pre-processing python-based script named *ms_replacR* has been implemented and is highly suggested to be used. You can find the source code here.
*ms_replacR* expects four arguments:
1. *output directory*, in which results and log files will be written
2. *genomic annotation*, in the GFF3 format
3. *methylation annotation*, a sequencer-produced modified base calls in the GFF3 format
4. *genomic sequence*, in fasta or fna file format

As a matter of fact MeStudio requires consistent formatting as far as the sequence identifiers are concerned but sometimes the annotation process can lead to genomic headers alterations. By default, *ms_replacR* changes pipe symbol with the underscore for what concerns “seqid” fields identity and makes sure that numbers and quality of the contigs reported on your file are good.

An example is reported below:

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

Moreover, depending on the annotator, we noticed that sometimes the order in which seqids are found, both in the GFF and fasta, differs.
*ms_replacR* also fixes this kind of anomaly in order to avoid biased results.


###### How do I run ms_replacR?

Before running *ms_replacR*, please make sure to have python3.8 (or above) properly installed on your computer. [Here](https://phoenixnap.com/kb/how-to-install-python-3-ubuntu) you find an "how to install python on Ubuntu" tutorial ;)

After the installation, in order to check the *ms_replacR* usage, you can run:

```
python3.8 ms_replacR.py --help

usage: ms_replacR.py [-h] [-out OUTPUTDIR] [-g GENOMIC] [-f FASTA] [-Me METHYLATION]

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
 ```

## MeStudio Core

MeStudio *core* is composed of five software components. These components match the nucleotide motifs to the genomic sequence and maps them to the corresponding category, which are extracted from the annotation file. Categories are defined as follows:
- protein-coding gene with accordant (sense) strand (**CDS**)
- discordant (antisense) strand (**nCDS**)
- regions that fall between annotated genes (true intergenic, **tIG**)
- regions upstream to the reading frame of a gene, with accordant strand (**US**)
 
The components, which in the paper we refer to as executables, are called internally by the pipeline version of MeStudio, in the following order: *mscheck*, *msmine*, *msfasta*, *msmatch*, *msx*. There is one additional executable, *msread*, which gets never called and serves only for debugging purpose for the wizard user. Following is a detailed explanation of what each executable does.

###### mscheck

*mscheck* requires four mandatory existing files, i.e., the genomic annotation, the genomic sequence, the methylated base calls, and a list containing the motifs to be scanned. It also expects a non-existing directory in which it will dump all the output. Through two optional arguments it is possible to specify the feature type that is to be extracted from the gff (default: gene) and the limit of nucleotides to include during the upstream category calculation. Note that 'gene' is the default feature type and it will ensure that all CDS fall within the reading frame of the gene, and that the intergenic regions will consist in regions that are between two neighboring and non-overlapping genes. It is highly recommended that if you use MeStudio core outside of a pipeline, you select a non existing name for the output argument in your current working directory. Among others, mscheck will produce a binary file called params.ms in the output directory. This file will serve as input for all the other executables that follow. Note that you cannot pass the directory tree to another user since params.ms stores the absolute paths. Run `msread params.ms` to see what's inside!

###### msmine

*msmine* is intended to "mine" the genomic annotation and compute the four categories reported in the introductory paragraph. In order to do so it calculates the length of each input chromosome, and then uses this information to calculate the range of the upstream and/or intergenic regions that are flanking the first and last feature type (hopefully gene, as we keep default). The rest of the data is taken directly from the annotation file, which however was rendered in a binary version by mscheck to minimize the number of type castings and dynamic memory allocations. msmine populates the output directory with binary files which resemble GFF3 tables and contain nCDS, tIG and upstream information which will be used downstream by the other executables.

###### msfasta

*msfasta* also requires the `params.ms` binary parameters file produced by *mscheck*, but it is only in charge of rendering the genomic sequence into a better computer-parsable format. Once again a binary file is written in the output directory.

###### msmatch

The current implementation of *msmatch* uses a naive matching algorithm to map motif sequences to the reference genome. Albeit the number of comparisons in the worst case is O(m*(n-m+1)), we assumed that the amount of input meaningful motifs seldom exceeds a few dozen. Therefore, we did not make use of more optimized transforms or algorithms for pattern finding, and opted for the most straightforward one. During the naive matching phase, each replicon or chromosome gets loaded into memory in the form of an array of chars (or null terminated C-string) one at a time as both strands are scanned for the presence of the motif sequences, which can hold ambiguity characters. After you run msmatch you will find as many new directories as there were motifs, having each directory containing a binary file which holds motif matching information (such as start, end, strand and exact pattern that was matched). This information, together with the pre-computed categories will be utilized by the last executable in the pipeline, msx, to produce the final output.

###### msx

The last-but-not-least executable is the one that checks whether the modified bases fall within patterns on the genomic sequence, and bins them by category. *msx* performs several nested loop calculations and outputs a number of GFF3 files which report the methylation status in every category. A detailed explanation of the output produced by msx is reported in the following paragraph.

###### MeStudio core output

*msx* produces GFF3 files that report the methylation status in every category. Following is a description of each of the nine fields:
1. SequenceID		Replicon/Chromosome identifier
2. Source		Normally this would report the algorithm or the procedure that generated this feature. MeStudio reports the category that holds the information (CDS, nCDS, true_intergenic or upstream)
3. Feature Type		Describes what the feature is (default: gene)
4. Feature Start	Start
5. Feature End		End
6. Score		Typically E-values for sequence similarity and P-values for predictions. MeStudio actually reports the distance from the beginning or end of the feature in its reading frame. When a motif is found but no methylation is present on the strand of interest, a 0 is placed instead of a positive integer
7. Strand		The strand on which the methylation lies on. If no methylation is reported, it displays the strand the feature lies on
8. Phase		The name of the (standard) methylated base found in within that motif. When a motif is found but no methylation is present on the strand of interest, a "." is placed instead of a letter
9. Attributes		MeStudio only reports the locus tag of the region in which the motif is found. For CDS and nCDS category this would be the locus tag of the gene that stretches between Feature Start and Feature End. For upstream category this reports the locus tag of the gene to which the region of interest is upstream to, whereas for a true_intergenic category it reports its rightmost gene.

Unless you're interested in your own parsing of the data, you normally wouldn't worry about the binaries and text files that MeStudio *core* produces, since its output is processed by *ms_analyzR*. Just so you know, however, we'd like to clarify what columns 6 and 8 report. Column 6 (Score) reports an unsigned integer that represents the distance of the methylated base from the start of the gene it associates with, in its correct reading frame. For example, if you're looking at a CDS GFF3, any methylated base that falls within a gene stretch on the positive strand, and both the gene and the methylation are on the "Plus" strand, you would see the distance from the end of such gene (Feature End). Instead, if the methylation is found on the reverse strand, and also the gene it was found on has its reading frame on the "Minus" strand, then you would see the distance calculated from the beginning of the gene (Feature Start), since it "makes sense" when in its reverse complement form. The distance calculation can be summed up in the following table:


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
MeStudio also implements a post-processing python-based script named *ms_analyzR* which has to be run on the four GFF3 MeStudio-deriving files. 
In addition, for comparative genomic analyses a “gene_presence_abscence.csv” file produced by [Roary](https://sanger-pathogens.github.io/Roary/) can be used to define the methylation level and patterns of core and dispensable genome fractions, as well as annotating the genes-coded proteins.
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
  -bed, --make_bed       Write in [OUT] tabular per-feature file ready for RCircos (OPTIONAL)
  -split, --split_by_chromosome
                        Rearrange your input GFFs for chromosomes (OPTIONAL)
```
The `-split` flag saves into the output directory the GFFs at “chromosome level” rather than “feature level”. Each GFF produced will be characterized not for feature (CDS, nCDS, true intergenic and upstream) but by chromosomes (or contigs), maintaining the MeStudio Core derived contents and layout. 

The `-bed` flag produces a BED file for each feature in which is reported: 
- the chrom column, with the name of each chromosome or contig
- start of the feature
- end of the feature
- the name of the ID found in that interval
- the number of methylations found for ID
- the protein product of the ID
 
As well as being significant, the information contained in BED files are directly related to an R script (see [src](/src/)) which plots the distribution of the methylation density for each feature analysed making use of the r-package [circlize](https://jokergoo.github.io/circlize_book/book/). 

## Results

`ms_analyzR` produces four different type of results: statistics, tabular, distributions and circular plot.

Below are reported these results.

*Statistics*
```
Total number of FOUND genes: 7159
Parsing -cds file..
	 Total number of METHYLATED genes: 4193 (58.57 %)
	 Total number of METHYLATIONS found: 8721
	 Max number of methylations found: 17 on gene JPHAALHC_06173
	 Roary annotation: JPHAALHC_06173 = hypothetical protein -> Core genome


Parsing -ncds file..
	 Total number of METHYLATED genes: 4193 (58.57 %)
	 Total number of METHYLATIONS found: 8726
	 Max number of methylations found: 17 on gene JPHAALHC_06173
	 Roary annotation: JPHAALHC_06173 = hypothetical protein -> Core genome


Parsing -inter file..
	 Total number of METHYLATED regions: 1590 (22.21 %)
	 Total number of METHYLATIONS found: 2500
	 Max number of methylations found: 14 on gene JPHAALHC_02514
	 Roary annotation: JPHAALHC_02514 = Acetoacetyl-CoA reductase -> Core genome


Parsing -ups file..
	 Total number of METHYLATED regions: 2719 (37.98 %)
	 Total number of METHYLATIONS found: 13761
	 Max number of methylations found: 84 up to gene JPHAALHC_03228
	 Roary annotation: JPHAALHC_03228 = Motility protein A -> Core genome
```

*Tabular*

| chrom | chromStart | chromEnd | name | score | protein |
|-------|------------|----------|------|-------|--------|
| 000000F_arrow |	1319 | 2620 |	JPHAALHC_00003 | 2 | putative zinc protease
| 000000F_arrow | 5122 | 6228 | JPHAALHC_00006 | 1 | hypothetical protein
| 000000F_arrow | 6357 | 7052 | JPHAALHC_00007 | 3 | 6-phosphogluconate phosphatase
| 000000F_arrow | 7237 | 7920 | JPHAALHC_00008 | 1 | HTH-type transcriptional regulator CmtR
| 000000F_arrow | 11274 | 11780 | JPHAALHC_00012 | 1 | hypothetical protein
| 000000F_arrow | 12999 | 16460 | JPHAALHC_00014 | 3 | Chromosome partition protein Smc
| ... | ... | ... | ... | ... | ... |

*Distributions*

![GCRDB_analyzed_cds_scatterplot](https://user-images.githubusercontent.com/97732011/161234235-57dd01c5-e8d3-40c0-aa8f-b5c3c253b9b5.png)
![GCRDB_analyzed_ncds_scatterplot](https://user-images.githubusercontent.com/97732011/161234359-26af5ea3-babf-4546-869e-b7991c910dae.png)
![GCRDB_analyzed_intergenic_scatterplot](https://user-images.githubusercontent.com/97732011/161234268-556c97dc-2356-449c-a34d-5a22c5873f6a.png)
![GCRDB_analyzed_upstream_scatterplot](https://user-images.githubusercontent.com/97732011/161234407-9d0695ff-bcb7-48c1-9363-9fb47b4fa6c5.png)

*Circular plot*

<img width="904" alt="Schermata 2022-04-01 alle 10 51 27" src="https://user-images.githubusercontent.com/97732011/161234104-7a509acd-e2f0-44be-a1d1-dd45a6c4ce44.png">

## Reference

Please, cite us from doi: https://doi.org/10.1101/2022.03.23.485463
