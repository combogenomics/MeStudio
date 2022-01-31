/*
===========================================================================================================================
* AUTHOR:		Christopher Riccardi, PhD student at Computational Biology Dept., University of Firenze
* PROJECT:		Extraction and analysis of methylation features from Pacific Biosciences SMRT reads using MeStudio
* CONTACT:		CHRISTOPHER.RICCARDI@UNIFI.IT
===========================================================================================================================
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#define BUF_SIZE 1024
#define SQID 32
#define LTAG 56
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define cwd getcwd
#define cd chdir
#define MEMORY_FAILURE(pointer) {                                   \
        if (pointer==NULL) {                                        \
            fprintf (stderr, "%s:%d: "                              \
                     " memory allocation failed. (Error: %p)\n",    \
                     __FILE__, __LINE__,                            \
                     strerror (errno));                             \
            exit (EXIT_FAILURE);                                    \
        }                                                           \
    }

typedef struct
{
	char* genomic_gff;		/* [M] */	
	char* sequencer_gff;	/* [M] */
	char* genomic_fasta;	/* [M] */
	char* motifs_list_txt;	/* [M] */
	char* output_dir;		/* [M] */
	char* genomic_bin;
	char* sequencer_bin;
	char* params_bin;
	char* matches_bin;
	char* genomic_fasta_bin;
	char* nCDS_bin;
	char* true_inter_bin;
	char* upstream_bin;
	int upstream_range;
	char type[128];				/* [O], stands for feature type (default: gene) */
} PARAMS;

char g_path_buffer[FILENAME_MAX + 1] = { 0x00 };
char m_path_buffer[FILENAME_MAX + 1] = { 0x00 };
char f_path_buffer[FILENAME_MAX + 1] = { 0x00 };
char o_path_buffer[FILENAME_MAX + 1] = { 0x00 };
char mo_path_buffer[FILENAME_MAX + 1] = { 0x00 };

void print_help()
{
	fprintf(stderr, "\nUsage: mscheck -g <str> -m <str> -f <str> -o <str> --mo <str> [--type <str> --up <int>]\n\n");

	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "\t-g\t\tGenomic annotation GFF3 file path\n");
	fprintf(stderr, "\t-m\t\tSequencer-produced modified-base call GFF3 file path\n");
	fprintf(stderr, "\t-f\t\tGenomic sequence (FASTA) file path\n");
	fprintf(stderr, "\t-o\t\tOutput working directory, must be a non-existing path\n");
	fprintf(stderr, "\t--mo\t\tNewline delimited file with motifs list\n\n");

	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "\t--type\t\tFeature type to parse Genomic GFF file with. Default: gene\n");
	fprintf(stderr, "\t--up\t\tLimit upstream nucleotides count up to <int>. Must be greater than 0\n");
}

int makedir(const char* dir)
{
	struct stat st = {0};
	if (stat(dir, &st) == -1)
	{
	    mkdir(dir, 0777);
	    return 1;
	}
	return 0;
}

void fstrip(char* str)
{	
	// strip after first space starting from C string beginning
	size_t s = strlen(str);
	size_t i = 0;
	char* p = str;
	while(*p!=' ')
	{
		++p;
		++i;
		if(i==s) return;
	}
	if(i==s) return;
	*p = 0x00;
}

void rstrip(char* str)
{	
	// shrink to first standard character ( remove CR (Carriage Return) - LF (Line Feed) etc)
	size_t i = strlen(str);
	do
	{
		--i;
		if(str[i]=='\n' ||
			str[i]=='\r'||
			str[i]=='\t'||
			str[i]=='\f'||
			str[i]=='\v'||
			str[i]==' ') str[i] = 0x00;
	}
	while(i);
}

int parse_args(int argc, char** argv, PARAMS* params)
{	
	// parse command line args
	/*
		Usage: mscheck -g <str> -m <str> -f <str> -o <str> --mo <str> [--type <str> --up <int>]
	*/

	// there are 5 mandatory parameters
	int MNDT_GENOMIC_GFF = 0;
	int MNDT_METHYL_GFF = 0;
	int MNDT_GENOMIC_FNA = 0;
	int MNDT_OUTPUT_DIR = 0;
	int MNDT_MOTIF_TXT = 0;

	params->upstream_range = -1;	// limit nucleotied number in upstream up to this amount. -1: do not limit
	memset(params->type, 0x00, 128 * sizeof(char));
	strcpy(params->type, "gene");		// default feature type is "gene"

	size_t i = 1;
	for(i=1;i<argc;++i)
	{
		if(strcmp(argv[i], "-g")==0)
		{
			if((i + 1) == argc) break;
			//params->genomic_gff = argv[i + 1];
			
			params->genomic_gff = realpath(argv[i + 1], g_path_buffer);
			if(params->genomic_gff==NULL)
			{
				fprintf(stderr, "[ERROR] Could not fetch file path (-g)\n");
				return 0;				
			}
			if(!MNDT_GENOMIC_GFF)
			{
				MNDT_GENOMIC_GFF = 1;
			}
			else
			{
				fprintf(stderr, "[ERROR] -g argument already provided. Please keep command line input neat\n");
				return 0;
			}
			++i;
		}
		else if(strcmp(argv[i], "-m")==0)	// modified bases GFF as produced by the sequencer
		{
			if((i + 1) == argc) break;
			//params->sequencer_gff = argv[i + 1];
			params->sequencer_gff = realpath(argv[i + 1], m_path_buffer);
			if(params->sequencer_gff==NULL)
			{
				fprintf(stderr, "[ERROR] Could not fetch file path (-m)\n");
				return 0;				
			}
			if(!MNDT_METHYL_GFF)
			{
				MNDT_METHYL_GFF = 1;
			}
			else
			{
				fprintf(stderr, "[ERROR] -m argument already provided. Please keep command line input neat\n");
				return 0;
			}
			++i;
		}
		else if(strcmp(argv[i], "-f")==0)
		{
			if((i + 1) == argc) break;
			//params->genomic_fasta = argv[i + 1];
			params->genomic_fasta = realpath(argv[i + 1], f_path_buffer);
			if(params->genomic_fasta==NULL)
			{
				fprintf(stderr, "[ERROR] Could not fetch file path (-f)\n");
				return 0;				
			}
			if(!MNDT_GENOMIC_FNA)
			{
				MNDT_GENOMIC_FNA = 1;
			}
			else
			{
				fprintf(stderr, "[ERROR] -f argument already provided. Please keep command line input neat\n");
				return 0;
			}
			++i;
		}
		else if(strcmp(argv[i], "-o")==0)
		{
			if((i + 1) == argc) break;
			params->output_dir = argv[i + 1];
			if(!MNDT_OUTPUT_DIR)
			{
				MNDT_OUTPUT_DIR = 1;
			}
			else
			{
				fprintf(stderr, "[ERROR] -o argument already provided. Please keep command line input neat\n");
				return 0;
			}
			++i;
		}
		else if(strcmp(argv[i], "--mo")==0)
		{
			if((i + 1) == argc) break;
			//params->genomic_fasta = argv[i + 1];
			params->motifs_list_txt = realpath(argv[i + 1], mo_path_buffer);
			if(params->motifs_list_txt==NULL)
			{
				fprintf(stderr, "[ERROR] Could not fetch file path (--mo)\n");
				return 0;				
			}
			if(!MNDT_MOTIF_TXT)
			{
				MNDT_MOTIF_TXT = 1;
			}
			else
			{
				fprintf(stderr, "[ERROR] --mo argument already provided. Please keep command line input neat\n");
				return 0;
			}
			++i;
		}
		else if(strcmp(argv[i], "--type")==0)
		{
			if((i + 1) == argc) break;
			strcpy(params->type, argv[i + 1]);
			++i;
		}
		else if(strcmp(argv[i], "--up")==0)
		{
			if((i + 1) == argc) break;
			char* p;
			params->upstream_range = (int)strtol(argv[i+1], &p, 10);
			if(params->upstream_range < 1)
			{
				params->upstream_range = -1;
			}
			++i;
		}
		else
		{
			fprintf(stderr, "Unrecognized option \"%s\"\n", argv[i]);
			return 0;
		}
	}

	// check mandatory input
	if(!MNDT_GENOMIC_GFF)
	{
		fprintf(stderr, "[ERROR] -g argument required\n");
		return 0;
	}
	if(!MNDT_METHYL_GFF)
	{
		fprintf(stderr, "[ERROR] -m argument required\n");
		return 0;
	}
	if(!MNDT_GENOMIC_FNA)
	{
		fprintf(stderr, "[ERROR] -f argument required\n");
		return 0;
	}
	if(!MNDT_OUTPUT_DIR)
	{
		fprintf(stderr, "[ERROR] -o argument required\n");
		return 0;
	}
	if(!MNDT_MOTIF_TXT)
	{
		fprintf(stderr, "[ERROR] --mo argument required\n");
		return 0;
	}
	// log upstream range
	if(params->upstream_range < 1)
	{
		fprintf(stderr, "Upstream range not set. Will include all nucleotides in upstream gff\n");
	}
	else
	{
		fprintf(stderr, "Upstream range set to include at most %d nucleotides before feature start\n", params->upstream_range);
	}


	fprintf(stderr, "Mandatory Input Parameters:\n");
	fprintf(stderr, "params->genomic_gff %s\n", params->genomic_gff);
	fprintf(stderr, "params->sequencer_gff %s\n", params->sequencer_gff);
	fprintf(stderr, "params->genomic_fasta %s\n", params->genomic_fasta);
	fprintf(stderr, "params->motifs_list_txt %s\n", params->motifs_list_txt);
	fprintf(stderr, "params->output_dir %s\n", params->output_dir);

	// create output working directory

	if(!makedir((const char*)params->output_dir)) 
	{
		fprintf(stderr, "Could not create directory.\n");
		fprintf(stderr, "Either path already exsists or you do not have sufficient privileges to run this program.\n");
		return 0;
	}

	char* out_tmp = realpath(params->output_dir, o_path_buffer);
	if(out_tmp==NULL)
	{
		fprintf(stderr, "[ERROR] Could not fetch created file path (-o)\n");
		return 0;				
	}
	params->output_dir = out_tmp;
	fprintf(stderr, "Created directory %s\n", params->output_dir);

	/* Update genomic bin paths */
	size_t genomic_bin_len = (1 + strlen(params->output_dir) + strlen("/genomic.ms")) * sizeof(char);
	params->genomic_bin = malloc(genomic_bin_len);
	MEMORY_FAILURE(params->genomic_bin);

	memset(params->genomic_bin, 0x00, genomic_bin_len);
	strcpy(params->genomic_bin, params->output_dir);
	strcat(params->genomic_bin, "/genomic.ms");
	fprintf(stderr, "params->genomic_bin %s\n", params->genomic_bin);

	/* Update sequencer bin paths */
	size_t sequencer_bin_len = (1 + strlen(params->output_dir) + strlen("/sequencer.ms")) * sizeof(char);
	params->sequencer_bin = malloc(sequencer_bin_len);
	MEMORY_FAILURE(params->sequencer_bin);

	memset(params->sequencer_bin, 0x00, sequencer_bin_len);
	strcpy(params->sequencer_bin, params->output_dir);
	strcat(params->sequencer_bin, "/sequencer.ms");
	fprintf(stderr, "params->sequencer_bin %s\n", params->sequencer_bin);

	/* Update params bin paths */
	size_t params_bin_len = (1 + strlen(params->output_dir) + strlen("/params.ms")) * sizeof(char);
	params->params_bin = malloc(params_bin_len);
	MEMORY_FAILURE(params->params_bin);

	memset(params->params_bin, 0x00, params_bin_len);
	strcpy(params->params_bin, params->output_dir);
	strcat(params->params_bin, "/params.ms");
	fprintf(stderr, "params->params_bin %s\n", params->params_bin);

	/* Update matches bin paths */
	size_t matches_bin_len = (1 + strlen(params->output_dir) + strlen("/matches.ms")) * sizeof(char);
	params->matches_bin = malloc(matches_bin_len);
	MEMORY_FAILURE(params->matches_bin);
	
	memset(params->matches_bin, 0x00, matches_bin_len);
	strcpy(params->matches_bin, params->output_dir);
	strcat(params->matches_bin, "/matches.ms");
	fprintf(stderr, "params->matches_bin %s\n", params->matches_bin);

	/* Update genomic_fasta_bin paths */
	size_t genomic_fasta_bin_len = (1 + strlen(params->output_dir) + strlen("/genomic_fasta.ms")) * sizeof(char);
	params->genomic_fasta_bin = malloc(genomic_fasta_bin_len);
	MEMORY_FAILURE(params->genomic_fasta_bin);
	
	memset(params->genomic_fasta_bin, 0x00, genomic_fasta_bin_len);
	strcpy(params->genomic_fasta_bin, params->output_dir);
	strcat(params->genomic_fasta_bin, "/genomic_fasta.ms");
	fprintf(stderr, "params->genomic_fasta_bin %s\n", params->genomic_fasta_bin);

	/* Update nCDS bin paths */
	size_t nCDS_bin_len = (1 + strlen(params->output_dir) + strlen("/nCDS.ms")) * sizeof(char);
	params->nCDS_bin = malloc(nCDS_bin_len);
	MEMORY_FAILURE(params->nCDS_bin);
	
	memset(params->nCDS_bin, 0x00, nCDS_bin_len);
	strcpy(params->nCDS_bin, params->output_dir);
	strcat(params->nCDS_bin, "/nCDS.ms");
	fprintf(stderr, "params->nCDS_bin %s\n", params->nCDS_bin);

	/* Update true intergenic bin paths */
	size_t true_inter_bin_len = (1 + strlen(params->output_dir) + strlen("/true_intergenic.ms")) * sizeof(char);
	params->true_inter_bin = malloc(true_inter_bin_len);
	MEMORY_FAILURE(params->true_inter_bin);
	
	memset(params->true_inter_bin, 0x00, true_inter_bin_len);
	strcpy(params->true_inter_bin, params->output_dir);
	strcat(params->true_inter_bin, "/true_intergenic.ms");
	fprintf(stderr, "params->true_inter_bin %s\n", params->true_inter_bin);

	/* Update upstream bin paths */
	size_t upstream_bin_len = (1 + strlen(params->output_dir) + strlen("/upstream.ms")) * sizeof(char);
	params->upstream_bin = malloc(upstream_bin_len);
	MEMORY_FAILURE(params->upstream_bin);
	
	memset(params->upstream_bin, 0x00, upstream_bin_len);
	strcpy(params->upstream_bin, params->output_dir);
	strcat(params->upstream_bin, "/upstream.ms");
	fprintf(stderr, "params->upstream_bin %s\n", params->upstream_bin);
	fprintf(stderr, "Using feature type: %s\n", params->type);
	return 1;
}

int reduce_gff(char* path, char* output, char* type, PARAMS* params)
{
	char buffer[BUF_SIZE];

	FILE* fptr;
	fptr = fopen(path, "r");
	if (fptr == NULL)
	{
		printf("[ERROR] Could not open GFF3 file\n");
		return 0;
	}

	int sequencer_gff = 0;
	if(strncmp(&output[strlen(output) - strlen("sequencer.ms")], "sequencer.ms", strlen("sequencer.ms"))==0) 
	{
		sequencer_gff = 1;
		fprintf(stderr, "Reducing sequencer-derived GFF ");
	}
	else fprintf(stderr, "Reducing general GFF ");
	
	FILE* pFile;
	pFile=fopen(output, "wb");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (%s)\n", output);
	}

	char seqid[SQID] = { 0x00 };
	char feature[SQID] = { 0x00 };
	char str_toi[SQID] = { 0x00 };
	unsigned long int start = 0;
	unsigned long int end = 0;
	unsigned long int mescore = 0;
	char strand = '.';

	char locus_tag[LTAG] = { 0x00 };
	unsigned long int bytes_written = 0;

	while (fgets(buffer, BUF_SIZE, fptr) != 0)
	{
		if(strncmp(buffer, "##FASTA", 7)==0) break;
		if(buffer[0]!='#')
		{	
			memset(seqid, 0x00, SQID * sizeof(char));
			memset(feature, 0x00, SQID * sizeof(char));
			memset(str_toi, 0x00, SQID * sizeof(char));
			strand = 0x00;
			memset(locus_tag, 0x00, LTAG * sizeof(char));
			
			start = 0;
			end = 0;

			/* Get seqid, update cursor till a tab is found */
			size_t i = 0;
			char* p = buffer;
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				seqid[i] = *p;
				++p;
				++i;
			}
			seqid[i] = 0x00;
			++p;
			i = 0;

			/* Get feature, update cursor till a tab is found */
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				++p;
			}
			++p;
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				feature[i] = *p;
				++p;
				++i;
			}
			feature[i] = 0x00;
			i = 0;
			++p;

			/* Get start, update cursor till a tab is found */
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				str_toi[i] = *p;
				++p;
				++i;
			}
			i = 0;
			char* ptr = p;
			start = (unsigned long int)strtol(str_toi, &ptr, 10);
			++p;
			memset(str_toi, 0x00, SQID * sizeof(char));

			/* Get end, update cursor till a tab is found */
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				str_toi[i] = *p;
				++p;
				++i;
			}
			i = 0;
			ptr = p;
			end = (unsigned long int)strtol(str_toi, &ptr, 10);
			++p;

			/* Get strand, update cursor. we don't know if score was written or left a dot */
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				++p;
			}
			++p;
			strand = *p;
			++p;
			/* same here */
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				++p;
			}
			++p;
			while(*p!='\n' && *p!='\t' && *p!='\r')
			{
				++p;
			}
			++p;
			if(strncmp(p, "ID=", 3)==0)	// GFF3 standard requires ID= keyword. If not present, write NA
			{
				if(strncmp(p, "ID=gene-", 7)==0)
				{
					p+=8;
					while(*p!='\n' && 
						*p!='\t' && 
						*p!='\r' && 
						*p!=',' && 
						*p!=';' && 
						*p!='|' && 
						*p!=' ')
					{
						locus_tag[i] = *p;
						++p;
						++i;
					}
					locus_tag[i] = 0x00;	//null terminate here
					--i;		
				}
				else
				{
					p+=3;
					while(*p!='\n' && 
					*p!='\t' && 
					*p!='\r' && 
					*p!=',' && 
					*p!=';' && 
					*p!='|' && 
					*p!=' ')
					{
						if(strncmp(p, "gene", 4)==0) break;
						if(strncmp(p, "_gene", 5)==0) break;
						if(strncmp(p, "-gene", 5)==0) break;
						locus_tag[i] = *p;
						++p;
						++i;
					}
					locus_tag[i] = 0x00;	//null terminate here
					--i;		
				}
			}
			else
			{
				strcpy(locus_tag, "NA");
			}

			// write binary file
			if(!sequencer_gff)	/* file not derived from sequencer */
			{
				if(strcmp(feature, type)==0)
				{
					fwrite(&seqid[0], SQID * sizeof(char), 1, pFile );
					fwrite(&feature[0], SQID * sizeof(char), 1, pFile );
					fwrite(&start, sizeof(unsigned long int), 1, pFile );
					fwrite(&end, sizeof(unsigned long int), 1, pFile );
					fwrite(&mescore, sizeof(unsigned long int), 1, pFile );
					fwrite(&strand, sizeof(char), 1, pFile );
					fwrite(&locus_tag[0], LTAG * sizeof(char), 1, pFile );
					bytes_written += SQID * sizeof(char);
					bytes_written += SQID * sizeof(char);
					bytes_written += sizeof(unsigned long int);
					bytes_written += sizeof(unsigned long int);
					bytes_written += sizeof(unsigned long int);
					bytes_written += sizeof(char);
					bytes_written += LTAG * sizeof(char);
				}						

			}
			else	/* file derived from sequencer */
			{
				//fprintf(stderr, "%s\t%s\t%ld\t%ld\t%c\t%s\n", seqid, feature, start, end, strand, locus_tag);
				fwrite(&seqid[0], SQID * sizeof(char), 1, pFile );
				fwrite(&feature[0], SQID * sizeof(char), 1, pFile );
				fwrite(&start, sizeof(unsigned long int), 1, pFile );
				fwrite(&end, sizeof(unsigned long int), 1, pFile );
				fwrite(&mescore, sizeof(unsigned long int), 1, pFile );
				fwrite(&strand, sizeof(char), 1, pFile );
				fwrite(&locus_tag[0], LTAG * sizeof(char), 1, pFile );
				bytes_written += SQID * sizeof(char);
				bytes_written += SQID * sizeof(char);
				bytes_written += sizeof(unsigned long int);
				bytes_written += sizeof(unsigned long int);
				bytes_written += sizeof(unsigned long int);
				bytes_written += sizeof(char);
				bytes_written += LTAG * sizeof(char);
			}
		}
	}
	fclose(pFile);
	fclose(fptr);
	fprintf(stderr, "... written %ld bytes\n", bytes_written);

	return 1;
}

void init_params(PARAMS* params)
{
	// initialize and populate params.ms file
    FILE* pFile;
    pFile = fopen(params->params_bin, "wb");
    if (pFile == NULL)
    {
		fprintf(stderr, "[ERROR] Could not point to params.ms file\n");
		return;
    }

    /* Write parameters to binary file using fwrite */
    size_t tmp = strlen(params->genomic_gff);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->genomic_gff[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->sequencer_gff);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->sequencer_gff[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->genomic_fasta);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->genomic_fasta[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->motifs_list_txt);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->motifs_list_txt[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->output_dir);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->output_dir[0], tmp * sizeof(char), 1, pFile );

	tmp = strlen(params->genomic_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->genomic_bin[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->sequencer_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->sequencer_bin[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->matches_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->matches_bin[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->genomic_fasta_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->genomic_fasta_bin[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->nCDS_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->nCDS_bin[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->true_inter_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->true_inter_bin[0], tmp * sizeof(char), 1, pFile );
	tmp = strlen(params->upstream_bin);
    fwrite(&tmp, sizeof(size_t), 1, pFile );
    fwrite(&params->upstream_bin[0], tmp * sizeof(char), 1, pFile );

    /* Write upstream range */
	int uprng = params->upstream_range;
	fwrite(&uprng, sizeof(int), 1, pFile );
    fclose(pFile);
}

void append_genome(PARAMS* params)
{
	// seqids need to take a different route, as they are variable through samples
   	FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    char header[SQID] = { 0x00 };
    fp = fopen(params->genomic_fasta, "r");
    if (fp == NULL)
    {
    	fprintf(stderr, "[ERROR] Could not point to file\n");
    	return;
    }

    size_t header_count = 0;
    while ((read = getline(&line, &len, fp)) != -1)
    {
        if(strlen(line)>0)
        {
            if(line[0]=='>') ++header_count;
        }
    }
    fseek(fp, 0, SEEK_SET);

    FILE* pFile;
    // open in append bytes mode
    pFile = fopen(params->params_bin, "a+b");
    if (pFile == NULL)
    {
		fprintf(stderr, "[ERROR] Could not point to file\n");
		return;
    }
    // write 1 * sizeof(size_t) bytes with number of seqids (herein defined as headers)
    fseek(pFile, 0, SEEK_END);
    fwrite(&header_count, sizeof(size_t), 1, pFile);

    int state = 0;
    size_t string_len = 0;

    while ((read = getline(&line, &len, fp)) != -1)
    {
    	rstrip(line);
        if(line[0] == '>' && state==0)
        {
            memcpy(header, &line[1], SQID);//strncpy(header, &line[1], SQID);
            if(strlen(header) == SQID)
            	header[strlen(header) - 1] = 0x00;	// trim trailing character if header == SQID
        }

        if (line[0]!='>')
        {
            state = 1;
            string_len += strlen(line);
        }

        if (line[0] == '>' && state)
        {
            fprintf(stderr, "%s: %ld\n", &header[0], string_len);
            fwrite(&header[0], SQID * sizeof(char), 1, pFile);	// &header[1] : skip leading '>' sign
            fwrite(&string_len, sizeof(size_t), 1, pFile);
            string_len = 0;
            memset(header, 0x00, SQID * sizeof(char));
            memcpy(header, &line[1], SQID);//strncpy(header, &line[1], SQID);
            if(strlen(header) == SQID)
            	header[strlen(header) - 1] = 0x00;
        }
    }
    fwrite(&header[0], SQID * sizeof(char), 1, pFile);
    fwrite(&string_len, sizeof(size_t), 1, pFile);
    fprintf(stderr, "%s: %ld\n", &header[0], string_len);
    string_len = 0;
    memset(header, 0x00, SQID * sizeof(char));

    fclose(fp);
    fclose(pFile);
    if (line)
    {
        free(line);
    }
}

int main(int argc, char** argv)
{
	PARAMS params;
	if(!parse_args(argc, argv, &params)) 
	{
		print_help();
		exit(1);
	}
	reduce_gff(params.genomic_gff, params.genomic_bin, params.type, &params);
	reduce_gff(params.sequencer_gff, params.sequencer_bin, params.type, &params);

	fprintf(stderr, "\tSequence summary\n");
	init_params(&params);
	append_genome(&params);

	free(params.genomic_bin);
	free(params.sequencer_bin); 
	free(params.params_bin);
	free(params.matches_bin);
	free(params.genomic_fasta_bin); 
	free(params.nCDS_bin);
	free(params.true_inter_bin); 
	free(params.upstream_bin);
	return 0;
}
