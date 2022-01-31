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
#define cwd get_cwd
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
	size_t headers_count;
	char** headers;
	size_t* sizes;
} PARAMS;

typedef struct
{
	char seqid[SQID];
	char feature[SQID];
	unsigned long int start;
	unsigned long int end;
	unsigned long int mscore;
	char strand;
	char locus_tag[LTAG];	
} GFFline;

void read_params_bin(PARAMS* params, char* path)	/* WARNING: slightly different from the others */
{
	FILE * pFile;
	long lSize;
	size_t result;
	pFile = fopen (path, "rb");
	if (pFile==NULL) 
	{
	fprintf(stderr, "[ERROR] Could not read params.ms\n"); 
	exit (1);
	}
	// obtain file size
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	if(!lSize) exit(1);
	size_t i = 0;
 	size_t size = 0;
 	int bool = 0;
 	char tmp[4096] = { 0x00 };

 	/* genomic_gff */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->genomic_gff = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->genomic_gff);
	
	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->genomic_gff, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* sequencer_gff */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->sequencer_gff = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->sequencer_gff);
	
	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->sequencer_gff, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* genomic_fasta */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->genomic_fasta = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->genomic_fasta);
	
	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->genomic_fasta, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* motifs_list_txt */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->motifs_list_txt = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->motifs_list_txt);
	
	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->motifs_list_txt, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* output_dir */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->output_dir = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->output_dir);
	
	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->output_dir, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* genomic_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->genomic_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->genomic_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->genomic_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* sequencer_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->sequencer_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->sequencer_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->sequencer_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* matches_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->matches_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->matches_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->matches_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* genomic_fasta_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->genomic_fasta_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->genomic_fasta_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->genomic_fasta_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* nCDS_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->nCDS_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->nCDS_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->nCDS_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* true_inter_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->true_inter_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->true_inter_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->true_inter_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* upstream_bin */
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->upstream_bin = malloc((1 + size) * sizeof(char));
	MEMORY_FAILURE(params->upstream_bin);

	result = fread (&tmp, size * sizeof(char), 1, pFile);
	strcpy(params->upstream_bin, tmp);
	memset(tmp, 0x00, size * sizeof(char));

 	/* upstream_range */
	result = fread (&bool, sizeof(int), 1, pFile);
	if(!result) return;
	params->upstream_range = bool;	

	/* STATIC PARAMS TERMINATED */
	/* READING VARIABLE PARAMS */

	size = 0;
	result = fread (&size, sizeof(size_t), 1, pFile);
	if(!result) return;
	params->headers_count = size;

	params->headers = malloc(params->headers_count * sizeof(char*));
	params->sizes = malloc(params->headers_count * sizeof(size_t));
	MEMORY_FAILURE(params->headers);
	MEMORY_FAILURE(params->sizes);

	i = 0;
	for(i=0;i<params->headers_count;++i)
	{
		params->headers[i] = malloc(SQID * sizeof(char));
		MEMORY_FAILURE(params->headers[i]);
		memset(params->headers[i], 0x00, SQID * sizeof(char));
	}
 	/* copy seq id strings */
	i = 0;
	for(i=0;i<params->headers_count; ++i)
	{
		size = 0;
		memset(tmp, 0x00, SQID * sizeof(char));	

		result = fread (&tmp, SQID * sizeof(char), 1, pFile);

		if(!result)
		{
			return;
		}
		strcpy(params->headers[i], tmp);
		memset(tmp, 0x00, size * sizeof(char));	
		result = fread (&size, sizeof(size_t), 1, pFile);
		params->sizes[i] = size;
	}
  fclose (pFile);
}

void print_file(char* path)
{
	FILE* fptr;
	fptr = fopen(path, "r");
	if(fptr==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file\n");
		exit(1);
	}
	char buffer[BUF_SIZE];
	while (fgets(buffer, BUF_SIZE, fptr) != 0)
	{
		if(strlen(buffer))
		{
			if(buffer[0]!='#')
			{
				fprintf(stderr, "\t%s", buffer);
			}
		}
	}
	fclose(fptr);
}

void write_nCDS(PARAMS* params)
{
	/* same as CDS but with inverted strands */
	unsigned long int number = 0;
	char character = 0x00;
	char string_32[SQID] = { 0x00 };
	char string_56[LTAG] = { 0x00 };
	int result = 0;
	FILE* fptr;
	fptr = fopen(params->genomic_bin, "rb");
	if(fptr == NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (argv[1])\n");
		exit(1);
	}

	FILE* pFile;
	pFile=fopen(params->nCDS_bin, "wb");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (%s)\n", params->nCDS_bin);
	}


	fprintf(stderr, "Writing non coding sequences binary\n");
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		fwrite(&string_32[0], SQID * sizeof(char), 1, pFile);

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		fwrite(&string_32[0], SQID * sizeof(char), 1, pFile );

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		fwrite(&number, sizeof(unsigned long int), 1, pFile );

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		fwrite(&number, sizeof(unsigned long int), 1, pFile );

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		fwrite(&number, sizeof(unsigned long int), 1, pFile );

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		if(character=='+')
		{
			character = '-';
			fwrite(&character, sizeof(char), 1, pFile );
		}
		else if(character=='-')
		{
			character = '+';
			fwrite(&character, sizeof(char), 1, pFile );
		}
		else fwrite(&character, sizeof(char), 1, pFile );
		

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		fwrite(&string_56[0], LTAG * sizeof(char), 1, pFile );
	}
	fclose(fptr);
	fclose(pFile);
}

void append_linear_forward(PARAMS* params, size_t i)
{
	// i is the index for the current replicon being analyzed

	// open output file for appending
	// open input file for reading
	// check number of features
	// detect end of last feature
	// treat end of feature, end of sequence and start of first feature in strand blind mode
	// limit up to first N=params->upstream_range nucleotides
	// write locus tags of features downstream to the calculated region
	// close files


	char replicon[SQID] = { 0x00 };
	strcpy(replicon, params->headers[i]);

	unsigned long int number = 0;
	char character = 0x00;
	char string_32[SQID] = { 0x00 };
	char string_56[LTAG] = { 0x00 };
	int result = 0;
	FILE* fptr;
	fptr = fopen(params->genomic_bin, "rb");
	if(fptr == NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (argv[1])\n");
		exit(1);
	}

	FILE* pFile;
	pFile=fopen(params->upstream_bin, "a+b");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (%s)\n", params->upstream_bin);
		exit(1);
	}

	// treat forward strand first
	unsigned long int latest_end = 0;
	unsigned long int first_start = 0;
	unsigned long int tmp_int = 0;
	int features_count = 0;
	int repl_on = 0;
	int boolean = 0;
	
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		// if is the right replicon switch is on
		if(strcmp(string_32, replicon)==0)
		{
			 repl_on = 1;
		}
		else  repl_on = 0;

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on) tmp_int = number;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on && character == '+')
		{
			latest_end = tmp_int;
		}
		if(character == '+') ++features_count;

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
	}

	if(!features_count) fprintf(stderr, "[WARNING] no features on forward strand found for %s\n", replicon);
	else
	{
		// rewind file pointer
		rewind(fptr);
		repl_on = 0;
		// find first forward strand feature
		while(1)
		{
			result = fread(&string_32, SQID * sizeof(char), 1, fptr);
			if(!result) break;
			// if is the right replicon switch is on
			if(strcmp(string_32, replicon)==0)
			{
				repl_on = 1;
			}
			else repl_on = 0;

			result = fread(&string_32, SQID * sizeof(char), 1, fptr );
			if(!result) break;
			
			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			if(repl_on) tmp_int = number;

			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;

			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;

			result = fread(&character, sizeof(char), 1, fptr );
			if(!result) break;
			if(repl_on && character == '+')
			{
				first_start = tmp_int;
				boolean = 1;
			}

			result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
			if(!result) break;
			if(boolean)
			{	
				character = '+';
				++latest_end;	// end + 1
				--first_start;	// start - 1
				unsigned long int uno = 1;

				// append upstream after begin of file
				fwrite(&replicon[0], SQID * sizeof(char), 1, pFile);
				fwrite(&string_32[0], SQID * sizeof(char), 1, pFile);
				fwrite(&uno, sizeof(unsigned long int), 1, pFile );
				fwrite(&first_start, sizeof(unsigned long int), 1, pFile );
				fwrite(&number, sizeof(unsigned long int), 1, pFile );
				fwrite(&character, sizeof(char), 1, pFile );
				fwrite(&string_56[0], LTAG * sizeof(char), 1, pFile );
				break;
			}
		}
		// track next features
		rewind(fptr);
		repl_on = 0;

		GFFline line1;
		GFFline line2;
		memset(line1.seqid, 0x00, SQID);
		memset(line1.feature, 0x00, SQID);
		memset(line1.locus_tag, 0x00, LTAG);
		line1.start = 0;
		line1.end = 0;
		line1.mscore = 0;
		line1.strand = 0x00;
		line2.start = 0;
		line2.end = 0;
		line2.mscore = 0;
		line2.strand = 0x00;
		memset(line2.seqid, 0x00, SQID);
		memset(line2.feature, 0x00, SQID);
		memset(line2.locus_tag, 0x00, LTAG);

		int skip = 0;
		boolean = 0;
		while(1)
		{
			result = fread(&string_32, SQID * sizeof(char), 1, fptr);
			if(!result) break;
			// if is the right replicon switch is on
			if(strcmp(string_32, replicon)==0)
			{
				 repl_on = 1;
			}
			else repl_on = 0;

			result = fread(&string_32, SQID * sizeof(char), 1, fptr );
			if(!result) break;
			if(repl_on && !skip) strcpy(line1.feature, string_32);
			if(repl_on && skip) strcpy(line2.feature, string_32);

			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			if(repl_on && !skip) line1.start = number;
			if(repl_on && skip) line2.start = number;

			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			if(repl_on && !skip) line1.end = number;
			if(repl_on && skip) line2.end = number;	

			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			if(repl_on && !skip) line1.mscore = number;
			if(repl_on && skip) line2.mscore = number;

			result = fread(&character, sizeof(char), 1, fptr );
			if(!result) break;
			if(repl_on && character == '+')
			{
				boolean = 1;
			}
			else if(repl_on && character == '-')
			{
				boolean = 0;
			}
			result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
			if(!result) break;
			if(repl_on && skip && boolean)
			{
				strcpy(line2.locus_tag, string_56);
				++line1.end;
				--line2.start;
				if(line2.start > line1.end)
				{
					// append output file
					fwrite(&replicon[0], SQID * sizeof(char), 1, pFile);
					fwrite(&line2.feature[0], SQID * sizeof(char), 1, pFile);
					fwrite(&line1.end, sizeof(unsigned long int), 1, pFile );
					fwrite(&line2.start, sizeof(unsigned long int), 1, pFile );
					fwrite(&line2.mscore, sizeof(unsigned long int), 1, pFile );
					fwrite(&character, sizeof(char), 1, pFile );
					fwrite(&line2.locus_tag[0], LTAG * sizeof(char), 1, pFile );					
				}
				--line1.end;	// reset to original
				++line2.start;	// reset to original

				// swap
				line1.start = line2.start;
				line1.end = line2.end;
				line1.mscore = line2.mscore;
				strcpy(line1.feature,line2.feature);
				strcpy(line1.locus_tag,line2.locus_tag);
			}
			if(repl_on && !skip && boolean)
			{
				strcpy(line1.locus_tag, string_56);
				skip = 1;
			}
		}
	}

	fclose(fptr);
	fclose(pFile);
}

void append_linear_reverse(PARAMS* params, size_t i)
{
	// i is the index for the current replicon being analyzed

	// open output file for appending
	// open input file for reading
	// check number of features
	// treat end of feature, end of sequence and start of first feature in strand blind mode
	// limit up to first N=params->upstream_range nucleotides
	// write locus tags of features downstream to the calculated region
	// detect end of last feature
	// close files


	char replicon[SQID] = { 0x00 };
	strcpy(replicon, params->headers[i]);

	unsigned long int number = 0;
	char character = 0x00;
	char string_32[SQID] = { 0x00 };
	char string_56[LTAG] = { 0x00 };
	int result = 0;
	FILE* fptr;
	fptr = fopen(params->genomic_bin, "rb");
	if(fptr == NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (argv[1])\n");
		exit(1);
	}

	FILE* pFile;
	pFile=fopen(params->upstream_bin, "a+b");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (%s)\n", params->upstream_bin);
		exit(1);
	}

	unsigned long int first_start = 0;
	int features_count = 0;
	int repl_on = 0;
	int boolean = 0;

	GFFline line1;
	GFFline line2;
	memset(line1.seqid, 0x00, SQID);
	memset(line1.feature, 0x00, SQID);
	memset(line1.locus_tag, 0x00, LTAG);
	line1.start = 0;
	line1.end = 0;
	line1.mscore = 0;
	line1.strand = 0x00;
	line2.start = 0;
	line2.end = 0;
	line2.mscore = 0;
	line2.strand = 0x00;
	memset(line2.seqid, 0x00, SQID);
	memset(line2.feature, 0x00, SQID);
	memset(line2.locus_tag, 0x00, LTAG);

	int skip = 0;
	boolean = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		// if is the right replicon switch is on
		if(strcmp(string_32, replicon)==0)
		{
			 repl_on = 1;
		}
		else repl_on = 0;

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) strcpy(line1.feature, string_32);
		if(repl_on && skip) strcpy(line2.feature, string_32);

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) line1.start = number;
		if(repl_on && skip) line2.start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) line1.end = number;
		if(repl_on && skip) line2.end = number;	

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) line1.mscore = number;
		if(repl_on && skip) line2.mscore = number;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on && character == '-')
		{
			boolean = 1;
		}
		else if(repl_on && character == '+')
		{
			boolean = 0;
		}
		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on && skip && boolean)
		{
			strcpy(line2.locus_tag, string_56);
			++line1.end;
			--line2.start;
			if(line2.start > line1.end)
			{
				// append output file
				fwrite(&replicon[0], SQID * sizeof(char), 1, pFile);
				fwrite(&line1.feature[0], SQID * sizeof(char), 1, pFile);
				fwrite(&line1.end, sizeof(unsigned long int), 1, pFile);
				fwrite(&line2.start, sizeof(unsigned long int), 1, pFile);
				fwrite(&line1.mscore, sizeof(unsigned long int), 1, pFile);
				fwrite(&character, sizeof(char), 1, pFile);
				fwrite(&line1.locus_tag[0], LTAG * sizeof(char), 1, pFile);			
			}
			--line1.end;	// reset to original
			++line2.start;	// reset to original

			// swap
			line1.start = line2.start;
			line1.end = line2.end;
			line1.mscore = line2.mscore;
			strcpy(line1.feature,line2.feature);
			strcpy(line1.locus_tag,line2.locus_tag);
		}
		if(repl_on && !skip && boolean)
		{
			strcpy(line1.locus_tag, string_56);
			skip = 1;
		}
	}

	// check the status of the last gene
	rewind(fptr);

	int SOF = 1;
	first_start = 0;
	line1.start = 0;
	line2.end = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		// if is the right replicon switch is on
		if(strcmp(string_32, replicon)==0)
		{
			 repl_on = 1;
		}
		else  repl_on = 0;

		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr);
		if(!result) break;
		if(repl_on) first_start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr);
		if(!result) break;
		if(repl_on) line2.end = number;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr);
		if(!result) break;
		

		result = fread(&character, sizeof(char), 1, fptr);
		if(!result) break;
		if(repl_on && character == '-')
		{
			if(SOF)
			{
				line1.start = first_start;
			}
			++features_count;
		}

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr);
		if(!result) break;
		if(line1.start!=0 && repl_on && SOF)
		{
			strcpy(line1.locus_tag, string_56);
			SOF = 0;
		}
		if(repl_on) strcpy(line2.locus_tag, string_56);
	}

	// append output file
	++line2.end;
	--line1.start;
	if(!features_count) fprintf(stderr, "[WARNING] no features on reverse strand found for %s\n", replicon);
	else
	{
		fwrite(&replicon[0], SQID * sizeof(char), 1, pFile);
		fwrite(&line2.feature[0], SQID * sizeof(char), 1, pFile);
		fwrite(&line2.end, sizeof(unsigned long int), 1, pFile);
		fwrite(&params->sizes[i], sizeof(unsigned long int), 1, pFile);
		fwrite(&line2.mscore, sizeof(unsigned long int), 1, pFile);
		fwrite(&character, sizeof(char), 1, pFile);
		fwrite(&line2.locus_tag[0], LTAG * sizeof(char), 1, pFile);		
	}

	fclose(fptr);
	fclose(pFile);
}

void write_true_intergenic(PARAMS* params, size_t i)
{
	char replicon[SQID] = { 0x00 };
	strcpy(replicon, params->headers[i]);

	unsigned long int number = 0;
	char character = 0x00;
	char string_32[SQID] = { 0x00 };
	char string_56[LTAG] = { 0x00 };
	int result = 0;
	FILE* fptr;
	fptr = fopen(params->genomic_bin, "rb");
	if(fptr == NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (argv[1])\n");
		exit(1);
	}

	FILE* pFile;
	pFile=fopen(params->true_inter_bin, "ab+");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (%s)\n", params->upstream_bin);
	}

	// treat forward strand first
	unsigned long int latest_end = 0;
	unsigned long int first_start = 0;
	unsigned long int tmp_int = 0;
	int features_count = 0;
	int repl_on = 0;
	int boolean = 0;

	GFFline line1;
	GFFline line2;
	memset(line1.seqid, 0x00, SQID);
	memset(line1.feature, 0x00, SQID);
	memset(line1.locus_tag, 0x00, LTAG);
	line1.start = 0;
	line1.end = 0;
	line1.mscore = 0;
	line1.strand = '.';
	line2.start = 0;
	line2.end = 0;
	line2.mscore = 0;
	line2.strand = '.';
	memset(line2.seqid, 0x00, SQID);
	memset(line2.feature, 0x00, SQID);
	memset(line2.locus_tag, 0x00, LTAG);
	int skip = 0;

	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		// if is the right replicon switch is on
		if(strcmp(string_32, replicon)==0)
		{
			repl_on = 1;
		}
		else repl_on = 0;

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		
		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on) tmp_int = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on)
		{
			first_start = tmp_int;
			boolean = 1;
		}

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		if(boolean)
		{	
			++latest_end;	// end + 1
			--first_start;	// start - 1
			unsigned long int uno = 1;

			// append upstream after begin of file
			fwrite(&replicon, SQID * sizeof(char), 1, pFile);
			fwrite(&string_32[0], SQID * sizeof(char), 1, pFile);
			fwrite(&uno, sizeof(unsigned long int), 1, pFile );
			fwrite(&first_start, sizeof(unsigned long int), 1, pFile );
			fwrite(&number, sizeof(unsigned long int), 1, pFile );
			fwrite(&character, sizeof(char), 1, pFile );
			fwrite(&string_56[0], LTAG * sizeof(char), 1, pFile );
			break;
		}
	}
	repl_on = 0;
	rewind(fptr);


	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		// if is the right replicon switch is on
		if(strcmp(string_32, replicon)==0)
		{
			 repl_on = 1;
		}
		else repl_on = 0;

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) strcpy(line1.feature, string_32);
		if(repl_on && skip) strcpy(line2.feature, string_32);

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) line1.start = number;
		if(repl_on && skip) line2.start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) line1.end = number;
		if(repl_on && skip) line2.end = number;	

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on && !skip) line1.mscore = number;
		if(repl_on && skip) line2.mscore = number;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on && skip)
		{
			strcpy(line2.locus_tag, string_56);
			++line1.end;
			--line2.start;
			if(line2.start > line1.end)
			{
				// append output file
				fwrite(&replicon, SQID * sizeof(char), 1, pFile);
				fwrite(&line2.feature[0], SQID * sizeof(char), 1, pFile);
				fwrite(&line1.end, sizeof(unsigned long int), 1, pFile );
				fwrite(&line2.start, sizeof(unsigned long int), 1, pFile );
				fwrite(&line2.mscore, sizeof(unsigned long int), 1, pFile );
				fwrite(&character, sizeof(char), 1, pFile );
				fwrite(&line2.locus_tag[0], LTAG * sizeof(char), 1, pFile );					
			}
			--line1.end;	// reset to original
			++line2.start;	// reset to original

			// swap
			line1.start = line2.start;
			line1.end = line2.end;
			line1.mscore = line2.mscore;
			strcpy(line1.feature,line2.feature);
			strcpy(line1.locus_tag,line2.locus_tag);
		}
		if(repl_on && !skip)
		{
			strcpy(line1.locus_tag, string_56);
			skip = 1;
		}
	}

	// check the status of the last feature
	rewind(fptr);

	int SOF = 1; 
	tmp_int = 0;
	latest_end = 0;
	first_start = 0;
	line1.start = 0;
	line2.end = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr);
		if(!result) break;
		// if is the right replicon switch is on
		if(strcmp(string_32, replicon)==0)
		{
			 repl_on = 1;
		}
		else  repl_on = 0;

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on) first_start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		if(repl_on) line2.end = number;
		

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		if(repl_on)
		{
			if(SOF)
			{
				line1.start = first_start;
			}
			++features_count;
		}

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		if(line1.start!=0 && repl_on && SOF)
		{
			strcpy(line1.locus_tag, string_56);
			SOF = 0;
		}
		if(repl_on) strcpy(line2.locus_tag, string_56);
	}

	// append output file
	++line2.end;
	--line1.start;
	if(!features_count) fprintf(stderr, "[WARNING] no features found for %s\n", replicon);
	else
	{
		fwrite(&replicon, SQID * sizeof(char), 1, pFile);
		fwrite(&line2.feature[0], SQID * sizeof(char), 1, pFile);
		fwrite(&line2.end, sizeof(unsigned long int), 1, pFile );
		fwrite(&params->sizes[i], sizeof(unsigned long int), 1, pFile );
		fwrite(&line2.mscore, sizeof(unsigned long int), 1, pFile );
		fwrite(&character, sizeof(char), 1, pFile );
		fwrite(&line2.locus_tag[0], LTAG * sizeof(char), 1, pFile );		
	}
	fclose(fptr);
	fclose(pFile);
}	

int main(int argc, char** argv)
{
	if(argc!=2)
	{
		fprintf(stderr, "Usage: %s params.ms (i.e., binary file produced by MeStudio pipeline)\n", argv[0]);
		exit(1);
	}
	PARAMS params;
	read_params_bin(&params, argv[1]);

	size_t i = 0;


	//core of the program
	write_nCDS(&params);


	fprintf(stderr, "Writing upstream sequences binary\n");

	//initialize file
	for(i=0;i<params.headers_count;++i)
	{
		append_linear_forward(&params, i);
		append_linear_reverse(&params, i);
	}

	fprintf(stderr, "Writing intergenic sequences binary\n");
	for(i=0;i<params.headers_count;++i)
	{
		write_true_intergenic(&params, i);
	}


	free(params.genomic_gff);
	free(params.sequencer_gff);
	free(params.genomic_fasta);
	free(params.motifs_list_txt);
	free(params.output_dir);
	free(params.genomic_bin);
	free(params.sequencer_bin);
	free(params.matches_bin);
	free(params.genomic_fasta_bin);
	free(params.nCDS_bin);
	free(params.true_inter_bin);
	free(params.upstream_bin);
	for(i=0;i<params.headers_count;++i)
	{
		free(params.headers[i]);
	}
	free(params.headers);
	free(params.sizes);

	return 0;
}