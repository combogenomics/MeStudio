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

void read_params_bin(PARAMS* params, char* path)
{
	FILE* pFile;
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
 	int bools = 0;
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
	result = fread (&bools, sizeof(int), 1, pFile);
	if(!result) return;
	params->upstream_range = bools;	

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

void write_binary(PARAMS* params)
{
	/* I/O variables */

	FILE* pFile;
	FILE* fptr;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

	pFile = fopen(params->genomic_fasta_bin, "wb");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (binary)\n");
		return;
	}
	fptr = fopen(params->genomic_fasta, "r");
	if(fptr==NULL)
	{
		fprintf(stderr, "[ERROR] Could not point to file (fasta: %s)\n", params->genomic_fasta);
		return;
	}

 	// initialize binary 

 	fwrite(&params->headers_count, sizeof(size_t), 1, pFile );


	size_t i = 0;
	for(i=0;i<params->headers_count;++i)
	{	
		// update binary
		fwrite(&params->headers[i][0], SQID * sizeof(char), 1, pFile );
		fwrite(&params->sizes[i], sizeof(size_t), 1, pFile );
		int printing = 0;
	    while ((read = getline(&line, &len, fptr)) != -1)
	    {
	        if(strlen(line)>0)
	        {
	          if(line[0]=='>')
	          {
	          	if(strncmp(&line[1], params->headers[i], strlen(params->headers[i]))==0)
	          	{
	          		printing = 1;
	          	}
	          	else
	          	{
	          		printing = 0;
	          	}
	          }
	          else
	          {
	          	if(printing)
	          	{
	          		// write to binary
	          		
					//if (line[read - 1] == '\n') line[read - 1] = 0x00;
					//if(read - 2 >= 0)
					//{
					//	if (line[read - 2] == '\r') line[read - 1] = 0x00;						
					//}
					rstrip(line);
					fwrite(&line[0], strlen(line) * sizeof(char), 1, pFile );
	          	}
	          }
	        }
	    }
	    fseek(fptr, 0, SEEK_SET);


	}


	fclose(pFile);
	fclose(fptr);
    if (line)
    {
        free(line);
    }
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
	write_binary(&params);	

	size_t i = 0;
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
