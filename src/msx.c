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
#include <limits.h>		// for ULONG_MAX
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
	int upstream_range;
	char* genomic_bin;
	char* sequencer_bin;
	char* params_bin;
	char* matches_bin;
	char* genomic_fasta_bin;
	char* nCDS_bin;
	char* true_inter_bin;
	char* upstream_bin;
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


void update_matches(char* motif_path, PARAMS* params)
{
	FILE* fptr;
	fptr = fopen(motif_path, "rb");
	if(fptr==NULL)
	{
		fprintf(stderr, "[ERROR]\n");
		return;
	}
	long line = SQID * sizeof(char) +
				SQID * sizeof(char) +
				3 * sizeof(unsigned long int) +
				sizeof(char) + 
				LTAG * sizeof(char);
	fseek(fptr, 0, SEEK_END);
	long match_size = ftell(fptr) / line;
	rewind(fptr);

	GFFline* gff_match = malloc(match_size * sizeof * gff_match);
	MEMORY_FAILURE(gff_match);
	/* initialize memory */
	size_t i = 0;
	for(i=0;i<match_size;++i)
	{
		memset(gff_match[i].seqid, 0x00, SQID * sizeof(char));
		memset(gff_match[i].feature, 0x00, SQID * sizeof(char));
		gff_match[i].start = 0;
		gff_match[i].end = 0;
		gff_match[i].mscore = 0;
		gff_match[i].strand = '.';
		memset(gff_match[i].locus_tag, 0x00, LTAG * sizeof(char));
	}

	i = 0;
	unsigned long int number = 0;
	char character = 0x00;
	char string_32[SQID] = { 0x00 };
	char string_56[LTAG] = { 0x00 };
	size_t result = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_match[i].seqid, string_32);

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_match[i].feature, string_32);

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_match[i].start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_match[i].end = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_match[i].mscore = number;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		gff_match[i].strand = character;

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_match[i].locus_tag, string_56);
		++i;
	}
	fclose(fptr);

	fptr = fopen(params->sequencer_bin, "rb");
	if(fptr==NULL)
	{
		fprintf(stderr, "[ERROR]\n");
		return;
	}
	fseek(fptr, 0, SEEK_END);
	long seq_size = ftell(fptr) / line;
	rewind(fptr);

	GFFline* gff_seq = malloc(seq_size * sizeof * gff_match);
	MEMORY_FAILURE(gff_seq);
	/* initialize memory */
	i = 0;
	for(i=0;i<seq_size;++i)
	{
		memset(gff_seq[i].seqid, 0x00, SQID * sizeof(char));
		memset(gff_seq[i].feature, 0x00, SQID * sizeof(char));
		gff_seq[i].start = 0;
		gff_seq[i].end = 0;
		gff_seq[i].mscore = 0;
		gff_seq[i].strand = '.';
		memset(gff_seq[i].locus_tag, 0x00, LTAG * sizeof(char));
	}

	i = 0;
	number = 0;
	character = 0x00;
	memset(string_32, 0x00, SQID * sizeof(char));
	memset(string_56, 0x00, LTAG * sizeof(char));
	result = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_seq[i].seqid, string_32);

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_seq[i].feature, string_32);

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_seq[i].start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_seq[i].end = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_seq[i].mscore = number;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		gff_seq[i].strand = character;

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_seq[i].locus_tag, string_56);
		++i;
	}
	fclose(fptr);

	/* motifs path gets updated. Mescore becomes the position of the methylated base. If not methylated mescore is 0 */
	FILE* pFile;
	pFile = fopen(motif_path, "wb");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR]\n");
		return;
	}

	size_t j = 0;
	int strcomp = 0;
	unsigned long int zero = 0;
	for(i=0;i<match_size;++i)
	{
		int isMe = 0;
		for(j=0;j<seq_size;++j)
		{
			strcomp = strcmp(gff_seq[j].seqid, gff_match[i].seqid);
			if(strcomp==0)
			{
				if(gff_seq[j].start >= gff_match[i].start &&
					gff_seq[j].end <= gff_match[i].end &&
					gff_seq[j].strand == gff_match[i].strand)
				{
					//printf("%s %s %ld %ld %ld %c %s\n", gff_match[i].seqid, gff_match[i].feature, gff_match[i].start, gff_match[i].end, gff_seq[j].start, gff_match[i].strand, gff_match[i].locus_tag);
					fwrite(&gff_match[i].seqid[0], SQID * sizeof(char), 1, pFile);
					fwrite(&gff_match[i].feature[0], SQID * sizeof(char), 1, pFile);
					fwrite(&gff_match[i].start, sizeof(unsigned long int), 1, pFile);
					fwrite(&gff_match[i].end, sizeof(unsigned long int), 1, pFile);
					fwrite(&gff_seq[j].start, sizeof(unsigned long int), 1, pFile);		// methylated base position
					fwrite(&gff_match[i].strand, sizeof(char), 1, pFile);
					fwrite(&gff_match[i].locus_tag[0], LTAG * sizeof(char), 1, pFile);
					isMe = 1;
				}
			}
		}
		if(!isMe)
		{
			//printf("%s %s %ld %ld %ld %c %s\n", gff_match[i].seqid, gff_match[i].feature, gff_match[i].start, gff_match[i].end, 0, gff_match[i].strand, gff_match[i].locus_tag);
			fwrite(&gff_match[i].seqid[0], SQID * sizeof(char), 1, pFile);
			fwrite(&gff_match[i].feature[0], SQID * sizeof(char), 1, pFile);
			fwrite(&gff_match[i].start, sizeof(unsigned long int), 1, pFile);
			fwrite(&gff_match[i].end, sizeof(unsigned long int), 1, pFile);
			fwrite(&zero, sizeof(unsigned long int), 1, pFile);		// methylated base position
			fwrite(&gff_match[i].strand, sizeof(char), 1, pFile);
			fwrite(&gff_match[i].locus_tag[0], LTAG * sizeof(char), 1, pFile);
		}
	}

	free(gff_match);
	free(gff_seq);
	fclose(pFile);
}


void cross(char* motif_path, char* feature_path, char* description, PARAMS* params)
{
	/* output file */
	char* output;	// like feature path, but having .gff instead of .ms (3 additional chars)
	size_t out_size = (strlen(motif_path) + strlen(description) + strlen(".gff")) * sizeof(char);	// reserve for 0x00
	output = malloc(out_size);
	MEMORY_FAILURE(output);
	memset(output, 0x00, out_size);
	strcpy(output, motif_path);
	output[strlen(output) - 3] = 0x00;
	strcat(output, "_");
	strcat(output, description);
	strcat(output, ".gff");

	FILE* pFile;
	pFile = fopen(output, "w");
	if(pFile==NULL)
	{
		fprintf(stderr, "[ERROR]\n");
		return;
	}

	FILE* fptr;
	fptr = fopen(motif_path, "rb");
	if(fptr==NULL)
	{
		fprintf(stderr, "[ERROR]\n");
		return;
	}
	long line = SQID * sizeof(char) +
				SQID * sizeof(char) +
				3 * sizeof(unsigned long int) +
				sizeof(char) + 
				LTAG * sizeof(char);
	fseek(fptr, 0, SEEK_END);
	long match_size = ftell(fptr) / line;
	rewind(fptr);

	GFFline* gff_match = malloc(match_size * sizeof * gff_match);
	MEMORY_FAILURE(gff_match);
	/* initialize memory */
	size_t i = 0;
	for(i=0;i<match_size;++i)
	{
		memset(gff_match[i].seqid, 0x00, SQID * sizeof(char));
		memset(gff_match[i].feature, 0x00, SQID * sizeof(char));
		gff_match[i].start = 0;
		gff_match[i].end = 0;
		gff_match[i].mscore = 0;
		gff_match[i].strand = '.';
		memset(gff_match[i].locus_tag, 0x00, LTAG * sizeof(char));
	}

	i = 0;
	unsigned long int number = 0;
	char character = 0x00;
	char string_32[SQID] = { 0x00 };
	char string_56[LTAG] = { 0x00 };
	size_t result = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_match[i].seqid, string_32);

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_match[i].feature, string_32);

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_match[i].start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_match[i].end = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_match[i].mscore = number;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		gff_match[i].strand = character;

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_match[i].locus_tag, string_56);
		++i;
	}
	fclose(fptr);


	fptr = fopen(feature_path, "rb");
	if(fptr==NULL)
	{
		fprintf(stderr, "[ERROR]\n");
		return;
	}
	fseek(fptr, 0, SEEK_END);
	long seq_size = ftell(fptr) / line;
	rewind(fptr);

	GFFline* gff_seq = malloc(seq_size * sizeof * gff_match);
	MEMORY_FAILURE(gff_seq);
	/* initialize memory */
	i = 0;
	for(i=0;i<seq_size;++i)
	{
		memset(gff_seq[i].seqid, 0x00, SQID * sizeof(char));
		memset(gff_seq[i].feature, 0x00, SQID * sizeof(char));
		gff_seq[i].start = 0;
		gff_seq[i].end = 0;
		gff_seq[i].mscore = 0;
		gff_seq[i].strand = '.';
		memset(gff_seq[i].locus_tag, 0x00, LTAG * sizeof(char));
	}

	i = 0;
	number = 0;
	character = 0x00;
	memset(string_32, 0x00, SQID * sizeof(char));
	memset(string_56, 0x00, LTAG * sizeof(char));
	result = 0;
	while(1)
	{
		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_seq[i].seqid, string_32);

		result = fread(&string_32, SQID * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_seq[i].feature, string_32);

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_seq[i].start = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_seq[i].end = number;

		result = fread(&number, sizeof(unsigned long int), 1, fptr );
		if(!result) break;
		gff_seq[i].mscore = number;

		result = fread(&character, sizeof(char), 1, fptr );
		if(!result) break;
		gff_seq[i].strand = character;

		result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
		if(!result) break;
		strcpy(gff_seq[i].locus_tag, string_56);
		++i;
	}
	fclose(fptr);

	/* if params->upstream_range is not defined (<0) set to max unsigned int */
	unsigned long int upstream = params->upstream_range;
	if(params->upstream_range < 0) 
	{
		upstream = ULONG_MAX;
	}

	/* motifs path gets updated. Mescore becomes the position of the methylated base. If not methylated mescore is 0 */
	size_t j = 0;
	int strcomp = 0;
	for(j=0;j<seq_size;++j)		//theoretically more populated than matches
	{
		for(i=0;i<match_size;++i)
		{
			strcomp = strcmp(gff_seq[j].seqid, gff_match[i].seqid);
			if(strcomp==0)
			{

				if(gff_match[i].start >= gff_seq[j].start &&
					gff_match[i].end <= gff_seq[j].end &&
					gff_seq[j].strand == gff_match[i].strand)
				{
					unsigned long int dist = 0;
					if(gff_match[i].mscore)
					{
						if(gff_seq[j].strand=='+')
						{
							dist = (gff_seq[j].end - gff_match[i].mscore) + 1;	// +1 bc is unsigned long int and cannot use -1 as default
						}						
						else if(gff_seq[j].strand=='-')
						{
							dist = (gff_match[i].mscore - gff_seq[j].start) + 1;	// +1 bc is unsigned long int and cannot use -1 as default
						}	
					}

					if(dist < upstream)
					{
						char c = '.';
						int g = gff_match[i].mscore - gff_match[i].start;

						if(g >= 0) c = gff_match[i].locus_tag[g];
						fprintf(pFile, 
							"%s\t%s\t%s\t%ld\t%ld\t%ld\t%c\t%c\t%s\n", 
							gff_seq[j].seqid,
							description, 
							gff_seq[j].feature, 
							gff_seq[j].start, 
							gff_seq[j].end, 
							dist, 
							gff_seq[j].strand,
							c,
							gff_seq[j].locus_tag);						
					}
				}
			}
		}
	}

	free(gff_match);
	free(gff_seq);
	free(output);
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

	FILE* fptr;
	fptr = fopen (params.matches_bin, "rb");
	if (fptr==NULL) 
	{
		fprintf(stderr, "[ERROR] Could not read matches.ms\n"); 
		exit(1);
	}

	size_t result = 0;
	size_t p_size = 0;
	while(1)
	{
		result = fread(&p_size, sizeof(size_t), 1, fptr);
		if(!result) break;
		char* motif_path = malloc((1 + p_size) * sizeof(char));
		MEMORY_FAILURE(motif_path);
		memset(motif_path, 0x00, (1 + p_size));
		result = fread(&motif_path[0], p_size * sizeof(char), 1, fptr);
		if(!result) break;


		fprintf(stderr, "Integrating methylation profiles for %s (will edit file)\n", motif_path);
		update_matches(motif_path, &params);


		fprintf(stderr, "Crossing %s\n", params.genomic_bin);
		cross(motif_path, params.genomic_bin, "CDS", &params);

		fprintf(stderr, "Crossing %s\n", params.nCDS_bin);
		cross(motif_path, params.nCDS_bin, "nCDS", &params);

		fprintf(stderr, "Crossing %s\n", params.true_inter_bin);
		cross(motif_path, params.true_inter_bin, "true_intergenic", &params);

		fprintf(stderr, "Crossing %s\n", params.upstream_bin);
		cross(motif_path, params.upstream_bin, "upstream", &params);

		free(motif_path);
	}

	fclose(fptr);

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