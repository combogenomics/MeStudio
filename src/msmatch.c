/*
===========================================================================================================================
* AUTHOR:		Christopher Riccardi, PhD student at Computational Biology Dept., University of Firenze
* PROJECT:		Extraction and analysis of methylation features from Pacific Biosciences SMRT reads using MeStudio
* CONTACT:		CHRISTOPHER.RICCARDI@UNIFI.IT
===========================================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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

int consensus(const char* s1, const char* s2)
{
	// compare and break as soon as at least
	// one match is found
	size_t j = 0;
	size_t i = 0;
	for(j=0; j <strlen(s1); ++j)
	{
		for(i=0; i<strlen(s2); ++i)
		{
			if(s1[j]==s2[i])
			{
				return 1;
				break;
			}
		}
	}
	return 0;
}

const char* hash(char c)
{
	//convert to upper and hash symbols
	switch(toupper(c))
	{
		case('A'):
		return "A";
		case('C'):
		return "C";
		case('G'):
		return "G";
		case('T'):
		return "TU";
		case('U'):
		return "UT";
		case('M'):
		return "AC";
		case('R'):
		return "AG";
		case('W'):
		return "ATU";
		case('S'):
		return "CG";
		case('Y'):
		return "CTU";
		case('K'):
		return "GTU";
		case('V'):
		return "ACG";
		case('H'):
		return "ACTU";
		case('D'):
		return "AGTU";
		case('B'):
		return "CGTU";
		case('N'):
		return "ACGTU";
		case('-'):
		return "-";
		default:
		return "-";
	}
}

int naiveMatch(char* pattern, char* string, const size_t size)
{
	/*	R: A or G
		N: any
		W: A or T
		Y: C or T
		B: T or G or C
		D: A or T or G	*/
	size_t i = 0;
	for (i=0; i < size; ++i)
	{
		//is O((n-m+1) m)
		if(!consensus(hash(pattern[i]),hash(string[i])) )
		{
			return 0;
		}
	}
	return 1;
}

void read_params_bin(PARAMS* params, char* path)
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
		//printf("Loaded %s : %ld\n", params->headers[i], params->sizes[i]);
	}



	//printf("Headers: %ld\n", params->headers_count);
	//char** headers;
	//size_t* sizes;


  // terminate
  fclose (pFile);
}

void revcmp(char* dest, char* src)
{
  char* d = &dest[0];
  char* p = &src[strlen(src) - 1];
  do
  {
    
    switch (*p)
    {
      case 'A':
      {
        *d = 'T';
      }break;
      case 'C':
      {
        *d = 'G';
      }break;
      case 'G':
      {
        *d = 'C';
      }break;
      case 'T':
      {
        *d = 'A';
      }break;
      default:
      {
        *d = *p;
      }
    }
    --p;
    ++d;
  }
  while(*p);
  *d = 0x00;
}

void cmpl(char* src)
{
  char* p = &src[strlen(src) - 1];
  do
  {
    switch (*p)
    {
      case 'A':
      {
        *p = 'T';
      }break;
      case 'C':
      {
        *p = 'G';
      }break;
      case 'G':
      {
        *p = 'C';
      }break;
      case 'T':
      {
        *p = 'A';
      }break;
      default:
      {
      }break;
    }
    --p;
  }
  while(*p);
}

void revert(char* dest, char* src)
{
	char* d = &dest[0];
	char* p = &src[strlen(src) - 1];
	do
	{
		*d = *p;
		--p;
		++d;
	}
	while(*p);
	*d = 0x00;
}

void msmatch(char* path, char* pattern, PARAMS* params)
{
	FILE* fptr;
	FILE* pFile;
	long lSize;
	size_t result;

	fptr = fopen (params->genomic_fasta_bin, "rb");
	if (fptr==NULL) 
	{
	fprintf(stderr, "[ERROR] Could not read genomic_fasta.ms\n"); 
	return;
	}

	pFile = fopen (path, "a+b");
	if (pFile==NULL) 
	{
	fprintf(stderr, "[ERROR] Could not open binary file\n"); 
	return;
	}

	// obtain file size
	fseek (fptr , 0 , SEEK_END);
	lSize = ftell (fptr);
	rewind (fptr);

	if(!lSize) return;

	size_t size = 0;
	size_t headers = 0;
	size_t nucleotides = 0;
	char seqid[SQID] = { 0x00 };
	char ref_str[SQID] = { 0x00 };
	char feature[SQID] = "match";

 	/* read headers count */
	result = fread (&size, sizeof(size_t), 1, fptr);
	if(!result) return;
	headers = size;
	if(!headers)
	{
		fprintf(stderr, "[ERROR] No seqids found\n");
		return;
	}
	size_t s = strlen(pattern);
	if(s>=LTAG)
	{
		fprintf(stderr, "[ERROR] Motif too long\n");
		exit(1);
	}
	char* tmprry = malloc((LTAG) * sizeof(char));
	MEMORY_FAILURE(tmprry);

	memset(tmprry, 0x00, (s + 1) * sizeof(char));
	unsigned long int mescore = 0;
	size_t i = 0;

	for(i=0;i<headers;++i)
	{
		long int j = 0;
		result = fread (&seqid, SQID * sizeof(char), 1, fptr);
		if(!result)
		{
			fprintf(stderr, "Binary file reading terminated\n");
			return;
		}	

		result = fread (&size, sizeof(size_t), 1, fptr);
		if(!result)
		{
			fprintf(stderr, "Binary file reading terminated\n");
			return;
		}	
		nucleotides = size;
		char* x = malloc((nucleotides + 1) * sizeof(char));
		MEMORY_FAILURE(x);
		memset(x, 0x00, (nucleotides + 1) * sizeof(char));
		result = fread (&x[0], nucleotides * sizeof(char), 1, fptr);
    if(!result)
    {
        fprintf(stderr, "Binary file reading terminated\n");
        return;
    }
		while(j < (nucleotides - s))
		{
			int pos = naiveMatch(pattern, &x[j], s);
			
			if(pos)
			{
				fwrite(&seqid[0], SQID * sizeof(char), 1, pFile);
				fwrite(&feature[0], SQID * sizeof(char), 1, pFile);
				unsigned long int start, end;
				start = 1 + (unsigned long int)j;
				end = start + (s - 1);
				fwrite(&start,  sizeof(unsigned long int), 1, pFile);
				fwrite(&end, sizeof(unsigned long int), 1, pFile);
				fwrite(&mescore, sizeof(unsigned long int), 1, pFile);
				char strand = '+';
				fwrite(&strand, sizeof(char), 1, pFile);
				strncpy(tmprry, &x[j], s);
				fwrite(&tmprry[0], LTAG * sizeof(char), 1, pFile);
				//printf("%s\t%s\t%ld\t%ld\t%ld\t%c\t%s\n", seqid, feature, start, end, mescore, strand, pattern);
			}
			++j;
		}

		// revert pattern, complement genome
		char rev_pattern[LTAG] = { 0x00 };
		revert(rev_pattern, pattern);
		fprintf(stderr, "Reverted pattern and complemented genome for double stranded matching\n");
		cmpl(x);


		j = 0;
		while(j < (nucleotides - s))
		{
			int pos = naiveMatch(rev_pattern, &x[j], s);
			
			if(pos)
			{
				fwrite(&seqid[0], SQID * sizeof(char), 1, pFile);
				fwrite(&feature[0], SQID * sizeof(char), 1, pFile);
				unsigned long int start, end;
				start = 1 + (unsigned long int)j;
				end = start + (s - 1);
				fwrite(&start,  sizeof(unsigned long int), 1, pFile);
				fwrite(&end, sizeof(unsigned long int), 1, pFile);
				fwrite(&mescore, sizeof(unsigned long int), 1, pFile );
				char strand = '-';
				fwrite(&strand, sizeof(char), 1, pFile );
				strncpy(tmprry, &x[j], s);
				fwrite(&tmprry[0], LTAG * sizeof(char), 1, pFile );
				start = 0;
				end = 0;
				//printf("%s\t%s\t%ld\t%ld\t%ld\t%c\t%s\n", seqid, feature, start, end, mescore, strand, pattern);
			}
			++j;
		}
		if(x) free(x);
		memset(seqid, 0x00, SQID * sizeof(char));	
		memset(ref_str, 0x00, SQID * sizeof(char));	
	}
	fclose(fptr);
	fclose(pFile);
	if(tmprry) free(tmprry);
}

void iterate_motifs(PARAMS* params)
{
	char buffer[LTAG] = {0x00};
		
	FILE* fptr;
	fptr = fopen(params->motifs_list_txt, "r");
	if (fptr == NULL)
	{
		printf("[ERROR] Could not open GFF3 file\n");
		return;
	}
	while (fgets(buffer, LTAG, fptr) != 0)
	{
		if(strlen(buffer) > 1)
		{
			if(buffer[0]!='#')
			{

				if(buffer[strlen(buffer) - 1] == '\n')
				{
					buffer[strlen(buffer) - 1] = 0x00;
				}
				if(buffer[strlen(buffer) - 1] == '\r')
				{
					buffer[strlen(buffer) - 1] = 0x00;
				}

				/* Update path with motif name */
				size_t motif_bin_len = (strlen(".ms") + 3 + strlen(params->output_dir) + (2 * strlen(buffer))) * sizeof(char);
				char* motif = malloc(motif_bin_len);
				MEMORY_FAILURE(motif);

				memset(motif, 0x00, motif_bin_len);
				strcpy(motif, params->output_dir);
				strcat(motif, "/");		// hence (strlen(".ms") + 3 + strlen(params->output_dir) + (2 * strlen(buffer))) * sizeof(char);
				strcat(motif, buffer);

				fprintf(stderr, "Creating directory %s\n", motif);
				makedir(motif);

				strcat(motif, "/");
				strcat(motif, buffer);
				strcat(motif, ".ms");
				fprintf(stderr, "Writing matches to %s\n", motif);

				fprintf(stderr, "Searching motif: %s\n", buffer);
				msmatch(motif, buffer, params);

				FILE* pFile;
				pFile = fopen (params->matches_bin, "ab");
				if (pFile==NULL) 
				{
					fprintf(stderr, "[ERROR] Could not open binary file\n"); 
					return;
				}

				size_t p_size = strlen(motif);
				fwrite(&p_size, sizeof(size_t), 1, pFile);
				fwrite(&motif[0], p_size * sizeof(char), 1, pFile);

				fclose(pFile);
				fprintf(stderr, "Matches were written to %s\n", motif);
				free(motif);
				memset(buffer, 0x00, LTAG *sizeof(char));
			}
		}
	}
	fclose(fptr);
}

int main(int argc, char* argv[])
{	
	if(argc!=2)
	{
		fprintf(stderr, "Usage: %s params.ms (i.e., binary file produced by MeStudio pipeline)\n", argv[0]);
		exit(1);
	}
	PARAMS params;
	read_params_bin(&params, argv[1]);
	iterate_motifs(&params);


	/* give memory back to the OS */
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
