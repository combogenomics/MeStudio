#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define BUF_SIZE 1024
#define SQID 32
#define LTAG 56
#include <errno.h>
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

void read_params_bin(PARAMS* params, char* path)
{
	FILE * pFile;
	long lSize;
	size_t result;
	pFile = fopen (path, "rb");
	if (pFile==NULL) 
	{
	fprintf(stderr, "[ERROR] Could not read file\n"); 
	exit (1);
	}
	// obtain file size
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);
	if(!lSize) exit(1);
	size_t i = 0;
 	size_t size = 0;
 	int booln = 0;
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
	result = fread (&booln, sizeof(int), 1, pFile);
	if(!result) return;
	params->upstream_range = booln;	

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

int main(int argc, char** argv)
{	
	if(argc!=2)
	{
		fprintf(stderr, "Usage: %s <.ms file> (i.e., any binary file produced by MeStudio pipeline)\n", argv[0]);
		exit(1);
	}
	if(strcmp(&(argv[1][strlen(argv[1]) - strlen("params.ms")]), "params.ms")==0)
	{
		PARAMS params;
		read_params_bin(&params, argv[1]);

		printf("Genomic GFF:\t\t%s\n", params.genomic_gff);
		printf("Sequencer GFF:\t\t%s\n", params.sequencer_gff);
		printf("Genomic FASTA:\t\t%s\n", params.genomic_fasta);
		printf("Motifs User List:\t%s\n", params.motifs_list_txt);
		printf("Output WorkDir:\t\t%s\n", params.output_dir);
		printf("Genomic GFF(bin):\t%s\n", params.genomic_bin);
		printf("Sequencer GFF(bin):\t%s\n", params.sequencer_bin);
		printf("Matches GFF(bin):\t%s\n", params.matches_bin);
		printf("Genomic FASTA(bin):\t%s\n", params.genomic_fasta_bin);
		printf("nCDS(bin):\t\t%s\n", params.nCDS_bin);
		printf("True Intergenic(bin):\t%s\n", params.true_inter_bin);
		printf("Upstream(bin):\t\t%s\n", params.upstream_bin);
		printf("Upstream range:\t\t%d\n", params.upstream_range);
		size_t i = 0;
		printf("\n");
		for(i=0;i<params.headers_count;++i)
		{
			printf("Seqid %s of size: %ld\n", params.headers[i], params.sizes[i]);
		}

		free(params.genomic_gff);
		free(params.sequencer_gff);
		free(params.genomic_fasta);
		free(params.output_dir);
		free(params.motifs_list_txt);
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

	}
	else if(strcmp(&(argv[1][strlen(argv[1]) - strlen("genomic_fasta.ms")]), "genomic_fasta.ms")==0)
	{
		FILE * fptr;;
		long lSize;
		size_t result;

		fptr = fopen (argv[1], "rb");
		if (fptr==NULL) 
		{
			fprintf(stderr, "[ERROR] Could not read genomic_fasta.ms\n"); 
			exit(1);
		}

		// obtain file size
		fseek (fptr , 0 , SEEK_END);
		lSize = ftell (fptr);
		rewind (fptr);

		if(!lSize) exit(1);

		size_t size = 0;
		size_t headers = 0;
		size_t nucleotides = 0;
		char string_32[SQID] = { 0x00 };

	 	/* read headers count */
		
		result = fread (&size, sizeof(size_t), 1, fptr);
		if(!result) exit(1);
		headers = size;
		if(!headers)
		{
			fprintf(stderr, "[ERROR] No seqids found\n");
			exit(1);
		}

		size_t i = 0;

		for(i=0;i<headers;++i)
		{
			long int j = 0;
			int s = 0;
			result = fread (&string_32, SQID * sizeof(char), 1, fptr);
			if(!result)
			{
				fprintf(stderr, "Binary file reading terminated\n");
				break;
			}	

			printf(">%s\n", string_32);
			result = fread (&size, sizeof(size_t), 1, fptr);
			if(!result)
			{
				fprintf(stderr, "Binary file reading terminated\n");
				break;
			}	
			nucleotides = size;

			char x;
			while(j < nucleotides)
			{
				result = fread (&x, sizeof(char), 1, fptr);
		        if(!result)
		        {
		            fprintf(stderr, "Binary file reading terminated\n");
		            break;
		        }
		        ++s;
		        if(s<=80) printf("%c", x);
		        else 
		        {
		        	printf("\n%c", x);
		        	s=1;
		        }
		        ++j;				
			}
			if(i<headers) printf("\n");
		}

		fclose(fptr);		
	}
	else if(strcmp(&(argv[1][strlen(argv[1]) - strlen("matches.ms")]), "matches.ms")==0)
	{
		// motif matches.ms binary file
		FILE * fptr;;
		long lSize;
		size_t result;

		fptr = fopen (argv[1], "rb");
		if (fptr==NULL) 
		{
			fprintf(stderr, "[ERROR] Could not read matches.ms\n"); 
			exit(1);
		}

		// obtain file size
		fseek (fptr , 0 , SEEK_END);
		lSize = ftell (fptr);
		rewind (fptr);

		if(!lSize) exit(1);
		size_t chars;
		
		while(1)
		{
			char* path;
			result = fread(&chars, sizeof(size_t), 1, fptr);
			if(!result) break;
			path = malloc((1 + chars) * sizeof(char));
			MEMORY_FAILURE(path);
			memset(path, 0x00, (1 + chars));
			result = fread(&path[0], chars * sizeof(char), 1, fptr);
			if(!result) break;
			printf("%s\n", path);
			free(path);	
		}
		fclose(fptr);
	}
	else	// generic pseudoGFF
	{
		unsigned long int number = 0;
		char character = 0x00;
		char string_32[SQID] = { 0x00 };
		char string_56[LTAG] = { 0x00 };
		int result = 0;
		FILE* fptr;
		fptr = fopen(argv[1], "rb");
		if(fptr == NULL)
		{
			fprintf(stderr, "[ERROR] Could not point to file (argv[1])\n");
			exit(1);
		}
		while(1)
		{
			result = fread(&string_32, SQID * sizeof(char), 1, fptr );
			if(!result) break;
			printf("%s\t", string_32);
			result = fread(&string_32, SQID * sizeof(char), 1, fptr );
			if(!result) break;
			printf("%s\t", string_32);
			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			printf("%ld\t", number);
			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			printf("%ld\t", number);
			result = fread(&number, sizeof(unsigned long int), 1, fptr );
			if(!result) break;
			printf("%ld\t", number);
			result = fread(&character, sizeof(char), 1, fptr );
			if(!result) break;
			printf("%c\t", character);
			result = fread(&string_56[0], LTAG * sizeof(char), 1, fptr );
			if(!result) break;
			printf("%s\n", string_56);
		}
	}
	return 0;
}