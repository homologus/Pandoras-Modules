/* The MIT License

   Copyright (c) 2014 Artisan Bioinformatics <info@artisanbio.info>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/


#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//
//	computes reverse complement of a string
//	takes ~5 minues (including reading time)  for a 
//	compressed file	with 100M reads
//

void rev_comp(char *s, char *rev_s, int size)
{
	int i;

	for(i=0; i<size; i++)
	{
		if(s[i]=='A') { rev_s[size-i-1]='T'; }
		else
		{
			if(s[i]=='C') { rev_s[size-i-1]='G'; }
			else
			{
				if(s[i]=='G') { rev_s[size-i-1]='C'; }
				else
				{
					if(s[i]=='T') { rev_s[size-i-1]='A'; }
					else
					{
						if(s[i]=='N') { rev_s[size-i-1]='N'; }
					}
				}
			}
		}
	}
	rev_s[size]='\0';
}



//
//	Main program
//

int main(int argc, char *argv[])
{
	gzFile fp1;
	kseq_t *seq1;
	char rev_s1[7000];
	int N=63, M=8, L;
	int i,j, l, pos, begin, strand=1; 
	char empty[300]="ZZZZZZZZZZZZZZZZZZ"; empty[M]='\0';
	char min[300]="ZZZZZZZZZZZZZZZZZZ"; min[M]='\0';

	if (argc == 1) 
	{
		fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
		return 1;
	}

	fp1 = gzopen(argv[1], "r");
	seq1 = kseq_init(fp1);
	while ((l = kseq_read(seq1)) >= 0) 
	{
		rev_comp(seq1->seq.s, rev_s1, strlen(seq1->seq.s));
		printf("seq: %s\n", seq1->seq.s);
		//printf("rev: %s\n", rev_s);

		pos=-1; 	// position of minimizer
		begin=0;	// position of begin
		strncpy(min,empty,M); min[M]='\0';
		L=strlen(seq1->seq.s);

		// this loop runs over the entire string - minus K-mer size
		for(i=0; i<=L-N; i++)
		{
			//printf("begin %i %i %i\n", pos, begin, i);
			if(i > pos)
			{
				// prints the substring and moves on to the next one
				if(min[0]!='Z')
				{
					//printf("here %i %i\n",i, pos);
					//printf("here1 %.*s %s %i %i\n", i-begin+N-1, seq1->seq.s+begin, min, begin, i-begin+N-1);
					printf("%.*s %s\n", i-begin+N-1, seq1->seq.s+begin, min);
				}

				strncpy(min,empty,M); min[M]='\0';
				for(j=0; j<=N-M; j++)
				{
					//printf("0min %s\n",min);
					if( strncmp(seq1->seq.s+i+j,min,M) < 0 )
					{
						pos=i+j;
						strncpy(min,seq1->seq.s+i+j,M);
					}
					//printf("1min %s\n",min);
					if( strncmp(rev_s1+L-i-j-M,min,M) < 0 )
					{
						pos=i+j;
						strncpy(min,rev_s1+L-i-j-M,M);
					}
					//printf("min %s\n",min);
				}
				begin=i;
			}
			else
			{
				//printf("compare %.*s %.*s %.*s\n",M, seq1->seq.s+i+j-1, M, rev_s1+L-i-j+1-M, M, min);
				strand=1;
				if(strncmp(seq1->seq.s+i+j-1,rev_s1+L-i-j+1-M,M)>0)
				{
					strand=-1;
				}

				if((strand==1) && ( strncmp(seq1->seq.s+i+j-1,min,M) < 0 ))
				{
					//printf("here2a %.*s %s %i\n", i-begin+N-1, seq1->seq.s+begin, min, begin);
					printf("%.*s %s\n", i-begin+N-1, seq1->seq.s+begin, min);
					begin=i; pos=i+N-M+1;
					strncpy(min,seq1->seq.s+i+j-1,M);
				}

				if((strand==-1) && ( strncmp(rev_s1+L-i-j+1-M, min,M) < 0 ))
				{
					//printf("here2b %.*s %s %i\n", i-begin+N-1, seq1->seq.s+begin, min, begin);
					printf("%.*s %s\n", i-begin+N-1, seq1->seq.s+begin, min);
					begin=i; pos=i+N-M+1;
					strncpy(min,rev_s1+L-i-j+1-M,M);
				}
			}
		}
				// prints the final substring
				if(min[0]!='Z')
				{
					//printf("here4 %.*s %s\n", i-begin+N, seq1->seq.s+begin, min);
					printf("%.*s %s\n", i-begin+N, seq1->seq.s+begin, min);
				}
	}

	printf("return value: %d\n", l);
	kseq_destroy(seq1);
	gzclose(fp1);
	return 0;
}
