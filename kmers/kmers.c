/* The MIT License

   Copyright (c) 2014 by Homolog.us <samanta@homolog.us>

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



#include <stdio.h>
#include <stdlib.h>


void count_kmer(int kmer)
{
	FILE* left;
	FILE* leftR;
	
	left=fopen("reads","r");
	leftR=fopen("reads-r","r");

	char Sequence[700];
	int iii=0, i, L;
	long int count[95369];
	__uint64_t v, vc;
	unsigned int x;
	__uint64_t mask=0;
	for(i=0; i<kmer; i++) { mask+=3; mask=mask<<2; }
	mask=mask>>2;;

	map<__uint64_t,int> kmer_seq;

	while(!feof(left))
	{
		fscanf(left,"%s\n",Sequence);

                v=0;
		L=strlen(Sequence);
		if(L>=kmer)
		{
                	for(i=0; i<kmer; i++)
                	{
                        	if(Sequence[i]=='A') { v+=0; }
                        	if(Sequence[i]=='T') { v+=3; }
                        	if(Sequence[i]=='C') { v+=1; }
                        	if(Sequence[i]=='G') { v+=2; }
                        	if(Sequence[i]=='N') { v+=0; }
                        	v=v<<2;
                	}
                	v=v>>2;
			//x=v%95369; count[x]++;
			kmer_seq[v]++;

                	for(i=kmer; i<L; i++)
                	{
                        	v=v<<2;
                        	v= v & mask;
                        	if(Sequence[i]=='A') { v+=0; }
                        	if(Sequence[i]=='T') { v+=3; }
                        	if(Sequence[i]=='C') { v+=1; }
                        	if(Sequence[i]=='G') { v+=2; }
                        	if(Sequence[i]=='N') { v+=0; }
				kmer_seq[v]++;
			//x=v%95369; count[x]++;
                	}
		}
	}

	while(!feof(leftR))
	{
		fscanf(leftR,"%s\n",Sequence);

                v=0;
		L=strlen(Sequence);
		if(L>=kmer)
		{
                	for(i=0; i<kmer; i++)
                	{
                        	if(Sequence[i]=='A') { v+=0; }
                        	if(Sequence[i]=='T') { v+=3; }
                        	if(Sequence[i]=='C') { v+=1; }
                        	if(Sequence[i]=='G') { v+=2; }
                        	if(Sequence[i]=='N') { v+=0; }
                        	v=v<<2;
                	}
                	v=v>>2;
			kmer_seq[v]++;
			//x=v%95369; count[x]++;

                	for(i=kmer; i<L; i++)
                	{
                        	v=v<<2;
                        	v= v & mask;
                        	if(Sequence[i]=='A') { v+=0; }
                        	if(Sequence[i]=='T') { v+=3; }
                        	if(Sequence[i]=='C') { v+=1; }
                        	if(Sequence[i]=='G') { v+=2; }
                        	if(Sequence[i]=='N') { v+=0; }
				kmer_seq[v]++;
			//x=v%95369; count[x]++;
                	}
		}
	}

	map<__uint64_t,int>::iterator curr, end;

       	for(curr=kmer_seq.begin(), end=kmer_seq.end(); curr !=end; curr++)
        {
               cout << (*curr).first << " " << (*curr).second << "\n";
        }

	//for(i=0; i<95369; i++) { printf("%i %i\n",i, count[i]); }
}


