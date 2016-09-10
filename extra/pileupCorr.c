#include "stdio.h"
#include "string.h"
char chromo[100];
char previousChromo[1000];
char coverage[100];
int pos;
int previousPos;
char refBase[10000];
char field3[100000];
char field4[100000];
FILE *in;
FILE *out;

int main(int argc,char** argv)
{
	in = fopen(argv[1],"r");
	out = fopen(argv[2],"w");
	strcpy(previousChromo,"");
	previousPos = 1;
	while(!feof(in))
	{
		fscanf(in,"%s",chromo);
		fscanf(in,"%d",&pos);
		fscanf(in,"%s",refBase);
		fscanf(in,"%s",coverage);
		fscanf(in,"%s",field3);
		fscanf(in,"%s",field4);
		if (feof(in)) break;
		if( strcmp(chromo,previousChromo) == 0)
		{
			if(pos>previousPos+1)
			{
				previousPos++;
				while(previousPos<pos)
				{
					if(previousPos!=0) fprintf(out,"%s\t%d\tN	0	0	0\n",chromo,previousPos);
					previousPos++;
					if (feof(in)) break;
				}
				fprintf(out,"%s\t%d\t%s\t%s\t%s\t%s\n",chromo,pos,refBase,coverage,field3, field4);
			}
			else 
			{
				fprintf(out,"%s\t%d\t%s\t%s\t%s\t%s\n",chromo,pos,refBase,coverage,field3, field4);
				previousPos = pos;
			}
		}
		else
		{
			previousPos=0;
			strcpy(previousChromo,chromo);
			if (pos==previousPos)  fprintf(out,"%s\t%d\t%s\t%s\t%s\t%s\n",chromo,pos,refBase,coverage,field3, field4);
			else
			{
				while(previousPos<pos)
				{
					if(previousPos!=0) fprintf(out,"%s\t%d\tN	0	0	0\n",chromo,previousPos);
					previousPos++;
					if (feof(in)) break;
				}
				fprintf(out,"%s\t%d\t%s\t%s\t%s\t%s\n",chromo,pos,refBase,coverage,field3, field4);
			}
		}
	}
}
	
