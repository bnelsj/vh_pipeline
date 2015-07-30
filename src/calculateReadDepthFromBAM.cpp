#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int maxNumSamples=160;
const int maxNumSV=10000000;
typedef struct SV{
	char chroName[10];
	int posOuterLeft;
	int posInnerLeft; 
	int posInnerRight;
	int posOuterRight;

	char SVName[20];
//	int sampleSup[maxNumSamples];
//	float sampleEditDist[maxNumSamples];
//	float averageEditDist;
	int totalSup;
	int minLen;
	int maxLen;
	bool yetGood;
}SV;

SV listSV[maxNumSV];
int numSV;

double expectedReadDepthGivenGC[101]; // the expected readdepth given GC for that individuals
int GenomeGC[25][300000000]; // for each [chro][pos] the window [pos-100, pos+200] the value (int) GC*100 sorted

int arrayReadDepth[25][300000000];

int GenomeLength[25];


int findChroArray(char *chroName)
{
	if (strcmp(chroName,"1")==0 || strcmp(chroName, "chr1")==0)
		return 1;
	if (strcmp(chroName,"2")==0 || strcmp(chroName, "chr2")==0)
		return 2;
	if (strcmp(chroName,"3")==0 || strcmp(chroName, "chr3")==0)
		return 3;
	if (strcmp(chroName,"4")==0 || strcmp(chroName, "chr4")==0)
		return 4;
	if (strcmp(chroName,"5")==0 || strcmp(chroName, "chr5")==0)
		return 5;
	if (strcmp(chroName,"6")==0 || strcmp(chroName, "chr6")==0)
		return 6;
	if (strcmp(chroName,"7")==0 || strcmp(chroName, "chr7")==0)
		return 7;
	if (strcmp(chroName,"8")==0 || strcmp(chroName, "chr8")==0)
		return 8;
	if (strcmp(chroName,"9")==0 || strcmp(chroName, "chr9")==0)
		return 9;
	if (strcmp(chroName,"10")==0 || strcmp(chroName, "chr10")==0)
		return 10;
	if (strcmp(chroName,"11")==0 || strcmp(chroName, "chr11")==0)
		return 11;
	if (strcmp(chroName,"12")==0 || strcmp(chroName, "chr12")==0)
		return 12;
	if (strcmp(chroName,"13")==0 || strcmp(chroName, "chr13")==0)
		return 13;
	if (strcmp(chroName,"14")==0 || strcmp(chroName, "chr14")==0)
		return 14;
	if (strcmp(chroName,"15")==0 || strcmp(chroName, "chr15")==0)
		return 15;
	if (strcmp(chroName,"16")==0 || strcmp(chroName, "chr16")==0)
		return 16;
	if (strcmp(chroName,"17")==0 || strcmp(chroName, "chr17")==0)
		return 17;
	if (strcmp(chroName,"18")==0 || strcmp(chroName, "chr18")==0)
		return 18;
	if (strcmp(chroName,"19")==0 || strcmp(chroName, "chr19")==0)
		return 19;
	if (strcmp(chroName,"20")==0 || strcmp(chroName, "chr20")==0)
		return 20;
	if (strcmp(chroName,"21")==0 || strcmp(chroName, "chr21")==0)
		return 21;
	if (strcmp(chroName,"22")==0 || strcmp(chroName, "chr22")==0)
		return 22;
	if (strcmp(chroName,"X")==0 || strcmp(chroName, "chrX")==0)
		return 23;
	if (strcmp(chroName,"Y")==0 || strcmp(chroName, "chrY")==0)
		return 24;

return 0;
}


int calculateMuStd(int SVid)
{

double mu=0;
double  std=0;
int chrId=0;
//chrId=1;
chrId=findChroArray(listSV[SVid].chroName);
int countLen=0;
	for (int count=listSV[SVid].posInnerLeft; count<listSV[SVid].posInnerRight; count++)	
	{
		//if (2*(double)arrayReadDepth[chrId][count]/(double)expectedReadDepthGivenGC[GenomeGC[chrId][count]]<11)
		{
			mu=mu+(2*(double)arrayReadDepth[chrId][count]/(double)expectedReadDepthGivenGC[GenomeGC[chrId][count]]);
			countLen++;
		}
	}
	
	mu=(double)mu/(double)countLen;
	
	for (int count=listSV[SVid].posInnerLeft; count<listSV[SVid].posInnerRight; count++)
	{
		//if (2*(double)arrayReadDepth[chrId][count]/(double)expectedReadDepthGivenGC[GenomeGC[chrId][count]]<11)
		{
			std=std+((2*(double)arrayReadDepth[chrId][count]/(double)expectedReadDepthGivenGC[GenomeGC[chrId][count]])-mu)*((2*(double)arrayReadDepth[chrId][count]/(double)expectedReadDepthGivenGC[GenomeGC[chrId][count]])-mu);
		}
	}
	std=(double)std/(double)countLen;
	//printf("%s\t%i\t%i\n", listSV[SVid].chroName, listSV[SVid].posInnerLeft, listSV[SVid].posInnerRight);
	printf("%s\t%i\t%i\tCNV:%lf\t%lf\t%lf\n",  listSV[SVid].chroName, listSV[SVid].posInnerLeft, listSV[SVid].posInnerRight, (double)mu, (double)sqrt(std), (double)sqrt(std)/(double)(sqrt(countLen)));

}


int initializeGenomeLen()
{
	char chroName[10];
	int len, id;
	FILE *fp=fopen("/net/eichler/vol23/projects/simons_genome_project/nobackups/fhormozd/VH_Initial_Calls/BAM_File_Depth/GC_Cal/ChroInfo","r");
	while(fscanf(fp,"%s\t%i\n", chroName, &len)!=EOF)
	{
		id=findChroArray(chroName);
		GenomeLength[id]=len;	
	}
}


int initializeReadDepth(FILE *fp)
{
char chroName[10];
int pos, count; 
int id=0;
long int tempCount[101];
	for (count=0; count<101; count++)
	{
		expectedReadDepthGivenGC[count]=0;
		tempCount[count]=0;
	}

	while(fscanf(fp, "%s\t%i\t%i\n", chroName, &pos, &count)!=EOF)
	{
		id = findChroArray(chroName);
		if (id>0)
		{
			arrayReadDepth[findChroArray(chroName)][pos]=count;
		}
		//if (pos%1000000==0)
		//	printf("%s\t%i\n", chroName, pos);
	}



	for (int count=1; count<23; count++)
	{
		for (int count2=0; count2<GenomeLength[count]; count2++)
		{
			if (arrayReadDepth[count][count2]>0)
			{
				expectedReadDepthGivenGC[GenomeGC[count][count2]]=expectedReadDepthGivenGC[GenomeGC[count][count2]]+arrayReadDepth[count][count2];
				tempCount[GenomeGC[count][count2]]++;
		//		printf("%i %i\n", GenomeGC[count][count2], tempCount[GenomeGC[count][count2]]);
			}
		}
	}	

	for (int count=0; count<101; count++)
	{
		expectedReadDepthGivenGC[count]=(double)expectedReadDepthGivenGC[count]/(double)tempCount[count];
//		printf("%i %f %i\n", count, expectedReadDepthGivenGC[count], tempCount[count]);
	}




}


int initializeGC(FILE *fp)
{
char chroName[10];
int pos1, pos2;
float gc;
int id;
while(fscanf(fp,"%s\t%i\t%i\t%f\n", chroName, &pos1, &pos2, &gc)!=EOF)
{
	id=findChroArray(chroName);
	if (id>0)
	{
		for (int count=pos1; count<pos2; count++)
		{
			GenomeGC[id][count]=GenomeGC[id][count]+(int)floorf(gc*100);	
			//printf("%i\n", GenomeGC[id][count]);
		} 
	//if (pos2%10000000==0)
	//	printf("%i\n", pos2);
	}
}

for (int count=1; count<25; count++)
{
	for (int count2=0; count2<300000000; count2++)
	{
		GenomeGC[count][count2]=GenomeGC[count][count2]/2;
	//	printf("%i %i\n", count2, GenomeGC[count][count2]);
	} 
}

}



int main( int argc, char *argv[]  )
{
	FILE *fp=fopen(argv[1],"r");


	initializeGenomeLen();	
	initializeGC(fp);
	
	
	FILE *fp2=fopen(argv[2],"r");

	initializeReadDepth(fp2);
	fclose(fp2);

//	FILE *fp3=fopen(argv[3],"r");
	FILE *fp4=fopen(argv[3],"r");	

/*	
	for (int count=0; count<160; count++)
	{
		fscanf(fp2,"%s\n", samplesName[count]);
	}
*/
	while(fscanf(fp4,"%s\t%i\t%i\t%s\n", listSV[numSV].chroName, &listSV[numSV].posInnerLeft, &listSV[numSV].posInnerRight, listSV[numSV].SVName)!=EOF)
	{
		numSV++;
	}

	for (int count=0; count<numSV; count++)
	{
		calculateMuStd(count);
	}

}
