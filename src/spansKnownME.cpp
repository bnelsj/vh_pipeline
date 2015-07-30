#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int maxNumSamples=160;

typedef struct SV{
	char chroName[20];
	int posOuterLeft;
	int posInnerLeft; 
	int posInnerRight;
	int posOuterRight;

	int sampleSup[maxNumSamples];
	float sampleEditDist[maxNumSamples];
	float averageEditDist;
	int totalSup;
	int minLen;
	int maxLen;
	bool yetGood;
	char SVName[20];
}SV;

SV SVTemp;

typedef struct MEI{

	char chroName[20];
	int pos1, pos2;

}MEI;
 

MEI listMEI[1000000];
int MEICount=0;

int main(int argv, char *argc[])
{
	FILE *fp=fopen(argc[1],"r");
	FILE *fp2=fopen(argc[2],"r");

	while(fscanf(fp,"%s\t%i\t%i\n", listMEI[MEICount].chroName, &listMEI[MEICount].pos1, &listMEI[MEICount].pos2)!=EOF)
	{
		MEICount++;
	}
	
		

	while(fscanf(fp2,"%s\t%i\t%i\t%s\t%i\t%f", SVTemp.chroName, &SVTemp.posInnerLeft, &SVTemp.posInnerRight, SVTemp.SVName, &SVTemp.totalSup, &SVTemp.averageEditDist)!=EOF)
	{
		for (int count=0; count<maxNumSamples; count++)
		{
			fscanf(fp2,"\t%i\t%f", &SVTemp.sampleSup[count], &SVTemp.sampleEditDist[count]);
		}
		fscanf(fp2, "\n");

		for (int count2=0; count2<MEICount; count2++)
		{
			if (strcmp(SVTemp.chroName, listMEI[count2].chroName)==0)
			{
				if (SVTemp.posInnerLeft<listMEI[count2].pos1 && listMEI[count2].pos1-SVTemp.posInnerLeft<150)
				{
					if (SVTemp.posInnerRight>listMEI[count2].pos2 && SVTemp.posInnerRight - listMEI[count2].pos2<150)
					{
						if (SVTemp.totalSup>20)
						printf("%s\t%i\t%i\t%s\n", SVTemp.chroName, SVTemp.posInnerLeft, SVTemp.posInnerRight, SVTemp.SVName);
					}
				}
			}
		}
	
	}


}
