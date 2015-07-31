#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int maxNumSamples=1000;
const int maxNumSV=10000000;
int numSV;
int totalNumSamples;
int numberSampleInput;
typedef struct SV{
	char chroName[10];
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
}SV; 

char samplesName[maxNumSamples][100]; 


SV listSV[maxNumSV];

int halfInsertSize=300;


bool matchs(int SV_Id1, int SV_Id2)
{// over 50% reciprical overlapp in region 
 // the breakpoints from both side <250bp
 //
int maxLeftPos, minRightPos;
if (strcmp(listSV[SV_Id1].chroName, listSV[SV_Id2].chroName)!=0)
	return false;
if (abs(listSV[SV_Id1].posInnerLeft - listSV[SV_Id2].posInnerLeft)>halfInsertSize || abs(listSV[SV_Id1].posInnerRight - listSV[SV_Id2].posInnerRight)>halfInsertSize)
	return false;
if (listSV[SV_Id1].posInnerLeft>listSV[SV_Id2].posInnerLeft)
{
	maxLeftPos=listSV[SV_Id1].posInnerLeft;
}else
{
	maxLeftPos=listSV[SV_Id2].posInnerLeft;
}

if (listSV[SV_Id1].posInnerRight>listSV[SV_Id2].posInnerRight)
{
	minRightPos=listSV[SV_Id2].posInnerRight;
}else
{
	minRightPos=listSV[SV_Id1].posInnerRight;
}



if ((minRightPos-maxLeftPos) < 0.25*(listSV[SV_Id1].posInnerRight-listSV[SV_Id1].posInnerLeft))
	return false;
if ((minRightPos-maxLeftPos) < 0.25*(listSV[SV_Id2].posInnerRight-listSV[SV_Id2].posInnerLeft))
	return false;

return true;

}

int compare(const void *a, const void *b)
{
	return (((SV*)a)->posInnerLeft - ((SV*)b)->posInnerLeft);
}
/*
int binarySearch(int minId, int maxId, int value)
{
	//int minId=0;
	//int maxId=numSV;
	printf("%i %i\n", minId, maxId);
	if (listSV[minId].posInnerLeft==value)
	{
		printf("L82, %i \n", listSV[minId].posInnerLeft);
		while(minId>0 && listSV[minId].posInnerLeft==value)
			minId--;
		return minId;
	}
	if (listSV[maxId].posInnerLeft==value)
	{	printf("L88, %i %i\n", maxId, listSV[maxId].posInnerLeft);
		while(maxId>0 && listSV[maxId].posInnerLeft==value)
			maxId--;
		return maxId;
	}	
	if (abs(minId-maxId)<2)
	{
		if (minId>0)
			minId--;
		return minId;
	}

	int midId=(minId+maxId)/2;
	printf("mid Id %i\n", midId);
	if (listSV[midId].posInnerLeft<value)
	{
		printf("L107\n");
		minId=midId;
		return (binarySearch(minId, maxId, value));
	}else if (listSV[midId].posInnerLeft>value)
	{
		maxId=midId;
		return (binarySearch(minId, maxId, value));
	}else if (listSV[midId].posInnerLeft==value)
	{
		while(midId>0 && listSV[midId].posInnerLeft==value)
			midId--;
		return midId;
	} 
	
}
	
*/

int mergeSV(int SV1, int SV2)
{
	if (listSV[SV1].totalSup>listSV[SV2].totalSup)
	{
		for (int count=0; count<totalNumSamples; count++)
		{
			listSV[SV1].sampleEditDist[count]=(float)(listSV[SV1].sampleEditDist[count]*listSV[SV1].sampleSup[count] + listSV[SV2].sampleEditDist[count]*listSV[SV2].sampleSup[count])/(float)(listSV[SV1].sampleSup[count] + listSV[SV2].sampleSup[count]);
			listSV[SV1].sampleSup[count]=listSV[SV1].sampleSup[count]+listSV[SV2].sampleSup[count];	
		}	
		
			listSV[SV1].averageEditDist=(float)(listSV[SV1].averageEditDist*listSV[SV1].totalSup + listSV[SV2].averageEditDist*listSV[SV2].totalSup)/(listSV[SV1].totalSup + listSV[SV2].totalSup);
			listSV[SV1].totalSup=listSV[SV1].totalSup+listSV[SV2].totalSup; 		
			listSV[SV2].yetGood=false;
	}else if (listSV[SV1].totalSup<listSV[SV2].totalSup)
	{
		for (int count=0; count<totalNumSamples; count++)
		{
			listSV[SV2].sampleEditDist[count]=(float)(listSV[SV1].sampleEditDist[count]*listSV[SV1].sampleSup[count] + listSV[SV2].sampleEditDist[count]*listSV[SV2].sampleSup[count])/(float)(listSV[SV1].sampleSup[count] + listSV[SV2].sampleSup[count]);
			listSV[SV2].sampleSup[count]=listSV[SV1].sampleSup[count]+listSV[SV2].sampleSup[count];	
		}
			listSV[SV1].yetGood=false;
			listSV[SV2].averageEditDist=(float)(listSV[SV1].averageEditDist*listSV[SV1].totalSup + listSV[SV2].averageEditDist*listSV[SV2].totalSup)/(listSV[SV1].totalSup + listSV[SV2].totalSup);
			listSV[SV2].totalSup=listSV[SV1].totalSup+listSV[SV2].totalSup; 		
	}else if (listSV[SV1].averageEditDist<listSV[SV2].averageEditDist)
	{
		for (int count=0; count<totalNumSamples; count++)
		{
			
			listSV[SV1].sampleEditDist[count]=(float)(listSV[SV1].sampleEditDist[count]*listSV[SV1].sampleSup[count] + listSV[SV2].sampleEditDist[count]*listSV[SV2].sampleSup[count])/(float)(listSV[SV1].sampleSup[count] + listSV[SV2].sampleSup[count]);
			listSV[SV1].sampleSup[count]=listSV[SV1].sampleSup[count]+listSV[SV2].sampleSup[count];	
		}
			listSV[SV2].yetGood=false;
			listSV[SV1].averageEditDist=(float)(listSV[SV1].averageEditDist*listSV[SV1].totalSup + listSV[SV2].averageEditDist*listSV[SV2].totalSup)/(listSV[SV1].totalSup + listSV[SV2].totalSup);
			listSV[SV1].totalSup=listSV[SV1].totalSup+listSV[SV2].totalSup; 		
	}else
	{	
		for (int count=0; count<totalNumSamples; count++)
		{	
		
			listSV[SV2].sampleEditDist[count]=(float)(listSV[SV1].sampleEditDist[count]*listSV[SV1].sampleSup[count] + listSV[SV2].sampleEditDist[count]*listSV[SV2].sampleSup[count])/(float)(listSV[SV1].sampleSup[count] + listSV[SV2].sampleSup[count]);
			listSV[SV2].sampleSup[count]=listSV[SV1].sampleSup[count]+listSV[SV2].sampleSup[count];	
			
		}
			
			listSV[SV1].yetGood=false;
			listSV[SV2].averageEditDist=(float)(listSV[SV1].averageEditDist*listSV[SV1].totalSup + listSV[SV2].averageEditDist*listSV[SV2].totalSup)/(listSV[SV1].totalSup + listSV[SV2].totalSup);
			listSV[SV2].totalSup=listSV[SV1].totalSup+listSV[SV2].totalSup; 		
	}

	for (int count=0; count<totalNumSamples; count++)
	{
		
			if (listSV[SV1].sampleSup[count]==0)
				listSV[SV1].sampleEditDist[count]=0;
			
			if (listSV[SV2].sampleSup[count]==0)
				listSV[SV2].sampleEditDist[count]=0;
	}

	

	
}

	
int findIdToStart(int count)
{
int id=count;
	while(id>0 && listSV[id].posInnerLeft>listSV[count].posInnerLeft-251)
	{
		id--;
	}
return id;
}

int findIdToEnd(int count)
{
int id=count;
	while(id<numSV && listSV[id].posInnerLeft<listSV[count].posInnerLeft+251)
	{
		id++;
	}
return id;
}


int main(int argv, char *argc[])
{
	FILE *fp=fopen(argc[1],"r"); // total SVs
	FILE *fp2=fopen(argc[2],"r"); // total samples file
	numSV=atoi(argc[3]);
	numberSampleInput=atoi(argc[4]); // number of samples in the input SV file
//	int totalSamples=20;
	float tempF;
	totalNumSamples=0;
	char sample[100];


	while(fscanf(fp2,"%s\n", samplesName[totalNumSamples])!=EOF)
	{
		totalNumSamples++;
	}



/*
	for (int count=0; count<160; count++)
	{
		fscanf(fp2,"%s\n", samplesName[count]);
		
	}
*/
//	printf("LLL\n");


//	chr1 Start_Outer:16510 Start_Inner:16795 End_Inner:16873 End_Outer:17132 SVtype:D sup:4 Sum_Weight:0 AvgEditDits:6.250000 Lib:SSC03070 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC03078 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC03092 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC03093 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05166 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05174 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05182 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05183 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05724 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05725 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC05726 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC12487 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC08058 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC08065 LibSup:2 LibHurScore:2 AvgEditDistInd:8 Lib:SSC08185 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC08186 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC10335 LibSup:2 LibHurScore:2 AvgEditDistInd:4.5 Lib:SSC10338 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC10358 LibSup:0 LibHurScore:0 AvgEditDistInd:0 Lib:SSC10379 LibSup:0 LibHurScore:0 AvgEditDistInd:0 minDelLen:6 maxDelLen:337

	
	for (int count=0; count<numSV; count++)
	{
		fscanf(fp, "Chr:%s Start_Outer:%i Start_Inner:%i End_Inner:%i End_Outer:%i SVtype:D sup:%i Sum_Weight:0 AvgEditDits:%f", listSV[count].chroName, &(listSV[count].posOuterLeft), &(listSV[count].posInnerLeft), &(listSV[count].posInnerRight), &(listSV[count].posOuterRight), &(listSV[count].totalSup), &(listSV[count].averageEditDist));
		//printf("%s %i\n", listSV[count].chroName, listSV[count].posOuterLeft);
		for (int count2=0; count2<numberSampleInput; count2++)
		{
			fscanf(fp," Lib:%s ", sample);
			//printf("%s\n", sample);
			for (int count2=0; count2<totalNumSamples; count2++)
			{
				if (strcmp(sample, samplesName[count2])==0)
				{
					fscanf(fp,"LibSup:%i LibHurScore:%f AvgEditDistInd:%f", &(listSV[count].sampleSup[count2]), &tempF, &(listSV[count].sampleEditDist[count2]));
				}
			}
		}
		listSV[count].yetGood=true;
			fscanf(fp," minDelLen:%i maxDelLen:%i\n", &(listSV[count].minLen), &(listSV[count].maxLen));
		
	}

	

//	printf("L113\n");


	qsort(listSV, numSV, sizeof(SV), compare);

	int maxSupSV=0;
	int maxSVId=0;
	int idStart, idEnd;
	bool additionalChange;
	do{
		additionalChange=false;
		for (int count=0; count<numSV; count++)
		{
			maxSupSV=0;
			maxSVId=-1;
	
			if (listSV[count].yetGood==true)
			{
				idStart=findIdToStart(count);
				idEnd=findIdToEnd(count);
			
				for (int id=idStart; id<idEnd; id++)
				{
				if (id!=count)
				{
					if (listSV[id].yetGood==true)
					{
					if (matchs(id, count)==true)
					{
						if (maxSupSV<listSV[id].totalSup)
						{
							maxSupSV=listSV[id].totalSup;
							maxSVId=id;
						}
					}
					}
				}
				}
			}
			if (maxSVId>-1)
			{
		//		printf("Merging %i %i\n", maxSVId, count); 
				mergeSV(maxSVId, count);
				additionalChange=true;
			}
		}
	//printf("LOOP\n");
	}while(additionalChange==true);

	int totalSV=0;
/*
	for (int count=0; count<numSV; count++)
	{
		if (listSV[count].yetGood==true)
		{
			printf("%s\t%i\t%i\t%i\n", listSV[count].chroName, listSV[count].posInnerLeft, listSV[count].posInnerRight, listSV[count].totalSup);
			totalSV++;
		}
	}
*/
	//printf("%i\n", totalSV);

	//printf("First run %i\n", listSV[binarySearch(0, numSV-1, 0)].posInnerLeft);
	//printf("Second run %i\n", listSV[binarySearch(0, numSV-1, 16056101)].posInnerLeft);
	//printf("thord run %i\n", listSV[binarySearch(0, numSV-1, 16340429)].posInnerLeft);
//	printf("L122\n");

/*	for (int count=0; count<numSV; count++)
	{
	printf("%i\n", count);
		for (int count2=count+1; count2<numSV; count2++)
		{
			if (matchs(count, count2)==true)
			{
				printf("%s %i %i %s %i %i\n", listSV[count].chroName, listSV[count].posInnerLeft, listSV[count].posInnerRight, listSV[count2].chroName, listSV[count2].posInnerLeft, listSV[count2].posInnerRight);
			}
		}
	}



//printf("%i\n", numSV);
*/

for (int count=0; count<numSV; count++)
{
	if ((listSV[count].yetGood==true))
	{
		printf("%s\t%i\t%i\t%s_%i\t%i\t%f", listSV[count].chroName, (listSV[count].posOuterLeft+listSV[count].posInnerLeft)/2, (listSV[count].posInnerRight + listSV[count].posOuterRight)/2,  listSV[count].chroName, count, listSV[count].totalSup, listSV[count].averageEditDist);
		for (int count2=0; count2<totalNumSamples; count2++)
		{
			printf(" %i %f", listSV[count].sampleSup[count2], listSV[count].sampleEditDist[count2]);
		}
		printf("\n");
	}	
}

/*
	for (int count=0; count<numSV; count++)
	{
		printf("%i\n", listSV[count].posInnerLeft);
	}
*/
}



