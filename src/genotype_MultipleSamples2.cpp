#include <stdio.h>
#include <stdlib.h>
#include <string.h>


const int maxNumSamples=1000;

typedef struct SV{
	char chroName[10];
	int posOuterLeft;
	int posInnerLeft;
	int posInnerRight;
	int posOuterRight;

	int sampleSup[maxNumSamples];
	float sampleEditDist[maxNumSamples];
	float wssdCNV_BWA[maxNumSamples]; // BWA based wssd 
	int genotypeArray[maxNumSamples];
	char SVName[20];
	float averageEditDist;
	int totalSup;
	int minLen;
	int maxLen;
	bool yetGood;
	//int windowsCovered;
	bool goodCall; 
	bool mrsFastReadDepthGoodCall;
	bool bwaReadDepthGoodCall;
	bool MEIgoodCall;
	int maxSup;
	
	float CNV1_BWA, CNV2_BWA;

}SV;


SV listSV[5000000];
int countSV;

char samplesName[maxNumSamples][100];
char samplesOtherName[maxNumSamples][100];
int isProbandOrSibling[maxNumSamples];// indicats if a smaple is probands ot sibling 
char wssdSamplesName[maxNumSamples][100];
int convertWSSDSampleNameToSV[maxNumSamples];// Array[WSSD]-> SV sample
int totalNumSamples;

int compare(const void *a, const void *b)
{
        return (strcmp(((SV*)a)->SVName, ((SV*)b)->SVName));
}










/*int findChroArray(char *chroName)
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
*/


















int binarySearch(int idStart, int idEnd, char *SVName)
{
int idMid=0;
	if (strcmp(listSV[idStart].SVName, SVName)==0)
		return idStart;
	if (strcmp(listSV[idEnd].SVName, SVName)==0)
		return idEnd;
	idMid=(idStart+idEnd)/2;
	if (strcmp(listSV[idMid].SVName, SVName)==0)
		return idMid;
	
	if (strcmp(listSV[idMid].SVName, SVName)<0)
	{
		return (binarySearch(idMid, idEnd, SVName));
	}else if (strcmp(listSV[idMid].SVName, SVName)>0)
	{
		return (binarySearch(idStart, idMid, SVName));
	}

}


int  assignBWA_WSSD_ReadDepthFrom(FILE *fp7, FILE *fp8)
{
	int SVIdArray[500000];
	char SVNameTemp[100];
	int counSVNameList=0;
	char fileName[1000];
	char sample[100];
	int sampleId;
	FILE *fp9;
	int pos1, pos2;
	char chroName[10];
	float tempF, cnvCount;
	while(fscanf(fp8, "%s\n", SVNameTemp)!=EOF)
	{
		SVIdArray[counSVNameList]=binarySearch(0, countSV, SVNameTemp);
		counSVNameList++;
	}


int countLoci2=0;
	while(fscanf(fp7,"%s\t%s\n", fileName, sample)!=EOF)
	{
countLoci2=0;
		for (int count=0; count<totalNumSamples; count++)
		{
			if (strcmp(sample, samplesName[count])==0)
				sampleId=count;
		}

		fp9=fopen(fileName, "r");
		
		//chr10	74306	74825	CNV:1.938747	0.410072	0.018000	
		
		for (int count=0; count<counSVNameList; count++)
		{
			fscanf(fp9,"%s\t%i\t%i\tCNV:%f\t%f\t%f\n", chroName, &pos1, &pos2, &(listSV[SVIdArray[count]].wssdCNV_BWA[sampleId]), &tempF, &tempF);
		}
		

//		for (int count=0; count<countSV; count++)
//		{
		/*while(fscanf(fp9,"%s\t%i\t%i\tCNV:%f\t%f\t%f\n", chroName, &pos1, &pos2, &cnvCount, &tempF, &tempF)!=EOF)
		{
		countLoci2++;
		printf("%i %i\n", pos1, countLoci2);
			for (int count2=0; count2<countSV; count2++)
			{
				if (strcmp(chroName, listSV[count2].chroName)==0 && pos1==listSV[count2].posInnerLeft && listSV[count2].posInnerRight==pos2)
				{
					listSV[count2].wssdCNV_BWA[sampleId]=cnvCount;
					continue;
				}
			}


		}

		
		fclose(fp9);
	*/
	}
	

}



int main(int argv, char *argc[])
{


	FILE *fp1=fopen(argc[1],"r");// Samples Name
	FILE *fp2=fopen(argc[2],"r");// the SV paired-end info 
//	FILE *fp3=fopen(argc[3],"r");// the number of windows
//	FILE *fp4=fopen(argc[4],"r");// wssd CNV (mrsFast/Peter version)	
	FILE *fp5=fopen(argc[3],"r");// Known ALu, L1 calls
	FILE *fp6=fopen(argc[4],"r");// probands name

	FILE *fp7=fopen(argc[5],"r");// read-depth BAM file using BWA mapping
	FILE *fp8=fopen(argc[6],"r");// the name of the SVs in the BWA readdepth7
	char tempStr[100];


	while(fscanf(fp1,"%s\n", samplesName[totalNumSamples])!=EOF)
	{
		totalNumSamples++;
	}

/*	for (int count=0; count<160; count++)
	{
	//	fscanf(fp1,"%s\t%s\t%s\t%s\n", samplesName[count], samplesOtherName[count], tempStr, tempStr);
		fscanf(fp1,"%s\n", samplesName[count]);
	}
*/



	while(fscanf(fp2,"%s\t%i\t%i\t%s\t%i\t%f", listSV[countSV].chroName, &listSV[countSV].posInnerLeft, &listSV[countSV].posInnerRight, listSV[countSV].SVName, &listSV[countSV].totalSup, &listSV[countSV].averageEditDist)!=EOF)
	{
		for (int count=0; count<totalNumSamples; count++)
		{
			fscanf(fp2, "\t%i\t%f", &listSV[countSV].sampleSup[count], &listSV[countSV].sampleEditDist[count]);
			if (listSV[countSV].sampleSup[count]>0)
				listSV[countSV].genotypeArray[count]=1;
		}
		fscanf(fp2,"\n");
		countSV++;
	}//	printf("L91\n");

	printf("L253\n");
	
	qsort(listSV, countSV, sizeof(SV), compare);
	
	printf("L257\n");

//chr10   74306   74825   chr10_5 0

	char chroName[10];
	int pos1, pos2, winCount;
	char SVName[20];
	char tempSampleName[100];

	while(fscanf(fp6,"%s\n", tempSampleName)!=EOF)
	{
		for (int count=0; count<totalNumSamples; count++)
		{
			if (strcmp(tempSampleName,samplesName[count])==0)
			{
				isProbandOrSibling[count]=1;
			}
		} 
	}

//contig  start   end     name    SSC05075        SSC06572        SSC06563        SSC07807        SSC07806        SSC05066        SSC07802        SSC05070        SSC09060     SSC10379        SSC09067        SSC07131        SSC04980        SSC07130        SSC04987        SSC04986        SSC08188        SSC08185        SSC09514    SSC08187 SSC08518


	int SVIdFound=0;
	/*
	while(fscanf(fp3,"%s\t%i\t%i\t%s\t%i\n", chroName, &pos1, &pos2, SVName, &winCount)!=EOF)
	{
		SVIdFound=binarySearch(0, countSV-1, SVName);
		//if (strcmp(SVName, listSV[count].SVName)==0)
		listSV[SVIdFound].windowsCovered=winCount;
	}*/

	//printf("L115\n");

/*
	fscanf(fp5,"%s\t%i\t%i\n", chroName, pos1, pos2)
	{
		for (int count=0; count<countSV; count++)
		{
			if ((listSV[countSV].posInnerRight-listSV[countSV].posInnerLeft>300 && listSV[countSV].posInnerRight-listSV[countSV].posInnerLeft<500)||(listSV[countSV].posInnerRight-listSV[countSV].posInnerLeft>5000 && listSV[countSV].posInnerRight-listSV[countSV].posInnerLeft<7000))
			{
				if (strcmp(chroName, listSV[countSV].chroName)==0)
				{
					if (listSV[countSV].posInnerRight>pos1 && listSV[countSV].posInnerLeft<pos1 && listSV[countSV].posInnerRight-100<pos1 && listSV[countSV].posInnerLeft+100>pos1)
					{
						listSV[SVIdFound].MEIgoodCall=true;						
					}
				}
			}
		}	
	}

*/

/*	fscanf(fp4,"contig\tstart\tend");
	for (int count=0; count<160; count++)
	{
		fscanf(fp4,"\t%s", wssdSamplesName[count]);
		for (int count2=0; count2<160; count2++)
		{
			if (strcmp(wssdSamplesName[count], samplesName[count2])==0)
				convertWSSDSampleNameToSV[count]=count2;
		}
	}    
	fscanf(fp4,"\n");
*/

	//printf("L131\n");	

/*	for (int count=0; count<countSV; count++)
	{
		fscanf(fp4,"%s\t%i\t%i\t%s", chroName, &pos1, &pos2, SVName);
	//	if (strcmp(SVName, listSV[count].SVName)!=0)
	//		printf("EROROROROR\n");
		SVIdFound=binarySearch(0, countSV-1, SVName);
		//printf("%i\t%s\t%s\n", SVIdFound, listSV[SVIdFound].SVName, SVName); 
		for (int count2=0; count2<160; count2++)
		{
			fscanf(fp4, "\t%f", &listSV[SVIdFound].wssdCNV[convertWSSDSampleNameToSV[count2]]);
		}	
		fscanf(fp4,"\n");
	

	
	}
*/

	while(fscanf(fp5, "%s\n", SVName)!=EOF)
	{
		SVIdFound=binarySearch(0, countSV-1, SVName);
	
		listSV[SVIdFound].MEIgoodCall=true;
		listSV[SVIdFound].goodCall=true;
	}

	printf("L150\n");

	assignBWA_WSSD_ReadDepthFrom(fp7,fp8);

	printf("L232\n");	
	int countTemp=0;
	//float sumCNV1=0;
	float sumCNV1_BWA=0;
	
	//float sumCNV2=0;
	float sumCNV2_BWA=0;
	
	int maxSup;

	for (int count=0; count<countSV; count++)
	{
	countTemp=0;
	//sumCNV1=0;
	sumCNV1_BWA=0;
	//sumCNV2=0;
	sumCNV2_BWA=0;
	maxSup=0;
		//if (listSV[count].windowsCovered>0)
		{
			for (int count2=0; count2<totalNumSamples; count2++)
			{
				if (listSV[count].sampleSup[count2]>0)
				{
					countTemp++;
					//sumCNV1=sumCNV1+listSV[count].wssdCNV[count2];
					sumCNV1_BWA=sumCNV1_BWA+listSV[count].wssdCNV_BWA[count2];
					if (listSV[count].sampleSup[count2]>maxSup)
						maxSup=listSV[count].sampleSup[count2];
				}
				else{
					//sumCNV2=sumCNV2+listSV[count].wssdCNV[count2];
					sumCNV2_BWA=sumCNV2_BWA+listSV[count].wssdCNV_BWA[count2];
				}
			}
			listSV[count].maxSup=maxSup;
		//printf("%s\t%i\t%i\t%i\t%f\t%f\t%i\t%i\t%s\n", listSV[count].chroName, listSV[count].posInnerLeft, listSV[count].posInnerRight, listSV[count].windowsCovered, (float)sumCNV1/(float)countTemp, (float)sumCNV2/(float)(160-countTemp), listSV[count].totalSup, maxSup, listSV[count].SVName);
		/*	if (listSV[count].windowsCovered>0 && (float)sumCNV1/(float)countTemp<1.4 && (float)sumCNV2/(float)(160-countTemp) - (float)sumCNV1/(float)countTemp > 0.6)
			{
				listSV[count].goodCall=true;
				listSV[count].mrsFastReadDepthGoodCall=true;
				for (int count2=0; count2<maxNumSamples; count2++)
				{
					if (listSV[count].wssdCNV[count2]<1.4 && listSV[count].sampleSup[count2]==0)
					{
						listSV[count].genotypeArray[count2]=1;
					}
					if (listSV[count].wssdCNV[count2]>1.6 && listSV[count].sampleSup[count2]>0)
					{
						listSV[count].genotypeArray[count2]=0;
					}
				}
				listSV[count].CNV1_mrsFast=(float)sumCNV1/(float)countTemp;
				listSV[count].CNV2_mrsFast=(float)sumCNV2/(float)(160-countTemp);
				listSV[count].CNV1_BWA=(float)sumCNV1_BWA/(float)countTemp;
				listSV[count].CNV2_BWA=(float)sumCNV2_BWA/(float)(160-countTemp);
				
			}else //if (listSV[count].windowsCovered==0 && (float)sumCNV1_BWA/(float)countTemp<1.4 && (float)sumCNV2_BWA/(float)(160-countTemp) - (float)sumCNV1_BWA/(float)countTemp > 0.6)
		*/	
		if ((float)sumCNV1_BWA/(float)countTemp<1.4 && (float)sumCNV2_BWA/(float)(totalNumSamples-countTemp) - (float)sumCNV1_BWA/(float)countTemp > 0.6)
			{
				
				listSV[count].bwaReadDepthGoodCall=true;
				listSV[count].goodCall=true;
				//listSV[count].realGoodCall=true;
				for (int count2=0; count2<maxNumSamples; count2++)
				{
					if (listSV[count].wssdCNV_BWA[count2]<1.4 && listSV[count].sampleSup[count2]==0)
					{
						listSV[count].genotypeArray[count2]=1;
					}
					if (listSV[count].wssdCNV_BWA[count2]>1.6 && listSV[count].sampleSup[count2]>0)
					{
						listSV[count].genotypeArray[count2]=0;
					}
				}
				listSV[count].CNV1_BWA=(float)sumCNV1_BWA/(float)countTemp;
				listSV[count].CNV2_BWA=(float)sumCNV2_BWA/(float)(totalNumSamples-countTemp);
			} /*else if ((float)sumCNV1_BWA/(float)countTemp<1)
			{
				listSV[count].bwaReadDepthGoodCall=true;
				listSV[count].goodCall=true;
			}*/
		}
	}



	for (int count=0; count<countSV; count++)
	{
		if (listSV[count].posInnerRight-listSV[count].posInnerLeft>750 && listSV[count].maxSup>5 && listSV[count].totalSup>10 && listSV[count].goodCall==false)
		//if (listSV[count].posInnerRight-listSV[count].posInnerLeft>750 && listSV[count].maxSup>5 && listSV[count].totalSup>10 && listSV[count].goodCall==false)
		{
		countTemp=0;
		sumCNV1_BWA=0;
		sumCNV2_BWA=0;
			listSV[count].goodCall=true;
			
			for (int count2=0; count2<totalNumSamples; count2++)
			{
				if (listSV[count].sampleSup[count2]>0)
				{
					countTemp++;
					//sumCNV1=sumCNV1+listSV[count].wssdCNV[count2];
					sumCNV1_BWA=sumCNV1_BWA+listSV[count].wssdCNV_BWA[count2];
					//if (listSV[count].sampleSup[count2]>maxSup)
					//	maxSup=listSV[count].sampleSup[count2];
				}
				else{
					//sumCNV2=sumCNV2+listSV[count].wssdCNV[count2];
					sumCNV2_BWA=sumCNV2_BWA+listSV[count].wssdCNV_BWA[count2];
				}
			}
		
			listSV[count].CNV1_BWA=(float)sumCNV1_BWA/(float)countTemp;
			listSV[count].CNV2_BWA=(float)sumCNV2_BWA/(float)(totalNumSamples-countTemp);

		}
			
	}


	int genotypeCount=0;
	int numProbandSib=0;
	int denovoSampleId=0;

	printf("Chro\tPos1\tPos2\tSVName\tTotalSup\tMaxSup");
	for (int count=0; count<totalNumSamples;count++)
	{
		printf("\t%s", samplesName[count]);		
	}
	printf("\n");
		

	for (int count=0; count<countSV; count++)
	{
	genotypeCount=0;
	numProbandSib=0;
	denovoSampleId=0;
		if (listSV[count].goodCall==true)
		{

			for (int count2=0; count2<totalNumSamples; count2++)
			{
				if (listSV[count].genotypeArray[count2]==1)
				{
					genotypeCount++;
				}
				if (listSV[count].genotypeArray[count2]==1 &&  isProbandOrSibling[count2]==1)
				{
					numProbandSib++;
					denovoSampleId=count2;
				}
			}
			for (int count2=0; count2<totalNumSamples; count2++)
			{
			//chr10	131095256	131096444	1	SSC10104	chr10_100117	
				if (listSV[count].genotypeArray[count2]==1)
				{
					if (!(listSV[count].bwaReadDepthGoodCall==false && listSV[count].MEIgoodCall==false && genotypeCount>totalNumSamples*0.9))
					{
						printf("%s\t%i\t%i\t1\t%s\t%s", listSV[count].chroName, listSV[count].posInnerLeft, listSV[count].posInnerRight, samplesName[count2], listSV[count].SVName);
						if (listSV[count].MEIgoodCall==true)
							printf("\tMEI");
						else if (listSV[count].bwaReadDepthGoodCall==true)
							printf("\tReadDepth");
						if (numProbandSib==1 && genotypeCount==1)
							printf("\tDeNovo\n");
						else printf("\n");
					}
				}
			}
		}
	}
	


}

