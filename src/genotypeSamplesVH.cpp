/*

This program was created at:  Wed Nov 25 10:16:59 2015
This program was created by:  Brad Nelson


Contact: bnelsj@uw.edu

Organization: Unviersity of Washington


The MIT License (MIT)

Copyright (c) <2015> <Brad Nelson>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*/

#include <string>
#include <iostream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include "split.h"
#include <fstream>
#include <map>
#include <algorithm>
#include "../tabixpp/tabix.hpp"

std::map<std::string, int> mei_names;

std::map<std::string, Tabix*> sampleRDPath;
std::vector<std::string> samples;

struct indivDat{
	double libSupport;
	double averageEdist;
	double readDepth;
	std::string genotype;
};

struct options{
	std::string file;
	std::string readDepthManifest;
	std::string MEI;
	int minSVLen;
	std::string probands;
}globalOpts;

static const char *optString = "hf:r:m:i:p:";

//-------------------------------   OPTIONS   --------------------------------
int parseOpts(int argc, char** argv)
    {
    int opt = 0;
    globalOpts.file = "NA";
    opt = getopt(argc, argv, optString);
    while(opt != -1){
	switch(opt){
		case 'h':
		{
		 std::cerr << "Useage" << std::endl;
		 break;
		}
		case 'f':
		{
		 globalOpts.file = optarg;
		 break;
		}
		case 'r':
		{
		 globalOpts.readDepthManifest = optarg;
		 break;
		}
		case 'm':
		{
		 globalOpts.MEI = optarg;
		 break;
		}
		case 'i':
		{
		globalOpts.minSVLen = atoi(optarg);
		break;
		}
		case 'p':
		{
		 globalOpts.probands = optarg;
		 break;
		}
		case '?':
		{
		 break;
		}
	}
 
  opt = getopt( argc, argv, optString ); 
   }
return 1;
}
//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
void sub()
{
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
std::string getInfoField(std::string &svName, std::vector<indivDat*> &dat)
{
	std::stringstream stream;
	stream << "\tSVTYPE=DEL;MEI=";
	if(mei_names.find(svName) != mei_names.end()) {
		stream << "1;";
	} else {
		stream << "0;";
	}

	int maxLibSup = 0;
	int nSampleSup = 0;
	int sumSupRD = 0;
	int sumNoSupRD = 0;

	int nSamples = dat.size();

	for(std::vector<indivDat*>::iterator j = dat.begin();
		j != dat.end(); j++) {
		if((*j)->libSupport > 0) {
			nSampleSup++;
			if((*j)->readDepth != -1) {
				sumSupRD += (*j)->readDepth;
			}

			if((*j)->libSupport > maxLibSup) {
				maxLibSup = (*j)->libSupport;
			}
		} else {
			if((*j)->readDepth != -1) {
				sumNoSupRD += (*j)->readDepth;
			}
		}
		// Genotype samples
		if((*j)->libSupport > 0 && (*j)->readDepth != -1 && (*j)->readDepth <= 1.6) {
			(*j)->genotype = "./1";
		} else if((*j)->readDepth != -1 && (*j)->readDepth < 1.4) {
			(*j)->genotype = "./1";
		} else {
			(*j)->genotype = "./0";
		}
	}

	float supRDMean = (float) sumSupRD/nSampleSup;
	float noSupRDMean = (float) sumNoSupRD/(nSamples - nSampleSup);

	int goodRDCall;
	if(supRDMean < 1.4 && noSupRDMean - supRDMean > 0.6)
		goodRDCall = 1;
	else
		goodRDCall = 0;

	stream << "GOODCNCALL=" << goodRDCall << ";SUPPORTCNMEAN=" << supRDMean << ";NOSUPPORTCNMEAN=" << noSupRDMean << ";MAXLIBSUPPORT=" << maxLibSup;
	stream << ";CIPOS=-10,0;CIEND=0,10";

	return stream.str();
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/
std::string indivDatGenotype(std::vector<indivDat*> &dat)
{
	std::stringstream stream;
	stream << "\tGT:LS:ED:CN";
	for(std::vector<indivDat*>::iterator j = dat.begin();
		j != dat.end(); j++) {
	  stream << "\t" << (*j)->genotype << ":" << (*j)->libSupport << ":" << (*j)->averageEdist << ":";
		if((*j)->readDepth == -1)
			stream << ".";
		else
			stream << (*j)->readDepth;
	}
	return stream.str();
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/

double getReadDepth(std::string &region, std::string &start, std::string &end, Tabix* tbx)
{
	std::string line;
	double rd = -1;
	tbx->setRegion(region);
	if(tbx->getNextLine(line)) {
		std::vector<std::string> linDat = split(line, "\t");
		if(start == linDat[1] && end == linDat[2]) {
			std::string rdStr = linDat[3].substr(4);
			rd = atof(rdStr.c_str());
		}
	}
	return rd;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  : vector of pointers to breakpoints and a bamtools RefVector

 Function does   : prints a vcf format

 Function returns: nada

*/


void printVCFHeader(vector<std::string> &samples){

  std::stringstream header;
  
  header << "##fileformat=VCFv4.2" << endl;
  header << "##source=VariationHunter" << endl;
  header << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
  header << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
  header << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
  header << "##INFO=<ID=POS,Number=2,Type=String,Description=\"POS and END\">" << endl;
  header << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">" << endl;
  header << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">" << endl;
  header << "##INFO=<ID=MEI,Number=1,Type=Integer,Description=\"SV intersects a known MEI\">" << endl;
  header << "##INFO=<ID=GOODCNCALL,Number=1,Type=Integer,Description=\"CN support for call based on allele balance\">" << endl;
  header << "##INFO=<ID=SUPPORTCNMEAN,Number=1,Type=Float,Description=\"Average CN of samples with read support for SV\">" << endl;
  header << "##INFO=<ID=NOSUPPORTCNMEAN,Number=1,Type=Float,Description=\"Average CN of samples without read support for SV\">" << endl;  
  header << "##INFO=<ID=MAXLIBSUPPORT,Number=1,Type=Integer,Description=\"Largest number of reads for a single sample that support the SV\">" << endl;
  header << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Hemizygous Genotype\">" << endl;
  header << "##FORMAT=<ID=LS,Number=1,Type=Integer,Description=\"Library reads supporting SV\">" << endl;
  header << "##FORMAT=<ID=ED,Number=1,Type=Float,Description=\"Average edit distance of supporting reads\">" << endl;
  header << "##FORMAT=<ID=CN,Number=1,Type=Float,Description=\"Copy number based on GC-corrected read depth\">" << endl;

  header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" ;

  for(vector<string>::iterator iz = samples.begin(); iz !=  samples.end(); iz++){
    header << "\t" << (*iz);
  }

  std::cout << header.str() << std::endl;
}

//------------------------------- SUBROUTINE --------------------------------
/*
 Function input  :

 Function does   :

 Function returns:

*/

void processLine(std::string &line)
{
	std::vector<std::string> lineDat;
	std::stringstream ss;

	lineDat = split(line, "\t");
	
	int svLen = atoi(lineDat[1].c_str()) - atoi(lineDat[2].c_str());

	ss << lineDat[0] << ":" << lineDat[1] << "-" << lineDat[2];
	std::string region = ss.str();

	std::vector<std::string> lineDat2 = split(lineDat[5], " ");
	std::vector<indivDat*> individuals;

	if((lineDat2.size() - 1) % 2 != 0) {
		std::cerr << "Error: Wrong number of columns in SV file" << std::endl;
		exit(1);
	}

	// Get read depth for calls > min size
	for(int i=1; i < lineDat2.size(); i+=2) {
		indivDat* indivSV = new indivDat;
		indivSV->libSupport = atof(lineDat2[i].c_str());
		indivSV->averageEdist = atof(lineDat2[i+1].c_str());

		if(sampleRDPath[samples[(i-1)/2]] != NULL && svLen > globalOpts.minSVLen) {
			indivSV->readDepth = getReadDepth(region, lineDat[1], lineDat[2], sampleRDPath[samples[(i-1)/2]]);
		} else {
			indivSV->readDepth = -1;
		}
		individuals.push_back(indivSV);
	}

	std::string svInfo = getInfoField(lineDat[3], individuals);

	std::cout << lineDat[0] << "\t" << lineDat[1] << "\t" << lineDat[3];
	std::cout << "\tN\t<DEL>\t.\t." << svInfo << ";SVLEN=" << svLen << indivDatGenotype(individuals) << std::endl;

	for(std::vector<indivDat*>::iterator j = individuals.begin();
		j != individuals.end(); j++) {
		delete (*j);
	}
}

//-------------------------------    MAIN     --------------------------------
/*
 Comments:
*/

int main( int argc, char** argv)
{
	globalOpts.minSVLen = -1;

	int parse = parseOpts(argc, argv);

	if (globalOpts.minSVLen == -1) {
		std::cerr << "Error: no min SV length specified for RD genotyping." << std::endl;
		exit(1);
	}

	// Read MEI file
	if (globalOpts.MEI.empty()) {
		std::cerr << "Error: no MEI file provided" << std::endl;
		exit(1);
	}

	std::ifstream mei;
	std::string line;
	mei.open(globalOpts.MEI.c_str()); 
	if(!mei.is_open()){
		std::cerr << "Error: Cannot open MEI file" << std::endl;
		exit(1);
	} else{
		while(mei.good()) {
			getline(mei, line);
			mei_names[line] = 1;
		}
	}

	mei.close();

	// Read depth manifest
	if (globalOpts.readDepthManifest.empty()) {
		std::cerr << "Error: no read depth manifest file provided" << std::endl;
		exit(1);
	}

	std::ifstream manifest;
	std::vector<std::string> sn_file_pair;
	manifest.open(globalOpts.readDepthManifest.c_str());
	if(!manifest.is_open()) {
		std::cerr << "Error: Cannot open read depth manifest file" << std::endl;
		exit(1);
	} else {
		while(manifest.good()) {
			getline(manifest, line);
			sn_file_pair = split(line, "\t");

			if(!sn_file_pair[0].empty()) {
				samples.push_back(sn_file_pair[0]);
			}

			// Get pointer map to tabix files
			if(!sn_file_pair[1].empty()) {
		   		Tabix* tmpTabix = new Tabix(sn_file_pair[1]);
		      	sampleRDPath[sn_file_pair[0]] = tmpTabix;
         } else {
				sampleRDPath[sn_file_pair[0]] = NULL;
			}
		}
	}
	manifest.close();
	

	// Read SV file

	if (globalOpts.file.empty()) {
		std::cerr << "Error: no file provided" << std::endl;
		exit(1);
	}
	std::ifstream infile;
	infile.open(globalOpts.file.c_str()); 
	if(!infile.is_open()){
		std::cerr << "Error: Cannot open file" << std::endl;
		exit(1);
	} else {
		printVCFHeader(samples);
		while(infile.good()) {
			getline(infile, line);
			if(line.empty()) {
				continue;
			}
			processLine(line);
		}
	}
	infile.close();

return 0;
}
