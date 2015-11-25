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

std::map<std::string, int> mei_names;

struct indivDat{
	double libSupport;
	double averageEdist;
};

struct options{
   std::string file;
	std::string readDepthManifest;
	std::string MEI;
	std::string probands;
}globalOpts;

static const char *optString = "hf:r:m:p:";

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
std::string indivDatGenotype(std::vector<indivDat*> &dat)
{
	std::stringstream stream;
	stream << "\tLS:ED";
	for(std::vector<indivDat*>::iterator j = dat.begin();
		j != dat.end(); j++) {
		stream << "\t" << (*j)->libSupport << ":" << (*j)->averageEdist;
	}
	return stream.str();
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
	lineDat = split(line, "\t");
	std::vector<std::string> lineDat2 = split(lineDat[5], " ");
	std::vector<indivDat*> individuals;

	if((lineDat2.size() - 1) % 2 != 0) {
		std::cerr << "Error: Wrong number of columns in SV file" << std::endl;
		exit(1);
	}
	std::cout << lineDat2.size() << std::endl;
	for(int i=1; i < lineDat2.size(); i+=2) {
		indivDat* tmp = new indivDat;
		tmp->libSupport = atof(lineDat2[i].c_str());
		tmp->averageEdist = atof(lineDat2[i+1].c_str());
		individuals.push_back(tmp);
	}

	if(mei_names.find(lineDat[3]) != mei_names.end()) {
		std::cout << "Found mei call " << line << std::endl;
		exit(0);
	}

	std::cout << indivDatGenotype(individuals) << std::endl;
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
int parse = parseOpts(argc, argv);

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
		while(infile.good()) {
			getline(infile, line);
			processLine(line);
		}
	}
	infile.close();

return 0;
}
