/* Copyright 2015 Salvatore Camiolo
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "completePileupAnalyzer.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <dirent.h>

DIR *inputFolder;
struct dirent *ent;

string pileupFolder;










completePileupAnalyzer::completePileupAnalyzer() {
	// TODO Auto-generated constructor stub

}




void completePileupAnalyzer::collectPolymorphysms(void)
{
	int a;
	ifstream pileupFile;
	int position, coverage, snp, indel;

	loadPileupFiles();
	//loading polymorphysms positions
	for(a=1;a<=numberOfFiles;a++)
	{
		cout << "Loading polymorphysms for chromosome " << fileInPileupFolder[a]<< endl;
        
		pileupFile.open((pileupFolder+fileInPileupFolder[a]).c_str());
		numberOfSnp[a]=0;
		numberOfIndel[a]=0;

		while(!pileupFile.eof())
		{
			pileupFile >> position >> coverage >> snp >> indel;
            
			//cout << position << "  " << coverage << "  " << snp << endl;
            
			if (snp==1)
			{
				numberOfSnp[a]++;
				snpPosition[a][numberOfSnp[a]]=position;
			}
			if (indel==1)
			{
				numberOfIndel[a]++;
				indelPosition[a][numberOfIndel[a]]=position;
			}
		}
		pileupFile.close();

	}
    

}

int completePileupAnalyzer::searchSnpInPileup(string chromosome,int position)
{
	int a,i;
	for(a=1;a<numberOfFiles;a++)
	{
		if (fileInPileupFolder[a]==chromosome) break;
	}


	for(i=0;i<numberOfSnp[a];i++)
	{
		if(snpPosition[a][i]==position) return 1;//{cout << "Troato!" << endl; getchar(); return 1;}
		if(snpPosition[a][i]>position) return 0;//{getchar(); return 0;}
	}
}


int completePileupAnalyzer::searchIndelInPileup(string chromosome,int position)
{
	int a,i;
	for(a=1;a<numberOfFiles;a++)
	{
		if (fileInPileupFolder[a]==chromosome) break;
	}

	for(i=0;i<numberOfIndel[a];i++)
	{
		//cout << indelPosition[a][i] << "  " << position << endl;
		if(indelPosition[a][i]==position) return 1;
		if(indelPosition[a][i]>position) return 0;
	}
}

int completePileupAnalyzer::loadPileupFiles(void)
{

			inputFolder = opendir (pileupFolder.c_str());
			if (inputFolder != NULL) {

			/* load all the files and directories within directory */

				numberOfFiles=0;
				while ((ent = readdir (inputFolder)) != NULL)
				{
					numberOfFiles++;
					fileInPileupFolder[numberOfFiles] = ent->d_name;
					if (fileInPileupFolder[numberOfFiles]=="." || fileInPileupFolder[numberOfFiles] == "..") numberOfFiles--;
				}
			closedir (inputFolder);
			} else {
			/* could not open directory */
			perror ("");
			return EXIT_FAILURE;
			}
			return 1;
}

void completePileupAnalyzer::setFolder(string folderName)
{
	pileupFolder = folderName;
	return;
}
completePileupAnalyzer::~completePileupAnalyzer() {
	// TODO Auto-generated destructor stub
}

