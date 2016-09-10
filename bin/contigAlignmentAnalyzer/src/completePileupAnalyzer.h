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

#ifndef COMPLETEPILEUPANALYZER_H_
#define COMPLETEPILEUPANALYZER_H_
#include <string>

using namespace std;




class completePileupAnalyzer {
public:
	completePileupAnalyzer();
	virtual ~completePileupAnalyzer();
	void setFolder(string);
	void collectPolymorphysms(void);
	int loadPileupFiles(void);
	int searchSnpInPileup(string,int);
	int searchIndelInPileup(string,int);
    
    int (*snpPosition)[10000000] = new int[100][10000000];
    
    int (*indelPosition)[10000000] = new int[100][10000000];
    
    int numberOfSnp[1000];
    int numberOfIndel[1000];
    
    int numberOfChromosomes;
    int numberOfFiles;
    string fileInPileupFolder[1000];
    
};

#endif /* COMPLETEPILEUPANALYZER_H_ */
