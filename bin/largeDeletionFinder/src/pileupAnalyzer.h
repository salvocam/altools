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


#include <string>
#include <dirent.h>

using namespace std;


#ifndef PILEUPANALYZER_H_
#define PILEUPANALYZER_H_

class pileupAnalyzer {
public:
	pileupAnalyzer();
	void setPileupFolderName(string);
	double getAverageCoverageInRange(string, int, int);
	virtual ~pileupAnalyzer();
	int loadFiles(void);
	void performSlidingAnalysis(string,int, int, int, int);

	string pileupFolderPath;

	int numberOfFiles;
	string fileInInputFolder[200];

	DIR *inputFolder;
	struct dirent *ent;
};

#endif /* PILEUPANALYZER_H_ */
