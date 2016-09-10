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

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>


using namespace std;

ifstream inputFile;
ifstream referenceFile;
ofstream outpuFile;

ofstream highCoverageFile;
ofstream highCoverageGEFile;
ofstream zeroCoverageFile;
ofstream zeroCoverageGEFile;
ofstream logRatioFile;
ofstream referenceCoverageFile;
ifstream cnvFile;
ifstream cnvSmoothedFile;
ifstream referenceCoverageFileRead;

ofstream losses;
ofstream gains;
ofstream losses_GE;
ofstream gains_GE;
ofstream losses_Smoothed;
ofstream gains_Smoothed;
ofstream losses_Smoothed_GE;
ofstream gains_Smoothed_GE;

string fileInInputFolder[100000];
string fileInReferenceFolder[10000];
string inputDirectory;
string referencePileupFolder;
string outputFolder;


char command[1000];
char filename[1000];
string linuxCommand;

int numberOfFiles,a;
int chromosomeToAnalyze;
int minimumReadsNumberLowCoverage;
int minimumReadsNumberHighCoverage;
int minimumAreaToOutput;
int verboseOutput;
int generateGEFiles;
int analyzedLines = 0;
int analyzedWindows = 0;
int analyzedReferenceLines = 0;
int deletionStart;
int deletionArea = 0;
int zeroCoverageInWindow = 0;
int zeroCoverageInReferenceWindow = 0;
int zeroCoverageWindows=0;
int	zeroCoverageWindowsReference =0;
int posInFile;
int distanteToCollapse;
int winSize;

float percentageCutoffLossAndGain;
float windowSize = 1;
float windowCoverage = 0;
float referenceWindowCoverage = 0;
float lossCutoff = 0;
float gainCutoff = 0;
float referenceAverageCoverage = 0;
float subjectAverageCoverage = 0;
float segmentCoverage[1000000];



DIR *inputFolder;
DIR *referencePileupDirectory;
struct dirent *ent;

int loadInputFiles(void);
float *coverage = new float[500000000];
float *referenceCoverage = new float[500000000];
int *sortedCoverage = new int[500000000];
int *sortedReferenceCoverage = new int[500000000];
void reportZeroCoverageAreas(void);
void reportHighDepthAreas(int coverageCutoff);
void sortCoverageValues(int left, int right);
void extractLossesAndGainsCutoff(void);
void generateLossesAndGainsOutput(void);
void refineTables(void);
int extractHighCoverageCutoff(int *coverageArray, int percentage, int linesAnalyzed);


int main(int argc, char** argv) {
    
    
    //get parameters from gui;
    inputDirectory = argv[1];
    minimumReadsNumberLowCoverage = atoi(argv[2]);
    minimumAreaToOutput = atoi(argv[3]);
    verboseOutput = atoi(argv[4]);
    generateGEFiles = atoi(argv[5]);
    outputFolder = argv[6];
    referencePileupFolder = argv[7];
    windowSize = atof(argv[8]);
    distanteToCollapse = atoi(argv[9]);
    
    
    winSize = atoi(argv[8]);
    //Collect file names in the pileup folder and put them in the fileInInputFolder[] string array
    loadInputFiles();
    
    //create the folder
    linuxCommand = "mkdir -p ";
    linuxCommand.append(outputFolder);
    system(linuxCommand.c_str());
    
    //create in the output folder the files logRatioFile and referenceCoveragefile
    // and write field names
    strcpy(filename,outputFolder.c_str());
    strcat(filename,"/logRatioFile");
    logRatioFile.open(filename);
    
    logRatioFile << "Position	Chromosome	logRatio" << endl;
    
    strcpy(filename,outputFolder.c_str());
    strcat(filename,"/referenceCoverageFile");
    referenceCoverageFile.open(filename);
    
    referenceCoverageFile << "Position	Chromosome	Coverage" << endl;
    
    //Create needed output files.
    //GE files if required. With fields names
    if(generateGEFiles==1)
    {
        zeroCoverageGEFile.open( (outputFolder+"ZeroCoverage_GE.txt").c_str() );
        losses_GE.open( (outputFolder+"Losses_GE.txt").c_str() );
        gains_GE.open( (outputFolder+"Gains_GE.txt").c_str() );
        losses_Smoothed_GE.open( (outputFolder+"Losses_smoothed_GE.txt").c_str() );
        gains_Smoothed_GE.open( (outputFolder+"Gains_smoothed_GE.txt").c_str() );
        
        zeroCoverageGEFile << "Scaffold\tChromosome Start   End Length" << endl;
        /*
        losses_GE << "Scaffold\tChromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
        gains_GE << "Scaffold\tChromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
        losses_Smoothed_GE << "Scaffold\tChromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
        gains_Smoothed_GE << "Scaffold\tChromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;*/
    }
    //All other files
    gains.open((outputFolder+"Gains.txt").c_str() );
    losses.open((outputFolder+"Losses.txt").c_str() );
    gains_Smoothed.open((outputFolder+"Gains_smoothed.txt").c_str() );
    losses_Smoothed.open((outputFolder+"Losses_smoothed.txt").c_str() );
    zeroCoverageFile.open((outputFolder+"ZeroCoverage.txt").c_str() );
    
    //Write field names on output files
    zeroCoverageFile << "Chromosome Start   End Length" << endl;
    losses << "Chromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
    gains << "Chromsosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
    losses_Smoothed << "Chromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
    gains_Smoothed << "Chromosome\tStart\tEnd\tLength\tCoverageRatio\tCoverageInReference\tCoverageInSubject" << endl;
    
    
    /* Analyze chromosomes reported in the pileup folder */
    chromosomeToAnalyze = 0;
    while(chromosomeToAnalyze<numberOfFiles)
    {
        chromosomeToAnalyze++;
        if(fileInInputFolder[chromosomeToAnalyze]!="." && fileInInputFolder[chromosomeToAnalyze]!="..")
        {
            inputFile.open((inputDirectory+fileInInputFolder[chromosomeToAnalyze]).c_str());
            
            referenceFile.open((referencePileupFolder+"/"+fileInInputFolder[chromosomeToAnalyze]).c_str());
            
            //analyze zero coverage areas and create coverage ratio informations
            reportZeroCoverageAreas();
            
            inputFile.close();
        }
    }
    
    
    
    logRatioFile.close();
    cout << "Extracting copy number variations......" << endl;
    cout.flush();
    
    strcpy(command,"Rscript DNAcopy.R ");
    
    strcat(command,outputFolder.c_str());
    cout << "Executing " << command << endl;
    cout.flush();
    system(command);
    
    generateLossesAndGainsOutput();
    
    cout << "The reference average Coverage is " << referenceAverageCoverage << endl;
    cout.flush();
    cout << "The subject coverage average is " << subjectAverageCoverage << endl;
    cout.flush();
    
    
    
    cout << "The coverage analysis was successfully performed. " << endl;
    cout << "Thanks for using Altools" << endl;
    cout << "Press any key to close this window" << endl;
    cout.flush();
    getchar();
    
    zeroCoverageFile.close();
    highCoverageFile.close();
    cnvFile.close();
    cnvSmoothedFile.close();
    losses.close();
    gains.close();
    losses_Smoothed.close();
    gains_Smoothed.close();
    
    if( generateGEFiles==1 )
    {
        zeroCoverageGEFile.close();
        highCoverageGEFile.close();
        losses_GE.close();
        gains_GE.close();
        losses_Smoothed_GE.close();
        gains_Smoothed_GE.close();
    }
    
    return 0;
}



void reportZeroCoverageAreas(void)
{
    
    
    int a;
    int position1, coverage1, snp1, indel1;
    int position2, coverage2, snp2, indel2;
    int zeroCoverageDetected;
    
    
    
    //Load subject and reference genome coverage data in memory
    cout << "Creating windows for subject and reference genomes for chromosome " <<fileInInputFolder[chromosomeToAnalyze] << endl;
    cout.flush();
    
    analyzedLines = 0;
    analyzedWindows = 0;
    
    while(!inputFile.eof())
    {
        windowCoverage = 0;
        referenceWindowCoverage = 0;
        
        if(!inputFile.eof() && !referenceFile.eof())
        {
            
            zeroCoverageInWindow = 0;
            zeroCoverageInReferenceWindow = 0;
            
            for(a=0;a<windowSize;a++)
            {
                analyzedLines++;
                
                if (fmod(analyzedLines ,3000000) == 0)
                {
                    cout << analyzedLines << " analyzed in chromosome "
                    << fileInInputFolder[chromosomeToAnalyze] << endl;
                    cout.flush();
                }
                
                
                if(!inputFile.eof()) inputFile >> position1 >> coverage1 >> snp1 >> indel1;
                else break;
                
                
                
                if(!referenceFile.eof()) referenceFile >> position2 >> coverage2 >> snp2 >> indel2;
                
                
                else break;
                
                
                
                windowCoverage = windowCoverage + coverage1;
                if (coverage1 == 0) zeroCoverageInWindow++;
                referenceWindowCoverage = referenceWindowCoverage + coverage2;
                if (coverage2 == 0) zeroCoverageInReferenceWindow++;
            }
            if(windowSize!=zeroCoverageInWindow) coverage[analyzedWindows] = windowCoverage / (windowSize - zeroCoverageInWindow);
            else coverage[analyzedWindows] = 0;
            
            if (windowSize!=zeroCoverageInReferenceWindow) referenceCoverage[analyzedWindows] = referenceWindowCoverage / (windowSize - zeroCoverageInReferenceWindow);
            else referenceCoverage[analyzedWindows] = 0;
        }
        else break;
        
        analyzedWindows++;
        
        
    }
    inputFile.close();
    referenceFile.close();
    
    
    //fillup zeroCoverage files
    cout << "Extracring zero coverage data......" << endl;
    cout.flush();
    
    
    zeroCoverageWindows = 0;
    zeroCoverageWindowsReference = 0;
    
    
    for (a=0;a<analyzedWindows;a++)
    {
        
        referenceCoverageFile << a*windowSize << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
        << referenceCoverage[a] << endl;
        
        if(coverage[a]==0) zeroCoverageWindows++;
        else subjectAverageCoverage = subjectAverageCoverage + coverage[a];
        
        if(referenceCoverage[a]==0) zeroCoverageWindowsReference++;
        else referenceAverageCoverage = referenceAverageCoverage + referenceCoverage[a];
        
        
        if(coverage[a]!=0 && referenceCoverage[a]!=0)
            logRatioFile << a*windowSize << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
            << (coverage[a]/referenceCoverage[a]) << endl;
        
        
        deletionArea = 0;
        deletionStart = a*windowSize;
        if(coverage[a]==0 && referenceCoverage[a]!=0)
        {
            deletionArea = windowSize;
            while(coverage[a+1]==0 && referenceCoverage[a]!=0)
            {
                deletionArea = deletionArea + windowSize;
                a=a+1;
                if(a==analyzedWindows) break;
            }
            
            if(deletionArea >= minimumAreaToOutput)
            {
                zeroCoverageFile << fileInInputFolder[chromosomeToAnalyze] << "\t"
                << deletionStart << "\t" << (int)(deletionStart+deletionArea) << "\t" << (int)deletionArea << endl;
                if (generateGEFiles == 1) zeroCoverageGEFile << "undefinedScaffold" << "\t" << fileInInputFolder[chromosomeToAnalyze] << "\t"
                    << deletionStart << "\t" << (int)(deletionStart+deletionArea+windowSize) << "\t" << (int)(deletionArea+windowSize) << endl;
            }
            
        }
    }
    
    referenceAverageCoverage = referenceAverageCoverage / (analyzedWindows-1 - zeroCoverageWindowsReference);
    subjectAverageCoverage = subjectAverageCoverage / (analyzedWindows-1 - zeroCoverageWindows);
    
    
    
    
    
    
}


int loadInputFiles(void)
{
    
    inputFolder = opendir (inputDirectory.c_str());
    referencePileupDirectory = opendir (referencePileupFolder.c_str());
    if (inputFolder != NULL) {
        
        /* load all the files and directories within directory */
        
        numberOfFiles=0;
        
        while ((ent = readdir (inputFolder)) != NULL)
        {
            numberOfFiles++;
            fileInInputFolder[numberOfFiles] = ent->d_name;
            fileInReferenceFolder[numberOfFiles] = ent->d_name;
        }
        closedir (inputFolder);
    } else {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }
    return 1;
}





void generateLossesAndGainsOutput(void)
{
    
    
    
    ofstream testOutput;
    
    
    int a;
    int start;
    int end;
    int numberOfMarks;
    float segmentMean;
    float referencePosition;
    float startReference;
    float refCov;
    float segmentCoverageInReference;
    float numPositions;
    
    string header;
    char chromo[200];
    char previousChromo[200];
    char field[200];
    char chromoReference[1000];
    char refCovChar[1000];
    char referencePositionChar[100];
    
    
    
    //testOutput.open("Test_output.txt");
    strcpy(filename,outputFolder.c_str());
    strcat(filename,"cnv.txt");
    cnvFile.open(filename);
    strcpy(filename,outputFolder.c_str());
    strcat(filename,"cnv_smoothed.txt");
    cnvSmoothedFile.open(filename);
    referenceCoverageFile.close();
    strcpy(filename,outputFolder.c_str());
    strcat(filename,"referenceCoverageFile");
    referenceCoverageFileRead.open(filename);
    
    
    
    
    
    //write losses and gain in each chromosome folder
    //read header from cnv file
    getline(cnvFile,header);
    getline(referenceCoverageFileRead,header);
    
    //get the first line
    
    for(a=0;a<2;a++)
    {
        cnvFile.get(field,100,' ');
        cnvFile.get();
    }
    
    cnvFile.get();
    cnvFile.get(chromo,100,'"');
    cnvFile.get();
    cnvFile.get();
    cnvFile.get(field,100,' ');
    cnvFile.get();
    start=atoi(field);
    cnvFile.get(field,100,' ');
    cnvFile.get();
    end=atoi(field);
    cnvFile.get(field,100,' ');
    cnvFile.get();
    numberOfMarks = atoi(field);
    cnvFile.get(field,100,'\n');
    cnvFile.get();
    segmentMean = atof(field);
    strcpy(previousChromo, chromo);
    
    
    
    
    //Calculate reference coverage within the retrieved segment
    referenceCoverageFileRead.seekg(0, ios::beg);
    getline(referenceCoverageFileRead,header);
    while(!referenceCoverageFileRead.eof())
    {
        referenceCoverageFileRead >> referencePosition;
        referenceCoverageFileRead >> chromoReference;
        referenceCoverageFileRead >> refCov;
        
        
        if(strcmp(chromoReference,chromo)==0)
        {
            
            //look for start position in reference coverage file
            while(referencePosition < start)
            {
                referenceCoverageFileRead >> referencePosition;
                referenceCoverageFileRead >> chromoReference;
                referenceCoverageFileRead >> refCov;
                if (referenceCoverageFileRead.eof())
                {
                    segmentCoverageInReference = 0;
                    break;
                }
            }
            numPositions = 0;
            segmentCoverageInReference = 0;
            while(referencePosition < end)
            {
                numPositions++;
                referenceCoverageFileRead >> referencePosition;
                referenceCoverageFileRead >> chromoReference;
                referenceCoverageFileRead >> refCov;
                segmentCoverageInReference = segmentCoverageInReference + refCov;
                if (referenceCoverageFileRead.eof())
                {
                    segmentCoverageInReference = 0;
                    break;
                }
            }
            segmentCoverageInReference = (segmentCoverageInReference/numPositions);
            break;
        }
        
    }
    
    //testOutput << chromo << "  " << start << " SegMean: " << "SegMean: " << segmentMean << "first: " << segmentCoverageInReference / referenceAverageCoverage + 0.5 << " Second: " << subjectAverageCoverage/segmentCoverageInReference << endl;
    segmentMean = segmentMean * (referenceAverageCoverage / subjectAverageCoverage);
    if((segmentMean)>= ( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage)) ) )
    {
        gains << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
        ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
								if (generateGEFiles==1)
                                    gains_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
    }
    
    
    if((segmentMean)<=  ( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage)) ) )
    {
        losses << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
        ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
								if (generateGEFiles==1)
                                    losses_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end +windowSize - start) << endl;
    }
    
    
    
    //Scan the rest of the file
    cout << "Scanning the segments file for Losses/Gains calling....." << endl;
    while(!cnvFile.eof())
    {
        
        for(a=0;a<2;a++)
        {
            cnvFile.get(field,100,' ');
            cnvFile.get();
        }
        if (cnvFile.eof()) break;
        cnvFile.get();
        cnvFile.get(chromo,100,'"');
        cnvFile.get();
        cnvFile.get();
        cnvFile.get(field,100,' ');
        cnvFile.get();
        start=atoi(field);
        cnvFile.get(field,100,' ');
        cnvFile.get();
        end=atoi(field);
        cnvFile.get(field,100,' ');
        cnvFile.get();
        numberOfMarks = atoi(field);
        cnvFile.get(field,100,'\n');
        cnvFile.get();
        segmentMean = atof(field);
        
        cout << "Retrieving coverage data for loss/gain in chromosome " << chromo << " in range " << start << "-" << end<<  endl;
        
        
        //Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
        referenceCoverageFileRead.seekg(0, ios::beg);
        getline(referenceCoverageFileRead,header);
        while(!referenceCoverageFileRead.eof())
        {
            referenceCoverageFileRead >> referencePosition;
            referenceCoverageFileRead >> chromoReference;
            referenceCoverageFileRead >> refCov;
            
            //look for chromosome in the reference coverage file
            
            if(strcmp(chromoReference,chromo)==0)
            {
                //cout << "trovato chromo " << chromoReference << endl;
                //look for start position in reference coverage file
                while(referencePosition < start)
                {
                    referenceCoverageFileRead >> referencePosition;
                    referenceCoverageFileRead >> chromoReference;
                    referenceCoverageFileRead >> refCov;
                    if (referenceCoverageFileRead.eof())
                    {
                        segmentCoverageInReference = 0;
                        break;
                    }
                }
                numPositions = 0;
                segmentCoverageInReference = 0;
                while(referencePosition < end)
                {
                    numPositions++;
                    referenceCoverageFileRead >> referencePosition;
                    referenceCoverageFileRead >> chromoReference;
                    referenceCoverageFileRead >> refCov;
                    segmentCoverageInReference = segmentCoverageInReference + refCov;
                    //cout << "refCov " << refCov << endl;
                    if (referenceCoverageFileRead.eof())
                    {
                        segmentCoverageInReference = 0;
                        break;
                    }
                }
                segmentCoverageInReference = segmentCoverageInReference/numPositions;
                break;
            }
            
        }
        
        
        if(!cnvFile.eof())
        {
            /*if(strcmp(previousChromo,chromo) != 0)
             {
             losses.close();
             gains.close();
             losses.open(filename_Losses);
             gains.open(filename_Gains);
             strcpy(previousChromo,chromo);
             if(generateGEFiles==1)
             {
             losses_GE.close();
             gains_GE.close();
             losses_GE.open(filename_Losses_GE);
             gains_GE.open(filename_Gains_GE);
             }
             }*/
            
            //testOutput << chromo << "  " << start << " SegMean: " << segmentMean << " first: " << segmentCoverageInReference / referenceAverageCoverage + 0.5 << " Second: " << subjectAverageCoverage/segmentCoverageInReference << endl;
            segmentMean = segmentMean * (referenceAverageCoverage / subjectAverageCoverage);
            if((segmentMean)>=  ( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage)) ) )
            {
                gains << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
                if (generateGEFiles==1)
                    gains_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end +windowSize) << "\t" << (int)(end +windowSize - start) << endl;
            }
            
            
            if((segmentMean)<=  ( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage)) ) )
            {
                losses << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
                if (generateGEFiles==1)
                    losses_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end + windowSize - start) << endl;
            }
        }
        
    }
    
    
    
    
    
    
    //write losses_Smoothed and gain from SMOOTHED file in each chromosome folder
    //read header from cnv file
    
    
    
    
    getline(cnvSmoothedFile, header);
    
    
    //get the first line
    
    for(a=0;a<2;a++)
    {
        cnvSmoothedFile.get(field,100,' ');
        cnvSmoothedFile.get();
    }
    
    cnvSmoothedFile.get();
    cnvSmoothedFile.get(chromo,100,'"');
    
    cnvSmoothedFile.get();
    cnvSmoothedFile.get();
    cnvSmoothedFile.get(field,100,' ');
    cnvSmoothedFile.get();
    start=atoi(field);
    cnvSmoothedFile.get(field,100,' ');
    cnvSmoothedFile.get();
    end=atoi(field);
    cnvSmoothedFile.get(field,100,' ');
    cnvSmoothedFile.get();
    numberOfMarks = atoi(field);
    cnvSmoothedFile.get(field,100,'\n');
    cnvSmoothedFile.get();
    segmentMean = atof(field);
    
				//Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
				referenceCoverageFileRead.seekg(0, ios::beg);
				getline(referenceCoverageFileRead,header);
				while(!referenceCoverageFileRead.eof())
                {
                    referenceCoverageFileRead >> referencePosition;
                    referenceCoverageFileRead >> chromoReference;
                    referenceCoverageFileRead >> refCov;
                    
                    //look for chromosome in the reference coverage file
                    
                    if(strcmp(chromoReference,chromo)==0)
                    {
                        
                        //look for start position in reference coverage file
                        while(referencePosition < start)
                        {
                            referenceCoverageFileRead >> referencePosition;
                            referenceCoverageFileRead >> chromoReference;
                            referenceCoverageFileRead >> refCov;
                            if (referenceCoverageFileRead.eof())
                            {
                                segmentCoverageInReference = 0;
                                break;
                            }
                        }
                        numPositions = 0;
                        segmentCoverageInReference = 0;
                        while(referencePosition < end)
                        {
                            numPositions++;
                            referenceCoverageFileRead >> referencePosition;
                            referenceCoverageFileRead >> chromoReference;
                            referenceCoverageFileRead >> refCov;
                            segmentCoverageInReference = segmentCoverageInReference + refCov;
                            //cout << "refCov " << refCov << endl;
                            if (referenceCoverageFileRead.eof())
                            {
                                segmentCoverageInReference = 0;
                                break;
                            }
                        }
                        segmentCoverageInReference = segmentCoverageInReference/numPositions;
                        break;
                    }
                    
                }
    
    segmentMean = segmentMean * (referenceAverageCoverage / subjectAverageCoverage);
    if((segmentMean)>=  ( (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage)) ) )
    {
        gains_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
        ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
        if (generateGEFiles==1)
            gains_Smoothed_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
    }
    
    
    if((segmentMean)<=  ( (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage)) ) )
    {
        losses_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
        ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
        if (generateGEFiles==1)
            losses_Smoothed_GE << "UndefinedScaffold" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
    }
    
    
    //Scan the rest of the file
    while(!cnvSmoothedFile.eof())
    {
        for(a=0;a<2;a++)
        {
            cnvSmoothedFile.get(field,100,' ');
            cnvSmoothedFile.get();
        }
        if (cnvSmoothedFile.eof()) break;
        cnvSmoothedFile.get();
        cnvSmoothedFile.get(chromo,100,'"');
        cnvSmoothedFile.get();
        cnvSmoothedFile.get();
        cnvSmoothedFile.get(field,100,' ');
        cnvSmoothedFile.get();
        start=atoi(field);
        cnvSmoothedFile.get(field,100,' ');
        cnvSmoothedFile.get();
        end=atoi(field);
        cnvSmoothedFile.get(field,100,' ');
        cnvSmoothedFile.get();
        numberOfMarks = atoi(field);
        cnvSmoothedFile.get(field,100,'\n');
        cnvSmoothedFile.get();
        segmentMean = atof(field);
        
        cout << "Retrieving coverage data for loss/gain in chromosome " << chromo << " in range " << start << "-" << end<< " in smoothed dataset" << endl;
        
        //Calculate reference coverage within the retrieved segment (normalized over the reference average coverage)
        referenceCoverageFileRead.seekg(0, ios::beg);
        getline(referenceCoverageFileRead,header);
        while(!referenceCoverageFileRead.eof())
        {
            referenceCoverageFileRead >> referencePosition;
            referenceCoverageFileRead >> chromoReference;
            referenceCoverageFileRead >> refCov;
            
            //look for chromosome in the reference coverage file
            
            if(strcmp(chromoReference,chromo)==0)
            {
                
                //look for start position in reference coverage file
                while(referencePosition < start)
                {
                    referenceCoverageFileRead >> referencePosition;
                    referenceCoverageFileRead >> chromoReference;
                    referenceCoverageFileRead >> refCov;
                    if (referenceCoverageFileRead.eof())
                    {
                        segmentCoverageInReference = 0;
                        break;
                    }
                }
                numPositions = 0;
                segmentCoverageInReference = 0;
                while(referencePosition < end)
                {
                    numPositions++;
                    referenceCoverageFileRead >> referencePosition;
                    referenceCoverageFileRead >> chromoReference;
                    referenceCoverageFileRead >> refCov;
                    segmentCoverageInReference = segmentCoverageInReference + refCov;
                    //cout << "refCov " << refCov << endl;
                    if (referenceCoverageFileRead.eof())
                    {
                        segmentCoverageInReference = 0;
                        break;
                    }
                }
                segmentCoverageInReference = segmentCoverageInReference/numPositions;
                break;
            }
            
        }
        if(!cnvSmoothedFile.eof())
        {
            
            /*if(strcmp(previousChromo,chromo) != 0)
             {
             losses_Smoothed.close();
             gains_Smoothed.close();
             losses_Smoothed.open(filename_Smoothed_Losses);
             gains_Smoothed.open(filename_Smoothed_Gains);
             strcpy(previousChromo,chromo);
             if(generateGEFiles==1)
             {
             losses_Smoothed_GE.close();
             gains_Smoothed_GE.close();
             losses_Smoothed_GE.open(filename_Smoothed_Losses_GE);
             gains_Smoothed_GE.open(filename_Smoothed_Gains_GE);
             }
             }*/
            
            segmentMean = segmentMean * (referenceAverageCoverage / subjectAverageCoverage);
            if(segmentMean>=   (segmentCoverageInReference / referenceAverageCoverage + 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage))  )
            {
                gains_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
                if (generateGEFiles==1)
                    gains_Smoothed_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
            }
            
            
            if(segmentMean<=   (segmentCoverageInReference / referenceAverageCoverage - 0.5)*(1/(segmentCoverageInReference / referenceAverageCoverage))  )
            {
                losses_Smoothed << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << "\t" << segmentMean  << "\t" <<
                ( segmentCoverageInReference/ referenceAverageCoverage) << "\t" << ((segmentCoverageInReference/referenceAverageCoverage)*segmentMean) << endl;
                if (generateGEFiles==1)
                    losses_Smoothed_GE << "UndefinedScaffold\t" << chromo << "\t" << start << "\t" << (int)(end+windowSize) << "\t" << (int)(end+windowSize - start) << endl;
            }
        }
        
    }
    
    gains.close();
    losses.close();
    gains_Smoothed.close();
    losses_Smoothed.close();
    if(generateGEFiles==1)
    {
        gains_GE.close();
        losses_GE.close();
        gains_Smoothed_GE.close();
        losses_Smoothed_GE.close();
    }
    refineTables();
}

void refineTables(void)
{
    
    string chromo_pre;
    int pos1_pre;
    int pos2_pre;
    double cov_pre;
    string chromo;
    int pos1;
    int pos2;
    double cov;
    string str;
    double max_cop;
    int initial_pos1;
    double refCov;
    double subCov;
    int length;
    int length_pre;
    int max_length;
    double cop_pre;
    double maxRefCov;
    double maxSubCov;
    
    
    string UndefinedScaffold;
    
    int windowSize;
    
    ifstream input;
    ofstream output;
    
    
    string fnames[] = {outputFolder+"Gains.txt",outputFolder+"Losses.txt",outputFolder+"Gains_smoothed.txt",outputFolder+"Losses_smoothed.txt",outputFolder+"Gains_GE.txt",outputFolder+"Losses_GE.txt",outputFolder+"Gains_smoothed_GE.txt",outputFolder+"Losses_smoothed_GE.txt"};
    
    
    //cout << "Valori " << distanteToCollapse << " " << winSize << endl;
    //getchar();
    
    for (a=0;a<4;a++)
    {
    begin:
        input.open(fnames[a].c_str());
        output.open((outputFolder+"refinedTable").c_str());
        output << "Chromosome\tStart\tEnd\tCovInSub\tCovInRef\tCopyRatio" << endl;
        
        //get fields names
        getline(input,str);
        
        getline(input,chromo,'\t');
        getline(input,str,'\t');
        pos1 = atoi(str.c_str());
        getline(input,str,'\t');
        pos2 = atoi(str.c_str());
        getline(input,str,'\t');
        length = atoi(str.c_str());
        getline(input,str,'\t');
        
        getline(input,str,'\t');
        refCov = atof(str.c_str());
        getline(input,str,'\n');
        subCov = atof(str.c_str());
        
        if(input.eof())
        {
            a++;
            input.close();
            output.close();
            if (a<4)goto begin;
            else break;
        }
        
        max_cop = subCov / refCov;
        maxSubCov = subCov;
        maxRefCov = refCov;
        initial_pos1 = pos1;
        pos2_pre = pos1;
        chromo_pre = chromo;
        max_length = length;
        
        
        while(!input.eof())
        {
            
            while(( (pos1 - distanteToCollapse*winSize) <= pos2_pre) && (chromo==chromo_pre))
            {
                cout << "pos1 " << pos1 << "pos1_mod " << (pos1 - distanteToCollapse*(int)windowSize) << "pos1_mod_prec " << (pos1 - distanteToCollapse*windowSize) << endl;
                pos1_pre = pos1;
                pos2_pre = pos2;
                chromo_pre = chromo;
                cop_pre = subCov / refCov;
                length_pre = length;
                
                
                getline(input,chromo,'\t');
                getline(input,str,'\t');
                pos1 = atoi(str.c_str());
                getline(input,str,'\t');
                pos2 = atoi(str.c_str());
                getline(input,str,'\t');
                length = atoi(str.c_str());
                getline(input,str,'\t');
                
                getline(input,str,'\t');
                refCov = atof(str.c_str());
                getline(input,str,'\n');
                subCov = atof(str.c_str());
                
                if(input.eof()) break;
                
                if ( ((pos1 - distanteToCollapse*winSize) <= pos2_pre) && (length >= max_length) && (chromo==chromo_pre))
                {
                    max_cop = subCov / refCov;
                    maxSubCov = subCov;
                    maxRefCov = refCov;
                    max_length = length;
                }
                //if (pos1 == pos2_pre)   length += length_pre;
            }
            
            if ((pos2_pre - initial_pos1) >= minimumAreaToOutput) cout << chromo_pre << "\t" << initial_pos1 << "\t" << (int)pos2_pre << "\t" << (pos2_pre - initial_pos1)  << "\t" << maxSubCov << "\t" << maxRefCov << "\t" << max_cop << endl;
            //getchar();
            output << chromo_pre << "\t" << initial_pos1 << "\t" << (int)pos2_pre << "\t" << (pos2_pre - initial_pos1)  << "\t" << maxSubCov << "\t" << maxRefCov << "\t" << max_cop << endl;
            pos2_pre = pos1;
            initial_pos1 = pos1;
            max_cop = subCov / refCov;
            maxSubCov = subCov;
            maxRefCov = refCov;
            length_pre = length;
            chromo_pre = chromo;
            max_length = length;
            
            
            
        }
        input.close();
        output.close();
        //replace original table
        strcpy(command,"mv ");
        strcat(command,(outputFolder+"refinedTable ").c_str() );
        strcat(command,fnames[a].c_str());
        system(command);
        
        
        
    }
    
    
    if(generateGEFiles==1)
    {
        for(a=4;a<8;a++)
        {
        begin2:
            
            input.open(fnames[a].c_str());
            
            output.open((outputFolder+"refinedTable").c_str());
            
            //get fields names
            getline(input,str);
            
            getline(input,UndefinedScaffold,'\t');
            getline(input,chromo,'\t');
            getline(input,str,'\t');
            pos1 = atoi(str.c_str());
            getline(input,str,'\t');
            pos2 = atoi(str.c_str());
            getline(input,str,'\n');
            length = atoi(str.c_str());
            if(input.eof())
            {
                a++;
                input.close();
                output.close();
                if(a<8)goto begin2;
                else break;
            }
            
            initial_pos1 = pos1;
            pos2_pre = pos1;
            chromo_pre = chromo;
            
            
            while(!input.eof())
            {
                
                while(( (pos1 - distanteToCollapse*winSize) <= pos2_pre) && (chromo==chromo_pre))
                {
                    
                    
                    
                    pos1_pre = pos1;
                    pos2_pre = pos2;
                    chromo_pre = chromo;
                    length_pre = length;
                    
                    getline(input,UndefinedScaffold,'\t');
                    getline(input,chromo,'\t');
                    getline(input,str,'\t');
                    pos1 = atoi(str.c_str());
                    getline(input,str,'\t');
                    pos2 = atoi(str.c_str());
                    getline(input,str,'\n');
                    length = atoi(str.c_str());
                    if(input.eof()) break;
                    
                    
                    
                    
                }
                
                //cout << chromo_pre << "\t" << (int)initial_pos1 << "\t" << (int)pos2_pre << "\t" << ((int)pos2_pre - (int)initial_pos1)  << "\t" << maxSubCov << "\t" << maxRefCov << "\t" << max_cop << endl;
                output << "UndefinedScaffold" << "\t" <<  chromo_pre << "\t" << initial_pos1 << "\t" << pos2_pre << "\t" << (pos2_pre - initial_pos1)  << endl;
                
                
                pos2_pre = pos1;
                initial_pos1 = pos1;
                length_pre = length;
                chromo_pre = chromo;
                
            }
            
            input.close();
            output.close();
            //replace original table
            strcpy(command,"mv ");
            strcat(command,(outputFolder+"refinedTable ").c_str() );
            strcat(command,fnames[a].c_str());
            system(command);
            
        }
    }
    
    
    
    
    
    
}
