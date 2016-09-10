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
#include <stdlib.h>
#include <string.h>



using namespace std;

//Variable declarations
string readsFolder;
string referenceFile;
string outputFolder;
string outputFile;
string numThreads;
string editDistance;
string additionalBwaFlags;
string additionalPileupFlags;
string file1;
string firstRead[1000];
string secondRead[1000];
string firstReadInPair[1000];
string secondReadInPair[1000];
string singleRead[1000];
string samFileName[1000];
string minimumAlnQual;
string minimumBaseQual;
ifstream temp;
ofstream logFile;


char command[1000];
char filename[1000];




int readsPaired;
int numberOfFirstRead = 0;
int numberOfSecondRead = 0;
int numberOfFilesInPair = 0;
int numberOfFilesInSingle = 0;
int numberOfSamFiles = 0;
int matePresent;
int a,b;

void checkReadsInFolder(void);
void alignReads(void);
void convertToBam(void);
void sortBam(void);
void mergeBam(void);
void indexReference(void);
void pileupReads(void);
void createFolders(void);
void moveFilesToFolders(void);


int main(int argc, char** argv) {

    readsFolder = argv[1];
    referenceFile = argv[2];
    outputFolder = argv[3];
    outputFile = argv[4];
    editDistance = argv[5];
    numThreads = argv[6];
    additionalBwaFlags = argv[7];
    readsPaired = atoi(argv[8]);
    minimumAlnQual = argv[9];
    minimumBaseQual = argv[10];
    additionalPileupFlags = argv[11];
    

    createFolders();
    indexReference();
    checkReadsInFolder();
    alignReads();
    convertToBam();
    mergeBam();
    sortBam();
    pileupReads();
    moveFilesToFolders();
    
    cout << "The alignment successfully terminated" << endl;
    cout << "Thanks for using Altools" << endl;
    cout << "Press any key to close this window" << endl;
    getchar();
    
    
}


void moveFilesToFolders(void)
{
    strcpy(command,"mv ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.sam ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/sam/");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"rm ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.sai");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mv ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.bam ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/bam/");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mv ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/*.bam ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/bam/");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mv ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,"*_pileup* ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/pileup/");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mv ");
    strcat(command,referenceFile.c_str());
    strcat(command,".* ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/indexedReference");
    logFile << command << endl;
    system(command);
    
    
}

void createFolders(void)
{
    logFile.open((outputFolder+"/logfile.txt").c_str());

    
    strcpy(command,"mkdir -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/sam");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mkdir  -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/bam");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mkdir -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/pileup");
    logFile << command << endl;
    system(command);
    
    strcpy(command,"mkdir -p ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/indexedReference");
    logFile << command << endl;
    system(command);
}

void pileupReads(void)
{
    cout << "Indexing reference for pileup " <<endl;
    cout.flush();
    strcpy(command,"./samtools faidx ");
    strcat(command,referenceFile.c_str());
    cout << "Executing: " << endl;
    cout << command << endl;
    cout.flush();
    logFile << command << endl;
    system(command);
    
    
    cout << "Performing pileup" << endl;
    strcpy(command,"./samtools mpileup -f ");
    strcat(command,referenceFile.c_str());
    strcat(command," -q ");
    strcat(command,minimumAlnQual.c_str());
    strcat(command," -Q ");
    strcat(command,minimumBaseQual.c_str());
    strcat(command," ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,outputFile.c_str());
    strcat(command,"_sorted.bam > ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,outputFile.c_str());
    strcat(command,"_pileup");
    cout << "Executing: " << endl;
    cout << command << endl;
    cout.flush();
    logFile << command << endl;
    system(command);
    
    
    //Correcting pileup
    cout << "Elaborating pileup " << endl;
    strcpy(command,"awk '{if($4!=0) print}' ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,outputFile.c_str());
    strcat(command,"_pileup >");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,"correctedPileup");
    cout << "Executing: " << endl;
    cout << command << endl;
    cout.flush();
    logFile << command << endl;
    system(command);

    strcpy(command,"./pileupCorr ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,"correctedPileup ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,outputFile.c_str());
    strcat(command,"_pileup");
    logFile << command << endl;
    system(command);
    strcpy(command,"rm ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,"correctedPileup ");
    cout << "Executing: " << endl;
    cout << command << endl;
    cout.flush();
    logFile << command << endl;
    system(command);

    
}

void indexReference(void)
{
    cout << "Indexing reference file....." << endl;
    cout.flush();
    strcpy(command,"./bwa index ");
    strcat(command,referenceFile.c_str());
    cout << "Executing: " << endl;
    cout << command << endl;
    cout.flush();
    logFile << command << endl;
    system(command);
    cout << "Reference successfullly indexed...." << endl;
    cout.flush();

}


void checkReadsInFolder(void)
{
    //Check the reads folder
    strcpy(command,"ls ");
    strcat(command,readsFolder.c_str());
    strcat(command,"/*.fastq >");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,"readsFiles");
    system(command);
    strcpy(filename,outputFolder.c_str());
    strcat(filename,"/");
    strcat(filename,"readsFiles");
    temp.open(filename);
    while(!temp.eof())
    {
        getline(temp,file1,'\n');
        if(!temp.eof())
        {
            //cout << file1.substr(file1.size()-8) << endl;
            if (file1.substr(file1.size()-8) == "_1.fastq")
            {
                numberOfFirstRead++;
                firstRead[numberOfFirstRead] = file1.substr(0,file1.size()-8);
            }
            if (file1.substr(file1.size()-8)== "_2.fastq")
            {
                numberOfSecondRead++;
                secondRead[numberOfSecondRead] = file1.substr(0,file1.size()-8);
            }
        }
        
    }
    
    
    
    
    //Split available reads between paired end and single end
    for (a=1;a<=numberOfFirstRead;a++)
    {
        matePresent = 0;
        for(b=1;b<=numberOfSecondRead;b++)
        {
            if (firstRead[a] == secondRead[b])
            {
                matePresent = 1;
                break;
            }
        }
        if (matePresent==1)
        {
            numberOfFilesInPair++;
            firstReadInPair[numberOfFilesInPair]=firstRead[a];
            secondReadInPair[numberOfFilesInPair]=secondRead[b];
        }
        else
        {
            numberOfFilesInSingle++;
            singleRead[numberOfFilesInSingle]= firstRead[a];
        }
        
    }
    
    temp.close();
}

void sortBam(void)
{
    cout << "Sorting bam file " << endl;
    cout.flush();
    
    strcpy(command,"./samtools sort ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,outputFile.c_str());
    strcat(command,".bam ");
    strcat(command,outputFolder.c_str());
    strcat(command,"/");
    strcat(command,outputFile.c_str());
    strcat(command,"_sorted");
    cout << "Executing: " << endl;
    cout << command << endl;
    cout.flush();
    logFile << command << endl;
    system(command);
 
}

void convertToBam(void)
{
    cout << "Converting sam file to bam format....." << endl;
    for(a=1;a<=numberOfSamFiles;a++)
    {
        strcpy(command,"./samtools view -bS ");
        strcat(command,samFileName[a].c_str());
        strcat(command,".sam >");
        strcat(command,samFileName[a].c_str());
        strcat(command,".bam");
        cout << "Executing: " << endl;
        cout << command << endl;
        cout.flush();
        logFile << command << endl;
        system(command);
        
    }
}

void mergeBam(void)
{
    cout << "Merging bam files " << endl;
    cout.flush();
    
    if(numberOfSamFiles>1)
    {
        strcpy(command,"./samtools merge ");
        strcat(command,outputFolder.c_str());
        strcat(command,"/");
        strcat(command,outputFile.c_str());
        strcat(command,".bam ");
        for(a=1;a<=numberOfSamFiles;a++)
        {
            strcat(command,samFileName[a].c_str());
            strcat(command,".bam ");
        }
        cout << "Executing: " << endl;
        cout << command << endl;
        cout.flush();
        logFile << command << endl;
        system(command);
    }
    else
    {
        strcpy(command,"cp ");
        strcat(command,samFileName[1].c_str());
        strcat(command,".bam ");
        strcat(command,outputFolder.c_str());
        strcat(command,"/");
        strcat(command,outputFile.c_str());
        strcat(command,".bam ");
        cout << "Executing: " << endl;
        cout << command << endl;
        cout.flush();
        logFile << command << endl;
        system(command);
    }
    
}



void alignReads(void)
{
    
    // PERFORM THE BWA ALN COMMAND
    cout << "Aligning reads " << endl;
    cout.flush();
    for(a=1;a<=numberOfFilesInPair;a++)
    {
        strcpy(command,"./bwa aln ");
        strcat(command,referenceFile.c_str());
        strcat(command," ");
        strcat(command,firstReadInPair[a].c_str());
        strcat(command,"_1.fastq ");
        strcat(command," -t ");
        strcat(command,numThreads.c_str());
        strcat(command," -n ");
        strcat(command,editDistance.c_str());
        strcat(command," ");
        if (additionalBwaFlags!="none")
        {
            strcat(command,additionalBwaFlags.c_str());
            strcat(command," ");
        }
        
        strcat(command," > ");
        strcat(command,firstReadInPair[a].c_str());
        strcat(command,"_1.sai");
        cout << "Executing: " << endl;
        cout << command << endl;
        cout.flush();
        logFile << command << endl;
        system(command);
        
        
        
        strcpy(command,"./bwa aln ");
        strcat(command,referenceFile.c_str());
        strcat(command," ");
        strcat(command,secondReadInPair[a].c_str());
        strcat(command,"_2.fastq ");
        strcat(command," -t ");
        strcat(command,numThreads.c_str());
        strcat(command," -n ");
        strcat(command,editDistance.c_str());
        strcat(command," ");
        if (additionalBwaFlags!="none")
        {
            strcat(command,additionalBwaFlags.c_str());
            strcat(command," ");
        }
        strcat(command," > ");
        strcat(command,secondReadInPair[a].c_str());
        strcat(command,"_2.sai");
        cout << "Executing: " << endl;
        cout << command << endl;
        cout.flush();
        logFile << command << endl;
        system(command);
        
    }
    
    for(a=1;a<=numberOfFilesInSingle;a++)
    {
        strcpy(command,"./bwa aln ");
        strcat(command,referenceFile.c_str());
        strcat(command," ");
        strcat(command,singleRead[a].c_str());
        strcat(command,"_1.fastq  ");
        strcat(command," -t ");
        strcat(command,numThreads.c_str());
        strcat(command," -n ");
        strcat(command,editDistance.c_str());
        strcat(command," ");
        if (additionalBwaFlags!="none")
        {
            strcat(command,additionalBwaFlags.c_str());
            strcat(command," ");
        }
        
        strcat(command," > ");
        strcat(command,singleRead[a].c_str());
        strcat(command,"_1.sai");
        cout << "Executing: " << endl;
        cout << command << endl;
        cout.flush();
        logFile << command << endl;
        system(command);
        

    }
    
    //PERFORM THE BWA SAMPE/SAMSE COMMAND
    if (readsPaired==1)
    {
        for(a=1;a<=numberOfFilesInPair;a++)
        {
            strcpy(command,"./bwa sampe ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.sai ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.fastq ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.fastq > ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,".sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = firstReadInPair[a];
            cout << "Executing: " << endl;
            cout << command << endl;
            cout.flush();
            logFile << command << endl;
            system(command);
            
        }
        for(a=1;a<=numberOfFilesInSingle;a++)
        {
            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.fastq >");
            strcat(command,singleRead[a].c_str());
            strcat(command,".sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = singleRead[a];
            cout << "Executing: " << endl;
            cout << command << endl;
            cout.flush();
            logFile << command << endl;
            system(command);
        
        }
    }
    else
    {
        for(a=1;a<=numberOfFilesInPair;a++)
        {
            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.fastq >");
            strcat(command,firstReadInPair[a].c_str());
            strcat(command,"_1.sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = firstReadInPair[a] + "_1";
            cout << "Executing: " << endl;
            cout << command << endl;
            cout.flush();
            logFile << command << endl;
            system(command);
            

            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.sai ");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.fastq >");
            strcat(command,secondReadInPair[a].c_str());
            strcat(command,"_2.sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = secondReadInPair[a] + "_2";
            cout << "Executing: " << endl;
            cout << command << endl;
            cout.flush();
            logFile << command << endl;
            system(command);
            
        }
        for(a=1;a<=numberOfFilesInSingle;a++)
        {
            strcpy(command,"./bwa samse ");
            strcat(command,referenceFile.c_str());
            strcat(command," ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.sai ");
            strcat(command,singleRead[a].c_str());
            strcat(command,"_1.fastq >");
            strcat(command,singleRead[a].c_str());
            strcat(command,".sam");
            numberOfSamFiles++;
            samFileName[numberOfSamFiles] = singleRead[a];
            cout << "Executing: " << endl;
            cout << command << endl;
            cout.flush();
            logFile << command << endl;
            system(command);
            
        }
    }

}




