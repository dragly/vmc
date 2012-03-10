// Stat stuff
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
// math stuff
#include <math.h>
// stream stuff
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#include "blocker.h"
#include "inih/cpp/INIReader.h"

//TODO make sure that this class actually does what it is supposed to

Blocker::Blocker()
{
}

void Blocker::loadConfiguration(INIReader *settings)
{
    m_settings = settings;
    m_nProcesses = atoi(settings->Get("Blocker","nProcesses", "1").c_str());
}

void Blocker::runBlocking() {
    // load block data
    struct stat result;
    int nLocal = 0;
    int nBlockData = 0;
    if(stat("blocks_rank0.dat", &result) == 0) {
        nLocal = result.st_size / sizeof(double);
        cout << nLocal << endl;
        nBlockData = nLocal * m_nProcesses;
    } else {
        cerr << "Trouble loading blocks_rank0.dat" << endl;
        exit(99);
    }

    double *mcResults = new double[nBlockData];
    for(int i = 0; i < m_nProcesses; i++) {
        ostringstream ost;
        ost << "blocks_rank" << i << ".dat";
        ifstream infile;
        infile.open(ost.str().c_str(), ios::in | ios::binary);
        infile.read((char*)&(mcResults[i*nLocal]), result.st_size);
        infile.close();
    }
    int nBlockSamples = 41;
    int minBlockSize = nBlockData / 50;
    int maxBlockSize = nBlockData / 10;
    int blockStepSize = (maxBlockSize - minBlockSize) / (nBlockSamples - 1);
    int blockSize = -1;
    double meanSigma[2];
    for(int i = 0; i < nBlockSamples; i++) {
        blockSize = minBlockSize + i * blockStepSize;
        blocking(mcResults, nBlockData, blockSize, meanSigma);
        double mean = meanSigma[0];
        double sigma = meanSigma[1];
        printf("BlockSize: %d    Mean: %.6f     Sigma: %.6f\n", blockSize, mean, sqrt(sigma/((nBlockData/blockSize)-1)));
    }
}

double Blocker::mean(double* values, double nValues) {
    double myMean = 0;
    for(int i = 0; i < nValues; i++) {
        myMean += values[i];
    }
    return myMean / nValues;
}

void Blocker::blocking(double *values, int nValues, int blockSize, double *result) {
    int nBlocks = nValues / blockSize;
    double totalBlockSum = 0;
    double totalBlockSumSquared = 0;
    for(int i = 0; i < nBlocks; i++) {
        double blockSum = 0;
        for(int j = 0; j < nValues; j++) {
            double val = values[i * blockSize + j];
            blockSum += val;
        }
        totalBlockSum += (blockSum / nValues);
        totalBlockSumSquared += (blockSum / nValues) * (blockSum / nValues);
    }
    result[0] = totalBlockSum / nBlocks;
    result[1] = (totalBlockSumSquared / nBlocks - result[0] * result[0]);
}
