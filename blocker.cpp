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
#include "inih/ini.h"
#include "config.h"

// TODO make sure that the Blocker class actually does what it is supposed to

Blocker::Blocker(Config* config_) :
    nProcesses(config_->nProcesses()),
    nBlockSamples(41),
    minBlockSize(1),
    maxBlockSize(10000),
    config(config_)
{
}

void Blocker::loadConfiguration(INIParser *settings)
{
    nBlockSamples = settings->GetDouble("Blocking", "nBlockSamples", nBlockSamples);
    minBlockSize = settings->GetDouble("Blocking", "minBlockSize", minBlockSize);
    maxBlockSize = settings->GetDouble("Blocking", "maxBlockSize", maxBlockSize);
    scratchDir = settings->GetString("Blocking", "scratchDir", "/scratch/blocking");
}

void Blocker::runBlocking() {
    // load block data
    struct stat result;
    int nLocal = 0;
    int nSamples = 0;

    ostringstream blockfile;
    ostringstream path;
    path << scratchDir << "/" << config->nParticles() << "p-omega" << config->omega();
    blockfile << path.str() << "/blocks_rank0.dat";
    if(stat(blockfile.str().c_str(), &result) == 0) {

    } else {
        cerr << "Trouble loading blocks_rank0.dat" << endl;
        exit(99);
    }

    nLocal = result.st_size / sizeof(double);
    nSamples = nLocal * nProcesses;

    double *monteCarloResults = new double[nSamples];
    for(int i = 0; i < nProcesses; i++) {
        ostringstream ost;
        std::cout << "Opening file " << i << std::endl;
        ost << path.str() << "/blocks_rank" << i << ".dat";
        ifstream infile;
        if(stat(ost.str().c_str(), &result) == 0) {
                infile.open(ost.str().c_str(), ios::in | ios::binary);
                infile.read((char*)&(monteCarloResults[i*nLocal]), result.st_size);
                infile.close();
            } else {
                std::cerr << "Could not find file " << ost.str() << std::endl;
            exit(953);
        }
    }
//    for(int i = 0; i < nSamples; i++) {
//        std::cout << monteCarloResults[i] << std::endl;
//    }
    int blockStepSize = (maxBlockSize - minBlockSize) / (nBlockSamples - 1);
    int blockSize = -1;
    double results[2];
    ofstream outfile;
    outfile.open("blockingdata.dat");
    for(int i = 0; i < nBlockSamples; i++) {
        blockSize = minBlockSize + i * blockStepSize;
        blocking(monteCarloResults, nSamples, blockSize, results);
        double mean = results[0];
        double standardError = results[1];
        outfile << blockSize << " " << mean << " " << standardError << std::endl;
    }
    outfile.close();
}

void Blocker::blocking(double *values, int nValues, int blockSize, double *result) {
    int nBlocks = nValues / blockSize;
    double totalBlockSum = 0;
    double totalBlockSumSquared = 0;
    for(int i = 0; i < nBlocks; i++) {
        double blockSum = 0;
        for(int j = 0; j < blockSize; j++) {
            double val = values[i * blockSize + j];
            blockSum += val;
        }
        totalBlockSum += (blockSum / blockSize);
        totalBlockSumSquared += (blockSum / blockSize) * (blockSum / blockSize);
    }
    result[0] = totalBlockSum / nBlocks;
    result[1] = sqrt(totalBlockSumSquared / nBlocks - result[0] * result[0]) / sqrt(nBlocks - 1);
}
