#ifndef BLOCKER_H
#define BLOCKER_H

class INIParser;
class Config;

#include <string>

class Blocker
{
public:
    Blocker(Config *config);
    void runBlocking();
    void blocking(double *values, int nValues, int blockSize, double *result);
    void loadConfiguration(INIParser *settings);
private:
    int nProcesses;
    int nBlockSamples;
    int minBlockSize;
    int maxBlockSize;
    std::string scratchDir;
    Config *config;
};

#endif // BLOCKER_H
