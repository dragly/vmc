#ifndef BLOCKER_H
#define BLOCKER_H

class INIParser;
class Config;

#include <string>


/*!
  * This class opens a blocking file and performs statistical analysis on that file.
  */
class Blocker
{
public:
    Blocker(Config *config);
    void runBlocking();
    void blocking(double *values, int nValues, int blockSize, double *result);
    void loadConfiguration(INIParser *settings);
private:
    int m_nProcesses;
    int nBlockSamples;
    int minBlockSize;
    int maxBlockSize;
    std::string scratchDir;
    Config *config;
};

#endif // BLOCKER_H
