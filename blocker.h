#ifndef BLOCKER_H
#define BLOCKER_H

class INIReader;

class Blocker
{
public:
    Blocker();
    void runBlocking();
    double mean(double *values, double nValues);
    void blocking(double *values, int nValues, int blockSize, double *result);
    void loadConfiguration(INIReader *settings);
private:
    int m_nProcesses;
    INIReader *m_settings;
};

#endif // BLOCKER_H
