#ifndef BLOCKER_H
#define BLOCKER_H

class INIParser;

class Blocker
{
public:
    Blocker();
    void runBlocking();
    double mean(double *values, double nValues);
    void blocking(double *values, int nValues, int blockSize, double *result);
    void loadConfiguration(INIParser *settings);
private:
    int m_nProcesses;
    INIParser *m_settings;
};

#endif // BLOCKER_H
