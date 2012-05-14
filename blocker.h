#ifndef BLOCKER_H
#define BLOCKER_H

class ini;

class Blocker
{
public:
    Blocker();
    void runBlocking();
    double mean(double *values, double nValues);
    void blocking(double *values, int nValues, int blockSize, double *result);
    void loadConfiguration(ini *settings);
private:
    int m_nProcesses;
    ini *m_settings;
};

#endif // BLOCKER_H
