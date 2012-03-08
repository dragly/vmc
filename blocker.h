#ifndef BLOCKER_H
#define BLOCKER_H

class Blocker
{
public:
    Blocker();
    void runBlocking();
    double mean(double *values, double nValues);
    void blocking(double *values, int nValues, int blockSize, double *result);
private:
    int m_rank;
    int m_nProcesses;
};

#endif // BLOCKER_H
