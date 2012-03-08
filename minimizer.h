#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <sstream>
#include <iomanip>
using namespace std;

/*!
  Runs with multiple parameters and attempts to find the lowest possible value of the energy by adjusting the parameters.
  */

class INIReader;

class Minimizer
{
public:
    Minimizer(int rank, int nProcesses);
    virtual void runMinimizer() = 0;
    void runBlocking();

    virtual void loadConfiguration(INIReader *settings) = 0;
    void writeBlockData();
    void blocking(double *values, int nValues, int blockSize, double *result);
    double mean(double *values, double nValues);
protected:
    // output file as global variable
    ofstream ofile;
    ofstream blockofile;
    int m_rank;
    int m_nProcesses;
    // monte carlo settings
    int nCycles;
    // energies to store in block files
    double *allEnergies;
private:
};

#endif // MINIMIZER_H
