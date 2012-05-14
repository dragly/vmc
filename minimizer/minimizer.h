#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <sstream>
#include <iomanip>
using namespace std;

#include "../config.h"

/*!
  Runs with multiple parameters and attempts to find the lowest possible value of the energy by adjusting the parameters.
  */

class ini;

class Minimizer
{
public:
    Minimizer(Config *config);
    virtual void runMinimizer() = 0;
    void runBlocking();

    virtual void loadConfiguration(ini *settings) {
        (void)settings;
    }
    void writeBlockData();
    void blocking(double *values, int nValues, int blockSize, double *result);
    double mean(double *values, double nValues);
protected:
    Config *m_config;
    int m_rank;
    int m_nProcesses;
    int m_nParticles;
    int m_nDimensions;
    // monte carlo settings
    int m_nCycles;
    // energies to store in block files
    double *m_allEnergies;
private:
};

#endif // MINIMIZER_H
