#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <sstream>
#include <iomanip>
using namespace std;

#include "../config.h"

/*!
  * \brief Runs with multiple parameters and attempts to find the lowest possible value of the energy by adjusting the parameters.
  */
class INIParser;

class Minimizer
{
public:
    Minimizer(Config *config);
    virtual void runMinimizer() = 0;
    void runBlocking();

    virtual void loadConfiguration(INIParser *settings) {
        (void)settings;
    }
    void blocking(double *values, int nValues, int blockSize, double *result);
    double mean(double *values, double nValues);
protected:
    Config *config;
    int m_rank;
    int m_m_nProcesses;
    int m_nParticles;
    int m_nDimensions;
    // monte carlo settings
    int m_nSamples;
    // energies to store in block files
    double *m_allEnergies;
private:
};

#endif // MINIMIZER_H
