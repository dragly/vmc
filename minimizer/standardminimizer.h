#ifndef STANDARDMINIMIZER_H
#define STANDARDMINIMIZER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "minimizer.h"

class INIParser;
class WaveFunction;
class Hamiltonian;
class MonteCarlo;

using namespace std;

/*!
  * \brief Scans the parametric space and saves the energies to file for plotting
  */
class StandardMinimizer : public Minimizer
{
public:
    StandardMinimizer(Config *config);
    void runMinimizer();
    void loadConfiguration(INIParser *settings);
private:

    //  Here we define global variables  used in various functions
    //  These can be changed by reading from file the different parameters
    int dimension; // three-dimensional system
    int nVariations;
    double alphaStart;
    double alphaEnd;
    double betaStart;
    double betaEnd;
    INIParser *m_settings;

    WaveFunction *m_wave;
    MonteCarlo *m_monteCarlo;
    Hamiltonian *m_hamiltonian;

};

#endif // STANDARDMINIMIZER_H
