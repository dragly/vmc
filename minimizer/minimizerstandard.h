#ifndef MINIMIZERSTANDARD_H
#define MINIMIZERSTANDARD_H

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

class MinimizerStandard : public Minimizer
{
public:
    MinimizerStandard(Config *config);
    void runMinimizer();
    void loadConfiguration(INIParser *settings);
private:

    //  Here we define global variables  used in various functions
    //  These can be changed by reading from file the different parameters
    int dimension; // three-dimensional system
    double charge;  //  we fix the charge to be that of the helium atom
    double stepLength;  //  we fix the brute force jump to 1 Bohr radius
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

#endif // MINIMIZERSTANDARD_H
