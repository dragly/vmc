#ifndef MINIMIZERSTANDARD_H
#define MINIMIZERSTANDARD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "minimizer.h"

class INIReader;
class WaveFunction;
class Hamiltonian;

using namespace std;

class MinimizerStandard : public Minimizer
{
public:
    MinimizerStandard(Config *config);
    void runMinimizer();
    void loadConfiguration(INIReader *settings);
private:

    //  Here we define global variables  used in various functions
    //  These can be changed by reading from file the different parameters
    int dimension; // three-dimensional system
    double charge;  //  we fix the charge to be that of the helium atom
    double stepLength;  //  we fix the brute force jump to 1 Bohr radius
    int m_nVariations;
    string hamiltonianClass;
    WaveFunction *m_wave;

    INIReader *m_settings;
    Hamiltonian *m_hamiltonian;

};

#endif // MINIMIZERSTANDARD_H
