#ifndef MINIMIZERSTANDARD_H
#define MINIMIZERSTANDARD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "minimizer.h"

class QSettings;

using namespace std;

class MinimizerStandard : public Minimizer
{
public:
    MinimizerStandard(int rank, int nProcesses);
    void runMinimizer();
    void loadConfiguration(QSettings *settings);
private:
    // output file as global variable
    ofstream ofile, blockofile;

    //  Here we define global variables  used in various functions
    //  These can be changed by reading from file the different parameters
    int my_rank;
    int numprocs;  //  these are the parameters used by MPI  to define which node and how many
    int dimension; // three-dimensional system
    double charge;  //  we fix the charge to be that of the helium atom
    double stepLength;  //  we fix the brute force jump to 1 Bohr radius
    int nParticles;  //  we fix also the number of electrons to be 2
    int nCycles;
    int maxVariations;

    QSettings *settings;
};

#endif // MINIMIZERSTANDARD_H
