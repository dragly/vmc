#ifndef MINIMIZERSTANDARD_H
#define MINIMIZERSTANDARD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "minimizer.h"

using namespace std;

class MinimizerStandard : public Minimizer
{
public:
    MinimizerStandard(int rank, int nProcesses);
    void run();
private:
    // output file as global variable
    ofstream ofile, blockofile;

    //  Here we define global variables  used in various functions
    //  These can be changed by reading from file the different parameters
    int dimension; // three-dimensional system
    double charge;  //  we fix the charge to be that of the helium atom
    int my_rank, numprocs;  //  these are the parameters used by MPI  to define which node and how many
    double step_length;  //  we fix the brute force jump to 1 Bohr radius
    int number_particles;  //  we fix also the number of electrons to be 2
};

#endif // MINIMIZERSTANDARD_H
