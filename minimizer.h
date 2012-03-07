#ifndef MINIMIZER_H
#define MINIMIZER_H
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

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

    virtual void loadConfiguration(INIReader *settings) = 0;
    void writeBlockData();
protected:
    // output file as global variable
    ofstream ofile;
    ofstream blockofile;
    int rank;
    int nProcesses;
    // monte carlo settings
    int nCycles;
    // energies to store in block files
    double *allEnergies;
private:
};

#endif // MINIMIZER_H
