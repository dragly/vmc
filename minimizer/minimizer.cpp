// Stat stuff
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
// math stuff
#include <math.h>
// stream stuff
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// own header
#include "minimizer.h"

// local includes


Minimizer::Minimizer(Config *config) :
    config(config),
    m_rank(config->rank()),
    m_nProcesses(config->nProcesses())
{
}

void Minimizer::writeBlockData() {
    ofstream blockofile;
    // Setting output file name for this rank:
    ostringstream ost;
    ost << "blocks_rank" << m_rank << ".dat";
    // Open file for writing:
    blockofile.open(ost.str().c_str(), ios::out | ios::binary);
    blockofile.write((char*)(m_allEnergies+1),
                     m_nSamples*sizeof(double));
    blockofile.close();
}
