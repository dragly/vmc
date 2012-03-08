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


Minimizer::Minimizer(int rank, int nProcesses) :
    m_rank(rank),
    m_nProcesses(nProcesses)
{
}

void Minimizer::writeBlockData() {
    // Setting output file name for this rank:
    ostringstream ost;
    ost << "blocks_rank" << m_rank << ".dat";
    // Open file for writing:
    blockofile.open(ost.str().c_str(), ios::out | ios::binary);
    blockofile.write((char*)(allEnergies+1),
                     nCycles*sizeof(double));
    blockofile.close();
}
