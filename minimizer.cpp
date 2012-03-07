#include "minimizer.h"

Minimizer::Minimizer(int rank, int nProcesses) :
    rank(rank),
    nProcesses(nProcesses)
{
}

void Minimizer::writeBlockData() {
    // Setting output file name for this rank:
    ostringstream ost;
    ost << "blocks_rank" << rank << ".dat";
    // Open file for writing:
    blockofile.open(ost.str().c_str(), ios::out | ios::binary);
    blockofile.write((char*)(allEnergies+1),
                     nCycles*sizeof(double));
    blockofile.close();
}
