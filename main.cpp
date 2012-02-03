// Variational Monte Carlo for atoms with up to two electrons

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "wavestandard.h"
#include "wavesimple.h"
#include "montecarlostandard.h"
#include "matrix.h"
#include "hamiltonianstandard.h"
#include "hamiltoniansimple.h"
#include "minimizerstandard.h"
using namespace  std;


// Begin of main program

int main(int argc, char* argv[])
{
    //  MPI initializations
    MPI_Init (0, 0);
    int rank, nProcesses;
    MPI_Comm_size (MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    cout << "Starting minimizer..." << endl;
    MinimizerStandard *minimizer = new MinimizerStandard(rank, nProcesses);
    cout << "Running minimizer..." << endl;
    minimizer->run();
    cout << "Minimizer done." << endl;

    // End MPI
    MPI_Finalize ();

    return 0;
}  //  end of main function

