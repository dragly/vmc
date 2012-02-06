// Variational Monte Carlo for atoms with up to two electrons
//#include <QSettings>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

#include "inih/cpp/INIReader.h"

#include "mainapplication.h"

#include "wavestandard.h"
#include "wavesimple.h"
#include "montecarlostandard.h"
#include "matrix.h"
#include "hamiltonianstandard.h"
#include "hamiltoniansimple.h"
#include "minimizerstandard.h"
#include "minimizer.h"
using namespace  std;

MainApplication::MainApplication(int* argc, char*** argv) :
    argc(argc),
    argv(argv)
{
    //  MPI initializations
    MPI_Init (argc, argv);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    loadConfiguration();
}

void MainApplication::loadConfiguration()
{
    settings = new INIReader("config.ini");
    if(settings->ParseError()) {
        cerr << "Warning: " << __PRETTY_FUNCTION__ << ": Could not load configuration file. Does it exist?" << endl;
    }

    cout << "MainApplication::loadConfiguration(): Starting minimizer..." << endl;
    minimizer = new MinimizerStandard(rank, nProcesses);
    minimizer->loadConfiguration(settings);
}

int MainApplication::runApplication()
{
    cout << "MinimizerQObject::runMinimizer(): called" << endl;
    minimizer->runMinimizer();

    // End MPI
    MPI_Finalize ();

    return 0;
}
