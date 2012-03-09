// Variational Monte Carlo for atoms with up to two electrons
//#include <QSettings>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#ifdef USE_MPI
// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif

#include "inih/cpp/INIReader.h"

#include "mainapplication.h"

#include "wavestandard.h"
#include "wavesimple.h"
#include "montecarlostandard.h"
#include "matrix.h"
#include "minimizerstandard.h"
#include "minimizer.h"
#include "blocker.h"
#include "densityplotter.h"
using namespace  std;

MainApplication::MainApplication(int* argc, char*** argv) :
    argc(argc),
    argv(argv)
{
#ifdef USE_MPI
    //  MPI initializations
    MPI_Init (argc, argv);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#else
    nProcesses = 1;
    rank = 0;
#endif
}

void MainApplication::loadConfiguration()
{
    settings = new INIReader("config.ini");
    if(settings->ParseError()) {
        cerr << "Warning: " << __PRETTY_FUNCTION__ << ": Could not load configuration file 'config.ini'. Does it exist?" << endl;
    }
    string modeString = settings->Get("General","mode","minimizer");
    if(modeString == "minimizer") {
        mode = MinimizerMode;
    } else if(modeString == "density") {
        mode = DensityMode;
    } else if(modeString == "blocking") {
        mode = BlockingMode;
    } else {
        cerr << __PRETTY_FUNCTION__ << ": Unknown mode '" << modeString << "'" << endl;
        exit(460);
    }
    if(rank == 0) {
        cout << __PRETTY_FUNCTION__ << ": Config loaded. Mode is " << modeString << endl;
    }
}

void MainApplication::runConfiguration()
{
    if(mode == MinimizerMode) {
        runMinimizer();
    } else if(mode== DensityMode) {
        runDensity();
    } else if(mode == BlockingMode) {
        runBlocking();
    } else {
        cerr << __PRETTY_FUNCTION__ << ": Unknown mode" << endl;
        exit(459);
    }
}

void MainApplication::runMinimizer()
{
    Minimizer *minimizer = new MinimizerStandard(rank, nProcesses);
    minimizer->loadConfiguration(settings);
    minimizer->runMinimizer();
}

void MainApplication::runBlocking()
{
    if(rank == 0) {
        Blocker* blocker = new Blocker();
        blocker->runBlocking();
    }
}

void MainApplication::runDensity()
{
    DensityPlotter *densityPlotter = new DensityPlotter(rank, nProcesses);
    densityPlotter->makePlot();
}

void MainApplication::finalize() {
#ifdef USE_MPI
    // End MPI
    MPI_Finalize ();
#endif
}
