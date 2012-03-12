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
#include "montecarlo/montecarlostandard.h"
#include "matrix.h"
#include "minimizerstandard.h"
#include "minimizer.h"
#include "blocker.h"
#include "densityplotter.h"
#include "config.h"
using namespace  std;

MainApplication::MainApplication(int* argc, char*** argv) :
    argc(argc),
    argv(argv)
{
#ifdef USE_MPI
    //  MPI initializations
    MPI_Init (argc, argv);
    MPI_Comm_size (MPI_COMM_WORLD, &m_nProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &m_rank);
#else
    m_nProcesses = 1;
    m_rank = 0;
#endif
}

void MainApplication::loadConfiguration()
{
    m_settings = new INIReader("config.ini");
    if(m_settings->ParseError()) {
        cerr << "Warning: " << __PRETTY_FUNCTION__ << ": Could not load configuration file 'config.ini'. Does it exist?" << endl;
    }
    m_nParticles = m_settings->GetInteger("General", "nParticles", 2);
    m_nDimensions = m_settings->GetInteger("General", "nDimensions", 2);

    string modeString = m_settings->Get("General","mode","minimizer");
    if(modeString == "minimizer") {
        m_mode = MinimizerMode;
    } else if(modeString == "density") {
        m_mode = DensityMode;
    } else if(modeString == "blocking") {
        m_mode = BlockingMode;
    } else {
        cerr << __PRETTY_FUNCTION__ << ": Unknown mode '" << modeString << "'" << endl;
        exit(460);
    }
    if(m_rank == 0) {
        cout << __PRETTY_FUNCTION__ << ": Config loaded. Mode is " << modeString << endl;
    }
    m_config = new Config(m_rank, m_nProcesses, m_nDimensions, m_nParticles);
}

void MainApplication::runConfiguration()
{
    if(m_mode == MinimizerMode) {
        runMinimizer();
    } else if(m_mode== DensityMode) {
        runDensity();
    } else if(m_mode == BlockingMode) {
        runBlocking();
    } else {
        cerr << __PRETTY_FUNCTION__ << ": Unknown mode" << endl;
        exit(459);
    }
}

void MainApplication::runMinimizer()
{
    Minimizer *minimizer = new MinimizerStandard(m_config);
    minimizer->loadConfiguration(m_settings);
    minimizer->runMinimizer();
}

void MainApplication::runBlocking()
{
//    if(m_rank == 0) {
    Blocker* blocker = new Blocker();
    blocker->loadConfiguration(m_settings);
    blocker->runBlocking();
//    }
}

void MainApplication::runDensity()
{
    DensityPlotter *densityPlotter = new DensityPlotter(m_config);
    densityPlotter->loadConfiguration(m_settings);
    densityPlotter->makePlot();
}

void MainApplication::finalize() {
#ifdef USE_MPI
    // End MPI
    MPI_Finalize ();
#endif
}
