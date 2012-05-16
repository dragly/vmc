// Variational Monte Carlo for atoms with up to two electrons
//#include <QSettings>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

//#include "inih/cpp/INIReader.h"
#include "inih/ini.h"

#include "mainapplication.h"

#include "wavefunction/wavesimple.h"
#include "montecarlo/montecarlostandard.h"
#include "matrix.h"
#include "minimizer/minimizerstandard.h"
#include "minimizer/minimizer.h"
#include "blocker.h"
#include "densityplotter.h"
#include "config.h"
using namespace  std;

MainApplication::MainApplication(int* argc, char*** argv) :
    argc(argc),
    argv(argv)
{
    //  MPI initializations
    MPI_Init (argc, argv);
    MPI_Comm_size (MPI_COMM_WORLD, &m_nProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &m_rank);
}

void MainApplication::loadConfiguration()
{
    // Make sure the config is loaded by one processor at the time
    for(int i = 0; i < m_nProcesses; i++) {
        if(i == m_rank) {
            std::cout << "Loading ini reader for rank " << m_rank << std::endl;
            m_settings = new INIParser("config.ini");

            if(!m_settings->Good()) {
                cerr << "Warning: " << __PRETTY_FUNCTION__ << ": Could not load configuration file 'config.ini'. Does it exist?" << endl;
            }

            m_config = new Config(m_rank, m_nProcesses);
            m_config->loadConfiguration(m_settings);

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
                cout << __PRETTY_FUNCTION__ << ": Config loaded. Mode is " << modeString << "." << endl;
                cout << __PRETTY_FUNCTION__ << ": Running with " << m_nProcesses << " process(es)." << endl;
                flush(cout);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void MainApplication::runConfiguration()
{
    loadConfiguration();
    if(m_mode == MinimizerMode) {
        runMinimizer();
    } else if(m_mode== DensityMode) {
        runDensity();
    } else if(m_mode == BlockingMode) {
        runBlocking();
    } else if(m_mode == DiffusionMode) {
        runDiffusion();
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

void MainApplication::runDiffusion()
{
    // DiffusionMonteCarlo *diffMC = new DiffusionMonteCarlo();
    // ...
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
    MPI_Finalize ();
}
