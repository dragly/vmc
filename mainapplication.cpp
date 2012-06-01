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
#include "montecarlo/standardmontecarlo.h"
#include "matrix.h"
#include "minimizer/standardminimizer.h"
#include "minimizer/minimizer.h"
#include "blocker.h"
#include "densityplotter.h"
#include "config.h"
#include "montecarlo/diffusionmontecarlo.h"
#include "minimizer/geneticminimizer.h"
#include "onerun/onerun.h"
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
    std::cout << "Loading ini reader for myRank " << m_rank << std::endl;
    m_settings = new INIParser("config.ini");

    if(!m_settings->Good()) {
        cerr << "Warning: " << __PRETTY_FUNCTION__ << ": Could not load configuration file 'config.ini'. Does it exist?" << endl;
    }

    m_config = new Config(m_rank, m_nProcesses);
    // Make sure the config is loaded by one processor at the time
    for(int i = 0; i < m_nProcesses; i++) {
        if(i == m_config->myRank()) {
            m_config->loadConfiguration(m_settings);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    string modeString = m_settings->Get("General","mode","minimizer");
    if(modeString == "minimizer") {
        m_mode = MinimizerMode;
    } else if(modeString == "density") {
        m_mode = DensityMode;
    } else if(modeString == "blocking") {
        m_mode = BlockingMode;
    } else if(modeString == "genetic") {
        m_mode = GeneticMode;
    } else if(modeString == "onerun") {
        m_mode = OneRunMode;
    } else if(modeString == "diffusion") {
        m_mode = DiffusionMode;
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

void MainApplication::runConfiguration()
{
    loadConfiguration();
    if(m_mode == MinimizerMode) {
        Minimizer *minimizer = new StandardMinimizer(m_config);
        // Make sure the config is loaded by one processor at the time
        for(int i = 0; i < m_nProcesses; i++) {
            if(i == m_config->myRank()) {
                minimizer->loadConfiguration(m_settings);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        minimizer->runMinimizer();
        delete minimizer;
    } else if(m_mode== DensityMode) {
        DensityPlotter *densityPlotter = new DensityPlotter(m_config);
        // Make sure the config is loaded by one processor at the time
        for(int i = 0; i < m_nProcesses; i++) {
            if(i == m_config->myRank()) {
                densityPlotter->loadConfiguration(m_settings);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        densityPlotter->makePlot();
        delete densityPlotter;
    } else if(m_mode == BlockingMode) {
        Blocker* blocker = new Blocker(m_config);

        // Make sure the config is loaded by one processor at the time
        for(int i = 0; i < m_nProcesses; i++) {
            if(i == m_config->myRank()) {
                blocker->loadConfiguration(m_settings);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        blocker->runBlocking();
        delete blocker;
    } else if(m_mode == DiffusionMode) {
        DiffusionMonteCarlo *diffusionMonteCarlo = new DiffusionMonteCarlo(m_config);

        // Make sure the config is loaded by one processor at the time
        for(int i = 0; i < m_nProcesses; i++) {
            if(i == m_config->myRank()) {
                diffusionMonteCarlo->loadConfiguration(m_settings);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        diffusionMonteCarlo->sample();
        double energy = diffusionMonteCarlo->energy();
        std::cout << "Diffusion monte carlo returned energy of " << energy << std::endl;
        delete diffusionMonteCarlo;
    } else if(m_mode == GeneticMode) {
        GeneticMinimizer *geneticMinimizer = new GeneticMinimizer(m_config);
        // Make sure the config is loaded by one processor at the time
        for(int i = 0; i < m_nProcesses; i++) {
            if(i == m_config->myRank()) {
                geneticMinimizer->loadConfiguration(m_settings);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        geneticMinimizer->runMinimizer();
        delete geneticMinimizer;
    } else if(m_mode == OneRunMode) {
        OneRun *oneRun = new OneRun(m_config);
        // Make sure the config is loaded by one processor at the time
        for(int i = 0; i < m_nProcesses; i++) {
            if(i == m_config->myRank()) {
                oneRun->loadConfiguration(m_settings);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        oneRun->run();
        delete oneRun;
    } else {
        cerr << __PRETTY_FUNCTION__ << ": Unknown mode" << endl;
        exit(459);
    }
}

void MainApplication::finalize() {
    MPI_Finalize ();
}
