//#include <QSettings>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#ifdef USE_MPI
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif
#include <stdio.h>
#include <stdlib.h>

// local stuff
#include "minimizerstandard.h"
#include "inih/cpp/INIReader.h"

#include "wavestandard.h"
#include "wavesimple.h"
#include "montecarlostandard.h"
#include "matrix.h"
#include "hamiltonian/hamiltonianstandard.h"
#include "hamiltonian/hamiltoniansimple.h"
#include "hamiltonian/hamiltonianideal.h"
#include "waveideal.h"

using namespace std;

MinimizerStandard::MinimizerStandard(Config *config) :
    Minimizer(config)
{
}

void MinimizerStandard::loadConfiguration(INIReader *settings)
{
    m_settings = settings;
    charge = atof(settings->Get("MinimizerStandard","charge", "1.0").c_str());
    stepLength = atof(settings->Get("MinimizerStandard","stepLength", "1.0").c_str());
    m_nCycles = atoi(settings->Get("MinimizerStandard","nCycles", "1000").c_str());
    m_nVariations = settings->GetInteger("MinimizerStandard","nVariations", 11);    //  default number of variations

    // Wave properties
    string waveClass = settings->Get("Wave","class", "WaveSimple");
    m_wave = WaveFunction::fromName(waveClass, m_config);
    if(m_wave == 0) {
        cerr << "Unknown wave class '" << waveClass << "'" << endl;
        exit(99);
    }
    m_wave->loadConfiguration(m_settings);

    // Hamiltonian
    m_hamiltonian = Hamiltonian::fromName(settings->Get("Hamiltonian","class", "HamiltonianSimple"), m_config, charge);
    if(m_hamiltonian == 0) {
        cerr << "Unknown hamiltonian class '" << hamiltonianClass << "'" << endl;
        exit(98);
    }
}

void MinimizerStandard::runMinimizer()
{
    int total_number_cycles, i;
    double *cumulative_e, *cumulative_e2;
    double *total_cumulative_e, *total_cumulative_e2;
    double timeStart;
    double timeEnd;
#ifndef USE_MPI
    timeStart = timeEnd = -1;
#endif
    double totalTime;
    double variance;
    double energy;
    double error;

    cout << "MinimizerStandard::runMinimizer(): called" << endl;
    //    WaveIdeal *wave = new WaveIdeal(nParticles, m_config->nDimensions());
    //    HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(nParticles, m_config->nDimensions(), charge);


#ifdef USE_MPI
    timeStart = MPI_Wtime();
#endif

    string outfilename;
    if (m_config->rank() == 0) {
        outfilename = "output.dat";
        ofile.open(outfilename.c_str());
    }

    total_cumulative_e = new double[m_nVariations+1];
    total_cumulative_e2 = new double[m_nVariations+1];
    cumulative_e = new double[m_nVariations+1];
    cumulative_e2 = new double[m_nVariations+1];

    //  initialize the arrays  by zeroing them
    for( i=1; i <= m_nVariations; i++){
        cumulative_e[i] = cumulative_e2[i]  = total_cumulative_e[i] = total_cumulative_e2[i]  = 0.0;
    }
#ifdef USE_MPI
    // broadcast the total number of  variations
    MPI_Bcast (&m_nVariations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&m_nCycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    total_number_cycles = m_nCycles*m_config->nProcesses();

    // array to store all energies for last variation of alpha
    m_allEnergies = new double[m_nCycles+1];

    double *energies = new double[2];

    //  Do the mc sampling  and accumulate data with MPI_Reduce
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(m_wave, m_hamiltonian, m_config->nParticles(), m_config->nDimensions(), charge, m_config->rank(), stepLength);

    double beta = 0.4;
    double alpha = 0.5*charge;
    // loop over variational parameters
    for (int variate=1; variate <= m_nVariations; variate++){
        m_wave->setParameters(alpha, beta);
        monteCarlo->sample(m_nCycles, energies, m_allEnergies);
        // update the energy average and its squared
        cumulative_e[variate] = energies[0];
        cumulative_e2[variate] = energies[1];
        alpha += 0.1;
    }
#ifdef USE_MPI
    //  Collect data in total averages
    for( i=1; i <= m_nVariations; i++){
        MPI_Reduce(&cumulative_e[i], &total_cumulative_e[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&cumulative_e2[i], &total_cumulative_e2[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    timeEnd = MPI_Wtime();
    totalTime = timeEnd-timeStart;
#else
    for( i=1; i <= m_nVariations; i++){
        total_cumulative_e[i] = cumulative_e[i];
        total_cumulative_e2[i] = cumulative_e2[i];
    }
    totalTime = -1;
#endif
    // Print out results
    if ( m_config->rank() == 0) {
        cout << "Time = " <<  totalTime  << " on number of processors: "  << m_config->nProcesses()  << endl;
        alpha = 0.5*charge;
        for( i=1; i <= m_nVariations; i++){
            energy = total_cumulative_e[i]/total_number_cycles;
            variance = total_cumulative_e2[i]/total_number_cycles-energy*energy;
            error=sqrt(variance/(total_number_cycles-1));
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << setw(15) << setprecision(8) << alpha;
            ofile << setw(15) << setprecision(8) << energy;
            ofile << setw(15) << setprecision(8) << variance;
            ofile << setw(15) << setprecision(8) << error << endl;
            alpha += 0.1;
        }
        ofile.close();  // close output file
    }
    writeBlockData();
    delete [] total_cumulative_e; delete [] total_cumulative_e2;
    delete [] cumulative_e;
    delete [] cumulative_e2;
    delete [] m_allEnergies;
}
