//#include <QSettings>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"
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

MinimizerStandard::MinimizerStandard(int rank, int nProcesses) :
    Minimizer(rank, nProcesses)
{
}

void MinimizerStandard::loadConfiguration(INIReader *settings)
{
    m_settings = settings;
    dimension = atoi(settings->Get("MinimizerStandard","dimension", "2").c_str());
    charge = atof(settings->Get("MinimizerStandard","charge", "1.0").c_str());
    stepLength = atof(settings->Get("MinimizerStandard","stepLength", "1.0").c_str());
    nParticles = atoi(settings->Get("MinimizerStandard","nParticles", "2").c_str());
    nCycles = atoi(settings->Get("MinimizerStandard","nCycles", "1000").c_str());
    maxVariations = atoi(settings->Get("MinimizerStandard","maxVariations", "11").c_str());    //  default number of variations
    waveClass = settings->Get("Wave","class", "WaveSimple");
    waveUseAnalyticalLaplace = atoi(settings->Get("Wave","useAnalyticalLaplace", "0").c_str());
    hamiltonianClass = settings->Get("Hamiltonian","class", "HamiltonianSimple");
}

void MinimizerStandard::runMinimizer()
{
    int total_number_cycles, i;
    double *cumulative_e, *cumulative_e2;
    double *total_cumulative_e, *total_cumulative_e2;
    double  timeStart;
    double timeEnd;
    double totalTime;
    double variance;
    double energy;
    double error;

    cout << "MinimizerStandard::runMinimizer(): called" << endl;
    //    WaveIdeal *wave = new WaveIdeal(nParticles, dimension);
    //    HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(nParticles, dimension, charge);
    WaveFunction *wave;
    if(waveClass == "WaveSimple") {
        WaveSimple *waveSimple = new WaveSimple(nParticles, dimension);
        wave = waveSimple;
    } else if(waveClass == "WaveIdeal") {
        WaveIdeal *waveIdeal = new WaveIdeal(nParticles, dimension);
        wave = waveIdeal;
    } else {
        cerr << "Unknown wave class!" << endl;
        exit(99);
    }
    wave->loadConfiguration(m_settings);

    Hamiltonian *hamiltonian;
    if(hamiltonianClass == "HamiltonianSimple") {
        HamiltonianSimple *hamiltonianSimple = new HamiltonianSimple(nParticles, dimension, charge);
        hamiltonian = hamiltonianSimple;
    } else if(hamiltonianClass == "HamiltonianIdeal") {
        HamiltonianIdeal *hamiltonianIdeal = new HamiltonianIdeal(nParticles, dimension, charge);
        hamiltonian = hamiltonianIdeal;
    } else {
        cerr << "Unknown hamiltonian class!" << endl;
        exit(98);
    }
    timeStart = MPI_Wtime();

    string outfilename;
    if (m_rank == 0) {
        outfilename = "output.dat";
        ofile.open(outfilename.c_str());
    }

    total_cumulative_e = new double[maxVariations+1];
    total_cumulative_e2 = new double[maxVariations+1];
    cumulative_e = new double[maxVariations+1];
    cumulative_e2 = new double[maxVariations+1];

    //  initialize the arrays  by zeroing them
    for( i=1; i <= maxVariations; i++){
        cumulative_e[i] = cumulative_e2[i]  = total_cumulative_e[i] = total_cumulative_e2[i]  = 0.0;
    }

    // broadcast the total number of  variations
    MPI_Bcast (&maxVariations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nCycles, 1, MPI_INT, 0, MPI_COMM_WORLD);

    total_number_cycles = nCycles*m_nProcesses;

    // array to store all energies for last variation of alpha
    allEnergies = new double[nCycles+1];

    double *energies = new double[2];

    //  Do the mc sampling  and accumulate data with MPI_Reduce
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(wave, hamiltonian, nParticles, dimension, charge, m_rank, stepLength);

    double beta = 0.4;
    double alpha = 0.5*charge;
    // loop over variational parameters
    for (int variate=1; variate <= maxVariations; variate++){
        wave->setParameters(alpha, beta);
        monteCarlo->sample(nCycles, energies, allEnergies);
        // update the energy average and its squared
        cumulative_e[variate] = energies[0];
        cumulative_e2[variate] = energies[1];
        alpha += 0.1;
    }
    //  Collect data in total averages
    for( i=1; i <= maxVariations; i++){
        MPI_Reduce(&cumulative_e[i], &total_cumulative_e[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&cumulative_e2[i], &total_cumulative_e2[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    timeEnd = MPI_Wtime();
    totalTime = timeEnd-timeStart;
    // Print out results
    if ( m_rank == 0) {
        cout << "Time = " <<  totalTime  << " on number of processors: "  << m_nProcesses  << endl;
        alpha = 0.5*charge;
        for( i=1; i <= maxVariations; i++){
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
    delete [] allEnergies;
}
