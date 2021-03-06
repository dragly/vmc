#include "standardminimizer.h"

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
#include <typeinfo>

// local stuff
#include "../inih/ini.h"

#include "../wavefunction/wavesimple.h"
#include "../montecarlo/standardmontecarlo.h"
#include "../matrix.h"
#include "../hamiltonian/hamiltoniansimple.h"
#include "../hamiltonian/hamiltonianideal.h"
#include "../wavefunction/waveideal.h"
#include "../wavefunction/waveslater.h"

using namespace std;

StandardMinimizer::StandardMinimizer(Config *config) :
    Minimizer(config),
    m_wave(config->wave()),
    m_monteCarlo(config->monteCarlo()),
    m_hamiltonian(config->hamiltonian())
{

}

void StandardMinimizer::loadConfiguration(INIParser *settings)
{
            m_settings = settings;
            m_nSamples = settings->GetDouble("MinimizerStandard","nCycles", 1000);
            nVariations = settings->GetDouble("MinimizerStandard","nVariations", 11);    //  default number of variations
            alphaStart = settings->GetDouble("MinimizerStandard","alphaStart", 0);    //  default number of variations
            alphaEnd= settings->GetDouble("MinimizerStandard","alphaEnd", 1);    //  default number of variations
            betaStart = settings->GetDouble("MinimizerStandard","betaStart", 0);    //  default number of variations
            betaEnd = settings->GetDouble("MinimizerStandard","betaEnd", 1);    //  default number of variations
}

void StandardMinimizer::runMinimizer()
{
    int total_number_cycles;
    mat cumulativeEnergy, cumulativeEnergySquared;
    mat totalCumulativeEnergy, totalCumulativeEnergySquared;

    mat parameter0Map;
    mat parameter1Map;
    double timeStart;
    double timeEnd;
    double totalTime;

    //    WaveIdeal *wave = new WaveIdeal(nParticles, m_config->nDimensions());
    //    HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(nParticles, m_config->nDimensions(), charge);


    timeStart = MPI_Wtime();

    // broadcast the total number of  variations
    MPI_Bcast (&nVariations, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&m_nSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //  initialize the arrays  by zeroing them
    totalCumulativeEnergy = zeros<mat>(nVariations,nVariations);
    totalCumulativeEnergySquared = zeros<mat>(nVariations,nVariations);
    cumulativeEnergy = zeros<mat>(nVariations,nVariations);
    cumulativeEnergySquared = zeros<mat>(nVariations,nVariations);

    parameter0Map = zeros<mat>(nVariations, nVariations);
    parameter1Map = zeros<mat>(nVariations, nVariations);

    total_number_cycles = m_nSamples*config->m_nProcesses();

    // array to store all energies for last variation of alpha


    //    double *energies = new double[2];
    double parameters[2];

    // loop over variational parameters
    double alphaStep = (alphaEnd - alphaStart) / (nVariations - 1);
    double betaStep = (betaEnd - betaStart) / (nVariations - 1);
    parameters[0] = alphaStart;
    for (int i=0; i < nVariations; i++){
        parameters[1] = betaStart;
        for (int j=0; j < nVariations; j++){
            std::cout << "Testing parameters " << parameters[0] << " " << parameters[1] << std::endl;
            std::cout << "with " << m_nSamples << " cycles" << std::endl;
            m_wave->setParameters(parameters);
            m_monteCarlo->sample(m_nSamples);
            // update the energy average and its squared
            cumulativeEnergy(i,j) = m_monteCarlo->energy();
            cumulativeEnergySquared(i,j) = m_monteCarlo->energySquared();
            std::cout << "Got energy of " << cumulativeEnergy(i,j) << std::endl;
            parameter0Map(i,j) = parameters[0];
            parameter1Map(i,j) = parameters[1];
            parameters[1] += betaStep;
        }
        parameters[0] += alphaStep;
    }
    //  Collect data in total averages
    for(int i=0; i < nVariations; i++){
        for(int j=0; j < nVariations; j++){
            double cumulativeEnergyDummy = cumulativeEnergy(i,j);
            double totalCumulativeEnergyDummy = totalCumulativeEnergy(i,j);
            MPI_Reduce(&cumulativeEnergyDummy, &totalCumulativeEnergyDummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if(m_rank == 0) {
                std::cout << totalCumulativeEnergyDummy << std::endl;
            }
            totalCumulativeEnergy(i,j) = totalCumulativeEnergyDummy / m_m_nProcesses;

            double cumulativeEnergySquaredDummy = cumulativeEnergySquared(i,j);
            double totalCumulativeEnergySquaredDummy = totalCumulativeEnergySquared(i,j);
            MPI_Reduce(&cumulativeEnergySquaredDummy, &totalCumulativeEnergySquaredDummy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            totalCumulativeEnergySquared(i,j) = totalCumulativeEnergySquaredDummy / m_m_nProcesses;
        }
    }
    timeEnd = MPI_Wtime();
    totalTime = timeEnd-timeStart;

    // Print out results
    if ( m_rank == 0) {

        // output file as global variable
        ofstream energyFile;
        ofstream parameters0File;
        ofstream parameters1File;
        ofstream errorFile;
        ofstream varianceFile;
        string outfilename;
        outfilename = "energies.dat";
        energyFile.open(outfilename.c_str());
        outfilename = "parameters0.dat";
        parameters0File.open(outfilename.c_str());
        outfilename = "parameters1.dat";
        parameters1File.open(outfilename.c_str());
        outfilename = "variances.dat";
        varianceFile.open(outfilename.c_str());
        outfilename = "errors.dat";
        errorFile.open(outfilename.c_str());
        cout << "Time = " <<  totalTime  << " on number of processors: "  << config->m_nProcesses()  << endl;
        for(int i=0; i < nVariations; i++){
            for(int j=0; j < nVariations; j++){
                double energy = totalCumulativeEnergy(i,j) ;
                double variance = totalCumulativeEnergySquared(i,j) - energy*energy;
                double error = sqrt(variance / (total_number_cycles-1));
                energyFile << setw(15) << setprecision(20) << energy << " ";
                parameters0File << setw(15) << setprecision(20) << parameter0Map(i,j) << " ";
                parameters1File << setw(15) << setprecision(20) << parameter1Map(i,j) << " ";
                varianceFile << setw(15) << setprecision(20) << variance << " ";
                errorFile << setw(15) << setprecision(20) << error << " ";
            }
            energyFile << std::endl;
            varianceFile << std::endl;
            errorFile << std::endl;
            parameters0File << std::endl;
            parameters1File << std::endl;
        }
        energyFile.close();  // close output file
    }
}
