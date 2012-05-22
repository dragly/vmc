#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

#include "densityplotter.h"

#include "wavefunction/wavefunction.h"
#include "inih/ini.h"
#include "config.h"
#include "montecarlo/standardmontecarlo.h"
#include "matrix.h"
#include "random.h"

using namespace std;

DensityPlotter::DensityPlotter(Config *config_) :
    config(config_),
    idum(config->idum())
{
    if(config->nDimensions() != 2) {
        cerr << "Density plots only implemented for two dimensions!" << endl;
        exit(928);
    }
    // allocate matrices which contain the position of the particles
    r_old = new vec2[config->nParticles()];
    r_new = new vec2[config->nParticles()];
}

DensityPlotter::~DensityPlotter()
{
    delete [] r_old;
    delete [] r_new;
}

void DensityPlotter::loadConfiguration(INIParser *settings)
{
    m_settings = settings;
    m_charge = atof(settings->Get("DensityPlotter", "charge", "1.0").c_str());
    m_stepLength = atof(settings->Get("DensityPlotter", "stepLength", "1.0").c_str());
    m_wave = WaveFunction::fromName(settings->Get("Wave", "class", "WaveSimple"), config);
    //    m_hamiltonian = Hamiltonian::fromName(settings->Get("Hamiltonian", "class", "HamiltonianSimple"), m_config, 1.0);
    m_nCycles = settings->GetInteger("DensityPlotter", "nCycles", 1000);
    double alpha = atof(settings->Get("DensityPlotter", "alpha", "1.0").c_str());
    double beta = atof(settings->Get("DensityPlotter", "beta", "1.0").c_str());
    double parameters[2];
    parameters[0] = alpha;
    parameters[1] = beta;

    m_wave->setParameters(parameters);
    aMax = atof(settings->Get("DensityPlotter", "aMax", "6.0").c_str());
    bMax = atof(settings->Get("DensityPlotter", "bMax", "6.0").c_str());
    aMin = atof(settings->Get("DensityPlotter", "aMin", "-6.0").c_str());
    bMin = atof(settings->Get("DensityPlotter", "bMin", "-6.0").c_str());
    aSteps = settings->GetInteger("DensityPlotter", "aSteps", 51);
    bSteps = settings->GetInteger("DensityPlotter", "bSteps", 51);
}

void DensityPlotter::divideSteps(int rank, int nProcesses, int totalSteps, StepConfig *stepConfig) {
    stepConfig->firstStep = (int)(rank * totalSteps / nProcesses);
    stepConfig->lastStep = (int)((rank + 1) * totalSteps / nProcesses) - 1;
    stepConfig->nSteps = stepConfig->lastStep - stepConfig->firstStep + 1;
}

void DensityPlotter::makePlot()
{
    MPI_Status status;
    double da = (aMax - aMin) / aSteps;
    double db = (bMax - bMin) / bSteps;
    double **probability = 0;
    double **myProbability = 0;
    if(config->rank() == 0) {
        probability = (double**)matrix(aSteps, bSteps, sizeof(double));
    }
    StepConfig myStepConfig;
    // split up the data for each process
    divideSteps(config->rank(), config->nProcesses(), aSteps, &myStepConfig);
    for(int i = 0; i < config->nProcesses(); i++) {
        if(config->rank() == i) {
            cout << myStepConfig.firstStep << " " << myStepConfig.lastStep << " " << myStepConfig.nSteps << endl;
            cout.flush();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // allocate my portion of the whole cake
    myProbability = (double**)matrix(myStepConfig.nSteps, bSteps, sizeof(double));
    // loop over positions
    for(int aStep = 0; aStep < myStepConfig.nSteps; aStep++) {
        for(int bStep = 0; bStep < bSteps; bStep++) {
            std::cout << "Calculating for a=" << aStep << " b=" << bStep << std::endl;
            r_new[0][0] = aMin + aStep * da + myStepConfig.firstStep * da;
            r_new[0][1] = bMin + bStep * db;
            double prob = 0;
            // loop over monte carlo cycles
            for (int cycle = 1; cycle <= m_nCycles; cycle++){
                // new positions for all particles
                for (int i = 1; i < config->nParticles(); i++) {
                    r_new[i][0] = aMin + (aMax - aMin) * ran2(idum);
                    r_new[i][1] = bMin + (bMax - bMin) * ran2(idum);
                }  //  end of loop over particles
                // compute probability
                prob += m_wave->evaluate(r_new) * m_wave->evaluate(r_new)/* * (r_new[0][0] * r_new[0][0]  + r_new[0][1] * r_new[0][1])*/;
            }   // end of loop over MC trials
            myProbability[aStep][bStep] = prob / m_nCycles;
        }
    }
    if(config->rank() == 0) {
        // collect all the data in the main matrix
        for(int i = 1; i < config->nProcesses(); i++) {
            StepConfig stepConfig;
            divideSteps(i, config->nProcesses(), aSteps, &stepConfig);
            MPI_Recv(probability[stepConfig.firstStep], stepConfig.nSteps * bSteps, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
        }
        memcpy(probability, myProbability, myStepConfig.nSteps * sizeof(double));
    } else {
        MPI_Send(myProbability[0], myStepConfig.nSteps * bSteps, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
    }

    if(config->rank() == 0) {
        plotFile.open("density.dat");
        params0File.open("density-grid-x.dat");
        params1File.open("density-grid-y.dat");
        for(int aStep = 0; aStep < aSteps; aStep++) {
            for(int bStep = 0; bStep < bSteps; bStep++) {
                plotFile << probability[aStep][bStep] << " ";
                params0File << aMin + aStep * da << " ";
                params1File << bMin + bStep * db << " ";
            }
            plotFile << "\n";
            params0File << "\n";
            params1File << "\n";
        }
        plotFile.close();
        params0File.close();
        params1File.close();
    }
}
