#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#ifdef USE_MPI
// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif

#include "densityplotter.h"

#include "wavefunction/wavefunction.h"
#include "inih/cpp/INIReader.h"
#include "config.h"
#include "montecarlo/montecarlostandard.h"
#include "matrix.h"
#include "random.h"

using namespace std;

DensityPlotter::DensityPlotter(Config *config) :
    m_config(config)
{
    if(config->nDimensions() != 2) {
        cerr << "Density plots only implemented for two dimensions!" << endl;
        exit(928);
    }
    // every node has its own seed for the random numbers
    idum = -1;
    // allocate matrices which contain the position of the particles
    r_old = new vec2[config->nParticles()];
    r_new = new vec2[config->nParticles()];
}

DensityPlotter::~DensityPlotter()
{
    delete [] r_old;
    delete [] r_new;
}

void DensityPlotter::loadConfiguration(INIReader *settings)
{
    m_settings = settings;
    m_charge = atof(settings->Get("DensityPlotter", "charge", "1.0").c_str());
    m_stepLength = atof(settings->Get("DensityPlotter", "stepLength", "1.0").c_str());
    m_wave = WaveFunction::fromName(settings->Get("Wave", "class", "WaveSimple"), m_config);
    //    m_hamiltonian = Hamiltonian::fromName(settings->Get("Hamiltonian", "class", "HamiltonianSimple"), m_config, 1.0);
    m_nCycles = settings->GetInteger("DensityPlotter", "nCycles", 1000);
    double alpha = atof(settings->Get("DensityPlotter", "alpha", "1.0").c_str());
    double beta = atof(settings->Get("DensityPlotter", "beta", "1.0").c_str());
    m_wave->setParameters(alpha, beta);
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
#ifdef USE_MPI
    MPI_Status status;
#endif
    double da = (aMax - aMin) / aSteps;
    double db = (bMax - bMin) / bSteps;
    double **probability = 0;
    double **myProbability = 0;
    if(m_config->rank() == 0) {
        plotFile.open("plot.dat");
        probability = (double**)matrix(aSteps, bSteps, sizeof(double));
    }
    StepConfig myStepConfig;
    // split up the data for each process
    divideSteps(m_config->rank(), m_config->nProcesses(), aSteps, &myStepConfig);
    for(int i = 0; i < m_config->nProcesses(); i++) {
        if(m_config->rank() == i) {
            cout << myStepConfig.firstStep << " " << myStepConfig.lastStep << " " << myStepConfig.nSteps << endl;
            cout.flush();
        }
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    // allocate my portion of the whole cake
    myProbability = (double**)matrix(myStepConfig.nSteps, bSteps, sizeof(double));
    // loop over positions
    for(int aStep = 0; aStep < myStepConfig.nSteps; aStep++) {
        for(int bStep = 0; bStep < bSteps; bStep++) {
            r_new[0][0] = aMin + aStep * da + myStepConfig.firstStep * da;
            r_new[0][1] = bMin + bStep * db;
            double prob = 0;
            // loop over monte carlo cycles
            for (int cycle = 1; cycle <= m_nCycles; cycle++){
                // new positions for all particles
                for (int i = 1; i < m_config->nParticles(); i++) {
                    r_new[i][0] = aMin + (aMax - aMin) * ran2(&idum);
                    r_new[i][1] = bMin + (bMax - bMin) * ran2(&idum);
                }  //  end of loop over particles
                // compute probability
                prob += m_wave->wave(r_new) * m_wave->wave(r_new) * (r_new[0][0] * r_new[0][0]  + r_new[0][1] * r_new[0][1]);
            }   // end of loop over MC trials
            myProbability[aStep][bStep] = prob / m_nCycles;
        }
    }
    if(m_config->rank() == 0) {
#ifdef USE_MPI
        // collect all the data in the main matrix
        for(int i = 1; i < m_config->nProcesses(); i++) {
            StepConfig stepConfig;
            divideSteps(i, m_config->nProcesses(), aSteps, &stepConfig);
            MPI_Recv(probability[stepConfig.firstStep], stepConfig.nSteps * bSteps, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
        }
#endif
        memcpy(probability, myProbability, myStepConfig.nSteps * sizeof(double));
    } else {
#ifdef USE_MPI
        MPI_Send(myProbability[0], myStepConfig.nSteps * bSteps, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
#endif
    }

    if(m_config->rank() == 0) {
        for(int aStep = 0; aStep < aSteps; aStep++) {
            for(int bStep = 0; bStep < bSteps; bStep++) {
                plotFile << probability[aStep][bStep] << "\t";
                //            printf("%.4f    ", probability[aStep][bStep]);
            }
            plotFile << "\n";
            //        printf("\n");
        }
        plotFile.close();
    }
}
