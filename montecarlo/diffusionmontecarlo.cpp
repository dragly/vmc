#include "diffusionmontecarlo.h"
#include "standardmontecarlo.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/waveslater.h"
#include "../random.h"
#include "metropolishastingsmontecarlo.h"
#include "../inih/ini.h"

#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>

// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

using std::cout;
using std::setprecision;

DiffusionMonteCarlo::DiffusionMonteCarlo(Config *config_) :
    MonteCarlo(config_),
    correlationStep(500),
    nWalkersMax(10000),
    nWalkersIdeal(500),
    nWalkersAlive(nWalkersIdeal),
    nSamples(10000),
    nThermalizationCycles(2000),
    initialMonteCarlo(config->monteCarlo()),
    timeStep(0.0001)
{
    parameters[0] = 1;
    parameters[1] = 1;
    config->wave()->setParameters(parameters);
    //    rNew = new vec2*[nWalkersMax];
    //    for(int i = 0; i < nWalkersMax; i++) {
    //        rNew[i] = new vec2[nParticles];
    //    }
    //    rOld = new vec2*[nWalkersMax];
    //    for(int i = 0; i < nWalkersMax; i++) {
    //        rOld[i] = new vec2[nParticles];
    //    }
    walkers = new DiffusionWalker*[nWalkersMax];
    for(int i = 0; i < nWalkersMax; i++) {
        walkers[i] = new DiffusionWalker(config, walkers, nWalkersMax);
    }
    updateWalkerParameters();
    //    waves = new WaveFunction*[nWalkersMax];
    //    for(int i = 0; i < nWalkersMax; i++) {
    //        waves[i] = wave->clone();
    //    }
    //    aliveOld = new bool[nWalkersMax];
    //    aliveNew = new bool[nWalkersMax];
    //    for(int i = 0; i < nWalkersMax; i++) {
    //        if(i < nWalkersAlive) {
    //            aliveOld[i] = true;
    //            aliveNew[i] = true;
    //        } else {
    //            aliveOld[i] = false;
    //            aliveNew[i] = false;
    //        }
    //    }
}

void DiffusionMonteCarlo::updateWalkerParameters()
{
    for(int i = 0; i < nWalkersMax; i++) {
        if(i < nWalkersAlive) {
            walkers[i]->setAliveNew(true);
            walkers[i]->setAliveOld(true);
        } else {
            walkers[i]->setAliveNew(false);
            walkers[i]->setAliveOld(false);
        }
        walkers[i]->setParameters(parameters);
        walkers[i]->setTimeStep(timeStep);
    }
    config->wave()->setParameters(parameters);
}

void DiffusionMonteCarlo::loadConfiguration(INIParser *settings) {
    MonteCarlo::loadConfiguration(settings);
    double alpha = settings->GetDouble("DiffusionMonteCarlo", "alpha", 1.0);
    double beta = settings->GetDouble("DiffusionMonteCarlo", "beta", 1.0);
    nWalkersMax = settings->GetDouble("DiffusionMonteCarlo", "nWalkersMax", nWalkersMax);
    nWalkersIdeal = settings->GetDouble("DiffusionMonteCarlo", "nWalkersIdeal", nWalkersIdeal);
    nWalkersAlive = nWalkersIdeal;
    correlationStep = settings->GetDouble("DiffusionMonteCarlo", "correlationStep", correlationStep);
    nSamples = settings->GetDouble("DiffusionMonteCarlo", "nSamples", nSamples);
    timeStep = settings->GetDouble("DiffusionMonteCarlo", "tau", timeStep);
    nThermalizationCycles = settings->GetDouble("DiffusionMonteCarlo", "nThermalizationCycles", nThermalizationCycles);
    parameters[0] = alpha;
    parameters[1] = beta;
    updateWalkerParameters();
    wave->setParameters(parameters);
}

void DiffusionMonteCarlo::sample() {
    sample(nSamples);
}


/*!
  * \note Implementation based on http://www.ornl.gov/~pk7/thesis/pkthnode21.html
  */
void DiffusionMonteCarlo::sample(int nSamplesLocal)
{
    if(config->myRank() == 0) {
        struct stat st;
        if(stat(scratchDir.c_str(), &st) != 0) {
            if(mkdir(scratchDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
                std::cerr << "Error creating directory " << scratchDir << std::endl;
                exit(948);
            }
        }
    }
    std::cout << "Running VMC to initialize DMC" << std::endl;
    // Initialize ensemble of walkers from VMC best guess
    initialMonteCarlo->setOutputEnergies(true);
    initialMonteCarlo->setRecordMoves(true, nWalkersAlive * nParticles, "/scratch/positions/vmc-positions.dat");
    initialMonteCarlo->setThermalizationEnabled(true);
    //    initialMonteCarlo->setThermalizationEnabled(false); // TODO set true
    initialMonteCarlo->sample(nWalkersAlive * correlationStep);
    double trialEnergy = initialMonteCarlo->energy();
    std::cout << "Done initializing. Initial trial energy was " << trialEnergy << std::endl;

    vec2 **moves = initialMonteCarlo->moves();
    for(int j = 0; j < nWalkersAlive; j++) {
        walkers[j]->initialize(moves[j]);
    }
    std::cout << "Done writing positions to file." << std::endl;

    ofstream energyFile;
    energyFile.open("dmc-energies.dat");

    setRecordMoves(true, 0, "/scratch/positions/dmc-positions.dat");

    int blockLength = 100;
    double energySum = 0;
    int nEnergySamples = 0;
    int blockSamples = 0;
    double totalEnergySum = 0;
    double nTotalEnergySamples = 0;
    int acceptances = 0;
    int rejections = 0;
    // For every cycle:
    for(int cycle = 0; cycle < nSamplesLocal; cycle++) {
        // For every walker (configuration)
        nWalkersAlive = 0;
        for(int i = 0; i < nWalkersMax; i++) {
            DiffusionWalker *walker = walkers[i];
            if(walker->aliveOld()) {
                walker->advance(trialEnergy);

                energySum += walker->energy();
                nEnergySamples += walker->changeInEnergySamples();
                totalEnergySum += walker->energy();
                nTotalEnergySamples += walker->changeInEnergySamples();
                acceptances += walker->acceptances();
                rejections += walker->rejections();
                nWalkersAlive++;
            } // END if walker alive
        } // END for each walker
        // Update state of every walker
        for(int i = 0; i < nWalkersMax; i++) {
            walkers[i]->setAliveOld(walkers[i]->aliveNew());
        }
        // Calculate the mean energy
        double meanEnergy = energySum / double(nEnergySamples);
        energyFile << cycle << " " << meanEnergy << " " << nWalkersAlive << std::endl;
        //        std::cout << "Alive walkers: " << nWalkersAlive << "\xd" << std::endl;
        // Repeat configuration moves for about 100 - 1000 steps
        if(cycle < nThermalizationCycles || !(cycle % blockLength)) {
            // Update trial energy ET to bring it closer to the current ensemble
            if(cycle < nThermalizationCycles) {
                trialEnergy = meanEnergy - log((double)nWalkersAlive / (double) nWalkersIdeal);
            } else {
                trialEnergy = meanEnergy;
            }
            std::cout << "DMC cycle " << cycle << ": Trial energy " << setw(14) << setprecision(10) << trialEnergy << ", average " << setw(14) << setprecision(10) << totalEnergySum / nTotalEnergySamples << ",  acceptance: " << setw(8) <<  setprecision(6) << acceptances / double(acceptances + rejections) << ", with " << nWalkersAlive << " walkers at cycle " << cycle << std::endl;
            // Renormalise the number of walkers to the target number by creating or deleting walkers
            //                        while(nWalkersAlive > nWalkersIdeal) {
            //                            int randomWalker = ran3(idumMC) * nWalkersMax;
            //                            if(walkers[randomWalker]->aliveNew()) {
            //                                walkers[randomWalker]->setAliveNew(false);
            //                                nWalkersAlive--;
            //                            }
            //                        }
            blockSamples++;
            energySum = 0;
            nEnergySamples = 0;
            acceptances = 0;
            rejections = 0;
            if(cycle == nThermalizationCycles) {
                totalEnergySum = 0;
                nTotalEnergySamples = 0;
            }
        }
        int testNum = nSamplesLocal * nParticles * nWalkersAlive / 1e9;
        int modulusA = 0;
        if(testNum > 0) {
            modulusA = cycle % testNum;
        }
        if(cycle > nThermalizationCycles && !modulusA) { // store 10^x positions
            for(int j = 0; j < nWalkersMax; j++) {
                if(walkers[j]->aliveOld()) {
                    for(int i = 0; i < nParticles; i++) {
                        if(walkers[j]->positionsNew()[i][0] == 0.0) {
                            std::cerr << "WTF?" << std::endl;
                        }
                        writePositionToFile(walkers[j]->positionsNew()[i]);
                    }
                }
            }
        }


    }
    energyFile.close();

    m_energy = trialEnergy;
}
