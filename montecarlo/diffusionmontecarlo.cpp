#include "diffusionmontecarlo.h"
#include "standardmontecarlo.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/waveslater.h"
#include "../random.h"
#include "metropolishastingsmontecarlo.h"
#include "../inih/ini.h"

#include <iostream>
#include <iomanip>

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
    parameters[0] = 0;
    parameters[1] = 0;
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
}

void DiffusionMonteCarlo::sample() {
    sample(nSamples);
}


/*!
  * \note Implementation based on http://www.ornl.gov/~pk7/thesis/pkthnode21.html
  */
void DiffusionMonteCarlo::sample(int nSamplesLocal)
{
    std::cout << "Running VMC to initialize DMC" << std::endl;
    // Initialize ensemble of walkers from VMC best guess
    initialMonteCarlo->setRecordMoves(true, nWalkersAlive * nParticles);
    initialMonteCarlo->setThermalizationEnabled(true);
//    initialMonteCarlo->setThermalizationEnabled(false); // TODO set true
    initialMonteCarlo->sample(nWalkersAlive * correlationStep);
    double trialEnergy = initialMonteCarlo->energy();
    std::cout << "Done initializing. Initial trial energy was " << trialEnergy << std::endl;

    ofstream positionFile;
    positionFile.open("dmc-positions-init.dat");
    vec2 **moves = initialMonteCarlo->moves();
    for(int j = 0; j < nWalkersAlive; j++) {
        walkers[j]->initialize(moves[j]);
        for(int i = 0; i < nParticles; i++) {
            positionFile << moves[j][i][0] << " " << moves[j][i][1] << std::endl;
        }
    }
    positionFile.close();
    std::cout << "Done writing positions to file." << std::endl;

    ofstream energyFile;
    energyFile.open("dmc-energies.dat");

    int blockLength = 100;
    double energySum = 0;
    int nEnergySamples = 0;
    int blockSamples = 0;
    double totalEnergySum = 0;
    double nTotalEnergySamples = 0;
    positionFile.open("dmc-positions-end.dat");
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
        double meanEnergy = energySum / nEnergySamples;
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
            std::cout << "Trial energy is now " << setprecision(20) << trialEnergy << ", average " << totalEnergySum / nTotalEnergySamples << ",  acceptance: " << acceptances / double(acceptances + rejections) << ", with " << nWalkersAlive << " walkers at cycle " << cycle << std::endl;
            // Renormalise the number of walkers to the target number by creating or deleting walkers
//                        while(nWalkersAlive > nWalkersIdeal) {
//                            int randomWalker = ran3(idumMC) * nWalkersMax;
//                            if(walkers[randomWalker]->aliveNew()) {
//                                walkers[randomWalker]->setAliveNew(false);
//                                nWalkersAlive--;
//                            }
//                        }


//            stringstream fileName;
//            fileName << "positions-" << blockSamples << ".dat";
//            scatterfile.open(fileName.str().c_str());
//            for(int j = 0; j < nWalkersMax; j++) {
//                if(walkers[j]->aliveOld()) {
//                    for(int i = 0; i < nParticles; i++) {
//                        scatterfile << walkers[j]->positionsNew()[i][0] << "\t" << walkers[j]->positionsNew()[i][1] << std::endl;
//                    }
//                }
//            }
//            scatterfile.close();
            blockSamples++;
            if(cycle > nThermalizationCycles) {
                for(int j = 0; j < nWalkersMax; j++) {
                    if(walkers[j]->aliveOld()) {
                        for(int i = 0; i < nParticles; i++) {
                            positionFile << walkers[j]->positionsNew()[i][0] << " " << walkers[j]->positionsNew()[i][1] << std::endl;
                        }
                    }
                }
            }
            energySum = 0;
            nEnergySamples = 0;
            acceptances = 0;
            rejections = 0;
            if(cycle == nThermalizationCycles) {
                totalEnergySum = 0;
                nTotalEnergySamples = 0;
            }
        }
        for(int i = 0; i < nWalkersMax; i++) {
            //            std::cout << "Walkers alive " << aliveOld[walker] << " " << aliveNew[walker] << std::endl;
            walkers[i]->progressToNextStep();
            //            std::cout << "Walkers alive after " << aliveOld[walker] << " " << aliveNew[walker] << std::endl;
        }
    }
    energyFile.close();

    positionFile.close();

    m_energy = trialEnergy;
}
