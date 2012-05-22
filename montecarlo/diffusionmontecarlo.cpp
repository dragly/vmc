#include "diffusionmontecarlo.h"
#include "standardmontecarlo.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/waveslater.h"
#include "../random.h"
#include "metropolishastingsmontecarlo.h"

DiffusionMonteCarlo::DiffusionMonteCarlo(Config *config) :
    MonteCarlo(config),
    tau(0.01)
{
    nWalkersMax = 10000;
    nWalkersIdeal = 500;
    nWalkersAlive = nWalkersIdeal;
    correlationStep = 500;

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
        if(i < nWalkersAlive) {
            walkers[i]->setAliveNew(true);
            walkers[i]->setAliveOld(true);
        }
    }
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
    std::cout << "Done constructing DMC class" << std::endl;
}

// Note: Implementation from http://www.ornl.gov/~pk7/thesis/pkthnode21.html

void DiffusionMonteCarlo::sample(int nCycles)
{
    // Initialize ensemble of walkers from VMC best guess
    MonteCarloMetropolisHastings *initialMonteCarlo = new MonteCarloMetropolisHastings(config);
    initialMonteCarlo->setRecordMoves(true, nWalkersAlive * nParticles);
    initialMonteCarlo->setThermalizationEnabled(true);
    initialMonteCarlo->sample(nWalkersAlive * correlationStep);
    double trialEnergy = initialMonteCarlo->energy();
    std::cout << "Initial trial energy was " << trialEnergy << std::endl;

    ofstream scatterfile;
    scatterfile.open("positions-init.dat");
    vec2 **moves = initialMonteCarlo->moves();
    for(int j = 0; j < nWalkersAlive; j++) {
        walkers[j]->initialize(moves[j]);
        if(j < nWalkersAlive) {
            for(int i = 0; i < nParticles; i++) {
                scatterfile << moves[j][i][0] << "\t" << moves[j][i][1] << std::endl;
            }
        }
    }
    scatterfile.close();

    int blockLength = 100;
    double energySum = 0;
    int nEnergySamples = 0;
    int blockSamples = 0;
    // For every cycle:
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // For every walker (configuration)
        nWalkersAlive = 0;
        for(int i = 0; i < nWalkersMax; i++) {
            DiffusionWalker *walker = walkers[i];
            if(walker->aliveOld()) {
                walker->advance(trialEnergy);
                energySum += walker->energy();
                nEnergySamples += walker->changeInEnergySamples();
                nWalkersAlive++;
            } // END if walker alive
        } // END for each walker
        //        std::cout << "Alive walkers: " << nWalkersAlive << "\xd" << std::endl;
        // Repeat configuration moves for about 100 - 1000 steps
        if(cycle < 2000 || !(cycle % blockLength)) {
            // Update trial energy ET to bring it closer to the current ensemble
            trialEnergy = energySum / nEnergySamples;
            std::cout << "Trial energy is now " << trialEnergy << " with " << nWalkersAlive << " walkers at cycle " << cycle << std::endl;
            // Renormalise the number of walkers to the target number by creating or deleting walkers
            //            while(nWalkersAlive > nWalkersIdeal) {
            //                int randomWalker = ran2(idum) * nWalkersMax;
            //                if(aliveNew[randomWalker]) {
            //                    aliveNew[randomWalker] = false;
            //                    nWalkersAlive--;
            //                }
            //            }

            energySum = 0;
            nEnergySamples = 0;

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
        }
        for(int i = 0; i < nWalkersMax; i++) {
            //            std::cout << "Walkers alive " << aliveOld[walker] << " " << aliveNew[walker] << std::endl;
            walkers[i]->progressToNextStep();
            //            std::cout << "Walkers alive after " << aliveOld[walker] << " " << aliveNew[walker] << std::endl;
        }
    }
    scatterfile.open("positions-end.dat");
    for(int j = 0; j < nWalkersMax; j++) {
        if(walkers[j]->aliveOld()) {
            for(int i = 0; i < nParticles; i++) {
                scatterfile << walkers[j]->positionsNew()[i][0] << "\t" << walkers[j]->positionsNew()[i][0] << std::endl;
            }
        }
    }
    scatterfile.close();

    m_energy = trialEnergy;
}
