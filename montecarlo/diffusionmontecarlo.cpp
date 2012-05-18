#include "diffusionmontecarlo.h"
#include "montecarlostandard.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/waveslater.h"
#include "../random.h"

DiffusionMonteCarlo::DiffusionMonteCarlo(Config *config) :
    MonteCarlo(config),
    tau(1)
{
    nTotalWalkers = 1000;
    correlationStep = 200;

    rNew = new vec2*[nTotalWalkers];
    for(int i = 0; i < nTotalWalkers; i++) {
        rNew[i] = new vec2[nParticles];
    }
    rOld = new vec2*[nTotalWalkers];
    for(int i = 0; i < nTotalWalkers; i++) {
        rOld[i] = new vec2[nParticles];
    }
    waves = new WaveFunction*[nTotalWalkers];
    for(int i = 0; i < nTotalWalkers; i++) {
        waves[i] = wave->clone();
    }
}

// Note: Implementation from http://www.ornl.gov/~pk7/thesis/pkthnode21.html

void DiffusionMonteCarlo::sample(int nCycles)
{
    // Initialize ensemble of walkers from VMC best guess
    MonteCarloStandard *monteCarloStandard = new MonteCarloStandard(config);
    monteCarloStandard->setRecordMoves(true, nTotalWalkers * nParticles);
    monteCarloStandard->sample(nTotalWalkers * correlationStep);

    ofstream scatterfile;
    scatterfile.open("positions-init.dat");
    for(int j = 0; j < nTotalWalkers; j++) {
        for(int i = 0; i < nParticles; i++) {
            vec2 **moves = monteCarloStandard->moves();
            rOld[j][i] = moves[j][i];
            rNew[j][i] = moves[j][i];
            scatterfile << moves[j][i][0] << "\t" << moves[j][i][1] << std::endl;
        }
    }
    scatterfile.close();

    // TODO: Terminalize all walkers

    vec *quantumForceNew = new vec[nTotalWalkers];
    vec *quantumForceOld = new vec[nTotalWalkers];

    for(int j = 0; j < nTotalWalkers; j++) {
        quantumForceNew[j] = zeros<vec>(nDimensions * nParticles);
        quantumForceOld[j] = zeros<vec>(nDimensions * nParticles);
        waves[j]->initialize(rNew[j]);
        waves[j]->gradient(rNew[j], 0, quantumForceNew[j]);
        quantumForceOld[j] = quantumForceNew[j];
    }
    // For every configuration:
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // For every walker (configuration)
        for(int j = 0; j < nTotalWalkers; j++) {
            // For every electron
            for(int i = 0; i < nParticles; i++) {
                waves[j]->gradient(rNew[j], i, quantumForceNew[j]);

                // Propose move (with quantum force)
                for(int k = 0; k < nDimensions; k++) {
                    // TODO per cartesian component tau?
                    rNew[j][i][k] = rOld[j][i][k] + tau * quantumForceNew[j][i * nDimensions + k] + simpleGaussRandom(idum);
//                    std::cout << "i ndim k " << i * nDimensions + k << " " << quantumForceNew->n_elem << std::endl;
//                    std::cout << "Quantum force " << j << " " << i << " " << k << " " << quantumForceNew[j][i * nDimensions + k] << std::endl;
//                    std::cout << "Old/New: " << rOld[j][i][k] << " " << rNew[j][i][k] << std::endl;
                }
                double ratio = waves[j]->ratio(rNew[j][i], i);
                // Apply fixed node approximation (keep sign or reject move)
                if(ratio > 0) {
                    // Compute weight function
                    double argSum = 0;
                    for(int k = 0; k < nDimensions; k++) {
                        double argument = rNew[j][i][k] - rOld[j][i][k] - tau * quantumForceOld[j][i * nDimensions + k];
                        argSum += argument * argument;
                    }
                    argSum *= 1 / 4.;
                    double greensRatio = pow((2 * M_PI * tau),(-3*nParticles / 2)) * exp(-(argSum) / 2 * tau);
                    double weight = ratio*ratio * greensRatio;
                    // Accept move according to Metropolis probability
                    if(weight > ran2(idum) - 0.5) {
                        waves[j]->acceptMove(i);
                        rOld[j][i] = rNew[j][i];
                    } else {
                        waves[j]->rejectMove();
                        rNew[j][i] = rOld[j][i];
                    }
                } else {
                    waves[j]->rejectMove();
                }
            }
            // Compute branching factor PB
            // Accumulate the energy and any observables weighted by PB
        }
    }
    scatterfile.open("positions-end.dat");
    for(int j = 0; j < nTotalWalkers; j++) {
        for(int i = 0; i < nParticles; i++) {
            scatterfile << rNew[j][i][0] << "\t" << rNew[j][i][1] << std::endl;
        }
    }
    scatterfile.close();
    // Repeat configuration moves for about 100 - 1000 steps
    // Update trial energy ET to bring it closer to the current ensemble
    // Renormalise the number of walkers to the target number by creating or deleting walkers
    // Repeat and rinse

}
