#include "diffusionmontecarlo.h"

DiffusionMonteCarlo::DiffusionMonteCarlo(Config *config)
    : MonteCarlo(config)
{
}

// Note: Implementation from http://www.ornl.gov/~pk7/thesis/pkthnode21.html

void DiffusionMonteCarlo::sample(int numberCycles)
{
    int nStartingPositions = 100;
    double sumEvaluation = 0;
    int nTotalWalkers = 1000;

    mat evaluations = zeros<mat>(nStartingPositions, nStartingPositions);
    // Initialize ensemble of particles from VMC best guess
    for(int i = 0; i < nStartingPositions; i++) {
        for(int i = 0; i < nStartingPositions; j++) {
            vec2 r;
            r[0] = -2 + 4 * i / (double)nStartingPositions;
            r[1] = -2 + 4 * j / (double)nStartingPositions;
            double evaluation = wave->evaluate(r);
            evaluations(i,j) = evaluation;
            sumEvaluation += evaluation;
        }
    }

    for(int i = 0; i < nStartingPositions; i++) {
        for(int i = 0; i < nStartingPositions; j++) {
            vec2 r;
            r[0] = -2 + 4 * i / (double)nStartingPositions;
            r[1] = -2 + 4 * j / (double)nStartingPositions;
            int nWalkersHere = evaluations(i,j) / sumEvaluation;
            for(int k = 0; k < nWalkersHere; k++) {
                // TODO: Spawn walker here
            }
        }
    }

    // TODO: Terminalize all walkers

    // For every configuration:
           // For every electron
                // Propose move (with quantum force)
                // Apply fixed node approximation (keep sign or reject move)
                // Compute weight function
                // TODO: Read up on and implement Green's function
                // Accept move according to Metropolis probability
            // Compute branching factor PB
            // Accumulate the energy and any observables weighted by PB
    // Repeat configuration moves for about 100 - 1000 steps
    // Update trial energy ET to bring it closer to the current ensemble
    // Renormalise the number of walkers to the target number by creating or deleting walkers
    // Repeat and rinse

}
