#include "montecarlometropolishastings.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"
#include "../config.h"

MonteCarloMetropolisHastings::MonteCarloMetropolisHastings(Config *config) :
    MonteCarlo(config),
    rank(config->rank()),
    stepLength(config->stepLength()),
    hamiltonian(config->hamiltonian())
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[ nParticles];
    rNew = new vec2[ nParticles];

    forceVectorNew = zeros<vec>(nParticles * nDimensions);
    forceVectorOld = zeros<vec>(nParticles * nDimensions);
}

MonteCarloMetropolisHastings::~MonteCarloMetropolisHastings()
{
}

void MonteCarloMetropolisHastings::quantumForce(vec2 rPosition[], vec &forceVector) {
    wave->gradient(rPosition, 0, forceVector); // TODO add particle number
    forceVector *= 2;
}

void MonteCarloMetropolisHastings::sample(int nCycles)
{
    m_energy = 0;
    m_energySquared = 0;
    terminalizationSum = 0;
    terminalizationNum = 1;
    double localEnergy = 0;
    double diffConstant = 1;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = stepLength*(ran2(idum)-0.5);
        }
    }
    wave->initialize(rOld);
//    wave->gradient(rOld, 0, waveGradientOld); // TODO add particle number
    quantumForce(rOld, forceVectorOld);
    // loop over monte carlo cycles
    for (int cycle = 1; cycle <= nCycles; cycle++){
        // new trial position
        for (int i = 0; i < nParticles; i++) {
            quantumForce(rOld, forceVectorNew);
            vec2 particleQuantumForce;
            particleQuantumForce(0) = forceVectorNew(i);
            particleQuantumForce(1) = forceVectorNew(i+1);
            rNew[i] = rOld[i] + diffConstant*particleQuantumForce*stepLength;
            for (int j=0; j < nDimensions; j++) {
                rNew[i][j] += simpleGaussRandom(idum) * sqrt(stepLength);
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < nParticles; k++) {
                if ( k != i) {
                    rNew[k] = rOld[k];
                }
            }
//            wfnew = wave->evaluate(rNew);
//            wave->gradient(rNew, i, waveGradientNew);
//            double argument = 0;
//            for( int j = 0; j < nDimensions; j++) {
//                forceVectorSum[j] = forceVectorNew[j] + forceVectorOld[j];
//                forceVectorDiff[j] = forceVectorOld[j] - forceVectorNew[j];
//                positionDiff[j] = rNew[i][j] - rOld[i][j];
//                argument += 0.5 * forceVectorSum[j] * (diffConstant * stepLength / 2 * forceVectorDiff[j] - positionDiff[j]);
//            }
//            double waveFrac = wfnew*wfnew/(wfold*wfold);
            // The Metropolis test is performed by moving one particle at the time
            double ratio = wave->ratio(rNew[i], i);
            if(ran2(idum) <= (ratio*ratio)) {
                rOld[i] = rNew[i];
                wave->acceptEvaluation(i);
//                std::cout << "Accepted" << std::endl;
            } else {
                rNew[i] = rOld[i]; // Move the particle back
                wave->refuseEvaluation();
            }
            localEnergy = hamiltonian->energy(wave, rOld);
            if(terminalized) {
                if(storeEnergies) {
                    m_allEnergies[cycle] = localEnergy;
                }
                //        }
                // update energies
                m_energy += localEnergy;
                m_energySquared += localEnergy*localEnergy;
            } else {
                checkTerminalization(localEnergy);
            }
        }  //  end of loop over particles
    }
    m_energy /= (nCycles * nParticles);
    m_energySquared /= (nCycles * nParticles);
}
