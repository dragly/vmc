#include "montecarlometropolishastings.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"
#include "../config.h"

MonteCarloMetropolisHastings::MonteCarloMetropolisHastings(Config *config) :
    MonteCarlo(config),
    rank(config->rank()),
    step_length(config->stepLength()),
    wave(config->wave()),
    hamiltonian(config->hamiltonian())
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[ nParticles];
    rNew = new vec2[ nParticles];
}

MonteCarloMetropolisHastings::~MonteCarloMetropolisHastings()
{
}

void MonteCarloMetropolisHastings::quantumForce(vec2 rPosition[], vec2 &forceVector) {
    double waveValue = config->wave()->evaluate(rPosition);
    config->wave()->gradient(rPosition, forceVector);
    forceVector = 2 * forceVector / waveValue;
}

void MonteCarloMetropolisHastings::sample(int nCycles)
{
    double wfnew = 0;
    double wfold = 0;
    m_energy = 0;
    m_energySquared = 0;
    terminalizationSum = 0;
    terminalizationNum = 1;
    double localEnergy = 0;
    double diffConstant = 1;
    // initialisations of variational parameters and energies
    m_energy = m_energySquared = 0; localEnergy=0;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = step_length*(ran2(idum)-0.5);
        }
    }
    wfold = wave->evaluate(rOld);
    wave->gradient(rOld, waveGradientOld);
    quantumForce(rOld, forceVectorOld);
    // loop over monte carlo cycles
    for (int cycle = 1; cycle <= nCycles; cycle++){
        // new trial position
        for (int i = 0; i < nParticles; i++) {
            quantumForce(rOld, forceVectorNew);
            rNew[i] = rOld[i] + diffConstant*forceVectorNew*step_length;
            for (int j=0; j < nDimensions; j++) {
                rNew[i][j] += simpleGaussRandom(idum);
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < nParticles; k++) {
                if ( k != i) {
                    rNew[k] = rOld[k];
                }
            }
            wfnew = wave->evaluate(rNew);
            wave->gradient(rNew, waveGradientNew);
            double argument = 0;
            for( int j = 0; j < nDimensions; j++) {
                forceVectorSum[j] = forceVectorNew[j] + forceVectorOld[j];
                forceVectorDiff[j] = forceVectorOld[j] - forceVectorNew[j];
                positionDiff[j] = rNew[i][j] - rOld[i][j];
                argument += 0.5 * forceVectorSum[j] * (diffConstant * step_length / 2 * forceVectorDiff[j] - positionDiff[j]);
            }
            double waveFrac = wfnew*wfnew/(wfold*wfold);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(idum) <= exp(argument) * waveFrac) {
                rOld[i]=rNew[i];
                waveGradientOld = waveGradientNew;
                wfold = wfnew;
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
