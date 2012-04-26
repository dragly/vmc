#include "montecarlometropolishastings.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"
#include "../config.h"

MonteCarloMetropolisHastings::MonteCarloMetropolisHastings(Config *config) :
    MonteCarlo(config),
    m_nParticles(config->nParticles()),
    m_nDimensions(config->nDimensions()),
    charge(config->charge()),
    rank(config->rank()),
    step_length(config->stepLength()),
    wave(config->wave()),
    hamiltonian(config->hamiltonian())
{
    // every node has its own seed for the random numbers
    idum = -1-rank;
    // allocate matrices which contain the position of the particles
    rOld = new vec2[ m_nParticles];
    rNew = new vec2[ m_nParticles];
}

MonteCarloMetropolisHastings::~MonteCarloMetropolisHastings()
{
}

void MonteCarloMetropolisHastings::quantumForce(vec2 rPosition[], vec2 &forceVector) {
    double waveValue = m_config->wave()->wave(rPosition);
    m_config->wave()->gradient(rPosition, forceVector);
    forceVector = 2 * forceVector / waveValue;
}

void MonteCarloMetropolisHastings::sample(int nCycles, double *energies, double *allEnergies)
{
    double wfnew = 0;
    double wfold = 0;
    double energy = 0;
    double energy2 = 0;
    double delta_e = 0;
    double diffConstant = 1;
    // initialisations of variational parameters and energies
    energy = energy2 = 0; delta_e=0;
    //  initial trial position, note calling with alpha
    for (int i = 0; i < m_nParticles; i++) {
        for (int j=0; j < m_nDimensions; j++) {
            rOld[i][j] = step_length*(ran2(&idum)-0.5);
        }
    }
    wfold = wave->wave(rOld);
    wave->gradient(rOld, waveGradientOld);
    quantumForce(rOld, forceVectorOld);
    // loop over monte carlo cycles
    for (int cycle = 1; cycle <= nCycles; cycle++){
        // new trial position
        for (int i = 0; i < m_nParticles; i++) {
            quantumForce(rOld, forceVectorNew);
            rNew[i] = rOld[i] + diffConstant*forceVectorNew*step_length;
            for (int j=0; j < m_nDimensions; j++) {
                rNew[i][j] += simpleGaussRandom(&idum);
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < m_nParticles; k++) {
                if ( k != i) {
                    rNew[k] = rOld[k];
                }
            }
            wfnew = wave->wave(rNew);
            wave->gradient(rNew, waveGradientNew);
            double argument = 0;
            for( int j = 0; j < m_nDimensions; j++) {
                forceVectorSum[j] = forceVectorNew[j] + forceVectorOld[j];
                forceVectorDiff[j] = forceVectorOld[j] - forceVectorNew[j];
                positionDiff[j] = rNew[i][j] - rOld[i][j];
                argument += 0.5 * forceVectorSum[j] * (diffConstant * step_length / 2 * forceVectorDiff[j] - positionDiff[j]);
            }
            double waveFrac = wfnew*wfnew/(wfold*wfold);
            // The Metropolis test is performed by moving one particle at the time
            if(ran2(&idum) <= exp(argument) * waveFrac) {
                rOld[i]=rNew[i];
                waveGradientOld = waveGradientNew;
                wfold = wfnew;
            }
        }  //  end of loop over particles
        // compute local energy
        delta_e = hamiltonian->energy(wave, rOld);
        // save all energies on last variate
        //        if(variate==max_variations){
        allEnergies[cycle] = delta_e;
        //        }
        // update energies
        energy += delta_e;
        energy2 += delta_e*delta_e;
    }   // end of loop over MC trials
    // return the energy and the energy squared
    energies[0] = energy;
    energies[1] = energy2;
}
