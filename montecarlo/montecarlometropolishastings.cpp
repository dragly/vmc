#include "montecarlometropolishastings.h"
#include <math.h>
#include "../matrix.h"
#include "../random.h"
#include "../utils.h"

MonteCarloMetropolisHastings::MonteCarloMetropolisHastings(WaveFunction *wave, Hamiltonian *hamiltonian, int nParticles, int nDimensions, double charge, int rank, double step_length) :
    MonteCarlo(wave, hamiltonian),
    m_nParticles(nParticles),
    m_nDimensions(nDimensions),
    charge(charge),
    rank(rank),
    step_length(step_length)
{
    // every node has its own seed for the random numbers
    idum = -1-rank;
    // allocate matrices which contain the position of the particles
    rOld = (double **) matrix( nParticles, nDimensions, sizeof(double));
    rNew = (double **) matrix( nParticles, nDimensions, sizeof(double));
    waveGradientOld = new double[nDimensions];
    waveGradientNew = new double[nDimensions];
    forceVectorNew = new double[m_nDimensions];
    forceVectorSum = new double[m_nDimensions];
    forceVectorDiff = new double[m_nDimensions];
    positionDiff = new double[m_nDimensions];
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = rNew[i][j] = 0;
        }
    }
}

MonteCarloMetropolisHastings::~MonteCarloMetropolisHastings()
{
    free_matrix((void **) rOld); // free memory
    free_matrix((void **) rNew); // free memory
    delete [] waveGradientOld ;
    delete [] waveGradientNew;
    delete [] forceVectorNew;
    delete [] forceVectorSum;
    delete [] forceVectorDiff;
    delete [] positionDiff;
}

void MonteCarloMetropolisHastings::quantumForce(double **rPosition, double *forceVector) {
    double waveValue = wave->wave(rPosition);
    wave->gradient(rPosition, forceVector);
    for(int j = 0; j < m_nDimensions; j++) {
        forceVector[j] = 2 * forceVector[j] / waveValue;
    }
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
            for (int j=0; j < m_nDimensions; j++) {
                rNew[i][j] = rOld[i][j] + diffConstant*forceVectorNew[j]*step_length + simpleGaussRandom(&idum);
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < m_nParticles; k++) {
                if ( k != i) {
                    for (int l=0; l < m_nDimensions; l++) {
                        rNew[k][l] = rOld[k][l];
                    }
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

            // The Metropolis test is performed by moving one particle at the time
            if(ran2(&idum) <= exp(argument) ) {
                for (int l=0; l < m_nDimensions; l++) {
                    rOld[i][l]=rNew[i][l];
                    waveGradientOld[l] = waveGradientNew[l];
                }
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
