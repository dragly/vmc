#ifndef MONTECARLOMETROPOLISHASTINGS_H
#define MONTECARLOMETROPOLISHASTINGS_H
#include <armadillo>
using namespace arma;

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction/wavefunction.h"

/*!
  * \brief Implements the Metropolis Hastings algortihm
  */
class MetropolisHastingsMonteCarlo : public MonteCarlo
{
public:
    MetropolisHastingsMonteCarlo(Config *config);
    void sample(int nCycles);

    ~MetropolisHastingsMonteCarlo();
private:
    int myRank;
//    vec waveGradientOld;
//    vec waveGradientNew;
    vec quantumForceNew;
    vec quantumForceOld;
//    vec forceVectorSum;
//    vec forceVectorDiff ;
    vec2 positionDiff ;
    Hamiltonian* hamiltonian;
    bool firstSample;
    double diffConstant;
};

#endif // MONTECARLOMETROPOLISHASTINGS_H
