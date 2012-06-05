#ifndef HAMILTONIANSIMPLE_H
#define HAMILTONIANSIMPLE_H
#define   ZERO       1.0E-10
#include <armadillo>
using namespace arma;
#include "hamiltonian.h"

/*!
  * \brief A simple Hamiltonian used during development and to test HamiltonianIdeal
  */
class HamiltonianSimple : public Hamiltonian
{
public:
    HamiltonianSimple(Config *config);
    double potentialEnergy(WaveFunction *wave, vec2 r[]);
    double kineticEnergy(WaveFunction *wave, vec2 r[]);
    void setAnalyticalKineticEnergy(bool val);
private:
};

#endif // HAMILTONIANSTANDARD_H
