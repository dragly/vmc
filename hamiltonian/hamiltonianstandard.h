#ifndef HAMILTONIANSTANDARD_H
#define HAMILTONIANSTANDARD_H
#include <armadillo>
using namespace arma;
#include "hamiltonian.h"
class HamiltonianStandard : public Hamiltonian
{
public:
    HamiltonianStandard(int number_particles, int dimension, double charge);
    double energy(WaveFunction *wave, vec2 *r, double wfold);
private:
    int number_particles;
    int dimension;
    double charge;
};

#endif // HAMILTONIANSTANDARD_H
