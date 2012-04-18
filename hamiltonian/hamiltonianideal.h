#ifndef HAMILTONIANIDEAL_H
#define HAMILTONIANIDEAL_H
#define   ZERO       1.0E-10
#include <armadillo>
using namespace arma;

#include "hamiltonian.h"

class HamiltonianIdeal : public Hamiltonian
{
public:
    HamiltonianIdeal(int number_particles, int dimension, double charge);
    double energy(WaveFunction *wave, const vec2 r[]);
private:
    int number_particles;
    int dimension;
    double charge;
    double kineticEnergy(WaveFunction *wave, const vec2 r[], double wfold);
};

#endif // HAMILTONIANIDEAL_H
