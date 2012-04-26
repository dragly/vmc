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
    double energy(WaveFunction *wave, vec2 r[]);
    ~HamiltonianIdeal();
private:
    int number_particles;
    int dimension;
    double charge;
//    double kineticEnergy(WaveFunction *wave, vec2 r[], double wfold);

    vec2 *r_plus;
    vec2 *r_minus;
};

#endif // HAMILTONIANIDEAL_H
