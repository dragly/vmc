#ifndef HAMILTONIANIDEAL_H
#define HAMILTONIANIDEAL_H
#define   ZERO       1.0E-10

#include "hamiltonian.h"

class HamiltonianIdeal : public Hamiltonian
{
public:
    HamiltonianIdeal(int number_particles, int dimension, double charge);
    double energy(WaveFunction *wave, double **r, double wfold);
private:
    int number_particles;
    int dimension;
    double charge;
    double kineticEnergy(WaveFunction *wave, double **r, double wfold);
};

#endif // HAMILTONIANIDEAL_H
