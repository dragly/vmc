#ifndef HAMILTONIANSIMPLE_H
#define HAMILTONIANSIMPLE_H
#include "hamiltonian.h"
class HamiltonianSimple : public Hamiltonian
{
public:
    HamiltonianSimple(int number_particles, int dimension, double charge);
    double energy(WaveFunction *wave, double **r, double alpha, double wfold);
private:
    int number_particles;
    int dimension;
    double charge;
};

#endif // HAMILTONIANSTANDARD_H
