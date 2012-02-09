#ifndef HAMILTONIANSIMPLE_H
#define HAMILTONIANSIMPLE_H
#define   ZERO       1.0E-10
#include "hamiltonian.h"
class HamiltonianSimple : public Hamiltonian
{
public:
    HamiltonianSimple(int number_particles, int dimension, double charge);
    double energy(WaveFunction *wave, double **r, double wfold);
    void setAnalyticalKineticEnergy(bool val);
private:
    int number_particles;
    int dimension;
    double charge;
    double kineticEnergy(WaveFunction *wave, double **r, double wfold);
    bool analyticalKineticEnergy;
};

#endif // HAMILTONIANSTANDARD_H
