#ifndef HAMILTONIANSIMPLE_H
#define HAMILTONIANSIMPLE_H
#define   ZERO       1.0E-10
#include <armadillo>
using namespace arma;
#include "hamiltonian.h"
class HamiltonianSimple : public Hamiltonian
{
public:
    HamiltonianSimple(int number_particles, int dimension, double charge);
    double energy(WaveFunction *wave, const vec2 r[]);
    void setAnalyticalKineticEnergy(bool val);
private:
    int number_particles;
    int dimension;
    double charge;
    double kineticEnergy(WaveFunction *wave, const vec2 r[], double wfold);
    double analyticalKineticEnergy(WaveFunction *wave, const vec2 r[], double wfold);
    bool useAnalyticalKineticEnergy;
};

#endif // HAMILTONIANSTANDARD_H
