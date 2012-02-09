#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

class WaveFunction;

class Hamiltonian
{
public:
    Hamiltonian();
    virtual double energy(WaveFunction *wave, double **r, double wfold) = 0;

};

#endif // HAMILTONIAN_H
