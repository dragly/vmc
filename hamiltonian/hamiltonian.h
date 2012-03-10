#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <string>

using namespace std;

class WaveFunction;
class Config;

class Hamiltonian
{
public:
    Hamiltonian();
    virtual double energy(WaveFunction *wave, double **r) = 0;
    static Hamiltonian *fromName(string hamiltonianClass, Config *config, double charge);

};

#endif // HAMILTONIAN_H
