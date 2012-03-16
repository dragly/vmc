#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <string>
#include <armadillo>
using namespace arma;

using namespace std;

class WaveFunction;
class Config;

class Hamiltonian
{
public:
    Hamiltonian();
    virtual double energy(WaveFunction *wave, vec2 *r) = 0;
    static Hamiltonian *fromName(string hamiltonianClass, Config *config, double charge);

};

#endif // HAMILTONIAN_H
