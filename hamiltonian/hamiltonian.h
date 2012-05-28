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
    Hamiltonian(Config *config);
    virtual double energy(WaveFunction *wave, vec2 r[]) {
        return potentialEnergy(wave, r) + kineticEnergy(wave, r);
    }
    virtual double potentialEnergy(WaveFunction *wave, vec2 r[]) = 0;
    virtual double kineticEnergy(WaveFunction *wave, vec2 r[]) = 0;

    static Hamiltonian *fromName(string hamiltonianClass, Config *config);
protected:
    int m_nParticles;
    int m_nDimensions;
    bool m_interactionEnabled;
};

#endif // HAMILTONIAN_H
