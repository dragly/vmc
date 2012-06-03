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
        double potEnergy = potentialEnergy(wave, r);
        double kinEnergy = kineticEnergy(wave, r);
        m_totalPotentialEnergy += potEnergy;
        m_totalKineticEnergy += kinEnergy;
        return potEnergy + kinEnergy;
    }
    virtual double potentialEnergy(WaveFunction *wave, vec2 r[]) = 0;
    virtual double kineticEnergy(WaveFunction *wave, vec2 r[]) = 0;
    virtual void resetTotalEnergies();
    virtual void outputTotals();
    virtual string totalsString();
    double totalPotentialEnergy () {
        return m_totalPotentialEnergy;
    }
    double totalKineticEnergy () {
        return m_totalKineticEnergy;
    }

    static Hamiltonian *fromName(string hamiltonianClass, Config *config);
protected:
    int m_nParticles;
    int m_nDimensions;
    bool m_interactionEnabled;
    double m_totalPotentialEnergy;
    double m_totalKineticEnergy;
};

#endif // HAMILTONIAN_H
