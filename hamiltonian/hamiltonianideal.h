#ifndef HAMILTONIANIDEAL_H
#define HAMILTONIANIDEAL_H
#define   ZERO       1.0E-10
#include <armadillo>
using namespace arma;

#include "hamiltonian.h"

class HamiltonianIdeal : public Hamiltonian
{
public:
    HamiltonianIdeal(Config *config);
    double kineticEnergy(WaveFunction *wave, vec2 r[]);
    double potentialEnergy(WaveFunction *wave, vec2 r[]);
    double externalPotentialEnergy(WaveFunction *wave, vec2 r[]);
    double interactionPotentialEnergy(WaveFunction *wave, vec2 r[]);

    double totalExternalPotentialEnergy() {
        return m_totalExternalPotentialEnergy;
    }
    double totalInteractionPotentialEnergy() {
        return m_totalInteractionPotentialEnergy;
    }
    void resetTotalEnergies();
    string totalsString();
    void outputTotals();

private:
    double omega;
    double m_totalExternalPotentialEnergy;
    double m_totalInteractionPotentialEnergy;
};

#endif // HAMILTONIANIDEAL_H
