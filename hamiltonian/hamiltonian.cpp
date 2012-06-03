#include "hamiltonian.h"
#include "hamiltoniansimple.h"
#include "hamiltonianideal.h"
#include "../config.h"

Hamiltonian::Hamiltonian(Config *config) :
    m_nParticles(config->nParticles()),
    m_nDimensions(config->nDimensions()),
    m_interactionEnabled(config->interactionEnabled())
{
}

void Hamiltonian::resetTotalEnergies()
{
    m_totalPotentialEnergy = 0;
    m_totalKineticEnergy = 0;
}

void Hamiltonian::outputTotals()
{
    double totalEnergy = totalKineticEnergy() + totalPotentialEnergy();
    std::cout << "Kinetic energy: " << totalKineticEnergy() / (totalEnergy) << " potential energy: " << totalPotentialEnergy() / (totalEnergy) << std::endl;
}

string Hamiltonian::totalsString()
{
    double totalEnergy = totalKineticEnergy() + totalPotentialEnergy();
    stringstream myString;
    myString << totalKineticEnergy() / (totalEnergy) << " " << totalPotentialEnergy() / (totalEnergy);
    return myString.str();
}


Hamiltonian *Hamiltonian::fromName(string hamiltonianClass, Config *config)
{
    if(hamiltonianClass == "HamiltonianSimple") {
        return new HamiltonianSimple(config);
    } else if(hamiltonianClass == "HamiltonianIdeal") {
        return new HamiltonianIdeal(config);
    } else {
        return 0;
    }
}
