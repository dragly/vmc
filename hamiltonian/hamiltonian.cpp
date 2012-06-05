#include "hamiltonian.h"
#include "hamiltoniansimple.h"
#include "hamiltonianideal.h"
#include "../config.h"

/*!
  * Constructs the Hamiltonian.
  */
Hamiltonian::Hamiltonian(Config *config) :
    m_nParticles(config->nParticles()),
    m_nDimensions(config->nDimensions()),
    m_interactionEnabled(config->interactionEnabled())
{
}

/*!
  * \brief Hamiltonian::resetTotalEnergies clears the total potential and kinetic energy counters
  */
void Hamiltonian::resetTotalEnergies()
{
    m_totalPotentialEnergy = 0;
    m_totalKineticEnergy = 0;
}

/*!
  * \brief Hamiltonian::outputTotals prints out a string with the relative amounts of total and potential energy
  */
void Hamiltonian::outputTotals()
{
    double totalEnergy = totalKineticEnergy() + totalPotentialEnergy();
    std::cout << "Kinetic energy: " << totalKineticEnergy() / (totalEnergy) << " potential energy: " << totalPotentialEnergy() / (totalEnergy) << std::endl;
}

/*!
  * \brief Hamiltonian::totalsString creates a string with total potential and kinetic energy for file printing
  */
string Hamiltonian::totalsString()
{
    double totalEnergy = totalKineticEnergy() + totalPotentialEnergy();
    stringstream myString;
    myString << totalKineticEnergy() / (totalEnergy) << " " << totalPotentialEnergy() / (totalEnergy);
    return myString.str();
}


/*!
  * \brief Hamiltonian::fromName returns an Hamiltonian class based on the class name
  */
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
