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
