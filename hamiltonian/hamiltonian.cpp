#include "hamiltonian.h"
#include "hamiltoniansimple.h"
#include "hamiltonianideal.h"
#include "../config.h"

Hamiltonian::Hamiltonian()
{
}

Hamiltonian *Hamiltonian::fromName(string hamiltonianClass, Config *config, double charge)
{
    if(hamiltonianClass == "HamiltonianSimple") {
        return new HamiltonianSimple(config->nParticles(), config->nDimensions(), charge);
    } else if(hamiltonianClass == "HamiltonianIdeal") {
        return new HamiltonianIdeal(config->nParticles(), config->nDimensions(), charge);
    } else {
        return 0;
    }
}
