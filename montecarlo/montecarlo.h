#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "../wavefunction.h"
#include "../hamiltonian/hamiltonian.h"

class MonteCarlo
{
public:
    MonteCarlo(Config *config);

    virtual void sample(int number_cycles, double *energies, double *allEnergies) = 0;
    virtual void loadConfiguration(INIReader *settings) {
        (void)settings;
    }
    static MonteCarlo *fromName(string monteCarloClass, Config *config);
protected:
    Config *m_config;
};

#endif // MONTECARLO_H
