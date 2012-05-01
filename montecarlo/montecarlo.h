#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "../wavefunction/wavefunction.h"
#include "../hamiltonian/hamiltonian.h"

class INIReader;

class MonteCarlo
{
public:
    MonteCarlo(Config *config);

    virtual void sample(int numberCycles, bool storeEnergies) = 0;
    virtual void loadConfiguration(INIReader *settings) {
        (void)settings;
    }
    static MonteCarlo *fromName(string monteCarloClass, Config *config);

    double energy() {
        return m_energy;
    }
    double energySquared() {
        return m_energySquared;
    }

protected:
    Config *m_config;
    double m_energy;
    double m_energySquared;
    double *m_allEnergies;
    long *idum;
};

#endif // MONTECARLO_H
