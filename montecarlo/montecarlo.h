#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "../wavefunction/wavefunction.h"
#include "../hamiltonian/hamiltonian.h"

class INIParser;

class MonteCarlo
{
public:
    MonteCarlo(Config *config);

    virtual void sample(int nCycles) = 0;
    virtual void loadConfiguration(INIParser *settings) {
        (void)settings;
    }
    static MonteCarlo *fromName(string monteCarloClass, Config *config);

    double energy() {
        return m_energy;
    }
    double energySquared() {
        return m_energySquared;
    }
    void setTerminalizationEnabled(bool arg) {
        terminalized = !arg;
    }

    void checkTerminalization(double localEnergy);
protected:
    Config *config;
    int nParticles;
    int nDimensions;
    double m_energy;
    double m_energySquared;
    double *m_allEnergies;
    long *idum;

    double terminalizationSum;
    int terminalizationNum;
    bool terminalized;
    double diffAverage;

    double prevTerminalizationAverage;

    int terminalizationTrials;

    WaveFunction* wave;
    Hamiltonian *hamiltonian;

    int cycle;

    bool storeEnergies;
};

#endif // MONTECARLO_H
