#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "../wavefunction/wavefunction.h"
#include "../hamiltonian/hamiltonian.h"

class INIParser;

class MonteCarlo
{
public:
    MonteCarlo(Config *config);
    ~MonteCarlo();

    virtual void sample(int nCycles) = 0;
    virtual void loadConfiguration(INIParser *settings);

    static MonteCarlo *fromName(string monteCarloClass, Config *config);

    double *allEnergies() {
        return m_allEnergies;
    }

    double energy() {
        return m_energy;
    }
    double energySquared() {
        return m_energySquared;
    }

    void setThermalizationEnabled(bool arg) {
        terminalized = !arg;
        terminalizationTrials = 0;
    }

    vec2 **moves() {
        return m_moves;
    }

    void setStoreEnergies(bool arg) {
        storeEnergies = arg;
    }

    void setOutputEnergies(bool arg) {
        outputEnergies = arg;
    }

    void checkTerminalization(double localEnergy);
    void setRecordMoves(bool arg, int nMoves);
    void recordMove(int i, int nthMove);
    void randomizePositions();

protected:
    Config *config;
    int nParticles;
    int nDimensions;

    double m_energy;
    double m_energySquared;
    double *m_allEnergies;
    long *idumMC;

    double terminalizationSum;
    int terminalizationNum;
    bool terminalized;
    double diffAverage;

    double prevTerminalizationAverage;

    int terminalizationTrials;

    WaveFunction* wave;
    Hamiltonian *hamiltonian;

    bool recordMoves;
    int nMoves;
    vec2 **m_moves;
    int move;

    int cycle;

    vec2 *rOld;
    vec2 *rNew;
    double stepLength;
    bool storeEnergies;
    double spawnRadius;

    bool sampleVariationalGradient;
    bool outputEnergies;
};

#endif // MONTECARLO_H
