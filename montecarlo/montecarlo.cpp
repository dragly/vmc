#include "montecarlo.h"

#include "montecarlostandard.h"
#include "montecarlometropolishastings.h"

MonteCarlo::MonteCarlo(Config *config) :
    config(config),
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    m_energy(0),
    m_energySquared(0),
    idum(config->idum()),
    terminalizationSum(0),
    terminalizationNum(0),
    terminalized(false),
    prevTerminalizationAverage(999999),
    terminalizationTrials(0),
    wave(config->wave()),
    hamiltonian(config->hamiltonian())
{
}

MonteCarlo* MonteCarlo::fromName(string monteCarloClass, Config *config)
{
    if(monteCarloClass == "MonteCarloStandard") {
        return new MonteCarloStandard(config);
    } else if(monteCarloClass == "MonteCarloMetropolisHastings") {
        return new MonteCarloMetropolisHastings(config);
    } else {
        return 0;
    }
}

void MonteCarlo::checkTerminalization(double localEnergy) {
    if(!(cycle % 500)) {
        double terminalizationAverage = terminalizationSum / terminalizationNum;
        diffAverage = fabs(terminalizationAverage - prevTerminalizationAverage);
        if(diffAverage < 2) {
            terminalized = true;
            cycle = 0;
        }
        prevTerminalizationAverage = terminalizationAverage;
        terminalizationSum = 0;
        terminalizationNum = 1;
        terminalizationTrials++;
    }
    terminalizationSum += localEnergy;
    terminalizationNum++;
}
