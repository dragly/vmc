#include "montecarlo.h"

#include "montecarlostandard.h"
#include "montecarlometropolishastings.h"

MonteCarlo::MonteCarlo(Config *config) :
    m_config(config),
    m_energy(0),
    m_energySquared(0),
    idum(config->idum())
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
