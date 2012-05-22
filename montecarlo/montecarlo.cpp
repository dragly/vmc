#include "montecarlo.h"

#include "standardmontecarlo.h"
#include "metropolishastingsmontecarlo.h"
#include "../random.h"

MonteCarlo::MonteCarlo(Config *config) :
    config(config),
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    m_energy(0),
    m_energySquared(0),
    idumMC(config->idum()),
    terminalizationSum(0),
    terminalizationNum(0),
    terminalized(false),
    prevTerminalizationAverage(999999),
    terminalizationTrials(0),
    wave(config->wave()),
    hamiltonian(config->hamiltonian()),
    recordMoves(false),
    nMoves(1),
    stepLength(config->stepLength())
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[nParticles];
    rNew = new vec2[nParticles];
    randomizePositions();
}

MonteCarlo::~MonteCarlo() {
    delete [] m_moves;
    delete [] rOld;
    delete [] rNew;
}

void MonteCarlo::randomizePositions() {
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = rNew[i][j] = stepLength*(ran2(idumMC)-0.5);
        }
    }
}

MonteCarlo* MonteCarlo::fromName(string monteCarloClass, Config *config)
{
    if(monteCarloClass == "MonteCarloStandard") {
        return new StandardMonteCarlo(config);
    } else if(monteCarloClass == "MonteCarloMetropolisHastings") {
        return new MonteCarloMetropolisHastings(config);
    } else {
        return 0;
    }
}

void MonteCarlo::setRecordMoves(bool arg, int nMoves) {
    this->recordMoves = arg;
    this->nMoves = nMoves;
    delete [] m_moves;
    m_moves = new vec2*[nMoves];
    for(int i = 0; i < nMoves; i++) {
        m_moves[i] = new vec2[nParticles];
    }
}

void MonteCarlo::checkTerminalization(double localEnergy) {
    if(!(cycle % 1000)) {
        double terminalizationAverage = terminalizationSum / terminalizationNum;
        diffAverage = fabs(terminalizationAverage - prevTerminalizationAverage);
        if(diffAverage < 2 && terminalizationTrials > 50) {
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
