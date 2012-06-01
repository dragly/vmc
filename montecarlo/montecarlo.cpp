#include "montecarlo.h"

#include "standardmontecarlo.h"
#include "metropolishastingsmontecarlo.h"
#include "../random.h"
#include "../inih/ini.h"

// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"

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
    stepLength(config->stepLength()),
    storeEnergies(false),
    spawnRadius(2)
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[nParticles];
    rNew = new vec2[nParticles];
    m_moves = new vec2*[1];
    randomizePositions();
}

MonteCarlo::~MonteCarlo() {
    delete [] m_moves;
    delete [] rOld;
    delete [] rNew;
}

void MonteCarlo::loadConfiguration(INIParser *settings)
{
    spawnRadius = settings->GetDouble("MonteCarlo", "spawnRadius", spawnRadius);
}

void MonteCarlo::randomizePositions() {
    std::cout << "Randomizing positions: " << std::endl;
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = rNew[i][j] = spawnRadius * gaussianDeviate(idumMC);
            std::cout << "rOld[" << i << "][" << j << "] = " << rOld[i][j] << std::endl;
        }
    }
}

MonteCarlo* MonteCarlo::fromName(string monteCarloClass, Config *config)
{
    if(monteCarloClass == "MonteCarloStandard") {
        return new StandardMonteCarlo(config);
    } else if(monteCarloClass == "MonteCarloMetropolisHastings") {
        return new MetropolisHastingsMonteCarlo(config);
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
        if(cycle >= 100000) {
            terminalized = true;
            std::cout << "Thermalized after " << cycle << " cycles." << std::endl;
            cycle = 0;
        }
        terminalizationTrials++;
    }
    terminalizationSum += localEnergy;
    terminalizationNum++;
}
