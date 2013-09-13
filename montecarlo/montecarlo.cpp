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
    spawnRadius(2),
    outputEnergies(false)
{
    // allocate matrices which contain the position of the particles
    rOld = new vec2[nParticles];
    rNew = new vec2[nParticles];
    m_moves = new vec2*[1];
    randomizePositions();
}

MonteCarlo::~MonteCarlo() {
    if(movesFile.is_open()) {
        movesFile.close();
    }
    delete [] m_moves;
    delete [] rOld;
    delete [] rNew;
}

void MonteCarlo::loadConfiguration(INIParser *settings)
{
    setSpawnRadius(settings->GetDouble("MonteCarlo", "spawnRadius", spawnRadius));
    scratchDir = settings->GetString("MonteCarlo", "scratchDir", "/scratch/positions");
}

void MonteCarlo::randomizePositions() {
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rOld[i][j] = rNew[i][j] = spawnRadius * gaussianDeviate(idumMC);
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

void MonteCarlo::setSpawnRadius(double arg)
{
    spawnRadius = arg;
    randomizePositions();
}

void MonteCarlo::setRecordMoves(bool arg, int nMoves, string fileName) {
    this->recordMoves = arg;
    if(nMoves > 0 && arg) {
        this->nMoves = nMoves;
        delete [] m_moves;
        m_moves = new vec2*[nMoves];
        for(int i = 0; i < nMoves; i++) {
            m_moves[i] = new vec2[nParticles];
        }
    }
    if(arg) {
        std::cout << "Opening moves file " << fileName << std::endl;
        movesFile.open(fileName.c_str(), ios::out | ios::binary);
        if(!movesFile.is_open()) {
            std::cerr << "Could not open file " << fileName << "!" << std::endl;
            exit(951);
        }
    }
}

void MonteCarlo::checkTerminalization(double localEnergy) {
    if(!(cycle % 10000)) {
        if(cycle >= 100000) {
            terminalized = true;
            hamiltonian->resetTotalEnergies();
            std::cout << "Thermalized after " << cycle << " cycles." << std::endl;
            cycle = 0;
            m_energy = 0;
        }
        terminalizationTrials++;
    }
    terminalizationSum += localEnergy;
    terminalizationNum++;
}

void MonteCarlo::writePositionToFile(vec2 &position) {
    if(movesFile.is_open()) {
        movesFile.write((char*)(&position[0]), nDimensions*sizeof(double));
    } else {
        std::cerr << "Movefile closed!" << std::endl;
        exit(945);
    }
}
