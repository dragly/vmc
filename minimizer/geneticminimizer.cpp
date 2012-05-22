#include "geneticminimizer.h"
#include "../inih/ini.h"
#include "../walker/evolutionarywalker.h"
#include "../random.h"

GeneticMinimizer::GeneticMinimizer(Config *config) :
    Minimizer(config),
    nSamples(100),
    nSamplesStart(100),
    nSamplesEnd(100),
    nCycles(100),
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    stepLength(config->stepLength()),
    monteCarlo(config->monteCarlo()),
    wave(config->wave()),
    hamiltonian(config->hamiltonian())
{
}

void GeneticMinimizer::loadConfiguration(INIParser *settings)
{
    nGenes = settings->GetInteger("GeneticMinimizer","nGenes", 2);
    nIndividuals = settings->GetInteger("GeneticMinimizer","nIndividuals", 16);
    nPopulations = settings->GetInteger("GeneticMinimizer","nPopulations", 2);
    nCycles = settings->GetInteger("GeneticMinimizer","nCycles", 20000);
    nSamplesStart = settings->GetInteger("GeneticMinimizer","nSamplesStart", 100);
    nSamplesEnd = settings->GetInteger("GeneticMinimizer","nSamplesEnd", 1000);
    rescaleCycles = settings->GetInteger("GeneticMinimizer","rescaleCycles", 2);
    double lowRescaleLimit = settings->GetDouble("GeneticMinimizer","lowRescaleLimit", 0.1);
    double highRescaleLimit = settings->GetDouble("GeneticMinimizer","highRescaleLimit", 2.0);
    setRescaleLimits(lowRescaleLimit, highRescaleLimit);
    setPopulationData(nGenes, nIndividuals, nPopulations);
}

void GeneticMinimizer::runMinimizer()
{
    ofstream dataFile;
    dataFile.open("minimizer-evolutionary.dat");
    energies = new vec[nPopulations];
    for(int i = 0; i < nPopulations; i++) {
        energies[i] = zeros<vec>(nIndividuals);
    }
    for(int acycle = 0; acycle < nCycles; acycle++) {
        allBestValue = INFINITY;
        evolve(1, nCycles / 10.);
        double energySum = 0;
        int nEnergySamples = 0;
        double alphaSum = 0;
        double betaSum = 0;
        for(int i = 0; i < nPopulations; i++) {
            for(uint j = 0; j < bestIndices[i].size() / 2; j++) {
                int index = bestIndices[i][j];
                double energy = energies[i][index];
                nEnergySamples++;
                energySum += energy;
                alphaSum += fabs(populations[i][index][0]);
                betaSum += fabs(populations[i][index][1]);
            }
        }
        // TODO output variance if possible
        dataFile << alphaSum / nEnergySamples << " ";
        dataFile << betaSum / nEnergySamples << " ";
        dataFile << energySum / nEnergySamples << " ";
        dataFile << nEnergySamples << std::endl;
        std::cout << cycle << ": Mean energy " <<  left << setw(10) << energySum / nEnergySamples << " using params " <<  left << setw(10) << alphaSum / nEnergySamples << " " <<  left << setw(10) << betaSum / nEnergySamples << " with best " <<  left << setw(10) << energies[allBestPopulationIndex][allBestIndex] << " using " <<  left << setw(10) << populations[allBestPopulationIndex][allBestIndex][0] << ", " << left << setw(10) <<  populations[allBestPopulationIndex][allBestIndex][1] << " @ scale " <<  left << setw(10) << scale << " samples: " << nSamples << " x " << nIndividuals << " x " << nPopulations << std::endl;
    }
    dataFile.close();

    delete [] energies;
}

double GeneticMinimizer::fitness(vec &coefficients, int population, int individual)
{
    double parameters[2];
    parameters[0] = fabs(coefficients[0]);
    parameters[1] = fabs(coefficients[1]);
    wave->setParameters(parameters);
    monteCarlo->setThermalizationEnabled(false);

    monteCarlo->sample(nSamples);
    double energy = monteCarlo->energy();

    if(population == 0 && individual == 0) {
        nSamples = nSamplesStart + (nSamplesEnd - nSamplesStart) * pow((double)cycle, 20) / pow((double)nCycles, 20);

        // reduce both scale limits to a tenth of their initial values
        double cycleFactor = (double)cycle / (double)nCycles;
        lowScaleLimit = initialLowScaleLimit - cycleFactor * initialLowScaleLimit * 9. / 10.;
        highScaleLimit = initialHighScaleLimit - cycleFactor * initialHighScaleLimit * 9. / 10.;
    }

    if(energy <= 0) {
        energy = INFINITY;
    }
    energies[population][individual] = energy;
    return energy;
}

