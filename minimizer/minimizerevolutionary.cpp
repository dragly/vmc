#include "minimizerevolutionary.h"
#include "../inih/ini.h"
#include "../walker/evolutionarywalker.h"
#include "../random.h"

MinimizerEvolutionary::MinimizerEvolutionary(Config *config) :
    Minimizer(config),
    nSamples(100),
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions()),
    stepLength(config->stepLength()),
    nCycles(10000)
{
    monteCarlo = config->monteCarlo();
    wave = config->wave();
    hamiltonian = config->hamiltonian();
    walker = new EvolutionaryWalker(config);

    vec2 *positions = new vec2[config->nParticles()];
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            positions[i][j] = stepLength*(ran2(idum)-0.5);
        }
    }
    walker->initialize(positions);
}

void MinimizerEvolutionary::runMinimizer()
{
    ofstream dataFile;
    dataFile.open("minimizer-evolutionary.dat");
    energies = new vec[nPopulations];
    for(int i = 0; i < nPopulations; i++) {
        energies[i] = zeros<vec>(nIndividuals);
    }
    for(int acycle = 0; acycle < nCycles; acycle++) {
        allBestValue = INFINITY;
        evolve(1,5);
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
        dataFile << energySum / nEnergySamples << std::endl;
        std::cout << cycle << ": Mean energy " <<  left << setw(12) << energySum / nEnergySamples << " using params " <<  left << setw(12) << alphaSum / nEnergySamples << " " <<  left << setw(12) << betaSum / nEnergySamples << " with best " <<  left << setw(12) << energies[allBestPopulationIndex][allBestIndex] << " using " <<  left << setw(12) << populations[allBestPopulationIndex][allBestIndex][0] << ", " << left << setw(12) <<  populations[allBestPopulationIndex][allBestIndex][1] << " @ scale " <<  left << setw(12) << scale << " samples: " << nSamples << " x " << nIndividuals << " x " << nPopulations << std::endl;
    }
    dataFile.close();

    delete [] energies;
}

void MinimizerEvolutionary::loadConfiguration(INIParser *settings)
{
    nIndividuals = settings->GetInteger("MinimizerEvolutionary","nIndividuals",16);
    nPopulations = settings->GetInteger("MinimizerEvolutionary","nPopulations",2);
    nGenes = settings->GetInteger("MinimizerEvolutionary","nGenes",2);
    constructor(nGenes, nPopulations, nIndividuals);
    wave = config->wave();
    hamiltonian = config->hamiltonian();
}

void MinimizerEvolutionary::startEvolution()
{
}

double MinimizerEvolutionary::fitness(vec &coefficients, int population, int individual)
{
    double parameters[2];
    parameters[0] = fabs(coefficients[0]);
    parameters[1] = fabs(coefficients[1]);
    wave->setParameters(parameters);
    monteCarlo->setTerminalizationEnabled(false);
    //    std::cout << "Sampling with " << nSamples << std::endl;
    //    double energySum = 0;
    //    int energySamples = 0;
    //    for(int i = 0; i < nSamples; i++) {
    //        walker->advance();
    //        energySum += walker->energy();
    //        energySamples += walker->changeInEnergySamples();
    //    }
    //    double energy = energySum / energySamples;
    monteCarlo->sample(nSamples);
    double energy = monteCarlo->energy();
    if(population == 0 && individual == 0) {
        //        nSamples *= 1.0001;
        nSamples = nSamplesStart + (nSamplesEnd - nSamplesStart) * pow(cycle, 3) / pow(nCycles, 3);
        lowScaleLimit -= 0.0005;
        highScaleLimit -= 0.0005;
    }
    // update the energy average and its squared
    //    std::cout << left << setw(12) << setprecision(10) << "Population" << population << ", individual" << individual << ", scale " << setw(12) << scale << ", testing parameters " << setw(12) << parameters[0] << ", " << setw(12) << setprecision(10) << parameters[1] << ", got energy " << setw(12) << setprecision(10) << energy << std::endl;
    //    std::cout << energy << std::endl;
    //    cumulativeEnergySquared(i,j) = m_monteCarlo->energySquared();
    //    std::cout << "Got energy of " << cumulativeEnergy(i,j) << std::endl;
    //    parameter0Map(i,j) = parameters[0];
    //    parameter1Map(i,j) = parameters[1];
    //    parameters[1] += betaStep;
    if(energy < 0) {
        energy = INFINITY;
    }
    energies[population][individual] = energy;
    if(energy == 0) {
        return INFINITY;
    }
    return energy;
}

