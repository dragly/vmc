#include "minimizerevolutionary.h"
#include "../inih/ini.h"

MinimizerEvolutionary::MinimizerEvolutionary(Config *config) :
    Minimizer(config),
    nSamples(100)
{
    monteCarlo = config->monteCarlo();
    wave = config->wave();
    hamiltonian = config->hamiltonian();
}

void MinimizerEvolutionary::runMinimizer()
{
    ofstream dataFile;
    dataFile.open("minimizer-evolutionary.dat");
    energies = new vec[nPopulations];
    for(int i = 0; i < nPopulations; i++) {
        energies[i] = zeros<vec>(nIndividuals);
    }
    for(int cycle = 0; cycle < 1000; cycle++) {
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
                alphaSum += populations[i][index][0];
                betaSum += populations[i][index][1];
            }
        }
        // TODO output variance if possible
        dataFile << alphaSum / nEnergySamples << "\t";
        dataFile << betaSum / nEnergySamples << "\t";
        dataFile << energySum / nEnergySamples << std::endl;
        std::cout << "Mean energy " <<  left << setw(12) << energySum / nEnergySamples << " using params " <<  left << setw(12) << alphaSum / nEnergySamples << " " <<  left << setw(12) << betaSum / nEnergySamples << " with best " <<  left << setw(12) << energies[allBestPopulationIndex][allBestIndex] << " using " <<  left << setw(12) << populations[allBestPopulationIndex][allBestIndex][0] << ", " << left << setw(12) <<  populations[allBestPopulationIndex][allBestIndex][1] << " @ scale " <<  left << setw(12) << scale << " samples: " << nSamples << " x " << nIndividuals << " x " << nPopulations << std::endl;
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
    parameters[0] = coefficients[0];
    parameters[1] = coefficients[1];
    wave->setParameters(parameters);
    monteCarlo->setTerminalizationEnabled(false);
//    std::cout << "Sampling with " << nSamples << std::endl;
    monteCarlo->sample(nSamples);
    if(population == individual == 0) {
//        nSamples *= 1.0001;
        nSamples += (nSamplesEnd / (double)nSamplesStart) / 100;
    }
    // update the energy average and its squared
    double energy = monteCarlo->energy();
//    std::cout << left << setw(12) << setprecision(10) << "Population" << population << ", individual" << individual << ", scale " << setw(12) << scale << ", testing parameters " << setw(12) << parameters[0] << ", " << setw(12) << setprecision(10) << parameters[1] << ", got energy " << setw(12) << setprecision(10) << energy << std::endl;
    //    std::cout << energy << std::endl;
    //    cumulativeEnergySquared(i,j) = m_monteCarlo->energySquared();
    //    std::cout << "Got energy of " << cumulativeEnergy(i,j) << std::endl;
    //    parameter0Map(i,j) = parameters[0];
    //    parameter1Map(i,j) = parameters[1];
    //    parameters[1] += betaStep;
    energies[population][individual] = energy;
    if(energy == 0) {
        return INFINITY;
    }
    return energy;
}

