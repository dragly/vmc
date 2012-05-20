#include "evolver.h"

#include <armadillo>
#include "../random.h"

using namespace arma;

Evolver::Evolver()
{
    constructor(32,32,10);
}

Evolver::Evolver(int nGenes, int nIndividuals, int nPopulations)
{
    constructor(nIndividuals, nGenes, nPopulations);
}

void Evolver::constructor(int nGenes, int nIndividuals, int nPopulations) {
    idum = -1;
    currentCycle = 0;
    this->nIndividuals = nIndividuals;
    this->nGenes = nGenes;
    this->nPopulations = nPopulations;

    lowScaleLimit = 1e-2;
    highScaleLimit = 1e5;

    populations = new vec*[nPopulations];
    values = new vec[nPopulations];
    bestIndices = new uvec[nPopulations];
    for(int i = 0; i < nPopulations; i++) {
        populations[i] = new vec[nIndividuals];
        values[i] = zeros<vec>(nIndividuals);
        bestIndices[i] = zeros<uvec>(nIndividuals);
        for(int j = 0; j < nIndividuals; j++) {
            populations[i][j] = randu<vec>(nGenes);
            values[i][j] = INFINITY;
        }
    }
    allBestIndex = -1;
    allBestPopulationIndex = -1;
    allBestValue = INFINITY;
    cyclesSinceLastImprovement = 0;
    cyclesSinceLastRescale = 0;

    std::cout << "Done constructing" << std::endl;
}

Evolver::~Evolver()
{
    delete [] populations;
    delete [] values;
}

void Evolver::updateBest()
{
    for(int i = 0; i < nPopulations; i++) {
        for(int j = 0; j < nIndividuals; j++) {
            vec &genes = populations[i][j];
            double value = fitness(genes);
            values[i][j] = value;
            if(value < allBestValue) {
                allBestValue = value;
                allBestIndex = j;
                allBestPopulationIndex = i;
                cyclesSinceLastImprovement = 0;
            }
        }
        bestIndices[i] = sort_index(values[i]);
    }
}

void Evolver::rescale()
{
    double expScale = lowScaleLimit + (highScaleLimit - lowScaleLimit) * ran2(&idum);
    scale = pow(10,expScale);
}

void Evolver::evolve(int nCycles, int populationMatchingPeriod)
{
    updateBest();
    for(int cycle = 0; cycle < nCycles; cycle++) {
        cyclesSinceLastImprovement++;
        cyclesSinceLastRescale++;
        // Mating of the first two quarters of individuals and add to the third quarter
        for(int i = 0; i < nPopulations; i++) {
            for(int j = 0; j < nIndividuals / 2; j += 2) {
                uint parent1Index = bestIndices[i][j];
                uint parent2Index = bestIndices[i][j + 1];
                uint childIndex = bestIndices[i][j / 2 + nIndividuals / 2];
//                std::cout << j << " " << j + 1 << " " << j / 2 + nIndividuals / 2 << std::endl;
                for(int k = 0; k < nGenes; k++) {
                    // a random selection of half of the genes comes from one parent
                    int parentIndex;
                    if(ran2(&idum) > 0.5) {
                        parentIndex = parent1Index;
                    } else {
                        parentIndex = parent2Index;
                    }
                    populations[i][childIndex][k] = populations[i][parentIndex][k];
                }
            }
        }

        // Mutate one gene in the second half of individuals
        for(int i = 0; i < nPopulations; i++) {
            for(int j = 0; j < nIndividuals / 2; j++) {
                uint individualIndex = bestIndices[i][nIndividuals / 2 + j];
                vec &genes = populations[i][individualIndex]; // note the use of reference!
                for(int k = 0; k < nGenes / 4; k++) {
                    int randomGene = ran2(&idum) * nGenes;
    //                double gauss = simpleGaussRandom(&idum);
                    double gauss = ran2(&idum);
                    genes[randomGene] += gauss * scale;
                }
            }
        }

        // Introduction of completely new individuals
        for(int i = 0; i < nPopulations; i++) {
            for(int j = 0; j < nIndividuals / 4; j++) {
                uint individualIndex = bestIndices[i][3 * nIndividuals / 4 + j];
                populations[i][individualIndex] = randu<vec>(nGenes) * scale;
            }
        }

        // Merge populations every populationMatching step. Adds the best from one population to the other population.
        if(!(cycle % populationMatchingPeriod)) {
            for(int i = 0; i < nPopulations / 2; i++) {
                int population1Index = (int)(ran2(&idum) * (nPopulations/2));
                int population2Index = (int)((nPopulations/2) + ran2(&idum) * (nPopulations/2));
                if(population1Index != population2Index) {
                    for(int j = 0; j < nIndividuals / 2; j++) {
                        uint bestIndex1 = bestIndices[population1Index][j];
                        uint bestIndex2 = bestIndices[population2Index][j];
                        uint worstIndex1 = bestIndices[population1Index][j + nIndividuals/2];
                        uint worstIndex2 = bestIndices[population2Index][j + nIndividuals/2];
                        populations[population1Index][worstIndex1] = populations[population2Index][bestIndex2];
                        populations[population2Index][worstIndex2] = populations[population1Index][bestIndex1];
                    }
                }
            }
        }

        // Testing
        updateBest();

        if(cyclesSinceLastImprovement > 50) {
            rescale();
            cyclesSinceLastImprovement = 20;
            cyclesSinceLastRescale = 0;
        } else if(cyclesSinceLastImprovement > 45) {
            scale = lastWorkingScale;
        } else if(cyclesSinceLastImprovement > 25) {
            scale *= 1.2;
        }
        if(cyclesSinceLastRescale > 100) {
            lastWorkingScale = scale;
            rescale();
        }


        std::cout << "All best value " << currentCycle << " " << allBestValue << " " << scale << std::endl;
        currentCycle++;
    }
//    std::cout << "All best genes " << populations[allBestPopulationIndex][allBestIndex] << std::endl;
    std::cout << "All best value " << allBestValue << std::endl;
    allBestGenes = populations[allBestPopulationIndex][allBestIndex];
}
