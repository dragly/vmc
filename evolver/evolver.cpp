#include "evolver.h"

#include <armadillo>
#include "../random.h"

using namespace arma;

Evolver::Evolver()
{
    setPopulationData(32,32,10);
}

Evolver::Evolver(int nGenes, int nIndividuals, int nPopulations)
{
    setPopulationData(nGenes, nIndividuals, nPopulations);
}

void Evolver::setPopulationData(int nGenes, int nIndividuals, int nPopulations) {
    idum = new long;
    *idum = -1*time(NULL);
    cycle = 0;
    rescaleCycles = 50;
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
    veryFirst = true;
    cycle = 0;
}

Evolver::~Evolver()
{
    delete [] populations;
    delete [] values;
    delete [] bestIndices;
}

void Evolver::updateBest()
{
    for(int i = 0; i < nPopulations; i++) {
        for(int j = 0; j < nIndividuals; j++) {
            vec &genes = populations[i][j];
            double value = fitness(genes, i, j);
//            std::cout << "Value: " << value << std::endl;
            if(isnan(fabs(value))) {
                value = INFINITY;
            }
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
    if(allBestIndex == -1 || allBestPopulationIndex == -1) {
        std::cerr << "No fitness function returned anything better than nan or inf. Evolver cannot continue." << std::endl;
        exit(961);
    }
}

void Evolver::rescale()
{
    double lowLimit = log(lowScaleLimit);
    double highLimit = log(highScaleLimit);
    double expScale = lowLimit + (highLimit - lowLimit) * ran2(idum);
    scale = pow(10,expScale);
}

void Evolver::evolve(int nSteps, int populationMatchingPeriod)
{
    if(cycle == 0) {
        updateBest();
//        std::cout << "Done first update" << std::endl;
    }
    for(int localCycle = 0; localCycle < nSteps; localCycle++) {
        cyclesSinceLastImprovement++;
        cyclesSinceLastRescale++;
        // Mating of the first two quarters of individuals and add to the third quarter
        for(int i = 0; i < nPopulations; i++) {
            for(int j = 0; j < nIndividuals / 2.; j += 2) {
                uint parent1Index = bestIndices[i][j];
                uint parent2Index = bestIndices[i][j + 1];
                uint childIndex = bestIndices[i][j / 2 + nIndividuals / 2];
//                std::cout << j << " " << j + 1 << " " << j / 2 + nIndividuals / 2 << std::endl;
                for(int k = 0; k < nGenes; k++) {
                    // a random selection of half of the genes comes from one parent
                    int parentIndex;
                    if(ran2(idum) > 0.5) {
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
            for(int j = nIndividuals * 2. / 4.; j < nIndividuals; j++) {
                uint individualIndex = bestIndices[i][j];
                vec &genes = populations[i][individualIndex]; // note the use of reference!
                for(int k = 0; k < nGenes; k++) {
                    int randomGene = ran2(idum) * nGenes;
    //                double gauss = simpleGaussRandom(idum);
                    double gauss = 2 * ran2(idum) - 1;
                    genes[randomGene] += gauss * scale;
                }
            }
        }

        // Introduction of completely new individuals, generated from the best but with random additions
        for(int i = 0; i < nPopulations; i++) {
            for(int j = 7. * nIndividuals / 8.; j < nIndividuals; j++) {
                uint individualIndex = bestIndices[i][j];
                for(int k = 0; k < nGenes; k++) {
                    double gauss = 2 * ran2(idum) - 1;
                    populations[i][individualIndex][k] = allBestGenes[k] + gauss * scale;
                }
            }
        }

        // Merge populations every populationMatching step. Adds the best from one population to the other population.
        if(!(cycle % populationMatchingPeriod)) {
            for(int i = 0; i < nPopulations / 2; i++) {
                int population1Index = (int)(ran2(idum) * (nPopulations/2));
                int population2Index = (int)((nPopulations/2) + ran2(idum) * (nPopulations/2));
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

        if(cyclesSinceLastImprovement > rescaleCycles) {
            rescale();
            cyclesSinceLastImprovement = rescaleCycles * 4. / 10.;
            cyclesSinceLastRescale = 0;
        } else if(cyclesSinceLastImprovement > rescaleCycles * 9. / 10.) {
            scale = lastWorkingScale;
        } else if(cyclesSinceLastImprovement > rescaleCycles / 2.) {
            scale *= 1.2;
        }
        if(cyclesSinceLastRescale > rescaleCycles * 2) {
            lastWorkingScale = scale;
            rescale();
        }


//        std::cout << "All best value @ cycle " << cycle << ": " << allBestValue << " with scale " << scale << std::endl;
        cycle++;
    }
    allBestGenes = populations[allBestPopulationIndex][allBestIndex];
}
