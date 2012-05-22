#include "functionevolver.h"

/*!
 * \brief FunctionEvolver::FunctionEvolver
 * \param nGenes
 * \param nIndividuals
 * \param nPopulations
 */
FunctionEvolver::FunctionEvolver(int nGenes, int nIndividuals, int nPopulations)
    : Evolver(nGenes, nIndividuals, nPopulations)
{
    x = linspace<vec>(-1, 1, 100);
    result = zeros<vec>(x.n_elem);
    for(uint j = 0; j < x.n_elem; j++) {
        result[j] = exp(- 2 * x[j] * x[j]) * 3 * exp(- 4 * x[j] * x[j]) * 3 * exp(- 4 * x[j] * x[j]) + 3 * sin(2 * x[j]);
    }
    fitnessResult = zeros<vec>(x.n_elem);
}

void FunctionEvolver::calculate(vec &coefficients) {
//    fitnessResult.zeros();
    fitnessResult.zeros();
    for(uint i = 0; i < coefficients.n_elem; i = i + 4) {
        for(uint j = 0; j < x.n_elem; j++) {
            //fitnessResult[j] = fitnessResult[j] + coefficients(i) * x[j] + coefficients(i + 1) * x[j] * x[j];
            fitnessResult[j] += coefficients(i) * exp( - fabs(coefficients(i + 1)) * x[j]  * x[j] );
            fitnessResult[j] += coefficients(i+2) * sin(coefficients(i + 3)   * x[j] );
        }
    }
}

double FunctionEvolver::fitness(vec &coefficients, int population, int individual)
{
    (void)population;
    (void)individual;
    calculate(coefficients);
    double diffSum = 0;
    for(uint i = 0; i < result.n_elem; i++) {
        diffSum += fabs(fitnessResult[i] - result[i]);
    }

    return diffSum;
}
