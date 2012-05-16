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
    x = linspace<vec>(0, 10, 100);
    result = sin(x);
    fitnessResult = zeros<vec>(result.n_elem);
}

void FunctionEvolver::calculate(vec &coefficients) {
    fitnessResult.zeros();
    for(uint i = 0; i < coefficients.n_elem; i = i + 4) {
        for(uint j = 0; j < x.n_elem; j++) {
            //fitnessResult[j] = fitnessResult[j] + coefficients(i) * x[j] + coefficients(i + 1) * x[j] * x[j];
            fitnessResult[j] += coefficients(i) * sin( coefficients(i + 1) * x[j] );
            fitnessResult[j] += coefficients(i+2) * cos( coefficients(i + 3) * x[j] );
        }
    }
}

double FunctionEvolver::fitness(vec &coefficients)
{
    calculate(coefficients);
    double diffSum = 0;
    for(uint i = 0; i < result.n_elem; i++) {
        diffSum += fabs(fitnessResult[i] - result[i]);
    }

    return diffSum;
}
