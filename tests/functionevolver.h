#ifndef FUNCTIONEVOLVER_H
#define FUNCTIONEVOLVER_H

#include <../evolver/evolver.h>

class FunctionEvolver : public Evolver
{
public:
    FunctionEvolver(int nGenes, int nIndividuals, int nPopulations);

    vec fitnessResult;
    vec x;
    vec result;
    void calculate(vec &coefficients);
private:
    double fitness(vec &coefficients);


    vec myFunction();

};

#endif // FUNCTIONEVOLVER_H
