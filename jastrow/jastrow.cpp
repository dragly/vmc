#include "jastrow.h"

#include "../config.h"

Jastrow::Jastrow(Config *config) :
    nParticles(config->nParticles())
{
    a = zeros<mat>(nParticles, nParticles);
    distancesOld = zeros<mat>(nParticles, nParticles);
    distancesNew = zeros<mat>(nParticles, nParticles);
    jastrowArgumentsOld = zeros<mat>(nParticles, nParticles);
    jastrowArgumentsNew = zeros<mat>(nParticles, nParticles);
    rOld = new vec2[nParticles];
    rNew = new vec2[nParticles];


    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nParticles; j++) {
            if((i < nParticles / 2 && j < nParticles / 2) || (i >= nParticles / 2 && j >= nParticles / 2)) {
                a(i,j) = 1. / 3;
            } else {
                a(i,j) = 1.;
            }
        }
    }
}

Jastrow::~Jastrow() {
    delete rOld;
    delete rNew;
}

/*!
  Note: Only sets values for the upper right triangle of the distance and Jastrow argument matrices.
  */
void Jastrow::calculateDistances(vec2 r[]) {
    for(int i = 0; i < nParticles; i++) {
        rOld[i] = r[i];
        rNew[i] = r[i];
        for(int j = i + 1; j < nParticles; j++) {
            vec2 diff = r[i] - r[j];
            distancesOld.at(i,j) = sqrt(dot(diff,diff));
            jastrowArgumentsOld.at(i,j) = argument(i,j,distancesOld);
        }
    }
}

void Jastrow::acceptEvaluation(int movedParticle)
{
    (void)movedParticle;
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = rOld[i];
    }
}

double Jastrow::argument(int i, int j, mat &distances) {
    return (a(i,j) * distances.at(i,j)) / (1 + m_beta * distances.at(i,j));
}

void Jastrow::gradient(const vec2 &r, vec2 &rGradient)
{

}

double Jastrow::evaluate(vec2 r[]) {
    double wf = 1;

    calculateDistances(r);
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            wf *= exp(jastrowArgumentsOld.at(i,j));
        }
    }
    return wf;
}

/*!
  Note: The distance and Jastrow argument matrices only have values in their upper right triangle.
  */
double Jastrow::ratio(vec2 &r, int particleNumber)
{
    rNew[particleNumber] = r;
    // we only need to update the elements in the matrix that are affected by the move of one particle
    distancesNew = distancesOld;
    double argumentChange = 0;
    for(int i = 0; i < particleNumber; i++) {
        vec2 diff = rOld[i] - r;
        distancesNew.at(i,particleNumber) = sqrt(dot(diff,diff));
        jastrowArgumentsNew.at(i,particleNumber) = argument(i,particleNumber,distancesNew);
        argumentChange += jastrowArgumentsNew.at(i,particleNumber) - jastrowArgumentsOld.at(i,particleNumber);
    }
    for(int j = particleNumber + 1; j < nParticles; j++) {
        vec2 diff = r - rOld[j];
        distancesNew.at(particleNumber,j) = sqrt(dot(diff,diff));
        jastrowArgumentsNew.at(particleNumber,j) = argument(particleNumber,j,distancesNew);
        argumentChange += jastrowArgumentsNew.at(particleNumber,j) - jastrowArgumentsOld.at(particleNumber,j);
    }
//    jastrowArgumentsNew = jastrowArgumentsOld;
//    calculateDistances(rNew);
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = i + 1; j < nParticles; j++) {
//            argumentChange += jastrowArgumentsOld.at(i,j) - jastrowArgumentsNew.at(i,j);
//        }
//    }
    return exp(argumentChange);
}

void Jastrow::setParameters(double *parameters)
{
    m_beta = parameters[1];
}
