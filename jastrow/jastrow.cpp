#include "jastrow.h"

#include "../config.h"

Jastrow::Jastrow(Config *config) :
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions())
{

    jastrowGradient = zeros<vec>(nParticles * nDimensions);
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
    delete [] rOld;
    delete [] rNew;
}

/*!
 * \brief Jastrow::initialize calculates the distances and jastrow arguments for an initial position
 * \param positions
 */
void Jastrow::initialize(vec2 positions[]) {
    calculateDistances(positions);
    distancesOld = distancesNew;
    jastrowArgumentsOld = jastrowArgumentsNew;
}

/*!
 * \brief Jastrow::calculateDistances sets the distance and argument matrix based on the given positions
 * \param r
 */
void Jastrow::calculateDistances(vec2 r[]) {
    vec2 diff;
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = r[i];
        for(int j = i + 1; j < nParticles; j++) {
            diff = (r[i] - r[j]);
            distancesNew.at(i,j) = norm(diff,2);
            jastrowArgumentsNew.at(i,j) = argument(i,j,distancesNew);
        }
    }
    distancesNew = symmatu(distancesNew); // reflect the upper triangle to the lower
    jastrowArgumentsNew = symmatu(jastrowArgumentsNew); // reflect the upper triangle to the lower
}

/*!
 * \brief Jastrow::acceptMove sets the old distance and argument matrices to their new values
 * \param movedParticle
 */
void Jastrow::acceptMove(int movedParticle)
{
    (void)movedParticle;
    for(int i = 0; i < nParticles; i++) {
        rOld[i] = rNew[i];
    }
    distancesOld = distancesNew;
    jastrowArgumentsOld = jastrowArgumentsNew;
}

/*!
 * \brief Jastrow::rejectMove sets the new distance and argument matrices to their old values
 */
void Jastrow::rejectMove()
{
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = rOld[i];
    }
    distancesNew = distancesOld;
    jastrowArgumentsNew = jastrowArgumentsOld;
}

double Jastrow::argument(int i, int j, mat &distances) {
    return (a(i,j) * distances(i,j)) / (1 + beta * distances(i,j));
}

/*!
 * Jastrow::laplacePartial calculates the Jastrow Laplace ratio, except for the gradient dot term
  * \warning Named "partial" because it does not include the square gradient term. (8.46 in Leirvåg's thesis) This must be included manually. This is done by WaveSlater in this code.
  */
double Jastrow::laplacePartial(vec2 r[]) {
    // Optimized by removing square gradient term. Reduced inclusive cost by 80 % !
    double laplaceSum = 0;
    for(int movedParticle = 0; movedParticle < nParticles; movedParticle++) {
        int p = movedParticle;
        for(int i = 0; i < nParticles; i++) {
            if(i != p) {
                double rpi = norm(r[p] - r[i],2);
//                double rpi = distancesOld.at(i,p);
                double numerator = (1 - beta * rpi);
                double denominator = (1 + beta * rpi);
                laplaceSum += a(i,p) * numerator / ( rpi * denominator * denominator * denominator);
            }
        }
    }

    return laplaceSum;
}

/*!
 * \brief Jastrow::gradient calculates the gradient analytically
 * \param r
 * \param movedParticlea
 * \param rGradient
 * \note Call ratio() or calculateDistances() before calling this function to update the distance matrix.
 */
void Jastrow::gradient(vec2 r[], vec &rGradient)
{
    rGradient.zeros();
    for(int p = 0; p < nParticles; p++) {
        for(int i = 0; i < nParticles; i++) {
            if(i != p) {
                // Optimized by writing as array operations instead of vector operations (18 % run time reduction)
                for(int j = 0; j < nDimensions; j++) {
                    rpiVec[j] = r[p][j] - r[i][j];
                }
                double rpi = norm(rpiVec,nDimensions);
                for(int j = 0; j < nDimensions; j++) {
                    double denominator = (1 + beta * rpi);
                    int gradientIndex = p*nDimensions + j;
                    rGradient(gradientIndex) += a(i,p) * rpiVec[j] / ( rpi * denominator * denominator);
                }
            }
        }
    }
}

/*!
 * \brief Jastrow::evaluate calculates the Jastrow part of the wave function
 * \param r
 * \return
 */
double Jastrow::evaluate(vec2 r[]) {
    double wf = 1;

    calculateDistances(r);
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            wf *= exp(jastrowArgumentsNew.at(i,j));
        }
    }
    return wf;
}

/*!
 * Jastrow::ratio calculates the Jastrow ratio
  * \note The distance and Jastrow argument matrices only have values in their upper right triangle.
  */
double Jastrow::ratio(vec2 &r, int movedParticle)
{
    rNew[movedParticle] = r;
//    // we only need to update the elements in the matrix that are affected by the move of one particle

    double argumentChange = 0;
    vec2 diff;
    for(int i = 0; i < nParticles; i++) {
        if(i != movedParticle) {
            int row = -1;
            int col = -1;
            if(i < movedParticle) {
                row = i;
                col = movedParticle;
            } else {
                row = movedParticle;
                col = i;
            }
            diff = rNew[i] - r;
            distancesNew.at(row,col) = sqrt(dot(diff,diff));
            jastrowArgumentsNew.at(row,col) = argument(row,col,distancesNew);
            argumentChange += jastrowArgumentsNew.at(row,col) - jastrowArgumentsOld.at(row,col);
        }
    }
    distancesNew = symmatu(distancesNew);
    jastrowArgumentsNew = symmatu(jastrowArgumentsNew);

    return exp(argumentChange);
//    return evaluate(rNew) / evaluate(rOld);
}

void Jastrow::setParameters(double *parameters)
{
    beta = parameters[1];
}
