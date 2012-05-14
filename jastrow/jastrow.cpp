#include "jastrow.h"

#include "../config.h"

Jastrow::Jastrow(Config *config) :
    nParticles(config->nParticles()),
    nDimensions(config->nDimensions())
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

void Jastrow::initialize(vec2 positions[]) {
    calculateDistances(positions);
    distancesOld = distancesNew;
    jastrowArgumentsOld = jastrowArgumentsNew;
}

/*!
  */
void Jastrow::calculateDistances(vec2 r[]) {
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = r[i];
        for(int j = i + 1; j < nParticles; j++) {
            vec2 diff = r[i] - r[j];
            distancesNew.at(i,j) = sqrt(dot(diff,diff));
            jastrowArgumentsNew.at(i,j) = argument(i,j,distancesNew);
        }
    }
    distancesNew = symmatu(distancesNew); // reflect the upper triangle to the lower
    jastrowArgumentsNew = symmatu(jastrowArgumentsNew); // reflect the upper triangle to the lower
}

void Jastrow::acceptEvaluation(int movedParticle)
{
    (void)movedParticle;
    for(int i = 0; i < nParticles; i++) {
        rOld[i] = rNew[i];
    }
    distancesOld = distancesNew;
    jastrowArgumentsOld = jastrowArgumentsNew;
}

void Jastrow::refuseEvaluation()
{
    for(int i = 0; i < nParticles; i++) {
        rNew[i] = rOld[i];
    }
    distancesNew = distancesOld;
    jastrowArgumentsNew = jastrowArgumentsOld;
}

double Jastrow::argument(int i, int j, mat &distances) {
    return (a(i,j) * distances.at(i,j)) / (1 + beta * distances.at(i,j));
}

double Jastrow::laplace(vec2 r[], int movedParticlea) {
    double laplaceSum = 0;

    vec jastrowGradient = zeros<vec>(nParticles * nDimensions);
    gradient(r, movedParticlea, jastrowGradient);

    laplaceSum += dot(jastrowGradient, jastrowGradient);
    for(int movedParticle = 0; movedParticle < nParticles; movedParticle++) {
        int p = movedParticle;
        for(int i = 0; i < nParticles; i++) {
            if(i != p) {
//                double rpi = norm(r[p] - r[i],2);
                double rpi = distancesNew.at(i,p);
                double denominator = (1 + beta * rpi);
                double numerator = (1 - beta * rpi);
                laplaceSum += a(p,i) * numerator / ( rpi * denominator * denominator * denominator);
            }
        }
    }

    return laplaceSum;
}

/*!
 * \brief Jastrow::gradient
 * \param r
 * \param movedParticlea
 * \param rGradient
 * Note: Call ratio() or calculateDistances() before calling this function to update the distance matrix.
 */
void Jastrow::gradient(vec2 r[], int movedParticlea, vec &rGradient)
{
//    for(int i = 0; i < nParticles; i++) {
//        if(r[i][0] != rNew[i][0] || r[i][1] != rNew[i][1]) {
//            std::cout << "Inequality" << i << "\n" << r[i] << std::endl << rNew[i] << std::endl << rOld[i] << std::endl;
//        }
//    }
    rGradient.zeros();
    for(int movedParticle = 0; movedParticle < nParticles; movedParticle++) {
        int p = movedParticle;
        for(int i = 0; i < nParticles; i++) {
            if(i != p) {
                vec2 rpiVec = r[movedParticle] - r[i];
//                if(norm(rpiVec,2) - distancesNew.at(i,p) > 1e-8) {
//                    std::cout << i << " " << p << " IP " << std::endl;
//                    std::cout << norm(rpiVec,2) << " " << distancesNew.at(i,p) << " " << distancesNew.at(p,i) << " " << distancesOld.at(i,p) << " " << distancesOld.at(p,i) << std::endl;
//                }
                double rpi = distancesNew.at(i,p);
                for(int j = 0; j < nDimensions; j++) {
                    double denominator = (1 + beta * distancesNew.at(i,p));
                    int gradientIndex = movedParticle*nDimensions + j;
                    rGradient(gradientIndex) += a(p,i) * rpiVec[j] / ( rpi * denominator * denominator);
                }
            }
        }
    }
}

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
  Note: The distance and Jastrow argument matrices only have values in their upper right triangle.
  */
double Jastrow::ratio(vec2 &r, int movedParticle)
{
//    std::cout << "Moved particle " << movedParticle << std::endl;
    rNew[movedParticle] = r;
//    for(int i = 0; i < nParticles; i++) {
//        std::cout << "rNew[" << i << "] = " << rNew[i] << std::endl;
//    }
    // we only need to update the elements in the matrix that are affected by the move of one particle
    double argumentChange = 0;
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
            vec2 diff = rNew[i] - r;
            distancesNew.at(row,col) = norm(diff,2);
            jastrowArgumentsNew.at(row,col) = argument(row,col,distancesNew);
            argumentChange += jastrowArgumentsNew.at(row,col) - jastrowArgumentsOld.at(row,col);
        }
    }
//    jastrowArgumentsNew = jastrowArgumentsOld;
//    calculateDistances(rNew);
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = i + 1; j < nParticles; j++) {
//            argumentChange += jastrowArgumentsOld.at(i,j) - jastrowArgumentsNew.at(i,j);
//        }
//    }
    distancesNew = symmatu(distancesNew);
    jastrowArgumentsNew = symmatu(jastrowArgumentsNew);

    return exp(argumentChange);
}

void Jastrow::setParameters(double *parameters)
{
    beta = parameters[1];
}
