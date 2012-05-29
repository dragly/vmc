#include "../wavefunction/wavesimple.h"
#include "../wavefunction/waveideal.h"
#include "../matrix.h"
#include "../config.h"
#include "../hamiltonian/hamiltonianideal.h"
#include "../montecarlo/standardmontecarlo.h"
#include "../montecarlo/metropolishastingsmontecarlo.h"
#include "../minimizer/geneticminimizer.h"
#include "../hermite.h"
#include "../orbital/orbital.h"
#include "../slater/slater.h"
#include "../wavefunction/waveslater.h"
#include "../jastrow/jastrow.h"
#include "../random.h"
#include "../montecarlo/diffusionmontecarlo.h"
#include "../montecarlo/evolutionarymontecarlo.h"
#include "functionevolver.h"

#include <QtCore/QString>
#include <QtTest/QtTest>

#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

class VmcTests : public QObject
{
    Q_OBJECT
    
public:
    VmcTests();

    // old tests
    void waveSimpleGradientTest(); // TODO - consider implementing this again

    // quick tests
    void orbitalGradientTest();
    void initTestCase();
    void orbitalTest();
    void cleanupTestCase();
    void waveSimpleLaplaceTest();
    void waveIdealLaplaceTest();
    void hermiteTest();
    void twoOrbitalsOneWavefunctionTest();
    void slaterFourParticleTest();
    void waveSlaterSixParticleTest();
    void slaterRatioTest();
    void orbitalLaplaceTest();
    void jastrowTest();
    void jastrowRatioTest();
    void slaterTest();
    void slaterInverse();
    void waveSlaterLaplaceTest();
    void waveSlaterGradientTest();
    // slow tests
    void fullIdealTest();
    void fullIdealHastingsTest();
    void fullIdealHastingsSlaterTest();
    void fullSlaterSixInteractionTest();
    void fullSlaterSixNoInteractionTest();

    // unfinished tests
    void evolverTest();
    void diffusionMonteCarloTest();
    void evolutionaryMonteCarloTest();
private slots:
    // quick tests
    // slow tests
    // unfinished tests
    void geneticMinimizerTest();

private:
    Config *oldConfig;
    WaveSimple *waveSimple;
    WaveIdeal *waveIdeal;
    HamiltonianIdeal *hamiltonianIdeal;
    vec2 *r_old;
    double charge;
};

VmcTests::VmcTests()
{
}

void VmcTests::initTestCase()
{
    oldConfig = new Config(0,1);
    // Set up waveSimple
    int nParticles = oldConfig->nParticles();
    cout << nParticles << endl;
    int nDimensions = oldConfig->nDimensions();
    charge = 1.0;
    r_old = new vec2[nParticles];
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            r_old[i][j] = 0.234 + i + 2*j;
        }
    }
    waveSimple = new WaveSimple(oldConfig);
    double parameters[2];
    parameters[0] = 2;
    parameters[1] = 1;
    waveSimple->setParameters(parameters);
    // Set up waveIdeal
    waveIdeal = new WaveIdeal(oldConfig);
    waveIdeal->setParameters(parameters);
    // Set up hamiltonianIdeal
    hamiltonianIdeal = new HamiltonianIdeal(oldConfig);
}

void VmcTests::cleanupTestCase()
{
}

void VmcTests::waveSimpleLaplaceTest()
{
    waveSimple->setUseAnalyticalLaplace(true);
    double analyticalLaplace = waveSimple->laplace(r_old);
    waveSimple->setUseAnalyticalLaplace(false);
    double numericalLaplace = waveSimple->laplaceNumerical(r_old);
    QVERIFY(fabs(analyticalLaplace - numericalLaplace) < 0.001);
}

void VmcTests::waveIdealLaplaceTest()
{
    waveIdeal->setUseAnalyticalLaplace(true);
    double analyticalLaplace = waveIdeal->laplace(r_old);
    waveIdeal->setUseAnalyticalLaplace(false);
    double numericalLaplace = waveIdeal->laplaceNumerical(r_old);
    QVERIFY(fabs(analyticalLaplace - numericalLaplace) < 0.001);
}

void VmcTests::waveSimpleGradientTest()
{
    cout << "Particles..." << endl;
    int nParticles = 1;
    WaveSimple *waveSimpleNew = new WaveSimple(oldConfig);
    vec2* rPositions = new vec2[nParticles];
    rPositions[0][0] = -1;
    rPositions[0][1] = 0.0;
    vec2 rGradient;
    double parameters[2];
    parameters[0] = 2;
    parameters[1] = 1;
    waveSimpleNew->setParameters(parameters);
    waveSimpleNew->gradient(rPositions, rGradient);
    QVERIFY(fabs(rGradient[0] - 0.735) < 0.001);
    QVERIFY(fabs(rGradient[1] - 0.000) < 0.0000001);
    cout << "Particles..." << endl;
}

void VmcTests::fullIdealTest()
{
    int nCycles = 500000;
    Config *config1 = new Config(1,1);
    config1->setNDimensions(2);
    config1->setNParticles(2);
    WaveIdeal* waveIdeal1 = new WaveIdeal(config1);
    config1->setWave(waveIdeal1);
    HamiltonianIdeal* hamiltonianIdeal1 = new HamiltonianIdeal(config1);
    config1->setHamiltonian(hamiltonianIdeal1);
    waveIdeal1->setUseAnalyticalLaplace(true);
    double parameters[2];
    parameters[0] = 1;
    parameters[1] = 0.4;
    waveIdeal1->setParameters(parameters);
    StandardMonteCarlo *monteCarlo = new StandardMonteCarlo(config1);
    double energy;
    //  Do the mc sampling
    monteCarlo->sample(nCycles);
    energy = monteCarlo->energy();
    //    std::cout << "Full ideal energy was " << energy << std::endl;
    QVERIFY(fabs(energy - 3.00034530284643397025) < 1e-2);
}

void VmcTests::slaterInverse() {
    Config *config1 = new Config(1,1);
    config1->setNDimensions(2);
    config1->setNParticles(4);
    double parameters[2];
    long idum = -1;
    Orbital **orbitals = new Orbital*[4];
    for(int i = 0; i < 4; i++) {
        orbitals[i] = new Orbital(0,i,config1);
    }
    vec2 r[4];
    for(int i = 0; i < config1->nParticles(); i++) {
        for(int j = 0; j < config1->nDimensions(); j++) {
            r[i].at(j) = ran2(&idum);
        }
    }

    for(int a = 1; a < 10; a++) {
        parameters[0] = a * 0.31;
        parameters[1] = a * 0.92;
        for(int i = 0; i < 4; i++) {
            orbitals[i]->setParameters(parameters);
        }
        Slater* slater1 = new Slater(config1, orbitals, true);
        Slater* slater2 = new Slater(config1, orbitals, true);
        slater1->initialize(r);
        slater2->initialize(r);
        //        std::cout << "Init:" << std::endl;
        //        std::cout << slater1->matrix() << std::endl;
        //        std::cout << slater2->matrix() << std::endl;
        for(int cycle = 0; cycle < 20; cycle++) {
            for(int movedParticle = 0; movedParticle < config1->nParticles(); movedParticle++) {
                for(int j = 0; j < config1->nDimensions(); j++) {
                    r[movedParticle].at(j) = ran2(&idum);
                }
                //                std::cout << "movedParticle:" << movedParticle << std::endl;
                slater1->ratio(r[movedParticle] , movedParticle);
                slater1->calculateInverse(movedParticle);
                mat analyticalInverse = slater1->inverse();
                slater2->initialize(r);
                mat numericalInverse = slater2->inverse();
                //                std::cout << "Matrix: " << std::endl;
                //                std::cout << slater1->matrix() << std::endl << slater2->matrix() << std::endl;
                //                std::cout << "Inverse:" << std::endl;
                //                std::cout << analyticalInverse << std::endl << numericalInverse << std::endl;
                for(int i = 0; i < 2; i++) {
                    for(int j = 0; j < 2; j++) {
                        QCOMPARE(analyticalInverse.at(i,j), numericalInverse.at(i,j));
                    }
                }
                slater1->acceptMove(movedParticle);
            }
        }
    }
}

void VmcTests::slaterRatioTest() {
    Config *config1 = new Config(1,1);
    config1->setNDimensions(2);
    config1->setNParticles(4);
    Orbital **orbitals = new Orbital*[4];
    for(int i = 0; i < 4; i++) {
        orbitals[i] = new Orbital(0,i,config1);
    }
    vec2 rOld[4];
    vec2 rNew[4];
    int movedParticle = 1;
    for(int i = 0; i < config1->nParticles(); i++) {
        for(int j = 0; j < config1->nDimensions(); j++) {
            rOld[i].at(j) = 0.2 * i * config1->nDimensions() + j * 0.1 + 0.24;
        }
        rNew[i] = rOld[i];
    }
    rNew[movedParticle].at(0) = 1;
    rNew[movedParticle].at(1) = 2;
    Slater* slater1 = new Slater(config1, orbitals, true);
    Slater* slater2 = new Slater(config1, orbitals, true);
    for(int a = 1; a < 10; a++) {
        double parameters[2];
        parameters[0] = a * 0.19;
        parameters[1] = a * 0.23;
        for(int i = 0; i < 4; i++) {
            orbitals[i]->setParameters(parameters);
        }
        double simpleRatio = slater1->determinant(rNew) / slater1->determinant(rOld);
        slater2->initialize(rOld);
        //    slater2->setPreviousMovedParticle(1);
        double fancyRatio = slater2->ratio(rNew[movedParticle], movedParticle);
        QCOMPARE(simpleRatio, fancyRatio);
    }
    delete slater2;
    delete slater1;
}

void VmcTests::jastrowRatioTest() {
    Config *config1 = new Config(1,1);
    config1->setNDimensions(2);
    config1->setNParticles(4);

    vec2 rOld[4];
    vec2 rNew[4];
    int movedParticle = 1;
    for(int i = 0; i < config1->nParticles(); i++) {
        for(int j = 0; j < config1->nDimensions(); j++) {
            rOld[i].at(j) = i * config1->nDimensions() + j;
        }
        rNew[i] = rOld[i];
    }
    rNew[movedParticle].at(0) = 1;
    rNew[movedParticle].at(1) = 2;

    Jastrow* jastrow1 = new Jastrow(config1);
    double simpleRatio = jastrow1->evaluate(rNew) / jastrow1->evaluate(rOld);

    Jastrow* jastrow2 = new Jastrow(config1);
    jastrow2->initialize(rOld);
    double fancyRatio = jastrow2->ratio(rNew[movedParticle], movedParticle);
    QCOMPARE(simpleRatio, fancyRatio);
}

void VmcTests::fullIdealHastingsTest()
{
    int nCycles = 1000000;
    Config *config1 = new Config(1,1);
    config1->setNDimensions(2);
    config1->setNParticles(2);
    WaveIdeal* waveIdeal1 = new WaveIdeal(config1);
    config1->setWave(waveIdeal1);
    HamiltonianIdeal* hamiltonianIdeal1 = new HamiltonianIdeal(config1);
    config1->setHamiltonian(hamiltonianIdeal1);
    waveIdeal1->setUseAnalyticalLaplace(true);
    double parameters[2];
    parameters[0] = 1;
    parameters[1] = 0.4;
    waveIdeal1->setParameters(parameters);
    MetropolisHastingsMonteCarlo *monteCarlo = new MetropolisHastingsMonteCarlo(config1);
    monteCarlo->setThermalizationEnabled(false);
    double energy;
    //  Do the mc sampling
    monteCarlo->sample(nCycles);
    energy = monteCarlo->energy();
        std::cout << "Full ideal Hastings energy was " << energy << std::endl;
    QVERIFY(fabs(energy - 3.000) < 1e-2);
}

void VmcTests::hermiteTest() {
    //    Hermite* hermite5 = new Hermite(5);

    QVERIFY(Hermite::evaluate(2,4) - (4*4*4 - 2) < 1e-20);
    QVERIFY(Hermite::evaluate(3,5) - (8*5*5*5 - 12) < 1e-20);
    QVERIFY(Hermite::evaluate(4,3) - (16*3*3*3*3 - 48*3*3 + 12) < 1e-20);
}

/*!
  Testing that two orbitals multiplied returns the same as the simple wave function class.
  */
void VmcTests::twoOrbitalsOneWavefunctionTest() {
    // two particles
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(2);
    vec2 *rpos = new vec2[config1->nParticles()];
    rpos[0][0] = 0.4;
    rpos[0][1] = 0.9;
    rpos[1][0] = 0.1;
    rpos[1][1] = 0.8;

    WaveSimple *waveSimple1 = new WaveSimple(config1);
    double parameters[2];
    parameters[0] = 2.0;
    parameters[1] = 1.0;
    waveSimple1->setParameters(parameters);

    Orbital *orbital1 = new Orbital(0,0,config1);
    Orbital *orbital2 = new Orbital(0,0,config1);
    orbital1->setParameters(parameters);
    orbital2->setParameters(parameters);

    double val1 = waveSimple1->evaluate(rpos);

    double val2 = orbital1->evaluate(rpos[0]) * orbital2->evaluate(rpos[1]);

    QCOMPARE(val1,val2);
}

/*!
  Testing multiple aspects of the slater class.
  */
void VmcTests::slaterTest() {
    // two particles
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(2);
    vec2 *rpos = new vec2[config1->nParticles()];
    rpos[0][0] = 0.4;
    rpos[0][1] = 0.9;
    rpos[1][0] = 0.1;
    rpos[1][1] = 0.8;


    WaveSimple *waveSimple1 = new WaveSimple(config1);
    double parameters[2];
    parameters[0] = 2.0;
    parameters[1] = 1.0;
    waveSimple1->setParameters(parameters);

    double val1 = waveSimple1->evaluate(rpos);

    Orbital *orbitals[1];
    orbitals[0] = new Orbital(0,0,config1);
    orbitals[0]->setParameters(parameters);

    Slater *slaterUp = new Slater(config1, orbitals, true);
    Slater *slaterDown = new Slater(config1, orbitals, false);

    double normFactorial = 1;
    for(int i = 2; i < config1->nParticles() / 2; i++) {
        normFactorial *= i;
    }

    double val2 =  1 / sqrt(normFactorial) * slaterUp->determinant(rpos) * slaterDown->determinant(rpos);

    QCOMPARE(val1,val2);
}

/*!
  Test the Metropolis-Hastings algorithm by using the Slater determinant in the wave function
  with two particles. This is compared to the energy found in Lars-Eivind Lervåg's thesis.
  */
void VmcTests::fullIdealHastingsSlaterTest()
{
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(2);
    WaveSlater *waveSlater1 = new WaveSlater(config1);
    HamiltonianIdeal *hamiltonian1 = new HamiltonianIdeal(config1);
    config1->setHamiltonian(hamiltonian1);
    config1->setWave(waveSlater1);
    double parameters[2];
    parameters[0] = 1.0;
    parameters[1] = 0.4;
    waveSlater1->setParameters(parameters);
    int nCycles = 100000;
    //    waveSlater1->setUseAnalyticalLaplace(false);
    MetropolisHastingsMonteCarlo *monteCarlo = new MetropolisHastingsMonteCarlo(config1);
    //    double *allEnergies = new double[nCycles+1];
    double energy;
    //  Do the mc sampling
    monteCarlo->sample(nCycles);
    energy = monteCarlo->energy();
    cout << "Two particle energy is " << fabs(energy) << endl;
    QVERIFY(fabs(energy - 3.000) < 1e-2);
}

/*!
  Tests whether the orbitals return their analytical equivalents
*/
void VmcTests::orbitalTest() {
    Config *config1 = new Config(0,1);
    double omega = 1;
    config1->setOmega(omega);
    Orbital *orbital11 = new Orbital(1,1,config1);
    Orbital *orbital23 = new Orbital(2,3,config1);
    double parameters[2];
    parameters[0] = 2.0;
    parameters[1] = 3.0;
    orbital11->setParameters(parameters);
    orbital23->setParameters(parameters);
    vec2 r;
    r[0] = 1.8;
    r[1] = 2.9;
    double A = 1;
    QCOMPARE(orbital11->evaluate(r), A * 2 * sqrt(omega * parameters[0]) * r[0] * 2 * sqrt(omega * parameters[0]) * r[1] * exp(- parameters[0] * omega * (r[0] * r[0] + r[1] * r[1]) / 2));
    double sqrtomegax = (sqrt(omega * parameters[0]) * r[0]);
    double sqrtomegay = (sqrt(omega * parameters[0]) * r[1]);
    QCOMPARE(orbital23->evaluate(r), A * (4 * sqrtomegax * sqrtomegax - 2) * (8 * sqrtomegay * sqrtomegay * sqrtomegay - 12) * exp(- parameters[0] * omega * (r[0] * r[0] + r[1] * r[1]) / 2));
}

void VmcTests::jastrowTest() {
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(4);
    double omega = 1;
    config1->setOmega(omega);
    config1->setOmega(omega);
    double parameters[2];
    parameters[0] = 2.0;
    parameters[1] = 3.0;
    vec2 *r = new vec2[config1->nParticles()];
    r[0][0] = 0.4;
    r[0][1] = 0.9;
    r[1][0] = 0.2;
    r[1][1] = 0.3;
    r[2][0] = 0.4;
    r[2][1] = 0.5;
    r[3][0] = 0.6;
    r[3][1] = 0.7;
    Jastrow *jastrow = new Jastrow(config1);
    jastrow->setParameters(parameters);
    double r01 = sqrt(dot(r[0] - r[1], r[0] - r[1]));
    double r02 = sqrt(dot(r[0] - r[2], r[0] - r[2]));
    double r03 = sqrt(dot(r[0] - r[3], r[0] - r[3]));
    double r12 = sqrt(dot(r[1] - r[2], r[1] - r[2]));
    double r13 = sqrt(dot(r[1] - r[3], r[1] - r[3]));
    double r23 = sqrt(dot(r[2] - r[3], r[2] - r[3]));
    double jastrowValue = exp(1./3. * r01 / (1 + parameters[1] * r01))
            * exp(1. * r02 / (1 + parameters[1] * r02))
            * exp(1. * r03 / (1 + parameters[1] * r03))
            * exp(1. * r12 / (1 + parameters[1] * r12))
            * exp(1. * r13 / (1 + parameters[1] * r13))
            * exp(1./3. * r23 / (1 + parameters[1] * r23));

    QCOMPARE(jastrow->evaluate(r), jastrowValue);
}

/*!
  Tests the slater gradient numerically against our implementation
  */
void VmcTests::waveSlaterGradientTest()
{
    int nParticles = 6;
    int nDimensions = 2;
    Config *config1 = new Config(0,1);
    config1->setNDimensions(nDimensions);
    config1->setNParticles(nParticles);
    double parameters[2];
    long idum = -1;
    for(int w = 1; w < 11; w++) {
        config1->setOmega(w * 0.1);
        if(w > 5) {
            config1->setInteractionEnabled(false);
        } else {
            config1->setInteractionEnabled(true);
        }
        WaveSlater *waveSlater1 = new WaveSlater(config1);
        vec2 *r = new vec2[nParticles];
        for(int i = 0; i < nParticles; i++) {
            r[i][0] = 4 * ran2(&idum) - 2;
            r[i][1] = 4 * ran2(&idum) - 2;
            //            r[i][0] = 0.34 + i * 0.2;
            //            r[i][1] = 0.12 + i * 0.42;
        }
        for(int p = 1; p < 11; p++) {
            parameters[0] = ran2(&idum);
            parameters[1] = ran2(&idum);
            waveSlater1->setParameters(parameters);
            waveSlater1->initialize(r);
            vec analyticalGradient = zeros<vec>(nParticles * nDimensions);
            vec numericalGradient = zeros<vec>(nParticles * nDimensions);
            waveSlater1->gradient(r, analyticalGradient);
            waveSlater1->gradientNumerical(r, numericalGradient);
            for(int i = 0; i < nDimensions * nParticles; i++) {
                //                std::cout << analyticalGradient[i] << " " << numericalGradient[i] << std::endl;
                //                std::cout << "Diff: " << fabs(analyticalGradient[i] - numericalGradient[i]) << std::endl;
                QVERIFY(fabs(analyticalGradient[i] - numericalGradient[i]) < 1e-2);
            }
        }
        delete waveSlater1;
    }
}

/*!
  Tests the slater laplace numerically against our implementation
  */
void VmcTests::waveSlaterLaplaceTest()
{
    int nParticles = 4;
    int nDimensions = 2;
    Config *config1 = new Config(0,1);
    config1->setNDimensions(nDimensions);
    config1->setNParticles(nParticles);
    double parameters[2];
    long idum = -1;
    vec2 *r = new vec2[nParticles];
    for(int i = 0; i < nParticles; i++) {
        r[i][0] = 4 * ran2(&idum) - 2;
        r[i][1] = 4 * ran2(&idum) - 2;
    }

    for(int w = 1; w < 11; w++) {
        config1->setOmega(w * 0.1);
        if(w > 5) {
            //            std::cout << "Interaction enabled" << std::endl;
            config1->setInteractionEnabled(true);
        } else {
            //            std::cout << "Interaction disabled" << std::endl;
            config1->setInteractionEnabled(false);
        }
        WaveSlater *waveSlater1 = new WaveSlater(config1);
        for(int p = 1; p < 10; p++) {
            parameters[0] = ran2(&idum);
            parameters[1] = ran2(&idum);
            if(p > 5) {
                config1->setInteractionEnabled(true);
            } else {
                config1->setInteractionEnabled(false);
            }
            waveSlater1->setParameters(parameters);
            waveSlater1->initialize(r);
            double analyticalLaplace = waveSlater1->laplace(r);
            double numericalLaplace = waveSlater1->laplaceNumerical(r);
            //            std::cout << analyticalLaplace << " " << numericalLaplace << std::endl;
            //            std::cout << "Diff: " << fabs(analyticalLaplace - numericalLaplace) << std::endl;
            QVERIFY(fabs(analyticalLaplace - numericalLaplace) < 1e-5);
        }
    }
}

/*!
  Test the slater wave function with four particles. We compare it to a result
  from multiplying a hand-crafted determinant and the jastrow factor (tested
  separately).
  */
void VmcTests::slaterFourParticleTest()
{
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(4);
    double omega = 1;
    config1->setOmega(omega);
    Orbital *orbital00 = new Orbital(0,0,config1);
    Orbital *orbital01 = new Orbital(0,1,config1);
    double parameters[2];
    parameters[0] = 2.0;
    parameters[1] = 3.0;
    orbital00->setParameters(parameters);
    orbital01->setParameters(parameters);
    Jastrow *jastrow = new Jastrow(config1);
    jastrow->setParameters(parameters);
    WaveSlater *waveSlater1 = new WaveSlater(config1);
    config1->setWave(waveSlater1);
    vec2 *r = new vec2[config1->nParticles()];
    r[0][0] = 0.4;
    r[0][1] = 0.9;
    r[1][0] = 0.2;
    r[1][1] = 0.3;
    r[2][0] = 0.4;
    r[2][1] = 0.5;
    r[3][0] = 0.6;
    r[3][1] = 0.7;
    waveSlater1->setParameters(parameters);
    double waveValue = (orbital00->evaluate(r[0]) * orbital01->evaluate(r[1]) - orbital01->evaluate(r[0]) * orbital00->evaluate(r[1]))
            * (orbital00->evaluate(r[2]) * orbital01->evaluate(r[3]) - orbital01->evaluate(r[2]) * orbital00->evaluate(r[3]))
            * jastrow->evaluate(r);
    QCOMPARE(waveSlater1->evaluate(r), waveValue);
}
// TODO: Test slater wave for four particles

/*!
  Test the slater wave function with multiple particles. There is no benchmark
  for this function, so we just make sure it runs and prints without error.
  */
void VmcTests::waveSlaterSixParticleTest()
{
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(6);
    config1->setHamiltonian(hamiltonianIdeal);
    WaveSlater *waveSlater1 = new WaveSlater(config1);
    config1->setWave(waveSlater1);
    vec2 *rpos = new vec2[config1->nParticles()];
    rpos[0][0] = 0.4;
    rpos[0][1] = 0.9;
    rpos[1][0] = 0.2;
    rpos[1][1] = 0.3;
    rpos[2][0] = 0.4;
    rpos[2][1] = 0.5;
    rpos[3][0] = 0.6;
    rpos[3][1] = 0.7;
    rpos[4][0] = 0.8;
    rpos[4][1] = 0.9;
    rpos[5][0] = 0.12;
    rpos[5][1] = 0.14;
    double parameters[2];
    parameters[0] = 1.0;
    parameters[1] = 0.4;
    waveSlater1->setParameters(parameters);
    //    cout << "Value for slaterSixParticleTest: " << waveSlater1->evaluate(rpos) << endl;
}

void VmcTests::orbitalGradientTest()
{
    int nParticles = 2;
    int nDimensions = 2;
    double hstep = 0.00001;
    Config *config1 = new Config(0,1);
    config1->setNParticles(nParticles);
    config1->setNDimensions(nDimensions);
    double parameters[2];
    for(int nx = 0; nx < 5; nx++) {
        for(int ny = 0; ny < 5; ny++) {
            for(int p = 1; p < 11; p++) {
                parameters[0] = 0.34 * p;
                parameters[1] = 0.12 * p;
                config1->setOmega(0.23 * p);
                //                            std::cout << std::endl << "Orbital " << nx  << " " << ny << std::endl;
                Orbital *orbital = new Orbital(nx,ny,config1);
                orbital->setParameters(parameters);

                vec2 rPlus;
                vec2 rMinus;
                vec2 rPos;
                vec2 simpleGradient;
                for(int i = 1; i < 10; i++) {
                    rPos[0] = i*0.12;
                    rPos[1] = i*0.05;
                    rPlus = rPos;
                    rMinus = rPos;
                    //                std::cout << rPos << std::endl;
                    vec2 analyticalGradient;
                    orbital->gradient(rPos, analyticalGradient);
                    for (int j = 0; j < nDimensions; j++) {
                        rPlus = rPos;
                        rMinus = rPos;
                        rPlus[j] = rPos[j] + hstep;
                        rMinus[j] = rPos[j] - hstep;
                        double wfminus = orbital->evaluate(rMinus);
                        double wfplus  = orbital->evaluate(rPlus);
                        simpleGradient[j] = (wfplus - wfminus)/(2*hstep);
                        //                        std::cout << analyticalGradient[j] << " " << simpleGradient[j] << std::endl;
                        QVERIFY(fabs(analyticalGradient[j] - simpleGradient[j]) < 1e-4);
                    }
                }
                delete orbital;
            }
        }
    }
}

void VmcTests::orbitalLaplaceTest()
{
    int nParticles = 2;
    int nDimensions = 2;
    double hstep = 0.00001;
    Config *config1 = new Config(0,1);
    config1->setNParticles(nParticles);
    config1->setNDimensions(nDimensions);
    double parameters[2];
    for(int nx = 0; nx < 5; nx++) {
        for(int ny = 0; ny < 5; ny++) {
            for(int p = 1; p < 11; p++) {
                parameters[0] = 0.34 * p;
                parameters[1] = 0.12 * p;
                config1->setOmega(0.23 * p);
                //                std::cout << std::endl << "Orbital " << nx  << " " << ny << std::endl;
                Orbital *orbital = new Orbital(nx,ny,config1);
                orbital->setParameters(parameters);

                vec2 rPlus;
                vec2 rMinus;
                vec2 rPos;
                for(int i = 1; i < 10; i++) {

                    rPos[0] = i*0.12;
                    rPos[1] = i*0.05;
                    rPlus = rPos;
                    rMinus = rPos;
                    //                std::cout << rPos << std::endl;
                    double analyticalLaplace = orbital->laplace(rPos);
                    double simpleLaplace = 0;
                    for (int j = 0; j < nDimensions; j++) {
                        rPlus = rPos;
                        rMinus = rPos;
                        rPlus[j] = rPos[j] + hstep;
                        rMinus[j] = rPos[j] - hstep;
                        double wfold = orbital->evaluate(rPos);
                        double wfminus = orbital->evaluate(rMinus);
                        double wfplus  = orbital->evaluate(rPlus);
                        double partialLaplace =  (wfminus+wfplus-2*wfold) / (hstep * hstep);
                        //                    std::cout << partialLaplace << std::endl;
                        simpleLaplace += partialLaplace;
                    }
                    //                    std::cout << analyticalLaplace << " " << simpleLaplace << std::endl;
                    //                    std::cout << "diff: " << fabs(analyticalLaplace - simpleLaplace) << std::endl;
                    QVERIFY(fabs(analyticalLaplace - simpleLaplace) < 1e-2);
                }
                delete orbital;
            }
        }
    }
}

/*!
  Test the Metropolis-Hastings algorithm by using the Slater determinant with 6 particles.
  This is benchmarked against the optimal result found in Lars-Eivind Lervåg's thesis.

  This test is slow, but gives a decent benchmark.
  */
void VmcTests::fullSlaterSixNoInteractionTest()
{
    Config *config1 = new Config(0,1);
    config1->setNDimensions(2);
    config1->setNParticles(6);
    config1->setInteractionEnabled(false);
    config1->setStepLength(0.01);
    HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(config1);
    config1->setHamiltonian(hamiltonian);
    WaveSlater *waveSlater1 = new WaveSlater(config1);
    waveSlater1->setUseAnalyticalGradient(true);
    waveSlater1->setUseAnalyticalLaplace(true);
    config1->setWave(waveSlater1);
    double parameters[2];
    parameters[0] = 1.0;
    parameters[1] = 1.0;
    waveSlater1->setParameters(parameters);
    int nCycles = 50000;
    //    waveSlater1->setUseAnalyticalLaplace(false);
    MetropolisHastingsMonteCarlo *monteCarlo1 = new MetropolisHastingsMonteCarlo(config1);
    //  Do the mc sampling
    monteCarlo1->sample(nCycles);
    double energy = monteCarlo1->energy();
    cout << "Six non-interacting energy was " << fabs(energy) << endl;
    QVERIFY(fabs(energy - 10) < 1e-2);
}

/*!
  Test the Metropolis-Hastings algorithm by using the Slater determinant with 6 particles.
  This is benchmarked against the optimal result found in Lars-Eivind Lervåg's thesis.

  This test is slow, but gives a decent benchmark.
  */
void VmcTests::fullSlaterSixInteractionTest()
{
    QBENCHMARK {
        Config *config1 = new Config(0,1);
        config1->setNParticles(6);
        config1->setNDimensions(2);
        config1->setStepLength(0.001);
        HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(config1);
        config1->setHamiltonian(hamiltonian);
        config1->setInteractionEnabled(true);
        int nCycles = 200000;
        //        int nCycles = 500; // using few cycles for profiling
        WaveSlater *waveSlater2 = new WaveSlater(config1);
        waveSlater2->setUseAnalyticalGradient(true);
        waveSlater2->setUseAnalyticalLaplace(true);
        double parameters[2];
        parameters[0] = 0.9303;
        parameters[1] = 0.5493;
        waveSlater2->setParameters(parameters);
        config1->setWave(waveSlater2);
        MetropolisHastingsMonteCarlo *monteCarlo2 = new MetropolisHastingsMonteCarlo(config1);
        monteCarlo2->setThermalizationEnabled(true);
        monteCarlo2->sample(nCycles);
        double energy = monteCarlo2->energy();
        cout << "Six interacting energy was " << fabs(energy) << endl;
        QVERIFY(fabs(energy - 20.190) < 1e-1);
    }
    //    std::cout << "Benchmark used to be 51.328 seconds @ hyperon" << std::endl;
    //    std::cout << "Benchmark used to be 6.835 seconds @ home after first optimizations" << std::endl;
    //    std::cout << "Benchmark used to be 4.045 seconds @ home after profiling optimizations" << std::endl;
}

void VmcTests::evolverTest()
{
    //    QBENCHMARK {
    double sum = 0;
    double sumsq = 0;
    for(int test = 0; test < 100; test++) {
        std::cout << test << std::endl;
        FunctionEvolver *evolver = new FunctionEvolver(16, 16, 2);
        evolver->setRescaleLimits(1e-10, 1e10);
        evolver->evolve(500,250);
        //        std::cout << evolver->fitnessResult << std::endl;
        evolver->calculate(evolver->allBestGenes);

        //        ofstream data;
        //        data.open("fitness.dat");
        //        for(uint i = 0; i < evolver->x.n_elem; i++) {
        //            data << evolver->x[i] << " " << evolver->result[i] << " " << evolver->fitnessResult[i] << std::endl;
        //        }
        //        data.close();

        for(uint i = 0; i < evolver->x.n_elem; i++) {
            //            QVERIFY(fabs(evolver->result[i] - evolver->fitnessResult[i]) < 1e-1);
            double val = fabs(evolver->result[i] - evolver->fitnessResult[i]);
            sum += val;
            sumsq += val*val;
        }

        delete evolver;
    }
    std::cout << "Sum: " << sum << " variance: " << (sum/1000. * sum/1000.) - sumsq/1000 << std::endl;
    //    }
}

void VmcTests::geneticMinimizerTest() {
    QBENCHMARK {
        Config *config = new Config(0, 1);
        config->setNParticles(2);
        config->setNDimensions(2);
        config->setOmega(1);
        config->setStepLength(2);
        WaveFunction *wave = new WaveSlater(config);
        double parameters[2];
        parameters[0] = 999;
        parameters[1] = 999;
        wave->setParameters(parameters);
        wave->setUseAnalyticalGradient(true);
        wave->setUseAnalyticalLaplace(true);
        config->setWave(wave);
        config->setHamiltonian(new HamiltonianIdeal(config));
        config->setMonteCarlo(new StandardMonteCarlo(config));
        //        config->setMonteCarlo(new MonteCarloMetropolisHastings(config));
        GeneticMinimizer *evolver = new GeneticMinimizer(config);
        evolver->setNCycles(1000);
        evolver->setNSamples(1000, 20000);
        evolver->setPopulationData(2, 16, 2);
        evolver->setRescaleLimits(0.1, 2);
        evolver->setRescaleCycles(2);
        //        evolver->setGeneLimits(0.3, 1.2);
        evolver->runMinimizer();

        delete evolver;
    }
}

void VmcTests::diffusionMonteCarloTest() {
    Config *config = new Config(0, 1);
    config->setNParticles(6);
    config->setNDimensions(2);
    config->setInteractionEnabled(true);
    config->setDiffusionConstant(0.5);
    config->setStepLength(0.001);
    config->setOmega(1.0);

    double parameters[2];
//        parameters[0] = 0.987;
//        parameters[1] = 0.4;
    parameters[0] = 0.9303;
    parameters[1] = 0.5493;
    //    parameters[0] = 0.92;
    //    parameters[1] = 0.565;

    WaveSlater *waveSlater = new WaveSlater(config);
    config->setWave(waveSlater);
    waveSlater->setParameters(parameters);
    waveSlater->setUseAnalyticalGradient(true);
    waveSlater->setUseAnalyticalLaplace(true);

    Hamiltonian *hamiltonianIdeal = new HamiltonianIdeal(config);
    config->setHamiltonian(hamiltonianIdeal);

    DiffusionMonteCarlo *diffusionMonteCarlo = new DiffusionMonteCarlo(config);
    diffusionMonteCarlo->setTimeStep(0.0001);
    diffusionMonteCarlo->sample(40000);
    double energy = diffusionMonteCarlo->energy();
    std::cout << "Diffusion monte carlo returned energy of " << energy << std::endl;
    QVERIFY(fabs(energy - 3.000) < 1e-2);
}

void VmcTests::evolutionaryMonteCarloTest() {
    Config *config = new Config(0, 1);
    config->setNParticles(2);
    config->setNDimensions(2);
    config->setInteractionEnabled(true);

    double parameters[2];
    parameters[0] = 0.987;
    parameters[1] = 0.398;
    //    parameters[0] = 0.928;
    //    parameters[1] = 0.55;

    WaveSlater *waveSlater = new WaveSlater(config);
    config->setWave(waveSlater);
    waveSlater->setParameters(parameters);
    waveSlater->setUseAnalyticalGradient(true);
    waveSlater->setUseAnalyticalLaplace(true);

    Hamiltonian *hamiltonianIdeal = new HamiltonianIdeal(config);
    config->setHamiltonian(hamiltonianIdeal);

    int nWalkers = 4;
    int nGenes = nWalkers * config->nParticles() * config->nDimensions();
    int nIndividuals = 16;
    int nPopulations = 1;
    EvolutionaryMonteCarlo *evolutionary = new EvolutionaryMonteCarlo(config, nGenes, nIndividuals, nPopulations);
    evolutionary->setRescaleLimits(-2, 10);
    evolutionary->setRescaleCycles(5);
    evolutionary->sample(100000);
    double energy = evolutionary->energy();
    std::cout << "Evolutionary monte carlo returned energy of " << energy << std::endl;
    QVERIFY(fabs(energy - 3.000) < 1e-2);
}

QTEST_APPLESS_MAIN(VmcTests)

#include "tst_vmctests.moc"
