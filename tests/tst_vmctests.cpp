#include <QtCore/QString>
#include <QtTest/QtTest>

#include <stdio.h>
#include <iostream>
#include <armadillo>

#include "../wavefunction/wavesimple.h"
#include "../wavefunction/waveideal.h"
#include "../matrix.h"
#include "../config.h"
#include "../hamiltonian/hamiltonianideal.h"
#include "../montecarlo/montecarlostandard.h"
#include "../montecarlo/montecarlometropolishastings.h"
#include "../minimizer/minimizerevolutionary.h"
#include "../hermite.h"
#include "../orbital/orbital.h"
#include "../slater/slater.h"
#include "../wavefunction/waveslater.h"
#include "../jastrow/jastrow.h"

#include "minimizerevolutionarytest.h"

using namespace std;
using namespace arma;

class VmcTests : public QObject
{
    Q_OBJECT
    
public:
    VmcTests();

    void waveSimpleGradientTest(); // TODO - consider implementing this again
    void initTestCase();
    void orbitalTest();
    void jastrowTest();
    void cleanupTestCase();
    void waveSimpleLaplaceTest();
    void waveIdealLaplaceTest();
    void fullIdealTest();
    void fullIdealHastingsTest();
    void hermiteTest();
    void minimizerEvolutionaryTest();
    void twoOrbitalsOneWavefunctionTest();
    void slaterTest();
    void fullIdealHastingsSlaterTest();
    void slaterFourParticleTest();
    void slaterSixParticleTest();
    void fullSlaterSixNoInteractionTest();
private slots:
    void fullSlaterSixInteractionTest();

private:
    Config *config;
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
    config = new Config(0,1);
    // Set up waveSimple
    int nParticles = config->nParticles();
    cout << nParticles << endl;
    int nDimensions = config->nDimensions();
    charge = 1.0;
    r_old = new vec2[nParticles];
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            r_old[i][j] = 0.234 + i + 2*j;
        }
    }
    waveSimple = new WaveSimple(config);
    double parameters[2];
    parameters[0] = 2;
    parameters[1] = 1;
    waveSimple->setParameters(parameters);
    // Set up waveIdeal
    waveIdeal = new WaveIdeal(config);
    waveIdeal->setParameters(parameters);
    // Set up hamiltonianIdeal
    hamiltonianIdeal = new HamiltonianIdeal(config);
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
    WaveSimple *waveSimpleNew = new WaveSimple(config);
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
    int nTotalCycles = nCycles;
    config->setWave(waveIdeal);
    config->setHamiltonian(hamiltonianIdeal);
    waveIdeal->setUseAnalyticalLaplace(true);
    double parameters[2];
    parameters[0] = 1;
    parameters[1] = 0.4;
    waveIdeal->setParameters(parameters);
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(config);
    double *allEnergies = new double[nCycles+1];
    double *energies = new double[2];
    //  Do the mc sampling
    monteCarlo->sample(nCycles, energies, allEnergies);
    QVERIFY(fabs(energies[0] / nTotalCycles - 3.00034530284643397025) < 1e-20);
}

void VmcTests::fullIdealHastingsTest()
{
    double parameters[2];
    parameters[0] = 1.0;
    parameters[1] = 0.4;
    int nCycles = 100000;
    int nTotalCycles = nCycles;
    config->setWave(waveIdeal);
    config->setHamiltonian(hamiltonianIdeal);
    waveIdeal->setUseAnalyticalLaplace(true);
    waveIdeal->setParameters(parameters);
    MonteCarloMetropolisHastings *monteCarlo = new MonteCarloMetropolisHastings(config);
    double *allEnergies = new double[nCycles+1];
    double *energies = new double[2];
    //  Do the mc sampling
    monteCarlo->sample(nCycles, energies, allEnergies);
    QVERIFY(fabs(energies[0] / nTotalCycles - 3.000) < 1e-2);
}

void VmcTests::hermiteTest() {
//    Hermite* hermite5 = new Hermite(5);

    QVERIFY(Hermite::evaluate(2,4) - (4*4*4 - 2) < 1e-20);
    QVERIFY(Hermite::evaluate(3,5) - (8*5*5*5 - 12) < 1e-20);
    QVERIFY(Hermite::evaluate(4,3) - (16*3*3*3*3 - 48*3*3 + 12) < 1e-20);
}

void VmcTests::minimizerEvolutionaryTest() {
    cout << "minimizerEvolutionaryTest not implemented" << endl;
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

    double val1 = waveSimple1->wave(rpos);

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

    double val1 = waveSimple1->wave(rpos);

    Slater *slater1 = new Slater(config1);

    Orbital *orbitals[1];
    orbitals[0] = new Orbital(0,0,config1);
    orbitals[0]->setParameters(parameters);

    double val2 = slater1->determinant(rpos, orbitals);

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
    config1->setHamiltonian(hamiltonianIdeal);
    WaveSlater *waveSlater1 = new WaveSlater(config1);
    config1->setWave(waveSlater1);
    double parameters[2];
    parameters[0] = 1.0;
    parameters[1] = 0.4;
    waveSlater1->setParameters(parameters);
    int nCycles = 100000;
    int nTotalCycles = nCycles;
//    waveSlater1->setUseAnalyticalLaplace(false);
    MonteCarloMetropolisHastings *monteCarlo = new MonteCarloMetropolisHastings(config1);
    double *allEnergies = new double[nCycles+1];
    double *energies = new double[2];
    //  Do the mc sampling
    monteCarlo->sample(nCycles, energies, allEnergies);
    cout << "Two particle energy is " << fabs(energies[0] / nTotalCycles) << endl;
    QVERIFY(fabs(energies[0] / nTotalCycles - 3.000) < 1e-2);
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
    QCOMPARE(orbital11->evaluate(r), A * 2 * sqrt(omega) * r[0] * 2 * sqrt(omega) * r[1] * exp(- parameters[0] * omega * (r[0] * r[0] + r[1] * r[1]) / 2));
    double sqrtomegax = (sqrt(omega) * r[0]);
    double sqrtomegay = (sqrt(omega) * r[1]);
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
    QCOMPARE(waveSlater1->wave(r), waveValue);
}
// TODO: Test slater wave for four particles

/*!
  Test the slater wave function with multiple particles. There is no benchmark
  for this function, so we just make sure it runs and prints without error.
  */
void VmcTests::slaterSixParticleTest()
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
    cout << "Value for slaterSixParticleTest: " << waveSlater1->wave(rpos) << endl;
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
    HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(config1);
    config1->setHamiltonian(hamiltonian);
    WaveSlater *waveSlater1 = new WaveSlater(config1);
    config1->setWave(waveSlater1);
    double parameters[2];
    parameters[0] = 1.0;
    parameters[1] = 1.0;
    waveSlater1->setParameters(parameters);
    int nCycles = 100000;
    int nTotalCycles = nCycles;
//    waveSlater1->setUseAnalyticalLaplace(false);
    MonteCarloStandard *monteCarlo1 = new MonteCarloStandard(config1);
    double *allEnergies = new double[nCycles+1];
    double *energies = new double[2];
    //  Do the mc sampling
    monteCarlo1->sample(nCycles, energies, allEnergies);
    cout << "Six non-interacting energy was " << fabs(energies[0] / nTotalCycles) << endl;
    QVERIFY(fabs(energies[0] / nTotalCycles - 10) < 1e-2);
}

/*!
  Test the Metropolis-Hastings algorithm by using the Slater determinant with 6 particles.
  This is benchmarked against the optimal result found in Lars-Eivind Lervåg's thesis.

  This test is slow, but gives a decent benchmark.
  */
void VmcTests::fullSlaterSixInteractionTest()
{
    QBENCHMARK{
        Config *config1 = new Config(0,1);
        config1->setNDimensions(2);
        config1->setNParticles(6);
        HamiltonianIdeal *hamiltonian = new HamiltonianIdeal(config1);
        config1->setHamiltonian(hamiltonian);
        config1->setInteractionEnabled(true);
        int nCycles = 5000;
        int nTotalCycles = nCycles;

        WaveSlater *waveSlater2 = new WaveSlater(config1);
        double parameters[2];
        parameters[0] = 0.92;
        parameters[1] = 0.565;
        waveSlater2->setParameters(parameters);
        config1->setWave(waveSlater2);
        MonteCarloStandard *monteCarlo2 = new MonteCarloStandard(config1);
        double *allEnergies = new double[nCycles+1];
        double *energies = new double[2];
        monteCarlo2->sample(nCycles, energies, allEnergies);
        cout << "Six interacting energy was " << fabs(energies[0] / nTotalCycles) << endl;
//        QVERIFY(fabs(energies[0] / nTotalCycles - 20.190) < 1e-2);
    }
}

QTEST_APPLESS_MAIN(VmcTests)

#include "tst_vmctests.moc"
