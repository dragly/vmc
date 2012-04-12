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

#include "minimizerevolutionarytest.h"

using namespace std;
using namespace arma;

class VmcTests : public QObject
{
    Q_OBJECT
    
public:
    VmcTests();

private slots:
    void initTestCase();
    void cleanupTestCase();
    void waveSimpleLaplaceTest();
    void waveIdealLaplaceTest();
    void fullIdealTest();
    void fullIdealHastingsTest();
    void waveSimpleGradientTest();
    void hermiteTest();
    void minimizerEvolutionaryTest();
    void twoOrbitalsOneWavefunction();

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
    waveSimple = new WaveSimple(nParticles,nDimensions);
    waveSimple->setParameters(2, 1);
    // Set up waveIdeal
    waveIdeal = new WaveIdeal(nParticles,nDimensions);
    waveIdeal->setParameters(2, 1);
    // Set up hamiltonianIdeal
    hamiltonianIdeal = new HamiltonianIdeal(nParticles, nDimensions, charge);
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

    int nParticles = 1;
    int nDimensions = 2;
    WaveSimple *waveSimpleNew = new WaveSimple(nParticles,nDimensions);
    vec2* rPositions = new vec2[nParticles];
    rPositions[0][0] = -1;
    rPositions[0][1] = 0.0;
    vec2 rGradient;
    waveSimpleNew->setParameters(2, 1);
    waveSimpleNew->gradient(rPositions, rGradient);
    QVERIFY(fabs(rGradient[0] - 0.735) < 0.001);
    QVERIFY(fabs(rGradient[1] - 0.000) < 0.0000001);
}

void VmcTests::fullIdealTest()
{
    double alpha = 1.0;
    double beta = 0.4;
    int nCycles = 500000;
    int nTotalCycles = nCycles;
    config->setWave(waveIdeal);
    config->setHamiltonian(hamiltonianIdeal);
    waveIdeal->setUseAnalyticalLaplace(true);
    waveIdeal->setParameters(alpha, beta);
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(config);
    double *allEnergies = new double[nCycles+1];
    double *energies = new double[2];
    //  Do the mc sampling
    monteCarlo->sample(nCycles, energies, allEnergies);
    QVERIFY(fabs(energies[0] / nTotalCycles - 3.00034530284643397025) < 1e-20);
}

void VmcTests::fullIdealHastingsTest()
{
    double alpha = 1.0;
    double beta = 0.4;
    int nCycles = 100000;
    int nTotalCycles = nCycles;
    config->setWave(waveIdeal);
    config->setHamiltonian(hamiltonianIdeal);
    waveIdeal->setUseAnalyticalLaplace(true);
    waveIdeal->setParameters(alpha, beta);
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
    cout << "test" << endl;
}

void VmcTests::twoOrbitalsOneWavefunction() {
    vec2 *rpos = new vec2[config->nParticles()];
    rpos[0][0] = 0.4;
    rpos[0][1] = 0.9;
    rpos[1][0] = 0.1;
    rpos[1][1] = 0.8;

    WaveSimple *waveSimple1 = new WaveSimple(config->nParticles(), config->nDimensions());
    waveSimple1->setParameters(1,1);

    Orbital *orbital1 = new Orbital(0,0,config);
    Orbital *orbital2 = new Orbital(0,0,config);
    orbital1->setParameters(1,1);
    orbital2->setParameters(1,1);

    double val1 = waveSimple1->wave(rpos);

    double val2 = orbital1->evaluate(rpos[0]) * orbital2->evaluate(rpos[1]);

    QCOMPARE(val1,val2);
}

QTEST_APPLESS_MAIN(VmcTests)

#include "tst_vmctests.moc"
