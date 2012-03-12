#include <QtCore/QString>
#include <QtTest/QtTest>

#include <stdio.h>

#include "../wavesimple.h"
#include "../waveideal.h"
#include "../matrix.h"
#include "../config.h"
#include "../hamiltonian/hamiltonianideal.h"
#include "../montecarlo/montecarlostandard.h"
#include "../montecarlo/montecarlometropolishastings.h"

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

private:
    Config *config;
    WaveSimple *waveSimple;
    WaveIdeal *waveIdeal;
    HamiltonianIdeal *hamiltonianIdeal;
    double **r_old;
    double charge;
};

VmcTests::VmcTests()
{
}

void VmcTests::initTestCase()
{
    config = new Config(0,1,2,2);
    // Set up waveSimple
    int nParticles = config->nParticles();
    int nDimensions = config->nDimensions();
    charge = 1.0;
    r_old = (double **) matrix( nParticles, nDimensions, sizeof(double));
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
    double** rPositions = (double **) matrix( nParticles, nDimensions, sizeof(double));
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < nDimensions; j++) {
            rPositions[i][j] = 1.0 - j;
        }
    }
    double rgrad[2];
    waveSimpleNew->gradient(r_old, rgrad);
    printf("Gradient: %.10f %.10f\n", rgrad[0], rgrad[1]);
    QVERIFY(fabs(rgrad[0] - 2) < 0.001);
}

void VmcTests::fullIdealTest()
{
    double alpha = 1.0;
    double beta = 0.4;
    double stepLength = 1.0;
    int nCycles = 500000;
    int nTotalCycles = nCycles;
    waveIdeal->setUseAnalyticalLaplace(true);
    waveIdeal->setParameters(alpha, beta);
    MonteCarloStandard *monteCarlo = new MonteCarloStandard(waveIdeal, hamiltonianIdeal, config->nParticles(), config->nDimensions(), charge, config->rank(), stepLength);
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
    double stepLength = 1.0;
    int nCycles = 500000;
    int nTotalCycles = nCycles;
    waveIdeal->setUseAnalyticalLaplace(true);
    waveIdeal->setParameters(alpha, beta);
    MonteCarloMetropolisHastings *monteCarlo = new MonteCarloMetropolisHastings(waveIdeal, hamiltonianIdeal, config->nParticles(), config->nDimensions(), charge, config->rank(), stepLength);
    double *allEnergies = new double[nCycles+1];
    double *energies = new double[2];
    //  Do the mc sampling
    monteCarlo->sample(nCycles, energies, allEnergies);
    QVERIFY(fabs(energies[0] / nTotalCycles - 3.00034530284643397025) < 1e-20);
}

QTEST_APPLESS_MAIN(VmcTests)

#include "tst_vmctests.moc"
