#include <QtCore/QString>
#include <QtTest/QtTest>

#include "../wavesimple.h"
#include "../waveideal.h"
#include "../matrix.h"

class VmcTests : public QObject
{
    Q_OBJECT
    
public:
    VmcTests();
    
private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void waveSimpleLaplaceTest();
    void waveIdealLaplaceTest();

private:
    WaveSimple *waveSimple;
    WaveIdeal *waveIdeal;
    double **r_old;
};

VmcTests::VmcTests()
{
}

void VmcTests::initTestCase()
{
    // Set up waveSimple
    int nParticles = 2;
    int dimensions = 2;
    r_old = (double **) matrix( nParticles, dimensions, sizeof(double));
    for (int i = 0; i < nParticles; i++) {
        for (int j=0; j < dimensions; j++) {
            r_old[i][j] = 0.234 + i + 2*j;
        }
    }
    waveSimple = new WaveSimple(nParticles,dimensions);
    waveSimple->setParameters(2, 1);
    // Set up waveIdeal
    waveIdeal = new WaveIdeal(nParticles,dimensions);
    waveIdeal->setParameters(2, 1);
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

QTEST_APPLESS_MAIN(VmcTests)

#include "tst_vmctests.moc"
