#ifndef DENSITYPLOTTER_H
#define DENSITYPLOTTER_H
#include <fstream>
#include <armadillo>
using namespace arma;

class WaveFunction;
class INIParser;
class Config;

using namespace std;
//class Hamiltonian;

/*!
  A class that uses the given wave function to plot the onebody density.
  */

class StepConfig {
public:
    int firstStep;
    int lastStep;
    int nSteps;
};

class DensityPlotter
{
    //Q_PROPERTY(double charge READ charge WRITE setCharge)
public:
    DensityPlotter(Config *config);
    void makePlot();
    void loadConfiguration(INIParser *settings);
    ~DensityPlotter();

    double charge() const {
        return m_charge;
    }

//public slots:
    void setCharge(double arg) {
        m_charge = arg;
    }

    void divideSteps(int rank, int nProcesses, int totalSteps, StepConfig *stepConfig);
private:
    Config *m_config;
    WaveFunction *m_wave;
    INIParser *m_settings;
    long idum;
    vec2 *r_old;
    vec2 *r_new;
    double m_charge;
    double m_stepLength;
    double aMin;
    double bMin;
    double aMax;
    double bMax;
    int aSteps;
    int bSteps;
    int m_nCycles;
    ofstream plotFile;
    ofstream params0File;
    ofstream params1File;
};

#endif // DENSITYPLOTTER_H
