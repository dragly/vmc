#ifndef DENSITYPLOTTER_H
#define DENSITYPLOTTER_H
#include <fstream>

class WaveFunction;
class INIReader;
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
    void loadConfiguration(INIReader *settings);
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
    INIReader *m_settings;
    long idum;
    double **r_old;
    double **r_new;
    double m_charge;
    double m_stepLength;
    int m_nCycles;
    ofstream plotFile;
};

#endif // DENSITYPLOTTER_H
