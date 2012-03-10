#ifndef DENSITYPLOTTER_H
#define DENSITYPLOTTER_H

class WaveFunction;
class INIReader;
class Config;
class Hamiltonian;

/*!
  A class that uses the given wave function to plot the onebody density.
  */

class DensityPlotter
{
public:
    DensityPlotter(Config *config);
    void makePlot();
    void loadConfiguration(INIReader *settings);
private:
    Config *m_config;
    WaveFunction *m_wave;
    Hamiltonian *m_hamiltonian;
    INIReader *m_settings;

    double m_charge;
    double m_stepLength;
    int m_nCycles;
};

#endif // DENSITYPLOTTER_H
