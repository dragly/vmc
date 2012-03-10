#ifndef DENSITYPLOTTER_H
#define DENSITYPLOTTER_H

class WaveFunction;
class INIReader;
class Config;

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
    INIReader *m_settings;
};

#endif // DENSITYPLOTTER_H
