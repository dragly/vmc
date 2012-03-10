#include "densityplotter.h"

#include "wavefunction.h"
#include "inih/cpp/INIReader.h"
#include "config.h"

DensityPlotter::DensityPlotter(Config *config) :
    m_config(config)
{

}

void DensityPlotter::makePlot()
{

}

void DensityPlotter::loadConfiguration(INIReader *settings)
{
    m_settings = settings;
    m_wave = WaveFunction::fromName(settings->Get("Wave", "class", "WaveSimple"), m_config);
}
