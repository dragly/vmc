#ifndef DENSITYPLOTTER_H
#define DENSITYPLOTTER_H

class WaveFunction;

/*!
  A class that uses the given wave function to plot the onebody density.
  */

class DensityPlotter
{
public:
    DensityPlotter(int rank, int nProcesses);
    void makePlot();
private:
    int m_rank;
    int m_nProcesses;
    WaveFunction *wave;
};

#endif // DENSITYPLOTTER_H
