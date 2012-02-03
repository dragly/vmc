#ifndef MINIMIZER_H
#define MINIMIZER_H

/*!
  Runs with multiple parameters and attempts to find the lowest possible value of the energy by adjusting the parameters.
  */

class Minimizer
{
public:
    Minimizer(int rank, int nProcesses);
    virtual void run() = 0;
private:
    int rank;
    int nProcesses;
};

#endif // MINIMIZER_H
