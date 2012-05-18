#ifndef MONTECARLOSTANDARD_H
#define MONTECARLOSTANDARD_H
#include <armadillo>
using namespace arma;

#include "montecarlo.h"
#include "../hamiltonian/hamiltonian.h"
#include "../wavefunction/wavefunction.h"
#include "../config.h"

class MonteCarloStandard : public MonteCarlo
{
public:
    MonteCarloStandard(Config* config);
    void setRecordMoves(bool arg, int nMoves = 0);

    void sample(int nCycles);

    vec2 **moves() {
        return m_moves;
    }

    ~MonteCarloStandard();
private:
    int rank;
    double stepLength;
    vec2 *rOld;
    vec2 *rNew;

    bool recordMoves;
    int nMoves;
    vec2 **m_moves;
};

#endif // MONTECARLOSTANDARD_H
