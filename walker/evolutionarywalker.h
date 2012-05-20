#ifndef EVOLUTIONARYWALKER_H
#define EVOLUTIONARYWALKER_H

#include "walker.h"

class EvolutionaryWalker : public Walker
{
public:
    EvolutionaryWalker(Config *config);

    void advance();

private:
    int nWalkersMax;

    double diffConstant;
    double tau;
    bool m_aliveNew;
    bool m_aliveOld;

};

#endif // EVOLUTIONARYWALKER_H
