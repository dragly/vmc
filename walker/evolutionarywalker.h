#ifndef EVOLUTIONARYWALKER_H
#define EVOLUTIONARYWALKER_H

#include "walker.h"

class EvolutionaryWalker : public Walker
{
public:
    EvolutionaryWalker(Config *config);

    void advance();

private:

};

#endif // EVOLUTIONARYWALKER_H
