#ifndef MINIMIZEREVOLUTIONARYTEST_H
#define MINIMIZEREVOLUTIONARYTEST_H

#include "../minimizer/minimizerevolutionary.h"

class MinimizerEvolutionaryTest : MinimizerEvolutionary
{
public:
    MinimizerEvolutionaryTest(Config *config);
    double value();
};

#endif // MINIMIZEREVOLUTIONARYTEST_H
