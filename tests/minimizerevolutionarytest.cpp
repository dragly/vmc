#include "minimizerevolutionarytest.h"

#include "../minimizer/minimizerevolutionary.h"

MinimizerEvolutionaryTest::MinimizerEvolutionaryTest(Config *config) : MinimizerEvolutionary(config)
{
}

double MinimizerEvolutionaryTest::value() {
    return 2;
}
