#include "minimizerevolutionarytest.h"

MinimizerEvolutionaryTest::MinimizerEvolutionaryTest(Config *config) : MinimizerEvolutionary(config)
{
}

double MinimizerEvolutionaryTest::value(vec *coefficients) {
    return 2;
}
