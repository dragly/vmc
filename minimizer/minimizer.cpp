// Stat stuff
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
// math stuff
#include <math.h>
// stream stuff
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// own header
#include "minimizer.h"

// local includes


Minimizer::Minimizer(Config *config) :
    config(config),
    m_rank(config->myRank()),
    m_m_nProcesses(config->m_nProcesses())
{
}
