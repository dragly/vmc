#include "minimizerqobject.h"

#include "minimizerstandard.h"

MinimizerQObject::MinimizerQObject(int rank, int nProcesses, QObject *parent) :
    QObject(parent),
    rank(rank),
    nProcesses(nProcesses)
{
    cout << "Starting minimizer..." << endl;
    minimizer = new MinimizerStandard(rank, nProcesses);
}

void MinimizerQObject::runMinimizer()
{
    cout << "MinimizerQObject::runMinimizer(): called" << endl;
    minimizer->runMinimizer();
}
