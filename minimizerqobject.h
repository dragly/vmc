#ifndef MINIMIZERQOBJECT_H
#define MINIMIZERQOBJECT_H

#include <QObject>

class Minimizer;

class MinimizerQObject : public QObject
{
    Q_OBJECT
public:
    explicit MinimizerQObject(int rank, int nProcesses, QObject *parent = 0);


signals:

public slots:
    void runMinimizer();

private:
    Minimizer *minimizer;

    int rank;
    int nProcesses;
};

#endif // MINIMIZERQOBJECT_H
