#ifndef HERMITE_H
#define HERMITE_H

class Hermite
{
public:
    Hermite(int degree);
    double evaluate(double x);
private:
    int m_degree;
};

#endif // HERMITE_H
