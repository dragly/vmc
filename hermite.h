#ifndef HERMITE_H
#define HERMITE_H

class Hermite
{
public:
    Hermite();
    static double evaluate(int degree, double x);
    static double derivative(int degree, double x);
    static double doubleDerivative(int degree, double x);
};

#endif // HERMITE_H
