#include "hermite.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

/*!
  Constructs a new Hermite polynomial of a given degree.

  @param degree Degree of the polynomial.
  */

Hermite::Hermite(int degree) :
    m_degree(degree)
{
    if(degree > 4 || degree < 0) {
        cerr << __PRETTY_FUNCTION__ << ": Hermite polynomial of unknown degree called. Degree was: " << m_degree << endl;
        exit(1);
    }
}

/*!
  Evaluate the polynomial at the given value based on the polynomial's degree.

  @param x Value at which to evaluate the polynomial.
  @returns Value of the polynomial at the given x
  */
double Hermite::evaluate(double x) {
    switch(m_degree) {
    case 0:
        return 1;
        break;
    case 1:
        return 2*x;
        break;
    case 2:
        return 4*x*x - 2;
        break;
    case 3:
        return 8*x*x*x - 12;
        break;
    case 4:
        return 16*x*x*x*x - 48*x*x + 12;
        break;
    default:
        cerr << __PRETTY_FUNCTION__ << ": Hermite polynomial of unknown degree called. Degree was: " << m_degree << endl;
        exit(1);
        return 0;
    }
}
