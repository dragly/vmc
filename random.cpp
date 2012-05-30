#include "random.h"

/*
** The function
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of
** L'Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
    int            j;
    long           k;
    static long    idum2 = 123456789;
    static long    iy=0;
    static long    iv[NTAB];
    double         temp;

    if(*idum <= 0) {
        if(-(*idum) < 1) *idum = 1;
        else             *idum = -(*idum);
        idum2 = (*idum);
        for(j = NTAB + 7; j >= 0; j--) {
            k     = (*idum)/IQ1;
            *idum = IA1*(*idum - k*IQ1) - k*IR1;
            if(*idum < 0) *idum +=  IM1;
            if(j < NTAB)  iv[j]  = *idum;
        }
        iy=iv[0];
    }
    k     = (*idum)/IQ1;
    *idum = IA1*(*idum - k*IQ1) - k*IR1;
    if(*idum < 0) *idum += IM1;
    k     = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if(idum2 < 0) idum2 += IM2;
    j     = iy/NDIV;
    iy    = iv[j] - idum2;
    iv[j] = *idum;
    if(iy < 1) iy += IMM1;
    if((temp = AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*
** The function
** ran3()
** returns a uniform random number deviate between 0.0 and 1.0. Set
** the idum to any negative value to initialize or reinitialize the
** sequence. Any large MBIG, and any small (but still large) MSEED
** can be substituted for the present values.
*/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
static int inext, inextp;
static long ma[56]; // value 56 is special, do not modify
static int iff = 0;
long mj, mk;
int i, ii, k;

if(*idum < 0 || iff == 0) { // initialization
  iff = 1;

  mj = MSEED - (*idum < 0 ? -*idum : *idum);
  mj %= MBIG;
  ma[55] = mj; // initialize ma[55]

  for(i = 1, mk = 1; i <= 54; i++) { // initialize rest of table
     ii = (21*i) % 55;
ma[ii] = mk;
mk = mj - mk;
if(mk < MZ) mk += MBIG;
mj = ma[ii];
  }

  for(k = 1; k <= 4; k++) { // randimize by "warming up" the generator
     for(i = 1; i <= 55; i++) {
ma[i] -= ma[1 + (i + 30) % 55];
if(ma[i] < MZ) ma[i] += MBIG;
}
  }

  inext = 0; // prepare indices for first generator number
  inextp = 31; // 31 is special
  *idum = 1;
}

if(++inext == 56) inext = 1;
if(++inextp == 56) inextp = 1;
mj = ma[inext] - ma[inextp];
if(mj < MZ) mj += MBIG;
ma[inext] = mj;
return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/*!
 * \brief simpleGaussRandom
 * \param idum
 * \return Number between -1 and 1 with a close to gaussian distribution
 */
double simpleGaussRandom(long *idum) {
    return (-1 + 2 * ran2(idum)) + (-1 + 2 * ran2(idum)) + (-1 + 2 * ran2(idum));
}

// End: function ran2()
