#ifndef __CaseStudy1__FourierFuncs__
#define __CaseStudy1__FourierFuncs__

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <complex>
#include <stdlib.h>

typedef std::complex<double> dcmplx;
#define PI  3.1415926535897932

//FourierLib includes static functions to do fourier related computations
class FourierFuncs{
public:
    static dcmplx* FFT(dcmplx*, int);
    static dcmplx* IFFT(dcmplx*, int);
    static void FFT_calculate(dcmplx*, int, int, dcmplx*, dcmplx*, dcmplx*);
    static dcmplx* FFT_get_twiddle_factors(int);
};

#endif
