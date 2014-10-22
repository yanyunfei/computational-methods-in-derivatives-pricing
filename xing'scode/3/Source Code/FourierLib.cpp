#include "FourierLib.h"

using namespace std;

// N must be a power of 2

dcmplx* FourierFuncs::FFT(dcmplx* x, int N)
{
    dcmplx* out = (dcmplx*) malloc(sizeof(dcmplx) * N);
    dcmplx* scratch = (dcmplx*) malloc(sizeof(dcmplx) * N);
    dcmplx* twiddles = FFT_get_twiddle_factors(N);
    FFT_calculate(x, N, 1, out, scratch, twiddles);
    free(twiddles);
    free(scratch);
    return out;
}


dcmplx* FourierFuncs::IFFT(dcmplx* x, int N)
{
    dcmplx* out = (dcmplx*) malloc(sizeof(dcmplx) * N);
    dcmplx* scratch = (dcmplx*) malloc(sizeof(dcmplx) * N);
    dcmplx* twiddles = FFT_get_twiddle_factors(N);
    FFT_calculate(x, N, 1, out, scratch, twiddles);
    free(twiddles);
    free(scratch);
    // Calculate IFFT via reciprocity property of DFT.
    int N2 = N/2;
    double tmp0, tmp1;
    out[0] = dcmplx(real(out[0])/N, imag(out[0])/N); out[N2] = dcmplx(real(out[N2])/N, imag(out[N2])/N);
    for(int i=1; i<N2; i++) {
        tmp0 = real(out[i])/N;
        tmp1 = imag(out[i])/N;
        out[i] = dcmplx(real(out[N-i])/N, imag(out[N-i])/N);
        out[N-i] = dcmplx(tmp0, tmp1);
    }
    return out;
}


void FourierFuncs::FFT_calculate(dcmplx* x, int N, int skip, dcmplx* X, dcmplx* D, dcmplx* twiddles) {
    dcmplx * E = D + N/2;
    int k;
    if (N == 1) {
        X[0] = x[0];
        return;
    }
    // use X as a scratch buffer
    FFT_calculate(x, N/2, skip*2, E, X, twiddles);
    FFT_calculate(x + skip, N/2, skip*2, D, X, twiddles);
    for(k = 0; k < N/2; k++) {
        // Multiply entries of D by the twiddle factors e^(-2*pi*i/N * k)
        D[k] = twiddles[k*skip] * D[k];
    }
    for(k = 0; k < N/2; k++) {
        X[k]       = E[k] + D[k];
        X[k + N/2] = E[k] - D[k];
    }
}


dcmplx* FourierFuncs::FFT_get_twiddle_factors(int N) {
    dcmplx* twiddles = (dcmplx*) malloc(sizeof(dcmplx) * N/2);
    int k;
    for (k = 0; k != N/2; ++k) {
        twiddles[k] = exp(dcmplx(0, -2.0*PI*k/N));
    }
    return twiddles;
}






