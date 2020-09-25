#ifndef FFT_h
#define FFT_h

#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

extern void fft(float *x_in, complex<double> *x_out, int N);

void fft_rec(complex<double> *x, int N);

#endif
