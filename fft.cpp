#include "fft.h"
#include <iostream>

using namespace std;


void fft(float *x_in, complex<double> *x_out, int N) {

	//cout << "HERE fft start" << endl;
	// Make copy of array and apply window
	for (int i = 0; i < N; i++) {
		x_out[i] = complex<double>(x_in[i], 0);
		x_out[i] *= 1; // Window
	}

	// Start recursion
	fft_rec(x_out, N);
	cout << "FFT Performed" << endl;
}

void fft_rec(complex<double> *x, int N) {
	// Check if it is splitted enough
	if (N <= 1) {
		return;
	}

	// Split even and odd
	complex<double> odd[N/2];
	complex<double> even[N/2];
	for (int i = 0; i < N / 2; i++) {
		even[i] = x[i*2];
		odd[i] = x[i*2+1];
	}

	// Split on tasks
	fft_rec(even, N/2);
	fft_rec(odd, N/2);


	// Calculate DFT
	for (int k = 0; k < N / 2; k++) {
		complex<double> t = exp(complex<double>(0, -2 * M_PI * k / N)) * odd[k];
		x[k] = even[k] + t;
		x[N / 2 + k] = even[k] - t;
	}
}
