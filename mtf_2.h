#ifndef MTF_2_H
#define MTF_2_H

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <complex.h>
#include <iomanip>
#include "fft.h"

using namespace cv;
using namespace std;

struct detected_circle
{
         float radius;
         float center_x;
         float center_y;
};

detected_circle detect_circle(struct detected_circle main_circle);
int obtain_distance_of_pixels_from_center(vector<float> &dist, vector<unsigned int> &intensity, struct detected_circle main_circle);
void print_vector_float(vector<float> &vector);
void print_vector_uint(vector<unsigned int> &vector);
void print_distance_intensity(vector<float> &dist, vector<unsigned int> &intensity);
void sort_pixel_by_ascending_distance(vector<float> &dist, vector<unsigned int> &intensity);
void image_binning_averaging(vector<float> &dist, vector<unsigned int> &intensity, unsigned int bin_size, vector<float> &averaged_distance, vector<unsigned int> &averaged_intensity);
void write_to_excel(vector<float> &dist, vector<float> &intensity);
void write_to_excel(vector<float> &dist, vector<unsigned int> &intensity);
void grouping_for_polynomial_fit(vector<float> &averaged_distance, vector<unsigned int> &averaged_intensity, unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &erf_cubic_fit_distance, vector<unsigned int> &erf_cubic_fit_intensity);
void least_square_polynomial_fit(float x[], int y[], unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &erf_cubic_fit_distance, vector<unsigned int> &erf_cubic_fit_intensity);
void grouping_for_derivative_of_cubic_fit(vector<float> &erf_cubic_fit_distance, vector<unsigned int> &erf_cubic_fit_intensity, unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &psf_cubic_fit_distance, vector<float> &psf_cubic_fit_intensity);
void derivative_of_cubic_fit(float x[], int y[], unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &psf_cubic_fit_distance, vector<float> &psf_cubic_fit_intensity);
void normalise_psf(vector<float> &psf_cubic_fit_distance, vector<float> &psf_cubic_fit_intensity);
void perform_fft(vector<float> &psf_cubic_fit_intensity, complex<double> *x_out);
void print_complex_array(complex<double> *x_out, unsigned int size);
void obtain_magnitude_of_fft_values(complex<double> *x_out, unsigned int size, vector<double> &fft_magnitude);
void normalize_fft(vector<double> &fft_magnitude);
void sort_fft_value(vector<double> &fft_magnitude);
void write_fft_to_excel(vector<double> &fft_magnitude, float frequency);
#endif
