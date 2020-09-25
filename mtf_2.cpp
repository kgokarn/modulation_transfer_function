#include "mtf_2.h"

#define ROI_RADIUS 10
int main(int argc, char** argv)
{

	//declarations
	struct detected_circle main_circle;

	vector<float> dist;
	vector<unsigned int> intensity;

	vector<float> averaged_distance;
	vector<unsigned int> averaged_intensity;

	vector<float> erf_cubic_fit_distance;
	vector<unsigned int> erf_cubic_fit_intensity;

	vector<float> psf_cubic_fit_distance;
	vector<float> psf_cubic_fit_intensity;

	vector<double> fft_magnitude;

	complex<double> *x_out;

	unsigned int bin_size[] = {38,11,11};
	unsigned int degree_of_polynomial = 3;
	unsigned int fit_size = 33;
	unsigned int erf_size = 0;
	unsigned int psf_size = 0;
	float frequency = 0.01;


	// function calls

	// Obtaining the center and radius of the main circle
	main_circle = detect_circle(main_circle);
	if(main_circle.radius != -1){
		cout << "circle radius " << main_circle.radius << endl;
		cout << "circle X " << main_circle.center_x << endl;
		cout << "circle Y " << main_circle.center_y << endl;
	}

	//ROI extraction
	obtain_distance_of_pixels_from_center(dist, intensity, main_circle);
	sort_pixel_by_ascending_distance(dist, intensity);

	//Preprocessing
	image_binning_averaging(dist, intensity, bin_size[0], averaged_distance, averaged_intensity);
	dist.clear();
	intensity.clear();
	sort_pixel_by_ascending_distance(averaged_distance, averaged_intensity);

	//ERF
	grouping_for_polynomial_fit(averaged_distance, averaged_intensity, fit_size, degree_of_polynomial, erf_cubic_fit_distance, erf_cubic_fit_intensity);
	erf_size = erf_cubic_fit_distance.size();
	averaged_distance.clear();
	averaged_intensity.clear();
	cout << "ERF obtained" << endl;


	//PSF
	grouping_for_derivative_of_cubic_fit(erf_cubic_fit_distance, erf_cubic_fit_intensity, fit_size, degree_of_polynomial, psf_cubic_fit_distance, psf_cubic_fit_intensity);
	erf_cubic_fit_distance.clear();
	erf_cubic_fit_intensity.clear();
	normalise_psf(psf_cubic_fit_distance, psf_cubic_fit_intensity);
	write_to_excel(psf_cubic_fit_distance, psf_cubic_fit_intensity);
	cout << "PSF obtained" << endl;

	//MTF
	x_out = new complex<double>[psf_cubic_fit_intensity.size()]();
	perform_fft(psf_cubic_fit_intensity, x_out);
	psf_size = psf_cubic_fit_distance.size();
	psf_cubic_fit_distance.clear();
	psf_cubic_fit_intensity.clear();
	//print_complex_array(x_out, psf_size);
	obtain_magnitude_of_fft_values(x_out, psf_size, fft_magnitude);
	normalize_fft(fft_magnitude);
	sort_fft_value(fft_magnitude);
	write_fft_to_excel(fft_magnitude, frequency);
	cout << "MTF obtained" << endl;

    return 0;
}

detected_circle detect_circle(struct detected_circle main_circle){
		Mat src, src_gray;

		/// Read the image
		src = imread( "Input4.tif", 1 );

		if( !src.data ){
		  cout << "detect_circle() : Image not found" << endl;
		  main_circle.radius = -1;
		  main_circle.center_x = -1;
		  main_circle.center_y = -1;
		  return main_circle;
		}

		/// Convert it to gray
		cvtColor( src, src_gray, CV_BGR2GRAY );

		/// Reduce the noise so we avoid false circle detection
		GaussianBlur( src_gray, src_gray, Size(9, 9), 2, 4 );



		vector<Vec3f> circles;


		/// Apply the Hough Transform to find the circles
		HoughCircles( src_gray, circles, CV_HOUGH_GRADIENT, 1, src_gray.rows/8, 70, 70, 314, 2000 );

		cout << "Number of Circles detected : " <<  circles.size() << endl;
		/// Draw the circles detected
		for( size_t i = 0; i < circles.size(); i++ )
		{
			Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));

			int radius = cvRound(circles[i][2]);

			main_circle.center_x = center.x;
			main_circle.center_y = center.y;
			main_circle.radius =  radius;
			// circle center
			circle( src, center, 3, Scalar(0,255,0), -1, 8, 0 );
			// circle outline
			circle( src, center, radius, Scalar(0,0,255), 3, 8, 0 );
			circle( src, center, radius-ROI_RADIUS, Scalar(0,0,255), 3, 8, 0 );
			circle( src, center, radius+ROI_RADIUS, Scalar(0,0,255), 3, 8, 0 );
		 }


	    // Show your results
//	    namedWindow( "Hough Circle Transform Demo", CV_WINDOW_AUTOSIZE );
//	    imshow( "Hough Circle Transform Demo", src );

//	    waitKey(0);
	    return main_circle;

}


// Calculating the distance of each pixel from the center and
// putting it into a vector/arraylist
int obtain_distance_of_pixels_from_center(vector<float> &dist, vector<unsigned int> &intensity, struct detected_circle main_circle){
	Mat src, src_gray;
	float temp_variable_distance;
	int temp_variable_intensity;

	/// Read the image
	src = imread( "Input4.tif", 1 );

	if( !src.data || main_circle.radius == -1 ){
		cout << "obtain_distance_from_center () : either image or circle is not found" << endl;
		return -1;
	  }

	/// Convert it to gray
	cvtColor( src, src_gray, CV_BGR2GRAY );

	/// Reduce the noise so we avoid false circle detection
	GaussianBlur( src_gray, src_gray, Size(9, 9), 2, 4 );

	for(int i = 0; i < src_gray.cols ; i++){
	  for( int j = 0; j < src_gray.rows ; j++){
	    if( (pow((j-main_circle.center_y),2) + pow((i-main_circle.center_x),2))  <= pow((main_circle.radius+ROI_RADIUS),2)){
		  if((pow((j-main_circle.center_y),2) + pow((i-main_circle.center_x),2))  >= pow((main_circle.radius-ROI_RADIUS),2)){
			  temp_variable_intensity = src_gray.at<uchar>(j, i);
			  temp_variable_distance = (sqrt((pow((j-main_circle.center_y),2) + pow((i-main_circle.center_x),2))));
			  dist.push_back(temp_variable_distance);
			  intensity.push_back(temp_variable_intensity);
		  }

	    }

	  }

	}
//	cout << "EXITING obtain_distance_of_pixels_from_center" << endl;
	return 0;
}

void print_vector_float(vector<float> &vector){
	for (unsigned int i=0; i < vector.size(); i++){
	       cout << vector[i] << endl;
	}
}

void print_vector_uint(vector<unsigned int> &vector){
	for (unsigned int i=0; i < vector.size(); i++){
	       cout << vector[i] << endl;
	}
}

void print_distance_intensity(vector<float> &dist, vector<unsigned int> &intensity){
	cout << "Distance -- Intensity" << endl;
	if(dist.size() == intensity.size()){
		cout << "Vector Size : " << dist.size()  << endl;
		for (unsigned int i=0; i < dist.size(); i++){
		       cout << dist[i] << " -- " << intensity[i] << endl;
		}
	}else{
		cout << "Vector size of intensity and distance do not match" << endl;
	}
}

void least_square_polynomial_fit(float x[], int y[], unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &erf_cubic_fit_distance, vector<unsigned int> &erf_cubic_fit_intensity){
	int i,j,k,n,N;
	unsigned int temp_intensity = 0;
	N = fit_size;
	n = degree_of_polynomial;
	double X[2*n+1];                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	for (i=0;i<2*n+1;i++)
	{
		X[i]=0;
		for (j=0;j<N;j++)
			X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	}
	double B[n+1][n+2],a[n+1];            //B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
	for (i=0;i<=n;i++)
		for (j=0;j<=n;j++)
			B[i][j]=X[i+j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
	double Y[n+1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	for (i=0;i<n+1;i++)
	{
		Y[i]=0;
		for (j=0;j<N;j++)
		Y[i]=Y[i]+pow(x[j],i)*y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}
	for (i=0;i<=n;i++)
		B[i][n+1]=Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
	n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations


	for (i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
		for (k=i+1;k<n;k++)
			if (B[i][i]<B[k][i])
				for (j=0;j<=n;j++)
				{
					double temp=B[i][j];
					B[i][j]=B[k][j];
					B[k][j]=temp;
				}

	for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
		for (k=i+1;k<n;k++)
			{
				double t=B[k][i]/B[i][i];
				for (j=0;j<=n;j++)
					B[k][j]=B[k][j]-t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
			}
	for (i=n-1;i>=0;i--)                //back-substitution
	{                        //x is an array whose values correspond to the values of x,y,z..
		a[i]=B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
		for (j=0;j<n;j++)
			if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
				a[i]=a[i]-B[i][j]*a[j];
		a[i]=a[i]/B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	}

//	for(unsigned int k = 0; k < fit_size; k++){
//		cubic_fit_distance.push_back(x[k]);
//		if(k == (fit_size/2)){
//			for (int m=0;m<n;m++){
//				temp_intensity += a[m]*pow(x[k],m);
//			}
//		}else{
//			temp_intensity = y[k];
//		}
//		cubic_fit_intensity.push_back(temp_intensity);
//		temp_intensity = 0;
//	}

	// replacement method
	unsigned int center_fit = fit_size/2;
	for (int m=0;m<n;m++){
		temp_intensity += a[m]*pow(x[center_fit],m);
	}
//	cout << "-----------------------------------" << endl;
	if(temp_intensity < 255){
		erf_cubic_fit_distance.push_back(x[center_fit]);
		erf_cubic_fit_intensity.push_back(temp_intensity);
	}
	return;
}


void sort_pixel_by_ascending_distance(vector<float> &dist, vector<unsigned int> &intensity){
	unsigned int temp_intensity;
	float temp_distance;
	cout << "Number of pixel in ROI : " << dist.size() << endl;
	if(dist.size() == intensity.size()){
		for (unsigned int i=0; i < dist.size(); i++){
			for (unsigned int j=0; j < dist.size()-1; j++){
				if(dist[j] > dist[j+1]){
 					temp_distance = dist[j];
					dist[j] =  dist[j+1];
					dist[j+1] =  temp_distance;

					temp_intensity = intensity[j];
					intensity[j] =  intensity[j+1];
					intensity[j+1] =  temp_intensity;
				}
			}
		}
	}
}

void image_binning_averaging(vector<float> &dist, vector<unsigned int> &intensity, unsigned int bin_size, vector<float> &averaged_distance, vector<unsigned int> &averaged_intensity){
	unsigned int i,j;
	float temp_distance;
	unsigned int temp_intensity;
	unsigned int count;
	if(dist.size() == intensity.size()){
		for(i = 0; i <  dist.size() ; i++){
			temp_distance = 0;
			temp_intensity = 0;
			count = 0;
			for(j = 0; j < bin_size; j++){
				if((i + j) <= dist.size()){
					temp_distance += dist[i+j];
					temp_intensity += intensity[i+j];
					count++;

				}else{
					//cout << "BREAK " << count << endl;
					break;
				}
			}
			if(count == bin_size){
				temp_distance = temp_distance/count;
				temp_intensity = temp_intensity/count;

				averaged_distance.push_back(temp_distance);
				averaged_intensity.push_back(temp_intensity);
				//cout << temp_distance << " : " << temp_intensity << " : " << count << " : "  << averaged_distance.size() << " : " << averaged_intensity.size()  <<  endl;
			}
		}
	}
}

void write_to_excel(vector<float> &dist, vector<float> &intensity){
	ofstream ex_file;
    ex_file.open("Values.csv");
    if(dist.size() == intensity.size()){
    		for (unsigned int i=0; i < dist.size(); i++){
    			ex_file << dist[i] << "," << intensity[i] << endl;
    		}
	}
    ex_file.close();

}

void write_to_excel(vector<float> &dist, vector<unsigned int> &intensity){
	ofstream ex_file;
    ex_file.open("Values.csv");
    if(dist.size() == intensity.size()){
    		for (unsigned int i=0; i < dist.size(); i++){
    			ex_file << dist[i] << "," << intensity[i] << endl;
    		}
	}
    ex_file.close();

}

void grouping_for_polynomial_fit(vector<float> &averaged_distance, vector<unsigned int> &averaged_intensity, unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &erf_cubic_fit_distance, vector<unsigned int> &erf_cubic_fit_intensity){
	float *distance_poly_fit = new float[fit_size];
	int *intensity_poly_fit = new int[fit_size];

	for(unsigned int i = 0; i < averaged_distance.size() ; i=i+fit_size){
		for(unsigned int j = 0; j < fit_size; j++){
			distance_poly_fit[j] = averaged_distance[i+j];
			intensity_poly_fit[j] = averaged_intensity[i+j];
		}
		if((i+fit_size) >= averaged_distance.size()){
			break;
		}
		least_square_polynomial_fit(distance_poly_fit, intensity_poly_fit, fit_size, degree_of_polynomial, erf_cubic_fit_distance, erf_cubic_fit_intensity);
	}


}

void grouping_for_derivative_of_cubic_fit(vector<float> &erf_cubic_fit_distance, vector<unsigned int> &erf_cubic_fit_intensity, unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &psf_cubic_fit_distance, vector<float> &psf_cubic_fit_intensity){
	float *distance_poly_fit = new float[fit_size];
	int *intensity_poly_fit = new int[fit_size];

	for(unsigned int i = 0; i < erf_cubic_fit_distance.size() ; i=i+fit_size){
		for(unsigned int j = 0; j < fit_size; j++){
			distance_poly_fit[j] = erf_cubic_fit_distance[i+j];
			intensity_poly_fit[j] = erf_cubic_fit_intensity[i+j];
		}
		if((i+fit_size) >= erf_cubic_fit_distance.size()){
			break;
		}
		derivative_of_cubic_fit(distance_poly_fit, intensity_poly_fit, fit_size, degree_of_polynomial, psf_cubic_fit_distance, psf_cubic_fit_intensity);
	}


}

void derivative_of_cubic_fit(float x[], int y[], unsigned int fit_size, unsigned int degree_of_polynomial, vector<float> &psf_cubic_fit_distance, vector<float> &psf_cubic_fit_intensity){
	int i,j,k,n,N;
	float	 temp_intensity = 0;
	N = fit_size;
	n = degree_of_polynomial;
	double X[2*n+1];                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	for (i=0;i<2*n+1;i++)
	{
		X[i]=0;
		for (j=0;j<N;j++)
			X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	}
	double B[n+1][n+2],a[n+1];            //B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
	for (i=0;i<=n;i++)
		for (j=0;j<=n;j++)
			B[i][j]=X[i+j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
	double Y[n+1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	for (i=0;i<n+1;i++)
	{
		Y[i]=0;
		for (j=0;j<N;j++)
		Y[i]=Y[i]+pow(x[j],i)*y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}
	for (i=0;i<=n;i++)
		B[i][n+1]=Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
	n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations


	for (i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
		for (k=i+1;k<n;k++)
			if (B[i][i]<B[k][i])
				for (j=0;j<=n;j++)
				{
					double temp=B[i][j];
					B[i][j]=B[k][j];
					B[k][j]=temp;
				}

	for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
		for (k=i+1;k<n;k++)
			{
				double t=B[k][i]/B[i][i];
				for (j=0;j<=n;j++)
					B[k][j]=B[k][j]-t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
			}
	for (i=n-1;i>=0;i--)                //back-substitution
	{                        //x is an array whose values correspond to the values of x,y,z..
		a[i]=B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
		for (j=0;j<n;j++)
			if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
				a[i]=a[i]-B[i][j]*a[j];
		a[i]=a[i]/B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	}

//	for(unsigned int k = 0; k < fit_size; k++){
//		cubic_fit_distance.push_back(x[k]);
//		if(k == (fit_size/2)){
//			for (int m=0;m<n;m++){
//				temp_intensity += a[m]*pow(x[k],m);
//			}
//		}else{
//			temp_intensity = y[k];
//		}
//		cubic_fit_intensity.push_back(temp_intensity);
//		temp_intensity = 0;
//	}

	// replacement method
	unsigned int center_fit = fit_size/2;
	for (int m=0;m<n;m++){
		// if m-1=-1 then m is zero hence a constant and its derivative is zero. thus ignored.
		if((m-1)>=0){
			temp_intensity = temp_intensity + m*a[m]*pow(x[center_fit],m-1);
			//cout << " Temp intensity : " << temp_intensity << " : " << (m-1) << " : " << a[m] << " : "<< x[center_fit] << " : " << pow(x[center_fit],m-1) << " : "<<  (m*a[m]*pow(x[center_fit],m-1)) << endl;
		}
	}
	//cout << "----------------------" << endl;
	psf_cubic_fit_distance.push_back(x[center_fit]);
	psf_cubic_fit_intensity.push_back((-1)*temp_intensity);
	return;
}


void normalise_psf(vector<float> &psf_cubic_fit_distance, vector<float> &psf_cubic_fit_intensity){
	float max_intensity = FLT_MIN;
	unsigned int i = 0;
	for(i = 0; i < psf_cubic_fit_intensity.size(); i++){
		if(psf_cubic_fit_intensity[i] > max_intensity){
			max_intensity = psf_cubic_fit_intensity[i];
		}
	}

	for(i = 0; i < psf_cubic_fit_intensity.size(); i++){
			psf_cubic_fit_intensity[i] = psf_cubic_fit_intensity[i]/max_intensity;
		}
}

void perform_fft(vector<float> &psf_cubic_fit_intensity, complex<double> *x_out){
	  // FFT of PSF
	  unsigned int size =  psf_cubic_fit_intensity.size();
	  float* psf_intensity = new float[size]();

	  for(unsigned int i =0; i < size; i++){
		  psf_intensity[i] = psf_cubic_fit_intensity[i];
		  //cout << "PSF Intensity : " << psf_intensity[i] << endl;
	  }
	  fft(psf_intensity,x_out,size);
}

void print_complex_array(complex<double> *x_out, unsigned int size){
	for(unsigned int i=0; i < size; i++){
		cout << x_out[i] << endl;
	}

}

void obtain_magnitude_of_fft_values(complex<double> *x_out, unsigned int size, vector<double> &fft_magnitude){
	for(unsigned int i = 0; i < size; i++){
		fft_magnitude.push_back(abs(x_out[i]));
	}
}

void normalize_fft(vector<double> &fft_magnitude){
	double max_intensity = FLT_MIN;
	unsigned int i = 0;
	for(i = 0; i < fft_magnitude.size(); i++){
		if(fft_magnitude[i] > max_intensity){
			max_intensity = fft_magnitude[i];
		}
	}

	for(i = 0; i < fft_magnitude.size(); i++){
		fft_magnitude[i] = fft_magnitude[i]/max_intensity;
	}
}

void sort_fft_value(vector<double> &fft_magnitude){
	double temp_fft_value = 0;
	for(unsigned int i = 0; i < fft_magnitude.size() ; i++){
		for(unsigned int j = 0; j < fft_magnitude.size()-1 ; j++){
			if(fft_magnitude[j] < fft_magnitude[j+1]){
				temp_fft_value = fft_magnitude[j+1];
				fft_magnitude[j+1] = fft_magnitude[j];
				fft_magnitude[j] = temp_fft_value;
			}
		}
	}
}

void write_fft_to_excel(vector<double> &fft_magnitude, float frequency){
	ofstream ex_file;
	ex_file.open("Values.csv");
	float interval = 0;
	for (unsigned int i=0; i < fft_magnitude.size(); i++){
		ex_file << interval << "," << fft_magnitude[i] << endl;
		interval += frequency;
	}

	ex_file.close();
}
