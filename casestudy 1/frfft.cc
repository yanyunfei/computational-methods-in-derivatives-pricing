#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <complex>

using namespace std;

typedef complex <double> dcmplx;

#define PI 3.1415926535897932

//
//FFT functions
//
dcmplx* fft ( dcmplx*, int );
dcmplx* ifft ( dcmplx*, int );
void fft_calculate ( dcmplx*, int, int, dcmplx*, dcmplx*, dcmplx* );
dcmplx* fft_get_twiddle_factors ( int );

//
//N mustbe power of 2
//
dcmplx* fft ( dcmplx* x, int N )
{
	dcmplx* out      = ( dcmplx* ) malloc ( sizeof ( dcmplx ) * N );
	dcmplx* scratch  = ( dcmplx* ) malloc ( sizeof ( dcmplx ) * N );
	dcmplx* twiddles = fft_get_twiddle_factors ( N );

	fft_calculate ( x, N, 1, out, scratch, twiddles );

	free ( twiddles );
	free ( scratch );
	return out;
}



dcmplx* ifft ( dcmplx* x, int N )
{
	dcmplx* out      = ( dcmplx* ) malloc ( sizeof ( dcmplx ) * N );
	dcmplx* scratch  = ( dcmplx* ) malloc ( sizeof ( dcmplx ) * N );
	dcmplx* twiddles = fft_get_twiddle_factors ( N );

	fft_calculate ( x, N, 1, out, scratch, twiddles );

	free ( twiddles );
	free ( scratch );
	
	//claculate ifft vi areciprocity property of dft
	
	int N2 = N / 2;
	double temp0, temp1;
	
	out [ 0 ]  = dcmplx ( real ( out [ 0 ] ) / N, imag ( out [ 0 ] ) / N );
	out [ N2 ] = dcmplx ( real ( out [ N2 ] ) / N, imag ( out [ N2 ] ) / N );

	int i;
	for ( i = 1; i < N2; i++ )
	{
		temp0 = real ( out [ i ] ) / N;
		temp1 = imag ( out [ i ] ) / N;

		out [ i ]  = dcmplx ( real ( out [ N - i ] ) / N, imag ( out [ N - i ] ) / N );
		out [ N - i ] = dcmplx ( temp0, temp1 );
	}

	return out;
}



void fft_calculate ( dcmplx* x, int N, int skip, dcmplx* bigx, dcmplx* D, dcmplx* twiddles )
{
	dcmplx* E = D + N / 2;

	if ( N == 1 )
	{
		bigx [ 0 ] = x [ 0 ];
		return;
	}

	//use x as a scratch buffer
	fft_calculate ( x, N / 2, skip * 2, E, bigx, twiddles );
	fft_calculate ( x + skip, N / 2, skip * 2, D, bigx, twiddles );

	int k;
	for ( k = 0; k < N / 2; k++ )
	{
		//mutiply entries of D by the twiddle factors e ^ ( -2 * PI * i / N * K )
		D [ k ] = twiddles [ k * skip ] * D [ k ];
	}

	for ( k = 0; k < N / 2; ++k )
	{
		bigx [ k ]         = E [ k ] + D [ k ];
		bigx [ k + N / 2 ] = E [ k ] - D [ k ];
	}
}



dcmplx* fft_get_twiddle_factors ( int N )
{
	dcmplx* twiddles = ( dcmplx* ) malloc ( sizeof ( dcmplx ) * N / 2 );
	
	int k;
	for ( k = 0; k != N / 2; ++k )
	{
		twiddles [ k ] = exp ( dcmplx ( 0, -2.0 * PI * k / N ) );
	}
	return twiddles;
}



int main ( int, char** )
{
	double r = 0.0025;
	double S0 = 1800;
	double q = 0.0203;
	double sigma = 0.3;
	double T = 0.5;
	double K = 900;

	double eta = 0.25;
	double alpha = -2;
	int    N = pow ( 2, 8 );
	double C = exp ( -r * T );
	double lambda = 0.1;
	double gamma = lambda * eta / 2 / PI;

	double nu [ N ];
	
	dcmplx input [ N ];
// get xi
	int i;
	for ( i = 0; i < N; i++ )
	{
		nu [ i ] = i * eta;
		if ( i == 0 )
		{
			input [ i ] =     eta / 2 * C / dcmplx ( alpha, nu [ i ] ) / dcmplx ( alpha + 1, nu [ i ] ) 
					* exp ( dcmplx ( 0, -1 ) * ( log ( K ) - lambda * N / 2 ) * nu [ i ] ) 
					* exp ( ( log ( S0 ) + ( r - q - 0.5 * sigma * sigma ) * T ) * dcmplx ( nu [ i ], - alpha - 1 ) 
					* dcmplx ( 0, 1 ) - 0.5 * sigma * sigma * T * dcmplx ( nu [ i ], -alpha - 1 ) * dcmplx ( nu [ i ], -alpha - 1 ) ); 
		}		
		else
			input [ i ] =     eta * C / dcmplx ( alpha, nu [ i ] ) / dcmplx ( alpha + 1, nu [ i ] ) 
					* exp ( dcmplx ( 0, -1 ) * ( log ( K ) - lambda * N / 2 ) * nu [ i ] ) 
					* exp ( ( log ( S0 ) + ( r - q - 0.5 * sigma * sigma ) * T ) * dcmplx ( nu [ i ], - alpha - 1 ) 
					* dcmplx ( 0, 1 ) - 0.5 * sigma * sigma * T * dcmplx ( nu [ i ], -alpha - 1 ) * dcmplx ( nu [ i ], -alpha - 1 ) );          
	}

	dcmplx y [ 2 * N ];
	dcmplx z [ 2 * N ];
// get yi
	for ( i = 0; i < ( 2 * N ); i++ )
	{
		if ( i < N )
		{
			y [ i ] = input [ i ] * exp ( -dcmplx ( 0, 1 ) * PI * gamma * pow ( i, 2 ) );
		}
		else
			y [ i ] = 0;
	}
// get zi
	for ( i = 0; i < ( 2 * N ); i++ )
	{
		if ( i < N )
		{
			z [ i ] = exp ( dcmplx ( 0, 1 ) * gamma * PI * pow ( i, 2 ) );
		}
		else
			z [ i ] = z [ 2 * N - 1 - i ];
	}
// use fft to get the ouput of yi and zi
	dcmplx* yhat = fft ( y, 2 * N );
	dcmplx* zhat = fft ( z, 2 * N );
// calculate the xi
	dcmplx xi [ 2 * N ];
	for ( i = 0; i < ( 2 * N ); i++ )
	{
		xi [ i ] = yhat [ i ] * zhat [ i ];
	}

	free ( yhat );
	free ( zhat );

	dcmplx* rawprice = ifft ( xi, 2 * N );

	printf ( "FRFFT\tK\n" );
	
	for ( i = N / 2 - 7; i < N / 2 + 8; i ++ )
	{
		printf ( "%f, %f\n",    ( rawprice [ i ] * exp ( -dcmplx ( 0, 1 ) * PI * gamma * pow ( i , 2 ) ) ).real () 
					* exp ( -alpha * ( log ( K ) - ( N / 2 - i ) * lambda ) ) / PI, exp ( log ( K ) - lambda * N / 2 + i * lambda ) );
	}

	free ( rawprice );

	return 0;
}

	
