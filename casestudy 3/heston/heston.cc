#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include<complex>
#include<gsl/gsl_multimin.h>
#include<algorithm>

#define length 100
#define PI     3.1415926535897932

using namespace std;
typedef complex <double> dcmplx;

// fft 
dcmplx* fft ( dcmplx*, int );
dcmplx* ifft ( dcmplx*, int );
dcmplx* fft_get_twiddle_factors ( int );
void    fft_calculate ( dcmplx*, int, int, dcmplx*, dcmplx*, dcmplx* );

// funtions I write

struct  Option        { double S; double K; double r; double T; double q; double bid; double ask; double S0; double origin; char type; };
void    die           ( const char *message );
void    datareading   ( struct Option *option );
dcmplx  hestoncf      ( dcmplx mu, double S0, double r, double q, double T, double sigma, double kappa, double theta, double rho, double v0 );
double  fftheston     ( struct Option option, double sigma, double kappa, double theta, double rho, double v0 );
double  myheston      ( const gsl_vector *v, void *params ); 
double  ffthestoncall ( double S0, double K, double T, double r, double q, double sigma, double kappa, double theta, double rho, double v0 );
void    hestonsurface ( double S0, double r, double q, double sigma, double kappa, double theta, double rho, double v0 );

//---------------------------------------------------------------------------
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

	fft_calculate ( x, N / 2, skip * 2, E, bigx, twiddles );
	fft_calculate ( x + skip, N / 2, skip * 2, D, bigx, twiddles );

	int k;
	for ( k = 0; k < N / 2; k++ )
	{
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

//-----------------------------------------------------------------------------------------
void die ( const char *message )
{
	perror ( message );
	exit ( 1 );
}

void datareading ( struct Option *option  )
{	
	char buffer [ 20 ];
	int i;
	int j;

	FILE *data;
	//---------- price
	data = fopen ( "price.txt", "r" );
	if( data == NULL )
		die ( "failed to open price.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].S = atof ( buffer );
		i++;
	}
	fclose ( data );

	//---------- r
	data = fopen ( "r.txt", "r" );
	if ( data == NULL )
		die ( "failed to open r.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].r = atof ( buffer );
		i++;
	}
	fclose ( data );

	//---------- q
	data = fopen ( "q.txt", "r" );
	if( data == NULL )
		die ( "failed to open q.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].q = atof ( buffer );
		i++;
	}
	fclose ( data );

	//---------- tau
	data = fopen ( "maturity.txt", "r" );
	if( data == NULL )
		die ( "failed to open maturity.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].T = atof ( buffer );
		i++;
	}
	fclose ( data );

	//B
	//---------- strike
	data = fopen ( "strike.txt", "r" );
	if( data == NULL )
		die ( "failed to open price.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].K = atof ( buffer );
		i++;
	}
	fclose ( data );

	//---------- bid
	data = fopen ( "bid.txt", "r" );
	if( data == NULL )
		die ( "failed to open price.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].bid = atof ( buffer );
		i++;
	}
	fclose ( data );

	//---------- ask
	data = fopen ( "ask.txt", "r" );
	if( data == NULL )
		die ( "failed to open price.txt" );

	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		option [ i ].ask = atof ( buffer );
		i++;
	}
	fclose ( data );
}

dcmplx hestoncf ( dcmplx mu, double S0, double r, double q, double T, double sigma, double kappa, double theta, double rho, double v0 )
{
	dcmplx factor = kappa - rho * sigma * mu * dcmplx ( 0.0, 1.0 ); 
	dcmplx gamma  = sqrt ( sigma * sigma * ( mu * mu + mu * dcmplx ( 0.0, 1.0 ) ) + factor * factor );
	dcmplx upper  = exp  ( dcmplx ( 0.0, 1.0 ) * mu * ( log ( S0 ) + ( r - q ) * T ) + kappa * theta * T * factor / sigma / sigma );
	dcmplx middle = exp  ( - ( mu * mu + mu * dcmplx ( 0.0, 1.0 ) ) * v0 / ( gamma * cosh ( gamma * T / 2.0 ) / sinh ( gamma * T / 2.0 ) + factor ) ) ;
	dcmplx lower  = pow  ( cosh ( gamma * T / 2.0 ) + factor / gamma * sinh ( gamma * T / 2.0 ), 2.0 * kappa * theta / sigma / sigma );
	return upper * middle / lower;
}

double fftheston ( struct Option option, double sigma, double kappa, double theta, double rho, double v0 )
{
	double eta    = 0.25;
	int    N      = pow ( 2, 8 );
	double C      = exp ( -option.r * option.T );
	double lambda = 2 * PI / N / eta;

	double S0 = option.S0;

	double alpha;
	if ( option.type == 'c' )
	{
		alpha = 1.5;
	}
	else
	{
		alpha = -1.5;
	}

	double nu [ N ];
	dcmplx *input = ( dcmplx * ) malloc ( N * sizeof ( dcmplx ) );

	int i;
	for ( i = 0; i < N; i++ )
	{
		nu [ i ] = i * eta;
		input [ i ] = eta * C / dcmplx ( alpha, nu [ i ] ) / dcmplx ( alpha + 1.0, nu [ i ] ) 
			* exp ( dcmplx ( 0, -1.0 ) * ( log ( option.K ) - PI / eta ) * nu [ i ] ) 
			* hestoncf ( dcmplx ( nu [ i ], - alpha - 1.0 ), S0, option.r, option.q, option.T, sigma, kappa, theta, rho, v0 );	
	}
	input [ 0 ] = input [ 0 ] / 2.0; // according to the FFT method, the cofficient of first element is 0.5.

	dcmplx *output = fft ( input, N );

	double value = output [ N / 2 ].real () * exp ( -alpha * ( log ( option.K ) - ( N / 2 - N / 2 ) * 2 * PI / N / eta ) ) / PI;

	free ( input  );
	free ( output );

	return value;
}

double myheston ( const gsl_vector *v, void *params ) 
{	
	double sigma = gsl_vector_get( v, 0 ); 
	double kappa = gsl_vector_get( v, 1 ); 
	double theta = gsl_vector_get( v, 2 ); 
	double rho   = gsl_vector_get( v, 3 ); 
	double v0    = gsl_vector_get( v, 4 ); 

	struct Option *option = ( struct Option * ) params;// trucnt for our use the data.actually, the params is a void point to te data structure

	double *s  = ( double * ) malloc ( length * sizeof ( double ) );

	int i, j;
	for ( j = 0; j < length; j++ )
	{
		if ( ( ( option [ j ].type == 'c' ) && ( option [ j ].S0 < option [ j ].K ) ) || ( ( option [ j ].type == 'p' ) && ( option [ j ].S0 > option [ j ].K ) ) )
		{

			s [ j ] = fftheston ( option [ j ], sigma, kappa, theta, rho, v0 );		

		}
		else
			s [ j ] = option [ j ].S;
	}

	double object = 0.0;
	for ( j = 0; j < length; j++ )
	{
		if ( 2 * kappa * theta <  sigma * sigma )
		{
/* spread method *///	object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] ) / abs ( option [ j ].ask - option [ j ].bid ) 
/* penalty       *///	+ ( s [ j ] - option [ j ].origin ) * ( s [ j ] - option [ j ].origin ) / abs ( option [ j ].ask - option [ j ].bid );
/* equal  method */	object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] ); 
/* penalty       *///	+ ( s [ j ] - option [ j ].origin ) * ( s [ j ] - option [ j ].origin );
		}
		else 
		{ 
/* spread method *///	object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] ) / abs ( option [ j ].ask - option [ j ].bid ); 
/* penalty       *///	object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] ) + ( s [ j ] - option [ j ].origin ) * ( s [ j ] - option [ j ].origin );
/* equal  method */	object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] );

		}
/* kappa theta   *///	object = object + ( max ( 0.0, -theta ) + max ( 0.0, -kappa ) ) * 10000000000.0;	
	}
	free ( s );
	return object / 50;

} 

double ffthestoncall ( double S0, double K, double T, double r, double q, double sigma, double kappa, double theta, double rho, double v0 )
{
	double alpha  = 1.5;
	double eta    = 0.25;
	int    N      = pow ( 2, 8 );
	double C      = exp ( - r * T );
	double lambda = 2 * PI / N / eta;
	double nu [ N ];

	dcmplx *input = ( dcmplx * ) malloc ( N * sizeof ( dcmplx ) );

	int i;
	for ( i = 0; i < N; i++ )
	{
		nu    [ i ] = i * eta;
		input [ i ] = eta * C / dcmplx ( alpha, nu [ i ] ) / dcmplx ( alpha + 1.0, nu [ i ] ) 
			* exp ( dcmplx ( 0, -1.0 ) * ( log ( K ) - PI / eta ) * nu [ i ] ) 
			* hestoncf ( dcmplx ( nu [ i ], - alpha - 1.0 ), S0, r, q, T, sigma,  kappa, theta, rho, v0 );	
	}
	input [ 0 ] = input [ 0 ] / 2.0; //the first coeficient of input is 0.5.

	dcmplx *output = fft ( input, N );

	double value = output [ N / 2 ].real () * exp ( -alpha * ( log ( K ) - ( N / 2 - N / 2 ) * 2 * PI / N / eta ) ) / PI;

	free ( input  );
	free ( output );

	return value;
}

void hestonsurface ( double S0, double r, double q, double sigma, double kappa, double theta, double rho, double v0 )
{
	int M = 100;
	int N = 100;

	double kmin   = 1800.0;
	double tmin   = 0;
	double deltak = 100.0 / M;
	double deltat = 1.0 / N;	// the numeritor has to be 1.0, because deltat is double!!!!!!!!!!!!

	double  *k = ( double *  ) malloc ( ( M + 1 ) * sizeof ( double ) );
	double  *t = ( double *  ) malloc ( ( N + 1 ) * sizeof ( double ) );

	int i, j;
	for ( i = 0; i < ( M + 1 ); i++ )
	{
		k [ i ] = kmin + i * deltak;
	}
	for ( i = 0; i < ( N + 1 ); i++ )
	{
		t [ i ] = tmin + i * deltat;
	}

	double **c = ( double ** ) malloc ( ( M + 1 ) * sizeof ( double * ) );
	for ( i = 0; i < M + 1; i++ )
	{
		c [ i ] = ( double * ) malloc ( ( N + 1 ) * sizeof ( double ) );
	}

	for ( i = 0; i < M + 1; i++ )
	{
		for ( j = 0; j < N + 1; j++ )
		{
			c [ i ] [ j ] = ffthestoncall ( S0, k [ i ], t [ j ], r, q, sigma, kappa, theta, rho, v0 );
			if ( c [ i ] [ j ] < 0 )
			{
				c [ i ] [ j ] = 0;
			}
		}
	}	

	double **vol = ( double ** ) malloc ( ( M + 1 ) * sizeof ( double * ) );
	for ( i = 0; i < M + 1; i++ )
	{
		vol [ i ] = ( double * ) malloc ( ( N + 1 ) * sizeof ( double ) );
	}

	for ( i = 1; i < M; i++ )
	{
		for ( j = 1; j < N; j ++ )
		{
			vol [ i ] [ j ] = ( c [ i ] [ j + 1 ] - c [ i ] [ j - 1 ] ) / 2 / deltat + q * c [ i ] [ j ];
			vol [ i ] [ j ] = vol [ i ] [ j ] + k [ i ] * ( r - q ) * ( c [ i + 1 ] [ j ] - c [ i - 1 ] [ j ] ) / 2 / deltak;
			vol [ i ] [ j ] = vol [ i ] [ j ] * 2 * deltak * deltak / k [ i ] / k [ i ];
			vol [ i ] [ j ] = vol [ i ] [ j ] / ( c [ i + 1 ] [ j ] - 2 * c [ i ] [ j ] + c [ i - 1 ] [ j ] );
			if ( vol [ i ] [ j ] < 0 )
				vol [ i ] [ j ] = 0;
			else
				vol [ i ] [ j ] = sqrt ( vol [ i ] [ j ] );
		}
	}

	for ( i = 0; i < M + 1; i++ )
	{
		vol [ i ] [ 0 ] = vol [ i ] [ 1 ];
		vol [ i ] [ N ] = vol [ i ] [ N - 1 ];
	}
	for ( j = 0; j < N + 1; j++ )
	{
		vol [ 0 ] [ j ] = vol [ 1 ] [ j ];
		vol [ M ] [ j ] = vol [ M - 1 ] [ j ];
	}

	cout << endl;

	for ( i = 0; i < M + 1; i++ )
	{
		for ( j = 0; j < N + 1; j++ )
		{
		}
	}


	FILE *data = fopen ( "surface.txt", "wb" );
	if( data == NULL )
		die ( "failed to open price.txt" );
	for ( i = 5; i < M - 4; i++ )
	{
		for ( j = 2; j < N - 1; j++ )
		{
			char buffer [ 20 ];
			sprintf ( buffer, "%f", vol [ i ] [ j ] );
			fprintf ( data, "%s\t", buffer );
		}
		fprintf ( data, "\n" );
	}
	fclose ( data );

	free ( k );
	free ( t );

	for ( i = 0; i < M + 1; i++ )
	{
		free ( c   [ i ] );
		free ( vol [ i ] );
	}
	free ( c );
	free ( vol );
	return;
}
//--------------------------------------------------------------------------------------------
int main ( int argc, char **argv ) 
{

	struct Option *option = ( struct Option * ) malloc ( length * sizeof ( struct Option ) );

	datareading ( option );

	double sigma = 1.0143;
	double kappa = 4.9549;
	double theta = 0.0562;
	double rho   = -0.6552;
	double v0    = 0.057;
	double S0    = 1861.34;	
/*
	double sigma = 1.4205;
	double kappa = 4.5408;
	double theta = 0.0454;
	double rho   = -0.7363;
	double v0    = 0.0235;
	double S0    = 1861.34;
*//*
	double sigma = 0.583266;
	double kappa = 3.224887;
	double theta = 0.038943;
	double rho   = -0.869719;
	double v0    = 0.017816;
	double S0    = 1861.34;
*/
/*	//	v0 = 0.0235	kappa = 4.5408	theta = 0.0454	sigma = 1.4205	rho = -0.7363 
	double sigma = 2.0;
	double kappa = 4.5;
	double theta = 0.015;
	double rho   = -1.0;
	double v0    = 0.015;
	double S0    = 1861.34;
*/

	int j;
	for ( j = 0; j < length; j++ )	
	{
		if ( j < ( length / 2 ) )
		{
			option [ j ].type = 'c';
		}
		else
		{
			option [ j ].type = 'p';
		}
		option [ j ].S0     = S0;

		if ( ( ( option [ j ].type == 'c' ) && ( S0 < option [ j ].K ) ) || ( ( option [ j ].type == 'p' ) && ( S0 > option [ j ].K ) ) )
		{

			option [ j ].origin = fftheston ( option [ j ], sigma, kappa, theta, rho, v0 );		

		}
		else
			option [ j ].origin = option [ j ].S;

		//		printf ( "%d\t%f\t%f\t%f\n", j + 1, option [ j ].origin, option [ j ].origin - option [ j ].S, option [ j ].S );
	}



	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; 

	gsl_multimin_fminimizer *s = NULL; 
	gsl_vector *ss, *x; 
	gsl_multimin_function minex_func; 

	int iter = 0;
	int status; 
	double size; 

	/* Starting point */ 
	x = gsl_vector_alloc ( 5 ); 
	gsl_vector_set ( x, 0, sigma ); 
	gsl_vector_set ( x, 1, kappa ); 
	gsl_vector_set ( x, 2, theta ); 
	gsl_vector_set ( x, 3, rho   ); 
	gsl_vector_set ( x, 4, v0    ); 

	/* Set initial step size s to 1 */ 
	ss = gsl_vector_alloc ( 5 ); 
	gsl_vector_set_all ( ss, 1.0 ); 

	/* Initialize method and iterate */ 
	minex_func.n = 5; 
	minex_func.f = myheston; 
	minex_func.params = option; 
	s = gsl_multimin_fminimizer_alloc ( T, 5 ); 
	gsl_multimin_fminimizer_set ( s, &minex_func, x, ss ); 

	do 
	{ 
		iter++; 
		status = gsl_multimin_fminimizer_iterate ( s ); 

		if ( status ) 
			break; 

		size = gsl_multimin_fminimizer_size ( s ); 
		status = gsl_multimin_test_size ( size, 1e-4 ); 

		if ( status == GSL_SUCCESS ) 
		{ 
			printf ( "converged to minimum at\n" ); 
		} 
		printf ("%5d %f\t%f\t%f\t%f\t%f f() = %7.3f size = %f\n", 
				iter, 
				gsl_vector_get ( s -> x, 0 ), 
				gsl_vector_get ( s -> x, 1 ),
				gsl_vector_get ( s -> x, 2 ),
				gsl_vector_get ( s -> x, 3 ),
				gsl_vector_get ( s -> x, 4 ),
				s -> fval, size ); 
	} 
	while ( status == GSL_CONTINUE && iter < 5000 );

	// for evaluate the results
	sigma = gsl_vector_get ( s -> x, 0 );
	kappa = gsl_vector_get ( s -> x, 1 );
	theta = gsl_vector_get ( s -> x, 2 );
	rho   = gsl_vector_get ( s -> x, 3 );
	v0    = gsl_vector_get ( s -> x, 4 );

	double *price  = ( double * ) malloc ( length * sizeof ( double ) );

	for ( j = 0; j < length; j++ )
	{

		price [ j ] = fftheston ( option [ j ], sigma, kappa, theta, rho, v0 );		

		printf ( "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", j + 1, price [ j ], option [ j ].S - price [ j ],
				option [ j ].S, option [ j ].K, option [ j ].r, option [ j ].q, option [ j ].T );
	}
	
	printf ( "%f\n", 2 * kappa * theta - sigma * sigma );

	double object = 0.0;
	for ( j = 0; j < length; j++ )
	{	
		// if we use different objective function, when we compare the results, we should use the same measure!!!!!!!!!
		object = object + ( option [ j ].S - price [ j ] ) * ( option [ j ].S - price [ j ] );
	}

	printf ( "%f\t%f\n", object / 100, S0 );

	double r = 0.0;
	double q = 0.0;
	for ( j = 0; j < length; j++ )
	{
		r = r + option [ j ].r;
		q = q + option [ j ].q;
	}
	r = r / length;
	q = q / length;

	hestonsurface ( S0, r, q, sigma, kappa, theta, rho, v0 );

	gsl_vector_free( x  ); 
	gsl_vector_free( ss ); 
	gsl_multimin_fminimizer_free ( s ); 

	free ( option );
	free ( price  );
	return status; 
} 


