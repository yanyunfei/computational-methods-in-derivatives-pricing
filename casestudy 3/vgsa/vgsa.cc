


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
struct Option        { double S; double K; double r; double T; double q; double bid; double ask; double S0; double origin; char type; };
void   die           ( const char *message );
void   datareading   ( struct Option *option );
dcmplx vgcf          ( dcmplx mu, double nu, double theta, double sigma );
dcmplx circf         ( dcmplx mu, double t, double y0, double kappa, double eta, double lambda );
dcmplx vgsacf        ( dcmplx mu, double S0, double y0, double r, double q, double t, double sigma, double theta, double kappa, double eta, double lambda, double nu );
double fftvgsa       ( struct Option option, double sigma, double theta, double kappa, double Eta, double Lambda, double Nu, double y0 );
double myvgsa        ( const gsl_vector *v, void *params ); 
double fftvgsacall   ( double S0, double K, double T, double r, double q, double sigma, double theta, double kappa, double eta, double lambda, double nu, double y0 );
void   vgsasurface   ( double S0, double r, double q, double sigma, double theta, double kappa, double eta, double lambda, double nu, double y0 );

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

dcmplx vgcf ( dcmplx mu, double nu, double theta, double sigma )
{
	return - log ( 1.0 - dcmplx ( 0.0, 1.0 ) * mu * theta * nu + sigma * sigma * nu * mu * mu / 2.0 ) / nu;
}

dcmplx circf ( dcmplx mu, double t, double y0, double kappa, double eta, double lambda )
{
	dcmplx gamma  = sqrt ( kappa * kappa - 2.0 * lambda * lambda * mu * dcmplx ( 0.0, 1.0 ) );
	dcmplx factor = gamma * t / 2.0;
	dcmplx upper  = exp ( kappa * kappa * t * eta / lambda / lambda ); 	
	dcmplx middle = exp ( 2.0 * mu * dcmplx ( 0.0, 1.0 ) * y0 / ( kappa + gamma * cosh ( factor ) / sinh ( factor ) ) );
	dcmplx lower  = pow ( cosh ( factor ) + kappa / gamma * sinh ( factor ) , 2.0 * kappa * eta / lambda / lambda );
	return upper * middle / lower;
}

dcmplx vgsacf ( dcmplx mu, double S0, double y0, double r, double q, double t, double sigma, double theta, double kappa, double eta, double lambda, double nu )
{
	dcmplx middle = exp ( mu * dcmplx ( 0.0, 1.0 ) * ( log ( S0 ) + ( r - q ) * t ) );
	dcmplx upper  = circf ( - dcmplx ( 0.0, 1.0 ) * vgcf ( mu, nu, theta, sigma ), t, 1.0 / nu, kappa, eta, lambda );
	dcmplx lower  = pow ( circf ( - dcmplx ( 0.0, 1.0 ) * vgcf ( dcmplx ( 0.0, - 1.0 ), nu, theta, sigma ), t, 1.0 / nu, kappa, eta, lambda ), dcmplx ( 0.0, 1.0 ) * mu );
	return middle * upper / lower;
}

double fftvgsa ( struct Option option, double sigma, double theta, double kappa, double Eta, double Lambda, double Nu, double y0 )
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
			* vgsacf ( dcmplx ( nu [ i ], - alpha - 1.0 ), S0, y0, option.r, option.q, option.T, sigma, theta, kappa, Eta, Lambda, Nu ); 
	}
	input [ 0 ] = input [ 0 ] / 2.0; //the first coeficient of xi is 0.5, that's why I have if here. this step is calculating the xi

	dcmplx *output = fft ( input, N );

	double value = output [ N / 2 ].real () * exp ( -alpha * ( log ( option.K ) - ( N / 2 - N / 2 ) * 2 * PI / N / eta ) ) / PI;

	free ( input  );
	free ( output );

	return value;
}

double myvgsa ( const gsl_vector *v, void *params ) 
{	
	double sigma = gsl_vector_get ( v, 0 ); 
	double theta = gsl_vector_get ( v, 1 ); 
	double kappa = gsl_vector_get ( v, 2 ); 
	double eta   = gsl_vector_get ( v, 3 ); 
	double lambda= gsl_vector_get ( v, 4 ); 
	double nu    = gsl_vector_get ( v, 5 ); 
	double y0 = 1.0 / nu;	

	struct Option *option = ( struct Option * ) params;// trucnt for our use the data.actually, the params is a void point to te data structure

	double *s  = ( double * ) malloc ( length * sizeof ( double ) );

	int j;
	for ( j = 0; j < length; j++ )
	{
		if ( ( ( option [ j ].type == 'c' ) && ( option [ j ].S0 < option [ j ].K ) ) || ( ( option [ j ].type == 'p' ) && ( option [ j ].S0 > option [ j ].K ) ) )
		{

			s [ j ] = fftvgsa ( option [ j ], sigma, theta, kappa, eta, lambda, nu, y0 );		

		}
		else
			s [ j ] = option [ j ].S;
	}

	double object = 0.0;
	for ( j = 0; j < length; j++ )
	{
//		object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] ) / abs ( option [ j ].ask - option [ j ].bid );
		object = object + ( option [ j ].S - s [ j ] ) * ( option [ j ].S - s [ j ] );
	}

//	object = object + max ( 0.0, theta ) * 10000000000.0;	
	free ( s );
	return object / 50;

} 

double  fftvgsacall   ( double S0, double K, double T, double r, double q, double sigma, double theta, double kappa, double Eta, double Lambda, double Nu, double y0 )
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
			* vgsacf ( dcmplx ( nu [ i ], - alpha - 1.0 ), S0, y0, r, q, T, sigma, theta, kappa, Eta, Lambda, Nu ); 
	}
	input [ 0 ] = input [ 0 ] / 2.0; //the first coeficient of xi is 0.5, that's why I have if here. this step is calculating the xi

	dcmplx *output = fft ( input, N );

	double value = output [ N / 2 ].real () * exp ( -alpha * ( log ( K ) - ( N / 2 - N / 2 ) * 2 * PI / N / eta ) ) / PI;

	free ( input  );
	free ( output );
	return value;
}

void vgsasurface ( double S0, double r, double q, double sigma, double theta, double kappa, double eta, double lambda, double nu, double y0 )
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
			c [ i ] [ j ] = fftvgsacall ( S0, k [ i ], t [ j ], r, q, sigma, theta, kappa, eta, lambda, nu, y0 );
			if ( c [ i ] [ j ] < 0 )	
				c [ i ] [ j ] = 0;
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
	for ( i = 0; i < M + 1; i++ )
	{
		for ( j = 0 ; j < N + 1; j++ )
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

	double sigma = 0.1022;
	double theta = -0.0761;
	double kappa = 8.1143;
	double eta   = 2.8060;
	double lambda= 10.3646;
	double nu    = 0.1819;
	double y0    = 1 / nu;
	double S0    = 1861.34;
/*
	double sigma = 0.066513;
	double theta = -0.054227;
	double kappa = 6.273433;
	double eta   = 3.014308;
	double lambda= 10.575303;
	double nu    = 0.237939;
	double y0    = 1 / nu;
	double S0    = 1861.34;
*
	double sigma = 0.0824;
	double theta = -0.0447;
	double kappa = 8.7843;
	double eta   = 2.6429;
	double lambda= 17.3831;
	double nu    = 0.2480;
	double y0    = 1 / nu;
	double S0    = 1861.34;

	double sigma = 0.1;
	double theta = -0.045;
	double kappa = 7.0;
	double eta   = 2.8;
	double lambda= 20;
	double nu    = 0.2;
	double y0    = 1 / nu;
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
	
		option [ j ].S0 = S0;
		
		if ( ( ( option [ j ].type == 'c' ) && ( S0 < option [ j ].K ) ) || ( ( option [ j ].type == 'p' ) && ( S0 > option [ j ].K ) ) )
		{
			option [ j ].origin = fftvgsa ( option [ j ], sigma, theta, kappa, eta, lambda, nu, y0 );		
		}
		else
		{	
			option [ j ].origin = option [ j ].S;
		}
		
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
	x = gsl_vector_alloc ( 6 ); 
	gsl_vector_set ( x, 0, sigma ); 
	gsl_vector_set ( x, 1, theta ); 
	gsl_vector_set ( x, 2, kappa ); 
	gsl_vector_set ( x, 3, eta   ); 
	gsl_vector_set ( x, 4, lambda); 
	gsl_vector_set ( x, 5, nu    );

	/* Set initial step size s to 1 */ 
	ss = gsl_vector_alloc ( 6 ); 
	gsl_vector_set_all ( ss, 3.0 ); 

	/* Initialize method and iterate */ 
	minex_func.n = 6; 
	minex_func.f = myvgsa; 
	minex_func.params = option; 
	s = gsl_multimin_fminimizer_alloc ( T, 6 ); 
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
		printf ("%5d %f\t%f\t%f\t%f\t%f\t%f f() = %7.3f size = %f\n", 
				iter, 
				gsl_vector_get ( s -> x, 0 ), 
				gsl_vector_get ( s -> x, 1 ),
				gsl_vector_get ( s -> x, 2 ),
				gsl_vector_get ( s -> x, 3 ),
				gsl_vector_get ( s -> x, 4 ),
				gsl_vector_get ( s -> x, 5 ),
				s -> fval, size ); 
	} 
	while ( status == GSL_CONTINUE && iter < 5000 );



	// for evaluate the results
	sigma = gsl_vector_get ( s -> x, 0 );
	theta = gsl_vector_get ( s -> x, 1 );
	kappa = gsl_vector_get ( s -> x, 2 );
	eta   = gsl_vector_get ( s -> x, 3 );
	lambda= gsl_vector_get ( s -> x, 4 );
	nu    = gsl_vector_get ( s -> x, 5 );
	y0    = 1 / nu;
	

	double *price  = ( double * ) malloc ( length * sizeof ( double ) );

	for ( j = 0; j < length; j++ )
	{
		price [ j ] = fftvgsa ( option [ j ], sigma, theta, kappa, eta, lambda, nu, y0 );
	
		printf ( "%d\t%10f\t%10f\t%10f\t%10f\t%10f\t%10f\t%10f\n", j + 1, price [ j ], option [ j ].S - price [ j ],
				option [ j ].S, option [ j ].K, option [ j ].r, option [ j ].q, option [ j ].T );
	}

	printf ( "%f\n", theta );
	
	double object = 0.0;
	for ( j = 0; j < length; j++ )
	{
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

	vgsasurface ( S0, r, q, sigma, theta, kappa, eta, lambda, nu, y0 );

	gsl_vector_free( x  ); 
	gsl_vector_free( ss ); 
	gsl_multimin_fminimizer_free ( s ); 

	free ( option );
	free ( price  );
	return status; 
} 
