#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<math.h>
#include<iostream>
#include<gsl/gsl_multimin.h>
#include<algorithm>
#include<ctime>

#define length 1589
#define PI     3.1415926535897932

using namespace std;

void   pricereading ( double *price );
double myparameter  ( const gsl_vector *vector, void *params );


void pricereading ( double *price )
{	
	char buffer [ 20 ];
	int i;
	int j;

	FILE *data;
	//---------- price
	data = fopen ( "price.txt", "r" );
	i = 0;
	while ( fgets ( buffer, sizeof ( buffer ), data ) != NULL )
	{
		price [ i ] = atof ( buffer );
		i++;
	}
	fclose ( data );
}


double myparameter ( const gsl_vector *vector, void *params ) 
{	
	// read the parameters
	double mu    = gsl_vector_get ( vector, 0 );
	double sigma = gsl_vector_get ( vector, 1 ); 
	double kappa = gsl_vector_get ( vector, 2 ); 
	double theta = gsl_vector_get ( vector, 3 ); 
	double rho   = gsl_vector_get ( vector, 4 ); 
	double v0    = gsl_vector_get ( vector, 5 ); 

	// read the price data
	double *y = ( double * ) params;

	// v is the post-v, v1 is the fore-v
	double v = v0;
	double v1;
	double likelihood = 0.0;
	double dt = 1.0 / 365.0;
	double A;
	double e;

	// p is the post-p, p1 is the fore-p
	double p  [ 4 ] = { 0.01, 0.0, 0.0, 0.01 };
	double p1 [ 4 ];
	double f  [ 4 ] = { 1, -0.5 * dt, 0, 1 - kappa * dt };
	double u  [ 4 ] = { sqrt ( v * dt ), 0, 0, sigma * sqrt ( v * dt ) };
	double q  [ 4 ] = { 1.0, rho, rho, 1.0 };

	// x is the post-x, x1 is the fore-x
	double h  [ 2 ] = { 1.0, 0 };
	double k  [ 2 ];
	double x  [ 2 ];
	double x1 [ 2 ];

	int i;
	for ( i = 1; i < length; i++ )
	{

		// get the model price vector
		x1 [ 0 ] = y [ i - 1 ] + ( mu - 0.5 * v ) * dt;
		x1 [ 1 ] = v + kappa * ( theta - v ) * dt;

		// refresh the elements in matrix u, it's depend ou our post-v, which will be modified at the end of the loop
		u [ 0 ] = sqrt ( v * dt );
		u [ 3 ] = sqrt ( v * dt ) * sigma;

		// compute the p1, p1 = F*P*F'+U*Q*U'
		p1 [ 0 ] = ( f [ 0 ] * p [ 0 ] + f [ 1 ] * p [ 2 ] ) * f [ 0 ] + ( f [ 0 ] * p [ 1 ] + f [ 1 ] * p [ 3 ] ) * f [ 1 ];
		p1 [ 1 ] = ( f [ 0 ] * p [ 0 ] + f [ 1 ] * p [ 2 ] ) * f [ 2 ] + ( f [ 0 ] * p [ 1 ] + f [ 1 ] * p [ 3 ] ) * f [ 3 ];
		p1 [ 2 ] = ( f [ 2 ] * p [ 0 ] + f [ 3 ] * p [ 2 ] ) * f [ 0 ] + ( f [ 2 ] * p [ 1 ] + f [ 3 ] * p [ 3 ] ) * f [ 1 ];
		p1 [ 3 ] = ( f [ 2 ] * p [ 0 ] + f [ 3 ] * p [ 2 ] ) * f [ 2 ] + ( f [ 2 ] * p [ 1 ] + f [ 3 ] * p [ 3 ] ) * f [ 3 ];

		p1 [ 0 ] = ( u [ 0 ] * q [ 0 ] + u [ 1 ] * q [ 2 ] ) * u [ 0 ] + ( u [ 0 ] * q [ 1 ] + u [ 1 ] * q [ 3 ] ) * u [ 1 ] + p1 [ 0 ];
		p1 [ 1 ] = ( u [ 0 ] * q [ 0 ] + u [ 1 ] * q [ 2 ] ) * u [ 2 ] + ( u [ 0 ] * q [ 1 ] + u [ 1 ] * q [ 3 ] ) * u [ 3 ] + p1 [ 1 ];
		p1 [ 2 ] = ( u [ 2 ] * q [ 0 ] + u [ 3 ] * q [ 2 ] ) * u [ 0 ] + ( u [ 2 ] * q [ 1 ] + u [ 3 ] * q [ 3 ] ) * u [ 1 ] + p1 [ 2 ];
		p1 [ 3 ] = ( u [ 2 ] * q [ 0 ] + u [ 3 ] * q [ 2 ] ) * u [ 2 ] + ( u [ 2 ] * q [ 1 ] + u [ 3 ] * q [ 3 ] ) * u [ 3 ] + p1 [ 3 ];

		// A = H*P1*H'+V*R*V',V = { 0, 0 }', which remove the second term.
		A = ( h [ 0 ] * p1 [ 0 ] + h [ 1 ] * p1 [ 2 ] ) * h [ 0 ] + ( h [ 0 ] * p1 [ 1 ] + h [ 1 ] * p1 [ 3 ] ) * h [ 1 ];

		// K = P1*H'/A
		k [ 0 ] = ( p1 [ 0 ] * h [ 0 ] + p1 [ 1 ] * h [ 1 ] ) / A;
		k [ 1 ] = ( p1 [ 2 ] * h [ 0 ] + p1 [ 3 ] * h [ 1 ] ) / A;

		// refresh post-v, which is the second element of x1
		v1 = x1 [ 1 ];
		// compute the error e = real price - model price = y [ i ] - y [ i ]',y [ i ]' = y [ i - 1 ] + ( mu - 0.5 * v1 ) * dt;
		e  = y [ i ] - y [ i - 1 ] - ( mu - 0.5 * v1 ) * dt;

		// eliminate the A == Nan;
		if( A == A )
		{
			likelihood += ( log ( A ) + e * e / A );
		}

		// add the penalty function, which means, theta, kappa, v0, sigma should not be nagative, also require that 2kappatheta>sigmasigma
		likelihood += ( max ( 0.0, -theta ) + max ( 0.0, -kappa ) + max ( 0.0, -v0 ) + max ( 0.0, -sigma ) 
				+ max ( 0.0, sigma * sigma - 2 * kappa * theta ) ) * 100000000; 
		
		// refresh the x
		x [ 0 ] = x1 [ 0 ] + k [ 0 ] * e;
		x [ 1 ] = x1 [ 1 ] + k [ 1 ] * e;

		// refresh v
		v = x [ 1 ];

		// post-p P = ( I - K * H ) * P1;
		p [ 0 ] = ( 1 - k [ 0 ] * h [ 0 ] ) * p1 [ 0 ] - k [ 0 ] * h [ 1 ] * p1 [ 2 ];
		p [ 1 ] = ( 1 - k [ 0 ] * h [ 0 ] ) * p1 [ 1 ] - k [ 0 ] * h [ 1 ] * p1 [ 3 ];
		p [ 2 ] = ( 1 - k [ 1 ] * h [ 1 ] ) * p1 [ 2 ] - k [ 1 ] * h [ 0 ] * p1 [ 0 ];
		p [ 3 ] = ( 1 - k [ 1 ] * h [ 1 ] ) * p1 [ 3 ] - k [ 1 ] * h [ 0 ] * p1 [ 1 ];

	}

	return likelihood;

} 


//--------------------------------------------------------------------------------------------
int main ( int argc, char **argv ) 
{

	double *y = ( double * ) malloc ( length * sizeof ( double ) );

	pricereading ( y );

	double mu    = 0.0125;
	double sigma = 0.3353;
	double kappa = 1.4440;
	double theta = 0.0184;
	double rho   = -0.5267;
	double v0    = 0.0225;

	/*	
		double mu    = 0.0125;
		double sigma = 0.3353;
		double kappa = 2.4440;
		double theta = 0.0124;
		double rho   = -0.5267;
		double v0    = 0.00001;

		double mu    = 0.011;
		double sigma = 0.4;
		double kappa = 1.0;
		double theta = 0.05;
		double rho   = -0.5;
		double v0    = 0.02;
		*/
	
	// record the time
	clock_t start, finish;
	start = clock ();

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; 

	gsl_multimin_fminimizer *s = NULL; 
	gsl_vector *ss, *x; 
	gsl_multimin_function minex_func; 

	int iter = 0;
	int status; 
	double size; 

	/* Starting point */ 
	x = gsl_vector_alloc ( 6 );
	gsl_vector_set ( x, 0, mu    ); 	
	gsl_vector_set ( x, 1, sigma ); 
	gsl_vector_set ( x, 2, kappa ); 
	gsl_vector_set ( x, 3, theta ); 
	gsl_vector_set ( x, 4, rho   ); 
	gsl_vector_set ( x, 5, v0    ); 

	/* Set initial step size s to 1 */ 
	ss = gsl_vector_alloc ( 6 ); 
	gsl_vector_set_all ( ss, 1.0 ); 

	/* Initialize method and iterate */ 
	minex_func.n = 6; 
	minex_func.f = myparameter; 
	minex_func.params = y; 
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
		printf ("%5d\t%f\t%f\t%f\t%f\t%f\t%f f() = %7.3f size = %f\n", 
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
	mu    =	gsl_vector_get ( s -> x, 0 );
	sigma = gsl_vector_get ( s -> x, 1 );
	kappa = gsl_vector_get ( s -> x, 2 );
	theta = gsl_vector_get ( s -> x, 3 );
	rho   = gsl_vector_get ( s -> x, 4 );
	v0    = gsl_vector_get ( s -> x, 5 );

	cout << endl << endl;
	cout << "Is 2 * kappa * theta > sigma * sigma?" << endl;
	if ( ( 2 * kappa * theta - sigma * sigma ) > 0 ) cout << "Yes! Good parameters!" << endl;
	else cout << "No! Poor job!" << endl;
	//cout << 2 * kappa * theta - sigma * sigma << endl;

	finish = clock ();
	
	cout << endl << endl;
	cout << "program operation time ( seconds ): " << ( ( float ) ( finish - start ) ) / CLOCKS_PER_SEC << endl << endl << endl;

	gsl_vector_free( x  ); 
	gsl_vector_free( ss ); 
	gsl_multimin_fminimizer_free ( s ); 

	free ( y  );

	return status; 
} 


