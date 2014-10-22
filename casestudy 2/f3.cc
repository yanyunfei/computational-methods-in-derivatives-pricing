#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *tridiagonalsolver ( double *, double *, double *, double *, int );

int main ( int agrc, char **argv )
{
	double S0 = 100;
	double K = 90;
	double r = 0.0025;
	double q = 0.0125;
	double T = 1;
	double sigma = 0.5;
	double deltasigma = 0.01;
	double smin = 0;
	double smax = 250;
	int    N = 1000;
	int    M = 250;
	double deltas = ( smax - smin ) / ( N );
	double deltat = T / ( M );

	// initialize grids for S and t
	double S [ N + 1 ];
	double t [ M + 1 ];

	int i;
	int j;
	for ( i = 0; i < N + 1; i++ )
	{
		S [ i ] = smin + i * deltas;
	}

	for ( i = 0; i < M + 1; i++ )
	{
		t [ i ] = i * deltat;
	}

	// initialize alpha and beta
	double alpha [ N + 1 ];
	double beta [ N + 1 ];
	for ( i = 0; i < N + 1; i++ )
	{
		alpha [ i ] = sigma * sigma * S [ i ] * S [ i ] * deltat / 2 / deltas / deltas;
		beta [ i ] = ( r - q ) * S [ i ] * deltat / 2 / deltas;
	}
	
	// initialize v
	double *v = ( double * ) malloc ( ( N - 1 ) * sizeof ( double ) );
	if ( v == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	for ( i = 0; i < N - 1; i++ )
	{
		* ( v + i ) = ( K > S [ i + 1 ] ) ? ( K - S [ i + 1 ] ) : 0;
	}

	// initialize the l, d, u
	double *l = ( double * ) malloc ( ( N - 1 ) * sizeof ( double ) );
	if ( l == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *d = ( double * ) malloc ( ( N - 1 ) * sizeof ( double ) );
	if ( d == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *u = ( double * ) malloc ( ( N - 1 ) * sizeof ( double ) );
	if ( u == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	for ( i = 0; i < N - 1; i++ )
	{
		* ( l + i ) = - alpha [ i + 1 ] + beta [ i + 1 ];
		* ( d + i ) = 1 + r * deltat + 2 * alpha [ i + 1 ];
		* ( u + i ) = - alpha [ i + 1 ] - beta [ i + 1 ];
	}

	// boundary condition initialization for elements l, d and u for matrix A
	* ( d + 0 ) = * ( d + 0 ) + * ( l + 0 ) + * ( l + 0 );
	* ( u + 0 ) = * ( u + 0 ) - * ( l + 0 );
	* ( d + N - 2 ) = * ( d + N - 2 ) + * ( u + N - 2 ) + * ( u + N - 2 );
	* ( l + N - 2 ) = * ( l + N - 2 ) - * ( u + N - 2 );  
	
	for ( i = 1; i < N - 1; i++ )
	{
		* ( d + i ) = * ( d + i ) - * ( l + i ) / * ( d + i - 1 ) * * ( u + i - 1 );
	}
	
	// repeat the tridiagonalsover to solve for the terminal v
	for ( i = 0; i < M; i++ )
//	{
		 tridiagonalsolver ( l, d, u, v, N - 1 );
//	}
	for ( i = 0; i < N - 1; i++ )
//	for ( i = 390; i < 405; i++ )
	{
		printf ( "%d\t\t%f\t%f\n", i + 1, S [ i + 1], * ( v + i ) );
	}
	
	free ( l );
	free ( d );
	free ( u );
	free ( v );
	return 0;
}
	
double *tridiagonalsolver ( double *l, double *d, double *u, double *v, int N )// N should equal to N - 1
{
	int i;
	for ( i = 1; i < N; i++ )
	{
		* ( v + i ) = * ( v + i ) - ( * ( l + i ) ) / ( * ( d + i - 1 ) ) * ( * ( v + i - 1 ) );
	}

	* ( v + N - 1 ) = * ( v + N - 1 ) / * ( d + N - 1 );
	
	for ( i = N - 2; i > -1; i-- )
	{
		* ( v + i ) = ( * ( v + i ) - ( * ( u + i ) ) * ( * ( v + i + 1 ) ) ) / ( * ( d + i ) );
	}

	return v;
}




