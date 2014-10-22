#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *pentadiagonalsolver ( double *, double *, double *, double *, double *, double *, int );

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
	double *alpha = ( double * ) malloc ( ( N + 1 ) * sizeof ( double ) );
	double *beta  = ( double * ) malloc ( ( N + 1 ) * sizeof ( double ) );
	for ( i = 0; i < N + 1; i++ )
	{
		* ( alpha + i ) = sigma * sigma * S [ i ] * S [ i ] * deltat / 24 / deltas / deltas;
		* ( beta  + i ) = ( r - q ) * S [ i ] * deltat / 12 / deltas;
	}
	
	// initialize v
	double *v = ( double * ) malloc ( ( N - 3 ) * sizeof ( double ) );
	if ( v == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	for ( i = 0; i < N - 3; i++ )
	{
		* ( v + i ) = ( K > S [ i + 2 ] ) ? ( K - S [ i + 2 ] ) : 0;
	}

	// initialize the k, l, d, u, m
	double *k = ( double * ) malloc ( ( N - 3 ) * sizeof ( double ) );
	if ( k == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *l = ( double * ) malloc ( ( N - 3 ) * sizeof ( double ) );
	if ( l == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *d = ( double * ) malloc ( ( N - 3 ) * sizeof ( double ) );
	if ( d == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *u = ( double * ) malloc ( ( N - 3 ) * sizeof ( double ) );
	if ( u == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *m = ( double * ) malloc ( ( N - 3 ) * sizeof ( double ) );
	if ( m == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	for ( i = 0; i < N - 3; i++ )
	{
		* ( k + i ) = alpha [ i + 2 ] - beta [ i + 2 ];
		* ( l + i ) = - 16 * alpha [ i + 2 ] + 8 * beta [ i + 2 ];
		* ( d + i ) = 1 + r * deltat + 30 * alpha [ i + 2 ];
		* ( u + i ) = - 16 * alpha [ i + 2 ] - 8 * beta [ i + 2 ];
		* ( m + i ) = beta [ i + 2 ] + alpha [ i + 2 ];
	}

	// boundary condition initialization for elements l, d and u for matrix A
	* ( d + 0 ) = * ( d + 0 ) + * ( l + 0 ) * 27 / 13 + * ( k + 0 ) * 42 / 13;
	* ( u + 0 ) = * ( u + 0 ) - * ( l + 0 ) * 15 / 13 - * ( k + 0 ) * 32 / 13;
	* ( m + 0 ) = * ( m + 0 ) + * ( l + 0 ) *  1 / 13 + * ( k + 0 ) *  3 / 13;
	
	* ( l + 1 ) = * ( l + 1 ) + * ( k + 1 ) * 27 / 13;
	* ( d + 1 ) = * ( d + 1 ) - * ( k + 1 ) * 15 / 13;
	* ( u + 1 ) = * ( u + 1 ) + * ( k + 1 ) *  1 / 13;
	
	* ( d + N - 4 ) = * ( d + N - 4 ) + * ( u + N - 4 ) * 27 / 13 + * ( m + N - 4 ) * 42 / 13;
	* ( l + N - 4 ) = * ( l + N - 4 ) - * ( u + N - 4 ) * 15 / 13 - * ( m + N - 4 ) * 32 / 13;
	* ( k + N - 4 ) = * ( k + N - 4 ) + * ( u + N - 4 ) *  1 / 13 + * ( m + N - 4 ) *  3 / 13;

	* ( u + N - 5 ) = * ( u + N - 5 ) + * ( m + N - 5 ) * 27 / 13;
	* ( d + N - 5 ) = * ( d + N - 5 ) - * ( m + N - 5 ) * 15 / 13;
	* ( l + N - 5 ) = * ( l + N - 5 ) + * ( m + N - 5 ) *  1 / 13;

//	for ( i = 0; i < N - 3; i++ )
//	{
//		printf ( "%d\t\t%f\t%f\t%f\t%f\t%f\t%f\n", i + 1, * ( k + i ), * ( l + i ), * ( d + i ), * ( u + i ), * ( m + i ), * ( v + i ) );
//	}

	// repeat the tridiagonalsover to solve for the terminal v
	for ( i = 0; i < M; i++ )
	{
		v = pentadiagonalsolver ( k, l, d, u, m, v, N - 3 );
	}
	
//	for ( i = 390; i < 405; i++ )
	for ( i = 0; i < 997; i++ )
	{
		printf ( "%d\t\t%f\t%f\n", i + 1,  S [ i + 2 ], * ( v + i ) );
	}
	
	free ( k );
	free ( l );
	free ( d );
	free ( u );
	free ( m );
	free ( v );
	free ( alpha );
	free ( beta );
	return 0;
}
	
double *pentadiagonalsolver ( double *k, double *l, double *d, double *u, double *m, double *v, int N )// N should equal to N - 1
{
	
	double *kk = ( double * ) malloc ( N * sizeof ( double ) );
	if ( kk == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *ll = ( double * ) malloc ( N * sizeof ( double ) );
	if ( ll == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *dd = ( double * ) malloc ( N * sizeof ( double ) );
	if ( dd == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *uu = ( double * ) malloc ( N * sizeof ( double ) );
	if ( uu == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *mm = ( double * ) malloc ( N * sizeof ( double ) );
	if ( mm == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	double *vv = ( double * ) malloc ( N * sizeof ( double ) );
	if ( vv == NULL )
	{
		perror ( "malloc returned NULL" );
		exit ( 1 );
	}

	int i;
	for ( i = 0; i < N; i++ )
	{
		* ( kk + i ) = * ( k + i );
		* ( ll + i ) = * ( l + i );
		* ( dd + i ) = * ( d + i );
		* ( uu + i ) = * ( u + i );
		* ( mm + i ) = * ( m + i );
		* ( vv + i ) = * ( v + i );
	}

//	for ( i = 0; i < N; i++ )
//	{
//		printf ( "%d\t\t%f\t%f\t%f\t%f\t%f\t%f\n", i + 1, * ( kk + i ), * ( ll + i ), * ( dd + i ), * ( uu + i ), * ( mm + i ), * ( vv + i ) );
//	}

//	printf ( "%f\n", *vv );
	for ( i = 1; i < N - 1; i++ )
	{
		* ( dd + i ) = * ( dd + i ) - * ( ll + i ) / * ( dd + i - 1 ) * * ( uu + i - 1 );
		* ( uu + i ) = * ( uu + i ) - * ( ll + i ) / * ( dd + i - 1 ) * * ( mm + i - 1 );
		* ( vv + i ) = * ( vv + i ) - * ( ll + i ) / * ( dd + i - 1 ) * * ( vv + i - 1 );
		
		* ( ll + i + 1 ) = * ( ll + i + 1 ) - * ( kk + i + 1 ) / * ( dd + i - 1 ) * * ( uu + i - 1 );
		* ( dd + i + 1 ) = * ( dd + i + 1 ) - * ( kk + i + 1 ) / * ( dd + i - 1 ) * * ( mm + i - 1 );
		* ( vv + i + 1 ) = * ( vv + i + 1 ) - * ( kk + i + 1 ) / * ( dd + i - 1 ) * * ( vv + i - 1 );	
	}
	
	* ( v + N - 1 ) = * ( vv + N - 1 ) / * ( dd + N - 1 );
	* ( v + N - 2 ) = ( * ( vv + N - 2 ) - * ( uu + N - 2 ) * * ( v + N - 1 ) ) / * ( dd + N - 2 );
	
	for ( i = N - 3; i > -1; i-- )
	{
		* ( v + i ) = ( * ( vv + i ) - * ( mm + i ) * * ( v + i + 2 ) - * ( uu + i ) * * ( v + i + 1 ) ) / ( * ( dd + i ) );
	}

//	for ( i = 0; i < N; i++ )
//	{
//		printf ( "%d\t\t%f\t%f\t%f\t%f\t%f\t%f\n", i + 1, * ( kk + i ), * ( ll + i ), * ( dd + i ), * ( uu + i ), * ( mm + i ), * ( vv + i ) );
//	}

	free ( kk );
	free ( ll );
	free ( dd );
	free ( uu );
	free ( mm );
	free ( vv );

	return v;
}




