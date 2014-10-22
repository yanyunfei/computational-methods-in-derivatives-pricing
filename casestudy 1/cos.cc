#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <complex>

using namespace std;

typedef complex <double> dcmplx;

#define PI 3.1415926535897932

int main ( int, char** )
{
	// characteristic function of log of stock price
	double r = 0.0025;
	double S0 = 1800;
	double q = 0.0203;
	double sigma = 0.3;
	double T = 0.5;
	double K [ 4 ] = { 900, 1100, 1300, 1500 };

	int    N = pow ( 2, 9 );
	double C = exp ( -r * T );

	double a [ 4 ] = { -2, -5, -10, -20 };
	double b [ 4 ] = {  2,  5,  10,  20 };
	double c [ 4 ] = { -2, -5, -10, -20 };
	double d [ 4 ] = {  0,  0,   0,   0 };


int l;
for ( l = 0; l < 4; l++ )
{
	printf ( "[ %f,%f ]", a [ l ], b [ l ] );
	printf ( "\n\n" );
	double chi [ N ];
	double phi [ N ];
	
	int i;
	for ( i = 0; i < N; i++ )
	{
		double coefupper = i * PI * ( d [ l ] - a [ l ] ) / ( b [ l ] - a [ l ] );
		double coeflower = i * PI * ( c [ l ] - a [ l ] ) / ( b [ l ] - a [ l ] );
		double coef = i * PI / ( b [ l ] - a [ l ] );
		
		chi [ i ] = ( 
				cos ( coefupper ) * exp ( d [ l ] ) - cos ( coeflower ) * exp ( c [ l ] ) 
			      + coef * sin ( coefupper ) * exp ( d [ l ] ) - coef * sin ( coeflower ) * exp ( c [ l ] ) 
			    )
		       	      / ( 1 + coef * coef );		
	}

	for ( i = 0; i < N; i++ )
	{
		if ( i == 0 )
		{
			phi [ i ] = d [ l ] - c [ l ];
		}
		else
		{
			double coefupper = i * PI * ( d [ l ] - a [ l ] ) / ( b [ l ] - a [ l ] );
			double coeflower = i * PI * ( c [ l ] - a [ l ] ) / ( b [ l ] - a [ l ] );
			double coef = i * PI / ( b [ l ] - a [ l ] );
	
			phi [ i ] = ( sin ( coefupper ) - sin ( coeflower ) ) / coef;
		}
	}

	int j;
	for ( j = 0; j < 4; j++ )
	{
		double V [ N ];
		for ( i = 0; i < N; i++ )
		{	
			V [ i ] = 2.0 / ( b [ l ] - a [ l ] ) * K [ j ] * ( -chi [ i ] + phi [ i ] );
		}
	
		dcmplx nu [ N ];
		for ( i = 0; i < N; i++ )
		{
			double coef = i * PI / ( b [ l ] - a [ l ] );

			nu [ i ] = exp ( 
					 ( log ( S0 / K [ j ] ) + ( r - q - 0.5 * sigma * sigma ) * T ) * dcmplx ( 0, 1 ) * coef
			   		 - 0.5 * sigma * sigma * T * coef * coef 
					 );
		}

		double value = 0.5 * nu [ 0 ].real () * V [ 0 ];
		for ( i = 1; i < N; i++ )
		{
			value = value + ( nu [ i ] * exp ( i * PI * dcmplx ( 0, -1 ) * a [ l ] / ( b [ l ] - a [ l ] ) ) ).real () * V [ i ];
		}
		
		printf ( "%f\t%f\n", value * C, K [ j ] );
	}
	printf ( "\n" );
}	
	return 0;
}

		
