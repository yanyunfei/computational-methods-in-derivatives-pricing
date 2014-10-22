#include "OptionPricing.h"

//-----------------------------------------------------------------------------
//Struct OptionParameterSet is used to wrap the required parameters for option
//-----------------------------------------------------------------------------

OptionParameterSet::OptionParameterSet() {};

OptionParameterSet::OptionParameterSet(double _S0, double _T, double _r, double _q, double _K): S0(_S0), T(_T), r(_r), q(_q), K(_K) {
    type = 'C';
    ask = 0.0;
    bid = 0.0;
    weight = 1;
}

//-----------------------------------------------------------------------------
//Struct ModelParameterSet is used to wrap the required parameters for pricing model
//-----------------------------------------------------------------------------

ModelParameterSet::ModelParameterSet() {};

ModelParameterSet::ModelParameterSet(double _sigma, double _kappa, double _theta, double _rou, double _nu, double _eta, double _lambda): sigma(_sigma), kappa(_kappa), theta(_theta), rou(_rou), nu(_nu), eta(_eta), lambda(_lambda) {};


//-----------------------------------------------------------------------------
//Class FTPricing - used to conduct FFT method
//-----------------------------------------------------------------------------

//Constructor
FTPricing::FTPricing(OptionParameterSet _op, ModelParameterSet _mp, dcmplx (*_CharacterFunc)(dcmplx, OptionParameterSet, ModelParameterSet)): op(_op), mp(_mp), CharacterFunc(_CharacterFunc) {
    IC = dcmplx(0.0,1.0);
    X = NULL;
};

//Function BuildXVector - compute X vector.
void FTPricing::BuildXVector(double lambda, double eta, double alpha, double k, int N){
    
    //Update X vector
    if (X != NULL)
        delete[] X;
    X = new dcmplx[N];
    
    double v;
    double C = exp(-op.r*op.T); //discount factor
    
    for (int i = 0; i<N; i++){
        v = i*eta;
        
        if (i == 0)
            X[i] = eta/2.0*C/(alpha + v*IC)/(alpha + v*IC + 1.0)*exp(-IC*v*(log(k)-lambda*N/2.0)) * CharacterFunc(v-(alpha+1)*IC, op, mp);
        else
            X[i] = eta*C/(alpha + v*IC)/(alpha + v*IC + 1.0) * exp( -IC*v*(log(k)-lambda*N/2.0) ) * CharacterFunc(v-(alpha+1)*IC, op, mp);
    }
    
    return;
}

//Function FastFT - conduct FastFT methods
//                - this method returns the price of the option
//                - this method also stores the prices for different k in the array prices
double FTPricing::FastFT(double* prices, double k, double eta, double alpha, int N){
    
 	double _lambda = 2.0*PI/N/eta;
    BuildXVector(_lambda, eta, alpha, k, N);
    
    dcmplx *FFT = new dcmplx[N];
    
    //Compute vector y
    FFT = FourierFuncs::FFT(X, N);
    
    //Computer price vector
    for (int i=0; i<N; i++)
        prices[i] = FFT[i].real() / PI * exp( -alpha*(log(k) - (N/2.0-i)*_lambda) );
    
    //clear
    delete[] FFT;
    
	return prices[N/2];
}

FTPricing::~FTPricing(){
    delete[] X;
    X = NULL;
}








