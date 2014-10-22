#include "OptionPricing.h"

//-----------------------------------------------------------------------------
//Struct ModelParameterSet is used to wrap the required parameters for pricing model
//-----------------------------------------------------------------------------

ModelParameterSet::ModelParameterSet(double _S0, double _T, double _r, double _q, double _sigma): S0(_S0), T(_T), r(_r), q(_q), sigma(_sigma) {};


//-----------------------------------------------------------------------------
//Class FTPricing - used to conduct FFT and FrFFT methods
//-----------------------------------------------------------------------------

//Constructor
FTPricing::FTPricing(ModelParameterSet _parameters, dcmplx (*_CharacterFunc)(dcmplx, ModelParameterSet)): p(_parameters), CharacterFunc(_CharacterFunc) {
    IC = dcmplx(0,1);
    X = NULL;
};

//Function BuildXVector - compute X vector.
void FTPricing::BuildXVector(double lambda, double eta, double alpha, double k, int N){
    
    //Update X vector
    if (X != NULL)
        delete[] X;
    X = new dcmplx[N];
    
    double v;
    double C = exp(-p.r*p.T); //discount factor
    
    for (int i = 0; i<N; i++){
        v = i*eta;
        
        if (i == 0)
            X[i] = eta/2.0*C/(alpha + v*IC)/(alpha + v*IC + 1.0)*exp(-IC*v*(log(k)-lambda*N/2.0)) * CharacterFunc(v-(alpha+1)*IC, p);
        else
            X[i] = eta*C/(alpha + v*IC)/(alpha + v*IC + 1.0) * exp( -IC*v*(log(k)-lambda*N/2.0) ) * CharacterFunc(v-(alpha+1)*IC, p);
    }
}

//Function FastFT - conduct FastFT methods
//                - this method returns the price of the option
//                - this method also stores the prices for different k in the array prices
double FTPricing::FastFT(double* prices, double k, double eta, double alpha, int N){
    
	double _lambda = 2*PI/N/eta;
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


//Function FrFFT - conduct FrFFT method
//               - this method returns the price of the option
//               - this method also stores the prices for different k in the array prices
double FTPricing::FrFFT(double* prices, double k, double eta, double alpha, int N, double lambda){
    
	double gamma = eta*lambda/2/PI;
    BuildXVector(lambda, eta, alpha, k, N);
    
    dcmplx *Fry = new dcmplx[2*N];
    dcmplx *Frz = new dcmplx[2*N];
    dcmplx *FrXi = new dcmplx[2*N];
    
    for (int i=0; i<2*N; i++){
        
        if (i<N){
            double i_square = i*i;
            Fry[i] = X[i] * exp(-IC*PI*gamma*i_square);
            Frz[i] = exp(IC*PI*gamma*i_square);
        }
        else{
            Fry[i] = 0;
            Frz[i] = exp(IC*PI*gamma*double(2*N-1-i)*double(2*N-1-i));
        }
    }
    
    Fry = FourierFuncs::FFT(Fry, 2*N);
    Frz = FourierFuncs::FFT(Frz, 2*N);
    
    for (int i=0; i<2*N; i++)
        FrXi[i] = Fry[i] * Frz[i];
    
    FrXi = FourierFuncs::IFFT(FrXi, 2*N);
    
    for (int i=0; i<N; i++){
        prices[i] = (FrXi[i] * exp(-IC*PI*gamma*double(i)*double(i)) ).real() / PI * exp( -alpha*(log(k)-(N/2.0-i)*lambda) );
    }
    
    //clear
    delete[] Fry;
    delete[] Frz;
    delete[] FrXi;

	return prices[N/2];
}

FTPricing::~FTPricing(){
    delete[] X;
    X = NULL;
}

//-----------------------------------------------------------------------------
//Class COSPricing - used to conduct COS methods
//-----------------------------------------------------------------------------

//Constructor
COSPricing::COSPricing(ModelParameterSet _parameters, dcmplx (*_CharacterFunc)(dcmplx, ModelParameterSet)): p(_parameters), CharacterFunc(_CharacterFunc) {};

double COSPricing::ComputeChi(double a, double b, double c, double d, double k){
    
	double theta1 = k*PI*(d-a)/(b-a);
	double theta2 = k*PI*(c-a)/(b-a);
	double alpha = k*PI/(b-a);
	double output = (cos(theta1)*exp(d)-cos(theta2)*exp(c)+alpha*sin(theta1)*exp(d)-alpha*sin(theta2)*exp(c))/(1+alpha*alpha);
	return output;
}

double COSPricing::ComputePhi(double a, double b, double c, double d, double k){
    
    double output;
	if (k == 0)
		output = d-c;
	else{
        
		double theta1 = k*PI*(d-a)/(b-a);
		double theta2 = k*PI*(c-a)/(b-a);
		double alpha = (b-a)/(k*PI);
		output = (sin(theta1)-sin(theta2))*alpha;
	}
    
	return output;
}

double COSPricing::COS(double k, int N, double a, double b){
		
    double L = 10;
    double c1 = (p.r-p.q)*p.T;
    double c2 = p.sigma*p.sigma*p.T;
    double c4 = 0;
    
    double lowbnd, upbnd;
    if (a == 0)
        lowbnd = c1-L*sqrt(c2+sqrt(c4)); //a
    else
        lowbnd = a;
    
    if (b == 0)
        upbnd = c1+L*sqrt(c2+sqrt(c4)); //b
    else
        upbnd = b;
    
    double sum = 0;
	
    for (int i=0; i<N; ++i){
        
        p.S0 = p.S0/k;
        double Value = 2.0/(upbnd-lowbnd)*k*(-ComputeChi(lowbnd, upbnd, lowbnd, 0, double(i)) + ComputePhi(lowbnd, upbnd, lowbnd, 0, double(i)));
        dcmplx ACoeff = exp(dcmplx(0,-1)*double(i)*PI*lowbnd/(upbnd-lowbnd)) * CharacterFunc(double(i)*PI/(upbnd-lowbnd), p);
        p.S0 = p.S0*k;
        
        if (i == 0)
            sum += ACoeff.real()*Value/2.0;
        else
            sum += ACoeff.real()*Value;
    }
    
    return exp(-p.r*p.T)*sum;
}




