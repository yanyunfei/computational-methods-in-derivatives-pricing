#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

typedef complex<double> dcmplx;
#define PI 3.14159265358979323846
#define IC dcmplx(0,1)

// FFT
dcmplx* FFT(dcmplx*, int);
dcmplx* IFFT(dcmplx*, int);
void FFT_calculate(dcmplx*, int, int, dcmplx*, dcmplx*, dcmplx*);
dcmplx* FFT_get_twiddle_factors(int);

double* FFTPricing(double*, double*, double*, double*, int, double, double); 

// Fractional FFT
double* FrFFTPricing(double*, double*, double*, double*, int, double, double, double); 

// COS
double CosChi(double, double, double, double, double);
double CosPhi(double, double, double, double, double);
double* COSPricing(double*, double*, double*, double*, int, double, double); 

// Character Function for Black-Scholes Model
dcmplx ChrFuncGBM(dcmplx, double, double*, double*); 

//---------------------------------------------------------------
// N must be a power of 2 
//---------------------------------------------------------------
dcmplx* FFT(dcmplx* x, int N) 
{
  dcmplx* out = (dcmplx*) malloc(sizeof(dcmplx) * N);
  dcmplx* scratch = (dcmplx*) malloc(sizeof(dcmplx) * N);
  dcmplx* twiddles = FFT_get_twiddle_factors(N);
  
  FFT_calculate(x, N, 1, out, scratch, twiddles);
  
  free(twiddles);
  free(scratch);
  return out;
}
//---------------------------------------------------------------
dcmplx* IFFT(dcmplx* x, int N) 
{
  dcmplx* out = (dcmplx*) malloc(sizeof(dcmplx) * N);
  dcmplx* scratch = (dcmplx*) malloc(sizeof(dcmplx) * N);
  dcmplx* twiddles = FFT_get_twiddle_factors(N);
  
  FFT_calculate(x, N, 1, out, scratch, twiddles);
  
  free(twiddles);
  free(scratch);

  // Calculate IFFT via reciprocity property of DFT. 
  int N2 = N/2;
  double tmp0, tmp1;
  out[0]   = dcmplx(real(out[0])/N, imag(out[0])/N);
  out[N2]  = dcmplx(real(out[N2])/N, imag(out[N2])/N);

  for(int i=1; i<N2; i++) {
    tmp0 = real(out[i])/N;      
    tmp1 = imag(out[i])/N;
    out[i] = dcmplx(real(out[N-i])/N, imag(out[N-i])/N);
    out[N-i] = dcmplx(tmp0, tmp1);
  }
  return out;
}
//---------------------------------------------------------------
void FFT_calculate(dcmplx* x, int N, int skip, dcmplx* X, dcmplx* D, dcmplx* twiddles) {
  dcmplx * E = D + N/2;
  int k;

  if (N == 1) {
    X[0] = x[0];
    return;
  }

  // use X as a scratch buffer 
  FFT_calculate(x, N/2, skip*2, E, X, twiddles);
  FFT_calculate(x + skip, N/2, skip*2, D, X, twiddles);

  for(k = 0; k < N/2; k++) {
    // Multiply entries of D by the twiddle factors e^(-2*pi*i/N * k) 
    D[k] = twiddles[k*skip] * D[k]; 
  }

  for(k = 0; k < N/2; k++) {
    X[k]       = E[k] + D[k];
    X[k + N/2] = E[k] - D[k];
  }
}
//---------------------------------------------------------------
dcmplx* FFT_get_twiddle_factors(int N) {
  dcmplx* twiddles = (dcmplx*) malloc(sizeof(dcmplx) * N/2);
  int k;
  for (k = 0; k != N/2; ++k) {
    twiddles[k] = exp(dcmplx(0, -2.0*PI*k/N));
  }
  return twiddles;
}

// FFT
double* FastFTPricing(double* K, double* PutPrice, double* common, double* p_mdl, int N, double fft_eta, double fft_alpha) //FastFT
{
	double S0 = common[0], T = common[1], r = common[2], q = common[3];
	double fft_lambda = 2*PI/N/fft_eta;
	double fft_const = exp(-r*T);
	dcmplx *FFTx = (dcmplx*) malloc(sizeof(dcmplx)*N);
	dcmplx *FFTy = (dcmplx*) malloc(sizeof(dcmplx)*N);
	double *output = (double*) malloc(sizeof(double)*N);
	double *fft_nu = (double*) malloc(sizeof(double)*N);	
	// begin iteration
	int i, j;
	for (j=0; j<6; j++) 
	{
		for (i=0; i<N; i++)
		{
			fft_nu[i] = i*fft_eta;
			
			if (i==0)
				FFTx[i]=fft_eta/2.0*fft_const/(fft_alpha+fft_nu[i]*IC)/(fft_alpha+fft_nu[i]*IC+1.0)
						*exp(-IC*(log(K[j])-fft_lambda*N/2.0)*fft_nu[i])*ChrFuncGBM(fft_nu[i]-(fft_alpha+1)*IC, S0, common, p_mdl);
			else
				FFTx[i]=fft_eta*fft_const/(fft_alpha+fft_nu[i]*IC)/(fft_alpha+fft_nu[i]*IC+1.0)
						*exp(-IC*(log(K[j])-fft_lambda*N/2.0)*fft_nu[i])*ChrFuncGBM(fft_nu[i]-(fft_alpha+1)*IC, S0, common, p_mdl);					
		}
		FFTy = FFT(FFTx, N);
		for (i=0; i<N; i++)
		{
			output[i] = FFTy[i].real()/PI*exp(-fft_alpha*(log(K[j])-(N/2.0-i)*fft_lambda));
		}
		PutPrice[j] = output[N/2];
	}
	return PutPrice;
}

// FrFFT
double* FrFFTPricing(double* K, double* PutPrice, double* common, double* p_mdl, int N, double fr_eta, double fr_alpha, double fr_lambda)
{
	double S0 = common[0], T = common[1], r = common[2], q = common[3];
	double fr_gamma = fr_eta*fr_lambda/2/PI;
	double fr_const = exp(-r*T);
	dcmplx *Frx = (dcmplx*) malloc(sizeof(dcmplx)*N);
	dcmplx *Fry = (dcmplx*) malloc(sizeof(dcmplx)*2*N);
	dcmplx *Frz = (dcmplx*) malloc(sizeof(dcmplx)*2*N);
	dcmplx *FryH = (dcmplx*) malloc(sizeof(dcmplx)*2*N);
	dcmplx *FrzH = (dcmplx*) malloc(sizeof(dcmplx)*2*N);
	dcmplx *FrXi = (dcmplx*) malloc(sizeof(dcmplx)*2*N);
	dcmplx *FrXiH = (dcmplx*) malloc(sizeof(dcmplx)*2*N);
	double *output = (double*) malloc(sizeof(double)*N);
	double *Frnu = (double*) malloc(sizeof(double)*N);
	// begin iteration
	int i, j;
	for (j=0; j<6; j++) 
	{
		for (i=0; i<N; i++)
		{
			Frnu[i] = i*fr_eta;

			if (i==0)
				Frx[i]=fr_eta/2.0*fr_const/(fr_alpha+Frnu[i]*IC)/(fr_alpha+Frnu[i]*IC+1.0)
						*exp(-IC*(log(K[j])-fr_lambda*N/2.0)*Frnu[i])*ChrFuncGBM(Frnu[i]-(fr_alpha+1)*IC, S0, common, p_mdl);
			else
				Frx[i]=fr_eta*fr_const/(fr_alpha+Frnu[i]*IC)/(fr_alpha+Frnu[i]*IC+1.0)
						*exp(-IC*(log(K[j])-fr_lambda*N/2.0)*Frnu[i])*ChrFuncGBM(Frnu[i]-(fr_alpha+1)*IC, S0, common, p_mdl);
		}
		for (i=0; i<2*N; i++)
		{
			if (i<N)
			{
				Fry[i] = Frx[i]*exp(-IC*PI*fr_gamma*double(i)*double(i));
				Frz[i] = exp(IC*PI*fr_gamma*double(i)*double(i));
			}
			else
			{
				Fry[i] = 0;
				Frz[i] = exp(IC*PI*fr_gamma*double(2*N-1-i)*double(2*N-1-i));
			}
		}
		FryH = FFT(Fry, 2*N);
		FrzH = FFT(Frz, 2*N);
		for (i=0; i<2*N; i++)
			FrXi[i] = FryH[i] * FrzH[i];
		FrXiH = IFFT(FrXi, 2*N);
		
		for (i=0; i<N; i++)
		{
			output[i] = (FrXiH[i]*exp(-IC*PI*fr_gamma*double(i)*double(i))).real()
						/PI*exp(-fr_alpha*(log(K[j])-(N/2.0-i)*fr_lambda));
		}
		PutPrice[j] = output[N/2];
	}
	return PutPrice;
}




//------------------------------------------------------------------
// Numerical Methods
double CosChi(double a, double b, double c, double d, double k)
{
	double theta1 = k*PI*(d-a)/(b-a);
	double theta2 = k*PI*(c-a)/(b-a);
	double alpha = k*PI/(b-a);
	double output = (cos(theta1)*exp(d)-cos(theta2)*exp(c)+alpha*sin(theta1)*exp(d)-alpha*sin(theta2)*exp(c))/(1+alpha*alpha);
	return output;
}

double CosPhi(double a, double b, double c, double d, double k)
{
	double output;
	if (k==0)
		output=d-c;
	else
	{
		double theta1 = k*PI*(d-a)/(b-a);
		double theta2 = k*PI*(c-a)/(b-a);
		double alpha = (b-a)/(k*PI);
		output = (sin(theta1)-sin(theta2))*alpha;	
	}
	return output;	
}

// COS
double* COSPricing(double* K, double* PutPrice, double* common, double* p_mdl, int N, double lowbnd, double uppbnd) //COS
{
	double S0 = common[0], T = common[1], r = common[2], q = common[3];
	int i, k;
	double *Value = (double*) malloc(sizeof(double)*N);
	dcmplx *ACoeff = (dcmplx*) malloc(sizeof(dcmplx)*N);
	double sum[6]={0, 0, 0, 0, 0, 0};
	for (i=0; i<6; i++) 
	{
		for (k=0; k<N; k++)
		{
			Value[k] = 2.0/(uppbnd-lowbnd)*K[i]*(-CosChi(lowbnd, uppbnd, lowbnd, 0, double(k))
				+CosPhi(lowbnd, uppbnd, lowbnd, 0, double(k)));
			
			ACoeff[k] = exp(dcmplx(0,-1)*double(k)*PI*lowbnd/(uppbnd-lowbnd)) 
								*ChrFuncGBM(double(k)*PI/(uppbnd-lowbnd), S0/K[i], common, p_mdl); 		
			
			if (k==0)
				sum[i] += ACoeff[k].real()*Value[k]/2.0;
			else
				sum[i] += ACoeff[k].real()*Value[k];
		}
	PutPrice[i] = exp(-r*T)*sum[i];
	}
	return PutPrice;
}

//------------------------------------------------------------------
// Character Functions for Black-Scholes Model
dcmplx ChrFuncGBM(dcmplx u, double S0, double* common, double* p_mdl)
{
	double T = common[1], r = common[2], q = common[3];
	double sigma = p_mdl[0];
	dcmplx im = u*(log(S0)+(r-q-sigma*sigma/2.0)*T);
	dcmplx re = -u*u*sigma*sigma*T/2.0;
	dcmplx output = exp(re+im*IC);
	return output;
}

int main()
{
	// Common Parameters
	double common[4] = {1800, 0.5, 0.0025, 0.0203};		// S0, T, r, q
	double K[4]={900, 1100, 1300, 1500};				// K
	double Sigma = 0.3;									// Sigma

	// FFT Parameters
	double EtaFFT = 0.25; 
	double AlphaFFT[4] = {-2, -5, -10, -20};
	int PowerFFT[4] = {8, 10, 12, 14};
	int NFFT = 1;

	// FrFFT Parameters
	double EtaFrFFT = 0.25;
	double AlphaFrFFT[4] = {-2, -5, -10, -20};
	int PowerFrFFT[4] = {7, 8, 9, 10};
	int NFrFFT = 1;
	double LamdaFrFFT = 0.1;

	// COS Parameters
	double LowerBound[4] = {-2, -5, -10, -20};
	double UpperBound[4] = {2, 5, 10, 20};

	// FFT Pricing
	double *PutPriceFFT = (double*) malloc(sizeof(K)); 
	PutPriceFFT = FastFTPricing(K, PutPriceFFT, common, &Sigma, int(pow(2, PowerFFT[0])), EtaFFT, AlphaFFT[0]);
	cout<<"FFT Pricing: "<<endl;
	cout<<endl<<"N = "<<pow(2, PowerFFT[0])<<", alpha = "<<AlphaFFT[0]<<", eta = "<<EtaFFT<<endl;
	cout<<"Price 1: "<<PutPriceFFT[0]<<", Price 2: "<<PutPriceFFT[1]<<", Price 3: "<<PutPriceFFT[2]<<", Price 4: "<<PutPriceFFT[3]<<endl;

	// FrFFT Pricing
	double *PutPriceFrFFT = (double*) malloc(sizeof(K));
	PutPriceFrFFT = FrFFTPricing(K, PutPriceFrFFT, common, &Sigma, int(pow(2, PowerFrFFT[0])), EtaFrFFT, AlphaFrFFT[0], LamdaFrFFT);
	cout<<"FrFFT Pricing: "<<endl;
	cout<<endl<<"N = "<<pow(2, PowerFrFFT[0])<<", alpha = "<<AlphaFrFFT[0]<<", eta = "<<EtaFrFFT<<", Lamda = "<<LamdaFrFFT<<endl;
	cout<<"Price 1: "<<PutPriceFrFFT[0]<<", Price 2: "<<PutPriceFrFFT[1]<<", Price 3: "<<PutPriceFrFFT[2]<<", Price 4: "<<PutPriceFrFFT[3]<<endl;

	// COS Pricing
	double *PutPriceCOS = (double*) malloc(sizeof(K));
	PutPriceCOS = COSPricing(K, PutPriceCOS, common, &Sigma, int(pow(2, PowerFFT[0])), LowerBound[0], UpperBound[0]);
	cout<<"COS Pricing: "<<endl;
	cout<<endl<<"N = "<<pow(2, PowerFFT[0])<<", LowerBound = "<<LowerBound[0]<<", UpperBound = "<<UpperBound[0]<<endl;
	cout<<"Price 1: "<<PutPriceCOS[0]<<", Price 2: "<<PutPriceCOS[1]<<", Price 3: "<<PutPriceCOS[2]<<", Price 4: "<<PutPriceCOS[3]<<endl;

	system("pause");

	return 0;
}
