#ifndef __CaseStudy1__Pricing_ByFourier__
#define __CaseStudy1__Pricing_ByFourier__

#include "FourierLib.h"

//-----------------------------------------------------------------------------------
//Struct ModelParameterSet is used to wrap the required parameters for pricing model
//-----------------------------------------------------------------------------------
struct ModelParameterSet{
    ModelParameterSet(double _S0, double _T, double _r, double _q, double _sigma);
    
    double S0;
    double T;
    double r;
    double q;
    double sigma;
};


//-----------------------------------------------------------------------------
//Class FTPricing - used to conduct FFT and FrFFT methods
//-----------------------------------------------------------------------------
class FTPricing{
public:
    FTPricing(ModelParameterSet _parameters, dcmplx (*_CharacterFunc)(dcmplx, ModelParameterSet));
    
    double FastFT(double* putPrices, double eta, double alpha, double k, int N);
    double FrFFT(double* putPrices, double eta, double alpha, double k, int N, double lambda);
    
    virtual ~FTPricing();

protected:
    void BuildXVector(double lambda, double eta, double alpha, double k, int N);
    
    ModelParameterSet p;
    dcmplx IC;
    
    dcmplx *X;
    dcmplx (*CharacterFunc)(dcmplx, ModelParameterSet);
    
};

//-----------------------------------------------------------------------------
//Class COSPricing - used to conduct COS methods
//-----------------------------------------------------------------------------
class COSPricing{
public:
    COSPricing(ModelParameterSet _parameters, dcmplx (*_CharacterFunc)(dcmplx, ModelParameterSet));
    double COS(double k, int N, double a, double b);
    virtual ~COSPricing(){};
    
protected:
    double ComputeChi(double a, double b, double c, double d, double k);
    double ComputePhi(double a, double b, double c, double d, double k);
    
    ModelParameterSet p;
    dcmplx (*CharacterFunc)(dcmplx, ModelParameterSet);
};

#endif /* defined(__CaseStudy1__Pricing_ByFourier__) */
