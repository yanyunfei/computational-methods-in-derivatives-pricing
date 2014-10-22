#ifndef __CaseStudy3__OptionPricing
#define __CaseStudy3__OptionPricing

#include "FourierLib.h"

//-----------------------------------------------------------------------------------
//Struct OptionParameterSet is used to wrap the required parameters of Option
//-----------------------------------------------------------------------------------
struct OptionParameterSet{
    OptionParameterSet();
    OptionParameterSet(double _S0, double _T, double _r, double _q, double _K);
    
    double S0;
    double T;
    double r;
    double q;
    double K;
    
    double barrier;
    
    char type;
    double ask;
    double bid;
    double mid;
    
    double weight;
    
    bool operator<(const OptionParameterSet &other) const {
        return (*this).T < other.T;
    }
    
    bool operator==(const OptionParameterSet &other) const {
        return (*this).T == other.T;
    }
    
    void operator=(const OptionParameterSet &other) {
        (*this).S0 = other.S0;
        (*this).T = other.T;
        (*this).r = other.r;
        (*this).q = other.q;
        (*this).K = other.K;
    }
};

//-----------------------------------------------------------------------------------
//Struct ModelParameterSet is used to wrap the required parameters for pricing model
//-----------------------------------------------------------------------------------

struct ModelParameterSet{
    ModelParameterSet();
    ModelParameterSet(double _sigma, double _kappa, double _theta, double _rou, double _nu, double _eta, double _lambda);
    
    double sigma;
    double kappa;
    double theta;
    double rou;
    double nu;
    double eta;
    double lambda;

};

//-----------------------------------------------------------------------------
//Class FTPricing - used to conduct FFT and FrFFT methods
//-----------------------------------------------------------------------------
class FTPricing{
public:
    FTPricing(OptionParameterSet _op, ModelParameterSet _mp, dcmplx (*_CharacterFunc)(dcmplx, OptionParameterSet, ModelParameterSet));
    
    double FastFT(double* prices, double k, double eta, double alpha, int N);
    
    OptionParameterSet op;
    ModelParameterSet mp;
    
    virtual ~FTPricing();

protected:
    virtual void BuildXVector(double lambda, double eta, double alpha, double k, int N);
    
    dcmplx IC;
    
    dcmplx *X;
    dcmplx (*CharacterFunc)(dcmplx, OptionParameterSet, ModelParameterSet);
    
};

#endif

