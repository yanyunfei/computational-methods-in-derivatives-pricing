#ifndef __CaseStudy3__CharacteristicFuncLib__
#define __CaseStudy3__CharacteristicFuncLib__

#include "OptionPricing.h"

//The Characteristic Function Library

//1. Heston stochastic volatility model
dcmplx cf_heston(dcmplx u, OptionParameterSet op, ModelParameterSet mp);

//2. VGSA
dcmplx cf_vgsa(dcmplx u, OptionParameterSet op, ModelParameterSet mp);
dcmplx vgsa_psi(dcmplx u, double sigma, double v, double theta);
dcmplx vgsa_phi(dcmplx u, double T, double y0, double k, double eta, double lambda);

//3. GBM
dcmplx cf_gbm(dcmplx u, OptionParameterSet op, ModelParameterSet mp);


#endif /* defined(__E4732CaseStudy1__CharacteristicFuncLib__) */
