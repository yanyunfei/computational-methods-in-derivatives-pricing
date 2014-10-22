#include "CharacteristicFuncLib.h"

//1. Geometric Brownian Motion
dcmplx GBMCharacterFunc(dcmplx u, ModelParameterSet p){
    
	dcmplx im = ( log(p.S0) + (p.r-p.q-p.sigma*p.sigma/2.0)*p.T ) * u;
	dcmplx re = -u*u*p.sigma*p.sigma*p.T/2.0;
	
    return exp(re+im*dcmplx(0,1));
}