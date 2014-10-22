#include "CharacteristicFuncLib.h"

//1. Heston stochastic volatility model
dcmplx cf_heston(dcmplx u, OptionParameterSet op, ModelParameterSet mp){
    
    dcmplx IC = dcmplx(0,1);
    
    double T = op.T, r = op.r, q = op.q;
    
	double sigma = mp.sigma, kappa = mp.kappa, theta = mp.theta, rou = mp.rou, nu = mp.nu;
    
	dcmplx gamma = sqrt(sigma*sigma*(u*u+IC*u)+(kappa-IC*rou*sigma*u)*(kappa-IC*rou*sigma*u));
	
    dcmplx numcomp1 = IC*u*log(op.S0) + IC*u*(r-q)*T;
    dcmplx numcomp2 = kappa*theta*T*(kappa-IC*rou*sigma*u)/sigma/sigma;
    dcmplx numcomp3 = -(u*u+IC*u)*nu/(gamma/tanh(gamma*T/2.0)+kappa-IC*rou*sigma*u);
	
    dcmplx dencomp1 = cosh(gamma*T/2.0)+(kappa-IC*rou*sigma*u)/gamma*sinh(gamma*T/2.0);
	dcmplx dencomp2 = 2.0*kappa*theta/sigma/sigma;
	dcmplx output = exp(numcomp1+numcomp2+numcomp3)/pow(dencomp1, dencomp2);
	
    return output;
}

//2. VGSA
dcmplx vgsa_psi(dcmplx u, double sigma, double v, double theta){
    
	dcmplx output = -1.0/v*log(1.0-dcmplx(0.0,1.0)*u*theta*v+sigma*sigma*v*u*u/2.0);
	return output;
}

dcmplx vgsa_phi(dcmplx u, double T, double y0, double k, double eta, double lambda){
    
	dcmplx gamma = sqrt(k*k-dcmplx(0.0,2.0)*lambda*lambda*u);
	dcmplx A = exp(k*k*eta*T/lambda/lambda)/pow(cosh(gamma*T/2.0)+k/gamma*sinh(gamma*T/2.0), 2.0*k*eta/lambda/lambda);
	dcmplx B = dcmplx(0,2)*u/(k + gamma/tanh(gamma*T/2.0) );
	return A*exp(B*y0);
}

dcmplx cf_vgsa(dcmplx u, OptionParameterSet op, ModelParameterSet mp){
    
    dcmplx IC = dcmplx(0.0,1.0);
    
	dcmplx comp1 = exp(IC*u*(log(op.S0)+(op.r-op.q)*op.T));
	dcmplx psiu = vgsa_psi(u, mp.sigma, mp.nu, mp.theta);
	dcmplx psii = vgsa_psi(dcmplx(0.0,-1.0), mp.sigma, mp.nu, mp.theta);
	dcmplx comp2 = vgsa_phi(dcmplx(0.0,-1.0)*psiu, op.T, 1.0/mp.nu, mp.kappa, mp.eta, mp.lambda);
	dcmplx comp3 = pow(vgsa_phi(dcmplx(0.0,-1.0)*psii, op.T, 1.0/mp.nu, mp.kappa, mp.eta, mp.lambda),dcmplx(0.0,1.0)*u);
	dcmplx output = comp1*comp2/comp3;
	return output;
}

//3. cf_gbm
dcmplx cf_gbm(dcmplx u, OptionParameterSet op, ModelParameterSet mp){
    
	dcmplx im = ( log(op.S0) + (op.r-op.q-mp.sigma*mp.sigma/2.0)*op.T ) * u;
	dcmplx re = -u*u*mp.sigma*mp.sigma*op.T/2.0;
	
    return exp(re+im*dcmplx(0,1));
}

