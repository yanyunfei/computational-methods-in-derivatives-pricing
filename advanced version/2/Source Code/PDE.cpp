#include "PDE.h"

BlackScholesPDE::BlackScholesPDE(VanillaOption* _option) : option(_option){}            // constructor

double BlackScholesPDE::diff_coeff(double t, double x) const {                          // we could use BSPDE.diff_coeff ( t, x ) to do some calculation
    
    double vol = option->sigma;
    return 0.5*vol*vol*x*x;
}

double BlackScholesPDE::conv_coeff(double t, double x) const {
    return (option->r - option->q)*x;
}

double BlackScholesPDE::zero_coeff(double t, double x) const {
    return -(option->r);
}

double BlackScholesPDE::init_cond(double x) const {
    return option->pay_off->operator()(x);                                              // why operator () ( x )?
}




