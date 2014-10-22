#ifndef __CaseStudy2__DerivativeLib__
#define __CaseStudy2__DerivativeLib__

#include <iostream>
#include "PayOff.h"

class VanillaOption {
public:
    //Parameters of vanilla option
    PayOff* pay_off;                                                                    // note, pay_off is a pointer which points to PayOff
    double K;
    double r;
    double q;
    double T;
    double sigma;
    
    VanillaOption();
    VanillaOption(double _K, double _r, double _q, double _T, double _sigma, PayOff* _pay_off);
};

#endif
