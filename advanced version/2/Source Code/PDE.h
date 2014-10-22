#ifndef __CaseStudy2__PDE__
#define __CaseStudy2__PDE__

#include <iostream>
#include "DerivativeLib.h"
#include <math.h>

class BlackScholesPDE {
        public:
                VanillaOption* option;                                  // note: it's a pointer, which means we could access the member of VanillaOption directly by ->

                BlackScholesPDE(VanillaOption* option);

                virtual double diff_coeff(double t, double x) const;
                virtual double conv_coeff(double t, double x) const;
                virtual double zero_coeff(double t, double x) const;

                virtual double init_cond(double x) const;
};


#endif
