#include "DerivativeLib.h"

VanillaOption::VanillaOption() {}

VanillaOption::VanillaOption(double _K, double _r , double _q, double _T, double _sigma, PayOff* _pay_off) : K(_K), r(_r), q(_q), T(_T), sigma(_sigma), pay_off(_pay_off) {}