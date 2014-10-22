#include "PayOff.h"

PayOff::PayOff(){}                                      // note: Payoff default constructor has no colon!

PayOffPut::PayOffPut( const double& _K){                // note: PayOffPut default constructor has no colon!
    K = _K;
}

double PayOffPut::operator()(const double &S) const{    // we could initialize a object as PayOff X ( K );
                                                        // and return the payof of the X by X ( S );
    return std::max(K-S, 0.0); // Standard European put
}
