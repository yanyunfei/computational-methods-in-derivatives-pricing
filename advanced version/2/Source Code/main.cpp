#include <iostream>
#include "PayOff.h"
#include "DerivativeLib.h"
#include "PDE.h"
#include "FDM.h"
#include <fstream>

int main(int argc, const char * argv[])
{
    //double S = 100;
    double K = 90;
    double r = 0.0025;
    double q = 0.0125;
    double T = 1.00;
    double sigma = 0.5;
    
    double S_max = 250;
    int N = 1000;
    int M = 250;
    
    PayOff* payoff_put = new PayOffPut(K);
    VanillaOption* put_option = new VanillaOption(K, r, q, T, sigma, payoff_put);
    
    BlackScholesPDE* bs_pde = new BlackScholesPDE(put_option);
    
    //Apply second order FDMImplicit
    FDMImplicit_Second fdm_second(S_max, N, T, M, bs_pde);
    
    std::ofstream fdm1_out("fdm_second.csv");
    fdm1_out.precision(8);
    fdm_second.compute(fdm1_out);
    
    //Apply higher order FDMImplicit
    FDMImplicit_Third fdm_third(S_max, N, T, M, bs_pde);
    
    std::ofstream fdm2_out("fdm_third.csv");
    fdm2_out.precision(8);
    fdm_third.compute(fdm2_out);

    return 0;
}

