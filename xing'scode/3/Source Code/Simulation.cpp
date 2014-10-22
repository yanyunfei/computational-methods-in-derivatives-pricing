#include "Simulation.h"

RandomNumGenerator::RandomNumGenerator() {
    rgen = gsl_rng_alloc(gsl_rng_taus);
};

double RandomNumGenerator::U(){
    return gsl_rng_uniform_pos(rgen);
}

double RandomNumGenerator::normal_box_muller(double m, double s){
    
    return m + gsl_ran_gaussian(rgen, s);
}

double RandomNumGenerator::normal_box_muller(){
    return normal_box_muller(0,1);
}

void RandomNumGenerator::normal_box_muller(std::vector<double>& normalVector, double m, double s){
    for (std::vector<double>::iterator iter = normalVector.begin(); iter != normalVector.end(); ++iter)
        (*iter) = normal_box_muller(m, s);
}

void RandomNumGenerator::correlatedNormal(std::vector<double>& normalVector1, std::vector<double>& normalVector2, double rho, double m1, double s1, double m2, double s2){
    normal_box_muller(normalVector1, 0, 1);
    normal_box_muller(normalVector2, 0, 1);
    
    for (std::vector<double>::iterator iter1 = normalVector1.begin(), iter2 = normalVector2.begin(); iter1 != normalVector1.end(); ++iter1, ++iter2){
        (*iter2) = ( rho*(*iter1) + sqrt(1-rho*rho)*(*iter2) )*s2 + m2;
        (*iter1) = m1 + (*iter1)*s1;
    }
}

double RandomNumGenerator::gamma(double alpha, double lambda){
    return gsl_ran_gamma(rgen, alpha, lambda);
}


