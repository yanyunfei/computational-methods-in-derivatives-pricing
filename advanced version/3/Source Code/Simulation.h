#ifndef __CaseStudy3__Simulation__
#define __CaseStudy3__Simulation__

#include <iostream>
#include <math.h>
#include <vector>
#include <gsl/gsl_randist.h>
#include <random>

class RandomNumGenerator{
public:
    RandomNumGenerator();
    
    double U();
    
    double normal_box_muller(double m, double s);
    double normal_box_muller();
    void normal_box_muller(std::vector<double>& normalVector, double m, double s);
    void correlatedNormal(std::vector<double>& normalVector1, std::vector<double>& normalVector2, double rho, double m1, double s1, double m2, double s2);
    
    double gamma(double alpha, double lambda);

private:
    gsl_rng * rgen;
};



#endif /* defined(__CaseStudy3__Simulation__) */

