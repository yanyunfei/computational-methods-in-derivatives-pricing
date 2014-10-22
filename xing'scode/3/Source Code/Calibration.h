#ifndef __CaseStudy3__Calibration__
#define __CaseStudy3__Calibration__

#include <iostream>
#include <vector>
#include <fstream>
#include "CharacteristicFuncLib.h"
#include <gsl/gsl_multimin.h>
#include "OptionPricing.h"
#include <algorithm>

struct RQPair{
    
    RQPair(double _R, double _Q, double _T) : R(_R), Q(_Q), T(_T) {};
    
    double R;
    double Q;
    double T;
};

class OptionLoader{
public:
    OptionLoader(std::string _file);
    
    //Load option data
    void LoadOption(int method);
    
    std::vector<OptionParameterSet> Options;
    
private:
    std::string file;
};

//Build Term Structure
class TermStructureBuilder{
public:
    TermStructureBuilder(double Smin, double Smax, int numOfSStep, double Tmin, double Tmax, int numOfTStep);
    
    void BuildRQ(std::vector<OptionParameterSet> data, std::vector<double>& R, std::vector<double>& Q, int method);
    
    std::vector<double> SGrid;
    std::vector<double> TGrid;
    
protected:
    
    double S_StepSize;
    double T_StepSize;
    
    //linear interpolation to build term structure
    void RQBuilder_linear(std::vector<OptionParameterSet> options, std::vector<double>& R, std::vector<double>& Q);
    //piecewise constant interpolation to build term structure
    void RQBuilder_step(std::vector<OptionParameterSet> options, std::vector<double>& R, std::vector<double>& Q);
    void RQBuilder_constant(std::vector<OptionParameterSet> options, std::vector<double>& R, std::vector<double>& Q);
};


double heston_objectiveFunct(const gsl_vector *v, void *params);
double vgsa_objectiveFunct(const gsl_vector *v, void *params);


class Calibrator{
public:
    Calibrator();
    
    std::vector<double> Calibrate(std::vector<double> initialVector, OptionParameterSet options[], double (*objectiveFunction)(const gsl_vector *v, void *params), std::ofstream &output);
    
private:
    
};

#endif /* defined(__CaseStudy3__Calibration__) */
