#ifndef __CaseStudy3__LocalVolSurface__
#define __CaseStudy3__LocalVolSurface__

#include <iostream>
#include <vector>
#include <fstream>
#include "CharacteristicFuncLib.h"
#include <gsl/gsl_multimin.h>
#include "OptionPricing.h"

class LocalVolSurfaceBuilder{
public:
    LocalVolSurfaceBuilder(std::vector<double> _SGrid, std::vector<double> _TGrid, std::vector<double> _S_R, std::vector<double> _S_Q, std::vector<OptionParameterSet> options);
    
    void Build( double S0, dcmplx (*_CharacterFunc)(dcmplx, OptionParameterSet, ModelParameterSet), ModelParameterSet mp, std::string title);
    std::vector< std::vector<double> > getLocalVolSurface();
    
    virtual ~LocalVolSurfaceBuilder();

private:
    std::vector<double> SGrid;
    std::vector<double> TGrid;
    
    std::vector<double> S_R;
    std::vector<double> S_Q;
    
    std::vector< std::vector<double> > *OptionPremiumSurface;
    std::vector< std::vector<double> > *LocalVolSurface;
};

#endif