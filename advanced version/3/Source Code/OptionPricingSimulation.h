#ifndef __CaseStudy3__OptionPricingSimulation__
#define __CaseStudy3__OptionPricingSimulation__

#include <iostream>
#include "Simulation.h"
#include "OptionPricing.h"
#include "CharacteristicFuncLib.h"
#include <fstream>

class OptionSimulation{
public:
    OptionSimulation(OptionParameterSet _op, ModelParameterSet _mp);
    
    virtual double Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q);
    
    virtual ~OptionSimulation();
    
protected:
    OptionParameterSet op;
    ModelParameterSet mp;
};

//Heston Simulation for UOC Pricing
class Heson_UOC: public OptionSimulation{
public:
    Heson_UOC(OptionParameterSet _op, ModelParameterSet _mp);
    
    virtual double Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q);
};

//VGSA Simulation for UOC Pricing
class VGSA_UOC: public OptionSimulation{
public:
    VGSA_UOC(OptionParameterSet _op, ModelParameterSet _mp);
    
    virtual double Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q);
};

//LocalVol Simulation for UOC Pricing
class LocalVol_UOC: public OptionSimulation{
public:
    LocalVol_UOC(OptionParameterSet _op, ModelParameterSet _mp);
    void SetLocalVol(std::vector< std::vector<double> > _localVol, std::vector<double> _TGrid, std::vector<double> _SGrid);
    
    virtual double Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q);
    
protected:
    std::vector< std::vector<double> > localVol;
    std::vector<double> Tgrid;
    std::vector<double> Sgrid;
    
    double getLocalVol(double T, double S);
};


#endif /* defined(__CaseStudy3__OptionPricingSimulation__) */
