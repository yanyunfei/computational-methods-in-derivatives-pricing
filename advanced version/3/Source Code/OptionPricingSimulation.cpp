#include "OptionPricingSimulation.h"

OptionSimulation::OptionSimulation(OptionParameterSet _op, ModelParameterSet _mp): op(_op), mp(_mp) {};
double OptionSimulation::Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q) { return 0.0;};
OptionSimulation::~OptionSimulation() {};

Heson_UOC::Heson_UOC(OptionParameterSet _op, ModelParameterSet _mp) : OptionSimulation(_op, _mp) {};

double Heson_UOC::Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q) {
    
    RandomNumGenerator rand;
    
    double result;
    
    double z1, z2;
    
    for (int i = 0; i < N; ++i) {
        
        double S = log( op.S0 );
        double V, V_pre = mp.nu;
        
        for (int j = 0; j < numStep; ++j) {
            
            z1 = rand.normal_box_muller();
            z2 = rand.normal_box_muller();
            
            V = V_pre + mp.kappa*(mp.theta - std::max(0.0,V_pre) )*dt + mp.sigma*sqrt( std::max(0.0,V_pre) )*sqrt(dt)*z1;
            
            S = S + (R[j] - Q[j] - 0.5*V_pre)*dt + sqrt( std::max(0.0,V_pre) )*sqrt(dt)*( z1*mp.rou + sqrt(1-mp.rou*mp.rou)*z2 );
            
            V_pre = V;
            
            if (exp(S) >= op.barrier){
                S = 0;
                break;
            }
        }
        
        result += std::max(0.0, exp(S) - op.K);
    }
    
    return result*exp(-R.back()*op.T)/(double)N;
}

VGSA_UOC::VGSA_UOC(OptionParameterSet _op, ModelParameterSet _mp) : OptionSimulation(_op, _mp) {};

double VGSA_UOC::Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q) {
    
    RandomNumGenerator rand;
    
    dcmplx IC = dcmplx(0.0,1.0);
    dcmplx psi = vgsa_psi(-IC, mp.sigma, mp.nu, mp.theta);
    
    double S, w;
    double x , y , y_pre, z, g;
    double result = 0.0;
    
    for (int i = 0; i < N; ++i) {
        
        S = log(op.S0);
        y_pre = 1.0/mp.nu;
        
        for (int j = 0; j < numStep; ++j) {
            
            w = log( vgsa_phi(-IC*psi, (j-1)*dt, 1/mp.nu, mp.kappa, mp.eta, mp.lambda).real() )
            - log( vgsa_phi(-IC*psi, j*dt, 1/mp.nu, mp.kappa, mp.eta, mp.lambda).real() );
            
            z = rand.normal_box_muller();
            
            y = y_pre + mp.kappa*(mp.eta - y_pre)*dt + mp.lambda*sqrt(y_pre*dt)*z + mp.lambda*mp.lambda/4.0*dt*(z*z - 1.0);
            y = (y>0) ? y : 0.0;
            
            g = rand.gamma(dt*(y + y_pre)/mp.nu/2.0, mp.nu);
            x = mp.theta*g + mp.sigma*sqrt(g)*rand.normal_box_muller();
            
            y_pre = y;
            
            S = S + (R[j] - Q[j])*dt + w + x;
            
            if ( exp(S) >= op.barrier){
                S = 0;
                break;
            }
        }
        
        result += std::max(0.0, exp(S) - op.K);
    }
    
    return result*exp(-R.back()*op.T)/(double)N;
}


LocalVol_UOC::LocalVol_UOC(OptionParameterSet _op, ModelParameterSet _mp) : OptionSimulation(_op, _mp) {};

void LocalVol_UOC::SetLocalVol(std::vector< std::vector<double> > _localVol, std::vector<double> _TGrid, std::vector<double> _SGrid){
    localVol = _localVol;
    Tgrid = _TGrid;
    Sgrid = _SGrid;
}

double LocalVol_UOC::getLocalVol(double T, double S){
    int iterT = 0;
    int iterS = 0;
    
    for ( ; iterT < Tgrid.size(); ++iterT) {
        if ( Tgrid[iterT] >= T )
            break;
    }
    
    for ( ; iterS < Sgrid.size(); ++iterS) {
        if ( Sgrid[iterS] >= S )
            break;
    }
    
    iterT = (iterT ==  Tgrid.size()) ? --iterT : iterT;
    iterS = (iterS ==  Sgrid.size()) ? --iterS : iterS;
    
    return localVol[iterT][iterS];
}

double LocalVol_UOC::Simulate(int N, double dt, int numStep, std::vector<double> R, std::vector<double> Q) {
    
    RandomNumGenerator rand;
    
    double result = 0.0;
    
    for (int i = 0; i < N; ++i) {
        
        double S = op.S0;
        
        for (int j = 0; j < numStep; ++j) {
            
            S = S + (R[j] - Q[j])*S*dt + getLocalVol( (dt*j), S) * S * sqrt(dt) * rand.normal_box_muller();

            if (S >= op.barrier){
                S = 0;
                break;
            }
        }
        
        result += std::max(0.0, S - op.K)*exp(-R.back()*op.T);
    }
    
    return result/(double)N;
}





