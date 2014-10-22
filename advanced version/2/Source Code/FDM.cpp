#include "FDM.h"

FDMBase::FDMBase(double _x_max , int _N, double _t_max, int _M, BlackScholesPDE* _pde) : x_max(_x_max), N(_N), t_max(_t_max), M(_M), pde(_pde){}


FDMImplicit_Second::FDMImplicit_Second(double _x_max , int _N, double _t_max, int _M, BlackScholesPDE* _pde) : FDMBase(_x_max, _N, _t_max, _M, _pde) {
        computeStepSize();
        setInitialCondition();                      // besides calling the constructor of base class FDMBase, it also calls the computerStepSize and setInitialCondition
}


void FDMImplicit_Second::computeStepSize(){     // in .h file, we have to specify that the computeStepSize is a virtual function, but in .cpp file, we don't have to do that!
        dx = x_max/(double)(N);
        N = N + 1;                                  // why we need N = N + 1????
        dt = t_max/(double)(M);
}

void FDMImplicit_Second::setInitialCondition(){

        double cur_spot = 0.0;
        old_result.resize(N, 0.0);                  // old_result, new_result and x_values are all vectors.
        new_result.resize(N, 0.0);
        x_values.resize(N, 0.0);

        for (int j=0; j<N; j++) {
                cur_spot = (double)(j)*dx;
                old_result[j] = pde->init_cond(cur_spot);// return option->pay_off->operator()(x); PayOffPut::operator()(const double& s ) const { return std:: max(K-s,0.0);}
        x_values[j] = cur_spot;
        }

        // Temporal settings
        prev_t = 0.0;
        cur_t = 0.0;
}

void FDMImplicit_Second::computeBoundaryCondition(){
        d[1] = d[1] + 2.0*lamda[1];
        mu[1] = mu[1] - lamda[1];
        lamda[N-1] = lamda[N-1] - mu[N-1];
        d[N-1] = d[N-1] + 2.0*mu[N-1];
}

//Implement Tridiagonal Solver
void FDMImplicit_Second::computeSpatialSlice(){

        lamda.resize(N, 0.0);                   // use N to define their length
        d.resize(N, 0.0);
        mu.resize(N, 0.0);

        //Construct Matrix A
        for (int j = 0; j < N; ++j) {
                alpha = (pde->diff_coeff(prev_t, x_values[j]) ) * dt/dx/dx;
                beta = (pde->conv_coeff(prev_t, x_values[j]) ) * dt/2.0/dx;

                lamda[j] = -alpha + beta;
                d[j] = 1.0 - (pde->zero_coeff(prev_t, x_values[j]) ) * dt + 2.0*alpha;
                mu[j] = -alpha - beta;
        }

        //Apply bounday conditions
        computeBoundaryCondition();

        //Solve Ax = b
        for (int j = 1; j < N; ++j) {
                d[j] = d[j] - lamda[j]/d[j-1]*mu[j-1];
                old_result[j] = old_result[j] - lamda[j]/d[j-1]*old_result[j-1];
        }

        new_result[N-1] = old_result[N-1]/d[N-1];

        for (int j = N-2; j >= 0; --j) {
                new_result[j] = (old_result[j] - mu[j]*new_result[j+1])/d[j];
        }
}

void FDMImplicit_Second::compute(std::ofstream &fdm_out){

        for (int j=0; j<N; j++) {
                fdm_out << "," << x_values[j];
        }
        fdm_out << std::endl;

        fdm_out << "T = 0";
        for (int j=0; j<N; j++) {
                fdm_out << "," << old_result[j];
        }
        fdm_out << std::endl;

        while(cur_t < t_max) {
                cur_t = prev_t + dt;
                fdm_out << "T = " << cur_t;

                computeSpatialSlice();

                for (int j=0; j<N; j++) {
                        fdm_out << "," << new_result[j];
                }
                fdm_out << std::endl;

                old_result = new_result;
                prev_t = cur_t;
        }

        std::cout.precision(8);
        std::cout << "When S = 100, Price of Put = " << old_result[(int)(100/dx)] << std::endl;

        fdm_out.close();
}

FDMImplicit_Third::FDMImplicit_Third(double _x_max , int _N, double _t_max, int _M, BlackScholesPDE* _pde) : FDMImplicit_Second(_x_max, _N, _t_max, _M, _pde) {}

void FDMImplicit_Third::computeBoundaryCondition(){
        d[2] = d[2] + 42.0/13.0*k[2] + 27.0/13.0*lamda[2];
        mu[2] = mu[2] - 15.0/13.0*lamda[2] - 32.0/13.0*k[2];
        v[2] = v[2] + lamda[2]/13.0 + 3.0/13.0*k[2];

        lamda[3] = lamda[3] + 27.0/13.0*k[3];
        d[3] = d[3] - 15.0/13.0*k[3];
        mu[3] = mu[3] + k[3]/13.0;

        lamda[N-3] = lamda[N-3] + v[N-3]/13.0;
        d[N-3] = d[N-3] - 15.0/13.0*v[N-3];
        mu[N-3] = mu[N-3] + 27.0/13.0*v[N-3];

        k[N-2] = k[N-2] + mu[N-2]/13.0 + 3.0/13.0*v[N-2];
        lamda[N-2] = lamda[N-2] - 15.0/13.0*mu[N-2] - 32.0/13.0*v[N-2];
        d[N-2] = d[N-2] + 27.0/13.0*mu[N-2] + 42.0/13.0*v[N-2];
}

//Implement Penta Solver
void FDMImplicit_Third::computeSpatialSlice(){

        lamda.resize(N, 0.0);
        d.resize(N, 0.0);
        mu.resize(N, 0.0);
        k.resize(N, 0.0);
        v.resize(N, 0.0);

        //Construct Matrix A
        for (int j = 0; j < N; ++j) {
                alpha = (pde->diff_coeff(prev_t, x_values[j]) )/12.0 * dt/dx/dx;
                beta = (pde->conv_coeff(prev_t, x_values[j]) )/12.0 * dt/dx;

                k[j] = alpha - beta;
                lamda[j] = -16.0*alpha + 8.0*beta;
                d[j] = 1.0 - (pde->zero_coeff(prev_t, x_values[j]) ) * dt + 30.0*alpha;
                mu[j] = -16.0*alpha - 8.0*beta;
                v[j] = alpha + beta;
        }

        //Apply bounday conditions
        computeBoundaryCondition();

        //Solve Ax = b
        for (int j = 1; j < N-1; ++j) {
                d[j] = d[j] - lamda[j]/d[j-1]*mu[j-1];
                mu[j] = mu[j] - lamda[j]/d[j-1]*v[j-1];
                old_result[j] = old_result[j] - lamda[j]/d[j-1]*old_result[j-1];

                lamda[j+1] = lamda[j+1] - k[j+1]/d[j-1]*mu[j-1];
                d[j+1] = d[j+1] - k[j+1]/d[j-1]*v[j-1];
                old_result[j+1] = old_result[j+1] - k[j+1]/d[j-1]*old_result[j-1];
        }

        d[N-1] = d[N-1] - lamda[N-1]/d[N-2]*mu[N-2];
        old_result[N-1] = old_result[N-1] - lamda[N-1]/d[N-2]*old_result[N-2];

        new_result[N-1] = old_result[N-1]/d[N-1];
        new_result[N-2] = (old_result[N-2] - mu[N-2]*new_result[N-1]) / d[N-2];

        for (int j = (int)N-3; j >= 0; --j) {
                new_result[j] = (old_result[j] - v[j]*new_result[j+2] - mu[j]*new_result[j+1]) / d[j];
        }

}



