#ifndef __CaseStudy2__FDM__
#define __CaseStudy2__FDM__

#include <iostream>
#include <vector>
#include "PDE.h"
#include <fstream>

class FDMBase {
        protected:
                BlackScholesPDE* pde;

                //X discretisation
                double x_max;      //Underlining [0.0 , x_max]
                int N;             //Number of points in x direction
                double dx;         //X step size
                std::vector<double> x_values; // Stores the coordinates of the x

                //Time discretisation
                double t_max;      //Time [0.0, t_max]
                int M;             //Number of points in t direction
                double dt;         //Time step size

                //Temporal directions
                double prev_t;     //Previous times
                double cur_t;      //Current times

                //Differencing coefficients
                double alpha, beta;

                //Storage
                std::vector<double> new_result;
                std::vector<double> old_result;

                //Constructor
                FDMBase(double _x_max , int _N, double _t_max, int _M, BlackScholesPDE* _pde);          // constructor only initialize x_max, N, t_max, M, pde!!!

                //
                virtual void computeStepSize() = 0;
                virtual void setInitialCondition() = 0;
                virtual void computeBoundaryCondition() = 0;
                virtual void computeSpatialSlice() = 0;

        public:
                //Carry out the actual time−stepping
                virtual void compute(std::ofstream &fdm_out) = 0;


};

class FDMImplicit_Second : public FDMBase{
        protected:
                std::vector<double> lamda, d, mu;

                virtual void computeStepSize();
                virtual void setInitialCondition();
                virtual void computeBoundaryCondition();
                virtual void computeSpatialSlice();

        public:
                FDMImplicit_Second(double _x_max , int _N, double _t_max, int _M, BlackScholesPDE* _pde);
                virtual void compute(std::ofstream &fdm_out);
};

class FDMImplicit_Third : public FDMImplicit_Second{
        protected:
                std::vector<double> k, v;

                virtual void computeBoundaryCondition();
                virtual void computeSpatialSlice();

        public:
                FDMImplicit_Third(double _x_max , int _N, double _t_max, int _M, BlackScholesPDE* _pde);
};


#endif
