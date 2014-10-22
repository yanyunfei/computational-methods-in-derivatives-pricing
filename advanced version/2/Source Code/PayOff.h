#ifndef __CaseStudy2__PayOff__
#define __CaseStudy2__PayOff__

#include <iostream>
#include <algorithm>

//Payoff Base class
class PayOff {
        public:
                PayOff();
                virtual double operator() (const double& S) const = 0;          //why this form?                putting a const here to show this function could not change the variables in class PayOff object!
                virtual ~PayOff() {};                                           //why virtual descructor?
};

//Put Payoff
class PayOffPut : public PayOff {
        public:
                PayOffPut(const double& K);
                virtual ~PayOffPut () {};                                               

                double operator() (const double& S) const;                      // overload ooperator function in base class

        private:
                double K; // Strike
};

#endif
