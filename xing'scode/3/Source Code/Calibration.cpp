#include "Calibration.h"

OptionLoader::OptionLoader(std::string _file): file(_file) {};

void OptionLoader::LoadOption(int method){
    std::ifstream myfile(file.c_str());
    
    std::string line;
    std::string cell;
    
    while(myfile){
        
        OptionParameterSet option;
        
        std::getline(myfile, line);
        std::stringstream lineStream(line);
        
        std::getline(lineStream, cell, ',');
        option.T = atof(cell.c_str())/365.0;
        
        std::getline(lineStream, cell, ',');
        option.type = cell[0];
        
        std::getline(lineStream, cell, ',');
        option.S0 = atof(cell.c_str());
        
        std::getline(lineStream, cell, ',');
        option.K = atof(cell.c_str());
        
        std::getline(lineStream, cell, ',');
        option.bid = atof(cell.c_str());
        
        std::getline(lineStream, cell, ',');
        option.ask = atof(cell.c_str());
        
        option.mid = (option.bid + option.ask)/2.0;
        
        std::getline(lineStream, cell, ',');
        option.r = atof(cell.c_str())/100.0;
        
        std::getline(lineStream, cell, ',');
        option.q = atof(cell.c_str())/100.0;
        
        if (method == 0)
            option.weight = 1.0;
        else
            option.weight = 1.0/fabs(option.ask - option.bid);
        
        //Choose Out of Money
        if ( option.type == 'C') {
            if ( option.S0 < option.K)
                Options.push_back(option);
        }
        else {
            if ( option.S0 > option.K)
                Options.push_back(option);
        }
    }
    
    if (method == 1) {
        double sum = 0.0;
        for (std::vector<OptionParameterSet>::iterator iter = Options.begin(); iter != Options.end(); ++iter)
            sum += iter->weight;
        
        for (std::vector<OptionParameterSet>::iterator iter = Options.begin(); iter != Options.end(); ++iter)
            iter->weight /= sum;
    }
    
    myfile.close();
}


TermStructureBuilder::TermStructureBuilder(double Smin, double Smax, int numOfSStep, double Tmin, double Tmax, int numOfTStep) {
    
    S_StepSize = (Smax - Smin)/(numOfSStep-1);
    for (int i = 0; i < numOfSStep; ++i)
        SGrid.push_back(Smin + i*S_StepSize);
    
    T_StepSize = (Tmax - Tmin)/(numOfTStep-1);
    for (int i = 0; i < numOfTStep; ++i)
        TGrid.push_back(Tmin + i*T_StepSize);
    
};
    
void TermStructureBuilder::BuildRQ(std::vector<OptionParameterSet> data, std::vector<double>& R, std::vector<double>& Q, int method){
    if (method == 0)
        RQBuilder_step(data, R, Q);
    
    if (method == 1)
        RQBuilder_linear(data, R, Q);
    
    if (method == 2)
        RQBuilder_constant(data, R, Q);
}

//linear interpolation
void TermStructureBuilder::RQBuilder_linear(std::vector<OptionParameterSet> options, std::vector<double>& R, std::vector<double>& Q) {
    
    std::sort(options.begin(), options.end());
    std::vector<RQPair> RQ;
    
    double tempT = -1.0;
    for (int i = 0; i < options.size(); ++i) {
        if ( fabs(tempT-options[i].T) > 1e-5 ){
            RQ.push_back( RQPair(options[i].r, options[i].q, options[i].T) );
            tempT = options[i].T;
        }
    }
    
    int i = 0;
    for (; i < TGrid.size(); ++i){
        if (TGrid[i] <= RQ[0].T) {
            R.push_back(RQ[0].R);
            Q.push_back(RQ[0].Q);
        }
        else
            break;
    }
    
    int count = 1;
    for (; i < TGrid.size(); ++i){
        if (TGrid[i] <= RQ[count].T) {
            R.push_back( (TGrid[i] - RQ[count-1].T)/T_StepSize * ( RQ[count].R - RQ[count-1].R ) + RQ[count-1].R  );
            Q.push_back( (TGrid[i] - RQ[count-1].T)/T_StepSize * ( RQ[count].Q - RQ[count-1].Q ) + RQ[count-1].Q  );
        }
        else{
            if (TGrid[i] > RQ.back().T)
                break;
            else
                count++;
        }
    }
    
    for (; i < TGrid.size(); ++i){
        R.push_back( RQ.back().R );
        Q.push_back( RQ.back().Q );
    }
    
    return;
}

//piecewise constant interpolation
void TermStructureBuilder::RQBuilder_step(std::vector<OptionParameterSet> options, std::vector<double>& R, std::vector<double>& Q) {
    
    std::sort(options.begin(), options.end());
    std::vector<RQPair> RQ;
    
    double tempT = -1.0;
    for (int i = 0; i < options.size(); ++i) {
        if ( fabs(tempT-options[i].T) > 1e-5 ){
            RQ.push_back( RQPair(options[i].r, options[i].q, options[i].T) );
            tempT = options[i].T;
        }
    }
    
    int count = 0;
    for (int i = 0; i < TGrid.size(); ++i){
        while ( count != (RQ.size()-1) ){
            if ( (TGrid[i] > RQ[count].T) )
                count++;
            else
                break;
        }
        
        R.push_back(RQ[count].R);
        Q.push_back(RQ[count].Q);
    }
    
    return;
}

//Constant
void TermStructureBuilder::RQBuilder_constant(std::vector<OptionParameterSet> options, std::vector<double>& R, std::vector<double>& Q) {
    
    double r = 0.03;
    double q = 0.0;
    
    for (int i = 0; i < TGrid.size(); ++i){
        R.push_back(r);
        Q.push_back(q);
    }
    
    return;
}


double heston_objectiveFunct(const gsl_vector *v, void *params) //Heston objective func in simplex
{
    int FFTStep = 512;
    
	OptionParameterSet *op = (OptionParameterSet *)params;
    
	double output = 0, temp;
	
    int i, NUM=0;
	
    //get the size of dataset
    while((*op).type == 'C' || (*op).type=='P'){
		op++;
		NUM++;
	}
    
    double sigma, kappa, theta, rou, nu;
    sigma = gsl_vector_get(v, 0);
	kappa = gsl_vector_get(v, 1);
	theta = gsl_vector_get(v, 2);
	rou = gsl_vector_get(v, 3);
	nu = gsl_vector_get(v, 4);
    
	op = (OptionParameterSet *)params;
    
    ModelParameterSet mp(sigma, kappa, theta, rou, nu, 0.0, 0.0);
    
	for (i=0; i<NUM; i++){
        
        FTPricing pricingAgent(op[i], mp, cf_heston);
        
        double noUse[FFTStep];
        double alpha = (op[i].type == 'C') ? 2.0:-2.0;
        double eta = 0.15;
        
        temp = pricingAgent.FastFT(noUse, op[i].K, eta, alpha, FFTStep);
        
        double a = (temp -  op[i].mid)*(temp -  op[i].mid);
        output +=  a * op[i].weight/NUM;
        
	}
    
	output += ( std::max(0.0,rou-1) + std::max(0.0, -rou-1) + std::max(0.0, -nu) + std::max(0.0, -theta) + std::max(0.0, -sigma) + std::max(0.0, sigma*sigma-2.0*kappa*theta) )*1e8;
    
	return output;
}

double vgsa_objectiveFunct(const gsl_vector *v, void *params) //VGSA objective func in simplex
{
    int FFTStep = 512;
    
	OptionParameterSet *op = (OptionParameterSet *)params;
    
	double output = 0, temp;
	
    int i, NUM=0;
	
    //get the size of dataset
    while((*op).type == 'C' || (*op).type=='P'){
		op++;
		NUM++;
	}
    
    double sigma, nu, theta, kappa, eta, lambda;
	sigma = gsl_vector_get(v, 0);
	nu = gsl_vector_get(v, 1);
	theta = gsl_vector_get(v, 2);
	kappa = gsl_vector_get(v, 3);
	eta = gsl_vector_get(v, 4);
	lambda = gsl_vector_get(v, 5);
    
	op = (OptionParameterSet *)params;
    
    ModelParameterSet mp(sigma, kappa, theta, 0.0, nu, eta, lambda);
    
	for (i=0; i<NUM; i++){
        
        FTPricing pricingAgent(op[i], mp, cf_vgsa);
        
        double noUse[FFTStep];
        double alpha = (op[i].type == 'C') ? 2:-2;
        double eta = 0.15;
        temp = pricingAgent.FastFT(noUse, op[i].K, eta, alpha, FFTStep);
        
        double a = (temp -  op[i].mid)*(temp -  op[i].mid);
        output +=  a * op[i].weight/NUM;
	}
    
    output += ( std::max(0.0, -sigma) + std::max(0.0, -eta) + std::max(0.0, -lambda) + std::max(0.0, -kappa) )*1e8;
    
	return output;
}

Calibrator::Calibrator() {};

std::vector<double> Calibrator::Calibrate(std::vector<double> initialVector, OptionParameterSet options[], double (*objectiveFunction)(const gsl_vector *v, void *params), std::ofstream &output){
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
    
	size_t iter = 0;
	int status;
	double size;
    
    std::vector<double> result;
	
    int numOfVariable = (int) initialVector.size();
    
    //Starting point
    x = gsl_vector_alloc(numOfVariable);
    
    for (int i = 0; i < numOfVariable; i++) {
        gsl_vector_set (x, i, initialVector[i]);
    }
    
    minex_func.f = objectiveFunction;
    ss = gsl_vector_alloc(numOfVariable);
    minex_func.n = numOfVariable;
    s = gsl_multimin_fminimizer_alloc (T, numOfVariable);
    
	//Set initial step sizes to 1
	gsl_vector_set_all (ss, 1);
    
	//Initialize method and iterate
	minex_func.params = options;
    
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	do{
        
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
            break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-3); //1e-3 for Heston, 1e-4 for VGSA for better result
        
		if (status == GSL_SUCCESS){
            
            std::cout << "A solution is found!" << std::endl;
            output << "Final Result:"<< std::endl;
            
            for (int i = 0; i < numOfVariable; ++i)
                result.push_back( gsl_vector_get (s->x, i) );
		}
        
        //Print current step info
        std::cout << "size = " << size << ", f() = " << s->fval << ": ";
        output << "size = " << size << ",f() = " << s->fval << ",";
        
        for (int i = 0; i < numOfVariable; ++i) {
            double val = gsl_vector_get (s->x, i);
            output << val << ",";
            std::cout << val << ", ";
        }
        
        output << std::endl;
        std::cout << std::endl;
        
        iter++;
	}
	while (status == GSL_CONTINUE && iter < 50000 );
	
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
    
    return result;
}













