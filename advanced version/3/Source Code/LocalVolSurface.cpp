#include "LocalVolSurface.h"

LocalVolSurfaceBuilder::LocalVolSurfaceBuilder(std::vector<double> _SGrid, std::vector<double> _TGrid, std::vector<double> _S_R, std::vector<double> _S_Q, std::vector<OptionParameterSet> options) : SGrid(_SGrid), TGrid(_TGrid), S_R(_S_R), S_Q(_S_Q){
    
    OptionPremiumSurface = new std::vector< std::vector<double> >( (int)TGrid.size(), std::vector<double>( (int)SGrid.size(), 0) );
    LocalVolSurface = new std::vector< std::vector<double> >( (int)TGrid.size(), std::vector<double>( (int)SGrid.size(), 0) );
    
}

void LocalVolSurfaceBuilder::Build(double S0, dcmplx (*_CharacterFunc)(dcmplx, OptionParameterSet, ModelParameterSet), ModelParameterSet mp, std::string title){
    
    int FFTStep = 512;
    double noUse[FFTStep];
    double alpha = 2;
    double eta = 0.15;
    
    double s_size = 1;
    double t_size = 1e-3;
    
    double C_C, C_U, C_D, C_L, C_R;
    
    std::cout << "Build Surfaces..." << std::endl;
    
    std::ofstream outputFile( title + "_LocalVolSurface.csv");
    std::ofstream outputFile2( title + "_CallSurface.csv");
    
    for (int i = 0; i < SGrid.size(); i++){
        
		for (int j = 0; j < TGrid.size(); j++){
            
            std::cout << "(" << i << "," << j << ")" << std::endl;
            
            //Build Parameter
            OptionParameterSet op(S0, TGrid[j], S_R[j], S_Q[j], SGrid[i]);
            FTPricing fft(op, mp, _CharacterFunc);
			(*OptionPremiumSurface)[j][i] = fft.FastFT(noUse, SGrid[i], eta, alpha, FFTStep); //C_T_K
			C_C = (*OptionPremiumSurface)[j][i];
			
            OptionParameterSet op2(S0, TGrid[j] - t_size, S_R[j], S_Q[j], SGrid[i]);
            fft.op = op2;
			C_L = fft.FastFT(noUse, SGrid[i], eta, alpha, FFTStep); //C_T-dt
			
            OptionParameterSet op3(S0, TGrid[j] + t_size, S_R[j], S_Q[j], SGrid[i]);
            fft.op = op3;
			C_R = fft.FastFT(noUse, SGrid[i], eta, alpha, FFTStep);  //C_T+dt
			
            fft.op = op;
			C_D = fft.FastFT(noUse, SGrid[i] - s_size, eta, alpha, FFTStep);  //C_T_K-dK
			C_U = fft.FastFT(noUse, SGrid[i] + s_size, eta, alpha, FFTStep);  //C_T_K+dK
            
            double num1 = (C_R-C_L)/2.0/t_size;
            double num2 = S_Q[j]*C_C;
            double num3 = SGrid[i]*( S_R[j]-S_Q[j] )*( C_U-C_D )/2.0/s_size;
            double num4 = (SGrid[i]*SGrid[i]*(C_U-2.0*C_C+C_D)/s_size/s_size);
            double result = sqrt( fabs( 2.0 * (num1+num2+num3)/num4 ) );
			
			(*LocalVolSurface)[j][i] = result;
            
            outputFile2 << (*OptionPremiumSurface)[j][i] << ",";
            outputFile << (*LocalVolSurface)[j][i] << ",";
		}
        
        outputFile << std::endl;
        outputFile2 << std::endl;
    }
    
    outputFile.close();
    outputFile2.close();
    
    std::cout << "Surfaces Built!" << std::endl;
}

std::vector< std::vector<double> > LocalVolSurfaceBuilder::getLocalVolSurface(){
    std::vector< std::vector<double> > output = (*LocalVolSurface);
    return output;
}

LocalVolSurfaceBuilder::~LocalVolSurfaceBuilder(){
    delete OptionPremiumSurface;
    delete LocalVolSurface;
}
