#include <iostream>
#include <algorithm>
#include "Calibration.h"
#include "LocalVolSurface.h"
#include "OptionPricingSimulation.h"

ModelParameterSet GetModelParaSet(std::vector<double> parameters);

void CaseStudy3(std::string file, std::string title);
void SubQuestion1(std::string file, int method, std::string title);


int main(int argc, const char * argv[]){
    
    //File Data_Day1 is the option data for 3/12
    std::string file1 = "/Users/Kiminian/Desktop/C++/CaseStudy3/CaseStudy3/Data_Day1.csv";
    CaseStudy3(file1, "Day1");
    
    //File Data_Day2 is the option data for 3/13
    std::string file2 = "/Users/Kiminian/Desktop/C++/CaseStudy3/CaseStudy3/Data_Day2.csv";
    CaseStudy3(file2, "Day2");
    
    return 0;
}

void CaseStudy3(std::string file, std::string title){
    
    //Equal Weight
    SubQuestion1(file, 0, title);
    //Inverse Weight
    SubQuestion1(file, 1, title);
    
}

void SubQuestion1(std::string file, int method, std::string title){
    
    //For heston model: sigma, kappa, theta, rou, nu
    std::vector<double> hestonInitialVector = {1.0143, 4.9549, 0.0562, -0.6552, 0.057};
    
    //For VGSA model: sigma, nu, theta, kappa, eta, lambda
    std::vector<double> vgsaInitialVector = {0.18, 0.15, -0.09, 7.6, 2.4, 9.6};
    
    std::cout << "Loading file " << file << "...";
    OptionLoader loader(file);
    loader.LoadOption(method);
    std::string output = (method == 0) ? "_EqualWeight" : "_InverseSpreadWeight";
    
    //Construct term structure
    std::vector<double> R, Q;
    TermStructureBuilder termStruct(1700, 2300, 61, 165.0/365.0, 365.0/365.0, 21);
    termStruct.BuildRQ(loader.Options, R, Q, 0);
     
    //Heston Model
    std::cout << "Calibrate Heston Model...";
    Calibrator HestonCalibrator;
    std::ofstream outputFile( ("Heston_Calibration_" + title + output + ".csv") );
    std::vector<double> Heston_ParameterSet = HestonCalibrator.Calibrate(hestonInitialVector, &(loader.Options[0]), heston_objectiveFunct, outputFile);
    outputFile.close();
    
    std::cout << "--Construct localvol surface!";
    ModelParameterSet heston_paras = GetModelParaSet(Heston_ParameterSet);
    LocalVolSurfaceBuilder heston_builder(termStruct.SGrid , termStruct.TGrid, R, Q, loader.Options);
    heston_builder.Build(loader.Options[0].S0, cf_heston, heston_paras, ( "Heston_" + title + output) );
    std::vector< std::vector<double> > heston_localvol = heston_builder.getLocalVolSurface();
    
    //VGSA model:
    std::cout << "Calibrate VGSA Model...";
    Calibrator VGSACalibrator;
    std::ofstream outputFile2( ("VGSA_Calibration_" + title + output + ".csv") );
    std::vector<double> VGSA_ParameterSet = VGSACalibrator.Calibrate(vgsaInitialVector, &(loader.Options[0]), vgsa_objectiveFunct, outputFile2);
    outputFile2.close();
    
    std::cout << "--Construct localvol surface!";
    ModelParameterSet vgsa_paras = GetModelParaSet(VGSA_ParameterSet);
    LocalVolSurfaceBuilder vgsa_builder(termStruct.SGrid , termStruct.TGrid, R, Q, loader.Options);
    vgsa_builder.Build(loader.Options[0].S0, cf_vgsa, vgsa_paras, ( "VGSA_" + title + output) );
    std::vector< std::vector<double> > vgsa_localvol = vgsa_builder.getLocalVolSurface();
    
    //Simulation
    std::ofstream outputFile3( ("Simulation_" + title + output + ".csv") );
    
    std::vector<double> R2, Q2;
    TermStructureBuilder termStruct2(1700, 2000, 50, 0, 0.5, 100);
    termStruct2.BuildRQ(loader.Options, R2, Q2, 0);
    
    OptionParameterSet op2(1845.59, 0.5, 0.0, 0.0, 1850);
    op2.barrier = 1950;
    
    std::cout << "Heston:" << std::endl; outputFile3 << "Heston:" << std::endl;
    Heson_UOC H_UOC(op2, heston_paras);
    double p = H_UOC.Simulate(10000, 0.5/100, 100, R2, Q2);
    std::cout << p << std::endl; outputFile3 << p << std::endl;
    
    std::cout << "VGSA:" << std::endl; outputFile3 << "VGSA:" << std::endl;
    VGSA_UOC V_UOC(op2, vgsa_paras);
    p = V_UOC.Simulate(10000, 0.5/100.0, 100, R2, Q2);
    std::cout << p << std::endl; outputFile3 << p << std::endl;
     
    std::cout << "LocalVol-Heston:" << std::endl; outputFile3 << "LocalVol-Heston:" << std::endl;
    LocalVol_UOC L_H_UOC(op2, heston_paras);
    L_H_UOC.SetLocalVol(heston_localvol, termStruct.TGrid, termStruct.SGrid);
    p = L_H_UOC.Simulate(10000, 0.5/100, 100, R2, Q2);
    std::cout << p << std::endl; outputFile3 << p << std::endl;
    
    std::cout << "LocalVol-VGSA:" << std::endl; outputFile3 << "LocalVol-VGSA:" << std::endl;
    LocalVol_UOC L_V_UOC(op2, vgsa_paras);
    L_V_UOC.SetLocalVol(vgsa_localvol, termStruct.TGrid, termStruct.SGrid);
    p = L_V_UOC.Simulate(10000, 0.5/100, 100, R2, Q2);
    std::cout << p << std::endl; outputFile3 << p << std::endl;
    outputFile3.close();
    
}


ModelParameterSet GetModelParaSet(std::vector<double> parameters){
    ModelParameterSet output;
    
    if (parameters.size() == 5) {
        output.sigma = parameters[0];
        output.kappa = parameters[1];
        output.theta = parameters[2];
        output.rou = parameters[3];
        output.nu = parameters[4];
    }
    else{
        output.sigma = parameters[0];
        output.nu = parameters[1];
        output.theta = parameters[2];
        output.kappa = parameters[3];
        output.eta = parameters[4];
        output.lambda = parameters[5];
    }
    
    return output;
}



