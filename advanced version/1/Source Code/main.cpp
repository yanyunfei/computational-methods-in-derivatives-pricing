#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "math.h"
#include "OptionPricing.h"
#include "CharacteristicFuncLib.h"

using namespace std;

//Function FTTest - Conduct tests for FFT and FrFFT. Results will be exported to csv file
void FTTest(ModelParameterSet modelParas, dcmplx (*charFunct)(dcmplx, ModelParameterSet), vector<double> K, vector<double> alpha, vector<int> N, double eta, double lambda);

//Function COSTest - Conduct tests for FFT and FrFFT. Results will be exported to csv file
void COSTest(ModelParameterSet modelParas, dcmplx (*charFunct)(dcmplx, ModelParameterSet), vector<double> K, vector<int> N, vector<double> ab);

//Function writeVectorToCSV - write a vector to csv file
template<class TYPE> void writeVectorToCSV(ofstream& file, vector<TYPE>& v, string name);
//Function ArrayToVector - covert array to vector
template<class TYPE> vector<TYPE> ArrayToVector(TYPE* in, int size);



int main(int argc, const char * argv[])
{
    //Model Parameters
    double S0 = 1800;
    double T = 0.5;
    double r = 0.0025;
    double q = 0.0203;
    double sigma = 0.3;
    ModelParameterSet modelParas(S0, T, r, q, sigma);
    
    //Set Test Cases
    double eta = 0.25;
    double lambda = 0.1; //only for FrFFT method
    
    double _K[] = {900, 1100, 1300, 1500};
    int _N[] = {(int)pow(2, 7), (int)pow(2, 8), (int)pow(2, 9), (int)pow(2, 10), (int)pow(2, 12), (int)pow(2, 14)};
    double _Alpha[] = {-2, -5, -10, -20};
    
    double _AB[] = {2 , 5, 10, 20};
    
    //Convert to vector
    vector<double> K = ArrayToVector(_K, sizeof(_K)/sizeof(_K[0]));
    vector<int> N = ArrayToVector(_N, sizeof(_N)/sizeof(_N[0]));
    vector<double> Alpha = ArrayToVector(_Alpha, sizeof(_Alpha)/sizeof(_Alpha[0]));
    vector<double> AB = ArrayToVector(_AB, sizeof(_AB)/sizeof(_AB[0]));

    
    //Function pointer to the desired Character Function
    //Later, it can be easily extented to other models
    dcmplx (*gbmCharFunct)(dcmplx, ModelParameterSet) = GBMCharacterFunc;
    
    //Conduct Test for FFT and FrFFT - results will be exported to FT_Test.csv
    FTTest(modelParas, gbmCharFunct, K, Alpha, N, eta, lambda);
    
    //Conduct Test for COS - results will be exported to COS_Test.csv
    COSTest(modelParas, gbmCharFunct, K, N, AB);

    return 0;
}


//Function FTTest - Conduct tests for FFT and FrFFT. Results will be exported to csv file
void FTTest(ModelParameterSet modelParas, dcmplx (*charFunct)(dcmplx, ModelParameterSet), vector<double> K, vector<double> alpha, vector<int> N, double eta, double lambda){
    
    //Track Execution Time
    clock_t timer;
    double sum1 = 0, sum2 = 0;
    double num = (double) alpha.size();
    
    //Output File Location.
    ofstream outputFile("FT_Test.csv");
    
    outputFile << "FFT and FrFFT Test V2" << endl << endl << ",";
    
    for (int i = 0 ; i < alpha.size(); ++i)
        outputFile << "  By FFT,  By FrFFT,";
    outputFile << endl;
    
    FTPricing ftPricingAgent(modelParas, charFunct);
    
    //For each K
    for (vector<double>::iterator k_iter = K.begin(); k_iter != K.end(); ++k_iter) {
        outputFile << "Set K = " << *k_iter << "," << endl << ",";
        writeVectorToCSV(outputFile, alpha, "Alpha = ");
        
        //For each N
        for (vector<int>::iterator n_iter = N.begin(); n_iter != N.end(); ++n_iter){
            sum1 = 0;
            sum2 = 0;
            outputFile << "N = " << *n_iter << ",";
            
            //For each alpha
            for (vector<double>::iterator a_iter = alpha.begin(); a_iter != alpha.end(); ++a_iter){
                //FFT
                double* prices = new double[*n_iter];
                timer = clock();
                outputFile << ftPricingAgent.FastFT(prices, *k_iter, eta, *a_iter, *n_iter) << ",";
                sum1 += (double) ( (clock() - timer)/(double)CLOCKS_PER_SEC );
                //FrFFT
                timer = clock();
                outputFile << ftPricingAgent.FrFFT(prices, *k_iter, eta, *a_iter, *n_iter, lambda) << ",";
                sum2 += (double) ( (clock() - timer)/(double)CLOCKS_PER_SEC );
                delete[] prices;
            }
            outputFile << ", N = " << *n_iter << ": FFT Avg Time = " << sum1/num*1000 << " ms. FrFFT Avg Time = " << sum2/num*1000 << " ms"<< endl;
        }
        outputFile << endl;
    }
    
    cout << "FT test is completed. CSV is generated!\n";
    
    outputFile.close();
}

//Function COSTest - Conduct tests for FFT and FrFFT. Results will be exported to csv file
void COSTest(ModelParameterSet modelParas, dcmplx (*charFunct)(dcmplx, ModelParameterSet), vector<double> K, vector<int> N, vector<double> ab){
    
    //Track Execution Time
    clock_t timer;
    double sum = 0;
    double num = (double) ab.size();
    
    //Output File Location.
    ofstream outputFile("COS_Test.csv");
    
    outputFile << "COS Test V2" << endl << endl;
    
    COSPricing cosPricingAgent(modelParas, charFunct);
    
    //For each K
    for (vector<double>::iterator k_iter = K.begin(); k_iter != K.end(); ++k_iter) {
        outputFile << "Set K = " << *k_iter << "," << endl << ",";
        
        for (int i = 0 ; i < ab.size(); ++i)
            outputFile << " [" << -ab[i] << " : " <<  ab[i] << "] ,";
        outputFile << endl;
        
        //For each N
        for (vector<int>::iterator n_iter = N.begin(); n_iter != N.end(); ++n_iter){
            sum = 0;
            outputFile << "N = " << *n_iter << ",";
            
            //For each [a,b]
            for (vector<double>::iterator ab_iter = ab.begin(); ab_iter != ab.end(); ++ab_iter){
                timer = clock();
                outputFile << cosPricingAgent.COS(*k_iter, *n_iter, -*ab_iter, *ab_iter) << ",";
                sum += (double) ( (clock() - timer)/(double)CLOCKS_PER_SEC );
            }
            
            outputFile << ",N = " << *n_iter << ": Time = " << sum/num*1000 << " ms" << endl;
        }
        
        outputFile << endl;
    }
    
    cout << "COS test is completed. CSV is generated!\n";
    
    outputFile.close();
}



//Helper Functions
template<class TYPE>
void writeVectorToCSV(ofstream& file, vector<TYPE>& v, string name){
    for (typename vector<TYPE>::iterator iter = v.begin(); iter != v.end(); ++iter){
        file << name << *iter << ", ,";
    }
    file << endl;
}

template<class TYPE>
vector<TYPE> ArrayToVector(TYPE* in, int size){
    return vector<TYPE>(in, in + size);
}




