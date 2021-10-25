#include <iostream>
#include <random>
#include "ESN.h"
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include "matplotlib-cpp/matplotlibcpp.h"

//BTCUSDT-1m-data.csv


using namespace std;
namespace plt = matplotlibcpp;


Matrix_wrapper read_dataset(string filename);

int main() {
    //plt::plot({1,3,2,4});
    //plt::show();
    //print_matrix(dataset.get_line(10));

    
    Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv");

    int Nr = 100;
    int Nu = 4;
    int Ny = 4;
    float lambda = 0.1;
    float nabla = 0.01;
    int i=0;

    ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
    Matrix_wrapper P = eye(Nr)*nabla;
    n.Wout = zeros(Ny,Nu); Matrix_wrapper y = zeros(4,1);
    while(i<dataset.n2){
        Matrix_wrapper candle = dataset.get_line(i);

        i++;
    }


    cout << "\nftt!\n";
    return 0;
}

Matrix_wrapper read_dataset(string filename){
    string line;
    ifstream myfile (filename);
    Matrix_wrapper dataset = Matrix_wrapper(nullptr,0,0); //container for the final dataset
    if( myfile.is_open() ){
        getline(myfile,line); // discard the first line as it contains intestation.
        while ( getline(myfile,line) ){ //get a line from the csv
            istringstream iss(line); //turn it into this thing
            string s; 
            string* ss = new string[12]; // each line is divided into 12 tokens (we are interested in token 1 to 5)
            int i=0;
            while( getline(iss, s, ',') ){
                ss[i] = s; //store the tokens in the apposite array
                i++;
            }
            float* numbers = new float[4]; // keep only ohlc values
            for(int i=0;i<5;i++){
                numbers[i] = stof( ss[i+1] );
            }
            Matrix_wrapper n = from_array(numbers,4); //dichiarazione va spostata fuori
            dataset = vstack(dataset,n); // put the collected line into the dataset
        }
        myfile.close();
    }
    return dataset;
}