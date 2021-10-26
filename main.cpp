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

    cout<< "reading dataset...";
    Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv");
    Matrix_wrapper dataset_n = normalize(dataset);
    cout<<"dataset read"<<endl;

    int Nr = 100;
    int Nu = 4;
    int Ny = 4;
    float lambda = 0.995;
    float nabla = 0.1;
    int i=0;

    vector<double> errors_norms {};

    ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
    n.Win = n.Win / 10;
    Matrix_wrapper P = eye(Nr+1)* (1/nabla);
    n.Wout = zeros(Ny,Nr+1); Matrix_wrapper y = zeros(4,1);
    while(i<1000){
        
        Matrix_wrapper u = dataset_n.get_line(i).transpose(); // shape(Nu,1)
        Matrix_wrapper d = dataset.get_line(i+1).transpose(); //shape(Ny,1)
        
        Matrix_wrapper x = vstack( n.compute_state( u ) , ones(1,1) ) ; //shape(Nr,1)
        //print_matrix(n.x);

        Matrix_wrapper y = n.compute_output( ); //shape(Ny,1)
        //print_matrix(y);

        Matrix_wrapper psi = d - y; // shape(Ny,1)
        //print_matrix(psi); //errore

        Matrix_wrapper zeta = P | x ; //shape(Nr,1)
        //print_matrix(zeta);

        float k_den = lambda +  ( x.transpose() | zeta ).to_float() ;
        //cout << k_den <<endl; 

        Matrix_wrapper k = zeta / k_den; //shape(Nr,1)
        //print_matrix(k);

        P = ( P - ( k | zeta.transpose() ) ) * (1/lambda) ;  //shape(Nr,Nr)
        //print_matrix(P);

        n.Wout = n.Wout + ( psi | k.transpose() ); // shape(Ny,Nr)
        //print_matrix(n.Wout);

        errors_norms.push_back( norm(psi) ) ;
        cout<< errors_norms[i] << " " ;
        i++;

    }

    cout<< endl;
    plt::plot(errors_norms);
    plt::show();    

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