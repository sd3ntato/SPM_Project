#include <random>
#include <iostream>
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;


void print_matrix(float** m, int n1, int n2){
    for (int i = 0; i < n1; i++)
    {
        cout << "[ ";
        for (int ii = 0; ii < n2; ii++)
        {
            cout<< m[i][ii]<<" , ";
        }
        cout<<"]\n";
    }
    cout<<"\n";
}

float** zeros(int n1, int n2){
    float **m = new float*[n1];
    for(int i=0;i<n1;i++){
        m[i] = new float[n2];
    }
    return m;
}

float** generate_random_sparse_matrix(int n1, int n2, float d){

    // matrix with all zeros
    float** m = zeros(n1,n2);
        
    // compute #elements to be filled out
    int n = ceil( n1 * n2 * d ) ; 

    // generate a list of n random values using uniform distribution on floats
    float* numbers = new float[n];
    int* positions_1 = new int[n];
    int* positions_2 = new int[n];
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd());
    uniform_real_distribution<> fdis(1.0, 2.0);
    uniform_int_distribution<> idis1(0, n1-1);
    uniform_int_distribution<> idis2(0, n2-2);
    for(int i=0; i<n; i++){
        numbers[i] = fdis(gen);
        positions_1[i] = idis1(gen); // non funziona quando ricapita lo stesso numero
        positions_2[i] = idis2(gen);
    }

    for(int i=0;i<n;i++){
        m[positions_1[i]][positions_2[i]] = numbers[i];
    }

    return m;
}

float** build_sparse_contractive_matrix(int n1, int n2){
    MatrixXd mat = (MatrixXd::Random(n1,n2).array() > 0.9).cast<double>() * MatrixXd::Random(n1,n2).array();
    auto rho_act = abs(mat.eigenvalues().array())[0];

    if(rho_act<0.1){
        cout<< "recurrent matrix creation failed bc spectral radius too small";
        throw;
    }

    mat = mat / rho_act;

    auto m = zeros(n1,n2);
    for(int i=0;i<n1;i++){
        for (int j = 0; j < n2; j++){
            m[i][j] = mat(i,j);
        }
    }
    return m;
}

