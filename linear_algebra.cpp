#include "linear_algebra.h"
#include <random>
#include <iostream>
#include <Eigen/Dense>
#include <assert.h> 


using namespace Eigen;
using namespace std;

Matrix_wrapper::Matrix_wrapper(float** m, int n1, int n2){
    this->m = m;
    this->n1 = n1;
    this->n2 = n2;
}

void print_matrix(Matrix_wrapper mat){
    float** m = mat.m;
    int n1 = mat.n1;
    int n2 = mat.n2;
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

Matrix_wrapper zeros(int n1, int n2){
    float **m = new float*[n1];
    for(int i=0;i<n1;i++){
        m[i] = new float[n2];
        for(int j=0;j<n2;j++){
            m[i][j]=0;
        }
    }
    Matrix_wrapper mat = Matrix_wrapper(m,n1,n2);
    return mat;
}

bool contains(int* l, int n, int t){
    for(int i=0;i<n;i++){
        if(l[i]==t){
            return true;
        }
    }
    return false;
}

Matrix_wrapper generate_random_sparse_matrix(int n1, int n2, float d, float interval_start, float interval_stop){

    // matrix with all zeros
    Matrix_wrapper mat = zeros(n1,n2);
    float** m = mat.m;
        
    // compute #elements to be filled out
    int n = ceil( n1 * n2 * d ) ; 

    // generate a list of n random values using uniform distribution on floats
    float* numbers = new float[n];
    int* positions = new int[n1 * n2];
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd());
    uniform_real_distribution<> fdis(interval_start, interval_stop);
    uniform_int_distribution<> idis(0, n1*n2-1);
    for(int i=0; i<n; i++){
        numbers[i] = fdis(gen);// decido che numero voglio inserire

        /* per la posizione del numero: non posso scegliere la stessa posizione due volte,
        // che significa che la lista positions non puo' avere duplicati.
        // modo molto stupido per inserire senza fare duplicati: prima di inserire controllo se 
        // il numero che sto inserendo l'ho gia messo. Se l'ho gia messo non va bene quindi ne prendo un altro
        */
        int p = idis(gen); // prendo un numero a caso
        while( contains(positions,n1*n2,p) ){ //controllo se gia sta nella lista
            p = idis(gen); // se ho verificato che quello che volevo inserire gia sta nella lista ne prendo un altro
        } // alla fine sono sicuro che not constains(positions,p) quindi posso inserire tranquillamente p nella lista
        positions[i] = p; 
    }

    // inserisco i numeri che ho scelto, nelle relative posizioni all'interno della matrice
    for(int i=0;i<n;i++){
        int p = positions[i];
        int row = (int) floor(p/n2);
        int col = p % n2;
        m[row][col] = numbers[i];
    }
    return mat;
}

Matrix_wrapper build_sparse_contractive_matrix(int n1, int n2){
    MatrixXd mat = (MatrixXd::Random(n1,n2).array() > 0.9).cast<double>() * MatrixXd::Random(n1,n2).array();
    auto rho_act = abs(mat.eigenvalues().array())[0];

    if(rho_act<0.1){
        cout<< "recurrent matrix creation failed bc spectral radius too small";
        throw;
    }

    mat = mat / rho_act;

    Matrix_wrapper mm = zeros(n1,n2);
    float** m = mm.m;
    for(int i=0;i<n1;i++){
        for (int j = 0; j < n2; j++){
            m[i][j] = mat(i,j);
        }
    }
    return mm;
}

Matrix_wrapper dot(Matrix_wrapper m1, Matrix_wrapper m2){
    if(m1.n2 != m2.n1){
        throw "invalid dimensions!";
    }

    Matrix_wrapper res = zeros(m1.n1, m2.n2);
    for(int i=0; i<m1.n1; i++){ // scorre le righe della prima matrice
        for(int j=0; j<m2.n2; j++){ // scorren le colonne della seconda matrice
            int s = 0; // accumulatore
            for(int cnt=0; cnt<m1.n2; cnt++){ //contatore per scorrere riga-colonna
                s = s + (m1.m[i][cnt] * m2.m[cnt][j] );
            }
            res.m[i][j] = s;
        }
    }

    return res;
}

