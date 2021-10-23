#include <random>
#include <iostream>
#include <Eigen/Dense>
#include <assert.h> 


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
        for(int j=0;j<n2;j++){
            m[i][j]=0;
        }
    }
    return m;
}

bool contains(int* l, int n, int t){
    for(int i=0;i<n;i++){
        if(l[i]==t){
            return true;
        }
    }
    return false;
}

float** generate_random_sparse_matrix(int n1, int n2, float d, float interval_start=-1, float interval_stop=1){

    // matrix with all zeros
    float** m = zeros(n1,n2);
        
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

