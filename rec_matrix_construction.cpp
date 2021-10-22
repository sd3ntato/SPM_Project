#include <iostream>
#include <Eigen/Dense>
 
using namespace Eigen;
using namespace std;

float** zeros(int n1, int n2){
    float **m = new float*[n1];
    for(int i=0;i<n1;i++){
        m[i] = new float[n2];
    }
    return m;
}

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

int main()
{
    int n1 = 100;
    int n2 = 100;
    auto m = zeros(n1,n2);

    MatrixXd mat = (MatrixXd::Random(n1,n2).array() > 0.9).cast<double>() * MatrixXd::Random(n1,n2).array();
    auto rho_act = abs(mat.eigenvalues().array())[0];

    cout << rho_act << endl;
    if(rho_act<0.1){
        cout<< "recurrent matrix creation failed bc spectral radius too small";
        return -1;
    }

    mat = mat / rho_act;

    for(int i=0;i<n1;i++){
        for (int j = 0; j < n2; j++){
            m[i][j] = mat(i,j);
        }
    }
    
    print_matrix(m,n1,n2);

    std::cout<< abs(mat.eigenvalues().array())[0] << std::endl;

}