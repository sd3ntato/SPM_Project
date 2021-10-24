#include <iostream>
#include <random>
#include "ESN.h"

using namespace std;

int main() {

    Matrix_wrapper m1 = zeros(2,3);
    m1.m[0][0] = 1.0; m1.m[0][1] = 2.0; m1.m[0][2] = 1.0;
    m1.m[1][0] = 2.0; m1.m[1][1] = 1.0; m1.m[1][2] = 1.0;
    Matrix_wrapper m2 = zeros(3,2);
    m2.m[0][0] = 1.0; m2.m[0][1] = 1.0;
    m2.m[1][0] = 1.0; m2.m[1][1] = 0.0;
    m2.m[2][0] = 0.0; m2.m[2][1] = 0.0;

    Matrix_wrapper m3 = dot(m1,m2);

    print_matrix(m1);
    print_matrix(m2);
    print_matrix(m3);

    Matrix_wrapper t1 = generate_random_sparse_matrix(300,300,0.7);
    Matrix_wrapper t2 = generate_random_sparse_matrix(300,300,0.7);
    Matrix_wrapper t3 = dot(t1,t2);

    cout << "\nftt!\n";
    return 0;
}