#include <iostream>
#include <random>
#include "ESN.h"

using namespace std;

int main() {

    Matrix_wrapper m1 = generate_random_sparse_matrix(3,1,1);
    Matrix_wrapper m2 = generate_random_sparse_matrix(3,1,1);

    Matrix_wrapper m3 = vstack(m1,m2);

    print_matrix(m1);
    print_matrix(m2);
    print_matrix(m3);


    cout << "\nftt!\n";
    return 0;
}