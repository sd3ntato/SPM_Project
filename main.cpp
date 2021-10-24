#include <iostream>
#include <random>
#include "ESN.h"

using namespace std;

int main() {

    Matrix_wrapper m1 = ones(3,3);
    Matrix_wrapper m2 = ones(3,3);
    Matrix_wrapper m3 = elementwise_sum(m1,m2);

    print_matrix(m3);



    cout << "\nftt!\n";
    return 0;
}