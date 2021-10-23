#include <iostream>
#include <random>
#include "ESN.h"

using namespace std;

int main() {

    ESN n =  ESN() ;

    print_matrix(n.W);

    print_matrix(n.Win);

    cout << "\nftt!\n";
    return 0;
}