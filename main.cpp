#include <iostream>
#include <random>
#include "ESN.h"
#include "funct.h"

using namespace std;

int main() {

    ESN n =  ESN() ;
    print_matrix(n.Win,n.Nr,n.Nu);

    cout << "\nftt!\n";
    return 0;
}