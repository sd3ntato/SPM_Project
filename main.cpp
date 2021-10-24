#include <iostream>
#include <random>
#include "ESN.h"
#include <math.h>

using namespace std;

int main() {

    auto m1 = ones(3,3);
    auto m2 = ones(3,3);

    //print_matrix( m1+m2 );
    //print_matrix( m1|m2 );


    auto n = ESN();
    for(int i=0;i<1000;i++){
        n.compute_output(ones(1,1));
    }
    print_matrix( n.compute_output(ones(1,1)) );

    cout << "\nftt!\n";
    return 0;
}