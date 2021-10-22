#include <iostream>
#include "funct.h"
 
using namespace std;

int main(int argc, char** argv){
    int n1, n2;
    n1 = stoi(argv[1]);
    n2 = stoi(argv[2]);

    float** m = build_sparse_contractive_matrix(n1,n2);
        
    print_matrix(m,n1,n2);
    return 0;
}