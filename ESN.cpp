#include "funct.h"
#include "ESN.h"

ESN::ESN(int Nr, int Ny, int Nu, float rho, float r_d) { // Constructor with parameters
    this->Nr = Nr;
    this->Ny = Ny;
    this->Nu = Nu;
    this->rho = rho;
    this->r_density = r_d;

    this->W = build_sparse_contractive_matrix(Nr,Nr);
    this->Win = generate_random_sparse_matrix(Nr,Nu,0.7);
}
