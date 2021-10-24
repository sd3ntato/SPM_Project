#include "ESN.h"

ESN::ESN(int Nr, int Ny, int Nu, float rho, float r_d) {
    this->Nr = Nr;
    this->Ny = Ny;
    this->Nu = Nu;
    this->rho = rho;
    this->r_density = r_d;

    this->W = build_sparse_contractive_matrix(Nr,Nr);
    this->Win = generate_random_sparse_matrix(Nr,Nu,0.7);

    this->x = zeros(Nr,1);
}

Matrix_wrapper ESN::compute_state(Matrix_wrapper u){
    u = vstack(u, ones(1,1));
    Matrix_wrapper z = dot( this->Win, u );
    return Matrix_wrapper(nullptr,0,0);
}
