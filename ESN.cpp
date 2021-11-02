#include "ESN.h"

ESN::ESN(int Nr, int Ny, int Nu, float rho, float r_d) {
    this->Nr = Nr;
    this->Ny = Ny;
    this->Nu = Nu;
    this->rho = rho;
    this->r_density = r_d;

    this->W = build_sparse_contractive_matrix(Nr,Nr);
    this->Win = generate_random_sparse_matrix(Nr,Nu+1,0.7);
    this->Wout = generate_random_sparse_matrix(Ny,Nr+1,0.7);

    this->x = zeros(Nr,1);
}

Matrix_wrapper ESN::compute_state(Matrix_wrapper u){
    u = vstack(u, ones(1,1));
    Matrix_wrapper z = (this->Win | u) + (this->W | this->x);
    Matrix_wrapper out = elementwise_tanh( z );
    free_matrices({this->x});
    this->x = out;
    free_matrices({z});
    return copy(this->x);
}

Matrix_wrapper ESN::compute_output(Matrix_wrapper u){
    return this->Wout | vstack( this->compute_state(u), ones(1,1) );
}

Matrix_wrapper ESN::compute_output(){
    return this->Wout | vstack( this->x, ones(1,1) );
}