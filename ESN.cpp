#include "ESN.h"
#include <iostream>

ESN::ESN(int Nr, int Ny, int Nu, float rho, float r_d)
{
  this->Nr = Nr;
  this->Ny = Ny;
  this->Nu = Nu;
  this->rho = rho;
  this->r_density = r_d;

  this->W = build_sparse_contractive_matrix(Nr, Nr);
  std::cout << "generated recurrent matrix \n"
            << std::flush;
  this->Win = generate_random_sparse_matrix(Nr, Nu + 1, 0.7);
  std::cout << "generated input matrix \n"
            << std::flush;
  this->Wout = generate_random_sparse_matrix(Ny, Nr + 1, 0.7);
  std::cout << "generated output matrix \n"
            << std::flush;

  this->x = zeros(Nr, 1);
}

Matrix_wrapper ESN::compute_state(Matrix_wrapper u)
{
  u = vstack(u, ones(1, 1));
  Matrix_wrapper p1 = (this->Win | u);
  Matrix_wrapper p2 = (this->W | this->x);
  Matrix_wrapper z = p1 + p2;
  Matrix_wrapper out = elementwise_tanh(z);
  free_matrices({this->x});
  this->x = out;
  free_matrices({z, u, p1, p2});
  return copy(this->x);
}

Matrix_wrapper ESN::compute_output(Matrix_wrapper u)
{
  Matrix_wrapper x = this->compute_state(u);
  Matrix_wrapper v = vstack(x, ones(1, 1));
  Matrix_wrapper out = this->Wout | v;
  free_matrices({x, v});
  return out;
}

Matrix_wrapper ESN::compute_output()
{
  Matrix_wrapper v = vstack(this->x, ones(1, 1));
  Matrix_wrapper out = this->Wout | v;
  free_matrices({v});
  return out;
}