#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef io_stream
#define io_stream
#include <iostream>
#endif


class ESN
{
public:
  ESN(int Nr = 100, int Ny = 1, int Nu = 1, float rho = 0.9, float r_d = 0.1);
  float rho;
  int Nr;              // Number recurrent units
  int Nu;              // Number inputs
  int Ny;              // Number outputs
  float r_density;     // density of recurrent matrix
  Matrix_wrapper Win;  // input-to-reservoir matrix
  Matrix_wrapper W;    //recurrent matrix
  Matrix_wrapper Wout; // readout
  Matrix_wrapper x;    //state

  Matrix_wrapper compute_state(Matrix_wrapper u);
  Matrix_wrapper compute_output(Matrix_wrapper u);
  Matrix_wrapper compute_output();
};
