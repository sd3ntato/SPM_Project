#ifndef io_stream
#define io_stream
#include <iostream>
#endif

#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#ifndef parallel_functs_h
#define parallel_functs_h
#include "parallel_functs.h"
#endif

#include "ESN.h" // do not need ifndef bc only included here
#include "matplotlibcpp.h"

//#include "matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace std;

Matrix_wrapper read_dataset(string filename, int n_samples);

int main()
{
  int n_samples = 1000;

  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv", n_samples);
  Matrix_wrapper dataset_n = normalize(dataset);
  cout << "dataset read" << endl;

  int Nr;
  string input;
  cout << "insert number of recurrent neurons: \n";
  cin >> input;
  Nr = stoi(input);

  int Nu = 4;
  int Ny = 4;
  float l = 0.995;
  float nabla = 0.1;
  int cnt = 0;

  vector<double> errors_norms{};

  ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
  float **W = n.W.m;
  float **Win = n.Win.m;
  float **Wout = zeros(Ny, Nr + 1).m;
  float **Wold = zeros(Ny, Nr + 1).m;
  float **P = (eye(Nr + 1) * (1 / nabla)).m;
  float **Pold = (eye(Nr + 1) * (1 / nabla)).m;

  float *x = zeros(1, Nr + 1).m[0];
  float *x_rec = zeros(1, Nr + 1).m[0];
  float *x_in = zeros(1, Nr + 1).m[0];
  float *x_old = zeros(1, Nr + 1).m[0];
  x_old[Nr] = 1.0;
  float *k = zeros(1, Nr + 1).m[0];
  float *z = zeros(1, Nr + 1).m[0];
  float *y = zeros(1, 4).m[0];
  float *u;
  float *d;

  float k_den;
  float s;

  Pool p(7);

  while (cnt < n_samples)
  {

    u = dataset_n.m[cnt];
    d = dataset.m[cnt + 1];

    // compute rec and in parts of state update
    parallel_matrix_dot_vector(0, Nr, 0, Nr, ref(p), W, x_old, x_rec);
    parallel_matrix_dot_vector(0, Nr, 0, Nu, ref(p), Win, u, x_in);
    p.barrier();

    // compute tanh(sum...)
    map1(0, Nr, ref(p), Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);
    p.barrier();

    // z = P|x
    parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, ref(p), P, x, z);
    p.barrier();

    // k_den = x.T | z , y = Wout|x
    parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, ref(p), Wout, x, y);
    p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});
    p.barrier();

    k_den += l;

    // k = z/k_den
    map1(0, Nr + 1, ref(p), Divide_by_const(), z, k_den, k);
    p.barrier();

    // Wold = ....
    map2(0, Ny, 0, Nr + 1, -1, ref(p), Compute_new_Wout(), Wout, Wold, d, y, k);
    p.barrier();

    // P = ...
    map2(0, Nr + 1, 0, Nr + 1, -1, ref(p), Compute_new_P(), P, Pold, k, z, l);
    p.barrier();

    s = 0;
    for (int i = 0; i < Ny; i++)
    {
      s += pow(d[i] - y[i], 2);
    }
    errors_norms.push_back(sqrt(s));

    cnt++;
  }

  cout << endl;
  errors_norms.erase(errors_norms.begin());
  plt::plot(errors_norms);
  plt::show();

  cout << "\nftt!\n";
  return 0;
}
