#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef io_stream
#define io_stream
#include <iostream>
#endif

#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#include "ESN.h"
#include "math.h"
#include "matplotlibcpp.h"
#include <thread>
#include <chrono>

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
    map2(0, Nr, 0, Nr, ref(p), Dot_task(), W, x_old, x_rec);
    map2(0, Nr, 0, Nu, ref(p), Dot_task(), Win, u, x_in);

    // compute tanh(sum...)
    map1(0, Nr, 100, ref(p), Comp_state_task(), Nu, x_rec, x_in, Win, x);

    // z = P|x
    map2(0, Nr + 1, 0, Nr + 1, ref(p), Dot_task(), P, x, z);

    // y = Wout|x
    map2(0, Ny, 0, Nr + 1, ref(p), Dot_task(), Wout, x, y);

    // k_den = x.T | z
    p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});
    p.await_no_tasks_todo();

    k_den += l;

    // k = z/k_den
    map1(0, Nr + 1, 100, ref(p), Divide_by_const(), z, k_den, k);

    // Wold = ....
    map2(0, Ny, 0, Nr + 1, ref(p), Compute_new_Wout(), Wout, Wold, d, y, k);

    // P = ...
    map2(0,Nr+1,0,Nr+1, ref(p), Compute_new_P(), P, Pold, k, z, l);

    // map
    for (int i = 0; i < Ny; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Wold[i][j] = Wout[i][j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Pold[i][j] = P[i][j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      x_old[i] = x[i];
    }

    s = 0;
    for (int i = 0; i < Ny; i++)
    {
      s += pow(d[i] - y[i], 2);
    }
    //cout << sqrt(s) << " " << flush;
    errors_norms.push_back(s);

    cnt++;
  }

  cout << endl;
  errors_norms.erase(errors_norms.begin());
  plt::plot(errors_norms);
  plt::show();

  cout << "\nftt!\n";
  return 0;
}
