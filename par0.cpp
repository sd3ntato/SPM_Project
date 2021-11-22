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

#ifndef ESN_h
#define ESN_h
#include "ESN.h"
#endif

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

  ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
  float **W = n.W.m;
  float **Win = n.Win.m;
  
  int c_line_size=128;

  for (int par_deg = 1; par_deg < 10; par_deg++)
  {
    vector<double> errors_norms = par_train(par_deg, c_line_size, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l,
                                            W, Win);

    errors_norms.erase(errors_norms.begin());
    plt::plot(errors_norms);
    plt::show();
  }

  cout << "\nftt!\n";
  return 0;
}
