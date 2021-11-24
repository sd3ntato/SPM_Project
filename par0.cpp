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

#ifndef utimer_cpp
#define utimer_cpp
#include "utimer.cpp"
#endif

#include <limits>

#include "matplotlibcpp.h"

//#include "matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace std;

int main()
{
  int n_samples = 300;
  int n_trials = 2;
  int c_line_size = 128;

  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv", n_samples);
  Matrix_wrapper dataset_n = normalize(dataset);
  cout << "dataset read" << endl;

  vector<double> times1 = compute_average_times(2000, n_samples, n_trials, c_line_size, dataset, dataset_n);
  vector<double> times2 = compute_average_times(1000, n_samples, n_trials, c_line_size, dataset, dataset_n);

  plt::named_plot("1000 neurons", times1, "-x"); // TODO error bar con varianza
  plt::named_plot("2000 neurons", times2, "-x");
  plt::legend();
  plt::show();

  cout << "\nftt!\n";
  return 0;
}
