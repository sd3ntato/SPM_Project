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

#include "seq_functs.h"
#include "utils.h"

#include <limits>

using namespace std;

int main()
{
  int n_samples = 100;
  int n_trials = 1;
  int c_line_size = 128;
  int max_par_degree = 10;

  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv", n_samples);
  Matrix_wrapper dataset_n = normalize(dataset);
  cout << "dataset read" << endl;

  vector<int> Nrs = {100, 200, 1000};
  vector<double> *times_for_each_Nr = new vector<double>[Nrs.size()];
  vector<double> *speedups_for_each_Nr = new vector<double>[Nrs.size()];
  vector<double> *scalabilities_for_each_Nr = new vector<double>[Nrs.size()];
  vector<double> *efficiencies_for_each_Nr = new vector<double>[Nrs.size()];
  for (int i = 0; i < Nrs.size(); i++)
  {

    ESN n = ESN(Nrs[i], 4, 4);
    float **W = n.W.m;
    float **Win = n.Win.m;
    double t0 = compute_sequential_time(Nrs[i], n_samples, n_trials, dataset, dataset_n, W, Win);
    times_for_each_Nr[i] = compute_average_times(Nrs[i], n_samples, n_trials, max_par_degree, c_line_size, dataset, dataset_n, W, Win);
    speedups_for_each_Nr[i] = compute_speedups(times_for_each_Nr[i], t0);
    scalabilities_for_each_Nr[i] = compute_scalabilities(times_for_each_Nr[i]);
    efficiencies_for_each_Nr[i] = compute_effieciencies(speedups_for_each_Nr[i]);
  }

  plot(Nrs, times_for_each_Nr, "times", n_trials);
  plot(Nrs, speedups_for_each_Nr, "speedups", n_trials);
  plot(Nrs, scalabilities_for_each_Nr, "scalabilities", n_trials);
  plot(Nrs, efficiencies_for_each_Nr, "efficiencies", n_trials);

  cout << "\nftt!\n";
  return 0;
}
