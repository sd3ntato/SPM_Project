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

#ifndef metrics_computation_h
#define metrics_computation_h
#include "metrics_computation.hpp"
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

#ifndef utils_h
#define utils_h
#include "utils.h"
#endif

#include <limits>
#include <fstream>

using namespace std;

int main()
{
  int n_samples = 100;
  int n_trials = 1;
  int max_par_degree = 8;

  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv", n_samples);
  Matrix_wrapper dataset_n = normalize(dataset);
  cout << "dataset read" << endl;

  vector<int> Nrs = {2000};

  prepare_vectors_for_statistics();

  for (int i = 0; i < Nrs.size(); i++)
  {
    declare_vectors_and_matrices();

    double t0 = compute_sequential_time(Nr, n_samples, n_trials, dataset, dataset_n, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);

    //without fastflow (my pool)
    compute_statistics("none", times_for_each_Nr, speedups_for_each_Nr, scalabilities_for_each_Nr, efficiencies_for_each_Nr);

    // with fastflow parallel_for
    compute_statistics("parfor", times_for_each_Nr_parfor, speedups_for_each_Nr_parfor, scalabilities_for_each_Nr_parfor, efficiencies_for_each_Nr_parfor);

    // with my implementation of pool in fastflow
    compute_statistics("ff_pool", times_for_each_Nr_ff_pool, speedups_for_each_Nr_ff_pool, scalabilities_for_each_Nr_ff_pool, efficiencies_for_each_Nr_ff_pool);

    // with mdf model + ff_pool
    //compute_statistics("mdf", times_for_each_Nr_mdf, speedups_for_each_Nr_mdf, scalabilities_for_each_Nr_mdf, efficiencies_for_each_Nr_mdf);
  }

  plot_comparison("time", times_for_each_Nr[0], times_for_each_Nr_parfor[0], times_for_each_Nr_ff_pool[0], n_trials);
  plot_comparison("speedup", speedups_for_each_Nr[0], speedups_for_each_Nr_parfor[0], speedups_for_each_Nr_ff_pool[0], n_trials);
  plot_comparison("scalability", scalabilities_for_each_Nr[0], scalabilities_for_each_Nr_parfor[0], scalabilities_for_each_Nr_ff_pool[0], n_trials);
  plot_comparison("efficiency", efficiencies_for_each_Nr[0], efficiencies_for_each_Nr_parfor[0], efficiencies_for_each_Nr_ff_pool[0], n_trials);

  ofstream out_file("out");
  if (out_file.is_open())
  {
    out_file << Nrs[0] << " neurons " << endl;

    out_file << "seq time: " << 

    out_file << "time_none" << " " << "time_parfor" << " " << "time_ff_pool" << endl;

    for (int i = 1; i < max_par_degree + 1; i++)
    {
      out_file << times_for_each_Nr[0][i] << " " << times_for_each_Nr_parfor[0][i] << " " << times_for_each_Nr_ff_pool[0][i] << " ";
      out_file << endl;
    }

  }

  /*
  plot(Nrs, times_for_each_Nr, "times", n_trials, false);
  plot(Nrs, speedups_for_each_Nr, "speedups", n_trials, false);
  plot(Nrs, scalabilities_for_each_Nr, "scalabilities", n_trials, false);
  plot(Nrs, efficiencies_for_each_Nr, "efficiencies", n_trials, false);

  plot(Nrs, times_for_each_Nr_ff, "times", n_trials, true);
  plot(Nrs, speedups_for_each_Nr_ff, "speedups", n_trials, true);
  plot(Nrs, scalabilities_for_each_Nr_ff, "scalabilities", n_trials, true);
  plot(Nrs, efficiencies_for_each_Nr_ff, "efficiencies", n_trials, true);
  */

  cout << "\nftt!\n";
  return 0;
}
