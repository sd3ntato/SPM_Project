#ifndef io_stream
#define io_stream
#include <iostream>
#endif

#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef metrics_computation_h
#define metrics_computation_h
#include "metrics_computation.hpp"
#endif

#ifndef utils_h
#define utils_h
#include "utils.h"
#endif

#include "seq_functs.h"
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

  for (int i = 0; i < Nrs.size(); i++)
  {
    prepare_vectors_for_statistics();
    declare_vectors_and_matrices();

    double t0 = compute_sequential_time(Nr, n_samples, n_trials, dataset, dataset_n, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);

    //without fastflow (my pool)
    compute_statistics("none", times_none, speedups_none, scalabilities_none, efficiencies_none);

    // with fastflow parallel_for
    compute_statistics("parfor", times_parfor, speedups_parfor, scalabilities_parfor, efficiencies_parfor);

    // with my implementation of pool in fastflow
    compute_statistics("ff_pool", times_ff_pool, speedups_ff_pool, scalabilities_ff_pool, efficiencies_ff_pool);

    // with mdf model + ff_pool
    //compute_statistics("mdf", times_mdf, speedups_mdf, scalabilities_mdf, efficiencies_mdf);

    do_plots();
    dump_statistics();
  }

  /*
  plot(Nrs, times, "times", n_trials, false);
  plot(Nrs, speedups, "speedups", n_trials, false);
  plot(Nrs, scalabilities, "scalabilities", n_trials, false);
  plot(Nrs, efficiencies, "efficiencies", n_trials, false);

  plot(Nrs, times_ff, "times", n_trials, true);
  plot(Nrs, speedups_ff, "speedups", n_trials, true);
  plot(Nrs, scalabilities_ff, "scalabilities", n_trials, true);
  plot(Nrs, efficiencies_ff, "efficiencies", n_trials, true);
  */

  cout << "\nftt!\n";
  return 0;
}
