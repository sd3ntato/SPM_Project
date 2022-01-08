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

  /* global parameters */
  int n_samples = 100; // number of training samples to be supplied
  int n_trials = 1;    // number of mesurements per congiguration
  int max_par_degree = 8;

  /* read dataset */
  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv", n_samples);
  Matrix_wrapper dataset_n = normalize(dataset); // in linear_algebra.cpp
  cout << "dataset read" << endl;

  /* numbers of neurons for each network at test. defines scale of the problem */
  vector<int> Nrs = {2000};

  /* conduct the experiments: 
  * for each network size instantiate the network and train it with different modalities,
  * then collect some stats on the training process,
  * In particular, will be collecting time, speedup, scalability, efficiency, each averaged over
  * n_trials trials.
  */
  for (int i = 0; i < Nrs.size(); i++)
  {
    prepare_vectors_for_statistics(); // prepares data structure in wich will put statistics
    declare_vectors_and_matrices();   // prepares data structure needed to define network and its functioning

    /* compute sequential time */
    double t0 = compute_sequential_time(Nr, n_samples, n_trials, dataset, dataset_n, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);

    //without fastflow (my pool)
    compute_statistics("none", times_none, speedups_none, scalabilities_none, efficiencies_none);
    dump_one("none_times", times_none, "none");

    // with fastflow parallel_for
    compute_statistics("parfor", times_parfor, speedups_parfor, scalabilities_parfor, efficiencies_parfor);
    dump_one("parfor_times", times_parfor, "parfor");

    // with my implementation of pool in fastflow
    compute_statistics("ff_pool", times_ff_pool, speedups_ff_pool, scalabilities_ff_pool, efficiencies_ff_pool);
    dump_one("ff_pool_times", times_ff_pool, "ff_pool");

    // with mdf model + ff_pool
    //compute_statistics("mdf", times_mdf, speedups_mdf, scalabilities_mdf, efficiencies_mdf);

    /* plot the collected stats and dump them into file "out" */
    //do_plots();
    //dump_statistics();
  }

  cout << "\ndone!\n";
  return 0;
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
