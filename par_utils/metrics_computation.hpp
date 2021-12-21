#include "matplotlibcpp.h"

#ifndef vector_h
#define vector_h
#include <vector>
#endif

#ifndef string_h
#define string_h
#include <string>
#endif

#ifndef parallel_functs_h
#define parallel_functs_h
#include "parallel_functs.h"
#endif

#ifndef utimer_cpp
#define utimer_cpp
#include "utimer.cpp"
#endif

using namespace std;

/* prepares data structures in wich statistics will be stored */
#define prepare_vectors_for_statistics() \
  vector<double> times_none;             \
  vector<double> speedups_none;          \
  vector<double> scalabilities_none;     \
  vector<double> efficiencies_none;      \
                                         \
  vector<double> times_parfor;           \
  vector<double> speedups_parfor;        \
  vector<double> scalabilities_parfor;   \
  vector<double> efficiencies_parfor;    \
                                         \
  vector<double> times_ff_pool;          \
  vector<double> speedups_ff_pool;       \
  vector<double> scalabilities_ff_pool;  \
  vector<double> efficiencies_ff_pool;   \
                                         \
  vector<double> times_mdf;              \
  vector<double> speedups_mdf;           \
  vector<double> scalabilities_mdf;      \
  vector<double> efficiencies_mdf;

/* directive param int compute_statistics is a string, it can be:
* - none for no ff
* - parfor for ff_parallel_for
* - ff_pool for my impl. of thread pool using ff
* - mdf for mdf model orchestrating the ff_pool object
*/
#define compute_statistics(directive, times, speedups, scalabilities, efficiencies)                                                                                   \
  times = compute_average_times(directive, Nr, n_samples, n_trials, max_par_degree, dataset, dataset_n, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y); \
  speedups = compute_speedups(times, t0);                                                                                                                             \
  scalabilities = compute_scalabilities(times);                                                                                                                       \
  efficiencies = compute_effieciencies(speedups);

#define dump_statistics()                                                                 \
  ofstream out_file("out");                                                               \
  if (out_file.is_open())                                                                 \
  {                                                                                       \
    out_file << Nr << " neurons " << endl;                                                \
                                                                                          \
    out_file << "seq time: " << t0 << endl;                                               \
                                                                                          \
    dump("time", times_none, times_parfor, times_ff_pool);                                \
    dump("speedup", speedups_none, speedups_parfor, speedups_ff_pool);                    \
    dump("scalability", scalabilities_none, scalabilities_parfor, scalabilities_ff_pool); \
    dump("efficiency", efficiencies_none, efficiencies_parfor, efficiencies_ff_pool);     \
  }

#define dump(what, stat_none, stat_parfor, stat_ff_pool)                                \
  out_file << what << "_none"                                                           \
           << " "                                                                       \
           << what << "_parfor"                                                         \
           << " "                                                                       \
           << what << "_ff_pool" << endl;                                               \
                                                                                        \
  for (int i = 1; i < max_par_degree + 1; i++)                                          \
  {                                                                                     \
    out_file << stat_none[i] << " " << stat_parfor[i] << " " << stat_ff_pool[i] << " "; \
    out_file << endl;                                                                   \
  }                                                                                     \
  out_file << endl;

#define do_plots()                                                                                           \
  plot_comparison("time", times_none, times_parfor, times_ff_pool, n_trials);                                \
  plot_comparison("speedup", speedups_none, speedups_parfor, speedups_ff_pool, n_trials);                    \
  plot_comparison("scalability", scalabilities_none, scalabilities_parfor, scalabilities_ff_pool, n_trials); \
  plot_comparison("efficiency", efficiencies_none, efficiencies_parfor, efficiencies_ff_pool, n_trials);

void plot_comparison(string what, vector<double> quantity_none, vector<double> quantity_parfor, vector<double> quantity_ff_pool, int n_trials)
{
  plt::figure_size(1200, 780);
  plt::named_plot("none", quantity_none, "-x");
  plt::named_plot("parfor", quantity_parfor, "-x");
  plt::named_plot("ff_pool", quantity_ff_pool, "-x");

  plt::xlabel("number of workers");
  plt::ylabel("average " + what + " on " + to_string(n_trials) + " executions");

  plt::title("comparison " + what);
  plt::legend();
  plt::save("imgs/comparison_" + what);

  //plt::show();
}

void plot(vector<int> Nrs, vector<double> *v, string s, int n_trials, bool ff)
{
  plt::figure_size(1200, 780);
  for (int i = 0; i < Nrs.size(); i++)
  {
    plt::named_plot(to_string(Nrs[i]) + " neurons", v[i], "-x");
  }
  plt::xlabel("number of workers");
  plt::ylabel("average " + s + " on " + to_string(n_trials) + " executions");
  if (ff)
  {
    plt::title(s + " vs. number of workers with ff");
    plt::legend();
    plt::save("imgs/" + s + "-nworkers_ff");
  }
  else
  {
    plt::title(s + " vs. number of workers");
    plt::legend();
    plt::save("imgs/" + s + "-nworkers");
  }
  //plt::show();
}

/************* FUNCTIONS TO COMPUTE OBSERVED METRICS *********************/

/* the first param in compute_avergage_times is a string, it can be:
* - none for no ff - uses my implementation of thread pool -
* - parfor for ff_parallel_for
* - ff_pool for my impl. of thread pool using ff
* - mdf for mdf model orchestrating the ff_pool object
*/
vector<double> compute_average_times(string ff, int Nr, int n_samples, int n_trials, int max_par_degree, Matrix_wrapper dataset, Matrix_wrapper dataset_n, float **W, float **Win,
                                     float **Wout, float **Wold, float **P, float **Pold,
                                     float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{
  // parameters initialization
  int Nu = 4;
  int Ny = 4;
  float l = 0.995;
  float nabla = 0.1;

  vector<double> times(max_par_degree + 1); // this is used to store the times taken to compute sequential training

  times[0] = std::numeric_limits<double>::quiet_NaN(); // just useful for later when plotting

  vector<double> err; // used to store error at each iteration

  // compute training n_trials times for each configuration of parallelism degree
  for (int par_deg = 1; par_deg < max_par_degree + 1; par_deg++)
  {
    vector<double> ts(n_trials); // this stores the training time for each trial, later used to compute a mean

    // repeat the training process n_trials times and collect the training times
    for (int i = 0; i < n_trials; i++)
    {
      utimer t(to_string(par_deg), &ts[i]); // store the training time into apposit placeholder

      // compute training 
      err = par_train(ff, par_deg, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);
      //plt::plot(err);
      //plt::title("err with " + to_string(Nr) + " neurons and par deg " + to_string(par_deg));
      //plt::save("imgs/" + to_string(Nr) + "-" + to_string(par_deg) + "-err");
      //plt::show();
    }

    //average the training times
    times[par_deg] = mean(ts);
  }
  return times;
}

vector<double> compute_speedups(vector<double> times, double t0)
{
  int n = times.size();
  vector<double> speedups(n);
  for (int i = 0; i < n; i++)
  {
    speedups[i] = t0 / times[i];
  }
  return speedups;
}

vector<double> compute_scalabilities(vector<double> times)
{
  int n = times.size();
  vector<double> scalabilities(n);
  for (int i = 0; i < n; i++)
  {
    scalabilities[i] = times[1] / times[i];
  }
  return scalabilities;
}

vector<double> compute_effieciencies(vector<double> speedups)
{
  int n = speedups.size();
  vector<double> efficiencies(n);
  for (int i = 0; i < n; i++)
  {
    efficiencies[i] = speedups[i] / i;
  }
  return efficiencies;
}