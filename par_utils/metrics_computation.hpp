#include <vector>
#include "matplotlibcpp.h"
#include <string>

#ifndef parallel_functs_h
#define parallel_functs_h
#include "parallel_functs.h"
#endif

#ifndef utimer_cpp
#define utimer_cpp
#include "utimer.cpp"
#endif

using namespace std;

#define prepare_vectors_for_statistics()                                              \
  vector<double> *times_for_each_Nr = new vector<double>[Nrs.size()];                 \
  vector<double> *speedups_for_each_Nr = new vector<double>[Nrs.size()];              \
  vector<double> *scalabilities_for_each_Nr = new vector<double>[Nrs.size()];         \
  vector<double> *efficiencies_for_each_Nr = new vector<double>[Nrs.size()];          \
                                                                                      \
  vector<double> *times_for_each_Nr_parfor = new vector<double>[Nrs.size()];          \
  vector<double> *speedups_for_each_Nr_parfor = new vector<double>[Nrs.size()];       \
  vector<double> *scalabilities_for_each_Nr_parfor = new vector<double>[Nrs.size()];  \
  vector<double> *efficiencies_for_each_Nr_parfor = new vector<double>[Nrs.size()];   \
                                                                                      \
  vector<double> *times_for_each_Nr_ff_pool = new vector<double>[Nrs.size()];         \
  vector<double> *speedups_for_each_Nr_ff_pool = new vector<double>[Nrs.size()];      \
  vector<double> *scalabilities_for_each_Nr_ff_pool = new vector<double>[Nrs.size()]; \
  vector<double> *efficiencies_for_each_Nr_ff_pool = new vector<double>[Nrs.size()];  \
                                                                                      \
  vector<double> *times_for_each_Nr_mdf = new vector<double>[Nrs.size()];             \
  vector<double> *speedups_for_each_Nr_mdf = new vector<double>[Nrs.size()];          \
  vector<double> *scalabilities_for_each_Nr_mdf = new vector<double>[Nrs.size()];     \
  vector<double> *efficiencies_for_each_Nr_mdf = new vector<double>[Nrs.size()];

#define compute_statistics(directive, times, speedups, scalabilities, efficiencies)                                                                                                   \
  times[i] = compute_average_times(directive, Nr, n_samples, n_trials, max_par_degree, c_line_size, dataset, dataset_n, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y); \
  speedups[i] = compute_speedups(times[i], t0);                                                                                                                                       \
  scalabilities[i] = compute_scalabilities(times[i]);                                                                                                                                 \
  efficiencies[i] = compute_effieciencies(speedups[i]);

/************* FUNCTIONS TO COMPUTE OBSERVED METRICS *********************/

vector<double> compute_average_times(string ff, int Nr, int n_samples, int n_trials, int max_par_degree, int c_line_size, Matrix_wrapper dataset, Matrix_wrapper dataset_n, float **W, float **Win,
                                     float **Wout, float **Wold, float **P, float **Pold,
                                     float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{
  int Nu = 4;
  int Ny = 4;
  float l = 0.995;
  float nabla = 0.1;

  vector<double> times(max_par_degree + 1);
  times[0] = std::numeric_limits<double>::quiet_NaN();
  vector<double> err;
  for (int par_deg = 1; par_deg < max_par_degree + 1; par_deg++)
  {
    {
      vector<double> ts(n_trials);
      for (int i = 0; i < n_trials; i++)
      {
        utimer t(to_string(par_deg), &ts[i]);
        err = par_train(ff, par_deg, c_line_size, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);
        //plt::plot(err);
        //plt::title("err with " + to_string(Nr) + " neurons and par deg " + to_string(par_deg));
        //plt::save("imgs/" + to_string(Nr) + "-" + to_string(par_deg) + "-err");
        //plt::show();
      }
      times[par_deg] = mean(ts);
    }
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