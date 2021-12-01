#include <vector>
#include "matplotlibcpp.h"
#include <string>

namespace plt = matplotlibcpp;

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

tuple<float **, float **, float **, float **> prepare_matrices(int Nr, int Ny, float **Wout, float **Wold, float **P, float **Pold, float nabla)
{
  Wout = zeros(Ny, Nr + 1, Wout);
  Wold = zeros(Ny, Nr + 1, Wold);
  P = resetP(Nr + 1, Nr + 1, P, nabla);
  Pold = resetP(Nr + 1, Nr + 1, P, nabla);
  return make_tuple(Wout, Wold, P, Pold);
}

tuple<float *, float *, float *, float *, float *, float *, float *> prepare_vectors(int Nr, float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{
  x = zeros(Nr + 1, x);
  x_rec = zeros(Nr + 1, x_rec);
  x_in = zeros(Nr + 1, x_in);
  x_old = zeros(Nr + 1, x_old);

  k = zeros(Nr + 1, k);
  z = zeros(Nr + 1, z);
  y = zeros(4, y);
  return make_tuple(x, x_rec, x_in, x_old, k, z, y);
}

float compute_error(float *d, float *y, int Ny)
{
  float s = 0;
  for (int i = 0; i < Ny; i++)
  {
    s += pow(d[i] - y[i], 2);
  }
  return s;
}