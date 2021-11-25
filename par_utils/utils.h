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

void plot(vector<int> Nrs, vector<double> *v, string s, int n_trials)
{
  plt::figure_size(1200, 780);
  for (int i = 0; i < Nrs.size(); i++)
  {
    plt::named_plot(to_string(Nrs[i]) + " neurons", v[i], "-x");
  }
  plt::xlabel("number of workers");
  plt::ylabel("average " + s + " on " + to_string(n_trials) + " executions");
  plt::title(s + " vs. number of workers");
  plt::legend();
  plt::save("imgs/" + s + "-nworkers");
  plt::show();
}