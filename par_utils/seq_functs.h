
#ifndef ESN_h
#define ESN_h
#include "ESN.h"
#endif

#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#include "math.h"

vector<double> seq_train(int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                         int Nr, int Nu, int Ny, float nabla, float l,
                         float **W, float **Win)
{
  float **Wout = zeros(Ny, Nr + 1).m;
  float **Wold = zeros(Ny, Nr + 1).m;
  float **P = (eye(Nr + 1) * (1 / nabla)).m;
  float **Pold = (eye(Nr + 1) * (1 / nabla)).m;

  float *x = zeros(1, Nr + 1).m[0];
  float *x_rec = zeros(1, Nr + 1).m[0];
  float *x_in = zeros(1, Nr + 1).m[0];
  float *x_old = zeros(1, Nr + 1).m[0];

  float *k = zeros(1, Nr + 1).m[0];
  float *z = zeros(1, Nr + 1).m[0];
  float *y = zeros(1, 4).m[0];

  vector<double> error_norms{};
  int i = 0;

  x[Nr] = 1.0;
  x_old[Nr] = 1.0;
  float *u;
  float *d;

  float k_den;
  float s;
  while (i < n_samples)
  {

    u = dataset_n.m[i];
    d = dataset.m[i + 1];

    for (int i = 0; i < Nr; i++)
    {
      x_rec[i] = dot(0, Nr, W[i], x_old);
    }

    for (int i = 0; i < Nr; i++)
    {
      x_in[i] = dot(0, Nu, Win[i], u);
    }

    for (int i = 0; i < Nr; i++)
    {
      x[i] = tanh(x_rec[i] + x_in[i] + Win[i][Nu]);
    }
    x[Nr] = 1.0;

    for (int i = 0; i < Nr + 1; i++)
    {
      z[i] = dot(0, Nr + 1, P[i], x);
    }

    for (int i = 0; i < Ny; i++)
    {
      y[i] = dot(0, Nr + 1, Wout[i], x);
    }

    k_den = l + dot(0, Nr + 1, x, z);

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      k[i] = z[i] / k_den;
    }

    // map
    for (int i = 0; i < Ny; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;
      }
    }

    // map
    for (int i = 0; i < Ny; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Wold[i][j] = Wout[i][j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Pold[i][j] = P[i][j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      x_old[i] = x[i];
    }

    s = 0;
    for (int i = 0; i < Ny; i++)
    {
      s += pow(d[i] - y[i], 2);
    }
    error_norms.push_back(sqrt(s));

    i++;
  }
  return error_norms;
}

double compute_sequential_time(int Nr, int n_samples, int n_trials, Matrix_wrapper dataset, Matrix_wrapper dataset_n, float **W, float **Win)
{
  int Nu = 4;
  int Ny = 4;
  float l = 0.995;
  float nabla = 0.1;

  vector<double> ts(n_trials);
  for (int i = 0; i < n_trials; i++)
  {
    utimer t("seq", &ts[i]);
    seq_train(n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win);
  }

  return mean(ts);
}