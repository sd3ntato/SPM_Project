#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#ifndef ESN_h
#define ESN_h
#include "ESN.h"
#endif

#include "math.h"

template <typename T, typename... Args>
void map1(int begin, int end, Pool &p, T task, Args... args)
{
  int n_points = end - begin;
  int diff = max((int)floor(n_points / p.n_workers), 128 / 4);
  if (diff == 0)
    diff = n_points;
  int start, stop;
  for (int i = begin; i < end; i += diff)
  {
    start = i;
    stop = min(i + diff, end);
    p.submit({new T(start, stop, args...)});
  }
}

template <typename T, typename... Args>
void map2(int begin, int end, int start, int stop, int diff, Pool &p, T task, Args... args)
{
  if (diff == -1)
  {
    int n_points = end - begin;
    if (n_points <= p.n_workers)
      diff = n_points;
    else
      diff = max((int)floor(n_points / p.n_workers), 128 / 4);
  }
  int i0, ii;
  for (int i = begin; i < end; i += diff)
  {
    i0 = i;
    ii = min(i + diff, end);
    //            horizontally   vertically
    //            |           |  |     |
    p.submit({new T(start, stop, i0, ii, args...)});
  }
}

void parallel_matrix_dot_vector(int row_start, int row_stop, int col_start, int col_stop, int diff, Pool &p, float **M, float *v, float *r)
{
  map2(row_start, row_stop, col_start, col_stop, diff, p, Multiple_Dot_task(), M, v, r);
  p.barrier();
}

vector<double> par_train(int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
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
  Pool p(par_degree);
  int cnt = 0;

  x[Nr] = 1.0;
  x_old[Nr] = 1.0;
  float *u;
  float *d;

  float k_den;
  float s;

  while (cnt < n_samples)
  {

    u = dataset_n.m[cnt];
    d = dataset.m[cnt + 1];

    // compute rec and in parts of state update
    parallel_matrix_dot_vector(0, Nr, 0, Nr, c_line_size / 4, ref(p), W, x_old, x_rec);
    parallel_matrix_dot_vector(0, Nr, 0, Nu, c_line_size / 4, ref(p), Win, u, x_in);
    p.barrier();

    // compute tanh(sum...)
    map1(0, Nr, ref(p), Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);
    p.barrier();

    // z = P|x
    parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, c_line_size / 4, ref(p), P, x, z);
    p.barrier();

    // k_den = x.T | z , y = Wout|x
    parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, c_line_size / 4, ref(p), Wout, x, y);
    p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});
    p.barrier();

    k_den += l;

    // k = z/k_den
    map1(0, Nr + 1, ref(p), Divide_by_const(), z, k_den, k);
    p.barrier();

    // Wold = .... ,  P = ...
    map2(0, Ny, 0, Nr + 1, -1, ref(p), Compute_new_Wout(), Wout, Wold, d, y, k);
    map2(0, Nr + 1, 0, Nr + 1, -1, ref(p), Compute_new_P(), P, Pold, k, z, l);
    p.barrier();

    s = 0;
    for (int i = 0; i < Ny; i++)
    {
      s += pow(d[i] - y[i], 2);
    }
    error_norms.push_back(sqrt(s));

    cnt++;
  }
  return error_norms;
}

#ifndef utimer_cpp
#define utimer_cpp
#include "utimer.cpp"
#endif
vector<double> compute_average_times(int Nr, int n_samples, int n_trials, int max_par_degree, int c_line_size, Matrix_wrapper dataset, Matrix_wrapper dataset_n)
{
  int Nu = 4;
  int Ny = 4;
  float l = 0.995;
  float nabla = 0.1;

  ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
  float **W = n.W.m;
  float **Win = n.Win.m;

  vector<double> times(max_par_degree + 1);
  times[0] = std::numeric_limits<double>::quiet_NaN();
  for (int par_deg = 1; par_deg < max_par_degree + 1; par_deg++)
  {
    {
      vector<double> ts(n_trials);
      for (int i = 0; i < n_trials; i++)
      {
        utimer t(to_string(par_deg), &ts[i]);
        par_train(par_deg, c_line_size, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win);
      }
      times[par_deg] = mean(ts);
    }
  }
  return times;
}