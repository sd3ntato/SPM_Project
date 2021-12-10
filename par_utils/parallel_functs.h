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

#ifndef utils_h
#define utils_h
#include "utils.h"
#endif

#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <ff/mdf.hpp>
#include "ff_Pool.hpp"

#ifndef utimer_cpp
#define utimer_cpp
#include "utimer.cpp"
#endif

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

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
                         float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                         float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  tie(Wout, Wold, P, Pold) = prepare_matrices(Nr, Ny, Wout, Wold, P, Pold, nabla);
  tie(x, x_rec, x_in, x_old, k, z, y) = prepare_vectors(Nr, x, x_rec, x_in, x_old, k, z, y);

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

void parallel_matrix_dot_vector_ff(int row_start, int row_stop, int col_start, int col_stop, ff::ParallelFor *p, float **M, float *v, float *r)
{
  auto dot_f = [&](const int i)
  {
    r[i] = dot(col_start, col_stop, M[i], v);
  };
  p->parallel_for(row_start, row_stop, dot_f);
}

vector<double> par_train_ff(int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                            int Nr, int Nu, int Ny, float nabla, float l,
                            float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                            float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  tie(Wout, Wold, P, Pold) = prepare_matrices(Nr, Ny, Wout, Wold, P, Pold, nabla);
  tie(x, x_rec, x_in, x_old, k, z, y) = prepare_vectors(Nr, x, x_rec, x_in, x_old, k, z, y);

  vector<double> error_norms{};
  ff::ParallelFor p(par_degree);
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
    // x_rec = W x_old, x_in = Win u
    parallel_matrix_dot_vector_ff(0, Nr, 0, Nr, &p, W, x_old, x_rec);
    parallel_matrix_dot_vector_ff(0, Nr, 0, Nu, &p, Win, u, x_in);

    // compute tanh(sum...)
    auto comp_state_f = [&](const int i)
    {
      x[i] = tanh(x_rec[i] + x_in[i] + Win[i][Nu]); // Win[i]-> b[i] TODO
      x_old[i] = x[i];
    };
    p.parallel_for(0, Nr, comp_state_f);

    // z = P|x
    parallel_matrix_dot_vector_ff(0, Nr + 1, 0, Nr + 1, &p, P, x, z);

    // k_den = x.T | z , y = Wout|x
    parallel_matrix_dot_vector_ff(0, Ny, 0, Nr + 1, &p, Wout, x, y);
    k_den = dot(0, Nr + 1, x, z);
    k_den += l;

    // k = z/k_den
    auto divide_by_const_f = [&](const int i)
    {
      k[i] = (z[i] / k_den);
    };
    p.parallel_for(0, Nr + 1, divide_by_const_f);

    // Wold = .... ,  P = ...
    auto compute_new_wout_f = [&](const int i)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];
        Wold[i][j] = Wout[i][j];
      }
    };
    p.parallel_for(0, Ny, compute_new_wout_f);

    auto compute_new_p_f = [&](const int i)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;
        Pold[i][j] = P[i][j];
      }
    };
    p.parallel_for(0, Nr + 1, compute_new_p_f);

    s = compute_error(d, y, Ny);
    error_norms.push_back(sqrt(s));

    cnt++;
  }
  return error_norms;
}

template <typename T>
struct Parameters
{
  int Nr, Nu, Ny, par_degree;
  float l;
  float **W, **Win, **Wout, **Wold, **P, **Pold;
  float *x, *x_rec, *x_in, *x_old, *u, *d, *k, *z, *y;
  T *mdf;
};

void comp_state_ff(int Nr, int Nu, float *x, float *x_rec, float *x_in, float **Win, float *x_old, ff::ParallelFor *p)
{
  auto comp_state_f = [&](const int i)
  {
    x[i] = tanh(x_rec[i] + x_in[i] + Win[i][Nu]); // Win[i]-> b[i] TODO
    x_old[i] = x[i];
  };
  p->parallel_for(0, Nr, comp_state_f);
}

void divide_by_const_ff(float *k, float *z, float k_den, int Nr, ff::ParallelFor *p)
{
  auto divide_by_const_f = [&](const int i)
  {
    k[i] = (z[i] / k_den);
  };
  p->parallel_for(0, Nr + 1, divide_by_const_f);
}

void dot_in_place(int start, int stop, float *v1, float *v2, float *x)
{
  float s = 0.0;
  for (int i = start; i < stop; i++)
  {
    s += v1[i] * v2[i];
  }
  *x = s;
}

void compute_new_Wout_ff(float **Wout, float *d, float *y, float *k, float **Wold, int Nr, int Ny, ff::ParallelFor *p)
{
  auto compute_new_wout_f = [&](const int i)
  {
    for (int j = 0; j < Nr + 1; j++)
    {
      Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];
      Wold[i][j] = Wout[i][j];
    }
  };
  p->parallel_for(0, Ny, compute_new_wout_f);
}

void compute_new_P_ff(float **P, float **Pold, float *k, float *z, float l, int Nr, ff::ParallelFor *p)
{
  auto compute_new_p_f = [&](const int i)
  {
    for (int j = 0; j < Nr + 1; j++)
    {
      P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;
      Pold[i][j] = P[i][j];
    }
  };
  p->parallel_for(0, Nr + 1, compute_new_p_f);
}

#define mdf_submit_matrix_dot_vector(mdf, Param, ptr_p, row_start, row_stop, col_start, col_stop, M, v, r)          \
  {                                                                                                                 \
    {                                                                                                               \
      Param.clear();                                                                                                \
      const ff::param_info _1 = {(uintptr_t)&row_start, ff::INPUT};                                                 \
      const ff::param_info _2 = {(uintptr_t)&row_stop, ff::INPUT};                                                  \
      const ff::param_info _3 = {(uintptr_t)&col_start, ff::INPUT};                                                 \
      const ff::param_info _4 = {(uintptr_t)&col_stop, ff::INPUT};                                                  \
      const ff::param_info _5 = {(uintptr_t)&ptr_p, ff::INPUT};                                                     \
      const ff::param_info _6 = {(uintptr_t)&M, ff::INPUT};                                                         \
      const ff::param_info _7 = {(uintptr_t)&v, ff::INPUT};                                                         \
      const ff::param_info _8 = {(uintptr_t)&r, ff::OUTPUT};                                                        \
      Param.push_back(_1);                                                                                          \
      Param.push_back(_2);                                                                                          \
      Param.push_back(_3);                                                                                          \
      Param.push_back(_4);                                                                                          \
      Param.push_back(_5);                                                                                          \
      Param.push_back(_6);                                                                                          \
      Param.push_back(_7);                                                                                          \
      Param.push_back(_8);                                                                                          \
      mdf->AddTask(Param, parallel_matrix_dot_vector_ff, row_start, row_stop, col_start, col_stop, ptr_p, M, v, r); \
    }                                                                                                               \
  }

#define mdf_submit_compute_state_state_task(mdf, Param, ptr_p, Nr, Nu, x, x_rec, x_in, Win, x_old) \
  {                                                                                                \
    {                                                                                              \
      Param.clear();                                                                               \
      const ff::param_info _1 = {(uintptr_t)&Nr, ff::INPUT};                                       \
      const ff::param_info _2 = {(uintptr_t)&Nu, ff::INPUT};                                       \
      const ff::param_info _3 = {(uintptr_t)&x, ff::OUTPUT};                                       \
      const ff::param_info _4 = {(uintptr_t)&x_rec, ff::INPUT};                                    \
      const ff::param_info _5 = {(uintptr_t)&x_in, ff::INPUT};                                     \
      const ff::param_info _6 = {(uintptr_t)&Win, ff::INPUT};                                      \
      const ff::param_info _7 = {(uintptr_t)&x_old, ff::OUTPUT};                                   \
      const ff::param_info _8 = {(uintptr_t)&ptr_p, ff::INPUT};                                    \
      Param.push_back(_1);                                                                         \
      Param.push_back(_2);                                                                         \
      Param.push_back(_3);                                                                         \
      Param.push_back(_4);                                                                         \
      Param.push_back(_5);                                                                         \
      Param.push_back(_6);                                                                         \
      Param.push_back(_7);                                                                         \
      Param.push_back(_8);                                                                         \
      mdf->AddTask(Param, comp_state_ff, Nr, Nu, x, x_rec, x_in, Win, x_old, ptr_p);               \
    }                                                                                              \
  }

#define mdf_submit_vector_dot_vector_task(mdf, Param, start, stop, v1, v2, r) \
  {                                                                           \
    {                                                                         \
      Param.clear();                                                          \
      const ff::param_info _1 = {(uintptr_t)&start, ff::INPUT};               \
      const ff::param_info _2 = {(uintptr_t)&stop, ff::INPUT};                \
      const ff::param_info _3 = {(uintptr_t)&v1, ff::INPUT};                  \
      const ff::param_info _4 = {(uintptr_t)&v2, ff::INPUT};                  \
      const ff::param_info _5 = {(uintptr_t)&r, ff::OUTPUT};                  \
      Param.push_back(_1);                                                    \
      Param.push_back(_2);                                                    \
      Param.push_back(_3);                                                    \
      Param.push_back(_4);                                                    \
      Param.push_back(_5);                                                    \
      mdf->AddTask(Param, dot_in_place, start, stop, v1, v2, &r);             \
    }                                                                         \
  }

#define mdf_submit_divide_by_const(mdf, Param, ptr_p, k, z, k_den, stop) \
  {                                                                      \
    {                                                                    \
      Param.clear();                                                     \
      const ff::param_info _1 = {(uintptr_t)&k, ff::OUTPUT};             \
      const ff::param_info _2 = {(uintptr_t)&z, ff::INPUT};              \
      const ff::param_info _3 = {(uintptr_t)&k_den, ff::INPUT};          \
      const ff::param_info _4 = {(uintptr_t)&stop, ff::INPUT};           \
      const ff::param_info _5 = {(uintptr_t)&ptr_p, ff::INPUT};          \
      Param.push_back(_1);                                               \
      Param.push_back(_2);                                               \
      Param.push_back(_3);                                               \
      Param.push_back(_4);                                               \
      Param.push_back(_5);                                               \
      mdf->AddTask(Param, divide_by_const_ff, k, z, k_den, stop, ptr_p); \
    }                                                                    \
  }

#define mdf_subumit_compute_new_wout(Wout, d, y, k, Wold, Nr, Ny, ptr_p)            \
  {                                                                                 \
    {                                                                               \
      Param.clear();                                                                \
      const ff::param_info _1 = {(uintptr_t)&Wout, ff::OUTPUT};                     \
      const ff::param_info _2 = {(uintptr_t)&d, ff::INPUT};                         \
      const ff::param_info _3 = {(uintptr_t)&y, ff::INPUT};                         \
      const ff::param_info _4 = {(uintptr_t)&k, ff::INPUT};                         \
      const ff::param_info _5 = {(uintptr_t)&Wold, ff::OUTPUT};                     \
      const ff::param_info _6 = {(uintptr_t)&Nr, ff::INPUT};                        \
      const ff::param_info _7 = {(uintptr_t)&Ny, ff::INPUT};                        \
      const ff::param_info _8 = {(uintptr_t)&ptr_p, ff::INPUT};                     \
      Param.push_back(_1);                                                          \
      Param.push_back(_2);                                                          \
      Param.push_back(_3);                                                          \
      Param.push_back(_4);                                                          \
      Param.push_back(_5);                                                          \
      Param.push_back(_6);                                                          \
      Param.push_back(_7);                                                          \
      Param.push_back(_8);                                                          \
      mdf->AddTask(Param, compute_new_Wout_ff, Wout, d, y, k, Wold, Nr, Ny, ptr_p); \
    }                                                                               \
  }

#define mdf_submit_compute_new_p(P, Pold, k, z, l, Nr, ptr_p)             \
  {                                                                       \
    {                                                                     \
      Param.clear();                                                      \
      const ff::param_info _1 = {(uintptr_t)&P, ff::OUTPUT};              \
      const ff::param_info _2 = {(uintptr_t)&Pold, ff::INPUT};            \
      const ff::param_info _3 = {(uintptr_t)&k, ff::INPUT};               \
      const ff::param_info _4 = {(uintptr_t)&z, ff::INPUT};               \
      const ff::param_info _5 = {(uintptr_t)&l, ff::INPUT};               \
      const ff::param_info _6 = {(uintptr_t)&Nr, ff::INPUT};              \
      const ff::param_info _7 = {(uintptr_t)&ptr_p, ff::INPUT};           \
      Param.push_back(_1);                                                \
      Param.push_back(_2);                                                \
      Param.push_back(_3);                                                \
      Param.push_back(_4);                                                \
      Param.push_back(_5);                                                \
      Param.push_back(_6);                                                \
      Param.push_back(_7);                                                \
      mdf->AddTask(Param, compute_new_P_ff, P, Pold, k, z, l, Nr, ptr_p); \
    }                                                                     \
  }

// istanzia DAG di una iterazione, sottopone le task. il risultato dell' operazione lo deposita in appsito puntatore
void taskGen(Parameters<ff::ff_mdf> *const Par)
{
  ff::ff_mdf *mdf = Par->mdf;

  int Nr = Par->Nr;
  int Nu = Par->Nu;
  int Ny = Par->Ny;
  int par_degree = Par->par_degree;

  float **W = Par->W;
  float **Win = Par->Win;
  float **Wout = Par->Wout;
  float **Wold = Par->Wold;
  float **P = Par->P;
  float **Pold = Par->Pold;

  float *x = Par->x;
  float *x_rec = Par->x_rec;
  float *x_in = Par->x_in;
  float *x_old = Par->x_old;
  float *u = Par->u;
  float *d = Par->d;
  float *k = Par->k;
  float *z = Par->z;
  float *y = Par->y;

  float l = Par->l;
  float k_den = 0.0;

  std::vector<ff::param_info> Param;
  const int c0 = 0;

  // x_rec = W | x_old
  ff::ParallelFor *ptr_p1 = new ff::ParallelFor(par_degree);
  mdf_submit_matrix_dot_vector(mdf, Param, ptr_p1, c0, Nr, c0, Nr, W, x_old, x_rec);
  //cout << "generated W xold" << endl << flush;

  // x_in = Win | u
  ff::ParallelFor *ptr_p2 = new ff::ParallelFor(par_degree);
  mdf_submit_matrix_dot_vector(mdf, Param, ptr_p2, c0, Nr, c0, Nu, Win, u, x_in);

  // compute tanh(sum...)
  ff::ParallelFor *ptr_p3 = new ff::ParallelFor(par_degree);
  mdf_submit_compute_state_state_task(mdf, Param, ptr_p3, Nr, Nu, x, x_rec, x_in, Win, x_old);

  // z = P|x
  ff::ParallelFor *ptr_p4 = new ff::ParallelFor(par_degree);
  mdf_submit_matrix_dot_vector(mdf, Param, ptr_p4, c0, Nr, c0, Nr, P, x, z);

  // k_den = x.T | z
  mdf_submit_vector_dot_vector_task(mdf, Param, c0, Nr + 1, x, z, k_den);

  //y = Wout|x
  ff::ParallelFor *ptr_p5 = new ff::ParallelFor(par_degree);
  mdf_submit_matrix_dot_vector(mdf, Param, ptr_p5, c0, Ny, c0, Nr + 1, Wout, x, y);

  // k = z/k_den
  ff::ParallelFor *ptr_p6 = new ff::ParallelFor(par_degree);
  mdf_submit_divide_by_const(mdf, Param, ptr_p6, k, z, k_den, Nr);

  // Wout = ...
  ff::ParallelFor *ptr_p7 = new ff::ParallelFor(par_degree);
  mdf_subumit_compute_new_wout(Wout, d, y, k, Wold, Nr, Ny, ptr_p7);

  // P = ...
  ff::ParallelFor *ptr_p8 = new ff::ParallelFor(par_degree);
  mdf_submit_compute_new_p(P, Pold, k, z, l, Nr, ptr_p8);
}

vector<double> par_train_mdf(int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                             int Nr, int Nu, int Ny, float nabla, float l,
                             float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                             float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  tie(Wout, Wold, P, Pold) = prepare_matrices(Nr, Ny, Wout, Wold, P, Pold, nabla);
  tie(x, x_rec, x_in, x_old, k, z, y) = prepare_vectors(Nr, x, x_rec, x_in, x_old, k, z, y);

  vector<double> error_norms{};
  int cnt = 0;

  x[Nr] = 1.0;
  x_old[Nr] = 1.0;
  float *u;
  float *d;

  float s;

  // data structure that contains data used by the DAG
  Parameters<ff::ff_mdf> Par;
  // fill up parameters int Par strurcture
  Par.Nr = Nr;
  Par.Nu = Nu;
  Par.Ny = Ny;
  Par.par_degree = par_degree;
  Par.l = l;
  Par.W = W;
  Par.Win = Win;
  Par.Wout = Wout;
  Par.Wold = Wold;
  Par.P = P;
  Par.Pold = Pold;
  Par.x = x;
  Par.x_rec = x_rec;
  Par.x_in = x_in;
  Par.x_old = x_old;
  Par.k = k;
  Par.z = z;
  Par.y = y;

  while (cnt < n_samples)
  {
    // mdf obj
    ff::ff_mdf mdf(taskGen, &Par, 8192, 3);

    // deposit input into parameters structure
    Par.mdf = &mdf;
    Par.u = dataset_n.m[cnt];
    Par.d = dataset.m[cnt + 1];

    // make the mdf run, it will use the given input to modify its data during execution
    mdf.run_and_wait_end(); // computes an iteration
    // mdf executes and deposits values into the float*s, so
    // after execution, e.g. y contains the output of the
    // network at the end of the iteration

    s = compute_error(Par.d, Par.y, Ny);
    error_norms.push_back(sqrt(s));

    cnt++;
  }
  return error_norms;
}

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
        if (ff == "parfor")
          err = par_train_ff(par_deg, c_line_size, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);
        if (ff == "mdf")
          err = par_train_mdf(par_deg, c_line_size, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);
        if (ff == "none")
          err = par_train(par_deg, c_line_size, n_samples, dataset, dataset_n, Nr, Nu, Ny, nabla, l, W, Win, Wout, Wold, P, Pold, x, x_rec, x_in, x_old, k, z, y);
        plt::plot(err);
        plt::title("err with " + to_string(Nr) + " neurons and par deg " + to_string(par_deg));
        plt::save("imgs/" + to_string(Nr) + "-" + to_string(par_deg) + "-err");
        plt::show();
      }
      times[par_deg] = mean(ts);
    }
  }
  return times;
}