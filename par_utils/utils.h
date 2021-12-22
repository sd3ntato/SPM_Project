#ifndef vector_h
#define vector_h
#include <vector>
#endif

#ifndef string_h
#define string_h
#include <string>
#endif

#ifndef ff_hpp
#define ff_hpp
#include <ff/ff.hpp>
#endif

#ifndef ff_parallel_for_h
#define ff_parallel_for_h
#include <ff/parallel_for.hpp>
#endif

#ifndef ff_pool_h
#define ff_pool_h
#include "ff_Pool.hpp"
#endif

#ifndef basic_calculations_h
#define basic_calculations_h
#include "basic_calculations.hpp"
#endif

namespace plt = matplotlibcpp;
using namespace std;

/************* PURE COMMON TEXT *********************/

/* prepares data structures needed to define the network and its functioning */
#define declare_vectors_and_matrices()    \
  int Nr = Nrs[i];                        \
  ESN n = ESN(Nr, 4, 4);                  \
  float **W = n.W.m;                      \
  float **Win = n.Win.m;                  \
  float **Wold = zeros(4, Nr + 1).m;      \
  float **Wout = zeros(4, Nr + 1).m;      \
  float **P = zeros(Nr + 1, Nr + 1).m;    \
  float **Pold = zeros(Nr + 1, Nr + 1).m; \
                                          \
  float *x = zeros(1, Nr + 1).m[0];       \
  float *x_rec = zeros(1, Nr + 1).m[0];   \
  float *x_in = zeros(1, Nr + 1).m[0];    \
  float *x_old = zeros(1, Nr + 1).m[0];   \
                                          \
  float *k = zeros(1, Nr + 1).m[0];       \
  float *z = zeros(1, Nr + 1).m[0];       \
  float *y = zeros(1, 4).m[0];

// resets the data structures needed to define the network and its functioning
#define init_train_loop()                                                                    \
  tie(Wout, Wold, P, Pold) = prepare_matrices(Nr, Ny, Wout, Wold, P, Pold, nabla);           \
  tie(x, x_rec, x_in, x_old, k, z, y) = prepare_vectors(Nr, x, x_rec, x_in, x_old, k, z, y); \
  vector<double> error_norms{};                                                              \
  int cnt = 0;                                                                               \
  x[Nr] = 1.0;                                                                               \
  x_old[Nr] = 1.0;                                                                           \
  float *u;                                                                                  \
  float *d;                                                                                  \
  float k_den;                                                                               \
  float s;

// loads input to the network and desired output for this iteration
#define init_train_iteration() \
  u = dataset_n.m[cnt];        \
  d = dataset.m[cnt + 1];

// computes error committed by the network this iteration, stores it into apposite placeholder and then continues training
#define end_train_iteration()     \
  s = compute_error(d, y, Ny);    \
  error_norms.push_back(sqrt(s)); \
  cnt++;

// resets matrices to restart training
tuple<float **, float **, float **, float **> prepare_matrices(int Nr, int Ny, float **Wout, float **Wold, float **P, float **Pold, float nabla)
{
  Wout = zeros(Ny, Nr + 1, Wout);
  Wold = zeros(Ny, Nr + 1, Wold);
  P = resetP(Nr + 1, Nr + 1, P, nabla);
  Pold = resetP(Nr + 1, Nr + 1, P, nabla);
  return make_tuple(Wout, Wold, P, Pold);
}

// resets vector to restart training
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

// d: desired output,
// y: actual output
float compute_error(float *d, float *y, int Ny)
{
  float s = 0;
  for (int i = 0; i < Ny; i++)
  {
    s += pow(d[i] - y[i], 2);
  }
  return s;
}

/************* FUNCTIONS FOR FF_PARFOR COMPUTATION *********************/

void parallel_matrix_dot_vector_ff(int row_start, int row_stop, int col_start, int col_stop, ff::ParallelFor *p, float **M, float *v, float *r)
{
  auto dot_f = [&](const int i)
  {
    dot_in_place(col_start, col_stop, M[i], v, &r[i]);
  };
  p->parallel_for(row_start, row_stop, dot_f);
}

void comp_state_ff(int Nr, int Nu, float *x, float *x_rec, float *x_in, float **Win, float *x_old, ff::ParallelFor *p)
{
  auto comp_state_f = [&](const int i)
  {
    comp_state_i(x, x_rec, x_in, x_in, Win, Nu, i);
  };
  p->parallel_for(0, Nr, comp_state_f);
}

void divide_by_const_ff(float *k, float *z, float k_den, int Nr, ff::ParallelFor *p)
{
  auto divide_by_const_f = [&](const int i)
  {
    divide_by_const(k, z, k_den, i);
  };
  p->parallel_for(0, Nr + 1, divide_by_const_f);
}

void compute_new_Wout_ff(float **Wout, float *d, float *y, float *k, float **Wold, int Nr, int Ny, ff::ParallelFor *p)
{
  auto compute_new_wout_f = [&](const int i)
  {
    compute_line_of_wout(0, Nr + 1, Wout, Wold, d, y, k, i)
  };
  p->parallel_for(0, Ny, compute_new_wout_f);
}

void compute_new_P_ff(float **P, float **Pold, float *k, float *z, float l, int Nr, ff::ParallelFor *p)
{
  auto compute_new_p_f = [&](const int i)
  {
    compute_line_of_P(0, Nr + 1, P, Pold, k, z, l, i);
  };
  p->parallel_for(0, Nr + 1, compute_new_p_f);
}

/************* MACROS FOR MDF CALCULATIONS *********************/

#define mdf_submit_matrix_dot_vector(mdf, Param, p, row_start, row_stop, col_start, col_stop, M, v, r)                    \
  {                                                                                                                       \
    {                                                                                                                     \
      Param.clear();                                                                                                      \
      const ff::param_info _1 = {(uintptr_t)&row_start, ff::INPUT};                                                       \
      const ff::param_info _2 = {(uintptr_t)&row_stop, ff::INPUT};                                                        \
      const ff::param_info _3 = {(uintptr_t)&col_start, ff::INPUT};                                                       \
      const ff::param_info _4 = {(uintptr_t)&col_stop, ff::INPUT};                                                        \
      const ff::param_info _5 = {(uintptr_t)&p, ff::INPUT};                                                               \
      const ff::param_info _6 = {(uintptr_t)&M, ff::INPUT};                                                               \
      const ff::param_info _7 = {(uintptr_t)&v, ff::INPUT};                                                               \
      const ff::param_info _8 = {(uintptr_t)&r, ff::OUTPUT};                                                              \
      Param.push_back(_1);                                                                                                \
      Param.push_back(_2);                                                                                                \
      Param.push_back(_3);                                                                                                \
      Param.push_back(_4);                                                                                                \
      Param.push_back(_5);                                                                                                \
      Param.push_back(_6);                                                                                                \
      Param.push_back(_7);                                                                                                \
      Param.push_back(_8);                                                                                                \
      mdf->AddTask(Param, ff_pool::parallel_matrix_dot_vector, p, row_start, row_stop, col_start, col_stop, -1, M, v, r); \
    }                                                                                                                     \
  }

#define mdf_submit_compute_state_state_task(mdf, Param, p, Nr, Nu, x, x_rec, x_in, Win, x_old) \
  {                                                                                            \
    {                                                                                          \
      Param.clear();                                                                           \
      const ff::param_info _1 = {(uintptr_t)&Nr, ff::INPUT};                                   \
      const ff::param_info _2 = {(uintptr_t)&Nu, ff::INPUT};                                   \
      const ff::param_info _3 = {(uintptr_t)&x, ff::OUTPUT};                                   \
      const ff::param_info _4 = {(uintptr_t)&x_rec, ff::INPUT};                                \
      const ff::param_info _5 = {(uintptr_t)&x_in, ff::INPUT};                                 \
      const ff::param_info _6 = {(uintptr_t)&Win, ff::INPUT};                                  \
      const ff::param_info _7 = {(uintptr_t)&x_old, ff::OUTPUT};                               \
      const ff::param_info _8 = {(uintptr_t)&p, ff::INPUT};                                    \
      Param.push_back(_1);                                                                     \
      Param.push_back(_2);                                                                     \
      Param.push_back(_3);                                                                     \
      Param.push_back(_4);                                                                     \
      Param.push_back(_5);                                                                     \
      Param.push_back(_6);                                                                     \
      Param.push_back(_7);                                                                     \
      Param.push_back(_8);                                                                     \
      mdf->AddTask(Param, ff_pool::comp_state, Nr, Nu, x, x_rec, x_in, Win, x_old, p);         \
    }                                                                                          \
  }

#define mdf_submit_comp_k_den(mdf, Param, p, start, stop, x, z, k_den_ptr, l)       \
  {                                                                                 \
    {                                                                               \
      Param.clear();                                                                \
      const ff::param_info _1 = {(uintptr_t)&start, ff::INPUT};                     \
      const ff::param_info _2 = {(uintptr_t)&stop, ff::INPUT};                      \
      const ff::param_info _3 = {(uintptr_t)&x, ff::INPUT};                         \
      const ff::param_info _4 = {(uintptr_t)&z, ff::INPUT};                         \
      const ff::param_info _5 = {(uintptr_t)&k_den_ptr, ff::OUTPUT};                \
      const ff::param_info _6 = {(uintptr_t)&l, ff::INPUT};                         \
      const ff::param_info _7 = {(uintptr_t)&p, ff::INPUT};                         \
      Param.push_back(_1);                                                          \
      Param.push_back(_2);                                                          \
      Param.push_back(_3);                                                          \
      Param.push_back(_4);                                                          \
      Param.push_back(_5);                                                          \
      Param.push_back(_6);                                                          \
      Param.push_back(_7);                                                          \
      mdf->AddTask(Param, ff_pool::comp_k_den, start, stop, x, z, k_den_ptr, l, p); \
    }                                                                               \
  }

#define mdf_submit_divide_by_const(mdf, Param, p, k, z, k_den_ptr, stop)    \
  {                                                                         \
    {                                                                       \
      Param.clear();                                                        \
      const ff::param_info _1 = {(uintptr_t)&k, ff::OUTPUT};                \
      const ff::param_info _2 = {(uintptr_t)&z, ff::INPUT};                 \
      const ff::param_info _3 = {(uintptr_t)&k_den_ptr, ff::INPUT};         \
      const ff::param_info _4 = {(uintptr_t)&stop, ff::INPUT};              \
      const ff::param_info _5 = {(uintptr_t)&p, ff::INPUT};                 \
      Param.push_back(_1);                                                  \
      Param.push_back(_2);                                                  \
      Param.push_back(_3);                                                  \
      Param.push_back(_4);                                                  \
      Param.push_back(_5);                                                  \
      mdf->AddTask(Param, ff_pool::div_by_const, k, z, k_den_ptr, stop, p); \
    }                                                                       \
  }

#define mdf_submit_compute_new_wout(Wout, d, y, k, Wold, Nr, Ny, p)                   \
  {                                                                                   \
    {                                                                                 \
      Param.clear();                                                                  \
      const ff::param_info _1 = {(uintptr_t)&Wout, ff::OUTPUT};                       \
      const ff::param_info _2 = {(uintptr_t)&d, ff::INPUT};                           \
      const ff::param_info _3 = {(uintptr_t)&y, ff::INPUT};                           \
      const ff::param_info _4 = {(uintptr_t)&k, ff::INPUT};                           \
      const ff::param_info _5 = {(uintptr_t)&Wold, ff::OUTPUT};                       \
      const ff::param_info _6 = {(uintptr_t)&Nr, ff::INPUT};                          \
      const ff::param_info _7 = {(uintptr_t)&Ny, ff::INPUT};                          \
      const ff::param_info _8 = {(uintptr_t)&p, ff::INPUT};                           \
      Param.push_back(_1);                                                            \
      Param.push_back(_2);                                                            \
      Param.push_back(_3);                                                            \
      Param.push_back(_4);                                                            \
      Param.push_back(_5);                                                            \
      Param.push_back(_6);                                                            \
      Param.push_back(_7);                                                            \
      Param.push_back(_8);                                                            \
      mdf->AddTask(Param, ff_pool::compute_new_wout, Wout, d, y, k, Wold, Nr, Ny, p); \
    }                                                                                 \
  }

#define mdf_submit_compute_new_p(P, Pold, k, z, l, Nr, p)                   \
  {                                                                         \
    {                                                                       \
      Param.clear();                                                        \
      const ff::param_info _1 = {(uintptr_t)&P, ff::OUTPUT};                \
      const ff::param_info _2 = {(uintptr_t)&Pold, ff::INPUT};              \
      const ff::param_info _3 = {(uintptr_t)&k, ff::INPUT};                 \
      const ff::param_info _4 = {(uintptr_t)&z, ff::INPUT};                 \
      const ff::param_info _5 = {(uintptr_t)&l, ff::INPUT};                 \
      const ff::param_info _6 = {(uintptr_t)&Nr, ff::INPUT};                \
      const ff::param_info _7 = {(uintptr_t)&p, ff::INPUT};                 \
      Param.push_back(_1);                                                  \
      Param.push_back(_2);                                                  \
      Param.push_back(_3);                                                  \
      Param.push_back(_4);                                                  \
      Param.push_back(_5);                                                  \
      Param.push_back(_6);                                                  \
      Param.push_back(_7);                                                  \
      mdf->AddTask(Param, ff_pool::compute_new_P, P, Pold, k, z, l, Nr, p); \
    }                                                                       \
  }

// fills up Par data strucure
#define pack_values(Par)         \
  {                              \
    Par.Nr = Nr;                 \
    Par.Nu = Nu;                 \
    Par.Ny = Ny;                 \
    Par.par_degree = par_degree; \
    Par.l = l;                   \
    Par.W = W;                   \
    Par.Win = Win;               \
    Par.Wout = Wout;             \
    Par.Wold = Wold;             \
    Par.P = P;                   \
    Par.Pold = Pold;             \
    Par.x = x;                   \
    Par.x_rec = x_rec;           \
    Par.x_in = x_in;             \
    Par.x_old = x_old;           \
    Par.k = k;                   \
    Par.z = z;                   \
    Par.y = y;                   \
  }

// gets values stored in Par data strucure
#define unpack_values(Par)        \
  {                               \
    Nr = Par->Nr;                 \
    Nu = Par->Nu;                 \
    Ny = Par->Ny;                 \
    par_degree = Par->par_degree; \
    W = Par->W;                   \
    Win = Par->Win;               \
    Wout = Par->Wout;             \
    Wold = Par->Wold;             \
    P = Par->P;                   \
    Pold = Par->Pold;             \
    x = Par->x;                   \
    x_rec = Par->x_rec;           \
    x_in = Par->x_in;             \
    x_old = Par->x_old;           \
    u = Par->u;                   \
    d = Par->d;                   \
    k = Par->k;                   \
    z = Par->z;                   \
    y = Par->y;                   \
    l = Par->l;                   \
    k_den = 0.0;                  \
  }
