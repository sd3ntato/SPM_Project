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

#ifndef ff_hpp
#define ff_hpp
#include <ff/ff.hpp>
#endif

#ifndef ff_parallel_for_h
#define ff_parallel_for_f
#include <ff/parallel_for.hpp>
#endif

#ifndef ff_mdf_h
#define ff_mdf_h
#include <ff/mdf.hpp>
#endif

#ifndef ff_pool_h
#define ff_pool_h
#include "ff_Pool.hpp"
#endif

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#include "math.h"

#define train_iteration_pool()                                                  \
  /* compute rec and in parts of state update */                                \
  p.parallel_matrix_dot_vector(0, Nr, 0, Nr, c_line_size / 4, W, x_old, x_rec); \
  p.parallel_matrix_dot_vector(0, Nr, 0, Nu, c_line_size / 4, Win, u, x_in);    \
  p.barrier();                                                                  \
                                                                                \
  /* compute tanh(sum...) */                                                    \
  p.map1(0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);             \
  p.barrier();                                                                  \
                                                                                \
  /* z = P|x */                                                                 \
  p.parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, c_line_size / 4, P, x, z); \
  p.barrier();                                                                  \
                                                                                \
  /* k_den = x.T | z , y = Wout|x */                                            \
  p.parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, c_line_size / 4, Wout, x, y);  \
  p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});                        \
  p.barrier();                                                                  \
                                                                                \
  k_den += l;                                                                   \
                                                                                \
  /* k = z/k_den */                                                             \
  p.map1(0, Nr + 1, Divide_by_const(), z, k_den, k);                            \
  p.barrier();                                                                  \
                                                                                \
  /* Wold = .... ,  P = ... */                                                  \
  p.map2(0, Ny, 0, Nr + 1, -1, Compute_new_Wout(), Wout, Wold, d, y, k);        \
  p.map2(0, Nr + 1, 0, Nr + 1, -1, Compute_new_P(), P, Pold, k, z, l);          \
  p.barrier();

#define train_iteration_parfor()                                    \
  /* compute rec and in parts of state update */                    \
  /* x_rec = W x_old, x_in = Win u */                               \
  parallel_matrix_dot_vector_ff(0, Nr, 0, Nr, &p, W, x_old, x_rec); \
  parallel_matrix_dot_vector_ff(0, Nr, 0, Nu, &p, Win, u, x_in);    \
                                                                    \
  /* compute tanh(sum...) */                                        \
  comp_state_ff(Nr, Nu, x, x_rec, x_in, Win, x_old, &p);            \
                                                                    \
  /* z = P|x */                                                     \
  parallel_matrix_dot_vector_ff(0, Nr + 1, 0, Nr + 1, &p, P, x, z); \
                                                                    \
  /* k_den = x.T | z , y = Wout|x */                                \
  parallel_matrix_dot_vector_ff(0, Ny, 0, Nr + 1, &p, Wout, x, y);  \
  k_den = dot(0, Nr + 1, x, z);                                     \
  k_den += l;                                                       \
                                                                    \
  /* k = z/k_den */                                                 \
  divide_by_const_ff(k, z, k_den, Nr, &p);                          \
                                                                    \
  /* Wold = .... ,  P = ... */                                      \
  compute_new_Wout_ff(Wout, d, y, k, Wold, Nr, Ny, &p);             \
  compute_new_P_ff(P, Pold, k, z, l, Nr, &p);

#define mdf_train_iteration()                                                             \
  /* mdf obj */                                                                           \
  ff::ff_mdf mdf(taskGen, &Par, 8192, 1);                                                 \
                                                                                          \
  /* deposit input into parameters structure */                                           \
  Par.mdf = &mdf;                                                                         \
  Par.u = dataset_n.m[cnt];                                                               \
  Par.d = dataset.m[cnt + 1];                                                             \
                                                                                          \
  /* make the mdf run, it will use the given input to modify its data during execution */ \
  mdf.run_and_wait_end(); /* computes an iteration */                                     \
  /* mdf executes and deposits values into the float*s, so */                             \
  /* after execution, e.g. y contains the output of the */                                \
  /* network at the end of the iteration */                                               \
                                                                                          \
  s = compute_error(Par.d, Par.y, Ny);                                                    \
  error_norms.push_back(sqrt(s));                                                         \
                                                                                          \
  cnt++;

#define ff_pool_train_iteration()                                               \
  /* compute rec and in parts of state update */                                \
  p.parallel_matrix_dot_vector(0, Nr, 0, Nr, c_line_size / 4, W, x_old, x_rec); \
  p.parallel_matrix_dot_vector(0, Nr, 0, Nu, c_line_size / 4, Win, u, x_in);    \
                                                                                \
  /* compute tanh(sum...)*/                                                     \
  p.map1(0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);             \
                                                                                \
  /* z = P|x*/                                                                  \
  p.parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, c_line_size / 4, P, x, z); \
                                                                                \
  /* k_den = x.T | z , y = Wout|x*/                                             \
  p.parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, c_line_size / 4, Wout, x, y);  \
  p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});                        \
                                                                                \
  k_den += l;                                                                   \
                                                                                \
  /* k = z/k_den*/                                                              \
  p.map1(0, Nr + 1, Divide_by_const(), z, k_den, k);                            \
                                                                                \
  /* Wold = .... ,  P = ...*/                                                   \
  p.map2(0, Ny, 0, Nr + 1, -1, Compute_new_Wout(), Wout, Wold, d, y, k);        \
  p.map2(0, Nr + 1, 0, Nr + 1, -1, Compute_new_P(), P, Pold, k, z, l);

template <typename T>
struct Parameters
{
  int Nr, Nu, Ny, par_degree;
  float l;
  float **W, **Win, **Wout, **Wold, **P, **Pold;
  float *x, *x_rec, *x_in, *x_old, *u, *d, *k, *z, *y;
  T *mdf;
  ff_pool* p;
};

// istanzia DAG di una iterazione, sottopone le task. il risultato dell' operazione lo deposita in appsito puntatore
void taskGen(Parameters<ff::ff_mdf> *const Par)
{
  ff::ff_mdf *mdf = Par->mdf;
  ff_pool* p = Par->p;

  int Nr, Nu, Ny, par_degree;
  float **W, **Win, **Wout, **Wold, **P, **Pold;
  float *x, *x_rec, *x_in, *x_old, *u, *d, *k, *z, *y, l, k_den;
  unpack_values(Par);

  std::vector<ff::param_info> Param;
  const int c0 = 0;

  // x_rec = W | x_old
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Nr, c0, Nr, W, x_old, x_rec);
  //cout << "generated W xold" << endl << flush;

  // x_in = Win | u
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Nr, c0, Nu, Win, u, x_in);

  // compute tanh(sum...)
  mdf_submit_compute_state_state_task(mdf, Param, p, Nr, Nu, x, x_rec, x_in, Win, x_old);

  // z = P|x
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Nr, c0, Nr, P, x, z);

  // k_den = x.T | z
  mdf_submit_vector_dot_vector_task(mdf, Param, c0, Nr + 1, x, z, k_den);

  //y = Wout|x
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Ny, c0, Nr + 1, Wout, x, y);

  // k = z/k_den
  mdf_submit_divide_by_const(mdf, Param, p, k, z, k_den, Nr);

  // Wout = ...
  mdf_submit_compute_new_wout(Wout, d, y, k, Wold, Nr, Ny, p);

  // P = ...
  mdf_submit_compute_new_p(P, Pold, k, z, l, Nr, p);
}

vector<double> par_train(string ff, int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                         int Nr, int Nu, int Ny, float nabla, float l,
                         float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                         float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  init_train_loop();

  if (ff == "none")
  {
    std::cout<<"doing with none"<<std::endl;
    Pool p(par_degree);
    while (cnt < n_samples)
    {
      init_train_iteration();
      train_iteration_pool();
      end_train_iteration();
    }
  }
  else if (ff == "parfor")
  {
    std::cout<<"doing with parfor"<<std::endl;
    ff::ParallelFor p(par_degree);
    while (cnt < n_samples)
    {
      init_train_iteration();
      train_iteration_parfor();
      end_train_iteration();
    }
  }
  else if (ff == "ff_pool")
  {
    std::cout<<"doing with ff_pool"<<std::endl;
    ff_pool p(par_degree);
    while (cnt < n_samples)
    {
      init_train_iteration();
      ff_pool_train_iteration();
      end_train_iteration();
    }
  }
  else if (ff == "mdf")
  {
    std::cout<<"doing with mdf"<<std::endl;
    // data structure that contains data used by the DAG
    Parameters<ff::ff_mdf> Par;
    
    // thread pool that will be orchestrated by the mdf
    ff_pool p(par_degree);
    Par.p = &p;

    // fill up other parameters in Par strurcture
    pack_values(Par);
    while (cnt < n_samples)
    {
      mdf_train_iteration();
    }
  }

  return error_norms;
}