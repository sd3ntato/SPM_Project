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

vector<double> par_train(int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                         int Nr, int Nu, int Ny, float nabla, float l,
                         float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                         float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  init_train_loop();

  Pool p(par_degree);

  while (cnt < n_samples)
  {

    init_train_iteration();

    // compute rec and in parts of state update
    p.parallel_matrix_dot_vector(0, Nr, 0, Nr, c_line_size / 4, W, x_old, x_rec);
    p.parallel_matrix_dot_vector(0, Nr, 0, Nu, c_line_size / 4, Win, u, x_in);
    p.barrier();

    // compute tanh(sum...)
    p.map1(0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);
    p.barrier();

    // z = P|x
    p.parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, c_line_size / 4, P, x, z);
    p.barrier();

    // k_den = x.T | z , y = Wout|x
    p.parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, c_line_size / 4, Wout, x, y);
    p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});
    p.barrier();

    k_den += l;

    // k = z/k_den
    p.map1(0, Nr + 1, Divide_by_const(), z, k_den, k);
    p.barrier();

    // Wold = .... ,  P = ...
    p.map2(0, Ny, 0, Nr + 1, -1, Compute_new_Wout(), Wout, Wold, d, y, k);
    p.map2(0, Nr + 1, 0, Nr + 1, -1, Compute_new_P(), P, Pold, k, z, l);
    p.barrier();

    end_train_iteration();
  }
  return error_norms;
}

vector<double> par_train_ff(int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                            int Nr, int Nu, int Ny, float nabla, float l,
                            float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                            float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  init_train_loop();

  ff::ParallelFor p(par_degree);

  while (cnt < n_samples)
  {

    init_train_iteration();

    // compute rec and in parts of state update
    // x_rec = W x_old, x_in = Win u
    parallel_matrix_dot_vector_ff(0, Nr, 0, Nr, &p, W, x_old, x_rec);
    parallel_matrix_dot_vector_ff(0, Nr, 0, Nu, &p, Win, u, x_in);

    // compute tanh(sum...)
    comp_state_ff(Nr, Nu, x, x_rec, x_in, Win, x_old, &p);

    // z = P|x
    parallel_matrix_dot_vector_ff(0, Nr + 1, 0, Nr + 1, &p, P, x, z);

    // k_den = x.T | z , y = Wout|x
    parallel_matrix_dot_vector_ff(0, Ny, 0, Nr + 1, &p, Wout, x, y);
    k_den = dot(0, Nr + 1, x, z);
    k_den += l;

    // k = z/k_den
    divide_by_const_ff(k, z, k_den, Nr, &p);

    // Wold = .... ,  P = ...
    compute_new_Wout_ff(Wout, d, y, k, Wold, Nr, Ny, &p);

    compute_new_P_ff(P, Pold, k, z, l, Nr, &p);

    end_train_iteration();
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

// istanzia DAG di una iterazione, sottopone le task. il risultato dell' operazione lo deposita in appsito puntatore
void taskGen(Parameters<ff::ff_mdf> *const Par)
{
  ff::ff_mdf *mdf = Par->mdf;

  int Nr, Nu, Ny, par_degree;
  float **W, **Win, **Wout, **Wold, **P, **Pold;
  float *x, *x_rec, *x_in, *x_old, *u, *d, *k, *z, *y, l, k_den;
  unpack_values(Par);

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

  init_train_loop();

  // data structure that contains data used by the DAG
  Parameters<ff::ff_mdf> Par;
  // fill up parameters int Par strurcture
  pack_values(Par);

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

vector<double> par_train_ffPool(int par_degree, int c_line_size, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                         int Nr, int Nu, int Ny, float nabla, float l,
                         float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                         float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  init_train_loop();

  ff_pool p(par_degree);

  while (cnt < n_samples)
  {

    init_train_iteration();

    // compute rec and in parts of state update
    p.parallel_matrix_dot_vector(0, Nr, 0, Nr, c_line_size / 4, W, x_old, x_rec);
    p.parallel_matrix_dot_vector(0, Nr, 0, Nu, c_line_size / 4, Win, u, x_in);

    // compute tanh(sum...)
    p.map1(0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);

    // z = P|x
    p.parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, c_line_size / 4, P, x, z);

    // k_den = x.T | z , y = Wout|x
    p.parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, c_line_size / 4, Wout, x, y);
    p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});

    k_den += l;

    // k = z/k_den
    p.map1(0, Nr + 1, Divide_by_const(), z, k_den, k);

    // Wold = .... ,  P = ...
    p.map2(0, Ny, 0, Nr + 1, -1, Compute_new_Wout(), Wout, Wold, d, y, k);
    p.map2(0, Nr + 1, 0, Nr + 1, -1, Compute_new_P(), P, Pold, k, z, l);

    end_train_iteration();
  }
  return error_norms;
}