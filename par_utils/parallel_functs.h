#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
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

#include "math.h"

/* computes a training iteration with a Pool object. Pool p is an instantiation of my
* own version of a thread pool. Note that when it needs to synch results it calls p.barrier()
*/
#define train_iteration_pool()                                      \
  /* compute rec and in parts of state update */                    \
  p.parallel_matrix_dot_vector(0, Nr, 0, Nr, -1, W, x_old, x_rec);  \
  p.parallel_matrix_dot_vector(0, Nr, 0, Nu, -1, Win, u, x_in);     \
  p.barrier();                                                      \
                                                                    \
  /* compute tanh(sum...) */                                        \
  p.map(0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);  \
  p.barrier();                                                      \
                                                                    \
  /* z = P|x */                                                     \
  p.parallel_matrix_dot_vector(0, Nr + 1, 0, Nr + 1, -1, P, x, z);  \
  p.barrier();                                                      \
                                                                    \
  /* k_den = x.T | z , y = Wout|x */                                \
  p.parallel_matrix_dot_vector(0, Ny, 0, Nr + 1, -1, Wout, x, y);   \
  p.submit({new Dot_task(0, Nr + 1, 0, &x, z, &k_den)});            \
  p.barrier();                                                      \
                                                                    \
  k_den += l;                                                       \
                                                                    \
  /* k = z/k_den */                                                 \
  p.map(0, Nr + 1, Divide_by_const(), z, k_den, k);                 \
  p.barrier();                                                      \
                                                                    \
  /* Wold = .... ,  P = ... */                                      \
  p.map(0, Ny, Compute_new_Wout(), 0, Nr + 1, Wout, Wold, d, y, k); \
  p.map(0, Nr + 1, Compute_new_P(), 0, Nr + 1, P, Pold, k, z, l);   \
  p.barrier();

/* computes a training iteration with an ff_Pool object. ff_Pool p is an instantiation of my
* version of a thread pool realized through an ff_Farm object and a protected queue.
* Note that this time I am not using barriers to synch because the calls are blocking this time
*/
#define ff_pool_train_iteration()                                             \
  /* compute rec and in parts of state update */                              \
  ff_pool::parallel_matrix_dot_vector(&p, 0, Nr, 0, Nr, -1, W, x_old, x_rec); \
  ff_pool::parallel_matrix_dot_vector(&p, 0, Nr, 0, Nu, -1, Win, u, x_in);    \
                                                                              \
  /* compute tanh(sum...)*/                                                   \
  ff_pool::comp_state(Nr, Nu, x, x_rec, x_in, Win, x_old, &p);                \
                                                                              \
  /* z = P|x*/                                                                \
  ff_pool::parallel_matrix_dot_vector(&p, 0, Nr + 1, 0, Nr + 1, -1, P, x, z); \
                                                                              \
  /* k_den = x.T | z , y = Wout|x*/                                           \
  ff_pool::parallel_matrix_dot_vector(&p, 0, Ny, 0, Nr + 1, -1, Wout, x, y);  \
  ff_pool::comp_k_den(0, Nr + 1, x, z, &k_den, l, &p);                        \
                                                                              \
  /* k = z/k_den*/                                                            \
  ff_pool::div_by_const(k, z, &k_den, Nr + 1, &p);                            \
                                                                              \
  /* Wold = .... ,  P = ...*/                                                 \
  ff_pool::compute_new_wout(Wout, d, y, k, Wold, Nr, Ny, &p);                 \
  ff_pool::compute_new_P(P, Pold, k, z, l, Nr, &p);

// computes a training iterarion with ff::parallel_for object. (object is p)
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

// computes a train iteration by the use of an mdf object that orchestrates an ff_pool object
#define mdf_train_iteration()                                                             \
  /* mdf obj */                                                                           \
  ff::ff_mdf mdf(taskGen, &Par, 8192, 3);                                                 \
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

/* This data structure contains 
* - Data structures needed to define network and its functioning
* - mdf object to wich tasks will be submitted
* - ff_pool object tha will be used by the mdf object to accomplish the tasks
*/
template <typename T>
struct Parameters
{
  int Nr, Nu, Ny, par_degree;
  float l;
  float **W, **Win, **Wout, **Wold, **P, **Pold;
  float *x, *x_rec, *x_in, *x_old, *u, *d, *k, *z, *y;
  T *mdf;
  ff_pool *p;
};

// instantiates DAG for an iteration, submits the tasks. The result of the operation is put into apposit placeholder inside Par d.s.
void taskGen(Parameters<ff::ff_mdf> *const Par)
{
  ff::ff_mdf *mdf = Par->mdf;
  ff_pool *p = Par->p;

  int Nr, Nu, Ny, par_degree;
  float **W, **Win, **Wout, **Wold, **P, **Pold;
  float *x, *x_rec, *x_in, *x_old, *u, *d, *k, *z, *y, l, k_den;
  unpack_values(Par);
  int Nrp1 = Nr + 1;

  std::vector<ff::param_info> Param;
  const int c0 = 0;

  // x_rec = W | x_old
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Nr, c0, Nr, W, x_old, x_rec);

  // x_in = Win | u
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Nr, c0, Nu, Win, u, x_in);

  // compute tanh(sum...)
  mdf_submit_compute_state_state_task(mdf, Param, p, Nr, Nu, x, x_rec, x_in, Win, x_old);

  // z = P|x
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Nrp1, c0, Nrp1, P, x, z);

  // k_den = x.T | z
  float *k_den_ptr = &k_den;
  mdf_submit_comp_k_den(mdf, Param, p, c0, Nrp1, x, z, k_den_ptr, l);

  //y = Wout|x
  mdf_submit_matrix_dot_vector(mdf, Param, p, c0, Ny, c0, Nrp1, Wout, x, y);

  // k = z/k_den
  mdf_submit_divide_by_const(mdf, Param, p, k, z, k_den_ptr, Nrp1);

  // Wout = ...
  mdf_submit_compute_new_wout(Wout, d, y, k, Wold, Nr, Ny, p);

  // P = ...
  mdf_submit_compute_new_p(P, Pold, k, z, l, Nr, p);
}

/* the first param in par_traing is a string, it can be:
* - none for no ff - uses my implementation of thread pool -
* - parfor for ff_parallel_for
* - ff_pool for my impl. of thread pool using ff
* - mdf for mdf model orchestrating the ff_pool object
*/
vector<double> par_train(string ff, int par_degree, int n_samples, Matrix_wrapper dataset, Matrix_wrapper dataset_n,
                         int Nr, int Nu, int Ny, float nabla, float l,
                         float **W, float **Win, float **Wout, float **Wold, float **P, float **Pold,
                         float *x, float *x_rec, float *x_in, float *x_old, float *k, float *z, float *y)
{

  // resets the data structures needed to define the network and its functioning
  init_train_loop();

  // differentiate between the various modalities:
  if (ff == "none") // no ff, uses my implementation of thread pool
  {
    std::cout << "doing with none" << std::endl;

    Pool p(par_degree);

    while (cnt < n_samples)
    {
      init_train_iteration();
      train_iteration_pool();
      end_train_iteration();
    }
  }
  else if (ff == "parfor") // uses ff:parallelFor
  {
    std::cout << "doing with parfor" << std::endl;

    ff::ParallelFor p(par_degree);

    while (cnt < n_samples)
    {
      init_train_iteration();
      train_iteration_parfor();
      end_train_iteration();
    }
  }
  else if (ff == "ff_pool") // uses my implementation of thread pool that in turn uses ff::Farm
  {
    std::cout << "doing with ff_pool" << std::endl;

    ff_pool p(par_degree);

    while (cnt < n_samples)
    {
      init_train_iteration();
      ff_pool_train_iteration();
      end_train_iteration();
    }
  }
  else if (ff == "mdf") // uses ff::mdf to orchestrate an ff_pool object
  {
    std::cout << "doing with mdf" << std::endl;

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