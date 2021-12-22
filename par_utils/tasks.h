#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef basic_calculations_h
#define basic_calculations_h
#include "basic_calculations.hpp"
#endif

#include "math.h"
#include <iostream>
#include <mutex>
#include <condition_variable>
using namespace std;

// This is the mother Task class. Its constructo instantiates tools for asyncronous wait on the completition of the task.
// New tasks can be designed by overloading the execute method and encapsulating her arguments into child class contructor.
// See below for example usage
struct Task
{
public:
  bool terminated; // indicates wheter the task has been executed
  std::mutex *mutex; // this is not cpyable, so I have to set it as a pointer and instantiate it in constructor.
  std::condition_variable *condition; // same as above mutex

  // this gets called whenever one constructs a child class instantiation
  Task()
  {
    mutex = new std::mutex;
    condition = new std::condition_variable;
    terminated = false;
  }

  // This method is supposed to incapsulate the business logic of the task.
  // It is supposed to work on members of the child class
  virtual void execute(){};
};

// To compute dot product between to vectors and deposit the result into apposite placeholder
struct Dot_task : public Task
{
private:
  int start, stop, i;
  float **M, *y, *res_ptr;

public:
  Dot_task() = default;
  Dot_task(int start, int stop, int i, float **M, float *y, float *r)
  {
    this->start = start;
    this->stop = stop;
    this->i = i;
    this->M = M;
    this->y = y;
    this->res_ptr = r;
  }
  void execute()
  {
    dot_in_place(start, stop, M[i], y, &res_ptr[i]);
  }
};

// To compute novel state of the network.
// State of the network is defined as a vector of scalar numbers.
// computes entries of the state vector whose indeces range from this->start to this->stop
struct Comp_state_task : public Task
{
private:
  int start, stop, Nu;
  float *x_rec, *x_in, *x, *x_old, **Win;

public:
  Comp_state_task() = default;
  Comp_state_task(int start, int stop, int Nu, float *x_rec, float *x_in, float **Win, float *x, float *x_old)
  {
    this->start = start;
    this->stop = stop;
    this->x_rec = x_rec;
    this->x_in = x_in;
    this->Nu = Nu;
    this->Win = Win;
    this->x = x;
    this->x_old = x_old;
  }
  void execute()
  {
    for (int i = start; i < stop; i++)
    {
      comp_state_i(x, x_rec, x_in, x_in, Win, Nu, i);
    }
  }
};

// To divide a vector of scalar numbers by a constant scalar number.
// Computes elements whose indeces range from this->start to this->stop
struct Divide_by_const : public Task
{
private:
  int start, stop;
  float *v, k, *r;

public:
  Divide_by_const() = default;
  Divide_by_const(int start, int stop, float *v, float k, float *r)
  {
    this->start = start;
    this->stop = stop;
    this->v = v;
    this->k = k;
    this->r = r;
  }
  void execute()
  {
    for (int i = start; i < stop; i++)
    {
      divide_by_const(r, v, k, i);
    }
  }
};

// To compute new entries of Wout and copy old values into apposit matrix.
// computes the entries one row at a time, wokring on row this->i0 to this->ii (i infinite)
struct Compute_new_Wout : public Task
{
private:
  int start, stop, i0, ii;
  float **Wout, **Wold, *d, *y, *k;

public:
  Compute_new_Wout() = default;
  Compute_new_Wout(int i0, int ii,int start, int stop, float **Wout, float **Wold, float *d, float *y, float *k)
  {
    this->i0 = i0;
    this->ii = ii;
    this->start = start;
    this->stop = stop;
    this->Wout = Wout;
    this->Wold = Wold;
    this->d = d;
    this->y = y;
    this->k = k;
  }
  void execute()
  {
    for (int i = i0; i < ii; i++)
    {
      compute_line_of_wout(start, stop, Wout, Wold, d, y, k, i)
    }
  }
};

// To compute new entries of P and copy old values into apposit matrix.
// computes the entries one row at a time, wokring on row this->i0 to this->ii (i infinite)
struct Compute_new_P : public Task
{
private:
  int start, stop, i0, ii;
  float **P, **Pold, *k, *z, l;

public:
  Compute_new_P() = default;
  Compute_new_P( int i0, int ii, int start, int stop, float **P, float **Pold, float *k, float *z, float l)
  {
    this->i0 = i0;
    this->ii = ii;
    this->start = start;
    this->stop = stop;
    this->P = P;
    this->Pold = Pold;
    this->k = k;
    this->z = z;
    this->l = l;
  }
  void execute()
  {
    for (int i = i0; i < ii; i++)
    {
      compute_line_of_P(start, stop, P, Pold, k, z, l, i);
    }
  }
};

// To compute part of a matrix-to-vector dot product.
// computes the entries one row at a time, wokring on row this->i0 to this->ii (i infinite)
struct Multiple_Dot_task : public Task
{
private:
  int start, stop, i0, ii;
  float **M, *x, *y;

public:
  Multiple_Dot_task() = default;
  Multiple_Dot_task(int i0, int ii, int start, int stop, float **M, float *x, float *y)
  {
    this->i0 = i0;
    this->ii = ii;
    this->start = start;
    this->stop = stop;
    this->M = M;
    this->x = x;
    this->y = y;
  }
  void execute()
  {
    for (int i = i0; i < ii; i++)
    {
      dot_in_place(start, stop, M[i], x, &y[i]);
    }
  }
};
