#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif
#include "math.h"

class Task
{
public:
  virtual void execute(){};
};

class Dot_task : public Task
{
private:
  int start;
  int stop;
  int i;
  float **M;
  float *y;
  float *res_ptr;

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
    res_ptr[i] = dot(start, stop, M[i], y);
  }
};

class Comp_state_task : public Task
{
private:
  int start, stop, Nu;
  float *x_rec, *x_in, *x, **Win;

public:
  Comp_state_task() = default;
  Comp_state_task(int start, int stop, int Nu, float *x_rec, float *x_in, float **Win, float *x)
  {
    this->start = start;
    this->stop = stop;
    this->x_rec = x_rec;
    this->x_in = x_in;
    this->Nu = Nu;
    this->Win = Win;
    this->x = x;
  }
  void execute()
  {
    for (int i = start; i < stop; i++)
    {
      x[i] = tanh(x_rec[i] + x_in[i] + Win[i][Nu]); // Win[i]-> b[i] TODO
    }
  }
};

class Divide_by_const : public Task
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
      r[i] = (v[i] / k);
    }
  }
};

class Compute_new_Wout : public Task
{
private:
  int start, stop, i;
  float **Wout, **Wold, *d, *y, *k;

public:
  Compute_new_Wout() = default;
  Compute_new_Wout(int start, int stop, int i, float **Wout, float **Wold, float *d, float *y, float *k)
  {
    this->start = start;
    this->stop = stop;
    this->i = i;
    this->Wout = Wout;
    this->Wold = Wold;
    this->d = d;
    this->y = y;
    this->k = k;
  }
  void execute()
  {
    for (int j = start; j < stop; j++)
    {
      Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];
    }
  }
};

class Compute_new_P : public Task
{
private:
  int start, stop, i;
  float **P, **Pold, *k, *z, l;

public:
  Compute_new_P() = default;
  Compute_new_P(int start, int stop, int i, float **P, float **Pold, float *k, float *z, float l)
  {
    this->start = start;
    this->stop = stop;
    this->i = i;
    this->P = P;
    this->Pold = Pold;
    this->k = k;
    this->z = z;
    this->l = l;
  }
  void execute()
  {
    for (int j = start; j < stop; j++)
    {
      P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;
    }
  }
};
