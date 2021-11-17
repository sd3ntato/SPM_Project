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
  float *v1;
  float *v2;
  float *res_ptr;

public:
  Dot_task(int s, int st, float *x, float *y, float *r)
  {
    start = s;
    stop = st;
    v1 = x;
    v2 = y;
    res_ptr = r;
  }
  void execute()
  {
    *res_ptr = dot(start, stop, v1, v2);
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
