#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#include "linear_algebra.h"

Dot_task::Dot_task(int s, int st, float *x, float *y, float *r)
{
  start = s;
  stop = st;
  v1 = x;
  v2 = y;
  res_ptr = r;
}
void Dot_task::execute()
{
  *res_ptr = dot(start, stop, v1, v2);
}

