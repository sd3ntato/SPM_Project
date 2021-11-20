#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#include "math.h"

template <typename T, typename... Args>
void map1(int begin, int end, int diff, Pool &p, T task, Args... args)
{
  int start, stop;
  for (int i = begin; i < end; i += diff)
  {
    start = i;
    stop = min(i + diff, end);
    p.submit({new T(start, stop, args...)});
  }
}

template <typename T, typename... Args>
void map2(int begin, int end, int start, int stop, Pool &p, T task, Args... args)
{
  for (int i = begin; i < end; i++)
  {
    p.submit({new T(start, stop, i, args...)});
  }
}

template <typename T, typename... Args>
void map3(int begin, int end, int start, int stop, int diff, Pool &p, T task, Args... args)
{
  int i0, ii;
  for (int i = begin; i < end; i += diff)
  {
    i0 = i;
    ii = min(i + diff, end);
    //cout << " submitting " << i0 << " " << ii << endl
    //     << flush;
    //            horizontally   vertically
    //            |           |  |     |
    p.submit({new T(start, stop, i0, ii, args...)});
  }
}

void parallel_matrix_dot_vector(int row_start, int row_stop, int col_start, int col_stop, Pool &p, float **M, float *v, float *r)
{
  int n_dots = row_stop - row_start;
  //int diff = max( (int)floor(n_dots / p.n_workers), 128/4 );
  //if (diff == 0)
  //  diff = n_dots;
  int diff = 128/4; // this does not work as I want for some reason
  map3(row_start, row_stop, col_start, col_stop, diff, p, Multiple_Dot_task(), M, v, r);
  p.barrier();
}