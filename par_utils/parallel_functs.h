#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

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

void parallel_matrix_dot_vector(int row_start, int row_stop, int col_start, int col_stop, Pool &p, float **M, float *v, float *r)
{
  map2(row_start, row_stop, col_start, col_stop, p, Dot_task(), M, v, r);
}