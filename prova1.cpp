#include <iostream>
#include <vector>

#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef ff_pool_h
#define ff_pool_h
#include "ff_Pool.hpp"
#endif

#ifndef pool_h
#define pool_h
#include "pool.h"
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

float dot(int start, int stop, float *v1, float *v2)
{
  float s = 0.0;
  for (int i = start; i < stop; i++)
  {
    s += v1[i] * v2[i];
  }
  return s;
}

float **zeros1(int n1, int n2)
{
  float **m = new float *[n1];
  for (int i = 0; i < n1; i++)
  {
    m[i] = new float[n2];
    for (int j = 0; j < n2; j++)
    {
      m[i][j] = 0;
    }
  }
  return m;
}

float *zeros1(int n1)
{
  float *m = new float[n1];
  for (int j = 0; j < n1; j++)
  {
    m[j] = 0;
  }
  return m;
}

int main(int argc, char *argv[])
{
  int nworkers = 5;

  ff_pool ffp(nworkers);

  float **M = zeros1(1000, 1000);
  float *v = zeros1(1000);
  float *r = zeros1(1000);

  ffp.parallel_matrix_dot_vector(0, 1000, 0, 1000, -1, M, v, r);
  
  Pool p(nworkers);
  p.parallel_matrix_dot_vector(0, 1000, 0, 1000, -1, M, v, r);

  return 0;
}

// g++ prova1.cpp  -o prova1 -I ./par_utils/ -I ./linear_algebra/ -I ./eigen/ -I ./spectra/include/ -I ./fastflow/ -pthread
