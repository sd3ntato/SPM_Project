#include <iostream>
#include "linear_algebra.h"
#include <memory>
#include <queue>
#include <thread>
#include "math.h"
#include "utimer.cpp"

#ifndef pool_h
#define pool_h
#include <pool.h>
#endif



using namespace std;


int main()
{
  int n_samples = 1000;
  Pool p(5);

  auto a1 = generate_random_sparse_matrix(n_samples, 1, 0.5);
  auto b1 = generate_random_sparse_matrix(n_samples, 1, 0.5);

  float *a = a1.transpose().m[0];
  float *b = b1.transpose().m[0];

  {
    utimer t("par");
    float c = parallel_dot(a, b, n_samples, ref(p));
    cout << c << endl;
  }

  {
    utimer t("seq");
    auto c = dot(0, n_samples, a, b);
    cout << c << endl;
  }

  {
    utimer t("seq1");
    auto c = a1.transpose() | b1;
    print_matrix(c);
  }

  p.terminate();

  /*
  float *results = zeros(1, n_samples).m[0];
  vector<Task *> tasks;

  for (int i = 0; i < n_samples; i++)
  {
    float *x = ones(1, 10).m[0];
    float *y = ones(1, 10).m[0];
    float *res_ptr = &results[i];

    Task *t = new Dot_task(0, 10, x, y, res_ptr); // will be deleted by worker thread after execution
    tasks.push_back(t);
  }

  Pool p(10);
  p.submit(tasks);
  p.await_no_tasks_todo();
  p.terminate();

  cout << endl;

  for (int i = 0; i < n_samples; i++)
  {
    cout << results[i] << " " << i << endl;
  }
  */

  return 0;
}