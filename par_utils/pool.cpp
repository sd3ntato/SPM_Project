#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#include <thread>
#include <queue>
#include "math.h"

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif


using namespace std;

void worker_fun(queue<Task *> &taskq, bool &stop, bool &no_tasks_todo)
{
  Task *task;
  while (!stop)
  {
    while (!taskq.empty())
    {
      task = taskq.front();
      taskq.pop();
      task->execute();
      delete task;
    }
    no_tasks_todo = 1;
  }
}

void join_threads(bool *&stops, int n, bool &global_stop, thread *&threads)
{
  while (!global_stop)
  {
  }
  for (int i = 0; i < n; i++)
  {
    stops[i] = 1;
  }
  for (int i = 0; i < n; i++)
  {
    threads[i].join();
  }
}

Pool::Pool(int n)
{
  n_workers = n;
  glob_stop = 0;
  taskqs = new queue<Task *>[n];
  no_tasks_todos = new bool[n];
  stops = new bool[n];
  threads = new thread[n];

  for (int i = 0; i < n; i++)
  {
    stops[i] = 0;
    no_tasks_todos[i] = 1;
    threads[i] = thread(worker_fun, ref(taskqs[i]), ref(stops[i]), ref(no_tasks_todos[i]));
  }

  thread join_deamon(join_threads, ref(stops), n, ref(glob_stop), ref(threads));
  join_deamon.detach();
}

void Pool::submit(vector<Task *> taskv)
{
  for (int i = 0; i < taskv.size(); i++)
  {
    taskqs[i % n_workers].push(taskv[i]);
  }
  for (int i = 0; i < n_workers; i++)
  {
    no_tasks_todos[i] = 0;
  }
}

void Pool::terminate()
{
  glob_stop = 1;
}

void Pool::await_no_tasks_todo()
{
  for (int i = 0; i < n_workers; i++)
  {
    while (!no_tasks_todos[i])
    {
    }
  }
  return;
}

int min(int a, int b)
{
  if (a < b)
  {
    return a;
  }
  return b;
}

float sum(float *a, int n)
{
  float s = 0.0;
  for (int i = 0; i < n; i++)
  {
    s += a[i];
  }
  return s;
}

float parallel_dot(float *a, float *b, int n, Pool &p)
{
  int k = 10; //size of cache
  float *results = new float[(int)ceil(n / k)];
  vector<Task *> taskv;
  for (int i = 0; i < n; i += k)
  {
    //cout<< i<< " "<< min(i + k, n)<< endl;
    taskv.push_back(new Dot_task(i, min(i + k, n), a, b, &results[int(i / k)]));
  }
  p.submit(taskv);
  p.await_no_tasks_todo();
  return sum(results, (int)ceil(n / k));
}
