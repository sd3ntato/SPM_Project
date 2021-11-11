#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#include <thread>
#include "math.h"


#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif


using namespace std;

void worker_fun(queue<Task *> &taskq, bool &stop)
{
  Task *task;
  while (!stop)
  {
    while (!taskq.empty()) // impelemtnare wait quando la coda e vuota senno fa attesa attiva 
    {
      task = taskq.pop();
      task->execute();
      delete task;
    }
  }
}

void join_threads(bool *&stops, int n, bool &global_stop, thread *&threads)
{
  while (!global_stop) // questo non puo proprio fare attesa attiva
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
  stops = new bool[n];
  threads = new thread[n];
  last_submitted = 0;

  for (int i = 0; i < n; i++)
  {
    stops[i] = 0;
    threads[i] = thread(worker_fun, ref(taskqs[i]), ref(stops[i]));
  }

  thread join_deamon(join_threads, ref(stops), n, ref(glob_stop), ref(threads));
  join_deamon.detach();
}

void Pool::submit(vector<Task *> taskv)
{
  for (int i = 0; i < taskv.size(); i++)
  {
    last_submitted = ( last_submitted + 1 ) % n_workers;
    taskqs[last_submitted].push(taskv[i]);
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
    while (! taskqs[i].empty())
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

float parallel_dot(float *a, float *b, int n, Pool &p, int k)
{
  float *results = new float[(int)ceil(n / k)];
  vector<Task *> taskv;
  for (int i = 0; i < n; i += k)
  {
    //cout<< i<< " "<< min(i + k, n)<< endl;
    taskv.push_back(new Dot_task(i, min(i + k, n), a, b, &results[int(i / k)]));
  }
  p.submit(taskv);
  p.await_no_tasks_todo();
  float s = sum(results, (int)ceil(n / k));
  delete [] results;
  return s;
}
