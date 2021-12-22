#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#include <thread>

#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef pool_utils_h
#define pool_utils_h
#include "pool_utils.hpp"
#endif

using namespace std;

void worker_fun(prot_queue<Task *> &taskq, bool &stop)
{
  Task *task;
  while (true)
  {
    task = taskq.pop(); // this blocks the thread until something inserts into the queue
    if (task)
    {
      task->execute();
      delete task;
    }
    else
      return; // nullptr = EOS
  }
}

class Pool
{
private:
  prot_queue<Task *> *taskqs; // array of queues of tasks
  bool *stops;
  thread *threads; // array of workers
  int last_submitted;
  mutex terminate_mutex;
  condition_variable terminate_condition;

public:
  int n_workers;
  Pool(int n)
  {
    n_workers = n;
    taskqs = new prot_queue<Task *>[n];
    stops = new bool[n];
    threads = new thread[n];
    last_submitted = 0;

    for (int i = 0; i < n; i++)
    {
      stops[i] = 0;
      threads[i] = thread(worker_fun, ref(taskqs[i]), ref(stops[i]));
    }
  }

  void submit(vector<Task *> taskv)
  {
    for (int i = 0; i < taskv.size(); i++)
    {
      last_submitted = (last_submitted + 1) % n_workers;
      taskqs[last_submitted].push(taskv[i]);
    }
  }

  void terminate()
  {
    for (int i = 0; i < n_workers; i++)
    {
      taskqs[i].push(nullptr);
    }
  }

  void barrier()
  {
    for (int i = 0; i < n_workers; i++)
    {
      taskqs[i].wait_empty();
      //assert( taskqs[i].empty() );
    }
    return;
  }

  ~Pool() { this->terminate(); }

  template <typename T, typename... Args>
  void map(int begin, int end, T task, Args... args)
  {
    int n_points = end - begin;
    int diff = max((int)floor(n_points / n_workers), 128 / 4);
    if (diff == 0)
      diff = n_points;
    int start, stop;
    for (int i = begin; i < end; i += diff)
    {
      start = i;
      stop = min(i + diff, end);
      this->submit({new T(start, stop, args...)});
    }
  }

  void parallel_matrix_dot_vector(int row_start, int row_stop, int col_start, int col_stop, int diff, float **M, float *v, float *r)
  {
    this->map(row_start, row_stop, Multiple_Dot_task(), col_start, col_stop, M, v, r);
    this->barrier();
  }
};