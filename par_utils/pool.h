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

void worker_fun(prot_queue<Task *> &taskq)
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
  thread *threads;            // array of workers
  int last_submitted;         // indicates last worker thread to wich a task has been submitted

public:
  int n_workers;
  Pool(int n)
  {
    n_workers = n;
    taskqs = new prot_queue<Task *>[n]; // one queue for each worker thread
    threads = new thread[n];
    last_submitted = 0;

    // instantiate worker threads
    for (int i = 0; i < n; i++)
      threads[i] = thread(worker_fun, ref(taskqs[i]));
  }

  // given a vector of tasks, pushes them all into the task queue and then waits for their termination
  void submit(vector<Task *> taskv)
  {
    for (int i = 0; i < taskv.size(); i++)
    {
      last_submitted = (last_submitted + 1) % n_workers; // update thread to wich I submitted the last task
      taskqs[last_submitted].push(taskv[i]);             // push the task into the queue of the desired thread
    }
  }

  // push EOS to each worker
  void terminate()
  {
    for (int i = 0; i < n_workers; i++)
    {
      taskqs[i].push(nullptr);
    }
  }

  // waits for all the taskqs to be empty (meaning that all submitted tasks up to here are terminated)
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

  // divides the indices from begin to end into a proper number of slices, then creates tasks of the given type (Task_t),
  // each working on one of these slices, and submits them to the pool.
  // blocking version in ff_pool.
  // args is arguments to the contructor of task object
  template <typename T, typename... Args>
  void map(int begin, int end, T task, Args... args)
  {
    int n_points = end - begin;

    // diff is the size of the slices in wich the interval begin - end is divided.
    // a slice is either long (end-begin)/nworkers or the size of a cache line divided by the size of a float
    // ( last one to avoid false sharing )
    int diff = max((int)floor(n_points / n_workers), 128 / 4);
    if (diff == 0)
      diff = n_points;

    // now that we know the size of the intervals
    int start, stop;
    for (int i = begin; i < end; i += diff) // note i+=diff
    {
      start = i;
      stop = min(i + diff, end);
      // the current interval goes from start to stop, obv stop cannot be > end (size of the entire slice of intervals to be iterated through)

      // create a task of the given type, that works on indeces of the given slice and immediatelly submit it to the pool
      this->submit({new T(start, stop, args...)});
    }
  }

  // explots the thread pool and its map function to compute dot product in parallel. Deposits the result into apposite vector
  void parallel_matrix_dot_vector(int row_start, int row_stop, int col_start, int col_stop, int diff, float **M, float *v, float *r)
  {
    this->map(row_start, row_stop, Multiple_Dot_task(), col_start, col_stop, M, v, r);
    this->barrier();
  }
};