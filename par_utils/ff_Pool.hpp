
#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#ifndef pool_utils_h
#define pool_utils_h
#include "pool_utils.hpp"
#endif

#ifndef ff_hpp
#define ff_hpp
#include <ff/ff.hpp>
#endif

#include <ff/farm.hpp>

#ifndef vector_h
#define vector_h
#include <vector>
#endif

using namespace ff;

// given a task it sets its terminated flag in thread safe manner
#define set_terminated(t)                        \
  {                                              \
    {                                            \
      unique_lock<std::mutex> lock(*(t->mutex)); \
      t->terminated = true;                      \
    }                                            \
    t->condition->notify_all();                  \
  }

// given a task it waits for its termination flag to be set
#define await_termination_of(task)                       \
  {                                                      \
    unique_lock<std::mutex> lock(*(task->mutex));        \
    task->condition->wait(lock, [=]                      \
                          { return task->terminated; }); \
  }

// the worker threads receveive tasks in their input channel
// and just execute them and set them as terminated
struct Worker : ff_node_t<Task>
{
  Task *svc(Task *t)
  {
    t->execute();
    set_terminated(t);
    return t;
  }
};

// the emitter node waits on his tasks-queue to be non-empty.
// When the queue is non empty it pushes the tasks to the worker threads
struct Emitter : ff_node_t<Task>
{
  prot_queue<Task *> *taskq;

  Emitter(prot_queue<Task *> *q)
  {
    taskq = q;
  }

  Task *svc(Task *)
  {
    Task *t;
    while (true)
    {
      t = taskq->pop(); // here emitter blocks until someone puts something in his queue
      if (t)            // not a nullptr: send out a task
        ff_send_out(t);
      else // when emitter pops a nullptr from his taskq its sends out EOS and terminates
        break;
    }
    return EOS;
  }
};

// function to build the farm that runs at the core of ff_pool object
ff_Farm<Task> *build_farm(prot_queue<Task *> *q, int nworkers)
{
  ff_Farm<Task> *farm = new ff_Farm<Task>([nworkers]()
                                          {
                                            std::vector<std::unique_ptr<ff_node>> Workers;
                                            for (int i = 0; i < nworkers; ++i)
                                              Workers.push_back(std::unique_ptr<ff_node_t<Task>>(new Worker));
                                            return Workers;
                                          }());

  Emitter *E = new Emitter(q);
  farm->add_emitter(*E);
  farm->remove_collector();

  if (farm->run_then_freeze() < 0)
    error("running farm");

  return farm;
}

/* ff_pool works by a combination of an ff_farm and a protected queue.
* The custom emitter in the ff_Farm waits for tasks to be pushed on the queue, 
* then just passes them to the worker threads
*/
struct ff_pool
{
  ff_Farm<Task> *farm;
  prot_queue<Task *> *taskq;
  int nworkers;

  ff_pool(int nworkers)
  {
    this->nworkers = nworkers;
    taskq = new prot_queue<Task *>;
    farm = build_farm(taskq, nworkers); // when farm is built the taskq is passed to the emitter thread
  }

  ~ff_pool()
  {
    this->terminate(); // makes the emitter push EOS to all workers
    farm->wait_freezing();
    delete farm;
    delete taskq;
  }

  // given a vector of tasks, pushes them all into the task queue and then waits for their termination
  void submit(vector<Task *> taskv)
  {
    for (int i = 0; i < taskv.size(); i++)
    {
      taskq->push(taskv[i]);
    }
    for (int i = 0; i < taskv.size(); i++)
    {
      await_termination_of(taskv[i]);
    } // when exiting the loop all tasks have been marked as terminated
  }

  // makes the emitter push EOS to all workers
  void terminate()
  {
    taskq->push(nullptr);
  }

  // divides the indices from begin to end into a proper number of slices, then creates tasks of the given type (Task_t),
  // each working on one of these slices, and submits them to the pool.
  // blocking version in ff_pool.
  // args is arguments to the contructor of task object
  template <typename Task_t, typename... Args>
  static void map1(ff_pool *p, int begin, int end, Task_t task, Args... args)
  {
    std::vector<Task *> tasks;
    int n_points = end - begin;

    // diff is the size of the slices in wich the interval begin - end is divided.
    // a slice is either long (end-begin)/nworkers or the size of a cache line divided by the size of a float 
    // ( last one to avoid false sharing )
    int diff = max((int)floor(n_points / p->nworkers), 128 / 4);

    if (diff == 0) 
      diff = n_points;

    // now that we know the size of the intervals
    int start, stop;
    for (int i = begin; i < end; i += diff) // note i+=diff
    {
      start = i;
      stop = min(i + diff, end);
      // the current interval goes from start to stop, obv stop cannot be > end (size of the entire slice of intervals to be iterated through)

      // create a task of the given type, that works on indeces of the given slice
      tasks.push_back(new Task_t(start, stop, args...));
    }

    //submit the tasks to the thread pool
    p->submit(tasks);
  }

  // identical to map1 except that it specifies two couples of indeces
  template <typename Task_t, typename... Args>
  static void map2(ff_pool *p, int begin, int end, int start, int stop, int diff, Task_t task, Args... args)
  {
    std::vector<Task *> tasks;
    
    // if diff is not given infer it from the number of points int the interval begin - end
    if (diff == -1)
    {
      int n_points = end - begin;
      if (n_points <= p->nworkers) // in case the interval is very small create a single task
        diff = n_points;
      else
        diff = max((int)floor(n_points / p->nworkers), 128 / 4);
    }

    int i0, ii;
    for (int i = begin; i < end; i += diff)
    {
      i0 = i;
      ii = min(i + diff, end);
      //                   horizontally   vertically
      //                   |          |  |     |
      tasks.push_back(new Task_t(start, stop, i0, ii, args...));
    }

    p->submit(tasks);
  }

  static void parallel_matrix_dot_vector(ff_pool *p, int row_start, int row_stop, int col_start, int col_stop, int diff, float **M, float *v, float *r)
  {
    ff_pool::map2(p, row_start, row_stop, col_start, col_stop, diff, Multiple_Dot_task(), M, v, r);
  }

  static void comp_state(int Nr, int Nu, float *x, float *x_rec, float *x_in, float **Win, float *x_old, ff_pool *p)
  {
    ff_pool::map1(p, 0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);
  }

  static void comp_k_den(int start, int stop, float *x, float *z, float *k_den, float l, ff_pool *p)
  {
    p->submit({new Dot_task(start, stop, 0, &x, z, k_den)});
    *k_den += l;
  }

  static void div_by_const(float *k, float *z, float *k_den, int stop, ff_pool *p)
  {
    ff_pool::map1(p, 0, stop, Divide_by_const(), z, *k_den, k);
  }

  static void compute_new_wout(float **Wout, float *d, float *y, float *k, float **Wold, int Nr, int Ny, ff_pool *p)
  {
    ff_pool::map2(p, 0, Ny, 0, Nr + 1, -1, Compute_new_Wout(), Wout, Wold, d, y, k);
  }

  static void compute_new_P(float **P, float **Pold, float *k, float *z, float l, int Nr, ff_pool *p)
  {
    ff_pool::map2(p, 0, Nr + 1, 0, Nr + 1, -1, Compute_new_P(), P, Pold, k, z, l);
  }
};