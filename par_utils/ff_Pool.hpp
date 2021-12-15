
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

#include <ff/ff.hpp>
#include <ff/farm.hpp>
#include <vector>

using namespace ff;

#define set_terminated(t)                        \
  {                                              \
    {                                            \
      unique_lock<std::mutex> lock(*(t->mutex)); \
      t->terminated = true;                      \
    }                                            \
    t->condition->notify_all();                  \
  }

#define await_termination_of(task)                       \
  {                                                      \
    unique_lock<std::mutex> lock(*(task->mutex));        \
    task->condition->wait(lock, [=]                      \
                          { return task->terminated; }); \
  }

struct Worker : ff_node_t<Task>
{
  Task *svc(Task *t)
  {
    t->execute();
    set_terminated(t);
    return t;
  }
};

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
      t = taskq->pop();
      if (t)
        ff_send_out(t);
      else
        break;
    }
    return EOS;
  }
};

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

struct ff_pool
{
  ff_Farm<Task> *farm;
  prot_queue<Task *> *taskq;
  int nworkers;

  ff_pool(int nworkers)
  {
    this->nworkers = nworkers;
    taskq = new prot_queue<Task *>;
    farm = build_farm(taskq, nworkers);
  }

  ~ff_pool()
  {
    this->terminate();
    farm->wait_freezing();
  }

  void submit(vector<Task *> taskv)
  {
    for (int i = 0; i < taskv.size(); i++)
    {
      taskq->push(taskv[i]);
    }
    for (int i = 0; i < taskv.size(); i++)
    {
      await_termination_of(taskv[i]);
    } //esco da qui che tutte le task che ho submittato sono terminate
  }

  void terminate()
  {
    taskq->push(nullptr);
  }

  template <typename T, typename... Args>
  static void map1(ff_pool *p, int begin, int end, T task, Args... args)
  {
    std::vector<Task *> tasks;
    int n_points = end - begin;
    int diff = max((int)floor(n_points / p->nworkers), 128 / 4);
    if (diff == 0)
      diff = n_points;
    int start, stop;
    for (int i = begin; i < end; i += diff)
    {
      start = i;
      stop = min(i + diff, end);
      tasks.push_back(new T(start, stop, args...));
    }
    p->submit(tasks);
  }

  template <typename T, typename... Args>
  static void map2(ff_pool *p, int begin, int end, int start, int stop, int diff, T task, Args... args)
  {
    int i0, ii;
    std::vector<Task *> tasks;
    if (diff == -1)
    {
      int n_points = end - begin;
      if (n_points <= p->nworkers)
        diff = n_points;
      else
        diff = max((int)floor(n_points / p->nworkers), 128 / 4);
    }
    for (int i = begin; i < end; i += diff)
    {
      i0 = i;
      ii = min(i + diff, end);
      //                   horizontally   vertically
      //                   |          |  |     |
      tasks.push_back(new T(start, stop, i0, ii, args...));
    }
    p->submit(tasks);
  }

  static void parallel_matrix_dot_vector(ff_pool *p, int row_start, int row_stop, int col_start, int col_stop, int diff, float **M, float *v, float *r)
  {
    ff_pool::map2(p, row_start, row_stop, col_start, col_stop, diff, Multiple_Dot_task(), M, v, r);
  }

  static void comp_state(int Nr, int Nu, float *x, float *x_rec, float *x_in, float **Win, float *x_old, ff_pool *p)
  {
    p->map1(p, 0, Nr, Comp_state_task(), Nu, x_rec, x_in, Win, x, x_old);
  }

  static void div_by_const(float *k, float *z, float k_den, int stop, ff_pool *p)
  {
    p->map1(p, 0, stop, Divide_by_const(), z, k_den, k);
  }

  static void compute_new_wout(float **Wout, float *d, float *y, float *k, float **Wold, int Nr, int Ny, ff_pool *p)
  {
    p->map2(p, 0, Ny, 0, Nr + 1, -1, Compute_new_Wout(), Wout, Wold, d, y, k);
  }

  static void compute_new_P(float **P, float **Pold, float *k, float *z, float l, int Nr, ff_pool *p)
  {
    p->map2(p, 0, Nr + 1, 0, Nr + 1, -1, Compute_new_P(), P, Pold, k, z, l);
  }
};