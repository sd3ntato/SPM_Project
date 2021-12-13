
#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
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

#define await_termination_of(task)                                                       \
  {                                                                                      \
    unique_lock<std::mutex> lock(*(task->mutex));                                        \
    task->condition->wait(lock, [=]                                                      \
                          {                                                              \
                            std::cout << "submit check" << task->terminated << std::endl \
                                      << std::flush;                                     \
                            return task->terminated;                                     \
                          });                                                            \
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

  ff_pool(int nworkers)
  {
    taskq = new prot_queue<Task *>;
    ;
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
};