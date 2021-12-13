#include <iostream>
#include <vector>
#include <ff/ff.hpp>
#include <ff/farm.hpp>
#include <mutex>
#include <condition_variable>
#include "prot_queue.h"
#include <chrono>
#include <thread>

float dot(int start, int stop, float *v1, float *v2)
{
  float s = 0.0;
  for (int i = start; i < stop; i++)
  {
    s += v1[i] * v2[i];
  }
  return s;
}

struct Task
{

public:
  bool terminated;
  std::mutex *mutex; // this is not cpyable, set it as a pointer and instantiate it in constructor, see difference struct and struct
  std::condition_variable *condition;
  Task()
  {
    mutex = new std::mutex;
    condition = new std::condition_variable;
    terminated = false;
  }
  virtual void execute(){};
};

struct Dot_task : public Task
{
private:
  int start;
  int stop;
  int i;
  float **M;
  float *y;
  float *res_ptr;

public:
  Dot_task() = default;
  Dot_task(int start, int stop, int i, float **M, float *y, float *r)
  {
    this->start = start;
    this->stop = stop;
    this->i = i;
    this->M = M;
    this->y = y;
    this->res_ptr = r;
  }
  void execute()
  {
    res_ptr[i] = dot(start, stop, M[i], y);
  }
};

float **zeros(int n1, int n2)
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

float *zeros(int n1)
{
  float *m = new float[n1];
  for (int j = 0; j < n1; j++)
  {
    m[j] = 0;
  }
  return m;
}

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
  } // it does nothing, just sends out tasks
};

struct Emitter : ff_node_t<Task>
{
  const long size = 10;
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

  ff_pool(prot_queue<Task *> *q, int nworkers)
  {
    taskq = q;
    farm = build_farm(q, nworkers);
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

  void terminate(){
    taskq->push(nullptr);
  }
};

int main(int argc, char *argv[])
{
  int nworkers = 5;
  prot_queue<Task *> *q = new prot_queue<Task *>;

  ff_pool p(q, nworkers);

  float **M = zeros(10, 10);
  float *v = zeros(10);
  float *r1 = zeros(10);
  float *r2 = zeros(10);
  float *r3 = zeros(10);

  std::vector<Task*> tasks = { new Dot_task(0, 10, 1, M, v, r1),new Dot_task(0, 10, 2, M, v, r2) ,new Dot_task(0, 10, 3, M, v, r3)};
  p.submit( tasks );

  return 0;
}

// g++ prova1.cpp  -o prova1 -I ./par_utils/ -I ./linear_algebra/ -I ./eigen/ -I ./spectra/include/ -I ./fastflow/ -pthread
