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

struct Worker : ff_node_t<Task>
{
  Task *svc(Task *t)
  {
    std::cout << t << std::endl;
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
    std::cout << "emitter built" << std::endl;
  }
  Task *svc(Task *)
  {
    Task *t;
    for (long i = 0; i <= size; ++i)
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

struct ff_pool
{
  prot_queue<Task *> *q;
  ff_Farm<Task>* farm;

  ff_pool(prot_queue<Task *> *q, int nworkers)
  {
    
  }
};

//void deamon_funct(prot_queue<Task *> *q, int nworkers)
//auto th = thread(deamon_funct,q,nworkers);

int
main(int argc, char *argv[])
{
  int nworkers = 5;
  prot_queue<Task *> *q = new prot_queue<Task *>;

  ff_Farm<Task> farm([nworkers]()
                     {
                       std::vector<std::unique_ptr<ff_node>> Workers;
                       for (int i = 0; i < nworkers; ++i)
                         Workers.push_back(std::unique_ptr<ff_node_t<Task>>(new Worker));
                       return Workers;
                     }());

  Emitter E(q);
  farm.add_emitter(E);
  farm.remove_collector();

  if (farm.run_then_freeze() < 0)
    error("running farm");

  q->push(new Task());
  q->push(new Task());
  q->push(new Task());
  q->push(nullptr);

  farm.wait_freezing();

  return 0;
}
