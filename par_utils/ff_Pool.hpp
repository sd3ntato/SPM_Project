
#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#include <ff/farm.hpp>
#include <ff/ff.hpp>
#include <vector>

struct Worker : ff::ff_node_t<Task>
{
  Task *svc(Task *task)
  {
    std::cout << "recieved task" << std::endl << std::flush;
    task->execute();
    {
      unique_lock<std::mutex> lock(*(task->mutex));
      task->terminated = true;
      std::cout << " task terminated " << std::endl << std::flush;
    }
    task->condition->notify_all();
    return task;
  }
};

struct Emitter : ff::ff_node_t<Task>
{
  prot_queue<Task *> *taskq;
  Emitter(prot_queue<Task *> *taskq)
  {
    this->taskq = taskq;
  }
  Task *svc(Task *t)
  {
    Task *task;
    while (true)
    {
      task = taskq->pop(); //this blocks until a task comes out
      std::cout << "popped task "<<task << std::endl << std::flush;
      if (task)
      {
        ff_send_out(task);
        std::cout << "sent task" << std::endl << std::flush;
      }
    }
    return t;
  }
};

class ff_Pool
{
private:
  ff::ff_farm *farm;
  prot_queue<Task *> taskq;

public:
  int nworkers;
  ff_Pool(int nworkers)
  {
    this->nworkers = nworkers;

    Emitter emitter(&taskq);

    std::vector<ff::ff_node *> Workers;
    for (int i = 0; i < nworkers; i++)
    {
      Workers.push_back(new Worker);
    }

    farm = new ff::ff_farm(Workers);
    farm->add_emitter(&emitter);
    farm->remove_collector();

    int e = farm->run_then_freeze();
    if(e<0){
      std::cout<<"error "<<e <<" running farm"<<std::endl<<std::flush;
    }
    else{
      std::cout<<"farm is running  "<<std::endl<<std::flush;
    }
  }

  void submit(vector<Task *> taskv)
  {
    for (int i = 0; i < taskv.size(); i++)
    {
      taskq.push(taskv[i]);
      std::cout << "pushed task "<<taskv[i] << std::endl << std::flush;
    }
    for (int i = 0; i < taskv.size(); i++)
    {
      unique_lock<std::mutex> lock(*(taskv[i]->mutex));
      taskv[i]->condition->wait(lock, [=]
                                { 
                                  std::cout << "submit check"<<taskv[i]->terminated << std::endl << std::flush;
                                  return taskv[i]->terminated; });
      std::cout << "submite task done" << std::endl << std::flush;
    } //esco da qui che tutte le task che ho submittato sono terminate
    std::cout << "submit done" << std::endl << std::flush;
  }
};