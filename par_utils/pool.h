#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#include <thread>

#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

using namespace std;

class Pool
{
private:
  int n_workers;
  queue<Task *> *taskqs; // array of queues of tasks
  bool *stops;
  thread *threads; // array of workers
  int last_submitted;
  mutex terminate_mutex;
  condition_variable terminate_condition;

public:
  Pool(int n);
  ~Pool() { this->terminate(); }
  void submit(vector<Task *> taskv);
  void terminate();
  void await_no_tasks_todo();
};

int min(int a, int b);
float sum(float *a, int n);
float parallel_dot(float *a, float *b, int n, Pool &p, int k);

template <typename T, typename... Args>
void mapp(int begin, int end, int diff, Pool &p, T task, Args... args)
{
  int start, stop;
  for (int i = begin; i < end; i += diff)
  {
    start = i;
    stop = min(i + diff, end);
    p.submit({new T(start, stop, args...)});
  }
  p.await_no_tasks_todo();
}