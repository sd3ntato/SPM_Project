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
  bool glob_stop;
  thread *threads; // array of workers
  int last_submitted;

public:
  Pool(int n);
  void submit(vector<Task *> taskv);
  void terminate();
  void await_no_tasks_todo();
};

int min(int a, int b);
float sum(float *a, int n);
float parallel_dot(float *a, float *b, int n, Pool &p, int k);
void worker_fun(queue<Task *> &taskq, bool &stop);
void join_threads(bool *&stops, int n, bool &global_stop, thread *&threads);
