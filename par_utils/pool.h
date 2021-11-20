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
  queue<Task *> *taskqs; // array of queues of tasks
  bool *stops;
  thread *threads; // array of workers
  int last_submitted;
  mutex terminate_mutex;
  condition_variable terminate_condition;

public:
  int n_workers;
  Pool(int n);
  ~Pool() { this->terminate(); }
  void submit(vector<Task *> taskv);
  void terminate();
  void barrier();
};

int min(int a, int b);
int max(int a, int b);
float sum(float *a, int n);
float parallel_dot(float *a, float *b, int n, Pool &p, int k);