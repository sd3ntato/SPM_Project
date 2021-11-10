#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif

#include <thread>
#include <queue>

using namespace std;

class Pool
{
private:
  int n_workers;
  queue<Task *> *taskqs; // array of queues of tasks
  bool *no_tasks_todos;
  bool *stops;
  bool glob_stop;
  thread *threads; // array of workers

public:
  Pool(int n);
  void submit(vector<Task *> taskv);

  void terminate();

  void await_no_tasks_todo();
};

int min(int a, int b);
float sum(float *a, int n);
float parallel_dot(float *a, float *b, int n, Pool &p);
void worker_fun(queue<Task *> &taskq, bool &stop, bool &no_tasks_todo);
void join_threads(bool *&stops, int n, bool &global_stop, thread *&threads);
