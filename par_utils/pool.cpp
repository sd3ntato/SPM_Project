#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#include <thread>
#include "math.h"


#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef tasks_h
#define tasks_h
#include <tasks.h>
#endif


using namespace std;

void worker_fun(queue<Task *> &taskq, bool& stop, condition_variable& thread_condition, mutex& thread_mutex)
{
  Task *task;
  while (true)
  {
    //wait on this thread's confition variable until queue is not empty or join deamon 
    //signal to stop
    unique_lock<mutex> lock(thread_mutex);
    thread_condition.wait(lock);

    while (!taskq.empty()) // can assure queue is not empty and then go and pop bc this thread is the only one that pops
    {
      task = taskq.pop();
      task->execute();
      delete task;
    }

    // if I have been woken up because i have to turn off I turn off
    if(stop) return;
  }
}

void join_threads(bool *&stops, int n, thread *&threads, mutex& terminate_mutex, condition_variable& terminate_condition, condition_variable*& thread_conditions)
{
  unique_lock<mutex> lock(terminate_mutex);
  terminate_condition.wait(lock);

  // when this thread get woken up it tells all workers to turn off
  for (int i = 0; i < n; i++)
  {
    stops[i] = 1;
    thread_conditions[i].notify_one();
  }
  for (int i = 0; i < n; i++)
  {
    threads[i].join();
  }
}

Pool::Pool(int n)
{
  n_workers = n;
  taskqs = new queue<Task *>[n];
  stops = new bool[n];
  threads = new thread[n];
  last_submitted = 0;
  thread_conditions = new condition_variable[n];
  thread_mutexes = new mutex[n];

  for (int i = 0; i < n; i++)
  {
    stops[i] = 0;
    threads[i] = thread(worker_fun, ref(taskqs[i]), ref(stops[i]), ref(thread_conditions[i]), ref(thread_mutexes[i]));
  }

  thread join_deamon(join_threads, ref(stops), n, ref(threads), ref(terminate_mutex), ref(terminate_condition), ref(thread_conditions));
  join_deamon.detach();
}

void Pool::submit(vector<Task *> taskv)
{
  for (int i = 0; i < taskv.size(); i++)
  {
    last_submitted = ( last_submitted + 1 ) % n_workers;
    taskqs[last_submitted].push(taskv[i]);
    thread_conditions[last_submitted].notify_one();
  }
}

void Pool::terminate()
{
  terminate_condition.notify_one();
}

void Pool::await_no_tasks_todo()
{
  for (int i = 0; i < n_workers; i++)
  {
    while (! taskqs[i].empty())
    {
    }
  }
  return;
}

int min(int a, int b)
{
  if (a < b)
  {
    return a;
  }
  return b;
}

float sum(float *a, int n)
{
  float s = 0.0;
  for (int i = 0; i < n; i++)
  {
    s += a[i];
  }
  return s;
}

float parallel_dot(float *a, float *b, int n, Pool &p, int k)
{
  float *results = new float[(int)ceil(n / k)];
  for (int i = 0; i < n; i += k)
  {
    p.submit({ new Dot_task(i, min(i + k, n), a, b, &results[int(i / k)]) });
  }
  p.await_no_tasks_todo();
  float s = sum(results, (int)ceil(n / k));
  delete [] results;
  return s;
}
