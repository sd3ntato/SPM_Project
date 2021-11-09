#include <iostream>
#include "linear_algebra.h"
#include <memory>
#include <queue>
#include <thread>
#include <chrono>

using namespace std;

float dot(int start, int stop, float *v1, float *v2)
{
  float s = 0.0;
  for (int i = start; i < stop; i++)
  {
    s += v1[i] * v2[i];
  }
  return s;
}

class Task
{
public:
  virtual void execute(){};
};

class Dot_task : public Task
{
private:
  int start;
  int stop;
  float *v1;
  float *v2;
  float *res_ptr;

public:
  Dot_task(int s, int st, float *x, float *y, float *r)
  {
    start = s;
    stop = st;
    v1 = x;
    v2 = y;
    res_ptr = r;
  }
  void execute()
  {
    *res_ptr = dot(start, stop, v1, v2);
  }
};

void worker_fun(queue<Task *> &taskq, bool &stop, bool &no_tasks_todo)
{
  Task *task;
  int i = 0;
  while (!stop)
  {
    while (!taskq.empty())
    {
      task = taskq.front();
      cout << "task " << i << endl;
      taskq.pop();
      task->execute();
      delete task;
      i++;
    }
    no_tasks_todo = 1;
  }
}

void join_threads(bool* &stops, int n, bool &global_stop, thread* & threads){
  while (!global_stop) {}
  for(int i=0;i<n;i++){
    stops[i] = 1;
  }
  for(int i=0;i<n;i++){
    threads[i].join();
  }
}

class Pool
{
private:
  int n_workers;
  queue<Task *> *taskqs; // array of queues of tasks
  bool* no_tasks_todos;
  bool* stops;
  bool glob_stop;
  thread *threads;      // array of workers

public:
  Pool(int n)
  {
    n_workers = n;
    glob_stop = 0;
    taskqs = new queue<Task *> [n];
    no_tasks_todos = new bool[n];
    stops = new bool[n];
    threads = new thread[n];

    for(int i=0; i<n;i++){
      stops[i] = 0; no_tasks_todos[i]=1;
      threads[i] = thread(worker_fun, ref(taskqs[i]), ref(stops[i]), ref(no_tasks_todos[i]) );
    }

    thread join_deamon(join_threads, ref(stops), n, ref(glob_stop), ref(threads) );
    join_deamon.detach();
  }

  void submit(vector<Task *> tasks){
    for(int i=0; i< tasks.size(); i++){
      taskqs[i%n_workers].push(tasks[i]); 
    }
    for(int i=0;i<n_workers;i++){
      no_tasks_todos[i]=0;
    }
  }

  void terminate(){
    glob_stop = 1;
  }

};

int main()
{
  int n_samples = 100;

  float *results = zeros(1, n_samples).m[0];
  vector<Task *> tasks;

  for (int i = 0; i < n_samples; i++)
  {
    float *x = ones(1, 10).m[0];
    float *y = ones(1, 10).m[0];
    float *res_ptr = &results[i];

    Task *t = new Dot_task(0, 10, x, y, res_ptr);
    tasks.push_back(t);
  }

  Pool p(10);
  p.submit(tasks);
  std::this_thread::sleep_for(std::chrono::seconds(1));
  p.terminate();

  cout << endl;

  for (int i = 0; i < n_samples; i++)
  {
    cout << results[i] << " " << i << endl;
  }

  return 0;
}