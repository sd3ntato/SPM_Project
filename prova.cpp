#include <iostream>
#include "linear_algebra.h"
#include <memory>
#include <queue>
#include <thread>
#include <chrono>
#include "math.h"

#define START(timename) auto timename = std::chrono::system_clock::now();
#define STOP(timename, elapsed) auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - timename).count();

class utimer
{
  std::chrono::system_clock::time_point start;
  std::chrono::system_clock::time_point stop;
  std::string message;
  using usecs = std::chrono::microseconds;
  using msecs = std::chrono::milliseconds;

private:
  long *us_elapsed;

public:
  utimer(const std::string m) : message(m), us_elapsed((long *)NULL)
  {
    start = std::chrono::system_clock::now();
  }

  utimer(const std::string m, long *us) : message(m), us_elapsed(us)
  {
    start = std::chrono::system_clock::now();
  }

  ~utimer()
  {
    stop = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    auto musec = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

    std::cout << message << " computed in " << musec << " usec " << std::endl;
    if (us_elapsed != NULL)
      (*us_elapsed) = musec;
  }
};

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
  while (!stop)
  {
    while (!taskq.empty())
    {
      task = taskq.front();
      taskq.pop();
      task->execute();
      delete task;
    }
    no_tasks_todo = 1;
  }
}

void join_threads(bool *&stops, int n, bool &global_stop, thread *&threads)
{
  while (!global_stop)
  {
  }
  for (int i = 0; i < n; i++)
  {
    stops[i] = 1;
  }
  for (int i = 0; i < n; i++)
  {
    threads[i].join();
  }
}

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
  Pool(int n)
  {
    n_workers = n;
    glob_stop = 0;
    taskqs = new queue<Task *>[n];
    no_tasks_todos = new bool[n];
    stops = new bool[n];
    threads = new thread[n];

    for (int i = 0; i < n; i++)
    {
      stops[i] = 0;
      no_tasks_todos[i] = 1;
      threads[i] = thread(worker_fun, ref(taskqs[i]), ref(stops[i]), ref(no_tasks_todos[i]));
    }

    thread join_deamon(join_threads, ref(stops), n, ref(glob_stop), ref(threads));
    join_deamon.detach();
  }

  void submit(vector<Task *> tasks)
  {
    for (int i = 0; i < tasks.size(); i++)
    {
      taskqs[i % n_workers].push(tasks[i]);
    }
    for (int i = 0; i < n_workers; i++)
    {
      no_tasks_todos[i] = 0;
    }
  }

  void terminate()
  {
    glob_stop = 1;
  }

  void await_no_tasks_todo()
  {
    for (int i = 0; i < n_workers; i++)
    {
      while (!no_tasks_todos[i])
      {
      }
    }
    return;
  }
};

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

float parallel_dot(float *a, float *b, int n, Pool &p)
{
  int k = 10; //size of cache
  float *results = new float[(int)ceil(n / k)];
  vector<Task *> tasks;
  for (int i = 0; i < n; i += k)
  {
    //cout<< i<< " "<< min(i + k, n)<< endl;
    tasks.push_back(new Dot_task(i, min(i + k, n), a, b, &results[int(i / k)]));
  }
  p.submit(tasks);
  p.await_no_tasks_todo();
  return sum(results, (int)ceil(n / k));
}

int main()
{
  int n_samples = 1000;
  Pool p(5);

  auto a1 = generate_random_sparse_matrix(n_samples, 1, 0.5);
  auto b1 = generate_random_sparse_matrix(n_samples, 1, 0.5);

  float *a = a1.transpose().m[0];
  float *b = b1.transpose().m[0];

  {
    utimer t("par");
    float c = parallel_dot(a, b, n_samples, ref(p));
    cout << c << endl;
  }

  {
    utimer t("seq");
    auto c = dot(0, n_samples, a, b);
    cout << c << endl;
  }

  {
    utimer t("seq1");
    auto c = a1.transpose() | b1;
    print_matrix(c);
  }

  p.terminate();

  /*
  float *results = zeros(1, n_samples).m[0];
  vector<Task *> tasks;

  for (int i = 0; i < n_samples; i++)
  {
    float *x = ones(1, 10).m[0];
    float *y = ones(1, 10).m[0];
    float *res_ptr = &results[i];

    Task *t = new Dot_task(0, 10, x, y, res_ptr); // will be deleted by worker thread after execution
    tasks.push_back(t);
  }

  Pool p(10);
  p.submit(tasks);
  p.await_no_tasks_todo();
  p.terminate();

  cout << endl;

  for (int i = 0; i < n_samples; i++)
  {
    cout << results[i] << " " << i << endl;
  }
  */

  return 0;
}