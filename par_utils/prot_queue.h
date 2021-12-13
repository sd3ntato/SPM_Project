#include <iostream>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <vector>
#include <cstddef>

using namespace std;

template <typename T>
class prot_queue
{
private:
  std::mutex d_mutex;
  std::condition_variable d_condition;
  std::deque<T> d_queue;

public:
  prot_queue() = default;
  void push(T const &value);
  T pop();
  bool empty();
  void wait_empty();
};

template <typename T>
inline void prot_queue<T>::push(T const &value)
{
  {
    unique_lock<std::mutex> lock(this->d_mutex);
    this->d_queue.push_front(value);
  }
  this->d_condition.notify_all();
}

template <typename T>
inline T prot_queue<T>::pop()
{
  T rc;
  {
    std::unique_lock<std::mutex> lock(this->d_mutex);
    this->d_condition.wait(lock, [=]
                           { return !this->d_queue.empty(); });
    rc = this->d_queue.back();
    this->d_queue.pop_back();
  }
  d_condition.notify_one();
  return rc;
}

template <typename T>
inline bool prot_queue<T>::empty()
{
  unique_lock<std::mutex> lock(this->d_mutex);
  return d_queue.empty();
  d_condition.notify_all();
}

template <typename T>
inline void prot_queue<T>::wait_empty()
{
  unique_lock<std::mutex> lock(d_mutex);
  d_condition.wait(lock, [=]
                   { return d_queue.empty(); });
}
