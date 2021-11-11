#include <iostream>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <vector>
#include <cstddef>

using namespace std;

template <typename T>
class queue
{
private:
  std::mutex d_mutex;
  std::condition_variable d_condition;
  std::deque<T> d_queue;

public:
  queue() = default;
  void push(T const &value);
  T pop();
  bool empty();
};

template <typename T>
inline void queue<T>::push(T const &value)
{
  {
    unique_lock<std::mutex> lock(this->d_mutex);
    this->d_queue.push_front(value);
  }
  this->d_condition.notify_all();
}

template <typename T>
inline T queue<T>::pop()
{
  std::unique_lock<std::mutex> lock(this->d_mutex);
  this->d_condition.wait(lock, [=]
                         { return !this->d_queue.empty(); });
  T rc(std::move(this->d_queue.back()));
  this->d_queue.pop_back();
  return rc;
}

template <typename T>
inline bool queue<T>::empty()
{
  unique_lock<std::mutex> lock(this->d_mutex);
  return d_queue.empty();
}
