#include <iostream>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <vector>
#include <chrono>
#include <cstddef>
#include <math.h>
#include <string>

#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

void queue<T>::push(T const &value)
{
  {
    std::unique_lock<std::mutex> lock(this->d_mutex);
    this->d_queue.push_front(value);
  }
  this->d_condition.notify_one();
}

T queue<T>::pop()
{
  std::unique_lock<std::mutex> lock(d_mutex);
  this->d_condition.wait(lock, [=]  { return !this->d_queue.empty(); });
  T rc(std::move(d_queue.back()));
  d_queue.pop_back();
  return rc;
}

bool queue<T>::empty()
{
  std::unique_lock<std::mutex> lock(this->d_mutex);
  this->d_condition.wait(lock);
  return this->d_queue.empty();
}
