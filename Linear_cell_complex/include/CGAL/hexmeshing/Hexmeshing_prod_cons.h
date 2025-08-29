// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//
#ifndef HEXMESHING_PROD_CONS_H
#define HEXMESHING_PROD_CONS_H

#include <queue>
#include <mutex>
#include <condition_variable>

namespace CGAL::internal::Hexmeshing {
  /**
   * @brief Thread-safe producer-consumer queue implementation
   * @tparam T The type of elements stored in the queue
   * 
   * This class implements a thread-safe queue using the producer-consumer pattern.
   * It provides synchronized access to a queue where multiple threads can safely
   * produce (push) and consume (wait and pop) items.
   */
  template <typename T>
  class ProdCons {
  public:
    /**
     * @brief Checks if the queue has any items
     * @return true if the queue is not empty, false otherwise
     */
    bool hasNext() const {
      return !queue.empty();
    }

    /**
     * @brief Waits for and retrieves the next item from the queue
     * @return The next item in the queue
     * 
     * If the queue is empty, this function blocks until an item becomes available.
     * When an item is available, it is removed from the queue and returned.
     */
    T waitNextItem(){
      std::unique_lock lock(m);

      while (!hasNext()){
        awake_signal.wait(lock);
      }

      auto item = std::move(queue.front());
      queue.pop();

      return item;
    }

    /**
     * @brief Pushes a new item into the queue
     * @param item The item to push (rvalue reference)
     * 
     * This overload pushes the item into the queue and notifies one waiting consumer.
     */
    void push(T&& item){
      std::unique_lock lock(m);
      queue.push(item);
      awake_signal.notify_one();
    }

    /**
     * @brief Pushes a new item into the queue
     * @param item The item to push (const reference)
     * 
     * This overload copies the item into the queue and notifies one waiting consumer.
     */
    void push(const T& item){
      std::unique_lock lock(m);
      queue.push(item);
      awake_signal.notify_one();
    }

  private:
    std::queue<T> queue;           ///< The underlying queue container
    std::mutex m;                  ///< Mutex for thread synchronization
    std::condition_variable awake_signal;  ///< Condition variable for thread notification
  };
}



#endif