/**
 * @file   typedefs_thread.h
 * @author Gernot Walzl
 * @date   2011-11-21
 */

#ifndef TYPEDEFS_THREAD_H
#define TYPEDEFS_THREAD_H

#include "smarter_ptr.h"
#include <chrono>
#include <functional>
#include <mutex>
#include <shared_mutex>
#include <thread>

typedef SHARED_PTR<std::thread> ThreadSPtr;

typedef std::recursive_mutex RecursiveMutex;
typedef std::scoped_lock<std::recursive_mutex> Lock;

typedef std::shared_mutex SharedMutex;
typedef std::shared_lock<std::shared_mutex> ReadLock;
typedef std::unique_lock<std::shared_mutex> WriteLock;

inline void thread_sleep(unsigned int sleep_ms) {
    std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
}

#endif /* TYPEDEFS_THREAD_H */
