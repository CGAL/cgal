// Copyright (C) 2000 Stephen Cleary
//
// This file can be redistributed and/or modified under the terms found
//  in "copyright.html"
// This software and its documentation is provided "as is" without express or
//  implied warranty, and with no claim as to its suitability for any purpose.
//
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_POOL_MUTEX_HPP
#define BOOST_POOL_MUTEX_HPP

// Extremely Light-Weight wrapper classes for OS thread synchronization

// Configuration: for now, we just choose between pthread or Win32 mutexes or none

#define BOOST_MUTEX_HELPER_NONE         0
#define BOOST_MUTEX_HELPER_WIN32        1
#define BOOST_MUTEX_HELPER_PTHREAD      2

#ifdef BOOST_NO_MT
  // No multithreading -> make locks into no-ops
  #define BOOST_MUTEX_HELPER BOOST_MUTEX_HELPER_NONE
#else
  #ifdef __WIN32__
    #define BOOST_MUTEX_HELPER BOOST_MUTEX_HELPER_WIN32
  #else
    #include <unistd.h>
    #ifdef _POSIX_THREADS
      #define BOOST_MUTEX_HELPER BOOST_MUTEX_HELPER_PTHREAD
    #endif
  #endif
#endif

#ifndef BOOST_MUTEX_HELPER
  #error Unable to determine platform mutex type; define BOOST_NO_MT to assume single-threaded
#endif


#ifdef __WIN32__
  #include <windows.h>
#endif
#ifdef _POSIX_THREADS
  #include <pthread.h>
#endif

namespace boost {

namespace details {
namespace pool {

#ifdef __WIN32__

class win32_mutex
{
  private:
    CRITICAL_SECTION mtx;

    win32_mutex(const win32_mutex &);
    void operator=(const win32_mutex &);

  public:
    win32_mutex()
    { InitializeCriticalSection(&mtx); }

    ~win32_mutex()
    { DeleteCriticalSection(&mtx); }

    void lock()
    { EnterCriticalSection(&mtx); }

    void unlock()
    { LeaveCriticalSection(&mtx); }
};

#endif // defined(__WIN32__)

#ifdef _POSIX_THREADS

class pthread_mutex
{
  private:
    pthread_mutex_t mtx;

    pthread_mutex(const pthread_mutex &);
    void operator=(const pthread_mutex &);

  public:
    pthread_mutex()
    { pthread_mutex_init(&mtx, 0); }

    ~pthread_mutex()
    { pthread_mutex_destroy(&mtx); }

    void lock()
    { pthread_mutex_lock(&mtx); }

    void unlock()
    { pthread_mutex_unlock(&mtx); }
};

#endif // defined(_POSIX_THREADS)

class null_mutex
{
  private:
    null_mutex(const null_mutex &);
    void operator=(const null_mutex &);

  public:
    null_mutex() { }

    static void lock() { }
    static void unlock() { }
};

#if BOOST_MUTEX_HELPER == BOOST_MUTEX_HELPER_NONE
  typedef null_mutex default_mutex;
#elif BOOST_MUTEX_HELPER == BOOST_MUTEX_HELPER_WIN32
  typedef win32_mutex default_mutex;
#elif BOOST_MUTEX_HELPER == BOOST_MUTEX_HELPER_PTHREAD
  typedef pthread_mutex default_mutex;
#endif

} // namespace pool
} // namespace details

} // namespace boost

#undef BOOST_MUTEX_HELPER_WIN32
#undef BOOST_MUTEX_HELPER_PTHREAD
#undef BOOST_MUTEX_HELPER_NONE
#undef BOOST_MUTEX_HELPER

#endif
