// Copyright (C) 2001-2003
// William E. Kempf
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  William E. Kempf makes no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.

#ifndef BOOST_THREAD_WEK070601_HPP
#define BOOST_THREAD_WEK070601_HPP

#include <boost/thread/detail/config.hpp>

#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/thread/mutex.hpp>
#include <list>
#include <memory>

#if defined(BOOST_HAS_PTHREADS)
#   include <pthread.h>
#   include <boost/thread/condition.hpp>
#elif defined(BOOST_HAS_MPTASKS)
#   include <Multiprocessing.h>
#endif

namespace boost {

struct xtime;

class BOOST_THREAD_DECL thread : private noncopyable
{
public:
    thread();
    explicit thread(const function0<void>& threadfunc);
    ~thread();

    bool operator==(const thread& other) const;
    bool operator!=(const thread& other) const;

    void join();

    static void sleep(const xtime& xt);
    static void yield();

private:
#if defined(BOOST_HAS_WINTHREADS)
    void* m_thread;
    unsigned int m_id;
#elif defined(BOOST_HAS_PTHREADS)
private:
    pthread_t m_thread;
#elif defined(BOOST_HAS_MPTASKS)
    MPQueueID m_pJoinQueueID;
    MPTaskID m_pTaskID;
#endif
    bool m_joinable;
};

class BOOST_THREAD_DECL thread_group : private noncopyable
{
public:
    thread_group();
    ~thread_group();

    thread* create_thread(const function0<void>& threadfunc);
    void add_thread(thread* thrd);
    void remove_thread(thread* thrd);
    void join_all();
        int size();

private:
    std::list<thread*> m_threads;
    mutex m_mutex;
};

} // namespace boost

// Change Log:
//    8 Feb 01  WEKEMPF Initial version.
//    1 Jun 01  WEKEMPF Added boost::thread initial implementation.
//    3 Jul 01  WEKEMPF Redesigned boost::thread to be noncopyable.

#endif // BOOST_THREAD_WEK070601_HPP
