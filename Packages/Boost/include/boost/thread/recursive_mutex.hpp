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

#ifndef BOOST_RECURSIVE_MUTEX_WEK070601_HPP
#define BOOST_RECURSIVE_MUTEX_WEK070601_HPP

#include <boost/thread/detail/config.hpp>

#include <boost/utility.hpp>
#include <boost/thread/detail/lock.hpp>

#if defined(BOOST_HAS_PTHREADS)
#   include <pthread.h>
#endif

#if defined(BOOST_HAS_MPTASKS)
#   include "scoped_critical_region.hpp"
#endif

namespace boost {

struct xtime;

class BOOST_THREAD_DECL recursive_mutex
    : private noncopyable
{
public:
    friend class detail::thread::lock_ops<recursive_mutex>;

    typedef detail::thread::scoped_lock<recursive_mutex> scoped_lock;

    recursive_mutex();
    ~recursive_mutex();

private:
#if (defined(BOOST_HAS_WINTHREADS) || defined(BOOST_HAS_MPTASKS))
    typedef std::size_t cv_state;
#elif defined(BOOST_HAS_PTHREADS)
    struct cv_state
    {
        long count;
        pthread_mutex_t* pmutex;
    };
#endif
    void do_lock();
    void do_unlock();
    void do_lock(cv_state& state);
    void do_unlock(cv_state& state);

#if defined(BOOST_HAS_WINTHREADS)
    void* m_mutex;
    bool m_critical_section;
    unsigned long m_count;
#elif defined(BOOST_HAS_PTHREADS)
    pthread_mutex_t m_mutex;
    unsigned m_count;
#   if !defined(BOOST_HAS_PTHREAD_MUTEXATTR_SETTYPE)
    pthread_cond_t m_unlocked;
    pthread_t m_thread_id;
    bool m_valid_id;
#   endif
#elif defined(BOOST_HAS_MPTASKS)
    threads::mac::detail::scoped_critical_region m_mutex;
    threads::mac::detail::scoped_critical_region m_mutex_mutex;
    std::size_t m_count;
#endif
};

class BOOST_THREAD_DECL recursive_try_mutex
    : private noncopyable
{
public:
    friend class detail::thread::lock_ops<recursive_try_mutex>;

    typedef detail::thread::scoped_lock<recursive_try_mutex> scoped_lock;
    typedef detail::thread::scoped_try_lock<
        recursive_try_mutex> scoped_try_lock;

    recursive_try_mutex();
    ~recursive_try_mutex();

private:
#if (defined(BOOST_HAS_WINTHREADS) || defined(BOOST_HAS_MPTASKS))
    typedef std::size_t cv_state;
#elif defined(BOOST_HAS_PTHREADS)
    struct cv_state
    {
        long count;
        pthread_mutex_t* pmutex;
    };
#endif
    void do_lock();
    bool do_trylock();
    void do_unlock();
    void do_lock(cv_state& state);
    void do_unlock(cv_state& state);

#if defined(BOOST_HAS_WINTHREADS)
    void* m_mutex;
    bool m_critical_section;
    unsigned long m_count;
#elif defined(BOOST_HAS_PTHREADS)
    pthread_mutex_t m_mutex;
    unsigned m_count;
#   if !defined(BOOST_HAS_PTHREAD_MUTEXATTR_SETTYPE)
    pthread_cond_t m_unlocked;
    pthread_t m_thread_id;
    bool m_valid_id;
#   endif
#elif defined(BOOST_HAS_MPTASKS)
    threads::mac::detail::scoped_critical_region m_mutex;
    threads::mac::detail::scoped_critical_region m_mutex_mutex;
    std::size_t m_count;
#endif
};

class BOOST_THREAD_DECL recursive_timed_mutex
    : private noncopyable
{
public:
    friend class detail::thread::lock_ops<recursive_timed_mutex>;

    typedef detail::thread::scoped_lock<recursive_timed_mutex> scoped_lock;
    typedef detail::thread::scoped_try_lock<
        recursive_timed_mutex> scoped_try_lock;
    typedef detail::thread::scoped_timed_lock<
        recursive_timed_mutex> scoped_timed_lock;

    recursive_timed_mutex();
    ~recursive_timed_mutex();

private:
#if (defined(BOOST_HAS_WINTHREADS) || defined(BOOST_HAS_MPTASKS))
    typedef std::size_t cv_state;
#elif defined(BOOST_HAS_PTHREADS)
    struct cv_state
    {
        long count;
        pthread_mutex_t* pmutex;
    };
#endif
    void do_lock();
    bool do_trylock();
    bool do_timedlock(const xtime& xt);
    void do_unlock();
    void do_lock(cv_state& state);
    void do_unlock(cv_state& state);

#if defined(BOOST_HAS_WINTHREADS)
    void* m_mutex;
    unsigned long m_count;
#elif defined(BOOST_HAS_PTHREADS)
    pthread_mutex_t m_mutex;
    pthread_cond_t m_unlocked;
    pthread_t m_thread_id;
    bool m_valid_id;
    unsigned m_count;
#elif defined(BOOST_HAS_MPTASKS)
    threads::mac::detail::scoped_critical_region m_mutex;
    threads::mac::detail::scoped_critical_region m_mutex_mutex;
    std::size_t m_count;
#endif
};

} // namespace boost

#endif // BOOST_RECURSIVE_MUTEX_WEK070601_HPP

// Change Log:
//    8 Feb 01  WEKEMPF Initial version.
//    1 Jun 01  WEKEMPF Modified to use xtime for time outs.  Factored out
//                      to three classes, mutex, try_mutex and timed_mutex.
//   11 Jun 01  WEKEMPF Modified to use PTHREAD_MUTEX_RECURSIVE if available.
//    3 Jan 03  WEKEMPF Modified for DLL implementation.
