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

#ifndef BOOST_TSS_WEK070601_HPP
#define BOOST_TSS_WEK070601_HPP

#include <boost/config.hpp>
// insist on threading support being available:
#include <boost/config/requires_threads.hpp>

#include <boost/utility.hpp>
#include <boost/thread/detail/config.hpp>

#if defined(BOOST_HAS_PTHREADS)
#   include <pthread.h>
#elif defined(BOOST_HAS_MPTASKS)
#   include <Multiprocessing.h>
#endif

namespace boost {

namespace detail {
class BOOST_THREAD_DECL tss : private noncopyable
{
public:
    tss(void (*cleanup)(void*)=0);
    ~tss();

    void* get() const;
    bool set(void* value);

private:
#if defined(BOOST_HAS_WINTHREADS)
    unsigned long m_key;
    void (*m_cleanup)(void*);
#elif defined(BOOST_HAS_PTHREADS)
    pthread_key_t m_key;
#elif defined(BOOST_HAS_MPTASKS)
    TaskStorageIndex m_key;
    void (*m_cleanup)(void*);
#endif
};

#if defined(BOOST_HAS_MPTASKS)
void thread_cleanup();
#endif
}

template <typename T>
class thread_specific_ptr : private noncopyable
{
public:
    thread_specific_ptr() : m_tss(&thread_specific_ptr<T>::cleanup) { }

    T* get() const { return static_cast<T*>(m_tss.get()); }
    T* operator->() const { return get(); }
    T& operator*() const { return *get(); }
    T* release() { T* temp = get(); m_tss.set(0); return temp; }
    void reset(T* p=0)
    {
        T* cur = get();
        if (cur == p) return;
        delete cur;
        m_tss.set(p);
    }

private:
    static void cleanup(void* p) { delete static_cast<T*>(p); }

    mutable detail::tss m_tss;
};

} // namespace boost

// Change Log:
//   6 Jun 01  WEKEMPF Initial version.

#endif // BOOST_TSS_WEK070601_HPP

