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

#include <boost/thread/detail/config.hpp>

#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/thread/exceptions.hpp>

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
    tss(boost::function1<void, void*>* pcleanup) {
        if (pcleanup == 0) throw boost::thread_resource_error();
        try
        {
            init(pcleanup);
        }
        catch (...)
        {
            delete pcleanup;
            throw boost::thread_resource_error();
        }
    }

    void* get() const;
    void set(void* value);
    void cleanup(void* p);

private:
    unsigned int m_slot; //This is a "pseudo-slot", not a native slot

    void init(boost::function1<void, void*>* pcleanup);
};

#if defined(BOOST_HAS_MPTASKS)
void thread_cleanup();
#endif

template <typename T>
struct tss_adapter
{
    template <typename F>
    tss_adapter(const F& cleanup) : m_cleanup(cleanup) { }
    void operator()(void* p) { m_cleanup(static_cast<T*>(p)); }
    boost::function1<void, T*> m_cleanup;
};

} // namespace detail

template <typename T>
class thread_specific_ptr : private noncopyable
{
public:
    thread_specific_ptr()
        : m_tss(new boost::function1<void, void*>(
                    boost::detail::tss_adapter<T>(
                        &thread_specific_ptr<T>::cleanup)))
    {
    }
    thread_specific_ptr(void (*clean)(T*))
        : m_tss(new boost::function1<void, void*>(
                    boost::detail::tss_adapter<T>(clean)))
    {
    }
    ~thread_specific_ptr() { reset(); }

    T* get() const { return static_cast<T*>(m_tss.get()); }
    T* operator->() const { return get(); }
    T& operator*() const { return *get(); }
    T* release() { T* temp = get(); if (temp) m_tss.set(0); return temp; }
    void reset(T* p=0)
    {
        T* cur = get();
        if (cur == p) return;
        m_tss.set(p);
        if (cur) m_tss.cleanup(cur);
    }

private:
    static void cleanup(T* p) { delete p; }
    detail::tss m_tss;
};

} // namespace boost

#endif //BOOST_TSS_WEK070601_HPP

// Change Log:
//   6 Jun 01  
//      WEKEMPF Initial version.
//  30 May 02  WEKEMPF 
//      Added interface to set specific cleanup handlers.
//      Removed TLS slot limits from most implementations.
//  22 Mar 04 GlassfordM for WEKEMPF
//      Fixed: thread_specific_ptr::reset() doesn't check error returned
//          by tss::set(); tss::set() now throws if it fails.
//      Fixed: calling thread_specific_ptr::reset() or 
//          thread_specific_ptr::release() causes double-delete: once on
//          reset()/release() and once on ~thread_specific_ptr().
