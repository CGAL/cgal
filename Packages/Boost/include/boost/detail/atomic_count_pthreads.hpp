#ifndef BOOST_DETAIL_ATOMIC_COUNT_PTHREADS_HPP_INCLUDED
#define BOOST_DETAIL_ATOMIC_COUNT_PTHREADS_HPP_INCLUDED

//
//  boost/detail/atomic_count_pthreads.hpp
//
//  Copyright (c) 2001, 2002 Peter Dimov and Multi Media Ltd.
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//

#include <pthread.h>

//
//  The generic pthread_mutex-based implementation sometimes leads to
//    inefficiencies. Example: a class with two atomic_count members
//    can get away with a single mutex.
//
//  Users can detect this situation by checking BOOST_AC_USE_PTHREADS.
//

namespace boost
{

namespace detail
{

class atomic_count
{
private:

    class scoped_lock
    {
    public:

        scoped_lock(pthread_mutex_t & m): m_(m)
        {
            pthread_mutex_lock(&m_);
        }

        ~scoped_lock()
        {
            pthread_mutex_unlock(&m_);
        }

    private:

        pthread_mutex_t & m_;
    };

public:

    explicit atomic_count(long v): value_(v)
    {
        pthread_mutex_init(&mutex_, 0);
    }

    ~atomic_count()
    {
        pthread_mutex_destroy(&mutex_);
    }

    void operator++()
    {
        scoped_lock lock(mutex_);
        ++value_;
    }

    long operator--()
    {
        scoped_lock lock(mutex_);
        return --value_;
    }

    operator long() const
    {
        scoped_lock lock(mutex_);
        return value_;
    }

private:

    atomic_count(atomic_count const &);
    atomic_count & operator=(atomic_count const &);

    mutable pthread_mutex_t mutex_;
    long value_;
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_ATOMIC_COUNT_PTHREADS_HPP_INCLUDED
