#ifndef BOOST_DETAIL_LWM_GCC_HPP_INCLUDED
#define BOOST_DETAIL_LWM_GCC_HPP_INCLUDED

//
//  boost/detail/lwm_gcc.hpp
//
//  lightweight_mutex for GNU libstdc++ v3
//
//  http://gcc.gnu.org/onlinedocs/porting/Thread-safety.html
//
//  Copyright (c) 2002 Peter Dimov and Multi Media Ltd.
//  Copyright (c) 2002 Lars Gullik Bjønnes <larsbj@lyx.org>
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include <bits/atomicity.h>
#include <sched.h>

namespace boost
{

namespace detail
{

class lightweight_mutex
{
private:

    _Atomic_word a_;

    lightweight_mutex(lightweight_mutex const &);
    lightweight_mutex & operator=(lightweight_mutex const &);

public:

    lightweight_mutex(): a_(0)
    {
    }

    class scoped_lock;
    friend class scoped_lock;

    class scoped_lock
    {
    private:

        lightweight_mutex & m_;

        scoped_lock(scoped_lock const &);
        scoped_lock & operator=(scoped_lock const &);

    public:

        explicit scoped_lock(lightweight_mutex & m): m_(m)
        {
            while( __exchange_and_add(&m_.a_, 1) )
            {
                __atomic_add(&m_.a_, -1);
                sched_yield();
            }
        }

        ~scoped_lock()
        {
            __atomic_add(&m_.a_, -1);
        }
    };
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_LWM_GCC_HPP_INCLUDED
