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
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
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

    lightweight_mutex(): a_(1)
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
            while( !__exchange_and_add(&m_.a_, -1) )
            {
                __atomic_add(&m_.a_, 1);
                sched_yield();
            }
        }

        ~scoped_lock()
        {
            __atomic_add(&m_.a_, 1);
        }
    };
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_LWM_GCC_HPP_INCLUDED
