#ifndef BOOST_DETAIL_LWM_LINUX_HPP_INCLUDED
#define BOOST_DETAIL_LWM_LINUX_HPP_INCLUDED

//
//  boost/detail/lwm_linux.hpp
//
//  Copyright (c) 2002 Peter Dimov and Multi Media Ltd.
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//

//
//  This implementation uses <asm/atomic.h>. This is a kernel header;
//  using kernel headers in a user program may cause a number of problems,
//  and not all flavors of Linux provide the atomic instructions.
//
//  This file is only provided because the performance of this implementation
//  is about 3.5 times higher than the pthreads version. Use at your own risk
//  (by defining BOOST_USE_ASM_ATOMIC_H.)
//

#include <asm/atomic.h>
#include <sched.h>

namespace boost
{

namespace detail
{

class lightweight_mutex
{
private:

    atomic_t a_;

    lightweight_mutex(lightweight_mutex const &);
    lightweight_mutex & operator=(lightweight_mutex const &);

public:

    lightweight_mutex()
    {
        atomic_t a = ATOMIC_INIT(1);
        a_ = a;
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
            while( !atomic_dec_and_test(&m_.a_) )
            {
                atomic_inc(&m_.a_);
                sched_yield();
            }
        }

        ~scoped_lock()
        {
            atomic_inc(&m_.a_);
        }
    };
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_LWM_LINUX_HPP_INCLUDED
