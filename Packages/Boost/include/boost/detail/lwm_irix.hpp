#ifndef BOOST_DETAIL_LWM_IRIX_HPP_INCLUDED
#define BOOST_DETAIL_LWM_IRIX_HPP_INCLUDED

//
//  boost/detail/lwm_irix.hpp
//
//  Copyright (c) 2002 Peter Dimov and Multi Media Ltd.
//  Copyright (c) 2002 Dan Gohman
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <sgidefs.h>
#include <mutex.h>
#include <sched.h>

namespace boost
{

namespace detail
{

class lightweight_mutex
{
private:

    __uint32_t l_;

    lightweight_mutex(lightweight_mutex const &);
    lightweight_mutex & operator=(lightweight_mutex const &);

public:

    lightweight_mutex(): l_(0)
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
            while( test_and_set32(&m_.l_, 1) )
            {
                sched_yield();
            }
        }

        ~scoped_lock()
        {
            m_.l_ = 0;
        }
    };
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_LWM_IRIX_HPP_INCLUDED
