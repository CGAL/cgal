#ifndef BOOST_DETAIL_LWM_WIN32_HPP_INCLUDED
#define BOOST_DETAIL_LWM_WIN32_HPP_INCLUDED

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  boost/detail/lwm_win32.hpp
//
//  Copyright (c) 2002, 2003 Peter Dimov
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//

#ifdef BOOST_USE_WINDOWS_H
#  include <windows.h>
#endif

#ifdef __BORLANDC__
# pragma warn -8027     // Functions containing while are not expanded inline
#endif

namespace boost
{

namespace detail
{

#ifndef BOOST_USE_WINDOWS_H

#ifdef _WIN64

// Intel 6.0 on Win64 version, posted by Tim Fenders to [boost-users]

extern "C" long_type __cdecl _InterlockedExchange(long volatile *, long);

#pragma intrinsic(_InterlockedExchange)

inline long InterlockedExchange(long volatile* lp, long l)
{
    return _InterlockedExchange(lp, l);
}

#else  // _WIN64

extern "C" __declspec(dllimport) long __stdcall InterlockedExchange(long volatile *, long);

#endif // _WIN64

extern "C" __declspec(dllimport) void __stdcall Sleep(unsigned long);

#endif // #ifndef BOOST_USE_WINDOWS_H

class lightweight_mutex
{
private:

    long l_;

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
            while( InterlockedExchange(&m_.l_, 1) )
            {
                // Note: changed to Sleep(1) from Sleep(0).
                // According to MSDN, Sleep(0) will never yield
                // to a lower-priority thread, whereas Sleep(1)
                // will. Performance seems not to be affected.

                Sleep(1);
            }
        }

        ~scoped_lock()
        {
            InterlockedExchange(&m_.l_, 0);

            // Note: adding a yield here will make
            // the spinlock more fair and will increase the overall
            // performance of some applications substantially in
            // high contention situations, but will penalize the
            // low contention / single thread case up to 5x
        }
    };
};

} // namespace detail

} // namespace boost

#ifdef __BORLANDC__
# pragma warn .8027     // Functions containing while are not expanded inline
#endif

#endif // #ifndef BOOST_DETAIL_LWM_WIN32_HPP_INCLUDED
