#ifndef BOOST_DETAIL_ATOMIC_COUNT_WIN32_HPP_INCLUDED
#define BOOST_DETAIL_ATOMIC_COUNT_WIN32_HPP_INCLUDED

// MS compatible compilers support #pragma once

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

//
//  boost/detail/atomic_count_win32.hpp
//
//  Copyright (c) 2001, 2002, 2003 Peter Dimov
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//

#ifdef BOOST_USE_WINDOWS_H
#  include <windows.h>
#endif

namespace boost
{

namespace detail
{

#ifndef BOOST_USE_WINDOWS_H

#ifdef _WIN64

// Intel 6.0 on Win64 version, posted by Tim Fenders to [boost-users]

extern "C" long_type __cdecl _InterlockedIncrement(long volatile *);
extern "C" long_type __cdecl _InterlockedDecrement(long volatile *);

#pragma intrinsic(_InterlockedIncrement)
#pragma intrinsic(_InterlockedDecrement)

inline long InterlockedIncrement(long volatile * lp)
{ 
    return _InterlockedIncrement(lp);
}

inline long InterlockedDecrement(long volatile* lp)
{ 
    return _InterlockedDecrement(lp);
}

#else  // _WIN64

extern "C" __declspec(dllimport) long __stdcall InterlockedIncrement(long volatile *);
extern "C" __declspec(dllimport) long __stdcall InterlockedDecrement(long volatile *);

#endif // _WIN64

#endif // #ifndef BOOST_USE_WINDOWS_H

class atomic_count
{
public:

    explicit atomic_count(long v): value_(v)
    {
    }

    long operator++()
    {
        // Some older <windows.h> versions do not accept volatile
        return InterlockedIncrement(const_cast<long*>(&value_));
    }

    long operator--()
    {
        return InterlockedDecrement(const_cast<long*>(&value_));
    }

    operator long() const
    {
        return value_;
    }

private:

    atomic_count(atomic_count const &);
    atomic_count & operator=(atomic_count const &);

    volatile long value_;
};

} // namespace detail

} // namespace boost

#endif // #ifndef BOOST_DETAIL_ATOMIC_COUNT_WIN32_HPP_INCLUDED
