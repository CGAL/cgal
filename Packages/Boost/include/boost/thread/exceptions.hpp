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

#ifndef BOOST_THREAD_EXCEPTIONS_PDM070801_H
#define BOOST_THREAD_EXCEPTIONS_PDM070801_H

#include <boost/config.hpp>
#include <boost/thread/detail/config.hpp>

//  pdm: Sorry, but this class is used all over the place & I end up
//       with recursive headers if I don't separate it
//  wek: Not sure why recursive headers would cause compilation problems
//       given the include guards, but regardless it makes sense to
//       seperate this out any way.

#include <stdexcept>

namespace boost {

class BOOST_THREAD_DECL lock_error : public std::logic_error
{
public:
    lock_error();
};

class BOOST_THREAD_DECL thread_resource_error : public std::runtime_error
{
public:
    thread_resource_error();
};

} // namespace boost

// Change log:
//    3 Jan 03  WEKEMPF Modified for DLL implementation.

#endif // BOOST_THREAD_CONFIG_PDM070801_H
