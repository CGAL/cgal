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

#ifndef BOOST_ONCE_WEK080101_HPP
#define BOOST_ONCE_WEK080101_HPP

#include <boost/thread/detail/config.hpp>

#if defined(BOOST_HAS_PTHREADS)
#   include <pthread.h>
#endif

namespace boost {

#if defined(BOOST_HAS_PTHREADS)

typedef pthread_once_t once_flag;
#define BOOST_ONCE_INIT PTHREAD_ONCE_INIT

#elif (defined(BOOST_HAS_WINTHREADS) || defined(BOOST_HAS_MPTASKS))

typedef long once_flag;
#define BOOST_ONCE_INIT 0

#endif

void BOOST_THREAD_DECL call_once(void (*func)(), once_flag& flag);

} // namespace boost

// Change Log:
//   1 Aug 01  WEKEMPF Initial version.

#endif // BOOST_ONCE_WEK080101_HPP
