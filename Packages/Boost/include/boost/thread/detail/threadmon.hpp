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

#include <boost/config.hpp>
// insist on threading support being available:
#include <boost/config/requires_threads.hpp>

#include <boost/thread/detail/config.hpp>

#ifdef BOOST_HAS_WINTHREADS

#include <boost/thread/detail/config.hpp>

extern "C" BOOST_THREAD_DECL int on_thread_exit(void (__cdecl * func)(void));

#endif // BOOST_HAS_WINTHREADS
