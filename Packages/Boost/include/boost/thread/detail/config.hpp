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

#ifndef BOOST_THREAD_CONFIG_WEK01032003_HPP
#define BOOST_THREAD_CONFIG_WEK01032003_HPP

#if defined(BOOST_HAS_WINTHREADS)
#   if defined(BOOST_THREAD_BUILD_DLL)
#       define BOOST_THREAD_DECL __declspec(dllexport)
#   else
#       define BOOST_THREAD_DECL __declspec(dllimport)
#   endif
#else
#   define BOOST_THREAD_DECL
#endif // BOOST_THREAD_SHARED_LIB

//
// Automatically link to the correct build variant where possible. 
// 
#if !defined(BOOST_ALL_NO_LIB) && !defined(BOOST_THREAD_NO_LIB) && !defined(BOOST_THREAD_BUILD_DLL)
//
// Set the name of our library, this will get undef'ed by auto_link.hpp
// once it's done with it:
//
#define BOOST_LIB_NAME boost_thread
//
// If we're importing code from a dll, then tell auto_link.hpp about it:
//
#  define BOOST_DYN_LINK
//
// And include the header that does the work:
//
#include <boost/config/auto_link.hpp>
#endif  // auto-linking disabled

#endif // BOOST_THREAD_CONFIG_WEK1032003_HPP
