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

#include <boost/config.hpp>

// insist on threading support being available:
#include <boost/config/requires_threads.hpp>

#if defined(BOOST_HAS_WINTHREADS)
#   if defined(BOOST_THREAD_BUILD_DLL) //Build dll
#       define BOOST_THREAD_DECL __declspec(dllexport)
#   elif defined(BOOST_THREAD_BUILD_LIB) //Build lib
#       define BOOST_THREAD_DECL
#   elif defined(BOOST_THREAD_USE_DLL) //Use dll
#       define BOOST_THREAD_DECL __declspec(dllimport)
#       define BOOST_DYN_LINK
#   elif defined(BOOST_THREAD_USE_LIB) //Use lib
#       define BOOST_THREAD_DECL
#   else //Use default
#       if defined(BOOST_MSVC) || defined(BOOST_INTEL_WIN)
            //For VC++, choose according to threading library setting
#           if defined(_DLL)
                //Threading library is dll: use Boost.Threads dll
#               define BOOST_THREAD_USE_DLL
#               define BOOST_THREAD_DECL __declspec(dllimport)
#               define BOOST_DYN_LINK
#           else
                //Threading library is lib: used Boost.Threads lib
#               define BOOST_THREAD_USE_LIB
#               define BOOST_THREAD_DECL
#           endif
#       else
            //For compilers not yet supporting auto-tss cleanup
            //with Boost.Threads lib, use Boost.Threads dll
#           define BOOST_THREAD_USE_DLL
#           define BOOST_THREAD_DECL __declspec(dllimport)
#           define BOOST_DYN_LINK
#       endif
#   endif
#else
#   define BOOST_THREAD_DECL
#   if defined(BOOST_THREAD_USE_LIB) //Use dll
#       define BOOST_THREAD_USE_DLL
#       define BOOST_DYN_LINK
#   elif defined(BOOST_THREAD_USE_DLL) //Use lib
#       define BOOST_THREAD_USE_LIB
#   else //Use default
#       define BOOST_THREAD_USE_LIB
#   endif
#endif // BOOST_HAS_WINTHREADS

//
// Automatically link to the correct build variant where possible.
//
#if !defined(BOOST_ALL_NO_LIB) && !defined(BOOST_THREAD_NO_LIB) && !defined(BOOST_THREAD_BUILD_DLL) && !defined(BOOST_THREAD_BUILD_LIB)
//
// Set the name of our library, this will get undef'ed by auto_link.hpp
// once it's done with it:
//
#if defined(BOOST_THREAD_LIB_NAME)
#    define BOOST_LIB_NAME BOOST_THREAD_LIB_NAME
#else
#    define BOOST_LIB_NAME boost_thread
#endif
//
// If we're importing code from a dll, then tell auto_link.hpp about it:
//
// And include the header that does the work:
//
#include <boost/config/auto_link.hpp>
#endif  // auto-linking disabled

#endif // BOOST_THREAD_CONFIG_WEK1032003_HPP
