//  (C) Copyright John Maddock 2001. 
//  (C) Copyright Jens Maurer 2001. 
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for most recent version.

//  config for libstdc++ v3
//  not much to go in here:

#ifdef __GLIBCXX__
#define BOOST_STDLIB "GNU libstdc++ version " BOOST_STRINGIZE(__GLIBCXX__)
#else
#define BOOST_STDLIB "GNU libstdc++ version " BOOST_STRINGIZE(__GLIBCPP__)
#endif

#if !defined(_GLIBCPP_USE_WCHAR_T) && !defined(_GLIBCXX_USE_WCHAR_T)
#  define BOOST_NO_CWCHAR
#  define BOOST_NO_CWCTYPE
#  define BOOST_NO_STD_WSTRING
#  define BOOST_NO_STD_WSTREAMBUF
#endif

#if defined(__osf__) && !defined(_REENTRANT) && defined(_GLIBCXX_HAVE_GTHR_DEFAULT)
// GCC 3.4 on Tru64 forces the definition of _REENTRANT when any std lib header
// file is included, therefore for consistency we define it here as well.
#  define _REENTRANT
#endif

#ifdef __GLIBCXX__ // gcc 3.4 and greater:
#  ifdef _GLIBCXX_HAVE_GTHR_DEFAULT
      // 
      // If the std lib has thread support turned on, then turn it on in Boost
      // as well.  We do this because some gcc-3.4 std lib headers define _REENTANT
      // while others do not...
      // 
#     define BOOST_HAS_THREADS
#  else
#     define BOOST_DISABLE_THREADS
#  endif
#endif

 
#if !defined(_GLIBCPP_USE_LONG_LONG) \
    && !defined(_GLIBCXX_USE_LONG_LONG)\
    && defined(BOOST_HAS_LONG_LONG)
// May have been set by compiler/*.hpp, but "long long" without library
// support is useless.
#  undef BOOST_HAS_LONG_LONG
#endif
