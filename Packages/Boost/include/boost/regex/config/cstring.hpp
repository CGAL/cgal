/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */
 
 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         boost/regex/config/cstring.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: regex narrow character string fixes.
  */

#ifndef BOOST_REGEX_CONFIG_CSTRING_HPP
#define BOOST_REGEX_CONFIG_CSTRING_HPP

#include <cstring>
#include <cctype>
#ifndef __sgi
#ifdef __KCC
#include <ios>
#endif
#include <boost/config.hpp>

namespace std{

#ifdef __BORLANDC__
#pragma option push -w-8008 -w-8066 -w-8004
#endif

#ifdef BOOST_NO_STDC_NAMESPACE

// Any function that is a macro is converted into an inline function:
#ifdef memcmp
inline int boost_memcmp(const void * p1, const void * p2, size_t s)
{ return memcmp(p1, p2, s); }
#undef memcmp
inline int memcmp(const void * p1, const void * p2, size_t s)
{ return boost_memcmp(p1, p2, s); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::memcmp;
#endif

#ifdef memcpy
inline void *boost_memcpy(void * p1, const void *p2, size_t s)
{ return memcpy(p1, p2, s); }
#undef memcpy
inline void *memcpy(void * p1, const void *p2, size_t s)
{ return boost_memcpy(p1, p2, s); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::memcpy;
#endif

#ifdef memmove
inline void *(memmove)(void *, const void *, size_t)
{ return memmove(p1,p2,s); }
#undef memmove
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::memmove;
#endif

#ifdef memset
inline void *(boost_memset)(void *p, int a, size_t b)
{ return memset(p,a,b); }
#undef memset
inline void *(memset)(void *p, int a, size_t b)
{ return boost_memset(p,a,b); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::memset;
#endif

#ifdef strcat
inline char *(boost_strcat)(char *p1, const char *p2)
{ return strcat(p1,p2); }
#undef strcat
inline char *(strcat)(char *p1, const char *p2)
{ return boost_strcat(p1,p2); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::strcat;
#endif

#ifdef strcmp
inline int (boost_strcmp)(const char *p1, const char *p2)
{ return strcmp(p1,p2); }
#undef strcmp
inline int (strcmp)(const char *p1, const char *p2)
{ return boost_strcmp(p1,p2); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::strcmp;
#endif

#ifdef strcoll
inline int (boost_strcoll) (const char *p1, const char *p2)
{ return strcoll(p1,p2); }
#undef strcoll
inline int (strcoll) (const char *p1, const char *p2)
{ return boost_strcoll(p1,p2); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::strcoll;
#endif

#ifdef strcpy
inline char *(boost_strcpy)(char *p1, const char *p2)
{ return strcpy(p1,p2); }
#undef strcpy
inline char *(strcpy)(char *p1, const char *p2)
{ return boost_strcpy(p1,p2); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::strcpy;
#endif

#ifdef strlen
inline size_t (boost_strlen)(const char *p)
{ return strlen(p); }
#undef strlen
inline size_t (strlen)(const char *p)
{ return boost_strlen(p); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::strlen;
#endif

#ifdef strxfrm
inline size_t (boost_strxfrm)(char *p1, const char *p2, size_t s)
{ return strxfrm(p1,p2,s); }
#undef strxfrm
inline size_t (strxfrm)(char *p1, const char *p2, size_t s)
{ return boost_strxfrm(p1,p2,s); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::strxfrm;
#endif

#ifdef isalnum
inline int (boost_isalnum)(int i)
{ return isalnum(i); }
#undef isalnum
inline int (isalnum)(int i)
{ return boost_isalnum(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isalnum;
#endif

#ifdef isalpha
inline int (boost_isalpha)(int i)
{ return isalpha(i); }
#undef isalpha
inline int (isalpha)(int i)
{ return boost_isalpha(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isalpha;
#endif

#ifdef iscntrl
inline int (boost_iscntrl)(int i)
{ return iscntrl(i); }
#undef iscntrl
inline int (iscntrl)(int i)
{ return boost_iscntrl(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::iscntrl;
#endif

#ifdef isdigit
inline int (boost_isdigit)(int i)
{ return isdigit(i); }
#undef isdigit
inline int (isdigit)(int i)
{ return boost_isdigit(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isdigit;
#endif

#ifdef isgraph
inline int (boost_isgraph)(int i)
{ return isgraph(i); }
#undef isgraph
inline int (isgraph)(int i)
{ return boost_isgraph(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isgraph;
#endif

#ifdef islower
inline int (boost_islower)(int i)
{ return islower(i); }
#undef islower
inline int (islower)(int i)
{ return boost_islower(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::islower;
#endif

#ifdef isprint
inline int (boost_isprint)(int i)
{ return isprint(i); }
#undef isprint
inline int (isprint)(int i)
{ return boost_isprint(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isprint;
#endif

#ifdef ispunct
inline int (boost_ispunct)(int i)
{ return ispunct(i); }
#undef ispunct
inline int (ispunct)(int i)
{ return boost_ispunct(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::ispunct;
#endif

#ifdef isspace
inline int (isspace)(int i)
{ return isspace(i); }
#undef isspace
inline int (boost_isspace)(int i)
{ return boost_isspace(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isspace;
#endif

#ifdef isupper
inline int (isupper)(int i)
{ return isupper(i); }
#undef isupper
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isupper;
#endif

#ifdef isxdigit
inline int (isxdigit)(int i)
{ return isxdigit(i); }
#undef isxdigit
inline int (boost_isxdigit)(int i)
{ return boost_isxdigit(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::isxdigit;
#endif

#ifdef tolower
inline int (boost_tolower)(int i)
{ return tolower(i); }
#undef tolower
inline int (tolower)(int i)
{ return boost_tolower(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::tolower;
#endif

#ifdef toupper
inline int (boost_toupper)(int i)
{ return toupper(i); }
#undef toupper
inline int (toupper)(int i)
{ return boost_toupper(i); }
#elif defined(BOOST_NO_STDC_NAMESPACE)
using ::toupper;
#endif

#else

#undef memcmp
#undef memcpy
#undef memmove
#undef memset
#undef strcat
#undef strcmp
#undef strcoll
#undef strcpy
#undef strlen
#undef strxfrm
#undef isalnum
#undef isalpha
#undef iscntrl
#undef isdigit
#undef isgraph
#undef islower
#undef isprint
#undef ispunct
#undef isspace
#undef isupper
#undef isxdigit
#undef tolower
#undef toupper

#endif


#ifdef __BORLANDC__
#pragma option pop
#endif

} // namespace std

#endif // __sgi

#endif

