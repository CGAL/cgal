/* boost integer_traits.hpp header file
 *
 * Copyright Jens Maurer 2000
 * Permission to use, copy, modify, sell, and distribute this software
 * is hereby granted without fee provided that the above copyright notice
 * appears in all copies and that both that copyright notice and this
 * permission notice appear in supporting documentation,
 *
 * Jens Maurer makes no representations about the suitability of this
 * software for any purpose. It is provided "as is" without express or
 * implied warranty.
 *
 * $Id$
 *
 * Idea by Beman Dawes, Ed Brey, Steve Cleary, and Nathan Myers
 */

//  See http://www.boost.org/libs/integer for documentation.


#ifndef BOOST_INTEGER_TRAITS_HPP
#define BOOST_INTEGER_TRAITS_HPP

#include <boost/config.hpp>
#include <boost/limits.hpp>

// These are an implementation detail and not part of the interface
#include <limits.h>
#if !defined(BOOST_NO_INTRINSIC_WCHAR_T) && !defined(BOOST_NO_CWCHAR)
#include <wchar.h>
#endif


namespace boost {
template<class T>
class integer_traits : public std::numeric_limits<T>
{
public:
  BOOST_STATIC_CONSTANT(bool, is_integral = false);
};

namespace detail {
template<class T, T min_val, T max_val>
class integer_traits_base
{
public:
  BOOST_STATIC_CONSTANT(bool, is_integral = true);
  BOOST_STATIC_CONSTANT(T, const_min = min_val);
  BOOST_STATIC_CONSTANT(T, const_max = max_val);
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
//  A definition is required even for integral static constants
template<class T, T min_val, T max_val>
const bool integer_traits_base<T, min_val, max_val>::is_integral;

template<class T, T min_val, T max_val>
const T integer_traits_base<T, min_val, max_val>::const_min;

template<class T, T min_val, T max_val>
const T integer_traits_base<T, min_val, max_val>::const_max;
#endif

} // namespace detail

template<>
class integer_traits<bool>
  : public std::numeric_limits<bool>,
    public detail::integer_traits_base<bool, false, true>
{ };

template<>
class integer_traits<char>
  : public std::numeric_limits<char>,
    public detail::integer_traits_base<char, CHAR_MIN, CHAR_MAX>
{ };

template<>
class integer_traits<signed char>
  : public std::numeric_limits<signed char>,
    public detail::integer_traits_base<signed char, SCHAR_MIN, SCHAR_MAX>
{ };

template<>
class integer_traits<unsigned char>
  : public std::numeric_limits<unsigned char>,
    public detail::integer_traits_base<unsigned char, 0, UCHAR_MAX>
{ };

#ifndef BOOST_NO_INTRINSIC_WCHAR_T
template<>
class integer_traits<wchar_t>
  : public std::numeric_limits<wchar_t>,
#if defined(WCHAR_MIN) && defined(WCHAR_MAX)
    public detail::integer_traits_base<wchar_t, WCHAR_MIN, WCHAR_MAX>
#elif defined(__BORLANDC__) || defined(__CYGWIN__) || defined(__MINGW32__) || (defined(__BEOS__) && defined(__GNUC__))
    // No WCHAR_MIN and WCHAR_MAX, whar_t is short and unsigned:
    public detail::integer_traits_base<wchar_t, 0, 0xffff>
#elif (defined(__sgi) && (!defined(__SGI_STL_PORT) || __SGI_STL_PORT < 0x400))\
    || (defined __APPLE__)\
    || (defined(__OpenBSD__) && defined(__GNUC__))\
    || (defined(__NetBSD__) && defined(__GNUC__))\
    || (defined(__FreeBSD__) && defined(__GNUC__))\
    || (defined(__hpux) && defined(__GNUC__) && (__GNUC__ == 3) && !defined(__SGI_STL_PORT))
    // No WCHAR_MIN and WCHAR_MAX, wchar_t has the same range as int.
    //  - SGI MIPSpro with native library
    //  - gcc 3.x on HP-UX
    //  - Mac OS X with native library
    //  - gcc on FreeBSD, OpenBSD and NetBSD
    public detail::integer_traits_base<wchar_t, INT_MIN, INT_MAX>
#elif defined(__hpux) && defined(__GNUC__) && (__GNUC__ == 2) && !defined(__SGI_STL_PORT)
    // No WCHAR_MIN and WCHAR_MAX, wchar_t has the same range as unsigned int.
    //  - gcc 2.95.x on HP-UX
    // (also, std::numeric_limits<wchar_t> appears to return the wrong values).
    public detail::integer_traits_base<wchar_t, 0, UINT_MAX>
#else
#error No WCHAR_MIN and WCHAR_MAX present, please adjust integer_traits<> for your compiler.
#endif
{ };
#endif // BOOST_NO_INTRINSIC_WCHAR_T

template<>
class integer_traits<short>
  : public std::numeric_limits<short>,
    public detail::integer_traits_base<short, SHRT_MIN, SHRT_MAX>
{ };

template<>
class integer_traits<unsigned short>
  : public std::numeric_limits<unsigned short>,
    public detail::integer_traits_base<unsigned short, 0, USHRT_MAX>
{ };

template<>
class integer_traits<int>
  : public std::numeric_limits<int>,
    public detail::integer_traits_base<int, INT_MIN, INT_MAX>
{ };

template<>
class integer_traits<unsigned int>
  : public std::numeric_limits<unsigned int>,
    public detail::integer_traits_base<unsigned int, 0, UINT_MAX>
{ };

template<>
class integer_traits<long>
  : public std::numeric_limits<long>,
    public detail::integer_traits_base<long, LONG_MIN, LONG_MAX>
{ };

template<>
class integer_traits<unsigned long>
  : public std::numeric_limits<unsigned long>,
    public detail::integer_traits_base<unsigned long, 0, ULONG_MAX>
{ };

#if !defined(BOOST_NO_INTEGRAL_INT64_T) && !defined(BOOST_NO_INT64_T)
#if defined(ULLONG_MAX) && defined(BOOST_HAS_LONG_LONG)

template<>
class integer_traits<long long>
  : public std::numeric_limits<long long>,
    public detail::integer_traits_base<long long, LLONG_MIN, LLONG_MAX>
{ };

template<>
class integer_traits<unsigned long long>
  : public std::numeric_limits<unsigned long long>,
    public detail::integer_traits_base<unsigned long long, 0, ULLONG_MAX>
{ };

#elif defined(ULONG_LONG_MAX) && defined(BOOST_HAS_LONG_LONG)

template<>
class integer_traits<long long>  : public std::numeric_limits<long long>,    public detail::integer_traits_base<long long, LONG_LONG_MIN, LONG_LONG_MAX>{ };
template<>
class integer_traits<unsigned long long>
  : public std::numeric_limits<unsigned long long>,
    public detail::integer_traits_base<unsigned long long, 0, ULONG_LONG_MAX>
{ };

#elif defined(ULONGLONG_MAX) && defined(BOOST_HAS_LONG_LONG)

template<>
class integer_traits<long long>
  : public std::numeric_limits<long long>,
    public detail::integer_traits_base<long long, LONGLONG_MIN, LONGLONG_MAX>
{ };

template<>
class integer_traits<unsigned long long>
  : public std::numeric_limits<unsigned long long>,
    public detail::integer_traits_base<unsigned long long, 0, ULONGLONG_MAX>
{ };

#elif defined(_LLONG_MAX) && defined(_C2) && defined(BOOST_HAS_LONG_LONG)

template<>
class integer_traits<long long>
  : public std::numeric_limits<long long>,
    public detail::integer_traits_base<long long, -_LLONG_MAX - _C2, _LLONG_MAX>
{ };

template<>
class integer_traits<unsigned long long>
  : public std::numeric_limits<unsigned long long>,
    public detail::integer_traits_base<unsigned long long, 0, _ULLONG_MAX>
{ };

#endif
#endif

} // namespace boost

#endif /* BOOST_INTEGER_TRAITS_HPP */



