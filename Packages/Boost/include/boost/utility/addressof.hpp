// Copyright (C) 2002 Brad King (brad.king@kitware.com) 
//                    Doug Gregor (gregod@cs.rpi.edu)
//                    Peter Dimov
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.

// For more information, see http://www.boost.org

#ifndef BOOST_UTILITY_ADDRESSOF_HPP
# define BOOST_UTILITY_ADDRESSOF_HPP

# include <boost/config.hpp>
# include <boost/detail/workaround.hpp>
# if BOOST_WORKAROUND(BOOST_MSVC, == 1300)
#  include <boost/type_traits/add_pointer.hpp>
# endif

namespace boost {

// Do not make addressof() inline. Breaks MSVC 7. (Peter Dimov)

// VC7 strips const from nested classes unless we add indirection here
# if BOOST_WORKAROUND(BOOST_MSVC, == 1300)
template <typename T> typename add_pointer<T>::type
# else
template <typename T> T*
# endif
addressof(T& v)
{
  return reinterpret_cast<T*>(
       &const_cast<char&>(reinterpret_cast<const volatile char &>(v)));
}

}

#endif // BOOST_UTILITY_ADDRESSOF_HPP
