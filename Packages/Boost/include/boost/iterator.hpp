//  interator.hpp workarounds for non-conforming standard libraries  ---------//

//  (C) Copyright Boost.org 2000. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

//  See http://www.boost.org/libs/utility for documentation.

//  Revision History
//  12 Jan 01 added <cstddef> for std::ptrdiff_t (Jens Maurer)
//  28 Jun 00 Workarounds to deal with known MSVC bugs (David Abrahams)
//  26 Jun 00 Initial version (Jeremy Siek)

#ifndef BOOST_ITERATOR_HPP
#define BOOST_ITERATOR_HPP

#include <iterator>
#include <cstddef>           // std::ptrdiff_t
#include <boost/config.hpp>

namespace boost
{
# if defined(BOOST_NO_STD_ITERATOR) && !defined(BOOST_MSVC_STD_ITERATOR)
  template <class Category, class T,
    class Distance = std::ptrdiff_t,
    class Pointer = T*, class Reference = T&>
  struct iterator
  {
    typedef T         value_type;
    typedef Distance  difference_type;
    typedef Pointer   pointer;
    typedef Reference reference;
    typedef Category  iterator_category;
  };
# else

  // declare iterator_base in namespace detail to work around MSVC bugs which
  // prevent derivation from an identically-named class in a different namespace.
  namespace detail {
   template <class Category, class T, class Distance, class Pointer, class Reference>
#  if !defined(BOOST_MSVC_STD_ITERATOR)
   struct iterator_base : std::iterator<Category, T, Distance, Pointer, Reference> {};
#  else
   struct iterator_base : std::iterator<Category, T, Distance>
   {
     typedef Reference reference;
     typedef Pointer pointer;
     typedef Distance difference_type;
   };
#  endif
  }

  template <class Category, class T, class Distance = std::ptrdiff_t,
            class Pointer = T*, class Reference = T&>
   struct iterator : detail::iterator_base<Category, T, Distance, Pointer, Reference> {};
# endif
} // namespace boost

#endif // BOOST_ITERATOR_HPP
