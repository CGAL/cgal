// Copyright David Abrahams 2004. Use, modification and distribution is
// subject to the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef IS_INCREMENTABLE_DWA200415_HPP
# define IS_INCREMENTABLE_DWA200415_HPP

# include <boost/type_traits/remove_cv.hpp>
# include <boost/mpl/bool.hpp>
# include <boost/detail/workaround.hpp>

namespace boost { namespace detail { 

// is_incrementable<T> metafunction
//
// Requires: Given x of type T&, if the expression ++x is well-formed
// it must have complete type; otherwise, it must neither be ambiguous
// nor violate access.

// This namespace ensures that ADL doesn't mess things up.
namespace is_incrementable_
{
  struct tag {};

  // any soaks up implicit conversions and makes the following
  // operator++ less-preferred than any other such operator which
  // might be found via ADL.
  struct any { template <class T> any(T const&); };
  tag operator++(any const&);

  // two check overloads help us identify which operator++ was picked
  char (& check(tag) )[2];
  
  template <class T>
  char check(T const&);
  

  template <class T>
  struct
# if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
  impl
# else 
  is_incrementable
# endif 
  {
      static typename remove_cv<T>::type& x;

      BOOST_STATIC_CONSTANT(
          bool
        , value = sizeof(is_incrementable_::check(++x)) == 1
      );

      typedef mpl::bool_<(
# if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
                             ::boost::detail::is_incrementable_::is_incrementable<T>::
# endif 
                             value)> type;
  };
}

# if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
template <class T>
struct is_incrementable : is_incrementable_::impl<T>
{
};
# else
using is_incrementable_::is_incrementable;
# endif 

}} // namespace boost::detail

#endif // IS_INCREMENTABLE_DWA200415_HPP
