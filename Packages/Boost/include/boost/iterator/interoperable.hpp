// (C) Copyright David Abrahams 2002.
// (C) Copyright Jeremy Siek    2002.
// (C) Copyright Thomas Witt    2002.
// Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef BOOST_INTEROPERABLE_23022003THW_HPP
# define BOOST_INTEROPERABLE_23022003THW_HPP

# include <boost/mpl/bool.hpp>
# include <boost/mpl/or.hpp>

# include <boost/type_traits/is_convertible.hpp>

# include <boost/iterator/detail/config_def.hpp> // must appear last

namespace boost
{

  //
  // Meta function that determines whether two
  // iterator types are considered interoperable.
  //
  // Two iterator types A,B are considered interoperable if either
  // A is convertible to B or vice versa.
  // This interoperability definition is in sync with the
  // standards requirements on constant/mutable container
  // iterators (23.1 [lib.container.requirements]).
  //
  // For compilers that don't support is_convertible 
  // is_interoperable gives false positives. See comments
  // on operator implementation for consequences.
  //
  template <typename A, typename B>
  struct is_interoperable
# ifdef BOOST_NO_STRICT_ITERATOR_INTEROPERABILITY
    : mpl::true_
# else
    : mpl::or_<
          is_convertible< A, B >
        , is_convertible< B, A > >
# endif
  { 
  };

} // namespace boost

# include <boost/iterator/detail/config_undef.hpp>

#endif // BOOST_INTEROPERABLE_23022003THW_HPP
