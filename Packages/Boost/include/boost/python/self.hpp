// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef SELF_DWA2002531_HPP
# define SELF_DWA2002531_HPP

# include <boost/python/detail/prefix.hpp>

namespace boost { namespace python {

#define BOOST_PYTHON_SELF_IS_CLASS

// Sink self_t into its own namespace so that we have a safe place to
// put the completely general operator templates which operate on
// it. It is possible to avoid this, but it turns out to be much more
// complicated and finally GCC 2.95.2 chokes on it.
namespace self_ns
{
# ifndef BOOST_PYTHON_SELF_IS_CLASS
  enum self_t { self };
# else 
  struct self_t {};
  extern BOOST_PYTHON_DECL self_t self;
# endif
}

using self_ns::self_t;
using self_ns::self;

}} // namespace boost::python

#endif // SELF_DWA2002531_HPP
