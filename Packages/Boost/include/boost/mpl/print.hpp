// Copyright David Abrahams 2003. Use, modification and distribution is
// subject to the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef PRINT_DWA20031231_HPP
# define PRINT_DWA20031231_HPP

#include <boost/mpl/identity.hpp>

namespace boost { namespace mpl
{
  
namespace aux
{
#if BOOST_MSVC
# pragma warning(push, 3)
// we only want one warning from MSVC, so turn off the other one
# pragma warning(disable: 4307)
#elif __MWERKS__
# pragma warn_hidevirtual on
   struct print_base { virtual void f() {} };
#endif

#if __EDG_VERSION__
  template <class T>
  struct dependent_unsigned
  {
      static const unsigned value = 1;
  };
#endif
}

  
template <class T>
struct print
  :  mpl::identity<T>
#if __MWERKS__
    , aux::print_base
#endif 
{
#if BOOST_MSVC
    enum { n = sizeof(T) + -1 };
#elif __MWERKS__
    void f(int);
#else 
    enum {
        n =
# if __EDG_VERSION__ 
        aux::dependent_unsigned<T>::value
# else 
           sizeof(T)
# endif 
           > -1, };
#endif 
};

}} // namespace boost::mpl

#if BOOST_MSVC
# pragma warning(pop)
#elif __MWERKS__
# pragma warn_hidevirtual reset
#endif

#endif // PRINT_DWA20031231_HPP
