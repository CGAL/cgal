#ifndef OTHER_DWA20020601_HPP
# define OTHER_DWA20020601_HPP

# include <boost/python/detail/prefix.hpp>
// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.

# if _MSC_VER+0 >= 1020
#  pragma once
# endif

# include <boost/config.hpp>

namespace boost { namespace python {

template<class T> struct other
{ 
    typedef T type;
};

# ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
namespace detail
{
  template<typename T>
  class is_other
  {
   public:
      BOOST_STATIC_CONSTANT(bool, value = false); 
  };

  template<typename T>
  class is_other<other<T> >
  {
   public:
      BOOST_STATIC_CONSTANT(bool, value = true);
  };

  template<typename T>
  class unwrap_other
  {
   public:
      typedef T type;
  };

  template<typename T>
  class unwrap_other<other<T> >
  {
   public:
      typedef T type;
  };
}
# else // no partial specialization

}} // namespace boost::python

#include <boost/type.hpp>

namespace boost { namespace python {

namespace detail
{
  typedef char (&yes_other_t)[1];
  typedef char (&no_other_t)[2];
      
  no_other_t is_other_test(...);

  template<typename T>
  yes_other_t is_other_test(type< other<T> >);

  template<bool wrapped>
  struct other_unwrapper
  {
      template <class T>
      struct apply
      {
          typedef T type;
      };
  };

  template<>
  struct other_unwrapper<true>
  {
      template <class T>
      struct apply
      {
          typedef typename T::type type;
      };
  };

  template<typename T>
  class is_other
  {
   public:
      BOOST_STATIC_CONSTANT(
          bool, value = (
              sizeof(detail::is_other_test(type<T>()))
              == sizeof(detail::yes_other_t)));
  };

  template <typename T>
  class unwrap_other
      : public detail::other_unwrapper<
      is_other<T>::value
  >::template apply<T>
  {};
}

# endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

}} // namespace boost::python

#endif // #ifndef OTHER_DWA20020601_HPP
