// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef REGISTERED_DWA2002710_HPP
# define REGISTERED_DWA2002710_HPP
# include <boost/python/type_id.hpp>
# include <boost/python/converter/registry.hpp>
# include <boost/python/converter/registrations.hpp>
# include <boost/type_traits/transform_traits.hpp>
# include <boost/type_traits/cv_traits.hpp>
# include <boost/detail/workaround.hpp>

namespace boost { namespace python { namespace converter { 

struct registration;

namespace detail
{
  template <class T>
  struct registered_base
  {
      static registration const& converters;
  };
}

template <class T>
struct registered
    : detail::registered_base<
        typename add_reference<
           typename add_cv<T>::type
        >::type
    >
{
};

# if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
    && !BOOST_WORKAROUND(BOOST_MSVC, BOOST_TESTED_AT(1310))
// collapses a few more types to the same static instance.  MSVC7.1
// fails to strip cv-qualification from array types in typeid.  For
// some reason we can't use this collapse there or array converters
// will not be found.
template <class T>
struct registered<T&>
  : registered<T> {};
# endif

//
// implementations
//
namespace detail
{
  template <class T>
  registration const& registered_base<T>::converters
     = registry::lookup(type_id<T>());
}
}}} // namespace boost::python::converter

#endif // REGISTERED_DWA2002710_HPP
