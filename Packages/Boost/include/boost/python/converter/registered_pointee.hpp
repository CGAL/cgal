// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef REGISTERED_POINTEE_DWA2002710_HPP
# define REGISTERED_POINTEE_DWA2002710_HPP
# include <boost/python/converter/registered.hpp>
# include <boost/python/converter/pointer_type_id.hpp>
# include <boost/python/converter/registry.hpp>
# include <boost/type_traits/transform_traits.hpp>
# include <boost/type_traits/cv_traits.hpp>

namespace boost { namespace python { namespace converter { 

struct registration;

# ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
template <class T>
struct registered_pointee
    : registered<
        typename remove_pointer<
           typename remove_cv<
              typename remove_reference<T>::type
           >::type
        >::type
    >
{
};
# else
namespace detail
{
  template <class T>
  struct registered_pointee_base
  {
      static registration const& converters;
  };
}

template <class T>
struct registered_pointee
    : detail::registered_pointee_base<
        typename add_reference<
           typename add_cv<T>::type
        >::type
    >
{
};

//
// implementations
//
namespace detail
{
  template <class T>
  registration const& registered_pointee_base<T>::converters
     = registry::lookup(pointer_type_id<T>());
}

# endif 
}}} // namespace boost::python::converter

#endif // REGISTERED_POINTEE_DWA2002710_HPP
