// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef BASE_TYPE_TRAITS_DWA2002614_HPP
# define BASE_TYPE_TRAITS_DWA2002614_HPP

# include <boost/python/detail/prefix.hpp>

namespace boost { namespace python { 

namespace detail
{
  struct unspecialized {};
}

// Derive from unspecialized so we can detect whether traits are
// specialized
template <class T> struct base_type_traits
  : detail::unspecialized
{};

template <>
struct base_type_traits<PyObject>
{
    typedef PyObject type;
};

template <>
struct base_type_traits<PyTypeObject>
{
    typedef PyObject type;
};

}} // namespace boost::python

#endif // BASE_TYPE_TRAITS_DWA2002614_HPP
